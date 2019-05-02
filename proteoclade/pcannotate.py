'''Annotates experimental files with taxa and genes.

Functions for user
---------
annotate_peptides
annotate_denovo
filter_taxa
'''

import sqlite3
from queue import Queue
from threading import Thread

import csv
csv.field_size_limit(1000000000)

from multiprocessing import cpu_count
import time #doesnt appear needed; try linter
import pickle
import os
import string

from .pcutilities import ncbi_check, _get_header_delim_reader, db_fetch_params, _db_connect, _enforce_worker_count
from .pcdb import _hashit, _digest
from .pcconstants import valid_quant_cols, valid_sample_cols, valid_score_cols, valid_pep_cols

def annotate_peptides(file, db, pctaxa, taxon_levels = None, worker_threads = None):
    '''Drives the taxonomic and gene annotation of peptide-containing files.

    Parameters
    ----------
    file: string
        csv or txt file containing wide-form, peptide entries
    db: string
        PCDB file containing digested peptides to match w/ experiment
    pctaxa: string
        PCTAXA file containing taxonomic mapping for species and above
    taxon_levels: None, string, or tuple
        Which taxa to annotate above the organism level (default None)
    worker_threads: None or integer
        Number of worker threads to use. (default None)

        if None: will use up to 6 threads.

    Notes
    -----
    Outputs csv or txt file with all data and appended taxonomic and gene annotations

    'annotated\_' + 'denovo_matched' + file
    '''  
    annotator = PeptideAnnotator(**locals())
    annotator.annotate()

def annotate_denovo(file, db, pctaxa, method = "dbconstrain",
                    taxon_levels = None, worker_threads = None):
    '''Drives the annotation of denovo/psm-containing files.

    Parameters
    ----------
    file: string
        csv or txt file containing long-form PSM entries
    db: string
        PCDB file containing digested peptides to match w/ experiment
    pctaxa: string
        PCTAXA file containing taxonomic mapping for species and above
    method: string
        "dbconstrain": serially checks PSM candidates against the PCDB

        "top": only looks at top scoring PSM candidate
        (default: "dbconstrain")
    taxon_levels: None, string, or tuple
        Which taxa to annotate above the organism level (default None)
    worker_threads: None or integer
        Number of worker threads to use. (default None)

        if None: will use up to 6 threads.

    Notes
    -----
    Output is csv or txt file with all data and appended taxonomic and gene annotations

    'denovo\_matched_\' + file

    'annotated\_' + 'denovo_matched' + file
    ''' 
    annotator = DeNovoAnnotator(**locals())
    annotator.annotate()


class AnnotationManager:
    '''Base class for annotations and de novo parsing.
    Only used for inheritance.
    '''
    
    def __init__(self, file, db, worker_threads = None, taxon_levels = None):

        #Common to targeted and denovo
        self.file_in = file
        assert os.path.exists(db), "DB needs to exist."
        self.db = db
        self.worker_threads = worker_threads
        self.reader_queue = Queue(30000) #Must be <= ~32k for MacOS
        self.writer_queue = Queue(30000)

        db_params = db_fetch_params(self.db)    
        if not db_params:
            raise Exception("Could not determine DB params. DB may be corrupted.")

        #Only need li_swap, but this is for remembering param order
        _minpep, _maxpep, _misscleave, _methcleave, self.li_swap, _rule, _date = db_params

        #Attributes specific to each *Annotator subclass
        self.type = None
        self.producer_args = None
        self.worker_args = None
        self.writer_args = None
        self.pctaxa = None #Dict containing taxon mapping
        self.taxon_levels = taxon_levels

        self.method = None #Denovo specific attribute

    def annotate(self):
        '''Function that manages the queues and threads.'''
        threads = _enforce_worker_count(self.worker_threads)

        #Start working
        p = Thread(target = self._produce, args = self.producer_args, daemon = True)
        p.start()

        workers = []
        for w in range(threads):
            t = Thread(target = self._worker, args = self.worker_args, daemon = True)
            workers.append(t)
            t.start()
        print(f'Started {threads} workers.')

        wri = Thread(target = self._writer, args = self.writer_args, daemon = True)
        wri.start()

        p.join()
        for w in workers:
            self.reader_queue.put(None)

        for w in workers:
            w.join()
        print("Annotation finished. Writing...")

        self.writer_queue.put(None)
        wri.join()

        print("Done.")

        if self.taxon_levels and self.type == "psm" and self.pctaxa:
            #At this point, the file has been exported with organisms and genes
            #Use annotate function to handle higher taxa
            PeptideAnnotator(f'denovo_matched_{os.path.basename(self.file_in)}', self.db, self.pctaxa, self.taxon_levels).annotate()


class PeptideAnnotator(AnnotationManager):
    '''Annotates peptide files. See annotate_peptide function docstrings.'''
    def __init__(self, file, db, pctaxa, taxon_levels = None, worker_threads = None):
        super().__init__(file, db, worker_threads)

        if taxon_levels:
            if type(taxon_levels) is str: #make this a tuple for iteration
                taxon_levels = (taxon_levels,)
            self.taxon_levels = ncbi_check(taxon_levels) #restrict to valid taxa
                
            #Load taxonomy dict into memory
            print('Loading taxa...')
            with open(pctaxa, 'rb') as tax_input:
                self.pc_taxdict = pickle.load(tax_input)
            print('Taxa loaded.')
        else:
            self.pc_taxdict = None
            
        self.type = "peptide"
        self.producer_args = (self.reader_queue, self.file_in)
        self.worker_args = (self.reader_queue, self.writer_queue, self.db, self.li_swap)
        self.writer_args = (self.writer_queue, self.file_in, self.pc_taxdict, self.taxon_levels)

    def _produce(self, q, file):
        '''Read from input file and put row/seq in queue.'''        
        pep_titles = valid_pep_cols
        with open(file) as input_file:
            header, _, file_reader = _get_header_delim_reader(file, input_file)
            
            possible_seqs = 0
            for pos, i in enumerate(header): #fix with .index method; or map(count())
                if i.lower() in pep_titles:
                    seq_position = pos
                    possible_seqs += 1
            assert possible_seqs == 1

            for row in file_reader:
                q.put((seq_position, row))

    def _worker(self, q_in, q_out, db, li_swap):
        '''Fetch information from PCDB'''
        c, conn = _db_connect(db)
        while True:
            item = q_in.get()
            
            if item is None:
                q_in.task_done()
                break
            
            seq_position, row = item
            pep_target = row[seq_position]
            
            #For strange formats:
            pep_target = ''.join(
                [x.upper() for x in pep_target if x in string.ascii_letters]
                )

            if li_swap:
                pep_target = pep_target.replace("L","I") # Check for LI here

            taxa_found = _fetch_taxa(c, pep_target, li_swap)
            if taxa_found:
                organisms, genes = taxa_found
            else:
                organisms, genes = ["NONE"], ["NONE"]

            #annotate even if nothing is found for those rows
            q_out.put((seq_position, row, organisms, genes))

            q_in.task_done()
            
    def _writer(self, q, file, pctaxa, taxa = None):
        '''Write annotated output.'''


        header, delim, _ = _get_header_delim_reader(file, None)

        if not all([x in header for x in ('organisms','genes')]):
            header.append('organisms')
            header.append('genes')
            include_organisms_genes = True
        else:
            include_organisms_genes = False                

        if taxa:
            for taxon in taxa:
                header.append(taxon)

        with open('annotated_' + file, 'w', newline = '') as output_file:

            output_writer = csv.writer(output_file, delimiter = delim)
            output_writer.writerow(header)

            counter = 1
            
            while True:
                if counter % 1000 == 0:
                    print(f'{counter} rows have been processed.')
                
                item = q.get()
                if item is None:
                    q.task_done()
                    break

                seq_position, row, organisms, genes = item
                organisms_write = '|'.join([str(x) for x in organisms])
                genes = '|'.join(genes)            

                if include_organisms_genes:
                    row += [organisms_write, genes]

                if taxa:
                    for taxon in taxa:
                        higher_taxa = []
                        for organism in organisms:
                            if organism in pctaxa:
                                higher_taxa.append(pctaxa[organism].get(taxon, f'{taxon}_NOT_FOUND'))
                            else:
                                higher_taxa.append('ORGANISM_NOT_FOUND')
                                
                        row.append('|'.join(higher_taxa))

                output_writer.writerow(row)
    ##        
    ##            try:
    ##                output_writer.writerow(row)
    ##            except Exception:
    ##                print(row)

                counter += 1
                q.task_done()

    
class DeNovoAnnotator(AnnotationManager):
    '''Annotates denovo/psm files. See annotate_denovo function docstrings.'''    
    def __init__(self, file, db, pctaxa, method = "dbconstrain",
                 taxon_levels = None, worker_threads = None):
        super().__init__(file, db, worker_threads, taxon_levels)
        
        self.method = method.lower()
        assert self.method in ("dbconstrain","top")
        
        #if pctaxa is made optional in the future
##        if self.taxon_levels:
##            assert self.pctaxa, 'Make sure to include the .pctaxa file.'

##        scol, qcol, pcol, sccol = self._file_check(self.file_in)
        cols = self._file_check(self.file_in)
        
        self.pctaxa = pctaxa

        self.type = "psm"
        self.producer_args = (self.reader_queue, self.file_in, *cols, self.li_swap)
        self.worker_args = (self.reader_queue, self.writer_queue, self.db,
                            self.li_swap, self.method)
        self.writer_args = (self.writer_queue, self.file_in)

    def _file_check(self, file):
        '''Tries to make sure file format is workable.'''       
        with open(file) as input_file:
            header, _, file_reader = _get_header_delim_reader(file, input_file)
            header_lower = [x.lower() for x in header]

            quant_possibles = 0
            quant_col = None
            for i in valid_quant_cols:
                quant_possibles += header_lower.count(i)
                if i in header_lower:
                    quant_col = header_lower.index(i)
            assert quant_possibles <= 1, \
            "File has more than one quant column."

            sample_possibles = 0
            for i in valid_sample_cols:
                sample_possibles += header_lower.count(i)
                if i in header_lower:
                    sample_col = header_lower.index(i)
            assert sample_possibles == 1, \
            f'File has {sample_possibles} sample columns and must be 1'

            pep_possibles = 0
            for i in valid_pep_cols:
                pep_possibles += header_lower.count(i)
                if i in header_lower:
                    pep_col = header_lower.index(i)
            assert pep_possibles == 1, \
            f'File has {pep_possibles} peptide columns and must be 1'

            #Not implemented
            score_col = None
##            score_possibles = 0
##            for i in valid_score_cols:
##                score_possibles += header_lower.count(i)
##                if i in header_lower:
##                    score_col = header_lower.index(i)
##            assert score_possibles <= 1, \
##            f'File has {score_possibles} score columns and must be <= 1'    

        #Index positions for columns. Score is not currently used    
        return sample_col, quant_col, pep_col, score_col

    def _produce(self, q, file, sample_col, quant_col, pep_col, score_col, li_swap):
        '''Read PSM file and put candidates in queue.'''
        #Assumed to be sorted by scan (any order) and then score (descending)
        #Queue will have all candidates for a single scan
        
        last_sample = None
        sample_holder = []

        #Debufstuff
        counter = 1

        with open(file) as input_file:
            header, _, file_reader = _get_header_delim_reader(file, input_file)
            
            for row in file_reader:
                sample = row[sample_col]

                counter += 1
                if counter % 1000 == 0:
                    print(f'{counter} rows have been produced.')

                if sample.count(':') != 1:
                    print("Skipping sample: ", sample, ". Needs to be formatted 'sample:scan'")
                    continue
                    
                if sample != last_sample and last_sample:
    ##                if score <= score_threshold: #Not implemented
                    q.put(sample_holder)
                    sample_holder = []
                    last_sample = sample
                elif sample != last_sample:
                    last_sample = sample

                if quant_col:
                    quant = float(row[quant_col])
                else:
                    quant = None

                peptide = row[pep_col]
                if li_swap:
                    peptide = peptide.replace("L","I")

                #Not implemented
                score = None
##                if score_col:
##                    score = int(row[score_col])
##                else:
##                    score = None
                    
                sample_holder.append((sample, quant, peptide, score))

            if sample_holder:
                q.put(sample_holder)

    def _worker(self, q_in, q_out, db, li_swap, method):
        '''Fetch taxonomic/gene information for candidates.'''
        c, conn = _db_connect(db)
        while True:
            item = q_in.get()

            if item is None:
                q_in.task_done()
                break

            if method == 'dbconstrain':
                for pos, data in enumerate(item):
                    taxa_found = _fetch_taxa(c, data[2], li_swap)
                    if taxa_found:
                        q_out.put((data, taxa_found, pos))
                        break #q.put stuff, then break
                    
            else: #method == "top"
                taxa_found = _fetch_taxa(c, item[0][2], li_swap)
                if taxa_found:
                    q_out.put((item[0], taxa_found, None)) #FIX THIS what is data referring to

            q_in.task_done()
  
    def _writer(self, q, file):
        '''Ouput annotations to file.
        Iterates through file twice to preserve memory with annotate_peptides.
        See base class for second pass.
        '''
        pep_dict = dict() #{Peptideseq : {sample: [summed_value, count]}}
        pep_taxa = dict() #{Peptideseq : [organisms, genes]}
        samples = []
        counter = 1
        
        while True:
            
##            if counter % 1000 == 0:
##                print(f'{counter} psms have been annotated.')
                      
            item = q.get()
            if item is None:
                q.task_done()
                break

            data, taxa_found, pos = item
            samplescan, quant, pepseq, score = data 
            if quant is None:
                quant = 0 #for addition  
            sample = samplescan.split(':')[0]

            if sample not in samples:
                samples.append(sample)

            #Add intensities if available; also add counts for PSMs    
            if pepseq not in pep_dict:
                pep_dict[pepseq] = {}
            if pepseq not in pep_taxa:
                pep_taxa[pepseq] = taxa_found
            if sample not in pep_dict[pepseq]:
                pep_dict[pepseq][sample] = [quant, 1]
            else:
                pep_dict[pepseq][sample][0] += quant
                pep_dict[pepseq][sample][1] += 1

            counter += 1
            q.task_done()
            
        ### Writing out data ###
        print("Writing out data.")
        
        file_out = f'denovo_matched_{os.path.basename(file)}'
        _, delim, _ = _get_header_delim_reader(file, None)
        
        header = ['peptide']
        header.extend(samples)
        header.extend([x + '_count' for x in samples])
        header.append('organisms')
        header.append('genes')

        with open(file_out, 'w', newline = '') as output_file:
            output_writer = csv.writer(output_file, delimiter = delim)
            output_writer.writerow(header)

            for peptide in pep_dict:
                sample_quants = []
                sample_counts = []
                
                for sample in samples:
                    samp_retrieve = pep_dict.get(peptide).get(sample)
                    if samp_retrieve:
                        sample_quants.append(samp_retrieve[0])
                        sample_counts.append(samp_retrieve[1])
                    else:
                        sample_quants.append(0)
                        sample_counts.append(0)

                organisms, genes = pep_taxa.get(peptide)
                organisms = [str(x) for x in organisms]

                row = [peptide]
                for i in (sample_quants, sample_counts):
                    row.extend(i)

                row.append('|'.join(organisms))
                row.append('|'.join(genes))

                output_writer.writerow(row)


def filter_taxa(file, taxon_levels, taxa, unique = False):
    '''Filters peptide files based on desired taxa.

    Parameters
    ----------
    file: string
        csv or txt file containing wide-form, peptide entries 
    taxon_levels: string, list, or tuple
        Taxonomic ranks to include in file search. Must be annotated
    taxa: string, list, or tuple
        Taxa to include in filter
    unique: bool
        Whether specified taxa must be unique in their given taxonomic rank

    Notes
    -----
    Output is csv or txt file pared down by filter specifications.

    'filtered\_' + file name
    '''
    
    if type(taxon_levels) not in (tuple, list):
        taxon_levels = [taxon_levels]
    if type(taxa) not in (tuple, list):
        taxa = [taxa]
    taxon_levels = [x.lower() for x in taxon_levels]
    taxa = [x.lower() for x in taxa]

    with open(file) as input_file:
        header, delim, reader = _get_header_delim_reader(file, input_file)
        lower_header = [x.lower() for x in header]
        
        for tax_level in taxon_levels:
            assert lower_header.count(tax_level) == 1, \
                   f'Taxon level {tax_level} was missing from your header.'

        output_file = open('filtered_' + file, 'w', newline = '')
        output_writer = csv.writer(output_file, delimiter = delim)
        output_writer.writerow(header)
    
        
        taxa_cols = [lower_header.index(x) for x in taxon_levels]

        for row in reader:
            keep_flag = False
            for col in taxa_cols:
                row_taxa = [x.lower() for x in row[col].split('|')]

                if unique is False and any([t in row_taxa for t in taxa]):
                    keep_flag = True
                elif (unique and any([t in row_taxa for t in taxa]) and
                      len(set(row_taxa)) == 1):
                    keep_flag = True
            if keep_flag:
                output_writer.writerow(row)
                

def _fetch_taxa(c, peptide, li_swap):
    '''Database query that finds organisms and genes from PCDB files.
    Unified to work for targeted and de novo searches.

    Parameters
    ----------
    c: sqlite3 cursor object
        sqlite3.connect().cursor()
    peptide: string
        Peptide sequence to search in database.
    li_swap: bool
        True/False for whether the database has converted all leucine to isoleucine

    Returns
    -------
    (organisms,genes): tuple
        Organisms as a list of integers and genes as a list of strings                        
    None: None
        No organisms found for peptide
    '''
    pep_target = ''.join(
        [x.upper() for x in peptide if x in string.ascii_letters]
        )
    organisms = []
    genes = []
    hashval = _hashit(pep_target)
    c.execute("SELECT ProtRowID FROM Reference WHERE HashPeptide = (?)", (hashval,))
    ref_results = c.fetchall()
    if ref_results:     
        for result in ref_results: #Check collisions
            search_result = result[0]
            c.execute("SELECT * FROM MainData WHERE rowID = (?)", (search_result,))
            main_results = c.fetchall()
            if main_results:
                if li_swap:
                    if pep_target.replace("L","I") in main_results[0][0].replace("L","I"):
                        genes.append(main_results[0][2])
                        organisms.append(main_results[0][1])
                else:
                     if pep_target in main_results[0][0]: 
                        genes.append(main_results[0][2])
                        organisms.append(main_results[0][1])                            

    if organisms and genes:
        return organisms, genes
    else:
        return None  


                
if __name__ == "__main__":
    pass #Reserved for later
