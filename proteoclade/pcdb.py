'''ProteoClade database (PCDB) creation module

Functions (for user)
---------
merge_fastas
create_pcdb

Functions for other submodules
---------
_hashit
'''
import sqlite3
import multiprocessing as mp
import csv
import os
import time
import datetime
import hashlib
import re
import string
import random

from .pcutilities import db_fetch_params, _get_todays_date, _db_connect, _enforce_worker_count
from .pcconstants import cleave_rules


def merge_fastas(merged_fasta_name, fasta_directory = 'fastas'):
    '''
    For merging fastas together in a directory.
    Used in preparation of a targeted database search, i.e. MaxQuant/Mascot.

    Parameters
    ----------
    merged_fasta_name: string
        Name of .fasta file that will result from merging other fastas
    fasta_directory: string
        Directory from which to read fastas (default 'fastas')

    Notes
    -----
    Creates a .fasta file containing all read fasta entries.
    '''
    
    def _check_header(header_line, incorrect_formats):
        '''
        Handles what happens if a bad header is found.
        These are the absolute minimum requirements.

        Parameters
        ----------
        header_line: string
            FASTA header line to check
        incorrect_formats: integer
            Number of times an incorrect format has been seen

        Returns
        -------
        incorrect_formats: integer
            Number of times an incorrect format has been seen
        '''       
        if not all([
            line.count('|') >= 2, #Normally would be == 2, but a few genes have | in their name
            line.count('OX=') == 1]):

            if incorrect_formats < 1:
                print('Warning: Incorrect format')
                print(line)
                print('These will be merged, but will not be added to PCDB.')
                input('Press the enter key to keep merging.')
                
            incorrect_formats += 1
            
        return incorrect_formats
            
    
    if not merged_fasta_name.endswith('.fasta'):
        merged_fasta_name += '.fasta'

    incorrect_formats = 0
    entry_count = 0
    with open(merged_fasta_name, 'w') as merged_file:
        
        for file in os.listdir(fasta_directory):
            if file.endswith('.fasta'):
                with open(os.path.join(fasta_directory, file)) as input_file:
                    for line in input_file:

                        if line.startswith('>'):
                            if entry_count and entry_count % 10000 == 0:
                                print(f'Merged {entry_count} entries.')
                            entry_count += 1
                            incorrect_formats = _check_header(line, incorrect_formats)
                        merged_file.write(line.strip() + '\n')
    
    print(f'Merging done. {incorrect_formats} incorrect formats were found.')

def create_pcdb(database_name, fasta_directory = 'fastas',
                min_length = 7, max_length = 55, missed_cleavages = 2,
                m_cleave = True, li_swap = True, rule = 'trypsin/p',
                temp_directory = None, worker_count = None, reverse = False):
    '''Creates the PCDB file which stores in silico digested peptides, genes, and organism info.

    Parameters
    ----------
        database_name: string
                Name of pcdb file; should be descriptive of what it contains
        fasta_directory: string
                Directory of FASTAs to use as input (default ‘fastas’)
        min_length: integer
                Minimum peptide amino acid count to include in database (default 7)
        max_length: integer
                Maximum peptide amino acid count to include in database (default 55)
        missed_cleavages: integer
                Number of times a protease is allowed to miss a cut site. (default 2)
        m_cleave: bool
                Whether or not N-terminal methionines are cleaved from proteins (default True)
        li_swap: bool
                Whether peptides stored will have leucines converted to isoleucines (default True)
        rule: string or tuple
                Protease rule for cutting sites (default ‘trypsin/p’).                

                if string: must be an enzyme option available in ProteoClade. See Appendix.

                if tuple: must be tuple of strings, (“regex_sites”,”terminus”) ex. (r”[RK]”, “C”).
                Use tuple for custom enzyme rules.
        temp_directory: None or string
                Directory for temporary database operations if space is a concern (default None)

                if None: uses working directory
        worker_count: None or integer
                Number of worker processes to use. Only set to experiment with performance. (default None)
                if None: determines processes up to a maximum of 6 to use. More processes does not necessarily increase performance.
        reverse: bool
            Whether to reverse protein sequences prior to digestion and storage. Used for FDR mitigation.
            
    Examples
    --------
        >>> create_pcdb("human_mouse_ref.pcdb") #creates a trypsin PCDB file using default settings
        >>> create_pcdb("bacteria_swissprot.pcdb", rule = "asp-n") #creates an AspN PCDB file
            
    Notes
    -----
        Creates a .pcdb SQLite file in the working directory.
'''

    fasta_files = os.listdir(fasta_directory)
    assert any([x.endswith('.fasta') for x in fasta_files]), \
       f"Your fasta folder {fasta_directory} contains no '.fasta' files"


    worker_count = _enforce_worker_count(worker_count) #guarantees 1 <= workers <= 6 for performance
    print("Using ", worker_count, " worker(s).")

    if temp_directory is None: #For indexing purposes - important if HD is full
        temp_directory = os.getcwd()

    if not database_name.endswith('.pcdb'): #for convention, and so people don't forget what the file is
        database_name += '.pcdb'
        
    if os.path.exists(database_name):
        raise Exception('Database already exists. Aborting.')


    db_time = _get_todays_date()
    db_params = (min_length, max_length, missed_cleavages, m_cleave, li_swap, rule, db_time)

    production_queue = mp.Queue(maxsize = 30000) #Must be <= ~32k for MacOS
    writing_queue = mp.Queue(maxsize = 30000)

    produce_p = mp.Process(target = _producer, args = (fasta_directory, production_queue), daemon = True)
    produce_p.start()

    workers = []
    for i in range(worker_count):
        worker_p = mp.Process(target = _worker, args = (production_queue, writing_queue, db_params, reverse), daemon = True)
        worker_p.start()
        workers.append(worker_p)

    db_writable_params = ['|'.join(x) if type(x) in (tuple, list) else x for x in db_params] #For inserting into dbcustom cleave rule
    writing_p = mp.Process(target = _writer, args = (writing_queue, database_name, worker_count, db_writable_params, temp_directory), daemon = True)
    writing_p.start()

    produce_p.join()
    print("Production done. (Task 1/3)")
    
    for i in range(worker_count):
        production_queue.put(None)
    for w in workers:
        w.join()
    print("Workers have consumed tasks. (Task 2/3)")

    writing_queue.put(None)
    writing_p.join()
    print("Writing process has finished. (Task 3/3)")
    print("All tasks have finished.")


##########Begin create_pcdb support functions##########
def _producer(fasta_directory, production_queue):
    '''Producer process to read fasta files and put organisms, genes, and sequences into header.

    Parameters
    ----------
    fasta_directory: string
        Folder to read fasta files from.
    production_queue: mp.Queue object
        Queue in which (gene, organism, sequence) are stored as a tuple.
    '''
    organisms_missed = 0 #perhaps remove later
    
    def fasta_entry_check(header, sequence):
        #Check organism
        if 'OX=' in header:
            organism = int(header.split('OX=')[1].split()[0])
        else:
            organisms_missed += 1
            return None #skip altogether if no organism info in line
            
        #Check gene; if no GN, use Uniprot Id
        if 'GN=' in header:
            gene = header.split('GN=')[1].split()[0]
        elif header.count('|') >= 2:
            #normally would be ==2, but there are some genes with | in them
            gene = header.split('|')[1]
        else:
            gene = "NOGENEIDFOUND"

        production_queue.put((gene, organism, sequence))
        
    #Parse fastas
    for file in os.listdir(fasta_directory):
        if file.endswith('.fasta'):
            with open(os.path.join(fasta_directory, file)) as input_file:
                print('Reading from file ', file)
                header = None
                seq = ''
                for line in input_file:
                    
                    if line.startswith('>') and seq:
                        fasta_entry_check(header, seq)
                        header = line.strip()
                        seq = ''
                    elif line.startswith('>'): #first entry of file
                        header = line.strip()
                    else:
                        seq += line.strip()

                fasta_entry_check(header, seq) #For last entry
                
                                            
def _worker(input_queue, output_queue, db_params, reverse):
    '''Worker process(es) that handle digestion and hashing of peptide sequences.

    Parameters
    ----------
    input_queue: mp.Queue object
        Retrieves entries from producer process.
    output_queue: mp.Queue object
        Holds (protein_seq, gene, organism, digest_results) for writer process.
    db_params: tuple
        (min_length, max_length, missed_cleavages, m_cleave, li_swap, rule, date) to keep track of db parameters
    reverse: bool
        Whether the protein sequences are to be reversed before digestion and storage
    '''
    *params, rule, date = db_params
    rule_to_use = _cleave_rule_determination(rule)
##    db_params = (min_length, max_length, missed_cleavages, m_cleave, li_swap, rule, db_time)
    
    while True:
        item = input_queue.get()
        if item is None:
            break

        gene, organism, protein_seq = item
        if reverse:
            protein_seq = protein_seq[::-1]
        
        digest_results = _digest(protein_seq, *params, rule_to_use, reverse)
        output_queue.put((protein_seq, gene, organism, digest_results))

def _writer(writing_queue, database, worker_count, db_params, temp_directory):
    '''Writer process for pcdb creation. Handles SQLite database setup and insertion/indexing of entries.

    Parameters
    ----------
    writing_queue: mp.Queue object
        Holds (protein_seq, gene, organism, digest_results) for inserting into database.
    database: string
        Name of database to connect to.
    worker_count: integer
        Number of workers to use for multithreaded SQLite indexing.
    db_params: tuple
        (min_length, max_length, missed_cleavages, m_cleave, li_swap, rule, date) to keep track of db parameters
    temp_directory: string
        Folder location to use for SQLite temporary use, specifically during indexing.
    '''
    #DB setup
    c, conn = _db_connect(database) #DB to write to

    pragmas = ['synchronous = OFF',
               'journal_mode = OFF',
               'locking_mode = EXCLUSIVE',
               'temp_store = 1']

    for pragma in pragmas:
        c.execute(f'PRAGMA {pragma}')
    
    c.execute("CREATE TABLE IF NOT EXISTS MainData (Protein TEXT, Organism INTEGER, Gene TEXT)")
    c.execute("CREATE TABLE IF NOT EXISTS Reference (HashPeptide INTEGER, ProtRowID INTEGER)")
    c.execute("BEGIN EXCLUSIVE TRANSACTION")

    #Queue processing and DB insertion
    db_start_time = time.time()
    counter = 0  
    while True:
        if counter % 10000 == 0:
            print(counter, " proteins have been digested and inserted into the database.")
        item = writing_queue.get() #Item is (PROTSEQ, Gene, organism, (hashresults))
        if item is None:
            break

        protseq, gene, organism, hashresults = item
        c.execute("INSERT INTO MainData (Protein, Organism, Gene) VALUES (?,?,?)", (protseq, organism, gene))
        lastid = c.lastrowid
        
        hashids = [(i, lastid) for i in hashresults]
        c.executemany("INSERT INTO Reference (HashPeptide, ProtRowID) VALUES (?,?)", hashids)        
        counter += 1

    insert_complete_time = time.time() - db_start_time
    print(f"Insertion took {insert_complete_time} seconds.")
    print("Done inserting entries, now indexing...")

    c.execute(f'PRAGMA threads = {worker_count}')
    c.execute(f'PRAGMA temp_store_directory = "{temp_directory}"') #Otherwise it ends up on default storage drive, usually C

    indexing_start_time = time.time()
    c.execute("CREATE INDEX HashIndex ON Reference (HashPeptide)")
    indexing_complete_time = time.time() - indexing_start_time
    print(f"Indexing took {indexing_complete_time} seconds.")

    #SQL formatting for DBParams
    db_param_info = (("min_length","INTEGER"),
                 ("max_length","INTEGER"),
                 ("missed_cleavages","INTEGER"),
                 ("m_cleave","BOOL"),
                 ("li_swap","BOOL"),
                 ("rule","TEXT"),
                 ("date","TEXT"))
    db_params_names = ', '.join([x[0] for x in db_param_info])
    db_params_insert = ', '.join([" ".join(x) for x in db_param_info])
    db_values_insert = ', '.join(['?' for x in db_param_info])
    
    c.execute(f"CREATE TABLE IF NOT EXISTS DBParams ({db_params_insert})")
    c.execute(f"INSERT INTO DBParams ({db_params_names}) VALUES ({db_values_insert})", db_params) 
    
    conn.commit()
    conn.close()

def _digest(sequence, min_length, max_length, missed_cleavages, m_cleave, li_swap, rule_to_use, reverse):
    '''Take a protein sequence and cut it into pieces, then hash for database insertion.

    Parameters
    ----------
    sequence: string
        Protein sequence in all capital letters to chop up.
    min_length: integer
        Minimum amino acid count of peptides to keep for database.
    max_length: integer
        Maximum amino acid count of peptides to keep for database.
    missed_cleavages: integer
        Number of times a protease is allowed to miss a specific site.
    m_cleave: bool
        Whether or not protein N-terminal methionines are removed.
    li_swap: bool
        Whether peptides will be stored with all leucines converted to isoleucines.
    rule_to_use: tuple
        Tuple of strings with tuple[0] being a regex expression for amino acid specificity and tuple[1] as either 'n' or 'c' cut direction.
    reverse: bool
        Whether protein sequence is to be reversed before storage
        
    Returns
    -------
    cut_set: set
        Set of integers (hashed peptides) that result from the cut rules used.
    '''
    site_specificity, cut_terminus = rule_to_use

    if reverse:
        if m_cleave and sequence[-1] == 'M':
            sequence = sequence[:-1]
    else:
        if m_cleave and sequence[0] == 'M':
            sequence = sequence[1:]

    #The following two blocks can be simplified to just use the re.match sites   

    cut_sites = []
    cut_peptides = []
    sites_matched = re.finditer(site_specificity, sequence)
    for site in sites_matched:
        cut_sites.append(site.start())        

    ###Find sites
    last_site = 0

    if cut_terminus.lower() == 'c':
        for i in cut_sites:
            cut_peptides.append(sequence[last_site:i+1])
            last_site = i + 1
        cut_peptides.append(sequence[last_site:])

    elif cut_terminus.lower() == 'n':
        for i in cut_sites:
            cut_peptides.append(sequence[last_site:i])
            last_site = i
        cut_peptides.append(sequence[last_site:])

    cut_and_missed = list(cut_peptides) #duplicate to add to for iteration
    
    missed_counter = 1
    while missed_counter <= missed_cleavages:
        
        missed_peptides = [
            ''.join(cut_peptides[i:i+1+missed_counter])
            for i in range(0,len(cut_peptides) - missed_counter)
            ] #subtract missed counter here to not duplicate c-terminal peps

        cut_and_missed += missed_peptides
        missed_counter += 1

#Diagnostic code for troubleshooting digest rules
##    with open(f'{random.randint(0,99999)}.txt','w') as w1:
##        for i in cut_and_missed:
##            if min_length <= len(i) <= max_length:
##                w1.write(i+ '\n')

    #Hash peptides
    if li_swap:
        cut_set = set([
            _hashit(i.replace("L","I"))
            for i in cut_and_missed
            if min_length <= len(i) <= max_length
            ])
    else:
        cut_set = set([
            _hashit(i)
            for i in cut_and_missed
            if min_length <= len(i) <= max_length
            ])

    return cut_set
##########End create_pcdb support functions##########

def _cleave_rule_determination(rule):     
    '''Handle whether the user chooses a built in digest rule or supplies their own.

    Parameters
    ----------
    rule: string or tuple
        if string: a built in rule, ex: "typsin/p"
        if tuple: a custom rule ("regexcutsites","terminus") ex: ("[RK]","c")

    Returns
    -------
    rules_to_use
        Tuple of strings containing (regex cutsites, terminus)
    '''
    if rule not in cleave_rules:
        assert type(rule) in (tuple, list), 'Cleave rule must be a tuple of strongs ("regexrule","terminus"), or be built in.'
        cutsites, direction = rule
        assert direction.lower() in ('c','n'), 'Second argument of cleave rule needs to be a valid protein terminus.'
        rules_to_use = rule #assume the user has put in a reg_ex string
    else:
        rules_to_use = cleave_rules.get(rule)

    return rules_to_use


def _hashit(sequence, digest_size = 8, order = "big"):
    '''A deterministic hash to convert peptide sequences to integers.

    Parameters
    ----------
    sequence: string
        Peptide sequence to be converted to an integer
    digest_size: integer
        Number of bytes (signed) to use for the returned integer (default 8)
    order: string
        Endianness of bytes. Recommended to leave at default. (default "big")

    Returns
    -------
    int_digested: integer
        Integer representing the peptide sequence converted by the hashing algorithm blake2b.
    '''
    hashed_value = hashlib.blake2b(bytes(sequence, encoding = "ASCII"), digest_size = digest_size)
    digested = hashed_value.digest()
    int_digested = int.from_bytes(digested, order, signed = True)
    return int_digested

if __name__ == "__main__":
    pass

    

    
    
    
    
