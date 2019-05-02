'''Utility functions that serve one of three purposes:
1. Interface with the web
2. Provide general utility to the user
3. Provide support to multiple submodules in ProteoClade

Functions for user
---------
download_uniprot
download_uniprot_batch
download_cRAP

download_taxonomy
load_taxonomy
ncbi_check

db_fetch_params
db_stats

Functions for other submodules
---------
ncbi_check
_get_header_delim_reader
_get_todays_date
_enforce_worker_count
_db_connect
'''

import os
from urllib import request
from urllib.parse import urlencode
from zipfile import ZipFile
import pickle
import sqlite3
import datetime
import csv
import multiprocessing as mp
import string

from .pcconstants import ncbi_ranks, uniprot_options


#---UNIPROT---
def download_uniprot(*targets, download_folder = 'fastas'):
    '''Download FASTA protein sequences from UniProt

    Parameters
    ----------
    targets: tuple
        One or more tuples containing an integer and a string (TaxonID, DBtype)
        TaxonID must be a NCBI-valid taxon identifier OR ‘all’ for all taxa (up to ~60GB)
        DBType must be one of: “s”, “sr”, “r”, “t”, “a”

        (SwissProt, SwissProt Reference, Reference, TrEMBL, or All, respectively)
    download_folder: string
        Folder which will contain downloaded FASTA files. Default: ‘fasta’ subdirectory
            
    Examples
    --------
    
    >>> download_uniprot((9606, ‘s’),(10090, ‘s’)) #downloads human and mouse SwissProt entries
    >>> download_uniprot((‘all’,’a’)) #downloads every entry in UniProt
    
    Notes
    -----
        Downloads UniProt-derived FASTA file(s) with specified parameters.
        Naming convention: taxonid_StartingEntryCount.fasta
        
    '''
    uniprot_restful_url = 'https://www.uniprot.org/uniprot/?'  
    uniprot_format = {"format":"FASTA"}

    for target in targets:
        assert type(target) in (tuple, list), f"Format for targets is container of containers: {target}"
    
    #Set up folder
    if not os.path.exists(download_folder):
        os.mkdir(download_folder)

    for target in targets:
        
        taxon, database_type = target #put this first in case of formatting issues

        #Check input to see if it is an integer
        if type(taxon) == str and taxon.lower() == 'all': #special case for uniprot; doesn't recognize root level = 1
            quotable_taxon = ''
        else:
            try:
                int(taxon)
                quotable_taxon = f'taxonomy:"{taxon}"'
            except ValueError:
                print("Not a valid taxon ID: ", taxon)
                continue

        #Modify query for rest API
        additional_mods = uniprot_options.get(database_type.upper())
        if additional_mods:
            additional_suffix = ' '.join([f'{key}:{val}' for key, val in additional_mods.items()])
            combined_string = ' '.join([quotable_taxon, additional_suffix])
        else:
            combined_string = quotable_taxon
        query = {"query":combined_string}

        #check header first for total number of values
        query_for_headers = uniprot_restful_url + urlencode(query)
        header_check = request.urlopen(query_for_headers)
        total_results = int(header_check.headers['X-Total-Results'])
        print(f'{total_results} total results for taxon {taxon}.')

        #Download in batches to prevent HTTP from dying on big files
        missed_downloads = []
        entries_per_batch = 1000000 #max entries to get in one http request
        for i in range(0, total_results, entries_per_batch):
            file_description = (f'Taxon: {taxon}, Starting Count: {i}')
            
            try:
##                query.update({'compress':'yes'}) Only if you want to unzip
                query.update({'offset':i, 'limit': entries_per_batch})
                query.update(uniprot_format)

                full_url = uniprot_restful_url + urlencode(query)
                target_file = os.path.join(download_folder, f'{taxon}_{i}.fasta')
                _download(full_url, target_file, file_description = file_description)
            except Exception: #null byte sent/server messing up ; change ValueError or http.client.IncompleteRead
                missed_downloads.append((full_url, target_file, file_description))

        while missed_downloads: #exhaustively fetch
            try:
                params = missed_downloads.pop()
                print("Attempting re-download: ", *params)
                _download(*params)
            except Exception: #Same, see as above
                missed_downloads.append(params)

    print("All downloads have finished.")

def download_uniprot_batch(file, download_folder = 'fastas'):
    '''Download UniProt entries from a tab delimited txt file

    Parameters
    ----------
    file: string
        tab-delimited txt file with 2 columns: [TaxonID,UniProtMods]

        Use the uniprot_options constant for available mods
    download_folder: string
        subdirectory to store downloads (default 'fastas')

    Notes
    -----
    Easier method to download larger numbers of taxa.
    '''
    with open(file) as input_file:
        for line in input_file:
            if len(line.split('\t')) == 2:  #"Must be a tab delimited file with 2 elements and no header: Taxon[tab]Reviewed."
                taxon, reviewed = line.strip().split('\t')
                download_uniprot((taxon, reviewed), download_folder = download_folder)
            else:
                print("Improper line found: ", line)

#---cRAP---
def download_cRAP(directory = 'fastas'):
    '''Downloads contaminant Repository for Affinity Purification,
    (cRAP) sequence database.

    Parameters
    ----------
    directory: string
        folder to download cRAP to (default 'fastas')

    Notes
    -----
    Supplies an edited fasta file from cRAP
    '''
    cRAP_url = "ftp://ftp.thegpm.org/fasta/crap/crap.fasta"
    
    if not os.path.exists(directory):
        os.mkdir(directory)

    cRAP_file = os.path.join(directory, 'cRAP.fasta')
    _download(cRAP_url, cRAP_file)
    print("cRAP downloaded.")
    _cRAP_edit(cRAP_file) #Produce edited file
    os.remove(cRAP_file)
    print("cRAP edited. Done.")

def _cRAP_edit(cRAP_file):
    """
    download_cRAP support function
    Editing the cRAP file to make sure it is obvious these are contaminants
    Handling the pipe symbol was necessary because it was omitted in some cRAP entries.
    """
    write_target = os.path.join(os.path.dirname(cRAP_file), "cRAP_edit.fasta")
    with open(cRAP_file) as file_in, open(write_target, 'w') as file_out:
        for line in file_in:
            if line.startswith('>'):
                if '|' in line:
                    uniprot_format = line.split('|')
                    gene = uniprot_format[1]
                    gene = "CONTAM_" + gene
                    file_out.write(f'{uniprot_format[0]}|{gene}| OS=CONTAM GN={gene} OX=-1\n')
                else: #One broken entry at end
                    file_out.write(f'>xx|CONTAM_{line[1:].strip()}| OS=CONTAM GN={line[1:].strip()} OX=-1\n')
            else:
                file_out.write(line)

#---NCBI---
def download_taxonomy(directory = 'taxonomy_downloads'):
    '''Download taxonomy mappings from NCBI.

    Parameters
    ----------
    directory: string
        Where to store temporary files downloaded from NCBI (default 'taxonomy_downloads')

    Notes
    -----
    Unzips taxonomy files from NCBI. 
    Will call _taxonomy_mapper to produce PCTAXA file in working directory. 
    PCTAXA file is a pickled dict of taxonomy mapping. 
    Naming will be Y-M-D formatted so you can remember when it was retrieved.
    '''
    taxonomy_url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip'
    
    if not os.path.exists(directory):
        os.mkdir(directory)

    taxonomy_file = os.path.join(directory, 'taxdmp.zip')
    _download(target = taxonomy_url, file_path = taxonomy_file, file_description = "NCBI Taxonomy")
    print("NCBI Download done.")
    _unzip(taxonomy_file, members = ("names.dmp", "nodes.dmp"))
    print("Unzipping done.")

    _taxonomy_mapper(directory_in = directory)


def _taxonomy_mapper(directory_in = 'taxonomy_downloads', directory_out = None):
    '''Generate taxonomy dict for lookups later on. Supports download_taxonomy.

    Parameters
    ----------
    directory_in: string
        Directory from which to read taxonomy information (default 'taxonomy_downloads')
    directory_out: string
        Deprecated - directory to deposit PCTAXA file

    Output
    ------
    PCTAXA file - pickled dict of taxonomy mapping.
    Naming will be Y-M-D formatted so you can remember when it was retrieved.
    '''
    
    #End file could be more memory efficient but works for now
    names = os.path.join(directory_in, 'names.dmp')
    nodes = os.path.join(directory_in, 'nodes.dmp')

    print("Reading taxon names.")
    name_dict = dict()
    with open(names) as input_file:
        for line in input_file:
            if 'scientific name' in line:
                line_split = line.split('\t')
                taxon = int(line_split[0])
                taxon_name = line_split[2]
                name_dict[taxon] = taxon_name

    print("Reading taxon nodes.")
    nodes_dict = dict()
    with open(nodes) as input_file:
        for line in input_file:
            line_split = line.split('\t')
            cur_id = int(line_split[0])
            parent_id = int(line_split[2])
            cur_taxon = line_split[4]
            nodes_dict[cur_id] = (parent_id, cur_taxon)

    print("Building dict...")

    taxa_dict = dict()
    for taxon in name_dict:
        taxon_dict = dict()
        taxon_dict['organism'] = name_dict.get(taxon)
        parent_id = None
        temp_taxon = taxon
        while parent_id != 1:
            parent_id, cur_taxon = nodes_dict.get(temp_taxon)
            cur_name = name_dict.get(temp_taxon)
            taxon_dict.update({cur_taxon:cur_name})
            temp_taxon = parent_id

        if 'species' in taxon_dict: #at or near terminal node            
            taxa_dict[taxon] = taxon_dict

    #Make sure to map contaminants to TaxonID == -1
    taxa_dict[-1] = {k: "CONTAMINANT" for k in ncbi_ranks}

    if not directory_out:
        directory_out = os.getcwd()

    todays_date = _get_todays_date()
    todays_file = todays_date + '.pctaxa'
    taxon_pickle = os.path.join(directory_out, todays_file)
    with open(taxon_pickle,'wb') as to_pickle:
        pickle.dump(taxa_dict, to_pickle)

    print("Finished taxonomy download. Assembled in file: ", todays_file)            
            
def load_taxonomy(file):
    '''Loads PCTAXA file into memory.

    Parameters
    ----------
    file: string
        A .pctaxa file created using the download_taxonomy
        function.

    Returns
    -------
    taxonomy dictionary: dict
        This dictionary contains all NCBI taxonomy mappings for an organism ID.

        dictionary[TaxID] = {'species': species, 'genus': genus, ... }

    Example
    -------
        >>> taxonomy = load_taxonomy('190101.pctaxa')
        >>> taxonomy[9606].get('species')
        Homo sapiens
    '''
    
    with open(file, 'rb') as input_file:
        return pickle.load(input_file)

def ncbi_check(taxa):
    '''Validates a list of taxa by making sure they are NCBI-valid ranks.

    Parameters
    ----------
    taxa: tuple or list
        List of taxonomic ranks, i.e. ('order','family')

    Returns
    -------
    list:
        List of only valid taxonomic ranks.
    '''
    
    return [x.lower() for x in taxa if x.lower() in ncbi_ranks]

def _unzip(file, members, output_path = None):
    '''Unzip files, primarily from NCBI downloads.'''
    if output_path is None:
        output_path = os.path.dirname(file)
    
    target_zipfile = ZipFile(file = file)
    for zip_member in members:
        target_zipfile.extract(zip_member, path = output_path)

#---Generic---
def _download(target, file_path, file_description = None):
    '''Downloads a file to a file_path, using a description if given one.'''
    if file_description:
        print(f"Download started for {file_description}.")
    else:
        print(f"Download started.")
    request_file = request.urlretrieve(target, file_path, reporthook = lambda *args: _fetch_progress(*args, description = file_description)) #make an intermediate function with lambda to take the 3 args and tack another on
    
def _fetch_progress(transferred_blocks, block_size, total_size, description):
    '''Prints file transfer progress'''
    size_transferred = transferred_blocks * block_size
    if size_transferred % 1000 == 0 and size_transferred != 0.0 and description:
        print(f'{size_transferred / 1e6} megabytes have been transferred for {description}.')

def _get_header_delim_reader(file, file_object):
    '''
    Returns header, delimiter, and reader from a file that will be opened
    with the csv module
    '''
    if file.endswith('.txt'):
        delim = '\t'
    elif file.endswith('.csv'):
        delim = ','
    else:
        raise Exception(
            'This file needs to be a csv or tab-delimited text file: ', file)

    if file_object:
        reader = csv.reader(file_object, delimiter = delim)
        header = next(reader)
    else: #Presumably not using the reader object here
        with open(file) as input_file:
            file_reader = csv.reader(input_file, delimiter = delim)
            header = next(file_reader)
            reader = None

    #Note that the reader is already advanced past the header
    return header, delim, reader

def _get_todays_date():
    '''Y-M-D formatted date, for PCTAXA and PCDB records.'''
    return datetime.datetime.today().strftime('%y%m%d')

def _enforce_worker_count(worker_count):
    '''
    Sets bottom and top limit for workers. Generally performance doesn't
    increase beyond workers in the range of 1:6.
    '''
    if worker_count is None or worker_count < 1:
        worker_count = max(1, min(mp.cpu_count() - 2, 6))
        return worker_count
    else:
        return worker_count
    

#---DB---
def db_fetch_params(db):
    '''Used in pcannotate.py but also may be useful for the user.
    Will retrieve the digest parameters specified when the database was created.

    Parameters
    ----------
    db: string
        PCDB (SQLite db) to connect to

    Returns
    -------
    results: tuple or None
        tuple of parameters
        (min_length, max_length, missed_cleavages, digest_rule, date_created).
        If parameters are not found, None is returned.
    '''
    c, conn = _db_connect(db)

    c.execute("SELECT * FROM DBParams")
    results = c.fetchone()
    conn.close()
    
    if results:
        return results
    else:
        return None

def db_stats(db):
    '''Prints database parameters from database for the user.

    Parameters
    ----------
    db: string
        .pcdb file created with create_pcdb; prints out stats.
    '''
    param_names = ["Min Length: ", "Max Length: ", "Missed Cleavages: ",
          "Digest Rule: ", "Date Created: "]
    param_vals = db_fetch_params(db)

    if param_vals:
        print_params = list(zip(param_names, [str(x) for x in param_vals]))
        print(printable_params)

def _db_connect(database):
    '''Handles connections to SQLite3 databases.
    Used in pcannotate.py and pcdb.py

    Parameters
    ----------
    database:string
        name of database to connect to

    Returns
    -------
    c, conn: tuple
        c: sqlite3.connect() cursor object for SQL statements
        conn: sqlite3.connect() object
    '''
    conn = sqlite3.connect(database)
    c = conn.cursor()
    return c, conn
    
if __name__ == "__main__":
    pass
