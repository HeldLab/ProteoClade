Getting Started
===============

Set Up Your Analysis
--------------------
Create a new folder to hold experimental files and results, and navigate to it in the command line of your operating system. This will be your **working directory**. ProteoClade will use this directory by default to look for sequence, taxonomy, and experimental files.

For experimental files, refer to the file format requirements for `targeted searches <targeted-database.html#file-format-requirements>`_ or `*de novo* searches <de-novo.html#file-format-requirements>`_.

Import ProteoClade
------------------
Common ProteoClade functions are most easily made available by importing everything from the package. Open a Python 3.6+ shell from the command line and then import::

    >>> from proteoclade import *

Assemble Required Files
-----------------------
ProteoClade needs to generate two files:

1. ProteoClade Database (PCDB)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to efficiently map experimental peptides to taxa and genes, a SQLite database (.pcdb) is generated using FASTA-formatted files. ProteoClade facilitates this by retrieving user-specified taxa from UniProt. 

Downloading FASTA files from UniProt can be accomplished using the function `download_uniprot <pcutilities.html#download-uniprot>`_::

    download_uniprot(*targets, download_folder = 'fastas')

Where targets are tuples containing (taxonomyID, database_type). The download folder will create a 'fastas' folder in your working directory if one does not already exist.

Example::

    >>> download_uniprot((9606,'sr'),(10090,'sr'))
    # Downloads the human and mouse Swiss-Prot reference proteomes

Alternatively, `download_uniprot_batch <pcutilities.html#download-uniprot-batch>`_ can be used to read a tab-delimited file to download multiple entries.

Optional step: once FASTA files have been downloaded, users may elect to include additional FASTA files by dragging and dropping them into the same directory as any that were downloaded by ProteoClade. User-supplied FASTA files are expected to have headers adhering loosely to the UniProt format. Only three fields are required::

    >db|UniqueIdentifier|EntryName OX=TaxonID GN=GeneName

Where TaxonID is the NCBI-assigned taxonomy identifier, and GeneName is the gene symbol. If no GeneName is found, the UniProt UniqueIdentifier is substituted.

Optional step: all FASTA files can be merged into a single file for database search engines using the function `merge_fastas <pcdb.html#merge-fastas>`_. This function will warn on the first incorrect header supplied but will not terminate. Note that incorrectly formatted headers will be discarded when a PCDB is created. 

Finally, the PCDB file can be created using the `create_pcdb <pcdb.html#create-pcdb>`_ function. ProteoClade enables the user to specify several parameters which should match the experimental, analytical, and informatics conditions of the experiment. An *in silico* digest is performed with these parameters so that experimentally-detected peptides can be matched. Please refer to the function's documentation for all of the parameters.

Example::

    >>> create_pcdb('human_mouse_swissprot')
	# Creates a pcdb file using trypsin and all default parameters

2. The ProteoClade Taxonomy Mapping File (PCTAXA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotating and taxa requires taxonomy mapping, which ProteoClade will source from NCBI. This is accomplished with the function `download_taxonomy <pcutilities.html#download-taxonomy>`_. Note that you should store the PCDB and PCTAXA files together, as taxonomy ID assignment may change over time. For ease of identification, contaminants are mapped with a TaxonID of -1 and a taxon name of “CONTAMINANT” at every level of the hierarchy, and this should be reflected in any contaminant FASTAs included (the function to download the cRAPome has this feature already built-in).

Example::

    >>> download_taxonomy()
	#Retrieves and builds the PCTAXA file in the working directory, named by the date it was retrieved