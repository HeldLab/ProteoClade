'''''Constants. Typically options or valid column names'''

#---For in silico digesting---
cleave_rules = {  
    #C terminal proteases
    'trypsin/p':(r'[RK]', 'c'),
    'trypsin!p':(r'[RK](?!P)', 'c'),
    'lys-c':(r'[K]', 'c'),
    
    #N terminal proteases
    'asp-n':(r'[D]', 'n'),
    'asp-nc':(r'[DC]', 'n'),
    'lys-n':(r'[K]', 'n')
    }

#---For UniProt FASTAS
uniprot_options = {"SR":{"reviewed":"yes", "keyword":"1185"},
                    "S":{"reviewed":"yes"},
                    "T":{"reviewed":"no"},
                    "R":{"keyword":"1185"},
                    "A":None}

#---NCBI---
ncbi_ranks = {'no rank', 'forma', 'species subgroup', 'family', 'infraclass',
              'parvorder', 'superphylum', 'cohort', 'subtribe', 'tribe',
              'species group', 'species', 'subkingdom', 'superkingdom', 'superfamily',
              'infraorder', 'class', 'subclass', 'subgenus', 'kingdom', 'genus',
              'superclass', 'subspecies', 'phylum', 'varietas', 'subfamily',
              'suborder', 'subphylum', 'superorder', 'order'}

#---For annotation pcannotate.py---
valid_quant_cols = ("area", "intensity")
valid_sample_cols = ("scan","sample")
valid_pep_cols = ("sequence", "pep_seq", "peptide")
valid_score_cols = ("alc (%)", "score")

#---Anything PC might have added/expected. Used for implicit rollup in pcquant.py---

basic_annotations = ("organisms", "genes")

proteoclade_cols = set().union(
                        basic_annotations,
                        ncbi_ranks,
                        valid_quant_cols,
                        valid_sample_cols,
                        valid_pep_cols,
                        valid_score_cols
                        )


