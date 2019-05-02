'''Rolls up peptides to the gene level, producing gene- and taxon specific outputs.

Functions for user
---------
roll_up
'''

import csv
import os

from .pcutilities import _get_header_delim_reader
from .pcconstants import proteoclade_cols
    
def roll_up(file, unique_taxon = 'species', inclusion_list = None,
            exclusion_list = None, samples = 'implicit', default_taxon = None,
            missing_values = True):
    '''Roll up peptides to the gene level, producing gene- and taxon specific outputs.

    Parameters
    ----------
    file: string
            .csv or tab .txt input file to roll peptides up to genes
    unique_taxon: string
        Taxon level which is used to determine uniqueness and identity (default “species”)
    inclusion_list: None or tuple
        Taxon members that must be matched for a peptide to be included, if specified (default None)
    exclusion_list: None or tuple
        Taxon members that must NOT be matched for a peptide to be included, if specified (default None)
    samples: string
        Method to use for finding samples in file, either “implicit” or “explicit” (default “implicit”)
    default_taxon: None or string
        Member of a taxonomy that will be assigned identity even if the peptide is shared with different members of the same taxonomic level (default None)
    missing_values: bool
        When adding samples, decide whether to include a gene if any samples have 0 or NaN in their quantitation (default True)

    Examples
    --------
        >>> roll_up("annotated_denovo_matched_experiment.txt") # species-specific gene rollup

    Notes
    -----
    Output is .csv or .txt file with unique taxon, unique genes, and samples
    '''

    print(f'Rolling up {file} to gene level.')

    with open(file) as input_file:
        header, delim, file_reader = _get_header_delim_reader(file, input_file)

        #Sample column
        if samples == 'implicit':
            #Something that proteoclade itself wouldn't work with
            sample_cols, samples = zip(*[(x,y) for (x,y) in enumerate(header) if y.lower() not in proteoclade_cols])
        elif samples == 'explicit':
            #Sample column starts with 'sample_'
            sample_cols, samples = zip(*[(x,y) for (x,y) in enumerate(header) if y.lower().startswith('sample_')])
        else:
            raise Exception('Invalid sample detection method.')        
        assert len(samples) == len(set(samples)), "Duplicate samples found."

        #Taxon column
        taxa_header_count = header.count(unique_taxon.lower())
        if taxa_header_count == 1:
            taxon_col = header.index(unique_taxon.lower())
        else:
            raise Exception(f'Found {taxa_header_count} taxa columns and need 1.')

        #Gene column
        assert header.count('genes') == 1, f'Genes column duplicated or missing.'
        gene_col = header.index('genes')

        assert inclusion_list == None or type(inclusion_list) in (tuple, list)
        assert exclusion_list == None or type(exclusion_list) in (tuple, list)

        if inclusion_list:
            inclusion_list = [x.lower() for x in inclusion_list]
        if exclusion_list:
            exclusion_list = [x.lower() for x in exclusion_list]

        sample_dict = dict() #much faster than dbm/shelve
        
        for row in file_reader:

            genes = row[gene_col].split('|')
            genes = [x.upper() for x in genes] #deal with species annotation issues, etc

            if len(set(genes)) > 1: #rollup to gene must be gene specific
                continue
            else:
                gene = genes.pop()

            #Set everything lowercase and then capitalize at end to handle input variations
            taxa = row[taxon_col].split('|')
            unique_taxa = set([x.lower() for x in taxa])

            if len(unique_taxa) > 1:
                if default_taxon and default_taxon.lower() in unique_taxa:
                    taxon_to_quant = default_taxon
                else:
                    continue #Not a unique taxon
            elif len(unique_taxa) == 0:
                continue
            else: #in the case of a single taxon in the taxon quant col
                taxon_to_quant = unique_taxa.pop()

            #Skip based on inclusion/exclusion rules
            if inclusion_list and taxon_to_quant not in inclusion_list:
                continue
            if exclusion_list and taxon_to_quant in exclusion_list:
                continue

            taxon_to_quant = taxon_to_quant.capitalize()

            taxon_gene = '|'.join((taxon_to_quant, gene))

            for col, sample in zip(sample_cols, samples):
                if taxon_gene not in sample_dict:
                    sample_dict[taxon_gene] = {}
                if sample not in sample_dict[taxon_gene]:
                    sample_dict[taxon_gene][sample] = 0

                try: #Skip non-numeric
                    float(row[col])
                except ValueError:
                    continue
                
                sample_dict[taxon_gene][sample] += float(row[col])
                
    
    ###Writing###
    file_out = f'rollup_{os.path.basename(file)}'
        
    with open(file_out, 'w', newline = '') as output_file:
        output_writer = csv.writer(output_file, delimiter = delim)
        
        output_header = ['Taxon','Gene']
        output_header.extend(samples)
        output_writer.writerow(output_header)

        for taxon_gene in sample_dict:

            row = taxon_gene.split('|')
            values = []
            assert len(row) == 2
            for sample in samples:
                sample_val = sample_dict[taxon_gene].get(sample, "NA")
                
                row.append(str(sample_val))
                values.append(sample_val)

            if not missing_values:
                if any([type(x) is str or int(x) == 0 for x in values]):
                    continue #skip row
                
            output_writer.writerow(row)

    print(f'Done rolling up file at {file_out}.')
    
if __name__ == '__main__':
    pass
