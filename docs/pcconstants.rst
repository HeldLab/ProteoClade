.. py:module:: proteoclade.pcconstants
.. py:currentmodule:: proteoclade.pcconstants

:py:mod:`pcconstants` Module
============================

* `cleave_rules`_
* `uniprot_options`_
* `ncbi_ranks`_
* `valid_quant_cols`_
* `valid_sample_cols`_
* `valid_pep_cols`_
* `valid_score_cols`_
* `basic_annotations`_
* `proteoclade_cols`_

cleave_rules
^^^^^^^^^^^^
A dictionary that maps a number of proteases to regular expressions.

uniprot_options
^^^^^^^^^^^^^^^
A dictionary that maps UniProt database types to their REST API equivalents.

ncbi_ranks
^^^^^^^^^^
A set of all valid taxonomic ranks.

valid_quant_cols
^^^^^^^^^^^^^^^^
A tuple of columns ProteoClade will assume are for quantification.

valid_sample_cols
^^^^^^^^^^^^^^^^^
A tuple of columns ProteoClade will assume are for samples.

valid_pep_cols
^^^^^^^^^^^^^^
A tuple of columns ProteoClade will assume are for peptide sequences.

valid_score_cols
^^^^^^^^^^^^^^^^
A tuple of columns ProteoClade will assume are for scores.

basic_annotations
^^^^^^^^^^^^^^^^^
A tuple ("organisms","genes") which ProteoClade will annotate by default.

proteoclade_cols
^^^^^^^^^^^^^^^^
A set of all columns ProteoClade uses and recognizes.