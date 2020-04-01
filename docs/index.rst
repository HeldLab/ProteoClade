
What is ProteoClade?
====================
ProteoClade is a Python library for **taxonomic-based annotation and quantification of bottom-up proteomics data**. It is designed to be user-friendly, and has been optimized for speed and storage requirements. 

ProteoClade helps you analyze two general categories of experiments: 

1. **Targeted Database Searches:** Experiments in which a limited number of species are defined ahead of time, such as those involving Patient-Derived Xenografts (PDXs) or host-pathogen interactions. Reference protein sequence databases are used for targeted searches (ex: using Mascot, MaxQuant).

2. **De Novo Searches:** Experiments in which the organisms are unspecified ahead of time or involve samples of high taxonomic complexity. Mass spectra are analyzed in the absence of a reference database (ex: using PEAKS, PepNovo).

ProteoClade scales from two organisms to every organism in UniProt. This documentation includes examples from each kind of proteomic workflow to familiarize the user with ProteoClade's features.

Check out and cite the publication at `PLOS Computational Biology <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007741>`_.

Mooradian AD, van der Post S, Naegle KM, Held JM (2020) ProteoClade: A taxonomic toolkit for multi-species and metaproteomic analysis. PLOS Computational Biology 16(3): e1007741.


.. toctree::
   :maxdepth: 3
   :caption: Introduction
   
   
   introduction
   


.. toctree::
   :maxdepth: 3
   :caption: User Guide
   
   getting-started
   targeted-database
   de-novo
   tutorial

.. toctree::
   :maxdepth: 3
   :caption: Reference
   
   pcutilities
   pcdb
   pcannotate
   pcquant
   pcconstants

   