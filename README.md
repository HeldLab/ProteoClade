# ProteoClade

ProteoClade is a Python library for **taxonomic-based annotation and quantification of bottom-up proteomics data**. It is designed to be user-friendly, and has been optimized for speed and storage requirements.

ProteoClade helps you analyze two general categories of experiments:

1. **_Targeted Database_ Searches:** Experiments in which a limited number of species are defined ahead of time, such as those involving Patient-Derived Xenografts (PDXs) or host-pathogen interactions. Reference protein sequence databases are used for targeted searches (ex: using Mascot, MaxQuant).

2. **_De Novo_ Searches:** Experiments in which the organisms are unspecified ahead of time or involve samples of high taxonomic complexity. Mass spectra are analyzed in the absence of a reference database (ex: using PEAKS, PepNovo).

ProteoClade scales from two organisms to every organism in UniProt. Please [refer to the complete documentation at proteoclade.readthedocs.io](https://proteoclade.readthedocs.io) for installation, a user's guide, and examples.

Citation: Mooradian AD, van der Post S, Naegle KM, Held JM (2020) ProteoClade: A taxonomic toolkit for multi-species and metaproteomic analysis. PLOS Computational Biology 16(3): e1007741.
