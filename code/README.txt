This README file was generated on 2023-01-11 by Bhavesh Borate.
It documents the contents of "code" folder for Cohen et al. manuscript: 
"First-in-human germline-targeting HIV nanoparticle vaccine induced broad and publicly-targeted helper T cell responses". 

NOTE: The current version of the code (which was used for the study) does not function on de-identified data.

Principal Investigator: William Schief, PhD
Institution: Scripps Research 
Email: schief@scripps.edu

CODE OVERVIEW

The code folder contains 3 .R files as follows:
1. Tcell-figures.R: 
  This code reads in cvd725_ics_2022DEC27.csv, cvd725_tfh_ics_2022DEC23.csv, and
  cvd725_association_TandBcell_2022DEC28.csv and plots Figures 1B, 1C, 1D 
  and 4 in the T-cell manuscript.

2. epitope_mapping_report.R:
	This code reads in cvd725_epitope_mapping_ics_2022DEC28.csv to map the 
	immunodominant T cell responses to eOD-GT8 and LumSyn (Figure 3).

3. UMAP-figures.R:
  This code reads in cvd725_cluster_ics_2022DEC27.csv and generates results 
  from the K-means clustering analysis (Figure 2).

