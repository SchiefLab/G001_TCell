This README file was generated on 2022-12-28 by Bhavesh Borate.
It documents the contents of "data" and "code" folders for Cohen et al. manuscript: 
"First-in-human germline-targeting HIV nanoparticle vaccine induced broad and publicly-targeted helper T cell responses". 

Principal Investigator: William Schief, PhD
Institution: Scripps Research 
Email: schief@scripps.edu

DATA & FILE OVERVIEW

The "data" folder contains 5 csv files as follows:
1. cvd725_ics_2022DEC27.csv: 
	This file contains data to plot the background-adjusted frequencies of cytokine-positive eOD-GT8 and LumSyn-specific CD4 (Figure 1B) and CD8 (Figure 1D) T cells 
	at week 10.
	
	Number of variables: 41
	Number of rows: 3359
	Variable List: study, pubid, visitno, visit_day, group,	treat, assay, antigen, parent, population, pop_table,	response_p, experiment_name, count_bg,	
	               parent_count_bg, drawdt, plate, runnum,	samp_ord, well_id, count, parent_count,	filter,	filter_reason, percent_cell, percent_cell_neg,	
		       percent_cell_net, percent_cell_net_trunc, filter_bg, response_adj_p, MIMOSA_response_prob, MIMOSA_N, Q2.5, Q25, Q50, Q75, Q97.5, 
		       response_fisher,	response_fisher_adj, response_MIMOSA, pop_grouping


2. cvd725_Tfh_ics_2022DEC23.csv:
	This file contains data to plot the frequencies of vaccine-specific cTfh-like cells measured in the PBMC at week 10 (Figure 1C)
	as the background-adjusted percent of IL-2+ or CD40L+ CXCR5+ CD4+ T cells out of total circulating CD4 T cells by stimulation.
	
	Number of variables: 12
	Number of rows: 180
	Variable List: study, pubid, week_visit, time_point, time_point_unit, group, treatment,	assay, antigen,	parent,	population, mag


3. cvd725_cluster_ics_2022DEC23.csv:
	This file contains data to generate results from the K-means clustering analysis (Figure 2) conducted on CD4 T cells that were positive for at least 2 cytokines or 
	activation markers in the eOD-GT8, LumSyn, negative control, and CMV stimulations by ICS assay.
	
	Number of variables: 23
	Number of rows: 31143
	Variable List: study, pubid, stim, replicate, num_pos_2_of_func, ifng, il2, tnfa, cd154, il17a, il4, gzb, cxcr5, pd-1, icos, cd45ra, ccr7, umap1, umap2, 
		       cluster, treat, visit_day, stim_treat


4. cvd725_epitope_mapping_ics_2022DEC28.csv:
	This file contains data to map the immunodominant T cell responses to eOD-GT8 and LumSyn (Figure 3).
	
	Number of variables: 48 
	Number of rows: 3519
	Variable List: study, pubid, peptide_id, peptide_grp, aa_sequence, start, end, antigen, population, parent, visitno, treat, response_p, specrole, count_bg, 
		       parent_count_bg, replicate, analysis_plan_id, exp_assay_id, group, plate, protocol, runnum, samp_ord, stdy_desc, well_id, count, parent_count, 
		       filter_1, filter_2, filter_3, filter_reason, percent_cell, percent_cell_neg, percent_cell_net, percent_cell_net_trunc, filter_bg, response_adj_p, 
		       parent_pop, mimosa_response_prob, mimosa_n, seed, response_fisher, response_fisher_adj, response_mimosa, parent_plot, peptide_group, 
		       ics_primary_filter


5. cvd725_association_TandBcell_2022DEC28.csv:
	This file contains data needed to associate the frequencies of T and B cell responses within and between peripheral and lymph node compartments (Figure 4). 
	
	Number of variables: 14
	Number of rows: 633
	Variable List: study, pubid, week_visit, group, treatment, time_point, time_point_unit, assay, parent, population, pop_name, mag, response, antigen

	
	


CODE OVERVIEW 

NOTE: The current version of the code (which was used for the study) does not function on de-identified data.

The "code" folder contains 3 .R files as follows:
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

