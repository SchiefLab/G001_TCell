# This code reads in cvd725_cluster_ics_2022DEC27.csv and generates results 
# from the K-means clustering analysis (Figure 2).

## ----install-packages, echo=FALSE, eval=TRUE, warning = FALSE, message = F------------------------------------------------------------------------------------------------

### Add additional packages needed here 
### Only works for CRAN packages (manually write library statement for Bioconductor or GitHub packages)
packages = c("MASS","exact2x2","Exact","tidyverse","knitr","kableExtra", 'DataSpaceR',"VISCfunctions","cowplot", "data.table",
             "ggplot2", "DataPackageR",   'Matrix', 'umap',  "here")
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    # VISCfunctions must be installed in separately (not in CRAN)
    if (x == 'VISCfunctions')
      stop('The package "VISCfunctions" must be installed through GitHub: https://github.com/FredHutch/VISCfunctions.git')
    
    install.packages(x, dependencies = TRUE,repos = "http://cran.us.r-project.org")
  }
  library(x, character.only = TRUE)
})

# Enforcing specific version of VISCfunctions
if (numeric_version(packageVersion('VISCfunctions')) < numeric_version('1.1'))
  stop('VISCfunctions must be at least version "1.1"')


# Create a theme
# This can be overloaded with other options
visc_theme <- theme(legend.position = "bottom", legend.margin = margin(unit = "cm"))

# VTN Color Palette
cbPalette <- c("#787873","#D92321","#1749FF","#0AB7C9","#FF6F1B","#810094","#FF5EBF","#8F8F8F")

# Setting up references depending on Word or PDf output
visc_ref <- function(ref_in) {
  ifelse(output_type == 'latex', 
          paste0('\\ref{', ref_in, '}'), 
          paste0('\\@ref(', ref_in, ')')) 
}

visc_clearpage <- function() {
  ifelse(output_type == 'latex', 
          '\\clearpage', 
         '#####')
} 


cluster_tab_org <- read.csv("N:/cavd/studies/cvd725/pdata/cohen_et_al_manuscript/data/cvd725_cluster_ics_2022DEC27.csv")

## ----marker-heatmap, fig.scap='Heatmap for marker positivity by cluster: the proportion of cells that are marker positive out of total number of cells in a cluster', fig.cap='Heatmap for marker positivity by cluster: the proportion of CD4+ T cells that are marker positive out of total number of cells in a cluster. Cell subset based on 2+ of 6 functional markers. Clusters determined using K-Means Clustering.', fig.height=7----

#look at main 4 rows (top) and then add in what's going on in the bottom row



fig2a <- ggplot(
  data = cluster_tab_org %>%
      filter(!Marker %in% c('CXCR3', 'CCR6', 'CD45RA', 'CCR7')) %>% 
  mutate(marker_factor = factor(Marker, levels = c("IFNg+","IL2+","TNFa+","CD154+",
                                                   "IL4+","IL17a+","GzB+","ICOS",
                                                   "CXCR5","PD-1", "EM", "CM", "TD", "Naive"),
                                labels = c(bquote("IFN-"*gamma*""),"IL-2",bquote("TNF-"*alpha*""),"CD40L",
                                                   bquote("IL"*'-'*"4"),bquote("IL"*'-'*"17a"),"GzB","ICOS",
                                                   "CXCR5","PD-1","EM","CM", "TEMRA", "Naive"),
                                ordered = T)),
  mapping = aes(x= factor(Cluster),
                y = fct_rev(marker_factor),
                fill = pos_prop * 100)) +
  geom_tile() +
  scale_fill_viridis_c('Marker\nPositivity\nRate',limits = c(0,100)) +
  xlab("Cluster") +  
  scale_y_discrete("Marker",labels = scales::parse_format())+
   theme(
                text=element_text(size=18),

        panel.spacing = unit(0, units = 'lines'))

  pdf(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig2a.pdf',
      width = 6, height = 10);
  
  fig2a
  
  dev.off()
  
  
  png(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig2a.png',
      res = 125,pointsize = 9,
      width = 800, height = 1200);
  
  fig2a
  
  dev.off()



## ----annotation-tab-------------------------------------------------------------------------------------------------------------------------------------------------------

#look at main 4 rows (top) and then add in what's going on in the bottom row


key_functions <- c("1" = "IFng+IL2+TNFa+154+ (4)",
                 "2" = "IFng-IL2+TNFa+154+ (3)",
                 "3" = "IFng-IL2+TNFa+154- (2)",
                 "4" = "IFng+IL2-TNFa+154+ (3)",
                 "5" = "IFng+IL2+TNFa+154+ (4)",
                 "6" = "IFng-IL2-TNFa+154- (1)",
                 "7" = "IFng-IL2-TNFa+154+ (2)",
                 "8" = "IFng+IL2-TNFa+154+ (3)",
                 "9" = "IFng-IL2+TNFa-154+ (2)",
                 "10" = "IFng-IL2+TNFa+154+ (3)")

other_features <- c("1" = "PD1+/EM",
                    "2" = "PD1+/EM",
                    "3" = "CM",
                    "4" = "GzB+/PD1+/EM",
                    "5" = "PD1+/CM/EM",
                    "6" = "IL17a+/PD1+/CM/EM",
                    "7" = "CM",
                    "8" = "PD1+/EM",
                    "9" = "PD1+/CM",
                    "10" = "CM")

annotation_tab <- bind_cols(c(1:10), key_functions, other_features)
names(annotation_tab) <- c("Cluster", "Key functions", "Other features")

kable(annotation_tab, col.names = c("Cluster", "Key functions", "Other features"),
      caption = "Cluster annotation, separated by key functions and other features. Numbers in parentheses indicate number of positive key functions out of 4.", escape = F, longtable = T, booktabs = T, linesep = "") %>% 
  kable_styling(font_size = 10, latex_options = c("hold_position", "repeat_header", "scale_down")) 



## ----percent-of-stim-in-cluster-------------------------------------------------------------------------------------------------------------------------------------------


# Getting CD4+ counts by pubid, visitno, Stim
all_samps_2_of_Func_n <- all_samps_2_of_Func %>%
  left_join(., dplyr::select(ics_adata, sample,treat) %>% distinct(), by = "sample") %>% 
  dplyr::count(sample, treat,visitno, Stim, name = 'pos_2_of_Func')



cluster_tab_stim <- bind_cols(
  all_samps_2_of_Func %>%
  left_join(., dplyr::select(ics_adata, sample,treat) %>% distinct(), by = "sample"),
  data.frame(Cluster = kmeans_out$cluster)
) %>%
  dplyr::count(treat,sample, visitno, Stim, Cluster, name = 'pos_cluster_2_of_Func') %>%
  # Filling missing clusters
  pivot_wider(names_from = 'Cluster',
              values_from = 'pos_cluster_2_of_Func',
              values_fill = 0) %>%
  pivot_longer(cols = -c(sample, treat,visitno, Stim), names_to = 'Cluster',
               values_to = 'pos_cluster_2_of_Func') %>%

  full_join(
    all_samps_2_of_Func_n,
    by = c("sample", "visitno", "Stim", "treat")
  ) %>%
  full_join(
    all_samps_n_counts ,
    by = c("sample",  "Stim")
  ) %>%
  mutate(
    percent_of_cd4 = 100 * pos_cluster_2_of_Func / cd4_pos,
    percent_of_2_of_Func = 100 * pos_cluster_2_of_Func / pos_2_of_Func,
    pubid = paste0('ID',
                   ifelse(as.numeric(factor(sample)) < 10, '0',''),
                   as.numeric(factor(sample))),
    Stim = factor(Stim),
    Cluster = factor(as.numeric(Cluster))
  ) %>%
  # getting BG subtracted
  group_by(sample, treat,visitno, Cluster) %>%
  mutate(
    percent_of_cd4_neg = percent_of_cd4[Stim == 'negctrl'],
    percent_of_cd4_net = percent_of_cd4 - percent_of_cd4_neg,
    `.groups` = "drop"
  ) %>%  mutate(treat = factor(treat ,
                              levels = c("DPBS sucrose",
                                         "20 µg eOD-GT8 60mer + AS01B",
                                         "100 µg eOD-GT8 60mer + AS01B"), ordered = T)) %>% 
  group_by(sample, Stim) %>% 
  mutate(
total_percent_net = sum(pmax(0, percent_of_cd4_net)),
rel_percent = 100 * pmax(0, percent_of_cd4_net) / total_percent_net
) %>% 
   mutate(treat_short = factor(case_when(treat == "DPBS sucrose" ~ "Placebo",
                                 treat == "20 µg eOD-GT8 60mer + AS01B" ~ "20 µg",
                                 treat == "100 µg eOD-GT8 60mer + AS01B" ~ "100 µg"), 
                               levels = c("Placebo","20 µg","100 µg"), ordered = T))
  







## ----percent-of-cd4-figs,  fig.scap=paste0("Percent positive cells, based on 2 of 6 function marker subset, for ", unique(cluster_tab_stim$Stim)[1:6], ", out of total CD4+ T-cells, at day 70 by treatment and cluster."), fig.cap=paste0("Percent positive cells, based on 2+ of 6 functional marker subset, for ", unique(cluster_tab_stim$Stim)[1:6], ", out of total CD4+ T-cells, at day 70 by treatment and cluster. Boxplots are superimposed on the distribution."), fig.height=6----

#do we want these by group or treatment?
#answers for cluster 1, how many 2 of 6 cells are there for each stim out of your total CD4 positive cells for each stim?


    
    fig2b <- cluster_tab_stim %>%
      filter(Stim %in% c("eOD-GT8", "LumSyn", "CMV")) %>%
      ggplot(aes(
        x = treat_short ,
        col = treat_short,
        y = pmax(0.001, percent_of_cd4_net)
      )) +
      geom_point(
        position = position_jitter(
          width = .1,
          height = 0,
          seed = 3241
        ),
        size = 2
      ) +
      geom_boxplot(
        fill = NA,
        lwd = .5,
        outlier.colour = NA,
        alpha = 1
      ) +
      facet_grid(Stim~ Cluster) +
      
      scale_x_discrete("") +
      scale_y_log10(
        paste0("% of CD4+ T cells"),
        breaks = c(.001, .01, .1, 1, 10, 100),
        labels = c(expression("" <= "0.001"),
                   '0.01', '0.1', "1", "10", "100")
      ) + theme_bw() +
      scale_color_manual(
        name = "",
        values = c("#D92321", "#1749FF", "#0AB7C9", "#FF6F1B"),
        limits = c(
          "Placebo",
          "20 µg",
          "100 µg"
        )
      ) + theme_bw()+
      theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
                text=element_text(size=18),
                axis.text.x = element_blank(),

        legend.position = "bottom",
        panel.spacing = unit(0, units = 'lines'))
      


pdf(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig2b.pdf',
    width = 10, height = 12);

fig2b
dev.off()


png(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig2b.png',
    res = 200,pointsize = 5,
    width = 1800, height = 2200);

fig2b

dev.off()



## ----percent-of-cd4-pooled,  fig.cap="Percent positive cells, based on 2 of 6 function marker subset, out of total CD4+ T-cells, at day 70 by antigen and cluster. Placebo recipients are not included, and vaccine recipients are pooled.", fig.scap="Percent positive cells, based on 2+ of 6 functional marker subset, out of total CD4+ T-cells, at day 70 by antigen and cluster.", fig.height=7.25----

#do we want these by group or treatment?
#answers for cluster 1, how many 2 of 6 cells are there for each stim out of your total CD4 positive cells for each stim?
plot_data <- cluster_tab_stim %>%
  mutate(Cluster = as.integer(Cluster)) %>%
  mutate(
    treat_pooled = case_when(
      treat != "DPBS sucrose" ~ "eOD-GT8 pooled",
      treat == "DPBS sucrose" ~ "DPBS sucrose"
    )
  ) %>%
  filter(treat_pooled == "eOD-GT8 pooled" & Stim != "negctrl") %>%
  filter(!(Stim %in% c("eOD-GT8-1", "eOD-GT8-2", "eOD-GT8-3"))) %>%
  full_join(annotation_tab, by = "Cluster") %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  mutate(annotation_plot = factor(
    paste0(Cluster, ": ", `Key functions`, "/", `Other features`),
    levels = c(
      "1: IFng+IL2+TNFa+154+ (4)/PD1+/EM",
      "2: IFng-IL2+TNFa+154+ (3)/PD1+/EM",
      "3: IFng-IL2+TNFa+154- (2)/CM",
      "4: IFng+IL2-TNFa+154+ (3)/GzB+/PD1+/EM",
      "5: IFng+IL2+TNFa+154+ (4)/PD1+/CM/EM",
      "6: IFng-IL2-TNFa+154- (1)/IL17a+/PD1+/CM/EM",
      "7: IFng-IL2-TNFa+154+ (2)/CM",
      "8: IFng+IL2-TNFa+154+ (3)/PD1+/EM",
      "9: IFng-IL2+TNFa-154+ (2)/PD1+/CM",
      "10: IFng-IL2+TNFa+154+ (3)/CM"
    )
  )) %>%
  mutate(annotation_plot_wrap = factor(
    paste0(Cluster, ": ", `Key functions`, "\n", `Other features`),
    levels = c(
      "1: IFng+IL2+TNFa+154+ (4)\nPD1+/EM",
      "2: IFng-IL2+TNFa+154+ (3)\nPD1+/EM",
      "3: IFng-IL2+TNFa+154- (2)\nCM",
      "4: IFng+IL2-TNFa+154+ (3)\nGzB+/PD1+/EM",
      "5: IFng+IL2+TNFa+154+ (4)\nPD1+/CM/EM",
      "6: IFng-IL2-TNFa+154- (1)\nIL17a+/PD1+/CM/EM",
      "7: IFng-IL2-TNFa+154+ (2)\nCM",
      "8: IFng+IL2-TNFa+154+ (3)\nPD1+/EM",
      "9: IFng-IL2+TNFa-154+ (2)\nPD1+/CM",
      "10: IFng-IL2+TNFa+154+ (3)\nCM"
    )
  ))



      ggplot(data = plot_data, aes(
        x = Stim ,
        col = Stim,
        y = pmax(0.001, percent_of_cd4_net)
      )) +
      geom_point(
        position = position_jitter(
          width = .1,
          height = 0,
          seed = 3241
        ),
        size = .85
      ) +
      geom_boxplot(
        fill = NA,
        lwd = .5,
        outlier.colour = NA,
        alpha = 1
      ) +
      facet_wrap(~ annotation_plot_wrap, nrow = 2) +
      
      scale_x_discrete("") +
      scale_y_log10(
        paste0("% of CD4+ T cells"),
        breaks = c(.001, .01, .1, 1, 10, 100),
        labels = c(expression("" <= "0.001"),
                   '0.01', '0.1', "1", "10", "100")
      ) + theme_bw() +
      scale_color_manual('', values = c("#0000FF","#1B9E77","#0AB7C9")) +
      theme(
        strip.text.y = element_text(size = 6),
                strip.text.x = element_text(size = 6),
        axis.text.x = element_text(
          size = 8,
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        legend.text = element_text(size = 7),
        legend.position = "bottom",
        panel.spacing = unit(0, units = 'lines'))
      
  






## ----rel-percent-of-cd4-pooled-bar,  fig.scap="Bar plot of relative percent positive cells, based on 2 of 6 function marker subset, out of total CD4+ T-cells, at day 70 by antigen and cluster.", fig.cap="Bar plot of relative percent positive cells, based on 2+ of 6 functional marker subset, out of total CD4+ T-cells, at day 70 by antigen and cluster. Placebo recipients are not included, and vaccine recipients are pooled. Error bars depicting 95\\% confidence intervals are displayed over each bar.", fig.height=7.25----


plot_data <- cluster_tab_stim %>% 
  mutate(treat_pooled = case_when(treat != "DPBS sucrose" ~ "eOD-GT8 pooled",
                                    treat == "DPBS sucrose" ~ "DPBS sucrose")) %>% 
  filter(treat_pooled == "eOD-GT8 pooled" & Stim != "negctrl") %>% 
  filter(!(Stim %in% c("eOD-GT8-1","eOD-GT8-2","eOD-GT8-3")))

summary_plot_dat <- Rmisc::summarySE(plot_data, measurevar="rel_percent", groupvars=c("Stim","Cluster")) %>% 
  mutate(Cluster = as.integer(Cluster)) %>% 
  full_join(annotation_tab, by = "Cluster") %>% 
  mutate(Cluster = as.factor(Cluster)) %>% 
  group_by(Stim) %>% 
  mutate(annotation_plot= factor(paste0(Cluster, ": ", `Key functions`, "/", `Other features`),
                                 levels = c("1: IFng+IL2+TNFa+154+ (4)/PD1+/EM",
                                            "2: IFng-IL2+TNFa+154+ (3)/PD1+/EM",
                                            "3: IFng-IL2+TNFa+154- (2)/CM",
                                            "4: IFng+IL2-TNFa+154+ (3)/GzB+/PD1+/EM",
                                            "5: IFng+IL2+TNFa+154+ (4)/PD1+/CM/EM",
                                            "6: IFng-IL2-TNFa+154- (1)/IL17a+/PD1+/CM/EM",
                                            "7: IFng-IL2-TNFa+154+ (2)/CM",
                                            "8: IFng+IL2-TNFa+154+ (3)/PD1+/EM",
                                            "9: IFng-IL2+TNFa-154+ (2)/PD1+/CM",
                                            "10: IFng-IL2+TNFa+154+ (3)/CM"))) %>% 
    mutate(annotation_plot_wrap= factor(paste0(Cluster, ": ", `Key functions`, "\n", `Other features`),
                                 levels = c("1: IFng+IL2+TNFa+154+ (4)\nPD1+/EM",
                                            "2: IFng-IL2+TNFa+154+ (3)\nPD1+/EM",
                                            "3: IFng-IL2+TNFa+154- (2)\nCM",
                                            "4: IFng+IL2-TNFa+154+ (3)\nGzB+/PD1+/EM",
                                            "5: IFng+IL2+TNFa+154+ (4)\nPD1+/CM/EM",
                                            "6: IFng-IL2-TNFa+154- (1)\nIL17a+/PD1+/CM/EM",
                                            "7: IFng-IL2-TNFa+154+ (2)\nCM",
                                            "8: IFng+IL2-TNFa+154+ (3)\nPD1+/EM",
                                            "9: IFng-IL2+TNFa-154+ (2)\nPD1+/CM",
                                            "10: IFng-IL2+TNFa+154+ (3)\nCM")))



fig2c <-  ggplot(data = summary_plot_dat, aes(fill = Stim,
        x = Stim ,
        col = Stim,
        y = pmax(1,rel_percent)
      )) +
          geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
    geom_errorbar(position=position_dodge(.9),colour = "black", width=.25, aes(ymin=pmax(1,rel_percent-ci), ymax=pmax(1,rel_percent+ci))) +
      facet_wrap(~ annotation_plot_wrap, nrow = 2) +
      scale_x_discrete("") +
      scale_y_continuous(
        paste0("Relative % of CD4+ T cells"),
        breaks = c(1, 10, 20, 30, 40, 50),
        labels = c(expression("" <= "1"),"10","20","30","40", "50")
      ) + theme_bw() +
      scale_color_manual('', values = c("#0000FF","#1B9E77","#0AB7C9")) +
         scale_fill_manual('', values = c("#0000FF","#1B9E77","#0AB7C9")) +
      theme(
        text=element_text(size=18),
        legend.text = element_text(size=24),
        axis.text.x = element_text(
          size = 14
        ),
        legend.position = "none",
        panel.spacing = unit(0, units = 'lines'))


pdf(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig2c.pdf',
    width = 14, height = 10);

fig2c

dev.off()


png(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig2c.png',
    res = 150,pointsize = 5,
    width = 2200, height = 1500);

fig2c

dev.off()



## ----load-data-cd8, include=FALSE, eval = T-------------------------------------------------------------------------------------------------------------------------------

load(file = '/networks/cavd/studies/cvd725/scratch/fi_cd8.RData')
load(file = '/networks/cavd/studies/cvd725/scratch/gates_cd8.RData')
load(file = '/networks/cavd/studies/cvd725/scratch/meta_data_cd8.RData')
load(file = '/networks/cavd/studies/cvd725/scratch/counts_dat_cd8.RData')

# data_path <- 'N:\\cavd\\studies\\cvd725\\pdata\\ICS_sharing'

# all_samps_gates_org <- readRDS(file.path(data_path,'gates_2_of_7_pos_03June2021.rds'))
all_samps_gates_org_cd8 <- gates_cd8 
  



# all_samps_meta_org <- readRDS(file.path(data_path,'meta_08June2021.rds'))
all_samps_meta_org_cd8 <- meta_data_cd8 %>% 
  mutate(antigen = str_remove(antigen, "-HxB2"),
         Stim = str_remove(Stim, "-HxB2")) %>% 
  mutate(antigen = case_when(antigen == "eOD-GT8-LumSyn" ~ "LumSyn",
                             antigen != "eOD-GT8-LumSyn" ~ antigen),
         Stim = case_when(Stim == "eOD-GT8-LumSyn" ~ "LumSyn",
                             Stim != "eOD-GT8-LumSyn" ~ Stim))




 # all_samps_counts_org <- readRDS(file.path(data_path,'UMAP_counts_03June2021.rds'))
all_samps_counts_org_cd8 <- full_join(
  meta_data_cd8, counts_dat_cd8, 
  by = c("Type", "ANALYSIS_PLAN_ID", "Stim", "name", "Comp", "Plate", 
         "PLATE ID", "PLATE NAME", "EXPERIMENT NAME", "TUBE NAME", 
         "Sample Order", "Replicate", "EXP_ASSAY_ID")) %>% 
  mutate(antigen = str_remove(antigen, "-HxB2"),
         Stim = str_remove(Stim, "-HxB2")) %>% 
  mutate(antigen = case_when(antigen == "eOD-GT8-LumSyn" ~ "LumSyn",
                             antigen != "eOD-GT8-LumSyn" ~ antigen),
         Stim = case_when(Stim == "eOD-GT8-LumSyn" ~ "LumSyn",
                             Stim != "eOD-GT8-LumSyn" ~ Stim)) %>% 
  full_join(ics %>% distinct(sample, treat, group), by = 'sample')


 method_here <- 'rlm'
             


## ----pops-of-int-cd8------------------------------------------------------------------------------------------------------------------------------------------------------


CD4_CUTOFF <- 10000
FISHER_RESPONSE_CUTOFF <- .00001
################



cd8_cytokine_names <-
  c( 'IFNg' = expression('IFN-' * gamma * '+'),
    'IL2' = expression('"IL-2+"'),
    'TNFa' = expression('TNF-' * alpha * '+'),
    'IFNg_OR_IL2_OR_TNFa' = expression(atop('IFN-' * gamma * '+ and/or IL-2+',
                                            'and/or TNF-' * alpha * '+')),
    'CD154' = expression('TNF-' * alpha * '+'),
    '2_of_Th1' = expression(atop('2+ of IFN-' * gamma * '+,', 'IL-2+, TNF-' * alpha * '+')),
    '2_of_TH1_CD154_IL17a_GzB' = expression("2 of TH1 CD154 IL17a GzB"))


cd8_cytokine_names_tab <-
  c('IFNg' = 'IFN-$\\gamma$+',
    'IL2' = 'IL-2+',
    'TNFa' = 'TNF-$\\alpha$+',
    'IFNg_OR_IL2' = 'IFN-$\\gamma$+ and/or IL-2+',
    'IFNg_OR_IL2_OR_TNFa' = 'IFN-$\\gamma$+ and/or IL-2+ and/or TNF-$\\alpha$+',
    '2_of_Th1' = '2+ of IFN-$\\gamma$+, IL-2+, TNF-$\\alpha$+',
    '2_of_TH1_CD154_IL17a_GzB' = '2 of TH1 CD154 IL17a GzB',
    "GzB" = "GzB+",
    'IL17a' = 'IL17-$\\alpha$+',
    'IL4' = 'IL-4+',
    "CD45RA" = "CD45RA",
    "CCR7" = "CCR7",
    "CXCR5" = "CXCR5")


parent_gate_cd8 <- "/S/Exclude/14-/Lv/L/3+/8+"


functional_markers_cd8 <- c("IFNg" = "IFNg+",
                            "IL2" = "IL2+",
                            "TNFa" = "TNFa+",
                            "CD154" = "8+154+",
                            "IL17a" = "IL17a+",
                            "IL4" = "IL4+")

additional_markers_cd8 <- c("GzB" = "GzB+",
                            "CXCR5" = 'CXCR5+',
                            "PD-1" = '8+PD1+',
                            "ICOS" = '8+ICOS+',
                             "CD45RA" = "CD45RA",
                        "CCR7" = "CCR7")



boolean_cd8_pops <- c('EM' = "EM",'CM' = "CM",'Naive' = "Naive",'TD' = "TD")

names(boolean_cd8_pops) <- boolean_cd8_pops

all_cd8_pops <- c(functional_markers_cd8, additional_markers_cd8, boolean_cd8_pops)


all_cd8_pops_path <- paste0(parent_gate_cd8, '/', all_cd8_pops)






## ----data-processing-cd8--------------------------------------------------------------------------------------------------------------------------------------------------



#merge on "name" - first part of gates rowname - create new variables "name" and "cell_id"
all_samps_cd8 <- full_join(
  all_samps_meta_org_cd8,
  all_samps_gates_org_cd8 %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rename_with(basename) %>% 
    rownames_to_column() %>% 
    mutate(name = str_sub(rowname, end = str_locate(rowname,'_')[,'start'] - 1)), 
  by = "name") %>% 
  filter(Stim != "sebctrl", !is.na(rowname))


#Getting number positive for 2 of 4
  all_samps_cd8$num_pos_2_of_Th1_GzB <- factor(rowSums(all_samps_cd8[, c("IFNg", "IL2", "TNFa")]) +
  all_samps_cd8$GzB * all_samps_cd8$IFNg_OR_IL2_OR_TNFa, levels = c(2, 3, 4), ordered = T)

  
  
  
# Creating marker vs. 2 of Th1, CD154, IL-17a, GzB results
percent_by_2_of_6_cd8 <- all_samps_cd8 %>%
  group_by(Stim) %>%
  dplyr::summarise(
    across(`IFNg`:`2_of_TH1_CD154_IL17a_GzB`, sum),
    parent_count = n(),
    `.groups` = "drop"
  ) %>% pivot_longer(
    cols = `IFNg`:`2_of_TH1_CD154_IL17a_GzB`,
    names_to = 'population',
    values_to = 'count'
  )

  
  



## ----adata-cd8------------------------------------------------------------------------------------------------------------------------------------------------------------


ics_adata_cd8 <- all_samps_counts_org_cd8 %>% 
  mutate(
    visitno = as.numeric(visitno),
    pop_plot = factor(population,
                      levels = names(cd8_cytokine_names),
                      labels = cd8_cytokine_names
    ),
    pop_tab = factor(population,
                     levels = names(cd8_cytokine_names_tab),
                     labels = cd8_cytokine_names_tab
    )
  ) %>% 
  mutate(
    visit_day = "Day 70",
    treat_plot = factor(treat, 
                        levels = c("DPBS sucrose", 
                                   "20 µg eOD-GT8 60mer + AS01B", 
                                   "100 µg eOD-GT8 60mer + AS01B"),
                        labels = c("DPBS sucrose", 
                                   "20 µg", 
                                   "100 µg"), ordered = T),
    treat = factor(treat ,
                   levels = c("DPBS sucrose",
                              "20 µg eOD-GT8 60mer + AS01B",
                              "100 µg eOD-GT8 60mer + AS01B"), ordered = T),
    percent_cell = count / parent_count * 100)




## ----summary-stats-cd8----------------------------------------------------------------------------------------------------------------------------------------------------

ics_stats_all_cd8 <- 
  bind_rows(ics_adata_cd8,
                           ics_adata_cd8 %>% 
                             mutate(antigen = 'Total')
                           ) %>%
  full_join(dplyr::select(all_samps_meta_org_cd8, name, antigen, sample), by = c("sample", "antigen")) %>% 
                             filter(!is.na(population)) %>% 
  group_by(group,antigen, population, pop_tab) %>% 
  dplyr::summarise(
    total = length(unique(sample)),
    stat_sum = sum(count),
    stat_mean_info = 10^mean(log10(pmax(count,1))),
    stat_med_info = 10^median(log10(pmax(count,1))),
    stat_info = paste0(stat_sum, ' / ', round_away_0(stat_mean_info, digits = 0), 
                       ' / ', round_away_0(stat_med_info, digits = 0)),
    `.groups` = "drop"
  ) %>% 
  mutate(visit_day = "Day 70") 

ics_stats_total_cd8 <- 
  bind_rows(ics_adata_cd8,
                           ics_adata_cd8 %>% 
                             mutate(antigen = 'Total')
                           ) %>%
  full_join(dplyr::select(all_samps_meta_org_cd8, name, antigen, sample), by = c("sample", "antigen")) %>% 
                             filter(!is.na(population)) %>% 
  group_by(antigen, population, pop_tab) %>% 
  dplyr::summarise(
    total = length(unique(sample)),
    stat_sum = sum(count),
    stat_mean_info = 10^mean(log10(pmax(count,1))),
    stat_med_info = 10^median(log10(pmax(count,1))),
    stat_info = paste0(stat_sum, ' / ', round_away_0(stat_mean_info, digits = 0), 
                       ' / ', round_away_0(stat_med_info, digits = 0)),
    `.groups` = "drop"
  ) %>% 
  dplyr::select(population,pop_tab, antigen, stat_info) %>% 
   filter(!(pop_tab %in% boolean_cd8_pops)) %>% 
  filter(!is.na(pop_tab)) %>% 
  mutate(Marker = pop_tab, `Peptide pool` = antigen)






## ----fi-dataset-cd8, cache=TRUE, cache.lazy = FALSE, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------

# need to subsample CMV 
# n_cmv <- sum(!all_samps$Stim %in% c('CMV', 'negctrl'))
n_cmv <- sum(all_samps_cd8$Stim == c('LumSyn'))
sample_ns <- all_samps_cd8 %>% count(Stim) %>% 
  mutate(n = if_else(Stim == 'CMV', n_cmv, n))

set.seed(51384152)
sample_index <- all_samps_cd8 %>% 
  full_join(sample_ns, by = 'Stim') %>% 
  mutate(index = 1:n()) %>%   
  group_by(Stim) %>% 
  sample_n(n) %>% 
  pull('index')

# Creating marker vs. 2 of Th1, CD154, IL-17a, GzB results
percent_by_2_of_6_subset_cd8 <- all_samps_cd8 %>%
  slice(sample_index) %>%
  group_by(Stim) %>%
  dplyr::summarise(
    across(`IFNg`:`2_of_TH1_CD154_IL17a_GzB`, sum),
    parent_count = n(),
    `.groups` = "drop"
  ) %>% pivot_longer(
    cols = `IFNg`:`2_of_TH1_CD154_IL17a_GzB`,
    names_to = 'population',
    values_to = 'count'
  )

all_samps_2_of_TH1_CD154_IL17a_GzB_cd8 <- all_samps_cd8 %>% 
  filter(`2_of_TH1_CD154_IL17a_GzB`, Stim != 'sebctrl')


# Getting CD8+ counts by pubid, visitno, Stim
all_samps_n_cd8 <- all_samps_cd8 %>% 
  dplyr::count(sample, visitno, Stim, name = 'cd8_pos')


all_samps_n_counts_cd8 <- all_samps_counts_org_cd8 %>% 
  group_by(sample,Stim) %>% 
  dplyr::summarize(cd8_pos = unique(parent_count), `.groups` = "drop")

# all_samps_fi_org <- readRDS(file.path(data_path,'fi_2_of_func_pos_03June2021.rds'))

all_samps_fi_org_cd8 <- fi_cd8
#continuous data that you put into UMAP and clustering, jimmy brought in later because of memory


all_samps_fi_pre_cd8 <-
  right_join(all_samps_meta_org_cd8,
             all_samps_fi_org_cd8 %>%
               mutate(name = str_sub(cell_id, 
                                     end = str_locate(cell_id,'_')[,'start'] - 1)), 
             by = "name") %>% 
  remove_rownames() %>% 
  # column_to_rownames('cell_id') %>% 
  pivot_longer(IFNg:CCR7,names_to = "channel", values_to = "value") %>% 
  group_by(channel) %>% 
  mutate(mean_val = mean(value),
         sd_val = sd(value)) %>% 
  mutate(mean_plus_3sig = mean_val + (1.5 * sd_val),
         mean_minus_3sig = mean_val - (1.5 * sd_val)) %>% 
  #truncate
  mutate(val_trunc = case_when(value > mean_minus_3sig & value < mean_plus_3sig ~ value,
                               value < mean_minus_3sig ~ mean_minus_3sig,
                               value > mean_plus_3sig ~ mean_plus_3sig)) %>% 
  #standardize
  mutate(val_std = (val_trunc - mean(val_trunc))/sd(val_trunc)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = channel, values_from = c(value,mean_val, sd_val, mean_plus_3sig, mean_minus_3sig, val_trunc, val_std)) 

all_samps_fi_cd8 = all_samps_fi_pre_cd8 %>% 
  slice(sample_index) %>% 
 #selecting only standardized values to run UMAP - actually not for now
  dplyr::select(value_IFNg:value_CCR7)
  
names(all_samps_fi_cd8) <- str_remove(names(all_samps_fi_cd8), "value_")

#rm(all_samps_fi_org)

# UMAP creation
 # umap.defaults
 # n_neighbors: 15
 # n_components: 2

 custom.config = umap.defaults
 custom.config$random_state = 7485613
 custom.config$n_neighbors = 70
 custom.config$min_dist = .2
 custom.config$metric = 'euclidean'
#
#
 full_umap_cd8 = umap(all_samps_fi_cd8,
                  config = custom.config,
                  verbose = TRUE)
# 
# kmeans clustering

n_clust <- 8

# Getting random seed
# sample(1:2^25, 1)
# 9030497

set.seed(9030497)
kmeans_out_cd8 <- kmeans(
  all_samps_fi_cd8,
  centers = n_clust,
  iter.max = 10000,
  #change nstart if not converging
  nstart = 100,
  algorithm = "Hartigan-Wong",
  trace = FALSE
)


#creating umap data

umap_data_cd8 <- bind_cols(
  all_samps_2_of_TH1_CD154_IL17a_GzB_cd8 %>% 
    slice(sample_index) %>% 
    select(sample, Stim, Replicate,num_pos_2_of_Th1_GzB),
  all_samps_fi_cd8 %>% 
  #want log10 value
    mutate(across(.fns = ~log10(pmax(1, .x)))),
  full_umap_cd8$layout %>% as.data.frame(),
  data.frame(Cluster = kmeans_out_cd8$cluster)
  ) %>% 
    full_join(., dplyr::select(Feinberg725_ICS_flowjo_data, sample, treat) %>% distinct(), by = "sample") %>% 

  mutate(
    Stim = factor(Stim),
    visit_day = "Day 70",
    Stim_treat = paste0(Stim, ' (', treat, ')')
  ) %>% 
  rename(`UMAP1` = 'V1', `UMAP2` = 'V2') 




## ----cluster-mag-tab-cd8, eval = T, results="asis", warning=kable_warnings------------------------------------------------------------------------------------------------
#of your cluster, how many are positive for each population? each column kind of summarizes what each cluster contains. annotation which is provided by the lab
cluster_tab_org_cd8 <- bind_cols(
  all_samps_2_of_TH1_CD154_IL17a_GzB_cd8 %>%
    slice(sample_index),
  data.frame(Cluster = kmeans_out_cd8$cluster)) %>% 

  pivot_longer(cols = IFNg:CCR7,
               names_to = 'Marker',
               values_to = 'Response') %>%
  mutate(Marker = factor(Marker, levels = names(all_cd8_pops),
                           labels = all_cd8_pops)) %>%
  group_by(Marker) %>%
  mutate(marker_pos_all = sum(Response)) %>%
  group_by(Cluster, Marker, marker_pos_all) %>%
  dplyr::summarise(
    cell_pos = sum(Response),
    cell_count = n(),
    `.groups` = "drop"
  ) %>%
  mutate(
    pos_prop = cell_pos / cell_count,
    pos_info = paste0(round_away_0(100 * pos_prop, 2),
                      '(', cell_pos, ')'),
    cluster_prop = cell_pos / marker_pos_all,
    cluster_pos_info = paste0(round_away_0(100 * cluster_prop, 2),
                              '(', cell_pos, ')'),
    Cluster_count = paste0(Cluster, ' (n=', cell_count, ')'),
    Marker_count = paste0(Marker, ' (n=', marker_pos_all, ')'),
    Grouping = case_when(
      !(Marker %in% c('EM','CM','Naive','TD') ) ~ 'Populations in UMAP',
      Marker %in% c('EM','CM','Naive','TD') ~ 'Memory Subsets'
    )
  )

marker_pos_tab_cd8 <- cluster_tab_org_cd8 %>%
  mutate(
    Marker = escape(Marker)
  ) %>%
  pivot_wider(id_cols = c(Grouping, Marker),
              names_from = Cluster_count,
              values_from = pos_info)

# 
marker_pos_tab_cd8 %>%
      mutate(Grouping = ifelse(duplicated(Grouping), NA, as.character(Grouping)) )%>%
  kable(
    format = output_type, longtable = FALSE, booktabs = T, linesep = "", escape = FALSE,
    caption.short = "Marker positivity by cluster: the proportion of cells that are marker positive out of total number of cells in a cluster",
    caption = 'Marker positivity by cluster: the proportion of CD8+ T cells that are marker positive out of total number of cells in a cluster. Cell subset based on 2+ of 6 functional markers. Clusters determined using K-Means Clustering, and clusters with cell counts < 10 (clusters 4 and 5) were removed.'
  ) %>%
  kable_styling(
    font_size = 6.5,
    latex_options = c("hold_position")
    )%>%
  # row_spec(which(marker_pos_tab$Marker %in% c('Ki67', 'TD')), hline_after = TRUE) %>%
  add_header_above(c(' ' = 2, 'Cluster Number (n)' = 8)) %>%
  landscape()



## ----marker-heatmap-cd8,eval = T, fig.scap='Heatmap for marker positivity by cluster: the proportion of cells that are marker positive out of total number of cells in a cluster', fig.cap='Heatmap for marker positivity by cluster: the proportion of CD8+ T cells that are marker positive out of total number of cells in a cluster. Cell subset based on 2+ of 6 functional markers. Clusters determined using K-Means Clustering.', fig.height=7----

#look at main 4 rows (top) and then add in what's going on in the bottom row


fig3a <- ggplot(
  data = cluster_tab_org_cd8 %>%
    filter(!Cluster %in% c(4, 5)) %>% 
      filter(!Marker %in% c('CXCR3', 'CCR6', 'CD45RA', 'CCR7')) %>% 
  mutate(marker_factor = factor(Marker, levels = c("IFNg+","IL2+","TNFa+","8+154+",
                                                   "IL4+","IL17a+","GzB+","8+ICOS+",
                                                   "CXCR5+","8+PD1+", "EM", "CM", "TD", "Naive"),
                                labels = c(bquote("IFN-"*gamma*""),"IL-2",bquote("TNF-"*alpha*""),"CD40L",
                                                   bquote("IL"*'-'*"4"),bquote("IL"*'-'*"17a"),"GzB","ICOS",
                                                   "CXCR5","PD-1","EM","CM", "TEMRA", "Naive"),
                                ordered = T)),
  mapping = aes(x= factor(Cluster),
                y = fct_rev(marker_factor),
                fill = pos_prop * 100)) +
  geom_tile() +
  scale_fill_viridis_c('Marker\nPositivity\nRate',limits = c(0,100)) +
  xlab("Cluster") +  
  scale_y_discrete("Marker",labels = scales::parse_format())+
   theme(
                text=element_text(size=18),

        panel.spacing = unit(0, units = 'lines'))


  pdf(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig3a.pdf',
      width = 6, height = 10);
  
  fig3a
  
  dev.off()
  
  
  png(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig3a.png',
      res = 125,pointsize = 9,
      width = 800, height = 1200);
  
  fig3a
  
  dev.off()





## ----annotation-tab-cd8, results = "asis", eval = T-----------------------------------------------------------------------------------------------------------------------

#look at main 4 rows (top) and then add in what's going on in the bottom row


key_functions_cd8 <- c("1" = "IFng+IL2+TNFa+154- (3)",
                 "2" = "IFng+IL2-TNFa+154- (2)",
                 "3" = "IFng-IL2+TNFa+154- (2)",
                 "6" = "IFng+IL2-TNFa+154- (2)",
                 "7" = "IFng+IL2-TNFa-154- (1)",
                 "8" = "IFng+IL2-TNFa-154- (1)"
                 )

other_features_cd8 <- c("1" = "EM",
                    "2" = "GzB+/EM",
                    "3" = "",
                    "6" = "GzB+/TD",
                    "7" = "IL17a+/ICOS/Naive",
                    "8" = "GzB+")

annotation_tab_cd8 <- bind_cols(c(1,2,3,6,7,8), key_functions_cd8, other_features_cd8)
names(annotation_tab_cd8) <- c("Cluster", "Key functions", "Other features")

# kable(annotation_tab, col.names = c("Cluster", "Key functions", "Other features"),
#       caption = "Cluster annotation, separated by key functions and other features. Numbers in parentheses indicate number of positive key functions out of 4.", escape = F, longtable = T, booktabs = T, linesep = "") %>%
#   kable_styling(font_size = 10, latex_options = c("hold_position", "repeat_header", "scale_down"))

 kable(annotation_tab_cd8,
    format = output_type, longtable = FALSE, booktabs = TRUE, linesep = "", escape = FALSE,
    caption = "Cluster annotation, separated by key functions and other features. Numbers in parentheses indicate number of positive key functions out of 4."
  ) %>%
  kable_styling(
    font_size = 6,
    latex_options = c("hold_position")
  )



## ----percent-of-stim-in-cluster-2-----------------------------------------------------------------------------------------------------------------------------------------


# Getting CD4+ counts by pubid, visitno, Stim
all_samps_2_of_Func_n_cd8 <- all_samps_2_of_TH1_CD154_IL17a_GzB_cd8  %>%
  left_join(., dplyr::select(ics_adata_cd8, sample,treat) %>% distinct(), by = "sample") %>% 
  dplyr::count(sample, treat,visitno, Stim, name = 'pos_2_of_Func')



cluster_tab_stim_cd8 <- bind_cols(
  all_samps_2_of_TH1_CD154_IL17a_GzB_cd8 %>%
    slice(sample_index)  %>%
  left_join(., dplyr::select(ics_adata_cd8, sample,treat) %>% distinct(), by = "sample"),
  data.frame(Cluster = kmeans_out_cd8$cluster))%>%
        filter(Cluster != 4 & Cluster != 5) %>% 

#     mutate(Cluster = case_when(Cluster %in% c(2, 4, 5, 7) ~ 2,
#                              Cluster == 1 ~ 1,
#                              Cluster == 3 ~ 3,
#                              Cluster == 6 ~ 4,
#                              Cluster == 8 ~ 5))
# ) %>%
  dplyr::count(treat,sample, visitno, Stim, Cluster, name = 'pos_cluster_2_of_Func') %>%
  # Filling missing clusters
  pivot_wider(names_from = 'Cluster',
              values_from = 'pos_cluster_2_of_Func',
              values_fill = 0) %>%
  pivot_longer(cols = -c(sample, treat,visitno, Stim), names_to = 'Cluster',
               values_to = 'pos_cluster_2_of_Func') %>%

  full_join(
    all_samps_2_of_Func_n_cd8,
    by = c("sample", "visitno", "Stim", "treat")
  ) %>%
  full_join(
    all_samps_n_counts_cd8 ,
    by = c("sample",  "Stim")
  ) %>%
  mutate(
    percent_of_cd8 = 100 * pos_cluster_2_of_Func / cd8_pos,
    percent_of_2_of_Func = 100 * pos_cluster_2_of_Func / pos_2_of_Func,
    pubid = paste0('ID',
                   ifelse(as.numeric(factor(sample)) < 10, '0',''),
                   as.numeric(factor(sample))),
    Stim = factor(Stim),
    Cluster = factor(as.numeric(Cluster))
  ) %>%
   group_by(sample, treat,visitno, Cluster) %>%
# getting BG subtracted
  group_by(sample, treat,visitno, Cluster) %>%
  mutate(
    percent_of_cd8_neg = mean(percent_of_cd8[Stim == 'negctrl']),
    percent_of_cd8_net = percent_of_cd8 - percent_of_cd8_neg,
    `.groups` = "drop"
  ) %>%  mutate(treat = factor(treat ,
                              levels = c("DPBS sucrose",
                                         "20 µg eOD-GT8 60mer + AS01B",
                                         "100 µg eOD-GT8 60mer + AS01B"), ordered = T)) %>%
  group_by(sample, Stim) %>%
  mutate(
total_percent_net = sum(pmax(0, percent_of_cd8_net)),
rel_percent = 100 * pmax(0, percent_of_cd8_net) / total_percent_net
) %>%
  filter(!is.na(treat))%>% 
   mutate(treat_short = factor(case_when(treat == "DPBS sucrose" ~ "Placebo",
                                 treat == "20 µg eOD-GT8 60mer + AS01B" ~ "20 µg",
                                 treat == "100 µg eOD-GT8 60mer + AS01B" ~ "100 µg"), 
                               levels = c("Placebo","20 µg","100 µg"), ordered = T))
  
  











## ----percent-of-cd8-figs--------------------------------------------------------------------------------------------------------------------------------------------------

#do we want these by group or treatment?
#answers for cluster 1, how many 2 of 6 cells are there for each stim out of your total CD4 positive cells for each stim?


    
    fig3b <- cluster_tab_stim_cd8 %>%
  filter(!is.na(Cluster)) %>% 
      filter(Stim %in% c("eOD-GT8", "LumSyn", "CMV")) %>%
      ggplot(aes(
        x = treat_short ,
        col = treat_short,
        y = pmax(0.001, percent_of_cd8_net)
      )) +
      geom_point(
        position = position_jitter(
          width = .1,
          height = 0,
          seed = 3241
        ),
        size = 2
      ) +
      geom_boxplot(
        fill = NA,
        lwd = .5,
        outlier.colour = NA,
        alpha = 1
      ) +
      facet_grid(Stim~ Cluster) +
      
      scale_x_discrete("") +
      scale_y_log10(
        paste0("% of CD8+ T cells"),
        breaks = c(.001, .01, .1, 1, 10, 100),
        labels = c(expression("" <= "0.001"),
                   '0.01', '0.1', "1", "10", "100")
      ) + theme_bw() +
      scale_color_manual(
        name = "",
        values = c("#D92321", "#1749FF", "#0AB7C9", "#FF6F1B"),
        limits = c(
          "Placebo",
          "20 µg",
          "100 µg"
        )
      ) + theme_bw()+
      theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
                text=element_text(size=18),
                axis.text.x = element_blank(),

        legend.position = "bottom",
        panel.spacing = unit(0, units = 'lines'))
      


pdf(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig3b.pdf',
    width = 10, height = 12);

fig3b
dev.off()


png(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig3b.png',
    res = 200,pointsize = 5,
    width = 1800, height = 2200);

fig3b

dev.off()



## ----rel-percent-of-cd8-pooled-bar,  fig.scap="Bar plot of relative percent positive cells, based on 2 of 6 function marker subset, out of total CD8+ T-cells, at day 70 by antigen and cluster.", fig.cap="Bar plot of relative percent positive cells, based on 2+ of 6 functional marker subset, out of total CD8+ T-cells, at day 70 by antigen and cluster. Placebo recipients are not included, and vaccine recipients are pooled. Error bars depicting 95\\% confidence intervals are displayed over each bar.", fig.height=7.25----
plot_data <- cluster_tab_stim_cd8 %>% 
  mutate(treat_pooled = case_when(treat != "DPBS sucrose" ~ "eOD-GT8 pooled",
                                    treat == "DPBS sucrose" ~ "DPBS sucrose")) %>% 
  filter(treat_pooled == "eOD-GT8 pooled" & Stim != "negctrl") %>% 
  filter(!(Stim %in% c("eOD-GT8-1","eOD-GT8-2","eOD-GT8-3"))) %>% 
  filter(rel_percent != "NaN")

summary_plot_dat_cd8 <- Rmisc::summarySE(plot_data, measurevar="rel_percent", groupvars=c("Stim","Cluster")) %>% 
  mutate(Cluster = as.integer(as.character(Cluster)) )%>% 
  full_join(annotation_tab_cd8, by = "Cluster") %>% 
  # mutate(Cluster = as.factor(Cluster)) %>% 
  group_by(Stim) %>% 
  # mutate(Cluster = as.factor(Cluster)) %>%
  mutate(annotation_plot = factor(
    paste0(Cluster, ": ", `Key functions`, "/", `Other features`),
    levels = c(
      "1: IFng+IL2+TNFa+154- (3)/EM",
                 "2: IFng+IL2-TNFa+154- (2)/GzB+/EM",
                 "3: IFng-IL2+TNFa+154- (2)/",
                 "6: IFng+IL2-TNFa+154- (2)/GzB+/TD",
                 "7: IFng+IL2-TNFa-154- (1)/IL17a+/ICOS/Naive",
      "8: IFng+IL2-TNFa-154- (1)/GzB+"
    )
  )) %>%
  mutate(annotation_plot_wrap = factor(
    paste0(Cluster, ": ", `Key functions`, "\n", `Other features`),
    levels = c(
      "1: IFng+IL2+TNFa+154- (3)\nEM",
                 "2: IFng+IL2-TNFa+154- (2)\nGzB+/EM",
                 "3: IFng-IL2+TNFa+154- (2)\n",
                 "6: IFng+IL2-TNFa+154- (2)\nGzB+/TD",
                 "7: IFng+IL2-TNFa-154- (1)\nIL17a+/ICOS/Naive",
      "8: IFng+IL2-TNFa-154- (1)\nGzB+"
    )
  ))


fig3c<- ggplot(data = summary_plot_dat_cd8, aes(fill = Stim,
        x = Stim ,
        col = Stim,
        y = pmax(1,rel_percent)
      )) +
          geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
    geom_errorbar(position=position_dodge(.9),colour = "black", width=.25, aes(ymin=pmax(1,rel_percent-ci), ymax=pmax(1,rel_percent+ci))) +
      facet_wrap(~ annotation_plot_wrap, nrow = 2) +
      scale_x_discrete("") +
      scale_y_continuous(
        paste0("Relative % of CD8+ T cells"),
        breaks = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
        labels = c(expression("" <= "1"),"10","20","30","40", "50", "60", "70", "80")
      ) + theme_bw() +
      scale_color_manual('', values = c("#0000FF","#1B9E77","#0AB7C9")) +
         scale_fill_manual('', values = c("#0000FF","#1B9E77","#0AB7C9")) +
   theme(
        text=element_text(size=18),
        legend.text = element_text(size=24),
        axis.text.x = element_text(
          size = 14
        ),
        legend.position = "none",
        panel.spacing = unit(0, units = 'lines'))
 
 
 


pdf(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig3c.pdf',
    width = 14, height = 10);

fig3c

dev.off()


png(file = 'Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig3c.png',
    res = 150,pointsize = 5,
    width = 2200, height = 1500);

fig3c

dev.off()


