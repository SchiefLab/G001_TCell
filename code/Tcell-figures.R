# This code reads in cvd725_ics_2022DEC27.csv, cvd725_tfh_ics_2022DEC23.csv, and
# cvd725_association_TandBcell_2022DEC28.csv and plots Figures 1B, 1C, 1D 
# and 4 in the T-cell manuscript.


### Add additional packages needed here 
### Only works for CRAN packages (manually write library statement for Bioconductor or GitHub packages)
packages = c("MASS","exact2x2","Exact","tidyverse","knitr","kableExtra", 'DataSpaceR',"VISCfunctions","cowplot", "data.table","GGally",
             "ggplot2",  "here")
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    # VISCfunctions must be installed in seperately (not in CRAN)
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

source(here('ICS', 'both_doses_PT_Report',"functions.R"))


tfh_ics <- read.csv("N:/cavd/studies/cvd725/pdata/cohen_et_al_manuscript/data/cvd725_tfh_ics_2022DEC23.csv")
ics <- read.csv("N:/cavd/studies/cvd725/pdata/cohen_et_al_manuscript/data/cvd725_ics_2022DEC27.csv")
TandBcell_assoc <- read.csv("N:/cavd/studies/cvd725/pdata/cohen_et_al_manuscript/data/cvd725_association_TandBcell_2022DEC28.csv")


## ----patient_info, message=FALSE, warning=TRUE, paged.print=TRUE----------------------------------------------------------------------------------------------------------
ppt_data <- ics %>% 
  group_by(visitno, treat) %>% 
  summarise(n = n_distinct(sample))


## ----response-rate-data---------------------------------------------------------------------------------------------------------------------------------------------------
response_ints <- ics %>% 
  group_by(treat,parent_plot,pop_grouping, population, antigen) %>% 
  summarise(npos_fisher = sum(response_fisher),
            npos_mimosa = sum(response_MIMOSA), 
            n = length(unique((sample))),
            fisher_lcl = binom::binom.confint(npos_fisher, n, methods = "wilson")$lower,
            fisher_rr = mean(response_fisher),
            fisher_ucl = binom::binom.confint(npos_fisher, n, methods = "wilson")$upper,
            mimosa_lcl = binom::binom.confint(npos_mimosa, n, methods = "wilson")$lower,
            mimosa_rr = mean(response_MIMOSA),
            mimosa_ucl = binom::binom.confint(npos_mimosa, n, methods = "wilson")$upper,
            cred_int2.5 = mean(Q2.5),
            #cred_int25 = mean(Q25),
            cred_int50 = mean(Q50),
           # cred_int75 = mean(Q75),
            cred_int97.5 = mean(Q97.5),
            all_zero_mimosa = mimosa_rr == 0)



## ----response-rate-testing, eval =F, message=FALSE------------------------------------------------------------------------------------------------------------------------
## 
## options(knitr.kable.NA = '', xtable.comment = FALSE)
## 
## 
## # rate test done using Barnard's test from the Exact package - will be incorporated into VISC functions later
## 
## treat_list <- c("DPBS sucrose","20 µg eOD-GT8 60mer + AS01B",
##                       "100 µg eOD-GT8 60mer + AS01B")
## 
## group_comparison_fun = function(study_in, response_var = "response_MIMOSA"){
##   study_data = copy(study_in) %>% arrange(treat)
##   study_groups = as.character(unique(study_data$treat))
##   if (any(is.na(study_groups))) stop("there is a NA group")
##   levels_here <- unique(study_data$treat)
##   setnames(study_data, response_var, "response_var")
##    map_df(head(study_groups, -1), function(i){
##     index_i = which(study_groups == i)
##     map_df(study_groups[-(1:index_i)], function(j){
##       test_data = subset(study_data, treat %in% c(i, j))
## 
##       if (dplyr::n_distinct(test_data$treat) != 2) stop("Something wrong with subsetting")
## 
##       comparison_here <- paste0(i, ' vs. ', j)
## 
##       test_data_summary = test_data %>% group_by(treat) %>%
##         dplyr::summarise(
##           n = n(),
##           response_frac = paste(sum(response_var), n(), sep = "/"),
##                   response_prop = eval(parse(text = response_frac)),
##                   response_pct = response_prop * 100,
##           npos = sum(response_var),
##             z = abs(qnorm(1 - (1 - .95 )/2)),
##             p = npos/n,
##             denom = 1 + z*z/n,
##             t1 = p + z*z/2/n,
##             t2 = z * sqrt(p*(1 - p)/n + z*z/4/n/n),
##             mean = p, lower = (t1 - t2)/denom, upper = (t1 + t2)/denom,
##           response_info = paste0(npos, '/', n, ' = ', round_away_0(response_prop * 100, 1), '%; (',
##                           round_away_0(lower * 100, 1), '%, ', round_away_0(upper * 100, 1),'%)'))
## 
##       response_props = test_data_summary$response_prop
## 
##       test_data$response_factor = factor(test_data$response_var, levels = c(0,1))
##       aa_table <- with(test_data, table(droplevels(treat), response_factor))
## 
##       if (nrow(aa_table) != 2 | ncol(aa_table) != 2) browser()
## 
##       if (all(response_props == 0) | all(response_props == 1)) exact_p = 1 else{
##           exact_p = exact.test(aa_table, method = 'Z-pooled', alternative = 'two.sided',
##                                to.plot = FALSE)$p.value
##         }
## 
## 
##       data.frame(
##         Comparison = comparison_here,
##         group1 = i,
##         group2 = j,
##         sample_sizes = paste0(test_data_summary$n,collapse = ' vs. '),
##         pos_responses = paste0(test_data_summary$npos,collapse = ' vs. '),
##         response_props = paste0(test_data_summary$response_prop,collapse = ' vs. '),
##         response_stats = paste0(test_data_summary$response_info,collapse = ' vs. '),
##         response_stats1 = test_data_summary$response_info[1],
##         response_stats2 = test_data_summary$response_info[2],
##         fisher_p = fisher.test(aa_table)$p.value,
##         exact_p = exact_p
##       )
##     })
##   })
## }
## 
## 
## rr_test <- subset(adata, !is.na(response_MIMOSA))  %>%
##   dplyr::group_by(antigen, parent_plot, pop_grouping,pop_table) %>%
##   do(data.frame(group_comparison_fun(., response_var = "response_MIMOSA"))) %>%
##   ungroup() %>%
##   mutate(p_tab = pretty_pvalues(exact_p, background = "yellow", bold = T, digits = 4)) %>%
##   arrange(Comparison,parent_plot, antigen) %>%
##   group_by(Comparison, parent_plot, antigen) %>%
##   mutate(parent_print = ifelse(duplicated(parent_plot), NA, as.character(parent_plot)),
##          # comp_print = ifelse(duplicated(parent_plot), NA, as.character(Comparison)),
##          stats_print = gsub("%", "\\\\%", response_stats),
##          pop_print =  gsub("IFN-g", "IFN-$\\\\gamma$", pop_table),
##          pop_print =  gsub("TNF-a", "TNF-$\\\\alpha$", pop_print),
## 
##          ant_print = ifelse(duplicated(parent_plot), NA, as.character(antigen)),
##          comp_print = ifelse(duplicated(parent_plot), NA,
##                              case_when(Comparison == "DPBS sucrose vs. 20 µg eOD-GT8 60mer + AS01B" ~ "DPBS sucrose vs. 20 µg",
##                                        Comparison == "DPBS sucrose vs. 100 µg eOD-GT8 60mer + AS01B" ~ "DPBS sucrose vs. 100 µg",
##                                        Comparison == "20 µg eOD-GT8 60mer + AS01B vs. 100 µg eOD-GT8 60mer + AS01B" ~ "20 µg vs. 100 µg"))) %>%
##   ungroup() %>%
##   dplyr::select( parent_print,ant_print, pop_print, comp_print, sample_sizes, pos_responses,
##                  stats_print, p_tab)
## 
## 
## 
##   kable(rr_test,caption = "Response rate testing between vaccine and placebo recipients and between vaccine groups two weeks after the second vaccination (day 70, visit 8). Testing was done using Barnard's test (two-sided, $\\alpha$ = 0.05) and p-values less than 0.05 are highlighted. Response rates (MIMOSA), proportions, and 95\\% Wilson confidence intervals (CI) are presented by T-cell subset, peptide pool, marker, and treatment.",
##         col.names = c("T-cell subset","Peptide Pool", "Marker", "Comparison", "N","Pos. Responses","Response Fraction = Response Rate; (95\\% Wilson score CI)", "P-value"), caption.short = "Response rate testing between vaccine and placebo recipients after the second vaccination.", format = output_type, escape = F, longtable = T, booktabs = T, linesep = "") %>%
##   kable_styling(font_size = 4.75, latex_options = c("hold_position", "repeat_header", "scale_down"))%>%
##       row_spec(which(!is.na(rr_test$ant_print))[-1] - 1, hline_after = TRUE)
##   # %>%
##   # row_spec(c(5,7,11,13,18,20,24,26,31,33,37,39,44,46,50,52,57,59,63), hline_after = T)
## 
## 


## ----rr-test-paired, eval = F ,results='asis', message = F, eval = F------------------------------------------------------------------------------------------------------
## options(knitr.kable.NA = '', xtable.comment = FALSE)
## 
## 
## rr_paired <- subset(adata,antigen %in% c("LumSyn","eOD-GT8") &
##                       treat !="DPBS sucrose" ) %>%
##   group_by(treat,parent_plot, pop_grouping,pop_table) %>%
##   dplyr::summarise(
##     comparison = paste(unique(antigen), collapse = " vs. "),
##     n = paste(c(length(response_MIMOSA[antigen == unique(antigen)[1]]),
##                 length(response_MIMOSA[antigen == unique(antigen)[2]])), collapse = " vs. "),
##     pos_responses = paste(c(sum(response_MIMOSA[antigen == unique(antigen)[1]]),
##                             sum(response_MIMOSA[antigen == unique(antigen)[2]])), collapse = " vs. "),
##     response_rates = paste0(round(((sum(response_MIMOSA[antigen == unique(antigen)[1]])/
##                                       length(response_MIMOSA[antigen == unique(antigen)[1]])))*100,1),"% vs. ",
##                             round(((sum(response_MIMOSA[antigen == unique(antigen)[2]])/
##                                       length(response_MIMOSA[antigen == unique(antigen)[2]])))*100,1), "%"),
##     p_value = if (all(response_MIMOSA == 0) | all(response_MIMOSA == 1)) 1 else
##       mcnemar_test_fun(sample, antigen, response_MIMOSA)
##   )   %>%
##   mutate(p_try = formatC(p_value, format = "e", digits = 2)) %>%
##   mutate(p_tab = pretty_pvalues(as.numeric(p_try),
##                                 background = "yellow", digits = 4, bold = T)) %>%
##   arrange(comparison, treat, parent_plot) %>%
##   group_by(treat, parent_plot) %>%
##   mutate(parent_print = ifelse(duplicated(parent_plot), NA, as.character(parent_plot)),
##          comp_print = ifelse(duplicated(parent_plot), NA, as.character(comparison)),
##          treat_print = ifelse(duplicated(parent_plot), NA, as.character(treat)),
##          rates_print = gsub("%", "\\\\%", response_rates),
##          pop_print =  gsub("IFN-g", "IFN-$\\\\gamma$", pop_table),
##          pop_print =  gsub("TNF-a", "TNF-$\\\\alpha$", pop_print),) %>%
##   ungroup() %>%
##   dplyr::select( parent_print, pop_print,treat_print, comp_print,n, pos_responses,
##                  rates_print, p_tab)
## 
## kable(rr_paired,caption = "Paired response rate testing between eOD-GT8 and LumSyn among vaccinees two weeks after the second vaccination (day 70, visit 8). Testing was done using McNemar's test for paired data (two-sided, $\\alpha$ = 0.05) and p-values less than 0.05 are highlighted. Response rates (MIMOSA) are presented by T-cell subset and marker.",
##       col.names = c("T-cell subset", "Marker","Treatment", "Comparison", "N","Pos. Responses","Response Rates", "P-value"), caption.short = "Paired response rate testing between eOD-GT8 and LumSyn.", format = output_type, escape = F, longtable = F, booktabs = T, linesep = "") %>%
##   kable_styling(font_size = 8, latex_options = c("hold_position", "repeat_header", "scale_down")) %>%
##   row_spec(which(!is.na(rr_paired$parent_print))[-1] - 1, hline_after = TRUE)
## 
## 
## 


## ----magnitude-testing-cxcr5, warning=FALSE-------------------------------------------------------------------------------------------------------------------------------

# Group magnitude testing
magnitude_results <-  ics %>%
  mutate(group_plot = factor(Treat_Short, 
         levels = c("Placebo", "20 µg", "100 µg"), ordered = T)) %>% 
  mutate(pop_plot = factor(case_when(population == "IL2+" ~ "CXCR5+ IL-2+",
                           population == "154+" ~ "CXCR5+ CD40L+"), 
         levels = c("CXCR5+ IL-2+","CXCR5+ CD40L+"), ordered = T)) %>% 
  filter(antigen != "Total Circulating") %>% 
  filter(assay == 'ICS')%>%
  group_by(pop_plot, antigen) %>% 
  group_modify(
    ~pairwise_test_cont(
      x =pmax(.$mag, 0), group = .$Treat_Short, method = 'wilcox', paired = FALSE, 
      alternative = 'two.sided', num_needed_for_test = 3, digits = 3, 
      verbose = FALSE
      ) %>% as.data.frame
    ) %>% 
   ungroup() %>% 
  mutate(
    MagnitudeTest = pretty_pvalues(
      MagnitudeTest, output_type = output_type, sig_alpha = .05,
      background = 'yellow',
      digits = 4
      )
    ) %>% 
  rename("Median (Range)" = Median_Min_Max, 'Mean (SD)' = Mean_SD) %>% 
  dplyr::select(-PerfectSeparation)

kable(magnitude_results,caption = "Response magnitude testing for eOD-GT8 and LumSyn two weeks after the second vaccination (day 70, visit 8). Testing was done using the Wilcoxon rank-sum test (two-sided, $\\alpha$ = 0.05) and p-values less than 0.05 are highlighted. ",
      col.names = c("Marker","Peptide Pool", "Comparison", "N","Median (Range)", "Mean (SD)", "P-value"), caption.short = "Response magnitude testing for eOD-GT8 and LumSyn.", format = output_type, escape = F, longtable = F, booktabs = T, linesep = "") %>% 
  kable_styling(font_size = 7, latex_options = c("hold_position", "repeat_header", "scale_down")) 




## ----fig-1, fig.height=7.6, eval =F---------------------------------------------------------------------------------------------------------------------------------------
## 
## 
## plot_data <- adata %>%
##   mutate(pop_print = factor(population,
##                             levels = c("IFN-g and/or IL-2",
##                                        "IFNg+",
##                                        "IL2+",
##                                        "TNFa+",
##                                        "154+",
##                                        "IL4+",
##                                        "IL17a+"),
##                             labels = c( bquote("IFN-"*gamma*" and"*'/'*"or IL-2"),
##                                         bquote("IFN-"*gamma*""),
##                                         "IL-2",
##                                         bquote("TNF-"*alpha*""),
##                                         "CD40L",
##                                         bquote("IL"*'-'*"4"),
##                                         bquote("IL"*'-'*"17a")), ordered = T),
##          parent_plot_print = factor(parent_plot,
##                                     levels = c("CD3+/CD4+","CD3+/CD8+"),
##                                     labels = c(bquote("CD3"*'+'*"CD4"*'+'),
##                                                bquote("CD3"*'+'*"CD8"*'+')))) %>%
##   filter(antigen %in% c("eOD-GT8","LumSyn")) %>%
##   filter(population %in% c("IFN-g and/or IL-2","IFNg+","IL2+","154+","TNFa+","IL17a+", "IL4+") & parent_plot == "CD3+/CD4+")
## 
## plot_data_rr <- plot_data %>%
##   group_by(treat_short, antigen,
##            parent_plot_print, pop_print) %>%
##   summarise(n = n_distinct(sample),
##             npos = sum(response_MIMOSA),
##             prop = paste0("frac(",npos, ",", n, ")"),
##             perc = round((npos/n) * 100, 0)) %>%
##   mutate(print = paste0(npos, "/", n, "\n", perc, "%")) %>%
##   unique() %>% arrange(treat_short)
## 
## 
##  fig1a <- ggplot(plot_data,
##                    aes(x = treat_short, y = pmax(0.001,percent_cell_net_trunc))) +
##      labs(title = "", tag = "A") +
## 
##   geom_jitter(data = plot_data %>% dplyr::filter(response_MIMOSA == 0),
##              size = 2.5, position = position_jitterdodge(jitter.width = .65, seed = 1),
##              col = cbPalette[1],
##              aes(shape = factor(response_MIMOSA)), show.legend = F) +
##   geom_jitter(data = plot_data %>% dplyr::filter(response_MIMOSA == 1),
##              size = 2.5, position = position_jitterdodge(jitter.width = .75, seed = 1),
##              aes(shape = factor(response_MIMOSA), col = factor(treat_short))) +
##   geom_boxplot(data = plot_data %>% dplyr::filter(response_MIMOSA == 1),
##                fill = NA, lwd = .5, outlier.colour = NA, position = position_dodge(width = .75),
##                aes(col = factor(treat_short))) +
##   scale_color_manual(name = "", values = c("#D92321","#1749FF","#0AB7C9","#FF6F1B"),
##                                           limits = c("DPBS Sucrose","20 µg","100 µg")) +
##   scale_shape_manual(name = "", breaks = c(0,1), labels = c("Non-responder", "Responder"),
##                      values = c(2, 19)) +
##   scale_x_discrete("") +
##     scale_y_log10("% of CD4+ T Cells", limits = c(0.001, 10),breaks = c(0.001,0.01,0.1, 1, 10),
##                   labels = c(expression("" <= 0.001),  "0.01","0.1", "1", "10")) +
##   facet_grid(antigen ~ pop_print, scales = "free_y", labeller = label_parsed) +
##   theme_bw() +
##   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
##         #strip.placement = "outside",
##                 text=element_text(size=18),
## 
## #         axis.text.y = element_text(size = 18),
## #          strip.text.x = element_text(size = 18),
## #         axis.text.x = element_text(size = 18),
## #         # strip.text.y = element_text(size = 5),
## #         # strip.background = element_rect(fill = "white"),
## # axis.text=element_text(size=18),
## # strip.text.y = element_text(size = 18),
##         legend.position = "none") +guides(col=FALSE)+
##    geom_text(data = plot_data_rr,
##   aes(label = print, y = 2, x = treat_short),
##   hjust = "center", vjust = 0, show.legend = F, size = 4.5, parse = F)
## 
## 
## 
## plot_data <- adata %>%
##   mutate(pop_print = factor(population,
##                             levels = c("IFN-g and/or IL-2",
##                                        "IFN-g and/or IL-2 and/or CD40L",
##                                        "IFNg+",
##                                        "IL2+",
##                                        "154+",
##                                        "TNFa+",
##                                        "IL4+",
##                                        "IL17a+"),
##                             labels = c( bquote("IFN-"*gamma*" and"*'/'*"or IL-2"),
##                                         bquote("IFN-"*gamma*"  and"*'/'*"or IL-2 and"*'/'*"or CD40L"),
##                                         bquote("IFN-"*gamma*""),
##                                         "IL-2",
##                                         "CD40L",
##                                         bquote("TNF-"*alpha*""),
##                                         bquote("IL"*'-'*"4"),
##                                         bquote("IL"*'-'*"17a")), ordered = T),
##          parent_plot_print = factor(parent_plot,
##                                     levels = c("CD3+/CD4+","CD3+/CD8+"),
##                                     labels = c(bquote("CD3"*'+'*"CD4"*'+'),
##                                                bquote("CD3"*'+'*"CD8"*'+')))) %>%
##   filter(antigen %in% c("eOD-GT8","LumSyn")) %>%
##   filter(population %in% c("IFN-g and/or IL-2","IFNg+","IL2+","TNFa+") & parent_plot == "CD3+/CD8+")
## 
## plot_data_rr <- plot_data %>%
##   group_by(treat_short, antigen,
##            parent_plot_print, pop_print) %>%
##   summarise(n = n_distinct(sample),
##             npos = sum(response_MIMOSA),
##             prop = paste0("frac(",npos, ",", n, ")"),
##             perc = round((npos/n) * 100, 0)) %>%
##   mutate(print = paste0(npos, "/", n, "\n", perc, "%")) %>%
##   unique() %>% arrange(treat_short)
## 
## 
## fig1c <- ggplot(plot_data,
##                    aes(x = treat_short, y = pmax(0.001,percent_cell_net_trunc))) +
##        labs(title = "", tag = "C") +
## 
##   geom_jitter(data = plot_data %>% dplyr::filter(response_MIMOSA == 0),
##              size = 2.5, position = position_jitterdodge(jitter.width = .65, seed = 1),
##              col = cbPalette[1],
##              aes(shape = factor(response_MIMOSA)), show.legend = F) +
##   geom_jitter(data = plot_data %>% dplyr::filter(response_MIMOSA == 1),
##              size = 2.5, position = position_jitterdodge(jitter.width = .75, seed = 1),
##              aes(shape = factor(response_MIMOSA), col = factor(treat_short))) +
##   geom_boxplot(data = plot_data %>% dplyr::filter(response_MIMOSA == 1),
##                fill = NA, lwd = .5, outlier.colour = NA, position = position_dodge(width = .75),
##                aes(col = factor(treat_short))) +
##   scale_color_manual(name = "", values = c("#D92321","#1749FF","#0AB7C9","#FF6F1B"),
##                      limits = c("DPBS Sucrose","20 µg","100 µg")) +
##   scale_shape_manual(name = "", breaks = c(0,1), labels = c("Non-responder", "Responder"),
##                      values = c(2, 19)) +
##   scale_x_discrete("") +
##     scale_y_log10("% of CD8+ T Cells", limits = c(0.001, 10),breaks = c(0.001,0.01,0.1, 1, 10),
##                   labels = c(expression("" <= 0.001),  "0.01","0.1", "1", "10")) +
##   facet_grid(antigen ~pop_print, scales = "free_y", labeller = label_parsed) +
##   theme_bw() +
##   theme(plot.margin = unit(c(-.5, 0, 0, 0), "cm"),
##         text=element_text(size=18),
##         #strip.placement = "outside",
## #         axis.text.y = element_text(size = 18),
## #          strip.text.x = element_text(size = 18),
##          axis.text.x = element_text(size = 15),
## #         # strip.text.y = element_text(size = 5),
## #         # strip.background = element_rect(fill = "white"),
## # axis.text=element_text(size=18),
## # strip.text.y = element_text(size = 18),
##         legend.position = "bottom") +guides(col=FALSE, shape = FALSE)+
##    geom_text(data = plot_data_rr,
##   aes(label = print, y = 1.5, x = treat_short),
##   hjust = "center", vjust = 0, show.legend = F, size = 4.5, parse = F)
## 
## 
## fig1b <- adata_bcell %>%
##   mutate(group_plot = factor(Treat_Short,
##          levels = c("Placebo", "20 µg", "100 µg"), ordered = T)) %>%
##   mutate(pop_plot = factor(case_when(population == "IL2+" ~ "CXCR5+ IL-2+",
##                            population == "154+" ~ "CXCR5+ CD40L+"),
##          levels = c("CXCR5+ IL-2+","CXCR5+ CD40L+"), ordered = T)) %>%
##   filter(antigen != "Total Circulating") %>%
##   filter(assay == 'ICS') %>%
## 
##   ggplot(aes(x = factor(group_plot),
##              # y = pmax(mag, 0.01),
##              y = pmax(mag, 0.001),
##              color = factor(group_plot))) +
##        labs(title = ""
##            # , tag = "B"
##             ) +
## 
##   geom_point(
##     position = position_jitter(height = 0, width = .25, seed = 3241), size = 2.5, shape = 1,
##   ) +
##   geom_boxplot(
##     fill = NA, lwd = .5, outlier.colour = NA
##   ) +
##   scale_color_manual(name = "", values = c("#D92321","#1749FF","#0AB7C9"),
##                      limits = c("Placebo","20 µg","100 µg")) +
##   scale_x_discrete("") +
##   # scale_y_continuous('Background-adjusted percent circulating-Tfh T cells expressing CD154 and IL2') +
##   scale_y_log10('% of CD4+ T Cells', limits = c(0.001, 0.1),
##                      breaks = c(.001,.003,0.01, .03, .1, .3, 1, 10, 100),
##                      labels = c(expression("" <= "0.001%"),"0.003%","0.001%", '0.03%',
##                                 '0.1%', '0.3%', '1%', '10%', '100%')
##   ) +theme_bw()+
##   facet_grid(antigen ~ pop_plot) +
##   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
##                 text=element_text(size=18),
## 
##         legend.position = "bottom") +guides(col=FALSE)
## #
## 
## 
## 
## #### editing y axis limits per kristen 10/11/22
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/ctfh001.pdf',
##     width = 8, height = 7);
## 
## fig1b
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/ctfh001.png',
##     res = 250,pointsize = 5,
##     width = 2000, height = 2000);
## 
## fig1b
## 
## dev.off()
## 
## ############
## 
## 
## 
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig1.pdf',
##     width = 19, height = 14);
## 
## cowplot::plot_grid(fig1a,
##                    cowplot::plot_grid(fig1b, fig1c, nrow = 1,rel_widths = c(3,5)), rel_heights = c(5, 4.5, 4.5),nrow = 2)
## 
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig1.png',
##     res = 250,pointsize = 5,
##     width = 4900, height = 4000);
## 
## cowplot::plot_grid(fig1a,
##                    cowplot::plot_grid(fig1b, fig1c, nrow = 1,rel_widths = c(3,5)), nrow = 2,
##                    rel_heights = c(5, 4.5, 4.5))
## 
## dev.off()
## 
## 


## ----fig-4a, eval =F------------------------------------------------------------------------------------------------------------------------------------------------------
## 
## plot_data <- adata %>%
##   mutate(pop_print = factor(population,
##                             levels = c("IFN-g and/or IL-2"),
##                             labels = c( bquote("IFN-"*gamma*" and"*'/'*"or IL-2")), ordered = T),
##          parent_plot_print = factor(parent_plot,
##                                     levels = c("CD3+/CD4+","CD3+/CD8+"),
##                                     labels = c(bquote("CD3"*'+'*"CD4"*'+'),
##                                                bquote("CD3"*'+'*"CD8"*'+')))) %>%
##   filter(population %in% c("IFN-g and/or IL-2") & parent_plot == "CD3+/CD4+") %>%
##   filter(antigen != "LumSyn")
## 
## plot_data_rr <- plot_data %>%
##   group_by(treat_short, antigen,
##            parent_plot_print, pop_print) %>%
##   summarise(n = n_distinct(sample),
##             npos = sum(response_MIMOSA),
##             prop = paste0("frac(",npos, ",", n, ")"),
##             perc = round((npos/n) * 100, 0)) %>%
##   mutate(print = paste0(npos, "/", n, "\n", perc, "%")) %>%
##   unique() %>% arrange(treat_short)
## 
## 
##  fig4a <- ggplot(plot_data,
##                    aes(x = treat_short, y = pmax(0.001,percent_cell_net_trunc))) +
##   geom_jitter(data = plot_data %>% dplyr::filter(response_MIMOSA == 0),
##              size = 2.5, position = position_jitterdodge(jitter.width = .65, seed = 1),
##              col = cbPalette[1],
##              aes(shape = factor(response_MIMOSA)), show.legend = F) +
##   geom_jitter(data = plot_data %>% dplyr::filter(response_MIMOSA == 1),
##              size = 2.5, position = position_jitterdodge(jitter.width = .75, seed = 1),
##              aes(shape = factor(response_MIMOSA), col = factor(treat_short))) +
##   geom_boxplot(data = plot_data %>% dplyr::filter(response_MIMOSA == 1),
##                fill = NA, lwd = .5, outlier.colour = NA, position = position_dodge(width = .75),
##                aes(col = factor(treat_short))) +
##   scale_color_manual(name = "", values = c("#D92321","#1749FF","#0AB7C9","#FF6F1B"),
##                                           limits = c("DPBS Sucrose","20 µg","100 µg")) +
##   scale_shape_manual(name = "", breaks = c(0,1), labels = c("Non-responder", "Responder"),
##                      values = c(2, 19)) +
##   scale_x_discrete("") +
##     scale_y_log10("% of CD4+ T Cells", limits = c(0.001, 10),breaks = c(0.001,0.01,0.1, 1, 10),
##                   labels = c(expression("" <= 0.001),  "0.01","0.1", "1", "10")) +
##   facet_grid(pop_print ~ antigen , scales = "free_y", labeller = label_parsed) +
##   theme_bw() +
##   theme(plot.margin = unit(c(1,1,1,1), "lines"),
##         #strip.placement = "outside",
##                 text=element_text(size=18),
## 
## #         axis.text.y = element_text(size = 18),
## #          strip.text.x = element_text(size = 18),
## #         axis.text.x = element_text(size = 18),
## #         # strip.text.y = element_text(size = 5),
## #         # strip.background = element_rect(fill = "white"),
## # axis.text=element_text(size=18),
## # strip.text.y = element_text(size = 18),
##         legend.position = "bottom") +guides(col=FALSE, shape = FALSE)+
##    geom_text(data = plot_data_rr,
##   aes(label = print, y = 3, x = treat_short),
##   hjust = "center", vjust = 0, show.legend = F, size = 4.5, parse = F)
## 
##  pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig4anew.pdf',
##     width = 11, height = 4);
## 
## fig4a
## dev.off()
## 
## 
##  png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig4anew.png',
##     res = 100,pointsize = 3,
##     width = 1100, height = 400);
## 
## fig4a
## dev.off()


## ----gc-tfh, eval =F------------------------------------------------------------------------------------------------------------------------------------------------------
## tfh_gc <- adata_bcell %>%
##   filter(assay == 'B cell') %>%
##   mutate(visit_plot = factor(Week_Visit, levels = c("Wk3 (V05)", "Wk11 (V09)"),
##                              labels = c("Post first-vaccination","Post second-vaccination"))) %>%
##   ggplot(aes(x = factor(Treat_Short,
##          levels = c("Placebo","20 µg","100 µg"), ordered = T),
##              y = pmax(mag, 0.01),
##              color = factor(Treat_Short,
##          levels = c("Placebo","20 µg","100 µg"), ordered = T))) +
##   geom_point(
##     position = position_jitter(height = 0, width = .25, seed = 3241), size = 2.5, shape = 1,
##   ) +
##   geom_boxplot(
##     fill = NA, lwd = .5, outlier.colour = NA
##   ) +
##  scale_color_manual(name = "", values = c("#D92321","#1749FF","#0AB7C9"),
##                      limits = c("Placebo","20 µg","100 µg")) +
##   scale_x_discrete("") +
##  # scale_y_continuous('Frequencies of % total GC Tfh from LN FNA') +
##   scale_y_log10('% GC Tfh of CD4+ T Cells',
##                 limits = c(0.01, 100),
##                 breaks = c(.01, .1, 1, 10, 100),
##                 labels = c(expression("" <= "0.01%"),
##                            '0.1%', '1%', '10%', '100%')
##   ) +
##   facet_wrap(. ~ visit_plot) + theme_bw()+
##   theme(                text=element_text(size=18),
## 
##     # axis.text.x = element_text(size = 7),
##     # axis.title.y = element_text(size = 9),
##     # axis.text.y = element_text(size = 7),
##     plot.margin = unit(c(.1,.1,.1,1), "lines"),
##     legend.position = "none"  )
## 
## pdf(file = '/home/cmahoney/Projects/Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/tfh_gc.pdf',
##     width = 7, height = 5.5);
## 
## tfh_gc
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/tfh_gc.png',
##     res = 100,pointsize = 3,
##     width = 700, height = 580);
## 
## tfh_gc
## dev.off()


## ----gc-tfh-magnitude-testing, warning=FALSE------------------------------------------------------------------------------------------------------------------------------

# Group magnitude testing
magnitude_results_gc_tfh <- adata_bcell %>%
  filter(assay == 'B cell') %>% 
  mutate(visit_plot = factor(Week_Visit, levels = c("Wk3 (V05)", "Wk11 (V09)"),
                             labels = c("Post first-vaccination","Post second-vaccination"))) %>%
  group_by(visit_plot) %>%
  group_modify(
    ~pairwise_test_cont(
      x = pmax(.$mag, 0), group = .$Treat_Short, method = 'wilcox', paired = FALSE, 
      alternative = 'two.sided', num_needed_for_test = 3, digits = 3, 
      verbose = FALSE
      ) %>% as.data.frame
    ) %>% 
   ungroup() %>% 
  mutate(
    MagnitudeTest = pretty_pvalues(
      MagnitudeTest, output_type = output_type, sig_alpha = .05,
      background = 'yellow',
      digits = 4
      )
    )  %>% 
  rename("Median (Range)" = Median_Min_Max, 'Mean (SD)' = Mean_SD) %>% 
  dplyr::select(-PerfectSeparation)

kable(magnitude_results_gc_tfh,caption = "Response magnitude testing for GC Tfh. Testing was done using the Wilcoxon rank-sum test (two-sided, $\\alpha$ = 0.05) and p-values less than 0.05 are highlighted. ",
      col.names = c("Visit", "Comparison", "N","Median (Range)", "Mean (SD)", "P-value"), caption.short = "Response magnitude testing for GC Thf.", format = output_type, escape = F, longtable = F, booktabs = T, linesep = "") %>% 
  kable_styling(font_size = 7, latex_options = c("hold_position", "repeat_header", "scale_down")) 




## ----gc-tfh-magnitude-testing-paired, warning=FALSE-----------------------------------------------------------------------------------------------------------------------

# Group magnitude testing
magnitude_results_gc_tfh_paired <- adata_bcell %>%
  filter(assay == 'B cell') %>% 
  mutate(visit_plot = factor(Week_Visit, levels = c("Wk3 (V05)", "Wk11 (V09)"),
                             labels = c("Post first-vaccination","Post second-vaccination"))) %>%
  arrange(pubid) %>% 
  group_by( Treat_Short) %>%
  group_modify(
    ~pairwise_test_cont(
      x = pmax(.$mag, 0), group = .$visit_plot, method = 'wilcox', paired = TRUE, 
      id = .$pubid,alternative = 'two.sided', num_needed_for_test = 3, digits = 3, 
      verbose = FALSE
      ) %>% as.data.frame
    ) %>% 
   ungroup() %>% 
  mutate(
    MagnitudeTest = pretty_pvalues(
      MagnitudeTest, output_type = output_type, sig_alpha = .05,
      background = 'yellow',
      digits = 4
      )
    )  %>% 
  rename("Median (Range)" = Median_Min_Max, 'Mean (SD)' = Mean_SD) %>% 
  dplyr::select(-PerfectSeparation)

kable(magnitude_results_gc_tfh_paired,caption = "Paired response magnitude testing for GC Tfh. Testing was done using the Wilcoxon signed-rank test for paired data (two-sided, $\\alpha$ = 0.05) and p-values less than 0.05 are highlighted. ",
      col.names = c("Treatment", "Comparison", "N","Median (Range)", "Mean (SD)", "P-value"), caption.short = "Paired response magnitude testing for GC Thf.", format = output_type, escape = F, longtable = F, booktabs = T, linesep = "") %>% 
  kable_styling(font_size = 7, latex_options = c("hold_position", "repeat_header", "scale_down")) 




## ----cor-tests, eval =F---------------------------------------------------------------------------------------------------------------------------------------------------
## 
## 
## corr_test_dat_low <- corr_data %>%
##           filter(Treat_Short != "Placebo" & !is.na(`B cell/Tot Tfh`) & Treat_Short == "20 µg") %>%
##           dplyr::select(Treat_Short, `B cell/Tot Tfh`, `B cell/GC`)
## 
## corr_test_dat_high <- corr_data %>%
##           filter(Treat_Short != "Placebo" & !is.na(`B cell/Tot Tfh`) & Treat_Short == "100 µg") %>%
##           dplyr::select(Treat_Short, `B cell/Tot Tfh`, `B cell/GC`)
## 
## corr_test_dat_both <- corr_data %>%
##           filter(Treat_Short != "Placebo" & !is.na(`B cell/Tot Tfh`) ) %>%
##           dplyr::select(Treat_Short, `B cell/Tot Tfh`, `B cell/GC`)
## 
## lowdose_cor_results <- data.frame(treat = "20 µg",
##            rho = round(stats::cor.test(corr_test_dat_low$`B cell/Tot Tfh`, corr_test_dat_low$`B cell/GC`, method = "spearman")$estimate, 3),
##            pval = pretty_pvalues(stats::cor.test(corr_test_dat_low$`B cell/Tot Tfh`, corr_test_dat_low$`B cell/GC`, method = "spearman")$p.value),
##            row.names = NULL)
## 
## highdose_cor_results <- data.frame(treat = "100 µg",
##            rho = round(stats::cor.test(corr_test_dat_high$`B cell/Tot Tfh`, corr_test_dat_high$`B cell/GC`, method = "spearman")$estimate, 3),
##            pval = pretty_pvalues(stats::cor.test(corr_test_dat_high$`B cell/Tot Tfh`, corr_test_dat_high$`B cell/GC`, method = "spearman")$p.value),
##            row.names = NULL)
## 
## overall_cor_results <- data.frame(treat = "Both",
##            rho = round(stats::cor.test(corr_test_dat_both$`B cell/Tot Tfh`, corr_test_dat_both$`B cell/GC`, method = "spearman")$estimate,3),
##            pval = pretty_pvalues(stats::cor.test(corr_test_dat_both$`B cell/Tot Tfh`, corr_test_dat_both$`B cell/GC`, method = "spearman")$p.value),
##            row.names = NULL)
## 
## cor_results = bind_rows(lowdose_cor_results, highdose_cor_results, overall_cor_results) %>%
##   mutate(cor_results = paste0("Corr: ", rho, " (p = ", pval, ")"),
##          cor_results_star = case_when(pval < 0.001 ~ paste0("Corr: ", rho, "***"),
##                                       pval >= 0.001 & pval < 0.01~ paste0("Corr: ", rho, "** "),
##                                       pval >= 0.01 & pval < 0.05 ~ paste0("Corr: ", rho, "*  ")))
## 
## 
## 


## ----cor-plot, fig.height=5, fig.width=6, fig.scap="Paper Figure 4 Panel C: Correlations", fig.cap= "Paper Figure 4 Panel C: Correlations, with scatterplots of magnitudes on the lower part, spearman correlation estimates on the upper part (overall and by treatment), and density plots on the diagonal.", warning=F, eval =F----
## 
## fig5_dat <-  corr_data %>%
##           filter(Treat_Short != "Placebo" & !is.na(`B cell/Tot Tfh`))
## 
## fig5 <- ggplot(data = fig5_dat) +
##   geom_jitter(aes(x = `B cell/Tot Tfh`, y = `B cell/GC`, color = factor(Treat_Short)),
##               size = 2, position = position_jitter(height = 0, width = 0),
##               show.legend = TRUE) +
##     scale_x_continuous('Tfh') +
##     scale_y_continuous("GC") +
##     theme_bw() +
##    scale_color_manual(name = "", values = c("#1749FF","#0AB7C9"),
##                      limits = c("20 µg","100 µg")) +
##     theme(legend.position = "bottom", legend.box = "vertical",
##          text = element_text(size = 18)) +
##    geom_text(data = filter(cor_results, treat == "Both"),
##             aes(label = cor_results_star, y = 45, x = 50),
##             hjust = 1, show.legend = F, size = 6, color = "black") +
##   geom_text(data = filter(cor_results, treat == "20 µg"),
##             aes(label = cor_results_star, y = 35, x = 50),
##             hjust = 1, show.legend = F, size = 6, color = "#1749FF") +
##   geom_text(data = filter(cor_results, treat == "100 µg"),
##             aes(label = cor_results_star, y = 25, x = 50),
##             hjust = 1, show.legend = F, size = 6, color = "#0AB7C9")
## 
## pdf(file = '/home/cmahoney/Projects/Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig5.pdf',
##     width = 7, height = 5.5);
## 
## fig5
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/Feinberg725Analysis/paper/Tcell_paper_figs/Tcell-figures/fig5.png',
##     res = 100,pointsize = 3,
##     width = 700, height = 580);
## 
## fig5
## dev.off()
## 


## ----cor-matrix, fig.height=7.25, fig.width=6, fig.scap="Paper Figure 4 Panel C: Correlations (no placebo)", fig.cap= "Paper Figure 4 Panel C: Correlations (no placebo), with scatterplots of magnitudes on the lower part, spearman correlation estimates on the upper part (overall and by treatment), and density plots on the diagonal.", warning=F, eval =F----
## 
## 
##   corr_plot_data <-  corr_data %>% filter(Treat_Short != 'Placebo') %>%
##   mutate('IFN-γ or IL-2+ \nVaccine CD4 T cell' = `ICS/IFNgorIL2/eOD-GT8` + `ICS/IFNgorIL2/LumSyn`) %>%
##           dplyr::select(Treat_Short,
##                 'GC Tfh' = 'B cell/Tot Tfh',
##                  'GC B cell' = 'B cell/GC',
##                  'eOD-GT8+ IgG \nGC B cells' = 'B cell/eOD V9',
##                  'VRC01+ IgG \nGC B cells' = 'B cell/VRC01 V09',
##                  'eOD-GT8+ IgG B cells' = 'B cell/eOD V8',
##                 'VRC01+ IgG B cells' = 'B cell/VRC01 V08',
##                  'IFN-γ or IL-2+ \neOD-GT8 CD4 T cell' = 'ICS/IFNgorIL2/eOD-GT8',
##                 #'bquote("IFN-"*gamma*" or IL-2") \neOD-GT8 CD4 T cell' = 'ICS/IFNgorIL2/eOD-GT8',
##                  'IFN-γ or IL-2+ \nLumSyn CD4 T cell' = 'ICS/IFNgorIL2/LumSyn',
##                 'IFN-γ or IL-2+ \nVaccine CD4 T cell')
## 
## 
##  corr_plot_data_simplified <-  corr_data %>% filter(Treat_Short != 'Placebo') %>%
##   mutate('IFN-γ or IL-2+ \nVaccine CD4 T cell' = `ICS/IFNgorIL2/eOD-GT8` + `ICS/IFNgorIL2/LumSyn`) %>%
##           dplyr::select(Treat_Short,
##                 'GC Tfh' = 'B cell/Tot Tfh',
##                 'Vaccine cTFH \n(CXCR5+ CD40L+)' = 'ICS/154+/TC',
##                 'Vaccine IFN-γ or IL-2 \nCD4 T cells'= 'IFN-γ or IL-2+ \nVaccine CD4 T cell',
## 
##                  'GC B cells' = 'B cell/GC',
##                  'eOD-GT8+ IgG \nGC B cells' = 'B cell/eOD V9',
##                  'eOD-GT8+ IgG B cells' = 'B cell/eOD V8')
## 
## p <- GGally::ggpairs(corr_plot_data,
## 
##         columns = 2:10, mapping = aes(color = Treat_Short),
##         lower = list(continuous =  wrap("points", size = .75)),
##         upper = list(continuous =  wrap("cor", method = "spearman", size = 3.5)),
##         xlab = 'Magnitude Values') +
##    # scale_color_manual(name = "", values = c("#1749FF","#0AB7C9"),
##    #                   limits = c("20 µg","100 µg")) +
##   theme_bw() +
##   theme(
##     text = element_text(size = 10))
## 
## 
## for(i in 1:p$nrow) {
##   for(j in 1:p$ncol){
##     p[i,j] <- p[i,j] +
##         scale_fill_manual(name = "", values = c("#1749FF","#0AB7C9"),
##                      limits = c("20 µg","100 µg")) +
##         scale_color_manual(name = "", values = c("#1749FF","#0AB7C9"),
##                      limits = c("20 µg","100 µg"))
##   }
## }
## 
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig6_july.pdf',
##     width = 12, height = 14);
## 
## p
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig6_july.png',
##     res = 100,pointsize = 3,
##     width = 1300, height = 1500);
## 
## p
## dev.off()
## 
## p <- GGally::ggpairs(corr_plot_data_simplified,
## 
##         columns = 2:7, mapping = aes(color = Treat_Short),
##         lower = list(continuous = wrap("points", size = .75)),
##         upper = list(continuous = wrap("cor", method = "spearman", size = 3.5)),
##         xlab = 'Magnitude Values') +
##    # scale_color_manual(name = "", values = c("#1749FF","#0AB7C9"),
##    #                   limits = c("20 µg","100 µg")) +
##   theme_bw() +
##   theme(
##     text = element_text(size = 10))
## 
## 
## for(i in 1:p$nrow) {
##   for(j in 1:p$ncol){
##     p[i,j] <- p[i,j] +
##         scale_fill_manual(name = "", values = c("#1749FF","#0AB7C9"),
##                      limits = c("20 µg","100 µg")) +
##         scale_color_manual(name = "", values = c("#1749FF","#0AB7C9"),
##                      limits = c("20 µg","100 µg"))
##   }
## }
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig6_simplified.pdf',
##     width = 12, height = 14);
## 
## p
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/fig6_simplified.png',
##     res = 100,pointsize = 3,
##     width = 800, height = 1000);
## 
## p
## dev.off()
## 
## 
## 


## ----corr-heatmap, eval =F------------------------------------------------------------------------------------------------------------------------------------------------
## 
## cors <- function(df) {
##   # turn all three matrices (r, n, and P into a data frame)
##   M <- Hmisc::rcorr(as.matrix(df), type="spearman")
##   # return the three data frames in a list return(Mdf)
##   Mdf <- map(M, ~data.frame(.x))
## }
## 
## formatted_cors <- function(df){
##   cors(df) %>%
##     map(~rownames_to_column(.x, var="measure1")) %>%
##     map(~pivot_longer(.x, -measure1, "measure2")) %>%
##     bind_rows(.id = "id") %>%
##     pivot_wider(names_from = id, values_from = value) %>%
##     mutate(sig_p = ifelse(P < .05, T, F),
##            p_if_sig = ifelse(P <.05, P, NA),
##            r_if_sig = ifelse(P <.05, r, NA))
## }
## 
## 
## ## correlation heatmap of combined dose groups
## 
##  corr_heat_data <-  corr_data %>% filter(Treat_Short != 'Placebo') %>%
##   mutate('IFN-γ or IL-2+ \nVaccine CD4 T cell' = `ICS/IFNgorIL2/eOD-GT8` + `ICS/IFNgorIL2/LumSyn`) %>%
##           dplyr::select(Treat_Short,
##                 `GC Tfh` = `B cell/Tot Tfh`,
##                 'Vaccine cTFH \n(CXCR5+ CD40L+)' = `ICS/154+/TC`,
##                 'Vaccine IFN-γ or IL-2 \nCD4 T cells'= 'IFN-γ or IL-2+ \nVaccine CD4 T cell',
## 
##                  `GC B cells` = `B cell/GC`,
##                  `eOD-GT8+ IgG \nGC B cells` = `B cell/eOD V9`,
##                  `eOD-GT8+ IgG \nB cells` = `B cell/eOD V8`)
## 
##  corr_heat_data_bothgroups <- corr_heat_data %>%
##    dplyr::select(-Treat_Short)
## 
##  corr_heat_data_lowdose <- corr_heat_data %>%
##    filter(Treat_Short == "20 µg") %>%
##    dplyr::select(-Treat_Short)
## 
##   corr_heat_data_highdose <- corr_heat_data %>%
##    filter(Treat_Short == "100 µg") %>%
##    dplyr::select(-Treat_Short)
## 
## 
##   ###both groups
## 
## corr_dat <- formatted_cors(corr_heat_data_bothgroups) %>%
##   mutate(measure2 = case_when(
##                 measure2 ==  "GC.Tfh"~"GC Tfh" ,
##                 measure2 ==   "Vaccine.cTFH...CXCR5..CD40L.."~"Vaccine cTFH \n(CXCR5+ CD40L+)",
##                 measure2 == "Vaccine.IFN.γ.or.IL.2..CD4.T.cells"~"Vaccine IFN-γ or IL-2 \nCD4 T cells",
## 
##                  measure2 ==   "GC.B.cells" ~ "GC B cells",
##                  measure2 ==   "eOD.GT8..IgG..GC.B.cells" ~ "eOD-GT8+ IgG \nGC B cells",
##                  measure2 ==   "eOD.GT8..IgG..B.cells" ~ "eOD-GT8+ IgG \nB cells")) %>%
##   filter(measure1!=measure2) %>%
##   mutate(
##     sig_label = ifelse(p_if_sig < 0.001, "***",
##                        ifelse(p_if_sig < 0.01, "**",
##                               ifelse(p_if_sig < 0.05, "*", ""))),
##     measure1 = factor(measure1),
##     measure2 = factor(measure2)) %>%
##   arrange(measure1)
## 
## # # drop lower estimates to create an 'upper' heatmap
## # for(i in unique(corr_dat$measure1)) {
## #   if(i %in% corr_dat$measure1){
## #     start <- max(which(corr_dat$measure1==i))
## #     corr_dat <- corr_dat[-which(corr_dat$measure2==i)[which(corr_dat$measure2==i) > start], ]
## #   }
## # }
## #write.csv(corr_dat, "../adata/corr_dat_version8.csv", row.names=FALSE)
## # BB: Sort measure1 and measure2 correctly to address reviewer response
## # First drop duplicated pairs
## uniq_measures <- unique(corr_dat$measure1) %>% as.character()
## corr_dat <- cbind(combn(uniq_measures, 2), rbind(uniq_measures, uniq_measures)) %>%
##   t() %>%
##   data.frame() %>%
##   rename(measure1 = uniq_measures,
##          measure2 = uniq_measures.1) %>%
##   left_join(corr_dat, by=c('measure1', 'measure2')) %>%
##   mutate(measure2_temp = measure2)
## corr_dat <- bind_rows(corr_dat %>% filter(measure1 == "eOD-GT8+ IgG \nB cells"), # This ensures correct format for manuscript figure: top left cor matrix
##                       corr_dat %>% filter(measure1 != "eOD-GT8+ IgG \nB cells") %>%
##                         mutate(measure2 = measure1,
##                                measure1 = measure2_temp)) %>%
##   mutate(
##     measure1 = fct_relevel(measure1,
##                            c("eOD-GT8+ IgG \nB cells", "Vaccine IFN-γ or IL-2 \nCD4 T cells",
##                              "Vaccine cTFH \n(CXCR5+ CD40L+)", "GC Tfh", "GC B cells",
##                              "eOD-GT8+ IgG \nGC B cells")),
##     measure2 = fct_relevel(measure2,
##                            c("eOD-GT8+ IgG \nB cells", "Vaccine IFN-γ or IL-2 \nCD4 T cells","Vaccine cTFH \n(CXCR5+ CD40L+)",
##                              "GC Tfh", "GC B cells", "eOD-GT8+ IgG \nGC B cells"))) %>%
##   arrange(measure1, measure2)
## #pdf('../figures/cov001_corr_heatmap_version8.pdf', width=10, height=10)
## stars <- ggplot(corr_dat, aes(measure1, measure2, col=r, label=sig_label)) +
##   geom_tile(col="black", fill=NA) +
##   geom_point(aes(size = abs(r)), shape=15) +
##   geom_text(color="black", size=5) +
##   labs(x = NULL, y = NULL, col = "Correlation") +
##   theme_minimal() +
##   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
##   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
##                          limits=c(-1,1), breaks=c(-1,0,1),
##                          labels=c("Negative","","Positive")) +
##   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
##   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
##   scale_size(range=c(1,8), guide=FALSE) +
##   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
##         plot.title = element_text(size=20, hjust = 0.5),
##         plot.caption = element_text(size=10, hjust = 0.5),
##         legend.text=element_text(size=12),
##         #  axis.text = element_text(size=16,hjust=0),
##         axis.title = element_text(size=20),
##         strip.text=element_text(size=11),
##         legend.position = "bottom",
##         legend.title=element_text(size=11),
##         strip.background = element_blank(),
##         panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(),
##         axis.text.x = element_text(angle = 90,size=11,hjust=0.5,vjust=1,color="black"),
##         axis.text.y = element_text(size=11,vjust=0.5,color="black"),
##         axis.ticks = element_blank(),
##         axis.line = element_blank())
## #dev.off()
## 
## 
## #pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_both.pdf',
## pdf(file = here("paper/Tcell_paper_figs/Tcell-figures/star_heatmap_both_UPDATED.pdf"),
##     width = 5, height = 5);
## 
## stars
## dev.off()
## 
## 
## #png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_both.png',
## png(file = here("paper/Tcell_paper_figs/Tcell-figures/star_heatmap_both_UPDATED.png"),
##     res = 100,pointsize = 3,
##     width = 800, height = 800);
## 
## stars
## dev.off()
## #
## #
## # #pdf('../figures/cov001_corr_heatmap_version8_n.pdf', width=10, height=10)
## # ggplot(corr_dat, aes(measure1, measure2, col=r, label=n)) +
## #   geom_tile(col="black", fill="white") +
## #   geom_point(aes(size = abs(r)), shape=15) +
## #   geom_text(color="black", size=4) +
## #   labs(x = NULL, y = NULL, col = "Correlation") +
## #   theme_minimal() +
## #   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
## #   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
## #                          limits=c(-1,1), breaks=c(-1,0,1),
## #                          labels=c("Negative","","Positive")) +
## #   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
## #   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
## #   scale_size(range=c(1,8), guide=FALSE) +
## #   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
## #         plot.title = element_text(size=20, hjust = 0.5),
## #         plot.caption = element_text(size=10, hjust = 0.5),
## #         legend.text=element_text(size=18),
## #         axis.text = element_text(size=16,hjust=0),
## #         axis.title = element_text(size=24),
## #         strip.text=element_text(size=16),
## #         legend.position = "bottom",
## #         legend.title=element_text(size=18),
## #         strip.background = element_blank(),
## #         panel.grid.major = element_blank(),
## #         panel.grid.minor = element_blank(),
## #         axis.text.x = element_text(angle = 90,size=16,hjust=0,color="black"),
## #         axis.text.y = element_text(size=16,vjust=0,color="black"),
## #         axis.ticks = element_blank(),
## #         axis.line = element_blank())
## # #dev.off()
## 
## ## Remove leading 0 from r:
## #corr_dat$rs <- sub("^0+", "", round(corr_dat$r,1))
## corr_dat$r1 <- round(corr_dat$r,2)
## 
## #pdf('../figures/cov001_corr_heatmap_version8_r.pdf', width=10, height=10)
## est <- ggplot(corr_dat, aes(measure1, measure2, col=r, label=r1)) +
##   geom_tile(col="black", fill="white") +
##   geom_point(aes(size = 2*abs(r)), shape=15) +
##   geom_text(color="black", size=4) +
##   labs(x = NULL, y = NULL, col = "Correlation") +
##   theme_minimal() +
##   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
##   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
##                          limits=c(-1,1), breaks=c(-1,0,1),
##                          labels=c("Negative","","Positive")) +
##   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
##   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
##   scale_size(range=c(1,8), guide=FALSE) +
##   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
##         plot.title = element_text(size=20, hjust = 0.5),
##         plot.caption = element_text(size=10, hjust = 0.5),
##         legend.text=element_text(size=18),
##         axis.text = element_text(size=16,hjust=0),
##         axis.title = element_text(size=24),
##         strip.text=element_text(size=16),
##         legend.position = "bottom",
##         legend.title=element_text(size=18),
##         strip.background = element_blank(),
##         panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(),
##         axis.text.x = element_text(angle = 90,size=16,hjust=0,color="black"),
##         axis.text.y = element_text(size=16,vjust=0,color="black"),
##         axis.ticks = element_blank(),
##         axis.line = element_blank())
## #dev.off()
## 
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/est_heatmap_both.pdf',
##     width = 10, height = 10);
## 
## est
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/pval_heatmap_both.png',
##     res = 100,pointsize = 3,
##     width = 800, height = 800);
## 
## est
## dev.off()


## ----corr-heat-low, eval =F-----------------------------------------------------------------------------------------------------------------------------------------------
## ## correlation heatmap of low dose group
## 
## 
##   ###low dose
## 
## 
## 
## corr_dat <- formatted_cors(corr_heat_data_lowdose) %>%
##   mutate(measure2 = case_when(
##                 measure2 ==  "GC.Tfh"~"GC Tfh" ,
##                 measure2 ==   "Vaccine.cTFH...CXCR5..CD40L.."~"Vaccine cTFH \n(CXCR5+ CD40L+)",
##                 measure2 == "Vaccine.IFN.γ.or.IL.2..CD4.T.cells"~"Vaccine IFN-γ or IL-2 \nCD4 T cells",
## 
##                  measure2 ==   "GC.B.cells" ~ "GC B cells",
##                  measure2 ==   "eOD.GT8..IgG..GC.B.cells" ~ "eOD-GT8+ IgG \nGC B cells",
##                  measure2 ==   "eOD.GT8..IgG..B.cells" ~ "eOD-GT8+ IgG \nB cells")) %>%
##   filter(measure1!=measure2) %>%
##   mutate(
##     sig_label = ifelse(p_if_sig < 0.001, "***",
##                        ifelse(p_if_sig < 0.01, "**",
##                               ifelse(p_if_sig < 0.05, "*", ""))),
##     measure1 = factor(measure1),
##     measure2 = factor(measure2)) %>%
##   arrange(measure1)
## 
## # # drop lower estimates to create an 'upper' heatmap
## # for(i in unique(corr_dat$measure1)) {
## #   if(i %in% corr_dat$measure1){
## #     start <- max(which(corr_dat$measure1==i))
## #     corr_dat <- corr_dat[-which(corr_dat$measure2==i)[which(corr_dat$measure2==i) > start], ]
## #   }
## # }
## #write.csv(corr_dat, "../adata/corr_dat_version8.csv", row.names=FALSE)
## # BB: Sort measure1 and measure2 correctly to address reviewer response
## # First drop duplicated pairs
## uniq_measures <- unique(corr_dat$measure1) %>% as.character()
## corr_dat <- cbind(combn(uniq_measures, 2), rbind(uniq_measures, uniq_measures)) %>%
##   t() %>%
##   data.frame() %>%
##   rename(measure1 = uniq_measures,
##          measure2 = uniq_measures.1) %>%
##   left_join(corr_dat, by=c('measure1', 'measure2')) %>%
##   mutate(measure2_temp = measure2)
## corr_dat <- bind_rows(corr_dat %>% filter(measure1 == "eOD-GT8+ IgG \nB cells"), # This ensures correct format for manuscript figure: top left cor matrix
##                       corr_dat %>% filter(measure1 != "eOD-GT8+ IgG \nB cells") %>%
##                         mutate(measure2 = measure1,
##                                measure1 = measure2_temp)) %>%
##   mutate(
##     measure1 = fct_relevel(measure1,
##                            c("eOD-GT8+ IgG \nB cells", "Vaccine IFN-γ or IL-2 \nCD4 T cells",
##                              "Vaccine cTFH \n(CXCR5+ CD40L+)", "GC Tfh", "GC B cells",
##                              "eOD-GT8+ IgG \nGC B cells")),
##     measure2 = fct_relevel(measure2,
##                            c("eOD-GT8+ IgG \nB cells", "Vaccine IFN-γ or IL-2 \nCD4 T cells","Vaccine cTFH \n(CXCR5+ CD40L+)",
##                              "GC Tfh", "GC B cells", "eOD-GT8+ IgG \nGC B cells"))) %>%
##   arrange(measure1, measure2)
## #pdf('../figures/cov001_corr_heatmap_version8.pdf', width=10, height=10)
## stars <- ggplot(corr_dat, aes(measure1, measure2, col=r, label=sig_label)) +
##   geom_tile(col="black", fill=NA) +
##   geom_point(aes(size = abs(r)), shape=15) +
##   geom_text(color="black", size=5) +
##   labs(x = NULL, y = NULL, col = "Correlation") +
##   theme_minimal() +
##   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
##   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
##                          limits=c(-1,1), breaks=c(-1,0,1),
##                          labels=c("Negative","","Positive")) +
##   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
##   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
##   scale_size(range=c(1,8), guide=FALSE) +
##  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
##         plot.title = element_text(size=20, hjust = 0.5),
##         plot.caption = element_text(size=10, hjust = 0.5),
##         legend.text=element_text(size=12),
##         #  axis.text = element_text(size=16,hjust=0),
##         axis.title = element_text(size=20),
##         strip.text=element_text(size=11),
##         legend.position = "bottom",
##         legend.title=element_text(size=11),
##         strip.background = element_blank(),
##         panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(),
##          axis.text.x = element_text(angle = 90,size=11,hjust=0.5,vjust=1,color="black"),
##         axis.text.y = element_text(size=11,vjust=0.5,color="black"),
##         axis.ticks = element_blank(),
##         axis.line = element_blank())
## #dev.off()
## 
## 
## #pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_low.pdf',
## pdf(file = here("/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_low_UPDATED.pdf"),
##     width = 5, height = 5);
## 
## stars
## dev.off()
## 
## 
## #png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_low.png',
## png(file = here("/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_low_UPDATED.png"),
##     res = 100,pointsize = 3,
##     width = 800, height = 800);
## 
## stars
## dev.off()
## #
## #
## # #pdf('../figures/cov001_corr_heatmap_version8_n.pdf', width=10, height=10)
## # ggplot(corr_dat, aes(measure1, measure2, col=r, label=n)) +
## #   geom_tile(col="black", fill="white") +
## #   geom_point(aes(size = abs(r)), shape=15) +
## #   geom_text(color="black", size=4) +
## #   labs(x = NULL, y = NULL, col = "Correlation") +
## #   theme_minimal() +
## #   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
## #   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
## #                          limits=c(-1,1), breaks=c(-1,0,1),
## #                          labels=c("Negative","","Positive")) +
## #   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
## #   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
## #   scale_size(range=c(1,8), guide=FALSE) +
## #   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
## #         plot.title = element_text(size=20, hjust = 0.5),
## #         plot.caption = element_text(size=10, hjust = 0.5),
## #         legend.text=element_text(size=18),
## #         axis.text = element_text(size=16,hjust=0),
## #         axis.title = element_text(size=24),
## #         strip.text=element_text(size=16),
## #         legend.position = "bottom",
## #         legend.title=element_text(size=18),
## #         strip.background = element_blank(),
## #         panel.grid.major = element_blank(),
## #         panel.grid.minor = element_blank(),
## #         axis.text.x = element_text(angle = 90,size=16,hjust=0,color="black"),
## #         axis.text.y = element_text(size=16,vjust=0,color="black"),
## #         axis.ticks = element_blank(),
## #         axis.line = element_blank())
## # #dev.off()
## 
## ## Remove leading 0 from r:
## #corr_dat$rs <- sub("^0+", "", round(corr_dat$r,1))
## corr_dat$r1 <- round(corr_dat$r,2)
## 
## #pdf('../figures/cov001_corr_heatmap_version8_r.pdf', width=10, height=10)
## est <- ggplot(corr_dat, aes(measure1, measure2, col=r, label=r1)) +
##   geom_tile(col="black", fill="white") +
##   geom_point(aes(size = abs(r)), shape=15) +
##   geom_text(color="black", size=4) +
##   labs(x = NULL, y = NULL, col = "Correlation") +
##   theme_minimal() +
##   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
##   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
##                          limits=c(-1,1), breaks=c(-1,0,1),
##                          labels=c("Negative","","Positive")) +
##   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
##   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
##   scale_size(range=c(1,8), guide=FALSE) +
##   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
##         plot.title = element_text(size=20, hjust = 0.5),
##         plot.caption = element_text(size=10, hjust = 0.5),
##         legend.text=element_text(size=18),
##         axis.text = element_text(size=16,hjust=0),
##         axis.title = element_text(size=24),
##         strip.text=element_text(size=16),
##         legend.position = "bottom",
##         legend.title=element_text(size=18),
##         strip.background = element_blank(),
##         panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(),
##         axis.text.x = element_text(angle = 90,size=16,hjust=0,color="black"),
##         axis.text.y = element_text(size=16,vjust=0,color="black"),
##         axis.ticks = element_blank(),
##         axis.line = element_blank())
## #dev.off()
## 
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/est_heatmap_low.pdf',
##     width = 8, height = 8);
## 
## est
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/pval_heatmap_low.png',
##     res = 100,pointsize = 3,
##     width = 800, height = 800);
## 
## est
## dev.off()
## 


## ----corr-heat-high, eval =F----------------------------------------------------------------------------------------------------------------------------------------------
## 
## 
## 
## 
##   ###high dose
## 
## corr_dat <- formatted_cors(corr_heat_data_highdose) %>%
##   mutate(measure2 = case_when(
##                 measure2 ==  "GC.Tfh"~"GC Tfh" ,
##                 measure2 ==   "Vaccine.cTFH...CXCR5..CD40L.."~"Vaccine cTFH \n(CXCR5+ CD40L+)",
##                 measure2 == "Vaccine.IFN.γ.or.IL.2..CD4.T.cells"~"Vaccine IFN-γ or IL-2 \nCD4 T cells",
## 
##                  measure2 ==   "GC.B.cells" ~ "GC B cells",
##                  measure2 ==   "eOD.GT8..IgG..GC.B.cells" ~ "eOD-GT8+ IgG \nGC B cells",
##                  measure2 ==   "eOD.GT8..IgG..B.cells" ~ "eOD-GT8+ IgG \nB cells")) %>%
##   filter(measure1!=measure2) %>%
##   mutate(
##     sig_label = ifelse(p_if_sig < 0.001, "***",
##                        ifelse(p_if_sig < 0.01, "**",
##                               ifelse(p_if_sig < 0.05, "*", ""))),
##     measure1 = factor(measure1),
##     measure2 = factor(measure2)) %>%
##   arrange(measure1)
## 
## # # drop lower estimates to create an 'upper' heatmap
## # for(i in unique(corr_dat$measure1)) {
## #   if(i %in% corr_dat$measure1){
## #     start <- max(which(corr_dat$measure1==i))
## #     corr_dat <- corr_dat[-which(corr_dat$measure2==i)[which(corr_dat$measure2==i) > start], ]
## #   }
## # }
## #write.csv(corr_dat, "../adata/corr_dat_version8.csv", row.names=FALSE)
## # BB: Sort measure1 and measure2 correctly to address reviewer response
## # First drop duplicated pairs
## uniq_measures <- unique(corr_dat$measure1) %>% as.character()
## corr_dat <- cbind(combn(uniq_measures, 2), rbind(uniq_measures, uniq_measures)) %>%
##   t() %>%
##   data.frame() %>%
##   rename(measure1 = uniq_measures,
##          measure2 = uniq_measures.1) %>%
##   left_join(corr_dat, by=c('measure1', 'measure2')) %>%
##   mutate(measure2_temp = measure2)
## corr_dat <- bind_rows(corr_dat %>% filter(measure1 == "eOD-GT8+ IgG \nB cells"), # This ensures correct format for manuscript figure: top left cor matrix
##                       corr_dat %>% filter(measure1 != "eOD-GT8+ IgG \nB cells") %>%
##                         mutate(measure2 = measure1,
##                                measure1 = measure2_temp)) %>%
##   mutate(
##     measure1 = fct_relevel(measure1,
##                            c("eOD-GT8+ IgG \nB cells", "Vaccine IFN-γ or IL-2 \nCD4 T cells",
##                              "Vaccine cTFH \n(CXCR5+ CD40L+)", "GC Tfh", "GC B cells",
##                              "eOD-GT8+ IgG \nGC B cells")),
##     measure2 = fct_relevel(measure2,
##                            c("eOD-GT8+ IgG \nB cells", "Vaccine IFN-γ or IL-2 \nCD4 T cells","Vaccine cTFH \n(CXCR5+ CD40L+)",
##                              "GC Tfh", "GC B cells", "eOD-GT8+ IgG \nGC B cells"))) %>%
##   arrange(measure1, measure2)
## #pdf('../figures/cov001_corr_heatmap_version8.pdf', width=10, height=10)
## stars <- ggplot(corr_dat, aes(measure1, measure2, col=r, label=sig_label)) +
##   geom_tile(col="black", fill=NA) +
##   geom_point(aes(size = abs(r)), shape=15) +
##   geom_text(color="black", size=5) +
##   labs(x = NULL, y = NULL, col = "Correlation") +
##   theme_minimal() +
##   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
##   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
##                          limits=c(-1,1), breaks=c(-1,0,1),
##                          labels=c("Negative","","Positive")) +
##   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
##   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
##   scale_size(range=c(1,8), guide=FALSE) +
##   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
##         plot.title = element_text(size=20, hjust = 0.5),
##         plot.caption = element_text(size=10, hjust = 0.5),
##         legend.text=element_text(size=12),
##         #  axis.text = element_text(size=16,hjust=0),
##         axis.title = element_text(size=20),
##         strip.text=element_text(size=11),
##         legend.position = "bottom",
##         legend.title=element_text(size=11),
##         strip.background = element_blank(),
##         panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(),
##         axis.text.x = element_text(angle = 90,size=11,hjust=0.5,vjust=1,color="black"),
##         axis.text.y = element_text(size=11,vjust=0.5,color="black"),
##         axis.ticks = element_blank(),
##         axis.line = element_blank())
## #dev.off()
## 
## 
## #pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_high.pdf',
## pdf(file = here("/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_high_UPDATED.pdf"),
##     width = 5, height = 5);
## 
## stars
## dev.off()
## 
## 
## #png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_high.png',
## png(file = here("/paper/Tcell_paper_figs/Tcell-figures/star_heatmap_high_UPDATED.png"),
##     res = 100,pointsize = 3,
##     width = 800, height = 800);
## 
## stars
## dev.off()
## #
## #
## # #pdf('../figures/cov001_corr_heatmap_version8_n.pdf', width=10, height=10)
## # ggplot(corr_dat, aes(measure1, measure2, col=r, label=n)) +
## #   geom_tile(col="black", fill="white") +
## #   geom_point(aes(size = abs(r)), shape=15) +
## #   geom_text(color="black", size=4) +
## #   labs(x = NULL, y = NULL, col = "Correlation") +
## #   theme_minimal() +
## #   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
## #   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
## #                          limits=c(-1,1), breaks=c(-1,0,1),
## #                          labels=c("Negative","","Positive")) +
## #   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
## #   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
## #   scale_size(range=c(1,8), guide=FALSE) +
## #   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
## #         plot.title = element_text(size=20, hjust = 0.5),
## #         plot.caption = element_text(size=10, hjust = 0.5),
## #         legend.text=element_text(size=18),
## #         axis.text = element_text(size=16,hjust=0),
## #         axis.title = element_text(size=24),
## #         strip.text=element_text(size=16),
## #         legend.position = "bottom",
## #         legend.title=element_text(size=18),
## #         strip.background = element_blank(),
## #         panel.grid.major = element_blank(),
## #         panel.grid.minor = element_blank(),
## #         axis.text.x = element_text(angle = 90,size=16,hjust=0,color="black"),
## #         axis.text.y = element_text(size=16,vjust=0,color="black"),
## #         axis.ticks = element_blank(),
## #         axis.line = element_blank())
## # #dev.off()
## 
## ## Remove leading 0 from r:
## #corr_dat$rs <- sub("^0+", "", round(corr_dat$r,1))
## corr_dat$r1 <- round(corr_dat$r,2)
## 
## #pdf('../figures/cov001_corr_heatmap_version8_r.pdf', width=10, height=10)
## est <- ggplot(corr_dat, aes(measure1, measure2, col=r, label=r1)) +
##   geom_tile(col="black", fill="white") +
##   geom_point(aes(size = abs(r)), shape=15) +
##   geom_text(color="black", size=4) +
##   labs(x = NULL, y = NULL, col = "Correlation") +
##   theme_minimal() +
##   #scale_color_gradient2(mid="#FBFEF9", low="royalblue3", high="red2", limits=c(-1,1)) +
##   scale_colour_gradientn(colours = c("royalblue2", "#FBFEF9", "red"),
##                          limits=c(-1,1), breaks=c(-1,0,1),
##                          labels=c("Negative","","Positive")) +
##   scale_x_discrete(expand=c(0,0), position="top", drop=FALSE) +
##   scale_y_discrete(expand=c(0,0), position="left", drop=FALSE) +
##   scale_size(range=c(1,8), guide=FALSE) +
##   theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "in"),
##         plot.title = element_text(size=20, hjust = 0.5),
##         plot.caption = element_text(size=10, hjust = 0.5),
##         legend.text=element_text(size=18),
##         axis.text = element_text(size=16,hjust=0),
##         axis.title = element_text(size=24),
##         strip.text=element_text(size=16),
##         legend.position = "bottom",
##         legend.title=element_text(size=18),
##         strip.background = element_blank(),
##         panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(),
##         axis.text.x = element_text(angle = 90,size=16,hjust=0,color="black"),
##         axis.text.y = element_text(size=16,vjust=0,color="black"),
##         axis.ticks = element_blank(),
##         axis.line = element_blank())
## #dev.off()
## 
## 
## pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/est_heatmap_high.pdf',
##     width = 8, height = 8);
## 
## est
## dev.off()
## 
## 
## png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/est_heatmap_high.png',
##     res = 100,pointsize = 3,
##     width = 800, height = 800);
## 
## est
## dev.off()

