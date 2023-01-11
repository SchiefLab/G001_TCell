# This code reads in cvd725_epitope_mapping_ics_2022DEC28.csv 
# to map the immunodominant T cell responses to eOD-GT8 and LumSyn (Figure 3).


## ----package-loading-and-options, include=FALSE---------------------------------------------------------------------------------------------------------------------------
packages_needed <- c("conflicted", "tidyverse", "knitr", "kableExtra", 
                      "VISCfunctions", "cowplot", 'reticulate', 'readxl',
                     'here','janitor')

install_load_cran_packages(packages_needed)

# Create a theme 
visc_theme <- theme_bw() + 
  theme(
    legend.position = "bottom", 
    legend.margin = margin(unit = "cm"),
    panel.grid.minor = element_blank()
    )

theme_set(visc_theme)

# Set group colors for this example.
# If you have group colors already in your project:
# source(here::here("docs", "group-colors.R"))
group_colors <- c(
  `No Response` = "#787873",
  `DPBS sucrose` = "#D92321",
  `20 µg eOD-GT8 60mer + AS01B` = "#1749FF",
  `100 µg eOD-GT8 60mer + AS01B` = '#0AB7C9'
)


## ----conflicted-preferences, message=FALSE--------------------------------------------------------------------------------------------------------------------------------

# filter() is used in both dplyr and stats, so need to set the preference to dplyr
conflict_prefer("filter", "dplyr")
conflict_prefer("fisher.test", "stats")


## ----constants------------------------------------------------------------------------------------------------------------------------------------------------------------

# Must have at least 3 allele present at 3 responders for HLA comparisons
HLA_PRESENT_NEEDED <- 4
HLA_RESPONSE_NEEDED <- 4




## ----data-package, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------
## 
## remotes::install_git("/Volumes/networks/cavd/studies/cvd725/pdata/Feinberg725.git", ref = "epm_ics_2021-april")
## 


## ----data-processing------------------------------------------------------------------------------------------------------------------------------------------------------

##Here is the full lum syn seq: MQIYEGKLTAEGLRFGIVASRFNHALVDRLVEGAIDAIVRHGGREEDITLVRVPGSWEIPVAAGELARKEDIDAVIAIGVLIRGATPHFDYIASEVSKGLADLSLELRKPITFGVITADTLEQAIERAGTKHGNKGWEAALSAIEMANLFKSLR


lum_names <- paste0('LumSyn ', 1:41)
HxB2_names <- paste0('HxB2.3-', 1:14)
all_names <- c(lum_names, HxB2_names)
all_names_dat <- tibble(antigen = all_names) %>% 
  mutate(peptide_group = case_when(
           antigen %>% str_detect('LumSyn') ~ 'LumSyn',
           antigen %>% str_detect('HxB2.3') ~ 'HxB2.3'
         ) %>% factor(levels = c('LumSyn','HxB2.3')))


mapping_data <- read.csv("N:/cavd/studies/cvd725/pdata/cohen_et_al_manuscript/data/cvd725_epitope_mapping_ics_2022DEC28.csv")

## ----data-out, eval=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------

#outputting file for Andrew to check mapping algorithm
ics_peptide_data_out <- full_join(
  Feinberg725_seq_positions %>%
    select(peptide_id, aa_sequence, start, end) %>%
    distinct(),
  ics_dat %>% 
    filter(!is.na(population)),
  by = "peptide_id") %>%
  select(peptide_group = peptide_group, peptide_id, aa_sequence, start, end,
         pubid, Parent = parent_plot, response = response_MIMOSA,
         treat)

## import pandas as pd

## from epitope_mapping_df import *

## 

## 

## lumsyn = 'MQIYEGKLTAEGLRFGIVASRFNHALVDRLVEGAIDAIVRHGGREEDITLVRVPGSWEIPVAAGELARKEDIDAVIAIGVLIRGATPHFDYIASEVSKGLADLSLELRKPITFGVITADTLEQAIERAGTKHGNKGWEAALSAIEMANLFKSLRGGSGG'

## 

## # r.ics_peptide_data_out pulls the ics_peptide_data_out object from the previous chunk

## df = r.ics_peptide_data_out

## df = df.query('response==1')

## df = df.rename({'start':'start',

##                   'aa_sequence':'seq',

##                   'end':'end'}, axis=1)

## 

## df = df.assign(RespID=['POS%3d' % i for i in range(df.shape[0])],

##                  protein='vaccine',

##                  start=df['start'] - 1)

## 

## df = df.assign(start_validated=df.apply(lambda r: r['seq'] == lumsyn[r['start']:r['end']], axis=1),

##                len_validated=df.apply(lambda r: r['end'] - r['start'] == len(r['seq']), axis=1))

## 

## def identify_epitopes(pos):

##     """Assign islands to each response"""

##     island_col = pos.groupby('pubid').apply(assignResponseIslands)

##     pos = pd.concat((pos, island_col), axis=1)

## 

##     ep = []

##     for (pubid, iid), gby in pos.groupby(['pubid', 'IslandID']):

##         fe = findEpitopes(gby, overlapRule, returnResponseIndexed=False, minSharedAA=0, minOverlap=8)

##         ep.append(fe)

##     ep = pd.concat(ep, axis=0)

##     ep = pd.merge(ep, pos, on=['pubid', 'IslandID', 'RespID'], how='left')

##     ep = ep.assign(EpID=ep.apply(lambda r: '-'.join(r[['IslandID', 'EpID']]), axis=1),

##                     EpRespSeq=ep.apply(sliceRespSeq, axis=1))

## 

##     vs = ep.groupby(['pubid', 'IslandID','EpID'])['EpRespSeq'].agg(encodeVariants)

##     vs.name = 'EpSeqVars'

##     ep = ep.join(vs, on=['pubid', 'IslandID', 'EpID'])

##     return ep

## 

## cd4 = df.query('Parent=="CD3+/CD4+"')

## cd8 = df.query('Parent=="CD3+/CD8+"')

## 

## cd4_ep = identify_epitopes(cd4)

## cd8_ep = identify_epitopes(cd8)

## 

## ep = pd.concat((cd4_ep, cd8_ep), axis=0)

## 

## 

## 


## ----seq-ics-data, warning=FALSE, message=FALSE---------------------------------------------------------------------------------------------------------------------------

# If coming from python chunk would be called py$ep.
# Just reading it in from repo
py$ep <- bind_rows(
  read_csv(here::here('ICS','epitope_mapping_report','epitope_mapping_results', 'epitopes.csv')) %>% 
    select(-`...1`) %>% 
    mutate(peptide_group = 'LumSyn'),
  read_csv(here::here('ICS','epitope_mapping_report','epitope_mapping_results', 'env_epitopes.csv')) %>% 
  select(-`...1`)
)

ics_epitope_seq_data <-
  full_join(
    py$ep %>% 
      distinct(EpSeq, EpStart, EpEnd, EpID, pubid, IslandID, Parent, peptide_group) %>% 
      mutate(parent = basename(Parent)),
    ics_dat %>% 
      distinct(pubid, treat, peptide_group, parent),
    by = c('pubid','peptide_group', 'parent')
) %>% 
  mutate() %>% 
  mutate(
    parent_plot = paste0('CD3+/', parent),
    peptide_group = peptide_group %>% factor(levels = c('LumSyn','HxB2.3')),
    # Need to impute some values for no 
    response = ifelse(is.na(EpStart), 0, 1),
    EpStart = case_when(
      is.na(EpStart) & peptide_group == 'LumSyn' ~ 88, 
      is.na(EpStart) & peptide_group == 'HxB2.3' ~ 300, 
      TRUE ~ EpStart),
    EpEnd = ifelse(is.na(EpEnd), EpStart + 1, EpEnd),
    )



## ----ics-summary-stats----------------------------------------------------------------------------------------------------------------------------------------------------
 response_fun <- function(xx) {
          wilson_est <- wilson_ci(xx, .95)
          paste0(sum(xx), "/", length(xx), " = ", stat_paste(wilson_est$mean * 
            100, wilson_est$lower * 100, wilson_est$upper * 
            100, digits = 1, trailing_zeros = TRUE, 
            suffix = '\\%'))
        }

ics_response_stats <- ics_dat %>% 
  filter(!is.na(response_MIMOSA),
        antigen %in% all_names) %>% 
  group_by(parent_plot, antigen_factor, peptide_group) %>% 
  summarise(
    total_info = response_fun(response_MIMOSA),
    low_info = response_fun(response_MIMOSA[group == 1]),
    high_info = response_fun(response_MIMOSA[group == 2]),
    `.groups` = "drop"
  )

# special pooled cases
#added one additional cases from kristen 8/22/22
pooled_cases <- ics_dat %>% 
  filter(!is.na(response_MIMOSA),
         antigen %in% all_names) %>% 
  pivot_wider(id_cols = c(pubid, parent_plot, group),
              names_from = antigen_factor, 
              values_from = response_MIMOSA) %>% 
  mutate(
    `LumSyn 5/6` = pmax(`LumSyn 5`, `LumSyn 6`),
    `LumSyn 22/23` = pmax(`LumSyn 22`, `LumSyn 23`),
    `LumSyn 24/25` = pmax(`LumSyn 24`, `LumSyn 25`),
    `LumSyn 28/29` = pmax(`LumSyn 28`, `LumSyn 29`),
    
    `HxB2.3-2/3` = pmax(`HxB2.3-2`, `HxB2.3-3`),
    `HxB2.3-12/13` = pmax( `HxB2.3-12`, `HxB2.3-13`)
  ) %>% 
    mutate(
    `LumSyn 22/23, 28/29` = pmax(`LumSyn 22/23`, `LumSyn 28/29`),
    `LumSyn 5/6, 22/23, 28/29` = pmax(`LumSyn 5/6`, `LumSyn 22/23`, `LumSyn 28/29`)
  ) %>% 
  pivot_longer(cols = `LumSyn 5/6`:`LumSyn 5/6, 22/23, 28/29`, 
               names_to = 'antigen',values_to = 'response_MIMOSA'
                ) %>% 
  select(pubid, parent_plot, group, antigen, response_MIMOSA) %>% 
  mutate(
    peptide_group = case_when(
      antigen %>% str_detect('LumSyn') ~ 'LumSyn',
      antigen %>% str_detect('HxB2.3') ~ 'HxB2.3'
    ) %>% factor(levels = c('LumSyn','HxB2.3')),
    antigen_factor = antigen %>% fct_inorder()
  ) %>% 
  filter(
    !is.na(response_MIMOSA),
    (antigen != 'LumSyn 24/25' & parent_plot == 'CD3+/CD4+') | 
      (antigen == 'LumSyn 24/25' & parent_plot == 'CD3+/CD8+')
         )


pooled_ics_response_stats <- pooled_cases %>% 
  group_by(parent_plot, antigen_factor, peptide_group) %>% 
  summarise(
    total_info = response_fun(response_MIMOSA),
    low_info = response_fun(response_MIMOSA[group == 1]),
    high_info = response_fun(response_MIMOSA[group == 2]),
    `.groups` = "drop"
  )




## ----breadth-info---------------------------------------------------------------------------------------------------------------------------------------------------------

ics_breadth <- ics_dat %>% 
  filter(!is.na(response_MIMOSA)) %>% 
  group_by(parent_plot, treat, peptide_group) %>% 
  summarise(
    breadth = sum(response_MIMOSA),
    n = n(),
    breadth_prop = breadth / n,
    `.groups` = "drop"
  )

epitope_breadth <- ics_epitope_seq_data  %>% 
  group_by(pubid, parent_plot, treat, peptide_group) %>% 
  summarise(
    breadth = EpID[response == 1] %>% unique %>% length,
    `.groups` = "drop"
  )

all_breadth <- 
  bind_rows(
    ics_breadth %>% 
      mutate(measure = 'Peptide Breadth'),
    epitope_breadth %>% 
      mutate(measure = 'Epitope Breadth')
  ) %>% 
  # adding 0s when missing
  pivot_wider(
    id_cols = c(pubid, measure, treat, peptide_group),
    names_from = parent_plot,
    values_from = breadth,
    values_fill = 0
  ) %>% 
  pivot_longer(
    cols = c(`CD3+/CD4+`, `CD3+/CD8+`),
    names_to = 'parent_plot',
    values_to = 'breadth'
  ) %>% 
  group_by(peptide_group) %>% 
  mutate(peptide_group_n = peptide_group %>% 
           fct_recode(`LumSyn (n=41)` = 'LumSyn',
                      `HxB2.3 (n=14)` = 'HxB2.3'),
         measure = measure %>% fct_inorder()) %>% 
  ungroup()

breadth_summary <- all_breadth %>% 
  group_by(measure, parent_plot, peptide_group, treat) %>% 
  summarise(
    mean_breadth = mean(breadth),
    ci_info = paste0(
      '(',
      paste0(
        round_away_0(pmax(0, mean_breadth + 
                       c(-1,1) * qnorm(.95)*(sd(breadth)/sqrt(10))),
                     2, trailing_zeros = TRUE),
        collapse = ', '),
      ')'),
    mean_breadth_info = paste0(
      round_away_0(mean_breadth, 2, trailing_zeros = TRUE),
      ' ', ci_info),
    med_info = paste0(
      round_away_0(median(breadth), 2), ' [',
      min(breadth) , ', ',
      max(breadth), ']'
    ),
    `.groups` = "drop"
  )





## ----hla-results----------------------------------------------------------------------------------------------------------------------------------------------------------


hla_dat_long <- Feinberg725_hla %>% 
  pivot_longer(!pubid,
               names_to = "allele", 
               values_to = "present") %>% 
  #any allele counts are present
  mutate(
    present = present > 0,
    class = case_when(
     allele %>% str_detect('^D') ~ 'II',
     TRUE ~ 'I'
    )
  )

hla_summary <- hla_dat_long %>% 
  group_by(allele, class) %>% 
  summarise(
    num_present = sum(present),
    num_info = paste0(num_present, ' (', round_away_0(num_present / n() * 100), '%)'),
    `.groups` = "drop"
  )



ics_hla_dat <- full_join(
  ics_dat %>%
    filter(antigen %in% all_names),
  hla_dat_long,
  by = 'pubid'
) %>% mutate(
  antigen = factor(antigen, levels = all_names),
)


ics_hla_results <- ics_hla_dat %>% 
  # Only want class 1 for CD8 and class 2 for CD4
  filter(
    (parent == 'CD4+' & class == 'II') |
      (parent == 'CD8+' & class == 'I')
  ) %>% 
  # Must have at least 5 allele present at 3 responders
  group_by(peptide_group, antigen, parent = parent_plot, allele) %>% 
  filter(
    sum(present) >= HLA_PRESENT_NEEDED,
    sum(present) <= n() - HLA_PRESENT_NEEDED,
    sum(response_MIMOSA == 1) >= HLA_RESPONSE_NEEDED,
    sum(response_MIMOSA == 0) >= HLA_RESPONSE_NEEDED
  ) %>% 
  summarise(
    present_response = sum(response_MIMOSA & present),
    not_present_response = sum(response_MIMOSA & !present),
    present_no_response = sum(!response_MIMOSA & present),
    not_present_no_response = sum(!response_MIMOSA & !present),
    pval = two_samp_bin_test(present, response_MIMOSA, method = 'fisher'),
    `.groups` = "drop"
  ) %>% 
  mutate(
    fdr_pval = p.adjust(pval, method = 'fdr')
  ) %>% 
  left_join(Feinberg725_seq_positions %>% 
              select(peptide_id, seq = aa_sequence, seq_start = start) %>% 
              distinct(),
            by = c('antigen' = 'peptide_id')) %>% 
  arrange(pval)





## ----useful-vals, echo = F, warning = F, message = F----------------------------------------------------------------------------------------------------------------------

total_n =  n_distinct(ICS_samps$sample)
total_ant = n_distinct(Feinberg725_ICS_epm$antigen)
total_treat = n_distinct(ICS_samps$treat)


group_total = ICS_samps %>% group_by(treat) %>% 
  dplyr::summarise(n = n_distinct(sample))
total_low <- group_total %>% filter(treat %>% str_detect('20')) %>% pluck('n')
total_high <- group_total %>% filter(treat %>% str_detect('100')) %>% pluck('n')

analysis_low <- Feinberg725_ICS_epm %>% 
  filter(treat == "20 µg eOD-GT8 60mer + AS01B") %>% 
  dplyr::summarise(n = n_distinct(sample))

analysis_high <- Feinberg725_ICS_epm %>% 
  filter(treat == "100 µg eOD-GT8 60mer + AS01B") %>% 
  dplyr::summarise(n = n_distinct(sample))



## ---- out.width = "400px",fig.pos="H", echo = F, fig.cap="Feinberg 725 study schema. Sample sizes per group are presented as vaccine/placebo."----------------------------

include_graphics(here::here( 'docs',"schema.png"))



## ----cohort-table, results = "asis", eval = T, message = F, warning=kable_warnings----------------------------------------------------------------------------------------

n_table <- ICS_samps %>% 
  mutate(
    filter_tab = filter_reason %>% 
      fct_collapse(Samples = c('no response', 'not filtered')) %>% 
      fct_relevel('Samples')
  ) %>% 
  count(`T-cell Subset` = parent, Group = peptide_group, Treatment = treat, filter_tab) %>% 
  pivot_wider(names_from = filter_tab, values_from = n, values_fill = 0)
  
# n_table <- ics_dat %>%
#   filter(antigen %in% lum_names) %>% 
#   distinct(pubid, parent_plot, treat) %>% 
#   count(parent_plot, Dose = treat) %>% 
#   pivot_wider(names_from = parent_plot, values_from = n)

kable(n_table,
      caption.short = "Sample size by T-cell Subset, peptide, treatment, and filter reason",
      caption = "Sample size by T-cell subset, peptide, treatment, and filter reason. One CD4+ sample for a Lumazine Synthase peptide pool (LumSyn 1) had low CD4+ cell count and was removed.", 
      format = output_type, escape = T, longtable = F, 
      booktabs = T, linesep = "") %>% 
  kable_styling(font_size = 8, latex_options = c("hold_position", "repeat_header")) %>% 
  collapse_rows(1:3, row_group_label_position = 'identity', 
                            headers_to_remove = 1:3, latex_hline = 'full')



## ----ics-response-tab, results="asis", warning=FALSE, message=FALSE-------------------------------------------------------------------------------------------------------


ics_response_stats %>%
  filter(peptide_group == 'LumSyn') %>% 
  select(
    `T-cell Subset` = parent_plot,
    # ` ` = peptide_group,
    `Peptide Pool` = antigen_factor,
    `All` = total_info,
    `Lose Dose` = low_info,
    `High Dose` = high_info
  ) %>% 
  kable(
    format = output_type, longtable = TRUE, booktabs = TRUE,
    linesep = "", escape = FALSE,
    caption.short = "CD4+ and CD8+ T-cell responses for Lumazine synthase by peptide and treatment",
    caption = "CD4+ and CD8+ T-cell responses for Lumazine synthase by peptide and treatment. Response rates (MIMOSA) and 95\\% Wilson confidence intervals (CI) are presented by T-cell subset, peptide, and treatment."
  ) %>%
  add_header_above(c(" " = 2, "Response Rate (95% CI)" = 3)) %>% 
  kable_styling(
    font_size = 6.25,
    # Note scale_down will overwrite font_size specifications
    latex_options = c("hold_position", "repeat_header"),
    repeat_header_method = 'replace'
  ) %>% 
  collapse_rows(columns = 1:2, row_group_label_position = 'identity', 
                latex_hline = 'full', 
                valign = 'top', longtable_clean_cut = TRUE)


## ----rr-lumazine-synthase-figs, fig.scap="Peptide pool and epitope breadth boxplots by treatment group", fig.cap= "Peptide pool and epitope breadth boxplots by treatment group.", fig.height=7.5, message = F----


bar_plot_dat_ind <- ics_dat %>% 
  filter(!is.na(response_MIMOSA),
        str_detect(antigen, "LumSyn")) %>% 
  group_by(parent_plot, antigen, peptide_group, treat) %>% 
  summarise(n = length(response_MIMOSA),
   
    resp = sum(response_MIMOSA)
    ) %>% 
  ungroup() %>% 
  mutate(prop = resp/n) %>% 
  mutate(peptide = str_sub(antigen, 8,9)) %>% 
  filter(str_detect(parent_plot, "CD4")) %>% 
  arrange(peptide) 



bar_plot_dat_pooled <- pooled_cases %>% 
  full_join(dplyr::select(ics_dat, pubid, treat) %>% distinct()) %>% 
  filter(str_detect(antigen, "LumSyn") & parent_plot == "CD3+/CD4+") %>% 
  group_by(parent_plot, antigen, peptide_group, treat) %>% 
  summarise(
  n = length(response_MIMOSA),
   
    resp = sum(response_MIMOSA)
  ) %>% 
   ungroup() %>% 
  mutate(prop = resp/n) %>% 
  mutate(peptide = str_sub(antigen, start = 8)) 

bar_plot_dat = bind_rows(bar_plot_dat_ind, bar_plot_dat_pooled) %>% 
  mutate(peptide_plot = factor(peptide, levels = c(1,2,3,4,5,
                                                   6:22,
                                                   23,24,
                                                   25, 26, 27, 28, 29:41, "5/6","22/23","24/25",
                                                   "28/29","22/23, 28/29","5/6, 22/23, 28/29"), ordered = T)) %>% 
  mutate(peptide_group = case_when(str_detect(peptide,"/") ~ 2,
                                   !str_detect(peptide,"/") ~ 1)) %>% 
  arrange(peptide_group, peptide_plot)




bar_plot_wide <- ggplot(data = bar_plot_dat, aes(fill = treat, x = peptide_plot, y = prop*100, color = treat)) +
  geom_bar(position="dodge", stat="identity") +  
  scale_fill_manual(name = "", values = group_colors[3:4]) +
  scale_color_manual(name = "", values = group_colors[3:4]) +
  scale_x_discrete("") +
  scale_y_continuous(
    "% of vaccinees with a \npositive response",
    limits = c(0,100),
    breaks = seq(0, 100, by = 20)
    ) +facet_grid(.~peptide_group,scales = "free_x", space = "free") +
     theme(
        axis.text.x = element_text(size = 12,angle = 50, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
strip.background = element_blank(), strip.text = element_blank(),
    axis.text.y = element_text(size = 12)
    ) 


  
  
    bar_plot_tall <- ggplot(data = bar_plot_dat, aes(fill = treat, x = fct_rev(peptide_plot), y = prop*100, color = treat)) +
  geom_bar(position="dodge", stat="identity") +  
  scale_fill_manual(name = "", values = group_colors[3:4]) +
      scale_color_manual(name = "", values = group_colors[3:4]) +

  scale_x_discrete("") +
  scale_y_continuous(
    "% of vaccinees with a \npositive response",
    limits = c(0,100),
    breaks = seq(0, 100, by = 20)
    ) + coord_flip() +facet_grid(peptide_group~.,scales = "free_y", space = "free") + 
  theme(
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12),
strip.background = element_blank(), strip.text = element_blank(),
    axis.text.y = element_text(size = 12)
    ) 


bar_plot_tall
  
    
# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_wide.pdf',
#     width = 12, height = 5.5);
# bar_plot_wide
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_wide.png',
#     res = 100,pointsize = 3,
#     width = 1200, height = 550);
# bar_plot_wide
# dev.off()
# 
# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_tall.pdf',
#     height = 13, width = 6);
# bar_plot_tall
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_tall.png',
#     res = 100,pointsize = 3,
#     height = 1300, width = 600);
# bar_plot_tall
# dev.off()


## ----rr-eOD-GT8-3-figs, fig.scap="Peptide pool and epitope breadth boxplots by treatment group", fig.cap= "Peptide pool and epitope breadth boxplots by treatment group.", fig.height=7.5, message = F----


bar_plot_dat_ind <- ics_dat %>% 
  filter(!is.na(response_MIMOSA),
        str_detect(antigen, "HxB2.3")) %>% 
  group_by(parent_plot, antigen, peptide_group, treat) %>% 
  summarise(n = length(response_MIMOSA),
   
    resp = sum(response_MIMOSA)
    ) %>% 
  ungroup() %>% 
  mutate(prop = resp/n) %>% 
  mutate(peptide = str_sub(antigen, 8,9)) %>% 
  filter(str_detect(parent_plot, "CD4")) %>% 
  arrange(peptide) 



bar_plot_dat_pooled <- pooled_cases %>% 
  full_join(dplyr::select(ics_dat, pubid, treat) %>% distinct()) %>% 
  filter(str_detect(antigen, "HxB2.3") & parent_plot == "CD3+/CD4+") %>% 
  group_by(parent_plot, antigen, peptide_group, treat) %>% 
  summarise(
  n = length(response_MIMOSA),
   
    resp = sum(response_MIMOSA)
  ) %>% 
   ungroup() %>% 
  mutate(prop = resp/n) %>% 
  mutate(peptide = str_sub(antigen, start = 8)) 

bar_plot_dat = bind_rows(bar_plot_dat_ind, bar_plot_dat_pooled) %>% 
  mutate(peptide_plot = factor(peptide, levels = c(1:14, "2/3","12/13"), ordered = T)) %>% 
  mutate(peptide_group = case_when(str_detect(peptide,"/") ~ 2,
                                   !str_detect(peptide,"/") ~ 1)) %>% 
  arrange(peptide_group, peptide_plot)




bar_plot_wide <- ggplot(data = bar_plot_dat, aes(fill = treat, x = peptide_plot, y = prop*100, color = treat)) +
  geom_bar(position="dodge", stat="identity") +  
  scale_fill_manual(name = "", values = group_colors[3:4]) +
  scale_color_manual(name = "", values = group_colors[3:4]) +
  scale_x_discrete("") +
  scale_y_continuous(
    "% of vaccinees with a \npositive response",
    limits = c(0,100),
    breaks = seq(0, 100, by = 20)
    ) +facet_grid(.~peptide_group,scales = "free_x", space = "free") +
     theme(
        axis.text.x = element_text(size = 12,angle = 50, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
strip.background = element_blank(), strip.text = element_blank(),
    axis.text.y = element_text(size = 12)
    ) 


  bar_plot_wide_resp_only <- ggplot(data = filter(bar_plot_dat,resp > 0),
                          aes(fill = treat, x = peptide_plot, y = prop*100, color = treat)) +
  geom_bar(position="dodge", stat="identity") +  
  scale_fill_manual(name = "", values = group_colors[3:4]) +
      scale_color_manual(name = "", values = group_colors[3:4]) +


  scale_x_discrete("") +
  scale_y_continuous(
    "% of vaccinees with a \npositive response",
    limits = c(0,100),
    breaks = seq(0, 100, by = 20)
    ) +facet_grid(.~peptide_group,scales = "free_x", space = "free") +
     theme(
        axis.text.x = element_text(size = 12,angle = 50, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
strip.background = element_blank(), strip.text = element_blank(),
    axis.text.y = element_text(size = 12)
    ) 

  
  
    bar_plot_tall_resp_only <- ggplot(data = filter(bar_plot_dat, resp > 0), 
                                      aes(fill = treat, x = fct_rev(peptide_plot), y = prop*100, color = treat)) +
  geom_bar(position="dodge", stat="identity") +  
  scale_fill_manual(name = "", values = group_colors[3:4]) +
      scale_color_manual(name = "", values = group_colors[3:4]) +

  scale_x_discrete("") +
  scale_y_continuous(
    "% of vaccinees with a \npositive response",
    limits = c(0,100),
    breaks = seq(0, 100, by = 20)
    ) + coord_flip() +facet_grid(peptide_group~.,scales = "free_y", space = "free") + 
  theme(
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12),
strip.background = element_blank(), strip.text = element_blank(),
    axis.text.y = element_text(size = 12)
    ) 


bar_plot_tall
  
    
# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_wide_eod3.pdf',
#     width = 10, height = 4.5);
# bar_plot_wide
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_wide_eod3.png',
#     res = 100,pointsize = 3,
#     width = 1000, height = 450);
# bar_plot_wide
# dev.off()
# 
# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_tall_eod3.pdf',
#     height = 10, width = 4.5);
# bar_plot_tall
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_tall_eod3.png',
#     res = 100,pointsize = 3,
#     height = 1000, width = 450);
# bar_plot_tall
# dev.off()

# 
# 
# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_wide_eod3_responly.pdf',
#     width = 8, height = 4.5);
# bar_plot_wide_resp_only
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_wide_eod3_responly.png',
#     res = 100,pointsize = 3,
#     width = 500, height = 200);
# bar_plot_wide_resp_only
# dev.off()
# 
# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_tall_eod3_responly.pdf',
#     height = 8, width = 4.5);
# bar_plot_tall_resp_only
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/bar_plot_tall_eod3_responly.png',
#     res = 100,pointsize = 3,
#     height = 500, width = 200);
# bar_plot_tall_resp_only
# dev.off()


## ----HxB2-ics-response-tab, results="asis", warning=FALSE, message=FALSE--------------------------------------------------------------------------------------------------


ics_response_stats %>%
  filter(peptide_group == 'HxB2.3') %>% 
  select(
    `T-cell Subset` = parent_plot,
    # ` ` = peptide_group,
    `Peptide Pool` = antigen_factor,
    `All` = total_info,
    `Lose Dose` = low_info,
    `High Dose` = high_info
  ) %>% 
  kable(
    format = output_type, longtable = TRUE, booktabs = TRUE,
    linesep = "", escape = FALSE,
    caption.short = "CD4+ and CD8+ T-cell responses for HxB2.3 by peptide and treatment",
    caption = "CD4+ and CD8+ T-cell responses for HxB2.3 by peptide and treatment. Response rates (MIMOSA) and 95\\% Wilson confidence intervals (CI) are presented by T-cell subset, peptide, and treatment."
  ) %>%
  add_header_above(c(" " = 2, "Response Rate (95% CI)" = 3)) %>% 
  kable_styling(
    font_size = 6.25,
    # Note scale_down will overwrite font_size specifications
    latex_options = c("hold_position", "repeat_header"),
  ) %>% 
  collapse_rows(columns = 1:2, row_group_label_position = 'identity', 
                latex_hline = 'full', 
                valign = 'top', longtable_clean_cut = TRUE)


## ----breadth-figs, fig.scap="Peptide pool and epitope breadth boxplots by treatment group", fig.cap= "Peptide pool and epitope breadth boxplots by treatment group.", fig.height=7.5----

plot1 <- all_breadth %>%
  filter(measure == 'Peptide Breadth') %>% 
  ggplot(aes(x = treat, y = breadth, color = treat)) +
  geom_point(
    position = position_jitter(width = .25, height = 0, seed = 3241),
    size = 1, show.legend = FALSE
    ) +
  geom_boxplot(
    fill = NA, lwd = .5, outlier.colour = NA, show.legend = FALSE
    ) +  
  scale_color_manual(name = "", values = group_colors) +
  scale_x_discrete("") +
  scale_y_continuous(
    "Peptide Breadth",
    breaks = seq(0, 20, by = 2)
    ) +
  facet_grid(peptide_group_n~parent_plot) +
  theme(
    axis.text.x = element_text(size = 7, angle = 30, vjust = 1, hjust = 1)
    )

plot2 <- all_breadth %>%
  filter(measure == 'Epitope Breadth') %>% 
  ggplot(aes(x = treat, y = breadth, color = treat)) +
  geom_point(
    position = position_jitter(width = .25, height = 0, seed = 3241),
    size = 1, show.legend = FALSE
    ) +
  geom_boxplot(
    fill = NA, lwd = .5, outlier.colour = NA, show.legend = FALSE
    ) +  
  scale_color_manual(name = "", values = group_colors) +
  scale_x_discrete("") +
  scale_y_continuous(
    "Epitope Sequence Breadth",
    breaks = seq(0, 20, by = 2)
    ) +
  facet_grid(peptide_group~parent_plot) +
  theme(
    axis.text.x = element_text(size = 7, angle = 30, vjust = 1, hjust = 1)
    )

plot_grid(plot1, plot2, nrow = 2, align = 'v')





## ----breadth-tab, results="asis", warning=kable_warnings------------------------------------------------------------------------------------------------------------------


breadth_summary %>%
  select(Measure = measure, `T Cell Subset` = parent_plot,
         Peptide = peptide_group, Treatment = treat,
         `Mean (95\\% CI)` = mean_breadth_info,  `Median (Range)` = med_info) %>% 
  kable(
    format = output_type, longtable = FALSE, booktabs = TRUE,
    linesep = "", escape = FALSE,
    caption.short = "Peptide pool and epitope breadth summary statistics",
    caption = "Peptide pool and epitope breadth summary statistics. Mean counts of positive responses are presented with 95\\% Wilson score confidence intervals. Median, min and max counts are also displayed."
  ) %>%
  kable_styling(
    font_size = 8,
    # Note scale_down will overwrite font_size specifications
    latex_options = c("hold_position", "scale_down", "repeat_header")
  ) %>% 
  collapse_rows(columns = 1:3, row_group_label_position = 'identity', valign = 'top', 
                latex_hline = 'full', headers_to_remove = 1:3)


## ----peptide-response, fig.scap="Lumazine synthase and HxB2.3 peptide response maps", fig.cap= "Lumazine synthase and HxB2.3 peptide response maps. Statistically significant responses are plotted using one horizontal line per pool. Each row (y-axis tick) indicates one participant and non-responders are included as blank rows. Due to the overlapping design of the peptides, recognition of a single epitope may be indicated by several peptide responses.", fig.height=7.5----

ics_peptide_data_plot <-  ics_peptide_data %>%
  mutate(
    pubid = fct_rev(pubid), 
    treat = fct_rev(treat)
  ) %>% 
  arrange(desc(response_MIMOSA), start) %>% 
  # filter(response_MIMOSA == 1) %>% 
  mutate(
    index = 1:n(),
    treat_pubid = fct_cross(pubid,treat),
    treat_pubid_num = as.numeric(treat_pubid)   
  ) %>% 
  group_by(pubid, parent_plot, peptide_group) %>% 
  mutate(y_val = treat_pubid_num + .08 * 1:n()) %>% 
  ungroup() %>% 
  pivot_longer(
    cols = c(start, end),
    names_to = 'start_or_end',
    values_to = 'position'
  )

ics_peptide_data_plot %>%   
  ggplot(aes(x = position, y = y_val, color = treat, group = index)) +
  geom_line(data = ics_peptide_data_plot %>% filter(response_MIMOSA == 1), 
            size = 1, alpha = .75) +  
  scale_color_manual(name = "", values = group_colors[4:3]) +
  scale_x_continuous("Peptide Position", 
                     breaks = c(1, 20, 40, 60, 80, 100, 120, 140, 159,
                                281, 300, 320, 339)) +
  scale_y_continuous("", breaks = unique(ics_peptide_data_plot$treat_pubid_num), 
                   labels =  unique(ics_peptide_data_plot$treat_pubid) %>%
                     str_sub(end = str_locate(unique(ics_peptide_data_plot$treat_pubid),
                                              ':')[,'start'] - 1)
                   ) +
  coord_cartesian(ylim = c(2,33)) +
  facet_grid(parent_plot ~ peptide_group, scales = 'free', space = 'free_x') +
  theme(
    axis.text.y = element_text(size = 7)
    ) +
  guides(color = guide_legend(reverse = TRUE))



## ----epitope-response, fig.scap="Minimal epitope map for lumazine synthase and HxB2.3", fig.cap= "Minimal epitope map for lumazine synthase and HxB2.3. The minimal set of epitopes able to explain each participant’s peptide responses were determined from the 15mer responses. Each row (y-axis tick) represents a participant and non-responders are included as blank rows.", fig.height=7.5----


ics_epitope_seq_data_plot <-  ics_epitope_seq_data %>%
  mutate(
    pubid = fct_rev(pubid), 
    treat = fct_rev(treat)
  ) %>% 
  arrange(EpStart) %>% 
  mutate(
    index = 1:n(),
    treat_pubid = fct_cross(pubid,treat),
    treat_pubid_num = as.numeric(treat_pubid)   
  ) %>% 
  group_by(pubid, parent_plot, peptide_group) %>% 
  mutate(y_val = treat_pubid_num + .08 * 1:n()) %>% 
  pivot_longer(
    cols = c(EpStart, EpEnd),
    names_to = 'start_or_end',
    values_to = 'position'
  )

ics_epitope_seq_data_plot %>%   
  ggplot(aes(x = position, y = y_val, color = treat, group = index)) +
  geom_line(data = ics_epitope_seq_data_plot %>% dplyr::filter(response == 1),
            size = 1, alpha = .75) +  
  scale_color_manual(name = "", values = group_colors[4:3]) +
  scale_x_continuous("Peptide Position", 
                     breaks = c(1, 20, 40, 60, 80, 100, 120, 140, 159,
                                285, 300, 320, 339)) +
  scale_y_continuous("", breaks = unique(ics_epitope_seq_data_plot$treat_pubid_num), 
                   labels =  unique(ics_epitope_seq_data_plot$treat_pubid) %>%
                     str_sub(end = str_locate(unique(ics_epitope_seq_data_plot$treat_pubid),
                                              ':')[,'start'] - 1)
                   ) +
  coord_cartesian(ylim = c(2,33)) +
  facet_grid(parent_plot ~ peptide_group, scales = 'free_x', space = 'free_x') +
  theme(
    axis.text.y = element_text(size = 7)
    ) +
  guides(color = guide_legend(reverse = TRUE))





## ----hotspot-response-tab, results="asis", warning=FALSE, message=FALSE---------------------------------------------------------------------------------------------------


pooled_ics_response_stats %>%
  select(
    `T-cell Subset` = parent_plot,
    ` ` = peptide_group,
    `Combined Peptide Pool` = antigen_factor,
    `All` = total_info,
    `Lose Dose` = low_info,
    `High Dose` = high_info
  ) %>% 
  kable(
    format = output_type, longtable = FALSE, booktabs = TRUE,
    linesep = "", escape = FALSE,
    caption.short = "CD4+ and CD8+ T-cell responses for hotspots by peptide and treatment",
  caption = paste0("CD4+ and CD8+ T-cell responses for hotspots by peptide and treatment. Hotspots correspond with positive peptide positions seen in figures ",  insert_ref('fig:peptide-response')," and ",  insert_ref('fig:epitope-response'), ". A hotspot is considered positive if any peptide pools within the hotspot are positive. Response rates and 95\\% Wilson confidence intervals (CI) are presented by T-cell subset, peptide, and treatment.")
  ) %>%
  column_spec(3, width = '2.15cm') %>% 
  add_header_above(c(" " = 3, "Response Rate (95% CI)" = 3)) %>% 
  kable_styling(
    font_size = 6.25,
    # Note scale_down will overwrite font_size specifications
    latex_options = c("hold_position", "repeat_header"),
  ) %>% 
  collapse_rows(columns = 1:3, row_group_label_position = 'identity', 
                latex_hline = 'full',
                valign = 'top', longtable_clean_cut = TRUE)



## ----hla-out, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
## 
##   write_csv(ics_hla_results, file = here::here(
##   'ICS','epitope_mapping_report',
##   'csv_output', 'hla_comparisons.csv'))
## 


## ----hla-info-tab, results = "asis", eval = T, message = F, warning=kable_warnings----------------------------------------------------------------------------------------

if (output_type == 'latex') {
  hla_info_tab <- hla_summary %>% 
    mutate(
      num_info_tab = case_when(
        num_present >= HLA_PRESENT_NEEDED & 
          num_present <= total_n - HLA_PRESENT_NEEDED ~ 
          cell_spec(num_info, bold = TRUE, format = output_type),
        TRUE ~ cell_spec(num_info, format = output_type)
      ),
      allele = escape(allele)
    )
} else {
    hla_info_tab <- hla_summary %>% 
    mutate(
      num_info_tab = case_when(
        num_present >= HLA_PRESENT_NEEDED & 
          num_present <= total_n - HLA_PRESENT_NEEDED ~ 
          paste0('**', num_info, '**'),
        TRUE ~ num_info
      ),
      allele = escape(allele)
    )
}

bind_cols(
  hla_info_tab %>% filter(class == 'I') %>% 
    select(`Alleles` = allele, 
           `Number Present (\\%)` = num_info_tab),
  hla_info_tab %>% filter(class == 'II') %>% 
    select(`Alleles ` = allele, 
           `Number Present (\\%) ` = num_info_tab)
) %>% 
  kable(
    format = output_type, longtable = FALSE, booktabs = TRUE,
    linesep = "", escape = FALSE,
    caption.short = "HLA information and allele presence",
    caption = paste0("HLA information and allele presence. Alleles with fewer than ", 
                     HLA_PRESENT_NEEDED, " participants having the allele present or fewer than ", HLA_PRESENT_NEEDED," participants having the allele not present are excluded from comparisons. Alleles considered for comparisons are bolded.")
  ) %>%
  column_spec(2, border_right = TRUE) %>% 
  # column_spec(5, width = '1.2cm') %>% 
  # column_spec(7:10, width = '2cm') %>% 
  add_header_above(c("Class I Alleles" = 2, "Class II Alleles" = 2)) %>%
  kable_styling(
    font_size = 7,
    # Note scale_down will overwrite font_size specifications
    latex_options = c("hold_position")
  )



## ----hla-tabs, results='asis', warning=kable_warnings---------------------------------------------------------------------------------------------------------------------

ics_hla_tab <- ics_hla_results %>% 
  filter(pval < 0.05) %>% 
  mutate(
    p_tab = pretty_pvalues(pval, digits = 4),
    p_adj_tab = pretty_pvalues(fdr_pval, digits = 4, background = 'yellow'),
    antigen = factor(antigen, levels = all_names)
  ) %>% 
  arrange(parent, peptide_group, antigen, pval) %>% 
  select(
    parent, Group = peptide_group, `Peptide Pool` = antigen,
    Sequence = seq, `Sequence Start` = seq_start, allele,
    `Present \\& Response` = present_response, 
    `Not Present \\& Response` = not_present_response,
    `Present \\& No Response` = present_no_response, 
    `Not Present \\& No Response` = not_present_no_response,
    `Unadjusted` = p_tab,
    `Adjusted` = p_adj_tab
  ) %>% 
  rename_with(str_to_title)


ics_hla_tab %>%
  kable(
    format = output_type, longtable = FALSE, booktabs = TRUE,
    linesep = "", escape = FALSE,
    caption.short = "Top HLA results comparing allele presence vs. peptide response",
    caption = paste0("Top HLA results comparing allele presence vs. peptide response. MHC class I (A, B, and C) are compared to CD8+ response, and MHC class II (DP, DM, DO, DQ, and DR) are compared to CD4+ response. Fisher's Exact test was computed for each peptide and allele combination, and results with a p value < 0.05 are displayed. Adjusted p values using the False Discovery Method (FDR) are also displayed, with p values < 0.05 highlighted. Alleles with fewer than ", HLA_PRESENT_NEEDED, " participants having the allele present or fewer than  ", HLA_PRESENT_NEEDED, " participants having the allele not present are excluded from comparisons.  Similarly, peptide pools with fewer than  ", HLA_RESPONSE_NEEDED, " participants having a positive response and fewer than ", HLA_RESPONSE_NEEDED, " participants having a negative response are excluded from comparisons.")
  ) %>%
  column_spec(3, width = '1.55cm') %>% 
  column_spec(5, width = '1.2cm') %>% 
  column_spec(7:10, width = '2cm') %>% 
  add_header_above(c(" " = 6, "Allele Present/Peptide Response" = 4, "P Value" = 2)) %>% 
  kable_styling(
    font_size = 8,
    # Note scale_down will overwrite font_size specifications
    latex_options = c("hold_position", "scale_down", "repeat_header")
  ) %>% 
  collapse_rows(columns = 1:4, row_group_label_position = 'identity', 
                latex_hline = 'full', headers_to_remove = 1:2,
                valign = 'top')



## ----setting-up-ind-plots-------------------------------------------------------------------------------------------------------------------------------------------------

pubs_to_run <- ics_peptide_data %>% 
  group_by(pubid, parent_plot) %>% 
  filter(any(response_MIMOSA == 1)) %>% 
  distinct(pubid, parent_plot) %>% 
  arrange(parent_plot, pubid)



## ----ind-plots-cd4, fig.width=8.5, fig.height=7.5, fig.scap=paste0('Individual Plots for ', escape(pubs_to_run$pubid), '(', pubs_to_run$parent_plot, ')'), fig.cap=paste0('Individual Plots for ', escape(pubs_to_run$pubid), '(', pubs_to_run$parent_plot, ')'), eval=FALSE----
## 
## 
## walk2(pubs_to_run$pubid,
##       pubs_to_run$parent_plot,
##       function(xx,yy){
## 
##   plot_data <- ics_peptide_data %>%
##     filter(pubid == xx,
##            parent_plot == yy,
##            !is.na(start)) %>%
##     mutate(peptide_factor = factor(peptide_id,
##                                    levels = unique(peptide_id),
##                                    ordered = T))
## 
## 
## 
##   plot1 <- plot_data %>%
##     ggplot() +
##     geom_segment(aes(x = start, xend = end,
##                      y = peptide_factor, yend = peptide_factor, color = factor(response_MIMOSA)), size = 2 ) +
##     geom_text(data = Feinberg725_seq_positions %>%
##                 right_join(plot_data %>%
##                              distinct(peptide_id, peptide_factor),
##                            by = 'peptide_id'),
##               aes(label = aa, x = position, y = peptide_factor),
##               nudge_y = -.3, size = 1.75) +
##     theme_bw() +
##     theme(panel.grid = element_blank()) +
##     scale_x_continuous("LumSyn Alignment Coordinate", breaks = seq(from = 0, to = 175, by = 10)) +
##     scale_y_discrete("") +
##     scale_color_discrete("MIMOSA Response", breaks = c(1, 0)) +
##     ggtitle(paste0(unique(plot_data$pubid)," (# peptides positive=", sum(plot_data$response_MIMOSA %>% na.omit()), ')')) +
##     facet_wrap(~parent_plot) +
##     theme(plot.title = element_text(size = 8), legend.position = 'bottom',
##           strip.background = element_rect(fill = 'white'))
## 
## 
## 
##   epitope_data_here <- ics_epitope_seq_data %>%
##     filter(pubid == xx,
##            parent_plot == yy) %>%
##     arrange(EpStart) %>%
##     mutate(EpSeq_factor = EpSeq %>% fct_inorder)
## 
##   if (nrow(epitope_data_here) == 0) {
##     print(plot1)
##   } else {
##     plot2 <- epitope_data_here %>%
##       ggplot() +
##       geom_segment(aes(x = EpStart, xend = EpEnd,
##                        y = EpSeq_factor, yend = EpSeq_factor),
##                    color = "#00BFC4", size = 2) +
##       scale_x_continuous("LumSyn Alignment Coordinate", breaks = seq(from = 0, to = 175, by = 10), limits = c(0,175)) +
##       scale_y_discrete("") +
##       ggtitle(paste0(unique(epitope_data_here$pubid)," (# epitopes positive=", nrow(epitope_data_here), ')')) +
##       theme_bw() +
##       theme(panel.grid = element_blank(),
##             plot.title = element_text(size = 8),
##             legend.position = 'bottom',
##             axis.text.y = element_text(size = 5))
## 
## 
##     print(cowplot::plot_grid(plot1, plot2, nrow = 2, rel_heights = c(1,.25),
##                              align = 'v'))
##   }
## 
## })
## 
## 
## 
## 


## ----pep-ind-plots, eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------
## 
## aaa <- expand_grid(
##   xx=levels(ics_peptide_data_plot$peptide_group),
##   yy=levels(ics_peptide_data_plot$parent_plot %>% factor))
## 
## walk2(
##   aaa$xx,
##   aaa$yy,
##   function(xx,yy){
## 
## tmp_plot_data <-  ics_peptide_data %>%
##       filter(peptide_group == xx,
##              parent_plot == yy) %>%
##   mutate(
##     pubid = fct_rev(pubid),
##     treat = fct_rev(treat)
##   ) %>%
##   arrange(desc(response_MIMOSA), start) %>%
##   # filter(response_MIMOSA == 1) %>%
##   mutate(
##     index = 1:n(),
##     treat_pubid = fct_cross(pubid,treat),
##     treat_pubid_num = as.numeric(treat_pubid)
##   ) %>%
##   group_by(pubid, parent_plot, peptide_group) %>%
##   mutate(y_val = treat_pubid_num + .08 * 1:n()) %>%
##   ungroup() %>%
##   pivot_longer(
##     cols = c(start, end),
##     names_to = 'start_or_end',
##     values_to = 'position'
##   )
## 
##     y_range <- range(tmp_plot_data$treat_pubid_num) + c(1,0)
## 
## 
##   aa = tmp_plot_data %>%
##   ggplot(aes(x = position, y = y_val, color = treat, group = index)) +
##   geom_line(data = tmp_plot_data %>% filter(response_MIMOSA == 1),
##             size = 1, alpha = .75) +
##   scale_color_manual(name = "", values = group_colors, drop = FALSE) +
##   scale_x_continuous("Peptide Position",
##                      breaks = c(1, 20, 40, 60, 80, 100, 120, 140, 159,
##                                 281, 300, 320, 339)) +
##   scale_y_continuous("", breaks = unique(tmp_plot_data$treat_pubid_num),
##                      labels =  unique(tmp_plot_data$treat_pubid) %>%
##                        str_sub(end = str_locate(unique(tmp_plot_data$treat_pubid),
##                                                 ':')[,'start'] - 1)
##   ) +
##     coord_cartesian(ylim = y_range,
##                     xlim = switch (xx,
##                                    'LumSyn' = c(1, 159),
##                                    'HxB2.3' = c(281, 339)
##                     )) +
##   facet_grid(parent_plot ~ peptide_group, scales = 'free_x', space = 'free_x') +
##   theme(
##     axis.text.y = element_text(size = 7)
##   ) +
##   guides(color = guide_legend(reverse = TRUE))
## 
## 
## pdf(file = paste0('~/Temp/Plots_for_Kristen/Peptide_', xx, '_',basename(yy) %>% str_replace('\\+',''),'.pdf'),
##     width = 7, height = 4.5)
## print(aa)
## dev.off()
## }
## )
## 
## 
## 
## 
## 
## 
## 
## aaa <- expand_grid(
##   xx=levels(ics_epitope_seq_data_plot$peptide_group),
##   yy=levels(ics_epitope_seq_data_plot$parent_plot %>% factor))
## 
## walk2(
##   aaa$xx,
##   aaa$yy,
##   function(xx,yy){
## 
## tmp_plot_data <-  ics_epitope_seq_data %>%
##       filter(peptide_group == xx,
##              parent_plot == yy) %>%
##   mutate(
##     pubid = fct_rev(pubid),
##     treat = fct_rev(treat)
##   ) %>%
##   arrange(EpStart) %>%
##   mutate(
##     index = 1:n(),
##     treat_pubid = fct_cross(pubid,treat),
##     treat_pubid_num = as.numeric(treat_pubid)
##   ) %>%
##   group_by(pubid, parent_plot, peptide_group) %>%
##   mutate(y_val = treat_pubid_num + .08 * 1:n()) %>%
##   pivot_longer(
##     cols = c(EpStart, EpEnd),
##     names_to = 'start_or_end',
##     values_to = 'position'
##   )
## 
## y_range <- range(tmp_plot_data$treat_pubid_num) + c(1,0)
## 
##   aa = tmp_plot_data %>%
##   filter(peptide_group == xx,
##          parent_plot == yy) %>%
##   ggplot(aes(x = position, y = y_val, color = treat, group = index)) +
##   geom_line(data = tmp_plot_data %>% filter(response == 1),
##             size = 1, alpha = .75) +
##   scale_color_manual(name = "", values = group_colors, drop = FALSE) +
##   scale_x_continuous("Peptide Position",
##                      breaks = c(1, 20, 40, 60, 80, 100, 120, 140, 159,
##                                 281, 300, 320, 339)) +
##   scale_y_continuous("", breaks = unique(tmp_plot_data$treat_pubid_num),
##                      labels =  unique(tmp_plot_data$treat_pubid) %>%
##                        str_sub(end = str_locate(unique(tmp_plot_data$treat_pubid),
##                                                 ':')[,'start'] - 1)
##   ) +
##     coord_cartesian(ylim = y_range,
##                     xlim = switch (xx,
##                                    'LumSyn' = c(1, 159),
##                                    'HxB2.3' = c(281, 339)
##                     )) +
##   facet_grid(parent_plot ~ peptide_group, scales = 'free_x', space = 'free_x') +
##   theme(
##     axis.text.y = element_text(size = 7)
##   ) +
##   guides(color = guide_legend(reverse = TRUE))
## 
## 
## pdf(file = paste0('~/Temp/Plots_for_Kristen/Epitope_', xx, '_',basename(yy) %>% str_replace('\\+',''),'.pdf'),
##     width = 7, height = 4.5)
## print(aa)
## dev.off()
## }
## )
## 
## 
## 
## 
## 


## ----peptide-response-pub, fig.scap="Lumazine synthase and HxB2.3 peptide response maps (with Hot Spots)", fig.cap= "Lumazine synthase and HxB2.3 peptide response maps (with Hot Spots). Statistically significant responses are plotted using one horizontal line per pool. Each row (y-axis tick) indicates one participant and non-responders are included as blank rows. Due to the overlapping design of the peptides, recognition of a single epitope may be indicated by several peptide responses.", fig.height=6.5----

ics_peptide_data_plot <-  ics_peptide_data %>%
  mutate(
    pubid = fct_rev(pubid), 
    treat = fct_rev(treat)
  ) %>% 
  arrange(desc(response_MIMOSA), start) %>% 
  # filter(response_MIMOSA == 1) %>% 
  mutate(
    index = 1:n(),
    treat_pubid = fct_cross(pubid,treat),
    treat_pubid_num = as.numeric(treat_pubid)   
  ) %>% 
  group_by(pubid, parent_plot, peptide_group) %>% 
  mutate(y_val = treat_pubid_num + .08 * 1:n()) %>% 
  ungroup() %>% 
  pivot_longer(
    cols = c(start, end),
    names_to = 'start_or_end',
    values_to = 'position'
  ) %>% 
  mutate(parent_plot = str_sub(parent, 1,3))

bg_color_data <- tibble(
  xmin = c(17, 85, 109,281, 321, 93),
  xmax = c(35, 103, 127,299, 335, 111),
  peptide_group = c('LumSyn','LumSyn',"LumSyn",'HxB2.3','HxB2.3','LumSyn') %>% 
    factor(levels = levels(ics_peptide_data_plot$peptide_group)),
  parent_plot = c(rep('CD4', 5), 'CD8'),
  fill = c( "#440154FF" ,"#2A788EFF" ,"#7AD151FF", "#414487FF", "#22A884FF", "#FDE725FF"),
   # fill = c( "black" ,"yellow" ,"red", "black", "red", "yellow"),

  hotspot_label = c("5/6", "22/23", "28/29", "2/3", "12/13", "24/25")
)

 hotspot_plot_withlabs <- ics_peptide_data_plot %>%   
  ggplot() +
  geom_line(data = ics_peptide_data_plot %>% filter(response_MIMOSA == 1), 
            aes(x = position, y = y_val, color = treat, group = index),
            size = 1, alpha = .75) +  
  geom_rect(data = bg_color_data,
            aes(xmin = xmin, xmax = xmax, fill = factor(hotspot_label, 

                                                        levels = c("5/6", "22/23", "28/29", "2/3", "12/13", "24/25"),
                                                        labels = c("5/6 (17-35)", "22/23 (85-103)", "28/29 (109-127)", "2/3 (281-299)", "12/13 (321-335)", "24/25 (93-111)"),
                                                        ordered = T)),
            ymin = -Inf, ymax = Inf, alpha = .25) +
 scale_color_manual(name = "", values = group_colors[4:3]) +
  # scales::viridis_pal()(6)[c(1,3,5,2,4)]+
 # scale_fill_discrete(name = "", labels = 'Hot Spot') +
     scale_fill_manual("Hot Spot\n (Peptide Position Range)",
                      # values = c( "#440154FF" ,"#2A788EFF" ,"#7AD151FF", "#414487FF", "#22A884FF", "#FDE725FF"))+
                       values = c("yellow", "red", "skyblue", "green", "purple", "pink"))+
   
  scale_x_continuous("Peptide Position", 
                     breaks = c(1, 20, 40, 60, 80, 100, 120, 140, 159,
                                281, 300, 320, 339)) +
  scale_y_continuous("Participant", breaks = unique(ics_peptide_data_plot$treat_pubid_num), 
                   labels =  unique(ics_peptide_data_plot$treat_pubid) %>%
                     str_sub(end = str_locate(unique(ics_peptide_data_plot$treat_pubid),
                                              ':')[,'start'] - 1)
                   ) +
  coord_cartesian(ylim = c(2,33)) +
  facet_grid(parent_plot ~ peptide_group, scales = 'free', space = 'free_x') +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    ) +
  guides(color = guide_legend(reverse = TRUE),
         fill = guide_legend(nrow=2,byrow=TRUE))


 
hotspot_plot_withlabs

# pdf(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/hotspot_plot_withlabs.pdf',
#     width = 11, height = 5.5);
# hotspot_plot_withlabs
# dev.off()
# 
# png(file = '/home/cmahoney/Projects/feinberg_july/paper/Tcell_paper_figs/Tcell-figures/hotspot_plot_withlabs.png',
#     res = 100,pointsize = 3,
#     width = 800, height = 550);
# hotspot_plot_withlabs
# dev.off()



## ----Software-Session-Information, results="asis", message=FALSE, warning=kable_warnings----------------------------------------------------------------------------------
# load in rmarkdown to capture version number
if (any(installed.packages()[,1] == 'rmarkdown')) suppressWarnings(library(rmarkdown))

my_session_info <- VISCfunctions::get_session_info()

kable(
  my_session_info$platform_table, 
  format = output_type, 
  booktabs = TRUE, 
  linesep = "", 
  caption = "Reproducibility software session information"
  ) %>% 
  kable_styling(font_size = 5,
                    latex_options = c("hold_position", "repeat_header")
)


## ----Software-Package-Version-Information, results="asis", warning=kable_warnings-----------------------------------------------------------------------------------------
kable(
  my_session_info$packages_table, 
  format = output_type, booktabs = TRUE, 
  linesep = "", 
  caption = "Reproducibility software package version information"
  ) %>% 
  kable_styling(font_size = 5,
                    latex_options = c("hold_position","repeat_header")
)

