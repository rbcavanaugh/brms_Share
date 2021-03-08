
############################################### Script to share ###########################################

# R version: 4.0.2. system info at bottom. 
library(tidyverse)
library(janitor)
library(brms)
library(performance)
library(patchwork)
library(RColorBrewer)
library(insight)
library(tidybayes)
library(ggridges)
library(here)
library(DHARMa)


############################################### DATA ###########################################
across <- read_csv(here("data", 'discourse_fake_data.csv')) # if you're not using the here package and r projects, I recommend you do so!!!!! otherwise, just run sd() from wherever you 

# fake data for treatment study
# pre post CIUs/words (proportion CIUs)
# severity
# research q: does severity moderate discourse improvements?

##################################### MODEL PROPORTION CIUS ############################################
across.binom <- brm(cius | trials(words) ~ time_point*severity.z + (time_point|participant),
                    data = across,
                    family = binomial(link = 'logit'), # proporiton data. critical to choose right distribution!!!!
                    iter = 3000,
                    control = list(adapt_delta = .9), # increase this to avoid divergent transitions
                    warmup = 1000, # these samples are discarded. 
                    chains = 2, # 3000 iterations - 1000 discarded = 2000x 2 chains = 4000 total samples. use 4 normally. 
                    cores = 2,
                    inits = 'random',
                    save_all_pars = TRUE) # you need this for some post prociessign sometimes. 

# in case this model doesn't run for you right now, you can load a pre-run one. 
#load(here('data', 'model1.RData'))

# mode checks
brms::pp_check(across.binom, nsamples = 200) # check the posterior distribution. dark line should be in lighter lines
summary(across.binom, prob = .9) # summary with 90% credible interval. 
conditional_effects(across.binom, probs = c(.05, .95)) # basic plots of marginal effects
# overdispersion test for binomial logit
preds = posterior_predict(across.binom)
preds = t(preds)
res = createDHARMa(simulatedResponse = preds,
                   fittedPredictedResponse = apply(preds, 1, median),
                   observedResponse = across.binom$data$cius, integerResponse = T)
testDispersion(res, alternative = 'greater') # some underdispersion oddly enough. ok for this example. 


################################## individual effect sizes #######################

# note! this only words right if you correctl specify random effect structures!
# need to include random interepts for participants and appropriate random slopes!
# if you use intercept only models, this will mischaracterize individual effects. 


# https://www.rensvandeschoot.com/tutorials/brms-started/
# to get individual estimates, we need to extract posterior samples. 
# first model data without cius
prediction_data2 <- across %>%
  select(participant, time_point, severity.z, words) %>%
  group_by(participant) %>%
  mutate(words = round(mean(words, na.rm = T)),0) # hold # of words constant in the prediction data. 

### fitted draws for existing data:
# 1. get posterior draws from the model
# 2. calculat ethe predicted percentCIUs at entry and exit for each draw
# 3. # take the difference
# 4. get the 90% highest density interval for each participant
# 5. done!
prop2 <- add_fitted_draws(across.binom, newdata = prediction_data2, re_formula = NULL, pred = 'value') %>% 
  ungroup() %>%
  mutate(.value = .value/words) %>%
  select(participant, time_point, value = .value, draw = .draw) %>%
  mutate(time_point = ifelse(time_point == 0, 'entry', 'exit')) %>%
  pivot_wider(names_from = time_point, values_from = value) %>%
  mutate(beta = exit - entry) %>%
  select(participant, beta) %>%
  ungroup() %>% 
  nest_by(participant) %>%
  summarize(qi = map(data, median_hdi, .width = .9)) %>%
  unnest(cols = c(qi))

# for plotting
# same stuff, except don't get the HDI, just keep all the draws and make individual density plots
ridges.plot.data <- add_fitted_draws(across.binom, newdata = prediction_data2, re_formula = NULL, pred = 'value') %>%
  ungroup() %>%
  mutate(.value = .value/words) %>%
  select(participant, time_point, value = .value, draw = .draw) %>%
  mutate(time_point = ifelse(time_point == 0, 'entry', 'exit')) %>%
  pivot_wider(names_from = time_point, values_from = value) %>%
  mutate(beta = exit - entry) %>%
  select(participant, beta) %>%
  ungroup() 

# to add severity back in
severity.df <- prediction_data2 %>%
  select(participant, severity.z)

# plot desnity for each participant. 
ridges.plot.data %>%
  left_join(severity.df, by = 'participant') %>%
  mutate(participant = as.factor(participant)) %>%
  group_by(participant) %>%
  ggplot(aes(x = beta, y = fct_reorder(participant, beta), fill = severity.z), color = 'grey') +
  geom_density_ridges(aes(height = ..density..)) +
  geom_vline(aes(xintercept = 0)) +
  scale_fill_viridis_c(direction = -1) +
  ylab('Participant') +
  xlab('Change in proportion CIU')



# sessionInfo()
# 
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] sjlabelled_1.1.6   DHARMa_0.3.3.0     here_0.1           knitr_1.30         ggridges_0.5.2     tidybayes_2.1.1   
# [7] insight_0.9.6      RColorBrewer_1.1-2 patchwork_1.0.1    performance_0.5.0  brms_2.13.5        Rcpp_1.0.5        
# [13] gghighlight_0.3.0  janitor_2.0.1      forcats_0.5.0      stringr_1.4.0      dplyr_1.0.2        purrr_0.3.4       
# [19] readr_1.3.1        tidyr_1.1.2        tibble_3.0.3       ggplot2_3.3.2      tidyverse_1.3.0   
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1         backports_1.1.10     plyr_1.8.6           igraph_1.2.5         splines_4.0.2       
# [6] svUnit_1.0.3         crosstalk_1.1.0.1    rstantools_2.1.1     inline_0.3.16        digest_0.6.25       
# [11] foreach_1.5.0        htmltools_0.5.0      rsconnect_0.8.16     fansi_0.4.1          magrittr_1.5        
# [16] modelr_0.1.8         RcppParallel_5.0.2   matrixStats_0.56.0   xts_0.12.1           prettyunits_1.1.1   
# [21] colorspace_1.4-1     blob_1.2.1           rvest_0.3.6          ggdist_2.2.0         haven_2.3.1         
# [26] xfun_0.18            callr_3.4.4          crayon_1.3.4         jsonlite_1.7.1       lme4_1.1-23         
# [31] zoo_1.8-8            iterators_1.0.12     glue_1.4.2           gtable_0.3.0         emmeans_1.5.0       
# [36] V8_3.2.0             distributional_0.2.0 pkgbuild_1.1.0       rstan_2.21.2         abind_1.4-5         
# [41] scales_1.1.1         mvtnorm_1.1-1        DBI_1.1.0            miniUI_0.1.1.1       viridisLite_0.3.0   
# [46] xtable_1.8-4         stats4_4.0.2         StanHeaders_2.21.0-6 DT_0.15              htmlwidgets_1.5.1   
# [51] httr_1.4.2           threejs_0.3.3        arrayhelpers_1.1-0   ellipsis_0.3.1       pkgconfig_2.0.3     
# [56] loo_2.3.1            farver_2.0.3         dbplyr_1.4.4         utf8_1.1.4           labeling_0.3        
# [61] tidyselect_1.1.0     rlang_0.4.7          reshape2_1.4.4       later_1.1.0.1        munsell_0.5.0       
# [66] cellranger_1.1.0     tools_4.0.2          cli_2.0.2            generics_0.0.2       broom_0.7.0         
# [71] evaluate_0.14        fastmap_1.0.1        processx_3.4.4       fs_1.5.0             packrat_0.5.0       
# [76] nlme_3.1-149         mime_0.9             xml2_1.3.2           compiler_4.0.2       bayesplot_1.7.2     
# [81] shinythemes_1.1.2    rstudioapi_0.11      curl_4.3             reprex_0.3.0         statmod_1.4.34      
# [86] stringi_1.5.3        ps_1.3.4             Brobdingnag_1.2-6    lattice_0.20-41      Matrix_1.2-18       
# [91] nloptr_1.2.2.2       markdown_1.1         shinyjs_2.0.0        vctrs_0.3.4          pillar_1.4.6        
# [96] lifecycle_0.2.0      bridgesampling_1.0-0 estimability_1.3     httpuv_1.5.4         R6_2.4.1            
# [101] promises_1.1.1       gridExtra_2.3        codetools_0.2-16     boot_1.3-25          colourpicker_1.0    
# [106] MASS_7.3-53          gtools_3.8.2         assertthat_0.2.1     rprojroot_1.3-2      withr_2.3.0         
# [111] shinystan_2.5.0      bayestestR_0.7.2     parallel_4.0.2       hms_0.5.3            grid_4.0.2          
# [116] minqa_1.2.4          coda_0.19-3          snakecase_0.11.0     shiny_1.5.0          lubridate_1.7.9     
# [121] base64enc_0.1-3      dygraphs_1.1.1.6   
