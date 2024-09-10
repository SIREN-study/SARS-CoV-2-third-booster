#
# Purpose: Runs the Cox proportional hazards models on the split SIREN data
#

here::i_am("r/cox_models.R")

library(dplyr)
library(finalfit)
library(here)
library(survival)

# load data
load(here("data/split_data.RData"))

dependent <- "Surv(tstart, tstop, event)"

# cox_model_1
explanatory_multi <- c("vaccine_short", "gender", "strata(occupation_setting)", "strata(agegr)", "household", "strata(region)", "cluster(trust_code)")
cox_model_1 <- siren_split |>
    coxphmulti(dependent, explanatory_multi)

# cox_model_2
explanatory_multi <- c("vaccine", "gender", "strata(occupation_setting)", "strata(agegr)", "household", "strata(region)", "cluster(trust_code)")
cox_model_2 <- siren_split |>
    coxphmulti(dependent, explanatory_multi)

# cox_model_3
explanatory_multi <- c("months_since_pos*vaccine_short", "gender", "strata(occupation_setting)", "strata(agegr)", "household", "strata(region)", "cluster(trust_code)")
cox_model_3 <- siren_split |>
    filter(months_since_pos != "No evidence of infection") |>
    coxphmulti(dependent, explanatory_multi)

# cox_model_4
explanatory_multi <- c("months_since_pos", "months_since_pos:vaccine_short", "gender", "strata(occupation_setting)", "strata(agegr)", "household", "strata(region)", "cluster(trust_code)")
cox_model_4 <- siren_split |>
    filter(months_since_pos != "No evidence of infection") |>
    coxphmulti(dependent, explanatory_multi)

cox_model_results <- grep("cox_model", ls(), value = TRUE)

save(list = cox_model_results, file = here("data/cox_models.RData"))