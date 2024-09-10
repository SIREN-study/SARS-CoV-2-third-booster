#
# Purpose: Runs the multi-state models for the SIREN data
#

here::i_am("r/msm_models.R")

library(dplyr)
library(future.callr)
library(future)
library(here)
library(msm)

plan(callr)

load(here("data/siren_processed.RData"))

# remove individuals with only one pcr result
study_ids <- siren_df_interim4 |>
    select(study_id, state) |>
    filter(state != 99) |>
    add_count(study_id) |>
    filter(n > 1) |>
    distinct(study_id)

# prepare two datasets with and without covid_symptom information
siren_df <- siren_df_interim4 |>
    filter(
        study_id %in% study_ids$study_id
    ) |>
    mutate(
        monthyear = fct_drop(monthyear),
        months_since_pos = fct_drop(months_since_pos)
    )

siren_df_sym <- siren_df_interim4 |>
    filter(
        study_id %in% study_ids$study_id
    ) |>
    mutate(
        covid_symptoms = if_else(covid_symptoms == "Other symptoms", "Asymptomatic", covid_symptoms),
        state = if_else(state == 2 & covid_symptoms == "Asymptomatic" & !is.na(covid_symptoms), 3, state),
        monthyear = fct_drop(monthyear),
        months_since_pos = fct_drop(months_since_pos)
    )

# prepare a constraint list for the models which include covid_symptoms
constraint_list <- list(
    `monthyearOct 2022` = c(1, 1, 2, 2),
    `monthyearNov 2022` = c(1, 1, 2, 2),
    `monthyearDec 2022` = c(1, 1, 2, 2),
    `monthyearJan 2023` = c(1, 1, 2, 2),
    `monthyearFeb 2023` = c(1, 1, 2, 2),
    `monthyearMar 2023` = c(1, 1, 2, 2),
    `regionEast of England` = c(1, 1),
    `regionLondon` = c(1, 1),
    `regionNorth East` = c(1, 1),
    `regionNorth West` = c(1, 1),
    `regionNorthern Ireland` = c(1, 1),
    `regionScotland` = c(1, 1),
    `regionSouth East` = c(1, 1),
    `regionSouth West` = c(1, 1),
    `regionWales` = c(1, 1),
    `regionWest Midlands` = c(1, 1),
    `regionYorkshire and the Humber` = c(1, 1),
    `agegr<25` = c(1, 1, 2, 2),
    `agegr25-34` = c(1, 1, 2, 2),
    `agegr35-44` = c(1, 1, 2, 2),
    `agegr55-64` = c(1, 1, 2, 2),
    `agegr65+` = c(1, 1, 2, 2),
    `householdLives_alone` = c(1, 1),
    `householdLives_with_others_with_child` = c(1, 1),
    `genderFemale` = c(1, 1),
    `occupation_settingPatient_facing_(non-clinical)` = c(1, 1),
    `occupation_settingAmbulance/Emergency_Department/Inpatient_Wards` = c(1, 1),
    `occupation_settingMaternity/Labour_Ward` = c(1, 1),
    `occupation_settingIntensive_Care` = c(1, 1),
    `occupation_settingTheatres` = c(1, 1),
    `occupation_settingOther` = c(1, 1),
    `occupation_settingOffice` = c(1, 1)
)

# MSM model 1
msm_model_1 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df,
    qmatrix = rbind(c(0, 0.1), c(0.1, 0)),
    covariates = list(
        "1-2" = ~ monthyear + vaccine_short + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ vaccine_short + agegr
    ),
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 40000)
))

# MSM model 2
msm_model_2 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df,
    qmatrix = rbind(c(0, 0.1), c(0.1, 0)),
    covariates = list(
        "1-2" = ~ monthyear + vaccine + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ vaccine + agegr
    ),
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 40000)
))

# MSM model 3
msm_model_3 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df,
    qmatrix = rbind(c(0, 0.1), c(0.1, 0)),
    covariates = list(
        "1-2" = ~ monthyear + (months_since_pos * vaccine_short) + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ (months_since_pos * vaccine_short) + agegr
    ),
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 40000)
))

# MSM model 4
msm_model_4 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df,
    qmatrix = rbind(c(0, 0.1), c(0.1, 0)),
    covariates = list(
        "1-2" = ~ monthyear + months_since_pos + months_since_pos:vaccine_short + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ months_since_pos + months_since_pos:vaccine_short + agegr
    ),
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 40000)
))

# MSM model 5
msm_model_5 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df_sym,
    qmatrix = rbind(c(0, 0.1, 0.1), c(0.1, 0, 0), c(0.1, 0, 0)),
    covariates = list(
        "1-2" = ~ monthyear + vaccine_short + region + agegr + household + gender + occupation_setting,
        "1-3" = ~ monthyear + vaccine_short + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ monthyear + vaccine_short + agegr,
        "3-1" = ~ monthyear + vaccine_short + agegr
    ),
    constraint = constraint_list,
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 70000)
))

# MSM model 6
msm_model_6 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df_sym,
    qmatrix = rbind(c(0, 0.1, 0.1), c(0.1, 0, 0), c(0.1, 0, 0)),
    covariates = list(
        "1-2" = ~ monthyear + vaccine + region + agegr + household + gender + occupation_setting,
        "1-3" = ~ monthyear + vaccine + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ monthyear + vaccine + agegr,
        "3-1" = ~ monthyear + vaccine + agegr
    ),
    constraint = constraint_list,
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 40000)
))

# MSM model 7
msm_model_7 <- future(msm::msm(
    formula = state ~ time,
    subject = study_id,
    data = siren_df_sym,
    qmatrix = rbind(c(0, 0.1, 0.1), c(0.1, 0, 0), c(0.1, 0, 0)),
    covariates = list(
        "1-2" = ~ monthyear + (months_since_pos * vaccine_short) + region + agegr + household + gender + occupation_setting,
        "1-3" = ~ monthyear + (months_since_pos * vaccine_short) + region + agegr + household + gender + occupation_setting,
        "2-1" = ~ monthyear + (months_since_pos * vaccine_short) + agegr,
        "3-1" = ~ monthyear + (months_since_pos * vaccine_short) + agegr
    ),
    constraint = constraint_list,
    censor = 99,
    control = list(maxit = 100000, reltol = 1e-16, fnscale = 40000)
))

msm_model_1 <- value(msm_model_1)
msm_model_2 <- value(msm_model_2)
msm_model_3 <- value(msm_model_3)
msm_model_4 <- value(msm_model_4)
msm_model_5 <- value(msm_model_5)
msm_model_6 <- value(msm_model_6)
msm_model_7 <- value(msm_model_7)

msm_model_results <- grep("msm_model", ls(), value = TRUE)

save(list = msm_model_results, file = here("data/msm_models.RData"))
