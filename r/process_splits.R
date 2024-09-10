#
# Purpose: Uses the survival library and tmerge function to split the processed data
# into time chunks for the Cox proportional hazards models
#

here::i_am("r/process_splits.R")

# libraries
library(readstata13)
library(dplyr)
library(tidyr)
library(forcats)
library(lubridate)
library(survival)
library(here)

# load data
load(here("data/processed_data.RData"))

end_date <- as_date("2023-03-31")

siren_cox <- siren_df_interim4 |>
    select(study_id, time, specimen_date, result, prev_var, last_pos, infection_date_1, infection_date_2, infection_date_3, vaccine_date4, gender, agegr, ethnicity, medical_group, staff_type, occupation_setting, patient_contact, imd, region, trust_code, covid_symptoms, household, time_since_pos, max_gap, true_naive) |>
    filter(result != 99, specimen_date <= end_date) |>
    mutate(
        event = case_when(
            specimen_date == infection_date_1 ~ 1,
            specimen_date == infection_date_2 ~ 2,
            specimen_date == infection_date_3 ~ 3,
            TRUE ~ 0
        ),
        vaccine = as.numeric(difftime(vaccine_date4, as_date("2022-09-12"), units = "weeks")),
        last_pos_wk = as.numeric(difftime(last_pos, as_date("2022-09-12"), units = "weeks")),
        ar = as.numeric(difftime(last_pos, as_date("2022-09-12"), units = "weeks")), # 6 weeks of dropout after infection_dates
        time = if_else(time == 0, 1 / 7, time),
        prev_var = fct_recode(prev_var, `No recorded previous infection` = "Naive")
    ) |>
    # filter to just the final event, with cohort re-entry
    arrange(-time) |>
    distinct(study_id, event, .keep_all = TRUE) |>
    arrange(study_id, time) |>
    filter(
        time > 0
    ) |>
    mutate(
        ar = if_else(is.na(ar), 0, ar),
        study_id = paste0(study_id, "_", event)
    )

# Split to generate time at risk
siren_split <- tmerge(
    data1 = siren_cox,
    data2 = siren_cox,
    id = study_id,
    tstop = time,
    event = event(time, event),
    # split the data by month
    monthyear = cumtdc(rep(2.71428571428571, count(siren_cox))),
    monthyear = cumtdc(rep(7.14285714285714, count(siren_cox))),
    monthyear = cumtdc(rep(11.4285714285714, count(siren_cox))),
    monthyear = cumtdc(rep(15.8571428571429, count(siren_cox))),
    monthyear = cumtdc(rep(20.2857142857143, count(siren_cox))),
    monthyear = cumtdc(rep(24.2857142857143, count(siren_cox))),
    #sep10 = tdc(rep(14.5, count(siren_cox))),
    eligible = tdc(ar),
    vax = tdc(vaccine),
    vax_l = cumtdc(vaccine),
    vax_l = cumtdc(vaccine + 8),
    vax_l = cumtdc(vaccine + 16),
    vax_m = cumtdc(vaccine),
    vax_m = cumtdc(vaccine + 12),
    msp = cumtdc(last_pos_wk + 104),
    msp = cumtdc(last_pos_wk + 52),
    msp = cumtdc(last_pos_wk + 26),
    msp = cumtdc(last_pos_wk)
) |>
    mutate(
        event = if_else(event %in% c(2, 3, 4), 1, event),
        vaccine_short = factor(vax, levels = c(0:1), labels = c("Waned third dose", "Fourth dose")),
        vaccine_2 = factor(vax_m, levels = c(0:2), labels = c("Waned third dose", "Fourth dose 0-3 months", "Fourth dose 3+ months")),
        vaccine = factor(vax_l, levels = c(0:3), labels = c("Waned third dose", "Fourth dose 0-2 months", "Fourth dose 2-4 months", "Fourth dose 4+ months")),
        msp = if_else(true_naive == 1 & msp == 0, 5, msp),
        months_since_pos = factor(msp, levels = c(0:5), labels = c("No evidence of infection", "0-6 months", "6-12 months", "1-2 years", "2+ years", "Confirmed naive")),
        months_since_pos = fct_relevel(months_since_pos, "2+ years"),
        monthyear = factor(monthyear,
            levels = c(0:7),
            labels = c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
        )
    ) |>
    filter(
        eligible == 1,
        tstop - tstart > 0.1
    )

# Save data
save(siren_split, siren_cox, file = here("data/split_data.RData"))
