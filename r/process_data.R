#
# Purpose: Imports and cleans the raw SIREN data and generates datasets
#

here::i_am("r/process_data.R")

library(dplyr)
library(forcats)
library(here)
library(janitor)
library(lubridate)
library(readstata13)
library(readxl)
library(tidyr)

# participant dataset
pt <- read.dta13("~/coviddata/PT_20230118.dta", convert.factors = TRUE, nonint.factors = TRUE) |>
    clean_names()

# pcr dataset
pcr <- read.dta13("~/coviddata/PCR_20230421_clean.dta") |>
    clean_names()

# follow up dataset
fu <- read.dta13("~/coviddata/FollowUp_20230421.dta") |>
    clean_names()

# antibody dataset
ab <- read.dta13("~/coviddata/AB_20230421.dta") |>
    clean_names()

# true naives dataset - manually updated to only contain those with no previous infection
true_naives <- read_xlsx(here("data/interim4_infection_history_updated.xlsx")) |>
    clean_names() |>
    filter(prev_var == "Naive") |>
    # replace missing values with 0
    mutate(naive = if_else(is.na(naive), 0, naive)) |>
    select(study_id, true_naive = naive)

# start and end dates
start_date <- as_date("2022-09-12")
end_date <- as_date("2023-03-31")

window_length <- 42
symptom_window <- 14

# generate a siren_cohort dataset with all participants who have valid information
siren_cohort <- pt |>
    filter(
        (vaccine_date1 < vaccine_date2 | is.na(vaccine_date2) | is.na(vaccine_date1)),
        (vaccine_date1 < vaccine_date3 | is.na(vaccine_date3) | is.na(vaccine_date1)),
        (vaccine_date2 < vaccine_date3 | is.na(vaccine_date3) | is.na(vaccine_date2)),
        !gender %in% c("Non_binary", "Prefer_not"),
        ethnicity != "Prefer_not",
        study_end_date > start_date
    ) |>
    mutate(
        agegr = factor(agegr, labels = c("<25", "25-34", "35-44", "45-54", "55-64", "65+")),
        agegr = fct_relevel(agegr, "45-54"),
        gender = factor(gender, labels = c("Male", "Female")),
        ethnicity = factor(ethnicity),
        region = factor(region),
        occupation_setting = if_else(is.na(occupation_setting_ext), occupation_setting, occupation_setting_ext),
        occupation_setting = if_else(is.na(occupation_setting), "Other", as.character(occupation_setting)),
        occupation_setting = factor(occupation_setting),
        occupation_setting = fct_relevel(occupation_setting, "Outpatient"),
        staff_type = as.double(staff_type),
        staff_type = if_else(staff_type == 11 & patient_contact == "No", 12, staff_type), # incorrectly coded in source data
        staff_type = factor(staff_type, levels = 1:12, labels = c(
            "Administrative/Executive (office based)", "Doctor", "Nursing", "Healthcare Assistant", "Midwife", "Healthcare Scientist", "Pharmacist",
            "Physiotherapist/Occupational Therapist/SALT", "Student (Medical/Nursing/Midwifery/Other)", "Estates/Porters/Security",
            "Other (non patient facing)", "Other (patient facing)"
        )),
        imd = factor(imd_5, labels = c("Most deprived (1)", "Deprivation 2", "Deprivation 3", "Deprivation 4", "Least deprived (5)")),
        household = factor(household),
        household = fct_relevel(household, "Lives_with_others_nochild")
    )

# calculate infection history using ab dataset
# anti-N positive confirms previous infection
# anti-S positive confirms previous infection and/or response to vaccination
antibody_data <- ab |>
    group_by(study_id) |>
    arrange(study_id, spec_date_ab) |> # sort by specimen date
    mutate(
        result_ab = case_when(
            n_result == "#+" ~ 1,
            s_result == "#+" & (spec_date_ab < vaccine_date1 | is.na(vaccine_date1)) ~ 1,
            result_sgss == "Positive" ~ 1,
            .default = 0
        ),
        result_ab = case_when(
            result_ab == 1 & lag(result_ab) == 0 & !is.na(lag(result_ab)) ~ 1, # new infection
            result_ab == 1 ~ 11, # could be legacy infection
            .default = 0
        ),
        prev_neg = if_else(result_ab == 1, lag(spec_date_ab), NA_Date_),
        result = result_ab,
        specimen_date = spec_date_ab,
        study_id = study_id
    ) |>
    ungroup()

# define a function to calculate the infection window
calculate_window <- function(df) {
    for (val in 1:7) {
        df <- df |>
            mutate(window_init = if_else(window_init == 9 & lag(window_init) == 1 & specimen_date <= lag(window_end), 0, window_init)) |>
            filter(window_init != 0)
    }
    return(df)
}

# last test (regardless of result) up to three months before the start date
last_test <- pcr |>
    bind_rows(antibody_data) |>
    filter(
        study_id %in% siren_cohort$study_id,
        specimen_date < start_date,
        specimen_date >= as_date("2021-06-01")
    ) |>
    group_by(study_id) |>
    arrange(specimen_date) |>
    mutate(
        # calculate the gap between successive tests
        gap = lead(specimen_date, default = start_date) - specimen_date
    ) |>
    summarise(max_gap = max(gap, na.rm = FALSE)) |> # max gap between tests
    ungroup() |>
    left_join(
        pcr |>
            bind_rows(antibody_data) |>
            filter(
                study_id %in% siren_cohort$study_id,
                specimen_date < start_date
            )
    ) |>
    group_by(study_id) |>
    arrange(specimen_date) |>
    slice_tail() |>
    ungroup() |>
    select(study_id, last_test = specimen_date, max_gap)

# generate infection windows for every infection
infection_window <- pcr |>
    bind_rows(antibody_data) |>
    filter(
        study_id %in% siren_cohort$study_id,
        result %in% c(1, 11)
    ) |>
    mutate(
        window_start = as_date(specimen_date),
        window_end = window_start + window_length
    ) |>
    group_by(study_id) |>
    arrange(study_id, specimen_date) |>
    # remove antibody positives when there is evidence of a prior PCR positive
    # in this case we cannot be certain that the antibody positive represents a new infection
    mutate(
        result = case_when(
            result_ab == 1 & result == 1 & lag(result) == 1 & !is.na(lag(result)) & lag(specimen_date) < prev_neg ~ 1, # negative in interim
            result_ab == 1 & result == 1 & lag(result) == 1 & !is.na(lag(result)) & lag(specimen_date) > prev_neg ~ 0, # unlikely a new positive
            result == 11 & is.na(lag(result)) ~ 1, # first record of infection
            result == 11 & !is.na(lag(result)) ~ 0, # legacy infection, unlikely a new positive
            .default = result
        )
    ) |>
    filter(result == 1) |>
    distinct(study_id, specimen_date, result, window_start, window_end) |> # remove duplicates
    # first window (infection)
    mutate(window_init = if_else(row_number() == 1, 1, 9)) |>
    calculate_window() |>
    # second window (reinfection)
    mutate(window_init = if_else(window_init == 9 & lag(window_init) == 1 & specimen_date > lag(window_end), 1, window_init)) |>
    calculate_window() |>
    # third window (reinfection)
    mutate(window_init = if_else(window_init == 9 & lag(window_init) == 1 & specimen_date > lag(window_end), 1, window_init)) |>
    calculate_window() |>
    # fourth window (reinfection)
    mutate(window_init = if_else(window_init == 9 & lag(window_init) == 1 & specimen_date > lag(window_end), 1, window_init)) |>
    calculate_window() |>
    # fifth window (reinfection)
    mutate(window_init = if_else(window_init == 9 & lag(window_init) == 1 & specimen_date > lag(window_end), 1, window_init)) |>
    calculate_window() |>
    # sixth window (reinfection)
    mutate(window_init = if_else(window_init == 9 & lag(window_init) == 1 & specimen_date > lag(window_end), 1, window_init)) |>
    calculate_window() |>
    ungroup()

infection_window_wide <- infection_window |>
    filter(
        specimen_date >= start_date - window_length,
        specimen_date <= end_date
    ) |>
    group_by(study_id) |>
    mutate(n = row_number()) |>
    select(study_id, window_start, window_end, n) |>
    pivot_wider(
        names_from = n,
        values_from = c(window_start, window_end)
    )

# prior variant periods
# label variant by episode_date
var_type <- infection_window |>
    transmute(
        study_id,
        wt = if_else(specimen_date <= lubridate::as_date("2020-12-15"), 1, 0),
        alpha = if_else(specimen_date >= lubridate::as_date("2020-12-16") & specimen_date <= lubridate::as_date("2021-04-30"), 1, 0),
        delta = if_else(specimen_date >= lubridate::as_date("2021-05-01") & specimen_date <= lubridate::as_date("2021-12-15"), 1, 0),
        ba1 = if_else(specimen_date >= lubridate::as_date("2021-12-16") & specimen_date <= lubridate::as_date("2022-03-21"), 1, 0),
        ba2 = if_else(specimen_date >= lubridate::as_date("2022-03-22") & specimen_date <= lubridate::as_date("2022-05-31"), 1, 0),
        ba45 = if_else(specimen_date >= lubridate::as_date("2022-06-01"), 1, 0)
    ) |>
    group_by(study_id) |>
    summarise(across(everything(), ~ 1 * (sum(.) != 0))) |> # to turn sum into a binary indicator
    ungroup()

# clean the pcr dataset
pcr_clean <- pcr |>
    filter(
        study_id %in% siren_cohort$study_id,
        !is.na(specimen_date),
        !is.na(result),
        specimen_date >= start_date - window_length,
        specimen_date <= end_date
    ) |>
    left_join(infection_window_wide, by = "study_id") |>
    mutate(
        specimen_date = as_date(specimen_date),
        group = case_when(
            specimen_date >= window_start_1 & specimen_date < window_end_1 ~ "w1",
            specimen_date >= window_start_2 & specimen_date < window_end_2 ~ "w2",
            specimen_date >= window_start_3 & specimen_date < window_end_3 ~ "w3"
        )
    ) |>
    arrange(study_id, specimen_date) |>
    group_by(study_id) |>
    mutate(
        state = case_when(
            result == 1 ~ 2, # result == 1 => positive
            # positive both sides
            lag(result) == 1 & result == 2 & lead(result) == 1 & lag(group) == lead(group) ~ 2,
            # positive with two gaps
            lag(result) == 1 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 1 & lag(group) == lead(group, n = 2) ~ 2,
            lag(result, n = 2) == 1 & lag(result) == 2 & result == 2 & lead(result) == 1 & lag(group, n = 2) == lead(group) ~ 2,
            # positive with three gaps
            lag(result, n = 3) == 1 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 1 & lag(group, n = 3) == lead(group) ~ 2,
            lag(result, n = 2) == 1 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 1 & lag(group, n = 2) == lead(group, n = 2) ~ 2,
            lag(result) == 1 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 1 & lag(group) == lead(group, n = 3) ~ 2,
            # positive with four gaps
            lag(result, n = 4) == 1 & lag(result, n = 3) == 2 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 1 & lag(group, n = 4) == lead(group) ~ 2,
            lag(result, n = 3) == 1 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 1 & lag(group, n = 3) == lead(group, n = 2) ~ 2,
            lag(result, n = 2) == 1 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 1 & lag(group, n = 2) == lead(group, n = 3) ~ 2,
            lag(result) == 1 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 2 & lead(result, n = 4) == 1 & lag(group) == lead(group, n = 4) ~ 2,
            # positive with five gaps
            lag(result, n = 5) == 1 & lag(result, n = 4) == 2 & lag(result, n = 3) == 2 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 1 & lag(group, n = 5) == lead(group) ~ 2,
            lag(result, n = 4) == 1 & lag(result, n = 3) == 2 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 1 & lag(group, n = 4) == lead(group, n = 2) ~ 2,
            lag(result, n = 3) == 1 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 1 & lag(group, n = 3) == lead(group, n = 3) ~ 2,
            lag(result, n = 2) == 1 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 2 & lead(result, n = 4) == 1 & lag(group, n = 2) == lead(group, n = 4) ~ 2,
            lag(result) == 1 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 2 & lead(result, n = 4) == 2 & lead(result, n = 5) == 1 & lag(group) == lead(group, n = 5) ~ 2,
            # positive with six gaps
            lag(result, n = 6) == 1 & lag(result, n = 5) == 2 & lag(result, n = 4) == 2 & lag(result, n = 3) == 2 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 1 & lag(group, n = 6) == lead(group) ~ 2,
            lag(result, n = 5) == 1 & lag(result, n = 4) == 2 & lag(result, n = 3) == 2 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 1 & lag(group, n = 5) == lead(group, n = 2) ~ 2,
            lag(result, n = 4) == 1 & lag(result, n = 3) == 2 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 1 & lag(group, n = 4) == lead(group, n = 3) ~ 2,
            lag(result, n = 3) == 1 & lag(result, n = 2) == 2 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 2 & lead(result, n = 4) == 1 & lag(group, n = 3) == lead(group, n = 4) ~ 2,
            lag(result, n = 2) == 1 & lag(result) == 2 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 2 & lead(result, n = 4) == 2 & lead(result, n = 5) == 1 & lag(group, n = 2) == lead(group, n = 5) ~ 2,
            lag(result) == 1 & result == 2 & lead(result) == 2 & lead(result, n = 2) == 2 & lead(result, n = 3) == 2 & lead(result, n = 4) == 2 & lead(result, n = 5) == 2 & lead(result, n = 6) == 1 & lag(group) == lead(group, n = 6) ~ 2,
            result == 2 ~ 1 # result == 2 => negative
        ),
        # time of first positive within the current episode
        episode_start = if_else(state == 2 & (lag(state) != 2 | is.na(lag(state))), specimen_date, NA_Date_),
        # time of first negative within the current episode
        episode_end = if_else(state == 2 & lead(state) == 1 & !is.na(lead(state)), lead(specimen_date), NA_Date_)
    ) |>
    group_by(study_id, state) |>
    fill(episode_start, .direction = "down") |>
    fill(episode_end, .direction = "up") |>
    ungroup() |>
    filter(!is.na(state))

# generate censored observations states (for the piecewise assumption)
pcr_pw <- expand_grid(
    siren_cohort |> filter(study_id %in% pcr_clean$study_id),
    # by not including the first month we have (effective) delayed entry into the cohort, at the date of the first PCR test
    specimen_date = c(
        "2022-10-01", "2022-11-01", "2022-12-01",
        "2023-01-01", "2023-02-01", "2023-03-01"
    )
) |>
    mutate(
        result = 99,
        state = 99,
        specimen_date = as_date(specimen_date)
    )

# last_pos dates to be added to the cleaned PCR data
lk_pos <- pcr |>
    bind_rows(infection_window) |>
    filter(
        study_id %in% siren_cohort$study_id,
        !is.na(specimen_date),
        !is.na(result),
        specimen_date < start_date,
        result == 1
    ) |>
    group_by(study_id) |>
    arrange(study_id, specimen_date) |>
    rename(prev_last_pos = specimen_date) |>
    # keep the most recent observation
    slice_tail(n = 1) |>
    ungroup() |>
    select(study_id, prev_last_pos)

# add the unobserved states and calculate the last_pos and time_since_pos dates
siren_df <- siren_cohort |>
    filter(study_id %in% pcr_clean$study_id) |>
    left_join(pcr_clean, by = "study_id") |>
    bind_rows(pcr_pw) |>
    arrange(study_id, specimen_date, result) |>
    distinct(study_id, specimen_date, .keep_all = TRUE) |>
    group_by(study_id) |>
    mutate(
        last_pos = lag(episode_start),
        # twice in case positive episodes carry across more than one month
        episode_start = if_else(state == 99, lag(episode_start), episode_start),
        episode_start = if_else(state == 99, lag(episode_start), episode_start),
        episode_end = if_else(state == 99, lag(episode_end), episode_end),
        episode_end = if_else(state == 99, lag(episode_end), episode_end)
    ) |>
    fill(last_pos, .direction = "down") |>
    ungroup() |>
    # add the pre-study last_pos dates
    left_join(lk_pos, by = "study_id") |>
    mutate(
        # update last_pos with the pre-study last_pos dates
        last_pos = if_else(is.na(last_pos), prev_last_pos, last_pos),
        # time since last positive in days and months
        time_since_pos = difftime(specimen_date, last_pos, units = "days"),
        week = lubridate::week(specimen_date),
        month = lubridate::month(specimen_date, label = TRUE),
        year = lubridate::year(specimen_date),
        monthyear = paste0(month, " ", year),

        # time of test, relative to the individual since start_date
        days = difftime(specimen_date, start_date, units = "days"), # day scale
        time = difftime(specimen_date, start_date, units = "weeks"), # week scale
        time = as.numeric(time),

        # waned dose 3 vaccination defined as >24 weeks (6 months) post 3rd dose
        eligible = case_when(
            vaccine_date3 + 168 <= specimen_date ~ 1,
            vaccine_date3 <= specimen_date ~ 2,
            TRUE ~ 0
        ),

        # vaccine regimen
        regimen = case_when(
            vaccine_name1 == "Pfizer-BioNTech" & vaccine_name2 == "Pfizer-BioNTech" & vaccine_name3 == "Pfizer-BioNTech" ~ "3PF",
            vaccine_name1 == "Pfizer-BioNTech" & vaccine_name2 == "Pfizer-BioNTech" & vaccine_name3 == "Moderna" ~ "2PF/1MO",
            vaccine_name1 == "Oxford-AstraZeneca" & vaccine_name2 == "Oxford-AstraZeneca" & vaccine_name3 == "Pfizer-BioNTech" ~ "2AZ/1PF",
            vaccine_name1 == "Oxford-AstraZeneca" & vaccine_name2 == "Oxford-AstraZeneca" & vaccine_name3 == "Moderna" ~ "2AZ/1MO",
            vaccine_name1 == "Oxford-AstraZeneca" & vaccine_name2 == "Oxford-AstraZeneca" & vaccine_name3 == "Oxford-AstraZeneca" ~ "3AZ",
            TRUE ~ NA_character_
        )
    )

# record first (and subsequent) ba.3/4 infections
# only amongst waned third dose individuals
infection_dates <- siren_df |>
    filter(
        eligible == 1,
        specimen_date >= start_date,
        !is.na(regimen),
        specimen_date == episode_start
    ) |>
    group_by(study_id) |>
    mutate(n = row_number()) |>
    ungroup() |>
    select(study_id, specimen_date, n) |>
    pivot_wider(id_cols = study_id, names_from = n, values_from = specimen_date, names_prefix = "infection_date_")

# symptom status
# must have experienced symptoms within 14 days of first positive pcr test to be considered symptomatic
# aymptomatic must have no reported symptoms and first positive pcr test within reporting period
symptoms <- fu |>
    filter(
        study_id %in% siren_df$study_id
    ) |>
    mutate(
        covid_symptoms = rowSums(across(cough_fu:itchy_red_patches_fu), na.rm = TRUE),
        covid_symptoms = if_else(covid_symptoms > 0 | !is.na(sx_onset_fu), 1, 0),
        covid_symptoms = if_else((!is.na(sorethroat_fu) | !is.na(cough_fu) | !is.na(fever_fu) | !is.na(anosmia_fu) | !is.na(dygeusia_fu)), 2, covid_symptoms),
        symptom_start = case_when(
            covid_symptoms > 0 & !is.na(sx_onset_fu) ~ sx_onset_fu - symptom_window,
            covid_symptoms > 0 & is.na(sx_onset_fu) ~ date_survey_fu - symptom_window,
            TRUE ~ rep_start
        ),
        symptom_end = case_when(
            covid_symptoms > 0 & !is.na(sx_onset_fu) ~ sx_onset_fu + symptom_window,
            covid_symptoms > 0 & is.na(sx_onset_fu) ~ date_survey_fu + symptom_window,
            TRUE ~ rep_end
        )
    ) |>
    full_join(
        infection_dates |>
            pivot_longer(-study_id, values_to = "specimen_date") |>
            select(-name) |>
            filter(!is.na(specimen_date)),
        by = "study_id",
        relationship = "many-to-many"
    ) |>
    filter(
        specimen_date > symptom_start,
        specimen_date < symptom_end
    ) |>
    arrange(study_id, specimen_date, -covid_symptoms) |>
    distinct(study_id, specimen_date, .keep_all = TRUE) |>
    select(study_id, specimen_date, covid_symptoms)

siren_df_interim4 <- siren_df |>
    filter(
        eligible == 1,
        !is.na(regimen),
        specimen_date >= start_date
    ) |>
    left_join(infection_dates, by = "study_id") |>
    left_join(var_type, by = "study_id") |>
    left_join(symptoms, by = c("study_id", "specimen_date")) |>
    left_join(last_test, by = "study_id") |> # add last test date
    # add covid symptoms for the entire episode by joining on last_pos
    left_join(symptoms, by = c("study_id", "last_pos" = "specimen_date")) |>
    left_join(true_naives, by = "study_id") |>
    mutate(
        monthyear = factor(monthyear, levels = c(
            "Sep 2022", "Oct 2022", "Nov 2022",
            "Dec 2022", "Jan 2023", "Feb 2023",
            "Mar 2023"
        )),
        monthyear = fct_relevel(monthyear, "Sep 2022"),
        state = if_else(state == 3, 1, state), # tbc
        prev_var = case_when(
            infection_date_1 < specimen_date & (infection_date_1 < episode_start | is.na(episode_start)) & !is.na(infection_date_1) ~ "BA.4/5",
            (wt == 1 | alpha == 1 | delta == 1) & (ba1 == 1 | ba2 == 1) ~ "BA.1/2 + WT/Alpha/Delta",
            (ba1 == 1 | ba2 == 1) ~ "BA.1/2",
            (wt == 1 | alpha == 1 | delta == 1) ~ "WT/Alpha/Delta",
            TRUE ~ "Naive"
        ),
        prev_var_long = case_when(
            infection_date_1 < specimen_date & (infection_date_1 < episode_start | is.na(episode_start)) & !is.na(infection_date_1) ~ "BA.4/5",
            (delta == 1) & (ba1 == 1 | ba2 == 1) ~ "BA.1/2 + Delta",
            (alpha == 1) & (ba1 == 1 | ba2 == 1) ~ "BA.1/2 + Alpha",
            (wt == 1) & (ba1 == 1 | ba2 == 1) ~ "BA.1/2 + WT",
            (ba1 == 1 | ba2 == 1) ~ "BA.1/2",
            delta == 1 ~ "Delta",
            alpha == 1 ~ "Alpha",
            (wt == 1 | alpha == 1) ~ "WT",
            TRUE ~ "Naive"
        ),
        vaccine = case_when(
            vaccine_date4 + 112 <= specimen_date ~ "Fourth dose 4+ months",
            vaccine_date4 + 56 <= specimen_date ~ "Fourth dose 2-4 months",
            vaccine_date4 <= specimen_date ~ "Fourth dose 0-2 months",
            TRUE ~ "Waned third dose"
        ),
        vaccine_short = case_when(
            vaccine_date4 <= specimen_date ~ "Fourth dose",
            TRUE ~ "Waned third dose"
        ),
        prev_var = factor(prev_var, levels = c("Naive", "WT/Alpha/Delta", "BA.1/2", "BA.1/2 + WT/Alpha/Delta", "BA.4/5")),
        prev_var_long = factor(prev_var_long, levels = c("Naive", "WT", "Alpha", "Delta", "BA.1/2", "BA.1/2 + WT", "BA.1/2 + Alpha", "BA.1/2 + Delta", "BA.4/5")),
        # create a variable which splits the time_since_pos variable into several groups
        months_since_pos = case_when(
            time_since_pos <= 182 ~ "0-6 months",
            time_since_pos <= 364 ~ "6-12 months",
            time_since_pos <= 728 ~ "1-2 years",
            time_since_pos > 728 ~ "2+ years",
            true_naive == 1 ~ "Confirmed naive",
            TRUE ~ "No evidence of infection"
        ),
        months_since_pos = factor(months_since_pos, levels = c("2+ years", "1-2 years", "6-12 months", "0-6 months", "No evidence of infection", "Confirmed naive")),
        vaccine = factor(vaccine, levels = c("Waned third dose", "Fourth dose 0-2 months", "Fourth dose 2-4 months", "Fourth dose 4+ months")),
        vaccine_2 = factor(vaccine_2, levels = c("Waned third dose", "Fourth dose 0-3 months", "Fourth dose 3+ months")),
        vaccine_short = factor(vaccine_short, levels = c("Waned third dose", "Fourth dose")),
        # add covid symptoms for the entire episode
        covid_symptoms = covid_symptoms.x,
        covid_symptoms = if_else(is.na(covid_symptoms) & state == 2 & !is.na(covid_symptoms.y) &
            (specimen_date != infection_date_1 | is.na(infection_date_1)) &
            (specimen_date != infection_date_2 | is.na(infection_date_2)), covid_symptoms.y, covid_symptoms),
        covid_symptoms = factor(covid_symptoms, levels = c(0, 1, 2), labels = c("Asymptomatic", "Other symptoms", "COVID symptoms")),
        covid_symptoms = if_else(covid_symptoms == "Other symptoms", "Asymptomatic", covid_symptoms)
    ) |>
    select(-covid_symptoms.x, -covid_symptoms.y)

# remove individuals with zero pcr results
study_ids <- siren_df_interim4 |>
    select(study_id, state) |>
    filter(state != 99) |>
    distinct()

siren_df_interim4 <- siren_df_interim4 |> filter(study_id %in% study_ids$study_id)

save(siren_cohort, siren_df_interim4, file = here("data/processed_data.RData"))
