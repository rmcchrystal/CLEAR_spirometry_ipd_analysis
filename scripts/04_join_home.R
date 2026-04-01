
source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import
a_imp_patients <- read_rds("processed_data/patients.rds") %>% distinct()
a_imp_home <- read_rds("processed_data/home.rds")

# Many-to-many checks
a_imp_patients %>%
  count(across(all_of(c("site", "spiro_id", "trn")))) %>%
  filter(n > 1)

# Join patients and home data
b_join <- a_imp_patients %>% select(-creation_date) %>% left_join(a_imp_home)

# Final tidying
b_tidy <- b_join %>% 
  rename(id = subjid) %>% 
  arrange(site, id) %>% 
  # select(-eot_criteria, -sot_criteria) %>%
  pivot_longer(
    c(fev1, fvc, pef, fev1_fvc, fef2575),
    names_to = "measure",
    values_to = "value"
  ) %>% 
  pivot_longer(
    fev1_unit:fef2575_unit,
    names_to = "measure_unit",
    values_to = "unit"
  ) %>% 
  mutate(
    measure_unit = gsub("_unit", "", measure_unit),
    unit = gsub("\\[|\\]", "", unit)
  ) %>% 
  filter(measure == measure_unit) %>% 
  select(-measure_unit)

# Remove NAs
b_tidy %>% reframe(across(everything(), ~ sum(is.na(.))))
b_remNA <- b_tidy %>% filter(!is.na(session_id))

# Format dates
b_fmt_dt <- b_remNA %>% 
  mutate(
    creation_date = gsub("/", ".", creation_date),
    date_ext = str_extract(creation_date, "^\\d{1,2}\\.\\d{1,2}\\.\\d{4}"),
    day_part = as.numeric(str_split_fixed(date_ext, "\\.", 3)[, 1]),
    month_part = as.numeric(str_split_fixed(date_ext, "\\.", 3)[, 2]),
    date_fmt = case_when(
      day_part > 12 ~ "dmy",
      month_part > 12 ~ "mdy",
      TRUE ~ "ambiguous"
    ),
    date = case_when(
      date_fmt == "dmy" ~ dmy(date_ext),
      date_fmt == "mdy" ~ mdy(date_ext),
      trn %in% c("aug22", "mar23", "sep21") ~ mdy(date_ext),
      TRUE ~ dmy(date_ext)
    )
  ) %>% 
  select(
    site:creation_date, date_ext, date, date_fmt, display_unit_id:unit,
    contains("criteria")
  )

# Checks
b_fmt_dt %>% filter(is.na(date)) # No missing dates
b_fmt_dt %>% count(id, spiro_id) %>% count(id) # 257 patients
b_fmt_dt %>% reframe(across(everything(), ~ sum(is.na(.)))) # Some NA measures

# Save before taking unique measurements
# write_csv(b_fmt_dt, "processed_data/home_full.csv")
# write_rds(b_fmt_dt, "processed_data/home_full.rds")

# Check NA measures
b_fmt_dt %>% filter(is.na(value)) %>% count(id) # Across 62 patients, remove

# Filter out NA measures
b_remNA2 <- b_fmt_dt %>% filter(!is.na(value))
b_remNA2 %>% count(id, spiro_id) %>% count(id) # 257 patients (same)

# Get unique measurements per patient (full set across transfers)
b_full_set <- b_remNA2 %>% 
  arrange(
    id, site, session_id, package_id, maneuver_id, measure, value, desc(date)
  ) %>%
  distinct(
    id, site, session_id, package_id, maneuver_id, measure, value,
    .keep_all = TRUE
  )

# 257 in final set
b_full_set %>% count(id)

# Get start and end criteria
b_full_set %>% 
  select(id, date, measure, value, contains("criteria")) %>% 
  write_csv("processed_data/home_criteria.csv")

# Import fixed dates and join with rest of home data
d_preCLEAR_fixed <- read_csv("outputs/preCLEAR/preCLEAR_df_fixed.csv") %>% 
  mutate(across(c(date, new_date), ~ as.Date(dmy(.)))) %>% 
  select(-v1)

d_dates_fixed <- b_full_set %>% 
  left_join(d_preCLEAR_fixed %>% select(id, date, new_date)) %>% 
  mutate(date = if_else(!is.na(new_date), new_date, date)) %>% 
  select(-new_date, -contains("criteria"))

# Save
write_csv(d_dates_fixed, "processed_data/home_distinct.csv")
write_rds(d_dates_fixed, "processed_data/home_distinct.rds")
