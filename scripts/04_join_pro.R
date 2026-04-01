
source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import
a_imp_patients <- read_rds("processed_data/patients.rds") %>% distinct()
a_imp_pro <- read_rds("processed_data/pro.rds")

# Many-to-many checks
a_imp_patients %>%
  count(across(all_of(c("site", "spiro_id", "trn")))) %>%
  filter(n > 1)

# Join patients and pro data
b_join <- a_imp_patients %>% select(-creation_date) %>% left_join(a_imp_pro)

# Final tidying
b_tidy <- b_join %>% 
  rename(id = subjid) %>% 
  arrange(site, id) %>% 
  pivot_longer(
    c(fev1, fvc, pef, fev1_fvc, fef2575, contains("zscore"), contains("predicted")),
    names_to = "measure",
    values_to = "value"
  ) %>% 
  pivot_longer(
    contains("unit"),
    names_to = "measure_unit",
    values_to = "unit"
  ) %>% 
  mutate(
    across(c(measure_unit, unit), ~ if_else(grepl("z|pred", measure), NA, .)),
    measure_unit = gsub("_unit", "", measure_unit),
    unit = gsub("\\[|\\]", "", unit)
  ) %>% 
  filter(measure == measure_unit  | is.na(measure_unit)) %>% 
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
    site:session_id, package_id, 
    creation_date, date_ext, date, date_fmt, calibration_ok:unit
  )

# Checks
b_fmt_dt %>% filter(is.na(date)) # No missing dates
b_fmt_dt %>% count(id, spiro_id) %>% count(id) # 278 patients
b_fmt_dt %>% reframe(across(everything(), ~ sum(is.na(.)))) # Some NA measures

# Save before taking unique measurements
# write_csv(b_fmt_dt, "processed_data/pro_full.csv")
# write_rds(b_fmt_dt, "processed_data/pro_full.rds")

# Check NA measures
b_fmt_dt %>% filter(is.na(value)) %>% count(id) # Across 138 patients, remove

# Filter out NA measures
b_remNA2 <- b_fmt_dt %>% filter(!is.na(value))
b_remNA2 %>% count(id, spiro_id) %>% count(id) # 278 patients (same)

# Get unique measurements per patient (full set across transfers)
b_full_set <- b_remNA2 %>% 
  arrange(
    id, site, session_id, package_id, measure, value, desc(date)
  ) %>%
  distinct(
    id, site, session_id, package_id, measure, value,
    .keep_all = TRUE
  )

# 278 in final set
b_full_set %>% count(id)

# Save
write_csv(b_full_set, "processed_data/pro_distinct.csv")
write_rds(b_full_set, "processed_data/pro_distinct.rds")
