
# Objective 4
# Evaluate the quality of lung function measurements for home and clinic
# Same spirometry data between proxy datasets, so only need to use one

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_denom <- readRDS("created_metadata/denom.rds")
a_imp_df <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_3d) %>% 
  unnest(pex_3d)

a_imp_home_all <- readRDS("processed_data/obj4_df.rds")
a_imp_clinic_all <- readRDS("processed_data/pro_distinct.rds")

# Prepare home
# Add measurement counts per date
# Identify two highest measures
b_home <- a_imp_home_all %>% 
  filter(measure %in% c("fev1", "fvc")) %>% 
  group_by(id, measure, date) %>% 
  reframe(
    n = n(),
    hi  = nth(value, 1, order_by = desc(value)),
    hi2 = nth(value, 2, order_by = desc(value), default = NA_real_)
  ) %>% 
  mutate(
    diff = abs(hi - hi2),
    win150 = if_else(diff <= 0.15, 1L, 0L),
    win200 = if_else(diff <= 0.20, 1L, 0L),
    win250 = if_else(diff <= 0.25, 1L, 0L),
    mt250 = if_else(diff > 0.25, 1L, 0L)
  ) %>% 
  arrange(id, date)

# Grade according to ATS criteria
c_ats <- b_home %>% 
  mutate(
    grade = case_when(
      n >= 3 & win150 == 1L ~ "A",
      n == 2 & win150 == 1L ~ "B",
      n >= 2 & win200 == 1L ~ "C",
      (n >= 2 & win250 == 1L) ~ "D",
      (n >= 2 & mt250 == 1L) | n == 1 ~ "E",
    )
  )

# Summarise overall
c_ats_sum <- c_ats %>% 
  group_by(measure, grade) %>% 
  reframe(n = n()) %>% 
  group_by(measure) %>% 
  mutate(
    total = sum(n),
    pcnt = round((n / total) * 100, 2)
  )

write_csv(c_ats_sum, "outputs/objs/04_quality/sum_qc_home_all.csv")

# Summarise by patient
c_ats_sum_patients <- c_ats %>% 
  group_by(id, measure, grade) %>% 
  reframe(n = n()) %>% 
  group_by(id, measure) %>% 
  mutate(
    total = sum(n),
    pcnt = round((n / total) * 100, 2)
  ) %>% 
  ungroup() %>% 
  select(id, grade, total, n, pcnt) %>% 
  mutate(value = paste0(n, " (", pcnt, ")"))

write_csv(c_ats_sum_patients, "outputs/objs/04_quality/sum_qc_home_patients.csv")

# Prepare clinic
c_clinic <- a_imp_clinic_all %>% 
  filter(measure %in% c("fev1", "fvc")) %>% 
  group_by(id, date, measure) %>%
  slice_max(value, n = 1) %>%
  ungroup() %>%
  select(id, date, measure, value, contains("reproducible")) %>%
  distinct() %>% 
  mutate(
    qc = case_when(
      measure == "fev1" & fev1_reproducible ~ 1,
      measure == "fev1" & !fev1_reproducible ~ 0,
      measure == "fvc" & fvc_reproducible ~ 1,
      measure == "fvc" & !fvc_reproducible ~ 0
    )
  )

# Summarise clinic based on reproducibility
c_clinic_sum <- c_clinic %>% 
  group_by(measure) %>% 
  reframe(
    n = n(),
    qc = sum(qc),
    qc_pcnt = round((qc / n) * 100, 2)
  )

# Save
write_csv(c_clinic_sum, "outputs/objs/04_quality/sum_qc_clinic.csv")
