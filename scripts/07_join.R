
# Join home, clinic and pex data
# Add proxy lung function where missing at start and end
# Calculate predicted FEV1/FVC, then % predicted
# 460 pex included
# 2 were removed without start, end or any treatment dates

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)
library(rspiro)

# Import datasets
a_imp_denom <- readRDS("created_metadata/denom_chars.rds")
a_imp_home <- readRDS("processed_data/home_distinct.rds")
a_imp_clinic <- readRDS("processed_data/pro_distinct.rds")

# Add a pex_id before any modifications
a_imp_pex <- readRDS("processed_data/pex_tidy.rds") %>% 
  select(-covid_tmpt) %>% 
  group_by(site, id) %>% 
  mutate(pex_id = row_number()) %>% 
  ungroup()

# Check and fix overlapping PEx ------------------------------------------------

# Has to be addressed first or else it creates difficulties with harmonising
# with spirometry data (overlapping PEx cycles) and adding proxy LF

# Identify overlap clusters
b_pex_clusters <- a_imp_pex %>%
  arrange(id, start, end) %>%
  group_by(id) %>%
  mutate(
    start_num = as.numeric(start),
    end_num   = as.numeric(end),
    end_running = cummax(end_num),
    new_episode = start_num > lag(end_running, default = -Inf),
    pex_cluster = cumsum(new_episode)
  ) %>%
  ungroup() %>%
  select(-start_num, -end_num, -end_running, -new_episode)

# Collapse overlapping clusters
b_collapse <- b_pex_clusters %>%
  group_by(id, pex_cluster) %>%
  filter(n() > 1) %>%
  reframe(
    site  = first(site),
    type  = first(type),
    start = min(start),
    end   = max(end),
    n_pex_merged = n(),
    pex_original = list(pex)
  )

# Remove original overlapping PEx
b_rm_overlap <- b_pex_clusters %>%
  anti_join(
    b_pex_clusters %>%
      group_by(id, pex_cluster) %>%
      filter(n() > 1) %>%
      ungroup() %>%
      select(id, pex),
    by = c("id", "pex")
  )

# Combine and re-index PEx
b_overlap_fix <- b_rm_overlap %>%
  select(site, id, type, start, end) %>%
  bind_rows(b_collapse) %>%
  arrange(id, start, end) %>%
  group_by(id) %>%
  mutate(pex_id = row_number()) %>%
  ungroup() %>%
  select(-c(pex_cluster:pex_original), -n_pex_merged)

# Compare before and after
# 460 -> 437 (reduced by 23 after collapsing)
a_imp_pex %>% reframe(n_pex = n())
b_overlap_fix %>% reframe(n_pex = n())

# Prepare and join spirometry data ---------------------------------------------

# Get relevant variables in home and clinic
# Take best measures on each date
# Fix handful of dates that were pre-CLEAR (device error)
c_rlv_home <- a_imp_home %>% 
  select(site:id, date, measure:unit) %>% 
  mutate(data = "home") %>% 
  select(-unit) %>% 
  mutate(
    date = case_when(
      id == "R07010" & grepl("2016", date) ~ update(date, year = 2019),
      id == "R21010" & grepl("2016", date) ~ update(date, year = 2023),
      TRUE ~ date
    )
  )

# Save for objective 4 of analysis (quality of home spirometry)
# saveRDS(c_rlv_home, "processed_data/obj4_df.rds")

c_home_best <- c_rlv_home %>% 
  group_by(id, date, measure) %>% 
  slice_max(value, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = measure,
    values_from = value
  )

c_rlv_clinic <- a_imp_clinic %>% 
  select(site:id, date, measure:unit) %>% 
  mutate(data = "clinic") %>% 
  select(-unit)

c_clinic_best <- c_rlv_clinic %>% 
  group_by(id, date, measure) %>% 
  slice_max(value, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = measure,
    values_from = value
  )

# Check there are no duplicates (none)
c_home_best %>% count(id, date) %>% filter(n > 1)
c_clinic_best %>% count(id, date) %>% filter(n > 1)

# Join home and clinic
# Remove any NA fev1 and fvc
c_join_spiro <- c_home_best %>% 
  full_join(c_clinic_best) %>% 
  arrange(id, date) %>% 
  select(site:pef, fev1p = fev1_predicted, fvcp = fvc_predicted)

# 56 missing fev1_fvc (non-issue, not used)
# c_join_spiro %>% reframe(across(everything(), ~ sum(is.na(.))))

# Join to denom
c_join_denom <- a_imp_denom %>% 
  select(site, id, status) %>% 
  left_join(c_join_spiro)

# Prepare and join PEx data ----------------------------------------------------

# Create PEx events explicit
d_pex_events <- b_overlap_fix %>%
  pivot_longer(
    cols = c(start, end),
    names_to = "tmpt",
    values_to = "date"
  ) %>%
  mutate(
    tmpt = case_match(
      tmpt,
      "start" ~ "pex_start",
      "end"   ~ "pex_end"
    )
  )

# 437 start and end dates as expected
d_pex_events %>% ungroup() %>% count(tmpt)

# Join PEx start and end dates to spirometry data
d_pex_join <- c_join_denom %>% 
  full_join(d_pex_events, by = c("site", "id", "date")) %>% 
  arrange(site, id, date)

# Classify timepoints and disease status
d_tmp <- d_pex_join %>% 
  select(-pex_id) %>%
  left_join(
    b_overlap_fix %>%
      select(site, id, pex_id, start, end),
    by = c("site", "id")
  ) %>% 
  mutate(in_pex = date >= start & date <= end)

d_pex_cls <- d_tmp %>%
  group_by(site, id, date, data) %>%
  reframe(
    pex_id = pex_id[in_pex][1],
    start = start[in_pex][1],
    end = end[in_pex][1],
    across(everything(), first)
  ) %>% 
  mutate(
    in_pex = !is.na(pex_id),
    is_start = in_pex & date == start,
    is_end   = in_pex & date == end,
    tmpt = case_when(
      is_start ~ "pex_start",
      is_end   ~ "pex_end",
      in_pex   ~ "during_pex",
      TRUE ~ "stable"
    ),
    pex_status = if_else(in_pex, "pex", "stable")
  ) %>%
  select(-start, -end)

# Check PEx counts - all 437 conserved
d_pex_cls %>%
  distinct(id, pex_id) %>%
  filter(!is.na(pex_id)) %>%
  nrow()

d_pex_cls %>% count(tmpt)

# Tidy dataset
d_pex_tidy <- d_pex_cls %>% 
  select(
    site:id, status, date:data, fef2575:fvcp, 
    pex_id, type:tmpt, disease_status = pex_status
  ) %>% 
  group_by(site, id) %>% 
  fill(status, .direction = "downup") %>% 
  group_by(site, id, pex_id) %>% 
  fill(type, .direction = "downup") %>% 
  ungroup() %>% 
  arrange(site, id, date)

# Add COVID timepoints ---------------------------------------------------------

# Start and end dates
zCOVID_pre <- as.Date("2020-03-12")
zCOVID_post <- as.Date("2021-11-01")

f_covid <- d_pex_tidy %>% 
  mutate(
    covid_tmpt = case_when(
      date < zCOVID_pre ~ "before",
      date >= zCOVID_pre & date < zCOVID_post ~ "during",
      date >= zCOVID_post ~ "after"
    )
  ) %>% 
  select(site:status, covid_tmpt, everything()) %>% 
  group_by(id) %>% 
  fill(status, .direction = "downup") %>% 
  ungroup()

# Final checks -----------------------------------------------------------------

# Check missingness
# 4 missing covid timepoints and dates had no lung function
f_covid %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything())

# All PEx conserved (437)
f_covid %>% 
  distinct(id, pex_id) %>% 
  distinct() %>% 
  filter(!is.na(pex_id)) %>% 
  nrow()

# Calculate % predicted for home spirometry ------------------------------------

# Join baseline characteristics
g_chars <- a_imp_denom %>% 
  select(site, id, age, sex, race, height, weight) %>% 
  left_join(f_covid)

# Calculate FEV1 and FVC % predicted for home spirometry
# Using rspiro package, pred_GLI
# Requires sex as factor, male as reference
# Height in metres
# Ethnicity as ordinal
g_pred <- g_chars %>% 
  mutate(
    sex = case_when(
      sex == 0L ~ 2,
      sex == 1L ~ 1,
      TRUE ~ sex
    ),
    race = case_match(
      race,
      "white" ~ 1,
      "black_aa" ~ 2,
      "asian" ~ 4,
      "other" ~ 5
    ),
    height = height/100,
    fev1p = case_when(
      data == "home" 
      & !is.na(fev1) & !is.na(age) & !is.na(sex) & !is.na(height)
      & !is.na(race) ~ rspiro::pred_GLI(age, height, sex, race, param = "FEV1"),
      TRUE ~ fev1p
    ),
    fvcp = case_when(
      data == "home" 
      & !is.na(fvc) & !is.na(age) & !is.na(sex) & !is.na(height)
      & !is.na(race) ~ rspiro::pred_GLI(age, height, sex, race, param = "FVC"),
      TRUE ~ fvcp
    ),
    fev1pp = (fev1/fev1p) * 100,
    fvcpp = (fvc/fvcp) * 100
  ) %>% 
  select(
    site, id, age, sex, race, height, weight, status, 
    covid_tmpt, date, data, fef2575, fev1, fev1p, fev1pp,
    fvc, fvcp, fvcpp, fev1_fvc, everything()
  ) %>% 
  arrange(site, id, date)

# Save
saveRDS(g_pred, "processed_data/all_joined_data.rds")
