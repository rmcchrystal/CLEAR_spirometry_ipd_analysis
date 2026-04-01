
# Objective 1 | Prepare data

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_denom <- readRDS("created_metadata/denom.rds")

a_imp_df <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig)

b_df <- a_imp_df

# Participants who did and did not have any home spirometry
z_has_home <- a_imp_denom %>% 
  select(id, status, arm) %>% 
  left_join(
    a_imp_df %>% 
      filter(data == "home") %>% 
      mutate(flag = 1) %>% 
      distinct(id, flag)
  ) %>% 
  mutate(flag = if_else(is.na(flag), 0, flag))

write_csv(z_has_home, "scratch_data/patients_with_home_spiro.csv")

# Same day regardless of disease status
c_day_all <- b_df %>%
  group_by(id, date) %>% 
  filter(length(unique(data)) > 1) %>% 
  ungroup()

# Same day during stable disease
d_day_stable <- b_df %>%
  group_by(id, date) %>% 
  filter(length(unique(data)) > 1 & tmpt == "stable") %>% 
  ungroup()

# Same day between start and end of PEx
e_day_pex <- b_df %>%
  group_by(id, date) %>% 
  filter(length(unique(data)) > 1 & tmpt != "stable") %>% 
  ungroup()

# Same week regardless of disease status ---------------------------------------

# Find home that's within 7 days of clinic
f_week_clinic <- b_df %>% 
  filter(data == "clinic") %>% 
  rename(dt_clinic = date) %>% 
  rename_with(~ paste0(.x, "_clinic"), -c(id, dt_clinic))

f_week_home <- b_df %>% 
  filter(data == "home") %>% 
  rename(dt_home = date) %>% 
  rename_with(~ paste0(.x, "_home"), -c(id, dt_home))

# Join pairs
f_week_pairs <- f_week_clinic %>% 
  inner_join(f_week_home, by = "id") %>% 
  filter(between(dt_home, dt_clinic - days(7), dt_clinic + days(7))) %>% 
  mutate(delta = abs(as.numeric(dt_home - dt_clinic))) %>% 
  group_by(id, dt_clinic) %>%
  slice_min(delta, n = 1, with_ties = FALSE) %>%
  ungroup()

# Tidy up
f_week_all <- f_week_pairs %>% 
  mutate(pair_id = row_number()) %>% 
  pivot_longer(
    cols = -c(id, dt_clinic, dt_home, delta, pair_id),
    names_to = c(".value", "source"),
    names_sep = "_(?=[^_]+$)"
  ) %>% 
  mutate(
    data = source,
    date = if_else(source == "clinic", dt_clinic, dt_home)
  ) %>% 
  arrange(id, pair_id, desc(data == "clinic")) %>% 
  select(-pair_id) %>% 
  distinct()

# Same week during stable disease ----------------------------------------------

# Find home that's within 7 days of clinic during stability
g_week_clinic_stable <- b_df %>% 
  filter(data == "clinic" & tmpt == "stable") %>% 
  rename(dt_clinic = date) %>% 
  rename_with(~ paste0(.x, "_clinic"), -c(id, dt_clinic))

g_week_home_stable <- b_df %>% 
  filter(data == "home" & tmpt == "stable") %>% 
  rename(dt_home = date) %>% 
  rename_with(~ paste0(.x, "_home"), -c(id, dt_home))

g_week_pairs_stable <- g_week_clinic_stable %>% 
  inner_join(g_week_home_stable, by = "id") %>% 
  filter(between(dt_home, dt_clinic - days(7), dt_clinic + days(7))) %>% 
  mutate(delta = abs(as.numeric(dt_home - dt_clinic))) %>% 
  group_by(id, dt_clinic) %>%
  slice_min(delta, n = 1, with_ties = FALSE) %>%
  ungroup()

g_week_all_stable <- g_week_pairs_stable %>% 
  mutate(pair_id = row_number()) %>% 
  pivot_longer(
    cols = -c(id, dt_clinic, dt_home, delta, pair_id),
    names_to = c(".value", "source"),
    names_sep = "_(?=[^_]+$)"
  ) %>% 
  mutate(
    data = source,
    date = if_else(source == "clinic", dt_clinic, dt_home)
  ) %>% 
  arrange(id, pair_id, desc(data == "clinic")) %>% 
  select(-pair_id) %>% 
  distinct()

# Same week during pex ---------------------------------------------------------

# Find home that's within 7 days of clinic during PEx
h_week_clinic_pex <- b_df %>% 
  filter(data == "clinic" & tmpt != "stable") %>% 
  rename(dt_clinic = date) %>% 
  rename_with(~ paste0(.x, "_clinic"), -c(id, dt_clinic))

h_week_home_pex <- b_df %>% 
  filter(data == "home" & tmpt != "stable") %>% 
  rename(dt_home = date) %>% 
  rename_with(~ paste0(.x, "_home"), -c(id, dt_home))

h_week_pairs_pex <- h_week_clinic_pex %>% 
  inner_join(h_week_home_pex, by = "id") %>% 
  filter(between(dt_home, dt_clinic - days(7), dt_clinic + days(7))) %>% 
  mutate(delta = abs(as.numeric(dt_home - dt_clinic))) %>% 
  group_by(id, dt_clinic) %>%
  slice_min(delta, n = 1, with_ties = FALSE) %>%
  ungroup()

h_week_all_pex <- h_week_pairs_pex %>% 
  mutate(pair_id = row_number()) %>% 
  pivot_longer(
    cols = -c(id, dt_clinic, dt_home, delta, pair_id),
    names_to = c(".value", "source"),
    names_sep = "_(?=[^_]+$)"
  ) %>% 
  mutate(
    data = source,
    date = if_else(source == "clinic", dt_clinic, dt_home)
  ) %>% 
  arrange(id, pair_id, desc(data == "clinic")) %>% 
  select(-pair_id) %>% 
  distinct()

# Save data --------------------------------------------------------------------

# List of datasets
i_lst <- list(
  same_day_all = c_day_all,
  same_day_stable = d_day_stable,
  same_day_pex = e_day_pex,
  same_week_all = f_week_all,
  same_week_stable = g_week_all_stable,
  same_week_pex = h_week_all_pex
)

# Save
saveRDS(i_lst, "processed_data/obj1_data.rds")
