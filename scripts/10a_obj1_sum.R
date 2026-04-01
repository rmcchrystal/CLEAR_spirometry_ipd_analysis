
# Objective 1
# Compare lung function measurements between home and clinic spirometry 
# performed on the same day and within the same week:
# a)	Irrespective of disease status
# b)	During clinically stable periods of bronchiectasis
# c)	Between the start and end of an exacerbation
# This script prepares plot data and generates summaries

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_lst <- readRDS("processed_data/obj1_data.rds")

# Add names as variables
b_lst <- lapply(names(a_imp_lst), function(i) {
  nm <- i
  a_imp_lst[[i]] %>% mutate(type = nm)
})

# Bind
c_df <- bind_rows(b_lst) %>% 
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    pex_cycle = cumsum(
      str_detect(type, "same_.*_pex") &
        !lag(str_detect(type, "same_.*_pex"), default = FALSE)
    ),
    pex_cycle = if_else(str_detect(type, "same_.*_pex"), pex_cycle, NA_integer_)
  ) %>%
  ungroup()

# PEx counts
pex_counts <- c_df %>%
  filter(str_detect(type, "same_.*_pex")) %>%
  distinct(type, tmpt = str_extract(type, "same_day|same_week"), id, pex_cycle) %>%
  filter(!is.na(pex_cycle)) %>%
  count(type, tmpt, name = "n_pex_total") %>% 
  mutate(
    tmpt = str_extract(type, "same_day|same_week"), 
    type = case_match(
      type, 
      "same_day_all" ~ "Irrespective of disease status", 
      "same_day_pex" ~ "Between PEx start and end", 
      "same_day_stable" ~ "During stable periods", 
      "same_week_all" ~ "Irrespective of disease status", 
      "same_week_pex" ~ "Between PEx start and end", 
      "same_week_stable" ~ "During stable periods" 
    ) 
  )

# Bland-Altman data
c_ba <- c_df %>% 
  select(
    id, data, date, fev1, fev1pp, fvc, fvcpp, fev1_fvc, fef2575, pef, type
  ) %>%
  mutate(tmpt = str_extract(type, "same_day|same_week")) %>% 
  pivot_longer(fev1:pef) %>% 
  distinct() %>% 
  group_by(id, date, type, tmpt, name) %>% 
  mutate(mu = sum(value, na.rm = TRUE) / 2) %>% 
  pivot_wider(names_from = data, values_from = value) %>% 
  group_by(id, type, tmpt, name) %>% 
  fill(clinic, .direction = "downup") %>% 
  ungroup() %>% 
  mutate(diff = clinic - home) %>% 
  group_by(type, tmpt, name) %>% 
  mutate(
    bias = mean(diff, na.rm = TRUE),
    stddev = sd(diff, na.rm = TRUE),
    hi = bias + 2 * stddev,
    lo = bias - 2 * stddev
  ) %>% 
  ungroup() %>% 
  filter(!is.na(home))

# Save
saveRDS(c_ba, "scratch_data/bland_altman_data.rds")

# Overall summary
c_sum_overall <- c_df %>% 
  select(
    id, data, date, fev1, fev1pp, fvc, fvcpp, fev1_fvc, fef2575, pef, type
  ) %>%
  pivot_longer(fev1:pef) %>% 
  distinct() %>% 
  filter(!is.na(value)) %>% 
  mutate(
    tmpt = str_extract(type, "same_day|same_week"), 
    type = case_match(
      type, 
      "same_day_all" ~ "Irrespective of disease status", 
      "same_day_pex" ~ "Between PEx start and end", 
      "same_day_stable" ~ "During stable periods", 
      "same_week_all" ~ "Irrespective of disease status", 
      "same_week_pex" ~ "Between PEx start and end", 
      "same_week_stable" ~ "During stable periods" 
    ) 
  ) %>% 
  group_by(data, name, type, tmpt) %>% 
  reframe(
    n = length(unique(id)), 
    mean = mean(value, na.rm = TRUE), 
    sd = sd(value, na.rm = TRUE), 
    se = (sd / sqrt(n)), 
    med = median(value, na.rm = TRUE), 
    min = min(value, na.rm = TRUE), 
    max = max(value, na.rm = TRUE), 
    lo = mean - 1.96 * se, 
    hi = mean + 1.96 * se 
  ) %>% 
  mutate(across(mean:hi, ~ round(., 2))) %>% 
  arrange(type, name, tmpt) %>%
  left_join(pex_counts, by = c("type", "tmpt")) %>% 
  select(data:n, n_pex = n_pex_total, everything())

# Save
write_csv(
  c_sum_overall,
  "outputs/objs/01_agreement/summaries/clinic_home_summaries.csv"
)

# Mean difference
c_sum_mudiff <- c_sum_overall %>% 
  select(data:se) %>% 
  pivot_wider(names_from = data, values_from = c(n, mean, sd, se)) %>% 
  group_by(name, type, tmpt) %>% 
  reframe(
    mean_diff = mean_clinic - mean_home,
    se_diff = sqrt(
      (sd_clinic^2 / n_clinic) + (sd_home^2 / n_home)
    ),
    lo_diff = mean_diff - 1.96 * se_diff,
    hi_diff = mean_diff + 1.96 * se_diff
  ) %>% 
  arrange(type, name, tmpt) %>% 
  mutate(across(mean_diff:hi_diff, ~ round(., 2)))

# Save
write_csv(
  c_sum_mudiff,
  "outputs/objs/01_agreement/summaries/mean_difference.csv"
)

# Summarise Bland-Altman
c_ba_sum <- c_ba %>%
  select(type, tmpt, name, bias, hi, lo) %>% 
  distinct() %>% 
  mutate(across(bias:lo, ~ round(., 4)))

# Save
write_csv(
  c_ba_sum,
  "outputs/objs/01_agreement/summaries/baltman_summaries.csv"
)

# Correlations -----------------------------------------------------------------

# Evidence of proportional bias for several timepoints/measurements
# Running for all

# Prepare data
d_cor_df <- c_ba %>% 
  mutate(
    type = case_match(
      type, 
      "same_day_all" ~ "Irrespective of disease status", 
      "same_day_pex" ~ "Between PEx start and end", 
      "same_day_stable" ~ "During stable periods", 
      "same_week_all" ~ "Irrespective of disease status", 
      "same_week_pex" ~ "Between PEx start and end", 
      "same_week_stable" ~ "During stable periods" 
    )
  )

d_cor <- d_cor_df %>% 
  group_by(type, tmpt, name) %>% 
  reframe(
    n = length(unique(id)),
    n_pairs = n(),
    cor_pearsons = cor(mu, diff, use = "complete.obs", method = "pearson"),
    cor_p = cor.test(mu, diff, method = "pearson")$p.value
  ) %>% 
  mutate(
    cor_pearsons = round(cor_pearsons, 2),
    cor_p = round(cor_p, 4)
  )

# Save
write_csv(d_cor, "outputs/objs/01_agreement/summaries/baltman_correlation.csv")

# Clinic vs home correlation
d_cor_spiro <- d_cor_df %>% 
  group_by(type, tmpt, name) %>% 
  reframe(
    n = length(unique(id)),
    n_pairs = n(),
    cor_pearsons = cor(clinic, home, use = "complete.obs", method = "pearson"),
    cor_p = cor.test(clinic, home, method = "pearson")$p.value
  ) %>% 
  mutate(
    cor_pearsons = round(cor_pearsons, 2),
    cor_p = round(cor_p, 4)
  )

write_csv(d_cor_spiro, "outputs/objs/01_agreement/summaries/spiro_correlation.csv")
