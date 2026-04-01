
# Combine variation summary statistics

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")

# Import summaries
a_imp_cng <- read_csv("outputs/objs/02_variation/summaries/change_from_baseline_all.csv")
a_imp_var <- read_csv("outputs/objs/02_variation/summaries/icc_cv_diff_summary.csv")

# Visit summary
b_prep <- a_imp_cng %>% 
  filter(!is.na(visit), popn == "all") %>% 
  select(measure, visit, n, mean, lo, hi) %>% 
  rename(
    chng_mean = mean,
    chng_lo = lo,
    chng_hi = hi
  ) %>% 
  left_join(
    a_imp_var %>% select(-c(n, sd_diff))
  )

b_tidy <- b_prep %>% 
  mutate(
    across(-c(measure, visit, n), ~ round(., 3)),
    chng = paste0(chng_mean, " [", chng_lo, ", ", chng_hi, "]"),
    diff = paste0(mean_diff, " [", lo_diff, ", ", hi_diff, "]"),
    icc = paste0(icc, " [", icc_lo, ", ", icc_hi, "]"),
    cv = paste0(cv_within, " [", cv_lo, ", ", cv_hi, "]")
  ) %>% 
  select(
    measure, visit, n, chng, diff, icc, cv
  )

# Save
write_csv(b_tidy, "outputs/objs/02_variation/summaries/all.csv")
