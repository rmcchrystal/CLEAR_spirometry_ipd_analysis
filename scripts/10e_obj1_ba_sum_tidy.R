
# Combine summary statistics for Bland-Altman analysis
# Empirical and regression outputs with correlations and p-values

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import summaries
a_ba_emp <- read_csv("outputs/objs/01_agreement/summaries/baltman_summaries.csv")
a_ba_cor <- read_csv("outputs/objs/01_agreement/summaries/baltman_correlation.csv")
a_ba_reg <- read_csv("outputs/objs/01_agreement/summaries/baltman_regression.csv")

# Tidy empirical summaries
b_emp <- a_ba_emp %>% 
  mutate(
    type_tidy = case_match(
      type,
      "same_week_all" ~ "± 7 days",
      "same_day_stable" ~ "Same day during stability",
      "same_day_pex" ~ "Same day during PEx",
      "same_week_pex" ~ "± 7 days during PEx"
    )
  ) %>% 
  filter(!is.na(type_tidy)) %>% 
  select(
    type = type_tidy, name, bias, hi, lo
  )

b_cor <- a_ba_cor %>% 
  mutate(
    type_tidy = case_when(
      type == "Irrespective of disease status" & tmpt == "same_week" ~ "± 7 days",
      type == "During stable periods" & tmpt == "same_day" ~ "Same day during stability",
      type == "Between PEx start and end" & tmpt == "same_day" ~ "Same day during PEx",
      type == "Between PEx start and end" & tmpt == "same_week" ~ "± 7 days during PEx"
    )
  ) %>% 
  filter(!is.na(type_tidy)) %>% 
  select(
    type = type_tidy, name, n, n_pairs, cor_pearsons, cor_p
  )

b_reg <- a_ba_reg %>% 
  mutate(
    reg_bias = intercept + slope * mu_max,
    reg_hi = reg_bias + 1.96 * sigma,
    reg_lo = reg_bias - 1.96 * sigma
  )

# Combine
b_cmbn <- b_emp %>% 
  left_join(b_cor) %>% 
  left_join(
    b_reg %>% 
      mutate(type = str_replace(type, "exacerbation", "PEx"))
  )

# Prepare for table
b_rdy <- b_cmbn %>%
  mutate(
    across(-c(type, name), ~ round(., 3)),
    across(c(cor_p, reg_p), ~ if_else(. < 0.01, "p < 0.01", paste0("p = ", .))),
    cor = paste0(cor_pearsons, " (", cor_p, ")"),
    mean_bias = paste0(bias, " [", hi, ", ", lo, "]"),
    reg_bias = paste0(reg_bias, " [", reg_hi, ", ", reg_lo, "]; ", reg_p)
  ) %>% 
  select(
    type, name, n, n_pairs, cor, mean_bias, reg_bias
  ) %>% 
  arrange(name)

# Save
write_csv(b_rdy, "outputs/objs/01_agreement/summaries/baltman_all.csv")
