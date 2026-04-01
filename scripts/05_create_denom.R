
# Creates an updated patient denominator
source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_denom <- read_csv("data/denom_old.csv")
a_imp_home <- readRDS("processed_data/home_distinct.rds")

# Get denom variables
b_denom_vars <- a_imp_denom %>%
  select(site:id, v1:status) %>% 
  mutate(across(v1:dropout_date, ~ as.Date(dmy(.))))

# Save
write_csv(b_denom_vars, "created_metadata/denom.csv")
write_rds(b_denom_vars, "created_metadata/denom.rds")
