
# Prepares participant-level denominator
source("Scripts/00_config.R")
library(tidyverse)
library(readxl)

# Get files and names
a_pt_files <- list.files("data/patients", full.names = TRUE, recursive = TRUE)
a_pt_files <- gsub("\\~\\$", "", a_pt_files)

a_pt_names <- list.files("data/patients", recursive = TRUE)
a_pt_names <- gsub("\\~\\$", "", a_pt_names)

# Import
a_pt_imp <- lapply(a_pt_files, read_excel)
names(a_pt_imp) <- a_pt_names

# Add transfer as variable
b_add_trn <- lapply(names(a_pt_imp), function(i) {
  transfer <- i
  transfer <- gsub("/.*", "", transfer)
  a_pt_imp[[i]] %>% 
    mutate(trn = transfer) %>% 
    select(-Transfer) %>% 
    rename_with(., tolower)
})

# Bind rows
b_bind <- bind_rows(b_add_trn) %>% 
  select(-note, -...4) %>% 
  mutate(
    site = substr(custom_id, 1, 3),
    site = gsub("R", "s", site),
    site = if_else(site == "s35", "s01", site)
  ) %>% 
  arrange(site, custom_id) %>% 
  rename(
    spiro_id = id,
    subjid = custom_id
  ) %>% 
  select(site, subjid, spiro_id, everything())

# Save
write_csv(b_bind, "processed_data/patients.csv")
write_rds(b_bind, "processed_data/patients.rds")
