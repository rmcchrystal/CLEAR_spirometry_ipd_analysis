
# Prepares clinic spirometry data
source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)
library(readxl)

# Get files and names
a_pro_files <- list.files("data/pro", full.names = TRUE, recursive = TRUE)
a_pro_files <- gsub("\\~\\$", "", a_pro_files)
a_pro_names <- list.files("data/pro", recursive = TRUE)
a_pro_names <- gsub("\\~\\$", "", a_pro_names)

# Import
a_pro_imp <- lapply(a_pro_files, read_excel)
names(a_pro_imp) <- a_pro_names

# Get variable lookup
a_imp_vars <- read_csv("created_metadata/defs_pro.csv") %>% 
  filter(relevant == 1) %>% 
  select(var)

# Add transfer as variable
b_add_trn <- lapply(names(a_pro_imp[sapply(a_pro_imp, nrow) > 0]), function(i) {
  nm <- i
  transfer <- gsub("/.*", "", nm)
  site <- str_extract(nm, regex_site)
  print(site)
  a_pro_imp[[i]] %>% 
    select(all_of(a_imp_vars$var)) %>% 
    rename_with(., tolower) %>% 
    mutate(
      site = site,
      trn = transfer
    )
})

b_bind %>% count(site) %>% print(n = 29)

# Bind rows
b_bind <- bind_rows(b_add_trn) %>% 
  mutate(
    site = case_when(
      grepl("SITE|^S", site) ~ tolower(site),
      TRUE ~ site
    ),
    site = case_when(
      grepl("site", site) ~ gsub("site", "s", site),
      TRUE ~ site
    ),
    site = case_match(
      site,
      "s 19" ~ "s19",
      "s003" ~ "s03",
      "s0221" ~ "s02",
      "s0320241017" ~ "s03",
      "s1920240417" ~ "s19",
      "s2120250129" ~ "s21",
      .default = site
    )
  ) %>% 
  rename(spiro_id = patient_id) %>% 
  select(trn, site, everything())

# Save
write_csv(b_bind, "processed_data/pro.csv")
write_rds(b_bind, "processed_data/pro.rds")
