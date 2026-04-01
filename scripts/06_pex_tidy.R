
# Tidy exacerbation data
source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import pex data
a_imp_pex <- read_csv("data/pex/pex_everything.csv")

# Import patient denom
a_imp_denom <- readRDS("created_metadata/denom.rds")

# Totals to 462 pex
glimpse(a_imp_pex)
a_imp_pex %>% count(type)

# Tidy dataset -----------------------------------------------------------------

# Simplify pex types
# Set NA for 1899 dates
# Format dates
b_pex_type <- a_imp_pex %>% 
  mutate(
    type_simple = case_when(
      type == "No exacerbation" ~ "none",
      grepl("Not new", type) ~ "continued_pex",
      TRUE ~ "pex"
    ),
    across(start:treat_end, ~ if_else(grepl("1899", .), NA, .)),
    across(start:treat_end, ~ as.Date(dmy(.)))
  )

b_pex_type %>% count(type_simple)

# Remove continous pex and none
b_keep_pex <- b_pex_type %>% filter(type_simple == "pex")

# Check
# 203 patients, 462 pex
b_keep_pex %>% count(id) 
b_keep_pex %>% count(type_simple)
b_keep_pex %>% count(type)

# Simplify pex categories
b_simplify <- b_keep_pex %>% 
  mutate(
    type = case_match(
      type,
      "Confirmed missed exacerbation" ~ "missed",
      "Fully qualifying exacerbation" ~ "full",
      "Other exacerbation" ~ "other",
      "Partially qualifying exacerbation" ~ "partial"
    )
  )

# Flag:
# Missing start or end dates
# Start or end dates that are treatment dates
# Patients missing all treatment dates
b_flags <- b_simplify %>% 
  select(-type_simple) %>% 
  mutate(
    na_start = if_else(is.na(start), 1L, 0L),
    na_end = if_else(is.na(treat_end), 1L, 0L),
    na_treat = if_else(is.na(treat1) & is.na(treat2) & is.na(treat3), 1L, 0L),
    start_is_treat = case_when(
      !is.na(start) & start %in% c(treat1, treat2, treat3) ~ 1L,
      TRUE ~ 0L
    )
  )

b_xtra_dt <- b_flags %>% 
  anti_join(b_flags %>% filter(na_start == 1L | na_end == 1L)) %>% 
  full_join(
    read_csv("scratch_data/pex_missing_start_end.csv") %>% 
      select(-note) %>% 
      mutate(across(start:treat_end, ~ as.Date(dmy(.))))
  ) %>% 
  arrange(site, id, start) %>% 
  mutate(
    na_start = if_else(is.na(start), 1L, 0L),
    na_end = if_else(is.na(treat_end), 1L, 0L),
    na_treat = if_else(is.na(treat1) & is.na(treat2) & is.na(treat3), 1L, 0L),
    start_is_treat = case_when(
      !is.na(start) & start %in% c(treat1, treat2, treat3) ~ 1L,
      TRUE ~ 0L
    )
  )

# Summarise flags according to pex type
sum_flags <- b_xtra_dt %>% 
  group_by(type) %>% 
  reframe(
    n = length(unique(id)),
    n_pex = n(),
    has_start_date = sum(na_start == 0),
    na_start_date = sum(na_start),
    has_end_date = sum(na_end == 0),
    na_end_date = sum(na_end),
    treat_dates = sum(na_treat == 0),
    no_treat_dates = sum(na_treat),
    start_treat_same = sum(start_is_treat)
  )

# write_csv(sum_flags, "outputs/pex_summaries/pex_date_flags.csv")

# Address missing dates --------------------------------------------------------

# Add proxy start dates where missing
# Keep start as is if missing (NA in new variable, coalesced at end)
# If missing, take the nearest treatment date - 3 days
# If no treatment dates available, take the end date - 14 days
# If the proxy start date crosses V1, use the nearest treatment date as is
# Otherwise leave it as NA
c_proxy_start <- a_imp_denom %>% 
  select(site, id, v1) %>% 
  inner_join(b_xtra_dt) %>% 
  mutate(
    start_proxy = case_when(
      !is.na(start) ~ NA,
      is.na(start) & !is.na(treat1) ~ treat1 - 3,
      is.na(start) & is.na(treat1) & !is.na(treat2) ~ treat2 - 3,
      is.na(start) & is.na(treat1) & is.na(treat2) & !is.na(treat3) ~ treat3 - 3,
      !is.na(start) & is.na(treat1) & is.na(treat2) & is.na(treat3) ~ treat_end - 14,
      TRUE ~ NA
    ),
    start_proxy = case_when(
      as.numeric(start_proxy - v1) <= 0 & !is.na(treat1) ~ treat1,
      as.numeric(start_proxy - v1) <= 0 & !is.na(treat3) ~ treat3,
      TRUE ~ start_proxy
    ),
    diff = as.numeric(start_proxy - v1),
    start_new = coalesce(start, start_proxy)
  ) %>% 
  select(-v1, -diff)

# Add proxy end dates where missing
# Keep end as is if missing (NA in new variable, coalesced at end)
# If no treatment dates available, take the start date + 14 days
# If treatment dates all missing, add 14 days to start (EMBARC guidelines)
# Otherwise, take latest treatment date as end date
c_proxy_end <- c_proxy_start %>% 
  mutate(
    end_proxy = case_when(
      !is.na(treat_end) ~ NA,
      !is.na(start_new) & is.na(treat_end) ~ start_new + 14,
      TRUE ~ NA
    ),
    end_new = coalesce(treat_end, end_proxy)
  )

# Tidy up dataset with new dates
c_tidy <- c_proxy_end %>% 
  select(site:id, type, start_new, end_new) %>% 
  rename(start = start_new, end = end_new)

# Two with no start, end or treatment dates, removed (2 missed, 2 other)
c_tidy %>%
  filter(is.na(start) & is.na(start)) %>%
  write_csv("outputs/pex_summaries/n2_pex_no_dates.csv")

# Remove the 2 PEx with no dates
c_rm <- c_tidy %>% filter(!is.na(start) & !is.na(end))
c_rm %>% count(id)

# Count PEx type
c_rm %>% count(type)

# Fix four start and end dates
# This was flagged at some point in the past during data checks
# Calculate duration
# Where a PEx that met EMBARC (either full or partial) had a duration of
# less than 2 days (symptoms had to have persisted for at least 48 hours),
# modify the end date to be start + 14 days (typical expected duration of 
# antibiotic therapy)
c_fix <- c_rm %>% 
  mutate(
    start = case_when(
      id == "R10007" & start == "2020-02-03" ~ as.Date("2020-02-06"),
      id == "R02011" & start == "2024-05-13" ~ as.Date("2024-05-16"),
      TRUE ~ start
    ),
    end = case_when(
      id == "R10007" & end == "2019-04-20" ~ as.Date("2020-02-18"),
      id == "R02011" & end == "2024-04-25" ~ as.Date("2024-05-25"),
      TRUE ~ end
    ),
    pex_dur = as.numeric(end - start),
    end = if_else(
      pex_dur <= 2 & type %in% c("full", "partial"), 
      start + 14, 
      end
    ),
    pex_dur = as.numeric(end - start)
  )

# Shortest is 1 day, longest is 94 days, median 16 days
summary(c_fix$pex_dur)

# Index each PEx
c_index <- c_fix %>% 
  group_by(id) %>% 
  mutate(pex = row_number())

# Summarise by PEx type and COVID timepoints -----------------------------------

# COVID timepoints
zCOVID_pre <- as.Date("2020-03-12")
zCOVID_post <- as.Date("2021-11-01")

# Categorise PEx in relation to COVID timepoints
# Identify timepoint start and end were at
# Then combine (if timepoints differ, consider it during COVID)
e_pex_covid <- c_index %>% 
  mutate(
    start_covid = case_when(
      start < zCOVID_pre ~ "before",
      start >= zCOVID_pre & start < zCOVID_post ~ "during",
      start >= zCOVID_post ~ "after"
    ),
    end_covid = case_when(
      end < zCOVID_pre ~ "before",
      end >= zCOVID_pre & end < zCOVID_post ~ "during",
      end >= zCOVID_post ~ "after"
    ),
    covid_tmpt = case_when(
      start_covid == "before" & end_covid == "before" ~ "before",
      start_covid == "during" & end_covid == "during" ~ "during",
      start_covid == "after" & end_covid == "after" ~ "after",
      TRUE ~ "during"
    )
  ) %>% 
  select(-c(start_covid, end_covid))

e_pex_covid %>% count(type, site)
e_pex_covid %>% count(type)
e_pex_covid %>% count(type, covid_tmpt)

# Save
saveRDS(e_pex_covid, "processed_data/pex_tidy.rds")
