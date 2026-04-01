
# Add proxy lung function where missing at PEx start and end

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_df <- readRDS("processed_data/all_joined_data.rds")

# Get PEx boundaries
b_pex_bounds <- a_imp_df %>%
  filter(!is.na(pex_id)) %>%
  group_by(site, type, id, pex_id) %>%
  reframe(
    start = min(date[tmpt == "pex_start"]),
    end = max(date[tmpt == "pex_end"])
  )

# Create spirometry lookup
b_spiro_lkp <- a_imp_df %>%
  select(
    site, id, date, data, 
    fev1, fev1p, fev1pp, 
    fvc, fvcp, fvcpp,
    fev1_fvc, fef2575, pef
  ) %>%
  filter(!is.na(fev1))

# Helper function to find proxy lung function
find_proxy <- function(target_date, df, lower_gap, upper_gap) {
  df %>%
    mutate(gap = as.numeric(date - target_date)) %>%
    filter(
      gap >= lower_gap,
      gap <= upper_gap
    ) %>%
    arrange(abs(gap)) %>%
    slice(1)
}

# Assign proxies at PEx episode-level for both rules
# +/- 3 days of start or end
# 7 days before start, 7 days after end
c_proxy <- b_pex_bounds %>% 
  rowwise() %>% 
  mutate(
    
    # Subset spirometry for a participant
    spiro_df = list(
      b_spiro_lkp %>% filter(site == .env$site, id == .env$id)
    ),
    
    # +/- 3 days
    start3 = list({
      df <- b_spiro_lkp %>% filter(site == .env$site, id == .env$id)
      obs <- df %>% filter(date == start)
      
      if (nrow(obs) > 0) {
        obs %>% mutate(flag = "observed", gap = 0)
      } else {
        find_proxy(start, df, -3, 3) %>%
          mutate(flag = "proxied")
      }
    }),
    
    end3 = list({
      df <- b_spiro_lkp %>% filter(site == .env$site, id == .env$id)
      obs <- df %>% filter(date == end)
      
      if (nrow(obs) > 0) {
        obs %>% mutate(flag = "observed", gap = 0)
      } else {
        find_proxy(end, df, -3, 3) %>%
          mutate(flag = "proxied")
      }
    }),
    
    # Assymetric 7 day rule (before = start, after = end)
    start7 = list({
      df <- b_spiro_lkp %>% filter(site == .env$site, id == .env$id)
      obs <- df %>% filter(date == start)
      
      if (nrow(obs) > 0) {
        obs %>% mutate(flag = "observed", gap = 0)
      } else {
        find_proxy(start, df, -7, 0) %>%
          mutate(flag = "proxied")
      }
    }),
    
    end7 = list({
      df <- b_spiro_lkp %>% filter(site == .env$site, id == .env$id)
      obs <- df %>% filter(date == end)
      
      if (nrow(obs) > 0) {
        obs %>% mutate(flag = "observed", gap = 0)
      } else {
        find_proxy(end, df, 0, 7) %>%
          mutate(flag = "proxied")
      }
    })
    
  ) %>% 
  ungroup()

# Summarise proxies
d_proxy_sum <- c_proxy %>%
  mutate(
    # ±3 rule
    start3_status = case_when(
      map_int(start3, nrow) == 0 ~ "missing",
      start3[[1]]$flag[1] == "observed" ~ "observed",
      TRUE ~ "proxied"
    ),
    end3_status = case_when(
      map_int(end3, nrow) == 0 ~ "missing",
      end3[[1]]$flag[1] == "observed" ~ "observed",
      TRUE ~ "proxied"
    ),
    
    # 7-day rule
    start7_status = case_when(
      map_int(start7, nrow) == 0 ~ "missing",
      start7[[1]]$flag[1] == "observed" ~ "observed",
      TRUE ~ "proxied"
    ),
    end7_status = case_when(
      map_int(end7, nrow) == 0 ~ "missing",
      end7[[1]]$flag[1] == "observed" ~ "observed",
      TRUE ~ "proxied"
    )
  )

# +/-3 day summary
d_proxy_sum %>%
  mutate(type = if_else(!type == "missed", "adj", type)) %>% 
  group_by(type) %>% 
  summarise(
    total_pex = n(),
    total_participants = n_distinct(id),
    
    start_observed_participants = n_distinct(id[start3_status == "observed"]),
    start_observed = sum(start3_status == "observed"),
    
    start_proxied_participants = n_distinct(id[start3_status == "proxied"]),
    start_proxied  = sum(start3_status == "proxied"),
    
    start_missing_participants = n_distinct(id[start3_status == "missing"]),
    start_missing  = sum(start3_status == "missing"),
    
    end_observed_participants = n_distinct(id[end3_status == "observed"]),
    end_observed = sum(end3_status == "observed"),
    
    end_proxied_participants = n_distinct(id[end3_status == "proxied"]),
    end_proxied  = sum(end3_status == "proxied"),
    
    end_missing_participants = n_distinct(id[end3_status == "missing"]),
    end_missing  = sum(end3_status == "missing"),
    
    both_participants = n_distinct(id[start3_status == "proxied" & end3_status == "proxied"]),
    both = sum(start3_status == "proxied" & end3_status == "proxied")
    
  ) %>% 
  pivot_longer(-type) %>% 
  print(n = Inf)

# Assymetric 7d summary
d_proxy_sum %>%
  mutate(type = if_else(!type == "missed", "adj", type)) %>% 
  group_by(type) %>% 
  summarise(
    total_pex = n(),
    total_participants = n_distinct(id),
    
    start_observed_participants = n_distinct(id[start7_status == "observed"]),
    start_observed = sum(start7_status == "observed"),
    
    start_proxied_participants = n_distinct(id[start7_status == "proxied"]),
    start_proxied  = sum(start7_status == "proxied"),
    
    start_missing_participants = n_distinct(id[start7_status == "missing"]),
    start_missing  = sum(start7_status == "missing"),
    
    end_observed_participants = n_distinct(id[end7_status == "observed"]),
    end_observed = sum(end7_status == "observed"),
    
    end_proxied_participants = n_distinct(id[end7_status == "proxied"]),
    end_proxied  = sum(end7_status == "proxied"),
    
    end_missing_participants = n_distinct(id[end7_status == "missing"]),
    end_missing  = sum(end7_status == "missing"),
    
    both_participants = n_distinct(id[start7_status == "proxied" & end7_status == "proxied"]),
    both = sum(start7_status == "proxied" & end7_status == "proxied")
  ) %>% 
  pivot_longer(-type) %>% 
  print(n = Inf)

# Prepare data with +/-3 proxies (or observed) ---------------------------------

# Keep adjudicated and missed PEx within 7d and +/- 3d
# Analysis using each set

# +/- 3d proxy
# Take start and end dates
e_pex_3d <- d_proxy_sum %>% 
  group_by(type) %>% 
  filter(
    start3_status %in% c("proxied", "observed") 
    & end3_status %in% c("proxied", "observed")
  ) %>% 
  ungroup() %>% 
  mutate(
    cmbn = map2(start3, end3, ~ {
      .x %>% 
        mutate(
          tmpt = "pex_start",
          flag = if_else(flag == "observed", flag, "proxied_3d")
        ) %>% 
        select(-gap) %>% 
        full_join(
          .y %>% 
            mutate(
              tmpt = "pex_end",
              flag = if_else(flag == "observed", flag, "proxied_3d")
            ) %>% 
            select(-gap)
        )
    })
  ) %>% 
  select(type, cmbn) %>% 
  unnest(cmbn) %>% 
  ungroup() %>% 
  arrange(site, id, date)

# Check for duplicates - none
e_pex_3d %>% group_by(site, id, date, data, fev1, fvc) %>% filter(n() > 1)

# Join with spirometry data
e_df_3d <- a_imp_df %>% 
  filter(!is.na(data)) %>% 
  select(-c(pex_id, type, tmpt, disease_status)) %>% 
  left_join(
    e_pex_3d %>% 
      select(
        site, id, date, data, fev1, fev1p, fev1pp, fvc, fvcp, fvcpp, fev1_fvc,
        fef2575, pef,
        type, flag, tmpt
      ) %>% 
      rename(proxy_flag = flag)
  ) %>% 
  group_by(site, id, date) %>% 
  fill(type, proxy_flag, tmpt, .direction = "downup") %>% 
  ungroup() %>% 
  arrange(site, id, date)

# Pull out start and end dates, calculate difference
# Save
e_pex_unq_3d <- e_df_3d %>% 
  distinct(site, id, date, tmpt) %>% 
  filter(!is.na(tmpt)) %>% 
  group_by(site, id, date) %>% 
  pivot_wider(names_from = "tmpt", values_from = "date") %>% 
  unnest() %>% 
  group_by(site, id) %>% 
  mutate(
    pex_diff = as.numeric(pex_end - pex_start),
    pex_id = row_number()
  ) %>% 
  ungroup()

# write_csv(e_pex_unq_3d, "Outputs/pex_summaries/3d_proxy_pex_dates.csv")

# Join PEx start and end dates and add flags identifying status
e_flag_3d <- e_df_3d %>% 
  left_join(
    e_pex_unq_3d %>% 
      rename(
        start = pex_start,
        end = pex_end
      ) %>% 
      select(-pex_diff)
  ) %>% 
  mutate(
    in_pex = date >= start & date <= end
  )

# Identify timepoint and disease status
e_id_3d <- e_flag_3d %>% 
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
      is_end ~ "pex_end",
      in_pex ~ "during_pex",
      TRUE ~ "stable"
    ),
    pex_status = if_else(in_pex, "pex", "stable")
  ) %>%
  select(-c(start, end, in_pex, is_start, is_end))

# Final tidying
e_tidy_3d <- e_id_3d %>% 
  group_by(site, id, pex_status) %>% 
  fill(type, .direction = "downup") %>% 
  ungroup() %>% 
  mutate(type = if_else(type != "missed", "full", type)) %>% 
  select(
    site, id, date, data, status, covid_tmpt, fev1, fev1p, fev1pp, 
    fvc, fvcp, fvcpp, fev1_fvc, fef2575, pef,
    tmpt, pex_status, pex_id, type, proxy_flag
  ) %>% 
  arrange(site, id, date)

# Final check that all 190 are present
e_tidy_3d %>% distinct(site, id, date, tmpt, type) %>% count(tmpt, type)

# Nest
e_nest_3d <- e_tidy_3d %>% group_by(site, id) %>% nest(.key = "pex_3d")

# Prepare data with 7d proxies (or observed) -----------------------------------

f_pex_7d <- d_proxy_sum %>% 
  group_by(type) %>% 
  filter(
    start7_status %in% c("proxied", "observed") 
    & end7_status %in% c("proxied", "observed")
  ) %>% 
  ungroup() %>% 
  mutate(
    cmbn = map2(start7, end7, ~ {
      .x %>% 
        mutate(
          tmpt = "pex_start",
          flag = if_else(flag == "observed", flag, "proxied_7d")
        ) %>% 
        select(-gap) %>% 
        full_join(
          .y %>% 
            mutate(
              tmpt = "pex_end",
              flag = if_else(flag == "observed", flag, "proxied_7d")
            ) %>% 
            select(-gap)
        )
    })
  ) %>% 
  select(type, cmbn) %>% 
  unnest(cmbn) %>% 
  ungroup() %>% 
  arrange(site, id, date)

# Check for duplicates - 10 rows, 5 PEx to be removed
f_pex_7d %>% 
  group_by(site, id, date, data, fev1, fvc) %>% 
  filter(n() > 1)

# Construct unique 7d intervals tibble
f_pex_unq_7d <- f_pex_7d %>% 
  distinct(site, id, date, tmpt) %>% 
  filter(!is.na(tmpt)) %>% 
  group_by(site, id, date) %>% 
  pivot_wider(names_from = "tmpt", values_from = "date") %>% 
  unnest() %>% 
  arrange(site, id, pex_start) %>% 
  group_by(site, id) %>% 
  mutate(
    collision = pex_start <= lag(pex_end)
  ) %>% 
  ungroup()

# Remove collisions - 195 PEx
f_pex_fix_7d <- f_pex_unq_7d %>%
  group_by(site, id) %>%
  mutate(
    collision = pex_start <= lag(pex_end)
  ) %>%
  ungroup() %>% 
  filter(
    !collision | is.na(collision),
    !is.na(pex_end)
  ) %>% 
  select(-collision)

# Rebuild tibble
f_pex_clean_7d <- f_pex_7d %>%
  left_join(
    f_pex_fix_7d %>%
      pivot_longer(
        cols = c(pex_start, pex_end),
        names_to = "tmpt",
        values_to = "date"
      ) %>%
      mutate(
        tmpt = recode(
          tmpt,
          pex_start = "pex_start",
          pex_end = "pex_end"
        ),
        flagged = "keep"
      ),
    by = c("site", "id", "date", "tmpt")
  ) %>% 
  rename(proxy_flag = flag) %>% 
  filter(!is.na(flagged)) %>% 
  select(-flagged)

# Join with spirometry data
f_df_7d <- a_imp_df %>% 
  filter(!is.na(data)) %>% 
  select(-c(pex_id, type, tmpt, disease_status)) %>% 
  arrange(site, id, date) %>% 
  left_join(f_pex_clean_7d) %>% 
  group_by(site, id, date) %>% 
  fill(type, proxy_flag, tmpt, .direction = "downup") %>% 
  ungroup() %>% 
  arrange(site, id, date)

# Pull out start and end dates, calculate difference
# Save
f_pex_unq_7d <- f_df_7d %>% 
  distinct(site, id, date, tmpt) %>% 
  filter(!is.na(tmpt)) %>% 
  group_by(site, id, date) %>% 
  pivot_wider(names_from = "tmpt", values_from = "date") %>% 
  unnest() %>% 
  group_by(site, id) %>% 
  mutate(
    pex_diff = as.numeric(pex_end - pex_start),
    pex_id = row_number()
  ) %>% 
  ungroup()

# write_csv(f_pex_unq_7d, "Outputs/pex_summaries/7d_proxy_pex_dates.csv")

# Join PEx start and end dates and add flags identifying status
date_map_7d <- f_df_7d %>%
  distinct(site, id, date) %>%
  left_join(f_pex_unq_7d) %>%
  filter(date >= pex_start & date <= pex_end) %>%
  select(site, id, date, pex_id, pex_start, pex_end)

f_flag_7d <- f_df_7d %>%
  left_join(
    date_map_7d %>%
      rename(
        start = pex_start,
        end   = pex_end
      ),
    by = c("site", "id", "date")
  ) %>%
  mutate(
    in_pex = !is.na(pex_id)
  )

# Identify timepoint and disease status
f_id_7d <- f_flag_7d %>% 
  mutate(
    in_pex = !is.na(pex_id),
    is_start = in_pex & date == start,
    is_end   = in_pex & date == end,
    tmpt = case_when(
      is_start ~ "pex_start",
      is_end   ~ "pex_end",
      in_pex   ~ "during_pex",
      TRUE     ~ "stable"
    ),
    pex_status = if_else(in_pex, "pex", "stable")
  ) %>%
  select(-c(start, end, in_pex, is_start, is_end))

# Final tidying
f_tidy_7d <- f_id_7d %>% 
  group_by(site, id, pex_status) %>% 
  fill(type, .direction = "downup") %>% 
  ungroup() %>% 
  mutate(type = if_else(type != "missed", "full", type)) %>% 
  select(
    site, id, date, data, status, covid_tmpt, fev1, fev1p, fev1pp, 
    fvc, fvcp, fvcpp, fev1_fvc,  fef2575, pef,
    tmpt, pex_status, pex_id, type, proxy_flag
  ) %>% 
  arrange(site, id, date)

# Final check that all 195 are present
f_tidy_7d %>% distinct(site, id, date, tmpt, type) %>% count(tmpt, type)

# Nest
f_nest_7d <- f_tidy_7d %>% group_by(site, id) %>% nest(.key = "pex_7d")

# Join and save ----------------------------------------------------------------

# Add original version (used in obj1 and obj2)
d_df <- e_nest_3d %>% 
  left_join(f_nest_7d) %>% 
  left_join(
    a_imp_df %>% 
      group_by(site, id) %>% 
      nest(.key = "pex_orig")
  ) %>% 
  select(site, id, pex_orig, pex_3d, pex_7d)

# Save
saveRDS(d_df, "Processed_data/analysis_data.rds")
