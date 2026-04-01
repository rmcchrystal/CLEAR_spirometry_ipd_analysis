
# Objective 3
# Assess patient adherence to weekly home spirometry throughout the study period
# Same spirometry data between proxy datasets, so only need to use one

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_denom <- readRDS("created_metadata/denom.rds")
a_imp_df <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_3d) %>% 
  unnest(pex_3d)

# Calculate participation time for each patient
# Get number of weeks of follow-up (guide for checking weekly adherence)
b_fu <- a_imp_denom %>% 
  select(id, v1:v5, status) %>% 
  pivot_longer(v1:v5) %>% 
  filter(!is.na(value)) %>% 
  group_by(id) %>% 
  mutate(
    v1 = value[name == "v1"],
    last = max(value),
    dur = as.numeric(last - v1),
    weeks = dur / 7
  ) %>% 
  ungroup()

# Get set of weeks from start to end of participation
b_seq_wks <- b_fu %>% 
  select(id, name, date = value) %>% 
  group_by(id) %>% 
  complete(
    date = seq(min(date), max(date), by = "day")
  ) %>% 
  mutate(
    week_start = floor_date(date, unit = "week", week_start = 1),
    week_idx = dense_rank(week_start)
  ) %>% 
  filter(week_idx <= 51) %>% 
  ungroup()

# Check weekly adherence -------------------------------------------------------

# Median number of total sessions
b_total_sum <- a_imp_df %>% 
  filter(data == "home") %>% 
  select(id, date, fev1) %>%
  count(id) %>% 
  reframe(
    med_sessions = median(n),
    min_sessions = min(n),
    max_sessions = max(n)
  )

write_csv(b_total_sum, "outputs/objs/03_adherence/median_total_adherence.csv")

# Prepare data
# Flag adherent weeks
b_weekly <- a_imp_df %>% 
  filter(data == "home") %>% 
  select(id, date, fev1) %>%
  full_join(b_seq_wks %>% select(id, date, week = week_idx)) %>% 
  arrange(id, date) %>% 
  group_by(id, week) %>% 
  mutate(
    adh = if_else(sum(is.na(fev1)) == length(week), 0L, 1L)
  ) %>% 
  ungroup() %>% 
  arrange(id, date)

# Get distinct rows
# Add visit timepoints
b_weekly_df <- b_weekly %>% 
  select(id, week, adh) %>%
  distinct() %>% 
  left_join(a_imp_denom %>% select(id, status)) %>% 
  group_split(status) %>% 
  purrr::map_dfr(~ {
    
    if (unique(.x$status) %in% c("completed", "completed (trt dropout)")) {
      .x %>%
        group_by(id) %>%
        complete(week = 1:52) %>%
        mutate(adh = replace_na(adh, 0L)) %>%
        ungroup()
      
    } else {
      .x
    }
    
  }) %>% 
  mutate(
    adh = replace_na(adh, 0L),
    visit = case_when(
      week == 1L ~ "V1",
      week == 2L ~ "V2",
      week == 8L ~ "V3",
      week == 26L ~ "V4",
      week == 52L ~ "V5"
    )
  ) %>% 
  ungroup()

# Summarise by week
b_weekly_sum <- b_weekly_df %>% 
  group_by(id) %>% 
  reframe(
    n_weeks = length(unique(week)),
    adh = sum(adh)
  ) %>% 
  mutate(adh_pcnt = round(adh / n_weeks * 100, 2)) %>% 
  ungroup() %>% 
  left_join(a_imp_denom %>% select(id, status))

b_weekly_sum_all <- b_weekly_sum %>% 
  reframe(
    med_adh = median(adh_pcnt),
    min_adh = min(adh_pcnt),
    max_adh = max(adh_pcnt)
  )

# Counts
sum(b_weekly_sum$adh_pcnt == 100) / 288
sum(b_weekly_sum$adh_pcnt < 100 & b_weekly_sum$adh_pcnt >= 79) / 288
sum(b_weekly_sum$adh_pcnt < 79 & b_weekly_sum$adh_pcnt >= 50) / 288
sum(b_weekly_sum$adh_pcnt < 50 & b_weekly_sum$adh_pcnt >= 21) / 288
sum(b_weekly_sum$adh_pcnt < 21) / 288

b_weekly_brackets <- tibble(
  adh_pcnt = c("100", "<100 - >=79", "<79 - >=50", "<50 - >=25", "<25"),
  count = c(7, 57, 92, 64, 68),
  pcnt = c(2.43, 19.79, 31.94, 22.22, 23.61)
)

# Save
write_csv(b_weekly_sum, "outputs/objs/03_adherence/weekly_adh_data.csv")
write_csv(b_weekly_sum_all, "outputs/objs/03_adherence/weekly_adh_sum.csv")
write_csv(b_weekly_brackets, "outputs/objs/03_adherence/weekly_adh_groups.csv")

# Summarise by visits with conditions
# V1 - strictly the week itself
# V2-V5 - the week before, during and after it (i.e. +/- 7 days)
b_adh_visits <- b_weekly_df %>% 
  arrange(id, week) %>% 
  group_by(id) %>% 
  mutate(
    pre = case_when(
      visit == "V1" ~ week,
      !is.na(visit) ~ week - 1L,
      TRUE ~ NA_integer_
    ),
    post = case_when(
      visit == "V1" ~ week,
      !is.na(visit) ~ week + 1L,
      TRUE ~ NA_integer_
    ),
    pre_adh  = lag(adh,  1),
    post_adh = lead(adh, 1)
  ) %>% 
  ungroup() %>% 
  filter(!is.na(visit)) %>% 
  mutate(
    visit_adh = case_when(
      visit == "V1" ~ adh,
      visit %in% c("V2", "V3", "V4") ~ pre_adh + adh + post_adh,
      visit == "V5" ~ pre_adh + adh
    ),
    rdy_adh = if_else(visit_adh != 0L, 1L, 0L)
  )

b_visit_sum <- b_adh_visits %>% 
  group_by(visit) %>% 
  reframe(
    n = length(unique(id)),
    adh = sum(rdy_adh, na.rm = TRUE)
  ) %>% 
  cbind(tibble(not_dropout = c(288, 258, 251, 241, 231))) %>% 
  mutate(adh_pcnt = round(adh / not_dropout * 100, 2)) %>% 
  select(visit, n, not_dropout, adh, adh_pcnt)

write_csv(b_visit_sum, "outputs/objs/03_adherence/visit_adh.csv")

# Plot of visit adherence
b_visit_plot <- b_visit_sum %>% 
  mutate(
    visit = case_match(
      visit,
      "V1" ~ "V1 (Baseline)",
      "V2" ~ "V2 (Week 2)",
      "V3" ~ "V3 (Week 8)",
      "V4" ~ "V4 (Week 26)",
      "V5" ~ "V5 (Week 52)"
    )
  ) %>% 
  ggplot(aes(x = visit, y = adh_pcnt, group = 1)) +
  geom_point(colour = "steelblue") +
  geom_line(colour = "steelblue") + 
  theme_classic() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Visit", y = "Adherence (% of enrolled weeks)")

b_visit_plot

ggsave(
  "outputs/objs/03_adherence/visit_adh_plot.png",
  b_visit_plot,
  width = 6,
  height = 4,
  units = "in"
)

# Number of dropouts at each visit
a_imp_denom %>% 
  filter(status == "dropout") %>% 
  select(id, status, date = dropout_date, v1:v5) %>% 
  pivot_longer(v1:v5) %>% 
  group_by(id) %>% 
  filter(value <= date) %>% 
  slice_max(value, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(id, last_vis = name) %>% 
  count(last_vis)
