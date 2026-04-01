
# Objective 5
# Examine the responsiveness of home spirometry to changes in lung function 
# across an exacerbation (i.e. pre- and post-exacerbation)
# Using 3d PEx proxy data

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)

# Import data
a_imp_denom <- readRDS("created_metadata/denom.rds")

a_imp_df <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_3d) %>% 
  unnest(pex_3d)

a_imp_bl <- readRDS("processed_data/obj5_bl_3d.rds")

# Create PEx dataset -----------------------------------------------------------

# Nest by patients who had exacerbations
# Pull exacerbation dates and home data
# Full-join PEx dates to add any missing (i.e. occurred when only a clinic
# measure was available)
# R12001 - some data clearly formatted incorrectly (2020), participant was
# enrolled in 2022, fix
b_nest <- a_imp_df %>%
  mutate(
    date = as.character(date),
    date = case_when(
      id == "R12001" & grepl("2020", date) ~ gsub("2020", "2022", date),
      TRUE ~ date
    ),
    date = as.Date(ymd(date))
  ) %>% 
  group_by(site, id) %>% 
  filter(!all(tmpt == "stable")) %>% 
  nest() %>% 
  mutate(
    home_data = map(data, function(df) {
      df %>% filter(data == "home")
    }),
    pex_dates = map(data, function(df) {
      df %>% 
        filter(tmpt %in% c("pex_start", "pex_end")) %>% 
        select(status, date, tmpt) %>% 
        distinct()
    }),
    add_pex = map2(home_data, pex_dates, function(home, pex) {
      home %>% 
        full_join(pex) %>% 
        arrange(date)
    })
  )

# Flag exacerbation timepoints
# Extract PEx timelines
b_flag <- b_nest %>% 
  mutate(
    flag_tmpts = map(add_pex, ~ flag_exacerbation_timepoints(.x, type = "3d")),
    pex_timelines = map(flag_tmpts, ~ .x %>% filter(!is.na(flagged_tmpts)))
  ) %>% 
  select(id, flag_tmpts, pex_timelines)

# Fix flagged tmpt name for when PEx ended
c_complete_dat <- b_flag %>% 
  select(id, pex_timelines) %>% 
  unnest(pex_timelines) %>% 
  ungroup() %>% 
  mutate(
    flagged_tmpts = case_when(
      flagged_tmpts == "pex" & tmpt == "stable" ~ "stable_again",
      TRUE ~ flagged_tmpts
    )
  )

# 114 participants
# Tidy timepoint names
c_complete_compl <- c_complete_dat %>% 
  group_by(id, pex_cycle, flagged_tmpts) %>% 
  mutate(
    flagged_tmpts = case_match(
      flagged_tmpts,
      "last_stable_pre" ~ "Pre-PEx stable",
      "pex" ~ "PEx",
      "stable_again" ~ "Post-PEx"
    ),
    flagged_tmpts = case_when(
      tmpt == "pex_start" ~ "PEx start",
      tmpt == "pex_end" ~ "PEx end",
      TRUE ~ flagged_tmpts
    )
  ) %>% 
  ungroup()

# Recovery to pre-PEx FEV1 -----------------------------------------------------

# Pivot by lung function
lf_measures <- c(
  "fev1", "fev1pp", "fvc", "fvcpp", "fev1_fvc", "fef2575", "pef"
)

d_lf_long <- c_complete_compl %>%
  pivot_longer(
    all_of(lf_measures),
    names_to = "measure",
    values_to = "value"
  )

# Extract timepoints
d_ext <- d_lf_long %>%
  group_by(id, pex_cycle, measure) %>%
  reframe(
    pre_val = value[flagged_tmpts == "Pre-PEx stable"][1],
    pre_dt = date[flagged_tmpts == "Pre-PEx stable"][1],
    pex_start_val = value[flagged_tmpts == "PEx start"][1],
    pex_start_dt = date[flagged_tmpts == "PEx start"][1],
    pex_end_val = value[flagged_tmpts == "PEx end"][1],
    pex_end_dt = date[flagged_tmpts == "PEx end"][1],
    pex_post_val = value[flagged_tmpts == "Post-PEx"][1],
    pex_post_dt = date[flagged_tmpts == "Post-PEx"][1]
  ) %>%
  mutate(
    diff_pct_start = abs(pre_val - pex_start_val) / pre_val,
    diff_pct = abs(pre_val - pex_end_val) / pre_val,
    diff_abs_start = abs(pre_val - pex_start_val),
    diff_abs = abs(pre_val - pex_end_val),
    rec_lt5_start = diff_pct_start <= 0.05,
    rec_lt5 = diff_pct <= 0.05,
    rec_lt100_start = diff_abs_start <= 0.1,
    rec_lt100 = diff_abs <= 0.1,
    recovery_time_pct = case_when(
      !rec_lt5_start & rec_lt5 ~ as.numeric(pex_post_dt - pre_dt),
      TRUE ~ NA
    ),
    recovery_time_abs = case_when(
      !rec_lt100_start & rec_lt100 ~ as.numeric(pex_post_dt - pre_dt),
      TRUE ~ NA
    )
  )

# Summarise lung function at each timepoint
d_sum <- d_ext %>% 
  select(
    id,
    measure,
    pre_val,
    pex_start_val,
    pex_end_val,
    pex_post_val
  ) %>%
  pivot_longer(
    cols = c(pre_val, pex_start_val, pex_end_val, pex_post_val),
    names_to = "timepoint",
    values_to = "value"
  ) %>%
  mutate(
    timepoint = recode(
      timepoint,
      pre_val = "Pre-PEx stable",
      pex_start_val = "PEx start",
      pex_end_val = "PEx end",
      pex_post_val = "Post-PEx"
    )
  ) %>%
  group_by(measure, timepoint) %>%
  reframe(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    lci = t.test(na.omit(value))$conf.int[1],
    uci = t.test(na.omit(value))$conf.int[2]
  ) %>%
  mutate(
    across(c(mean, sd, lci, uci), ~ round(., 2))
  )

# Summarise by relative and absolute measures
d_sum_lt5 <- d_ext %>%
  group_by(measure) %>%
  reframe(
    n_subj = n_distinct(id),
    n_pex = n(),
    n_drop_start = sum(!rec_lt5_start, na.rm = TRUE),
    n_rec = sum(!rec_lt5_start & rec_lt5, na.rm = TRUE),
    n_notrec = n_drop_start - n_rec,
    n_subj_drop = n_distinct(id[!rec_lt5_start]),
    n_subj_rec = n_distinct(id[!rec_lt5_start & rec_lt5]),
    n_subj_notrec = n_subj_drop - n_subj_rec,
    med_time = median(recovery_time_pct[!rec_lt5_start & rec_lt5], na.rm = TRUE),
    min_time = min(recovery_time_pct[!rec_lt5_start & rec_lt5], na.rm = TRUE),
    max_time = max(recovery_time_pct[!rec_lt5_start & rec_lt5], na.rm = TRUE)
  )

d_sum_lt100 <- d_ext %>%
  group_by(measure) %>%
  reframe(
    n_subj = n_distinct(id),
    n_pex = n(),
    n_drop_start = sum(!rec_lt100_start, na.rm = TRUE),
    n_rec = sum(!rec_lt100_start & rec_lt100, na.rm = TRUE),
    n_notrec = n_drop_start - n_rec,
    n_subj_drop = n_distinct(id[!rec_lt100_start]),
    n_subj_rec = n_distinct(id[!rec_lt100_start & rec_lt100]),
    n_subj_notrec = n_subj_drop - n_subj_rec,
    med_time = median(recovery_time_pct[!rec_lt100_start & rec_lt100], na.rm = TRUE),
    min_time = min(recovery_time_pct[!rec_lt100_start & rec_lt100], na.rm = TRUE),
    max_time = max(recovery_time_pct[!rec_lt100_start & rec_lt100], na.rm = TRUE)
  )

# Save
write_csv(d_sum, "outputs/objs/05_responsiveness/summaries/proxies/sum_lf_all_3d.csv")
write_csv(d_sum_lt5, "outputs/objs/05_responsiveness/summaries/proxies/sum_recovery_lt5_3d.csv")
write_csv(d_sum_lt100, "outputs/objs/05_responsiveness/summaries/proxies/sum_recovery_lt100_3d.csv")

# Difference from baseline (pre-PEx stability)
# e_diff_df <- d_ext %>%
#   mutate(
#     pre_start_diff = pre_val - pex_start_val,
#     pre_end_diff = pre_val - pex_end_val,
#     pre_start_pct = abs(pre_start_diff) / pre_val,
#     pre_end_pct = abs(pre_end_diff) / pre_val
#   )
# 
# e_diff_pct_sum <- e_diff_df %>%
#   pivot_longer(
#     c(pre_start_pct, pre_end_pct),
#     names_to = "timepoint",
#     values_to = "diff"
#   ) %>%
#   group_by(measure, timepoint) %>%
#   reframe(
#     mean = mean(diff, na.rm = TRUE),
#     sd = sd(diff, na.rm = TRUE),
#     lci = t.test(na.omit(diff))$conf.int[1],
#     uci = t.test(na.omit(diff))$conf.int[2],
#     p = t.test(na.omit(diff))$p.value
#   ) %>%
#   mutate(
#     across(c(mean, sd, lci, uci), ~ round(., 2)),
#     p = if_else(p < 0.01, "<0.01", as.character(round(p, 3)))
#   )
# 
# e_diff_abs_sum <- e_diff_df %>%
#   pivot_longer(
#     c(pre_start_diff, pre_end_diff),
#     names_to = "timepoint",
#     values_to = "diff"
#   ) %>%
#   group_by(measure, timepoint) %>%
#   reframe(
#     mean = mean(diff, na.rm = TRUE),
#     sd = sd(diff, na.rm = TRUE),
#     lci = t.test(na.omit(diff))$conf.int[1],
#     uci = t.test(na.omit(diff))$conf.int[2],
#     p = t.test(na.omit(diff))$p.value
#   ) %>%
#   mutate(
#     across(c(mean, sd, lci, uci), ~ round(., 2)),
#     p = if_else(p < 0.01, "<0.01", as.character(round(p, 3)))
#   )
# 
# # Save
# write_csv(e_diff_pct_sum, "outputs/objs/05_responsiveness/summaries/proxies/sum_diff_pct_3d.csv")
# write_csv(e_diff_abs_sum, "outputs/objs/05_responsiveness/summaries/proxies/sum_diff_ml_3d.csv")
# 
# # Plot differences
# y_scales <- list(
#   
#   "FEV1 (L/sec)" = scale_y_continuous(
#     limits = c(-0.6, 0.6),
#     n.breaks = 5
#   ),
#   
#   "FEV1 % predicted" = scale_y_continuous(
#     limits = c(-20, 20),
#     n.breaks = 5
#   ),
#   
#   "FVC (L)" = scale_y_continuous(
#     limits = c(-1, 1),
#     n.breaks = 5
#   ),
#   
#   "FVC % predicted" = scale_y_continuous(
#     limits = c(-40, 40),
#     n.breaks = 5
#   ),
#   
#   "FEV1/FVC" = scale_y_continuous(
#     limits = c(-20, 20),
#     n.breaks = 5
#   ),
#   
#   "FEF2575 (%)" = scale_y_continuous(
#     limits = c(-60, 60),
#     n.breaks = 5
#   ),
#   
#   "PEF (L)" = scale_y_continuous(
#     limits = c(-150, 150),
#     n.breaks = 5
#   )
#   
# )
# 
# d_plot <- e_diff_df %>%
#   pivot_longer(
#     c(pre_start_diff, pre_end_diff),
#     names_to = "timepoint",
#     values_to = "diff"
#   ) %>%
#   mutate(
#     timepoint = recode(
#       timepoint,
#       pre_start_diff = "PEx start",
#       pre_end_diff = "PEx end"
#     ),
#     timepoint = factor(
#       timepoint,
#       levels = c("PEx start", "PEx end")
#     ),
#     measure = case_match(
#       measure,
#       "fef2575" ~ "FEF2575 (%)",
#       "fev1" ~ "FEV1 (L/sec)",
#       "fev1_fvc" ~ "FEV1/FVC",
#       "fev1pp" ~ "FEV1 % predicted",
#       "fvc" ~ "FVC (L)",
#       "fvcpp" ~ "FVC % predicted",
#       "pef" ~ "PEF (L)"
#     ),
#     measure = factor(
#       measure,
#       levels = c(
#         "FEV1 (L/sec)",
#         "FEV1 % predicted",
#         "FVC (L)",
#         "FVC % predicted",
#         "FEV1/FVC",
#         "FEF2575 (%)",
#         "PEF (L)"
#       )
#     )
#   ) %>%
#   ggplot(aes(x = timepoint, y = diff)) +
#   geom_boxplot(width = 0.15, outliers = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
#   facet_wrap(~measure, scales = "free_y", ncol = 4) +
#   ggh4x::facetted_pos_scales(
#     y = y_scales
#   ) +
#   theme_classic() +
#   labs(x = "Timepoint", y = "Difference from pre-PEx stable measure")
# 
# ggsave(
#   "outputs/objs/05_responsiveness/plots/pex_diff_from_baseline_3d.png",
#   d_plot,
#   width = 8,
#   height = 4,
#   units = "in"
# )

# Examine delayed recovery among those who did not recover ---------------------

# Isolate PEx episodes
f_notrec_lt5 <- d_ext %>% 
  filter(
    (!rec_lt5_start & !rec_lt5) | (!rec_lt5_start & is.na(rec_lt5)),
    measure == "fev1"
  )

f_notrec_lt100 <- d_ext %>% 
  filter(
    (!rec_lt100_start & !rec_lt100) | (!rec_lt100_start & is.na(rec_lt100)),
    measure == "fev1"
  )

# Identify next PEx start date
g_nextpex <- c_complete_compl %>% 
  filter(flagged_tmpts == "PEx start") %>% 
  select(id, pex_cycle, next_start = date) %>% 
  mutate(pex_cycle = pex_cycle - 1)

# Join next start to timeline
g_followup <- c_complete_compl %>% 
  left_join(g_nextpex, by = c("id","pex_cycle"))

# Keep only post-PEx measurements
g_post <- g_followup %>% 
  group_by(id, pex_cycle) %>% 
  mutate(
    pex_end_dt = date[flagged_tmpts == "PEx end"][1],
    pre_fev1 = fev1[flagged_tmpts == "Pre-PEx stable"][1]
  ) %>% 
  ungroup() %>% 
  filter(
    date > pex_end_dt,
    is.na(next_start) | date < next_start
  )

# Identify prolonged recovery
g_rec <- g_post %>% 
  mutate(
    diff_pct = abs(pre_fev1 - fev1) / pre_fev1,
    diff_ml = abs(pre_fev1 - fev1),
    rec_lt5 = diff_pct <= 0.05,
    rec_lt100 = diff_ml <= 0.1
  )

# Determine first recovery time
# Set those exceeding 1-year FU as NA, did not recover by end of FU
g_first_rec <- g_rec %>% 
  group_by(id, pex_cycle) %>% 
  reframe(
    rec_date_lt5 = date[rec_lt5][1],
    rec_date_lt100 = date[rec_lt100][1],
    pex_end_dt = first(pex_end_dt)
  ) %>% 
  mutate(
    rec_time_lt5 = as.numeric(rec_date_lt5 - pex_end_dt),
    rec_time_lt100 = as.numeric(rec_date_lt100 - pex_end_dt),
    across(
      c(rec_time_lt5, rec_time_lt100),
      ~ case_when(. >= 365 ~ NA, TRUE ~ .)
    )
  )

# Restrict to participants who did not recover within analysis window
g_lt5_notrec <- f_notrec_lt5 %>% 
  select(id, pex_cycle) %>% 
  left_join(g_first_rec, by = c("id","pex_cycle"))

g_lt100_notrec <- f_notrec_lt100 %>% 
  select(id, pex_cycle) %>% 
  left_join(g_first_rec, by = c("id","pex_cycle"))

# Summarise
# <=5% threshold
g_sum_lt5_eventual <- g_lt5_notrec %>% 
  distinct() %>% 
  reframe(
    n_subj = n_distinct(id),
    n_pex = n(),
    n_eventual_rec = sum(!is.na(rec_date_lt5)),
    n_never_rec = sum(is.na(rec_date_lt5)),
    n_subj_eventual_rec = n_distinct(id[!is.na(rec_date_lt5)]),
    n_subj_never_rec = n_distinct(id[is.na(rec_date_lt5)]),
    med_time = median(rec_time_lt100[!is.na(rec_date_lt5)], na.rm = TRUE),
    min_time = min(rec_time_lt100[!is.na(rec_date_lt5)], na.rm = TRUE),
    max_time = max(rec_time_lt100[!is.na(rec_date_lt5)], na.rm = TRUE)
  )

# <=100ml threshold
g_sum_lt100_eventual <- g_lt100_notrec %>% 
  distinct() %>% 
  reframe(
    n_subj = n_distinct(id),
    n_pex = n(),
    n_eventual_rec = sum(!is.na(rec_date_lt100)),
    n_never_rec = sum(is.na(rec_date_lt100)),
    n_subj_eventual_rec = n_distinct(id[!is.na(rec_date_lt100)]),
    n_subj_never_rec = n_distinct(id[is.na(rec_date_lt100)]),
    med_time = median(rec_time_lt100[!is.na(rec_date_lt100)], na.rm = TRUE),
    min_time = min(rec_time_lt100[!is.na(rec_date_lt100)], na.rm = TRUE),
    max_time = max(rec_time_lt100[!is.na(rec_date_lt100)], na.rm = TRUE)
  )

# # Save
write_csv(
  g_sum_lt5_eventual,
  "outputs/objs/05_responsiveness/summaries/proxies/sum_eventual_recovery_lt5_3d.csv"
)

write_csv(
  g_sum_lt100_eventual,
  "outputs/objs/05_responsiveness/summaries/proxies/sum_eventual_recovery_lt100_3d.csv"
)

# Reasons for non-recovery among remaining participants and PEx ----------------

# Get dropout and last visit dates
h_denom <- a_imp_denom %>% 
  select(id, v1, v2, v3, v4, v5, dropout_date) %>% 
  pivot_longer(contains("v")) %>% 
  group_by(id) %>% 
  mutate(last_visit = max(value)) %>% 
  select(id, dropout_date, last_visit) %>% 
  distinct()

# 41 participants, 54 PEx 
h_notrec_pct <- g_lt5_notrec %>%
  filter(is.na(rec_date_lt5)) %>% 
  select(id, pex_cycle, rec_date_lt5, pex_end_dt)

# 36 participants, 48 PEx 
h_notrec_ml <- g_lt100_notrec %>% 
  filter(is.na(rec_date_lt100)) %>% 
  select(id, pex_cycle, rec_date_lt100, pex_end_dt)

# Identify last measurement before follow-up ends
h_last_measure <- c_complete_compl %>% 
  group_by(id, pex_cycle) %>% 
  reframe(
    last_meas_dt = max(date, na.rm = TRUE)
  )

# Join data together
# Join dropout dates
h_notrec_pct <- h_notrec_pct %>% 
  left_join(h_last_measure, by = c("id", "pex_cycle")) %>% 
  left_join(g_nextpex, by = c("id", "pex_cycle")) %>% 
  left_join(h_denom, by = "id")

h_notrec_ml <- h_notrec_ml %>% 
  left_join(h_last_measure, by = c("id","pex_cycle")) %>% 
  left_join(g_nextpex, by = c("id","pex_cycle")) %>% 
  left_join(h_denom, by = "id")

# Flag reasons recovery never recovered
h_notrec_pct <- h_notrec_pct %>% 
  distinct() %>% 
  mutate(
    reason = case_when(
      !is.na(dropout_date) & as.numeric(dropout_date - last_meas_dt) <= 7 ~ "Dropped out",
      is.na(next_start) ~ "No additional data",
      !is.na(next_start) ~ "Next PEx occurred",
      TRUE ~ "No additional data"
    )
  )

h_notrec_ml <- h_notrec_ml %>% 
  distinct() %>% 
  mutate(
    reason = case_when(
      !is.na(dropout_date) & as.numeric(dropout_date - last_meas_dt) <= 7 ~ "Dropped out",
      is.na(next_start) ~ "No additional data",
      !is.na(next_start) ~ "Next PEx occurred",
      TRUE ~ "No additional data"
    )
  )

h_reason_sum_pct <- h_notrec_pct %>% 
  group_by(reason) %>% 
  reframe(
    n_pex = n(),
    n_subj = n_distinct(id)
  )

h_reason_sum_ml <- h_notrec_ml %>% 
  group_by(reason) %>% 
  reframe(
    n_pex = n(),
    n_subj = n_distinct(id)
  )

write_csv(
  h_reason_sum_pct,
  "outputs/objs/05_responsiveness/summaries/proxies/sum_nonrecovery_reasons_lt5_3d.csv"
)

write_csv(
  h_reason_sum_ml,
  "outputs/objs/05_responsiveness/summaries/proxies/sum_nonrecovery_reasons_lt100_3d.csv"
)
