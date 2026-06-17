
# Objective 2
# Explore variation over time in lung function measurements between 
# home and clinic spirometry performed on the same day and within the same week:
# a.) Irrespective of disease status (i.e. stable or during an exacerbation)
# b.) During clinically stable periods of bronchiectasis
# c.) Between the start and end of an exacerbation
# Uses same data as obj1

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")
library(tidyverse)
library(irr)
library(lme4)

# Function to check for NAs
check_na <- function(df) {
  df %>% reframe(across(everything(), ~ sum(is.na(.))))
}

# Import data
a_imp_lst <- readRDS("processed_data/obj1_data.rds")

# Import patient denominator
a_imp_denom <- readRDS("created_metadata/denom.rds") %>% select(id, v1:v5)

# Add names as variables
b_lst <- lapply(names(a_imp_lst), function(i) {
  nm <- i
  a_imp_lst[[i]] %>% mutate(type = nm)
})

# Bind
c_df <- bind_rows(b_lst)

# Prepare ----------------------------------------------------------------------

# Long data
d_long <- c_df %>% 
  select(
    id, date, data, fev1, fev1pp, fvc, fvcpp, fev1_fvc, fef2575, pef, type
  ) %>% 
  pivot_longer(fev1:pef) %>% 
  distinct() %>% 
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
  filter(!grepl("PEx", type))

# Join visit timepoints
d_visits <- d_long %>% 
  left_join(
    a_imp_denom %>% 
      pivot_longer(v1:v5) %>% 
      rename(
        visit = name, 
        date = value
      )
  ) %>% 
  mutate(
    visit = gsub("v", "V", visit),
    visit = factor(visit, levels = c("V1", "V2", "V3", "V4", "V5"))
  )

# Setup dataset
d_df <- d_visits %>%
  filter(!is.na(visit)) %>%
  group_by(data, name, type, tmpt, visit) %>% 
  reframe(
    n    = length(unique(id)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    se   = sd/sqrt(n),
    cv   = sd/mean
  ) %>% 
  mutate(
    nm = case_match(
      name,
      "fev1"     ~ "FEV1 (L/sec)",
      "fvc"      ~ "FVC (L)",
      "fev1pp"   ~ "FEV1 % predicted",
      "fvcpp"    ~ "FVC % predicted",
      "fev1_fvc" ~ "FEV1/FVC",
      "fef2575"  ~ "FEF25-75 (%)",
      "pef"      ~ "PEF (L/min)"
    ),
    data = case_match(
      data,
      "clinic" ~ "Clinic",
      "home"   ~ "Home"
    )
  )

d_sum <- d_df %>% 
  mutate(across(mean:cv, ~ round(., 2))) %>% 
  arrange(tmpt, type, name, data)

# write_csv(
#   d_sum,
#   "outputs/objs/02_variation/summaries/clinic_home_summaries.csv"
# )

# Plots ------------------------------------------------------------------------

# Whisker plots
e_plot <- d_df %>% 
  group_by(tmpt, name) %>% 
  nest() %>% 
  mutate(
    data_plot = map2(data, name, ~ {
      
      p <- .x %>% 
        mutate(
          type = factor(
            type,
            levels = c(
              "Irrespective of disease status",
              "During stable periods",
              "Between PEx start and end"
            )
          )
        ) %>% 
        ggplot(aes(x = visit, y = mean, colour = data)) +
        geom_point(position = position_dodge(width = .6), size = 2) +
        geom_errorbar(
          aes(ymin = mean - se, ymax = mean + se),
          width = 0.2,
          position = position_dodge(width = .6)
        ) +
        facet_wrap(~type, scales = "free") +
        theme_classic() +
        labs(
          x = "Visit",
          y = unique(.x$nm),
          colour = "Spirometer"
        )
      
      if (.y == "fev1") {
        p + scale_y_continuous(limits = c(2, 2.5), n.breaks = 6)
      } else if (.y == "fvc") {
        p + scale_y_continuous(limits = c(2.9, 3.5), n.breaks = 6)
      } else if (.y == "fev1pp") {
        p + scale_y_continuous(limits = c(75, 92), n.breaks = 6)
      } else if (.y == "fvcpp") {
        p + scale_y_continuous(limits = c(85, 100), n.breaks = 6)
      } else if (.y == "fef2575") {
        p + scale_y_continuous(limits = c(70, 130), n.breaks = 6)
      } else if (.y == "pef") {
        p + scale_y_continuous(limits = c(300, 400), n.breaks = 6)
      }
    }),
    plot_name = paste0(name, "_", tmpt),
    plot_path = paste0(
      "outputs/objs/02_variation/plots/",
      plot_name, ".png"
    )
  )

# walk2(
#   e_plot$plot_path,
#   e_plot$data_plot,
#   ~ ggsave(.x, .y, height = 4, width = 6, units = "in")
# )

# Change from baseline in weekly spirometry ------------------------------------

# Save for obj5
d_visits %>%
  filter(
    visit == "V1" 
    & name == "fev1"
  ) %>%
  na.omit() %>%
  saveRDS("processed_data/obj5_bl.rds")

# Baseline values
g_bl <- d_visits %>% 
  filter(
    visit == "V1",
    name %in% c("fev1", "fvc", "fev1pp", "fvcpp", "fef2575", "pef")
  ) %>% 
  pivot_wider(names_from = data, values_from = value) %>% 
  transmute(id, name, home_bl = home) %>% 
  distinct() %>% 
  drop_na()

g_home <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig) %>% 
  filter(data == "home") %>% 
  select(id, date, fev1, fvc, fev1pp, fvcpp, fef2575, pef, disease_status) %>% 
  pivot_longer(
    cols = c(fev1, fvc, fev1pp, fvcpp, fef2575, pef),
    names_to = "name",
    values_to = "value"
  ) %>% 
  left_join(
    d_visits %>% distinct(id, date, visit),
    by = c("id","date")
  ) %>% 
  fill(visit, .direction = "down") %>% 
  group_by(id, name) %>%
  arrange(date) %>% 
  mutate(
    first    = first(date),
    week_num = as.integer(floor((date - first) / 7))
  ) %>% 
  ungroup() %>% 
  left_join(g_bl, by = c("id","name")) %>% 
  mutate(
    chng_bl = home_bl - value
  )

g_scenarios <- crossing(
  measure = c("fev1", "fvc", "fev1pp", "fvcpp", "fef2575", "pef"),
  popn    = c("all", "stable"),
  tscale  = c("weekly", "visit")
)

g_nested <- g_scenarios %>% 
  mutate(
    data = map2(measure, popn, ~ {
      
      df <- g_home %>% filter(name == .x)
      
      if (.y == "stable") {
        df <- df %>% filter(disease_status == "stable")
      }
      
      df %>% filter(week_num <= 52)
    })
  ) %>% 
  mutate(
    summary = map2(data, tscale, ~ {
      
      if (.y == "weekly") {
        
        .x %>%
          group_by(week_num) %>%
          reframe(
            n    = length(unique(id)),
            mean = mean(chng_bl, na.rm = TRUE),
            sd   = sd(chng_bl, na.rm = TRUE),
            se   = sd/sqrt(n),
            lo   = mean - 1.96 * se,
            hi   = mean + 1.96 * se
          ) %>%
          mutate(
            across(mean:hi, ~ if_else(week_num == 0, 0, .))
          )
        
      } else {
        
        .x %>%
          group_by(visit) %>%
          reframe(
            n    = length(unique(id)),
            mean = mean(chng_bl, na.rm = TRUE),
            sd   = sd(chng_bl, na.rm = TRUE),
            se   = sd/sqrt(n),
            lo   = mean - 1.96 * se,
            hi   = mean + 1.96 * se
          ) %>%
          mutate(
            across(mean:hi, ~ if_else(visit == "V1", 0, .))
          )
      }
    })
  )

g_nested <- g_nested %>% 
  mutate(
    plot = pmap(
      list(summary, measure, tscale),
      function(.x, .m, .t) {
        
        measure_label <- case_match(
          .m,
          "fev1"    ~ "FEV1 (L/sec)",
          "fvc"     ~ "FVC (L)",
          "fev1pp"  ~ "FEV1 % predicted",
          "fvcpp"   ~ "FVC % predicted",
          "fef2575" ~ "FEF25-75 (%)",
          "pef"     ~ "PEF (L/min)"
        )
        
        if (.t == "weekly") {
          
          p <- ggplot(.x, aes(x = week_num, y = mean)) +
            geom_line() +
            geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme_classic() +
            labs(
              x = "Follow-up (weeks)",
              y = paste("Change from baseline", measure_label)
            )
          
        } else {
          
          p <- ggplot(.x, aes(x = as.numeric(visit), y = mean)) +
            geom_point(size = 2) +
            geom_errorbar(
              aes(ymin = lo, ymax = hi),
              width = 0.2
            ) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme_classic() +
            labs(
              x = "Visit",
              y = paste("Change from baseline", measure_label)
            )
        }
        
        if (.m == "fev1") {
          p + scale_y_continuous(limits = c(-0.2, 0.4), n.breaks = 6)
        } else if (.m == "fvc") {
          p + scale_y_continuous(limits = c(-0.2, 0.5), n.breaks = 6)
        } else if (.m == "fev1pp") {
          p + scale_y_continuous(limits = c(-10, 15), n.breaks = 6)
        } else if (.m == "fvcpp") {
          p + scale_y_continuous(limits = c(-10, 12), n.breaks = 6)
        } else {
          p
        }
      }
    ),
    file_stub = paste0(
      measure, "_",
      tscale, "_change_from_baseline",
      if_else(popn == "stable","_stable","")
    ),
    csv_path = paste0(
      "outputs/objs/02_variation/summaries/",
      file_stub, ".csv"
    ),
    plot_path = paste0(
      "outputs/objs/02_variation/plots/",
      file_stub, ".png"
    )
  )

# Combine visit-based and follow-up based change plots
g_combined <- g_nested %>%
  select(measure, popn, tscale, plot) %>%
  pivot_wider(
    names_from = tscale,
    values_from = plot
  ) %>%
  mutate(
    combined_plot = map2(
      visit, weekly,
      ~ cowplot::plot_grid(
        .x, .y,
        labels = c("A", "B"),
        ncol = 2,
        align = "v",
        rel_heights = c(1, 1)
      )
    ),
    file_stub = paste0(
      measure, "_change_from_baseline_combined",
      if_else(popn == "stable", "_stable", "")
    ),
    plot_path = paste0(
      "outputs/objs/02_variation/plots/",
      file_stub, ".png"
    )
  ) %>% 
  filter(popn == "all")

walk2(
  g_combined$plot_path,
  g_combined$combined_plot,
  ~ ggsave(.x, .y, height = 4, width = 8, units = "in")
)

g_all_summaries <- g_nested %>%
  select(measure, popn, tscale, summary) %>%
  unnest(summary) %>%
  mutate(
    measure = case_match(
      measure,
      "fev1"    ~ "FEV1 (L/sec)",
      "fvc"     ~ "FVC (L)",
      "fev1pp"  ~ "FEV1 % predicted",
      "fvcpp"   ~ "FVC % predicted",
      "fef2575" ~ "FEF25-75 (%)",
      "pef"     ~ "PEF (L/min)"
    )
  )

write_csv(
  g_all_summaries,
  "outputs/objs/02_variation/summaries/change_from_baseline_all.csv"
)

# walk2(g_nested$summary, g_nested$csv_path, write_csv)
# 
# walk2(
#   g_nested$plot_path,
#   g_nested$plot,
#   ~ ggsave(.x, .y, height = 4, width = 4, units = "in")
# )

# ICC --------------------------------------------------------------------------

# Define iteration setup
f_setup <- tribble(
  ~measure,
  "fev1",
  "fvc",
  "fev1pp",
  "fvcpp",
  "fef2575",
  "pef"
)

f_home <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig) %>% 
  filter(data == "home") %>% 
  select(id, date, status, fev1, fvc, fev1pp, fvcpp, fef2575, pef)

f_clinic <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig) %>% 
  filter(data == "clinic") %>% 
  select(id, date, status, fev1, fvc, fev1pp, fvcpp, fef2575, pef)

f_base <- d_visits %>%
  select(id, visit, date) %>%
  distinct() %>%
  filter(!is.na(visit)) %>%
  left_join(f_home, by = c("id", "date")) %>%
  rename_with(
    ~ paste0(.x, "_home"), 
    c(fev1, fvc, fev1pp, fvcpp, fef2575, pef)
  ) %>%
  left_join(f_clinic, by = c("id", "date")) %>%
  rename_with(
    ~ paste0(.x, "_clinic"), 
    c(fev1, fvc, fev1pp, fvcpp, fef2575, pef)
  )

f_cross <- f_setup %>%
  crossing(visit = unique(f_base$visit)) %>%
  mutate(
    data = map2(measure, visit, ~ {
      home_col   <- paste0(.x,"_home")
      clinic_col <- paste0(.x,"_clinic")
      f_base %>%
        filter(visit == .y) %>%
        select(
          id,
          home = all_of(home_col),
          clinic = all_of(clinic_col)
        ) %>%
        filter(!is.na(home), !is.na(clinic)) %>% 
        left_join(
          readRDS("processed_data/analysis_data.rds") %>% 
            select(site, id, pex_orig) %>% 
            unnest(pex_orig) %>% 
            ungroup() %>% 
            select(id, status)
        )
    })
  )

# Helper function to fit linear mixed effects model for ICC
fit_icc_lmm <- function(df) {
  
  n <- length(unique(df$id))
  if(n < 5) return(NULL)
  
  diff_vec <- df$clinic - df$home
  
  long <- df %>%
    pivot_longer(
      cols = c(home,clinic),
      names_to = "method",
      values_to = "value"
    ) %>%
    mutate(
      method = factor(method, levels = c("clinic", "home"))
    )
  
  # Same specification for ICC and CV
  # Method as fixed effect, random participant-level intercepts
  icc_fit <- lmer(value ~ method + (1 | id), data = long, REML=TRUE)
  
  vc <- as.data.frame(VarCorr(icc_fit))
  
  var_id <- vc %>% 
    filter(grp == "id") %>% 
    pull(vcov)
  
  var_resid <- sigma(icc_fit)^2
  
  icc_lmm <- var_id / (var_id + var_resid)
  
  icc_boot <- tryCatch(
    suppressMessages(
      bootMer(
        icc_fit,
        FUN = function(fit) {
          
          vc <- as.data.frame(VarCorr(fit))
          
          var_id <- vc %>% 
            filter(grp == "id") %>% 
            pull(vcov)
          
          var_resid <- sigma(fit)^2
          
          var_id / (var_id + var_resid)
        },
        nsim = 1000,
        type = "parametric",
        use.u = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  icc_ci <- if(is.null(icc_boot)) {
    c(NA_real_, NA_real_)
  } else {
    quantile(icc_boot$t[,1], probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  fit <- lmer(log(value) ~ method + (1 | id), data = long, REML=TRUE)
  sigma_w <- sigma(fit)
  
  ci_sigma <- suppressMessages(
    confint(fit, parm = ".sigma", method = "profile")
  )
  
  sigma_lo <- ci_sigma[".sigma", "2.5 %"]
  sigma_hi <- ci_sigma[".sigma", "97.5 %"]
  
  cv_within <- sqrt(exp(sigma_w^2) - 1)
  cv_lo     <- sqrt(exp(sigma_lo^2) - 1)
  cv_hi     <- sqrt(exp(sigma_hi^2) - 1)
  
  tibble(
    n = n,
    icc = icc_lmm,
    icc_lo = icc_ci[1],
    icc_hi = icc_ci[2],
    cv_within = 100 * cv_within,
    cv_lo = 100 * cv_lo,
    cv_hi = 100 * cv_hi,
    mean_diff = mean(diff_vec, na.rm = TRUE),
    sd_diff = sd(diff_vec, na.rm = TRUE),
    mean_diff_lo = mean(diff_vec, na.rm = TRUE) - 1.96 * sd(diff_vec, na.rm = TRUE) / sqrt(n),
    mean_diff_hi = mean(diff_vec, na.rm = TRUE) + 1.96 * sd(diff_vec, na.rm = TRUE) / sqrt(n),
    loa_lower = mean(diff_vec, na.rm = TRUE) - 1.96 * sd(diff_vec, na.rm = TRUE),
    loa_upper = mean(diff_vec, na.rm = TRUE) + 1.96 * sd(diff_vec, na.rm = TRUE),
    random_intercept = "participant",
    random_slope = "not fitted"
  )
}

set.seed(1234)

f_icc_nest <- f_cross %>%
  mutate(
    results = map(data, ~ {
      fit_icc_lmm(.x)
    })
  )

f_icc_summary <- f_icc_nest %>%
  select(measure, visit, results) %>%
  unnest(results) %>%
  mutate(
    measure = case_match(
      measure,
      "fev1"    ~ "FEV1 (L/sec)",
      "fvc"     ~ "FVC (L)",
      "fev1pp"  ~ "FEV1 % predicted",
      "fvcpp"   ~ "FVC % predicted",
      "fef2575" ~ "FEF25-75 (%)",
      "pef"     ~ "PEF (L)"
    ),
    across(contains("cv"), ~ round(., 3)),
    across(starts_with("icc"), ~ round(., 3))
  ) %>%
  arrange(measure, visit)

write_csv(
  f_icc_summary,
  "outputs/objs/02_variation/summaries/icc_cv_diff_summary.csv"
)

# Completers-only
# f_icc_summary_compl <- f_icc_nest %>% 
#   select(measure, visit, results_compl) %>%
#   unnest(results_compl) %>% 
#   ungroup() %>% 
#   mutate(
#     measure = case_match(
#       measure,
#       "fev1"    ~ "FEV1 (L/sec)",
#       "fvc"     ~ "FVC (L)",
#       "fev1pp"  ~ "FEV1 % predicted",
#       "fvcpp"   ~ "FVC % predicted",
#       "fef2575" ~ "FEF25-75 (%)",
#       "pef"     ~ "PEF (L/min)"
#     ),
#     across(contains("cv"), ~ round(., 3))
#   ) %>%
#   arrange(measure, visit)

# Combine summaries
# f_icc_summary_all <- f_icc_summary %>% 
#   mutate(type = "all") %>% 
#   full_join(f_icc_summary_compl %>% mutate(type = "completers")) %>% 
#   arrange(measure, visit, type) %>% 
#   select(type, everything())

# Join with full population summary
# f_icc_sum_all <- f_icc_summary %>% 
#   # pivot_longer(-c(measure, visit, n)) %>% 
#   mutate(popn = "FAS") %>% 
#   full_join(
#     f_icc_summary_compl %>% 
#       # pivot_longer(-c(measure, visit, n)) %>% 
#       mutate(popn = "Completers only") 
#   ) %>% 
#   arrange(measure, visit, popn) %>% 
#   select(popn, everything())

# write_csv(
#   f_icc_sum_all,
#   "outputs/objs/02_variation/summaries/icc_cv_diff_summary_vs_compl.csv"
# )

# Plot
# f_icc_plot <- f_icc_summary %>%
#   mutate(
#     measure = factor(
#       measure,
#       levels = c(
#         "FEV1 (L/sec)", "FEV1 % predicted", "FVC (L)", "FVC % predicted",
#         "FEF25-75 (%)", "PEF (L/min)"
#       )
#     ),
#     popn = factor(popn, levels = c("FAS", "Completers only")),
#     visit = case_match(
#       visit,
#       "V1" ~ "V1 (Baseline)",
#       "V2" ~ "V2 (Week 2)",
#       "V3" ~ "V3 (Week 8)",
#       "V4" ~ "V4 (Week 26)",
#       "V5" ~ "V5 (Week 52)"
#     )
#   ) %>%
#   group_by(measure) %>%
#   nest() %>%
#   mutate(
#     plot_nm = paste0(
#       "plot_",
#       str_trim(gsub("%|\\-|\\(|\\)|/| ", "", measure))
#     ),
#     plot = map2(data, measure, ~ {
#       .x %>%
#         ggplot(aes(x = visit, y = cv_within, colour = popn)) +
#         geom_point(size = 0.5) +
#         geom_linerange(aes(ymin = cv_lo, ymax = cv_hi)) +
#         theme_classic() +
#         labs(
#           x = "Visit",
#           y = "Coefficient of variation (%)",
#           colour = "Analysis set"
#         ) +
#         theme(
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "bottom"
#         ) +
#         ggtitle(label = .y)
#     })
#   )
# 
# walk2(
#   .x = f_icc_plot$plot,
#   .y = f_icc_plot$plot_nm,
#   ~ ggsave(
#     paste0("outputs/objs/02_variation/plots/cv_compare/", .y, ".png"),
#     .x,
#     width = 4,
#     height = 3,
#     units = "in"
#   )
# )

# Annualised rate of change ----------------------------------------------------

# Supplementary analysis
# Baseline visit dates
h_v1 <- d_visits %>% 
  filter(visit == "V1") %>% 
  select(id, bl_date = date) %>% 
  distinct()

# Baseline paired home and clinic measurements
h_bl <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig) %>% 
  filter(data %in% c("clinic", "home")) %>% 
  select(id, date, data, fev1, fvc, fev1pp, fvcpp, fef2575, pef, disease_status) %>% 
  inner_join(
    h_v1 %>% rename(date = bl_date),
    by = c("id", "date")
  ) %>% 
  rename(bl_date = date) %>% 
  pivot_longer(
    cols = c(fev1, fvc, fev1pp, fvcpp, fef2575, pef),
    names_to = "measure",
    values_to = "value"
  ) %>% 
  select(id, bl_date, data, measure, value) %>% 
  distinct() %>% 
  pivot_wider(names_from = data, values_from = value) %>% 
  filter(
    !is.na(home),
    !is.na(clinic)
  ) %>% 
  select(id, measure, bl_date, home_bl = home, clinic_bl = clinic)

# Check baseline denominator
h_bl_check <- h_bl %>% 
  count(measure, name = "n_bl")

# Longitudinal home and clinic spirometry
h_long <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig) %>% 
  filter(data %in% c("clinic", "home")) %>% 
  select(id, date, data, fev1, fvc, fev1pp, fvcpp, fef2575, pef, disease_status) %>% 
  pivot_longer(
    cols = c(fev1, fvc, fev1pp, fvcpp, fef2575, pef),
    names_to = "measure",
    values_to = "value"
  ) %>% 
  inner_join(
    h_bl %>% select(id, measure, bl_date),
    by = c("id", "measure")
  ) %>% 
  filter(!is.na(value)) %>% 
  mutate(
    time_yr = as.numeric(date - bl_date) / 365.25,
    data = factor(data, levels = c("clinic", "home"))
  ) %>% 
  filter(
    time_yr >= 0,
    time_yr <= 1
  )

# Define scenarios
h_scenarios <- crossing(
  measure = c("fev1", "fvc", "fev1pp", "fvcpp", "fef2575", "pef"),
  popn    = c("all", "stable")
)

# Helper function to estimate annualised slopes
fit_slope_lmm <- function(df) {
  
  n <- length(unique(df$id))
  if(n < 5) return(NULL)
  
  if(length(unique(df$data)) < 2) return(NULL)
  
  fit <- lmer(value ~ time_yr * data + (1 | id), data = df, REML=FALSE)
  
  b <- fixef(fit)
  v <- vcov(fit)
  
  clinic_slope <- b["time_yr"]
  home_slope   <- b["time_yr"] + b["time_yr:datahome"]
  slope_diff   <- -b["time_yr:datahome"]
  
  clinic_se <- sqrt(v["time_yr", "time_yr"])
  
  home_se <- sqrt(
    v["time_yr", "time_yr"] +
      v["time_yr:datahome", "time_yr:datahome"] +
      2 * v["time_yr", "time_yr:datahome"]
  )
  
  slope_diff_se <- sqrt(v["time_yr:datahome", "time_yr:datahome"])
  
  tibble(
    n = n_distinct(df$id),
    n_obs = nrow(df),
    n_clinic = length(unique(df$id[df$data == "clinic"])),
    n_home = length(unique(df$id[df$data == "home"])),
    clinic_slope = clinic_slope,
    clinic_lo = clinic_slope - 1.96 * clinic_se,
    clinic_hi = clinic_slope + 1.96 * clinic_se,
    home_slope = home_slope,
    home_lo = home_slope - 1.96 * home_se,
    home_hi = home_slope + 1.96 * home_se,
    slope_diff = slope_diff,
    slope_diff_lo = slope_diff - 1.96 * slope_diff_se,
    slope_diff_hi = slope_diff + 1.96 * slope_diff_se,
    slope_diff_p = 2 * pnorm(abs(slope_diff / slope_diff_se), lower.tail = FALSE)
  )
}

# Fit annualised slope models
h_slope <- h_scenarios %>% 
  mutate(
    data = map2(measure, popn, ~ {
      
      df <- h_long %>% filter(measure == .x)
      
      if(.y == "stable") {
        df <- df %>% filter(disease_status == "stable")
      }
      
      df
    }),
    results = map(data, fit_slope_lmm)
  )

# Summarise annualised slopes
h_slope_summary <- h_slope %>% 
  select(measure, popn, results) %>% 
  unnest(results) %>% 
  mutate(
    measure = case_match(
      measure,
      "fev1"    ~ "FEV1 (L/sec)",
      "fvc"     ~ "FVC (L)",
      "fev1pp"  ~ "FEV1 % predicted",
      "fvcpp"   ~ "FVC % predicted",
      "fef2575" ~ "FEF25-75 (%)",
      "pef"     ~ "PEF (L/min)"
    ),
    across(
      c(clinic_slope:slope_diff_hi),
      ~ round(., 3)
    ),
    slope_diff_p = round(slope_diff_p, 3)
  ) %>% 
  arrange(measure, popn)

write_csv(
  h_slope_summary,
  "outputs/objs/02_variation/summaries/annualised_slope_summary.csv"
)
