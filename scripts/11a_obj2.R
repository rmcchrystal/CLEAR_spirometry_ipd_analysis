
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
      "pef"      ~ "PEF (L)"
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
          "pef"     ~ "PEF (L)"
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
      "pef"     ~ "PEF (L)"
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
  select(id, date, fev1, fvc, fev1pp, fvcpp, fef2575, pef)

f_clinic <- readRDS("processed_data/analysis_data.rds") %>% 
  select(site, id, pex_orig) %>% 
  unnest(pex_orig) %>% 
  filter(data == "clinic") %>% 
  select(id, date, fev1, fvc, fev1pp, fvcpp, fef2575, pef)

f_base <- d_visits %>%
  select(id,visit,date) %>%
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
        filter(!is.na(home), !is.na(clinic))
    })
  )

f_icc_nest <- f_cross %>%
  mutate(
    results = map(data, ~ {
      
      df <- .x
      if(nrow(df)<5) return(NULL)
      n <- length(unique(df$id))
      
      icc_fit <- icc(
        df %>% select(clinic,home),
        model = "twoway",
        type = "agreement",
        unit = "single"
      )
      
      diff_vec <- df$clinic - df$home
      
      long <- df %>%
        pivot_longer(
          cols=c(home,clinic),
          names_to="method",
          values_to="value"
        )
    
      fit <- lmer(log(value)~1+(1|id),data=long,REML=TRUE)
      sigma_w <- sigma(fit)
      cv_hat  <- sqrt(exp(sigma_w^2)-1)
      
      ci_sigma <- suppressMessages(
        confint(fit,parm=".sigma",method="profile")
      )
      
      tibble(
        n          = n,
        icc        = icc_fit$value,
        icc_lo     = icc_fit$lbound,
        icc_hi     = icc_fit$ubound,
        cv_within  = cv_hat,
        cv_lo      = sqrt(exp(ci_sigma[1]^2)-1),
        cv_hi      = sqrt(exp(ci_sigma[2]^2)-1),
        mean_diff  = mean(diff_vec),
        sd_diff    = sd(diff_vec),
        lo_diff    = mean(diff_vec) - 1.96 * sd(diff_vec) / sqrt(n),
        hi_diff    = mean(diff_vec) + 1.96 * sd(diff_vec) / sqrt(n)
      )
      
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
    across(contains("cv"), ~ round(. * 100, 3))
  ) %>%
  arrange(measure, visit)

write_csv(
  f_icc_summary,
  "outputs/objs/02_variation/summaries/icc_cv_diff_summary.csv"
)
