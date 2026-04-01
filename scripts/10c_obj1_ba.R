
# Objective 1
# Generate Bland-Altman plots (descriptive version)

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")

library(tidyverse)
z_filedir <- "outputs/objs/01_agreement/summaries"

# Import materials
a_imp_ba <- readRDS("scratch_data/bland_altman_data.rds")
a_imp_ba_cor <- read_csv(paste0(z_filedir, "/baltman_correlation.csv"))

# Define layout
facet_levels <- c(
  "± 7 days",
  "Same day during stability",
  "Same day during exacerbation",
  "± 7 days during exacerbation"
)

facet_lookup <- tribble(
  ~type_pattern,      ~facet_label,
  "same_week_all",    "± 7 days",
  "same_day_stable",  "Same day during stability",
  "same_day_pex",     "Same day during exacerbation",
  "same_week_pex",    "± 7 days during exacerbation",
)

# Prepare correlation tibble
b_ba_cor <- a_imp_ba_cor %>%
  mutate(across(cor_pearsons:cor_p, ~ round(., 2))) %>%
  rename(r = cor_pearsons, p = cor_p) %>% 
  mutate(
    type = case_when(
      type == "Irrespective of disease status" & tmpt == "same_week" ~ "same_week_all",
      type == "During stable periods"          & tmpt == "same_day"  ~ "same_day_stable",
      type == "Between PEx start and end"      & tmpt == "same_day"  ~ "same_day_pex",
      type == "Between PEx start and end" & tmpt == "same_week" ~ "± 7 days during PEx",
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(type))

# Define measures
measure_specs <- tribble(
  ~name,        ~xlab,                                              ~ylab,                                                   ~xlims,        ~ylims,        ~breaks,
  "fev1",       "Mean of home and clinic FEV1 (L/s)",               "Difference between home and clinic FEV1 (L/s)",         c(0, 6),        c(-4,4),       5,
  "fvc",        "Mean of home and clinic FVC (L)",                  "Difference between home and clinic FVC (L)",            c(0, 8),        c(-4,6),       5,
  "fev1pp",     "Mean of home and clinic FEV1 % predicted",         "Difference between home and clinic FEV1 % predicted",   c(0, 200),      c(-100,100),   5,
  "fvcpp",      "Mean of home and clinic FVC % predicted",          "Difference between home and clinic FVC % predicted",    c(0, 200),      c(-100,100),   5,
  "fev1_fvc",   "Mean of home and clinic FEV1/FVC",                 "Difference between home and clinic FEV1/FVC",           c(40 ,130),     c(-60,60),     5,
  "fef2575",    "Mean of home and clinic FEF25-75 (%)",             "Difference between home and clinic FEF25-75 (%)",       c(0, 350),      c(-200, 200),  5,
  "pef",        "Mean of home and clinic PEF (L)",                  "Difference between home and clinic PEF (L)",            c(0, 900),      c(-500, 500),  5
)

# Create scenario grid
scenario_grid <- measure_specs

# Nest and plot
ba_nested <- scenario_grid %>%
  mutate(
    
    # ---- Filter + recode data ----
    data = map(name, ~ {
      
      a_imp_ba %>%
        filter(
          name == .x,
          type %in% facet_lookup$type_pattern
        ) %>%
        left_join(facet_lookup, by = c("type" = "type_pattern")) %>%
        mutate(
          type = factor(facet_label, levels = facet_levels)
        )
    }),
    
    # ---- Compute BA features ----
    features = map2(data, name, ~ {
      
      .x %>%
        distinct(type, bias, hi, lo) %>%
        left_join(
          b_ba_cor %>%
            filter(name == .y) %>%
            filter(type %in% facet_lookup$type_pattern) %>%
            left_join(facet_lookup, by = c("type" = "type_pattern")) %>%
            select(type = facet_label, r, p) %>% 
            mutate(type = factor(type, levels = facet_levels)),
          by = "type"
        ) %>%
        mutate(
          across(bias:lo, ~ round(., 2)),
          p_str = if_else(p < 0.001, "<0.01", as.character(p)),
          lbl_bias = sprintf("Mean = %.2f", bias),
          lbl_hi   = sprintf("+1.96SD = %.2f", hi),
          lbl_lo   = sprintf("-1.96SD = %.2f", lo),
          lbl_cor_p = if_else(
            p_str == "<0.01",
            sprintf("r = %s; p %s", r, p_str),
            sprintf("r = %s; p = %s", r, p_str)
          )
        )
    }),
    
    # ---- Build plot ----
    plot = pmap(
      list(data, features, xlab, ylab, xlims, ylims, breaks),
      function(df, feat, xlab, ylab, xlims, ylims, nbreaks) {
        
        base_plot <- ggplot(df, aes(x = mu, y = diff)) +
          geom_point(alpha = 0.2, size = 0.75) +
          geom_hline(
            data = feat,
            aes(yintercept = bias),
            colour = "forestgreen",
            linewidth = 0.5
          ) +
          geom_hline(
            data = feat,
            aes(yintercept = hi),
            colour = "firebrick1",
            linetype = "dashed",
            linewidth = 0.5
          ) +
          geom_hline(
            data = feat,
            aes(yintercept = lo),
            colour = "firebrick1",
            linetype = "dashed",
            linewidth = 0.5
          ) +
          facet_wrap(~type, scales = "fixed", axes = "all") +
          scale_x_continuous(limits = xlims, n.breaks = nbreaks) +
          scale_y_continuous(limits = ylims, n.breaks = nbreaks) +
          theme_classic() +
          labs(x = xlab, y = ylab)
        
        base_plot +
          geom_text(
            data = feat,
            aes(x = Inf, y = bias, label = lbl_bias),
            colour = "forestgreen",
            hjust = 1.05, vjust = -0.6,
            size = 3, fontface = "bold",
            inherit.aes = FALSE
          ) +
          geom_text(
            data = feat,
            aes(x = Inf, y = hi, label = lbl_hi),
            colour = "firebrick1",
            hjust = 1.05, vjust = -0.6,
            size = 3, fontface = "bold",
            inherit.aes = FALSE
          ) +
          geom_text(
            data = feat,
            aes(x = Inf, y = lo, label = lbl_lo),
            colour = "firebrick1",
            hjust = 1.05, vjust = -0.6,
            size = 3, fontface = "bold",
            inherit.aes = FALSE
          ) +
          geom_text(
            data = feat,
            aes(x = Inf, y = Inf, label = lbl_cor_p),
            hjust = 1.05, vjust = 1.5,
            size = 3,
            inherit.aes = FALSE
          )
      }
    ),
    
    plot_path = paste0(
      "outputs/objs/01_agreement/plots/baltman_",
      name, ".png"
    )
  )

# Save
walk2(
  ba_nested$plot_path,
  ba_nested$plot,
  ~ ggsave(.x, .y, height = 4, width = 10, units = "in")
)

# Regression-based Bland-Altman ------------------------------------------------

# Repeat empirical structure
ba_nested <- scenario_grid %>%
  mutate(
    
    # ---- Filter + recode data ----
    data = map(name, ~ {
      
      a_imp_ba %>%
        filter(
          name == .x,
          type %in% facet_lookup$type_pattern
        ) %>%
        left_join(facet_lookup, by = c("type" = "type_pattern")) %>%
        mutate(
          type = factor(facet_label, levels = facet_levels)
        )
    }),
    
    # ---- Compute BA features ----
    features = map2(data, name, ~ {
      
      # Fit regression per facet
      reg <- .x %>%
        group_by(type) %>%
        group_modify(~{
          
          fit <- lm(diff ~ mu, data = .x)
          
          coef_summary <- summary(fit)$coefficients
          
          tibble(
            intercept = coef_summary["(Intercept)", "Estimate"],
            slope     = coef_summary["mu", "Estimate"],
            reg_p     = coef_summary["mu", "Pr(>|t|)"],  # <- key addition
            sigma     = summary(fit)$sigma,
            mu_min    = min(.x$mu, na.rm = TRUE),
            mu_max    = max(.x$mu, na.rm = TRUE)
          )
        })
      
      # Generate prediction lines
      reg_lines <- reg %>%
        rowwise() %>%
        mutate(
          mu_seq = list(seq(mu_min, mu_max, length.out = 100))
        ) %>%
        unnest(mu_seq) %>%
        mutate(
          bias = intercept + slope * mu_seq,
          hi   = bias + 1.96 * sigma,
          lo   = bias - 1.96 * sigma
        )
      
      # Prepare labels
      labels <- reg %>%
        mutate(
          mu_lab = mu_max,
          bias_lab = intercept + slope * mu_lab,
          hi_lab = bias_lab + 1.96 * sigma,
          lo_lab = bias_lab - 1.96 * sigma
        ) %>%
        left_join(
          b_ba_cor %>%
            filter(name == .y) %>%
            filter(type %in% facet_lookup$type_pattern) %>%
            left_join(facet_lookup, by = c("type" = "type_pattern")) %>%
            select(type = facet_label, r, p) %>% 
            mutate(
              type = factor(type, levels = facet_levels)
            ),
          by = "type"
        ) %>%
        mutate(
          p_str = if_else(p < 0.001, "<0.01", as.character(p)),
          lbl_bias = sprintf("Mean = %.2f", bias_lab),
          lbl_hi   = sprintf("+1.96SD = %.2f", hi_lab),
          lbl_lo   = sprintf("-1.96SD = %.2f", lo_lab),
          lbl_cor_p = if_else(
            p_str == "<0.01",
            sprintf("r = %s; p %s", r, p_str),
            sprintf("r = %s; p = %s", r, p_str)
          )
        )
      
      reg_lines %>%
        left_join(labels, by = "type") %>%
        mutate(
          lbl_bias = sprintf("Mean = %.2f", bias),
          lbl_hi   = sprintf("+1.96SD = %.2f", hi),
          lbl_lo   = sprintf("-1.96SD = %.2f", lo)
        )
      
      return(list(lines = reg_lines, labels = labels, reg = reg))
      
    }),
    
    # ---- Build plot ----
    plot = pmap(
      list(data, features, xlab, ylab, xlims, ylims, breaks),
      function(df, feat, xlab, ylab, xlims, ylims, nbreaks) {
        
        lines <- feat$lines
        labs  <- feat$labels
        
        base_plot <- ggplot(df, aes(mu, diff)) +
          geom_point(alpha = 0.2, size = 0.75) +
          geom_line(
            data = lines,
            aes(mu_seq, bias),
            colour = "forestgreen",
            linewidth = 0.5
          ) +
          geom_line(
            data = lines,
            aes(mu_seq, hi),
            colour = "firebrick1",
            linetype = "dashed",
            linewidth = 0.5
          ) +
          geom_line(
            data = lines,
            aes(mu_seq, lo),
            colour = "firebrick1",
            linetype = "dashed",
            linewidth = 0.5
          ) +
          facet_wrap(~type, scales = "fixed", axes = "all") +
          scale_x_continuous(limits = xlims, n.breaks = nbreaks) +
          scale_y_continuous(limits = ylims, n.breaks = nbreaks) +
          theme_classic() +
          labs(x = xlab, y = ylab)
        
        base_plot +
          geom_text(
            data = labs,
            aes(x = Inf, y = Inf, label = lbl_cor_p),
            hjust = 1.05, vjust = 1.5,
            size = 3,
            inherit.aes = FALSE
          )
      }
    ),
    
    plot_path = paste0(
      "outputs/objs/01_agreement/plots/baltman_regression_",
      name, ".png"
    )
  )

# Save
walk2(
  ba_nested$plot_path,
  ba_nested$plot,
  ~ ggsave(.x, .y, height = 4, width = 10, units = "in")
)

# Save regression outputs
ba_reg_sum <- ba_nested %>%
  select(name, features) %>%
  mutate(reg = map(features, "reg")) %>%
  select(-features) %>%
  unnest(reg)

write_csv(
  ba_reg_sum,
  "outputs/objs/01_agreement/summaries/baltman_regression.csv"
)
