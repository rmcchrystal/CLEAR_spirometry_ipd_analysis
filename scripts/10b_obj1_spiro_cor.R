
# Objective 1
# Generate spirometry correlation plots

source("Scripts/00_config.R")
source("Scripts/00_functions_regex.R")

library(tidyverse)

# Import directory
z_filedir <- "outputs/objs/01_agreement/summaries"

# Import materials
a_imp_df <- readRDS("scratch_data/bland_altman_data.rds")
a_imp_cor <- read_csv(paste0(z_filedir, "/spiro_correlation.csv"))

# Define labels
facet_levels <- c(
  "± 7 days",
  "Same day during stability",
  "Same day during PEx"
)

facet_lookup <- tribble(
  ~type_pattern,      ~facet_label,
  "same_week_all",    "± 7 days",
  "same_day_stable",  "Same day during stability",
  "same_day_pex",     "Same day during PEx"
)

# Define measures
measure_specs <- tribble(
  ~name,          ~xlab,                                   ~ylab,                                   ~limits,        ~breaks, ~label_x,
  "fev1",         "Home spirometry FEV1 (L/s)",            "Clinic spirometry FEV1 (L/s)",          c(0, 6),        4,       0.9,
  "fvc",          "Home spirometry FVC (L)",               "Clinic spirometry FVC (L)",             c(0, 8),        5,       1.1,
  "fev1pp",       "Home spirometry FEV1 % predicted",      "Clinic spirometry FEV1 % predicted",    c(0, 200),      5,       30,
  "fvcpp",        "Home spirometry FVC % predicted",       "Clinic spirometry FVC % predicted",     c(0, 200),      5,       30,
  "fev1_fvc",     "Home spirometry FEV1/FVC",              "Clinic spirometry FEV1/FVC",            c(25, 100),     4,       35
)

# Prepare correlation tibble
b_cor <- a_imp_cor %>%
  mutate(across(cor_pearsons:cor_p, ~ round(., 2))) %>%
  rename(r = cor_pearsons, p = cor_p) %>% 
  mutate(
    type = case_when(
      type == "Irrespective of disease status" & tmpt == "same_week" ~ "same_week_all",
      type == "During stable periods"          & tmpt == "same_day"  ~ "same_day_stable",
      type == "Between PEx start and end"      & tmpt == "same_day"  ~ "same_day_pex",
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(type))

# Create a scenario grid
scenario_grid <- measure_specs

# Nest data and build plots
plots_nested <- scenario_grid %>%
  mutate(
    
    # ---- Filter + recode data ----
    data = map(name, ~ {
      
      a_imp_df %>%
        filter(
          name == .x,
          type %in% facet_lookup$type_pattern
        ) %>%
        left_join(facet_lookup, by = c("type" = "type_pattern")) %>%
        mutate(
          type = factor(facet_label, levels = facet_levels)
        )
    }),
    
    # ---- Merge correlation features (properly grouped per facet) ----
    features = map2(data, name, ~ {
      
      cor_sub <- b_cor %>%
        filter(
          name == .y,
          type %in% facet_lookup$type_pattern
        ) %>%
        left_join(facet_lookup, by = c("type" = "type_pattern")) %>%
        mutate(
          type = factor(facet_label, levels = facet_levels),
          p_str = if_else(p < 0.001, "<0.01", as.character(p)),
          lbl_cor_p = if_else(
            p_str == "<0.01",
            sprintf("r = %s; p %s", r, p_str),
            sprintf("r = %s; p = %s", r, p_str)
          )
        ) %>%
        distinct(type, lbl_cor_p)
      
      cor_sub
    }),
    
    # ---- Build plot ----
    plot = pmap(
      list(data, features, xlab, ylab, limits, breaks, label_x),
      function(df, feat, xlab, ylab, lims, nbreaks, label_x) {
        
        base_plot <- ggplot(df, aes(x = home, y = clinic)) +
          geom_point(alpha = 0.2, size = 0.75) +
          geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed",
            colour = "red"
          ) +
          facet_wrap(~type, scales = "fixed", axes = "all") +
          scale_y_continuous(limits = lims, n.breaks = nbreaks) +
          scale_x_continuous(limits = lims, n.breaks = nbreaks) +
          theme_classic() +
          labs(x = xlab, y = ylab)
        
        base_plot +
          geom_text(
            data = feat,
            aes(x = label_x, y = Inf, label = lbl_cor_p),
            vjust = 1.5,
            size = 3,
            inherit.aes = FALSE
          )
      }
    ),
    
    plot_path = paste0(
      "outputs/objs/01_agreement/plots/spirocor_",
      name, ".png"
    )
  )

# Save plots in one pass
walk2(
  plots_nested$plot_path,
  plots_nested$plot,
  ~ ggsave(.x, .y, height = 4, width = 10, units = "in")
)
