
regex_site <- "site[0-9]+|s[0-9]+|SITE[0-9]+|Site[0-9]+|Site [0-9]+|S[0-9]+"

# Function to flag exacerbation timepoints
flag_exacerbation_timepoints <- function(df, type = "3d") {
  
  # define recovery window
  if (type == "3d") {
    recovery_gap <- 3
  } else if (type == "7d") {
    recovery_gap <- 7
  } else {
    stop("type must be '3d' or '7d'")
  }
  
  df <- df %>%
    arrange(date) %>%
    mutate(
      flagged_tmpts = NA_character_,
      pex_cycle = NA_integer_
    )
  
  # valid pex starts/ends
  pex_starts <- df %>%
    filter(tmpt == "pex_start", !is.na(date)) %>%
    pull(date)
  
  pex_ends <- df %>%
    filter(tmpt == "pex_end", !is.na(date)) %>%
    pull(date)
  
  # pair starts with the next valid end
  pex_pairs <- tibble(
    pex_start = pex_starts,
    pex_end = sapply(pex_starts, function(s) {
      ends_after <- pex_ends[pex_ends > s]
      if (length(ends_after) == 0) NA else min(ends_after)
    })
  ) %>%
    filter(!is.na(pex_end)) %>%
    mutate(pex_cycle = row_number())
  
  for (i in seq_len(nrow(pex_pairs))) {
    
    s <- pex_pairs$pex_start[i]
    e <- pex_pairs$pex_end[i]
    recovery_date <- e + recovery_gap
    cycle_id <- pex_pairs$pex_cycle[i]
    
    next_pex_start <- if (i < nrow(pex_pairs)) {
      pex_pairs$pex_start[i + 1]
    } else {
      as.Date(Inf)
    }
    
    # a) last stable before pex_start
    last_stable_date <- df %>%
      filter(tmpt == "stable", date < s) %>%
      slice_tail(n = 1) %>%
      pull(date)
    
    if (length(last_stable_date) == 1) {
      df <- df %>%
        mutate(
          flagged_tmpts = if_else(
            date == last_stable_date,
            "last_stable_pre",
            flagged_tmpts
          ),
          pex_cycle = if_else(
            date == last_stable_date,
            cycle_id,
            pex_cycle
          )
        )
    }
    
    # b) pex window INCLUDING recovery gap
    df <- df %>%
      mutate(
        flagged_tmpts = if_else(
          date >= s & date < recovery_date,
          "pex",
          flagged_tmpts
        ),
        pex_cycle = if_else(
          date >= s & date < recovery_date,
          cycle_id,
          pex_cycle
        )
      )
    
    # c) FIRST stable_again only
    first_stable_again_date <- df %>%
      filter(
        tmpt == "stable",
        date >= recovery_date
      ) %>%
      slice_head(n = 1) %>%
      pull(date)
    
    if (length(first_stable_again_date) == 1) {
      
      df <- df %>%
        mutate(
          flagged_tmpts = if_else(
            date >= first_stable_again_date &
              date < next_pex_start,
            "stable_again",
            flagged_tmpts
          ),
          pex_cycle = if_else(
            date >= first_stable_again_date &
              date < next_pex_start,
            cycle_id,
            pex_cycle
          )
        )
    }
  }
  
  df
}