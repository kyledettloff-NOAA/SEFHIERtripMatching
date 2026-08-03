# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script: SEFHIER Trip Matching & Threshold Optimization
# Repository: SEFHIERtripMatching
# Authors: Michelle Masi, Kyle Dettloff
# Dependencies: tidyverse
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

# 1. Data Acquisition ----------------------------------------------------------
base_url <- "https://github.com/kyledettloff-NOAA/SEFHIERtripMatching/raw/main/"

# Read in data directly from GitHub
logbooks_list <- readRDS(url(paste0(base_url, "fake_logbookdata.rds")))
surveys_list  <- readRDS(url(paste0(base_url, "fake_surveydata.rds")))

# combine lists into dataframes and assign unique row numbers
log_df  <- bind_rows(logbooks_list) %>% rename_with(~paste0("Log_", .x)) %>%
  mutate(Log_Logbook_RowID = row_number(),
         # Convert time to minutes since midnight for proper subtraction
         Log_Mins  = as.numeric(Log_TIME) %/% 100 * 60 + as.numeric(Log_TIME) %% 100)
surv_df <- bind_rows(surveys_list)  %>% rename_with(~paste0("Surv_", .x)) %>%
  mutate(Surv_Survey_RowID = row_number(),
         # Convert time to minutes since midnight for proper subtraction
         Surv_Mins = as.numeric(Surv_TIME) %/% 100 * 60 + as.numeric(Surv_TIME) %% 100)

# 2. Match Logbook and Survey Data by Date -------------------------------------
message("Joining datasets and calculating similarity scores...")

matched_pool <- inner_join(
  log_df, surv_df, 
  by = c("Log_Full_Date" = "Surv_Full_Date"),
  relationship = "many-to-many"
) %>%
  mutate(
    # Exponential Similarity
    Anglers_Sim = exp(log(0.8) * abs(as.numeric(Log_Num_Anglers) - as.numeric(Surv_Num_Anglers))),
    Hours_Sim =   exp(log(0.8) * abs(as.numeric(Log_Hours_Fished) - as.numeric(Surv_Hours_Fished))),
    
    # Continuous Similarity: 1 / (1 + abs(diff))
    Time_Sim    = 1 / (1 + abs(Log_Mins - Surv_Mins) / 60),
    
    # Binary Exact-Match Features (1 = Match, 0 = Disagreement)
    Site_Sim    = as.numeric(as.character(Log_State) == as.character(Surv_State) & as.character(Log_County) == as.character(Surv_County))
  ) %>%
  # assign 0 similarity if data is missing (NA)
  mutate(across(ends_with("_Sim"), ~replace_na(.x, 0)))

# 3. Ground Truth & Evaluation Set ---------------------------------------------
# identify true matches
true_matches <- matched_pool %>%
  # filter to matching vessel numbers
  filter(Log_Vessel_Official_Num == Surv_Vessel_Official_Num) %>%
  # keep record with highest time similarity after prioritizing same site when multiple matches for same vessel on same date
  arrange(Surv_Survey_RowID, desc(Site_Sim), desc(Time_Sim), desc(Anglers_Sim), desc(Hours_Sim)) %>%
  slice_head(n = 1, by = Surv_Survey_RowID) %>%
  mutate(is_match = 1)

eval_df <- matched_pool %>%
  # require site to match for optimization of other thresholds
  filter(Site_Sim == 1) %>%
  left_join(select(true_matches, Log_Logbook_RowID, Surv_Survey_RowID, is_match), 
            by = c("Log_Logbook_RowID", "Surv_Survey_RowID")) %>%
  mutate(is_match = replace_na(is_match, 0))

# 4. Grid Search Optimization --------------------------------------------------
# Pre-extract atomic vectors to eliminate data frame creation overhead during the optimization
surv_ids <- eval_df$Surv_Survey_RowID
a_sim    <- eval_df$Anglers_Sim
t_sim    <- eval_df$Time_Sim
h_sim    <- eval_df$Hours_Sim
is_m     <- eval_df$is_match
tm       <- nrow(true_matches)

calc_f1 <- function(ang, tim, hrs) {
  # filter candidates by thresholds
  idx <- which(a_sim >= ang & t_sim >= tim & h_sim >= hrs)
  
  if (length(idx) == 0) {
    return(c(f1_score = 0, n_matches = 0))
  }
  
  surv_sub <- surv_ids[idx]
  t_sub    <- t_sim[idx]
  a_sub    <- a_sim[idx]
  h_sub    <- h_sim[idx]
  m_sub    <- is_m[idx]
  
  # de-duplicate surviving candidates
  ord      <- order(surv_sub, -t_sub, -a_sub, -h_sub, na.last = TRUE)
  surv_ord <- surv_sub[ord]
  keep     <- !duplicated(surv_ord)
  m_final  <- m_sub[ord][keep]
  
  # keep number of matches for plotting
  total_matches <- length(m_final)
  
  tp <- sum(m_final == 1)
  fp <- sum(m_final == 0)
  fn <- tm - tp
  
  if (tp == 0) {
    return(c(f1_score = 0, n_matches = total_matches))
  }
  
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)
  f1        <- 2 * (precision * recall) / (precision + recall)
  
  return(c(f1_score = f1, n_matches = total_matches))
}

# create parameter grid to optimize over
threshold_grid <- expand.grid(
  t_anglers = seq(0, 1, by = 0.2),
  t_time    = seq(0, 1, by = 0.01),
  t_hours   = seq(0, 1, by = 0.01)
)

# calculate f1 scores
message("Optimizing thresholds...")

results <- threshold_grid %>%
  bind_cols(
    pmap_dfr(list(threshold_grid$t_anglers, threshold_grid$t_time, threshold_grid$t_hours), 
             ~calc_f1(..1, ..2, ..3))
  )

opt <- results %>% arrange(desc(f1_score), desc(t_anglers), desc(t_hours), desc(t_time)) %>% slice(1)

# Apply optimal thresholds to obtain the matched set
full_matches <- eval_df %>%
  filter(Anglers_Sim >= opt$t_anglers, Time_Sim >= opt$t_time, Hours_Sim >= opt$t_hours) %>%
  # apply same tiebreaking logic as used in the optimization step
  arrange(Surv_Survey_RowID, desc(Time_Sim), desc(Anglers_Sim), desc(Hours_Sim)) %>%
  slice_head(n = 1, by = Surv_Survey_RowID)

# 5. Visualize Optimization ----------------------------------------------------
plot_data <- results %>% 
  filter(round(t_anglers / 0.2, 5) %% 1 == 0) %>%
  rename("Angler Threshold" = t_anglers)

# F1 score
p1 <- ggplot(plot_data, aes(x = t_time, y = t_hours, fill = f1_score)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
    # center contrast around 40th %ile score
    midpoint = quantile(plot_data$f1_score, probs = 0.4, na.rm = TRUE),
    name = "F1 Score"
  ) +
  facet_wrap(~`Angler Threshold`, labeller = label_both, ncol = 1) +
  geom_point(data = opt %>% rename("Angler Threshold" = t_anglers),
             aes(x = t_time, y = t_hours), 
             color = "black", shape = 8, size = 3.5, stroke = 1.5) +
  theme_bw(base_size = 12) + 
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    legend.key.width = unit(2, "cm")
  ) +
  labs(
    x = "Time Similarity Threshold",
    y = "Hours Fished Similarity Threshold",
    caption = paste0("Global Optimal F1: ", round(opt$f1_score, 3))
  )

# find "core" F1 score cutoff within optimal angler threshold
optimal_f1_cutoff <- results %>%
  filter(t_anglers == opt$t_anglers) %>%
  arrange(f1_score) %>%
  mutate(
    rank_f1 = row_number(),
    score   = (f1_score - min(f1_score)) / (max(f1_score) - min(f1_score)) - 
      (rank_f1 - 1) / (n() - 1)
  ) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  pull(f1_score)

# filter top F1 results restricted to optimal angler threshold
top_f1_data <- plot_data %>%
  filter(`Angler Threshold` == opt$t_anglers, f1_score >= optimal_f1_cutoff) %>%
  mutate(match_ratio = n_matches / tm,
         percent_bias = (1 / match_ratio - 1) * 100)

# number of matches relative to actual within core F1 range
p2 <- ggplot() +
  # Base layer: draw all tiles in light grey for the cut-out combinations
  geom_tile(
    data = plot_data %>% filter(`Angler Threshold` == opt$t_anglers),
    aes(x = t_time, y = t_hours),
    fill = "grey90"
  ) +
  # Top layer: draw tiles for the retained top F1 combinations
  geom_tile(
    data = top_f1_data,
    aes(x = t_time, y = t_hours, fill = match_ratio)
  ) +
  scale_fill_gradient2(
    low = "#d7191c",
    mid = "white",
    high = "#2c7bb6",
    midpoint = 1,
    limits = c(min(0.5, min(top_f1_data$match_ratio, na.rm = TRUE)),
               max(1.5, max(top_f1_data$match_ratio, na.rm = TRUE))),
    breaks = c(0.5, 1, 1.5),
    oob = scales::squish,
    name = "True Match Ratio"
  ) +
  facet_wrap(~`Angler Threshold`, labeller = label_both, ncol = 1) +
  geom_point(
    data = opt %>% rename("Angler Threshold" = t_anglers),
    aes(x = t_time, y = t_hours), 
    color = "black", shape = 8, size = 3.5, stroke = 1.5
  ) +
  theme_bw(base_size = 13) + 
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    legend.key.width = unit(2, "cm")
  ) +
  labs(
    x = "Time Similarity Threshold",
    y = "Hours Fished Similarity Threshold",
    caption = paste0("Matches at Optimal F1 Threshold: ", opt$n_matches, " (Ratio: ", round(opt$n_matches / tm, 2), ")")
  )

# distribution of bias introduced
p3 <- ggplot(top_f1_data, aes(x = percent_bias)) +
  geom_histogram(
    aes(
      y = after_stat(count / sum(count)), 
      fill = after_stat(1 / (x / 100 + 1))
    ), 
    binwidth = 5, 
    boundary = 0,
    color = "gray30", 
    alpha = 0.9
  ) +
  # Secondary reference lines at +/- 5%
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "gray40", linewidth = 0.5) +
  # Reference line at 0% bias
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#d7191c",
    mid = "white",
    high = "#2c7bb6",
    midpoint = 1,
    limits = c(
      min(0.5, min(top_f1_data$match_ratio, na.rm = TRUE)),
      max(1.5, max(top_f1_data$match_ratio, na.rm = TRUE))
    ),
    breaks = c(0.5, 1, 1.5),
    oob = scales::squish,
    name = "True Match Ratio"
  ) +
  scale_x_continuous(
    n.breaks = 10,
    labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%")
  ) +
  scale_y_continuous(
    n.breaks = 8,
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold"),
    legend.position = "top",
    legend.key.width = unit(2, "cm")
  ) +
  labs(
    x = "Bias in Estimates of Total",
    y = "Proportion of Threshold Combinations",
    caption = paste0(
      "Combinations at F1 >= ", round(optimal_f1_cutoff, 3), "\n",
      round(mean(abs(top_f1_data$percent_bias) <= 5, na.rm = TRUE) * 100, 1),
      "% of combinations within ±5% bias"
    )
  )

# 6. Final Outputs and Diagnostics -------------------------------------------------------------
cat("\n--- Optimal Thresholds ---\n")
print(opt)

print(p1)

# Calculate error rates
tp <- sum(full_matches$is_match == 1)
fp <- sum(full_matches$is_match == 0)

fmr  <- fp / nrow(full_matches)  # False Match Rate (1 - Precision)
fnmr <- (tm - tp) / tm           # False Non-Match Rate (1 - Recall)

cat(sprintf("\nMatches: %d | FMR: %.2f%% | FNMR: %.2f%%\n", nrow(full_matches), fmr * 100, fnmr * 100))

print(p2)
print(p3)
