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

### Diagnostic plots
# Candidate Pool Sizes Plot (Fig. 1) -----------------------------------
plot_fig1 <- eval_df %>%
  group_by(Surv_Survey_RowID) %>%
  summarize(
    Candidate_Count = n(),
    Is_True_Match_Present = if_else(any(is_match == 1), "True Match Found", "No Match in Pool")
  )

p1 <- ggplot(plot_fig1, aes(x = Candidate_Count, fill = Is_True_Match_Present)) +
  geom_histogram(binwidth = 5, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("True Match Found" = "black", "No Match in Pool" = "gray75"),
    name = NULL 
  ) +
  labs(
    x = "Candidate Pool Size (Logbooks per Survey)",
    y = "Number of Surveys"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "inside",
    legend.position.inside = c(0.70, 0.80),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    legend.key.size = unit(1.2, "lines"),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Signal-to-Noise Overlap Plot (Fig. 2) ---------------------------------
plot_fig2 <- eval_df %>%
  select(is_match, Anglers = Anglers_Sim, Hours = Hours_Sim, `Trip Time` = Time_Sim) %>%
  pivot_longer(
    cols = c(Anglers, Hours, `Trip Time`),
    names_to = "Variable",
    values_to = "SimilarityScore"
  ) %>%
  mutate(
    Variable = case_when(
      Variable == "Anglers" ~ "Number of Anglers",
      Variable == "Hours"   ~ "Hours Fished",
      .default = Variable
    ),
    Variable = factor(Variable, levels = c("Number of Anglers", "Hours Fished", "Trip Time")),
    Category = if_else(is_match == 1, "True Match", "Potential Mismatch")
  )

p2 <- ggplot(plot_fig2, aes(x = SimilarityScore, fill = Category)) +
  geom_histogram(
    aes(y = after_stat(count / ave(count, PANEL, group, FUN = sum))),
    position = position_dodge(width = 0.05),
    binwidth = 0.05,
    boundary = 0,
    closed = "left",
    alpha = 1.0,
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~Variable, scales = "fixed", ncol = 1, axes = "all_x") +
  scale_y_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1, by = 0.20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = c("True Match" = "black", "Potential Mismatch" = "gray65")) +
  labs(x = "Similarity Score", y = "Proportion", fill = "") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 14, margin = margin(b = 10)),
    strip.background = element_blank(),
    axis.text.x = element_text(color = "black", size = 14, face = "plain"),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.key.size = unit(1.2, "lines")
  )

print(p1)
print(p2)

# 4. Grid Search Optimization --------------------------------------------------
# Vectorize data inputs for faster optimization
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
opt_matches <- eval_df %>%
  filter(Anglers_Sim >= opt$t_anglers, Time_Sim >= opt$t_time, Hours_Sim >= opt$t_hours) %>%
  # apply same tiebreaking logic as used in the optimization step
  arrange(Surv_Survey_RowID, desc(Time_Sim), desc(Anglers_Sim), desc(Hours_Sim)) %>%
  slice_head(n = 1, by = Surv_Survey_RowID)

# 5. Visualize Optimization ----------------------------------------------------
plot_data <- results %>% 
  filter(round(t_anglers / 0.2, 5) %% 1 == 0) %>%
  rename("Angler Threshold" = t_anglers)

# F1 score
p3 <- ggplot(plot_data, aes(x = t_time, y = t_hours, fill = f1_score)) +
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
p4 <- ggplot() +
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
    name = "True Match Ratio"
  ) +
  facet_wrap(~`Angler Threshold`, labeller = label_both, ncol = 1) +
  geom_point(
    data = opt %>% rename("Angler Threshold" = t_anglers),
    aes(x = t_time, y = t_hours), 
    color = "black", shape = 8, size = 3.5, stroke = 1.5
  ) +
  theme_bw(base_size = 14) + 
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
p5 <- ggplot(top_f1_data, aes(x = percent_bias)) +
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
  theme_bw(base_size = 14) +
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

print(p3)

# Calculate error rates
tp <- sum(opt_matches$is_match == 1)
fp <- sum(opt_matches$is_match == 0)

fmr  <- fp / nrow(opt_matches)  # False Match Rate (1 - Precision)
fnmr <- (tm - tp) / tm           # False Non-Match Rate (1 - Recall)

cat(sprintf("\nMatches: %d | FMR: %.2f%% | FNMR: %.2f%%\n", nrow(opt_matches), fmr * 100, fnmr * 100))

print(p4)
print(p5)
