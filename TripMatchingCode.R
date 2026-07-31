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
  
  if (length(idx) == 0) return(0)
  
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
  
  tp <- sum(m_final == 1)
  fp <- sum(m_final == 0)
  fn <- tm - tp
  
  if (tp == 0) return(0)
  
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)
  
  return(2 * (precision * recall) / (precision + recall))
}

# create parameter grid to optimize over
threshold_grid <- expand.grid(
  t_anglers = seq(0, 1, by = 0.2),
  t_time    = seq(0, 1, by = 0.05),
  t_hours   = seq(0, 1, by = 0.05)
)

# calculate f1 scores
message("Optimizing thresholds...")

results <- threshold_grid %>%
  mutate(f1_score = pmap_dbl(list(t_anglers, t_time, t_hours), 
                             ~calc_f1(..1, ..2, ..3)))

opt <- results %>% arrange(desc(f1_score), desc(t_anglers), desc(t_hours), desc(t_time)) %>% slice(1)

# Apply optimal thresholds to obtain the matched set
full_matches <- eval_df %>%
  filter(Anglers_Sim >= opt$t_anglers, Time_Sim >= opt$t_time, Hours_Sim >= opt$t_hours) %>%
  # apply same tiebreaking logic as used in the optimization step
  arrange(Surv_Survey_RowID, desc(Time_Sim), desc(Anglers_Sim), desc(Hours_Sim)) %>% slice_head(n = 1, by = Surv_Survey_RowID)

# 5. Visualize Optimization ----------------------------------------------------
plot_data <- results %>% 
  filter(round(t_anglers / 0.2, 5) %% 1 == 0) %>%
  rename("Angler Threshold" = t_anglers)

ggplot(plot_data, aes(x = t_time, y = t_hours, fill = f1_score)) +
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

# 6. Final Outputs and Diagnostics -------------------------------------------------------------
cat("\n--- Optimal Thresholds ---\n")
print(opt)

# Calculate error rates
tp <- sum(full_matches$is_match == 1)
fp <- sum(full_matches$is_match == 0)

fmr  <- fp / nrow(full_matches)  # False Match Rate (1 - Precision)
fnmr <- (tm - tp) / tm           # False Non-Match Rate (1 - Recall)

cat(sprintf("\nMatches: %d | FMR: %.2f%% | FNMR: %.2f%%\n", nrow(full_matches), fmr * 100, fnmr * 100))
