# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script: SEFHIER Trip Matching & Threshold Optimization
# Repository: SEFHIERtripMatching
# Authors: Michelle Masi, Kyle Dettloff
# Dependencies: tidyverse, stringdist
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(stringdist)

# 1. Data Acquisition ----------------------------------------------------------
base_url <- "https://github.com/kyledettloff-NOAA/SEFHIERtripMatching/raw/main/"

# Read in data directly from GitHub
logbooks_list <- readRDS(url(paste0(base_url, "fake_logbookdata.rds")))
surveys_list  <- readRDS(url(paste0(base_url, "fake_surveydata.rds")))

# combine lists into dataframes and assign unique row numbers
log_df  <- bind_rows(logbooks_list) %>% rename_with(~paste0("Log_", .x)) %>% mutate(Log_Logbook_RowID = row_number())
surv_df <- bind_rows(surveys_list)  %>% rename_with(~paste0("Surv_", .x)) %>% mutate(Surv_Survey_RowID = row_number())

# 2. Match Logbook and Survey Data by Date -------------------------------------
message("Joining datasets and calculating similarity scores...")

matched_pool <- inner_join(
  log_df, surv_df, 
  by = c("Log_Full_Date" = "Surv_Full_Date"),
  relationship = "many-to-many"
) %>%
  mutate(
    # Character Similarity (Jaro-Winkler)
    VslNum_Sim    = stringsim(as.character(Log_Vessel_Official_Num), as.character(Surv_Vessel_Official_Num), method = "jw"),
    VslName_Sim   = stringsim(as.character(Log_Vessel_Name), as.character(Surv_Vessel_Name), method = "jw"),
    
    # Continuous Similarity: 1 / (1 + abs(diff))
    # assign 0 similarity if data is missing (NA)
    Anglers_Sim = 1 / (1 + abs(as.numeric(Log_Num_Anglers) - as.numeric(Surv_Num_Anglers))),
    Hours_Sim   = 1 / (1 + abs(as.numeric(Log_Hours_Fished) - as.numeric(Surv_Hours_Fished))),
    Time_Sim    = 1 / (1 + abs(as.numeric(Log_TIME) - as.numeric(Surv_TIME))),
    
    # Binary Exact-Match Features (1 = Match, 0 = Disagreement)
    Caught_Sim  = as.numeric(as.character(Log_Anything_Caught_Flag) == as.character(Surv_Anything_Caught_Flag)),
    Site_Sim    = as.numeric(as.character(Log_State) == as.character(Surv_State) & as.character(Log_County) == as.character(Surv_County))
  ) %>%
  mutate(across(ends_with("_Sim"), ~replace_na(.x, 0)))

# 3. Ground Truth & Evaluation Set ---------------------------------------------
# set threshold
ves_sim_thres <- 1
# identify true matches
true_matches <- matched_pool %>%
  # filter to high vessel similarities and exclude vessel name matches when "UNNAMED"
  filter(VslNum_Sim >= ves_sim_thres | VslName_Sim >= ves_sim_thres & Log_Vessel_Name != "UNNAMED" & Surv_Vessel_Name != "UNNAMED") %>%
  # keep record with highest time similarity after prioritizing same site when multiple matches for same vessel on same date
  group_by(Surv_Survey_RowID) %>% arrange(desc(Site_Sim), desc(Time_Sim)) %>% slice(1) %>% %>% ungroup() %>%
  mutate(is_match = 1)

eval_df <- matched_pool %>%
  left_join(select(true_matches, Log_Logbook_RowID, Surv_Survey_RowID, is_match), 
            by = c("Log_Logbook_RowID", "Surv_Survey_RowID")) %>%
  mutate(is_match = replace_na(is_match, 0))

# 4. Grid Search Optimization --------------------------------------------------
is_m       <- eval_df$is_match == 1
a_sim      <- eval_df$Anglers_Sim
t_sim      <- eval_df$Time_Sim
h_sim      <- eval_df$Hours_Sim
site_sim   <- eval_df$Site_Sim
caught_sim <- eval_df$Caught_Sim

calc_f1 <- function(ang, tim, hrs) {
  pred <- (site_sim == 1) & (caught_sim == 1) & (a_sim >= ang) & (t_sim >= tim) & (h_sim >= hrs)
  pred[is.na(pred)] <- FALSE 
  
  tp <- sum(pred & is_m)
  fp <- sum(pred & !is_m)
  fn <- sum(!pred & is_m)
  
  if (tp == 0) return(0)
  
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)
  
  return(2 * (precision * recall) / (precision + recall))
}

# create parameter grid to optimize over
threshold_grid <- expand.grid(
  t_anglers = seq(0, 1, by = 0.2),
  t_time    = seq(min(true_matches$Time_Sim, na.rm = TRUE), 
                  max(true_matches$Time_Sim, na.rm = TRUE), by = 0.05),
  t_hours   = seq(min(true_matches$Hours_Sim, na.rm = TRUE), 
                  max(true_matches$Hours_Sim, na.rm = TRUE), by = 0.05)
)

# calculate f1 scores
message("Optimizing thresholds...")
results <- threshold_grid %>%
  mutate(f1_score = pmap_dbl(list(t_anglers, t_time, t_hours), 
                             ~calc_f1(..1, ..2, ..3)))

opt <- results %>% slice_max(f1_score, n = 1, with_ties = FALSE)

# 5. Visualize Optimization ----------------------------------------------------
plot_data <- results %>% 
  filter(t_anglers <= max(opt$t_anglers, 0.5)) %>% 
  rename("Angler Threshold" = t_anglers)

ggplot(plot_data, aes(x = t_time, y = t_hours, fill = f1_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "F1 Score") +
  facet_wrap(~`Angler Threshold`, labeller = label_both, ncol = 1) +
  geom_point(data = opt %>% rename("Angler Threshold" = t_anglers),
             aes(x = t_time, y = t_hours), 
             color = "cyan", shape = 8, size = 3, stroke = 1.5) +
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

# 6. Final Outputs -------------------------------------------------------------
cat("\n--- Final Model Results ---\n")
print(opt)
