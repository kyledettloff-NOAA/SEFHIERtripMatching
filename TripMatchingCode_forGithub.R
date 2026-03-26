#This code matches Survey Data to Logbook Data with and without a unique ID, using Jaro-Winkler Similarity, then uses dplyr piping to filter on common variables when the unique ID is withheld


#Libraries needed ----
library(tidyverse)
library(stringdist)
library(yardstick) #for comparing similarity score distributions to ID threshold values


#Get Data ----
#delete this after data is on Git repo ----
Michelles_path <- "C:/Users/michelle.masi/Documents/SEFHIER/R code/Validation Survey data and analyses/Matching Logbooks to Validation Survey Intercepts/include all permit types/"
# !! Change to your path !!
Path <- Michelles_path
#---------------------------

##KEEP THIS AFTER UPLOADING TO GIT ----
#read in RDS list files (fake data) from git
url <- "https://github.com/kyledettloff-NOAA/SEFHIERtripMatching"
#get Survey list
download.file(url, "fake_surveydata.rds", mode = "wb")
Surveys.list <- readRDS("fake_surveydata.rds")
#get Logbook list
download.file(url, "fake_logbookdata.rds", mode = "wb")
Logbooks.list <- readRDS("fake_logbookdata.rds")



#Matching Logbook to Survey (run this for the (monthly) survey and logbook lists ----
##IMPORTANT NOTE: Data parsed into monthly chunks for loop efficiency (otherwise it'd take the for-loop weeks to loop through 49k logbooks x 1835 surveys ----
##Loop through logbooks list for each survey (by month) ----

#create a temp list to store data from for-loop (note: use temp list vs DF for efficiency in looping through large amounts of data)
results_list <- list()

#time for loop run time
start_time <- Sys.time()

# Nested for loop~~~~~~~~~~
#iterate along each element in list i
for (i in seq_along(Logbooks.list)) { 
  
    #iterate along each element in list j
    for (j in seq_along(Surveys.list)) { 
      
        #iterate along each row in element i of list i
        for (row_i_index in 1:nrow(Logbooks.list[[i]])) { 
          
          #iterate along each row in element j of list j
          for (row_j_index in 1:nrow(Surveys.list[[j]])) { 
            
            ##only iterate through j for each row in i IF:
            #the list element numbers match (since element numbers are in order of months - Jan through Dec, in both lists) AND
            #Only calc similarity scores and save the results IF the Dates (YEAR-MON-DAY) are equal between row index i of Logbooks DF and row index j of the Survey DF
            if (i == j && Logbooks.list[[i]]$Full_Date[[row_i_index]] == Surveys.list[[i]]$Full_Date[[row_j_index]]) {
            
              #print out which logbook DF and row iteration we're on, during the loop
              Logbook_ListElement_ofi <- Logbooks.list[[i]]
              Survey_ListElement_ofj <- Surveys.list[[j]]
              print(paste0("Logbook DF Number: ", i, 
                           " ; Logbook Row Number: ", Logbooks.list[[i]]$Logbook_RowID[[row_i_index]],
                           " ; Survey DF Number: ", j,
                           " ; Survey Row Number: ", Surveys.list[[i]]$Survey_RowID[[row_j_index]])) 
                  
            
              #Apply a similarity formula for each pair of logbook vs survey rows, along the loop~~~~~~~~~
              #calc difference between dates
                diff_dates <- difftime(Logbooks.list[[i]]$Full_Date[[row_i_index]],
                                       Surveys.list[[j]]$Full_Date[[row_j_index]], units = "days")
                Dateresult <- as.numeric(diff_dates)
              
              #jw distance scores on characters
                Countyresult <- stringsim(as.character(Logbooks.list[[i]]$County[[row_i_index]]), 
                                          as.character(Surveys.list[[j]]$County[[row_j_index]]), 
                                          method = "jw")
                Stateresult <- stringsim(as.character(Logbooks.list[[i]]$State[[row_i_index]]), 
                                         as.character(Surveys.list[[j]]$State[[row_j_index]]), 
                                         method = "jw")
                VslNumberresult <- stringsim(as.character(Logbooks.list[[i]]$Vessel_Official_Num[[row_i_index]]), 
                                             as.character(Surveys.list[[j]]$Vessel_Official_Num[[row_j_index]]), 
                                             method = "jw")
                VslNameresult <- stringsim(as.character(Logbooks.list[[i]]$Vessel_Name[[row_i_index]]), 
                                           as.character(Surveys.list[[j]]$Vessel_Name[[row_j_index]]),
                                           method = "jw")
              
             
                 
                # 2. Continuous Variables (Calc the inverse distance similarity)
                # Formula: 1 / (1 + abs(val1 - val2))
                
                #ignores the ratio and looks only at the raw difference, where a larger number compared to 0 will result in the lower similiarity than a small number compared to 0
                #The formula \(1/(1+|x-y|)\) : mathematically safe because: the absolute difference \(|x-y|\) can never be less than zero.
                #The constant 1 ensures that even if \(x\) and \(y\) are identical (making \(|x-y|=0\)), the denominator becomes \(1+0=1\).
                #The denominator will always be greater than or equal to 1, making it impossible for it to ever reach zero
                
                
                # NumAnglers
                NumAnglersresult <- 1 / (1 + abs(as.numeric(Logbooks.list[[i]]$Num_Anglers[[row_i_index]]) - 
                                                   as.numeric(Surveys.list[[j]]$Num_Anglers[[row_j_index]])))
                
                # AnythingCaught
                AnythingCaughtresult <- 1 / (1 + abs(as.numeric(Logbooks.list[[i]]$Anything_Caught_Flag[[row_i_index]]) - 
                                                       as.numeric(Surveys.list[[j]]$Anything_Caught_Flag[[row_j_index]])))
                
                # HoursFished
                HoursFishedresult <- 1 / (1 + abs(as.numeric(Logbooks.list[[i]]$Hours_Fished[[row_i_index]]) - 
                                                    as.numeric(Surveys.list[[j]]$Hours_Fished[[row_j_index]])))
                
                # Time (Inverse distance similarity, no buffer)
                Timeresult <- 1 / (1 + abs(as.numeric(Logbooks.list[[i]]$TIME[[row_i_index]]) - 
                                             as.numeric(Surveys.list[[j]]$TIME[[row_j_index]])))    
        
              # Create a temporary data frame for this iteration's results~~~~~~~~~~~~~
              iteration_df <- data.frame(Logbooks.list[[i]]$Logbook_RowID[[row_i_index]], #ID the logbook row #
                                         Surveys.list[[j]]$Survey_RowID[[row_j_index]],   #ID the survey row #
                                         
                                         Logbooks.list[[i]]$Full_Date[[row_i_index]],  #add date to confirm dates match
                                         Surveys.list[[j]]$Full_Date[[row_j_index]],
                                         
                                         Logbooks.list[[i]]$County[[row_i_index]], #logbook county
                                         Surveys.list[[j]]$County[[row_j_index]],  #survey county
                                         
                                         Logbooks.list[[i]]$State[[row_i_index]], #logbook state
                                         Surveys.list[[j]]$State[[row_j_index]],  #survey state
                                         
                                         Logbooks.list[[i]]$Num_Anglers[[row_i_index]], #logbook number anglers
                                         Surveys.list[[j]]$Num_Anglers[[row_j_index]],   #survey number of anglers
                                         
                                         Logbooks.list[[i]]$Vessel_Official_Num[[row_i_index]],  #logbook vessel official number
                                         Surveys.list[[j]]$Vessel_Official_Num[[row_j_index]],   #survey vessel official number
                                         
                                         Logbooks.list[[i]]$Vessel_Name[[row_i_index]],
                                         Surveys.list[[j]]$Vessel_Name[[row_j_index]],
                                         
                                         Logbooks.list[[i]]$Anything_Caught_Flag[[row_i_index]],
                                         Surveys.list[[j]]$Anything_Caught_Flag[[row_j_index]],
                                         
                                         Logbooks.list[[i]]$Hours_Fished[[row_i_index]],
                                         Surveys.list[[j]]$Hours_Fished[[row_j_index]],
                                         
                                         Logbooks.list[[i]]$TIME[[row_i_index]],
                                         Surveys.list[[j]]$TIME[[row_j_index]],
                                         
                                         Dateresult, Countyresult,Stateresult,NumAnglersresult,
                                         VslNumberresult,VslNameresult,AnythingCaughtresult,
                                         HoursFishedresult,Timeresult) 
              
                # Append the iteration dataframe to the list with a unique name~~~~~~~~~~~~
                results_list[[paste0("Logbook.DF.Iteration_", i, "_",
                                     "Logbook.Row_", row_i_index, "_",
                                     "Survey.DF.Iteration_",j, "_",
                                     "Survey.Row_", row_j_index)]] <- iteration_df 
                
        }
      }
    } 
  }
}

#calculate total run time of for loop
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)



# Combine all temporary data frames into a single data frame

#####Date Matched DF----
DateMatched_df_withScores <- dplyr::bind_rows(results_list)

  #rename cols
  names(DateMatched_df_withScores) <- c("Logbook_RowID","Survey_RowID",
                       "LogbookDate","SurveyDate",
                       "LogbookCounty", "SurveyCounty",
                       "LogbookState","SurveyState",
                       "LogbookNumAnglers","SurveyNumAnglers",
                       "LogbookVslNumber","SurveyVslNumber",
                       "LogbookVslName","SurveyVslName",
                       "LogbookAnythingCaught","SurveyAnythingCaught",
                       "LogbookHoursFished","SurveyHoursFished",
                       "LogbookTime","SurveyTime",
                       "Dateresult","CountyResult","StateResult","NumAnglersResult",
                       "VslNumberresult","VslNameesult","AnythingCaughtresult",
                       "HoursFishedresult","Timeresult")



View(DateMatched_df_withScores)

####Save Date Matched DF----
#write.csv(DateMatched_df_withScores, paste0(Path, Outputs, "/DateMatched_df_withScores.csv"))  

#check how many survey rows were match to a logbook in the final DF, based on date alone matching between logbook and survey record (since the for-loop only inserts them into final DF IF the dates match)
NumberUniqueSurveys_MatchedAll <- length(unique(DateMatched_df_withScores$Survey_RowID))  

#check how many survey rows were match to a logbook in the final DF, based on date alone matching between logbook and survey record (since the for-loop only inserts them into final DF IF the dates match)
nrow(DateMatched_df_withScores)  


###True Match DF; set Vsl ID Threshold to 1, meaning vessel # or name have to match exactly ----
VslIDThreshold <- 1  #all pairs have to have a matching VSL name and number

#now filter by using Threshold Approach
TrueMatch_DF <- DateMatched_df_withScores %>% 
  #filter using the VslIDThreshold (assuming either the number or name has to be within some min "error" range, here 1.0 threshold to eliminate those that aren't close)
  filter(VslNumberresult >= VslIDThreshold | VslNameesult >= VslIDThreshold) %>%
    #however, need to exclude those vessel name matches when vessel name is "UNNAMED"
      filter(VslNumberresult >= VslIDThreshold & !LogbookVslName == "UNNAMED" & !SurveyVslName == "UNNAMED")

#check how many matches in the final DF
nrow(TrueMatch_DF) 

#########Save final_Logbook_RowsToKeep_MatchedCounty DF
#write.csv(TrueMatch_DF, paste0(Path, Outputs, "/TrueMatch_DF.csv"))



#Histogram showing the distribution of record linkage scores for unique matches among all possible matches----

##create a list of cols of data to plot
fields_toPlot <- grep(c("result"), names(TrueMatch_DF), value = TRUE, ignore.case = TRUE)

##select just cols of interest from TrueMatchDF
TrueMatchDF_short <- TrueMatch_DF %>%
  select(all_of(fields_toPlot)) %>%
  select(-Dateresult) # Assuming Dateresult is removed as in your code

##Reshape the data
TrueMatchDF_long <- TrueMatchDF_short %>%
  pivot_longer(cols = everything(), 
               names_to = "variable",
               values_to = "value") %>%
  ##Explicitly define the high reliability fields based on your requirement
  mutate(
    reliability_group = case_when(
      variable %in% c("AnythingCaughtresult", "CountyResult", "VslNumberresult", "StateResult") ~ "High Reliability (Single Score)",
      TRUE ~ "Lower Reliability (Score Range)"
    )
  )

##Summarize the data first to get counts for geom_col
TrueMatchDF_summary <- TrueMatchDF_long %>%
  group_by(variable, value, reliability_group) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # CRITICAL STEP: Sort the summary data frame by reliability_group before plotting
  arrange(reliability_group)

##Define the new, cleaner facet labels using a named vector
facet_labels <- c(
  "AnythingCaughtresult" = "Anything Caught",
  "CountyResult" = "County",
  "StateResult" = "State",
  "VslNumberresult" = "Vessel Number",
  "HoursFishedresult" = "Hours Fished",
  "NumAnglersResult" = "Number of Anglers",
  "Timeresult" = "Trip Time"
)

#plot
ggplot(TrueMatchDF_summary %>% filter(reliability_group == "Lower Reliability (Score Range)"),
       aes(x = value, weight = count)) +
  # Use geom_col for pre-counted data and increase the width
  geom_histogram(color = "white", bins = 20) +
  labs(x = "Similarity Score",
       y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ fct_inorder(variable), ncol = 1, scales = "fixed", labeller = as_labeller(facet_labels)) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")
  )

#save plot to tiff
# ggsave(paste0(Path, Outputs,"/Histograms- Distribution of Similiarity Scores for TrueMatch_DF.tiff"), 
#        units="in", width=15, height=10, dpi=600)#, compression = 'lzw')


#Combined Threshold Evaluation (Grid Search) ----
#This involves testing every reasonable combination of thresholds for your three variables and determining which combination yields the optimal overall performance (e.g., highest F1 score)

##Prepare Data, add a column to the true match DF to flag is_match ="1"
DateMatched_df_withScores_matchFlagAdded <- DateMatched_df_withScores %>%
  #For some reason there are rows with the same logbook rowID matched to the same SurveyRowID, so need to eliminate exact duplciates 
  distinct(Logbook_RowID, Survey_RowID, .keep_all = TRUE) %>% #.keep_all keeps all other columns associated with the first unique instance.

  #Now, create the column first as "0"
  mutate(is_match = "0") %>% 
  left_join(
    TrueMatch_DF %>% 
      select(Logbook_RowID, Survey_RowID) %>% 
      mutate(is_match_temp = "1"), # Use a temporary name
    by = c("Logbook_RowID", "Survey_RowID")
  ) %>%
  # Update is_match where a match was found (is_match_temp is "1")
  mutate(is_match = if_else(!is.na(is_match_temp), "1", is_match)) %>%
  # Remove the temp column
  select(-is_match_temp) 

#reduce DF to just needed cols
cols_to_retain <- c("NumAnglersResult","HoursFishedresult","Timeresult","is_match" )
DateMatched_df_withScores_matchFlagAdded_short <- DateMatched_df_withScores_matchFlagAdded %>% 
  select(all_of(cols_to_retain))
  
#Define the evaluation dataset
df_eval <- DateMatched_df_withScores_matchFlagAdded_short 


##Define the thresholds to test ---- 
#use the min and max scores for hours fished, number of anglers, and time linking variables

##first limit the DF to just true matches
DateMatched_df_withScores_matchFlagAdded_short_truematches <- DateMatched_df_withScores_matchFlagAdded_short %>%
  filter(is_match == "1")

##now calc the min and max scores for each of these linking variables, to define thresholds
SimiliarityScores_OnlyMinMaxScores <- tibble(row_id = c("NumAnglers_SimScore",
                                           "HoursFished_SimScore",
                                           "TripTime_SimScore"),
                                
                                           Min = c(min(DateMatched_df_withScores_matchFlagAdded_short_truematches$NumAnglersResult),
                                           min(DateMatched_df_withScores_matchFlagAdded_short_truematches$HoursFishedresult),
                                           min(DateMatched_df_withScores_matchFlagAdded_short_truematches$Timeresult)),
                                           
                                           Max = c(max(DateMatched_df_withScores_matchFlagAdded_short_truematches$NumAnglersResult),
                                           max(DateMatched_df_withScores_matchFlagAdded_short_truematches$HoursFishedresult),
                                           max(DateMatched_df_withScores_matchFlagAdded_short_truematches$Timeresult)))

##define the range of thresholds to test - test in increments of 0.05 (adjust as needed for precision).
##anglers
thresholds_anglers <- seq(
  from = 0,
  to = 1,
  by = 0.2
)

##time
thresholds_time <- seq(
  from = SimiliarityScores_OnlyMinMaxScores %>% filter(row_id == "HoursFished_SimScore") %>% pull(Min),
  to = SimiliarityScores_OnlyMinMaxScores %>% filter(row_id == "HoursFished_SimScore") %>% pull(Max),
  by = 0.05
)

##hoursFished
thresholds_hoursfished <- seq(
  from = SimiliarityScores_OnlyMinMaxScores %>% filter(row_id == "TripTime_SimScore") %>% pull(Min),
  to = SimiliarityScores_OnlyMinMaxScores %>% filter(row_id == "TripTime_SimScore") %>% pull(Max),
  by = 0.05
)

###Create testing grid - by generating every possible combination of these thresholds ----
threshold_grid <- expand.grid(
  thresh_anglers = thresholds_anglers,
  thresh_time = thresholds_time,
  thresh_hoursfished = thresholds_hoursfished
)

# --- Function to apply thresholds and calculate performance metric (F1 Score) ---
calculate_f1 <- function(ang_t, time_t, hours_t, data) {
  
  #Generate predictions based on thresholds across the FULL dataset
  eval_data <- data %>%
    mutate(
      #Ensure truth is a factor with levels "1" (Match) and "0" (Non-match)
      is_match = factor(is_match, levels = c("1", "0")),
      #Create an estimate column based on the thresholds
      pred_match = factor(ifelse(
        NumAnglersResult >= ang_t & 
          Timeresult >= time_t & 
          HoursFishedresult >= hours_t, "1", "0"), 
        levels = c("1", "0"))
    )
  
  ##Calculate F1 using yardstick package
  #event_level = "first" ensures "1" is treated as the positive class
  f1_val <- eval_data %>%
    yardstick::f_meas(truth = is_match, estimate = pred_match, event_level = "first") %>%
    pull(.estimate)
  
  return(ifelse(is.na(f1_val), 0, f1_val))
}

###Iterate over the grid ----
results <- threshold_grid %>%
  rowwise() %>%
  mutate(
    f1_score = calculate_f1(thresh_anglers, thresh_time, thresh_hoursfished, df_eval)
  ) %>%
  ungroup()

#Find the combination that produced the highest F1 score
optimal_combination <- results %>%
  arrange(desc(f1_score)) %>%
  slice(1)

print("Optimal Threshold Combination:")
print(optimal_combination)

# A tibble: 1 × 4
# thresh_anglers thresh_time thresh_hoursfished f1_score
# <dbl>       <dbl>              <dbl>    <dbl>
# 0.533      0.0435              0.501    0.225

###Plot optimal threshold combination----
##use ggplot to create heat-map to visualize every combination of the threshold combinations to show optimal combination
#Identify the optimal point for highlighting
opt_point <- optimal_combination %>% slice_max(f1_score, n = 1, with_ties = FALSE)

#Create the faceted heat map
# We facet by 'thresh_anglers' to see how the other two interact at different levels
ggplot(results %>% filter(thresh_time <= 0.5) %>% rename("Angler Threshold" = thresh_anglers),
       aes(x = thresh_time, y = thresh_hoursfished, fill = f1_score)) +
  geom_tile() +
  # Use a colorblind-friendly and perceptually uniform scale
  scale_fill_viridis_c(option = "magma", name = "F1 Score") +
  # Create panels for each Angler Threshold level
  facet_wrap(~`Angler Threshold`, labeller = label_both, ncol = 1) +
  # Highlight the optimal combination with a distinct border or point
  geom_point(data = opt_point %>% rename("Angler Threshold" = thresh_anglers),
             aes(x = thresh_time, y = thresh_hoursfished), 
             color = "cyan", shape = 8, size = 3, stroke = 1.5) +
  # Scientific styling
  theme_bw(base_size = 12) + 
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "top",
    legend.key.width = unit(2, "cm")
  ) +
  labs(
    x = "Time Similarity Threshold",
    y = "Hours Fished Similarity Threshold",
    caption = paste0("Optimal F1: ", round(opt_point$f1_score, 3))
  )

#Save as high-resolution TIFF or PDF for publication
#ggsave(paste0(Path, Outputs, "/Threshold_Optimization_Figure.tiff"), width = 10, height = 8, dpi = 300)

#Print the total number of combinations tested - for manuscript
n_total <- nrow(results)
print(paste("Total combinations tested:", n_total))

###Extract Min and Max F1 Scores ----
# na.rm = TRUE is important if some combinations resulted in NA scores (though the function was set to return 0)
Min_F1 <- min(results$f1_score, na.rm = TRUE)
Max_F1 <- max(results$f1_score, na.rm = TRUE)

###Extract Optimal Configuration Values ---_
# We use the opt_point dataframe generated previously
Optimal_F1_Score <- opt_point$f1_score
Optimal_Anglers_Thresh <- opt_point$thresh_anglers
Optimal_Time_Thresh <- opt_point$thresh_time
Optimal_HoursFished_Thresh <- opt_point$thresh_hoursfished

###Extract Near-Optimal Threshold (e.g., F1 > 0.8) ----
# Define a subjective threshold for "near-optimal" if needed for discussion
Near_Optimal_Threshold <- 0.8 # Example: Adjust this based on your actual results
# Count how many combinations were above this
Num_Near_Optimal <- results %>% 
  filter(f1_score > Near_Optimal_Threshold) %>% 
  nrow()

#Print statements to insert into manuscript ----
cat("--- Values to Insert ---\n")
cat(sprintf("[Insert Min F1]: %.3f\n", Min_F1))
cat(sprintf("[Insert Max F1]: %.3f\n", Max_F1))
cat(sprintf("[Insert Score]: %.3f\n", Optimal_F1_Score))
cat(sprintf("[Value] for anglers: %.2f\n", Optimal_Anglers_Thresh))
cat(sprintf("[Value] for trip timing: %.2f\n", Optimal_Time_Thresh))
cat(sprintf("[Value] for hours fished: %.2f\n", Optimal_HoursFished_Thresh))
cat(sprintf("[Value] for near-optimal threshold (F1 > [Value]): %.2f (found in %d combinations)\n", Near_Optimal_Threshold, Num_Near_Optimal))



#MATCHED DF (matched on common fields; ignoring unique ID); using the DF optimal_combination for the thresholds ----

#Apply threshold values based on Distribution of Similarity Scores Histograms (Histograms- Distribution of Similarity Scores for Matched Fields_Vessel ID Threshold.jpeg)
Site_threshold <- 1 #intuitively - sites need to match so this should be a 1
AnythingCaughtThreshold <- 1
Timethreshold <- optimal_combination$thresh_time
HoursFishedThreshold <- optimal_combination$thresh_hoursfished
NumAnglersThreshold <- optimal_combination$thresh_anglers


#Use threshold approach to filter out mismatches
fullMatch_DF <- DateMatched_df_withScores %>% 
  #filter out only those with the same Site (County and State)
  filter(CountyResult >= Site_threshold & StateResult >= Site_threshold) %>%
  #now use Anything Caught threshold to weed out some unlikely matches (where a threshold of 1 means they identical, both 1 or both 0, among the survey and logbook DFs)
  filter(AnythingCaughtresult == AnythingCaughtThreshold) %>%
  #now filter by TimeScore
  filter(Timeresult >= Timethreshold) %>%
  #now filter by hours fished threshold
  filter(HoursFishedresult >= HoursFishedThreshold) %>%
  #now only keep rows where anglers are equal (per accsp method)
  filter(NumAnglersResult >= NumAnglersThreshold) 


#check how many matches in the final DF
nrow(fullMatch_DF) 

#########Save final_Logbook_RowsToKeep_MatchedCounty DF
#write.csv(fullMatch_DF, paste0(Path, Outputs, "/fullMatch_DF.csv"))


#check how many rows are true matches, using vessel ID threshold
fullMatch_DF_IncludeVesselIDThreshold_check <- fullMatch_DF %>%
  filter(VslNumberresult >= VslIDThreshold | VslNameesult >= VslIDThreshold) %>%
  #however, need to exclude those vessel name matches when vessel name is "UNNAMED"
  filter(VslNumberresult >= VslIDThreshold & !LogbookVslName == "UNNAMED" & !SurveyVslName == "UNNAMED") %>%
  #For some reason there are 2 rows with the same logbook rowID matched to the same SurveyRowID, so need to eliminate duplciates 
  distinct(Logbook_RowID, Survey_RowID, .keep_all = TRUE) # .keep_all keeps all other columns associated with the first unique instance.

nrow(fullMatch_DF_IncludeVesselIDThreshold_check) 
