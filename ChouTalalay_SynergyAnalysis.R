# Purpose: This script processes data from high-throughput screening assays for drug combination synergy using the Chou-Talalay method.
# Filter and group the results data frame by relevant columns and calculate the standard deviation, replicate number, and mean value for each group.
# Set the names of the drug combinations and create output directories.
# Normalize drug concentrations, normalize values between zero and 1, and calculate Fa, Fu, and FaFu.
# Sort the data frame by Drug1 and Drug2, and count the number of unique Drug1 and Drug2 concentrations.
# Plot the individual dose responses and eliminate doses that are outside the linear range of the median.
# Normalize data between 1 and 0 and separate and process the data for each drug.
# Adjust the MedianEffect range if necessary, fit median effect curve, calculate Dx1 and Dx2, and calculate the CI.
# Assign the synergism "ranks" based on the CI values and create a function to look up the synergy "values" based on the CI values.
# Determine the mean CI and standard error of the mean CI.
# Save the synergy data to CSV files.
# Export the data to generate plots.

# Load required packages and install missing ones
load_required_packages <- function() {
  list.of.packages <- c("dplyr", "reshape", "scales", "reshape", "data.table", "stringr", "tidyr", "tidyext")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
}

#### IF RUNNING FROM RSTUDIO ######
# Set the working directory to the directory of the current file 
# set_working_directory <- function() {
#   current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#   setwd(current_working_dir)
# }

data_location = "ExampleData"

# Retrieve the list of experiments
get_experiments <- function(data_location) {
  experiments = list.dirs(data_location)
  experiments = gsub(paste0(data_location, "/"), "", experiments)
  experiments = experiments[-grep(data_location, experiments)]
  return(experiments)
}

# Create output directories for the results and data
create_output_directories <- function() {
  dir.create("Output")
  dir.create("Output/ByExperiment")
  dir.create("Output/Data")
  dir.create("Output/Data/CI_Data")
  dir.create("Output/Data/Exported_Data")
}

# Load data for each experiment
load_experiment_data <- function(data_location, exp) {
  results <- fread(paste0(data_location, "/", exp, "/responses.csv"), stringsAsFactors = FALSE, fill = TRUE)
  results$BlockId = as.character(results$BlockId)
  results$Row = as.character(results$Row)
  results$Col = as.character(results$Col)
  
  metadata <- fread(paste0(data_location,"/", exp, "/metadata.csv"), stringsAsFactors = FALSE, fill = TRUE)
  metadata$BlockId = as.character(metadata$BlockId)
  return(list(results, metadata))
}

# Parse rows and cols from metadata
parse_rows_and_cols <- function(metadata) {
  rows = as.data.frame(str_split_fixed(metadata$RowConcs, pattern = ",", n = 10), stringsAsFactors = FALSE)
  rows = mutate_all(rows, function(x) as.numeric(as.character(x)))
  colnames(rows) = 1:10
  rows$Drug = metadata$RowName
  rows$BlockId = metadata$BlockId
  rows = rows %>% pivot_longer(as.character(1:10), names_to = "Row")
  
  colnames(rows) = c("Name1", "BlockId", "Row", "Drug1")
  rows$Row = as.character(rows$Row)
  
  cols = as.data.frame(str_split_fixed(metadata$ColConcs, pattern = ",", n = 10), stringsAsFactors = FALSE)
  cols = mutate_all(cols, function(x) as.numeric(as.character(x)))
  colnames(cols) = 1:10
  cols$Drug = metadata$ColName
  cols$BlockId = metadata$BlockId
  cols = cols %>% pivot_longer(as.character(1:10), names_to = "Row")
  colnames(cols) = c("Name2", "BlockId", "Col", "Drug2")
  cols$Col = as.character(cols$Col)
  cols = cols %>% filter(!is.na(Drug2))
  
  return(list(rows, cols))
}

# Merge results, rows, and cols
merge_data <- function(results, rows, cols) {
  results = data.table(results)
  rows = data.table(rows)
  cols = data.table(cols)
  merge = merge(results, rows, by = c("BlockId", "Row"), all.x = TRUE)
  merge = merge(merge, cols, by = c("BlockId", "Col"), all.x = TRUE)
  return(merge)
}

# Organize drug combos
organize_drug_combos <- function(results) {
  drugcombos = unique(results %>% select(BlockId, Name1, Name2)) %>% arrange(BlockId)
  rownames(drugcombos) = 1:nrow(drugcombos)
  
  count1 = length(unique(drugcombos$Name1))
  count2 = length(unique(drugcombos$Name2))
  if (count2 < count1) {
    blockids = drugcombos$BlockId
    drugcombos = drugcombos %>% select(-BlockId)
    colnames(drugcombos) = rev(colnames(drugcombos))
    drugcombos$BlockId = blockids
    colnames(results) = gsub(pattern = "1", "XX", colnames(results))
    colnames(results) = gsub(pattern = "2", "1", colnames(results))
    colnames(results) = gsub(pattern = "XX", "2", colnames(results))
  }
  return(drugcombos)
}


### MAIN SCRIPT ###
load_required_packages()
set_working_directory()
experiments <- get_experiments(data_location = "ExampleData/")
create_output_directories()
dose_cutoff <- 5
data_location = "data" # path to where data is stored

exp = experiments[1] # example for debugging
for (exp in experiments) {
  print(paste("############ Analyzing folder......", exp, "   ############ "))
  
  # save output location by experiment
  output_location = paste0("Output/ByExperiment/", exp, "/")
  dir.create(output_location)
  
  data <- load_experiment_data(data_location = "ExampleData", exp)
  results <- data[[1]]
  metadata <- data[[2]]
  
  parsed_data <- parse_rows_and_cols(metadata)
  rows <- parsed_data[[1]]
  cols <- parsed_data[[2]]
  
  results <- merge_data(results, rows, cols)
  drugcombos <- organize_drug_combos(results)
  
  list_of_drug_combo_ids = unique(drugcombos$BlockId)
  
  
  id = list_of_drug_combo_ids[1] # example for debugging
  for (id in list_of_drug_combo_ids) {
    
    # Filter and group the results data frame by BlockId, Col, Row, Name1, Name2, Drug1, and Drug2
    # Calculate the standard deviation, replicate number, and mean value for each group
    compiled_plates <- results %>%
      filter(BlockId == id) %>%
      arrange() %>%
      group_by(BlockId, Col, Row, Name1, Name2, Drug1, Drug2) %>%
      summarise(sd = sd(Value), Replicate = n(), Value = mean(Value))
    
    # Set the names of the drug combinations
    savename1 <- compiled_plates$Name1[1]
    savename2 <- compiled_plates$Name2[1]
    
    # Print the drug combination being analyzed
    print(paste(exp, "...Analyzing Combination.......", savename1, " + ", savename2))
    
    # Append "_1" and "_2" to the drug names
    compiled_plates$Name1 <- paste0(compiled_plates$Name1, "_1")
    compiled_plates$Name2 <- paste0(compiled_plates$Name2, "_2")
    
    # Create output directories
    combo_output <- paste0(output_location, savename1, "_", savename2, "/")
    dir.create(combo_output)
    
    # Round drug concentrations to 6 decimal places
    compiled_plates$Drug1 <- round(compiled_plates$Drug1, digits = 6)
    compiled_plates$Drug2 <- round(compiled_plates$Drug2, digits = 6)
    
    # Normalize to average no-dose condition
    max_value <- compiled_plates %>%
      filter(Drug1 == 0, Drug2 == 0) %>%
      summarise(Avg = mean(Value)) %>%
      select(Avg)
    max_value <- max_value$Avg
    
    # Normalize values between zero and 1 and prevent log transformation issues by adding 0.01
    compiled_plates <- compiled_plates %>%
      group_by(Name1, Name2) %>%
      mutate(Value = .99 * (0.01 + (Value - min(Value)) / (max_value - min(Value)))) %>%
      mutate(Value = ifelse(Value > 1, .99, Value))
    
    # Calculate Fa, Fu, and FaFu
    compiled_plates$Fa <- 1 - compiled_plates$Value          # Fa = 1 - proliferation value
    compiled_plates$Fu <- compiled_plates$Value              # Fu = proliferation value
    compiled_plates$FaFu <- compiled_plates$Fa / compiled_plates$Fu  # Fa/Fu
    
    # Sort the data frame by Drug1 and Drug2
    compiled_plates <- compiled_plates %>% arrange(Drug1, Drug2)
    
    # Count the number of unique Drug1 and Drug2 concentrations
    count1 = compiled_plates %>% group_by(Drug1) %>% summarise(count = n())
    count1$log = log10(count1$Drug1)
    
    count2 = compiled_plates %>% group_by(Drug2) %>% summarise(count = n())
    count2$log = log10(count2$Drug2)
    
    # Plot the individual dose responses
    individual_doses = rbind(compiled_plates %>% filter(Drug1 == 0) %>% mutate(Name = Name2, Drug = Drug2, DrugNum = 2),
                             compiled_plates %>% filter(Drug2 == 0) %>% mutate(Name = Name1, Drug = Drug1, DrugNum = 1))
    
    name1 = unique(individual_doses$Name1)
    name2 = unique(individual_doses$Name2)
    
    
    # Eliminate doses that are outside inflection point of the IC50 (aka will be outside the median effect line)  --------
    average_dose <- individual_doses %>%
      group_by(Drug, Name) %>%
      mutate(Average = Value) %>%
      arrange(Name, Drug) %>%
      ungroup() %>%
      group_by(Name) %>%
      mutate(min = Drug[2]) %>%
      data.frame()
    
    # Normalize data between 1 and 0
    drugs <- unique(average_dose$Name)
    
    # Separate and process the data for each drug
    average_dose1 <- average_dose %>%
      ungroup() %>%
      filter(Name == drugs[1]) %>%
      mutate(Dist = abs(Average - 0.5))
    
    average_dose2 <- average_dose %>%
      ungroup() %>%
      filter(Name == drugs[2]) %>%
      mutate(Dist = abs(Average - 0.5))
    
    # Combine the processed data for both drugs
    average_dose <- rbind(average_dose1, average_dose2)
    
    # Sort the data by name and drug, group by name
    average_dose <- average_dose %>%
      arrange(Name, Drug) %>%
      group_by(Name) %>%
      
      # Find the inflection point for each name
      mutate(Inflection = ifelse(Dist == min(Dist), row_number(), NA)) %>%
      mutate(Inflection = mean(Inflection, na.rm = TRUE)) %>%
      
      # Calculate the distance between each row and the inflection point
      mutate(Neighbor = Inflection - row_number()) %>%
      
      # Determine which doses are within the linear range
      mutate(MedianEffect = ifelse(abs(Neighbor) < dose_cutoff / 2, TRUE, FALSE)) %>%
      
      # Count the number of doses within the linear range
      mutate(count = length(which(MedianEffect)))
    
    # Adjust the MedianEffect range if necessary
    for (c in 1:10) {
      average_dose <- average_dose %>%
        # Use case_when to modify MedianEffect
        mutate(MedianEffect = case_when(
          # If count is less than dose_cutoff and Neighbor is less than c, set MedianEffect to TRUE
          count < dose_cutoff & abs(Neighbor) < c ~ TRUE,  # both tests: group A
          # If count is greater than or equal to dose_cutoff - 1, don't modify MedianEffect
          count > dose_cutoff - 1 ~ MedianEffect
        ),
        # Recalculate the count of doses within the linear range
        count = length(which(MedianEffect)))
    }
    
    # Pull out which is the zero dose point (since this always needs to be included)
    zeros <- average_dose %>%
      filter(Drug == 0) %>%
      select(-min)
    
    average_dose <- average_dose %>%
      filter(MedianEffect == TRUE) %>%
      select(-min)
    
    # Keep ones that are less than the previous
    if (nrow(average_dose) > 0) {
      # Calculate the dose range for each drug for each Name
      dose_range = average_dose %>% group_by(Name) %>%
        filter(!Drug == 0) %>%
        mutate(max = max(Drug), min = min(Drug)) %>%
        select(Name, min, max) %>% unique()
      
      # Add the information about the no-drug conditons the average_dose data frame, then remove any duplicate rows
      average_dose = rbind(zeros, average_dose) %>% unique() %>% arrange(Name, Drug)
      
      # Merge the average_dose data frame with the dose range data frame to add min and max columns to each row,
      # then modify the Drug column by setting it to min/2 if it is equal to zero
      average_dose = merge(average_dose, dose_range, all.x = TRUE) %>%
        mutate(Drug = ifelse(Drug == 0, min/2, Drug))
      
      # Create two data frames (range1 and range2) that contain the dose range data for each drug for Name1 and Name2
      range1 = dose_range %>% filter(Name == name1)
      range2 = dose_range %>% filter(Name == name2)
      colnames(range1) = paste0(colnames(range1),1)
      colnames(range2) = paste0(colnames(range2),2)
      dose_range2 = cbind(range1, range2)
      
      # Merge dose_range2 with individual_doses
      individual_doses = merge(individual_doses, dose_range2) 
      
      # Filter individual_doses to only include rows where Drug1 is either zero or a unique dose for Name1, and
      # Drug2 is either zero or a unique dose for Name2
      doses1 = average_dose %>% filter(Name == name1) %>% pull(Drug) %>% sort() %>% unique()
      doses2 = average_dose %>% filter(Name == name2) %>% pull(Drug) %>% sort() %>% unique()
      individual_doses = individual_doses %>% filter(Drug1 %in% c(0,doses1[-1]) & Drug2 %in% c(0,doses2[-1]))
      
      # Modify Drug1 and Drug2 columns by setting them to min1/2 and min2/2 respectively if they are equal to zero,
      # then arrange individual_doses by Name, Drug1, and Drug2
      individual_doses = individual_doses %>%
        mutate(Drug1 = ifelse(Drug1 == 0, min1/2, Drug1)) %>%
        mutate(Drug2 = ifelse(Drug2 == 0, min2/2, Drug2)) %>%
        arrange(Name, Drug1, Drug2)
      
      
      compiled_plates = merge(compiled_plates, dose_range2)  %>%
        filter(Drug1 %in% c(0,doses1[-1]) & Drug2 %in% c(0,doses2[-1])) %>% 
        mutate(Drug1 = ifelse(Drug1 == 0, min1/2, Drug1)) %>%
        mutate(Drug2 = ifelse(Drug2 == 0, min2/2, Drug2)) %>%
        arrange(Drug1, Drug2) 
      
      
      compiled_plates = compiled_plates %>% select(-Row, -Col)
      
      doses1 = compiled_plates %>% select(Name1, Drug1) %>% unique() %>% arrange(Drug1) %>%
        mutate(Group1 = row_number())
      doses2 = compiled_plates %>% select(Name2, Drug2) %>% unique() %>% arrange(Drug2) %>%
        mutate(Group2 = row_number())
      
      compiled_plates = merge(compiled_plates, doses1)
      compiled_plates = merge(compiled_plates, doses2)
      
      compiled_plates$Group1 = factor(compiled_plates$Group1, levels = 1:10)
      compiled_plates$Group2 = factor(compiled_plates$Group2, levels = 1:10)
      compiled_plates = compiled_plates %>% arrange(Drug1, Drug2)
      
      compiled_plates_avg = compiled_plates %>% group_by(Name1, Drug1, Group1,
                                                         Name2, Drug2, Group2) %>%
        summarise(sd = sd(Value)/sqrt(n()), Mean = mean(Value))
      
      
      compiled_plates_avg = compiled_plates_avg %>% arrange(Group1, Group2)
      
      
      # Calculate Combination Index from linear Median Effect Equation ----------
      # Skip if there is an error or non-linearity, etc
      # Helper function to calculate synergy description based on CI values
      syn_desc <- c("Strong synergism", "Synergism", "Additive", "Antagonism", "Strong antagonism")
      
      syn_symbols <- function(ci) {
        syn_val <- c(0.1, 0.5, 1, 1.45, 3)
        
        symb <- ci
        for (i in 1:length(syn_val)) {
          symb[ci > syn_val[i]] <- syn_desc[i]
        }
        return(symb)
      }
      
      # Check for lm.fit error
      check_for_error <- function(lm_call) {
        possible_error <- tryCatch(lm_call, error = function(e) e)
        return(!grepl("lm.fit", as.character(possible_error$call)[1]))
      }
      
      # Perform regression and calculate Dm
      regression_and_dm <- function(data, dependent_var, independent_var) {
        if (check_for_error(lm(log10(data[[dependent_var]]) ~ log10(data[[independent_var]])))) {
          model <- lm(log10(data[[dependent_var]]) ~ log10(data[[independent_var]]))
          r_squared <- summary(model)$r.squared
          coefficients <- coef(model)
          df <- data.frame(t(coefficients))
          colnames(df) <- c("b", "m")
          df$Dm <- 10^(-df$b / df$m)
          return(list("model" = df, "r_squared" = r_squared))
        } else {
          return(NULL)
        }
      }
      
      # Get just the data for Drug1 alone (no combination)
      d1 =  compiled_plates %>%
        filter(Group2 == 1) %>%
        filter(!Drug1 == min(Drug1))
      
      # Get just the data for Drug2 alone (no combination)
      d2 = compiled_plates %>%
        filter(Group1 == 1) %>%
        filter(!Drug2 == min(Drug2))
      
      # Main part of the code
      # Perform regression and Dm calculation for Drug1
      drug1_result <- regression_and_dm(d1, "FaFu", "Drug1")
      if (!is.null(drug1_result)) {
        d1 <- drug1_result$model
        r2_1 <- drug1_result$r_squared
        
        # Perform regression and Dm calculation for Drug2
        drug2_result <- regression_and_dm(d2, "FaFu", "Drug2")
        if (!is.null(drug2_result)) {
          d2 <- drug2_result$model
          r2_2 <- drug2_result$r_squared
          
          # Create R2_plot dataframe
          R2_plot <- data.frame(Name = c(name1, name2),
                                Drug = c(min(individual_doses$Drug1), min(individual_doses$Drug2)), 
                                R2 = c(r2_1, r2_2),
                                R2_text = paste("R^2 =", round(c(r2_1, r2_2), digits = 2)),
                                FaFu = max(individual_doses$FaFu))
          
          # Update individual_doses dataframe
          individual_doses <- individual_doses %>% mutate(Drug = ifelse(Name == name1 & Drug == 0, min1/2,Drug)) %>%
            mutate(Drug = ifelse(Name == name2 & Drug == 0, min2/2,Drug)) %>%
            filter(!(Drug1 == min1/2 & Drug2 == min2/2))
          
          # Calculate Dx1 and Dx2, and create cell_dat dataframe 
          cell_dat <- compiled_plates %>% group_by(Drug1, Drug2, Name1, Name2) %>%
            filter(!Drug1 == min1/2, !Drug2 == min2/2) %>% 
            summarize(Fa = mean(Fa), Fu = mean(Fu))
          
          # Dx = Dm[fa/fu]^1/m
          # Update cell_dat dataframe
          cell_dat <- cell_dat %>% filter(!is.na(Fa),
                                          !is.na(Drug1),
                                          !is.na(Drug2)) %>%
            mutate(Dx1 = d1$Dm * (Fa/Fu)^(1/d1$m),
                   Dx2 = d2$Dm * (Fa/Fu)^(1/d2$m)) %>%
            mutate(I1 = Drug1/Dx1,
                   I2 = Drug2/Dx2,
                   DRI = Dx1/Drug1 + Dx2/Drug2,
                   DRI1 = Dx1/Drug1,
                   DRI2 = Dx2/Drug2) %>%
            mutate(CI = I1 + I2)
          
          # Create isobol_dat dataframe
          isobol_dat <- cell_dat %>% filter(!is.na(CI))
          isobol_dat$CIorg <- isobol_dat$CI
          
          # Assign synergy descriptions
          isobol_dat$CI_desc <- factor(syn_symbols(isobol_dat$CI), levels = syn_desc)
          
          # Truncate CI values and calculate mean CI and its standard error
          if (max(isobol_dat$CI) > 3) {
            isobol_dat$CI[which(isobol_dat$CI > 3)] <- 3
          }
          
          ci_max <- ifelse(max(isobol_dat$CI) > 2, max(isobol_dat$CI), 2)
          mean_CI <- round(mean(isobol_dat$CI), 3)
          mean_CI_se <- round(sd(isobol_dat$CI) / sqrt(length(isobol_dat$CI)), 3)
          
          # Save data to CSV file
          if (nrow(isobol_dat) > 0) {
            isobol_dat <- isobol_dat %>% data.frame() %>% mutate(Experiment = exp) %>%
              mutate(rsquared1 = r2_1,
                     rsquared2 = r2_2)
            
            write.csv(x = isobol_dat, paste0(combo_output, Sys.Date(), " ", savename1, "_", savename2, "_SynergyData.csv"))
            write.csv(x = isobol_dat, paste0("Output/Data/CI_data/", Sys.Date(), " ", exp, "-", savename1, "_", savename2, "_SynergyData.csv"))
            
            
            export_data = compiled_plates %>%
              mutate(Drug1 = ifelse(Drug1 == min1/2, 0, Drug1)) %>%
              mutate(Drug2 = ifelse(Drug2 == min2/2, 0, Drug2))
            
            export_data_avg = export_data %>%
              group_by(Name1, Name2, Drug1, Drug2) %>%
              summarise(sd = sd(Value)/sqrt(n()), Average = mean(Value), 
                        Fa = 1-Average,
                        Fu = Average,
                        FaFu = Fa/Fu) %>%
              ungroup() %>% 
              dplyr::select(-Name1, -Name2)
            
            drugnames =  c(savename1, savename2)
            drugnames = gsub(" ", "_", drugnames)
            colnames(export_data_avg)[1:2] = drugnames
            
            colnames(export_data_avg) = gsub(" ", "-", colnames(export_data_avg))
            export_data_avg = export_data_avg %>% data.frame() %>% mutate(Experiment = exp)
            
            write.csv(x = export_data_avg, paste0(combo_output, Sys.Date(),savename1, "_", savename2,"_SummaryData.csv"))
            write.csv(x = export_data_avg, paste0("Output/Data/Exported_Data/", Sys.Date(), " ", exp, "-",savename1, "_", savename2,"_SummaryData.csv"))
          }
        }
      }
    }
  }
}



