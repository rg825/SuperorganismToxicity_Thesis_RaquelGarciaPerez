# This script presents the complete analysis used to assess colony-level responses to chronic cadmium 
# exposure in Bombus terrestris audax as a model pollinator species.

# Setting up the working directory:
setwd("~/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB")

# Cleaning R environment:
rm(list = ls()) 
graphics.off()

# Loading the packages required for the complete data analysis:
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(zoo)
library(lmerTest)
library(ggeffects)
library(av)
library(terra)
library(reticulate)
library(jpeg)
library(ggforce)

# 1.1 HYPOTHESIS 1: DATA LOADING AND WRANGLING ---------------------------------------------------------

# Defining all different colony IDs:
colonies <- c("C6", "C14", "C1", "C2", "C8", "C10", "C4", "C5", "C7", "C12", "C9", "C11")

# Loading the file containing temperature measurements for all colonies:
abiotic_data <- read.csv("temperature_humidity.csv", header = TRUE)

# Transforming the timestamp column into date-time format:
abiotic_data$Time <- mdy_hms(abiotic_data$Time)

# Filtering the dataset to include only measurements collected during the experimental exposure 
# period:

start_date <- as.Date("2026-06-17")     # Start day of the experiment (Day = 0).

abiotic_data <- abiotic_data %>% filter(Time >= ymd_hm("2026-06-17 14:30"),
                                        Time <= ymd_hm("2026-07-17 12:00"))

# Reformatting, creating useful columns and assigning colonies to sensors:
abiotic_data <- abiotic_data %>% 
  mutate(Date = as.Date(Time),
         Day = as.integer(Date - start_date),  # Calculating the experimental day
         Date = format(Time, "%d/%m/%Y"),
         Hour = format(Time, "%H:%M"),
         Colony = recode(as.character(Sensor),
                         "50" = "C6",
                         "51" = "C14",
                         "52" = "C1",
                         "53" = "C2",
                         "54" = "C8",
                         "55" = "C10",
                         "56" = "C4",
                         "57" = "C5",
                         "58" = "C7",
                         "59" = "C12",
                         "60" = "C11",
                         "61" = "C9"))

# Assigning the Cadmium concentrations to each colony:
abiotic_data$Concentration <- NA
abiotic_data$Concentration[abiotic_data$Colony %in% c("C6","C14")] <- "0 mg/L"
abiotic_data$Concentration[abiotic_data$Colony %in% c("C1","C2")] <- "0.02 mg/L"
abiotic_data$Concentration[abiotic_data$Colony %in% c("C8","C10")] <- "0.2 mg/L"
abiotic_data$Concentration[abiotic_data$Colony %in% c("C4","C5")] <- "1 mg/L"
abiotic_data$Concentration[abiotic_data$Colony %in% c("C7","C12")] <- "2 mg/L"
abiotic_data$Concentration[abiotic_data$Colony %in% c("C9","C11")] <- "20 mg/L"

# Keeping only columns of interest:
abiotic_data <- abiotic_data %>% select(Colony,
                                        Concentration,
                                        Date,
                                        Day,
                                        Hour,
                                        Temperature,
                                        Time)

# Calculating daily summary statistics for temperature of all colonies:
daily_summary <- abiotic_data %>%
  group_by(Colony, Concentration, Date, Day) %>%
  summarise(mean_temp = mean(Temperature),
            sd_temp = sd(Temperature),
            min_temp = min(Temperature),
            max_temp = max(Temperature),
            .groups = "drop") %>%
  arrange(Day, Concentration)



# 1.2 HYPOTHESIS 1: INSPECTING ANOMALIES -----------------------------------------------------------------

# Sometimes sensors would fall inside the colony, which could be translated into big shifts in temperature.
# To examine this possible noise, I visualise the daily standard deviation to check if there is any big
# spread of data for any sensor any day.

ggplot(daily_summary,
       aes(Day, sd_temp)) +
  geom_line() +
  facet_wrap(~Colony)

# Colonies C2, C8 and C12 show unusually high daily variability. I inspect each case independently:

C2_Day26 <- abiotic_data %>% filter(Colony == "C2", 
                                    Day == 26)

ggplot(C2_Day26,                              # Temperature shows a decline between 11:00 - 14:00.
       aes(x = Time,                          # The sensor was directly taken from inside the colony, 
           y = Temperature)) +                # when it fell.
  geom_line() +
  geom_point() +
  scale_x_datetime(date_breaks = "1 hour",
                   date_labels = "%H:%M") +
  theme_bw()

C12_Day14 <- abiotic_data %>% filter(Colony == "C12", 
                                     Day == 14)

ggplot(C12_Day14,                             # Temperature shows a rise from 07:00 - 12:00.
       aes(x = Time,                          # Due to the high risk of bee escape, the sensor 
           y = Temperature)) +                # remained inside the colony after falling for 
  geom_line() +                               # several hours until it could be repositioned 
  geom_point() +                              # safely.
  scale_x_datetime(date_breaks = "1 hour",
                   date_labels = "%H:%M") +
  theme_bw()

# The anomalies detected in C2 and C12 matched known periods when the sensors had fallen into  
# the colonies (day 26 and 14 for colonies 2 and 12, respectively). Therefore, the affected 
# time intervals were removed from the dataset to minimise measurement noise.

# Removing data associated with the time intervals when the sensors fell:
clean_abiotic_data <- abiotic_data %>%
  filter(!(Colony == "C2" & 
             Time >= ymd_hms("2026-07-13 11:00:00") &
             Time <= ymd_hms("2026-07-13 14:00:00"))) %>%
  filter(!(Colony == "C12" & 
             Time >= ymd_hms("2026-07-01 07:00:00") &
             Time <= ymd_hms("2026-07-01 12:00:00")))

# Inspecting C8, it looks like the sensor stops registering data half-way the experiment. 
# Let´s find out which day it stopped collecting data:

unique(clean_abiotic_data$Date[clean_abiotic_data$Colony == "C8"])  # 28/06/2026 was the last  
# day of the sensor collecting  
# data.                                         

C8_data <- clean_abiotic_data %>% filter(Colony == "C8")

# Plotting the mean daily temperature of C8 throughout the experimental exposure period:

ggplot(daily_summary %>%
         filter(Colony == "C8"),
       aes(x = Day,
           y = mean_temp)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  labs(x = "Experimental day",
       y = "Mean daily temperature (°C)") +
  theme_classic()

# We can confirm that due to sensor failure, temperature data from colony C8 was not collected  
# during the second half of the experiment. In this context, we can adopt different approaches to 
# sort the missing data:

# OPTION 1: Extrapolating the remaining 19 days of C8 from only the first 11 days of data.
# However, nest temperature is influenced by many factors (e.g. room temperature, colony size, 
# bee activity, ventilation...). These factors can change throughout the experiment, so I can´t 
# guarantee that a function fitted to the first half of the experiment would accurately represent
# the second half.

# OPTION 2: Predicting the remaining 19 days of C8 from the 19 last days of C10. C8 and C10 were 
# exposed to the same Cd concentration, housed in the same room, and may have experienced similar
# environmental conditions. I can check this approach by fitting a linear regression and checking
# R squared value alongside the resulting slope.

temp_C8_C10 <- daily_summary %>%                         # Extracting the mean daily temperature
  filter(Colony %in% c("C8", "C10"),        # of C8 and C10 during the first 11 days 
         Day <= 11) %>%                     # of exposure.
  select(Colony, Day, mean_temp) %>%
  pivot_wider(names_from = Colony,
              values_from = mean_temp)

# Plotting changes in mean daily temperature for both colonies:
pA <- ggplot(temp_C8_C10,                          
             aes(Day)) +
  geom_line(aes(y = C8, colour = "C8"), linewidth = 1) +
  geom_line(aes(y = C10, colour = "C10"), linewidth = 1) +
  geom_point(aes(y = C8, colour = "C8"), size = 2) +
  geom_point(aes(y = C10, colour = "C10"), size = 2) +
  labs(x = "Experimental day",
       y = "Average nest temperature (°C)",
       colour = "Colony ID") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15))

pA

# Both colonies appear to show opposite temperature trends. To examine this relationship in 
# more detail, I fit a linear regression model:

model_C8 <- lm(C8 ~ C10, data = temp_C8_C10) 

summary(model_C8)   # C10 explained only 25% of the variation in C8 (R² = 0.253). 
# Moreover, the slope was negative, suggesting that as C10 gets warmer,
# C8 becomes cooler. This pattern is not biologically feasible for two 
# colonies housed separately but under the same environmental conditions.
# These results indicate that C10 is a poor predictor of C8.

# Plotting the linear regression model:
pB <- ggplot(temp_C8_C10, aes(x = C10, y = C8)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm",
              se = TRUE,
              colour = "black",
              alpha = 0.15,
              linewidth = 1.2) +
  annotate("text",
           x = Inf, y = Inf,
           label = expression(R^2 == 0.253 ~ "," ~ p == 0.096),
           hjust = 1.2,
           vjust = 1.5,
           size = 5) +
  labs(x = "C10 average nest temperature (°C)",
       y = "C8 average nest temperature (°C)") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15))

pB

# As no reliable method was available to predict the missing C8 temperature values, the observed 
# data for C8 was retained unchanged. Merging both plots into a single one:

pA <- pA + theme(plot.margin = margin(t = 15, r = 20, b = 15, l = 15))

pB <- pB + theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15))

plots_combined_1 <- (pA | pB) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.85),
        aspect.ratio = 0.9)

plots_combined_1

# Organising the columns in the final dataset containing temperature data:
abiotic_final <- clean_abiotic_data  %>%
  group_by(Colony, Concentration, Date, Day) %>%
  summarise(mean_temp = mean(Temperature),
            sd_temp = sd(Temperature),
            min_temp = min(Temperature),
            max_temp = max(Temperature),
            .groups = "drop") %>%
  arrange(Day, Colony)

# Looking at overall descriptive statistics in temperature dataset:
abiotic_final %>% summarise(N = n(),
                            Mean = mean(mean_temp, na.rm = TRUE),
                            SD = sd(mean_temp, na.rm = TRUE),
                            Median = median(mean_temp, na.rm = TRUE),
                            Min = min(mean_temp, na.rm = TRUE),
                            Max = max(mean_temp, na.rm = TRUE),
                            Q1 = quantile(mean_temp, 0.25, na.rm = TRUE),
                            Q3 = quantile(mean_temp, 0.75, na.rm = TRUE))

# Looking at descriptive statistics per Cd treatment:
abiotic_final %>% group_by(Concentration) %>% 
  summarise(N = n(),
            Mean = mean(mean_temp, na.rm = TRUE),
            SD = sd(mean_temp, na.rm = TRUE),
            Min = min(mean_temp, na.rm = TRUE),
            Max = max(mean_temp, na.rm = TRUE),
            .groups = "drop")








# 2. HYPOTHESIS 2: DATA LOADING AND WRANGLING -----------------------------------------------------------------

# Loading the file containing pollen consumption data for all colonies:
pollen <- read.csv("pollen_consumption.csv", header = TRUE, skip = 1)

# Keeping only the columns required for the analysis:
pollen <- pollen[, c("Date", "Day", colonies)]

# I replaced pollen balls either when they were unlikely to last until the following day or, at 
# most, after 3–4 days to minimise the risk of contamination. Therefore, I standardised consumption 
# and calculated it as a 3-day rate to enable comparisons across colonies over time. I retained daily 
# values for the first two experimental days, as 3-day rate could not be calculated.

pollen_rate <- data.frame(Date = pollen$Date,
                          Day  = pollen$Day,
                          lapply(pollen[colonies], function(x)
                            replace(round(zoo::rollmean(x, k = 3, align = "right", fill = NA), 3),
                                    1:2,
                                    x[1:2])))

# Converting table to long format:
pollen_final <- pollen_rate %>% pivot_longer(cols = C6:C11,
                                             names_to = "Colony",
                                             values_to = "Pollen_consumption_rate")

# Assigning the Cadmium (Cd) concentrations to each colony:
pollen_final$Concentration <- NA
pollen_final$Concentration[pollen_final$Colony %in% c("C6","C14")] <- "0 mg/L"
pollen_final$Concentration[pollen_final$Colony %in% c("C1","C2")] <- "0.02 mg/L"
pollen_final$Concentration[pollen_final$Colony %in% c("C8","C10")] <- "0.2 mg/L"
pollen_final$Concentration[pollen_final$Colony %in% c("C4","C5")] <- "1 mg/L"
pollen_final$Concentration[pollen_final$Colony %in% c("C7","C12")] <- "2 mg/L"
pollen_final$Concentration[pollen_final$Colony %in% c("C9","C11")] <- "20 mg/L"

# Organising the columns in the final dataset:
pollen_final <- pollen_final %>% 
  select(Colony,
         Concentration,
         Date,
         Day,
         Pollen_consumption_rate)

# Calculating overall summary statistics:
pollen_final %>% summarise(n = sum(!is.na(Pollen_consumption_rate)),
                           mean = mean(Pollen_consumption_rate, na.rm = TRUE),
                           sd = sd(Pollen_consumption_rate, na.rm = TRUE),
                           median = median(Pollen_consumption_rate, na.rm = TRUE),
                           min = min(Pollen_consumption_rate, na.rm = TRUE),
                           max = max(Pollen_consumption_rate, na.rm = TRUE),
                           se = sd / sqrt(n))

# Calculating summary statistics per Cd treatment:
pollen_final %>% group_by(Concentration) %>%
  summarise(n = sum(!is.na(Pollen_consumption_rate)),
            mean = mean(Pollen_consumption_rate, na.rm = TRUE),
            sd = sd(Pollen_consumption_rate, na.rm = TRUE),
            median = median(Pollen_consumption_rate, na.rm = TRUE),
            min = min(Pollen_consumption_rate, na.rm = TRUE),
            max = max(Pollen_consumption_rate, na.rm = TRUE),
            se = sd / sqrt(n),
            .groups = "drop")



# 3.1 HYPOTHESIS 3: COLONY RECONSTRUCTION ------------------------------------------------------------

# AIM: Reconstruct the structure of bumblebee colonies by comparing changes in pixel brightness
# across multiple frames extracted from colony video recordings.

# CHALLENGE: bumblebees moving around the colony can mask structures of interest (cocoons). 

# APPROACH: We will use changes in pixel brightness across multiple frames extracted from colony 
# video recordings to distinguish static nest structures (e.g. cocoons, which are generally bright) 
# from mobile elements (bumblebees, which are mainly dark). Identifying which pixels belong to 
# static structures and which belong to moving bees allows us to reconstruct bumblebee colonies 
# by selecting pixels individually from the frames where they are associated with static structures, 
# minimising the influence of moving bumblebees and producing a cleaner reconstruction of the colony.

# Note that this script was used to reconstruct all 12 colonies across seven sampling days,
# with the final sampling day consisting of photographs taken at the end of the experiment.
# However, only six randomly selected reconstructions were ultimately used to test hypothesis 3.

# Loading all bumblebee colony videos subject to analysis. They were stored in individual folders 
# representing each filming day:

day_folder <- "17-06-26"     # Indicating the recording day subject to analysis (eg. 17-06-26).
date_file <- gsub("-", "", day_folder)    # Extracting folder date to name files in the future.

videos <- list.files(day_folder,               # Finding all video files from the selected day.
                     pattern = "\\.MOV$",      # Example: 17-06-26/C1_170626.MOV
                     full.names = TRUE)        # Each folder contains 12 recordings, each belonging 
                                               # to one of the 12 study colonies. 

# At the end of the experiment, colonies were frozen, and bumblebees were removed to take a final 
# photo of the last stage of each colony. This made the nest structures of all colonies completely 
# visible and clear, requiring only a final photo to document their final stage. I account for this 
# situation in this code so the script can work regardless of whether the input is colony videos
# or photographs:

use_photos <- length(videos) == 0

if(use_photos){
  files <- list.files(day_folder,
                      pattern = "\\.png$",
                      full.names = TRUE)
} else {
  files <- videos
}

# Creating a folder where to save all outputs from the analysis.
dir.create("Results", showWarnings = FALSE)

# Image settings for all figures (high resolution to preserve details for future image processing 
# and analysis):
img_width  <- 5000   # pixels
img_height <- 5000   # Pixels
img_res    <- 300    # Resolution of the output image (DPI - Dots Per Inch)

# Automating the colony reconstruction procedure for all colonies in the same recording day using 
# a loop:

for(file in files){    # There are 12 files for each recording day, so the loop will run 12 times.
  
  # ORGANISATION OF THE RESULTS DIRECTORY -------------------------------------------------------------
  
  # Creating a central results directory, with individual subfolders for each study colony.
  
  colony <- sub("_\\d{6}\\.(MOV|png)$",     # Extracting the name of the colony subject to study 
                "",                         # from the name of the files. 
                basename(file))             # Example: 17-06-26/C1_170626.MOV --> C1
  
  # Creating descriptive filenames for the analysis outputs:
  original_filename <- paste0("Original_", colony, "_", date_file)
  stretch_filename  <- paste0("LinearContrastStretch_", colony, "_", date_file)
  
  colony_dir <- file.path("Results",colony)             # Defining where analysis results will be
  day_results_dir <- file.path(colony_dir,day_folder)   # stored and creating subfolders for
  dir.create(day_results_dir, recursive = TRUE)         # the study colonies.
  
  if(use_photos){
    
    clean_colony <- rast(file)     # Read the clean colony image directly.
    
    cat("Loaded clean colony image for", colony, "...\n")   # Printing messages to track code progress.
    
  } else {      # The final photo of the colonies does not need colony reconstruction, as they do not
    # have any moving element masking our structures of interest. Therefore, I create a
    # unique loop for the recording days where colonies were documented through videos:
    
    # EXTRACTION OF FRAMES FROM VIDEOS ------------------------------------------------------------------
    
    cat("Processing:", colony, "\n")    # Printing progress. Which colony is currently being processed?
    
    # Creating a subfolder within each colony directory to store the video frames to be analysed:
    frames_dir <- file.path(day_results_dir, "frames")                 
    dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)  
    
    # Loading the video for the selected colony and recording day + extracting frames at 5 fps 
    # (5 frames per second):
    av_video_images(file,                      
                    destdir = frames_dir,      
                    format = "png",
                    fps = 5)
    
    # Creating a list with all extracted frames:
    frame_files <- list.files(frames_dir,             
                              pattern = "\\.png$",
                              full.names = TRUE)
    
    # Each image is made of pixels organised into rows and columns. Since my goal is to compare the 
    # same pixel positions across all frames and identify which areas change over time (e.g. moving 
    # bees), I need to read the frames and store them as a multi-layer raster object. By stacking all
    # frames together, R can examine each pixel position across frames and assess whether their brightness 
    # remains stable or changes significantly due to bumblebee movement.
    
    frames <- rast(frame_files)   # Reading all extracted frames and saving them in a raster stack. 
    
    # Each image becomes a layer. Because the PNG frames are in RGB colour, each image is not saved
    # as a single layer but three, each one corresponding to one of the 3 RGB different channels (Red, 
    # Green and Blue). Example: Image 1 --> R1 layer, G1 layer, B1 layer (3 layers per frame). 
    
    # Layers are ordered like this: R1, G1, B1, R2, G2, B2, R3, G3, B3... Since my goal is to compare 
    # the same colour channel across time, we need to separate them:
    
    r <- frames[[seq(1, nlyr(frames), 3)]]  # Selecting all layers corresponding to the red channel.
    g <- frames[[seq(2, nlyr(frames), 3)]]  # Selecting all layers corresponding to the green channel.
    b <- frames[[seq(3, nlyr(frames), 3)]]  # Selecting all layers corresponding to the blue channel.
    
    cat("Extracted all frames from", colony, "...\n")   # Printing messages to track code progress.
    
    # COLONY RECONSTRUCTION: 95TH PERCENTILE PROJECTION ----------------------------------------------
    
    # The question now remains: which brightness value do we use as a threshold to disentangle static 
    # nest structures from moving elements? We could choose the maximum intensity value for each pixel, 
    # however, some bees have small patches of bright hairs that can produce high-intensity pixel values. 
    # Choosing the maximum brightness value would add noise to the reconstruction, resulting in a blurry 
    # nest reconstruction because of the highly reflective bumblebee hairs. Therefore, we choose the 
    # brightness value according to the 95th percentile (value below which 95% of observations fall) 
    # to provide a good balance between reducing the influence of bees and preserving the brighter 
    # underlying nest structures, resulting in a cleaner and more accurate reconstruction overall.
    
    # Calculating the 95th percentile brightness value for each pixel in each colour channel separately:
    
    cat("Calculating the 95th percentile projection of", colony, "...\n") 
    
    p95_r <- app(r, function(x) quantile(x, 0.95))      
    p95_g <- app(g, function(x) quantile(x, 0.95))
    p95_b <- app(b, function(x) quantile(x, 0.95))
    
    # Combining the three colour channels into a reconstructed RGB image:
    clean_colony <- c(p95_r, p95_g, p95_b)  
    
  }
  
  # The following code below is applied regardless of whether the input was a photograph (taken on the
  # final day of the experiment) or a video reconstructed using the 95th percentile projection.
  
  # Preparing the file where the reconstructed colony in RGB image will be saved: 
  png(file.path(day_results_dir,
                paste0(original_filename, ".png")),
      width = img_width,
      height = img_height,
      res = img_res)
  
  plotRGB(clean_colony) # Plotting the clean RGB image.
  
  dev.off()
  
  # Prepare the file where the reconstructed colony with a linear contrast stretch will be saved.
  # Linear contrast stretching rescales the brightness of images increasing visual contrast and 
  # making structures easier to distinguish:
  
  png(file.path(day_results_dir,
                paste0(stretch_filename, ".png")),
      width = img_width,
      height = img_height,
      res = img_res)
  
  if(!use_photos) {
    cat("Completed the 95th percentile projection for", colony, ".\n")}
  
  cat("Applying linear contrast stretch to", colony, "...\n")
  
  plotRGB(clean_colony, stretch = "lin")    # Plotting the clean image with a linear contrast stretch.
  
  dev.off()
  
  writeRaster(clean_colony,                                    # Saving the raster object to load it again
              file.path(day_results_dir,                       # in the future easily without having to 
                        paste0(original_filename, ".tif")),    # recalculate the 95th percentile again.
              overwrite = TRUE)
  
  cat("Completed applying the linear contrast stretch for", colony, ".\n")
  cat(colony, "has been successfully processed.\n")
  
}







# 3.2 HYPOTHESIS 3: COMBINING ALL MANUAL ANNOTATIONS INTO ONE SINGLE DATASET -----------------------

# Defining colonies subject to study:

colonies <- c("C6","C14","C1","C2","C8","C10", "C4","C5","C7","C12","C9","C11")

dates <- c("17-06-26", "19-06-26", "22-06-26", "25-06-26", "02-07-26", "09-07-26", "16-07-26")

# Function to import and format csv files containing Fiji measurements:

import_fiji_measurements <- function(colony, date_folder, confidence, object){
  
  date_file <- gsub("-", "", date_folder)     # Date used inside the filename (170626)
  date_label <- gsub("-", "/", date_folder)   # Date displayed in the dataframe (17/06/26)
  
  # All Fiji measurement files are stored in the "Manual annotations" directory:
  file <- file.path("Manual annotations",
                    paste0(confidence, "_", object, "_", colony, "_", date_file, ".csv"))
  
  if (!file.exists(file)) {
    warning("File not found: ", file) 
    return(NULL)
  }
  
  data <- read.csv(file) %>% rename(ID = X.1,
                                    Mean_brightness = Mean,
                                    Min_brightness = Min,
                                    Max_brightness = Max,
                                    SD_brightness = StdDev,
                                    Circularity = `Circ.`) %>%
                             mutate(ID = paste(confidence,
                                               ifelse(object == "Cocoons", "Cocoon", "Honey Pot"),
                                               colony,
                                               date_file,
                                               ID,
                                               sep = "_"),
                                    Image = paste0("CLAHE_", colony, "_", date_file, ".jpeg"),
                                    Range_brightness = Max_brightness - Min_brightness,
                                    CV_brightness = SD_brightness / Mean_brightness,
                                    Aspect_ratio = Major / Minor,
                                    Colony = colony,
                                    Date = date_label,
                                    Group = factor(ifelse(object == "Cocoons", "Cocoon", "Honey Pot"),
                                                   levels = c("Cocoon", "Honey Pot")),
                                    Confidence = factor(ifelse(confidence == "HC", "High confidence", "Low confidence"),
                                                        levels = c("High confidence", "Low confidence"))) %>%
                              select(Image, Colony, Date, ID, Area, Mean_brightness, Min_brightness, 
                                     Max_brightness, SD_brightness, Range_brightness, CV_brightness,
                                     Circularity, Aspect_ratio, X, Y, Group, Confidence)
  return(data)
}


# Reading everything:

objects <- c("Cocoons", "HoneyPots")
confidences <- c("HC", "LC")

all_measurements <- list()

for(date in dates){
  for(colony in colonies){
    for(object in objects){
      for(confidence in confidences){
        
        all_measurements[[length(all_measurements) + 1]] <- import_fiji_measurements(colony,
                                                                                     date,
                                                                                     confidence,
                                                                                     object)
        
      }
    }
  }
}

# Merging everything into one single dataframe:

all_measurements <- bind_rows(all_measurements)

write.csv(all_measurements,
          file = "manual_measurements.csv",
          row.names = FALSE)






# 3.3 HYPOTHESIS 3: VALIDATING A BLOB DETECTOR BASED ON THE DETERMINANT OF THE HESSIAN ALGORITHM -----

# AIM: Design a blob detector algorithm to automate the detection of cocoons as accurate as 
# possible in colony nests whose structures are visible.

# APPROACH: Build a blob detector based on the Determinant of the Hessian (DoH). This image-processing 
# algorithm uses 4 different parameters to identify blob-like structures, which are described below:

# 1. Gaussian blur (sigma; minimum and maximum): controls the amount of image blur applied 
# during blob detection. While excessive blur can hide small structures, a moderate amount 
# of blur reduces image noise and makes rounded objects easier to distinguish. The DoH 
# algorithm applies different levels of Gaussian blur by defining a range between a minimum 
# and a maximum sigma value, smoothing the image by averaging the brightness of neighbouring 
# pixels and allowing blob-like structures of different sizes to be detected more reliably.

# 2. Threshold: represents the minimum score assigned by the algorithm to a pixel location 
# to be considered a potential blob centre. This score is based on how the brightness of 
# each pixel position changes in different perpendicular directions. If brightness changes 
# in all directions, as expected at the centre of a rounded object when moving towards its 
# darker boundaries, the pixel receives a high score and can be identified as the centre of 
# a potential blob.

# Overlap: determines whether overlapping blob detections should be retained or discarded.

# THE GROUND TRUTH: To build a reliable blob detector and evaluate its accuracy in automating cocoon 
# detection, it is important to use a robust set of manual annotations (in this context, 
# cocoons and honey pots) as the ground truth against which to compare the performance of the 
# automated detector. 

# Loading the measurements of manual annotations (i.e cocoons and honey pots) from Fiji:

manual_measurements <- read.csv("manual_measurements.csv")    # This file is obtained from the previous 
                                                              # section of code.

table(manual_measurements$Group)   # Checking number of manual annotated objects measured.

# Calculating sigma for my manual annotations, which maintains a mathematical relationship with the radius 
# and diameter of rounded objects. As sigma increases, small blobs become less detectable, 
# whereas larger blobs remain detectable:
manual <- manual_measurements %>% mutate(Radius = sqrt(Area / pi),
                                         Diameter = Radius * 2,
                                         DoH_sigma = Diameter / (2 * sqrt(2)),   
                                         Group = factor(Group,
                                                        levels = c("Cocoon", "Honey Pot"),
                                                        labels = c("Cocoon", "Honey_Pot")))

# Images that will be used during parameter optimization:
images <- unique(manual$Image)


# BUILDING THE DETECTOR FUNCTION: I first load the Python library scikit-image, which provides several 
# advanced algorithms for detecting image features.

feature <- import("skimage.feature")   # Contains blob and other feature detection algorithms.
cv2 <- import("cv2")                   # Reads images and loads them into Python.

# Specifying the file path to the images where cocoons and honey pots were annotated: 
annotated_images <- path.expand("~/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB/Annotated_Images")

# Building and wrapping the blob detection algorithm in a function for later use:

blob_detector <- function(image,        # Image subject to be analysed. 
                          min_sigma,    # Smallest blob size considered. Small sigma = small blobs.
                          max_sigma,    # Largest blob size considered. Large sigma = large blobs.
                          threshold,    # Blob score. 
                          overlap){     # Overlappig blobs.
  
  # Reading the images as a single grayscale channel (CLAHE images are already grayscale).
  img <- cv2$imread(image, cv2$IMREAD_GRAYSCALE)  
  
  # Many image processing algorithms assume that pixel brightness values are normalized to 0-1:
  img <- img / 255   
  
  # Applying the Determinant of the Hessian (DoH) algorithm to detect blob-like structures:
  blobs <- feature$blob_doh(img,                     
                            min_sigma=min_sigma,
                            max_sigma=max_sigma,
                            threshold=threshold,
                            overlap=overlap)
  
  # Converting Python object into a R data frame:
  blobs <- as.data.frame(py_to_r(blobs))      
  
  # Accounting for situations when the detector does not find any blob:
  if(nrow(blobs) == 0 || ncol(blobs) == 0){blobs <- data.frame(Y = numeric(0),
                                                               X = numeric(0),
                                                               DoH_Sigma = numeric(0))} 
  
  else {colnames(blobs) <- c("Y","X","DoH_Sigma")}
  
  # The function will return a dataframe containing the coordinates of the centre of the blobs that 
  # were detected and the Gaussian blur scale at which each blob was detected:
  blobs
  
} 



# BUILDING THE EVALUATION FUNCTION: How do we evaluate the performance of the blob detector? Using the 
# manual annotations from Fiji as the ground truth and comparing if the centre of each detected blobs 
# falls within the radius of the closest manually annotated object.

# To that end, we build a distance matrix to record and compare the distance between every annotated object 
# and every detected blob to identify matches using the Euclidean distance between their centres. 
# Creating the distance matrix function for later use:

distance_matrix <- function(manual, blobs){
  
  # Extracting the coordinates of the manually annotated objects and the detected blobs.
  d <- as.matrix(dist(rbind(manual[, c("X","Y")],    
                            blobs[, c("X","Y")])))
  
  n_manual <- nrow(manual)
  
  D <- d[1:n_manual, 
         (n_manual + 1):nrow(d),
         drop = FALSE]          # Keeping the result as a matrix, even if it has only one row or column.
  
  return(D)  # Returning the distances between every annotated object and every detected blob.
}


# Now, we build the evaluation function of the blob detector. The idea is to identify the 
# closest detected blob for each manually annotated object and compare the Euclidean distance between 
# their centres. If the centre of a detected blob falls within 50% of the radius of a manual annotation, 
# the detection is considered successful (i.e. a True Positive, TP). I used 50% of the radius to ensure 
# that the matching criterion was sufficiently strict, such that all matches were as centred as possible 
# and provided a reliable estimate of object size for subsequent analyses. The evaluation will
# be built in a way that each detected blob can be matched to only one annotated object, preventing 
# duplication. All manual annotations that are not associated with any automatically detected blob, will 
# be considered as "missed" and classified as a False Negative (FN). Similarly, any unmatched detected 
# blob will be classified as False Positives (FP), as they would not correspond to a real cocoon or
# honey pot.

# Building and wrapping the evaluation of the blob detector in a function for later use:
evaluate_detector <- function(manual, blobs){
  
  # If the detector does not detect any blob, no distances are calculated:
  if(nrow(blobs)==0){
    
    return(list(TP = 0,
                FP = 0,
                FN = nrow(manual),
                Precision = 0,
                Sensitivity = 0,
                F1 = 0,
                detected = rep(FALSE, nrow(manual)),
                matched_blob = rep(NA, nrow(manual)),
                blobs = blobs))}
  
  # Distance matrix:
  D <- distance_matrix(manual, blobs)
  
  # As we match blobs with annotated objects, I will register them as "used" to avoid assigning 
  # them several times to several detections:
  used_blob <- rep(FALSE, nrow(blobs))
  
  # I will also record which blob matches with each annotated object and whether a manual annotation has 
  # been detected or missed:
  matched_blob <- rep(NA_integer_, nrow(manual))
  detected <- rep(FALSE, nrow(manual))
  
  TP <- 0 # This is the counter of true positives starting at 0.
  
  # Loop over every manual annotated object:
  for(i in seq_len(nrow(manual))){
    
    order_blob <- order(D[i,])   # I assess blobs from the closest to the furthest.
    
    for(j in order_blob){
      
      if(used_blob[j]) next  # I first check if the blob has been already matched. 
                             # If so, the blob is skipped and the algorithm moves on to the 
                             # next closest blob. This way, I avoid duplicated matches.
      
      if(D[i,j] <= 0.50 * manual$Radius[i]){   # Once the closest unmatched blob has been found, 
                                               # I check whether the its centre falls within the 
                                               # 50% of the radius of the annotated object.
                                               # Strict without losing many true detections. 
        
        # If the minimum distance falls inside the radius of the manual annotated object, I increase 
        # the number of true detections.
        TP <- TP+1                
        
        detected[i] <- TRUE       # Marking the annotated object as detected.
        
        matched_blob[i] <- j      # Recording which blob matched the annotated object.
        
        used_blob[j] <- TRUE      # The blob is recorded as "used" so it cannot be used for
                                  # another match.
        
        # When a suitable blob has been found, there is no need to examine the remaining blobs 
        # for the current annotated object:
        break                    
        
      }
      
    }
    
  }
  
  # For each combination of parameters I calculate False Negatives: manually annotated objects that 
  # were not matched to any detected blob:
  FN <- nrow(manual) - TP   
  
  # I also calculate False Positives (FP): blobs that were never assigned to a manual annotation.
  FP <- sum(!used_blob)  
  
  # I also calculate the precision of the blob detector: Among all blobs detected, how many 
  # corresponded to real manual annotations?
  Precision <- ifelse(TP + FP == 0, 0,   
                      TP/(TP + FP))   
  
  # I also calculate the sensitivity of the blob detector: Among all real annotated objectss, how many were
  # detected?
  Sensitivity <- ifelse(TP + FN == 0, 0,  
                        TP/(TP + FN))      
  
  # I now compute the F1 score, which combines Precision and Sensitivity into a single metric.
  # It penalizes detectors that perform poorly in either Precision or Sensitivity:
  F1 <- ifelse(Precision + Sensitivity == 0, 0,           
               2 * Precision * Sensitivity / (Precision + Sensitivity))   
  
  # Listing all useful parameters:
  list(TP = TP,  # True Positives.
       FP = FP,  # False Positives.
       FN = FN,  # False Negatives.
       Precision = Precision,  # Precision.
       Sensitivity = Sensitivity,  # Sensitivity.
       F1 = F1,  # F1 score.
       detected = detected,  # Manually annotated objects that were detected.
       matched_blob = matched_blob,   # Blobs that were matched.
       blobs = blobs)  # Original data frame of blobs.
  
}


# Now I build for later use the evaluation function for all annotated images using the same four 
# parameters:

evaluate_all_images <- function(manual,             # All manual annotated objects from all analysed images.
                                annotated_images ,  # Folder containing the images.
                                min_sigma,          # Parameter to test.
                                max_sigma,          # Parameter to test.
                                threshold,          # Parameter to test.
                                overlap){           # Parameter to test.
  
  TP <- FP <- FN <- 0
  
  # Analysing each image independently:
  for(im in unique(manual$Image)){
    
    manual_image <- manual %>% filter(Image == im)
    
    # Running the detector in the current image:
    blobs <- blob_detector(image = file.path(annotated_images, im),
                           min_sigma = min_sigma,
                           max_sigma = max_sigma,
                           threshold = threshold,
                           overlap = overlap)
    
    # Comparing detected blobs with manual annotated objects:
    res <- evaluate_detector(manual_image, blobs)
    
    TP <- TP + res$TP  # Counting True Positives.
    FN <- FN + res$FN  # Counting False Negatives.
    FP <- FP + res$FP  # Counting False Positives.
    
  }
  
  # Calculating overall precision (all images together).
  Precision <- ifelse(TP + FP == 0, 0,   
                      TP/(TP + FP))
  
  # Calculating overall sensitivity (all images together).
  Sensitivity <- ifelse(TP + FN == 0, 0,  
                        TP/(TP + FN))
  
  # Calculating overall F1 score (all images together).
  F1 <- ifelse(Precision + Sensitivity == 0, 0,      
               2 * Precision * Sensitivity / (Precision + Sensitivity))
  
  # Returning the results:
  data.frame(min_sigma,
             max_sigma,
             threshold,
             overlap,
             TP,
             FP,
             FN,
             Precision,
             Sensitivity,
             F1)
  
}


# SEARCHING FOR THE BEST BLOB DETECTOR: The DoH blob detector depends on four parameters: min_sigma, 
# max_sigma, threshold and overlap. Different combinations of these parameters produce detectors with 
# different levels of performance. To identify the best detector, we test multiple parameter combinations 
# and evaluate each one using the F1 score.

# First, we inspect the distribution of DoH sigma values measured for the manual annotated objects to see 
# which ranges for min_sigma and max_sigma we should include in the parameters combinations:

hist(manual$DoH_sigma,
     breaks = 20,
     main = "Distribution of DoH sigma values",
     xlab = "DoH sigma",
     col = "grey80",
     border = "white")

quantile(manual$DoH_sigma,
         probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

# I select the range for min_sigma and max_sigma based on the sigma distribution seen in my annotated 
# objects. Most of them have sigma values between approximately 20 and 40, so the blob detector search
# will focus on combinations of sigma within this interval:

# Total combinations to be tested: 450 (5 x 5 x 6 x 3):

grid <- expand.grid(min_sigma = seq(20, 28, 2),    # Min_sigma values: 20, 22, 24, 26, 28
                    max_sigma = seq(32, 40, 2),    # Max_sigma values: 32, 34, 36, 38, 40
                    threshold = seq(0.005, 0.03, by = 0.005),
                    overlap = c(0.3, 0.4, 0.5))        

results <- data.frame()

for(i in seq_len(nrow(grid))){
  
  cat("Combination", i, "of", nrow(grid), "\n")  # Printing messages to track code progress.
  
  pars <- grid[i,]  # Extracting one parameter combination:
  
  # Evaluating the blob detector using the selected combination of parameters:
  score <- evaluate_all_images(manual,
                               annotated_images,
                               pars$min_sigma,
                               pars$max_sigma,
                               pars$threshold,
                               pars$overlap)
  
  results <- rbind(results,score)  # Storing the results.
  
}

# Choosing the best detector: I select the blob detector that has the highest F1 score, that is, 
# the one with the best precision and sensitivity).

best <- results[which.max(results$F1), ]

print(best)

# Saving every tested parameter combination as a CSV:
write.csv(results, "results.csv", row.names = FALSE)

# As running all 450 parameter combinations is time-consuming and requires high computational power, I 
# provide the complete results here for reference:

results <- read.csv("results.csv")         
best <- results[which.max(results$F1), ]

# APPLICATION OF THE BEST BLOB DETECTOR IN ALL IMAGES SUBJECT TO STUDY: The best blob detector is now 
# applied to all study images, including those without manually annotated objects. For each image, the 
# detector records the coordinates and the sigma of every detected blob and export the results as CSV files.

# Specifying where the images subject to be analysed are stored (the same ones where I first annotated 
# manually cocoons and honey pots):

all_images_folder <- path.expand("~/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB/Annotated_Images")

all_images <- list.files(all_images_folder,         # Listing all JPEG images.
                         pattern = "\\.jpeg$",
                         full.names = FALSE)

# Creating the folder and subfolder where the detected blobs for each image are going to be 
# exported:
output_folder <- "Detected_ROIs"
dir.create(output_folder, showWarnings = FALSE)
dir.create(output_folder, showWarnings = FALSE)    # Creating the folder if it doesn´t exist.

for(im in all_images){
  
  cat("Processing:", im, "\n")      # Printing message to track code progress.
  
  # Applying the best blob detector to the current image using the parameters that achieved the
  # highest F1 score during the optimization:
  blobs <- blob_detector(image = file.path(all_images_folder, im),
                         min_sigma = best$min_sigma,
                         max_sigma = best$max_sigma,
                         threshold = best$threshold,
                         overlap = best$overlap)
  
  # Extracting the coordinates of each detected blob together with the amount of blur at which
  # they were detected (DoH_Sigma):
  rois <- data.frame(X = blobs$X,
                     Y = blobs$Y,
                     DoH_Sigma = blobs$DoH_Sigma)
  
  # Specifying the output filename:
  output_name <- paste0("Detected_ROIs_",
                        tools::file_path_sans_ext(im),
                        ".csv")
  
  # Exporting the detected blobs as a CSV file:
  write.csv(rois,
            file = file.path(output_folder, output_name),
            row.names = FALSE)
  
}

# Now that all automated detections for each colony have been exported, they can be imported
# into Fiji to measure their areas. I developed the macro "Automatic_blob_measurement.ijm" to 
# automatically obtain the measurements of interest from all blobs detected in the six selected CLAHE 
# images. At this point of the script, this macro needs to be run in the software Fiji. 

# These automated measurements will be then paired and compared to their corresponding manual
# measurements to assess the final accuracy of the algorithm. 

# BUILDING THE VALIDATION DATASET: To validate the selected blob detector, I need to link each 
# manually annotated cocoon from the six original images to its corresponding automated detection.
# Here, I match all manually measured objects from the six CLAHE images with their corresponding
# automated blobs and extract the scales at which those blobs were detected.

results <- read.csv("results.csv")    # Reading all results from the 450 parameter combinations. 

best <- results[which.max(results$F1), ]     # Selecting the best detector. 

matched_measurements <- data.frame()

for(im in images){     # This loop will run for each of the six annotated images.
  
  # Keeping only the manual annotations corresponding to the current image:
  manual_image <- manual %>% filter(Image == im)
  
  # Running the best blob detector on the current annotated image:
  blobs <- blob_detector(image = file.path(annotated_images, im),
                         min_sigma = best$min_sigma,
                         max_sigma = best$max_sigma,
                         threshold = best$threshold,
                         overlap = best$overlap)
  
  # Matching each detected blob to the manual annotations:
  res <- evaluate_detector(manual_image, blobs)
  
  # Recording whether each manual annotation was successfully detected and, if so,
  # the corresponding blob ID and the Gaussian scale (sigma) at which it was detected:
  temp <- manual_image %>% mutate(Detected = res$detected,
                                  Blob_ID = res$matched_blob,
                                  Detector_sigma = blobs$DoH_Sigma[Blob_ID])
  
  matched_measurements <- bind_rows(matched_measurements, temp)
}

# I filter the results and focus on the nest structures of interest: the cocoons. 

matched_measurements <- matched_measurements %>% filter(Group == "Cocoon", Detected) %>%
                                                 select(Image,
                                                        Colony,
                                                        Date,
                                                        Group,
                                                        manual_area = Area,
                                                        manual_sigma = DoH_sigma,
                                                        Blob_number = Blob_ID,
                                                        Detector_sigma)

# Exporting the comparison dataset:
write.csv(matched_measurements,
          file = "comparison_dataset.csv",
          row.names = FALSE)


# This dataset will later allow me to link the manual measurements to the automated measurements
# obtained by running the Fiji macro described above.




# 3.4 HYPOTHESIS 3: DATA WRANGLING OF BOTH MANUAL AND AUTOMATED MEASUREMENTS -----------------------------

# Loading the automated blob measurements extracted in Fiji (including blob area) for all
# colonies and sampling days:

automated_measurements <- read.csv("automated_measurements.csv")

# Adding useful columns to bind with the dataset containing the manual measurements:

automated_measurements <- automated_measurements %>%
  mutate(Colony = sub("CLAHE_(C[0-9]+)_.*", "\\1", Image),    
         Date = sub("CLAHE_C[0-9]+_([0-9]{6})", "\\1", Image), 
         Date = format(as.Date(Date, format = "%d%m%y"),
                       format = "%d/%m/%y")) %>%
  rename(detected_area = Area)   


# Loading the file containing manual measurements from the 6 six randomly selected photos: 
comparison_dataset <- read.csv("comparison_dataset.csv")

# Binding both datasets:
validation_dataset <- comparison_dataset %>% 
  left_join(automated_measurements %>%
              select(Colony,
                     Date,
                     Blob_number,
                     detected_area),
            by = c("Colony", "Date", "Blob_number")) %>%
  # Calculating the difference between manual and automated measurements:
  mutate(Error = detected_area - manual_area,
         Absolute_error = abs(Error),
         Percent_error = 100 * Absolute_error / manual_area)










# 4. STATISTICAL MODELLING TO TEST THE THREE HYPOTHESIS ----------------------------------------------

# HYPOTHESIS 1:

model_1 <-lmer(mean_temp ~ Concentration * Day + (1 | Colony),
               data = abiotic_final)

summary(model_1)

# Plotting predicted changes in nest temperature over time across Cd treatments:
pred_model1 <- ggpredict(model_1, terms = c("Day [0:30]", "Concentration"))

p1 <- ggplot(pred_model1,
             aes(x = x,
                 y = predicted,
                 colour = group,
                 group = group,
                 linetype = group)) +
  geom_line(linewidth = 1.25) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                     limits = c(0, 30),
                     expand = c(0, 0)) + 
  scale_y_continuous(limits = c(20.5, 24.5),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("0 mg/L" = "black",
                                 "0.02 mg/L" = "#D55E00",
                                 "0.2 mg/L" = "#009E73",
                                 "1 mg/L" = "#0072B2",
                                 "2 mg/L" = "#CC79A7",
                                 "20 mg/L" = "#F0E450")) +
  scale_linetype_manual(values = c("0 mg/L" = "dotted",
                                   "0.02 mg/L" = "solid",
                                   "0.2 mg/L" = "solid",
                                   "1 mg/L" = "solid",
                                   "2 mg/L" = "solid",
                                   "20 mg/L" = "solid")) +
  labs(title = "Predicted changes in nest temperature over time across Cd treatments",
       x = "Experimental day",
       y = "Predicted mean daily temperature (ºC)",
       colour = "Cd concentration",
       linetype = "Cd concentration") +
  theme_classic() +
  theme(plot.title = element_text(size = 11.5, face = "bold", hjust = 0.9, margin = margin(b = 35)),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10),
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15)) +
  guides(colour = guide_legend(title.position = "top",
                               title.hjust = 0.5,
                               nrow = 1,
                               byrow = TRUE),
         linetype = guide_legend(title.position = "top",
                                 title.hjust = 0.5,
                                 nrow = 1,
                                 byrow = TRUE))

p1

# HYPOTHESIS 2:

model_2 <-lmer(Pollen_consumption_rate ~ Concentration * Day + (1 | Colony),        
               data = pollen_final)

summary(model_2)   # This model assumes that every treatment can fed on pollen at different rates 
                   # over time.

# Plotting predicted changes in pollen consumption rate over time across Cd treatments:
pred_model2 <- ggpredict(model_2,
                         terms = c("Day", "Concentration"))

p2 <- ggplot(pred_model2,
             aes(x = x,
                 y = predicted,
                 colour = group,
                 group = group,
                 linetype = group)) +
  geom_line(linewidth = 1.25) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                     limits = c(0, 30),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 35),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("0 mg/L" = "black",
                                 "0.02 mg/L" = "#D55E00",
                                 "0.2 mg/L" = "#009E73",
                                 "1 mg/L" = "#0072B2",
                                 "2 mg/L" = "#CC79A7",
                                 "20 mg/L" = "#F0E450")) +
  scale_linetype_manual(values = c("0 mg/L" = "dotted",
                                   "0.02 mg/L" = "solid",
                                   "0.2 mg/L" = "solid",
                                   "1 mg/L" = "solid",
                                   "2 mg/L" = "solid",
                                   "20 mg/L"   = "solid")) +
  labs(title = "Predicted changes in pollen consumption rate over time across Cd treatments",
       x = "Experimental day",
       y = "Predicted pollen consumption rate (g/day)",
       colour = "Cd concentration",
       linetype = "Cd concentration") +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 1.2, margin = margin(b = 15)),
        axis.title.x = element_text(, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10),
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15)) +
  guides(colour = guide_legend(title.position = "top",
                               title.hjust = 0.5,
                               nrow = 1,
                               byrow = TRUE),
         fill = "none")

p2

# Merging the two plots from the first two hypotheses together:

p1 <- p1 + theme(axis.title.x = element_text(size = 13),
                 axis.title.y = element_text(size = 13),
                 aspect.ratio = 0.7,
                 plot.margin = margin(5, 25, 5, 5))

p2 <- p2 + theme(axis.title.x = element_text(size = 13),
                 axis.title.y = element_text(size = 13),
                 aspect.ratio = 0.7,
                 plot.margin = margin(5, 5, 5, 25))

combined_plot_2 <- (p1 | p2) + 
                   plot_layout(widths = c(1, 1), guides = "collect") +
                   plot_annotation(tag_levels = "A") &
                   theme(legend.position = "bottom", 
                         legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0),
                         legend.margin = margin(t = 0, r = 0, b = 110, l = 0),
                         plot.tag = element_text(face = "bold", size = 14),
                         plot.tag.position = c(0.02, 1),
                         plot.title = element_blank())

combined_plot_2


# HYPOTHESIS 3:

# I first assess how closely the blob detector reproduces manual measurements by comparing
# automated and manual estimates of cocoon size --> Does the automated method assign larger 
# areas to larger manually annotated cocoons?

cor.test(validation_dataset$manual_area,         # r = 0.64. Larger manually annotated cocoons 
         validation_dataset$detected_area)       # were also estimated as larger by the blob detector.
                                                 # However, the correlation is not very strong, 
                                                 # it is moderate. 

# Are the automated measurements systematically different from the manual ones? I assess 
# whether automated cocoon size estimation is statistically different from the manual ones 
# using a paired t-test, since each automated measurement corresponds to one manual measurement:

t.test(validation_dataset$detected_area,       # The automated estimates differ significantly 
       validation_dataset$manual_area,         # from the manual annotations on average by 
       paired = TRUE)                          # 2001 pixels² (the blob detector overestimates
                                               # cocoon size).

# Checking that the differences between the paired measurements are approximately normally 
# distributed (assumption of paired t-tests).

hist(validation_dataset$Error,
     breaks = 20,
     main = "Distribution of paired differences",
     xlab = "Detected area - Manual area")

# Plotting differences between both methods in cocoon size estimation:

validation_long <- validation_dataset %>% 
  select(manual_area, detected_area) %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Area") %>%
  mutate(Method = recode(Method, manual_area = "Manual annotations",
                         detected_area = "Blob detector"))


ggplot(validation_long,
       aes(x = Method,
           y = Area,
           fill = Method)) +
  geom_boxplot(width = 0.5,
               alpha = 0.6,
               outlier.shape = NA) +   
  geom_jitter(width = 0.12,
              size = 2,
              alpha = 0.35,
              colour = "black") +
  scale_fill_manual(values = c("Manual annotations" = "#0072B2",
                               "Blob detector" = "#D55E00")) +
  theme_classic() +
  labs(x = NULL,
       y = "Cocoon area (pixels²)",
       fill = "Method") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 1, margin = margin(b = 15)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 7)),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15)) 








