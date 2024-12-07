memory.limit(size = 16000) # Set memory limit to 16GB 


library("readxl")
library("tidyverse")
RegTable <- read.csv('/Users/lichaohui/Desktop/calculation/construction/RegTable.csv')
head(RegTable)

# Load necessary libraries
library(readxl)
library(tidyverse)
library(zoo)

# Ensure yearid and regionid are treated as factors
RegTable$yearid <- as.factor(RegTable$yearid)
RegTable$regionid <- as.factor(RegTable$regionid)

# Loop through each region and perform independent OLS regression
unique_regions <- unique(RegTable$regionid)

# Create an empty list to store data for all regions
all_regions_predictions <- list()

for (region in unique_regions) {
  
  # Filter data for the current region
  region_data <- RegTable %>% filter(regionid == region)
  
  # Handle missing values by linear interpolation for population
  region_data$x <- na.approx(region_data$x, na.rm = FALSE)
  
  # Linear OLS regression for the current region (carbon footprint ~ population)
  model <- lm(y ~ x, data = region_data)
  
  # Predict future carbon emissions and confidence intervals for future years
  # Future years are years 29-56
  future_data <- region_data %>% filter(as.numeric(as.character(yearid)) > 28)
  future_data$predicted_carbon <- predict(model, newdata = future_data)
  
  # 95% prediction intervals for future data
  prediction_with_interval <- predict(model, newdata = future_data, interval = "prediction", level = 0.95)
  future_data$lower_bound <- prediction_with_interval[, "lwr"]
  future_data$upper_bound <- prediction_with_interval[, "upr"]
  
  # Add predictions and confidence intervals back into the original dataset for historical + future data
  region_data$predicted_carbon <- predict(model, newdata = region_data)
  region_data$lower_bound <- NA  # Set NA for confidence intervals for historical data
  region_data$upper_bound <- NA  # Set NA for confidence intervals for historical data
  
  # Combine historical and future data for the region
  all_data_region <- bind_rows(region_data, future_data)
  
  # Save results for each region in the list
  all_regions_predictions[[region]] <- all_data_region
  
  # Plot the actual, fitted, and predicted values for each region
  plot <- ggplot(all_data_region, aes(x = as.numeric(as.character(yearid)))) +
    geom_ribbon(data = future_data, aes(ymin = lower_bound, ymax = upper_bound), fill = "blue", alpha = 0.2) +  # Prediction interval for future data
    geom_line(aes(y = predicted_carbon), color = "red", linetype = "solid") +  # Predicted carbon footprint (both historical and future)
    geom_line(aes(y = y), color = "black", linetype = "dashed") +  # Actual historical carbon footprint as lines
    labs(title = paste("Carbon Footprint Prediction for Region", region),
         x = "Year",
         y = "Carbon Footprint") +
    theme_minimal()
  
  # Display the plot
  print(plot)
  
# Combine all regions' data for further analysis
combined_data <- bind_rows(all_regions_predictions)
