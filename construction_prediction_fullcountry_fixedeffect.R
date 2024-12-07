memory.limit(size = 16000) # Set memory limit to 16GB 

library("readxl")
library("tidyverse")
RegTable <- read.csv('/Users/lichaohui/Desktop/calculation/construction/RegTable.csv')
head(RegTable)

library(ggplot2)
library(dplyr)
library(openxlsx)

# Check the structure of the data
head(RegTable)

# Split historical and future data
historical_data <- filter(RegTable, yearid <= 28)  # Historical data (1995-2022)
future_data <- filter(RegTable, yearid > 28)  # Future data (2023-2050)

# Linear Regression Model using historical data
# Predict future carbon footprint based on population (x) and regionid as fixed effect
model1 <- lm(y ~ x + factor(regionid), data = historical_data)

# Make predictions for the historical data and future data
historical_data$fitted_y <- predict(model1, newdata = historical_data)
future_data$predicted_y <- predict(model1, newdata = future_data)

# Calculate confidence intervals for the future predictions
prediction_with_interval <- predict(model1, newdata = future_data, interval = "confidence", level = 0.95)
future_data$lower_bound <- prediction_with_interval[, "lwr"]
future_data$upper_bound <- prediction_with_interval[, "upr"]

# Combine the historical data and predicted future data
combined_data <- bind_rows(historical_data, future_data)

# Plotting each region from year 1-56 (1995-2050) with confidence intervals and fitted historical lines
for (region_id in unique(combined_data$regionid)) {
  region_data <- filter(combined_data, regionid == region_id)
  
  # Plot
  p <- ggplot(region_data, aes(x = yearid)) +
    geom_line(aes(y = y), color = "black", na.rm = TRUE, linetype = "dashed") +  # Actual historical carbon footprint
    geom_line(aes(y = fitted_y), color = "red", na.rm = TRUE) +  # Fitted historical carbon footprint
    geom_line(aes(y = predicted_y), color = "blue", na.rm = TRUE) +  # Predicted future carbon footprint
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "lightblue", alpha = 0.5) +  # Confidence intervals for future predictions
    labs(title = paste("Carbon Footprint Prediction for Region", region_id),
         x = "Year", y = "Carbon Footprint") +
    theme_minimal()
  
  print(p)  # Display plot for each region
  
  # Optionally save the plot
  ggsave(paste0("region_", region_id, "_carbon_prediction_with_fitted_CI.png"), plot = p)
}

# Export the predicted future data along with confidence intervals to an Excel file
wb <- createWorkbook()
addWorksheet(wb, "Predicted_Carbon_Data")
writeData(wb, sheet = "Predicted_Carbon_Data", future_data)
saveWorkbook(wb, "Predicted_Carbon_Data.xlsx", overwrite = TRUE)
