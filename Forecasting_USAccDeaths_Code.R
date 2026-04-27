# Load library for time series statistical testing
library(tseries)

# Load library for forecasting models and Auto ARIMA
library(forecast)

# Load built-in US accidental deaths dataset
data = USAccDeaths

# View dataset in tabular format
View(data)

# Plot original time series data to observe trend and seasonality
plot(data)

# Check dimensions / total observations of dataset
dim(data)   # dimension of the data

# Apply log10 transformation to stabilize variance
log_data = log10(data)   # Special type of scaling done

# Plot transformed data
plot(log_data)

# Set graph layout to display two plots side by side
par(mfrow = c(1,2))

# Plot Autocorrelation Function to inspect lag relationships
acf(log_data)

# Plot Partial Autocorrelation Function
pacf(log_data)

# Apply first differencing to remove trend
diff1 <- diff(log_data, differences = 1)

# Plot first differenced series
plot.ts(diff1)

# Set graph layout for two plots
par(mfrow = c(1,2))

# Plot ACF of first differenced data
acf(diff1, lag.max = 50)

# Plot PACF of first differenced data
pacf(diff1)

# Apply seasonal differencing with lag 12 (monthly seasonality)
diff2 <- diff(diff1, lag = 12)

# Plot seasonally differenced data
plot(diff2)

# Reset graph layout to single plot
par(mfrow = c(1,1))

# Plot ACF after seasonal differencing
acf(diff2, lag.max = 50)

# Plot PACF after seasonal differencing
pacf(diff2)

# Automatically fit best seasonal ARIMA model
# seasonal = TRUE enables seasonal model search
# stepwise = FALSE performs full model search
auto_model = auto.arima(log_data, seasonal = TRUE, stepwise = FALSE)

# Display detailed summary of selected ARIMA model
summary(auto_model)

# Show residual errors of fitted model
auto_model$residuals

# Forecast next 12 future periods
forecast_value = forecast(auto_model, h = 12)

# Print forecasted values
print(forecast_value)

# Extract non-seasonal MA coefficient (theta)
theta <- coef(auto_model)["ma1"]

# Extract seasonal MA coefficient (Theta)
Theta <- coef(auto_model)["sma1"]

# Print MA coefficient
cat("Non-seasonal MA (theta):", theta, "\n")

# Print seasonal MA coefficient
cat("Seasonal MA (Theta):", Theta, "\n")

# Convert forecast values back to original scale from log scale
original_forecast <- 10^(forecast_value$mean)

# Print actual forecasted deaths count
print(original_forecast)

# Display heading for accuracy metrics
cat("\n--- Corrected Accuracy (Original Scale) ---\n")

# Calculate model accuracy on original scale
accuracy(10^fitted(auto_model), data)
