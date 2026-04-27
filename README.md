# Time Series Forecasting using Regression Based Methods & Auto ARIMA (USAccDeaths Dataset)

## 📌 Project Overview
This project demonstrates **Time Series Forecasting using Regression Based Methods** and **Auto ARIMA Modeling** using the built-in **USAccDeaths** dataset in **R Programming**. The dataset contains monthly accidental deaths in the United States from **1973 to 1978**. It is widely used for understanding seasonal time series forecasting.

The project covers data visualization, log transformation, trend removal, seasonal differencing, ACF/PACF diagnostics, Auto ARIMA model selection, parameter extraction, forecasting future values, and model accuracy evaluation.

## 👨‍💼 Author
**Biswaditya07 (Biswaditya Saha)**  
MBA in Business Analytics

## 🛠️ Tools & Technologies Used
- R Programming  
- Built-in R Dataset: `USAccDeaths`  
- Libraries Used:
  - `tseries`
  - `forecast`

## 📂 Files Included
| File Name | Description |
|-----------|-------------|
| `USAccDeaths Forecasting Script` | R script for complete forecasting analysis using USAccDeaths dataset |

## 📊 Dataset Information
- Dataset Name: `USAccDeaths`
- Frequency: Monthly
- Period: 1973 - 1978
- Total Observations: 72
- Data Type: Time Series
- Category: Monthly accidental deaths in the United States

## 📈 Project Workflow

### 1️⃣ Load Dataset
The built-in dataset is loaded for analysis.

```r
data = USAccDeaths
View(data)
plot(data)
dim(data)
````

Initial visualization helps identify trend and seasonality.

### 2️⃣ Log Transformation

Variance is stabilized using logarithmic scaling.

```r id="46wlh1"
log_data = log10(data)
plot(log_data)
```

This reduces fluctuations and prepares the data for ARIMA modeling.

### 3️⃣ ACF & PACF Analysis

Autocorrelation and Partial Autocorrelation plots help inspect lag relationships.

```r id="22rr9x"
acf(log_data)
pacf(log_data)
```

These plots indicate the need for differencing.

### 4️⃣ First Differencing

Trend is removed using first-order differencing.

```r id="nq1xvh"
diff1 <- diff(log_data, differences = 1)
plot.ts(diff1)
```

### 5️⃣ Seasonal Differencing

Monthly seasonality is removed using lag 12 differencing.

```r id="j80k34"
diff2 <- diff(diff1, lag = 12)
plot(diff2)
```

This removes repeating yearly seasonal patterns.

### 6️⃣ Recheck ACF & PACF

After differencing, ACF and PACF are re-evaluated.

```r id="jlwmv2"
acf(diff2, lag.max = 50)
pacf(diff2)
```

Used for confirming stationarity and model structure.

### 7️⃣ Auto ARIMA Model Selection

Best seasonal ARIMA model is selected automatically.

```r id="g4j0va"
auto_model = auto.arima(log_data, seasonal = TRUE, stepwise = FALSE)
summary(auto_model)
```

* `seasonal = TRUE` enables seasonal model search
* `stepwise = FALSE` performs complete model optimization

### 8️⃣ Forecast Future Values

```r id="myagx0"
forecast_value = forecast(auto_model, h = 12)
print(forecast_value)
```

Forecasts the next **12 months** of accidental deaths.

### 9️⃣ Extract Model Parameters

```r id="6wdrza"
theta <- coef(auto_model)["ma1"]
Theta <- coef(auto_model)["sma1"]
```

* `theta` = Non-seasonal Moving Average coefficient
* `Theta` = Seasonal Moving Average coefficient

### 🔟 Convert Forecast to Original Scale

```r id="3x70qo"
original_forecast <- 10^(forecast_value$mean)
print(original_forecast)
```

Since log transformation was applied, forecasts are converted back to original values.

### 1️⃣1️⃣ Accuracy Evaluation

```r id="d2z5sn"
accuracy(10^fitted(auto_model), data)
```

Measures model performance using metrics such as:

* RMSE
* MAE
* MAPE

## 📊 Business Applications

* Public Safety Forecasting
* Accident Risk Trend Analysis
* Healthcare Resource Planning
* Emergency Response Preparation
* Seasonal Risk Management
* Government Policy Planning
* Time Series Demand Forecasting

## 🔍 Key Learnings

* Seasonal Time Series Analysis
* Variance Stabilization using Log Transform
* Trend & Seasonal Differencing
* ACF/PACF Interpretation
* Auto ARIMA Model Selection
* Forecasting Future Values
* Coefficient Interpretation
* Accuracy Measurement in R

## 💻 Sample Code

```r id="4w1m1f"
library(tseries)
library(forecast)

data = USAccDeaths
plot(data)

log_data = log10(data)

auto_model = auto.arima(log_data, seasonal = TRUE)

forecast(auto_model, h = 12)
```

## 🚀 Conclusion

This project provides complete practical exposure to **seasonal time series forecasting** using the USAccDeaths dataset. It demonstrates how raw monthly data can be transformed, modeled, forecasted, and evaluated using Auto ARIMA techniques for real-world analytical decision-making.

## ⭐ Support

If you like this project, give it a **star ⭐** on GitHub and connect for more analytics projects.

```
```
