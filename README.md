# Technical Olympic S.A. Stock Price Forecasting Project

## Overview

This project focuses on a comprehensive time series analysis and forecasting of the closing stock prices of **Technical Olympic S.A.** Using historical data, the analysis involves:

- Time series decomposition to uncover underlying trends, seasonal patterns, and residuals.
- Polynomial regression models to explore the fit and select the best degree of polynomial.
- ARIMA modeling, following stationarity tests and parameter selection through the Box-Jenkins method.
- Model diagnostics and residual checks to validate the goodness of fit.
- Forecasting future stock prices based on the best-performing ARIMA model.

The objective is to gain insights into stock price movements and forecast future closing values.

## Table of Contents

- [Project Structure](#project-structure)
- [Data](#data)
- [Dependencies](#dependencies)
- [Methodology](#methodology)
  - [1. Data Loading and Preprocessing](#1-data-loading-and-preprocessing)
  - [2. Time Series Decomposition](#2-time-series-decomposition)
  - [3. Box-Cox Transformation](#3-box-cox-transformation)
  - [4. Polynomial Fitting](#4-polynomial-fitting)
  - [5. ARIMA Modeling](#5-arima-modeling)
  - [6. Forecasting](#6-forecasting)
- [Results](#results)
- [Usage](#usage)
- [Future Work](#future-work)
- [License](#license)

## Project Structure

The repository contains the following:

- `data/`: Directory for the historical stock data.
- `code.R`: Main R script for analysis and forecasting.
- `README.md`: This file, providing an overview of the project.
- `plots/`: Directory to store generated plots and visualizations.
- `forecast_results.csv`: File containing the predicted values.

## Data

The dataset used in this project represents the historical stock closing prices of **Technical Olympic S.A.**, including the corresponding dates.

- **Source**: [Data.csv](data/data.csv)
- **Columns**: 
  - `date`: Date of the stock price record.
  - `value`: Closing stock price on the given date.

## Dependencies

The project uses several R libraries for time series analysis and visualization. Install the required packages as follows:

```R
install.packages(c("knitr", "stats", "MARSS", "forecast", "datasets", "tseries"))
```

## Methodology

### 1. Data Loading and Preprocessing

The dataset is loaded into R, and the `date` column is converted into `Date` format. A `ts` object is created from the stock prices for time series analysis.

```R
data = read.csv('data/data.csv')
data$date = as.Date(data$date, format="%d/%m/%y")
value_ts = ts(data$value, frequency = 25, start = c(1, 1))
```

### 2. Time Series Decomposition

The time series is decomposed into three components:
- **Trend**
- **Seasonality**
- **Residuals**

```R
value_dec = decompose(value_ts)
plot(value_dec)
```

### 3. Box-Cox Transformation

A Box-Cox transformation is applied to stabilize the variance in the data:

```R
lambda = BoxCox.lambda(value_ts)
log_ts = BoxCox(value_ts, lambda)
```

### 4. Polynomial Fitting

We explored different polynomial degrees (up to 4) to model the stock prices. The best polynomial degree was selected based on **Adjusted R-squared**, **AIC**, and **BIC** metrics.

```R
fit_poln = lm(value_ts ~ poly(t, 3))
```

### 5. ARIMA Modeling

After making the time series stationary, ARIMA models were fitted. The process involved:
- **ADF** and **KPSS** tests to check stationarity.
- Auto-correlation (ACF) and partial auto-correlation (PACF) plots for identifying ARIMA parameters.
- Diagnostic checks for residuals.

```R
fit_arima = arima(value_ts, order = c(1, 1, 11), method = 'CSS')
checkresiduals(fit_arima)
```

### 6. Forecasting

Using the best ARIMA model, forecasts were generated for the next 7 periods.

```R
frcast = forecast(fit_arima, h = 7)
plot(frcast)
```

## Results

The model accurately captured the underlying patterns in the stock prices. Diagnostic checks confirmed the adequacy of the ARIMA(1,1,11) model, and forecasts for future values were provided. 

Key Results:
- **ADF test** indicated stationarity after differencing.
- **ARIMA(1,1,11)** model was selected based on AIC, BIC, and residual diagnostics.
- Forecasts for the next 7 periods were generated, showing potential future stock prices.

## Usage

To run this project on your local machine:
1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/stock-forecasting.git
    cd stock-forecasting
    ```
2. Open the `code.R` file and run the analysis in R.

- **Model Comparison**: Evaluate performance across different models.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

[notebook](notebook.html)
