# Load required packages
library(knitr)     # For creating tables and reports
library(stats)     # Base R statistical functions
library(MARSS)     # Multivariate time series analysis
library(forecast)  # For ARIMA modeling and forecasting
library(datasets)  # Access to built-in datasets in R
library(tseries)   # Time series analysis and testing

# Load dataset from CSV file and format date column
data = read.csv('data.csv')
data$date = as.Date(data$date, format = "%d/%m/%y")  # Convert date column to Date type
data = data.frame('date' = data$date, 'value' = data$value)  # Create new dataframe with date and value

# Calculate the range of stock prices (used for plotting)
r = range(data$value)[2] - range(data$value)[1]

# Create a time series object for the stock prices
f = 25  # Frequency of the time series (e.g., 25 trading days in a month)
a = 1   # Start of the time series
b = 1   # Arbitrary starting point for the series
value_ts = ts(data$value, frequency = f, start = a)

# Plot the raw time series
plot.ts(value_ts, main = "Time series plot", ylab = 'value', col = '#69b3a2', bty = 'n')

# Decompose the time series into trend, seasonality, and residuals
par(mfrow = c(1, 2))  # Set up plotting space for two plots side by side
value_dec = decompose(value_ts)  # Decompose the time series
plot(value_dec, yax.flip = TRUE, col = '#69b3a2', bty = 'n')  # Plot the decomposition

# Apply Box-Cox transformation to stabilize variance
lambda = BoxCox.lambda(value_ts)  # Find optimal lambda for transformation
log_ts = BoxCox(value_ts, lambda)  # Apply Box-Cox transformation

# Decompose the transformed time series
par(mfrow = c(1, 2))  # Prepare for side-by-side plotting
value_dec = decompose(log_ts)  # Decompose the transformed series
plot(value_dec, yax.flip = TRUE, col = '#69b3a2', bty = 'n')  # Plot the decomposition

# Polynomial fitting to explore the degree of polynomial model
n = 4  # Maximum degree of the polynomial
par(mfrow = c(2, 2))  # 2x2 plot grid
t = 1:length(value_ts)  # Time index for polynomial fitting

for (i in 1:n) {
  # Fit polynomial of degree i
  fit_poln = lm(value_ts ~ poly(t, i))
  
  # Print performance metrics
  print(paste('Degree of polynomial: ', i))
  print(paste('adj R squared: ', summary(fit_poln)$adj.r.squared))
  print(paste('AIC: ', AIC(fit_poln)))
  print(paste('BIC: ', BIC(fit_poln)))
  
  # Plot observed vs fitted values for polynomial model
  plot(value_ts, ylab = 'value', col = '#69b3a2', bty = 'n')
  lines(ts(fit_poln$fitted.values, frequency = f, start = a), type = "l", col = "red", lwd = 1, pch = 10)
}

# Select best polynomial model (degree 3 in this case) and analyze residuals
par(mfrow = c(1, 1))  # Reset plotting space
fit_poln = lm(value_ts ~ poly(t, 3))  # Fit a degree 3 polynomial
value_pol = ts(value_ts - fit_poln$fitted.values, frequency = f, start = a)  # Calculate residuals
plot.ts(value_pol, col = 'orange', main = 'Residuals', ylab = 'res', bty = 'n', ylim = c(-r/2, r/2))  # Plot residuals

# Perform Augmented Dickey-Fuller test for stationarity
adf.test(value_ts)  # Null hypothesis: data is non-stationary (p < 0.05 to reject)

# Perform KPSS test for stationarity
kpss.test(value_ts)  # Null hypothesis: data is stationary (p > 0.05 to fail to reject)

# Differencing to eliminate trend and seasonality
l = 1  # Lag for differencing to remove seasonality
d = ndiffs(value_ts, test = "adf")  # Find number of differences needed to remove trend

# Apply differencing and plot residuals
value_dif = diff(value_ts, lag = l, differences = d)
plot.ts(value_dif, main = 'Residuals', ylab = 'res', col = 'orange', ylim = c(-r/2, r/2), bty = 'n')

# Plot ACF and PACF for identifying ARIMA parameters
par(mfrow = c(1, 2))
plot(acf(value_ts, plot = FALSE), main = 'Auto-Correlation Function Estimation', col = 'orange', bty = 'n')
q = 11  # Tentative MA order based on PACF
plot(pacf(value_ts, plot = FALSE), main = 'Partial Auto-Correlation Function Estimation', col = 'orange', bty = 'n')

# Fit AR(1) model and plot the fitted values and residuals
p = 1  # AR order
par(mfrow = c(1, 2))  # Set up plotting space for two plots
fit_ar = arima(value_ts, order = c(p, 0, 0), method = 'CSS')  # Fit AR(1) model
value_ar = fitted(fit_ar)  # Fitted values from AR model
plot.ts(value_ts, ylab = 'value', main = 'AR(1) model', col = '#69b3a2', bty = 'n')  # Plot original series
lines(value_ar, col = 'red')  # Overlay AR(1) fitted values
legend('bottomright', legend = c('observed', 'fitted'), col = c('#69b3a2', 'red'), pch = c('-', '-'), bty = 'n')

# Plot residuals of AR(1) model
plot.ts(value_ts - value_ar, yax.flip = TRUE, col = 'orange', main = 'Residuals', ylim = c(-r/2, r/2), ylab = 'res', bty = 'n')

# Fit MA(11) model and plot the fitted values and residuals
par(mfrow = c(1, 2))
fit_ma = arima(value_ts, order = c(0, 0, q), method = 'CSS')  # Fit MA(11) model
value_ma = fitted(fit_ma)  # Fitted values from MA model
plot.ts(value_ts, ylab = 'value', main = 'MA(11) model', col = '#69b3a2', bty = 'n')  # Plot original series
lines(value_ma, col = 'red')  # Overlay MA(11) fitted values
legend('bottomright', legend = c('observed', 'fitted'), col = c('#69b3a2', 'red'), pch = c('-', '-'), bty = 'n')

# Plot residuals of MA(11) model
plot.ts(value_ts - value_ma, yax.flip = TRUE, col = 'orange', main = 'Residuals', ylab = 'res', ylim = c(-r/2, r/2), bty = 'n')

# Fit final ARIMA(1,1,11) model
par(mfrow = c(1, 2))  # Set up for two plots
fit_arima = arima(value_ts, order = c(p, d, q), method = 'CSS')  # Fit ARIMA model
value_arima = fitted(fit_arima)  # Fitted values from ARIMA model
plot.ts(value_ts, ylab = 'value', main = 'ARIMA(1,1,11) model', col = '#69b3a2', bty = 'n')  # Plot original series
lines(value_arima, col = 'red')  # Overlay ARIMA fitted values
legend('bottomright', legend = c('observed', 'fitted'), col = c('#69b3a2', 'red'), pch = c('-', '-'), bty = 'n')

# Plot residuals of ARIMA model
plot.ts(value_ts - value_arima, yax.flip = TRUE, col = 'orange', main = 'Residuals', ylab = 'res', ylim = c(-r/2, r/2), bty = 'n')

# Perform diagnostic checks on ARIMA residuals
checkresiduals(fit_arima)  # Ensure residuals are uncorrelated (p > 0.05)

# Forecast future values using ARIMA model
frcast = forecast(fit_arima, h = 7)  # Forecast next 7 periods
kable(frcast)  # Display forecasted values

# Plot forecasted values
plot(frcast, main = "Prediction", ylab = 'value', col = '#69b3a2', bty = 'n')

# Assess accuracy of the fitted ARIMA model
kable(accuracy(frcast$fitted, data$value))  # Compare fitted values to actual values for accuracy metrics
