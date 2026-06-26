source("cv-hsbvar-nll.R")
source("hsbvar-bea.R")

# ---- Load and align S&P and BTC on common trading dates ----

sp_raw <- read.csv("data/sp.csv")
sp_raw$Date <- as.Date(sp_raw$Date, "%m/%d/%Y")
sp_raw <- sp_raw[order(sp_raw$Date), ]
sp_raw <- sp_raw[!duplicated(sp_raw$Date), ]
stopifnot(all(sp_raw$Value > 0))

btc_raw <- read.csv("data/btc_1d_data_2018_to_2025.csv")
btc_raw$Date <- as.Date(substr(btc_raw$Open.time, 1, 10))
btc_raw <- btc_raw[order(btc_raw$Date), ]
btc_raw <- btc_raw[!duplicated(btc_raw$Date), ]
stopifnot(all(btc_raw$Close > 0))

# Log returns for each series
sp_ret <- data.frame(
  Date = sp_raw$Date[-1],
  r_sp = diff(log(sp_raw$Value))
)
btc_ret <- data.frame(
  Date  = btc_raw$Date[-1],
  r_btc = diff(log(btc_raw$Close))
)

# Inner join on date so both series cover the same days
merged <- merge(sp_ret, btc_ret, by = "Date")
merged <- merged[order(merged$Date), ]
dates <- merged$Date

# ---- Step 2: Remove conditional mean (demeaning) ----
r_sp <- merged$r_sp - mean(merged$r_sp)
r_btc <- merged$r_btc - mean(merged$r_btc)

cat(sprintf(
  "Series length after alignment: %d observations (%s to %s)\n",
  nrow(merged), min(dates), max(dates)
))

# ---- Step 3: Flag outliers (do not remove) ----
flag_outliers <- function(r, label, dates) {
  mad_s <- median(abs(r - median(r))) / 0.6745
  idx <- which(abs(r) > 10 * mad_s)
  if (length(idx) > 0) {
    cat(sprintf("Flagged outlier candidates in %s (>10 robust SDs):\n", label))
    print(data.frame(date = dates[idx], r = r[idx]))
  } else {
    cat(sprintf("No outlier candidates in %s.\n", label))
  }
}
flag_outliers(r_sp, "S&P 500", dates)
flag_outliers(r_btc, "BTC", dates)

# ---- Build n x 2 return matrix ----
Y <- cbind(r_sp, r_btc)
n <- nrow(Y)
plag <- 1L

# ---- CV to select lambda ----
# lambda_grid <- seq(0.01, 0.5, by = 0.01)
# cat("Running CV (H-SBVAR NLL) to select lambda...\n")
# cv_fit <- cv_hsbvar_nll(Y, p = plag, lambda = lambda_grid,
#                         val_spacing = plag + 1L, verbose = FALSE,
#                         max_iter = 5000)
# best_lambda <- cv_fit$best$lambda
# cat(sprintf("CV selected lambda = %.4f\n", best_lambda))

# ---- Stage 1: H-SBVAR with CV-selected lambda ----
cat("Fitting H-SBVAR (stage 1)...\n")
fit <- hsbvar(Y, p = plag, lambda = 0.01, max_iter = 5000, verbose = FALSE)
cp_s1 <- fit$cp
cat(sprintf("Stage 1 changepoints: %d\n", length(cp_s1)))
cat("Stage 1 dates: ")
print(dates[cp_s1])

# ---- Stage 2: BEA pruning ----
cat("Running BEA pruning (stage 2)...\n")
d <- ncol(Y)
bea_fit <- hsbvar_bea(fit,
  Y = Y, p = plag,
  omega_n = (d * d + d * plag) * log(n) / 2
)
cp_s2 <- bea_fit$cp
cp_dates_s1 <- dates[cp_s1]
cp_dates <- dates[cp_s2]
cat(sprintf("Stage 2 (BEA) changepoints: %d\n", length(cp_s2)))
cat("Stage 2 dates: ")
print(cp_dates)

# ---- Plots ----
years <- seq(
  as.integer(format(min(dates), "%Y")),
  as.integer(format(max(dates), "%Y"))
)
year_ticks <- as.Date(paste0(years, "-01-01"))
add_year_axis <- function() axis.Date(1, at = year_ticks, format = "%Y")

pdf("data/sp-btc_hsbvar.pdf", width = 10, height = 13)
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

# plot_cv_hsbvar_nll(cv_fit)

plot(dates, r_sp,
  type = "l", col = "steelblue", lwd = 0.6,
  xaxt = "n", xlab = "", ylab = "Log return",
  main = "S&P 500 daily log returns"
)
add_year_axis()
abline(v = cp_dates_s1, lty = 2, col = "gray60")
abline(v = cp_dates, lty = 2, col = "red")

plot(dates, r_btc,
  type = "l", col = "darkorange", lwd = 0.6,
  xaxt = "n", xlab = "", ylab = "Log return",
  main = "BTC/USD daily log returns"
)
add_year_axis()
abline(v = cp_dates_s1, lty = 2, col = "gray60")
abline(v = cp_dates, lty = 2, col = "red")

dev.off()
cat("Plots saved to data/sp-btc_hsbvar.pdf\n")
