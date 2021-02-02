
#' Builds index given a fixed income term structure
#'
#' @param cc_hist fixed income curve history
#' @param index_ini index base value
#' @param dates date range for which the calculation will be executed
#' @param target_matur_in_days synthetic bond maturity in days
#' @param bond_freq coupon payment frequency (times per year)
#' @param base n_days in year
#' @param slippage
#' @return list 
#' @export

cc_curve_index <- function(cc_hist, index_ini = 1000, dates = NULL, target_matur_in_days = 1095, 
                           bond_freq = 2, base = 365, slippage = 5, uniroot_interval = c(0, 0.3)){

  cc_curve <- na.omit(cc_hist[, as.numeric(colnames(cc_hist)) <= target_matur_in_days]) # omits NAs and slices columns with maturity less than target
  if(!is.null(dates)){ # If dates is passed, dataframe is sliced with passed dates
    cc_curve <- cc_curve[paste(dates, collapse = '/')]
  }
  cc_dates <- index(cc_curve) # Date vector
  ini_date <- cc_dates[1] # First date in date vector
  last_date <- tail(cc_dates, 1) # Last date in date vector
  target_mat <- target_matur_in_days/base # Target maturity as a proportion of number of days in year (base)

  per <- 1/bond_freq # Coupon payment period per year
  cpn_days <- round(seq(per, target_mat, per) * 365) # Days in which a coupon will be paid 
  cpn_times <- cpn_days/base # Days in which a coupon will be paid as a proportion of 1
  ldates <- length(cc_dates) # Number of dates in date vector
  days <- as.numeric(colnames(cc_curve)) # Maturity of cc_hist columns in days
  index_val <- rep(0, ldates) # Vector of length ldates which will be populated with index values
  index_val[1] <- index_ini # Assign index_ini to first position in index_val

  ytm <- function(x, cpn_times, cpn_rates){ # 
    df <- 1/((1 + cpn_rates) ** cpn_times) # Discount factor - effective rate
    lcpn <- length(cpn_times)
    cf <- rep(x, lcpn) * c(cpn_times[1], diff(cpn_times))
    cf[lcpn] <- cf[lcpn] + 1
    obj <- sum(cf * df) - 1
    return(obj)
  }

  rates <- as.vector(cc_curve[1])/100 # Takes values of first line in cc_curve as vector / 100
  cpn_rates <- approxExtrap(x = days, y = rates, xout = cpn_days)$y # Extrapolates rates to all coupon payment maturities
  par_rate <- uniroot(ytm, interval = uniroot_interval, cpn_times, cpn_rates)$root # Estimates the coupon rate for which a synthetic bond built with all cpn_rates as yields of cash flow payments will be priced at par
  mac_dur <- rep(0, length(cc_dates)) # Macdur vector with length cc_dates

  for(i in 2:ldates){
    delta_days <- as.numeric(cc_dates[i] - cc_dates[i-1]) # Days difference between date which is being populated and previous
    cpn_daysi <- cpn_days - delta_days # Days until coupon payments of i

    rates <- as.vector(cc_curve[i])/100 # Rates of the ith row on cc_curve
    cpn_ratesi <- approxExtrap (x = days, y = rates, xout = cpn_daysi)$y # Rate interpolation for cpn_daysi days to next coupon payment
    cpn_timesi <- cpn_daysi/base # Days until next coupon payment divided by base (multiplier)
    lcpn <- length(cpn_times) # Vector of length cpn_times
    cf <- rep(par_rate, lcpn) * c(cpn_times[1], diff(cpn_times)) # ¿?
    cf[lcpn] <- cf[lcpn] + 1 
    df0 <- 1/((1 + cpn_rates) ** cpn_times) # Discount factors for all coupon payments at i - 1
    bond_pr0 <- sum(cf * df0) # Bond's price at i - 1

    df <- 1/((1 + cpn_ratesi) ** cpn_timesi) # Discount factors for all coupon payments at i
    bond_pr <- sum(cf * df) # Bond's price at i
    bond_ret <- bond_pr/bond_pr0 # Bond's return i / i - 1

    bond_exec_pr <- bond_pr + bond_pr * slippage/10000 # Bond's price adjusted for slippage (for transaction costs)
    tc <- index_val[i-1] * (bond_exec_pr - bond_pr) # Transaction costs
    index_val[i] <- index_val[i-1] * bond_ret - tc # Previous index value multiplied by the return and adjusted by transaction costs

    tdf <- cpn_timesi * df # Discount factor adjusted by cpn_timesi period
    mac_dur[i] <- round(-sum(cf * tdf)/bond_pr, 3) # Macaulay Duration for time i

    cpn_rates <- approxExtrap(x = days, y = rates, xout = cpn_days)$y
    par_rate <- uniroot(ytm, interval = uniroot_interval, cpn_times, cpn_rates)$root
  }
  index_series <- xts(index_val, order.by = index(cc_curve))
  mac_dur <- xts(mac_dur, order.by = index(cc_curve))
  return(list(index_series=index_series, mac_dur=mac_dur[-1]))
}
