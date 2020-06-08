par_rate_index <- function(par_rate_hist, index_ini = 1000, target_mat = 1, bond_freq = 2, base = 365, cc_mat = 0.9){

  rate_dates <- index(par_rate_hist)
  ini_date <- rate_dates[1]
  last_date <- tail(rate_dates, 1)

  per <- 1/bond_freq
  if (target_mat<cc_mat){
    cpn_days <- round(target_mat*365)
  }else{
    cpn_days <- round(seq(per, target_mat, per) * 365)
  }
  cpn_times <- cpn_days/base
  ldates <- length(rate_dates)
  index_val <- rep(0, ldates)
  index_val[1] <- index_ini

  par_rate <- as.numeric(par_rate_hist[1])
  mac_dur <- rep(0, ldates)

  for(i in 2:ldates){
    delta_days <- as.numeric(rate_dates[i] - rate_dates[i-1])
    cpn_daysi <- cpn_days - delta_days

    cpn_timesi <- cpn_daysi/base
    par_ratei <- as.numeric(par_rate_hist[i])
    lcpn <- length(cpn_times)
    if (lcpn>1){
      cf <- rep(par_rate, lcpn) * c(cpn_times[1], diff(cpn_times))
      cf[lcpn] <- cf[lcpn] + 1
      df0 <- 1/((1 + par_rate) ** cpn_times)
      bond_pr0 <- sum(cf * df0)
      df <- 1/((1 + par_ratei) ** cpn_timesi)
      bond_pr <- sum(cf * df)
      bond_ret <- bond_pr/bond_pr0
      index_val[i] <- index_val[i-1] * bond_ret
    }else{
      cf <- 1
      df <- bond_pr <- 1/((1 + par_ratei) ** cpn_timesi)
      cc_ret <- ((1 + par_rate) ** cpn_times)/((1 + par_ratei) ** cpn_timesi)
      index_val[i] <- index_val[i-1] * cc_ret
    }


    tdf <- cpn_timesi * df
    mac_dur[i] <- round(-sum(cf * tdf)/bond_pr, 3)
    par_rate <- par_ratei
  }
  index_series <- xts(index_val, order.by = rate_dates)
  mac_dur <- xts(mac_dur, order.by = rate_dates)
  return(list(index_series=index_series, mac_dur=mac_dur[-1]))
}
