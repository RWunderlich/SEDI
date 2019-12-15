# R functions to compute SEDI and it 95% confidence interval
# (Hit Rate, False Alarm Ratio, TSS, Odds Ratio, and ORSS can also be computed)
# by Rainer Wunderlich
# References
# Wunderlich, R. F., Lin, Y. P., Anthony, J., & Petway, J. R. (2019). Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribution models. Nature Conservation, 35, 97-116.

# Original paper introducing SEDI
# Ferro, C.A.T. & Stephensons, D.B. (2011), # Extremal Dependence Indices: Improved Verification Measures for
# Deterministic Forecasts of Rare Binary Events, Weather & Forecasting, 26, 699-713.

# PLEASE NOTE:
# To avoid undefined values, individual zeros are substituted by infinitely small values (1e-9) and
# all returned results indicate whether they represent precise valuse or upper/lower approximations.

# Approximations
minusInf <- function(a = a, b = b, c = c, d = d) {
  return(-1e+09)
}

plusZero <- function(a = a, b = b, c = c, d = d) {
  return(1e-09)
}

# helper
# prevalence or base rate
BaseRate <- function(a = a, b = b, c = c, d = d) {
  return ((a + c) / (a + b + c + d))
}

# helper
# sensitivity or hit rate
HitRate <- function(a = a, b = b, c = c, d = d) {
  return (a / (a + c))
}

# 1 minus hit rate
OneMHR <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return(1 - (a / (a + c)))
}

# log of hit rate
logHR <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(HitRate(a = a, b = b, c = c, d = d)))
}

# helper
# false alarm
FalseAlarm <- function(a = a, b = b, c = c, d = d) {
  return (b / (b + d))
}

# 1 minus false alarm
OneMFA <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return (1 - (b / (b + d)))
}

# log of false alarm
logFA <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(FalseAlarm(a = a, b = b, c = c, d = d)))
}

# not required for SEDI but useful
# odds ratio
# (please cite:
# Stephenson, D.B., 2000. Use of the “odds ratio” for diagnosing forecast skill. Weather and Forecasting, 15(2), pp.221-232.)
OddsRatio <- function(a = a, b = b, c = c, d = d) {
  return (
    (HitRate(a = a, b = b, c = c, d = d) * (OneMFA(a = a, b = b, c = c, d = d)))
    /
      (FalseAlarm(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d)))
  )
}

# not required for SEDI but useful
# odds ratio skill score
# (please cite:
# Stephenson, D.B., 2000. Use of the “odds ratio” for diagnosing forecast skill.Weather and Forecasting, 15(2), pp.221-232.)
ORSS <- function(a = a, b = b, c = c, d = d) {
  return ((HitRate(a = a, b = b, c = c, d = d) - FalseAlarm(a = a, b = b, c = c, d = d))
          /
            (HitRate(a = a, b = b, c = c, d = d) + FalseAlarm(a = a, b = b, c = c, d = d) - 2*(HitRate(a = a, b = b, c = c, d = d)*FalseAlarm(a = a, b = b, c = c, d = d))))
}

# not required but useful
# true skill statistic
# (please cite:
# Allouche, O., Tsoar, A. and Kadmon, R., 2006. Assessing the accuracy of species distribution models:
# prevalence, kappa and the true skill statistic (TSS). Journal of applied ecology, 43(6), pp.1223-1232.)
TSS <- function(a = a, b = b, c = c, d = d) {
  return (HitRate(a = a, b = b, c = c, d = d) - FalseAlarm(a = a, b = b, c = c, d = d))
}

# helper
# denominator of confidence interval formula
CIdeno <- function(a = a, b = b, c = c, d = d) {
  return ( 2 * abs(
    (((OneMHR(a = a, b = b, c = c, d = d)) * (OneMFA(a = a, b = b, c = c, d = d)) + HitRate(a = a, b = b, c = c, d = d) * FalseAlarm(a = a, b = b, c = c, d = d))/((OneMHR(a = a, b = b, c = c, d = d)) * (OneMFA(a = a, b = b, c = c, d = d)))) *
      (log( FalseAlarm(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d)))) + ((2 * HitRate(a = a, b = b, c = c, d = d)) / (OneMHR(a = a, b = b, c = c, d = d))) * (log(HitRate(a = a, b = b, c = c, d = d) * (OneMFA(a = a, b = b, c = c, d = d))))
  )  )
}

# helper
# divisor of confidence interval formula
CIdivi <- function(a = a, b = b, c = c, d = d) {
  return (HitRate(a = a, b = b, c = c, d = d) * (log(FalseAlarm(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d)))   +  log(HitRate(a = a, b = b, c = c, d = d) * (OneMFA(a = a, b = b, c = c, d = d))))^2)
}

# helper
# trailing sq term of confidence interval formula
CIsqtail <- function(a = a, b = b, c = c, d = d) {
  return (sqrt((HitRate(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d))) / (a + c)))
}

# main function
# symmetric extremal dependence index
sedi <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(x = (max(a,b,c,d) <= 1) && sum(a,b,c,d) <= 1)) { # no value can be computed
    return(list(" (Almost) empty confusion matrix", NaN))  
  }
  
  if (isTRUE(x = (min(a,b,c,d) < 1))) { # catch all matrices which contain at least one zero
    if (isTRUE(x = (sum(a,b) < 1) &&  (sum(c,d) > 1) )) {
      return(list("Always predicting absence:", 0))
    }
    if (isTRUE(x = (sum(c,d) < 1) &&  (sum(a,b) > 1) )) {
      return(list("Always predicting presence:", 0))
    }
    if (isTRUE(x = (sum(b,c) < 1) &&  (sum(a,d) > 1) )) {
      return(list("Perfect prediction:", 1))
    }
    if (isTRUE(x = (sum(a,d) < 1) &&  (sum(b,c) > 1) )) {
      return(list("Opposite prediction:", 0))
    }
  
    if (isTRUE(x = (a < 1 && min(b,c,d) >=1))) { # only true presences equal to zero
      return(list("Undefined, but smaller than",
                  (logFA(a = a+1e-09, b = b, c = c, d = d) - logHR(a = a+1e-09, b = b, c = c, d = d) - log(OneMFA(a = a+1e-09, b = b, c = c, d = d)) + log(OneMHR(a = a+1e-09, b = b, c = c, d = d)))
                  /
                  (logFA(a = a+1e-09, b = b, c = c, d = d) + logHR(a = a+1e-09, b = b, c = c, d = d) + log(OneMFA(a = a+1e-09, b = b, c = c, d = d)) + log(OneMHR(a = a+1e-09, b = b, c = c, d = d)))))  
    }
  
    if (isTRUE(x = (d < 1 && min(a,b,c) >=1))) { # only true absences (also pseudo-absences, incl. background) equal to zero
      return(list("Undefined, but smaller than",
                  (logFA(a = a, b = b, c = c, d = d+1e-09) - logHR(a = a, b = b, c = c, d = d+1e-09) - log(OneMFA(a = a, b = b, c = c, d = d+1e-09)) + log(OneMHR(a = a, b = b, c = c, d = d+1e-09)))
                  /
                  (logFA(a = a, b = b, c = c, d = d+1e-09) + logHR(a = a, b = b, c = c, d = d+1e-09) + log(OneMFA(a = a, b = b, c = c, d = d+1e-09)) + log(OneMHR(a = a, b = b, c = c, d = d+1e-09)))))  
    }
  
    if (isTRUE(x = (b < 1 && min(a,c,d) >=1))) { # only zero commission errors
      return(list("Undefined, but greater than",
                  (logFA(a = a, b = b+1e-09, c = c, d = d) - logHR(a = a, b = b+1e-09, c = c, d = d) - log(OneMFA(a = a, b = b+1e-09, c = c, d = d)) + log(OneMHR(a = a, b = b+1e-09, c = c, d = d)))
                  /
                  (logFA(a = a, b = b+1e-09, c = c, d = d) + logHR(a = a, b = b+1e-09, c = c, d = d) + log(OneMFA(a = a, b = b+1e-09, c = c, d = d)) + log(OneMHR(a = a, b = b+1e-09, c = c, d = d)))))  
    }
  
    if (isTRUE(x = (c < 1 && min(a,b,d) >=1))) { # only zero omission errors
      return(list("Undefined, but greater than",
                  (logFA(a = a, b = b, c = c+1e-09, d = d) - logHR(a = a, b = b, c = c+1e-09, d = d) - log(OneMFA(a = a, b = b, c = c+1e-09, d = d)) + log(OneMHR(a = a, b = b, c = c+1e-09, d = d)))
                  /
                  (logFA(a = a, b = b, c = c+1e-09, d = d) + logHR(a = a, b = b, c = c+1e-09, d = d) + log(OneMFA(a = a, b = b, c = c+1e-09, d = d)) + log(OneMHR(a = a, b = b, c = c+1e-09, d = d)))))  
    }
  }
  else { # all regular cases with reasonable values across the confusion matrix
    return (list("SEDI equal to",
                 (logFA(a = a, b = b, c = c, d = d) - logHR(a = a, b = b, c = c, d = d) - log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
                 /
                 (logFA(a = a, b = b, c = c, d = d) + logHR(a = a, b = b, c = c, d = d) + log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
            ))  
  }
}
# confidence Interval 1
# standard error for SEDI
sediSE <- function(a = a, b = b, c = c, d = d) {
  return ((CIdeno(a = a, b = b, c = c, d = d) / CIdivi(a = a, b = b, c = c, d = d)) * CIsqtail(a = a, b = b, c = c, d = d) * 1.95996)
}

# confidence Interval 2
# return sedi and its 95% confidence interval
sediCI <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(min(a,b,c,d) >= 1))  {
  return (list("SEDI_lcl, SEDI, SEDI_ucl",sedi(a = a, b = b, c = c, d = d)[[2]] - sediSE(a = a, b = b, c = c, d = d),
            sedi(a = a, b = b, c = c, d = d)[[2]],
            sedi(a = a, b = b, c = c, d = d)[[2]] + sediSE(a = a, b = b, c = c, d = d)))
  }
  else{
    return(list("Confidence interval cannot be determined", sedi(a = a, b = b, c = c, d = d)))
  }
}
