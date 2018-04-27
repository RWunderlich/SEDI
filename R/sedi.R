# R functions to compute SEDI and it 95% CI (also Hit Rate, False Alarm Ratio, TSS, Odds Ratio, and ORSS)
# by Rainer Wunderlich
# based on Ferro, C.A.T. & Stephensons, D.B. (2011),
# Extremal Dependence Indices: Improved Verification Measures for Deterministic Forecasts of Rare Binary Events,
# Weather & Forecasting, 26, 699-713
# DOI:10.1175/WAF-D-10-05030.1

# PLEASE NOTE:
# To avoid undefined values, individual zeros are substituted by 1 and all returned results indicated whether they represent precise valuse or upper/lower approximations
# in addition zeroes in divisions and minusInf are approximated by 1e-09 and -1e09, respectively.

# Approximations
minusInf <- function(a = a, b = b, c = c, d = d) {
  return(-1e+09)
}

plusZero <- function(a = a, b = b, c = c, d = d) {
  return(1e-09)
}

# helper
# base rate
BaseRate <- function(a = a, b = b, c = c, d = d) {
  return ((a + c) / (a + b + c + d))
}

# helper
# hit rate
HitRate <- function(a = a, b = b, c = c, d = d) {
  return (a / (a + c))
}

OneMHR <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return(1 - (a / (a + c)))
}

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

OneMFA <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return (1 - (b / (b + d)))
}

logFA <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(FalseAlarm(a = a, b = b, c = c, d = d)))
}

# not required but useful
# odds ratio
OddsRatio <- function(a = a, b = b, c = c, d = d) {
  return (
    (HitRate(a = a, b = b, c = c, d = d) * (OneMFA(a = a, b = b, c = c, d = d)))
    /
      (FalseAlarm(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d)))
  )
}

# not required but useful
# ORSS
ORSS <- function(a = a, b = b, c = c, d = d) {
  return ((HitRate(a = a, b = b, c = c, d = d) - FalseAlarm(a = a, b = b, c = c, d = d))
          /
            (HitRate(a = a, b = b, c = c, d = d) + FalseAlarm(a = a, b = b, c = c, d = d) - 2*(HitRate(a = a, b = b, c = c, d = d)*FalseAlarm(a = a, b = b, c = c, d = d))))
}

# not required but useful
# TSS
TSS <- function(a = a, b = b, c = c, d = d) {
  return (HitRate(a = a, b = b, c = c, d = d) - FalseAlarm(a = a, b = b, c = c, d = d))
}

# helper
# denominator CI
CIdeno <- function(a = a, b = b, c = c, d = d) {
  return ( 2 * abs(
    (((OneMHR(a = a, b = b, c = c, d = d)) * (OneMFA(a = a, b = b, c = c, d = d)) + HitRate(a = a, b = b, c = c, d = d) * FalseAlarm(a = a, b = b, c = c, d = d))/((OneMHR(a = a, b = b, c = c, d = d)) * (OneMFA(a = a, b = b, c = c, d = d)))) *
      (log( FalseAlarm(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d)))) + ((2 * HitRate(a = a, b = b, c = c, d = d)) / (OneMHR(a = a, b = b, c = c, d = d))) * (log(HitRate(a = a, b = b, c = c, d = d) * (OneMFA(a = a, b = b, c = c, d = d))))
  )  )
}

# helper
# divisor CI
CIdivi <- function(a = a, b = b, c = c, d = d) {
  return (HitRate(a = a, b = b, c = c, d = d) * (log(FalseAlarm(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d)))   +  log(HitRate(a = a, b = b, c = c, d = d) * (OneMFA(a = a, b = b, c = c, d = d))))^2)
}

# helper
# trailing sq term
CIsqtail <- function(a = a, b = b, c = c, d = d) {
  return (sqrt((HitRate(a = a, b = b, c = c, d = d) * (OneMHR(a = a, b = b, c = c, d = d))) / (a + c)))
}

# main function
# symmetric extremal dependence index
sedi <- function(a = a, b = b, c = c, d = d) {
  if (isTRUE(x = (max(a,b,c,d) <= 1) && sum(a,b,c,d) <= 1)) {
    return(list(" (Almost) empty confusion matrix", NaN))  
  }
  
  if (isTRUE(x = (min(a,b,c,d) < 1))) {
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
      return(list("Inverse prediction:", 0))
    }
  
    if (isTRUE(x = (a < 1 && min(b,c,d) >=1))) {
      return(list("Undefined, but smaller than",
                  (logFA(a = a+1, b = b, c = c, d = d) - logHR(a = a+1, b = b, c = c, d = d) - log(OneMFA(a = a+1, b = b, c = c, d = d)) + log(OneMHR(a = a+1, b = b, c = c, d = d)))
                  /
                  (logFA(a = a+1, b = b, c = c, d = d) + logHR(a = a+1, b = b, c = c, d = d) + log(OneMFA(a = a+1, b = b, c = c, d = d)) + log(OneMHR(a = a+1, b = b, c = c, d = d)))))  
    }
  
    if (isTRUE(x = (d < 1 && min(a,b,c) >=1))) {
      return(list("Undefined, but smaller than",
                  (logFA(a = a, b = b, c = c, d = d+1) - logHR(a = a, b = b, c = c, d = d+1) - log(OneMFA(a = a, b = b, c = c, d = d+1)) + log(OneMHR(a = a, b = b, c = c, d = d+1)))
                  /
                  (logFA(a = a, b = b, c = c, d = d+1) + logHR(a = a, b = b, c = c, d = d+1) + log(OneMFA(a = a, b = b, c = c, d = d+1)) + log(OneMHR(a = a, b = b, c = c, d = d+1)))))  
    }
  
    if (isTRUE(x = (b < 1 && min(a,c,d) >=1))) {
      return(list("Undefined, but greater than",
                  (logFA(a = a, b = b+1, c = c, d = d) - logHR(a = a, b = b+1, c = c, d = d) - log(OneMFA(a = a, b = b+1, c = c, d = d)) + log(OneMHR(a = a, b = b+1, c = c, d = d)))
                  /
                  (logFA(a = a, b = b+1, c = c, d = d) + logHR(a = a, b = b+1, c = c, d = d) + log(OneMFA(a = a, b = b+1, c = c, d = d)) + log(OneMHR(a = a, b = b+1, c = c, d = d)))))  
    }
  
    if (isTRUE(x = (c < 1 && min(a,b,d) >=1))) {
      return(list("Undefined, but greater than",
                  (logFA(a = a, b = b, c = c+1, d = d) - logHR(a = a, b = b, c = c+1, d = d) - log(OneMFA(a = a, b = b, c = c+1, d = d)) + log(OneMHR(a = a, b = b, c = c+1, d = d)))
                  /
                  (logFA(a = a, b = b, c = c+1, d = d) + logHR(a = a, b = b, c = c+1, d = d) + log(OneMFA(a = a, b = b, c = c+1, d = d)) + log(OneMHR(a = a, b = b, c = c+1, d = d)))))  
    }
  }
  else {
    return (list("SEDI equal to",
                 (logFA(a = a, b = b, c = c, d = d) - logHR(a = a, b = b, c = c, d = d) - log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
                 /
                 (logFA(a = a, b = b, c = c, d = d) + logHR(a = a, b = b, c = c, d = d) + log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
            ))  
  }
}
# confidence Interval 1
# SE for SEDI
sediSE <- function(a = a, b = b, c = c, d = d) {
  return ((CIdeno(a = a, b = b, c = c, d = d) / CIdivi(a = a, b = b, c = c, d = d)) * CIsqtail(a = a, b = b, c = c, d = d) * 1.95996)
}

# confidence Interval 2
# LCL sedi UCL 95%
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
