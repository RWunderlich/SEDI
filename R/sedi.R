# R functions to compute SEDI and it 95% CI (also Hit Rate, False Alarm Ratio, TSS, Odds Ratio, and ORSS)
# by Rainer Wunderlich
# based on Ferro, C.A.T. & Stephensons, D.B. (2011),
# Extremal Dependence Indices: Improved Verification Measures for Deterministic Forecasts of Rare Binary Events,
# Weather & Forecasting, 26, 699-713
# DOI:10.1175/WAF-D-10-05030.1

# PLEASE NOTE:
# To avoid undefined values, I substitute some terms/values with -1e+09 or 1e-09

# helper
# base rate
BaseRate <- function(a, b, c, d) {
  return ((a + c) / (a + b + c + d))
}

# helper
# hit rate
HitRate <- function(a, c) {
  return (a / (a + c))
}

OneMHR <- function(a, c) {
  if (HitRate(a, c) == 1) {
    return(1e-09)
  }
  else
  return(1 - (a / (a + c)))
}

logHR <- function(a, c) {
  if (HitRate(a, c) == 0) {
    return(-1e+09)
  }
  else
    return(log(HitRate(a, c)))
}

# helper
# false alarm
FalseAlarm <- function(b, d) {
  return (b / (b + d))
}

OneMFA <- function(b, d) {
  if (FalseAlarm(b, d) == 1) {
    return(1e-09)
  }
  else
  return (1 - (b / (b + d)))
}

logFA <- function(b, d) {
  if (FalseAlarm(b, d) == 0) {
    return(-1e+09)
  }
  else
    return(log(FalseAlarm(b, d)))
}

# not required but useful
# odds ratio
OddsRatio <- function(a, b, c, d) {
  return (
          (HitRate(a, c) * (OneMFA(b, d)))
          /
            (FalseAlarm(b, d) * (OneMHR(a, c)))
          )
}

# not required but useful
# ORSS
ORSS <- function(a, b, c, d) {
  return ((a*d - b*c) / (a*d + b*c))
}

# not required but useful
# TSS
TSS <- function(a, b, c, d) {
  return (HitRate(a = a, c = c) - FalseAlarm(b = b, d = d))
}

# helper
# denominator CI
CIdeno <- function(a, b, c, d) {
  return ( 2 * abs( (((OneMHR(a, c)) * (OneMFA(b, d)) + HitRate(a, c) * FalseAlarm(b, d))/((OneMHR(a, c)) * (OneMFA(b, d)))) *
                      (log( FalseAlarm(b, d) * (OneMHR(a, c)))) + ((2 * HitRate(a, c)) / (OneMHR(a, c))) * (log(HitRate(a, c) * (OneMFA(b, d))))  )  )
}

# helper
# divisor CI
CIdivi <- function(a, b, c, d) {
  return (HitRate(a, c) * (log(FalseAlarm(b, d) * (OneMHR(a, c)))   +  log(HitRate(a, c) * (OneMFA(b, d))))^2)
}

# helper
# trailing sq term
CIsqtail <- function(a, b, c, d) {
  return (sqrt((HitRate(a, c) * (OneMHR(a, c))) / (a + c)))
}

# main function
# symmetric extremal dependence index
sedi <- function(a, b, c, d) {
  return (
          (logFA(b, d) - logHR(a, c) - log(OneMFA(b, d)) + log(OneMHR(a, c)))
          /
          (logFA(b, d) + logHR(a, c) + log(OneMFA(b, d)) + log(OneMHR(a, c)))
          )
}

# confidence Interval 1
# SE for SEDI
sediSE <- function(a, b, c, d) {
  return ((CIdeno(a = a, b = b, c = c, d = d) / CIdivi(a = a, b = b, c = c, d = d)) * CIsqtail(a = a, c = c) * 1.95996)
}

# confidence Interval 2
# LCL sedi UCL 95%
sediCI <- function(a, b, c, d) {
  return (c(sedi(a = a, b = b, c = c, d = d) - sediSE(a = a, b = b, c = c, d = d),
            sedi(a = a, b = b, c = c, d = d),
            sedi(a = a, b = b, c = c, d = d) + sediSE(a = a, b = b, c = c, d = d)))
}
