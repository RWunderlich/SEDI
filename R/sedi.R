# R functions to compute SEDI and it 95% CI
# by Rainer Wunderlich
# based on Ferro, C.A.T. & Stephensons, D.B. (2011), Weather & Forecasting, 26, 699-713

# helper
# hit rate
HitRate <- function(a, c) {
  return (a / (a + c))
}

# helper
# false alarm
FalseAlarm <- function(b, d) {
  return (b / (b + d))
}

# helper
# odds ratio
OddsRatio <- function(a, b, c, d) {
  return (
          (HitRate(a, c) * (1 - FalseAlarm(b, d)))
          /
            (FalseAlarm(b, d) * (1 - HitRate(a, c)))
          )
}

# helper
# denominator CI
CIdeno <- function(a, b, c, d) {
  return ( 2 * abs( (((1 - HitRate(a, c)) * (1 - FalseAlarm(b, d)) + HitRate(a, c) * FalseAlarm(b, d))/((1 - HitRate(a, c)) * (1 - FalseAlarm(b, d)))) *
                      (log( FalseAlarm(b, d) * (1 - HitRate(a, c)))) + ((2 * HitRate(a, c)) / (1 - HitRate(a, c))) * (log(HitRate(a, c) * (1 - FalseAlarm(b, d))))  )  )
}

# helper
# divisor CI
CIdivi <- function(a, b, c, d) {
  return (HitRate(a, c) * (log(FalseAlarm(b, d) * (1 - HitRate(a, c)))   +  log(HitRate(a, c) * (1 - FalseAlarm(b, d))))^2)
}

# helper
# trailing sq term
CIsqtail <- function(a, b, c, d) {
  return (sqrt((HitRate(a, c) * (1 - HitRate(a, c))) / (a + c)))
}



# main function
# symmetric extremal dependence index
sedi <- function(a, b, c, d) {
  return (
          (log(FalseAlarm(b, d)) - log(HitRate(a, c)) - log(1 - FalseAlarm(b, d)) + log(1 - HitRate(a, c)))
          /
          (log(FalseAlarm(b, d)) + log(HitRate(a, c)) + log(1 - FalseAlarm(b, d)) + log(1 - HitRate(a, c)))
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
