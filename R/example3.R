# Example to illustrate the utility of SEDI for evaluation of species distribution models
# by Rainer Wunderlich
# based on Ferro, C.A.T. & Stephensons, D.B. (2011),
# Extremal Dependence Indices: Improved Verification Measures for Deterministic Forecasts of Rare Binary Events,
# Weather & Forecasting, 26, 699-713
# DOI:10.1175/WAF-D-10-05030.1

source("/home/affu/Desktop/Code/R.Code/sedi.R")

# we assume a true prevalence of 0.25

csum <- round( (1:20)^1.25 * 1000) + 248
bplusc <- round( 0.75 * csum)
a <- round( 0.25 * csum)
b <- csum - a
c <- rep(1, 20)
d <- rep(1, 20)


Cases <- c()

for (i in 1:20) {
  Cases <- c(Cases, paste0("#", as.character(i)))
}

BR_1 <- BaseRate(a = a, b = b, c = c, d = d)

H_1 <- HitRate(a = a, c = c)

F_1 <- FalseAlarm(b = b, d = d)

TSS_1 <- TSS(a = a, b = b, c = c, d = d)

ORSS_1 <- ORSS(a = a, b = b, c = c, d = d)

SEDI_1 <- sedi(a = a, b = b, c = c, d = d)[[2]] # We can think about marking those cases where either b or c are zero somehow
SEDI_1_W <- sedi(a = a, b = b, c = c, d = d)[[1]]

results <- data.frame("Case" = Cases,
                      "Correctly Predicted Absences" = d,
                      "H" = H_1,
                      "F" = F_1,
                      "TSS" = TSS_1,
                      "ORSS" = ORSS_1,
                      "SEDI" = SEDI_1)

mycol <- c("purple", "darkorange", "darkgreen", "cyan")

# create empty plot
plot(1, type = "n", xlim = log(c(csum[1]+2, csum[20]+2), 10), ylim = c(-1.0, 1.0),
     ylab = "Value of TSS | ORSS | SEDI | H",
     xlab = expression("log"[10]*"(a + b + c + d)"),
     family ="arial") #log 10 for display

# plot hit rates
lines(x = log(csum[1:20]+2, 10), y = results[1:20, 3], col = mycol[4], lty = "dashed", lwd = 2) # remember log 10 HIT RATE 1


# for each metric plot 2 lines
for (i in 1:3){
  lines(x =  log(csum[1:20]+2, 10), y = results[1:20, 4+i], col = mycol[i], lty = "dashed", lwd = 2) # remember log 10
}

legend("bottomright", legend = c("TSS", "ORSS", "SEDI", "H"), col=mycol, pch = "—", lwd = 2) # optional legend
legend("bottomleft", legend = "A", col="black", lty = "dashed", lwd = 2) # optional legend