# Example 2 to illustrate the utility of SEDI for evaluation of species distribution models
# Hit rate constant at two levels; n2 <- (1:20)^2; d <- n2 * 100; a <- n2 * 50; b <- round(a/0.75 - a)
# by Rainer Wunderlich
# based on Ferro, C.A.T. & Stephensons, D.B. (2011),
# Extremal Dependence Indices: Improved Verification Measures for Deterministic Forecasts of Rare Binary Events,
# Weather & Forecasting, 26, 699-713
# DOI:10.1175/WAF-D-10-05030.1

source("/home/affu/Desktop/Code/R.Code/sedi.R")

d = d <- round( (1:20)^1.25 * 1000)

a <- c(175:189, rep(190, 5))

c <- 200 - a

b1 <- rev(c)
b2 <- 3*rev(c)

csum <- a + b1 + c + d
csum2 <- a + b2 + c + d



Cases <- c()

for (i in 1:40) {
  Cases <- c(Cases, paste0("#", as.character(i)))
}

BR_1 <- BaseRate(a = a, b = b1, c = c, d = d)
BR_2 <- BaseRate(a = a, b = b2, c = c, d = d)

H_1 <- HitRate(a = a, c = c)
H_2 <- HitRate(a = a, c = c)

F_1 <- FalseAlarm(b = b1, d = d)
F_2 <- FalseAlarm(b = b2, d = d)

TSS_1 <- TSS(a = a, b = b1, c = c, d = d)
TSS_2 <- TSS(a = a, b = b2, c = c, d = d)

ORSS_1 <- ORSS(a = a, b = b1, c = c, d = d)
ORSS_2 <- ORSS(a = a, b = b2, c = c, d = d)

SEDI_1 <- sedi(a = a, b = b1, c = c, d = d)[[2]]
SEDI_1_W <- sedi(a = a, b = b1, c = c, d = d)[[1]]

SEDI_2 <- sedi(a = a, b = b2, c = c, d = d)[[2]]
SEDI_2_W <- sedi(a = a, b = b2, c = c, d = d)[[1]]

results <- data.frame("Case" = Cases,
                      "Correctly Predicted Absences" = rep(d, 2),
                      "H" = c(H_1, H_2),
                      "F" = c(F_1, F_2),
                      "TSS" = c(TSS_1, TSS_2),
                      "ORSS" = c(ORSS_1, ORSS_2),
                      "SEDI" = c(SEDI_1, SEDI_2))

mycol <- c("purple", "darkorange", "darkgreen", "cyan")

# create empty plot
plot(1, type = "n", xlim = log(c(d[1], d[20]), 10), ylim = c(0.80, 1.0),
     ylab = "Value of TSS | ORSS | SEDI | H",
     xlab = expression("log"[10]*"(d)"),
     family ="arial") #log 10 for display

# plot hit rates
lines(x = log(d[1:20], 10), y = results[1:20, 3], col = mycol[4], lty = "dashed", lwd = 2) # remember log 10 HIT RATE 1
lines(x = log(d[1:20], 10), y = results[21:40, 3], col = mycol[4], lty = "dotted", lwd = 2) # remember log 10 HIT RATE 1



# for each metric plot 2 lines
for (i in 1:3){
  lines(x = log(results[1:20, 2], 10), y = results[1:20, 4+i], col = mycol[i], lty = "dashed", lwd = 2) # remember log 10
  lines(x = log(results[21:40, 2], 10), y = results[21:40, 4+i], col = mycol[i], lty = "dotted", lwd = 2) # remember log 10
}

legend("bottomright", legend = c("TSS", "ORSS", "SEDI", "H"), col=mycol, pch = "—", lwd = 2) # optional legend
legend("bottomleft", legend = c("C1", "C2"), col="black", lty = c("dashed", "dotted"), lwd = 2) # optional legend
