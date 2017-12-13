### Health equity plot
## One point per scenario

## x axis is gain in LE or QALY

## y axis is gain in LE for Q5 minus gain in LE for Q1
## "absolute equity reduction" (inequity reduction?)

## exclude these
scx <- match(c("Base case","Attenders keep attending","Attenders keep attending, target non-attenders","CVD incidence declining","Baseline QRisk uncertainty +- 20%"), scn)

scnwrap <- c("Base case","Invite high BP", "Attend age 50-74", "Attend age 40-80", "Attend age 50-80",
         "Baseline uptake +30%","Uptake +30% in\nmost deprived","Uptake of\nsmokers +30%","+30% uptake for\nQRisk > 20",
         "Target non-attenders", "Attenders keep\nattending", "Attenders keep attending,\ntarget non-attenders", 
         "Statin prescription\nx 2.5", "AHT prescription\nx 2.5", "Smoking referral\nx 2.5", "Weight referral\nx 2.5",
         "All treatments\nx 2.5", "Higher treatment,\nhigher uptake,\ninvite high BP","CVD incidence\declining",
         "Baseline QRisk\nuncertainty +- 20%")[-scx]

scnnowrap <- c("Base case",
               "Include people with hypertension",
               "Starting age 50 years",
               "Upper age of eligibility 79 years",
               "Starting age 50 years, upper eligibility age 79 years",
               "Uptake increased by 30% for everyone",
               "Uptake increased by 30% for the most deprived quintile",
               "Uptake increased by 30% for smokers",
               "Uptake increased by 30% for those at high risk of CVD",
               "Increase likelihood of offer to previous non-attenders",
               "Attenders keep attending",
               "Attenders keep attending, target non-attenders", 
               "Statin treatment x 2.5",
               "Anti-hypertensive treatment x 2.5",
               "Smoking cessation x 2.5",
               "Weight loss programme x 2.5",
               "All treatments x 2.5",
               "All treatments x 2.5, eligibility extended, uptake increased",
               "CVD incidence declining",
               "Baseline QRisk uncertainty +- 20%"
               )[-scx]

cols <- rep(c("black","black","red","blue","black"), c(1, 4, 7, 6, 2))[-scx]

x <- M[-scx,"IQALY","All"]
y <- M[-scx,"IQALY","Most deprived"] - M[-scx,"IQALY","Least deprived"]
xse <- SE[-scx,"IQALY","All"]
xup <- x + qnorm(0.975)*xse
xlo <- x - qnorm(0.975)*xse
yse <- sqrt(SE[-scx,"IQALY","Most deprived"]^2 + SE[-scx,"IQALY","Least deprived"]^2)
yup <- y + qnorm(0.975)*yse
ylo <- y - qnorm(0.975)*yse

library(devEMF)

emf("../paper/equity_big.emf", width=14, height=14)

par(mar=c(5,6.5,5,2))
plot(x, y,  xlim=c(-0.5,18), ylim=c(-7,19), type="n", axes=FALSE, cex.lab=1.5,
     xlab="Overall effectiveness (days of quality-adjusted life gained per person)",
     ylab="Change in health equity (gain in days of quality adjusted life for most – least deprived fifth)")
lim <- par("usr")
rect(lim[1], lim[3], lim[2], lim[4], col="gray92", border="gray92")
axis(1, at=seq(-1, 17, by=1), cex.lab=2, cex.axis=1.3)
axis(2, at=seq(-7, 19, by=1), cex.lab=2, cex.axis=1.3)
abline(v=seq(-1, 17, by=1), col="white")
abline(h=seq(-7, 19, by=1), col="white")
abline(h=0, v=0, col="black", lwd=2)

points(x[x<2], y[x<2], pch=19, cex=1, col=cols[x<2])
points(x[x>2], y[x>2], pch=1, cex=7, col=cols[x>2])
text(x[x>2], y[x>2], seq_along(x[x>2]), font=2, cex=2, col=cols[x>2])
text(1.5, 3, "(see second\nmagnified plot)", cex=1.5)
segments(xup, y, xlo, y, col=cols)
segments(x, yup, x, ylo, col=cols)
inds <- c(10,14,15); ni <- length(inds)
text(0, 17 - (1:ni)*1.5, paste0(seq_along(x[x>2]), ": ", scnnowrap[inds]),
     pos=4, font=2, cex=1.5, col=cols[x>2])

dev.off()


emf("../paper/equity_small.emf", width=14, height=14)

par(mar=c(5,6.5,5,2))
plot(x, y,  xlim=c(-0.7, 2.7), ylim=c(-1.5, 3), type="n", axes=FALSE, cex.lab=1.5, 
     xlab="Overall effectiveness (days of quality-adjusted life gained per person)",
     ylab="Change in health equity (gain in days of quality adjusted life for most – least deprived fifth)")
lim <- par("usr")
rect(lim[1], lim[3], lim[2], lim[4], col="gray92", border="gray92")
axis(1, at=seq(-1.5, 3.0, by=0.5), cex.axis=1.3)
axis(2, at=seq(-5, 5, by=0.5), cex.axis=1.3)
abline(v=seq(-1.5, 3.0, by=0.5), col="white")
abline(h=seq(-5, 5, by=1), col="white")
abline(h=0, v=0, col="black", lwd=2)
points(x, y, pch=1, cex=7, lwd=2, col=cols)
inds <- seq(along=x[x<2]); ni <- length(inds)
tx <- inds
text(x[x<2], y[x<2], tx, font=2, cex=2, col=cols)
segments(xup, y, xlo, y, lwd=1.5, col=cols)
segments(x, yup, x, ylo, lwd=1.5, col=cols)
text(1.0, 3 - c(1:9, 10:ni+0.7)*0.3,
     paste0(inds, ": ", scnnowrap[x<2]), pos=4, font=2, cex=1.2, col=cols[x<2])

dev.off()

## looks like increasing overall effects also reduces inequality.
