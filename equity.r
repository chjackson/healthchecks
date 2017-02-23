### Health equity plot
## One point per scenario

## x axis is gain in LE or QALY

## y axis is gain in LE for Q5 minus gain in LE for Q1
## "absolute equity reduction" (inequity reduction?)

## exclude these
scx <- match(c("Base case","Attenders keep attending"), scn)

scnwrap <- c("Base case","Invite high BP", "Attend age 50-74", "Attend age 40-80", "Attend age 50-80",
         "Baseline uptake +30%","Uptake +30% in\nmost deprived","Uptake of\nsmokers +30%","+30% uptake for\nQRisk > 20",
         "Target non-attenders", "Attenders keep\nattending", "Attenders keep attending,\ntarget non-attenders", 
         "Statin prescription\nx 2.5", "AHT prescription\nx 2.5", "Smoking referral\nx 2.5", "Weight referral\nx 2.5",
         "All treatments\nx 2.5", "Higher treatment,\nhigher uptake,\ninvite high BP")[-scx]

x <- M[-scx,"IQALY","All"]
y <- M[-scx,"IQALY","Most deprived"] - M[-scx,"IQALY","Least deprived"]
xse <- SE[-scx,"IQALY","All"]
xup <- x + qnorm(0.975)*xse
xlo <- x - qnorm(0.975)*xse
yse <- sqrt(SE[-scx,"IQALY","Most deprived"]^2 + SE[-scx,"IQALY","Least deprived"]^2)
yup <- y + qnorm(0.975)*yse
ylo <- y - qnorm(0.975)*yse

library(devEMF)
emf("equity_big.emf")
plot(x, y,  xlim=c(-0.5,17), ylim=c(-2,12), pch=19,
     xlab="QALY gain for population",
     ylab="QALY gain for most deprived fifth - gain for least deprived fifth")
abline(h=0, v=0, col="gray", lty=2)
text(x[x>2], y[x>2], scnwrap[x>2], cex=1, pos=4)
text(1, 2, "(see second\nmagnified plot)")
segments(xup, y, xlo, y)
segments(x, yup, x, ylo)
dev.off()

emf("equity_small.emf")
plot(x, y,  xlim=c(-0.7,1.5), ylim=c(-1.2,2), pch=19, 
     xlab="QALY gain for population",
     ylab="QALY gain for most deprived fifth - gain for least deprived fifth")
abline(h=0, v=0, col="gray")
text(x, y+0.05, scnwrap, cex=1, pos=4)
segments(xup, y, xlo, y)
segments(x, yup, x, ylo)
dev.off()

## looks like increasing overall effects also reduces inequality.
## TODO standard errors as crosses 
