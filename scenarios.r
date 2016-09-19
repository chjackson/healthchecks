## Produce graphs and tables of scenario results, given model output from scenarios.py

nsc <- 20; nouts=6; npops <- 4
nsc <- 3
M <- array(as.matrix(read.csv("results/scenoff_mean.csv",header=FALSE)), dim=c(nsc, nouts, npops))*365.25
S <- array(as.matrix(read.csv("results/scenoff_SD.csv",header=FALSE)), dim=c(nsc, nouts, npops))*365.25
N <- as.matrix(read.csv("results/scenoff_N.csv",header=FALSE))

outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY")
mn <- c("All", "Eligible", "Attending", "Any treatment", "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")

scn <- c("Base case","BP","CVD","CVD+BP",
         "Diabetes","Diabetes+BP","Diabetes+CVD","Diabetes+CVD+BP",
         "Attend from age 50","Baseline uptake 0.85 (from 0.48)","Uptake (+10)/+20% in most deprived","Uptake of smokers +50%","+20% uptake for QRisk > 15","Statin prescription x 2.5","Attenders keep attending","Target non-attenders", "No extra statin effect","Sudden death p=0.1","Sudden death; background reduction","Sudden death p=0.4")

scnp <- c("Base case","BP","CVD","CVD+BP",
         "Diabetes","Diabetes+BP","Diabetes+CVD","Diabetes+CVD+BP",
         "Attend from age 50","Baseline 0.85 (from 0.48)","(+10)/+20% in most deprived","Smokers +50%","+20% for QRisk > 15","Statin prescription x 2.5","Attenders keep attending","Target non-attenders", "No extra statin effect","p=0.1","p=0.1, background reduction","p=0.4")

scn <- scnp <- c("Base","Off 0.4", "Off 0.4 red 0.1")
dimnames(M) <- dimnames(S) <- list(scn, outn, mn[1:4])
dimnames(N) <- list(scn, mn)

N

SE <- L <- U <- S
for (j in 1:nouts){
    SE[,j,] <- S[,j,] / sqrt(N[,1:4])
    L[,j,] <- M[,j,] - qnorm(0.975)*SE[,j,]
    U[,j,] <- M[,j,] + qnorm(0.975)*SE[,j,]
}

library(devEMF)
emf("scenarios.emf")

## Incremental QALY for total population
j <- "IQALY";
y <- rev(c(1, 2, 3,  5, 6, 7, 8,9, 10, 11, 12, 14:21, 23))
x0 <- -15; xm <- 45
par(mar=c(3.2,0,0,0), mgp=c(2,1,0))
plot(M[,j,1], y, pch=19, xlim=c(x0, xm), ylim=c(1, 24), axes=FALSE, bty="n", xlab="QALY(days) gained from HC", ylab="")
lim <- par("usr")
rect(0, lim[3], lim[2], lim[4], col="gray92", border="gray92")
at <- seq(0,xm,by=5)
axis(1, at=at); abline(v=at, col="white")

k <- 1
abline(v=M[1,j,k], col="gray", lwd=2)
points(M[,j,k], y, pch=19)
segments(L[,j,k], y, U[,j,k], y)
k <- 4
abline(v=M[1,j,k], col="pink", lwd=2)
points(M[,j,k], y, pch=19, col="red")
segments(L[,j,k], y, U[,j,k], y, col="red")

text(rep(c(x0+0.5, x0, x0+0.5,x0,x0+0.5,x0), c(3,4,4,1,7,1)), rev(y), rev(scnp), pos=4,
     cex=rep(c(0.7, 1, 0.7, 1, 0.7, 1), c(3,4,4,1,7,1)))
text(x0, 22, "Include people on registers", pos=4)
text(x0, 13, "Increase uptake", pos=4)
text(x0, 4, "Sudden death", pos=4)
text(M[1,j,1], 24, "WHOLE POPULATION")
text(M[1,j,4], 24, "TREATED POPULATION", col="red")

dev.off()

## Table of N
P <- N[,2:9] / N[,1]
P[,1:7] <- 100 * P[,1:7]
perc <- round(P, 2)
scen <- rep(1:nrow(perc), ncol(perc))
options(width=180)
perc
write.csv(perc, quote=FALSE)
write.csv(perc, "scenarios_N.csv", quote=FALSE)
1

## QALY gained per HC
## mean qaly over pop * npop / number of HCs
## = mean QALY / average no of HCs per person
## doesn't include error in number of HCs

library(devEMF)
emf("scenarios-perhc.emf")

par(mar=c(3.2,0,0,0), mgp=c(2,1,0))
Mnh <- M[,3,1] / P[,"Number of HCs"]
Lnh <- L[,3,1] / P[,"Number of HCs"]
Unh <- U[,3,1] / P[,"Number of HCs"]
x0 <- -2; xm <- 8
plot(Mnh, y, pch=19, xlim=c(x0, xm), ylim=c(1, 24), axes=FALSE, bty="n", xlab="QALY(days) gained from HC", ylab="")
lim <- par("usr")
rect(0, lim[3], lim[2], lim[4], col="gray92", border="gray92")
at <- seq(0,xm,by=1)
axis(1, at=at); abline(v=seq(0,xm,by=0.5), col="white")
abline(v=Mnh[1], col="lightblue", lwd=2)
points(Mnh, y, col="blue", pch=19)
segments(Lnh, y, Unh, col="blue")
text(Mnh[1], 24, "QALY GAIN PER HC", col="blue")
text(x0, 22, "Include people on registers", pos=4)
text(x0, 13, "Increase uptake", pos=4)
text(x0, 4, "Sudden death", pos=4)
text(rep(c(x0+0.5, x0, x0+0.5,x0,x0+0.5,x0), c(3,4,4,1,7,1)), rev(y), rev(scnp), pos=4,
     cex=rep(c(0.7, 1, 0.7, 1, 0.7, 1), c(3,4,4,1,7,1)))
text(8, y, round(P[,8], 1), pos=2)
text(8, 24, "Mean HCs per person", pos=2)

dev.off()

## TODO interpretation
## why worse gain per hc if increase baseline uptake?  checking more well people?
##









### NEW OFFER SCENARIOS 

## Incremental QALY for total population
j <- "IQALY";
y <- rev(c(1, 2, 3))
x0 <- -5; xm <- 30
par(mar=c(3.2,0,0,0), mgp=c(2,1,0))
plot(M[,j,1], y, pch=19, xlim=c(x0, xm), ylim=c(1, 24), axes=FALSE, bty="n", xlab="QALY(days) gained from HC", ylab="")
lim <- par("usr")
rect(0, lim[3], lim[2], lim[4], col="gray92", border="gray92")
at <- seq(0,xm,by=5)
axis(1, at=at); abline(v=at, col="white")
k <- 1
abline(v=M[1,j,k], col="gray", lwd=2)
points(M[,j,k], y, pch=19)
segments(L[,j,k], y, U[,j,k], y)
k <- 4
abline(v=M[1,j,k], col="pink", lwd=2)
points(M[,j,k], y, pch=19, col="red")
segments(L[,j,k], y, U[,j,k], y, col="red")
text(x0, y, scnp)

res <- cbind(M[,"IQALY","All"], L[,"IQALY","All"], U[,"IQALY","All"])
round(res, 2)

Attached new "flow diagrams" document including the reduced weight/smoking treatment.

I've also just run a scenario with increasing the probability (p) that a HC is offered for people who did not accept their previous offer of a HC.  Last week I had showed scenarios where we increased the accept rate given an offer, which made less sense as "targeting" non-attenders.   Optionally we can also reduce the offer rate (q) for people who previously attended. 

QALY gain (days) with computational 95% CI
Base case (p=q=0.197): 3.71 (2.86, 4.57) 
p = 0.4:               4.38 (3.49, 5.26) 
p = 0.4, q=0.1         3.20 (2.35, 4.06) 

With p=0.4, the proportion attending at least once and the proportion treated at least once both go up by 2%.   With the additional q=0.1, the attendance rate goes up similarly, but  the proportion getting any treatment goes down by 0.7%, presumably because some higher risk people miss out if they don't get multiple checks.

TODO more HCs

Other points we might discuss today.

* Exactly what we will be presenting to Nick (and Simon?) 

* How to get data on relative incidence and mortality for CVD events, and exactly what data we need
