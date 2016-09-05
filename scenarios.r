## Produce graphs and tables of scenario results, given model output from scenarios.py 

nsc <- 20; nouts=6; npops <- 4
M <- array(as.matrix(read.csv("results/scen_mean.csv",header=FALSE)), dim=c(nsc, nouts, npops))*365.25
S <- array(as.matrix(read.csv("results/scen_SD.csv",header=FALSE)), dim=c(nsc, nouts, npops))*365.25
N <- as.matrix(read.csv("results/scen_N.csv",header=FALSE))

outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY")
mn <- c("All", "Eligible", "Attending", "Any treatment", "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC")

scn <- c("Base case","BP","CVD","CVD+BP",
         "Diabetes","Diabetes+BP","Diabetes+CVD","Diabetes+CVD+BP",
         "Attend from age 50","Baseline uptake 0.85 (from 0.48)","Uptake (+10)/+20% in most deprived","Uptake of smokers +50%","+20% uptake for QRisk > 15","Statin prescription x 2.5","Attenders keep attending","Target non-attenders", "No extra statin effect","Sudden death","Sudden death, background reduction")

scnp <- c("Base case","BP","CVD","CVD+BP",
         "Diabetes","Diabetes+BP","Diabetes+CVD","Diabetes+CVD+BP",
         "Attend from age 50","Baseline 0.85 (from 0.48)","(+10)/+20% in most deprived","Smokers +50%","+20% for QRisk > 15","Statin prescription x 2.5","Attenders keep attending","Target non-attenders", "No extra statin effect","p=0.1","p=0.1, background reduction","p=0.4")

dimnames(M) <- dimnames(S) <- list(scn, outn, mn[1:4])
dimnames(N) <- list(scn, mn)


SE <- L <- U <- S
for (j in 1:nouts){
    SE[,j,] <- S[,j,] / sqrt(N[,1:4])
    L[,j,] <- M[,j,] - qnorm(0.975)*SE[,j,]
    U[,j,] <- M[,j,] + qnorm(0.975)*SE[,j,]
}

library(devEMF)
emf("scenarios.emf")

## Incremental QALY for total population
j <- 3;
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

k <- 2; d <- 0.2
#abline(v=M[1,j,k], col="gray", lwd=2)
#points(M[,j,k], y+d, pch=19, col="red")
#segments(L[,j,k], y+d, U[,j,k], y+d, col="red")
k <- 3; d <- 0.1
#abline(v=M[1,j,k], col="gray", lwd=2)
#points(M[,j,k], y+d, pch=19, col="blue")
#segments(L[,j,k], y+d, U[,j,k], y+d, col="blue")
k <- 4
abline(v=M[1,j,k], col="gray", lwd=2)
points(M[,j,k], y, pch=19, col="red")
segments(L[,j,k], y, U[,j,k], y, col="red")

# text(x0, y, scn, pos=4)
text(rep(c(x0+0.5, x0, x0+0.5,x0,x0+0.5,x0), c(3,4,4,1,7,1)), rev(y), rev(scnp), pos=4,
     cex=rep(c(0.7, 1, 0.7, 1, 0.7, 1), c(3,4,4,1,7,1)))
text(x0, 22, "Include people on registers", pos=4)
text(x0, 13, "Increase uptake", pos=4)
text(x0, 4, "Sudden death", pos=4)
text(M[1,j,1], 24, "WHOLE POPULATION")
text(M[1,j,4], 24, "TREATED POPULATION", col="red")

dev.off()

## Table of N
perc <- round(100 * N[,2:8] / N[,1], 2)
scen <- rep(1:nrow(perc), ncol(perc))
pop <- factor(rep(1:ncol(perc), each=nrow(perc)), labels=c("EL","ATT","TRT","STAT","AHT","SMK","WT"))
plot(perc, scen, type="n")
text(perc, scen, pop)
options(width=150)
perc
write.csv(perc, "scenarios_N.csv", quote=FALSE)
1
