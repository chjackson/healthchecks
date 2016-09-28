library(tidyverse)
## Produce graphs and tables of scenario results, given model output from scenarios.py

nsc <- 12; nouts <- 6; npops <- 4

# Source block below for single run with no comp unc - use for rough results 
aname <- "scenall"
M <- array(as.matrix(read.csv(sprintf("results/%s_mean.csv",aname),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
S <- array(as.matrix(read.csv(sprintf("results/%s_SD.csv",aname),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
N <- as.matrix(read.csv(sprintf("results/%s_N.csv",aname),header=FALSE))
# NT <- as.matrix(read.csv(sprintf("results/%s_NT.csv",aname),header=FALSE))
outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY")
mn <- c("All", "Eligible", "Attending", "Any treatment", "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
scn <- c("Base case","BP", "Diabetes","Diabetes+BP", "Attend from age 50","Baseline uptake 0.85 (from 0.48)","Uptake (+10)/+20% in most deprived","Uptake of smokers +50%","+20% uptake for QRisk > 15","Statin prescription x 2.5", "Attenders keep attending", "Target non-attenders")
dimnames(M) <- dimnames(S) <- list(scn, outn, mn[1:npops])
dimnames(N) <- list(scn, mn)
# dimnames(NT) <- list(scn, mn[5:8])
SE <- L <- U <- S
for (j in 1:nouts){
    SE[,j,] <- S[,j,] / sqrt(N[,1:4])
#    SE[,j,] <- S[,j,] / sqrt(cbind(N[,1:4], NT))
    L[,j,] <- M[,j,] - qnorm(0.975)*SE[,j,]
    U[,j,] <- M[,j,] + qnorm(0.975)*SE[,j,]
}

M[,3,]

# Source for big pops / subset trick run - use this if want comp unc
# source("scenarios_combine.r")

## Graph of incremental QALYs (not for now) 
library(devEMF)
emf("scenarios.emf")
## 
j <- "IQALY";
y <- nsc:1 # rev(c(1, 2, 3,  5, 6, 7, 8,9, 10, 11, 12, 14:21, 23))
x0 <- -20; xm <- 45
par(mar=c(3.2,0,0,0), mgp=c(2,1,0))
plot(M[,j,1], y, pch=19, xlim=c(x0, xm), ylim=c(1, nsc+1), axes=FALSE, bty="n", xlab="QALY(days) gained from HC", ylab="")
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
text(x0, y, scn, pos=4)
text(x0, 22, "Include people on registers", pos=4)
text(x0, 13, "Increase uptake", pos=4)
text(x0, 4, "Sudden death", pos=4)
text(M[1,j,1], 24, "WHOLE POPULATION")
text(M[1,j,4], 24, "TREATED POPULATION", col="red")
dev.off()

## Table of N (not for now) 
P <- N[,2:9] / N[,1]
P[,1:7] <- 100 * P[,1:7]
perc <- round(P, 2)
scen <- rep(1:nrow(perc), ncol(perc))
options(width=180)
perc
write.csv(perc, quote=FALSE)
write.csv(perc, "scenarios_N.csv", quote=FALSE)
1

## Graph of QALY gained per HC
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









### TODO: add uncertainty to these tables. 

#sh1 <- short[short$scen=="Base case",] %>% 
#  gather(varname, outcome, Q1:deaths_oth_c) %>% 
#  spread(year, outcome)
#short[short$trt=="sc" & short$scen=="Base case", ] 
#      c("scen","trt","year","n","deaths_s","deaths_c","pd","pdl","pdu")]

## Unclear how to interpret effect in treated population for scenarios
## format as HC, control, diff for key outputs, reduction in short term mortality (only break down further if need to explain)



### FORM TABLE OF LONG TERM AND SELECTED SHORT TERM OUTCOMES 

procraw <- read.csv("results/scenall_nomortred/scenall_process.csv", header=FALSE,
                 col.names=c("popsize","eligible","attending","att_q","att_chol","att_hdl","qhi","qhi_presc_stat","qhi_stat","qlo","qlo_presc_stat","qlo_stat","nonatt","nonatt_q","nonatt_chol","nonatt_h","att_sbp","att_dbp","att_nhbp","qhi_presc_aht","qhi_aht","qlo_presc_aht","qlo_aht","nonatt_sbp","nonatt_dbp", "att_bmi","hbmi","hbmi_presc_wr","hbmi_wr","lobmi","lobmi_presc_wr","lobmi_wr","nonatt_bmi","att_smk","att_sc","att_quit1","att_nsmk"))
scn <- c("Base case","BP", "Diabetes","Diabetes+BP", "Attend from age 50","Baseline uptake 0.85 (from 0.48)","Uptake (+10)/+20% in most deprived","Uptake of smokers +50%","+20% uptake for QRisk > 15","Statin prescription x 2.5", "Attenders keep attending", "Target non-attenders")
nsc <- length(scn)

proc <- data.frame(as.vector(t(procraw))); names(proc) <- "num"
proc$scen <- factor(rep(scn, each=ncol(procraw)))
options(scipen=1e+7)
denomid <- c(NA, 1, 1, 
             NA, NA, NA,
             3, 7, 8,
             3, 10, 11, 
             2, NA, NA, NA, 
             NA, NA, 3,
             7, 20, 10, 22,
             NA, NA,
             NA, 3, 27, 28, 3, 30, 31,
             NA,
             3, 34, 35, 3)
proc$denom <- proc$num[rep(denomid,nsc)]
proc <- within(proc, perc <- 100*num/denom)

pnames <- c("Population", "Eligible at least once", "Attending at least once",
            "QRisk (attenders; over all HCs)","Cholesterol (attenders; over all HCs)","HDL (attenders; over all HCs)",
            "QRisk>20 (first HC)", "Prescribed statins (Q>20)","Statins (Q>20)",
            "QRisk<20 (first HC)", "Prescribed statins (Q<20)","Statins (Q<20)",
            "Non-attending", "QRisk (non-attenders)","Cholesterol (non-attenders)","HDL (non-attenders)",
            "SBP (attenders; all HCs)", "DBP (attenders; all HCs)", "SBP>140 (first HC)",
            "Prescribed AHT (Q>20)", "AHT (Q>20)",  "Prescribed AHT (Q<20)", "AHT (Q<20)",
            "SBP (non-attenders; all HCs)", "DBP (non-attenders; all HCs)",
            "BMI (attenders; all HCs)", "BMI>30 (first HC)", "Referred to WR", "WR", 
            "BMI<30 (first HC)", "Referred to weight reduction (BMI<30)", "Weight reduction (BMI<30)",
            "BMI(non-attenders)",
            "Smokers (first HC)", "Referred to smoking cessation", "Quit by 1 year due to cessation service", "Non-smokers")

proc$names <- rep(pnames, nsc)
proc <- proc[,c("scen","names","num","denom","perc")]
proc <- within(proc, seperc <- 100*sqrt((num/denom)*(1 - num/denom)/denom)) ## SEs around percentages

short <- read.csv("results/scenall_nomortred/scenall_short.csv", header=FALSE,
                  col.names=c("trt",
                              "Q1", "Q", "CH1", "CH", "HDL1", "HDL", "SBP1", "SBP", "DBP1", "DBP", "BMI1", "BMI", "SM_t", "SM_c",
                              "Q1SD", "QSD", "CH1SD", "CHSD", "HDL1SD", "HDLSD", "SBP1SD", "SBPSD", "DBP1SD", "DBPSD", "BMI1SD", "BMISD",
                              "cvd_events_t", "cvd_events_c", "lc_events_t", "lc_events_c", "deaths_t", "deaths_c", "deaths_cvd_t", "deaths_cvd_c", "deaths_lc_t", "deaths_lc_c", "deaths_dem_t", "deaths_dem_c", "deaths_oth_t", "deaths_oth_c"))   ## nsc blocks of 4 trt x 4 time point blocks 

nt <- 4 
short$scen <- factor(rep(scn, each=4*4))
short$year <- rep(c(0,1,5,10), nsc*4)
nstat <- with(procraw, qhi_stat + qlo_stat) 
naht <- with(procraw, qhi_aht + qlo_aht) 
nwr <- with(procraw, hbmi_wr + lobmi_wr) 
nsc <- with(procraw, att_sc)
nn <- cbind(nstat, naht, nwr, nsc)
short$n <- as.vector(t(nn[,rep(1:4, each=nt)]))
short <- within(short,{
    pd <- 100*(deaths_t - deaths_c)/deaths_c
    dlog <- log(deaths_t / deaths_c)
    selog <- sqrt(1/deaths_t - 2/n + 1/deaths_c)
    pdl <- 100*(exp(dlog - qnorm(0.975)*selog) - 1)
    pdu <- 100*(exp(dlog + qnorm(0.975)*selog) - 1)
})

j <- "IQALY"
resm <- round(M[,j,],1)
P <- N[,2:9] / N[,1]
P[,1:7] <- 100 * P[,1:7]
perc <- round(P, 2)
resn <- M[,"IQALY",1] / P[,"Number of HCs"]
foo <- short[short$trt=="statins" & short$year %in% c(5,10), c("scen","year","pd")]
foo <- short[short$year %in% c(5,10), c("scen","trt","year","pd")]
foo$scen <- factor(foo$scen, levels=scn)
foo2 <- spread(foo, year, pd)
resshort <- cbind(
    proc[proc$names=="QRisk>20 (first HC)","perc"],
    subset(foo2, trt=="statins")[,c("5","10")],
    proc[proc$names=="SBP>140 (first HC)","perc"],
    subset(foo2, trt=="aht")[,c("5","10")],
    proc[proc$names=="BMI>30 (first HC)","perc"],
    subset(foo2, trt=="wr")[,c("5","10")],
    proc[proc$names=="Smokers (first HC)","perc"],
    subset(foo2, trt=="sc")[,c("5","10")]
)
                        
res <- t(cbind(resm,perc,resn,resshort))

rownames(res) <- c("QALY gain (days; all)",
                   "QALY gain (days; eligible at least once)",
                   "QALY gain (days; attending at least once)",
                   "QALY gain (days; offered treatment via HC)",
                   "% eligible at least once",
                   "% attending at least once", 
                   "% offered treatment via HC" ,
                   "% taking statins via HC" ,
                   "% taking AHT via HC" ,
                   "% smoking cessation via HC" ,
                   "% weight management via HC",
                   "Mean number of HCs per person",
                   "QALY gain (days; per HC)",
                   "P(Qrisk>20) for first HC attenders",
                   "reduction in deaths; 5 years after statins",
                   "reduction in deaths; 10 years after statins",
                   "P(SBP>140) for first HC attenders",
                   "reduction in deaths; 5 years after AHT",
                   "reduction in deaths; 10 years after AHT",
                   "P(BMI>30) for first HC attenders",
                   "reduction in deaths; 5 years after WM",
                   "reduction in deaths; 10 years after WM",
                   "P(smoker) for first HC attenders",
                   "reduction in deaths; 5 years after SC",
                   "reduction in deaths; 10 years after SC"
                   )

## what is actually missing
## all outcomes for all treatments 

res <- res[c(
    "QALY gain (days; all)",
    "Mean number of HCs per person",
    "QALY gain (days; per HC)",
    "% eligible at least once",
    "QALY gain (days; eligible at least once)",
    "% attending at least once",
    "QALY gain (days; attending at least once)",
    "% offered treatment via HC" ,
    "QALY gain (days; offered treatment via HC)",
    "P(Qrisk>20) for first HC attenders",
    "% taking statins via HC" ,
    "reduction in deaths; 5 years after statins",
    "reduction in deaths; 10 years after statins",

    "P(SBP>140) for first HC attenders",
    "% taking AHT via HC",
    "reduction in deaths; 5 years after AHT",
    "reduction in deaths; 10 years after AHT",

    "P(BMI>30) for first HC attenders",
    "% weight management via HC",
    "reduction in deaths; 5 years after WM",
    "reduction in deaths; 10 years after WM",

    "P(smoker) for first HC attenders",
    "% smoking cessation via HC" ,
    "reduction in deaths; 5 years after SC",
    "reduction in deaths; 10 years after SC"),]

resf <- res; colnames(resf) <- NULL
resf <- round(resf,2)
percid <- grep("%|reduction|P\\(", rownames(res))
resf[percid,] <- paste(resf[percid,], "%", sep="")
colnames(resf) <- scn
resf

write.csv(resf, file="hc_scenarios_results2.csv", quote=FALSE)




#### BELOW HERE FOR TESTING

short <- read.csv("results/batch/scenall_short_1.csv",
                  col.names = c("trt",
                              "Q1", "Q", "CH1", "CH", "HDL1", "HDL", "SBP1", "SBP", "DBP1", "DBP", "BMI1", "BMI", "SM_t", "SM_c",
                              "Q1SD", "QSD", "CH1SD", "CHSD", "HDL1SD", "HDLSD", "SBP1SD", "SBPSD", "DBP1SD", "DBPSD", "BMI1SD", "BMISD",
                              "cvd_events_t", "cvd_events_c", "lc_events_t", "lc_events_c", "deaths_t", "deaths_c", "deaths_cvd_t", "deaths_cvd_c", "deaths_lc_t", "deaths_lc_c", "deaths_dem_t", "deaths_dem_c", "deaths_oth_t", "deaths_oth_c"), header=FALSE)
short$scen <- factor(rep(scn, each=4*4))
short$year <- rep(c(0,1,5,10), nsc*4)


short[short$trt=="sc" & short$scen=="Base case", ]

read.csv("data/lifetable_ihd.csv", skip=2)[40:60,] # men/women
read.csv("data/lifetable_lungcancer.csv",skip=2)[40:60,]
