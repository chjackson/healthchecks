### COMBINE RESULTS FROM MULTIPLE SMALL POPS INTO ONE BIG ONE 

library(tidyverse)

nruns <- 100
aname <- "results/oct16/scen"
nsc <- 16; 

#nruns <- 1000
#aname <- "~/scratch/hc/healthchecks/results/unc/unc"
#nsc <- 1; 

### TODO REMOVE lobmi_wr in future runs and edit pnms

nouts <- 6; npops <- 8
outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY")

mn <- c("All", "Eligible", "Attending", "Any treatment", "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")

scn <- c("Base case","BP", "Attend age 50-74", "Attend age 40-80", "Attend age 50-80", "Baseline uptake 0.85 (from 0.48)","Uptake 10 - 20% more in most deprived","Uptake of smokers +30%","+20% uptake for QRisk > 15","Statin prescription x 2.5", "AHT prescription x 2.5", "Smoking referral x 2.5", "Weight referral x 2.5", "All treatments x 2.5", "Attenders keep attending", "Target non-attenders")[1:nsc]

M <- S <- N <- NT <- Npop <- Tot <- Sumsq <- P <- SH <- vector(nruns, mode="list")
Marr <- array(dim=c(nruns, nsc, nouts, npops)) # different format, for statistical uncertainty
Parr <- array(dim=c(nruns, nsc, npops)) # different format, for statistical uncertainty

for (i in 1:nruns){
    M[[i]] <- array(as.matrix(read.csv(sprintf("%s_mean_%s.csv",aname,i),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
    Marr[i,,,] <- array(as.matrix(read.csv(sprintf("%s_mean_%s.csv",aname,i),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
    S[[i]] <- array(as.matrix(read.csv(sprintf("%s_SD_%s.csv",aname,i),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
    N[[i]] <- as.matrix(read.csv(sprintf("%s_N_%s.csv",aname,i),header=FALSE))
    NT[[i]] <- as.matrix(read.csv(sprintf("%s_NT_%s.csv",aname,i),header=FALSE))
    if (nsc==1){
        Parr[i,1,] <- N[[i]][,2:9,drop=FALSE] / N[[i]][,1]
    }
    dimnames(N[[i]]) <- list(scn, mn)
    Npop[[i]] <- N[[i]][,1:npops,drop=FALSE]
    Npop[[i]][,5:8] <-  NT[[i]]
    Tot[[i]] <- Sumsq[[i]] <- array(dim=dim(M[[i]]))
    for (j in 1:nouts){
        Tot[[i]][,j,] <- M[[i]][,j,] * Npop[[i]]
        Sumsq[[i]][,j,] <- S[[i]][,j,] * Npop[[i]]
    }
}

Nall <- Reduce("+", N)
Npopall <- Reduce("+", Npop)
NTall <- Reduce("+", NT)
Tot <- Reduce("+", Tot)
Sumsq <- Reduce("+", Sumsq)

Mall <- Sall <- array(dim=dim(Tot))
for (j in 1:nouts){ 
    Mall[,j,] <- Tot[,j,] / Npopall
    Sall[,j,] <- Sumsq[,j,] / Npopall
}

M <- Mall; S <- Sall; Npop <- Npopall
dimnames(M) <- dimnames(S) <- list(scn, outn, mn[1:npops])
dimnames(Npop) <- list(scn, mn[1:8])
SE <- L <- U <- S
for (j in 1:nouts){
    SE[,j,] <- S[,j,] / sqrt(Npop)
    L[,j,] <- M[,j,] - qnorm(0.975)*SE[,j,]
    U[,j,] <- M[,j,] + qnorm(0.975)*SE[,j,]
}
M[,3,]
SE[,3,]
SE[,3,] / M[,3,] ## SE is 5% of mean after 12 runs

if (nsc==1) { 
### Statistical standard error (posterior standard deviation)
SD <- sqrt(apply(Marr, c(2,3,4), var))
## Bias-corrected version, v little difference from npop 25k x npar 1000 run 
# SDA <- sqrt(SD^2 - SE^2)
## todo comp SE 
SEbig <- sqrt((SD^2 + SE^2)/nruns)
(secomp <- cbind(SD[1,3,], SEbig[1,3,])) # comp SE is 3% of the stat SD, satisfying bugs manual criterion, so neglect comp SE
secomp[,2] / secomp[,1]
PSD <- sqrt(apply(Parr, c(2,3), var))
}


mn <- c("Q1", "Q", "CH1", "CH", "HDL1", "HDL", "SBP1", "SBP", "DBP1", "DBP", "BMI1", "BMI")
sn <- c("Q1SD", "QSD", "CH1SD", "CHSD", "HDL1SD", "HDLSD", "SBP1SD", "SBPSD", "DBP1SD", "DBPSD", "BMI1SD", "BMISD")
nn <- c("n","SM_t","SM_c","cvd_events_t", "cvd_events_c", "lc_events_t", "lc_events_c", "deaths_t", "deaths_c", "deaths_cvd_t", "deaths_cvd_c", "deaths_lc_t", "deaths_lc_c", "deaths_dem_t", "deaths_dem_c", "deaths_oth_t", "deaths_oth_c")

P <- Pdenom <- Ptot <- SH <- SHtot <- SHsumsq <- vector(nruns, mode="list")

SHarr <- array(dim=c(nruns, 16, 50))
Prarr <- array(dim=c(nruns, 37))

for (i in 1:nruns){
    P[[i]] <- read.csv(sprintf("%s_process_%s.csv",aname,i), header=FALSE,
                       col.names=c("popsize","eligible","attending","att_q","att_chol","att_hdl","qhi","qhi_presc_stat","qhi_stat","qlo","qlo_presc_stat","qlo_stat","nonatt","nonatt_q","nonatt_chol","nonatt_h","att_sbp","att_dbp","att_nhbp","qhi_presc_aht","qhi_aht","qlo_presc_aht","qlo_aht","nonatt_sbp","nonatt_dbp", "att_bmi","hbmi","hbmi_presc_wr","hbmi_wr","lobmi","lobmi_presc_wr","lobmi_wr","nonatt_bmi","att_smk","att_sc","att_quit1","att_nsmk"))
    SH[[i]] <- read.csv(sprintf("%s_short_%s.csv",aname,i), header=FALSE,
                        col.names=c("trt",
                                    "Q1", "Q", "CH1", "CH", "HDL1", "HDL", "SBP1", "SBP", "DBP1", "DBP", "BMI1", "BMI", "SM_t", "SM_c",
                                    "Q1SD", "QSD", "CH1SD", "CHSD", "HDL1SD", "HDLSD", "SBP1SD", "SBPSD", "DBP1SD", "DBPSD", "BMI1SD", "BMISD",
                                    "cvd_events_t", "cvd_events_c", "lc_events_t", "lc_events_c", "deaths_t", "deaths_c", "deaths_cvd_t", "deaths_cvd_c", "deaths_lc_t", "deaths_lc_c", "deaths_dem_t", "deaths_dem_c", "deaths_oth_t", "deaths_oth_c"))
    SH[[i]]$year <- rep(c(0,1,5,10), nsc*4)

    cs <- c("att_q","att_chol","att_hdl", # denom: total number of HCs.
            "nonatt_q","nonatt_chol","nonatt_h", # denom: number not attending
            "att_sbp","att_dbp",
            "nonatt_sbp","nonatt_dbp",
            "att_bmi","nonatt_bmi"
            )
    N[[i]] <- cbind(N[[i]], "Not attending" = N[[i]][,"All"] - N[[i]][,"Attending"])
    Pdenom[[i]] <- N[[i]][,c(rep("Number of HCs",3),rep("Not attending",3),
                             rep("Number of HCs",2),rep("Not attending",2),
                             rep("Number of HCs",1),rep("Not attending",1))]
    Ptot[[i]] <- P[[i]][,cs]*Pdenom[[i]]

    ## number getting treatment at first HC and adherent (stat,aht,wr), 
    NTrep <- NT[[i]][,rep(c(1,2,4), each=4),drop=FALSE]
    ## number referred to SC  (NT stores quitters)
    Psmrep <- matrix(P[[i]][,"att_smk"], ncol=1)[,rep(1,4),drop=FALSE]
    SH[[i]]$n <- as.vector(t(cbind(NTrep, Psmrep)))
    SHtot[[i]] <- SH[[i]][,mn]*SH[[i]]$n
    SHsumsq[[i]] <- SH[[i]][,sn]*SH[[i]]$n

    ncounts <- c("cvd_events","lc_events","deaths","deaths_cvd","deaths_lc","deaths_dem","deaths_oth","SM")
    for (nam in ncounts){ 
        ds <- SH[[i]][,paste0(nam,"_t")]
        dc <- SH[[i]][,paste0(nam,"_c")]
        SH[[i]][,paste0(nam,"_pd")] <- 100*(ds - dc)/dc
    }
    if(nsc==1){
        SHarr[i,,] <- as.matrix(select(SH[[i]], -trt))
        Prarr[i,] <- as.matrix(P[[i]])
    }
}
if (nsc==1)
    dimnames(SHarr)[2:3] <- dimnames(select(SH[[1]], -trt))


ns <- c("popsize","eligible","attending","qhi","qhi_presc_stat","qhi_stat","qlo","qlo_presc_stat","qlo_stat","nonatt","att_nhbp","qhi_presc_aht","qhi_aht","qlo_presc_aht","qlo_aht","hbmi","hbmi_presc_wr","hbmi_wr","lobmi","lobmi_presc_wr","lobmi_wr","att_smk","att_sc","att_quit1","att_nsmk")
Pns <- lapply(P, function(x)x[,ns])
Pns <- Reduce("+", Pns)
Ptot <- Reduce("+", Ptot) / Reduce("+", Pdenom)
Pall <- array(dim=dim(P[[1]]), dimnames=dimnames(P[[1]]))
Pall[,ns] <- as.matrix(Pns)
Pall[,cs] <- as.matrix(Ptot)
## Pall is equivalent of procraw for usage in scenarios.r 

## Stat SE around process outcomes
## todo apply to percentages, but this needs percs to be calculated for each i
## Same with short term outcomes.
## Wouldn't all this be better in Python 
if (nsc==1) PrSD <- apply(Prarr, 2, sd)

pnms <- c("Population", "Eligible at least once", "Attending at least once",
            "QRisk (attenders; over all HCs)","Cholesterol (attenders; over all HCs)","HDL (attenders; over all HCs)",
            "QRisk>20 (first HC)", "Prescribed statins (Q>20)","Statins (Q>20)",
            "QRisk<20 (first HC)", "Prescribed statins (Q<20)","Statins (Q<20)",
            "Non-attending", "QRisk (non-attenders)","Cholesterol (non-attenders)","HDL (non-attenders)",
            "SBP (attenders; all HCs)", "DBP (attenders; all HCs)", "SBP>140 (first HC)",
            "Prescribed AHT (Q>20)", "AHT (Q>20)",  "Prescribed AHT (Q<20)", "AHT (Q<20)",
            "SBP (non-attenders; all HCs)", "DBP (non-attenders; all HCs)",
            "BMI (attenders; all HCs)", "BMI>30 (first HC)", "Referred to weight management", "Weight management", 
            "BMI<30 (first HC)", "Referred to weight management (BMI<30)", "Weight management (BMI<30)", 
            "BMI(non-attenders)",
            "Smokers (first HC)", "Referred to smoking cessation", "Quit by 1 year due to cessation service", "Non-smokers")

denoms <- read.table("denoms.txt", as.is=TRUE, header=TRUE)
Pall <- as.data.frame(Pall)
Pall$scen <- scn

proc <- 
  gather(as.data.frame(Pall), outcome, num, -scen) %>%
  group_by(scen) %>%
  mutate(denomname = denoms$denom[match(outcome, denoms$num)],
         denom = num[match(denomname, outcome)],
         perc=100*num/denom       ,
         seperc = 100*sqrt((num/denom)*(1 - num/denom)/denom),
         numf = ifelse(is.na(denom), round(num,2), paste0(num,"/",denom)),
         percf = ifelse(is.na(denom), NA, paste0(round(perc, 2),
                                                 "% (",round(seperc,2), ")")), 
         pnames = pnms)

SHnn <- lapply(SH, function(x)x[,nn])
SHnn <- Reduce("+", SHnn)
SHtot <- Reduce("+", SHtot)
SHmean <- SHtot/SHnn$n
SHsumsq <- Reduce("+", SHsumsq)
SHsd <- SHsumsq/SHnn$n
## short term results for big population : equivalent of short in scenarios.r
short <- cbind(SHmean, SHsd, SHnn)
short$scen <- factor(rep(scn, each=4*4))
short$year <- rep(c(0,1,5,10), nsc*4)
short$trt <- rep(c("statins","aht","wr","sc"), each=4)

## Percentages and computational uncertainty intervals around short term count outcomes 
### var (exp(x) ) = exp(mu)^2 var(x)
### se(exp(x)) = exp(mu) se(x) 

cuncp <- function(short, nam, denom=short$n){
    ds <- short[,paste0(nam,"_t")]
    dc <- short[,paste0(nam,"_c")]
    pd <- 100*(ds - dc)/dc
    dlog <- log(ds / dc)
    selog <- sqrt(1/ds - 2/short$n + 1/dc)
    pdse <- 100*exp(dlog)*selog
    pdl <- 100*(exp(dlog - qnorm(0.975)*selog) - 1)
    pdu <- 100*(exp(dlog + qnorm(0.975)*selog) - 1)
    short[,paste0(nam,"_pd")] <- pd
    short[,paste0(nam,"_pdse")] <- pdse
    short[,paste0(nam,"_pdl")] <- pdl
    short[,paste0(nam,"_pdu")] <- pdu
    short
}
nn <- c("cvd_events","lc_events","deaths","deaths_cvd","deaths_lc","deaths_dem","deaths_oth","SM")
for (i in nn)
    short <- short %>% cuncp(i)

## Stat standard errors around short term count outcomes 
if(nsc==1){
    SHmean <- apply(SHarr[,,paste0(nn, "_pd")], c(2,3), mean)
    SHSD <- apply(SHarr[,,paste0(nn, "_pd")], c(2,3), sd)
}
## compare with computational: e.g. 0.8, 0.5 for 25000*100, but with 1000, should be negligible comp unc 

### Though for events - isn't it poisson not binomial

#SH[[1]][SH[[1]]$trt=="wr" & SH[[1]]$year==10, ]
#short %>% filter(trt=="sc" & scen=="Base case") %>% select(n, SM_t, SM_c, SM_pd, SM_pdl, SM_pdu)
