### COMBINE RESULTS FROM MULTIPLE SMALL POPS INTO ONE BIG ONE 

library(tidyverse)
nruns <- 50
aname <- "base/scenall"
nsc <- 1; nouts <- 6; npops <- 8

outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY")
mn <- c("All", "Eligible", "Attending", "Any treatment", "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
scn <- c("Base case","BP", "Diabetes","Diabetes+BP", "Attend from age 50","Baseline uptake 0.85 (from 0.48)","Uptake (+10)/+20% in most deprived","Uptake of smokers +30%","+20% uptake for QRisk > 15","Statin prescription x 2.5", "Attenders keep attending", "Target non-attenders")[1:nsc]

M <- S <- N <- NT <- Npop <- Tot <- Sumsq <- P <- SH <- vector(nruns, mode="list")

for (i in 1:nruns){
    M[[i]] <- array(as.matrix(read.csv(sprintf("results/%s_mean_%s.csv",aname,i),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
    S[[i]] <- array(as.matrix(read.csv(sprintf("results/%s_SD_%s.csv",aname,i),header=FALSE)), dim=c(nsc, nouts, npops))*365.25
    N[[i]] <- as.matrix(read.csv(sprintf("results/%s_N_%s.csv",aname,i),header=FALSE))
    NT[[i]] <- as.matrix(read.csv(sprintf("results/%s_NT_%s.csv",aname,i),header=FALSE))
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


mn <- c("Q1", "Q", "CH1", "CH", "HDL1", "HDL", "SBP1", "SBP", "DBP1", "DBP", "BMI1", "BMI")
sn <- c("Q1SD", "QSD", "CH1SD", "CHSD", "HDL1SD", "HDLSD", "SBP1SD", "SBPSD", "DBP1SD", "DBPSD", "BMI1SD", "BMISD")
nn <- c("n","SM_t","SM_c","cvd_events_t", "cvd_events_c", "lc_events_t", "lc_events_c", "deaths_t", "deaths_c", "deaths_cvd_t", "deaths_cvd_c", "deaths_lc_t", "deaths_lc_c", "deaths_dem_t", "deaths_dem_c", "deaths_oth_t", "deaths_oth_c")

P <- Pdenom <- Ptot <- SH <- SHtot <- SHsumsq <- vector(nruns, mode="list")
for (i in 1:nruns){
    P[[i]] <- read.csv(sprintf("results/%s_process_%s.csv",aname,i), header=FALSE,
                       col.names=c("popsize","eligible","attending","att_q","att_chol","att_hdl","qhi","qhi_presc_stat","qhi_stat","qlo","qlo_presc_stat","qlo_stat","nonatt","nonatt_q","nonatt_chol","nonatt_h","att_sbp","att_dbp","att_nhbp","qhi_presc_aht","qhi_aht","qlo_presc_aht","qlo_aht","nonatt_sbp","nonatt_dbp", "att_bmi","hbmi","hbmi_presc_wr","hbmi_wr","lobmi","lobmi_presc_wr","lobmi_wr","nonatt_bmi","att_smk","att_sc","att_quit1","att_nsmk"))
    SH[[i]] <- read.csv(sprintf("results/%s_short_%s.csv",aname,i), header=FALSE,
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

    ## number getting treatment at first HC and adherent (stat,aht,wr), or number referrd to SC.
    sdenoms <- with(P[[i]], as.data.frame(cbind(statins = qhi_stat+qlo_stat, aht=qhi_aht+qlo_aht, wr=hbmi_presc_wr, sc=att_sc)))
    SH[[i]]$n <- rep(as.vector(t(sdenoms)), each=4)
    SHtot[[i]] <- SH[[i]][,mn]*SH[[i]]$n
    SHsumsq[[i]] <- SH[[i]][,sn]*SH[[i]]$n
}


ns <- c("popsize","eligible","attending","qhi","qhi_presc_stat","qhi_stat","qlo","qlo_presc_stat","qlo_stat","nonatt","att_nhbp","qhi_presc_aht","qhi_aht","qlo_presc_aht","qlo_aht","hbmi","hbmi_presc_wr","hbmi_wr","lobmi","lobmi_presc_wr","lobmi_wr","att_smk","att_sc","att_quit1","att_nsmk")
Pns <- lapply(P, function(x)x[,ns])
Pns <- Reduce("+", Pns)
Ptot <- Reduce("+", Ptot) / Reduce("+", Pdenom)
Pall <- array(dim=dim(P[[1]]), dimnames=dimnames(P[[1]]))
Pall[,ns] <- as.matrix(Pns)
Pall[,cs] <- as.matrix(Ptot)
## Pall is equivalent of procraw for usage in scenarios.r 


pnms <- c("Population", "Eligible at least once", "Attending at least once",
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

denoms <- read.table("denoms.txt", as.is=TRUE, header=TRUE)
Pall <- as.data.frame(Pall)
Pall$scen <- scn
proc <- 
  gather(as.data.frame(Pall), outcome, num, -scen) %>%
  group_by(scen) %>%
  mutate(denomname = denoms$denom[match(outcome, denoms$num)],
         denom = num[match(denomname, outcome)],
         perc=100*num/denom,
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

## Computational uncertainty intervals around short term count outcomes 
cuncp <- function(short, nam, denom=short$n){
    ds <- short[,paste0(nam,"_t")]
    dc <- short[,paste0(nam,"_c")]
    pd <- 100*(ds - dc)/dc
    dlog <- log(ds / dc)
    selog <- sqrt(1/ds - 2/short$n + 1/dc)
    pdl <- 100*(exp(dlog - qnorm(0.975)*selog) - 1)
    pdu <- 100*(exp(dlog + qnorm(0.975)*selog) - 1)
    short[,paste0(nam,"_pd")] <- pd
    short[,paste0(nam,"_pdl")] <- pdl
    short[,paste0(nam,"_pdu")] <- pdu
    short
}
nn <- c("cvd_events","lc_events","deaths","deaths_cvd","deaths_lc","deaths_dem","deaths_oth","SM")
for (i in nn)
    short <- short %>% cuncp(i)

###SHnn[SH[[1]]$trt=="statins" & SH[[1]]$year==10,
###      c("n","deaths_s","deaths_c","pd","pdl","pdu")]

short %>% filter(trt=="sc" & scen=="Base case") %>% select(n, SM_t, SM_c, SM_pd, SM_pdl, SM_pdu)
