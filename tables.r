library(tidyverse)
library(forcats)

############################################################
###  TABLES OF SCENARIOS 
############################################################

## Scenarios
scn <- c("Base case","Invite high BP", "Attend age 50-74", "Attend age 40-80", "Attend age 50-80",
         "Baseline uptake +30%","Uptake +30% in most deprived","Uptake of smokers +30%","+30% uptake for QRisk > 20",
         "Target non-attenders", "Attenders keep attending", "Attenders keep attending, target non-attenders", 
         "Statin prescription x 2.5", "AHT prescription x 2.5", "Smoking referral x 2.5", "Weight referral x 2.5",
         "All treatments x 2.5", "Higher treatment, higher uptake and invite BP",
         "CVD incidence declining", "Baseline QRisk uncertainty +- 20%")
# scn <- scn[1:13] # c("Base case","Invite high BP")
nsc <- length(scn)

## Model outputs 
evn <- c("IHD","STR","DEM","LC") # counts of disease events
outn.ev <- paste(rep(evn, each=3), rep(c(80,100), each=length(evn)*3), paste(rep(c("_con","_hc","_i"), length(evn)),sep=""),sep="")
dn <- c("D")                     # counts of deaths 
outn.d <- paste(rep(dn, each=3), rep(c(75,80), each=length(dn)*3), paste(rep(c("_con","_hc","_i"), length(dn)),sep=""),sep="")
outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY","IQALYperhc", "ILYperhc", outn.ev,outn.d)
nouts <- length(outn)

## Population subsets for which results are calculated
pops <- c("All", "Eligible", "Attending", "Any treatment", "Most deprived", "Least deprived")
npops <- length(pops)
popsall <- c(pops, 
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)",
          "Attending and eligible for any treatment (Q20)", "Attending and eligible for any treatment (Q10)",
          "Offered statins via HC", "Offered AHT via HC", "Offered SC via HC", "Offered WR via HC", 
          "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
npopsall <- length(popsall)
postn <- c("QRisk","SBP","DBP")
npostn <- length(postn)

## Read results with or without statistical uncertainty
statunc <- TRUE
# statunc <- FALSE

## Read and combine aggregate results from "nruns" batches 
## rearranging Python np.savetxt CSV output into properly-indexed arrays 
fname <- if (statunc) "~/scratch/hc/healthchecks/results/paperrev/scen" else "~/scratch/hc/healthchecks/results/paperrev/scen"
resm <- as.matrix(read.table(sprintf("%s_mean.csv",fname), colClasses="numeric", sep=",", header=FALSE))
## SD of output over population
ress <- as.matrix(read.table(sprintf("%s_SD.csv",fname), colClasses="numeric", sep=",", header=FALSE))
n <- as.matrix(read.table(sprintf("%s_N.csv",fname), colClasses="numeric", sep=",", header=FALSE))[,(1:npopsall)+1]  # nruns*nsc  x  npopsall
runid <- resm[,1]
nruns <- if (statunc) length(unique(runid)) else 40  # ~/scratch/hc/healthchecks/
resm <- resm[,-1]; ress <- ress[,-1]
Mrep <- array(resm, dim=c(nsc, nruns, nouts, npops))
Srep <- array(ress, dim=c(nsc, nruns, nouts, npops))
dimnames(Mrep) <- dimnames(Srep) <- list(scn, 1:nruns, outn, pops)
Nrep <- array(n, dim=c(nsc, nruns, npopsall))
dimnames(Nrep) <- list(scn, 1:nruns, popsall)
Nrep4D <- array(Nrep, dim=c(nsc, nruns, 1, npopsall))[,,rep(1,nouts),] # replicate to match dimensions of Srep
dimnames(Nrep4D) <- list(scn, 1:nruns, outn, popsall)
n <- apply(Nrep4D, c(1,3,4), sum) # total pop size over batches
Npop <- n[,,1:npops] 
Nall <- n[,1,]
## convert events averted from proportion to "per million population"
Mrep[,,c(outn.d, outn.ev),] <- Mrep[,,c(outn.d, outn.ev),]*1000000
Srep[,,c(outn.d, outn.ev),] <- Srep[,,c(outn.d, outn.ev),]*1000000

## examine some results 
Mrep["Base case",,"IQALY","All"]
t(apply(Mrep[,,"IQALY","All"], 1, quantile, c(0.025, 0.5, 0.975)))
t(apply(Mrep[,,"ILY","All"], 1, quantile, c(0.025, 0.5, 0.975)))
t(apply(Mrep[,,"IHD80_i","All"], 1, quantile, c(0.025, 0.5, 0.975)))

Mrep["CVD incidence declining",,"IQALY","All"]
Mrep["Baseline QRisk uncertainty +- 20%",,"IQALY","All"]


if (!statunc) {
    Msum <- apply(Mrep*Nrep4D[,,,1:npops], c(1,3,4), sum) # total 
    M <- Msum / Npop # mean after combining all batches
    Sumsq <- apply(Srep^2*Nrep4D[,,,1:npops], c(1,3,4), sum) # total squared error summed over batches 
    Sall <- sqrt(Sumsq / Npop) # SD after combining all batches
    SE <- Sall / sqrt(Npop)    # Monte Carlo SE of mean 
} else if (statunc) {  ### Results with statistical uncertainty
    M <- apply(Mrep, c(1,3,4), mean) # mean of means over all batches
    n <- apply(Nrep, c(1,3), mean) # mean within-simulation pop size
    n <- array(n, dim=c(nsc,1,npopsall))[,rep(1,nouts),,drop=FALSE] # match dims of M, S
    dimnames(n) <- list(scn, outn, popsall)
    nmu <- n[,,1:npops]
    sigsq <- apply(Mrep, c(1,3,4), var)
    LCL <- apply(Mrep, c(1,3,4), function(x)quantile(x, 0.025))
    UCL <- apply(Mrep, c(1,3,4), function(x)quantile(x, 0.975))
    ## Smooth the confidence limits for incremental event counts
    ## (otherwise will come out as multiples of 40, as counts from runs of population size 25000 scaled up to population size 1 million 
    qfn <- function(x, p){
        set.seed(1)
        quantile(x + rnorm(length(x), 0, sd(x)/10), p)
    }
    LCLN <- apply(Mrep, c(1,3,4), function(x)qfn(x, 0.025))
    UCLN <- apply(Mrep, c(1,3,4), function(x)qfn(x, 0.975))
    tausq <- apply(Srep^2, c(1,3,4), mean)
    ## Monte Carlo standard error around the mean
    semu <- sqrt(sigsq / nruns + tausq/(nruns*nmu))
    ## Statistical uncertainty / posterior standard deviation 
    SE <- sqrt(sigsq) # standard MC est is biased
    ## Monte Carlo standard error around the variance (SD^2)
    sevar <- sqrt(2/(nruns-1))*(sigsq + tausq/nmu)
    ## Monte Carlo standard error around the SD (using delta method)
    sesd <- 0.5*sevar/SE
    ## Bias-corrected estimates, may not work for smaller population sizes, and big between-person variability (see O'Hagan et al, Health Econ 2007, p1014)
    sigsqA <- sigsq - tausq/nmu # ...as this could be negative
    SEA <- sqrt(sigsqA)
    sevarA <- sqrt(2*((sigsqA + tausq/nmu)^2 / (nruns-1) + tausq^2 /(nruns*nmu^2*(nmu-1))))
    seSEA <- 0.5*sevarA/SEA
    LCLA <- M - qnorm(0.975)*SEA
    UCLA <- M + qnorm(0.975)*SEA
    seratio <- semu / SE
    res <- array(c(M, SE, SEA, semu, seratio, LCL, UCL, LCLN, UCLN, LCLA, UCLA), dim=c(dim(M), 11),
                 dimnames=c(dimnames(M), list(c("M","SE","SEA","semu","seratio","LCL","UCL","LCLN","UCLN","LCLA","UCLA"))))
    Nall <- n[,1,]
}


res[,"IQALY","All",]
res[,"IQALY","All",c("SE","SEA","LCL","UCL","LCLA","UCLA")] # bias-corrected SD estimator doesn't work unfortunately in the case where it would be useful (i.e. where the lower confidence limit for health gains from one scenario is implausibly negative)
### Computational error is about 10-15% of statistical error.
res[,,"All","seratio"]
1

rows.ev <- c("IHD80_i",  "STR80_i", "DEM80_i",  "LC80_i", 
          "IHD100_i", "STR100_i", "DEM100_i", "LC100_i",
          "D75_i", "D80_i")   
rows.ly <- c("IQALY", "ILY","IQALYperhc", "ILYperhc")
resev <- round(t(M[,rows.ev,"All"]), 0)
resev <- array(as.character(resev), dim=dim(resev), dimnames=dimnames(resev))
resly <- round(t(M[,rows.ly,"All"]), 1)
resly <- array(as.character(resly), dim=dim(resly), dimnames=dimnames(resly))
res_el <- t(M[,c("IQALY","ILY"),"Eligible"]); rownames(res_el) <- c("IQALY_el","ILY_el")
res_att <- t(M[,c("IQALY","ILY"),"Attending"]); rownames(res_att) <- c("IQALY_att","ILY_att")
res_dep <- t(M[,c("IQALY","ILY"),"Most deprived"]); rownames(res_dep) <- c("IQALY_dep","ILY_dep")
res_ndep <- t(M[,c("IQALY","ILY"),"Least deprived"]); rownames(res_ndep) <- c("IQALY_ndep","ILY_ndep")
ressubpop <- round(rbind(res_el, res_att, res_dep, res_ndep), 1) 
ressubpop <- array(as.character(ressubpop), dim=dim(ressubpop), dimnames=dimnames(ressubpop))
resm <- rbind(resev, resly, ressubpop)
resn <- t(Nall)
resperc <- 100*resn/resn[1,1]
resperc["Number of HCs",] <- resperc["Number of HCs",]/100
resperc <- round(resperc,1)
resperc <- array(as.character(resperc), dim=dim(resperc), dimnames=dimnames(resperc))

## Paste standard errors or CIs for life / event count gains 
sem <- t(SE[,rows.ly,"All"])
se_el <- t(SE[,c("IQALY","ILY"),"Eligible"]); rownames(se_el) <- c("IQALY_el","ILY_el")
se_att <- t(SE[,c("IQALY","ILY"),"Attending"]); rownames(se_att) <- c("IQALY_att","ILY_att")
se_dep <- t(SE[,c("IQALY","ILY"),"Most deprived"]); rownames(se_dep) <- c("IQALY_dep","ILY_dep")
se_ndep <- t(SE[,c("IQALY","ILY"),"Least deprived"]); rownames(se_ndep) <- c("IQALY_ndep","ILY_ndep")
sem <- round(rbind(sem, se_el, se_att, se_dep, se_ndep), 1) 
sem <- array(as.character(sem), dim=dim(sem), dimnames=dimnames(sem))

lclm <- t(LCL[,rows.ly,"All"])
lcl_el <- t(LCL[,c("IQALY","ILY"),"Eligible"]); rownames(lcl_el) <- c("IQALY_el","ILY_el")
lcl_att <- t(LCL[,c("IQALY","ILY"),"Attending"]); rownames(lcl_att) <- c("IQALY_att","ILY_att")
lcl_dep <- t(LCL[,c("IQALY","ILY"),"Most deprived"]); rownames(lcl_dep) <- c("IQALY_dep","ILY_dep")
lcl_ndep <- t(LCL[,c("IQALY","ILY"),"Least deprived"]); rownames(lcl_ndep) <- c("IQALY_ndep","ILY_ndep")
lclm <- round(rbind(lclm, lcl_el, lcl_att, lcl_dep, lcl_ndep), 1) 
lclm <- array(as.character(lclm), dim=dim(lclm), dimnames=dimnames(lclm))

uclm <- t(UCL[,rows.ly,"All"])
ucl_el <- t(UCL[,c("IQALY","ILY"),"Eligible"]); rownames(ucl_el) <- c("IQALY_el","ILY_el")
ucl_att <- t(UCL[,c("IQALY","ILY"),"Attending"]); rownames(ucl_att) <- c("IQALY_att","ILY_att")
ucl_dep <- t(UCL[,c("IQALY","ILY"),"Most deprived"]); rownames(ucl_dep) <- c("IQALY_dep","ILY_dep")
ucl_ndep <- t(UCL[,c("IQALY","ILY"),"Least deprived"]); rownames(ucl_ndep) <- c("IQALY_ndep","ILY_ndep")
uclm <- round(rbind(uclm, ucl_el, ucl_att, ucl_dep, ucl_ndep), 1) 
uclm <- array(as.character(uclm), dim=dim(uclm), dimnames=dimnames(uclm))

rows.ly <- rownames(uclm)
resm.ly <- array(paste0(resm[rows.ly,], " (", lclm, ", ", uclm, ")"),
                 dim=dim(resm[rows.ly,]),  dimnames=dimnames(resm[rows.ly,]))
se_ev <- round(t(SE[,rows.ev,"All"]), 0)
lcl_ev <- round(t(LCLN[,rows.ev,"All"]), 0)
ucl_ev <- round(t(UCLN[,rows.ev,"All"]), 0)
resm.ev <- array(paste0(resm[rows.ev,], " (", lcl_ev, ", ", ucl_ev, ")"),
                 dim=dim(resm[rows.ev,]),  dimnames=dimnames(resm[rows.ev,]))

heads <- c("Eligibility and uptake", "Treatment",
           "Cases prevented by age 80",
           "Cases prevented",
           "Deaths prevented",
           "Change in QALY", "Change in lifetime")
reshead <- array(dim=c(length(heads), nsc))
rownames(reshead) <- heads

rows <- c("Eligibility and uptake",
          "Eligible","Attending","Number of HCs",
          "Treatment",
          "Eligible for any treatment (Q10)",
          "Attending and eligible for any treatment (Q10)",
          "Any treatment",
          "Offered statins via HC", "Offered AHT via HC", "Offered WR via HC", "Offered SC via HC", 
          "Statins via HC","Antihypertensives via HC","Weight management via HC","Smoking cessation via HC",
          "Cases prevented by age 80",
          "IHD80_i",  "STR80_i", "DEM80_i",  "LC80_i", 
          "Cases prevented",
          "IHD100_i", "STR100_i", "DEM100_i", "LC100_i",
          "Deaths prevented",
          "D75_i", "D80_i",
          "Change in QALY",
          "IQALY","IQALY_el","IQALY_att","IQALYperhc","IQALY_dep","IQALY_ndep",
          "Change in lifetime",
          "ILY","ILY_el","ILY_att","ILYperhc","ILY_dep","ILY_ndep")

res <- rbind(resperc, resm.ev, resm.ly, reshead)[rows,]
options(width=180)
res[,1:5]
res[,6:10]
res[,11:12]
res[,13:18]

## text file to paste into Word and convert to table. 
write.table(res, "~/hc/paper/tab_res.txt", quote=FALSE, sep="\t", row.names=TRUE, na="")

### Absolute numbers of outcomes under base case (in paper text but not table)
round(M["Base case",c("IHD100_con","STR100_con","DEM100_con","LC100_con","D80_con"), "All"], -3)/1000


## Post-HC outputs (in paper text but not table)
resp <- as.matrix(read.table(sprintf("%s_post.csv",fname), colClasses="numeric", sep=",", header=FALSE))[,-1]
Prep <- array(resp, dim=c(nsc, nruns, 2, npostn))
PM <- Prep[,,1,]; PS <- Prep[,,2,]
PM <- apply(PM, c(1,3), mean); rownames(PM) <- scn
Psumsq <- apply(PS^2 * Nrep4D[,,1:3,"Attending"], c(1,3), sum)
PS <- sqrt(Psumsq / Npop[,1:3,"Attending"])
PSE <- PS / sqrt(Npop[,1:3,"Attending"])
PM["Base case",]   ### Chang results 
PSE["Base case",]


## Expected value of partial perfect information for each parameter
# library(earth) ## earth method doesn't work well here. perhaps since many low EVPPIs?
library(mgcv)
library(mvtnorm)
B <- 500
source("parnames.r")

### uptake rates and logeff not uncertain 

pars <- as.matrix(read.table("~/scratch/hc/healthchecks/results/paperjun/unc_pars.csv", colClasses="numeric",header=TRUE,sep=","))
## Note the above file might need to be edited by hand so the column names are placed in the first row - due to parallel processing, run number 1 might not have been the first one to complete and write to this file!
npars <- ncol(pars)
outs <- t(Mrep[c("Base case","Higher treatment, higher uptake and invite BP","Statin prescription x 2.5","AHT prescription x 2.5","Smoking referral x 2.5","Weight referral x 2.5"),, "IQALY","All"]) # headline output 
# outs <- matrix(Mrep[c("Base case"),, "IQALY","All"], ncol=1) # headline output 
pevppi <- ps <- matrix(nrow=npars, ncol=ncol(outs)) # could extend to different outputs 
dimnames(pevppi) <- list(colnames(pars), colnames(outs))
calc.se <- TRUE
for(i in 1:npars){
    print(i)
    x <- pars[,i]
    for (j in 1:ncol(outs)) {
        y <- outs[,j]
        mod <- gam(y ~ s(x, bs="cr"))
        pevppi[i,j] <- var(mod$fitted) / var(y)
        if (calc.se) { 
            P <- predict.gam(mod, type="lpmatrix")
            frep <- rmvnorm(B, mod$fitted, P %*% mod$Vp %*% t(P), method="svd")
            pevppirep <- apply(frep, 1, var) / var(y)  ## EVPPI as prop of EVPI
            ps[i,j] <- sd(pevppirep) # SD of EVPPI estimate
        }
    }
#        fitted(earth(x, y))) / var(y)  ## EVPPI as prop of EVPI
}
pevppi[pevppi<1e-04] <- 0
eres <- data.frame(var=colnames(pars), pevppi=pevppi)
eres %>% arrange(desc(pevppi[,1]))

if (calc.se) { 
    eres <- data.frame(var=colnames(pars), pevppi=pevppi, ps=ps)
    eres %>% arrange(desc(pevppi[,1])) # most SEs same mag as error, so ests practically zero
}

## Uncertainty analysis showed that the parameters contributing most to the uncertainty in the results were the initial adherence to statin prescription and the annual dropout rate from statins, each explaining 10-20% of the variability in the expected QALY gains for the base case relative to no health checks, and in the scenarios where extra statin treatment is given.

## why aren't more of them bigger?
## var(y) is mostly statistical error, MC error 5%
## Must be HSE, ELSA sampling uncertainty - hard to quantify contributions of these, we're at the boundaries of methodology
## Some other results may look funny, but the ests are all within a SE of zero, so can't interpret these.   Only statins compliance is valuable to learn







############################################################





############################################################

## BASELINE CHARACTERISTICS
## First col: all pop (40-45).
## Second, people who go on to attend at least one HC: chars measured at baseline, not first HC time
## Produced from simulation results at first time.
## Done with one run of the model, as big as will fit in memory (250k).

pops <- c("All", "Eligible", "Attending", "Any treatment", "Most deprived", "Least deprived")
npops <- length(pops)
popsall <- c(pops, 
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)",
          "Attending and eligible for any treatment (Q20)", "Attending and eligible for any treatment (Q10)",
          "Offered statins via HC", "Offered AHT via HC", "Offered SC via HC", "Offered WR via HC", 
          "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
npopsall <- length(popsall)

N <- read.csv("results/paperfeb/base_N_0.csv", header=FALSE, row.names=popsall)
B <- read.csv("results/paperfeb/base_baseline_0.csv", header=FALSE, col.names=c("all","elig","hc"))
vnames <-
    rbind(cbind("Gender",c("Male","Female")),
          cbind("Ethnicity", c("0","White","Indian","Pakistani","Bangladeshi","OtherAsian","Caribbean","African","Chinese","Other")),
          cbind("Deprivation",c("0","1_Lowest","2","3","4","5_Highest")),
          cbind("Education",1:4),
          cbind("Age",c("mean","q25","q75")),
          cbind("QRisk",c("mean","q25","q75")),
          cbind("QRisk20","n"),
          cbind("Qrisk1020","n"),
          cbind("SBP",c("mean","q25","q75")),
          cbind("DBP",c("mean","q25","q75")),
          cbind("BP>140/90","n"),
          cbind("TreatedHypertension","n"),
          cbind("Cholesterol",c("mean","q25","q75")),
          cbind("TC/HDL",c("mean","q25","q75")),
          cbind("BMI",c("mean","q25","q75")),
          cbind("BMI30","n"),
          cbind("Hba1c",c("mean","q25","q75")),
          cbind("HBA65","n"),
          cbind("Diabetes","n"),
          cbind("Smoking",c("Never","Ex","Current<10","Current10-19","Current20+"))
          )
ndec <- 1
colnames(vnames) <- c("var","val")
B <- cbind(B, vnames)

B <- B %>%
  filter(!(val==0)) %>%
  filter(!(var=="Education" & val==2)) %>% 
    mutate(all=round(all,ndec),
           elig=round(elig,ndec),
           hc=round(hc,ndec)) %>% 
    mutate(merge=1:n()) %>%
    mutate(merge = ifelse(val %in% c("Bangladeshi","OtherAsian","Chinese","Other"), "OtherEth", merge)) %>%
    mutate(merge = factor(merge, levels=unique(merge))) %>% 
    group_by(merge) %>%
    summarise(all=sum(all), elig=sum(elig), hc=sum(hc), var=last(var), val=last(val)) %>%
    mutate(merge=1:n()) %>%
    mutate(merge = ifelse(val %in% c("Current<10","Current10-19","Current20+"), "Current", merge)) %>%
    mutate(merge = factor(merge, levels=unique(merge))) %>% 
    group_by(merge) %>%
    summarise(all=sum(all), elig=sum(elig), hc=sum(hc), var=last(var), val=last(val)) %>%
    mutate(order = 1:n()) %>%
    as.data.frame()
inds <- B$order[match(c("Caribbean", "African", "Other"), B$val)]  # put Other ethnicity at the end
B$order[inds] <- sort(B$order[inds])
levels(B$val)[levels(B$val)=="Current20+"] <- "Current"

BN <- data.frame(val="N", all=Npop["Base case","All"], elig=Npop["Base case","Eligible"], hc=Npop["Base case","Attending"])
BN[,-1] <- 1000000 * BN[,-1] / BN[,"all"]
Bmean <- filter(B, val=="mean") %>% 
    mutate(meanall=all, meanel=elig, meanhc=hc) %>%
    select(var, meanall, meanel, meanhc, order)
Bq25 <- filter(B, val=="q25") %>%
    mutate(q25all=all, q25el=elig, q25hc=hc) %>%
    select(var, q25all, q25el, q25hc)
Bq75 <- filter(B, val=="q75") %>%
    mutate(q75all=all, q75el=elig, q75hc=hc) %>%
    select(var, q75all, q75el, q75hc)
Bcont <- left_join(left_join(Bmean, Bq25), Bq75) %>%
  mutate(all=paste0(meanall, " (", q25all, ", ", q75all, ")")) %>%
  mutate(elig=paste0(meanel, " (", q25el, ", ", q75el, ")")) %>%
  mutate(hc=paste0(meanhc, " (", q25all, ", ", q75hc, ")")) %>%
  mutate(val=var) %>% 
  select(var, val, all, elig, hc, order) 
Bn <- filter(B, val=="n") %>%
    mutate(val=var,
           all=paste0(round(100*all/N["All",],ndec), "%"),
           elig=paste0(round(100*elig/N["Eligible",],ndec), "%"),
           hc=paste0(round(100*hc/N["Attending",],ndec), "%")) %>% 
    select(var, val, all, elig, hc, order)
catvars <- c("Gender", "Ethnicity", "Deprivation", "Education", "Smoking")
Bcat <- filter(B, var %in% catvars) %>%
    mutate(all=paste0(round(100*all/N["All",],ndec),"%"),
           elig=paste0(round(100*elig/N["Eligible",],ndec),"%"),
           hc=paste0(round(100*hc/N["Attending",],ndec), "%")) %>% 
    select(var, val, all, elig, hc, order)
Bhead <- Bcat %>% group_by(var) %>% summarise(order=min(order) - 0.5) %>%
    mutate(val=var, all="", elig="", hc="") %>%
    select(var, val, all, elig, hc, order)
Bres <- rbind(Bn, Bcont, Bcat, Bhead) %>% 
    arrange(order) %>%
    select(val, all, elig, hc)
Bres <- rbind(BN, Bres)
Bres

write.table(Bres, "tab_base.txt", quote=FALSE, sep="\t", row.names=FALSE, na="")

### what is standard error around percentages?
### n=250k 

r <- filter(B, var %in% c("n", catvars))$all
n <- 250000
p <- r/n
pse <- sqrt(p*(1 - p)/ n)
round(100*cbind(p, pse, p-1.96*pse, p+1.96*pse), 2)
Nall

# These are only really accurate to 2 sf 


