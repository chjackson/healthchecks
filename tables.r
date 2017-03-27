library(tidyverse)
library(forcats)

############################################################
###  TABLES OF SCENARIOS 
############################################################

nruns <- 40
npoprun <- 25000
aname <- "results/papermar/scen"
n <- npoprun*nruns
statunc <- TRUE

scn <- c("Base case","Invite high BP", "Attend age 50-74", "Attend age 40-80", "Attend age 50-80",
         "Baseline uptake +30%","Uptake +30% in most deprived","Uptake of smokers +30%","+30% uptake for QRisk > 20",
         "Target non-attenders", "Attenders keep attending", "Attenders keep attending, target non-attenders", 
         "Statin prescription x 2.5", "AHT prescription x 2.5", "Smoking referral x 2.5", "Weight referral x 2.5",
         "All treatments x 2.5", "Higher treatment, higher uptake and invite BP")
nsc <- length(scn)

evn <- c("IHD","STR","DEM","LC")
outn.ev <- paste(rep(evn, each=3), rep(c(80,100), each=length(evn)*3), paste(rep(c("_con","_hc","_i"), length(evn)),sep=""),sep="")
dn <- c("D")
outn.d <- paste(rep(dn, each=3), rep(c(75,80), each=length(dn)*3), paste(rep(c("_con","_hc","_i"), length(dn)),sep=""),sep="")
outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY","IQALYperhc", "ILYperhc", outn.ev,outn.d)
nouts <- length(outn)

# for N
popn <- c("All", "Eligible", "Attending", "Any treatment", "Most deprived", "Least deprived")
npops <- length(popn)
npopn <- c(popn, 
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)",
          "Attending and eligible for any treatment (Q20)", "Attending and eligible for any treatment (Q10)",
          "Offered statins via HC", "Offered AHT via HC", "Offered SC via HC", "Offered WR via HC", 
          "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
npopsn <- length(npopn)
postn <- c("QRisk","SBP","DBP")
npostn <- length(postn)

## Results without statistical uncertainty

if (!statunc) { 
    M <- S <- N <- Npop <- Tot <- Sumsq <- vector(nruns, mode="list")
    P <- PM <- PS <- Ptot <- Psumsq <- vector(nruns, mode="list")

    for (i in 1:nruns){
        M[[i]] <- array(as.matrix(read.csv(sprintf("%s_mean_%s.csv",aname,i),header=FALSE)),
                        dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
        S[[i]] <- array(as.matrix(read.csv(sprintf("%s_SD_%s.csv",aname,i),header=FALSE)),
                        dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
        N[[i]] <- as.matrix(read.csv(sprintf("%s_N_%s.csv",aname,i),header=FALSE))
        dimnames(N[[i]]) <- list(scn, npopn)
        P[[i]] <- array(as.matrix(read.csv(sprintf("%s_post_%s.csv",aname,i),header=FALSE)),
                        dim=c(nsc, 2, npostn), dimnames=list(scn, c("Mean","SD"), postn))
        PM[[i]] <- P[[i]][,1,]; PS[[i]] <- P[[i]][,2,]
        Npop[[i]] <- N[[i]][,1:npops,drop=FALSE]
        Tot[[i]] <- Sumsq[[i]] <- array(dim=dim(M[[i]]))
        for (j in 1:nouts){
            Tot[[i]][,j,] <- M[[i]][,j,] * Npop[[i]]
            Sumsq[[i]][,j,] <- S[[i]][,j,] * Npop[[i]]
        }
        Ptot[[i]] <- P[[i]][,1,] * Npop[[i]][,"Attending"]
        Psumsq[[i]] <- P[[i]][,2,] * Npop[[i]][,"Attending"]
    }

    Nall <- Reduce("+", N)
    Npopall <- Reduce("+", Npop)
    Tot <- Reduce("+", Tot)
    Sumsq <- Reduce("+", Sumsq)
    Ptot <- Reduce("+", Ptot)
    Psumsq <- Reduce("+", Psumsq)

    Mall <- Sall <- array(dim=dim(Tot))
    for (j in 1:nouts){ 
        Mall[,j,] <- Tot[,j,] / Npopall
        Sall[,j,] <- Sumsq[,j,] / Npopall
    }
    Pall <- Ptot / Npopall[,"Attending"]
    PSall <- Psumsq / Npopall[,"Attending"]
    
    M <- Mall; S <- Sall; Npop <- Npopall; P <- Pall; PS <- PSall
    dimnames(M) <- dimnames(S) <- list(scn, outn, popn)
    dimnames(Npop) <- list(scn, popn)
    SE <- S
    for (j in 1:nouts){
        SE[,j,] <- S[,j,] / sqrt(Npop)
    }
    PSE <- PS / sqrt(Npop[,"Attending"])
}

### Results with statistical uncertainty

if (statunc) { 
    npsa <- 1000 # number of runs with different parameter values. 100 took 7051 sec on HPC

    resm <- as.matrix(read.table("results/papermar/unc_mean.csv", colClasses="numeric", sep=",", header=FALSE))
    ress <- as.matrix(read.table("results/papermar/unc_SD.csv", colClasses="numeric", sep=",", header=FALSE))
    str(resm) # npsa*nsc  x  npops*nouts + 1.
    runid <- resm[,1]; resm <- resm[,-1]; ress <- ress[,-1]
    n <- as.matrix(read.table("results/papermar/unc_N.csv", colClasses="numeric", sep=",", header=FALSE))[,(1:npops)+1]
    str(n) # npsa*nsc  x  npops
    Mrep <- array(resm, dim=c(nsc, npsa, nouts, npops))
    Srep <- array(ress, dim=c(nsc, npsa, nouts, npops))
    dimnames(Mrep) <- dimnames(Srep) <- list(scn[1:nsc], 1:npsa, outn, popn)
    Nrep <- array(n, dim=c(nsc, npsa, npops))
    dimnames(Nrep) <- list(scn[1:nsc], 1:npsa, popn)

### for count outcomes, uses normal approx to dist of proportion
    n <- apply(Nrep, c(1,3), mean) # mean within-simulation pop size
    n <- array(n, dim=c(nsc,1,npops))[,rep(1,nouts),,drop=FALSE] # match dims of M, S
    mu <- apply(Mrep, c(1,3,4), mean)
    sigsq <- apply(Mrep, c(1,3,4), var)
    tausq <- apply(Srep^2, c(1,3,4), mean)
    ## Monte Carlo standard error around the mean
    semu <- sqrt(sigsq / npsa + tausq/(npsa*n))
    ## Statistical uncertainty
    sd <- sqrt(sigsq) # standard MC est is biased
    ## Monte Carlo standard error around the variance (SD^2)
    sevar <- sqrt(2/(npsa-1))*(sigsq + tausq/n)
    ## Monte Carlo standard error around the SD (using delta method
    sesd <- 0.5*sevar/sd
    ## Bias-corrected estimates, may not work for smaller n
    sigsqA <- sigsq - tausq/n # ...as this could be negative
    sdA <- sqrt(sigsqA)
    sevarA <- sqrt(2*((sigsqA + tausq/n)^2 / (npsa-1) + tausq^2 /(npsa*n^2*(n-1))))
    sesdA <- 0.5*sevarA/sdA
    res <- array(c(mu, sd, semu), dim=c(dim(mu), 3),
                 dimnames=c(dimnames(mu), list(c("mu","sd","semu"))))
    M <- mu
    SE <- sd 
}

## convert events averted from proportion to "per million population"
M[,c(outn.d, outn.ev),] <- M[,c(outn.d, outn.ev),]*1000000
SE[,c(outn.d, outn.ev),] <- SE[,c(outn.d, outn.ev),]*1000000
rows.ev <- c("IHD80_i",  "STR80_i", "DEM80_i",  "LC80_i", 
          "IHD100_i", "STR100_i", "DEM100_i", "LC100_i",
          "D75_i", "D80_i")   
rows.ly <- c("IQALY", "ILY","IQALYperhc", "ILYperhc")
resm <- t(M[,c(rows.ev,rows.ly),"All"])
res_el <- t(M[,c("IQALY","ILY"),"Eligible"]); rownames(res_el) <- c("IQALY_el","ILY_el")
res_att <- t(M[,c("IQALY","ILY"),"Attending"]); rownames(res_att) <- c("IQALY_att","ILY_att")
res_dep <- t(M[,c("IQALY","ILY"),"Most deprived"]); rownames(res_dep) <- c("IQALY_dep","ILY_dep")
res_ndep <- t(M[,c("IQALY","ILY"),"Least deprived"]); rownames(res_ndep) <- c("IQALY_ndep","ILY_ndep")
resm <- rbind(resm, res_el, res_att, res_dep, res_ndep)
resm <- round(resm, 2)
resm <- array(as.character(resm), dim=dim(resm), dimnames=dimnames(resm))
resn <- t(Nall)
resperc <- 100*resn/resn[1,1]
resperc["Number of HCs",] <- resperc["Number of HCs",]/100
resperc <- round(resperc,3)
resperc <- array(as.character(resperc), dim=dim(resperc), dimnames=dimnames(resperc))

## Paste standard errors for life / event count gains 
sem <- t(SE[,rows.ly,"All"])
se_el <- t(SE[,c("IQALY","ILY"),"Eligible"]); rownames(se_el) <- c("IQALY_el","ILY_el")
se_att <- t(SE[,c("IQALY","ILY"),"Attending"]); rownames(se_att) <- c("IQALY_att","ILY_att")
se_dep <- t(SE[,c("IQALY","ILY"),"Most deprived"]); rownames(se_dep) <- c("IQALY_dep","ILY_dep")
se_ndep <- t(SE[,c("IQALY","ILY"),"Least deprived"]); rownames(se_ndep) <- c("IQALY_ndep","ILY_ndep")
sem <- formatC(rbind(sem, se_el, se_att, se_dep, se_ndep), digits=2)
rows.ly <- rownames(sem)
resm.ly <- array(paste0(resm[rows.ly,], " (", sem, ")"),
                 dim=dim(resm[rows.ly,]),  dimnames=dimnames(resm[rows.ly,]))
se_ev <- round(t(SE[,rows.ev,"All"]), 0)
resm.ev <- array(paste0(resm[rows.ev,], " (", se_ev, ")"),
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
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)",
          "Attending and eligible for any treatment (Q20)", "Attending and eligible for any treatment (Q10)",
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
res[,c(6:9)]
res[,10:12]
res[,13:18]

## text file to paste into Word and convert to table. 
write.table(res, "tab_res.txt", quote=FALSE, sep="\t", row.names=FALSE, na="")

### Chang results 
P["Base case",]
PSE["Base case",]

## Compare results with stat unc line by line to current results
## Point ests a bit different due to nonlinearity in model
## was 1092 now 1059
## 511 now 449 for stroke
## 5-10% increases in QALY gains
## stat SD is 10-20 times comp SE 
## notice large stat SDs for event counts 


## Expected value of partial perfect information for each parameter
library(earth)
source("parnames.r")

pars <- as.matrix(read.table("results/papermar/unctest_pars.csv", colClasses="numeric",sep=","))
colnames(pars) <- parnames
pars <- pars[,setdiff(colnames(pars), fixedpars)]
pars <- pars[rep(1:5, 1000/5), ] # replicate to 1000 for the moment and add noise, until proper run done 
pars <- pars + rnorm(1000, 0, 0.01)
npars <- ncol(pars)
## headline output 
y <- as.numeric(Mrep["Base case",,"IQALY","All"])
evpi <- sqrt(var(y))
pevppi <- matrix(nrow=npars, ncol=1)
for(i in 1:npars){
    x <- pars[,i]
    pevppi[i] <- var(fitted(earth(x, y))) / var(y)  ## EVPPI as prop of EVPI
}
pevppi[pevppi<1e-05] <- 0
eres <- data.frame(var=colnames(pars), pevppi=pevppi)
eres %>% arrange(desc(pevppi))
# 
### how accurate is this with 1000?
### SEs around regression estimates

## mod <- earth(X, Y, nfold=10, ncross=30, varmod.method="const", Get.leverages=TRUE)
## se <- sqrt(mod$varmod$model.var)
## B <- 1000
## evppi.rep <- numeric(B)
## for(i in 1:B)
##     evppi.rep[i] <- var(rnorm(length(fitted(mod)), fitted(mod), se)) 
## pe[i,j] <- sd(evppi.rep) / var(Y)

1


############################################################

## BASELINE CHARACTERISTICS
## First col: all pop (40-45).
## Second, people who go on to attend at least one HC: chars measured at baseline, not first HC time
## Produced from simulation results at first time.
## Done with one run of the model, as big as will fit in memory (250k).

popn <- c("All", "Eligible", "Attending", "Any treatment", "Most deprived", "Least deprived")
npops <- length(popn)
npopn <- c(popn, 
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)",
          "Attending and eligible for any treatment (Q20)", "Attending and eligible for any treatment (Q10)",
          "Offered statins via HC", "Offered AHT via HC", "Offered SC via HC", "Offered WR via HC", 
          "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
npopsn <- length(npopn)

N <- read.csv("results/paperfeb/base_N_0.csv", header=FALSE, row.names=npopn)
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


