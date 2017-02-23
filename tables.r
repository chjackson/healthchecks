### TODO
### Statistical uncertainty for at least base case, ideally all scenarios 

library(tidyverse)
library(forcats)

############################################################
###  TABLES OF SCENARIOS 
############################################################

nruns <- 40
npoprun <- 25000
aname <- "results/paperfeb/scen"
n <- npoprun*nruns

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

M <- S <- N <- Npop <- Tot <- Sumsq <- vector(nruns, mode="list")
P <- PM <- PS <- Ptot <- Psumsq <- vector(nruns, mode="list")
Marr <- array(dim=c(nruns, nsc, nouts, npops)) # different format, for statistical uncertainty
Percarr <- array(dim=c(nruns, nsc, npops)) # different format, for statistical uncertainty

for (i in 1:nruns){
    M[[i]] <- array(as.matrix(read.csv(sprintf("%s_mean_%s.csv",aname,i),header=FALSE)),
                    dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
    Marr[i,,,] <- M[[i]]
    S[[i]] <- array(as.matrix(read.csv(sprintf("%s_SD_%s.csv",aname,i),header=FALSE)),
                    dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
    N[[i]] <- as.matrix(read.csv(sprintf("%s_N_%s.csv",aname,i),header=FALSE))
    dimnames(N[[i]]) <- list(scn, npopn)
    P[[i]] <- array(as.matrix(read.csv(sprintf("%s_post_%s.csv",aname,i),header=FALSE)),
                    dim=c(nsc, 2, npostn), dimnames=list(scn, c("Mean","SD"), postn))
    PM[[i]] <- P[[i]][,1,]; PS[[i]] <- P[[i]][,2,]
    if (nsc==1){
        Percarr[i,1,] <- N[[i]][,2:9,drop=FALSE] / N[[i]][,1]
    }
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

M[,c(outn.d, outn.ev),] <- M[,c(outn.d, outn.ev),]*n
S[,c(outn.d, outn.ev),] <- S[,c(outn.d, outn.ev),]*n
SE <- S
for (j in 1:nouts){
    SE[,j,] <- S[,j,] / sqrt(Npop)
}
PSE <- PS / sqrt(Npop[,"Attending"])

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
res[,c(1,6:9)]
res[,10:12]
res[,13:18]

## text file to paste into Word and convert to table. 
write.table(res, "tab_res.txt", quote=FALSE, sep="\t", row.names=FALSE, na="")

### Chang results 
P["Base case",]
PSE["Base case",]



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
Nall <- Npop["Base case","All"]
Nel <- Npop["Base case","Eligible"]
Natt <- Npop["Base case","Attending"]
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
B <- cbind(B, vnames) %>%
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

BN <- data.frame(val="N", all=Nall, elig=Nel, hc=Natt)
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
           all=paste0(round(100*all/Nall,ndec), "%"),
           elig=paste0(round(100*elig/Nel,ndec), "%"),
           hc=paste0(round(100*hc/Natt,ndec), "%")) %>% 
    select(var, val, all, elig, hc, order)
catvars <- c("Gender", "Ethnicity", "Deprivation", "Education", "Smoking")
Bcat <- filter(B, var %in% catvars) %>%
    mutate(all=paste0(round(100*all/Nall,ndec),"%"),
           elig=paste0(round(100*elig/Nel,ndec),"%"),
           hc=paste0(round(100*hc/Natt,ndec), "%")) %>% 
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
