
############################################################

## TODO
## Do baseline table from a 250k pop base case run .   1mil ideal but faffy 
## Full unc on baseline 

## IQRs for baseline, but can argue not needed 

## Bullet points for both






############################################################

### TABLES FOR PAPER, DESIGNED BY OLI

## TABLE 2: BASELINE CHARACTERISTICS
## First col: all pop (40-45).    Second, people who go on to attend at least one HC: chars measured at baseline, not first HC time
## Produce from simulation results at first time.
## Check matches values from HSE 

## Does the sim record all this?
## self.eth, SES, etc, educ, age, qrisk, SBP, DBP, blah
## do once. 

## DO THIS WITH ONE RUN OF THE MODEL, AS BIG AS WILL FIT IN MEMORY.
## 250k then multiply by 4.
library(tidyverse)
library(forcats)


mn <- c("All", "Eligible", "Attending", "Any treatment", "Deprived", "Eligible for treatment (Q20)", "Eligible for treatment (Q10)", "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
N <- read.csv("results/paper/base_N_0.csv", header=FALSE, row.names=mn)
Nall <- N["All",]
Natt <- N["Attending",]
B <- read.csv("results/baseline.csv", header=FALSE, col.names=c("all","hc"))
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
colnames(vnames) <- c("var","val")
B <- cbind(B, vnames) %>%
  filter(!(val==0)) %>%
  filter(!(var=="Education" & val==2)) %>% 
    mutate(all=round(all, 2),
           hc=round(hc,2)) %>% 
    mutate(merge=1:n()) %>%
    mutate(merge = ifelse(val %in% c("Bangladeshi","OtherAsian","Chinese","Other"), "OtherEth", merge)) %>%
    mutate(merge = factor(merge, levels=unique(merge))) %>% 
    group_by(merge) %>%
    summarise(all=sum(all), hc=sum(hc), var=last(var), val=last(val)) %>%
    mutate(merge=1:n()) %>%
    mutate(merge = ifelse(val %in% c("Current<10","Current10-19","Current20+"), "Current", merge)) %>%
    mutate(merge = factor(merge, levels=unique(merge))) %>% 
    group_by(merge) %>%
    summarise(all=sum(all), hc=sum(hc), var=last(var), val=last(val)) %>%
    mutate(order = 1:n()) %>%
    as.data.frame()
inds <- B$order[match(c("Caribbean", "African", "Other"), B$val)]  # put Other ethnicity at the end
B$order[inds] <- sort(B$order[inds])
levels(B$val)[levels(B$val)=="Current20+"] <- "Current"

BN <- data.frame(val="N", all=Nall, hc=Natt)
Bmean <- filter(B, val=="mean") %>% 
    mutate(meanall=all, meanhc=hc) %>%
    select(var, meanall, meanhc, order)
Bq25 <- filter(B, val=="q25") %>%
    mutate(q25all=all, q25hc=hc) %>%
    select(var, q25all, q25hc)
Bq75 <- filter(B, val=="q75") %>%
    mutate(q75all=all, q75hc=hc) %>%
    select(var, q75all, q75hc)
Bcont <- left_join(left_join(Bmean, Bq25), Bq75) %>%
  mutate(all=paste0(meanall, " (", q25all, ", ", q75all, ")")) %>%
  mutate(hc=paste0(meanhc, " (", q25all, ", ", q75hc, ")")) %>%
  mutate(val=var) %>% 
  select(var, val, all, hc, order) 
Bn <- filter(B, val=="n") %>%
    mutate(val=var,
           all=paste0(round(100*all/Nall,2), "%"),
           hc=paste0(round(100*hc/Natt,2), "%")) %>% 
    select(var, val, all, hc, order)
catvars <- c("Gender", "Ethnicity", "Deprivation", "Education", "Smoking")
Bcat <- filter(B, var %in% catvars) %>%
    mutate(all=paste0(round(100*all/Nall,2),"%"),
           hc=paste0(round(100*hc/Natt,2), "%")) %>% 
    select(var, val, all, hc, order)
Bhead <- Bcat %>% group_by(var) %>% summarise(order=min(order) - 0.5) %>%
    mutate(val=var, all="", hc="") %>%
    select(var, val, all, hc, order)
Bres <- rbind(Bn, Bcont, Bcat, Bhead) %>% 
    arrange(order) %>%
    select(val, all, hc)
Bres <- rbind(BN, Bres)
Bres

write.table(Bres, "tab_base.txt", quote=FALSE, sep="\t", row.names=FALSE, na="")

## TODO

## IQR 
# Compute summaries using 250000, but state pop size used for results is 1000000 - the summaries will be the same within presented precision.

## dplyr learning todo
## how to use mutate to change a subset of elements without using ifelse
Bres %>% within(val[val=="Other"] <- "foo")



############################################################

## Each of next four tables has same rows in four row sets: 

## Eligibility and uptake
#Number eligible for HC (at any point)
#Number have one or more HC
#Mean HC per head of population

## Treatment
## Cases prevented
## Deaths prevented
## Change in survival




############################################################

### Table A1 (probably goes in the appendix or may be come T2): Eligiblity and uptake for HC under ‘business as usual’


### Table 3a: Change in eligibility and outcomes where eligibility criteria change relative to ‘standard’ programme (n=1,000,000)


### Table 3b (4): Change in uptake and outcomes where uptake of the healthcheck programme changes


### Table 3c (5): Change in uptake and outcomes with increased uptake of treatments offered by the health check programme changes


## do all four at once
## we dont need short term outcomes 

nruns <- 5
aname <- "results/paper/scen"
n <- 200000*nruns

## for M 
popn <- c("All", "Eligible", "Attending", "Any treatment", "Most deprived")
npops <- length(popn)

scn <- c("Base case","Invite high BP", "Attend age 50-74", "Attend age 40-80", "Attend age 50-80", "Baseline uptake 0.85 (from 0.48)","Uptake 10 - 20% more in most deprived","Uptake of smokers +30%","+20% uptake for QRisk > 15", "Target non-attenders","Statin prescription x 2.5", "AHT prescription x 2.5", "Smoking referral x 2.5", "Weight referral x 2.5", "All treatments x 2.5")
nsc <- length(scn)

evn <- c("IHD","STR","DEM","LC")
outn.ev <- paste(rep(evn, each=3), rep(c(80,105), each=length(evn)*3), paste(rep(c("_con","_hc","_i"), length(evn)),sep=""),sep="")
dn <- c("D")
outn.d <- paste(rep(dn, each=3), rep(c(75,80), each=length(dn)*3), paste(rep(c("_con","_hc","_i"), length(dn)),sep=""),sep="")
outn <- c("QALY_con","QALY_hc","IQALY","LY_con","LY_hc","ILY","IQALYperhc", "ILYperhc", outn.ev,outn.d)
nouts <- length(outn)

## for N 
npopn <- c(popn, 
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)",
          "Statins via HC", "Antihypertensives via HC", "Smoking cessation via HC", "Weight management via HC","Number of HCs")
npopsn <- length(npopn)


M <- S <- N <- Npop <- Tot <- Sumsq <- vector(nruns, mode="list")
Marr <- array(dim=c(nruns, nsc, nouts, npops)) # different format, for statistical uncertainty
Parr <- array(dim=c(nruns, nsc, npops)) # different format, for statistical uncertainty

for (i in 1:nruns){
    M[[i]] <- array(as.matrix(read.csv(sprintf("%s_mean_%s.csv",aname,i),header=FALSE)),
                    dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
    Marr[i,,,] <- M[[i]]
    S[[i]] <- array(as.matrix(read.csv(sprintf("%s_SD_%s.csv",aname,i),header=FALSE)),
                    dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
    N[[i]] <- as.matrix(read.csv(sprintf("%s_N_%s.csv",aname,i),header=FALSE))
    if (nsc==1){
        Parr[i,1,] <- N[[i]][,2:9,drop=FALSE] / N[[i]][,1]
    }
    dimnames(N[[i]]) <- list(scn, npopn)
    Npop[[i]] <- N[[i]][,1:npops,drop=FALSE]
    Tot[[i]] <- Sumsq[[i]] <- array(dim=dim(M[[i]]))
    for (j in 1:nouts){
        Tot[[i]][,j,] <- M[[i]][,j,] * Npop[[i]]
        Sumsq[[i]][,j,] <- S[[i]][,j,] * Npop[[i]]
    }
}

Nall <- Reduce("+", N)
Npopall <- Reduce("+", Npop)
Tot <- Reduce("+", Tot)
Sumsq <- Reduce("+", Sumsq)

Mall <- Sall <- array(dim=dim(Tot))
for (j in 1:nouts){ 
    Mall[,j,] <- Tot[,j,] / Npopall
    Sall[,j,] <- Sumsq[,j,] / Npopall
}

M <- Mall; S <- Sall; Npop <- Npopall
dimnames(M) <- dimnames(S) <- list(scn, outn, popn)
dimnames(Npop) <- list(scn, popn)
SE <- S
for (j in 1:nouts){
    SE[,j,] <- S[,j,] / sqrt(Npop)
}

M[,c(outn.d, outn.ev),] <- M[,c(outn.d, outn.ev),]*n
S[,c(outn.d, outn.ev),] <- S[,c(outn.d, outn.ev),]*n

rows.ev <- c("IHD80_i",  "STR80_i", "DEM80_i",  "LC80_i", 
          "IHD105_i", "STR105_i", "DEM105_i", "LC105_i",
          "D75_i", "D80_i")   
rows.ly <- c("IQALY", "ILY","IQALYperhc", "ILYperhc")
resm <- t(M[,c(rows.ev,rows.ly),"All"])
res_el <- t(M[,c("IQALY","ILY"),"Eligible"]); rownames(res_el) <- c("IQALY_el","ILY_el")
res_att <- t(M[,c("IQALY","ILY"),"Attending"]); rownames(res_att) <- c("IQALY_att","ILY_att")
res_dep <- t(M[,c("IQALY","ILY"),"Most deprived"]); rownames(res_dep) <- c("IQALY_dep","ILY_dep")
resm <- rbind(resm, res_el, res_att, res_dep)
resm <- round(resm, 2)
resm <- array(as.character(resm), dim=dim(resm), dimnames=dimnames(resm))
resn <- t(Nall)
resperc <- 100*resn/resn[1,1]
resperc["Number of HCs",] <- resperc["Number of HCs",]/100
resperc <- round(resperc,1)
resperc <- array(as.character(resperc), dim=dim(resperc), dimnames=dimnames(resperc))

## Add SE to each of above?  Looks small, so can present assuming negligible.
## QALY gain from narrowing eligibility is still implausible though. 

## Paste standard errors on the end
## SE for difference in two binomial event counts
rows.evb <- c("IHD80",  "STR80", "DEM80",  "LC80", "IHD105", "STR105", "DEM105", "LC105", "D75", "D80")
rows.evc <- paste0(rows.evb, "_con")
rows.evh <- paste0(rows.evb, "_hc")
pc <- t(M[,c(rows.evc),"All"])/n
ph <- t(M[,c(rows.evh),"All"])/n
sediff <- round(sqrt(n*(pc*(1-pc) + ph*(1-ph))),0) # se of diff in y. 
## SE for continuous LY and QALY gains 
## NO THIS IS WRONG, THEY'RE NOT INDEPENDENT

sem <- t(SE[,rows.ly,"All"])
se_el <- t(SE[,c("IQALY","ILY"),"Eligible"]); rownames(se_el) <- c("IQALY_el","ILY_el")
se_att <- t(SE[,c("IQALY","ILY"),"Attending"]); rownames(se_att) <- c("IQALY_att","ILY_att")
se_dep <- t(SE[,c("IQALY","ILY"),"Most deprived"]); rownames(se_dep) <- c("IQALY_dep","ILY_dep")
sem <- formatC(rbind(sem, se_el, se_att, se_dep), digits=2)
rows.ly <- rownames(sem)
resm.ev <- resm[rows.ev,]
resm.ly <- array(paste0(resm[rows.ly,], " (", rbind(sem), ")"), dim=dim(resm[rows.ly,]),  dimnames=dimnames(resm[rows.ly,]))

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
          "Eligible for any treatment (Q20)", "Eligible for any treatment (Q10)", "Any treatment",
          "Statins via HC","Antihypertensives via HC","Weight management via HC","Smoking cessation via HC",
          "Cases prevented by age 80",
          "IHD80_i",  "STR80_i", "DEM80_i",  "LC80_i", 
          "Cases prevented",
          "IHD105_i", "STR105_i", "DEM105_i", "LC105_i",
          "Deaths prevented",
          "D75_i", "D80_i",
          "Change in QALY",
          "IQALY","IQALY_el","IQALY_att","IQALYperhc","IQALY_dep",
          "Change in lifetime",
          "ILY","ILY_el","ILY_att","ILYperhc","ILY_dep")

res <- rbind(resperc, resm.ev, resm.ly, reshead)[rows,]
options(width=180)
res[,1:5]
res[,6:10]
res[,11:15]

write.table(res, "tab_res.txt", quote=FALSE, sep="\t", row.names=FALSE, na="")

SE[,c("IQALY","ILY"),"All"]



M[1,3,]

## CHECK OF BASE CASE 
## NUMBERS OK.  Note eligible higher then attending, lower than HC eligible as conditions on HC el

## CASES
## 1500 CVD cases out of a million by age 80
## 105 because cases delayed rather than prevented.  Will die of something eventually.  Exactly zero surprising though, worth investigating
## Dementia - more cases as CVD delayed 
## Deaths, 200, 300 by 75, 80 instead of 72 by 70 out of half pop size. OK 

## High number of dementia cases prevented under all scenarios is fishy.
## doesn't condition on currently alive FIXME 


## QALYs
## 3 days in base case.  was 4.3 for open.
## ?? 6 days for age 35-40 in excel results.  Fishy , todo run both again to compare 
## 5.3 days for a single run of 250000.   Suspect indexing fucked.
## Much narrower increase for eligible, attend.  Fishy. 
## OK M[0,0,2] contains correct QALY gain on exit, and this is propagated to res, with one run. 
## Is it the aggregation over runs that's going wrong? 
## No, looks ok from 2 run 2500 test (papertest) 
## has something else changed since the december results?  something when I looked at events 
## randseed thing in the baseline cvd prevalence? -- shall we rerun without this? 
## YES - can reproduce excel results when revert to old seed behaviour.
## So stored QALYs seem correct 

## I might have some explaining to do about QALY gains among those attending at least once.  1.5 times as high now, when they were 2.5 times as high before.   "Any treatment" ok though - 3-4 times "attenders" figure. 

## New results on deprivation: gains 50% higher for most deprived 

## TODO CHECKS ON SCENARIOS

## BP ?? IQALY is this really 3 on top of base case?  looks like compared to no hc... look at others first.

## Invite high BP is good now.  TODO check open cohort results with new seed thing 



## No Attend 50-74 looks wrong.   If ILY is negative, then IQALY should be too. 
## or is this plausible?  Never actually studied LYs before.
## Under new policy of not eligible for HC till age 50, have reduced LY, but increased QALY
## screen later, don't pick up disease until later

## how can that happen?  less time spent alive, but also less time spent with disease 

## which ones have negative ILY: OK only the attendance ones.  
## MC error causes negative values for small pops. 
## Disease prevention leads to bigger QALY gains where small LY gains, OK makes sense
## but switching sign still odd.

## SCREEN FROM 40      :  SCREEN FROM 50 
## treat earlier (45)  : treat later (55) 
## disease later (70)  :  disease earlier (60)
## die later     (80)  : die earlier (70)

## 40                  : 30 
## QALY since age 40
## 30 + 0.6*10         :  20 + 0.6+10

# QALY = LY*(1-p + 0.6*p) =  LY*(1 - 0.4*p)
# where p is prop of time spent with disease, in (0,1) 

# QALY2 - QALY1  =     LY2*(1 - 0.4*p2) -  LY1*(1 - 0.4*p1)
#                      LY2*[0.6 to 1]   - LY1*[0.6 to 1]
## QALY2 / QALY1 = (LY2/LY1) * (1 - 0.4p2)/(1 - 0.4p1) 
## Suppose LY2 / LY1 < 1, so life lost in scenario 2.
## Can QALY2/QALY1 be > 1,  so QALYs gained in scenario 2. 
## Yes it can if (1 - 0.4p2)/(1 - 0.4p1) >> 1.
## Can this happen?
## (1 - 0.4p2) >> (1 - 0.4p1)
## p2 << p1 
## Yes if p1 >> p2
## so much lower prop of lifetime spent with disease in scenario 2. 
## die more quickly, shorter LY, but less time with disease 
## huh? 



## it's random seeding surely.
## Seeds should ensure that person shouldn't live longer if screened later.
## difference should only arise after treatment
## why is np.random.seed not seeded with self.randseed?
## sometimes it is, sometimes it's called with 0 or fixed i value
## fixed i at start of each timestep.
## but should this be i + self.randseed?   arguably. shouldn't affect increments though

## only called with default (=None=seed from OS or clock, ie not reproducible), when instructed from program args

## ADDED TO MAIN
## two self.randseeds and a rnd.seed (0) in InitialisePopulation
## added seeds before treatment
## seeds before death already there, so people should have same mortality in different scenarios if they have same disease.  or always die sooner if older 

## todo search for all random gen functions and ensure seed set
# np.random   and rnd. 
# 

## Now all qaly gains about 5 times what they were with fewer seed settings
## What if do more seed settings will they change even more?!

## are absolute QALYs similar

i <- 1
aname <- "results/papertest/scen"
Mt <- array(as.matrix(read.csv(sprintf("%s_mean_%s.csv",aname,i),header=FALSE)),
            dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
St <- array(as.matrix(read.csv(sprintf("%s_SD_%s.csv",aname,i),header=FALSE)),
                    dim=c(nsc, nouts, npops), dimnames=list(scn, outn, popn))
N <- as.matrix(read.csv(sprintf("%s_N_%s.csv",aname,i),header=FALSE))

Mt[1:3,c("QALY_con","QALY_hc","IQALY"),1]

## Aha so if run once 

## Eh?  It's 4 from scenario 1.   So what's with big res 

St[2,c("QALY_con","QALY_hc"),1] / sqrt(N[1,1])

## no, similar control outcome, but better HC outcome.  similar SEs. 
## better outcomes for all HC scenarios 
