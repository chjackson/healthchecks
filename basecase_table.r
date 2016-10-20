### TABLES ONLY NEEDED FOR BASE CASE, not other scenarios 

### PROCESS MEASURES 

source("scenarios_combine.r")

procf <- as.data.frame(proc) %>%
  subset(scen=="Base case") %>%  select(outcome, numf, percf)
procf$outcome <- pnms

write.table(procf, file="process.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.csv(procf, file="", quote=FALSE, na="", row.names=FALSE)






### SHORT TERM OUTCOMES 

sh1 <- short[short$scen=="Base case",]

cnames <- c("Q","CH","HDL","SBP","DBP","BMI")
diffs <- sh1[,paste0(cnames,"1")] - sh1[,cnames]
colnames(diffs) <- paste0(cnames,"diff")
bnames <- c("SM","cvd_events","lc_events","deaths","deaths_cvd","deaths_lc","deaths_dem","deaths_oth")
percs <- 100*(sh1[,paste0(bnames,"_t")] - sh1[,paste0(bnames,"_c")]) / sh1[,paste0(bnames,"_c")]
colnames(percs) <- paste0(bnames, "pdiff")
diffs <- round(diffs,2)
percs <- round(percs,1)
cf <- array(paste0(as.matrix(round(sh1[,cnames],2)), ",",
                  as.matrix(round(sh1[,paste0(cnames,"1")],2)),
                  " (", as.matrix(diffs), ")"), dim=dim(diffs))
cf <- as.data.frame(cf); colnames(cf) <- cnames 
bf <- array(paste0(as.matrix(round(sh1[,paste0(bnames,"_c")],2)), ",",
                  as.matrix(round(sh1[,paste0(bnames,"_t")],2)),
                  " (", as.matrix(percs), "%)"), dim=dim(percs))
bf <- as.data.frame(bf); colnames(bf) <- bnames 
sh1$trt <- factor(sh1$trt, levels=unique(sh1$trt))

res <-
  cbind(sh1[,c("trt","year")], cf, bf) %>% 
  gather(varname, outcome, Q:deaths_oth) %>% 
  spread(year, outcome) %>%
  mutate(varname = factor(varname, levels=c(cnames, bnames))) %>%
  arrange(trt, varname) %>%
  filter((trt=="statins" & varname %in% c("Q","CH","HDL","cvd_events","deaths_cvd","deaths"))|
         (trt=="aht" & varname %in% c("Q","SBP","DBP","cvd_events","deaths_cvd","deaths")) | 
         (trt=="wr" & varname %in% c("Q","BMI","SBP","cvd_events","deaths_cvd","deaths")) | 
         (trt=="sc" & varname %in% c("Q","SM","cvd_events","lc_events","deaths_cvd","deaths_lc","deaths")) 
         )
levels(res$varname) <- c("Mean QRisk","Mean cholesterol","Mean HDL","Mean SBP","Mean DBP","Mean BMI","Smokers", "CVD events","Lung cancer events","Deaths","Deaths (CVD)","Deaths (lung cancer)", "Deaths (dementia)", "Deaths (other)")

a <- factor(c("a","c","b"))
levels(a) <- c("A","B","C")
res$trt <- as.character(res$trt)
res$varname <- as.character(res$varname)
res2 <- select(res, -trt)
write.table(
    rbind(
        c("Statins",rep(NA,4)), res2[res$trt=="statins",],
        c("Antihypertensives",rep(NA,4)), res2[res$trt=="aht",],
        c("Weight management",rep(NA,4)), res2[res$trt=="wr",],
        c("Smoking cessation",rep(NA,4)), res2[res$trt=="sc",]),
    sep="\t", quote=FALSE, na="", row.names=FALSE
)

## mort reduction lowers death rate by 0.3%
## and overall qaly gains by 0.2 days



## number of smokers going down and then up

## methods appendix + short document results




## check SC prop number 
#0.001 of population referred to SC at any time 
#0.18 are smokers at first HC
#0.036 of smokers referred to SC 
#then should be at least 0.0065 referred to SC 
#  sc = (H1.SmokingCessation.sum(axis=1)>0)
# Aha, this is the indicator for those who quit after one year
# SmokingCessation_Offered is those offered SC 
# other pars are those offered and adherent, so we're being consistent. 

## inviting more smokers should be bigger effect
## RR from 0.75 to 1, so more like 30% increased rate of accepting offer for smokers 
## Increases P(smoker at first HC) from 18.8 to 20.5%
## 0.11 to 0.14 % of population quitting after SC, change of 0.003%.  30% increase
## which maps to only 0.1 days of QALY. 

## doubling baseline uptake rate should increase attendance rate by more
## currently 66% to 71% : 1.07 increase in ever attending
## Aha it's because under current assumptions, once people become eligible, they're highly likely to attend anyway 
## though this will increase attendance for those eligible for a short time 
## 91% eligible at least once
## up_HC_takeup goes from 0.48 to 0.85, "best practice" from Tower Hamlets 
## people with already > 0.5 chance of attending don't have doubled chance
## e.g. age 60 + have 1.5, 1.6 the odds 
## e.g. 0.48 * OR 1.5 = 0.58 
## e.g. 0.85 * OR 1.5 = 0.89 more like 1.5 increase 
## Qrisk > 15 0.63 to 0.91, 1.3 


### FIXME

