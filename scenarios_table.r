## results from 50 runs of 25000, including subset trick 

## TODO

## replace SEs by stat unc for base case - true stat unc known at least for base case 
## stat unc measure includes some MC error already 
## So just need to replace SE by S
## ie calculate true S from within scenarios_combine using o'hagan method

## output is true mean_i  +  indep within run error_i
## variance of output is sum of these variances:
## statistical error + MC error 

## where do these come from.  SEs only presented currently for QALY gains and short term death reductions 


### BASE CASE AS FIRST SCENARIO WITH UNCERTAINTY, using 25k x 1000 run with different pars.  Comp unc is 3% of stat unc, so neglect.

source("scenarios_combine.r") ## TODO clean up workspace here 

j <- "IQALY"
Mbase <- M
resm <- round(M[,j,],1); ress <- round(SD[,3,],1)
resm <- paste0(resm, " (", ress, ")")
P <- Nall[,2:9] / Nall[,1]
P[1:7] <- 100 * P[1:7]
Pperc <- P[1:7]
perc <- paste0(round(Pperc, 1), "% (", signif(PSD[1:7]*100,1), ")")
nhc <- paste0(round(P[8], 1), " (", signif(PSD[8], 1), ")")
Mhcbase <- M[,"IQALY",1] / P["Number of HCs"]
resn <- round(Mhcbase, 1)
resns <- round(SD[,3,1] / P["Number of HCs"], 1)
resn <- paste0(resn, " (", resns, ")")
## no sds yet for process table 
short <- short %>%
  mutate(pd = paste0(round(deaths_pd,1), "% (", round(SHSD[,"deaths_pd"],1), ")")) %>%
  mutate(pd = ifelse(is.nan(deaths_pd), round(deaths_pd,1), pd))
foo <- short[short$year %in% c(5,10), c("scen","trt","year","pd")]
foo <- spread(foo, year, pd)
resshort <- unlist(c( ## can this be cleaner? reordered later 
    paste0(round(proc$perc[proc$outcome=="qhi"],1), "%"),
    subset(foo, trt=="statins")[,c("5","10")],
    paste0(round(proc$perc[proc$outcome=="att_nhbp"],1), "%"),
    subset(foo, trt=="aht")[,c("5","10")],
    paste0(round(proc$perc[proc$outcome=="hbmi"],1), "%"),
    subset(foo, trt=="wr")[,c("5","10")],
    paste0(round(proc$perc[proc$outcome=="att_smk"],1), "%"),
    subset(foo, trt=="sc")[,c("5","10")]
))

(res.base <- c(resm,perc,"Number of HCs"=nhc,resn,resshort))
Mbase[,"IQALY",]
length(res.base)





### ALL OTHER SCENARIOS

source("scenarios_combine.r")

j <- "IQALY"
resm <- round(M[,j,],1); ress <- signif(SE[,j,],1)
resm <- array(paste0(resm, " (", ress, ")"), dim=dim(resm))

P <- Nall[,2:9] / Nall[,1]
P[,1:7] <- 100 * P[,1:7]
Pperc <- P[,1:7]
perc <- array(paste0(round(Pperc, 1), "%"), dim=dim(Pperc), dimnames=dimnames(Pperc))
nhc <- round(P[,8], 1)
Mhc <- M[,"IQALY",1] / P[,"Number of HCs"]
resn <- round(Mhc, 1)
resns <- round(SE[,"IQALY",1] / P[,"Number of HCs"], 1)
resn <- paste0(resn, " (", resns, ")")

### GAINS VS NO HC AS WELL AS VS BASE 
Mabs <- M[-1,"IQALY",] + Mbase[rep(1,nsc-1),"IQALY",]
Mhcabs <- Mhc[-1] + Mhcbase
resmabs <- rbind(t(Mabs[,1:4]), Mhcabs)
resmabs <- round(resmabs[c(1,5,2,3,4),], 1)
resmabs <- resmabs[,c(1:4, 5:8, 15, 9:13)]

write.table(resmabs, file="hc_scenarios_vsnohc.txt", sep="\t", quote=FALSE)
## TODO gains per hc 

short <- short %>%
  mutate(pd = paste0(round(deaths_pd,1), "% (", round(deaths_pdse,1), ")")) %>%
  mutate(pd = ifelse(is.nan(deaths_pd), round(deaths_pd,1), pd))

foo <- short[short$year %in% c(5,10), c("scen","trt","year","pd")]
foo <- spread(foo, year, pd)
resshort <- cbind( ## can this be cleaner? reordered later 
    paste0(round(proc$perc[proc$outcome=="qhi"],1), "%"),
    subset(foo, trt=="statins")[,c("5","10")],
    paste0(round(proc$perc[proc$outcome=="att_nhbp"],1), "%"),
    subset(foo, trt=="aht")[,c("5","10")],
    paste0(round(proc$perc[proc$outcome=="hbmi"],1), "%"),
    subset(foo, trt=="wr")[,c("5","10")],
    paste0(round(proc$perc[proc$outcome=="att_smk"],1), "%"),
    subset(foo, trt=="sc")[,c("5","10")]
)

res <- t(cbind(resm, perc, "Number of HCs"=nhc,resn,resshort))

res[,1] <- res.base

rownames(res) <- c("qall", "qel","qatt","qoff","qstat","qaht","qsc","qwr",
                   "pel","patt","poff","pstat","paht","psc","pwr","nhc","qhc",
                   "hiq","d5stat","d10stat","hibp","d5aht","d10aht",
                   "hbmi","d5wr","d10wr","smk","d5sc","d10sc"
                   )
res <- res[c("qall","nhc","qhc",
             "pel","qel","patt","qatt","poff","qoff",
             "hiq","pstat","d5stat","d10stat","qstat",
             "hibp","paht","d5aht","d10aht","qaht",
             "hbmi","pwr","d5wr","d10wr","qwr",
             "smk","psc","d5sc","d10sc","qsc"
             ),]

rownames(res) <- c(
    "QALD gain (days; all)",
    "Mean number of HCs per person",
    "QALD gain (days; per HC)",
    "% eligible at least once",
    "QALD gain (days; eligible at least once)",
    "% attending at least once",
    "QALD gain (days; attending at least once)",
    "% offered treatment via HC" ,
    "QALD gain (days; offered treatment via HC)",
    "P(Qrisk>20) for first HC attenders",
    "% taking statins via HC" ,
    "reduction in deaths; 5 years after statins",
    "reduction in deaths; 10 years after statins",
    "Long term QALD gain from HC, statins population",
   
    "P(SBP>140) for first HC attenders",
    "% taking AHT via HC",
    "reduction in deaths; 5 years after AHT",
    "reduction in deaths; 10 years after AHT",
    "Long term QALD gain from HC, AHT population",

    "P(BMI>30) for first HC attenders",
    "% weight management via HC",
    "reduction in deaths; 5 years after WM",
    "reduction in deaths; 10 years after WM",
    "Long term QALD gain from HC, WM population",

    "P(smoker) for first HC attenders",
    "% smoking cessation via HC" ,
    "reduction in deaths; 5 years after SC",
    "reduction in deaths; 10 years after SC",
    "Long term QALD gain from HC, SC population"
)

write.table(res[, 1:5], file="hc_scenarios_results1.txt", sep="\t", quote=FALSE)
write.table(res[, c(6:9,16)], file="hc_scenarios_results2.txt", sep="\t", quote=FALSE)
write.table(res[, c(10:14)], file="hc_scenarios_results3.txt", sep="\t", quote=FALSE)

