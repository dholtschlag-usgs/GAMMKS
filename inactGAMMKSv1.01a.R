# This unified version was developed to produce more robust estimates by bootstrapping
# Clear memor
rm(list=ls());
path  <- 'C:/Home/Projects/LoadEstimation/Data/';
load(file=paste(path,"sUniqStaNumberParmCode.RData",sep=""));
colNameVec  <- NULL;
# Load coda library for read.and.check function
library("coda", lib.loc="C:/Users/dholtsch/Documents/R/win-library/3.0")
#
Analyst     <- "dholtschlag"
runDateTime <- Sys.time();
titleAnal   <- "Test Sin Smooth With CC Basis"
#
# Specify minimum number of observations to estimate the more complex equation
### >>>User Input<<<
nObsThres <- read.and.check("Minimum number of observations for extended form of equation. Default [100]: ",
                            what=numeric(), default=100);
# if (as.numeric(nObsUser)>25){nObsThres <- as.numeric(nObsUser)}
#
# Filename substring to ID run
fName <- 'WeightsAssignedTestUni01b'; # readline("Enter string to help ID analysis in file name: ")
# Setup output file
metaFile    <- paste(fName,nObsThres,sep="");
meta <- file(paste(path,metaFile,".txt",sep=""),open="wt")
cat("TITLE: ",titleAnal," by ",Analyst," on ",as.character(as.POSIXlt(runDateTime)),"\n",
    file = meta, sep = "")
#
biasCorr <- "c"
# Interactive user input for bias correction method
# biasCorr <- readline("Bias Correction Method: Enter [C/c] for Cohn-Finney or [S/s] for Smearing. Default [C]: ");
if (toupper(biasCorr)=="C"){
  cat("The Cohn-Finney bias correction method was specified.\n",file=meta);
  biasCorrect <- "CohnFinney";
} else if (toupper(biasCorr)=="S"){
  cat("The Smearing bias correction method was specified.\n",file=meta);
  biasCorrect <- "Smearing";
} else {
  biasCorrect <- "CohnFinney";
  cat("The default Cohn-Finney bias correction method was assigned:\n",file=meta);
}
# Load the measurement data set
# import wqtruth by water year data
file <- 'wqtruth_wyear.txt'; 
# Create sample dataframe (sdfrm)
targetWyearLoads    <- read.table(paste(path,file,sep=""),sep='\t',header=TRUE);
# 
#
# Load libraries used in all the analysis
library("mgcv",       lib.loc="C:/Program Files/R/R-3.0.1/library");
library("lubridate",  lib.loc="C:/Users/dholtsch/Documents/R/win-library/3.0");
library("XML",        lib.loc="C:/Users/dholtsch/Documents/R/win-library/3.0");
library("USGSwsBase", lib.loc="C:/Users/dholtsch/Documents/R/win-library/3.0");
library("MCMCglmm",   lib.loc="C:/Users/dholtsch/Documents/R/win-library/3.0");
# Source the finneyGM function
source("finneyGM.R");
# source("fGAMMKSnObsv1.00a.R");
#
# biasCorrect <- "Smearing"
# biasCorrect <- "CohnFinney"
# cat("Bias Correction: ",biasCorrect,"\n",file=meta);
#
#

i <- read.and.check("Enter index of station/constituent[1:22]: ",what=numeric(),
                    answer.in=c(seq(1,22)),default=1)
# for (i in 1:22) {
ichar <- as.character(i);
# Select site, constituent, and model formula
switch(ichar,
       "1" = {StaNumber  <- "01463500";   StaName <- "Delaware River at Trenton NJ";
              parmCode   <- "SSC";       parmName <- "Suspended-Sediment";
              ParmUnits  <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "2" = {StaNumber <- "01463500C"; StaName <- "Delaware River at Trenton NJ";
              parmCode <- "SC"; parmName <- "Specific Conductance"; 
              ParmUnits <- "microSeimens per centimeter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "3" = {StaNumber <- "01646500"; StaName <- "Potomac River near Wash., DC Little Falls Pump Station";
              parmCode <- "SC"; parmName <- "Specific Conductance";  
              ParmUnits <- "microSeimens per centimeter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "4" = {StaNumber <- "03150000";  StaName <- "Muskingum River at McConnelsville OH";
              parmCode <- "NO23";       parmName <- "Nitrate plus Nitrite Nitrogen";    
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "5" = {StaNumber <- "03150000";  StaName <- "Muskingum River at McConnelsville OH"; 
              parmCode  <- "TN";        parmName <- "Total Nitrogen";                   
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "6" = {StaNumber <- "03150000";  StaName <- "Muskingum River at McConnelsville OH"; 
              parmCode  <- "TP";         parmName <- "Total Phosphorus";                 
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "7" = {StaNumber <- "03271601";  StaName <- "Great Miami River below Miamisburg OH";
              parmCode  <- "NO23";       parmName <- "Nitrate plus Nitrite Nitrogen";    
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "8" = {StaNumber <- "03271601";  StaName <- "Great Miami River below Miamisburg OH";
              parmCode  <- "TN";         parmName <- "Total Nitrogen";                   
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "9" = {StaNumber <- "03271601";  StaName <- "Great Miami River below Miamisburg OH";
              parmCode  <- "TP";         parmName <- "Total Phosphorus";                 
              ParmUnits <- "milligrams per liter";
              GAMformula  <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              LINformula  <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
       },
       "10"= {StaNumber <- "03428200";  StaName <- "WEST FORK STONES RIVER AT MURFREESBORO,TN";
              parmCode  <- "SC";   parmName <- "Specific Conductance";
              ParmUnits <- "microSeimens per centimeter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "11"= {StaNumber <- "04193500";  StaName <- "Maumee River at Waterville OH";
              parmCode  <- "TP";        parmName <- "Total Phosphorus";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "12"= {StaNumber <- "04193500S"; StaName <- "Maumee River at Waterville OH"; 
              parmCode  <- "SSC";  parmName     <- "Suspended-Sediment"; 
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "13"= {StaNumber <- "04197100";  StaName <- "Honey Creek at Melmore OH";
              parmCode  <- "NO23"; parmName <- "Nitrate plus Nitrite Nitrogen";    
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "14"= {StaNumber <- "04197100";  StaName <- "Honey Creek at Melmore OH";
              parmCode  <- "TN";   parmName <- "Total Nitrogen";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "15"= {StaNumber <- "04197100";  StaName <- "Honey Creek at Melmore OH";
              parmCode  <- "TP";   parmName <- "Total Phosphorus";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "16"= {StaNumber <- "04197170";  StaName <- "Honey Creek at Melmore OH";
              parmCode  <- "TP";   parmName <- "Total Phosphorus";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "17"= {StaNumber <- "05057000";  StaName <- "SHEYENNE RIVER NR COOPERSTOWN, ND";
              parmCode  <- "SC";   parmName     <- "Specific Conductance";
              ParmUnits <- "microSeimens per centimeter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "18"= {StaNumber <- "05325000";  StaName <- "MINNESOTA RIVER AT MANKATO, MN";
              parmCode  <- "SSC";  parmName <- "Suspended-Sediment";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "19"= {StaNumber <- "05474000";  StaName <- "Skunk River at Augusta, IA";
              parmCode  <- "SSC";       parmName <- "Suspended-Sediment";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "20"= {StaNumber <- "05554500";  StaName <- "VERMILION RIVER AT PONTIAC, IL";
              parmCode  <- "NO23";       parmName     <- "Nitrate plus Nitrite Nitrogen";
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "21"= {StaNumber <- "09163500";  StaName <- "COLORADO RIVER NEAR COLORADO-UTAH STATE LINE";
              parmCode  <- "SC";       parmName <- "Specific Conductance";
              ParmUnits <- "microSeimens per centimeter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       },
       "22"= {StaNumber <- "11447650";  StaName <- "SACRAMENTO R A FREEPORT CA";
              parmCode  <- "SSC";      parmName <- "Suspended-Sediment"; 
              ParmUnits <- "milligrams per liter";
              GAMformulaS <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
              GAMformulaC <- as.formula('lnWq_value ~ s(logDflow,sinDOY, bs="ts") + s(anom_1day_30day, bs="ts") + s(dectime, bs="ts")')
              LINformulaS <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
              LINformulaC <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + anom_1day_30day + fourier(dectime) + dectime");
       }
)
paste("Index",i,"selected for",parmName,"at",StaNumber,StaName,sep=" ")
### 
strat <- read.and.check("Enter strategy code: L80 H8070 S U: ", what=character(),
                        answer.in=c("L80","H8070","S","U"),default="L80")
cat("Strategy is",strat)
freqy <- read.and.check("Enter frequency: 6 12 24 52: ", what=numeric(),
                        answer.in=c(6,12,24,52),default=52)
cat("Frequency is",freqy)
frlgh <- read.and.check("Enter furlough code: 00 E40 M40 T40 ", what=character(),
                  answer.in=c("00","E40","M40","T40"),default="00")
cat("Furlough is",frlgh)
rpttn <- read.and.check("Enter repetition[1:10]: ", what=numeric(),
                        answer.in=c(seq(1,10)),default=1)
cat("Repetition is",rpttn)

# ID and write site and constituent to meta file
cat("Analysis: ",ichar," for ",file = meta)
siteName <- paste(StaNumber,StaName);
cat("site: ",siteName," and ", sep="",file = meta)
cnstName <- paste(parmCode,parmName);
cat("parameter: ",cnstName,"\n", sep="", file = meta)
# 
# Write GAM terms to meta file
cat("GAM Formula: ",as.character(attr(terms(GAMformulaC),"variables"))[-1],"\n",file=meta)
# Extract GAM terms for analysis
GAMExpTerms <- as.character(attr(terms(GAMformulaC), "variables"))[-(1:2)];
GAMnTerms   <- length(GAMExpTerms);
#
# Write LIN terms to meta file
cat("LIN Formula: ",as.character(attr(terms(LINformulaC),"variables"))[-1],"\n",file=meta)
# Extract LIN terms for analysis
LINExpTerms <- as.character(attr(terms(LINformulaC), "variables"))[-(1:2)];
LINnTerms   <- length(LINExpTerms);
#
selNdx        <- which(targetWyearLoads$wq_station_id == StaNumber & 
                         targetWyearLoads$parm_code == parmCode)
selWyearLoads <- targetWyearLoads[selNdx,];
#
# biasCorrect <- "Smearing"
# biasCorrect <- "CohnFinney"
# cat("Bias Correction Method: ",biasCorrect,"\n",file=meta);
#
pathSub <- paste("StaNumber/",StaNumber,"/",parmCode,"/",sep="");
load(paste(path,pathSub,"sdfrm.RData",sep=""));
load(paste(path,pathSub,"cdfrm.RData",sep=""));
#
cat("Running GAMMKS v1.00...\n")   
# Indentify unique strategies, frequencies, furloughs, and repetitions
uniqStrategy     <- unique(sdfrm$Strategy);
uniqFrequency    <- unique(sdfrm$Frequency);
uniqFurlough     <- unique(sdfrm$Furlough); # ndx <- which(uniqFurlough=="00"); uniqFurlough[ndx] = "N00";
uniqRepetition   <- unique(sdfrm$Repetition);
#
# Determine number of days in the continuous (daily) dataframe 
nDays                <- nrow(cdfrm);
popMeanLogDflow      <- mean(cdfrm$logDflow)
# 
# subset the sample data for estimation
sdfrmTestSet     <- subset(sdfrm,subset=c(Strategy==strat & Frequency==freqy & 
                                            Furlough==frlgh & Repetition==rpttn));
# Drop wq_code 
sdfrmTestSet     <- subset(sdfrmTestSet,select=-c(wq_code));
# Extract just the Date field from the daily value dataframe 
cdfrmDate         <- subset(cdfrm,select=Date); 
#
# Form dataframe of daily values with only sampled values
mdfrm             <- merge(cdfrmDate,subset(sdfrmTestSet,select=c(Date,lnWq_value)),all=TRUE);
# Form matrix of sampled values
SAMPmeasLnConcDaily  <- as.matrix(subset(mdfrm,select="lnWq_value")); 
#
# determine the sample size for model estimation
nObs     <- nrow(sdfrmTestSet);
#
# Enter the number of bootstrap samples
nBoot    <- read.and.check("Enter number of bootstap samples: ",what=numeric(),
                           lower=1,upper=10000,default=500)
cat("User specified",nBoot,"bootstrap samples.")

samMeanlogDflow <- rep(NA,nBoot)
#
# Allocate matrix to contain all load estimates for GAMMKS
GAMMKSpredLoad <- matrix(NA,nBoot,22);
#
# Allocate Scenario matrix for analysis summary
modelStatSummary           <- matrix(NA,nBoot,12);
hObs                       <- rep(NA,nBoot)
colnames(modelStatSummary) <- c("Strategy","Frequency","Furlough","Repetition",
                                "GAMadjR2","GAMrmse","GAMmodelDF","GAMresidDF",
                                "LINadjR2","LINrmse","LINmodelDF","LINresidDF");
#
modelStatSummary[,1] <- strat;
modelStatSummary[,2] <- as.character(sprintf("%02d",freqy));
modelStatSummary[,3] <- frlgh;
modelStatSummary[,4] <- rpttn;
# 
if (nObs >= nObsThres ){
  GAMformula <- GAMformulaC;
  LINformula <- LINformulaC;
}  else {
  GAMformula <- GAMformulaS;
  LINformula <- LINformulaS;
}
# 
# Avoid duplicating fields in merge by intersecting column names in the sampled and continuous dframes
# Get column names from the sample data
sColNames   <- colnames(sdfrmTestSet);
# get column names from the continuous data
cColNames   <- colnames(cdfrm);
# merge the sample subset with the continuous data by common columns to get anom_1day_30day
mdfrmTestSet <- merge(sdfrmTestSet,cdfrm,by=c(intersect(sColNames,cColNames)))
# Select only the Date and lnWq_value 
mdfrmTestSet <- subset(mdfrmTestSet, select=c("Date","lnWq_value","sinDOY",
                                              "cosDOY","anom_1day_30day","dectime","logDflow","dflow"));
# Sort by flow magnitude
# mdfrmTestSet <- mdfrmTestSet[order(mdfrmTestSet$logDflow),];
# colnames(mdfrmTestSet)  <- c("Date","sampleLnQW");
sampMeanLogDflow <- mean(mdfrmTestSet$logDflow)
# Allocate load estimation by water-year matrix
wySpan                      <- min(cdfrm$wyear):max(cdfrm$wyear);
GAMMpredLoadWYear           <- matrix(NA,nBoot,length(wySpan));
colnames(GAMMpredLoadWYear) <- c(wySpan);
#
# Time execution of the run
cat('Start loop over all data subsets\n')
time.start <- Sys.time();
# Allocate matrices to contain logs predicted and sampled concentrations
GAMpredLnConcDaily  <- matrix(NA,nDays,nBoot);
LINpredLnConcDaily  <- matrix(NA,nDays,nBoot);
# SAMPmeasLnConcDaily <- matrix(NA,nDays,nBoot);
GAMMpredLnConcDaily <- matrix(NA,nDays,nBoot);
GAMMpredConcDaily   <- matrix(NA,nDays,nBoot);
GAMMpredLoadDaily   <- matrix(NA,nDays,nBoot);
samMeanLogDflow     <- matrix(NA,nBoot,1);
#
#
# # Specify standard formula for the GAM (Notice double quotes inside single quotes)
# GAMformula      <- as.formula('lnWq_value ~ s(logDflow,bs="ts") + s(sinDOY, bs="ts") + s(dectime, bs="ts")');
# Extract the explanatory terms
# GAMExpTerms     <- as.character(attr(terms(GAMformula), "variables"))[-(1:2)];
# Count the number of smooth explanatory terms
GAMnTerms       <- length(GAMExpTerms);
#
# Initialize GAM parameter matrices
# GAM effective degrees of freedom for each smoothed term
GAMsTermEdf      <- matrix(NA,nrow=nBoot,ncol=GAMnTerms); colnames(GAMsTermEdf)  <- GAMExpTerms;
# GAM Reference degrees of freedom
GAMsTermRdf      <- matrix(NA,nrow=nBoot,ncol=GAMnTerms); colnames(GAMsTermRdf)  <- GAMExpTerms;
# GAM p-values for smooth terms
GAMsTermPval    <- matrix(NA,nrow=nBoot,ncol=GAMnTerms);  colnames(GAMsTermPval) <- GAMExpTerms;
#
# # Specify formula for linear model
# LINformula      <- as.formula("lnWq_value ~ logDflow + I(logDflow^2) + fourier(dectime) + dectime");
# Extract the explanatory terms
LINExpTerms     <- as.character(attr(terms(LINformulaC), "variables"))[-(1:2)];
# Count the explanatory terms
LINnTerms       <- length(LINExpTerms);
# Find the position of the fourier term in the explanatory variable list
fourierTerm      <- which("fourier(dectime)"==LINExpTerms);
#
if (!is.na(fourierTerm)) {
  LINExpTerms   <- c(LINExpTerms[1:fourierTerm-1],"fourier(dectime)sin(k=1)","fourier(dectime)cos(k=1)",
                     LINExpTerms[(fourierTerm+1):LINnTerms])
}
# Prepend the intercept term 
LINExpTerms   <- c("(Intercept)",LINExpTerms); 
# Update the number of terms to reflect the intercept and two individual fourier terms (sin,cos) 
LINnTerms     <- length(LINExpTerms);  
# Initialize Parameters Estimate matrices
LINParmEst    <- matrix(NA,nrow=nBoot,ncol=LINnTerms); colnames(LINParmEst)  <- LINExpTerms;
LINParmStd    <- matrix(NA,nrow=nBoot,ncol=LINnTerms); colnames(LINParmStd)  <- LINExpTerms;
LINParmPval   <- matrix(NA,nrow=nBoot,ncol=LINnTerms); colnames(LINParmPval) <- LINExpTerms;
#
h            <- 0;
# 
logTF        <- FALSE
for (iInt in 1:nBoot){
#  for (iInt in 1:2         ){
  h <- h + 1;
  if (iInt>1) logTF <- TRUE
  # beginning index
  # begNdx   <- (iInt-1)*minSamples + 1;
  # ending index is the 
#   if( (iInt+1)*minSamples < nObs) { endNdx <- iInt*minSamples } else { endNdx <- nObs }
#   hObs[h] <- length(begNdx:endNdx)
#   meanlogDflow[iInt] <- mean(mdfrmTestSet$logDflow[begNdx:endNdx])
#   cat('For index ',iInt,' range ',begNdx,':',endNdx,', the mean logflow is ',meanlogDflow[iInt],':',sep="")   
  # Standard Model
  ### Estimate Gaussian Model Form
  # Select a sample 
  ndxBootStrap <- sample(seq(1,nrow(mdfrmTestSet),by=1),nrow(mdfrmTestSet),replace=logTF)
  samMeanlogDflow[iInt] <- mean(mdfrmTestSet[ndxBootStrap,"logDflow"])
  GAMStdModel  <- gam(GAMformula,family=gaussian(link="identity"),gamma=1.4,
                     weights=NULL,data=mdfrmTestSet[ndxBootStrap,]);
  sumGAMStdModel    <- summary(GAMStdModel);
  # print(sumGAMStdModel)
  #
  # Summarize smooth-term specific info from 
  GAMsTermEdf[h,]   <- sumGAMStdModel$edf;      
  # GAM Reference degrees of freedom
  GAMsTermRdf[h,]   <- sumGAMStdModel$s.table[,2];
  # GAM p-values for smooth terms
  GAMsTermPval[h,]    <- sumGAMStdModel$s.pv;
  #
  # Store r2 for model in modelStatSummary matrix
  modelStatSummary[h,"GAMadjR2"]    <- sumGAMStdModel$r.sq;         # Adjusted r2, actually
  modelStatSummary[h,"GAMrmse"]     <- sqrt(GAMStdModel$sig2);
  modelStatSummary[h,"GAMmodelDF"]  <- sum(sumGAMStdModel$edf)+1;
  modelStatSummary[h,"GAMresidDF"]  <- sumGAMStdModel$residual.df;
  #
  ## Linear Model
  LINFulModel       <- lm(LINformula, weights=NULL,data=mdfrmTestSet[ndxBootStrap,]);
  # summary(LINFulModel);
  # 
  LINSelModel       <- step(LINFulModel,trace=0);
  sumLINSelModel    <- summary(LINSelModel);
  # 
  ndxLINSelModel                <- match(names(LINSelModel$coefficients),LINExpTerms);
  # Remove NAs from indices
  # ndxLINSelModel                <- ndxLINSelModel[!is.na(ndxLINSelModel)];
  LINParmEst[ h,ndxLINSelModel] <- coefficients(sumLINSelModel)[,1];
  LINParmStd[ h,ndxLINSelModel] <- coefficients(sumLINSelModel)[,2];
  LINParmPval[h,ndxLINSelModel] <- coefficients(sumLINSelModel)[,4];
  
  # Get model terms
  LINSelModel.trm    <- as.character(attr(terms(model.frame(LINSelModel)), "variables"))[-(1:2)]
  ndxLnQ             <- match("logDflow",LINSelModel.trm);
  ndxLnQ2            <- match("I(logDflow^2)",LINSelModel.trm);
  sumLINSelModel     <- summary(LINSelModel);
  #
  # Determine whether or not linear and quadratic terms were estimate
  modelStatSummary[h,"LINadjR2"]    <- sumLINSelModel$r.squared;         
  modelStatSummary[h,"LINrmse"]     <- sumLINSelModel$sigma;
  modelStatSummary[h,"LINmodelDF"]  <- sumLINSelModel$df[1];
  modelStatSummary[h,"LINresidDF"]  <-    LINSelModel$df.residual;
  # 
  # Print summary of selected LINgau model
  cat(h,modelStatSummary[h,"Strategy"],modelStatSummary[h,"Frequency"],modelStatSummary[h,"Furlough"],modelStatSummary[h,"Repetition"],
      "GAMr2:",format(as.numeric(modelStatSummary[h,"GAMadjR2"]),digits=4),
      "GAMrmse:",format(as.numeric(modelStatSummary[h,"GAMrmse"]),digits=4),
      "GAMmodelDF:",format(as.numeric(modelStatSummary[h,"GAMmodelDF"]),digits=4),
      "LINr2:",format(as.numeric(modelStatSummary[h,"LINadjR2"]),digits=4),
      "LINrmse:",format(as.numeric(modelStatSummary[h,"LINrmse"]),digits=4),
      "LINmodelDF:",format(as.numeric(modelStatSummary[h,"LINmodelDF"]),digits=4),
      "nObs",nObs,"\n");
  #
  # Store the original values of cdfrm$dectime -> decTime 
  decTime                 <- cdfrm$dectime;
  fracDate                <- cdfrm$dectime - floor(cdfrm$dectime);
  yearDate                <- cdfrm$dectime - fracDate;
  minYear                 <- year(min(mdfrmTestSet$Date));
  ndxLTmin                <- which(yearDate < minYear);
  yearDate[ndxLTmin]      <- minYear;
  cdfrm$dectime[ndxLTmin] <- yearDate[ndxLTmin] + fracDate[ndxLTmin];
  #
  maxYear                 <- year(max(mdfrmTestSet$Date));
  ndxGTmax                <- which(yearDate > maxYear);
  yearDate[ndxGTmax]      <- maxYear;
  cdfrm$dectime[ndxGTmax] <- yearDate[ndxGTmax] + fracDate[ndxGTmax];
  # Populate 
  newData                 <- with(cdfrm,data.frame(logDflow,sinDOY,cosDOY,anom_1day_30day,dectime));
  # Restore the original values of cdfrm$dectime <- decTime
  cdfrm$dectime           <- decTime;
  #
  # Compute GAM model estimates and standard errors
  GAMlnConcDaily      <- predict(GAMStdModel,newdata=newData,type="response",se.fit=TRUE);
  # Compute LIN model estimates and standard errors
  LINlnConcDaily      <- predict(LINSelModel,newdata=newData,type="response",se.fit=TRUE);
  # GAM daily estimate
  GAMpredLnConcDaily[,h]  <- GAMlnConcDaily$fit;
  #
  #cat("GAMpredLnConcDaily",GAMpredLnConcDaily[1:12,h],"\n")   
  # LIN daily estimate will either OLS or Bayesian model predictions
  LINpredLnConcDaily[,h]  <- LINlnConcDaily$fit;
  # cat("LINpredLnConcDaily",LINpredLnConcDaily[1:12,h],"\n")
  # If the Bayesian model exists, use it. 
  if (exists("LINBayModel")){
    newBayData             <- as.matrix(cbind(rep(1,times=nrow(cdfrm)),cdfrm$logDflow,
                                              cdfrm$logDflow^2,cdfrm$sinDOY,cdfrm$cosDOY));
    # Extract parameters from model
    betaLINBayModel <- as.matrix(sumLINBayModel$solutions[,1],ncol=1);
    LINpredLnConcDaily[,h] <- newBayData %*% betaLINBayModel;
    # Once it has been used, delete it
    rm("LINBayModel","newBayData")
  }
  # Weighted Average
  GAMMpredLnConcDaily[,h] <- (GAMlnConcDaily$fit/GAMlnConcDaily$se.fit^2 + 
                                LINpredLnConcDaily[,h]/LINlnConcDaily$se.fit^2) /
    (1/(GAMlnConcDaily$se.fit^2) + 1/(LINlnConcDaily$se.fit^2));
  # Weighted Variance of Individual estimates
  GAMMpredLnVariDaily <- 1 / (1/LINlnConcDaily$se.fit + 1/LINlnConcDaily$se.fit^2);
  ## Variance of the mean
  GAMMpredLnVariModel <- 1 / (1/sumLINSelModel$sigma^2  + 1/GAMStdModel$sig2);
  # 
  # Compute the GAM-gaussian model estmate of variance for individual values
  #         vm             <- GAMMpredLnVariDaily/GAMMpredLnVariModel;
  #         minDailyVar    <- apply(cbind(GAMlnConcDaily$se.fit^2,LINlnConcDaily$se.fit^2),1,"min");
  #         vm             <- minDailyVar/min(sumLINSelModel$sigma^2,GAMStdModel$sig2);
  #         meanDailyVar     <- apply(cbind(GAMlnConcDaily$se.fit^2,LINlnConcDaily$se.fit^2),1,"mean");
  #         vm               <- meanDailyVar/mean(c(sumLINSelModel$sigma^2,GAMStdModel$sig2));
  #
  # Pooled residual degrees of freedom
  GAMMresidDF    <- (sumGAMStdModel$residual.df/GAMStdModel$sig2 + 
                       sumLINSelModel$df[2]/sumLINSelModel$sigma^2) /
    (1/GAMStdModel$sig2 + 1/(sumLINSelModel$sigma^2));
  #
  if (biasCorrect == "CohnFinney"){
    # Allocate and compute v0;        
    X      <- as.matrix(LINSelModel$model[-1]);
    Xaug   <- cbind(rep(1,nrow(X)),X);
    v0     <- matrix(NA,nrow=nrow(Xaug),ncol=1);
    InvXtX <- (t(Xaug) %*% Xaug)^-1;
    for (u in 1:nrow(X)){
      v0[u] <- t(Xaug[u,]) %*% InvXtX %*% Xaug[u,];
      
    }
    sigma2y_x <- sumLINSelModel$sigma^2;
    v0      <- LINlnConcDaily$se.fit^2/sigma2y_x;
    RSS     <- sum(residuals(LINSelModel)^2);
    m       <- LINSelModel$df.residual;
    n       <- sum(sumLINSelModel$df[1:2]);
    rterm   <- m*RSS / ((2 * m + n * v0)*m + RSS);
    #
    # Compute the argument of Finney's function
    argm    <- (GAMMresidDF+1)/(2*GAMMresidDF)*
      ((1-v0)*GAMMpredLnVariModel);
    #         argm         <- (GAMMresidDF+1)/(2*GAMMresidDF) * 
    #                         ((1-v0)*min(sumLINSelModel$sigma^2,GAMStdModel$sig2));
    #         argm         <- (GAMMresidDF+1)/(2*GAMMresidDF) * 
    #                         ((1-v0)*mean(c(sumLINSelModel$sigma^2,GAMStdModel$sig2)));
    # Compute Finney's GM
    gm           <- finneyGM(GAMMresidDF,argm)
    # Simple inverse transform of predicted log concentrations
    # DWAgauStdPredConc   <- exp(DWAgauStdlnPred);
    # Finney-adjusted inverse transform of predicted log concentrations
    GAMMpredConcDaily[,h] <- exp(GAMMpredLnConcDaily[,h]) * gm;  # should this be changed?
    #GAMMpredConcDaily[,h] <- exp(GAMMpredLnConcDaily[,h] + rterm);
  } else if (biasCorrect == "Smearing") {
    compSig2 <- 2/(1/GAMStdModel$sig2 + 1/sumLINSelModel$sigma^2);
    GAMMpredConcDaily[,h] <- exp(GAMMpredLnConcDaily[,h]+compSig2/2);
  } else {
    cat("Bias correction technique *",biasCorrect,"* not recognized! Stopping.");
    stop;
  }
  # 
  # Compute the load using the adjusted concentrations 
  # The coef 2.44657555 converts mg/L * ft3/s to kg/day
  GAMMpredLoadDaily[,h] <- GAMMpredConcDaily[,h] * cdfrm$dflow * 2.44657555;
  # 
  # form dataframe for summary statistics
  dfLoadWY              <- data.frame(Load=GAMMpredLoadDaily[,h],wyear=as.factor(cdfrm$wyear));
  # store summary statistics in temporary dataframe
  ldfrm                 <- aggregate(Load ~ wyear, data=dfLoadWY,sum);
  # wyear loads are expressed in metric tons
  GAMMpredLoadWYear[h,]     <- ldfrm$Load;
  # Pause the action every 5 time steps
  # if (!(h %% 5)){
  # readline("Pausing execution. Press any key to continue. ")
  # }
  #           }
  #         }
  #       }
  #     }
}
# # # 
# ndxF00 <- which(modelStatSummary[,"Furlough"]=="00"); modelStatSummary[ndxF00,"Furlough"] <- "N00";
# time.end <- Sys.time();
# cat("Finished equation development in",format(time.end-time.start,digits=3),"\n" );
# #
# # Compute Kalman filter estimates 
# #
# # Compute innovation matrix: Note LnPredConc is full, LnSampConc is sparse (NA filled)
INNOcompLnConcDaily   <- log(GAMMpredConcDaily) - rep(SAMPmeasLnConcDaily,nBoot); 
# # Allocate big Kalman Smooth Matrix
GAMMKSpredLnConcDaily <- matrix(NA,nDays,nBoot);
GAMMKFpredLnConcDaily <- matrix(NA,nDays,nBoot);

GAMMKSpredLoadDaily   <- matrix(NA,nDays,nBoot);
time.start <- Sys.time();
# 
cat("Starting the Kalman filter/smoother for all scenarios.\n")
# Run daily values for all scenarios
for (h in 1:nBoot){
  # Initialize Kalman filter variables for each scenario
  xp    <- rep(0,nDays);  # State estimate of x(t) at x(t+)
  xm    <- rep(0,nDays);  # State estimate of x(t) at x(t-)
  yf    <- rep(0,nDays);  # Kalman filter estimate
  ys    <- rep(0,nDays);  # Kalman smoother estimate
  Pp    <- rep(0,nDays);  # Posterior estimate of state covariance
  Pm    <- rep(0,nDays);  # Prior estimate of the state covariance
  #
  # Initialize the Kalman filter
  xp[1] <- 0;               # Error at time zero
  Pp[1] <- 1;               # State covariance
  Rww   <- 1e-2;            # Process variance
  Rvv   <- 0.001;           # Measurement variance
  #
  # Kalman filter
  infoDecay <- 0.95;        # 
  for (k in 2:nDays){
    xm[k] <- infoDecay*xp[k-1];
    Pm[k] <- infoDecay^2 * Pp[k-1] + Rww;
    if(!is.na(INNOcompLnConcDaily[k])){
      # Apply measurement update
      # cat("Measurement update at k= ",k,"\n")
      Ree    <- Pm[k] + Rvv;
      Kgain  <- Pm[k] / Ree;
      xp[k]  <- xm[k] + Kgain * (INNOcompLnConcDaily[k,h] - xm[k]);
      Pp[k]  <- (1-Kgain)^2 * Pm[k] + Kgain^2 * Rvv;
    } 
    else {
      xp[k]  <- xm[k];
      Pp[k]  <- Pm[k];
    }
  }
  #
  # RTS Smoother 
  # Initialize smooth state and state covariance
  xs  <- xp;
  Ps  <- Pp;
  # 
  # Run smoother
  for (k in seq(nDays-1,1,-1)){
    a     <- Pp[k] * infoDecay  / Pm[k+1];
    xs[k] <- xp[k] + a*(xs[k+1] - xm[k+1]);
    Ps[k] <- Pp[k] + a*(Ps[k+1] - Pm[k+1])*a;
  }
  # Compute filtered ln conc daily estimates
  GAMMKFpredLnConcDaily[,h]  <- GAMMpredLnConcDaily[,h]-xp;
  # Compute smoothed ln conc daily estimates
  GAMMKSpredLnConcDaily[,h]  <- GAMMpredLnConcDaily[,h]-xs;
  # Compute smoothed load estimates
  GAMMKSpredLoadDaily[,h]    <- exp(GAMMKSpredLnConcDaily[,h]) * 
    cdfrm$dflow * 2.44657555;
}
#
# Allocate KSmoLoadWY for all sampling scenarios and water years
GAMMKSpredLoadWYear  <- matrix(NA,nBoot,length(wySpan));
# Sum KSmoLoads over water years
for (h in 1:nBoot){
  #    if (nObs[h]>=nObsThres){
  dfLoadWY$Load           <- GAMMKSpredLoadDaily[,h];
  tmp                     <- aggregate(Load ~ wyear, 
                                       data=dfLoadWY, sum);
  GAMMKSpredLoadWYear[h,] <- tmp$Load;
  #    }
}
time.end <- Sys.time();
cat("Time for Kalman filtering and smoothing was",format(time.end-time.start,digits=4),"\n");
#
# # Store results in a data frame
if (!exists("dfResult")){
  dfResult <- as.data.frame(modelStatSummary[,c("Strategy","Frequency","Furlough","Repetition")]);
}
dfResult[,paste(parmCode,"_",StaNumber,sep="")] <- rowMeans(GAMMKSpredLoadWYear);

plot(ecdf(log10(dfResult[1:nrow(dfResult),paste(parmCode,"_",StaNumber,sep="")])),
        main=paste("Cumulative Distribution Function of",parmName,"Load Estimates Based on",nBoot,
        "Bootstrap Samples \nat",StaNumber,StaName,"for Strategy",strat,"Frequency",
                   freqy,"Furlough",frlgh,"Repetition",rpttn),
        xlab="log10(load)",cex=0.5,col="blue",cex.main=0.75)
abline(v=log10(dfResult[1,paste(parmCode,"_",StaNumber,sep="")]),pch=16,col="red")
abline(v=mean(log10(selWyearLoads[,"wyload_kg"])),col="black",lty="dashed")
legend("topleft",legend=c("Bootstrap","Mean Bootstrap","Measured"),
       pch=c(16,NA,NA),col=c("blue","red","black"),lty=c(NA,"solid","dashed"))
#
# You don't know selWyearLoads
plot(samMeanlogDflow/popMeanLogDflow,log10(dfResult[,5]),
     pch=16,col="blue",cex=0.75,log="")
abline(v=samMeanlogDflow[1]/popMeanLogDflow,col="red",lty="dashed")
y1    <- log10(dfResult[,paste(parmCode,"_",StaNumber,sep="")])
x1    <- samMeanlogDflow/popMeanLogDflow
lmObj <- step(lm(y1 ~ x1))
abline(reg=lmObj,col="green")
# samMeanlogDflow is the average
predVal <- data.frame(x1=samMeanlogDflow[1]/popMeanLogDflow)
yhat    <- predict(lmObj,newdata=predVal)
# estY    <- 10^(log10(dfResult[1,5])/yhat)
cat("The load estimate is ",10^yhat," for a measured load of ",mean(selWyearLoads[,8]))
# 
# close(meta)
# write.table(dfResult,file=paste(path,"wgtNull01GAMMKS.csv",sep=""),sep=",")
#
denBLoad <- density(log10(dfResult[1:nrow(dfResult),
                    paste(parmCode,"_",StaNumber,sep="")]))
par(las=1)
plot(denBLoad,
     col="blue",type="l",cex.main=0.8,
     main=paste("Probability Density of",parmName,"Load Estimates Based on",nBoot,"Bootstrap Samples",
                "\nat",StaNumber,StaName,"for Strategy",strat,"Frequency",
                freqy,"Furlough",frlgh,"Repetition",rpttn),
     xlab="log10(load)")
# The first sample is constructed w/o replacement
abline(v=log10(dfResult[1,paste(parmCode,"_",StaNumber,sep="")]),pch=16,col="red")
abline(v=mean(log10(selWyearLoads[,"wyload_kg"])),col="black",lty="dashed")
legend("topleft",legend=c("Bootstrap","All samples","Measured"),
       pch=c(NA,NA,NA),col=c("blue","red","black"),lty=c("solid","solid","dashed"))
#
# Find indices of 0.025 and 0.975 quantiles
ndx025 <- max(which(denBLoad$y<=0.025))
save(denBLoad,parmCode,StaNumber,strat,freqy,frlgh,rpttn,nBoot,
     file=paste(parmCode,"_",StaNumber,"str",strat,"fre",freqy,"fur",frlgh,
                        "rep",rpttn,"nBoot",nBoot,"d.RData",sep=""))
