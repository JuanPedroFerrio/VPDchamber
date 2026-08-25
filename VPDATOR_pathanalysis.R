library(tidyverse)
# library(quantreg) # quantile regression
library(readxl) #para leer xls
# library(lubridate) #handle dates
# library(imputeTS) #gapfilling in time series
library(lavaan) # path analysis
library(semPlot) # plot Structural Equation Models
library(cowplot) # to combine panes

citation("lavaan")
citation("semPlot")
citation("tidyverse")

#### SET DIR ####
## Setting the working directory to the path of your current open file
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
# stores mainpath
mainpath<-getwd()


#### Function esat ####
# calculates saturated VP (mbar) from T (Celsius)
esat<-function(T=25){ # esat in mbars, from T in celsius
  esat<-6.13753*exp((T*(18.564-T/254.4))/(T+255.57))
  return(esat)
}

#### ALLCOMBINE #####

ALL<-read.csv("dataset.csv",header=TRUE)
ALL<-dplyr::rename(ALL,SWP=WP)
ALL<-filter(ALL,Sp!="coc")
ALL<-filter(ALL,SWP>-0.5)

ALL<-mutate(ALL,
            vp=esat(Tc)*HR/100,
            E=log(E),
            # # iWUE=A/gs,
            gs=log(gs),
            WUE=log(WUE)
            )
            
ALL<-summarise(ALL,
               .by=c(Sp,Timehour),
               across(c(PAR, vp, Tc, Tleaf, SWP, VPDleaf,gs,E,A,WUE),
                      mean,na.rm=TRUE))

ALL<-mutate(ALL,
             across(c(PAR, vp, Tc, Tleaf, SWP, VPDleaf,gs,E,A), scale))



#### OPTION1 NO gs ####

ENVmodel<-
  'E~VPDleaf+SWP+Tc
  A~PAR+VPDleaf+SWP+Tc
WUE~A+E
'

##### free model= fit separate paths for each species #####
results_free<-lavaan::sem(ENVmodel,
                  data=ALL,
                  group="Sp")

summary(results_free, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)

##### equal model= forces same path for the two species #####
results_equal<-sem(ENVmodel,data=ALL,
                   group="Sp",
                   group.equal="regressions")

summary(results_equal, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)

##### check if global differences between species #####
anova(results_free, results_equal)

##### check separate effects for each varaible #####

###### template~=free ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### E~VPDleaf ######
ENVmodel_var<-
  'E~c(VPDleafe,VPDleafe)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### E~SWP ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWPe,SWPe)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### E~Tc ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWPe,SWPe)*SWP+c(Tce,Tce)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### A~PAR ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PARa,PARa)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### A~VPDleaf ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleafa,VPDleafa)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### A~SWP ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWPa,SWPa)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### A~Tc ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tca,Tca)*Tc
WUE~c(A1,A2)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### WUE~A ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(Aw,Aw)*A+c(E1,E2)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

###### WUE~E ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PAR+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
WUE~c(A1,A2)*A+c(Ew,Ew)*E
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
  group.equal="regressions")
anova(results_free, results_var)

ENVorden <- matrix(c(
  NA,   NA,    NA,     NA,     NA,    "Tc",    NA,
  NA,   NA,      NA,       NA,       NA,    NA,    NA,
  NA,   NA,      NA,       NA,       NA,    NA,    NA,
  NA,   NA,     "VPDleaf", NA,       NA,     NA,    NA,
  NA,   NA,      NA,       NA,       NA,    NA,    "E",
  "SWP",NA,      NA,       NA,        NA,     NA,   NA,
  NA,   NA,      NA,       NA,      NA,     NA,    NA,
  NA,   NA,      NA,       NA,        NA,    NA,    NA,
  NA,   NA,      NA,       NA,        NA,    NA,    NA,
  NA,   "PAR",      NA,       NA,        "A",    NA,    NA,
  NA,   NA,      NA,       NA,        NA,    NA,    NA,
  NA,   NA,      NA,       NA,        NA,    NA,    NA,
  NA,   NA,      NA,       NA,        NA,   "WUE",   NA
), nrow = 13, byrow = TRUE)


tiff(paste("env_ILEFAG4_nogs.tiff"), width=5000, height=2500,res=600,units="px",compression="lzw")

semPaths(
  results_free,
  what = "std",
  whatLabels = "std",
  style = "lisrel",
  layout = ENVorden,
  reorder = FALSE,
  
  panelGroups = TRUE,   # ← THIS creates two panels
  
  # Node appearance
  sizeMan = 6,
  nCharNodes = 0,
  
  # Clean look
  intercepts = FALSE,
  residuals = FALSE,
  exoCov = FALSE,
  
  # Path appearance
  edge.label.cex = 0.8,
  edge.width = 1.2,
  curve = 0.2,
  curvePivot = FALSE,
  curvature=2,
  curveAdjacent=TRUE,
  
  # Colors
  color = "LightGrey",
  theme = "colorblind")

dev.off()

# 
# 
# #### OPTION2 gs ####
# # 
# ENVmodel<-
#   'Tleaf~Tc+E+PAR
#   VPDleaf~Tleaf+vp
#   gs~VPDleaf+SWP+PAR
#   A~PAR+Tleaf+gs
#   E~gs+VPDleaf
# WUE~A+E
# '
# 
# ##### free model= fit separate paths for each species #####
# results_free<-lavaan::sem(ENVmodel,
#                   data=ALL,
#                   group="Sp")
# 
# summary(results_free, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)
# 
# ##### equal model= forces same path for the two species #####
# results_equal<-lavaan::sem(ENVmodel,data=ALL,
#                    group="Sp",
#                    group.equal="regressions")
# 
# summary(results_equal, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)
# 
# ##### check if global differences between species #####
# anova(results_free, results_equal)
# 
# ##### check separate effects for each varaible #####
# 
# ###### template~=free ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### VPDleaf~Tc ######
# ENVmodel_var<-
#   'VPDleaf~c(Tca,Tca)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### gs~VPDleaf ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleafa,VPDleafa)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### gs~SWP ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWPa,SWPa)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### gs~PAR ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PARa,PARa)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### A~PAR ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PARa,PARa)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### A~Tc ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tca,Tca)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### A~gs ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gsa,gsa)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### E~VPDleaf ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleafa,VPDleafa)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### E~gs ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gsa,gsa)*gs
# WUE~c(A1,A2)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### WUE~A ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(Aa,Aa)*A+c(E1,E2)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ###### WUE~E ######
# ENVmodel_var<-
#   'VPDleaf~c(Tc1v,Tc2v)*Tc
#   gs~c(VPDleaf1g,VPDleaf2g)*VPDleaf+c(SWP1,SWP2)*SWP+c(PAR1g,PAR2g)*PAR
#   A~c(PAR1a,PAR2a)*PAR+c(Tc1a,Tc2a)*Tc+c(gs1a,gs2a)*gs
# E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(gs1e,gs2e)*gs
# WUE~c(A1,A2)*A+c(Ea,Ea)*E
# '
# results_var <- lavaan::sem(
#   ENVmodel_var,
#   data = ALL,
#   group = "Sp",
#   group.equal="regressions")
# anova(results_free, results_var)
# 
# ENVorden <- matrix(c(
#   NA,NA,    NA,     NA,     "SWP",    NA,    NA,
#   NA,   NA,      NA,       NA,       NA,    NA,    NA,
#   NA,   NA,      NA,       NA,       NA,    NA,    NA,
#   "vp",   NA,     "VPDleaf", NA,       NA,     "gs",    NA,
#   NA,   NA,      NA,       NA,       NA,    NA,    NA,
#   NA,   NA,      NA,       NA,        NA,     NA,   NA,
#   NA,   NA,      NA,       NA,        NA,    NA,    NA,
#   NA,   "PAR",      NA,       NA,        "A",    NA,    NA,
#   NA,   NA,      NA,       NA,     NA,    NA,    "E",
#   NA,   NA,        "Tleaf",       NA,        NA,    NA,    NA,
#   NA,   NA,      NA,       NA,        NA,    NA,    NA,
#   NA,   "Tc",      NA,       NA,        NA,    NA,    NA,
#   NA,   NA,      NA,       NA,        NA,   "WUE",   NA
# ), nrow = 13, byrow = TRUE)
# 
# tiff(paste("env_ILEFAG4_gs.tiff"), width=5000, height=2500,res=600,units="px",compression="lzw")
# 
# semPaths(
#   results_free,
#   what = "std",
#   whatLabels = "std",
#   style = "lisrel",
#   layout = ENVorden,
#   reorder = FALSE,
# 
#   panelGroups = TRUE,   # ← THIS creates two panels
# 
#   # Node appearance
#   sizeMan = 6,
#   nCharNodes = 0,
# 
#   # Clean look
#   intercepts = FALSE,
#   residuals = FALSE,
#   exoCov = FALSE,
# 
#   # Path appearance
#   edge.label.cex = 0.8,
#   edge.width = 1.2,
#   curve = 0.2,
#   curvePivot = FALSE,
#   curvature=2,
#   curveAdjacent=TRUE,
# 
#   # Colors
#   color = "LightGrey",
#   theme = "colorblind")
# 
# dev.off()
