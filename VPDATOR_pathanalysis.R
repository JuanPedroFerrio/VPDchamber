library(tidyverse)
# library(quantreg) # quantile regression
library(readxl) #para leer xls
# library(lubridate) #handle dates
# library(imputeTS) #gapfilling in time series
library(lavaan) # path analysis
library(semPlot) # plot Structural Equation Models
library(cowplot) # to combine panes
# library(lme4)
# library(piecewiseSEM)

citation("lavaan")
# citation("semPlot")
citation("tidyverse")
# citation("piecewiseSEM")

#### SET DIR ####
## Setting the working directory to the path of your current open file
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
# stores mainpath
mainpath<-getwd()


#### ALLCOMBINE #####

ALL<-read.csv("dataset.csv",header=TRUE)
ALL<-dplyr::rename(ALL,SWP=WP)
ALL<-filter(ALL,Sp!="coc")
# ALL<-filter(ALL,SWP>-0.5)

ALL<-dplyr::mutate(ALL,
            E=log(E),
            # # iWUE=A/gs,
            gs=log(gs),
            WUE=log(WUE)
            )
            
ALL<-dplyr::summarise(ALL,
               .by=c(Sp,IDplant,Timehour),
               dplyr::across(c(PARc, vpc, Tc, Tleaf, SWP, VPDleaf,gs,E,A,WUE),
                      mean,na.rm=TRUE))

ALL<-dplyr::mutate(ALL,
             dplyr::across(c(PARc, vpc, Tc, Tleaf, SWP, VPDleaf,gs,E,A),
                           ~ as.numeric(scale(.x))))

ALL$Sp <- factor(ALL$Sp)
ALL$IDplant <- factor(ALL$IDplant)

#### NEW MODEL EXCLUDING WUE AND ACCOUNTING FOR IDPLANT ####

ENVmodel<-
  'E~VPDleaf+SWP+Tc
  A~PARc+VPDleaf+SWP+Tc
'

##### free model= fit separate paths for each species #####
results_nocluster<-lavaan::sem(ENVmodel,
                          data=ALL,
                          group="Sp")

lavInspect(results_nocluster, "converged")

summary(results_nocluster, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)

##### free model= fit separate paths for each species #####
results_free<-lavaan::sem(ENVmodel,
                  data=ALL,
                  group="Sp",
                  cluster="IDplant")

lavInspect(results_free, "converged")

summary(results_free, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)

##### residual covariance #####
EA_cov <- lavaan::parameterEstimates(
  results_free,
  standardized = TRUE,
  ci = TRUE
) %>%
  dplyr::filter(
    op == "~~",
    (lhs == "E" & rhs == "A") |
      (lhs == "A" & rhs == "E")
  )

EA_cov

##### equal model= forces same path for the two species #####
results_equal<-sem(ENVmodel,data=ALL,
                   group="Sp",
                   group.equal="regressions",
                   cluster="IDplant")

summary(results_equal, standardized=TRUE, rsquare=TRUE,fit.measures = TRUE)

##### check if global differences between species #####
anova(results_free, results_equal)

##### check separate effects for each varaible #####

###### template~=free ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### E~VPDleaf ######
ENVmodel_var<-
  'E~c(VPDleafe,VPDleafe)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### E~VPDleaf without clustering ######
ENVmodel_var<-
  'E~c(VPDleafe,VPDleafe)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp")
anova(results_nocluster, results_var)

###### E~SWP ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWPe,SWPe)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### E~SWP without clustering ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWPe,SWPe)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp")
anova(results_nocluster, results_var)

###### E~Tc ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWPe,SWPe)*SWP+c(Tce,Tce)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### A~PARc ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PARa,PARa)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### A~VPDleaf ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleafa,VPDleafa)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### A~SWP ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWPa,SWPa)*SWP+c(Tc1a,Tc2a)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

###### A~Tc ######
ENVmodel_var<-
  'E~c(VPDleaf1e,VPDleaf2e)*VPDleaf+c(SWP1e,SWP2e)*SWP+c(Tc1e,Tc2e)*Tc
  A~c(PAR1a,PAR2a)*PARc+c(VPDleaf1a,VPDleaf2a)*VPDleaf+c(SWP1a,SWP2a)*SWP+c(Tca,Tca)*Tc
'
results_var <- lavaan::sem(
  ENVmodel_var,
  data = ALL,
  group = "Sp",
cluster="IDplant")
anova(results_free, results_var)

ENVorden <- matrix(c(
  NA,   NA,    NA,     NA,     NA,    "Tc",    NA,
  NA,   NA,      NA,       NA,       NA,    NA,    NA,
  NA,   NA,      NA,       NA,       NA,    NA,    NA,
  NA,   NA,     "VPDleaf", NA,       NA,     NA,    NA,
  NA,   NA,      NA,       NA,       NA,    NA,    "E",
  NA,   NA,      NA,       NA,        NA,     NA,   NA,
  "SWP",NA,      NA,       NA,        NA,     NA,   NA,
  NA,   NA,      NA,       NA,        NA,    NA,    NA,
  NA,   NA,      NA,       NA,        NA,    NA,    NA,
  NA,   NA,      NA,       NA,       "A",    NA,    NA,
  NA,   NA,     NA,       NA,        NA,    NA,    NA,
  NA,   "PARc",  NA,       NA,        NA,    NA,    NA,
  NA,   NA,      NA,       NA,        NA,    NA,   NA
), nrow = 13, byrow = TRUE)


tiff(paste("env_ILEFAG4_nogsnoWUE.tiff"), width=5000, height=2500,res=600,units="px",compression="lzw")

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
