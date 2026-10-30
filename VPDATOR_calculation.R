library(tidyverse)
library(quantreg) # quantile regression
library(ggforce) # expanded features for ggplot
library(readxl) #para leer xls
library(lubridate) #handle dates
library(imputeTS) #gapfilling in time series
library(cowplot) # to combine panes
library(drc) # regressions
library(aomisc) # selfstarting functions
library(latex2exp) # convierte expresiones de latex en plotmath con la función "TeX", 
library(DHARMa) #para evaluar calidad residuos
library(lmerTest) #para sacar anova lmer con pval
library(emmeans) #comparación medias

citation("tidyverse")
citation("imputeTS")
citation("lubridate")
citation("ggforce")
citation("quantreg")

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

#### FUNCTION DATACIRAS ####
DATACIRASsum <- function(filename){
  # column names and types from Sheet1 (CIRAScombi)
  coltypes1<-c("skip","date","skip","date","numeric",rep("skip",2),
               rep("numeric",6),rep("skip",2),
               rep("numeric",2),"text",rep("numeric",15),rep("skip",2),
               rep("numeric",21))
  colnames1<-c("TimeMIN","ExcelTime","Timelag_sec","CO2r","CO2a",
               "CO2d","H2Or","H2Oa","H2Od","Tirga","Pr","Status","Position",
               paste("CO2_",seq(1:6),sep=""),"REFCO2_7",
               paste("H2O_",seq(1:6),sep=""),"REFH2O_7",
               "PAR_ext","HR_ext","T_ext",
               paste("Tc_",seq(1:6),sep=""),
               paste("HR_",seq(1:6),sep=""),
               paste("PAR_",seq(1:6),sep=""))
  #
  DataCIRAS<-read_excel(filename,sheet=1,skip=1,
             col_names=colnames1,col_types=coltypes1)

  # convert status text column to boolean TRUE=ok FALSE=any message
  # this way the flagged lines can be identified in summary
  DataCIRAS$Status<-is.na(DataCIRAS$Status)
  #
  # Summarise to 1 min and position, i.e. one value per measurement
  # this makes the interpolation cleaner
  DataCIRAS<-DataCIRAS %>%
    group_by(TimeMIN,Position) %>%
    summarise_all(mean,na.rm=TRUE)
  #
  # Interpolates ref. values (REFCO2_7,REFH2O_7)
  DataCIRAS$REFCO2_7<-na_interpolation(DataCIRAS$REFCO2_7,option="linear")
  DataCIRAS$REFH2O_7<-na_interpolation(DataCIRAS$REFH2O_7,option="linear")
  
  # Keeps only when position > 0 and < 7 (measuring chambers = positions 1:6)
  DataCIRAS<- subset(DataCIRAS,Position>0 & Position<7)
  
  #
  # Creates "cycle" variable, idenfified as the time when the change in position is <0 (i.e. from 6 to 1)
  # This will be used to summarise by cycle. Instead of using a fixed time,
  # so we maximise the time resolution to the cycle duration
  DataCIRAS<-DataCIRAS %>% as_tibble() %>% dplyr::mutate(Cycle = cumsum( c(1, diff(Position)) < 0 ) )
  
  # Caculates differential values (dC, dH)
  # for CO2 we use REFCO2_7 (input air measured through analyzer IRGA), 
  # because CO2 is more sensitive to the lack of regular matching
  # for H2O we use H2Or (input air measured through reference IRGA)
  # because H2O is less sensitive to matching, and we observe a larger
  # memory effect for H2O measured in REFH2O_7
  DataCIRAS<-dplyr::mutate(DataCIRAS,
                    dC_1=CO2_1-REFCO2_7,dC_2=CO2_2-REFCO2_7,dC_3=CO2_3-REFCO2_7,
                    dC_4=CO2_4-REFCO2_7,dC_5=CO2_5-REFCO2_7,dC_6=CO2_6-REFCO2_7,
                    dH_1=H2O_1-H2Or,dH_2=H2O_2-H2Or,dH_3=H2O_3-H2Or,
                    dH_4=H2O_4-H2Or,dH_5=H2O_5-H2Or,dH_6=H2O_6-H2Or)
  
  # forces NAs in environmental variables when the given chamber is not measured
  # to do this, multiplies and divides the variable by "dC_x"
  # if dC_ is NA then variable will be NA, if not will remain the same
  # By doing this, the means for each of these variables will correspond to the
  # measurement time of each chamber after summarising by cycle
  # because the periods without measurement will be omitted
  #
  DataCIRAS<-dplyr::mutate(DataCIRAS,
                    Tc_1=Tc_1*dC_1/dC_1,Tc_2=Tc_2*dC_2/dC_2,Tc_3=Tc_3*dC_3/dC_3,
                    Tc_4=Tc_4*dC_4/dC_4,Tc_5=Tc_5*dC_5/dC_5,Tc_6=Tc_6*dC_6/dC_6,
                    HR_1=HR_1*dC_1/dC_1,HR_2=HR_2*dC_2/dC_2,HR_3=HR_3*dC_3/dC_3,
                    HR_4=HR_4*dC_4/dC_4,HR_5=HR_5*dC_5/dC_5,HR_6=HR_6*dC_6/dC_6,
                    PAR_1=PAR_1*dC_1/dC_1,PAR_2=PAR_2*dC_2/dC_2,PAR_3=PAR_3*dC_3/dC_3,
                    PAR_4=PAR_4*dC_4/dC_4,PAR_5=PAR_5*dC_5/dC_5,PAR_6=PAR_6*dC_6/dC_6)
 
  # adds VP ext and chamber
  DataCIRAS<-dplyr::mutate(DataCIRAS,
                    e_1=esat(Tc_1)*HR_1/100,e_2=esat(Tc_2)*HR_2/100,
                    e_3=esat(Tc_3)*HR_3/100,e_4=esat(Tc_4)*HR_4/100,
                    e_5=esat(Tc_5)*HR_5/100,e_6=esat(Tc_6)*HR_6/100,
                    e_ext=esat(T_ext)*HR_ext/100)
  # adds vpd ext and chamber
  DataCIRAS<-dplyr::mutate(DataCIRAS,
                    VPD_1=esat(Tc_1)-e_1,VPD_2=esat(Tc_2)-e_2,
                    VPD_3=esat(Tc_3)-e_3,VPD_4=esat(Tc_4)-e_4,
                    VPD_5=esat(Tc_5)-e_5,VPD_6=esat(Tc_6)-e_6,
                    VPD_ext=esat(T_ext)-e_ext)
  
  # Summarise by cycle, i.e. one value per measurement cycle
  DataCIRAS<-DataCIRAS %>%
    group_by(Cycle) %>%
    summarise_all(mean,na.rm=TRUE)

  # adds column with Time 10MIN,but NOT summarise to 10 min
  # used to combine with  other sets of variables (CR1000-FLOW,ZL6-Water Potential)
  DataCIRAS$Time10MIN<-round_date(DataCIRAS$TimeMIN,unit="10 minutes")
  #
  return(DataCIRAS)
  }

#### FUNCTION DATACR1000 ####
DATACR1000sum <- function(filename){
  # column names and types from Sheet2 (CR1000)
  coltypes2<-c("skip","skip","date",rep("numeric",28),
               rep("skip",10))
  colnames2<-c("TimeCR1000","RECORD","BattV","PTemp_C",
               paste("Flow_",seq(1:6),sep=""),
               "RTemp_C",
               paste("T_tc",seq(1:17),sep=""),"T_Tc19")
  #
  DataCR1000<- read_excel(filename,sheet=2,skip=4,
                            col_names=colnames2,col_types=coltypes2)
  #
  # adds column with Time 10MIN AND summarise to 10 min
  DataCR1000$Time10MIN<-round_date(DataCR1000$TimeCR1000,unit="10 minutes")
  DataCR1000<-DataCR1000 %>%
    group_by(Time10MIN) %>%
    summarise_all(mean,na.rm=TRUE)
  return(DataCR1000)
  }
#

#### FUNCTION DATAZL6 ####
DATAZ6sum <- function(filename){
  # column names and types from Sheet1 (CIRAScombi)
  coltypes2<-c("date",rep("numeric",12),rep("skip",3),"numeric")
  colnames2<-c("Time10MIN",
               paste(c("WP_","TS_"),sort(c(seq(1:6),seq(1:6))),sep=""),
               "PkPa")
  #
  DataZ6<- read_excel(filename,sheet=1,skip=3,
                          col_names=colnames2,col_types=coltypes2)
  return(DataZ6)
}
#

####   BOUNDARY LAYER FUNCTION      ####

BOUNDARY<-function(wind=0.2,input=0.0004,Area=TRUE){
  
  # LeafArea= leaf area, in m2
  # d= Characteristic Dimension, in m
  # wind = wind speed, in m s-1
  
  # if Area=TRUE input = LeafArea, in m2
  # if Area=FALSE input = char. dimension (d) in m
  if (Area==TRUE) {
    LeafArea<-input # convert LA to m2
    # sensible heat conductance for an elliptical leaf,
    # after Ball et al 1988; gbH=1/(3.8*LeafArea^0.25*wind^-0.5)
    # rb = 1.78*rbH
    rb<-1.78*3.8*LeafArea^0.25*wind^(-0.5)
    } 
  else {
    d<-input # comvert d to m
    ### boundary layer from characteristic dimension and wind speed ###
    DiffW<-2.42*10^-5 # diffusivity of water in air (m2 s-1), after Jones 1992 appendix 2
    bT<-1000*(2*DiffW*((wind/d)^-0.5)/0.00662) # boundary layer thickness, in mm, after Jones 1992, p.63
    gb<-1000000*DiffW/bT # gb in mm s-1, from bT and diffussivity
    gb<-gb*0.04 #gb in mol m2 s-1
    # boundary resistance from char dimension d, in m2 s mol-1
    rb<-1/gb
    }

  return(rb)
  #  
} #function ends


#### TleafBarbourFUNCTION ####


TleafBarbour<-function(H=124,Tair=39,eout=45,P=980,E=0.0047,rb=2,maxiter=1000){

  
  ###### UNIVERSAL CONSTANTS ######
  
  #
  
  # H= Incident radiation absorbed by the leaf, in W m-2
  # eleaf= initial estimate of sat. vapour pressure in the leaf, in mbar
  # eout= vapour presssure in the cuvette, in mbar
  # E= transpiration rate, in mol m-2 s-1
  # rb = boundary layer resistance, in m2 s mol-1
  # P, pressure in mbar
  # maxiter, maximum number of iterations to converge
  
  # rho (Stefan-Boltzmann, W m-2 K-4)
  Rho<-5.67*10^-8
  # epsilon-Leaf = emissivity
  Emiss<-0.98
  # Calor espec?fico aire, por mol
  Cp<-29.2 # en J mol-1 K-1 # Barbour et al. 2000
  # Cp<-1.012 # en kJ kg-1 K-1 # CIRAS3 manual
  # Calor Latente de evaporaci?n; L (J mol-1)
  # LatHeat<-44012 # in Barbour et al. 2000
  LatHeat<-45064.3-(Tair*42.9) # formula in CIRAS 3 manual
  
  ###### INPUT VALUES #####
  
  # boundary layer resistance to water vapour (rb) and heat transfer (rbH) defined from chamber, in m2 s mol-1
  # rb is given as input parameter in the function, e.g. rb=0.4 for 2.5cm2 chamber in CIRAS, 
  # if outside depends on wind speed, leaf size and morphology
  # Here we calculate rbH from rb, with different ratios depending on stomatal distribution
  #
  rbH<-rb/1.78 # we use here the value for hypostomatous leaves from Ball et al.1988

  # gbH = 1 /rbH, conductance = inverse of resistance
  gbH<-1/rbH
  
  # Initial estimates of Tleaf and gtotal, using formulae in CIRAS 3 manual
  # A.6 Tleaf from Parkinson 1983 in CIRAS 3 manual
  #     Tl<-Tc_+(H-LatHeat*E)/(((0.93*Ma*Cp)/rb)+(4*Rho*(Tc_+273)^3))
  Tleaf=Tair+(H-LatHeat*E)/(Cp/rbH+4*Rho*(Tair+273)^3)
  # calculate eleaf from Tleaf, assuming saturation
  eleaf=esat(Tleaf)
  # (A.9) gtotal from E in mmol m-2 s-1
  gtotal=(E*(P-(eleaf+eout)/2)/(eleaf-eout))

  # Initialise convergence variables
  gold<-gtotal*1.1
  Told<-Tleaf*1.1
  # *1.1 forces "previous" values of convergence variables to be
  # different from "present" in the first iteration

  iteration<-0
  convergence<-FALSE
  while(mean(convergence,na.rm=TRUE)<1 & iteration<=maxiter)
    # this way it can check convergence for vector inputs
    # if convergence is FALSE in any cases mean <1
    # if convergence is TRUE in all cases mean=1
    # so it will run the iterations until ALL cases converge
  {

    ##### DeltaT after Barbour et al 2000 #####
    
    # total resistance (rt), 
    #
    rtotal <- 1/(gtotal)
    
    # radiative heat conductance (gr)
    # 
    gr <- (4*Rho*Emiss*(273+Tair)^3)/Cp
    
    # combined resistance to sensible and 
    # radiative heat transfer in parallel
    # r*bh (m2 s mol-1) = =1/(gbH+gr)
    # after dePury&Farquhar,unpublished, in Barbour et al. 2000
    #
    rxbh <- 1/(gbH+gr)
    
    # slope of the curve relating temperature 
    # to saturated vapour pressure s_mbar_C
    # after Postl and Bolh?r-Nordenkampf (1983), in Barbour et al. 2000
    #
    s_mbar_C <-6.13753*(((Tair+255.57)*(18.564-2*Tair/254.4)-Tair*(18.564-Tair/254.4))/
                          (Tair+255.57)^2)*exp(Tair*(18.564-Tair/254.4)/(Tair+255.57))
    # epsilon
    # s*L/(Cp*P) after Cowan 1977
    #
    Epsilon <- (s_mbar_C*LatHeat)/(Cp*P)
    
    # Molar leaf to air gradient
    Dmol <-(eleaf-eout)/P
    # DeltaT
    # modelled from Barbour et al. 2000
    #
    DeltaT <- (rxbh*(H*rtotal-LatHeat*Dmol)/
                 (Cp*(rtotal+Epsilon*rxbh)))
    Tleaf<-Tair+DeltaT
    
    # recalculates eleaf from Tleaf
    eleaf<- esat(Tleaf)
    
    # recalculates gtotal with new eleaf
    # (A.9) E from gtotal in mol m-2 s-1
    gtotal<-(E*(P-(eleaf+eout)/2)/(eleaf-eout))
   
    # convergence criteria, equal previous and current values with 6 digits
    convergence<-(round(gold,4)==round(gtotal,4))&(round(Tleaf,2)==round(Told,2))
    
    gold<-gtotal
    Told<-Tleaf
    iteration<-iteration+1
    
  } #while ends
  
  # keeps NAs for Tleaf in non-converged cases
  Tleaf<-case_when(convergence==FALSE ~ NA,
            TRUE ~ 1)*Tleaf
            
  return(Tleaf)
  #  
} #function ends


#### Function CIRAS ####
CIRAS <- function(Dataset=tibble(Time10MIN=1,P=980,eout=45,deltaH=12,
                  CO2out=338,deltaC=-20,
                  Flow=30,Tair=39,PAR=800,
                  d=1,LeafArea=2,Atotal=600,
                  Sp="ilex"),ID=1,wind=0.2){
  # INPUTS
  # eout, ein, P      mbar
  # CO2in deltaC      ppm
  # Flow              LPM
  # Atotal Leafarea   cm2
  # Tair Tleaf        celsius
  # Q                 PPFD micromol m-2 s-1
  # wind              m s-1
  #
  ##### CONSTANTS #####
  # rho (Stefan-Boltzmann, W m-2 K-4)
  Rho<-5.67*10^-8
  # epsilon-Leaf = emissivity
  Emiss<-0.98
  # Air specific heat, por mol
  Cp<-29.2 # in J mol-1 K-1 # Barbour et al. 2000
  # Evaporation latent heat; L (J mol-1)
  # LatHeat<-44012 # in Barbour et al. 2000
  Dataset<-dplyr::mutate(Dataset,LatHeat=45064.3-(Tair*42.9)) 
  # Sets LHE "light harvesting efficiency"
  # associated to crown arquitecture and self-shading
  # Esteso et al. 2006, Annals For. Sci
  Dataset<-dplyr::mutate(Dataset,
            LHE=case_when(startsWith(Sp,"ile") ~ 0.62,
                          startsWith(Sp,"coc") ~ 0.61,
                          startsWith(Sp,"fag") ~ 0.69,
                          TRUE ~ 1.00))
  # LHE<-1
#
  ##### INPUT DATA CONVERSION #####
  # Incident radiation absorbed by the leaf
  # (H=0.5 *PAR converting photon flux to energy
  # * 0.5 ratio of IR to visible 
  # for sun radiation, from Jones 1983, in Barbour et al. 2000
  # * LHE "light harvesting efficiency"
  # associated to crown arquitecture and self-shading
  # Esteso et al. 2006, Annals For. Sci
  Dataset<-dplyr::mutate(Dataset,H=PAR*0.5*0.5*LHE)

  # convert Total Leaf Area cm2 to m2
  Dataset<-dplyr::mutate(Dataset,Atotal=Atotal/10000)
  Dataset<-dplyr::mutate(Dataset,LeafArea=LeafArea/10000) # Leaf area
  # convert volume flow L/min to L/s
  Dataset<-dplyr::mutate(Dataset,V0=Flow/60)
  # ein = eout-deltaH
  Dataset<-dplyr::mutate(Dataset,ein=eout-deltaH)
  # CO2in = CO2out-deltaC
  Dataset<-dplyr::mutate(Dataset,CO2in=CO2out-deltaC)
  # CALCULATIONS
  # (A.1) volume flow in L/s to mass flow in mmol m-2 s-1)
  Dataset<-dplyr::mutate(Dataset,W=V0*(1/22.414)/Atotal)
  # (A.5) E in mol m-2 s-1 from eout, ein and W
  Dataset<-dplyr::mutate(Dataset,E=(W*(eout-ein))/(P-eout))
  # A.14 in CIRAS MANUAL
  Dataset<-dplyr::mutate(Dataset,A=-((deltaC*W)+(CO2out*E)))
  
  ###### boundary layer resistance from Leaf Area, hypostomatous, following Ball et al 1988 ######
  Dataset<-dplyr::mutate(Dataset,rb=BOUNDARY(wind,LeafArea)) # in m2 s mol-1
  ###### sensible heat conductance ######
  Dataset<-dplyr::mutate(Dataset,rbH=rb/1.78) # in m2 s mol-1
  
  ##### Tleaf folowing Barbour et al. #####
  Dataset<-dplyr::mutate(Dataset,Tleaf=TleafBarbour(H,Tair,eout,P,E,rb))
 
  ##### calculate eleaf from Tleaf, assuming saturation #####
  Dataset<-dplyr::mutate(Dataset,eleaf=esat(Tleaf))
  
  # (A.9) gtotal from E in mmol m-2 s-1
  Dataset<-dplyr::mutate(Dataset,gtotal=(E*(P-(eleaf+eout)/2)/(eleaf-eout)))
  # gs from gtotal and rb
  Dataset<-dplyr::mutate(Dataset,gs=1/(1/gtotal-rb))
  # E, gtotal and gs to mmol m-2 s-1
  Dataset<-dplyr::mutate(Dataset,E=E*1000,gtotal=gtotal*1000,gs=gs*1000)
  # gs as eleaf-eout
  Dataset<-dplyr::mutate(Dataset,VPDleaf=eleaf-eout)
  # Clean output
  # output filters
  output<-dplyr::filter(Dataset,gtotal>0 & gs<1500 & gs>0 & Flow > 20 & Flow <50)
  #output<-Dataset
  output<-dplyr::select(output,c(Time10MIN,Sp,E,A,Tleaf,VPDleaf,gtotal,gs))
  # rename columns to add ID
  colnames(output)[2:8]<-paste(colnames(output)[2:8],ID,sep="")
  # return(Dataset)
  return(output)
}

#### Function CIRAS2 ####
CIRAS2 <- function(Dataset=tibble(Time10MIN=1,P=980,eout=45,edil=45,deltaH=12,
                                 CO2out=338,deltaC=-20,
                                 Flow=30,Tair=39,PAR=800,
                                 d=1,LeafArea=2,Atotal=600,
                                 Sp="ilex"),ID=1,wind=0.2){
  
  # Special version of CIRAS used
  # to refine calculations using humidity sensor in camera to calculate eout
  # instead of the IRGA value for eout
  # in this case, we use the IRGA value (edil) to correct for the dilution effect
  # of transpiration on CO2
  # but the humidity sensor value to calculate the E and gs (eout)
  #
  # INPUTS
  # eout, edil, ein, P      mbar
  # CO2in deltaC      ppm
  # Flow              LPM
  # Atotal Leafarea   cm2
  # Tair Tleaf        celsius
  # Q                 PPFD micromol m-2 s-1
  # wind              m s-1
  #
  # CONSTANTS
  # rho (Stefan-Boltzmann, W m-2 K-4)
  Rho<-5.67*10^-8
  # epsilon-Leaf = emissivity
  Emiss<-0.98
  # Air specific heat, por mol
  Cp<-29.2 # in J mol-1 K-1 # Barbour et al. 2000
  # Evaporation latent heat; L (J mol-1)
  # LatHeat<-44012 # in Barbour et al. 2000
  Dataset<-dplyr::mutate(Dataset,LatHeat=45064.3-(Tair*42.9)) 
  # Sets LHE "light harvesting efficiency"
  # associated to crown arquitecture and self-shading
  # Esteso et al. 2006, Annals For. Sci
  Dataset<-dplyr::mutate(Dataset,
                  LHE=case_when(startsWith(Sp,"ile") ~ 0.62,
                                startsWith(Sp,"coc") ~ 0.61,
                                startsWith(Sp,"fag") ~ 0.69,
                                TRUE ~ 1.00))
  # LHE<-1
  #
  # INPUT DATA CONVERSION
  # Incident radiation absorbed by the leaf
  # (H=0.5 *PAR converting photon flux to energy
  # * 0.5 ratio of IR to visible 
  # for sun radiation, from Jones 1983, in Barbour et al. 2000
  # * LHE "light harvesting efficiency"
  # associated to crown arquitecture and self-shading
  # Esteso et al. 2006, Annals For. Sci
  Dataset<-dplyr::mutate(Dataset,H=PAR*0.5*0.5*LHE)
  
  # convert Total Leaf Area cm2 to m2
  Dataset<-dplyr::mutate(Dataset,Atotal=Atotal/10000)
  Dataset<-dplyr::mutate(Dataset,LeafArea=LeafArea/10000) # Leaf area
  # convert volume flow L/min to L/s
  Dataset<-dplyr::mutate(Dataset,V0=Flow/60)
  # ein = eout-deltaH
  Dataset<-dplyr::mutate(Dataset,ein=eout-deltaH)
  # CO2in = CO2out-deltaC
  Dataset<-dplyr::mutate(Dataset,CO2in=CO2out-deltaC)
  # CALCULATIONS
  # (A.1) volume flow in L/s to mass flow in mmol m-2 s-1)
  Dataset<-dplyr::mutate(Dataset,W=V0*(1/22.414)/Atotal)
  ##### (A.5) E in mol m-2 s-1 from eout, ein and W #####
  ##### USED FOR CALCULATION OF E AND GS #####
  Dataset<-dplyr::mutate(Dataset,E=(W*(eout-ein))/(P-eout))
  ##### (A.5) Edil in mol m-2 s-1 from edil, ein and W #####
  ##### Edil USED ONLY FOR DILUTION EFFECT ON CO2 #####
  Dataset<-dplyr::mutate(Dataset,Edil=(W*(edil-ein))/(P-eout))
  # A.14 in CIRAS MANUAL
  # NOTE THAT HERE WE USE Edil to account for dilution effect based on IRGA
  # Using eout from sensor would give a wrong estimate of the dilution effect due to memory effects
  Dataset<-dplyr::mutate(Dataset,A=-((deltaC*W)+(CO2out*Edil)))
  
  # sensible heat conductance for an elliptical leaf, 
  # # after Ball et al 1988; gbH=1/(3.8*LeafArea^0.25*wind^-0.5)
  #Dataset<-dplyr::mutate(Dataset,rbH=3.8*LeafArea^0.25*wind^(-0.5))
  # # boundary layer resistance from LA, hypostomatous, following Ball et al 1988
  #Dataset<-dplyr::mutate(Dataset,rb=1.78*rbH) # in m2 s mol-1
  #
  # boundary layer resistance from Leaf Area, hypostomatous, following Ball et al 1988
  Dataset<-dplyr::mutate(Dataset,rb=BOUNDARY(wind,LeafArea)) # in m2 s mol-1
  # sensible heat conductance 
  Dataset<-dplyr::mutate(Dataset,rbH=rb/1.78) # in m2 s mol-1
  # A.6 Tleaf from Parkinson 1983 in CIRAS 3 manual
  #     Tl<-Tc_+(H-LatHeat*E)/(((0.93*Ma*Cp)/rb)+(4*Rho*(Tc_+273)^3))
  #Dataset<-dplyr::mutate(Dataset,Tleaf=Tair+(H-LatHeat*E)/(Cp/rbH+4*Rho*(Tair+273)^3))
  Dataset<-dplyr::mutate(Dataset,Tleaf=TleafBarbour(H,Tair,eout,P,E,rb))
  
  # calculate eleaf from Tleaf, assuming saturation
  Dataset<-dplyr::mutate(Dataset,eleaf=esat(Tleaf))
  
  # ##### CIRAS TEST BEGIN #####
  # # modified formulae to replicate conditions in chamber CIRAS, used to test calculations
  # Cp<-1.012 # kJ kg-1 K-1 = J g K-1
  # Ma<-28.96 # air molar mass g mol-1
  # Dataset<-dplyr::mutate(Dataset,LatHeat=45064.3-(Tair*42.9)) 
  # Dataset<-dplyr::mutate(Dataset,H=Q*0.14) 
  # rb<-0.4 # m2 s mol-1 
  # rbH<-rb/1.78
  # Dataset<-dplyr::mutate(Dataset,Tleaf=Tair+(H-LatHeat*E)/((Ma*Cp)/rbH+4.639+0.5834*Tair)) 
  # Dataset<-dplyr::mutate(Dataset,eleaf=6.1365*exp((Tleaf*17.502)/(Tleaf+240.97))) 
  # ##### CIRAS TEST END #####
  
  # (A.9) gtotal from E in mmol m-2 s-1
  Dataset<-dplyr::mutate(Dataset,gtotal=(E*(P-(eleaf+eout)/2)/(eleaf-eout)))
  # gs from gtotal and rb
  Dataset<-dplyr::mutate(Dataset,gs=1/(1/gtotal-rb))
  # E, gtotal and gs to mmol m-2 s-1
  Dataset<-dplyr::mutate(Dataset,E=E*1000,gtotal=gtotal*1000,gs=gs*1000)
  # gs as eleaf-eout
  Dataset<-dplyr::mutate(Dataset,VPDleaf=eleaf-eout)
  # Clean output
  # output filters
  output<-dplyr::filter(Dataset,gtotal>0 & gs<1500 & gs>0 & Flow > 20 & Flow <50)
  #output<-Dataset
  output<-dplyr::select(output,c(Time10MIN,Sp,E,A,Tleaf,VPDleaf,gtotal,gs))
  # rename columns to add ID
  colnames(output)[2:8]<-paste(colnames(output)[2:8],ID,sep="")
  # return(Dataset)
  return(output)
}


#### READ DATA ####
#
# GETDATA
files<-str_subset(list.files(),"xlsx")
#
files

# starts dataframe with first file, adds one column with round
i<-1
DataCIRAS<-data.frame(round=1,DATACIRASsum(files[i]))
DataCR1000<-data.frame(round=1,DATACR1000sum(files[i]))
i<-2
for (i in 2:6){
DataCIRAS<-rbind(DataCIRAS,data.frame(round=i,DATACIRASsum(files[i])))
DataCR1000<-rbind(DataCR1000,data.frame(round=i,DATACR1000sum(files[i])))
}

#### combines CIRAS & CR1000 data ####
Alldata<-left_join(DataCIRAS,DataCR1000,by=c("round","Time10MIN"))

# add WP and Tsoil (ZL6)
i<-7
DataZ6<-data.frame(DATAZ6sum(files[i]))

for (i in 8:16){
  DataZ6<-rbind(DataZ6,data.frame(DATAZ6sum(files[i])))
}

# rearrange columns in DataZ6
DataZ6<-cbind(DataZ6[,c(1:2,4,6,8,10,12)],DataZ6[,c(3,5,7,9,11,13)])

#### combines CIRAS+CR1000 & ZL6 ####
Alldata<-left_join(Alldata,DataZ6,by="Time10MIN")

# Fills initial data for WP (before drought treatment, field capacity)
Alldata<-Alldata %>% mutate_at(vars(starts_with("WP_")),na_interpolation)

# For Chamber 1, MAKE 'NA' in  in the periods we detected failure of the fan
# (low flow measured coinciding with erratic values of dC_, 
# therefore we discarded an electrical problem of the sensor)
# After checking the fan, we found out that the fan blade was uncoupled
Alldata<-Alldata %>% dplyr::mutate(Flow_1 = case_when(
                                      Flow_1 < 27 ~ as.numeric(NA), 
                                      TRUE   ~ Flow_1))
#
# transfers "NA" from Flow_1 to other _1 variables
Alldata<-Alldata %>% dplyr::mutate(CO2_1=CO2_1*Flow_1/Flow_1,
                            H2O_1=H2O_1*Flow_1/Flow_1,
                            Tc_1=Tc_1*Flow_1/Flow_1,
                            HR_1=HR_1*Flow_1/Flow_1,
                            PAR_1=PAR_1*Flow_1/Flow_1,
                            dC_1=dC_1*Flow_1/Flow_1,
                            dH_1=dH_1*Flow_1/Flow_1,
                            e_1=e_1*Flow_1/Flow_1,
                            VPD_1=VPD_1*Flow_1/Flow_1)

# For Chamber 3, forces NAs in the period when the tube was 
# out of the chamber (i.e. not measuring CO2 inside chamber)
Alldata<-Alldata %>% dplyr::mutate(Flow_3=case_when((Time10MIN>"2024-08-14 14:30" 
                       & Time10MIN<"2024-08-14 7:00") ~ as.numeric(NA),
                       TRUE ~ Flow_3))

# transfers "NA" from Flow_3 to other _3 variables
Alldata<-Alldata %>% dplyr::mutate(CO2_3=CO2_3*Flow_3/Flow_3,
                            H2O_3=H2O_3*Flow_3/Flow_3,
                            Tc_3=Tc_3*Flow_3/Flow_3,
                            HR_3=HR_3*Flow_3/Flow_3,
                            PAR_3=PAR_3*Flow_3/Flow_3,
                            dC_3=dC_3*Flow_3/Flow_3,
                            dH_3=dH_3*Flow_3/Flow_3,
                            e_3=e_3*Flow_3/Flow_3,
                            VPD_3=VPD_3*Flow_3/Flow_3)

#### FUNCTION Extremflow ####
# used to remove extreme values of flow
Extremeflow<-function(x){
           case_when ((x < 27 | x > 42) ~as.numeric(NA),
                      TRUE ~ x)
}

# makes NAs in extreme Flow values due to wrong measurement
# i.e. not affecting the other variables
# (e.g. sensor malfunction or high battery voltage)
Alldata<-Alldata %>% dplyr::mutate(Flow_1=Extremeflow(Flow_1),
                   Flow_2=Extremeflow(Flow_2),
                   Flow_3=Extremeflow(Flow_3),
                   Flow_4=Extremeflow(Flow_4),
                   Flow_5=Extremeflow(Flow_5),
                   Flow_6=Extremeflow(Flow_6))

# Replaces NAs in Flow with mean Flow
Alldata<-Alldata %>% dplyr::mutate(Flow_1=na_mean(Flow_1),
                      Flow_2=na_mean(Flow_2),
                      Flow_3=na_mean(Flow_3),
                      Flow_4=na_mean(Flow_4),
                      Flow_5=na_mean(Flow_5),
                      Flow_6=na_mean(Flow_6),
                      )

# Creates new variable Flow_mean
Alldata<-Alldata %>% dplyr::mutate(Flow_mean=rowMeans(dplyr::across(Flow_1:Flow_6),
                                           na.rm=TRUE))

#### export csv files ####
write.csv(DataCIRAS,"DataCIRAS.csv",row.names=FALSE)
write.csv(DataCR1000,"DataCR1000.csv",row.names=FALSE)
write.csv(DataZ6,"DataZ6.csv",row.names=FALSE)
write.csv(Alldata,"Alldata.csv",row.names=FALSE)

#### CALCULATIONS ####

##### mean flow, hobo for E, CIRAS for dilution #####
# CALCULATION Use mean Flow value dplyr::across chambers,
# Use hobo sensors for dH (for E calculation)
# Use CIRAS input for CO2 dilution correction

Alldata<-read.csv("Alldata.csv")

Alldatacalc3d<-dplyr::filter(Alldata,PAR_ext>100)
i<-1
for (i in 1:6){
  datain1<-Alldatacalc3d %>% dplyr::select(round,Time10MIN,Pr)
  # datain2<-Alldatacalc3d %>% select(
  #   # starts_with(c("H2O_","dH_","CO2_","dC_"))
  #   # & ends_with(paste(i)))
  # e_chamber is used as H2O chamber
  datain2<-Alldatacalc3d %>% dplyr::select(
    starts_with(c("e_","H2O_","dH_","CO2_","dC_")) 
    & ends_with(paste(i)))
  # recalculates dH_ as e_chamber - H2Or=(H2O_-dH)
  datain2[,3]<-datain2[,1] -(datain2[,2]-datain2[,3]) 
  # datain2[,1]="e_", datain2[,2]="H2O_", datain2[,3]="dH_"
  #
  datain3<-Alldatacalc3d %>% dplyr::select(Flow_mean)
  datain4<-Alldatacalc3d %>% dplyr::select(
    starts_with(c("Tc_","PAR_")) 
    & ends_with(paste(i)))
  datain<-cbind(datain1,datain2,datain3,datain4)
  Area<-read.csv("LEAFAREA.csv")
  Area<-Area[,c(1,1+i,7+i,13+i,19+i)]
  datain<-left_join(datain,Area[,c(1,c(3,4,5),2)],by="round")
  colnames(datain)<-c("round","Time10MIN","P","eout","edil","deltaH",
                      "CO2out","deltaC","Flow",
                      "Tair","PAR","d","LeafArea","Atotal","Sp")
  
  dataout<-CIRAS2(datain[,c(2:15)],ID=i,wind=0.2)
  Alldatacalc3d<-left_join(Alldatacalc3d,dataout,by="Time10MIN")
}

# rearrange dataout columns
Spcols<-Alldatacalc3d %>% dplyr::select(starts_with("Sp"))
Ecols<-Alldatacalc3d %>% dplyr::select(c(E1,E2,E3,E4,E5,E6))
Acols<-Alldatacalc3d %>% dplyr::select(starts_with("A"))
Tleafcols<-Alldatacalc3d %>% dplyr::select(starts_with("Tleaf"))
VPDleafcols<-Alldatacalc3d %>% dplyr::select(starts_with("VPDleaf"))
gtotalcols<-Alldatacalc3d %>% dplyr::select(starts_with("gtotal"))
gscols<-Alldatacalc3d %>% dplyr::select(starts_with("gs"))
#
Alldatacalc3d <-cbind(Alldatacalc3d[,c(1:119)],
                      Spcols,Ecols,Acols,Tleafcols,
                      VPDleafcols,gtotalcols,gscols)
#
write.csv(Alldatacalc3d,"Alldatacalc3d.csv",row.names=FALSE)

Alldatacalclean<-Alldatacalc3d[,c(1:2,77,14,15,80,119,30:32,76,120:125,33:50,70:75,107:112,126:161)]

# Starts "gathered" dataset with Sp (IDplant), Chamber=1:6, Sp=c(ile,fag,coc)
Spcols<-Alldatacalclean[,1:17] %>% gather(key="Chamber",value="IDplant",starts_with("Sp"))
Spcols<-Spcols %>% dplyr::mutate(Chamber=as.numeric(str_sub(Chamber,-1)))
Spcols<-Spcols %>% dplyr::mutate(Sp=str_sub(IDplant,1,3))

# Adds gathered data for the rest of variables
varlist<-c("Tc","HR","PAR","VPD_","WP","E","A","Tleaf","VPDleaf","gtotal","gs")

for (i in 1:11) {
  gathercols<-Alldatacalclean[,12:83] %>% dplyr::select(starts_with(varlist[i])) %>% gather()
  gathercols<-dplyr::select(gathercols,starts_with("value"))
  colnames(gathercols)<-varlist[i]
  Spcols<-cbind(Spcols,gathercols)
}

# renames VPD_ as VPDc (inside chamber)
Spcols3d<-Spcols %>% dplyr::rename(VPDc=VPD_)
Spcols3d<-data.frame(Option="3d",Spcols3d)

# Cleaned dataset
Combicols<-na.omit(Spcols3d)

# For consistency, rename "chamber" variables 
# with subscript "c"
Combicols<-dplyr::rename(Combicols,
                         HRc=HR,
                         PARc=PAR)

write.csv(Combicols,"CombicolsD.csv",row.names=FALSE)

#### 1-hour summary ####
Combicols<-read.csv("CombicolsD.csv",header=TRUE)
Combicols$Time10MIN<-as_datetime(Combicols$Time10MIN)
Combicols$Timehour<-round_date(Combicols$Time10MIN,unit="1 hour")
CombiSp<-Combicols %>%
  group_by(Option,Sp,IDplant,Chamber,Timehour) %>%
  summarise_all(mean,na.rm=TRUE)

CombiSp<-dplyr::mutate(CombiSp,
                WPgroup=factor(case_when(WP>=-100 ~"(0, -0.1 MPa)",
                                  WP<(-100) & WP>=-500 ~"[-0.1, -0.5)",
                                  WP<(-500) & WP>=-1500 ~"[-0.5 -1.5)",
                                  WP<(-1500) & WP>=-3500~"[-1.5 -3.5)")),
                VPDgroup=factor(case_when(VPDleaf>=10 & VPDleaf<20 ~"VPD 1-2 kPa",
                                          VPDleaf>=30 & VPDleaf<40 ~"VPD 3-4 kPa",
                                          VPDleaf>=50 & VPDleaf<90~"VPD 5-9 kPa")),
                WUE=A/E)


#####  Qfag #####
dataset1<-dplyr::filter(CombiSp,Sp=="fag" & Option=="3d" 
                 # & WP>=-1500 
                 & PARc>500
                 & HRc<80
                 & WUE>0.01
                 & WUE<12
                 )
#####  Qile #####
dataset2<-dplyr::filter(CombiSp,Sp=="ile"& Option=="3d"
                 # & WP>=-1500 
                 & PARc>500
                 & HRc<80
                 & WUE>0.01
                 & WUE<12
                 )
##### Qcoc #####
dataset3<-dplyr::filter(CombiSp,Sp=="coc" & Option=="3d" 
                 # & WP>=-1500 
                 & PARc>500
                 & HRc<80
                 & WUE>0.01
                 & WUE<12
)
dataset<-rbind(dataset1,dataset2,dataset3)
# add vp_air inside chamber
dataset<-dplyr::mutate(dataset,
                vpc=esat(Tc)*HRc/100)

write.csv(dataset,"dataset.csv",row.names=FALSE)

#### SUPP. GRAPHS ####


dataset<-read.csv("dataset.csv",header=TRUE)

dataset$Sp<-factor(dataset$Sp,levels=c("fag","ile","coc"))
Sp.labs <- c("Q. faginea","Q. rotundifolia","Q. coccifera")
names(Sp.labs) <- c("fag","ile","coc")

WPlinePalette <- c("#753386", "#377EB8", "#4DAF4A", "#E4B61D")
WPfillPalette <- c("#B991C4", "#91B9DC", "#A0D79E", "#F2DC82")			

WPbreaks<-c("(0, -0.1 MPa)","[-0.1, -0.5)","[-0.5 -1.5)","[-1.5 -3.5)")

# TeX("\\textit{R}$_{s}$ ($mu$mol m$^{-2}$ s$^{-1}$)"

#### VPDc vs VPDleaf ####
tiff(paste("VPDleaf_VPDair.tiff"), width=6000, height=3000,res=600,units="px",compression="lzw")

refline<-data.frame(VPDc=c(0,12))

ggplot(data=dplyr::filter(dataset),aes(x=VPDc/10,y=VPDleaf/10,color=WPgroup))+
  geom_point(aes(shape=WPgroup,fill=WPgroup),alpha=1,size=2,stroke=NA)+
  geom_line(data=refline,aes(x=VPDc,y=VPDc),color="darkgrey",lty="dashed")+
  scale_y_continuous(name=TeX("VPD$_{leaf}$ (kPa)"),breaks=c(0,3,6,9,12),limits=c(0,12))+
  scale_x_continuous(name=TeX("VPD$_{chamber}$ (kPa)"),breaks=c(0,3,6,9,12),limits=c(0,12))+
  geom_smooth(method = "lm", se = FALSE)+
  scale_shape_manual(name="SWP (MPa)",values=c(21,17,15,18),breaks=WPbreaks)+
  scale_colour_manual(name="SWP (MPa)",values=WPlinePalette,breaks=WPbreaks)+
  scale_fill_manual(name="SWP (MPa)",values=WPfillPalette,breaks=WPbreaks)+
  # geom_line(data=Extrapolate_cols,aes(x=VPDleaf/10,y=gs),
  #           colour="blue",,linetype="dashed",linewidth=1.5)+
  facet_wrap(facets=vars(Sp),labeller=labeller(Sp=Sp.labs))+
  # guides(shape="none")+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=12,face="italic"),
    axis.text.x = element_text(size = 12,color="black"),
    axis.text.y = element_text(size = 12,color="black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size=10,hjust=1),
    legend.title= element_text(size=12))
dev.off()

# create quartiles for Tc
fag<-(dplyr::filter(dataset,Sp=="fag" & WP>=-100))
ile<-(dplyr::filter(dataset,Sp=="ile" & WP>=-100))
coc<-(dplyr::filter(dataset,Sp=="coc" & WP>=-100))
fag<-dplyr::mutate(fag,
              Tc_lev=factor(case_when(Tc<=quantile(Tc,0.25)~"Q1",
                                      Tc>quantile(Tc,0.25) & Tc<=quantile(Tc,0.5)~"Q2",
                                      Tc>quantile(Tc,0.5) & Tc<=quantile(Tc,0.75)~"Q3",
                                      Tc>quantile(Tc,0.25)~"Q4")))
ile<-dplyr::mutate(ile,
            Tc_lev=factor(case_when(Tc<=quantile(Tc,0.25)~"Q1",
                                    Tc>quantile(Tc,0.25) & Tc<=quantile(Tc,0.5)~"Q2",
                                    Tc>quantile(Tc,0.5) & Tc<=quantile(Tc,0.75)~"Q3",
                                    Tc>quantile(Tc,0.25)~"Q4")))
coc<-dplyr::mutate(coc,
            Tc_lev=factor(case_when(Tc<=quantile(Tc,0.25)~"Q1",
                                    Tc>quantile(Tc,0.25) & Tc<=quantile(Tc,0.5)~"Q2",
                                    Tc>quantile(Tc,0.5) & Tc<=quantile(Tc,0.75)~"Q3",
                                    Tc>quantile(Tc,0.25)~"Q4")))

# compare curves by Tc level
VPDfag<-drm(gs ~ VPDleaf,Tc_lev,#weights=1/Incubotron$std.error,
            data=fag,fct=L.3(names=c("slope","gsmax","P50")))
VPDfag
compParm(VPDfag,"slope","-")
compParm(VPDfag,"gsmax","-")
compParm(VPDfag,"P50","-")

VPDile<-drm(gs ~ VPDleaf,Tc_lev,#weights=1/Incubotron$std.error,
            data=ile,fct=L.3(names=c("slope","gsmax","P50")))
VPDile
compParm(VPDile,"slope","-")
compParm(VPDile,"gsmax","-")
compParm(VPDile,"P50","-")

VPDcoc<-drm(gs ~ VPDleaf,Tc_lev,#weights=1/Incubotron$std.error,
            data=coc,fct=L.3(names=c("slope","gsmax","P50")))
VPDcoc
compParm(VPDcoc,"slope","-")
compParm(VPDcoc,"gsmax","-")
compParm(VPDcoc,"P50","-")

# combine dataset
dataset_tcvp<-rbind(fag,ile,coc)

# summary data by Tc level
dataTclev<-dataset_tcvp %>%
  group_by(Sp,Tc_lev) %>%
  dplyr::summarise(Tcmin=min(Tc,na.rm=TRUE),
            Tcmax=max(Tc,na.rm=TRUE),
            PARmean=mean(PARc,na.rm=TRUE),
            Prmean=mean(Pr,na.rm=TRUE),
            Tleafmean=mean(Tleaf,na.rm=TRUE),
            vpmean=mean(vpc,na.rm=TRUE),
            VPDleafmean=mean(VPDleaf,na.rm=TRUE),
            gsmean=mean(gs,na.rm=TRUE))

dataTclev<-dplyr::mutate(dataTclev,
                  Tcrange=paste(Tc_lev," (",as.character(round(Tcmin,digits=1)),
                          " ",as.character(round(Tcmax,digits=1)),")",
                          sep=""))

write.csv(dataTclev,file="Tcquartiles.csv")

summaryq<-cbind(data.frame(x=6,y=c(1325,1250,1175,1100)),
                        dataTclev[,c(1:2,11)])
summaryq$Sp<-factor(summaryq$Sp,levels=c("fag","ile","coc"))

Rangetitle<-data.frame(x=5.5,y=1400,Sp=c("fag","ile","coc"))
Rangetitle$Sp<-factor(Rangetitle$Sp,levels=c("fag","ile","coc"))

tiff(paste("gs_VPD_Tc.tiff"), width=6000, height=2500,res=600,units="px",compression="lzw")
ggplot(data=dplyr::filter(dataset_tcvp),aes(x=VPDleaf/10,y=gs))+
  geom_point(aes(color=Tc_lev),alpha=1,shape=21,size=1)+
  geom_text(data=Rangetitle,aes(x=x,y=y,label=TeX("\\textit{T}$_{c}$ range (ºC)")),size=4,show.legend=FALSE,hjust=0)+
  geom_text(data=summaryq,aes(x=x,y=y,label=Tcrange,color=Tc_lev),size=4,show.legend=FALSE,hjust=0)+
  scale_y_continuous(name=TeX("\\textit{g}$_{s}$ (mmol m$^{-2}$ s$^{-1}$)"),limits=c(0,1500))+
  scale_x_continuous(name=TeX("VPD$_{leaf}$ (kPa)"),breaks=c(0,3,6,9,12),limits=c(0,12))+
  geom_smooth(aes(color=Tc_lev),method = "drm", method.args = list(fct = L.3()), se = FALSE)+
  scale_colour_manual(name=TeX("quantile \\textit{T}$_{c}$"),values=WPlinePalette)+
  # geom_line(data=Extrapolate_cols,aes(x=VPDleaf/10,y=gs),
  #           colour="blue",,linetype="dashed",linewidth=1.5)+
  facet_wrap(facets=vars(Sp),labeller=labeller(Sp=Sp.labs))+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=12,face="italic"),
    axis.text.x = element_text(size = 12,color="black"),
    axis.text.y = element_text(size = 12,color="black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size=10,hjust=1),
    legend.title= element_text(size=12),
    legend.position="none")
dev.off()

#### gs vs VPD x WP linear ####

##### faginea #####
dat_fag <- dataset |>
  dplyr::filter(
    Sp == "fag",
    base::is.finite(gs),
    gs > 0,
    base::is.finite(VPDleaf),
    base::is.finite(WP),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    IDplant = base::factor(IDplant),
    
    # Positive soil-drought index:
    # 0 = wet; larger values = drier soil
    Dsoil = -WP,
    
    # Centre predictors to improve interpretation and reduce
    # correlation between main effects and the interaction
    VPD_c = VPDleaf -
      base::mean(VPDleaf, na.rm = TRUE),
    
    Dsoil_c = Dsoil -
      base::mean(Dsoil, na.rm = TRUE),
    
    log_gs = base::log(gs)
  ) |>
  base::droplevels()

dat_fag |>
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_plants = dplyr::n_distinct(IDplant),
    min_gs = base::min(gs),
    max_gs = base::max(gs),
    VPD_min = base::min(VPDleaf),
    VPD_max = base::max(VPDleaf),
    WP_min = base::min(WP),
    WP_max = base::max(WP)
  )

mod_VPD <- lme4::lmer(
  log_gs ~ VPD_c +
    (1 | IDplant),
  
  data = dat_fag,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_additive <- lme4::lmer(
  log_gs ~ VPD_c + Dsoil_c +
    (1 | IDplant),
  
  data = dat_fag,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_interaction <- lme4::lmer(
  log_gs ~ VPD_c * Dsoil_c +
    (1 | IDplant),
  
  data = dat_fag,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

stats::anova(
  mod_VPD,
  mod_additive,
  mod_interaction
)

stats::AIC(
  mod_VPD,
  mod_additive,
  mod_interaction
)

lme4::fixef(mod_interaction)

#DHARMA
simulationOutput <- simulateResiduals(fittedModel = mod_interaction, plot = F)
plot(simulationOutput)

##### ilex #####

dat_ile <- dataset |>
  dplyr::filter(
    Sp == "ile",
    base::is.finite(gs),
    gs > 0,
    base::is.finite(VPDleaf),
    base::is.finite(WP),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    IDplant = base::factor(IDplant),
    
    # Positive soil-drought index:
    # 0 = wet; larger values = drier soil
    Dsoil = -WP,
    
    # Centre predictors to improve interpretation and reduce
    # correlation between main effects and the interaction
    VPD_c = VPDleaf -
      base::mean(VPDleaf, na.rm = TRUE),
    
    Dsoil_c = Dsoil -
      base::mean(Dsoil, na.rm = TRUE),
    
    log_gs = base::log(gs)
  ) |>
  base::droplevels()

dat_ile |>
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_plants = dplyr::n_distinct(IDplant),
    min_gs = base::min(gs),
    max_gs = base::max(gs),
    VPD_min = base::min(VPDleaf),
    VPD_max = base::max(VPDleaf),
    WP_min = base::min(WP),
    WP_max = base::max(WP)
  )

mod_VPD <- lme4::lmer(
  log_gs ~ VPD_c +
    (1 | IDplant),
  
  data = dat_ile,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_additive <- lme4::lmer(
  log_gs ~ VPD_c + Dsoil_c +
    (1 | IDplant),
  
  data = dat_ile,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_interaction <- lme4::lmer(
  log_gs ~ VPD_c * Dsoil_c +
    (1 | IDplant),
  
  data = dat_ile,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

stats::anova(
  mod_VPD,
  mod_additive,
  mod_interaction
)

stats::AIC(
  mod_VPD,
  mod_additive,
  mod_interaction
)

lme4::fixef(mod_interaction)

#DHARMA
simulationOutput <- simulateResiduals(fittedModel = mod_interaction, plot = F)
plot(simulationOutput)

#### An vs VPD x WP linear ####

##### faginea #####
dat_fag <- dataset |>
  dplyr::filter(
    Sp == "fag",
    base::is.finite(A),
    A > 0,
    base::is.finite(VPDleaf),
    base::is.finite(WP),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    IDplant = base::factor(IDplant),
    
    # Positive soil-drought index:
    # 0 = wet; larger values = drier soil
    Dsoil = -WP,
    
    # Centre predictors to improve interpretation and reduce
    # correlation between main effects and the interaction
    VPD_c = VPDleaf -
      base::mean(VPDleaf, na.rm = TRUE),
    
    Dsoil_c = Dsoil -
      base::mean(Dsoil, na.rm = TRUE),
    
    log_A = base::log(A)
  ) |>
  base::droplevels()

dat_fag |>
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_plants = dplyr::n_distinct(IDplant),
    min_A = base::min(A),
    max_A = base::max(A),
    VPD_min = base::min(VPDleaf),
    VPD_max = base::max(VPDleaf),
    WP_min = base::min(WP),
    WP_max = base::max(WP)
  )

mod_VPD <- lme4::lmer(
  A ~ VPD_c +
    (1 | IDplant),
  
  data = dat_fag,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_additive <- lme4::lmer(
  A ~ VPD_c + Dsoil_c +
    (1 | IDplant),
  
  data = dat_fag,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_interaction <- lme4::lmer(
  A ~ VPD_c * Dsoil_c +
    (1 | IDplant),
  
  data = dat_fag,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

stats::anova(
  mod_VPD,
  mod_additive,
  mod_interaction
)

stats::AIC(
  mod_VPD,
  mod_additive,
  mod_interaction
)

lme4::fixef(mod_interaction)

#DHARMA
simulationOutput <- simulateResiduals(fittedModel = mod_interaction, plot = F)
plot(simulationOutput)

##### ilex #####

dat_ile <- dataset |>
  dplyr::filter(
    Sp == "ile",
    base::is.finite(A),
    A > 0,
    base::is.finite(VPDleaf),
    base::is.finite(WP),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    IDplant = base::factor(IDplant),
    
    # Positive soil-drought index:
    # 0 = wet; larger values = drier soil
    Dsoil = -WP,
    
    # Centre predictors to improve interpretation and reduce
    # correlation between main effects and the interaction
    VPD_c = VPDleaf -
      base::mean(VPDleaf, na.rm = TRUE),
    
    Dsoil_c = Dsoil -
      base::mean(Dsoil, na.rm = TRUE),
    
    log_A = base::log(A)
  ) |>
  base::droplevels()

dat_ile |>
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_plants = dplyr::n_distinct(IDplant),
    min_A = base::min(A),
    max_A = base::max(A),
    VPD_min = base::min(VPDleaf),
    VPD_max = base::max(VPDleaf),
    WP_min = base::min(WP),
    WP_max = base::max(WP)
  )

mod_VPD <- lme4::lmer(
  A ~ VPD_c +
    (1 | IDplant),
  
  data = dat_ile,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_additive <- lme4::lmer(
  A ~ VPD_c + Dsoil_c +
    (1 | IDplant),
  
  data = dat_ile,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

mod_interaction <- lme4::lmer(
  A ~ VPD_c * Dsoil_c +
    (1 | IDplant),
  
  data = dat_ile,
  REML = FALSE,
  
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

stats::anova(
  mod_VPD,
  mod_additive,
  mod_interaction
)

stats::AIC(
  mod_VPD,
  mod_additive,
  mod_interaction
)

lme4::fixef(mod_interaction)

#DHARMA
simulationOutput <- simulateResiduals(fittedModel = mod_interaction, plot = F)
plot(simulationOutput)


#### TEST Tleaf ####

Leafcombi<-read.csv("CombicolsD.csv",header=TRUE)
Alldata<-read.csv("Alldata.csv",header=TRUE)

Leaftemp<-dplyr::select(Alldata,c(Time10MIN,starts_with("T_tc")))
colnames(Leaftemp)
summary(Leaftemp)

i<-1
# rearrange thermocouple data by chamber (3 Tcs per chamber)
for (i in 1:3) {
  Leaftempsub<-Leaftemp[,c(1,i+seq(1,16,by=3))]
  colnames(Leaftempsub)<-c("Time10MIN",1,2,3,4,5,6)
  Leaftempsub<-gather(Leaftempsub,key="Chamber",value="T_Tc",-Time10MIN)
  Leaftempsub<-filter(Leaftempsub,T_Tc>0)
  colnames(Leaftempsub)<-c("Time10MIN","Chamber",paste("T_Tc",i,sep=""))
  Leaftempsub$Chamber<-as.numeric(Leaftempsub$Chamber)
  # remove outliers
  Leafcombi<-left_join(Leafcombi,Leaftempsub,by=c("Time10MIN","Chamber"))
}

Leafcombi<-dplyr::mutate(Leafcombi,
                         T_Tcmean=rowMeans(dplyr::across(T_Tc1:T_Tc3),na.rm=TRUE))

Leafcombi<-dplyr::mutate(Leafcombi,
                         PARlevel=case_when(PARc>=500 ~ "PAR>=500",
                                            PARc<500 ~ "PAR<500"),
                         WPgroup=factor(case_when(WP>=-100 ~"(0, -0.1 MPa)",
                                                  WP<(-100) & WP>=-500 ~"[-0.1, -0.5)",
                                                  WP<(-500) & WP>=-1500 ~"[-0.5 -1.5)",
                                                  WP<(-1500) & WP>=-3500~"[-1.5 -3.5)")),
                         VPDgroup=factor(case_when(VPDleaf>=10 & VPDleaf<20 ~"VPD 1-2 kPa",
                                                   VPDleaf>=30 & VPDleaf<40 ~"VPD 3-4 kPa",
                                                   VPDleaf>=50 & VPDleaf<90~"VPD 5-9 kPa")),
                         WUE=A/E)

faginea<-Leafcombi%>% 
  filter(Sp=="fag") %>% dplyr::mutate(Erel=E/quantile(E,0.999,na.rm=TRUE))

ilex<-Leafcombi %>% 
  filter(Sp=="ile") %>% dplyr::mutate(Erel=E/quantile(E,0.999,na.rm=TRUE))

coccifera<-Leafcombi %>% 
  filter(Sp=="coc") %>% dplyr::mutate(Erel=E/quantile(E,0.999,na.rm=TRUE))

Leafcombi<-rbind(faginea,ilex,coccifera)

Leafcombi$Sp<-factor(Leafcombi$Sp,levels=c("fag","ile","coc"))
Leafcombi$PARlevel<-factor(Leafcombi$PARlevel)

Leafcombilong<-pivot_longer(Leafcombi,c(T_Tc1:T_Tc3),
                            names_to="Thermocouple",values_to="T_Tc")

Leafcombilong<-dplyr::mutate(Leafcombilong,
                             VPDleafTc=esat(T_Tc)-esat(Tc)*HRc/100)

write.csv(Leafcombilong,"Leafcombilong.csv",row.names=FALSE)

# remove round 4 due to electrical instability
Leafcombi2<-filter(Leafcombi,round!=4)     

# remove outliers
Leafcombi3<-Leafcombi2 %>%
  group_by(Sp) %>%
  dplyr::summarise(T_Tcdev_loq=quantile(T_Tcmean-Tc,0.001,na.rm=TRUE),
                   T_Tcdev_hiq=quantile(T_Tcmean-Tc,0.999,na.rm=TRUE))

Leafcombi<-left_join(Leafcombi2,Leafcombi3,by="Sp")

Leafcombi<-filter(Leafcombi,(T_Tcmean-Tc)>T_Tcdev_loq & (T_Tcmean-Tc)<T_Tcdev_hiq)

subset1<-Leafcombi %>% filter(IDplant=="coc01") %>% dplyr::mutate(T_Tc=T_Tc3)
subset2<-Leafcombi %>% filter(IDplant=="coc01") %>% dplyr::mutate(T_Tc=T_Tc2)
subset3<-Leafcombi %>% filter(IDplant=="coc05") %>% dplyr::mutate(T_Tc=T_Tc2)
subset4<-Leafcombi %>% filter(IDplant=="coc05") %>% dplyr::mutate(T_Tc=T_Tc3)

coccifera<-rbind(subset1,subset2,subset3,subset4)

subset1<-Leafcombi %>% filter(IDplant=="fag01") %>% dplyr::mutate(T_Tc=T_Tc2)
subset2<-Leafcombi %>% filter(IDplant=="fag01") %>% dplyr::mutate(T_Tc=T_Tc3)
subset3<-Leafcombi %>% filter(IDplant=="fag03") %>% dplyr::mutate(T_Tc=T_Tc2)

faginea<-rbind(subset1,subset2,subset3)

subset1<-Leafcombi %>% filter(IDplant=="ilex06") %>% dplyr::mutate(T_Tc=T_Tc2)
subset2<-Leafcombi %>% filter(IDplant=="ilex07") %>% dplyr::mutate(T_Tc=T_Tc3)
subset3<-Leafcombi %>% filter(IDplant=="ilex09") %>% dplyr::mutate(T_Tc=T_Tc1)

ilex<-rbind(subset1,subset2,subset3)

Leaftemp<-rbind(faginea,ilex,coccifera)

Leaftemp<-dplyr::mutate(Leaftemp,
                        VPDleafTc=esat(T_Tc)-esat(Tc)*HRc/100)

# PAR.labs <- c("PAR <500", "PAR >=500")
names(PAR.labs) <- c("PAR<500", "PAR>=500")
Sp.labs <- c("Q. faginea","Q. rotundifolia","Q. coccifera")
names(Sp.labs) <- c("fag","ile","coc")

write.csv(Leaftemp,"Leaftemp_Tcfilter.csv",row.names=FALSE)

#### Graph Tair_Tchamber ####
tiff(paste("Tair_Tchamber_E.tiff"), width=4975, height=3900,res=600,units="px",compression="lzw")
ggplot(data=filter(Leafcombi,round> 0), aes(x=T_ext,y=Tc,colour=E))+
  geom_point(shape=16,alpha=1,size=0.75)+
  geom_line(aes(x=Tc,y=Tc),colour="black",lty="dashed",size=1)+
  geom_smooth(method="lm")+
  facet_grid(PARlevel~Sp,
             labeller=labeller(PARlevel=PAR.labs,Sp=Sp.labs))+
  scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),
                         breaks=c(5,10,15),limits=c(0,16),palette="Spectral",direction=1)+
  scale_y_continuous(name="Chamber temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
  scale_x_continuous(name="Ambient temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10,face="italic"),
    strip.text.y=element_text(size=10),
    axis.text.x = element_text(size = 10,color="black"),
    axis.text.y = element_text(size = 10,color="black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size=10,hjust=0.5),
    legend.title= element_text(size=12,vjust=1),
    legend.position="top")
dev.off()

#### Graph Tchamber Tleafmod ####
tiff(paste("Tchamber_Tleafmod_byPAR_E.tiff"), width=4975, height=3900,res=600,units="px",compression="lzw")
ggplot(data=filter(Leafcombi,round> 0), aes(x=Tc,y=Tleaf,colour=E))+
  geom_point(shape=16,alpha=1,size=0.75)+
  geom_line(aes(x=Tc,y=Tc),colour="black",lty="dashed",size=1)+
  geom_smooth(method="lm",se=TRUE)+
  facet_grid(PARlevel~Sp,
             labeller=labeller(PARlevel=PAR.labs,Sp=Sp.labs))+
  scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),
                         breaks=c(5,10,15),limits=c(0,16),palette="Spectral",direction=1)+
  scale_y_continuous(name="Modelled leaf temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
  scale_x_continuous(name="Chamber temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10,face="italic"),
    strip.text.y=element_text(size=10),
    axis.text.x = element_text(size = 10,color="black"),
    axis.text.y = element_text(size = 10,color="black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size=10,hjust=0.5),
    legend.title= element_text(size=12,vjust=1),
    legend.position="top")
dev.off()

#### Graph Tchamber Tleaf TC ####
tiff(paste("Tleafmod_T_Tc_byPAR_E.tiff"), width=4975, height=3900,res=600,units="px",compression="lzw")
ggplot(data=filter(Leaftemp,round> 0), aes(x=Tleaf,y=T_Tc,colour=E))+
  geom_point(data=Leafcombilong,aes(x=Tleaf,y=T_Tc),colour="gray",alpha=0.25,size=0.5)+
  geom_point(shape=16,alpha=1,size=0.75)+
  geom_line(aes(x=Tc,y=Tc),colour="black",lty="dashed",size=1)+
  geom_smooth(method="lm")+
  facet_grid(PARlevel~Sp,
             labeller=labeller(PARlevel=PAR.labs,Sp=Sp.labs))+
  # facet_grid(~Sp,
  #            labeller=labeller(Sp=Sp.labs))+
  scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),
                         breaks=c(5,10,15),,limits=c(0,16),palette="Spectral",direction=1)+
  # scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),palette="Spectral",direction=1)+
  scale_y_continuous(name="Thermocouple leaf temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
  scale_x_continuous(name="Modelled leaf temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10,face="italic"),
    strip.text.y=element_text(size=10),
    axis.text.x = element_text(size = 10,color="black"),
    axis.text.y = element_text(size = 10,color="black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size=10,hjust=0.5),
    legend.title= element_text(size=12,vjust=1),
    legend.position="top")
dev.off()


#### Linear models Tleaf ####

##### function asteriskator #####

Asteriskator <- function(pvalue=0.049){
  asterisk<-case_when(pvalue< 0.001 ~"***",
                      pvalue< 0.01 ~"**",
                      pvalue< 0.05 ~"*",
                      TRUE ~"ns")
  return(asterisk)
}

##### rearrange factor PARlevel to make easier interpretation of day values #####
Leaftemp$PARlevel<-factor(Leaftemp$PARlevel,
                          levels=c("PAR>=500","PAR<500"))

##### faginea #####
subset<-filter(Leaftemp,
               round> 0 & Sp=="fag")

model<-lm(T_Tc ~ Tleaf*PARlevel,
          data=subset)

sumod<-summary(model)
cimod<-confint(model)
RSS <- c(crossprod(model$residuals))
MSE<-RSS/length(model$residuals)
RMSE<-sqrt(MSE)

fag<-data_frame(N=length(sumod$residuals),
                Intercept=paste(round(model$coefficients[1],digits=2),
                                "(",round(cimod[1],digits=2),
                                " ",round(cimod[5],digits=2),")",
                                Asteriskator(sumod$coefficients[1,4]),
                                sep=""),
                Tleaf=paste(round(model$coefficients[2],digits=2),
                            "(",round(cimod[2],digits=2),
                            " ",round(cimod[6],digits=2),")",
                            Asteriskator(sumod$coefficients[1,4]),
                            sep=""),
                PARlevel=paste(round(model$coefficients[3],digits=2),
                               "(",round(cimod[3],digits=2),
                               " ",round(cimod[7],digits=2),")",
                               Asteriskator(sumod$coefficients[3,4]),
                               sep=""),
                TleafxPARlevel=paste(round(model$coefficients[4],digits=2),
                                     "(",round(cimod[4],digits=2),
                                     " ",round(cimod[8],digits=2),")",
                                     Asteriskator(sumod$coefficients[4,4]),
                                     sep=""),
                R2adj=round(sumod$adj.r.squared,digits=3),
                RMSE=round(RMSE,digits=2))

##### ilex #####
subset<-filter(Leaftemp,
               round> 0 & Sp=="ile")

model<-lm(T_Tc ~ Tleaf*PARlevel,
          data=subset)

sumod<-summary(model)
cimod<-confint(model)
RSS <- c(crossprod(model$residuals))
MSE<-RSS/length(model$residuals)
RMSE<-sqrt(MSE)

ile<-data_frame(N=length(sumod$residuals),
                Intercept=paste(round(model$coefficients[1],digits=2),
                                "(",round(cimod[1],digits=2),
                                " ",round(cimod[5],digits=2),")",
                                Asteriskator(sumod$coefficients[1,4]),
                                sep=""),
                Tleaf=paste(round(model$coefficients[2],digits=2),
                            "(",round(cimod[2],digits=2),
                            " ",round(cimod[6],digits=2),")",
                            Asteriskator(sumod$coefficients[1,4]),
                            sep=""),
                PARlevel=paste(round(model$coefficients[3],digits=2),
                               "(",round(cimod[3],digits=2),
                               " ",round(cimod[7],digits=2),")",
                               Asteriskator(sumod$coefficients[3,4]),
                               sep=""),
                TleafxPARlevel=paste(round(model$coefficients[4],digits=2),
                                     "(",round(cimod[4],digits=2),
                                     " ",round(cimod[8],digits=2),")",
                                     Asteriskator(sumod$coefficients[4,4]),
                                     sep=""),
                R2adj=round(sumod$adj.r.squared,digits=3),
                RMSE=round(RMSE,digits=2))

##### coccifera #####
subset<-filter(Leaftemp,
               round> 0 & Sp=="coc")

model<-lm(T_Tc ~ Tleaf*PARlevel,
          data=subset)

sumod<-summary(model)
cimod<-confint(model)
RSS <- c(crossprod(model$residuals))
MSE<-RSS/length(model$residuals)
RMSE<-sqrt(MSE)

coc<-data_frame(N=length(sumod$residuals),
                Intercept=paste(round(model$coefficients[1],digits=2),
                                "(",round(cimod[1],digits=2),
                                " ",round(cimod[5],digits=2),")",
                                Asteriskator(sumod$coefficients[1,4]),
                                sep=""),
                Tleaf=paste(round(model$coefficients[2],digits=2),
                            "(",round(cimod[2],digits=2),
                            " ",round(cimod[6],digits=2),")",
                            Asteriskator(sumod$coefficients[1,4]),
                            sep=""),
                PARlevel=paste(round(model$coefficients[3],digits=2),
                               "(",round(cimod[3],digits=2),
                               " ",round(cimod[7],digits=2),")",
                               Asteriskator(sumod$coefficients[3,4]),
                               sep=""),
                TleafxPARlevel=paste(round(model$coefficients[4],digits=2),
                                     "(",round(cimod[4],digits=2),
                                     " ",round(cimod[8],digits=2),")",
                                     Asteriskator(sumod$coefficients[4,4]),
                                     sep=""),
                R2adj=round(sumod$adj.r.squared,digits=3),
                RMSE=round(RMSE,digits=2))

predTleaf<-rbind(fag,ile,coc)

write.csv(predTleaf,"predTleaf.csv")

#### LEAF SUBSAMPLING TEST ####


# LEAF-AREA SUBSAMPLING ANALYSIS WITH PERCENTILE-BASED LEAF SELECTION

#
# Objective:
#   Estimate the number of leaves that must be measured to obtain an
#   acceptable estimate of total tree leaf area.
#
# Reference total leaf area:
#   Total number of leaves counted on the tree *
#   mean area of ALL valid measured leaves.
#
# Sampling pool:
#   Only leaves between selected percentiles are eligible for random sampling.
#   By default, leaves below the 5th percentile and above the 95th percentile
#   are excluded from subsampling.
#
# Important:
#   Excluded leaves remain part of the full reference pool. Consequently,
#   error at the maximum sampling size may remain above zero because of
#   systematic selection bias.
#
# Sampling:
#   Simple random sampling WITHOUT replacement.
#




# 0. GENERAL SETTINGS


base::set.seed(123)

# Number of repeated random samples for every sample size.
# Use 1000 while testing and 5000-10000 for final results.
B <- 1000L

# Percentile interval defining leaves eligible for sampling.
#
# 0.05 and 0.95 retain the central 90% of the measured leaves.
lower_sampling_percentile <- 0.05
upper_sampling_percentile <- 0.95

# Smallest sample size to assess
minimum_sample_size <- 1L

# Accuracy target:
# select the smallest sample size for which the 95th percentile
# of absolute relative error is <= this value.
target_error_percent <- 5

# Equivalent probability target:
# at least 95% of random samples should remain within the error threshold.
target_probability <- 0.95

# Output prefix
output_prefix <- "leaf_area_subsampling"



# 1. INPUT FILES

#
# Each file should contain the leaf-area measurements for one tree.
# Replace these names with the actual file names.
#
# Example CSV structure:
#
# LeafArea
# 4.25
# 3.87
# 5.16
#


input_files <- c(
  faginea = "leafarea_faginea.csv",
  ilex = "leafarea_ilex.csv",
  coccifera = "leafarea_coccifera.csv"
)

# Name of the leaf-area column in all three files
area_column <- "LeafArea"



# 2. TOTAL NUMBER OF LEAVES COUNTED ON EACH TREE

#
# Replace these illustrative values with the total number of leaves counted
# on each complete tree.
#
# The number of counted leaves does not have to equal the number of leaves
# whose individual areas were measured.
#


total_leaf_counts <- c(
  faginea = 1000,
  ilex = 1500,
  coccifera = 1200
)



# 3. CHECK REQUIRED PACKAGES


required_packages <- c(
  "dplyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !base::vapply(
    required_packages,
    base::requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (base::length(missing_packages) > 0L) {
  
  base::stop(
    "Install the following packages before running the script: ",
    base::paste(
      missing_packages,
      collapse = ", "
    ),
    "\nFor example:\ninstall.packages(c(",
    base::paste(
      base::shQuote(missing_packages),
      collapse = ", "
    ),
    "))"
  )
}



# 4. CHECK GENERAL SETTINGS


if (
  lower_sampling_percentile < 0 ||
  upper_sampling_percentile > 1 ||
  lower_sampling_percentile >= upper_sampling_percentile
) {
  
  base::stop(
    "Percentile limits must satisfy: ",
    "0 <= lower_sampling_percentile < ",
    "upper_sampling_percentile <= 1."
  )
}

if (
  target_probability <= 0 ||
  target_probability > 1
) {
  
  base::stop(
    "target_probability must be greater than 0 and no greater than 1."
  )
}

if (B < 100L) {
  
  base::warning(
    "B is below 100. Use at least 1000 iterations for final inference."
  )
}

if (!base::identical(
  base::names(input_files),
  base::names(total_leaf_counts)
)) {
  
  base::stop(
    "input_files and total_leaf_counts must have the same names ",
    "and occur in the same order."
  )
}



# 5. READ THE THREE DATASETS


missing_files <- input_files[
  !base::file.exists(input_files)
]

if (base::length(missing_files) > 0L) {
  
  base::stop(
    "The following input files were not found:\n",
    base::paste(
      missing_files,
      collapse = "\n"
    ),
    "\nWorking directory: ",
    base::getwd()
  )
}

leaf_datasets <- base::lapply(
  input_files,
  function(current_file) {
    
    utils::read.csv(
      file = current_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      na.strings = c(
        "",
        "NA",
        "NaN"
      )
    )
  }
)



# 6. FUNCTION FOR ONE TREE


subsample_leaf_area <- function(
    data,
    species_name,
    total_leaf_count,
    area_column = "LeafArea",
    B = 5000L,
    minimum_sample_size = 1L,
    lower_percentile = 0.05,
    upper_percentile = 0.95,
    seed = 123L
) {
  
  base::set.seed(seed)
  

  # 6.1 Check inputs

  
  if (!area_column %in% base::names(data)) {
    
    base::stop(
      "Column '",
      area_column,
      "' was not found for ",
      species_name,
      ". Available columns: ",
      base::paste(
        base::names(data),
        collapse = ", "
      )
    )
  }
  
  if (
    !base::is.numeric(total_leaf_count) ||
    base::length(total_leaf_count) != 1L ||
    !base::is.finite(total_leaf_count) ||
    total_leaf_count <= 0
  ) {
    
    base::stop(
      "total_leaf_count must be a single positive number for ",
      species_name,
      "."
    )
  }
  

  # 6.2 Complete measured leaf pool

  
  leaf_area_all <- base::as.numeric(
    data[[area_column]]
  )
  
  leaf_area_all <- leaf_area_all[
    base::is.finite(leaf_area_all) &
      leaf_area_all > 0
  ]
  
  N_full_pool <- base::length(
    leaf_area_all
  )
  
  if (N_full_pool < 2L) {
    
    base::stop(
      "Fewer than two valid positive leaf-area observations remain for ",
      species_name,
      "."
    )
  }
  

  # 6.3 Percentile thresholds

  
  percentile_limits <- stats::quantile(
    leaf_area_all,
    probs = c(
      lower_percentile,
      upper_percentile
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  lower_area_limit <- percentile_limits[1]
  upper_area_limit <- percentile_limits[2]
  

  # 6.4 Restricted pool eligible for sampling

  
  sampling_eligible <- (
    leaf_area_all >= lower_area_limit &
      leaf_area_all <= upper_area_limit
  )
  
  leaf_area_sampling <- leaf_area_all[
    sampling_eligible
  ]
  
  N_sampling_pool <- base::length(
    leaf_area_sampling
  )
  
  if (N_sampling_pool < 2L) {
    
    base::stop(
      "Fewer than two leaves remain in the sampling pool for ",
      species_name,
      ". Change the percentile thresholds."
    )
  }
  
  if (minimum_sample_size > N_sampling_pool) {
    
    base::stop(
      "minimum_sample_size exceeds the eligible sampling-pool size for ",
      species_name,
      "."
    )
  }
  

  # 6.5 Reference values calculated from ALL leaves

  
  reference_mean <- base::mean(
    leaf_area_all
  )
  
  reference_median <- stats::median(
    leaf_area_all
  )
  
  reference_sd <- stats::sd(
    leaf_area_all
  )
  
  reference_CV_percent <- 100 *
    reference_sd /
    reference_mean
  
  reference_total_area <-
    total_leaf_count *
    reference_mean
  

  # 6.6 Characteristics of the percentile-restricted pool

  
  sampling_pool_mean <- base::mean(
    leaf_area_sampling
  )
  
  sampling_pool_median <- stats::median(
    leaf_area_sampling
  )
  
  sampling_pool_sd <- stats::sd(
    leaf_area_sampling
  )
  
  sampling_pool_CV_percent <- 100 *
    sampling_pool_sd /
    sampling_pool_mean
  
  # Bias remaining if every eligible leaf were measured
  selection_bias_area <-
    sampling_pool_mean -
    reference_mean
  
  selection_bias_percent <- 100 *
    selection_bias_area /
    reference_mean
  
  maximum_sampling_total_estimate <-
    total_leaf_count *
    sampling_pool_mean
  
  selection_bias_total_area <-
    maximum_sampling_total_estimate -
    reference_total_area
  
  excluded_leaf_count <-
    N_full_pool -
    N_sampling_pool
  
  excluded_leaf_percent <- 100 *
    excluded_leaf_count /
    N_full_pool
  

  # 6.7 Descriptive output

  
  base::cat(
    "\n============================================================\n",
    "LEAF-AREA SAMPLING POOLS: ",
    species_name,
    "\n============================================================\n",
    "Full measured pool: ",
    N_full_pool,
    " leaves\n",
    "Eligible sampling pool: ",
    N_sampling_pool,
    " leaves\n",
    "Excluded leaves: ",
    excluded_leaf_count,
    " (",
    base::round(
      excluded_leaf_percent,
      2
    ),
    "%)\n",
    "Lower percentile: ",
    lower_percentile,
    "\n",
    "Upper percentile: ",
    upper_percentile,
    "\n",
    "Lower leaf-area threshold: ",
    base::round(
      lower_area_limit,
      4
    ),
    "\n",
    "Upper leaf-area threshold: ",
    base::round(
      upper_area_limit,
      4
    ),
    "\n",
    "Mean area, full pool: ",
    base::round(
      reference_mean,
      4
    ),
    "\n",
    "Mean area, sampling pool: ",
    base::round(
      sampling_pool_mean,
      4
    ),
    "\n",
    "Maximum-sample selection bias: ",
    base::round(
      selection_bias_percent,
      3
    ),
    "%\n",
    sep = ""
  )
  

  # 6.8 Sample sizes to assess

  
  sample_sizes <- base::seq.int(
    from = minimum_sample_size,
    to = N_sampling_pool
  )
  
  summary_list <- base::vector(
    mode = "list",
    length = base::length(
      sample_sizes
    )
  )
  

  # 6.9 Repeated finite-pool sampling without replacement

  
  for (
    ii in base::seq_along(
      sample_sizes
    )
  ) {
    
    n_sample <- sample_sizes[ii]
    
    if (n_sample == N_sampling_pool) {
      
      # Measuring all eligible leaves gives the exact restricted-pool mean.
      sample_means <- base::rep(
        sampling_pool_mean,
        B
      )
      
    } else {
      
      sample_means <- base::vapply(
        X = base::seq_len(B),
        FUN = function(bb) {
          
          base::mean(
            base::sample(
              x = leaf_area_sampling,
              size = n_sample,
              replace = FALSE
            )
          )
        },
        FUN.VALUE = numeric(1)
      )
    }
    
    # Estimate complete-tree leaf area from the sampled mean
    estimated_total_areas <-
      total_leaf_count *
      sample_means
    
    signed_mean_errors <-
      sample_means -
      reference_mean
    
    relative_errors_percent <- 100 *
      signed_mean_errors /
      reference_mean
    
    absolute_errors_percent <-
      base::abs(
        relative_errors_percent
      )
    
    total_area_errors <-
      estimated_total_areas -
      reference_total_area
    
    # -------------------------------------------------------------------------
    # Summarise results for this sample size
    # -------------------------------------------------------------------------
    
    summary_list[[ii]] <- base::data.frame(
      Species = species_name,
      Sample_size = n_sample,
      
      Full_pool_size =
        N_full_pool,
      
      Sampling_pool_size =
        N_sampling_pool,
      
      Sampling_fraction_percent =
        100 *
        n_sample /
        N_sampling_pool,
      
      Mean_sample_area =
        base::mean(
          sample_means
        ),
      
      Reference_mean_area =
        reference_mean,
      
      Sampling_pool_mean =
        sampling_pool_mean,
      
      Mean_bias_area =
        base::mean(
          signed_mean_errors
        ),
      
      Relative_bias_percent =
        base::mean(
          relative_errors_percent
        ),
      
      Mean_absolute_error_percent =
        base::mean(
          absolute_errors_percent
        ),
      
      Median_absolute_error_percent =
        stats::median(
          absolute_errors_percent
        ),
      
      P75_absolute_error_percent =
        stats::quantile(
          absolute_errors_percent,
          probs = 0.75,
          names = FALSE
        ),
      
      P90_absolute_error_percent =
        stats::quantile(
          absolute_errors_percent,
          probs = 0.90,
          names = FALSE
        ),
      
      P95_absolute_error_percent =
        stats::quantile(
          absolute_errors_percent,
          probs = 0.95,
          names = FALSE
        ),
      
      P99_absolute_error_percent =
        stats::quantile(
          absolute_errors_percent,
          probs = 0.99,
          names = FALSE
        ),
      
      RMSE_percent =
        base::sqrt(
          base::mean(
            relative_errors_percent^2
          )
        ),
      
      Probability_within_2.5_percent =
        base::mean(
          absolute_errors_percent <= 2.5
        ),
      
      Probability_within_5_percent =
        base::mean(
          absolute_errors_percent <= 5
        ),
      
      Probability_within_7.5_percent =
        base::mean(
          absolute_errors_percent <= 7.5
        ),
      
      Probability_within_10_percent =
        base::mean(
          absolute_errors_percent <= 10
        ),
      
      Lower_95_mean_area =
        stats::quantile(
          sample_means,
          probs = 0.025,
          names = FALSE
        ),
      
      Upper_95_mean_area =
        stats::quantile(
          sample_means,
          probs = 0.975,
          names = FALSE
        ),
      
      Mean_estimated_total_area =
        base::mean(
          estimated_total_areas
        ),
      
      Reference_total_area =
        reference_total_area,
      
      Mean_total_area_error =
        base::mean(
          total_area_errors
        ),
      
      Lower_95_total_area =
        stats::quantile(
          estimated_total_areas,
          probs = 0.025,
          names = FALSE
        ),
      
      Upper_95_total_area =
        stats::quantile(
          estimated_total_areas,
          probs = 0.975,
          names = FALSE
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  error_summary <- dplyr::bind_rows(
    summary_list
  )
  

  # 6.10 Reference and pool information

  
  reference_summary <- base::data.frame(
    Species = species_name,
    
    Full_leaf_pool =
      N_full_pool,
    
    Sampling_leaf_pool =
      N_sampling_pool,
    
    Excluded_leaf_count =
      excluded_leaf_count,
    
    Excluded_leaf_percent =
      excluded_leaf_percent,
    
    Lower_percentile =
      lower_percentile,
    
    Upper_percentile =
      upper_percentile,
    
    Lower_area_limit =
      lower_area_limit,
    
    Upper_area_limit =
      upper_area_limit,
    
    Total_leaves_counted =
      total_leaf_count,
    
    Reference_mean_all_leaves =
      reference_mean,
    
    Reference_median_all_leaves =
      reference_median,
    
    Reference_SD_all_leaves =
      reference_sd,
    
    Reference_CV_percent =
      reference_CV_percent,
    
    Sampling_pool_mean =
      sampling_pool_mean,
    
    Sampling_pool_median =
      sampling_pool_median,
    
    Sampling_pool_SD =
      sampling_pool_sd,
    
    Sampling_pool_CV_percent =
      sampling_pool_CV_percent,
    
    Selection_bias_area =
      selection_bias_area,
    
    Selection_bias_percent =
      selection_bias_percent,
    
    Reference_total_area =
      reference_total_area,
    
    Maximum_sample_total_area_estimate =
      maximum_sampling_total_estimate,
    
    Selection_bias_total_area =
      selection_bias_total_area,
    
    stringsAsFactors = FALSE
  )
  

  # 6.11 Record every measured leaf and sampling eligibility

  
  leaf_classification <- base::data.frame(
    Species = species_name,
    Leaf_index =
      base::seq_along(
        leaf_area_all
      ),
    
    Leaf_area =
      leaf_area_all,
    
    Lower_area_limit =
      lower_area_limit,
    
    Upper_area_limit =
      upper_area_limit,
    
    Eligible_for_sampling =
      sampling_eligible,
    
    Selection_class =
      base::ifelse(
        leaf_area_all < lower_area_limit,
        "Below lower percentile",
        base::ifelse(
          leaf_area_all > upper_area_limit,
          "Above upper percentile",
          "Eligible sampling range"
        )
      ),
    
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      summary =
        error_summary,
      
      reference =
        reference_summary,
      
      leaf_classification =
        leaf_classification,
      
      full_pool =
        leaf_area_all,
      
      sampling_pool =
        leaf_area_sampling
    )
  )
}



# 7. RUN ANALYSIS FOR ALL THREE TREES


analysis_results <- base::lapply(
  base::names(
    leaf_datasets
  ),
  function(species_name) {
    
    subsample_leaf_area(
      data =
        leaf_datasets[[species_name]],
      
      species_name =
        species_name,
      
      total_leaf_count =
        total_leaf_counts[[species_name]],
      
      area_column =
        area_column,
      
      B =
        B,
      
      minimum_sample_size =
        minimum_sample_size,
      
      lower_percentile =
        lower_sampling_percentile,
      
      upper_percentile =
        upper_sampling_percentile,
      
      seed =
        123L
    )
  }
)

base::names(analysis_results) <-
  base::names(
    leaf_datasets
  )



# 8. COMBINE RESULTS


reference_table <- dplyr::bind_rows(
  base::lapply(
    analysis_results,
    function(x) {
      x$reference
    }
  )
)

error_summary <- dplyr::bind_rows(
  base::lapply(
    analysis_results,
    function(x) {
      x$summary
    }
  )
)

leaf_classification <- dplyr::bind_rows(
  base::lapply(
    analysis_results,
    function(x) {
      x$leaf_classification
    }
  )
)

base::cat(
  "\n============================================================\n",
  "REFERENCE AND SAMPLING-POOL VALUES\n",
  "============================================================\n",
  sep = ""
)

base::print(
  reference_table
)



# 9. FIND RECOMMENDED SAMPLE SIZE

#
# Criterion:
#   Smallest n for which at least 95% of samples have an absolute error
#   no greater than target_error_percent.
#
# This is equivalent to:
#   P95 absolute error <= target_error_percent.
#


recommended_sample_size <-
  error_summary |>
  dplyr::group_by(
    Species
  ) |>
  dplyr::filter(
    P95_absolute_error_percent <=
      target_error_percent
  ) |>
  dplyr::slice_min(
    order_by = Sample_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    Species,
    
    Recommended_sample_size =
      Sample_size,
    
    Sampling_fraction_percent,
    
    P95_absolute_error_percent,
    
    Mean_absolute_error_percent,
    
    Relative_bias_percent,
    
    Probability_within_5_percent,
    
    Reference_total_area,
    
    Mean_estimated_total_area
  )



# 10. IDENTIFY TREES FOR WHICH THE TARGET CANNOT BE ACHIEVED


species_with_recommendation <-
  base::unique(
    recommended_sample_size$Species
  )

species_without_recommendation <-
  base::setdiff(
    base::names(
      leaf_datasets
    ),
    species_with_recommendation
  )

if (
  base::length(
    species_without_recommendation
  ) > 0L
) {
  
  base::warning(
    "No sample size met the ",
    target_error_percent,
    "% error target for: ",
    base::paste(
      species_without_recommendation,
      collapse = ", "
    ),
    ". The systematic bias caused by percentile exclusion may exceed ",
    "the selected accuracy threshold."
  )
}

base::cat(
  "\n============================================================\n",
  "RECOMMENDED SAMPLE SIZE\n",
  "============================================================\n",
  "Criterion: at least ",
  target_probability * 100,
  "% of random samples within ±",
  target_error_percent,
  "% of the complete-pool reference.\n",
  sep = ""
)

base::print(
  recommended_sample_size
)



# 11. RECOMMENDATIONS FOR SEVERAL ERROR THRESHOLDS


error_thresholds <- c(
  2.5,
  5,
  7.5,
  10,
  15
)

recommendation_table <- dplyr::bind_rows(
  base::lapply(
    error_thresholds,
    function(current_threshold) {
      
      current_result <-
        error_summary |>
        dplyr::group_by(
          Species
        ) |>
        dplyr::filter(
          P95_absolute_error_percent <=
            current_threshold
        ) |>
        dplyr::slice_min(
          order_by = Sample_size,
          n = 1,
          with_ties = FALSE
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          Error_threshold_percent =
            current_threshold
        ) |>
        dplyr::select(
          Species,
          
          Error_threshold_percent,
          
          Recommended_sample_size =
            Sample_size,
          
          Sampling_fraction_percent,
          
          P95_absolute_error_percent,
          
          Relative_bias_percent
        )
      
      current_result
    }
  )
)

base::cat(
  "\nRecommended sample sizes by error threshold:\n"
)

base::print(
  recommendation_table
)



# 12. PLOT ERROR AGAINST SAMPLE SIZE


error_plot <- ggplot2::ggplot(
  data = error_summary,
  mapping = ggplot2::aes(
    x = Sample_size,
    colour = Species
  )
) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      y = Median_absolute_error_percent
    ),
    linewidth = 0.8
  ) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      y = P95_absolute_error_percent
    ),
    linewidth = 1,
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept =
      target_error_percent,
    linetype = "dotted",
    colour = "black"
  ) +
  ggplot2::labs(
    x = "Number of leaves measured",
    
    y = paste0(
      "Absolute error relative to all leaves (%)"
    ),
    
    colour = "Species",
    
    caption = base::paste0(
      "Solid lines: median absolute error; ",
      "dashed lines: 95th percentile of absolute error; ",
      "horizontal line: ",
      target_error_percent,
      "% target."
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    
    legend.position =
      "top"
  )

base::print(
  error_plot
)



# 13. PLOT PROBABILITY OF MEETING THE ERROR TARGET


probability_plot <- ggplot2::ggplot(
  data = error_summary,
  mapping = ggplot2::aes(
    x = Sample_size,
    y = Probability_within_5_percent,
    colour = Species
  )
) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_hline(
    yintercept =
      target_probability,
    
    linetype =
      "dashed",
    
    colour =
      "black"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    
    breaks = base::seq(
      0,
      1,
      by = 0.1
    ),
    
    labels = function(x) {
      base::paste0(
        base::round(
          100 * x
        ),
        "%"
      )
    }
  ) +
  ggplot2::labs(
    x = "Number of leaves measured",
    
    y = paste0(
      "Probability of absolute error ≤ ",
      target_error_percent,
      "%"
    ),
    
    colour = "Species"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    
    legend.position =
      "top"
  )

base::print(
  probability_plot
)



# 14. PLOT FULL AND RESTRICTED LEAF-AREA DISTRIBUTIONS


distribution_plot <- ggplot2::ggplot(
  data = leaf_classification,
  mapping = ggplot2::aes(
    x = Leaf_area,
    fill = Selection_class
  )
) +
  ggplot2::geom_histogram(
    bins = 30,
    colour = "white",
    alpha = 0.8
  ) +
  ggplot2::facet_wrap(
    facets = ggplot2::vars(
      Species
    ),
    
    scales = "free"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Below lower percentile" =
        "#0072B2",
      
      "Eligible sampling range" =
        "#999999",
      
      "Above upper percentile" =
        "#D55E00"
    )
  ) +
  ggplot2::labs(
    x = "Individual leaf area",
    y = "Number of leaves",
    fill = "Sampling eligibility"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    
    legend.position =
      "top"
  )

base::print(
  distribution_plot
)



# 15. EXPORT TABLES


utils::write.csv(
  reference_table,
  
  file = base::paste0(
    output_prefix,
    "_reference_and_selection_bias.csv"
  ),
  
  row.names = FALSE
)

utils::write.csv(
  error_summary,
  
  file = base::paste0(
    output_prefix,
    "_error_by_sample_size.csv"
  ),
  
  row.names = FALSE
)

utils::write.csv(
  recommended_sample_size,
  
  file = base::paste0(
    output_prefix,
    "_recommended_sample_size.csv"
  ),
  
  row.names = FALSE
)

utils::write.csv(
  recommendation_table,
  
  file = base::paste0(
    output_prefix,
    "_sample_sizes_by_error_threshold.csv"
  ),
  
  row.names = FALSE
)

utils::write.csv(
  leaf_classification,
  
  file = base::paste0(
    output_prefix,
    "_leaf_classification.csv"
  ),
  
  row.names = FALSE
)



# 16. EXPORT FIGURES


ggplot2::ggsave(
  filename = base::paste0(
    output_prefix,
    "_error_curves.tiff"
  ),
  
  plot =
    error_plot,
  
  device =
    "tiff",
  
  width =
    180,
  
  height =
    140,
  
  units =
    "mm",
  
  dpi =
    600,
  
  compression =
    "lzw"
)

ggplot2::ggsave(
  filename = base::paste0(
    output_prefix,
    "_target_probability.tiff"
  ),
  
  plot =
    probability_plot,
  
  device =
    "tiff",
  
  width =
    180,
  
  height =
    140,
  
  units =
    "mm",
  
  dpi =
    600,
  
  compression =
    "lzw"
)

ggplot2::ggsave(
  filename = base::paste0(
    output_prefix,
    "_leaf_area_distributions.tiff"
  ),
  
  plot =
    distribution_plot,
  
  device =
    "tiff",
  
  width =
    200,
  
  height =
    140,
  
  units =
    "mm",
  
  dpi =
    600,
  
  compression =
    "lzw"
)



# 17. SAVE COMPLETE R OBJECTS


base::saveRDS(
  object =
    analysis_results,
  
  file = base::paste0(
    output_prefix,
    "_complete_results.rds"
  )
)



# 18. SAVE SESSION INFORMATION


utils::capture.output(
  utils::sessionInfo(),
  
  file = base::paste0(
    output_prefix,
    "_sessionInfo.txt"
  )
)



# 19. FINAL SUMMARY


base::cat(
  "\n============================================================\n",
  "FINAL SUMMARY\n",
  "============================================================\n",
  sep = ""
)

base::cat(
  "\nReference values and systematic selection bias:\n"
)

base::print(
  dplyr::select(
    reference_table,
    
    Species,
    
    Full_leaf_pool,
    
    Sampling_leaf_pool,
    
    Excluded_leaf_percent,
    
    Reference_mean_all_leaves,
    
    Sampling_pool_mean,
    
    Selection_bias_percent,
    
    Reference_total_area
  )
)

base::cat(
  "\nRecommended sample sizes:\n"
)

base::print(
  recommended_sample_size
)

base::cat(
  "\nAnalysis completed.\n"
)