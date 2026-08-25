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
  output<-filter(Dataset,gtotal>0 & gs<1500 & gs>0 & Flow > 20 & Flow <50)
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
  output<-filter(Dataset,gtotal>0 & gs<1500 & gs>0 & Flow > 20 & Flow <50)
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

Alldatacalc3d<-filter(Alldata,PAR_ext>100)
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
dataset1<-filter(CombiSp,Sp=="fag" & Option=="3d" 
                 # & WP>=-1500 
                 & PARc>500
                 & HRc<80
                 & WUE>0.01
                 & WUE<12
                 )
#####  Qile #####
dataset2<-filter(CombiSp,Sp=="ile"& Option=="3d"
                 # & WP>=-1500 
                 & PARc>500
                 & HRc<80
                 & WUE>0.01
                 & WUE<12
                 )
##### Qcoc #####
dataset3<-filter(CombiSp,Sp=="coc" & Option=="3d" 
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

ggplot(data=filter(dataset),aes(x=VPDc/10,y=VPDleaf/10,color=WPgroup))+
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

# create quartiles for Tc and vpc
fag<-(filter(dataset,Sp=="fag" & WP>=-100))
ile<-(filter(dataset,Sp=="ile" & WP>=-100))
coc<-(filter(dataset,Sp=="coc" & WP>=-100))
fag<-dplyr::mutate(fag,
              Tc_lev=factor(case_when(Tc<=quantile(Tc,0.25)~"Q1",
                                      Tc>quantile(Tc,0.25) & Tc<=quantile(Tc,0.5)~"Q2",
                                      Tc>quantile(Tc,0.5) & Tc<=quantile(Tc,0.75)~"Q3",
                                      Tc>quantile(Tc,0.25)~"Q4")),
              vp_lev=factor(case_when(vpc<=quantile(vpc,0.25)~"Q1",
                                      vpc>quantile(vpc,0.25) & vpc<=quantile(vpc,0.5)~"Q2",
                                      vpc>quantile(vpc,0.5) & vpc<=quantile(vpc,0.75)~"Q3",
                                      vpc>quantile(vpc,0.25)~"Q4")))
ile<-dplyr::mutate(ile,
            Tc_lev=factor(case_when(Tc<=quantile(Tc,0.25)~"Q1",
                                    Tc>quantile(Tc,0.25) & Tc<=quantile(Tc,0.5)~"Q2",
                                    Tc>quantile(Tc,0.5) & Tc<=quantile(Tc,0.75)~"Q3",
                                    Tc>quantile(Tc,0.25)~"Q4")),
            vp_lev=factor(case_when(vpc<=quantile(vpc,0.25)~"Q1",
                                    vpc>quantile(vpc,0.25) & vpc<=quantile(vpc,0.5)~"Q2",
                                    vpc>quantile(vpc,0.5) & vpc<=quantile(vpc,0.75)~"Q3",
                                    vpc>quantile(vpc,0.25)~"Q4")))
coc<-dplyr::mutate(coc,
            Tc_lev=factor(case_when(Tc<=quantile(Tc,0.25)~"Q1",
                                    Tc>quantile(Tc,0.25) & Tc<=quantile(Tc,0.5)~"Q2",
                                    Tc>quantile(Tc,0.5) & Tc<=quantile(Tc,0.75)~"Q3",
                                    Tc>quantile(Tc,0.25)~"Q4")),
             vp_lev=factor(case_when(vpc<=quantile(vpc,0.25)~"Q1",
                                    vpc>quantile(vpc,0.25) & vpc<=quantile(vpc,0.5)~"Q2",
                                    vpc>quantile(vpc,0.5) & vpc<=quantile(vpc,0.75)~"Q3",
                                    vpc>quantile(vpc,0.25)~"Q4")))
dataset_tcvp<-rbind(fag,ile,coc)

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

datavplev<-dataset_tcvp %>%
  group_by(Sp,vp_lev) %>%
  dplyr::summarise(vpmin=min(vpc,na.rm=TRUE),
                   vpmax=max(vpc,na.rm=TRUE),
                   PARmean=mean(PARc,na.rm=TRUE),
                   Prmean=mean(Pr,na.rm=TRUE),
                   Tleafmean=mean(Tleaf,na.rm=TRUE),
                   Tcmean=mean(Tc,na.rm=TRUE),
                   VPDleafmean=mean(VPDleaf,na.rm=TRUE),
                   gsmean=mean(gs,na.rm=TRUE),
                   vpminchar=as.character(round(vpmin/10,digits=1)),
                   vpmaxchar=as.character(round(vpmax/10,digits=1)),
                   vpminchar=case_when(
                     str_length(vpminchar)==1 ~ paste(vpminchar,".0",sep=""),
                     TRUE ~ vpminchar),
                   vpmaxchar=case_when(
                     str_length(vpmaxchar)==1 ~ paste(vpmaxchar,".0",sep=""),
                     TRUE ~ vpmaxchar))

datavplev<-dplyr::mutate(datavplev,
                  vprange_kPa=paste(vp_lev," (",vpminchar,
                                    " ",vpmaxchar,")",sep=""))

datavplev<-dplyr::select(datavplev,!ends_with("char")) 

write.csv(datavplev,file="vpquartiles.csv")

summaryq<-cbind(data.frame(x=6,y=c(1325,1250,1175,1100)),
                        dataTclev[,c(1:2,11)],datavplev[,c(2,11)])
summaryq$Sp<-factor(summaryq$Sp,levels=c("fag","ile","coc"))

Rangetitle<-data.frame(x=5.5,y=1400,Sp=c("fag","ile","coc"))
Rangetitle$Sp<-factor(Rangetitle$Sp,levels=c("fag","ile","coc"))

tiff(paste("gs_VPD_Tc_vpc.tiff"), width=6000, height=5000,res=600,units="px",compression="lzw")
pA<-ggplot(data=filter(dataset_tcvp),aes(x=VPDleaf/10,y=gs))+
  geom_point(aes(color=Tc_lev),alpha=1,shape=21,size=1)+
  # geom_text(data=Rangetitle,aes(x=x,y=y,label=TeX("\\textit{T}$_{c}$ range (ºC)")),size=4,show.legend=FALSE,hjust=0)+
  # geom_text(data=summaryq,aes(x=x,y=y,label=Tcrange,color=Tc_lev),size=4,show.legend=FALSE,hjust=0)+
  # scale_y_continuous(name=TeX("\\textit{g}$_{s}$ (mmol m$^{-2}$ s$^{-1}$)"),limits=c(0,1500))+
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
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12,color="black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size=10,hjust=1),
    legend.title= element_text(size=12),
    legend.position="none")

pB<-ggplot(data=filter(dataset_tcvp),aes(x=VPDleaf/10,y=gs))+
  geom_point(alpha=1,shape=21,size=1)+
  geom_point(aes(color=vp_lev),alpha=1,shape=21,size=1)+
  geom_text(data=Rangetitle,aes(x=x,y=y,label=TeX("\\textit{vp}$_{c}$ range (kPa)")),size=4,show.legend=FALSE,hjust=0)+
  geom_text(data=summaryq,aes(x=x,y=y,label=vprange_kPa,color=vp_lev),size=4,show.legend=FALSE,hjust=0)+
  # geom_text(data=labels[1:2,],aes(x=x,y=y,label=label),size=8)+
  scale_y_continuous(name=TeX("\\textit{g}$_{s}$ (mmol m$^{-2}$ s$^{-1}$)"),limits=c(0,1500))+
  scale_x_continuous(name=TeX("VPD$_{leaf}$ (kPa)"),breaks=c(0,3,6,9,12),limits=c(0,12))+
  geom_smooth(aes(color=vp_lev),method = "drm", method.args = list(fct = L.3()), se = FALSE)+
  scale_colour_manual(name=TeX("quantile \\textit{vp}$_{c}$"),values=WPlinePalette)+
  # geom_line(data=Extrapolate_cols,aes(x=VPDleaf/10,y=gs),
  #           colour="blue",,linetype="dashed",linewidth=1.5)+
  facet_wrap(facets=vars(Sp),labeller=labeller(Sp=Sp.labs))+
  theme_bw()+
  theme(
    panel.grid.minor=element_blank(),
    strip.text.x=element_blank(),
    axis.text.x = element_text(size = 12,color="black"),
    axis.text.y = element_text(size = 12,color="black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size=10,hjust=1),
    legend.title= element_text(size=12),
    legend.position="none")

cowplot::plot_grid(pA, 
                   pB,ncol=1,
                                      align='v'
                   ,
                   rel_heights=c(0.32,
                                 0.35)
                   )
dev.off()

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

Leafcombi<-dplyr::mutate(Leafcombi,T_Tcmean=rowMeans(dplyr::across(T_Tc1:T_Tc3),na.rm=TRUE))

write.csv(Leafcombi,"Leafcombi.csv",row.names=FALSE)

Leafcombi<-read.csv("Leafcombi.csv",header=TRUE)

subset1<-Leafcombi %>% filter(IDplant=="coc01") %>% dplyr::mutate(T_Tc=T_Tc3)
subset2<-Leafcombi %>% filter(IDplant=="coc01") %>% dplyr::mutate(T_Tc=T_Tc2)
subset3<-Leafcombi %>% filter(IDplant=="coc03") %>% dplyr::mutate(T_Tc=T_Tc2)
subset4<-Leafcombi %>% filter(IDplant=="coc05") %>% dplyr::mutate(T_Tc=T_Tc2)
subset5<-Leafcombi %>% filter(IDplant=="coc05") %>% dplyr::mutate(T_Tc=T_Tc3)

coccifera<-rbind(subset1,subset2,subset3,subset4,subset5)

subset1<-Leafcombi %>% filter(IDplant=="fag01") %>% dplyr::mutate(T_Tc=T_Tc2)
subset2<-Leafcombi %>% filter(IDplant=="fag01") %>% dplyr::mutate(T_Tc=T_Tc3)
subset3<-Leafcombi %>% filter(IDplant=="fag03") %>% dplyr::mutate(T_Tc=T_Tc2)
subset4<-Leafcombi %>% filter(IDplant=="fag07") %>% dplyr::mutate(T_Tc=T_Tc3)

faginea<-rbind(subset1,subset2,subset3,subset4)

subset1<-Leafcombi %>% filter(IDplant=="ilex06") %>% dplyr::mutate(T_Tc=T_Tc2)
subset2<-Leafcombi %>% filter(IDplant=="ilex07") %>% dplyr::mutate(T_Tc=T_Tc1)
subset3<-Leafcombi %>% filter(IDplant=="ilex07") %>% dplyr::mutate(T_Tc=T_Tc2)
subset4<-Leafcombi %>% filter(IDplant=="ilex07") %>% dplyr::mutate(T_Tc=T_Tc3)
subset5<-Leafcombi %>% filter(IDplant=="ilex09") %>% dplyr::mutate(T_Tc=T_Tc1)

ilex<-rbind(subset1,subset2,subset3,subset4,subset5)

Leaftemp<-rbind(faginea,ilex,coccifera)

Leaftemp<-dplyr::mutate(Leaftemp,
                 PARlevel=case_when(PARc>=500 ~ "PAR>=500",
                                  PARc<500 ~ "PAR<500"))

Leaftemp$Sp<-factor(Leaftemp$Sp)
levels(Leaftemp$Sp)<-c("fag","ile","coc")
PAR.labs <- c("PAR <500", "PAR >=500")
names(PAR.labs) <- c("PAR<500", "PAR>=500")
Sp.labs <- c("Q. faginea","Q. rotundifolia","Q. coccifera")
names(Sp.labs) <- c("fag","ile","coc")

#### Graph Tair_Tchamber ####
tiff(paste("Tair_Tchamber.tiff"), width=4975, height=3900,res=600,units="px",compression="lzw")
ggplot(data=filter(Leaftemp,round> 0), aes(x=T_ext,y=Tc,colour=E))+
  geom_point(shape=16,alpha=1,size=0.5)+
  geom_line(aes(x=Tc,y=Tc),colour="black",lty="dashed",size=1)+
  geom_smooth(method="lm")+
  facet_grid(PARlevel~Sp,
             labeller=labeller(PARlevel=PAR.labs,Sp=Sp.labs))+
  scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),palette="Spectral",direction=1)+
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
tiff(paste("Tchamber_Tleafmod_byPAR.tiff"), width=4975, height=3900,res=600,units="px",compression="lzw")
ggplot(data=filter(Leaftemp,round> 0), aes(x=Tc,y=Tleaf,colour=E))+
  geom_point(shape=16,alpha=1,size=0.5)+
  geom_line(aes(x=Tc,y=Tc),colour="black",lty="dashed",size=1)+
  geom_smooth(method="lm")+
  facet_grid(PARlevel~Sp,
             labeller=labeller(PARlevel=PAR.labs,Sp=Sp.labs))+
  scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),palette="Spectral",direction=1)+
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
tiff(paste("Tchamber_Tleaf_Tc_byPAR.tiff"), width=4975, height=3900,res=600,units="px",compression="lzw")
ggplot(data=filter(Leaftemp,round> 0), aes(x=Tc,y=T_Tc,colour=E))+
  geom_point(shape=16,alpha=1,size=0.5)+
  geom_line(aes(x=Tc,y=Tc),colour="black",lty="dashed",size=1)+
  geom_smooth(method="lm")+
  facet_grid(PARlevel~Sp,
             labeller=labeller(PARlevel=PAR.labs,Sp=Sp.labs))+
  scale_colour_distiller(name=TeX("\\textit{E} (mmol m$^{-2}$ s$^{-1}$)"),palette="Spectral",direction=1)+
  scale_y_continuous(name="Thermocouple leaf temperature (ºC)",limits=c(5,55),breaks=c(10,20,30,40,50))+
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