# VPDchamber
Codes, functions and data for calculation of gas exchange in whole-plant chambers


# a) Input data

DataCIRAS.csv
  Extracted data from first sheet in excel files VPD...xlsx, Includes data from CIRAS3DC, PAR, and T/RH sensors
  
DataCR1000.csv
  Extracted data from second sheet in excel files VPD...xlsx, Includes data fron logger CR1000 (Flow sensors and Thermocouples)
  
DataZ6.csv
  Extracted data from first sheet in excel files z6...xlsx, Data from soil water potential probes (Teros 21)
  
LEAFAREA.csv
  Input data for each tree: ID, characteristic dimension of the leaf (d,cm), mean leaf area (cm2), total leaf area per tree (cm2)
  
Alldata.csv
  Combined dataset including all input data from different sources, grouped by 10 min slots


# b) Outputs

Alldatacalc3d.csv
  Initial output file including E, A, Tleaf, VPDleaf, gtotal and gs calculations, grouped by 10 min. Each chamber has its own column
  
CombicolsD.csv
  Rearranged dataset with chamber, plantID and species as factors, grouped by 10 min. One single column per variable
  
dataset.csv
  Final working dataset used for figures and calculations, hourly means of the key output variables


# c) Codes

VPDATOR_calculation.R
  Main R script for data arrangement, statistics and calculations, including the functions used to determine boundary layer conductance, modelled leaf temperature and chamber gas exchange
  
VPDATOR_pathanalysis.R
  Additional R script for path analysis
  
