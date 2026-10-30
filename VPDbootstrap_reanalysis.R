library(dplyr)
library(minpack.lm)
library(rstudioapi)
library(ggplot2)

#### SET DIR ####
## Setting the working directory to the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
# stores mainpath
mainpath<-getwd()

#### TAB1 PARAMETER BOOTSTRAP ####
##### P12,P88,PEAK,SLOPE #####
# redo bootstrap including PlantID #

library(readxl); library(dplyr); library(purrr); library(tidyr)
library(minpack.lm); library(multcompView)
set.seed(123); B <- 2000; alpha <- 0.05; ngrid <- 5000

#####  READ DATA #####
raw <- read.csv("dataset.csv",header=TRUE)
dat <- raw %>%
  dplyr::transmute(species=Sp, IDplant=IDplant, VPD=VPDleaf,
                   gs=gs, E=E, AN=A, WUE=WUE, PARc=PARc, WPgroup=WPgroup) %>%
  dplyr::filter(PARc > 500, WPgroup == "(0, -0.1 MPa)") %>%
  dplyr::mutate(species=dplyr::recode(species, fag="faginea", ile="rotundifolia",
                                      coc="coccifera", .default=as.character(species)),
                species=factor(species, levels=c("faginea","rotundifolia","coccifera")))

# verify
print(dat %>% dplyr::count(species, IDplant, name="n_obs"))
print(dat %>% group_by(species) %>% dplyr::summarise(n_plants=n_distinct(IDplant), n_obs=n()))

##### MODELS #####
logistic3 <- function(x,a,x0,b) a/(1+exp((x-x0)/b))
exp_offset <- function(x,a,b,c) c+a*exp(b*x)
gauss_peak <- function(x,A,mu,sigma,c) c+A*exp(-0.5*((x-mu)/sigma)^2)

first_down_crossing <- function(x,y,target){
  z<-y-target; idx<-which(z[-length(z)]>=0 & z[-1]<=0)
  if(!length(idx)) return(NA_real_); i<-idx[1]
  x1<-x[i]; x2<-x[i+1]; y1<-y[i]; y2<-y[i+1]
  if(y2==y1) return(x1)
  x1+(target-y1)*(x2-x1)/(y2-y1)
}

derived_decline <- function(fit,model,xmin,xmax,ngrid=5000){
  xg<-seq(xmin,xmax,length.out=ngrid); cf<-coef(fit)
  if(model=="logistic") {yg<-logistic3(xg,cf["a"],cf["x0"],cf["b"]); ymin_asym<-0}
  else {yg<-exp_offset(xg,cf["a"],cf["b"],cf["c"]); ymin_asym<-cf["c"];
  if(cf["b"]>=0) return(c(P88=NA,P12=NA,slope_max=NA))}
  ymax<-max(yg,na.rm=TRUE); amp<-ymax-ymin_asym
  if(!is.finite(amp)||amp<=0) return(c(P88=NA,P12=NA,slope_max=NA))
  P88<-first_down_crossing(xg,yg,ymin_asym+0.88*amp)
  P12<-first_down_crossing(xg,yg,ymin_asym+0.12*amp)
  slope_max<-max(abs(diff(yg)/diff(xg)),na.rm=TRUE)
  c(P88=P88,P12=P12,slope_max=slope_max)
}

fit_logistic <- function(d,yvar){
  d<-d %>% dplyr::filter(is.finite(VPD),is.finite(.data[[yvar]])); y<-d[[yvar]]; x<-d$VPD
  nlsLM(as.formula(paste0(yvar," ~ logistic3(VPD,a,x0,b)")),data=d,
        start=list(a=max(y),x0=median(x),b=max(diff(range(x))/5,.05)),
        lower=c(a=0,x0=min(x),b=1e-5),upper=c(a=Inf,x0=max(x),b=Inf),
        control=nls.lm.control(maxiter=1000))
}
fit_exp <- function(d,yvar="WUE"){
  d<-d %>% dplyr::filter(is.finite(VPD),is.finite(.data[[yvar]])); y<-d[[yvar]]; x<-d$VPD
  c0<-max(0,min(y)*.8)
  nlsLM(as.formula(paste0(yvar," ~ exp_offset(VPD,a,b,c)")),data=d,
        start=list(a=max(y)-c0,b=-.3,c=c0),lower=c(a=0,b=-Inf,c=-Inf),
        upper=c(a=Inf,b=0,c=Inf),control=nls.lm.control(maxiter=1000))
}
fit_peak <- function(d,yvar="E"){
  d<-d %>% dplyr::filter(is.finite(VPD),is.finite(.data[[yvar]])); y<-d[[yvar]]; x<-d$VPD
  mu0<-x[which.max(y)][1]; c0<-min(y)
  nlsLM(as.formula(paste0(yvar," ~ gauss_peak(VPD,A,mu,sigma,c)")),data=d,
        start=list(A=max(y)-c0,mu=mu0,sigma=max(sd(x),.2),c=c0),
        lower=c(A=0,mu=min(x),sigma=1e-5,c=-Inf),
        upper=c(A=Inf,mu=max(x),sigma=Inf,c=Inf),control=nls.lm.control(maxiter=1000))
}
extract_metrics <- function(d,variable){
  xmin<-min(d$VPD,na.rm=TRUE); xmax<-max(d$VPD,na.rm=TRUE)
  if(variable %in% c("gs","AN")) {m<-derived_decline(fit_logistic(d,variable),"logistic",xmin,xmax,ngrid); return(m)}
  if(variable=="WUE") {m<-derived_decline(fit_exp(d),"exp",xmin,xmax,ngrid); return(c(P12=m["P12"],slope_max=m["slope_max"]))}
  if(variable=="E") {fit<-fit_peak(d); mu<-unname(coef(fit)["mu"]); if(mu<xmin||mu>xmax) mu<-NA_real_; return(c(VPDpeak=mu))}
}

##### Cluster bootstrap: IDplant #####
resample_plants <- function(d){
  plants<-unique(d$IDplant); draws<-sample(plants,length(plants),replace=TRUE)
  map2_dfr(draws,seq_along(draws),function(id,k){
    d %>% dplyr::filter(IDplant==id) %>% dplyr::mutate(boot_cluster=paste0(id,"__draw",k))
  })
}
bootstrap_one <- function(d,variable,B=10000){
  d<-d %>% dplyr::filter(is.finite(VPD),is.finite(.data[[variable]]))
  np<-n_distinct(d$IDplant); if(np<2) stop("Menos de 2 plantas válidas")
  out<-vector("list",B)
  for(b in seq_len(B)){
    db<-resample_plants(d)
    out[[b]]<-tryCatch(extract_metrics(db,variable),error=function(e) NULL)
  }
  good<-out[!vapply(out,is.null,logical(1))]
  bind_rows(lapply(good,as.list)) %>% dplyr::mutate(iter=row_number())
}

##### global & bootstrap x plant #####
variables<-c("gs","E","AN","WUE"); species_names<-levels(dat$species)
boot<-list(); central<-list()
for(v in variables) for(sp in species_names){
  d<-dat %>% dplyr::filter(species==sp); key<-paste(v,sp,sep="__")
  central[[key]]<-tryCatch(extract_metrics(d,v),error=function(e) NULL)
  boot[[key]]<-bootstrap_one(d,v,B)
}
metrics_wanted<-list(gs=c("P88","P12","slope_max"),E="VPDpeak",
                     AN=c("P88","P12","slope_max"),WUE=c("P12","slope_max"))
summary_rows<-list()
for(v in names(metrics_wanted)) for(sp in species_names){
  key<-paste(v,sp,sep="__")
  for(m in metrics_wanted[[v]]){
    z<-boot[[key]][[m]]; z<-z[is.finite(z)]
    est<-if(!is.null(central[[key]])) unname(central[[key]][m]) else NA_real_
    summary_rows[[length(summary_rows)+1]]<-tibble(variable=v,species=sp,metric=m,
                                                   estimate=est,CI_low=quantile(z,.025,na.rm=TRUE),CI_high=quantile(z,.975,na.rm=TRUE),
                                                   n_boot_valid=length(z))
  }
}
summary_tbl<-bind_rows(summary_rows)

write.csv(summary_tbl,"summary_tbl.csv",row.names=FALSE)

##### species comparison #####
pairwise_boot <- function(variable,metric,sp1,sp2){
  z1<-boot[[paste(variable,sp1,sep="__")]][[metric]]; z1<-z1[is.finite(z1)]
  z2<-boot[[paste(variable,sp2,sep="__")]][[metric]]; z2<-z2[is.finite(z2)]
  n<-min(length(z1),length(z2)); if(n<100) return(tibble(diff=NA,low=NA,high=NA,p=NA))
  dd<-sample(z1,n,replace=FALSE)-sample(z2,n,replace=FALSE)
  tibble(diff=median(dd),low=quantile(dd,.025),high=quantile(dd,.975),
         p=min(1,2*min(mean(dd<=0),mean(dd>=0))))
}
pairs<-combn(species_names,2,simplify=FALSE); pw_rows<-list()
for(v in names(metrics_wanted)) for(m in metrics_wanted[[v]]) for(pp in pairs){
  tmp<-pairwise_boot(v,m,pp[1],pp[2])
  pw_rows[[length(pw_rows)+1]]<-tmp %>% dplyr::mutate(variable=v,metric=m,sp1=pp[1],sp2=pp[2],.before=1)
}
pairwise_tbl<-bind_rows(pw_rows)

write.csv(pairwise_tbl,"pairwise_tbl.csv",row.names=FALSE)

##### letters CLD #####
get_letters <- function(vv,mm){
  pw<-pairwise_tbl %>% dplyr::filter(variable==vv,metric==mm)
  pv<-setNames(pw$p,paste(pw$sp1,pw$sp2,sep="-")); sig<-pv<alpha
  L<-multcompView::multcompLetters(sig)$Letters
  tibble(species=names(L),letter=unname(L))
}

for(v in variables) for(sp in species_names){
  letters_tbl<-purrr::map_dfr(names(metrics_wanted),function(v)
    purrr::map_dfr(metrics_wanted[[v]],function(m) get_letters(v,m) %>% dplyr::mutate(variable=v,metric=m,.before=1)))
  final_tbl<-summary_tbl %>% left_join(letters_tbl,by=c("variable","metric","species")) %>%
    dplyr::mutate(result=sprintf("%.2f (%.2f–%.2f) %s",estimate,CI_low,CI_high,letter))
  write.csv(final_tbl,"Table1_plant_bootstrap_results.csv",row.names=FALSE)
  write.csv(pairwise_tbl,"Table1_plant_bootstrap_pairwise.csv",row.names=FALSE)
  print(final_tbl); print(pairwise_tbl)
}

#### WHOLE-CURVE BOOTSTRAP ####

#### gs-VPD CURVE COMPARISONS ####
# THREE-PARAMETER LOGISTIC CURVES WITH PLANT-LEVEL WILD CLUSTER BOOTSTRAP

#
# Null hypothesis:
#   All species follow one common gs-VPD response curve.
#
# Alternative hypothesis:
#   Each species follows a separate gs-VPD response curve.
#
# Test statistic:
#   Delta RSS = RSS_pooled - RSS_species-specific
#
# Repeated observations are accounted for using a wild cluster bootstrap.
# One Rademacher multiplier (-1 or +1) is assigned to all observations
# belonging to the same plant.
#
# Required packages:
#   dplyr
#   minpack.lm
#
# Optional package for the final figure:
#   ggplot2
#

# 0. GENERAL SETTINGS
base::set.seed(123)

# Use B = 100 while testing the script.
# Use B >= 1000, preferably 2000, for the final analysis.
B <- 1000L

input_file <- "dataset.csv"

# Exact label of the soil-water-potential category to retain
selected_WPgroup <- "(0, -0.1 MPa)"

# Minimum PAR used to define non-light-limiting conditions
PAR_threshold <- 500

# Minimum acceptable fraction of successful bootstrap iterations
minimum_convergence_rate <- 0.90

# Output files
global_output_file <-"global_gs_VPD_curve_test.csv"
pairwise_output_file <-"pairwise_gs_VPD_curve_tests.csv"
adjusted_pairwise_output_file <-"pairwise_gs_VPD_curve_tests_adjusted.csv"
parameter_output_file <-"observed_gs_logistic3_curve_parameters.csv"
global_bootstrap_output_file <-"global_gs_bootstrap_distribution.csv"
prediction_output_file <-"fitted_gs_VPD_curves.csv"
figure_output_file <-"fitted_gs_VPD_curves.tiff"

# 1. CHECK REQUIRED PACKAGES


required_packages <- c(
  "dplyr",
  "minpack.lm"
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
    "Install the following packages before running the analysis: ",
    base::paste(missing_packages, collapse = ", "),
    "\nFor example:\n",
    "install.packages(c(",
    base::paste(
      base::shQuote(missing_packages),
      collapse = ", "
    ),
    "))"
  )
}



# 2. READ THE ORIGINAL DATA


if (!base::file.exists(input_file)) {
  
  base::stop(
    "The input file was not found in the working directory: ",
    input_file,
    "\nCurrent working directory: ",
    base::getwd()
  )
}

raw <- utils::read.csv(
  file = input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)



# 3. CHECK REQUIRED COLUMNS


required_columns <- c(
  "Sp",
  "IDplant",
  "VPDleaf",
  "gs",
  "E",
  "A",
  "WUE",
  "PARc",
  "WPgroup"
)

missing_columns <- base::setdiff(
  required_columns,
  base::names(raw)
)

if (base::length(missing_columns) > 0L) {
  
  base::stop(
    "The following required columns are absent from the dataset: ",
    base::paste(
      missing_columns,
      collapse = ", "
    )
  )
}



# 4. INSPECT AND VERIFY WPgroup


base::cat(
  "\n============================================================\n",
  "AVAILABLE WPgroup VALUES\n",
  "============================================================\n",
  sep = ""
)

WPgroup_counts <- dplyr::count(
  raw,
  WPgroup,
  sort = TRUE,
  name = "n_observations"
)

base::print(WPgroup_counts)

if (!selected_WPgroup %in% raw$WPgroup) {
  
  base::stop(
    "\nThe requested WPgroup label was not found: ",
    selected_WPgroup,
    "\nCheck spelling, spaces, parentheses, and decimal notation."
  )
}



# 5. DATA PREPARATION


dat <- raw |>
  dplyr::transmute(
    species = Sp,
    IDplant = IDplant,
    VPD = VPDleaf,
    gs = gs,
    E = E,
    AN = A,
    WUE = WUE,
    PARc = PARc,
    WPgroup = WPgroup
  ) |>
  dplyr::filter(
    PARc > PAR_threshold,
    WPgroup == selected_WPgroup,
    base::is.finite(VPD),
    base::is.finite(gs),
    !base::is.na(species),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    species = dplyr::recode(
      base::as.character(species),
      fag = "faginea",
      ile = "rotundifolia",
      coc = "coccifera",
      .default = base::as.character(species)
    ),
    species = base::factor(
      species,
      levels = c(
        "faginea",
        "rotundifolia",
        "coccifera"
      )
    )
  ) |>
  dplyr::filter(
    !base::is.na(species)
  ) |>
  dplyr::mutate(
    # This identifier remains unique even if plant codes are repeated
    # among species.
    PlantCluster = base::interaction(
      species,
      IDplant,
      drop = TRUE,
      sep = "_"
    )
  ) |>
  base::droplevels()



# 6. VALIDATE THE FILTERED DATA


if (base::nrow(dat) == 0L) {
  
  base::stop(
    "No observations remain after filtering. ",
    "Check PARc, WPgroup, VPDleaf, and gs."
  )
}

data_summary <- dat |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    n_plants = dplyr::n_distinct(PlantCluster),
    n_obs = dplyr::n(),
    n_VPD_values = dplyr::n_distinct(VPD),
    min_VPD = base::min(VPD, na.rm = TRUE),
    max_VPD = base::max(VPD, na.rm = TRUE),
    min_gs = base::min(gs, na.rm = TRUE),
    max_gs = base::max(gs, na.rm = TRUE),
    .groups = "drop"
  )

base::cat(
  "\n============================================================\n",
  "FILTERED-DATA SUMMARY\n",
  "============================================================\n",
  sep = ""
)

base::print(data_summary)

plant_summary <- dplyr::count(
  dat,
  species,
  IDplant,
  PlantCluster,
  name = "n_obs"
)

base::cat("\nObservations per plant:\n")
base::print(plant_summary)

expected_species <- c(
  "faginea",
  "rotundifolia",
  "coccifera"
)

observed_species <- base::unique(
  base::as.character(dat$species)
)

absent_species <- base::setdiff(
  expected_species,
  observed_species
)

if (base::length(absent_species) > 0L) {
  
  base::stop(
    "The following species have no observations after filtering: ",
    base::paste(absent_species, collapse = ", ")
  )
}

if (base::any(data_summary$n_plants < 2L)) {
  
  base::stop(
    "At least one species has fewer than two plant clusters."
  )
}

if (base::any(data_summary$n_obs < 6L)) {
  
  base::stop(
    "At least one species has fewer than six valid observations."
  )
}

if (base::any(data_summary$n_VPD_values < 4L)) {
  
  base::stop(
    "At least one species has fewer than four distinct VPD values."
  )
}



# 7. THREE-PARAMETER LOGISTIC FUNCTION

#
# a  = upper asymptote, approximately maximum gs
# x0 = VPD at the midpoint of the decline
# b  = response-width parameter
#
# With b > 0, gs declines as VPD increases.
#


logistic3 <- function(x, a, x0, b) {
  
  a / (
    1 + base::exp(
      (x - x0) / b
    )
  )
}



# 8. FIT ONE LOGISTIC CURVE

#
# Several starting-value combinations are tested. The converged fit with
# the smallest RSS is retained.
#


fit_one_curve <- function(d) {
  
  d <- d |>
    dplyr::filter(
      base::is.finite(VPD),
      base::is.finite(gs)
    )
  
  if (base::nrow(d) == 0L) {
    
    base::stop(
      "fit_one_curve() received an empty dataset."
    )
  }
  
  if (base::nrow(d) < 6L) {
    
    base::stop(
      "Too few observations to fit a three-parameter curve: ",
      base::nrow(d)
    )
  }
  
  if (dplyr::n_distinct(d$VPD) < 4L) {
    
    base::stop(
      "Too few distinct VPD values to fit a logistic curve."
    )
  }
  
  x <- d$VPD
  y <- d$gs
  
  x_min <- base::min(
    x,
    na.rm = TRUE
  )
  
  x_max <- base::max(
    x,
    na.rm = TRUE
  )
  
  x_range <- x_max - x_min
  
  y_max <- base::max(
    y,
    na.rm = TRUE
  )
  
  values_to_check <- c(
    x_min,
    x_max,
    x_range,
    y_max
  )
  
  if (!base::all(
    base::is.finite(values_to_check)
  )) {
    
    base::stop(
      "The fitting data generated non-finite limits or starting values."
    )
  }
  
  if (x_range <= 0) {
    
    base::stop(
      "VPD has no usable variation."
    )
  }
  
  # Ensure a positive initial upper asymptote
  a_base <- base::max(
    y_max,
    1e-4
  )
  
  x0_candidates <- base::as.numeric(
    stats::quantile(
      x,
      probs = c(
        0.25,
        0.35,
        0.50,
        0.65,
        0.75
      ),
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  b_candidates <- c(
    base::max(x_range / 20, 0.01),
    base::max(x_range / 10, 0.02),
    base::max(x_range / 5, 0.05),
    base::max(x_range / 2, 0.10)
  )
  
  a_candidates <- c(
    a_base,
    a_base * 1.25,
    a_base * 1.75,
    a_base * 2
  )
  
  start_grid <- base::expand.grid(
    a = a_candidates,
    x0 = x0_candidates,
    b = b_candidates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  best_fit <- NULL
  best_RSS <- Inf
  
  for (
    jj in base::seq_len(
      base::nrow(start_grid)
    )
  ) {
    
    current_fit <- base::tryCatch(
      
      minpack.lm::nlsLM(
        formula = gs ~ logistic3(
          VPD,
          a,
          x0,
          b
        ),
        data = d,
        start = list(
          a = start_grid$a[jj],
          x0 = start_grid$x0[jj],
          b = start_grid$b[jj]
        ),
        lower = c(
          a = 0,
          x0 = x_min,
          b = 1e-5
        ),
        upper = c(
          a = Inf,
          x0 = x_max,
          b = Inf
        ),
        control = minpack.lm::nls.lm.control(
          maxiter = 1000,
          ftol = 1e-10,
          ptol = 1e-10
        )
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    if (!base::is.null(current_fit)) {
      
      current_residuals <- stats::residuals(
        current_fit
      )
      
      current_RSS <- base::sum(
        current_residuals^2
      )
      
      if (
        base::is.finite(current_RSS) &&
        current_RSS < best_RSS
      ) {
        
        best_fit <- current_fit
        best_RSS <- current_RSS
      }
    }
  }
  
  if (base::is.null(best_fit)) {
    
    base::stop(
      "The logistic curve failed for all starting-value combinations."
    )
  }
  
  list(
    fit = best_fit,
    RSS = best_RSS,
    coefficients = stats::coef(best_fit),
    n_obs = base::nrow(d)
  )
}



# 9. POOLED AND SPECIES-SPECIFIC MODELS


fit_pooled <- function(d) {
  
  if (base::nrow(d) == 0L) {
    
    base::stop(
      "fit_pooled() received an empty dataset."
    )
  }
  
  fit_one_curve(d)
}


fit_separate <- function(d) {
  
  d <- base::droplevels(d)
  
  species_data <- base::split(
    d,
    d$species,
    drop = TRUE
  )
  
  group_sizes <- base::vapply(
    species_data,
    base::nrow,
    FUN.VALUE = integer(1)
  )
  
  if (base::length(species_data) < 2L) {
    
    base::stop(
      "Species-specific fitting requires at least two non-empty species."
    )
  }
  
  if (base::any(group_sizes == 0L)) {
    
    base::stop(
      "An empty species group was passed to fit_separate()."
    )
  }
  
  fits <- base::lapply(
    species_data,
    fit_one_curve
  )
  
  separate_RSS <- base::sum(
    base::vapply(
      fits,
      function(z) {
        z$RSS
      },
      FUN.VALUE = numeric(1)
    )
  )
  
  list(
    fits = fits,
    RSS = separate_RSS
  )
}



# 10. FIT THE OBSERVED GLOBAL MODELS


H0 <- fit_pooled(dat)
H1 <- fit_separate(dat)

Delta_RSS_obs <- H0$RSS - H1$RSS

base::cat(
  "\n============================================================\n",
  "OBSERVED GLOBAL CURVE COMPARISON\n",
  "============================================================\n",
  sep = ""
)

base::cat(
  "RSS, pooled curve            = ",
  H0$RSS,
  "\n",
  sep = ""
)

base::cat(
  "RSS, species-specific curves = ",
  H1$RSS,
  "\n",
  sep = ""
)

base::cat(
  "Observed Delta RSS           = ",
  Delta_RSS_obs,
  "\n",
  sep = ""
)



# 11. EXTRACT OBSERVED CURVE PARAMETERS


pooled_parameters <- base::data.frame(
  Model = "Pooled",
  Species = "All species",
  a = base::unname(
    H0$coefficients["a"]
  ),
  x0 = base::unname(
    H0$coefficients["x0"]
  ),
  b = base::unname(
    H0$coefficients["b"]
  ),
  RSS = H0$RSS,
  stringsAsFactors = FALSE
)

species_parameters <- dplyr::bind_rows(
  base::lapply(
    base::names(H1$fits),
    function(sp) {
      
      current_fit <- H1$fits[[sp]]
      
      base::data.frame(
        Model = "Species-specific",
        Species = sp,
        a = base::unname(
          current_fit$coefficients["a"]
        ),
        x0 = base::unname(
          current_fit$coefficients["x0"]
        ),
        b = base::unname(
          current_fit$coefficients["b"]
        ),
        RSS = current_fit$RSS,
        stringsAsFactors = FALSE
      )
    }
  )
)

observed_parameters <- dplyr::bind_rows(
  pooled_parameters,
  species_parameters
)

base::cat("\nObserved logistic parameters:\n")
base::print(observed_parameters)

utils::write.csv(
  observed_parameters,
  file = parameter_output_file,
  row.names = FALSE
)



# 12. PREPARE RESIDUALS UNDER THE GLOBAL NULL MODEL


dat$fit_H0 <- base::as.numeric(
  stats::predict(
    H0$fit,
    newdata = dat
  )
)

dat$res_H0 <- dat$gs - dat$fit_H0

# Centre residuals under the null model
dat$res_H0 <- dat$res_H0 -
  base::mean(
    dat$res_H0,
    na.rm = TRUE
  )

plants <- base::unique(
  base::as.character(
    dat$PlantCluster
  )
)



# 13. GLOBAL WILD CLUSTER BOOTSTRAP


boot_delta <- base::rep(
  NA_real_,
  B
)

boot_error <- base::rep(
  "",
  B
)

for (bb in base::seq_len(B)) {
  
  plant_sign <- base::sample(
    x = c(-1, 1),
    size = base::length(plants),
    replace = TRUE
  )
  
  base::names(plant_sign) <- plants
  
  db <- dat
  
  cluster_multiplier <- plant_sign[
    base::as.character(
      db$PlantCluster
    )
  ]
  
  if (base::any(
    base::is.na(cluster_multiplier)
  )) {
    
    base::stop(
      "At least one plant cluster did not receive a bootstrap multiplier."
    )
  }
  
  db$gs <- db$fit_H0 +
    cluster_multiplier *
    db$res_H0
  
  boot_delta[bb] <- base::tryCatch(
    
    {
      m0 <- fit_pooled(db)
      m1 <- fit_separate(db)
      
      m0$RSS - m1$RSS
    },
    
    error = function(e) {
      
      boot_error[bb] <<-
        base::conditionMessage(e)
      
      NA_real_
    }
  )
  
  if (bb %% 100L == 0L) {
    
    base::cat(
      "\rGlobal bootstrap replicate ",
      bb,
      " of ",
      B,
      sep = ""
    )
  }
}

base::cat("\n")

valid_global <- base::is.finite(
  boot_delta
)

boot_delta_valid <- boot_delta[
  valid_global
]

B_valid_global <- base::length(
  boot_delta_valid
)

if (B_valid_global == 0L) {
  
  global_error_table <- base::sort(
    base::table(
      boot_error[
        boot_error != ""
      ]
    ),
    decreasing = TRUE
  )
  
  base::print(global_error_table)
  
  base::stop(
    "All global bootstrap replicates failed."
  )
}

p_global <- (
  1 +
    base::sum(
      boot_delta_valid >= Delta_RSS_obs
    )
) / (
  1 + B_valid_global
)

global_convergence_rate <-
  B_valid_global / B

base::cat(
  "\n============================================================\n",
  "GLOBAL WILD CLUSTER BOOTSTRAP\n",
  "============================================================\n",
  sep = ""
)

base::cat(
  "Requested bootstrap replicates = ",
  B,
  "\n",
  sep = ""
)

base::cat(
  "Valid bootstrap replicates     = ",
  B_valid_global,
  "\n",
  sep = ""
)

base::cat(
  "Convergence rate               = ",
  base::round(
    global_convergence_rate,
    4
  ),
  "\n",
  sep = ""
)

base::cat(
  "Global cluster-bootstrap P     = ",
  p_global,
  "\n",
  sep = ""
)

if (
  global_convergence_rate <
  minimum_convergence_rate
) {
  
  base::warning(
    "Fewer than ",
    minimum_convergence_rate * 100,
    "% of the global bootstrap fits converged. ",
    "Interpret the bootstrap result cautiously."
  )
}



# 14. SUMMARISE GLOBAL BOOTSTRAP ERRORS


global_error_table <- base::sort(
  base::table(
    boot_error[
      boot_error != ""
    ]
  ),
  decreasing = TRUE
)

if (base::length(global_error_table) > 0L) {
  
  base::cat(
    "\nMost frequent global-bootstrap errors:\n"
  )
  
  base::print(
    utils::head(
      global_error_table,
      n = 10L
    )
  )
}



# 15. EXPORT GLOBAL RESULTS


global_bootstrap_distribution <-
  base::data.frame(
    replicate = base::which(
      valid_global
    ),
    Delta_RSS = boot_delta_valid
  )

utils::write.csv(
  global_bootstrap_distribution,
  file = global_bootstrap_output_file,
  row.names = FALSE
)

global_results <- base::data.frame(
  Comparison = "All species",
  RSS_pooled = H0$RSS,
  RSS_species_specific = H1$RSS,
  Delta_RSS = Delta_RSS_obs,
  P_bootstrap = p_global,
  B_requested = B,
  B_valid = B_valid_global,
  convergence_rate =
    global_convergence_rate,
  stringsAsFactors = FALSE
)

utils::write.csv(
  global_results,
  file = global_output_file,
  row.names = FALSE
)



# 16. PAIRWISE CURVE-COMPARISON FUNCTION


pairwise_curve_test <- function(
    data,
    sp1,
    sp2,
    B = 1000L,
    seed = 123L
) {
  
  base::set.seed(seed)
  
  d <- data |>
    dplyr::filter(
      species %in% c(
        sp1,
        sp2
      )
    ) |>
    base::droplevels()
  
  species_counts <- dplyr::count(
    d,
    species,
    name = "n_obs"
  )
  
  if (base::nrow(species_counts) != 2L) {
    
    base::stop(
      "One or both requested species are absent after filtering: ",
      sp1,
      " versus ",
      sp2
    )
  }
  
  if (base::any(
    species_counts$n_obs < 6L
  )) {
    
    base::stop(
      "At least one species has fewer than six observations: ",
      sp1,
      " versus ",
      sp2
    )
  }
  
  plant_counts <- d |>
    dplyr::distinct(
      species,
      PlantCluster
    ) |>
    dplyr::count(
      species,
      name = "n_plants"
    )
  
  if (base::any(
    plant_counts$n_plants < 2L
  )) {
    
    base::stop(
      "At least one species has fewer than two plant clusters: ",
      sp1,
      " versus ",
      sp2
    )
  }
  
  # Observed pairwise models
  m0 <- fit_pooled(d)
  m1 <- fit_separate(d)
  
  observed_delta <- m0$RSS - m1$RSS
  
  # Null-model fitted values and residuals
  d$fit_H0 <- base::as.numeric(
    stats::predict(
      m0$fit,
      newdata = d
    )
  )
  
  d$res_H0 <- d$gs - d$fit_H0
  
  d$res_H0 <- d$res_H0 -
    base::mean(
      d$res_H0,
      na.rm = TRUE
    )
  
  pair_plants <- base::unique(
    base::as.character(
      d$PlantCluster
    )
  )
  
  bootstrap_delta <- base::rep(
    NA_real_,
    B
  )
  
  bootstrap_error <- base::rep(
    "",
    B
  )
  
  for (
    bb in base::seq_len(B)
  ) {
    
    plant_sign <- base::sample(
      x = c(-1, 1),
      size = base::length(
        pair_plants
      ),
      replace = TRUE
    )
    
    base::names(plant_sign) <-
      pair_plants
    
    db <- d
    
    cluster_multiplier <-
      plant_sign[
        base::as.character(
          db$PlantCluster
        )
      ]
    
    if (base::any(
      base::is.na(cluster_multiplier)
    )) {
      
      base::stop(
        "At least one plant cluster did not receive a bootstrap multiplier."
      )
    }
    
    db$gs <- db$fit_H0 +
      cluster_multiplier *
      db$res_H0
    
    bootstrap_delta[bb] <-
      base::tryCatch(
        
        {
          q0 <- fit_pooled(db)
          q1 <- fit_separate(db)
          
          q0$RSS - q1$RSS
        },
        
        error = function(e) {
          
          bootstrap_error[bb] <<-
            base::conditionMessage(e)
          
          NA_real_
        }
      )
  }
  
  valid_bootstrap <-
    base::is.finite(
      bootstrap_delta
    )
  
  bootstrap_delta_valid <-
    bootstrap_delta[
      valid_bootstrap
    ]
  
  B_valid <- base::length(
    bootstrap_delta_valid
  )
  
  if (B_valid == 0L) {
    
    error_table <- base::sort(
      base::table(
        bootstrap_error[
          bootstrap_error != ""
        ]
      ),
      decreasing = TRUE
    )
    
    base::print(error_table)
    
    base::stop(
      "All bootstrap replicates failed for ",
      sp1,
      " versus ",
      sp2
    )
  }
  
  p_value <- (
    1 +
      base::sum(
        bootstrap_delta_valid >=
          observed_delta
      )
  ) / (
    1 + B_valid
  )
  
  convergence_rate <- B_valid / B
  
  if (
    convergence_rate <
    minimum_convergence_rate
  ) {
    
    base::warning(
      "Fewer than ",
      minimum_convergence_rate * 100,
      "% of the bootstrap fits converged for ",
      sp1,
      " versus ",
      sp2,
      "."
    )
  }
  
  error_table <- base::sort(
    base::table(
      bootstrap_error[
        bootstrap_error != ""
      ]
    ),
    decreasing = TRUE
  )
  
  most_common_error <- if (
    base::length(error_table) > 0L
  ) {
    
    base::names(error_table)[1]
    
  } else {
    
    NA_character_
  }
  
  base::data.frame(
    species1 = sp1,
    species2 = sp2,
    RSS_pooled = m0$RSS,
    RSS_separate = m1$RSS,
    Delta_RSS = observed_delta,
    P_bootstrap = p_value,
    B_requested = B,
    B_valid = B_valid,
    convergence_rate =
      convergence_rate,
    most_common_error =
      most_common_error,
    stringsAsFactors = FALSE
  )
}



# 17. RUN ALL PAIRWISE COMPARISONS


pairwise_results <- dplyr::bind_rows(
  
  pairwise_curve_test(
    data = dat,
    sp1 = "faginea",
    sp2 = "rotundifolia",
    B = B,
    seed = 101L
  ),
  
  pairwise_curve_test(
    data = dat,
    sp1 = "faginea",
    sp2 = "coccifera",
    B = B,
    seed = 102L
  ),
  
  pairwise_curve_test(
    data = dat,
    sp1 = "rotundifolia",
    sp2 = "coccifera",
    B = B,
    seed = 103L
  )
)

base::cat(
  "\n============================================================\n",
  "PAIRWISE CURVE COMPARISONS\n",
  "============================================================\n",
  sep = ""
)

base::print(pairwise_results)

utils::write.csv(
  pairwise_results,
  file = pairwise_output_file,
  row.names = FALSE
)



# 18. MULTIPLE-COMPARISON ADJUSTMENT


pairwise_results_adjusted <-
  pairwise_results |>
  dplyr::mutate(
    P_Holm = stats::p.adjust(
      P_bootstrap,
      method = "holm"
    ),
    P_FDR = stats::p.adjust(
      P_bootstrap,
      method = "BH"
    )
  )

base::cat(
  "\nPairwise results with adjusted P-values:\n"
)

base::print(
  pairwise_results_adjusted
)

utils::write.csv(
  pairwise_results_adjusted,
  file =
    adjusted_pairwise_output_file,
  row.names = FALSE
)



# 19. GENERATE FITTED CURVES


prediction_grid <- base::data.frame(
  VPD = base::seq(
    from = base::min(
      dat$VPD,
      na.rm = TRUE
    ),
    to = base::max(
      dat$VPD,
      na.rm = TRUE
    ),
    length.out = 300L
  )
)

pooled_prediction <-
  prediction_grid |>
  dplyr::mutate(
    species = "Pooled",
    gs_predicted = base::as.numeric(
      stats::predict(
        H0$fit,
        newdata = prediction_grid
      )
    )
  )

species_predictions <- dplyr::bind_rows(
  base::lapply(
    base::names(H1$fits),
    function(sp) {
      
      current_fit <-
        H1$fits[[sp]]$fit
      
      prediction_grid |>
        dplyr::mutate(
          species = sp,
          gs_predicted =
            base::as.numeric(
              stats::predict(
                current_fit,
                newdata =
                  prediction_grid
              )
            )
        )
    }
  )
)

curve_predictions <- dplyr::bind_rows(
  pooled_prediction,
  species_predictions
)

utils::write.csv(
  curve_predictions,
  file = prediction_output_file,
  row.names = FALSE
)



# 20. OPTIONAL FIGURE WITH ggplot2


if (
  base::requireNamespace(
    "ggplot2",
    quietly = TRUE
  )
) {
  
  species_colours <- c(
    faginea = "#D55E00",
    rotundifolia = "#0072B2",
    coccifera = "#009E73",
    Pooled = "#000000"
  )
  
  plot_data <- dat |>
    dplyr::mutate(
      species = base::as.character(
        species
      )
    )
  
  curve_plot <- ggplot2::ggplot(
    data = plot_data,
    mapping = ggplot2::aes(
      x = VPD,
      y = gs,
      colour = species
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.25,
      size = 1
    ) +
    ggplot2::geom_line(
      data = curve_predictions,
      mapping = ggplot2::aes(
        x = VPD,
        y = gs_predicted,
        colour = species,
        linetype = species
      ),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = species_colours,
      breaks = c(
        "faginea",
        "rotundifolia",
        "coccifera",
        "Pooled"
      ),
      labels = c(
        "Q. faginea",
        "Q. rotundifolia",
        "Q. coccifera",
        "Pooled"
      )
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        faginea = "solid",
        rotundifolia = "solid",
        coccifera = "solid",
        Pooled = "dashed"
      ),
      breaks = c(
        "faginea",
        "rotundifolia",
        "coccifera",
        "Pooled"
      ),
      labels = c(
        "Q. faginea",
        "Q. rotundifolia",
        "Q. coccifera",
        "Pooled"
      )
    ) +
    ggplot2::labs(
      x = expression(
        VPD[leaf]
      ),
      y = expression(
        g[s]
      ),
      colour = "Curve",
      linetype = "Curve"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor =
        ggplot2::element_blank(),
      legend.position = "top"
    )
  
  base::print(curve_plot)
  
  ggplot2::ggsave(
    filename = figure_output_file,
    plot = curve_plot,
    device = "tiff",
    width = 180,
    height = 140,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
} else {
  
  base::message(
    "Package ggplot2 is not installed. ",
    "The statistical analysis was completed, ",
    "but the TIFF figure was not created."
  )
}



# 21. FINAL SUMMARY


base::cat(
  "\n============================================================\n",
  "FINAL ANALYSIS SUMMARY\n",
  "============================================================\n",
  sep = ""
)

base::cat("\nGlobal comparison:\n")
base::print(global_results)

base::cat(
  "\nPairwise comparisons with Holm correction:\n"
)

base::print(
  dplyr::select(
    pairwise_results_adjusted,
    species1,
    species2,
    Delta_RSS,
    P_bootstrap,
    P_Holm,
    B_valid,
    convergence_rate
  )
)

base::cat(
  "\nObserved logistic parameters:\n"
)

base::print(observed_parameters)

base::cat(
  "\nOutput files created:\n",
  "- ", global_output_file, "\n",
  "- ", pairwise_output_file, "\n",
  "- ", adjusted_pairwise_output_file, "\n",
  "- ", parameter_output_file, "\n",
  "- ", global_bootstrap_output_file, "\n",
  "- ", prediction_output_file, "\n",
  sep = ""
)

if (
  base::file.exists(
    figure_output_file
  )
) {
  
  base::cat(
    "- ",
    figure_output_file,
    "\n",
    sep = ""
  )
}



# 22. SAVE SESSION INFORMATION


utils::capture.output(
  utils::sessionInfo(),
  file = "sessionInfo_gs_VPD_bootstrap.txt"
)

base::cat(
  "- sessionInfo_gs_VPD_bootstrap.txt\n"
)

base::cat(
  "\nAnalysis completed.\n"
)


### AN-VPD CURVE COMPARISONS ####
# THREE-PARAMETER LOGISTIC CURVES WITH PLANT-LEVEL WILD CLUSTER BOOTSTRAP

# Null: all species share one common curve.
# Alternative: each species has its own curve.
# Test statistic: Delta RSS = RSS_pooled - RSS_species-specific.
# One Rademacher multiplier (-1 or +1) is assigned to all residuals from
# the same plant in each bootstrap replicate.


base::set.seed(123)
B <- 1000L                    # use 100 while testing; >=1000 for final analysis
input_file <- "dataset.csv"
selected_WPgroup <- "(0, -0.1 MPa)"
PAR_threshold <- 500
minimum_convergence_rate <- 0.90

# Variable-specific output files
global_output_file <- "global_AN_VPD_curve_test.csv"
pairwise_output_file <- "pairwise_AN_VPD_curve_tests.csv"
adjusted_pairwise_output_file <- "pairwise_AN_VPD_curve_tests_adjusted.csv"
parameter_output_file <- "observed_AN_logistic3_curve_parameters.csv"
global_bootstrap_output_file <- "global_AN_bootstrap_distribution.csv"
prediction_output_file <- "fitted_AN_VPD_curves.csv"
figure_output_file <- "fitted_AN_VPD_curves.tiff"
session_output_file <- "sessionInfo_AN_VPD_bootstrap.txt"


# 1. CHECK REQUIRED PACKAGES

required_packages <- c("dplyr", "minpack.lm")
missing_packages <- required_packages[
  !base::vapply(required_packages, base::requireNamespace,
                quietly = TRUE, FUN.VALUE = logical(1))
]
if (base::length(missing_packages) > 0L) {
  base::stop("Install before running: ",
             base::paste(missing_packages, collapse = ", "))
}


# 2. READ AND CHECK DATA

if (!base::file.exists(input_file)) {
  base::stop("Input file not found: ", input_file,
             "\nWorking directory: ", base::getwd())
}
raw <- utils::read.csv(input_file, header = TRUE,
                       stringsAsFactors = FALSE,
                       na.strings = c("", "NA", "NaN"))
required_columns <- c("Sp", "IDplant", "VPDleaf", "A",
                      "PARc", "WPgroup")
missing_columns <- base::setdiff(required_columns, base::names(raw))
if (base::length(missing_columns) > 0L) {
  base::stop("Missing columns: ", base::paste(missing_columns, collapse = ", "))
}

base::cat("\nAVAILABLE WPgroup VALUES\n")
WPgroup_counts <- dplyr::count(raw, WPgroup, sort = TRUE,
                               name = "n_observations")
base::print(WPgroup_counts)
if (!selected_WPgroup %in% raw$WPgroup) {
  base::stop("Requested WPgroup not found: ", selected_WPgroup)
}


# 3. DATA PREPARATION

dat <- raw |>
  dplyr::transmute(
    species = Sp,
    IDplant = IDplant,
    VPD = VPDleaf,
    response = A,
    PARc = PARc,
    WPgroup = WPgroup
  ) |>
  dplyr::filter(
    PARc > PAR_threshold,
    WPgroup == selected_WPgroup,
    base::is.finite(VPD),
    base::is.finite(response),
    response > 0,
    !base::is.na(species),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    species = dplyr::recode(
      base::as.character(species),
      fag = "faginea", ile = "rotundifolia", coc = "coccifera",
      .default = base::as.character(species)
    ),
    species = base::factor(species,
                           levels = c("faginea", "rotundifolia", "coccifera")),
    PlantCluster = base::interaction(species, IDplant,
                                     drop = TRUE, sep = "_")
  ) |>
  dplyr::filter(!base::is.na(species)) |>
  base::droplevels()


# 4. VALIDATE FILTERED DATA

if (base::nrow(dat) == 0L) base::stop("No observations remain after filtering.")
data_summary <- dat |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    n_plants = dplyr::n_distinct(PlantCluster),
    n_obs = dplyr::n(),
    n_VPD_values = dplyr::n_distinct(VPD),
    min_VPD = base::min(VPD, na.rm = TRUE),
    max_VPD = base::max(VPD, na.rm = TRUE),
    min_response = base::min(response, na.rm = TRUE),
    max_response = base::max(response, na.rm = TRUE),
    .groups = "drop"
  )
base::cat("\nFILTERED-DATA SUMMARY: AN\n")
base::print(data_summary)
base::print(dplyr::count(dat, species, IDplant, PlantCluster, name = "n_obs"))

expected_species <- c("faginea", "rotundifolia", "coccifera")
absent_species <- base::setdiff(expected_species,
                                base::unique(base::as.character(dat$species)))
if (base::length(absent_species) > 0L) {
  base::stop("Species absent after filtering: ",
             base::paste(absent_species, collapse = ", "))
}
if (base::any(data_summary$n_plants < 2L))
  base::stop("At least one species has fewer than two plant clusters.")
if (base::any(data_summary$n_obs < 6L))
  base::stop("At least one species has too few observations.")
if (base::any(data_summary$n_VPD_values < 4L))
  base::stop("At least one species has too few distinct VPD values.")


# 5. RESPONSE FUNCTION AND ROBUST SINGLE-CURVE FIT

logistic3 <- function(x, a, x0, b) {
  a / (1 + base::exp((x - x0) / b))
}

fit_one_curve <- function(d) {
  d <- d |> dplyr::filter(base::is.finite(VPD), base::is.finite(response))
  if (base::nrow(d) < 6L || dplyr::n_distinct(d$VPD) < 4L)
    base::stop("Insufficient data for the AN logistic curve.")
  x <- d$VPD; y <- d$response
  x_min <- base::min(x); x_max <- base::max(x); x_range <- x_max - x_min
  y_max <- base::max(y)
  if (!base::all(base::is.finite(c(x_min, x_max, x_range, y_max))) || x_range <= 0)
    base::stop("Invalid data range for AN curve.")
  start_grid <- base::expand.grid(
    a = c(base::max(y_max, 1e-4), base::max(y_max, 1e-4) * 1.25,
          base::max(y_max, 1e-4) * 1.75, base::max(y_max, 1e-4) * 2),
    x0 = base::as.numeric(stats::quantile(x,
                                          probs = c(0.25, 0.35, 0.50, 0.65, 0.75), names = FALSE)),
    b = c(base::max(x_range / 20, 0.01), base::max(x_range / 10, 0.02),
          base::max(x_range / 5, 0.05), base::max(x_range / 2, 0.10)),
    KEEP.OUT.ATTRS = FALSE
  )
  best_fit <- NULL; best_RSS <- Inf
  for (jj in base::seq_len(base::nrow(start_grid))) {
    z <- base::tryCatch(minpack.lm::nlsLM(
      response ~ logistic3(VPD, a, x0, b), data = d,
      start = base::as.list(start_grid[jj, ]),
      lower = c(a = 0, x0 = x_min, b = 1e-5),
      upper = c(a = Inf, x0 = x_max, b = Inf),
      control = minpack.lm::nls.lm.control(maxiter = 1000,
                                           ftol = 1e-10, ptol = 1e-10)
    ), error = function(e) NULL)
    if (!base::is.null(z)) {
      rss <- base::sum(stats::residuals(z)^2)
      if (base::is.finite(rss) && rss < best_RSS) {
        best_fit <- z; best_RSS <- rss
      }
    }
  }
  if (base::is.null(best_fit)) base::stop("All AN logistic fits failed.")
  list(fit = best_fit, RSS = best_RSS,
       coefficients = stats::coef(best_fit), n_obs = base::nrow(d))
}


# 6. POOLED AND SPECIES-SPECIFIC MODELS

fit_pooled <- function(d) {
  if (base::nrow(d) == 0L) base::stop("fit_pooled() received empty data.")
  fit_one_curve(d)
}
fit_separate <- function(d) {
  d <- base::droplevels(d)
  species_data <- base::split(d, d$species, drop = TRUE)
  group_sizes <- base::vapply(species_data, base::nrow, integer(1))
  if (base::length(species_data) < 2L || base::any(group_sizes == 0L))
    base::stop("Species-specific fitting requires >=2 non-empty groups.")
  fits <- base::lapply(species_data, fit_one_curve)
  list(fits = fits,
       RSS = base::sum(base::vapply(fits, function(z) z$RSS, numeric(1))))
}


# 7. OBSERVED GLOBAL TEST AND PARAMETERS

H0 <- fit_pooled(dat)
H1 <- fit_separate(dat)
Delta_RSS_obs <- H0$RSS - H1$RSS
base::cat("\nOBSERVED GLOBAL AN-VPD CURVE COMPARISON\n")
base::cat("RSS pooled = ", H0$RSS, "\n", sep = "")
base::cat("RSS species-specific = ", H1$RSS, "\n", sep = "")
base::cat("Observed Delta RSS = ", Delta_RSS_obs, "\n", sep = "")

pooled_parameters <- base::data.frame(
  Model = "Pooled", Species = "All species",
  base::as.list(H0$coefficients), RSS = H0$RSS,
  check.names = FALSE, stringsAsFactors = FALSE
)
species_parameters <- dplyr::bind_rows(base::lapply(base::names(H1$fits), function(sp) {
  z <- H1$fits[[sp]]
  base::data.frame(Model = "Species-specific", Species = sp,
                   base::as.list(z$coefficients), RSS = z$RSS,
                   check.names = FALSE, stringsAsFactors = FALSE)
}))
observed_parameters <- dplyr::bind_rows(pooled_parameters, species_parameters)
base::print(observed_parameters)
utils::write.csv(observed_parameters, parameter_output_file, row.names = FALSE)


# 8. GLOBAL WILD CLUSTER BOOTSTRAP

dat$fit_H0 <- base::as.numeric(stats::predict(H0$fit, newdata = dat))
dat$res_H0 <- dat$response - dat$fit_H0
dat$res_H0 <- dat$res_H0 - base::mean(dat$res_H0, na.rm = TRUE)
plants <- base::unique(base::as.character(dat$PlantCluster))
boot_delta <- base::rep(NA_real_, B)
boot_error <- base::rep("", B)

for (bb in base::seq_len(B)) {
  plant_sign <- base::sample(c(-1, 1), base::length(plants), replace = TRUE)
  base::names(plant_sign) <- plants
  db <- dat
  multiplier <- plant_sign[base::as.character(db$PlantCluster)]
  if (base::any(base::is.na(multiplier)))
    base::stop("A plant cluster did not receive a multiplier.")
  db$response <- db$fit_H0 + multiplier * db$res_H0
  boot_delta[bb] <- base::tryCatch({
    m0 <- fit_pooled(db); m1 <- fit_separate(db); m0$RSS - m1$RSS
  }, error = function(e) {
    boot_error[bb] <<- base::conditionMessage(e); NA_real_
  })
  if (bb %% 100L == 0L)
    base::cat("\rGlobal bootstrap replicate ", bb, " of ", B, sep = "")
}
base::cat("\n")
valid_global <- base::is.finite(boot_delta)
boot_delta_valid <- boot_delta[valid_global]
B_valid_global <- base::length(boot_delta_valid)
if (B_valid_global == 0L) {
  base::print(base::sort(base::table(boot_error[boot_error != ""]), decreasing = TRUE))
  base::stop("All global bootstrap replicates failed.")
}
p_global <- (1 + base::sum(boot_delta_valid >= Delta_RSS_obs)) /
  (1 + B_valid_global)
global_convergence_rate <- B_valid_global / B
base::cat("Valid global replicates = ", B_valid_global, "\n", sep = "")
base::cat("Global bootstrap P = ", p_global, "\n", sep = "")
if (global_convergence_rate < minimum_convergence_rate)
  base::warning("Global convergence rate below threshold: ",
                base::round(global_convergence_rate, 3))

global_error_table <- base::sort(
  base::table(boot_error[boot_error != ""]), decreasing = TRUE)
if (base::length(global_error_table) > 0L) {
  base::cat("\nMost frequent global-bootstrap errors:\n")
  base::print(utils::head(global_error_table, 10L))
}

utils::write.csv(
  base::data.frame(replicate = base::which(valid_global),
                   Delta_RSS = boot_delta_valid),
  global_bootstrap_output_file, row.names = FALSE
)
global_results <- base::data.frame(
  Variable = "AN", Comparison = "All species",
  RSS_pooled = H0$RSS, RSS_species_specific = H1$RSS,
  Delta_RSS = Delta_RSS_obs, P_bootstrap = p_global,
  B_requested = B, B_valid = B_valid_global,
  convergence_rate = global_convergence_rate,
  stringsAsFactors = FALSE
)
utils::write.csv(global_results, global_output_file, row.names = FALSE)


# 9. PAIRWISE CURVE-COMPARISON FUNCTION

pairwise_curve_test <- function(data, sp1, sp2, B = 1000L, seed = 123L) {
  base::set.seed(seed)
  d <- data |>
    dplyr::filter(species %in% c(sp1, sp2)) |>
    base::droplevels()
  species_counts <- dplyr::count(d, species, name = "n_obs")
  if (base::nrow(species_counts) != 2L)
    base::stop("A requested species is absent: ", sp1, " vs ", sp2)
  if (base::any(species_counts$n_obs < 6L))
    base::stop("Too few observations: ", sp1, " vs ", sp2)
  plant_counts <- d |>
    dplyr::distinct(species, PlantCluster) |>
    dplyr::count(species, name = "n_plants")
  if (base::any(plant_counts$n_plants < 2L))
    base::stop("Too few plant clusters: ", sp1, " vs ", sp2)
  
  m0 <- fit_pooled(d); m1 <- fit_separate(d)
  observed_delta <- m0$RSS - m1$RSS
  d$fit_H0 <- base::as.numeric(stats::predict(m0$fit, newdata = d))
  d$res_H0 <- d$response - d$fit_H0
  d$res_H0 <- d$res_H0 - base::mean(d$res_H0, na.rm = TRUE)
  pair_plants <- base::unique(base::as.character(d$PlantCluster))
  bootstrap_delta <- base::rep(NA_real_, B)
  bootstrap_error <- base::rep("", B)
  
  for (bb in base::seq_len(B)) {
    sign <- base::sample(c(-1, 1), base::length(pair_plants), replace = TRUE)
    base::names(sign) <- pair_plants
    db <- d
    multiplier <- sign[base::as.character(db$PlantCluster)]
    db$response <- db$fit_H0 + multiplier * db$res_H0
    bootstrap_delta[bb] <- base::tryCatch({
      q0 <- fit_pooled(db); q1 <- fit_separate(db); q0$RSS - q1$RSS
    }, error = function(e) {
      bootstrap_error[bb] <<- base::conditionMessage(e); NA_real_
    })
  }
  
  valid <- base::is.finite(bootstrap_delta)
  boot_valid <- bootstrap_delta[valid]
  B_valid <- base::length(boot_valid)
  if (B_valid == 0L) base::stop("All replicates failed: ", sp1, " vs ", sp2)
  p <- (1 + base::sum(boot_valid >= observed_delta)) / (1 + B_valid)
  convergence_rate <- B_valid / B
  error_table <- base::sort(base::table(bootstrap_error[bootstrap_error != ""]),
                            decreasing = TRUE)
  common_error <- if (base::length(error_table) > 0L)
    base::names(error_table)[1] else NA_character_
  base::data.frame(
    Variable = "AN", species1 = sp1, species2 = sp2,
    RSS_pooled = m0$RSS, RSS_separate = m1$RSS,
    Delta_RSS = observed_delta, P_bootstrap = p,
    B_requested = B, B_valid = B_valid,
    convergence_rate = convergence_rate,
    most_common_error = common_error,
    stringsAsFactors = FALSE
  )
}


# 10. RUN PAIRWISE COMPARISONS AND ADJUST P-VALUES

pairwise_results <- dplyr::bind_rows(
  pairwise_curve_test(dat, "faginea", "rotundifolia", B, 101L),
  pairwise_curve_test(dat, "faginea", "coccifera", B, 102L),
  pairwise_curve_test(dat, "rotundifolia", "coccifera", B, 103L)
)
utils::write.csv(pairwise_results, pairwise_output_file, row.names = FALSE)
pairwise_results_adjusted <- pairwise_results |>
  dplyr::mutate(
    P_Holm = stats::p.adjust(P_bootstrap, method = "holm"),
    P_FDR = stats::p.adjust(P_bootstrap, method = "BH")
  )
base::print(pairwise_results_adjusted)
utils::write.csv(pairwise_results_adjusted,
                 adjusted_pairwise_output_file, row.names = FALSE)


# 11. PREDICTIONS AND OPTIONAL FIGURE

prediction_grid <- base::data.frame(
  VPD = base::seq(base::min(dat$VPD, na.rm = TRUE),
                  base::max(dat$VPD, na.rm = TRUE), length.out = 300L)
)
pooled_prediction <- prediction_grid |>
  dplyr::mutate(species = "Pooled",
                response_predicted = base::as.numeric(
                  stats::predict(H0$fit, newdata = prediction_grid)))
species_predictions <- dplyr::bind_rows(base::lapply(base::names(H1$fits), function(sp) {
  prediction_grid |>
    dplyr::mutate(species = sp,
                  response_predicted = base::as.numeric(
                    stats::predict(H1$fits[[sp]]$fit,
                                   newdata = prediction_grid)))
}))
curve_predictions <- dplyr::bind_rows(pooled_prediction, species_predictions)
utils::write.csv(curve_predictions, prediction_output_file, row.names = FALSE)

if (base::requireNamespace("ggplot2", quietly = TRUE)) {
  species_colours <- c(faginea = "#D55E00", rotundifolia = "#0072B2",
                       coccifera = "#009E73", Pooled = "#000000")
  plot_data <- dat |>
    dplyr::mutate(species = base::as.character(species))
  curve_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = VPD, y = response, colour = species)
  ) +
    ggplot2::geom_point(alpha = 0.25, size = 1) +
    ggplot2::geom_line(
      data = curve_predictions,
      ggplot2::aes(x = VPD, y = response_predicted,
                   colour = species, linetype = species),
      linewidth = 1, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = species_colours,
      breaks = c("faginea", "rotundifolia", "coccifera", "Pooled"),
      labels = c("Q. faginea", "Q. rotundifolia", "Q. coccifera", "Pooled")
    ) +
    ggplot2::scale_linetype_manual(
      values = c(faginea = "solid", rotundifolia = "solid",
                 coccifera = "solid", Pooled = "dashed"),
      breaks = c("faginea", "rotundifolia", "coccifera", "Pooled"),
      labels = c("Q. faginea", "Q. rotundifolia", "Q. coccifera", "Pooled")
    ) +
    ggplot2::labs(x = expression(VPD[leaf]), y = "Net assimilation (AN)",
                  colour = "Curve", linetype = "Curve") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
  base::print(curve_plot)
  ggplot2::ggsave(figure_output_file, curve_plot, device = "tiff",
                  width = 180, height = 140, units = "mm", dpi = 600,
                  compression = "lzw")
}


# 12. FINAL SUMMARY AND SESSION INFORMATION

base::cat("\nGLOBAL RESULT\n"); base::print(global_results)
base::cat("\nPAIRWISE RESULTS WITH HOLM CORRECTION\n")
base::print(dplyr::select(pairwise_results_adjusted,
                          species1, species2, Delta_RSS,
                          P_bootstrap, P_Holm, B_valid, convergence_rate))
base::cat("\nOBSERVED logistic3 PARAMETERS\n"); base::print(observed_parameters)
utils::capture.output(utils::sessionInfo(), file = session_output_file)
base::cat("\nAnalysis completed for AN.\n")


#### E-VPD CURVE COMPARISONS ####
# FOUR-PARAMETER GAUSSIAN PEAK CURVES WITH PLANT-LEVEL WILD CLUSTER BOOTSTRAP

# Null: all species share one common curve.
# Alternative: each species has its own curve.
# Test statistic: Delta RSS = RSS_pooled - RSS_species-specific.
# One Rademacher multiplier (-1 or +1) is assigned to all residuals from
# the same plant in each bootstrap replicate.


base::set.seed(123)
B <- 1000L                    # use 100 while testing; >=1000 for final analysis
input_file <- "dataset.csv"
selected_WPgroup <- "(0, -0.1 MPa)"
PAR_threshold <- 500
minimum_convergence_rate <- 0.90

# Variable-specific output files
global_output_file <- "global_E_VPD_curve_test.csv"
pairwise_output_file <- "pairwise_E_VPD_curve_tests.csv"
adjusted_pairwise_output_file <- "pairwise_E_VPD_curve_tests_adjusted.csv"
parameter_output_file <- "observed_E_gauss_peak_curve_parameters.csv"
global_bootstrap_output_file <- "global_E_bootstrap_distribution.csv"
prediction_output_file <- "fitted_E_VPD_curves.csv"
figure_output_file <- "fitted_E_VPD_curves.tiff"
session_output_file <- "sessionInfo_E_VPD_bootstrap.txt"


# 1. CHECK REQUIRED PACKAGES

required_packages <- c("dplyr", "minpack.lm")
missing_packages <- required_packages[
  !base::vapply(required_packages, base::requireNamespace,
                quietly = TRUE, FUN.VALUE = logical(1))
]
if (base::length(missing_packages) > 0L) {
  base::stop("Install before running: ",
             base::paste(missing_packages, collapse = ", "))
}


# 2. READ AND CHECK DATA

if (!base::file.exists(input_file)) {
  base::stop("Input file not found: ", input_file,
             "\nWorking directory: ", base::getwd())
}
raw <- utils::read.csv(input_file, header = TRUE,
                       stringsAsFactors = FALSE,
                       na.strings = c("", "NA", "NaN"))
required_columns <- c("Sp", "IDplant", "VPDleaf", "E",
                      "PARc", "WPgroup")
missing_columns <- base::setdiff(required_columns, base::names(raw))
if (base::length(missing_columns) > 0L) {
  base::stop("Missing columns: ", base::paste(missing_columns, collapse = ", "))
}

base::cat("\nAVAILABLE WPgroup VALUES\n")
WPgroup_counts <- dplyr::count(raw, WPgroup, sort = TRUE,
                               name = "n_observations")
base::print(WPgroup_counts)
if (!selected_WPgroup %in% raw$WPgroup) {
  base::stop("Requested WPgroup not found: ", selected_WPgroup)
}


# 3. DATA PREPARATION

dat <- raw |>
  dplyr::transmute(
    species = Sp,
    IDplant = IDplant,
    VPD = VPDleaf,
    response = E,
    PARc = PARc,
    WPgroup = WPgroup
  ) |>
  dplyr::filter(
    PARc > PAR_threshold,
    WPgroup == selected_WPgroup,
    base::is.finite(VPD),
    base::is.finite(response),
    response >= 0,
    !base::is.na(species),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    species = dplyr::recode(
      base::as.character(species),
      fag = "faginea", ile = "rotundifolia", coc = "coccifera",
      .default = base::as.character(species)
    ),
    species = base::factor(species,
      levels = c("faginea", "rotundifolia", "coccifera")),
    PlantCluster = base::interaction(species, IDplant,
                                     drop = TRUE, sep = "_")
  ) |>
  dplyr::filter(!base::is.na(species)) |>
  base::droplevels()


# 4. VALIDATE FILTERED DATA

if (base::nrow(dat) == 0L) base::stop("No observations remain after filtering.")
data_summary <- dat |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    n_plants = dplyr::n_distinct(PlantCluster),
    n_obs = dplyr::n(),
    n_VPD_values = dplyr::n_distinct(VPD),
    min_VPD = base::min(VPD, na.rm = TRUE),
    max_VPD = base::max(VPD, na.rm = TRUE),
    min_response = base::min(response, na.rm = TRUE),
    max_response = base::max(response, na.rm = TRUE),
    .groups = "drop"
  )
base::cat("\nFILTERED-DATA SUMMARY: E\n")
base::print(data_summary)
base::print(dplyr::count(dat, species, IDplant, PlantCluster, name = "n_obs"))

expected_species <- c("faginea", "rotundifolia", "coccifera")
absent_species <- base::setdiff(expected_species,
                               base::unique(base::as.character(dat$species)))
if (base::length(absent_species) > 0L) {
  base::stop("Species absent after filtering: ",
             base::paste(absent_species, collapse = ", "))
}
if (base::any(data_summary$n_plants < 2L))
  base::stop("At least one species has fewer than two plant clusters.")
if (base::any(data_summary$n_obs < 8L))
  base::stop("At least one species has too few observations.")
if (base::any(data_summary$n_VPD_values < 5L))
  base::stop("At least one species has too few distinct VPD values.")


# 5. RESPONSE FUNCTION AND ROBUST SINGLE-CURVE FIT

gauss_peak <- function(x, A, mu, sigma, c) {
  c + A * base::exp(-0.5 * ((x - mu) / sigma)^2)
}

fit_one_curve <- function(d) {
  d <- d |> dplyr::filter(base::is.finite(VPD), base::is.finite(response))
  if (base::nrow(d) < 8L || dplyr::n_distinct(d$VPD) < 5L)
    base::stop("Insufficient data for the E Gaussian peak curve.")
  x <- d$VPD; y <- d$response
  x_min <- base::min(x); x_max <- base::max(x); x_range <- x_max - x_min
  y_min <- base::min(y); y_max <- base::max(y); amp <- base::max(y_max - y_min, 1e-4)
  mu_data <- x[base::which.max(y)][1]
  if (!base::all(base::is.finite(c(x_min, x_max, x_range, y_min, y_max, mu_data))) || x_range <= 0)
    base::stop("Invalid data range for E curve.")
  start_grid <- base::expand.grid(
    A = c(amp * 0.75, amp, amp * 1.5, amp * 2),
    mu = base::unique(c(mu_data,
      base::as.numeric(stats::quantile(x, c(0.35, 0.50, 0.65), names = FALSE)))),
    sigma = c(base::max(x_range / 10, 0.02), base::max(x_range / 5, 0.05),
              base::max(x_range / 3, 0.10), base::max(x_range / 2, 0.20)),
    c = c(base::max(0, y_min), base::max(0, stats::quantile(y, 0.05, names = FALSE))),
    KEEP.OUT.ATTRS = FALSE
  )
  best_fit <- NULL; best_RSS <- Inf
  for (jj in base::seq_len(base::nrow(start_grid))) {
    z <- base::tryCatch(minpack.lm::nlsLM(
      response ~ gauss_peak(VPD, A, mu, sigma, c), data = d,
      start = base::as.list(start_grid[jj, ]),
      lower = c(A = 0, mu = x_min, sigma = 1e-5, c = 0),
      upper = c(A = Inf, mu = x_max, sigma = Inf, c = Inf),
      control = minpack.lm::nls.lm.control(maxiter = 1000,
                                           ftol = 1e-10, ptol = 1e-10)
    ), error = function(e) NULL)
    if (!base::is.null(z)) {
      rss <- base::sum(stats::residuals(z)^2)
      if (base::is.finite(rss) && rss < best_RSS) {
        best_fit <- z; best_RSS <- rss
      }
    }
  }
  if (base::is.null(best_fit)) base::stop("All E Gaussian fits failed.")
  list(fit = best_fit, RSS = best_RSS,
       coefficients = stats::coef(best_fit), n_obs = base::nrow(d))
}


# 6. POOLED AND SPECIES-SPECIFIC MODELS

fit_pooled <- function(d) {
  if (base::nrow(d) == 0L) base::stop("fit_pooled() received empty data.")
  fit_one_curve(d)
}
fit_separate <- function(d) {
  d <- base::droplevels(d)
  species_data <- base::split(d, d$species, drop = TRUE)
  group_sizes <- base::vapply(species_data, base::nrow, integer(1))
  if (base::length(species_data) < 2L || base::any(group_sizes == 0L))
    base::stop("Species-specific fitting requires >=2 non-empty groups.")
  fits <- base::lapply(species_data, fit_one_curve)
  list(fits = fits,
       RSS = base::sum(base::vapply(fits, function(z) z$RSS, numeric(1))))
}


# 7. OBSERVED GLOBAL TEST AND PARAMETERS

H0 <- fit_pooled(dat)
H1 <- fit_separate(dat)
Delta_RSS_obs <- H0$RSS - H1$RSS
base::cat("\nOBSERVED GLOBAL E-VPD CURVE COMPARISON\n")
base::cat("RSS pooled = ", H0$RSS, "\n", sep = "")
base::cat("RSS species-specific = ", H1$RSS, "\n", sep = "")
base::cat("Observed Delta RSS = ", Delta_RSS_obs, "\n", sep = "")

pooled_parameters <- base::data.frame(
  Model = "Pooled", Species = "All species",
  base::as.list(H0$coefficients), RSS = H0$RSS,
  check.names = FALSE, stringsAsFactors = FALSE
)
species_parameters <- dplyr::bind_rows(base::lapply(base::names(H1$fits), function(sp) {
  z <- H1$fits[[sp]]
  base::data.frame(Model = "Species-specific", Species = sp,
                   base::as.list(z$coefficients), RSS = z$RSS,
                   check.names = FALSE, stringsAsFactors = FALSE)
}))
observed_parameters <- dplyr::bind_rows(pooled_parameters, species_parameters)
base::print(observed_parameters)
utils::write.csv(observed_parameters, parameter_output_file, row.names = FALSE)


# 8. GLOBAL WILD CLUSTER BOOTSTRAP

dat$fit_H0 <- base::as.numeric(stats::predict(H0$fit, newdata = dat))
dat$res_H0 <- dat$response - dat$fit_H0
dat$res_H0 <- dat$res_H0 - base::mean(dat$res_H0, na.rm = TRUE)
plants <- base::unique(base::as.character(dat$PlantCluster))
boot_delta <- base::rep(NA_real_, B)
boot_error <- base::rep("", B)

for (bb in base::seq_len(B)) {
  plant_sign <- base::sample(c(-1, 1), base::length(plants), replace = TRUE)
  base::names(plant_sign) <- plants
  db <- dat
  multiplier <- plant_sign[base::as.character(db$PlantCluster)]
  if (base::any(base::is.na(multiplier)))
    base::stop("A plant cluster did not receive a multiplier.")
  db$response <- db$fit_H0 + multiplier * db$res_H0
  boot_delta[bb] <- base::tryCatch({
    m0 <- fit_pooled(db); m1 <- fit_separate(db); m0$RSS - m1$RSS
  }, error = function(e) {
    boot_error[bb] <<- base::conditionMessage(e); NA_real_
  })
  if (bb %% 100L == 0L)
    base::cat("\rGlobal bootstrap replicate ", bb, " of ", B, sep = "")
}
base::cat("\n")
valid_global <- base::is.finite(boot_delta)
boot_delta_valid <- boot_delta[valid_global]
B_valid_global <- base::length(boot_delta_valid)
if (B_valid_global == 0L) {
  base::print(base::sort(base::table(boot_error[boot_error != ""]), decreasing = TRUE))
  base::stop("All global bootstrap replicates failed.")
}
p_global <- (1 + base::sum(boot_delta_valid >= Delta_RSS_obs)) /
            (1 + B_valid_global)
global_convergence_rate <- B_valid_global / B
base::cat("Valid global replicates = ", B_valid_global, "\n", sep = "")
base::cat("Global bootstrap P = ", p_global, "\n", sep = "")
if (global_convergence_rate < minimum_convergence_rate)
  base::warning("Global convergence rate below threshold: ",
                base::round(global_convergence_rate, 3))

global_error_table <- base::sort(
  base::table(boot_error[boot_error != ""]), decreasing = TRUE)
if (base::length(global_error_table) > 0L) {
  base::cat("\nMost frequent global-bootstrap errors:\n")
  base::print(utils::head(global_error_table, 10L))
}

utils::write.csv(
  base::data.frame(replicate = base::which(valid_global),
                   Delta_RSS = boot_delta_valid),
  global_bootstrap_output_file, row.names = FALSE
)
global_results <- base::data.frame(
  Variable = "E", Comparison = "All species",
  RSS_pooled = H0$RSS, RSS_species_specific = H1$RSS,
  Delta_RSS = Delta_RSS_obs, P_bootstrap = p_global,
  B_requested = B, B_valid = B_valid_global,
  convergence_rate = global_convergence_rate,
  stringsAsFactors = FALSE
)
utils::write.csv(global_results, global_output_file, row.names = FALSE)


# 9. PAIRWISE CURVE-COMPARISON FUNCTION

pairwise_curve_test <- function(data, sp1, sp2, B = 1000L, seed = 123L) {
  base::set.seed(seed)
  d <- data |>
    dplyr::filter(species %in% c(sp1, sp2)) |>
    base::droplevels()
  species_counts <- dplyr::count(d, species, name = "n_obs")
  if (base::nrow(species_counts) != 2L)
    base::stop("A requested species is absent: ", sp1, " vs ", sp2)
  if (base::any(species_counts$n_obs < 8L))
    base::stop("Too few observations: ", sp1, " vs ", sp2)
  plant_counts <- d |>
    dplyr::distinct(species, PlantCluster) |>
    dplyr::count(species, name = "n_plants")
  if (base::any(plant_counts$n_plants < 2L))
    base::stop("Too few plant clusters: ", sp1, " vs ", sp2)

  m0 <- fit_pooled(d); m1 <- fit_separate(d)
  observed_delta <- m0$RSS - m1$RSS
  d$fit_H0 <- base::as.numeric(stats::predict(m0$fit, newdata = d))
  d$res_H0 <- d$response - d$fit_H0
  d$res_H0 <- d$res_H0 - base::mean(d$res_H0, na.rm = TRUE)
  pair_plants <- base::unique(base::as.character(d$PlantCluster))
  bootstrap_delta <- base::rep(NA_real_, B)
  bootstrap_error <- base::rep("", B)

  for (bb in base::seq_len(B)) {
    sign <- base::sample(c(-1, 1), base::length(pair_plants), replace = TRUE)
    base::names(sign) <- pair_plants
    db <- d
    multiplier <- sign[base::as.character(db$PlantCluster)]
    db$response <- db$fit_H0 + multiplier * db$res_H0
    bootstrap_delta[bb] <- base::tryCatch({
      q0 <- fit_pooled(db); q1 <- fit_separate(db); q0$RSS - q1$RSS
    }, error = function(e) {
      bootstrap_error[bb] <<- base::conditionMessage(e); NA_real_
    })
  }

  valid <- base::is.finite(bootstrap_delta)
  boot_valid <- bootstrap_delta[valid]
  B_valid <- base::length(boot_valid)
  if (B_valid == 0L) base::stop("All replicates failed: ", sp1, " vs ", sp2)
  p <- (1 + base::sum(boot_valid >= observed_delta)) / (1 + B_valid)
  convergence_rate <- B_valid / B
  error_table <- base::sort(base::table(bootstrap_error[bootstrap_error != ""]),
                            decreasing = TRUE)
  common_error <- if (base::length(error_table) > 0L)
    base::names(error_table)[1] else NA_character_
  base::data.frame(
    Variable = "E", species1 = sp1, species2 = sp2,
    RSS_pooled = m0$RSS, RSS_separate = m1$RSS,
    Delta_RSS = observed_delta, P_bootstrap = p,
    B_requested = B, B_valid = B_valid,
    convergence_rate = convergence_rate,
    most_common_error = common_error,
    stringsAsFactors = FALSE
  )
}


# 10. RUN PAIRWISE COMPARISONS AND ADJUST P-VALUES

pairwise_results <- dplyr::bind_rows(
  pairwise_curve_test(dat, "faginea", "rotundifolia", B, 101L),
  pairwise_curve_test(dat, "faginea", "coccifera", B, 102L),
  pairwise_curve_test(dat, "rotundifolia", "coccifera", B, 103L)
)
utils::write.csv(pairwise_results, pairwise_output_file, row.names = FALSE)
pairwise_results_adjusted <- pairwise_results |>
  dplyr::mutate(
    P_Holm = stats::p.adjust(P_bootstrap, method = "holm"),
    P_FDR = stats::p.adjust(P_bootstrap, method = "BH")
  )
base::print(pairwise_results_adjusted)
utils::write.csv(pairwise_results_adjusted,
                 adjusted_pairwise_output_file, row.names = FALSE)


# 11. PREDICTIONS AND OPTIONAL FIGURE

prediction_grid <- base::data.frame(
  VPD = base::seq(base::min(dat$VPD, na.rm = TRUE),
                  base::max(dat$VPD, na.rm = TRUE), length.out = 300L)
)
pooled_prediction <- prediction_grid |>
  dplyr::mutate(species = "Pooled",
                response_predicted = base::as.numeric(
                  stats::predict(H0$fit, newdata = prediction_grid)))
species_predictions <- dplyr::bind_rows(base::lapply(base::names(H1$fits), function(sp) {
  prediction_grid |>
    dplyr::mutate(species = sp,
                  response_predicted = base::as.numeric(
                    stats::predict(H1$fits[[sp]]$fit,
                                   newdata = prediction_grid)))
}))
curve_predictions <- dplyr::bind_rows(pooled_prediction, species_predictions)
utils::write.csv(curve_predictions, prediction_output_file, row.names = FALSE)

if (base::requireNamespace("ggplot2", quietly = TRUE)) {
  species_colours <- c(faginea = "#D55E00", rotundifolia = "#0072B2",
                       coccifera = "#009E73", Pooled = "#000000")
  plot_data <- dat |>
    dplyr::mutate(species = base::as.character(species))
  curve_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = VPD, y = response, colour = species)
  ) +
    ggplot2::geom_point(alpha = 0.25, size = 1) +
    ggplot2::geom_line(
      data = curve_predictions,
      ggplot2::aes(x = VPD, y = response_predicted,
                   colour = species, linetype = species),
      linewidth = 1, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = species_colours,
      breaks = c("faginea", "rotundifolia", "coccifera", "Pooled"),
      labels = c("Q. faginea", "Q. rotundifolia", "Q. coccifera", "Pooled")
    ) +
    ggplot2::scale_linetype_manual(
      values = c(faginea = "solid", rotundifolia = "solid",
                 coccifera = "solid", Pooled = "dashed"),
      breaks = c("faginea", "rotundifolia", "coccifera", "Pooled"),
      labels = c("Q. faginea", "Q. rotundifolia", "Q. coccifera", "Pooled")
    ) +
    ggplot2::labs(x = expression(VPD[leaf]), y = "Transpiration (E)",
                  colour = "Curve", linetype = "Curve") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
  base::print(curve_plot)
  ggplot2::ggsave(figure_output_file, curve_plot, device = "tiff",
                  width = 180, height = 140, units = "mm", dpi = 600,
                  compression = "lzw")
}


# 12. FINAL SUMMARY AND SESSION INFORMATION

base::cat("\nGLOBAL RESULT\n"); base::print(global_results)
base::cat("\nPAIRWISE RESULTS WITH HOLM CORRECTION\n")
base::print(dplyr::select(pairwise_results_adjusted,
                          species1, species2, Delta_RSS,
                          P_bootstrap, P_Holm, B_valid, convergence_rate))
base::cat("\nOBSERVED gauss_peak PARAMETERS\n"); base::print(observed_parameters)
utils::capture.output(utils::sessionInfo(), file = session_output_file)
base::cat("\nAnalysis completed for E.\n")

#### WUE-VPD CURVE COMPARISONS ####
# THREE-PARAMETER EXPONENTIAL-OFFSET CURVES WITH PLANT-LEVEL WILD CLUSTER BOOTSTRAP

# Null: all species share one common curve.
# Alternative: each species has its own curve.
# Test statistic: Delta RSS = RSS_pooled - RSS_species-specific.
# One Rademacher multiplier (-1 or +1) is assigned to all residuals from
# the same plant in each bootstrap replicate.


base::set.seed(123)
B <- 1000L                    # use 100 while testing; >=1000 for final analysis
input_file <- "dataset.csv"
selected_WPgroup <- "(0, -0.1 MPa)"
PAR_threshold <- 500
minimum_convergence_rate <- 0.90

# Variable-specific output files
global_output_file <- "global_WUE_VPD_curve_test.csv"
pairwise_output_file <- "pairwise_WUE_VPD_curve_tests.csv"
adjusted_pairwise_output_file <- "pairwise_WUE_VPD_curve_tests_adjusted.csv"
parameter_output_file <- "observed_WUE_exp_offset_curve_parameters.csv"
global_bootstrap_output_file <- "global_WUE_bootstrap_distribution.csv"
prediction_output_file <- "fitted_WUE_VPD_curves.csv"
figure_output_file <- "fitted_WUE_VPD_curves.tiff"
session_output_file <- "sessionInfo_WUE_VPD_bootstrap.txt"


# 1. CHECK REQUIRED PACKAGES

required_packages <- c("dplyr", "minpack.lm")
missing_packages <- required_packages[
  !base::vapply(required_packages, base::requireNamespace,
                quietly = TRUE, FUN.VALUE = logical(1))
]
if (base::length(missing_packages) > 0L) {
  base::stop("Install before running: ",
             base::paste(missing_packages, collapse = ", "))
}


# 2. READ AND CHECK DATA

if (!base::file.exists(input_file)) {
  base::stop("Input file not found: ", input_file,
             "\nWorking directory: ", base::getwd())
}
raw <- utils::read.csv(input_file, header = TRUE,
                       stringsAsFactors = FALSE,
                       na.strings = c("", "NA", "NaN"))
required_columns <- c("Sp", "IDplant", "VPDleaf", "WUE",
                      "PARc", "WPgroup")
missing_columns <- base::setdiff(required_columns, base::names(raw))
if (base::length(missing_columns) > 0L) {
  base::stop("Missing columns: ", base::paste(missing_columns, collapse = ", "))
}

base::cat("\nAVAILABLE WPgroup VALUES\n")
WPgroup_counts <- dplyr::count(raw, WPgroup, sort = TRUE,
                               name = "n_observations")
base::print(WPgroup_counts)
if (!selected_WPgroup %in% raw$WPgroup) {
  base::stop("Requested WPgroup not found: ", selected_WPgroup)
}


# 3. DATA PREPARATION

dat <- raw |>
  dplyr::transmute(
    species = Sp,
    IDplant = IDplant,
    VPD = VPDleaf,
    response = WUE,
    PARc = PARc,
    WPgroup = WPgroup
  ) |>
  dplyr::filter(
    PARc > PAR_threshold,
    WPgroup == selected_WPgroup,
    base::is.finite(VPD),
    base::is.finite(response),
    response > 0,
    !base::is.na(species),
    !base::is.na(IDplant)
  ) |>
  dplyr::mutate(
    species = dplyr::recode(
      base::as.character(species),
      fag = "faginea", ile = "rotundifolia", coc = "coccifera",
      .default = base::as.character(species)
    ),
    species = base::factor(species,
                           levels = c("faginea", "rotundifolia", "coccifera")),
    PlantCluster = base::interaction(species, IDplant,
                                     drop = TRUE, sep = "_")
  ) |>
  dplyr::filter(!base::is.na(species)) |>
  base::droplevels()


# 4. VALIDATE FILTERED DATA

if (base::nrow(dat) == 0L) base::stop("No observations remain after filtering.")
data_summary <- dat |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    n_plants = dplyr::n_distinct(PlantCluster),
    n_obs = dplyr::n(),
    n_VPD_values = dplyr::n_distinct(VPD),
    min_VPD = base::min(VPD, na.rm = TRUE),
    max_VPD = base::max(VPD, na.rm = TRUE),
    min_response = base::min(response, na.rm = TRUE),
    max_response = base::max(response, na.rm = TRUE),
    .groups = "drop"
  )
base::cat("\nFILTERED-DATA SUMMARY: WUE\n")
base::print(data_summary)
base::print(dplyr::count(dat, species, IDplant, PlantCluster, name = "n_obs"))

expected_species <- c("faginea", "rotundifolia", "coccifera")
absent_species <- base::setdiff(expected_species,
                                base::unique(base::as.character(dat$species)))
if (base::length(absent_species) > 0L) {
  base::stop("Species absent after filtering: ",
             base::paste(absent_species, collapse = ", "))
}
if (base::any(data_summary$n_plants < 2L))
  base::stop("At least one species has fewer than two plant clusters.")
if (base::any(data_summary$n_obs < 6L))
  base::stop("At least one species has too few observations.")
if (base::any(data_summary$n_VPD_values < 4L))
  base::stop("At least one species has too few distinct VPD values.")


# 5. RESPONSE FUNCTION AND ROBUST SINGLE-CURVE FIT

exp_offset <- function(x, a, b, c) {
  c + a * base::exp(b * x)
}

fit_one_curve <- function(d) {
  d <- d |> dplyr::filter(base::is.finite(VPD), base::is.finite(response))
  if (base::nrow(d) < 6L || dplyr::n_distinct(d$VPD) < 4L)
    base::stop("Insufficient data for the WUE exponential-offset curve.")
  x <- d$VPD; y <- d$response
  x_min <- base::min(x); x_max <- base::max(x); x_range <- x_max - x_min
  y_min <- base::min(y); y_max <- base::max(y); amp <- base::max(y_max - y_min, 1e-4)
  if (!base::all(base::is.finite(c(x_min, x_max, x_range, y_min, y_max))) || x_range <= 0)
    base::stop("Invalid data range for WUE curve.")
  # b is constrained below zero because the specified model is exponential decay.
  start_grid <- base::expand.grid(
    a = c(amp * 0.5, amp, amp * 1.5, amp * 2),
    b = c(-0.01, -0.025, -0.05, -0.10, -0.20),
    c = c(base::max(0, y_min), base::max(0, stats::quantile(y, 0.05, names = FALSE))),
    KEEP.OUT.ATTRS = FALSE
  )
  best_fit <- NULL; best_RSS <- Inf
  for (jj in base::seq_len(base::nrow(start_grid))) {
    z <- base::tryCatch(minpack.lm::nlsLM(
      response ~ exp_offset(VPD, a, b, c), data = d,
      start = base::as.list(start_grid[jj, ]),
      lower = c(a = 0, b = -Inf, c = 0),
      upper = c(a = Inf, b = -1e-8, c = Inf),
      control = minpack.lm::nls.lm.control(maxiter = 1000,
                                           ftol = 1e-10, ptol = 1e-10)
    ), error = function(e) NULL)
    if (!base::is.null(z)) {
      rss <- base::sum(stats::residuals(z)^2)
      if (base::is.finite(rss) && rss < best_RSS) {
        best_fit <- z; best_RSS <- rss
      }
    }
  }
  if (base::is.null(best_fit)) base::stop("All WUE exponential fits failed.")
  list(fit = best_fit, RSS = best_RSS,
       coefficients = stats::coef(best_fit), n_obs = base::nrow(d))
}


# 6. POOLED AND SPECIES-SPECIFIC MODELS

fit_pooled <- function(d) {
  if (base::nrow(d) == 0L) base::stop("fit_pooled() received empty data.")
  fit_one_curve(d)
}
fit_separate <- function(d) {
  d <- base::droplevels(d)
  species_data <- base::split(d, d$species, drop = TRUE)
  group_sizes <- base::vapply(species_data, base::nrow, integer(1))
  if (base::length(species_data) < 2L || base::any(group_sizes == 0L))
    base::stop("Species-specific fitting requires >=2 non-empty groups.")
  fits <- base::lapply(species_data, fit_one_curve)
  list(fits = fits,
       RSS = base::sum(base::vapply(fits, function(z) z$RSS, numeric(1))))
}


# 7. OBSERVED GLOBAL TEST AND PARAMETERS

H0 <- fit_pooled(dat)
H1 <- fit_separate(dat)
Delta_RSS_obs <- H0$RSS - H1$RSS
base::cat("\nOBSERVED GLOBAL WUE-VPD CURVE COMPARISON\n")
base::cat("RSS pooled = ", H0$RSS, "\n", sep = "")
base::cat("RSS species-specific = ", H1$RSS, "\n", sep = "")
base::cat("Observed Delta RSS = ", Delta_RSS_obs, "\n", sep = "")

pooled_parameters <- base::data.frame(
  Model = "Pooled", Species = "All species",
  base::as.list(H0$coefficients), RSS = H0$RSS,
  check.names = FALSE, stringsAsFactors = FALSE
)
species_parameters <- dplyr::bind_rows(base::lapply(base::names(H1$fits), function(sp) {
  z <- H1$fits[[sp]]
  base::data.frame(Model = "Species-specific", Species = sp,
                   base::as.list(z$coefficients), RSS = z$RSS,
                   check.names = FALSE, stringsAsFactors = FALSE)
}))
observed_parameters <- dplyr::bind_rows(pooled_parameters, species_parameters)
base::print(observed_parameters)
utils::write.csv(observed_parameters, parameter_output_file, row.names = FALSE)


# 8. GLOBAL WILD CLUSTER BOOTSTRAP

dat$fit_H0 <- base::as.numeric(stats::predict(H0$fit, newdata = dat))
dat$res_H0 <- dat$response - dat$fit_H0
dat$res_H0 <- dat$res_H0 - base::mean(dat$res_H0, na.rm = TRUE)
plants <- base::unique(base::as.character(dat$PlantCluster))
boot_delta <- base::rep(NA_real_, B)
boot_error <- base::rep("", B)

for (bb in base::seq_len(B)) {
  plant_sign <- base::sample(c(-1, 1), base::length(plants), replace = TRUE)
  base::names(plant_sign) <- plants
  db <- dat
  multiplier <- plant_sign[base::as.character(db$PlantCluster)]
  if (base::any(base::is.na(multiplier)))
    base::stop("A plant cluster did not receive a multiplier.")
  db$response <- db$fit_H0 + multiplier * db$res_H0
  boot_delta[bb] <- base::tryCatch({
    m0 <- fit_pooled(db); m1 <- fit_separate(db); m0$RSS - m1$RSS
  }, error = function(e) {
    boot_error[bb] <<- base::conditionMessage(e); NA_real_
  })
  if (bb %% 100L == 0L)
    base::cat("\rGlobal bootstrap replicate ", bb, " of ", B, sep = "")
}
base::cat("\n")
valid_global <- base::is.finite(boot_delta)
boot_delta_valid <- boot_delta[valid_global]
B_valid_global <- base::length(boot_delta_valid)
if (B_valid_global == 0L) {
  base::print(base::sort(base::table(boot_error[boot_error != ""]), decreasing = TRUE))
  base::stop("All global bootstrap replicates failed.")
}
p_global <- (1 + base::sum(boot_delta_valid >= Delta_RSS_obs)) /
  (1 + B_valid_global)
global_convergence_rate <- B_valid_global / B
base::cat("Valid global replicates = ", B_valid_global, "\n", sep = "")
base::cat("Global bootstrap P = ", p_global, "\n", sep = "")
if (global_convergence_rate < minimum_convergence_rate)
  base::warning("Global convergence rate below threshold: ",
                base::round(global_convergence_rate, 3))

global_error_table <- base::sort(
  base::table(boot_error[boot_error != ""]), decreasing = TRUE)
if (base::length(global_error_table) > 0L) {
  base::cat("\nMost frequent global-bootstrap errors:\n")
  base::print(utils::head(global_error_table, 10L))
}

utils::write.csv(
  base::data.frame(replicate = base::which(valid_global),
                   Delta_RSS = boot_delta_valid),
  global_bootstrap_output_file, row.names = FALSE
)
global_results <- base::data.frame(
  Variable = "WUE", Comparison = "All species",
  RSS_pooled = H0$RSS, RSS_species_specific = H1$RSS,
  Delta_RSS = Delta_RSS_obs, P_bootstrap = p_global,
  B_requested = B, B_valid = B_valid_global,
  convergence_rate = global_convergence_rate,
  stringsAsFactors = FALSE
)
utils::write.csv(global_results, global_output_file, row.names = FALSE)


# 9. PAIRWISE CURVE-COMPARISON FUNCTION

pairwise_curve_test <- function(data, sp1, sp2, B = 1000L, seed = 123L) {
  base::set.seed(seed)
  d <- data |>
    dplyr::filter(species %in% c(sp1, sp2)) |>
    base::droplevels()
  species_counts <- dplyr::count(d, species, name = "n_obs")
  if (base::nrow(species_counts) != 2L)
    base::stop("A requested species is absent: ", sp1, " vs ", sp2)
  if (base::any(species_counts$n_obs < 6L))
    base::stop("Too few observations: ", sp1, " vs ", sp2)
  plant_counts <- d |>
    dplyr::distinct(species, PlantCluster) |>
    dplyr::count(species, name = "n_plants")
  if (base::any(plant_counts$n_plants < 2L))
    base::stop("Too few plant clusters: ", sp1, " vs ", sp2)
  
  m0 <- fit_pooled(d); m1 <- fit_separate(d)
  observed_delta <- m0$RSS - m1$RSS
  d$fit_H0 <- base::as.numeric(stats::predict(m0$fit, newdata = d))
  d$res_H0 <- d$response - d$fit_H0
  d$res_H0 <- d$res_H0 - base::mean(d$res_H0, na.rm = TRUE)
  pair_plants <- base::unique(base::as.character(d$PlantCluster))
  bootstrap_delta <- base::rep(NA_real_, B)
  bootstrap_error <- base::rep("", B)
  
  for (bb in base::seq_len(B)) {
    sign <- base::sample(c(-1, 1), base::length(pair_plants), replace = TRUE)
    base::names(sign) <- pair_plants
    db <- d
    multiplier <- sign[base::as.character(db$PlantCluster)]
    db$response <- db$fit_H0 + multiplier * db$res_H0
    bootstrap_delta[bb] <- base::tryCatch({
      q0 <- fit_pooled(db); q1 <- fit_separate(db); q0$RSS - q1$RSS
    }, error = function(e) {
      bootstrap_error[bb] <<- base::conditionMessage(e); NA_real_
    })
  }
  
  valid <- base::is.finite(bootstrap_delta)
  boot_valid <- bootstrap_delta[valid]
  B_valid <- base::length(boot_valid)
  if (B_valid == 0L) base::stop("All replicates failed: ", sp1, " vs ", sp2)
  p <- (1 + base::sum(boot_valid >= observed_delta)) / (1 + B_valid)
  convergence_rate <- B_valid / B
  error_table <- base::sort(base::table(bootstrap_error[bootstrap_error != ""]),
                            decreasing = TRUE)
  common_error <- if (base::length(error_table) > 0L)
    base::names(error_table)[1] else NA_character_
  base::data.frame(
    Variable = "WUE", species1 = sp1, species2 = sp2,
    RSS_pooled = m0$RSS, RSS_separate = m1$RSS,
    Delta_RSS = observed_delta, P_bootstrap = p,
    B_requested = B, B_valid = B_valid,
    convergence_rate = convergence_rate,
    most_common_error = common_error,
    stringsAsFactors = FALSE
  )
}


# 10. RUN PAIRWISE COMPARISONS AND ADJUST P-VALUES

pairwise_results <- dplyr::bind_rows(
  pairwise_curve_test(dat, "faginea", "rotundifolia", B, 101L),
  pairwise_curve_test(dat, "faginea", "coccifera", B, 102L),
  pairwise_curve_test(dat, "rotundifolia", "coccifera", B, 103L)
)
utils::write.csv(pairwise_results, pairwise_output_file, row.names = FALSE)
pairwise_results_adjusted <- pairwise_results |>
  dplyr::mutate(
    P_Holm = stats::p.adjust(P_bootstrap, method = "holm"),
    P_FDR = stats::p.adjust(P_bootstrap, method = "BH")
  )
base::print(pairwise_results_adjusted)
utils::write.csv(pairwise_results_adjusted,
                 adjusted_pairwise_output_file, row.names = FALSE)


# 11. PREDICTIONS AND OPTIONAL FIGURE

prediction_grid <- base::data.frame(
  VPD = base::seq(base::min(dat$VPD, na.rm = TRUE),
                  base::max(dat$VPD, na.rm = TRUE), length.out = 300L)
)
pooled_prediction <- prediction_grid |>
  dplyr::mutate(species = "Pooled",
                response_predicted = base::as.numeric(
                  stats::predict(H0$fit, newdata = prediction_grid)))
species_predictions <- dplyr::bind_rows(base::lapply(base::names(H1$fits), function(sp) {
  prediction_grid |>
    dplyr::mutate(species = sp,
                  response_predicted = base::as.numeric(
                    stats::predict(H1$fits[[sp]]$fit,
                                   newdata = prediction_grid)))
}))
curve_predictions <- dplyr::bind_rows(pooled_prediction, species_predictions)
utils::write.csv(curve_predictions, prediction_output_file, row.names = FALSE)

if (base::requireNamespace("ggplot2", quietly = TRUE)) {
  species_colours <- c(faginea = "#D55E00", rotundifolia = "#0072B2",
                       coccifera = "#009E73", Pooled = "#000000")
  plot_data <- dat |>
    dplyr::mutate(species = base::as.character(species))
  curve_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = VPD, y = response, colour = species)
  ) +
    ggplot2::geom_point(alpha = 0.25, size = 1) +
    ggplot2::geom_line(
      data = curve_predictions,
      ggplot2::aes(x = VPD, y = response_predicted,
                   colour = species, linetype = species),
      linewidth = 1, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = species_colours,
      breaks = c("faginea", "rotundifolia", "coccifera", "Pooled"),
      labels = c("Q. faginea", "Q. rotundifolia", "Q. coccifera", "Pooled")
    ) +
    ggplot2::scale_linetype_manual(
      values = c(faginea = "solid", rotundifolia = "solid",
                 coccifera = "solid", Pooled = "dashed"),
      breaks = c("faginea", "rotundifolia", "coccifera", "Pooled"),
      labels = c("Q. faginea", "Q. rotundifolia", "Q. coccifera", "Pooled")
    ) +
    ggplot2::labs(x = expression(VPD[leaf]), y = "Water-use efficiency (WUE)",
                  colour = "Curve", linetype = "Curve") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
  base::print(curve_plot)
  ggplot2::ggsave(figure_output_file, curve_plot, device = "tiff",
                  width = 180, height = 140, units = "mm", dpi = 600,
                  compression = "lzw")
}


# 12. FINAL SUMMARY AND SESSION INFORMATION

base::cat("\nGLOBAL RESULT\n"); base::print(global_results)
base::cat("\nPAIRWISE RESULTS WITH HOLM CORRECTION\n")
base::print(dplyr::select(pairwise_results_adjusted,
                          species1, species2, Delta_RSS,
                          P_bootstrap, P_Holm, B_valid, convergence_rate))
base::cat("\nOBSERVED exp_offset PARAMETERS\n"); base::print(observed_parameters)
utils::capture.output(utils::sessionInfo(), file = session_output_file)
base::cat("\nAnalysis completed for WUE.\n")
