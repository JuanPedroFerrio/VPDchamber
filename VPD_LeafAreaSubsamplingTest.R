library(dplyr)
library(rstudioapi)
library(ggplot2)

#### SET DIR ####
## Setting the working directory to the path of your current open file
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
# stores mainpath
mainpath<-getwd()

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
B <- 5000L

# Percentile interval defining leaves eligible for sampling.
#
# 0.05 and 0.95 retain the central 90% of the measured leaves.
lower_sampling_percentile <- 0.00
upper_sampling_percentile <- 1

# Smallest sample size to assess
minimum_sample_size <- 1L

# Accuracy target:
# select the smallest sample size for which the 95th percentile
# of absolute relative error is <= this value.
target_error_percent <- 5

# Equivalent probability target:
# at least 50% of random samples should remain within the error threshold.
target_probability <- 0.5


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
  faginea = "VPD_Leafarea_fag.csv",
  rotundifolia = "VPD_Leafarea_rot.csv",
  coccifera = "VPD_Leafarea_coc.csv"
)

# Name of the leaf-area column in all three files
area_column <- "LA_cm2"


# 2. TOTAL NUMBER OF LEAVES COUNTED ON EACH TREE

#
# Replace these illustrative values with the total number of leaves counted
# on each complete tree.
#
# The number of counted leaves does not have to equal the number of leaves
# whose individual areas were measured.
#


total_leaf_counts <- c(
  faginea = 157,
  rotundifolia = 182,
  coccifera = 169
)

# coc	186.2068966	251.2616974	361.3259669
# fag	43.90651085	72.48391974	122.8387097
# rot	147.099768	309.7284314	385.1626016


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
    
    # Summarise results for this sample size

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
      
      Lower_90_mean_area =
        stats::quantile(
          sample_means,
          probs = 0.05,
          names = FALSE
        ),
      
      Upper_90_mean_area =
        stats::quantile(
          sample_means,
          probs = 0.95,
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
#   Smallest n for which at least 75% of samples have an absolute error
#   no greater than target_error_percent.
#
# This is equivalent to:
#   P75 absolute error <= target_error_percent.
#


recommended_sample_size <-
  error_summary |>
  dplyr::group_by(
    Species
  ) |>
  dplyr::filter(
    P75_absolute_error_percent <=
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
    
    Median_absolute_error_percent,
    
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
          Median_absolute_error_percent <=
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
          
          Median_absolute_error_percent,
          
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
    ), linewidth = 0.8) +
  # ggplot2::geom_line(
  #   mapping = ggplot2::aes(
  #     y = P75_absolute_error_percent
  #   ),
  #   linewidth = 1,
  #   linetype = "dashed"
  # ) +
  ggplot2::geom_vline(
    xintercept =
      10,
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
      "vertical line: ",
      "N= 10"
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

