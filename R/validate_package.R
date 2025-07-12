#!/usr/bin/env Rscript
#' Validation Script: R vs Stata Comparison
#' 
#' This script runs both R and Stata implementations of mhtexp2 and compares results.
#' Supports both generated test data and user-provided datasets.

# Required packages loaded via DESCRIPTION

#' Find Stata executable automatically
#' @export
find_stata <- function() {
  # Common Stata locations on different systems (newest versions first)
  possible_paths <- c(
    # macOS locations - Stata 19 first
    "/Applications/Stata/Stata19.app/Contents/MacOS/stata-mp",
    "/Applications/Stata/Stata19.app/Contents/MacOS/stata-se", 
    "/Applications/Stata/Stata19.app/Contents/MacOS/stata-ic",
    "/Applications/StataNow/Stata19.app/Contents/MacOS/stata-mp",
    "/Applications/StataNow/Stata19.app/Contents/MacOS/stata-se",
    "/Applications/StataNow/Stata19.app/Contents/MacOS/stata-ic",
    # macOS Stata 18, 17, 16
    "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp",
    "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se", 
    "/Applications/Stata/StataIC.app/Contents/MacOS/stata-ic",
    "/Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp",
    "/Applications/StataNow/StataSE.app/Contents/MacOS/stata-se",
    "/Applications/StataNow/StataIC.app/Contents/MacOS/stata-ic",
    # Windows locations
    "C:/Program Files/Stata19/StataMP-64.exe",
    "C:/Program Files/Stata18/StataMP-64.exe",
    "C:/Program Files/Stata17/StataMP-64.exe",
    "C:/Program Files/Stata16/StataMP-64.exe", 
    "C:/Program Files/Stata15/StataMP-64.exe",
    # Linux locations
    "/usr/local/stata19/stata-mp",
    "/usr/local/stata18/stata-mp",
    "/usr/local/stata17/stata-mp",
    "/usr/local/stata16/stata-mp",
    "/usr/local/stata15/stata-mp"
  )
  
  # Check if stata is in PATH
  stata_in_path <- Sys.which("stata")
  if (stata_in_path != "") {
    return(stata_in_path)
  }
  
  # Check common installation paths
  for (path in possible_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }
  
  # If not found, return NULL
  return(NULL)
}

#' Main validation function
#' 
#' @export
#' @param data_source "generate" to create test data, or path to CSV file
#' @param n_obs Number of observations (if generating data)
#' @param combo Comparison type: "pairwise" or "treatmentcontrol" 
#' @param bootstrap Number of bootstrap replications
#' @param studentized Whether to studentize test statistics
#' @param transitivity_check Whether to apply transitivity corrections
#' @param seed Random seed for reproducibility (NULL for random)
#' @param stata_path Path to Stata executable
#' @param verbose Whether to print detailed output
validate_mhtexp2 <- function(data_source = "generate",
                            n_obs = 1000,
                            combo = "pairwise", 
                            bootstrap = 5000,
                            studentized = TRUE,
                            transitivity_check = TRUE,
                            seed = NULL,
                            stata_path = NULL,
                            verbose = TRUE) {
  
  if (verbose) cat("=== MHTEXP2 R vs STATA VALIDATION ===\n\n")
  
  # Auto-detect Stata if path not provided
  if (is.null(stata_path)) {
    stata_path <- find_stata()
    if (is.null(stata_path)) {
      stop("Stata not found. Please install Stata or provide stata_path parameter.")
    }
    if (verbose) cat("Found Stata at:", stata_path, "\n")
  }
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    if (verbose) cat("Random seed set to:", seed, "\n")
  }
  
  # 1. Load or generate data
  if (data_source == "generate") {
    if (verbose) cat("Generating test data with", n_obs, "observations...\n")
    test_data <- generate_test_data(n_obs)
  } else {
    if (verbose) cat("Loading data from:", data_source, "\n")
    test_data <- read.csv(data_source)
    n_obs <- nrow(test_data)
  }
  
  # Validate data structure
  required_cols <- c("y1", "y2", "treat", "x1", "x2", "subgroup")
  missing_cols <- setdiff(required_cols, colnames(test_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) {
    cat("Data summary:\n")
    cat("  Observations:", n_obs, "\n")
    cat("  Treatments:", paste(sort(unique(test_data$treat)), collapse = ", "), "\n")
    cat("  Subgroups:", paste(sort(unique(test_data$subgroup)), collapse = ", "), "\n")
    cat("  Y1 mean:", round(mean(test_data$y1), 4), "\n")
    cat("  Y2 mean:", round(mean(test_data$y2), 4), "\n\n")
  }
  
  # 2. Run R implementation
  if (verbose) cat("Running R implementation...\n")
  
  # Reset seed for R run
  if (!is.null(seed)) set.seed(seed)
  
  start_time <- Sys.time()
  r_results <- mhtexp2_r(
    Y = test_data[, c("y1", "y2")],
    treatment = test_data$treat,
    controls = test_data[, c("x1", "x2")],
    subgroup = test_data$subgroup,
    combo = combo,
    bootstrap = bootstrap,
    studentized = studentized,
    transitivity_check = transitivity_check
  )
  r_time <- as.numeric(Sys.time() - start_time)
  
  if (verbose) cat("R execution time:", round(r_time, 2), "seconds\n")
  
  # 3. Run Stata implementation
  if (verbose) cat("Running Stata implementation...\n")
  
  stata_results <- run_stata_mhtexp2(
    data = test_data,
    combo = combo,
    bootstrap = bootstrap,
    studentized = studentized,
    transitivity_check = transitivity_check,
    seed = seed,
    stata_path = stata_path,
    verbose = verbose
  )
  
  # 4. Compare results
  if (verbose) cat("Comparing results...\n")
  comparison <- compare_results(r_results$output, stata_results, verbose = verbose)
  
  # 5. Simple summary
  if (verbose) {
    cat("\nComparison complete:", nrow(comparison), "hypotheses tested\n")
    cat("Use View(result$comparison) to see the full comparison table\n")
  }
  
  return(list(
    comparison = comparison,
    r_results = r_results,
    stata_results = stata_results,
    data = test_data,
    params = list(
      combo = combo,
      bootstrap = bootstrap,
      studentized = studentized,
      transitivity_check = transitivity_check,
      seed = seed,
      n_obs = n_obs
    )
  ))
}

#' Generate test data for validation
generate_test_data <- function(n = 1000) {
  data.frame(
    treat = sample(0:2, n, replace = TRUE),
    subgroup = sample(1:2, n, replace = TRUE),
    x1 = rnorm(n, mean = 0, sd = 1),
    x2 = rnorm(n, mean = 0, sd = 1),
    y1 = rnorm(n, mean = 0, sd = 1) + 0.1 * rep(c(0, 0.5, 1), length.out = n),
    y2 = rnorm(n, mean = 0, sd = 1) + 0.1 * rep(c(0, 0.3, 0.8), length.out = n)
  )
}

#' Run Stata mhtexp2 and return results
run_stata_mhtexp2 <- function(data, combo, bootstrap, studentized, transitivity_check, 
                             seed = NULL, stata_path, verbose = TRUE) {
  
  # Write data to temporary file
  temp_data_file <- "temp_validation_data.csv"
  write.csv(data, temp_data_file, row.names = FALSE)
  
  # Create Stata script
  stata_script <- paste0(
    "clear all\n",
    "set more off\n",
    "adopath + \"", getwd(), "\"\n",  # Add current directory to ado path
    if (!is.null(seed)) paste0("set seed ", seed, "\n") else "",
    "import delimited \"", temp_data_file, "\"\n",
    "mhtexp2 y1 y2, treatment(treat) controls(x1 x2) subgroup(subgroup) ",
    "combo(", combo, ") bootstrap(", bootstrap, ") ",
    "studentized(", as.numeric(studentized), ") ",
    "transitivitycheck(", as.numeric(transitivity_check), ")\n",
    "mata:\n",
    "result_matrix = st_matrix(\"results\")\n",
    "fh = fopen(\"temp_stata_results.csv\", \"w\")\n",
    "fput(fh, \"outcome,subgroup,t1,t2,coefficient,Remark3_2,Thm3_1,Remark3_8,Bonf,Holm\")\n",
    "for (i = 1; i <= rows(result_matrix); i++) {\n",
    "  row_str = strofreal(result_matrix[i,1])\n",
    "  for (j = 2; j <= cols(result_matrix); j++) {\n",
    "    row_str = row_str + \",\" + strofreal(result_matrix[i,j])\n",
    "  }\n",
    "  fput(fh, row_str)\n",
    "}\n",
    "fclose(fh)\n",
    "end\n",
    "exit\n"
  )
  
  # Write and execute Stata script
  writeLines(stata_script, "temp_stata_script.do")
  
  start_time <- Sys.time()
  system_result <- system(paste(stata_path, "-b do temp_stata_script.do"), 
                         ignore.stdout = !verbose, ignore.stderr = !verbose)
  stata_time <- as.numeric(Sys.time() - start_time)
  
  if (verbose) cat("Stata execution time:", round(stata_time, 2), "seconds\n")
  
  # Read results
  if (file.exists("temp_stata_results.csv")) {
    stata_results <- read.csv("temp_stata_results.csv")
    
    # Add comparison column to match R output
    stata_results$comparison <- with(stata_results, 
      ifelse(t1 == 0 & t2 == 1, 1,
             ifelse(t1 == 0 & t2 == 2, 2,
                    ifelse(t1 == 1 & t2 == 2, 3, NA)))
    )
    
    # Clean up temporary files
    file.remove(c("temp_validation_data.csv", "temp_stata_script.do", "temp_stata_results.csv"))
    if (file.exists("temp_stata_script.log")) file.remove("temp_stata_script.log")
    
    return(stata_results)
  } else {
    stop("Stata execution failed - no results file generated")
  }
}

#' Compare R and Stata results
compare_results <- function(r_output, stata_output, verbose = TRUE) {
  
  # Match results by (outcome, subgroup, t1, t2)
  r_sorted <- r_output[order(r_output$outcome, r_output$subgroup, r_output$t1, r_output$t2), ]
  stata_sorted <- stata_output[order(stata_output$outcome, stata_output$subgroup, stata_output$t1, stata_output$t2), ]
  
  if (nrow(r_sorted) != nrow(stata_sorted)) {
    warning("Different number of hypotheses: R=", nrow(r_sorted), ", Stata=", nrow(stata_sorted))
  }
  
  # Create comparison table
  n_rows <- min(nrow(r_sorted), nrow(stata_sorted))
  comparison <- data.frame(
    hypothesis = 1:n_rows,
    outcome = r_sorted$outcome[1:n_rows],
    subgroup = r_sorted$subgroup[1:n_rows],
    t1 = r_sorted$t1[1:n_rows],
    t2 = r_sorted$t2[1:n_rows],
    
    # R results
    coef_R = r_sorted$coefficient[1:n_rows],
    Remark3_2_R = r_sorted$Remark3_2[1:n_rows],
    Thm3_1_R = r_sorted$Thm3_1[1:n_rows],
    Remark3_8_R = r_sorted$Remark3_8[1:n_rows],
    
    # Stata results
    coef_Stata = stata_sorted$coefficient[1:n_rows],
    Remark3_2_Stata = stata_sorted$Remark3_2[1:n_rows],
    Thm3_1_Stata = stata_sorted$Thm3_1[1:n_rows],
    Remark3_8_Stata = stata_sorted$Remark3_8[1:n_rows],
    
    # Differences
    diff_coef = r_sorted$coefficient[1:n_rows] - stata_sorted$coefficient[1:n_rows],
    diff_Remark3_2 = r_sorted$Remark3_2[1:n_rows] - stata_sorted$Remark3_2[1:n_rows],
    diff_Thm3_1 = r_sorted$Thm3_1[1:n_rows] - stata_sorted$Thm3_1[1:n_rows],
    diff_Remark3_8 = r_sorted$Remark3_8[1:n_rows] - stata_sorted$Remark3_8[1:n_rows]
  )
  
  
  return(comparison)
}

# Example usage:
# validate_mhtexp2()
# validate_mhtexp2(n_obs=500, bootstrap=500)
# validate_mhtexp2(data_source="mydata.csv", combo="treatmentcontrol")