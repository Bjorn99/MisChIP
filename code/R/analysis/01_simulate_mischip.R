# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# MisChIP: Global Parameter Definitions
# Phase1, Task 2: Define all simulation parameters

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


library(tidyverse)


# >>> 1. SIMULATION PARAMETERS <<<


# Basic genome parameters

GENOME_SIZE <- 1e6                # 1Mb simulated genome
BIN_SIZE <- 200                   # 200 bp bins
N_BINS <- GENOME_SIZE / BIN_SIZE  # 5,000 bins total

# Reproducibility
RANDOM_SEED <- 42


# >>> 2. BIOLOGICAL SCENARIO PARAMETERS <<<


# Parameter lists for each protein type
CTCF_PARAMS <- list(
  name = "CTCF",
  type = "transcription_factor",

  # Peak characteristics
  n_peaks = 100,                 # ~1 peak per 10kb
  peak_width_mean = 400,         # Narrow peak (bp)
  peak_width_sd = 100,           # Variation in width
  enrichment_mean = 25,          # Strong enrichment (fold over background)
  enrichment_sd = 10,            # Variation in binding strength

  # Background characteristics
  background_lambda = 5,         # Average reads per bin in background
  dispersion = 0.3,              # Lower dispersion (Poisson-like)

  # Spatial limitations
  min_peak_distance = 800        # At least 2X peak width apart
)

H3K27ME3_PARAMS <- list(
  name = "H2K27me3",
  type = "histone_modification",

  # Peak characteristics (broad domains)
  n_peaks = 25,                  # Fewer domains
  peak_width_mean = 5000,        # Broad domains (bp)
  peak_width_sd = 2000,          # High variation in domain size
  enrichment_mean = 5,           # Weaker enrichment (spread out)
  enrichment_sd = 2,             # Less variation

  # Background characteristics
  background_lambda = 8,         # Higher background (antibody noise)
  dispersion = 0.5,              # Higher dispersion (more variance)

  # Spatial limitations
  min_peak_distance = 1000       # Domains well-separated
)


# >>> 3. TECHNICAL BIAS PARAMETERS <<<


TECH_PARAMS <- list(
  # GC bias parameters
  gc_bias_strength = 0.2,        # +/-20% effect on read counts
  gc_bias_peak = 0.5,            # Optimal GC content (50%)
  gc_bias_width = 0.2,           # Width of GC preference curve


  # Mappability parameters
  mappability_fraction_low = 0.05, # 5% of genome has low mappability
  mappability_penalty = 0.3,       # 70% reduction in reads
  mappability_region_size = 1000,  # Size of low-mappability regions

  # Global noise
  noise_inflation_factor = 1.1     # 10% technical noise
)


# >>> 4. ANALYSIS PARAMETERS <<<


ANALYSIS_PARAMS <- list(
  # Peak calling method parameters
  fold_enrichment_threshold = 4,   # 4-fold enrichment
  poisson_pvalue_threshold = 0.01,
  fdr_threshold = 0.05,            # Multiple testing: FDR < 5%

  # Local background window
  local_bg_window_small = 1000,    # For sharp peaks (5 bins each side)
  local_bg_window_large = 10000,   # For broad peaks (50 bins each side)

  # Evaluation metrics
  peak_calling_resolution = BIN_SIZE # Call peaks at bin resolution
)


# >>> 5. PARAMETER VALIDATION <<<


validate_parameters <- function(params) {
  # Check that peak width doesn't exceed genome size
  max_width <- params$peak_width_mean + 3 * params$peak_width_sd
  if (max_width > GENOME_SIZE / 10) {
    warning(sprintf("%s: Maximum peak width (~%.0f bp) is > 10%% of genome size",
                    params$name, max_width))
  }

  # Check that all peaks can fit
  total_peak_space <- params$n_peaks * params$min_peak_distance
  if (total_peak_space > GENOME_SIZE * 0.8) {
    stop(sprintf("%s: Cannot fit %d peaks with minimum distance %d bp in genome",
                 params$name, params$n_peaks, params$min_peak_distance))
  }

  # Check dispersion parameter
  if (params$dispersion <= 0 || params$dispersion > 2) {
    warning(sprintf("%s: Unusual dispersion parameter %.2f",
                    params$name, params$dispersion))
  }

  return(TRUE)

}

# Validate both parameter
validate_parameters(CTCF_PARAMS)
validate_parameters(H3K27ME3_PARAMS)


# >>> 6. HELPER FUNCTIONS FOR PARAMETER ACCESS <<<


# Function to get all parameters for a scenario
get_scenario_params <- function(scenario = c("CTCF", "H3K27me3")) {
  scenario <- match.arg(scenario)

  # Combine biological and technical parameters
  if (scenario == "CTCF") {
    params <- CTCF_PARAMS
  } else {
    params <- H3K27ME3_PARAMS
  }

# Add technical parameters
params$tech <- TECH_PARAMS

# Add global parameters
params$genome_size <- GENOME_SIZE
params$bin_size <- BIN_SIZE
params$n_bins <- N_BINS

return(params)

}

# Function to print parameter summary
print_parameter_summary <- function() {
  cat(">>> MisChIP Parameter Summary <<<\n\n")

  cat("Genome Parameters:\n")
  cat(sprintf("  Genome size: %.0f bp (%.1f kb)\n", GENOME_SIZE, GENOME_SIZE/1000))
  cat(sprintf("  Bin size: %d bp\n", BIN_SIZE))
  cat(sprintf("  Total bins: %.0f\n\n", N_BINS))

  scenarios <- list(CTCF_PARAMS, H3K27ME3_PARAMS)

  for (params in scenarios) {
    cat(sprintf("%s Parameters:\n", params$name))
    cat(sprintf("  Type: %s\n", params$type))
    cat(sprintf("  Number of peaks: %d\n", params$n_peaks))
    cat(sprintf("  Peak width: %.0f ± %.0f bp\n",
                params$peak_width_mean, params$peak_width_sd))
    cat(sprintf("  Enrichment: %.1f ± %.1f fold\n",
                params$enrichment_mean, params$enrichment_sd))
    cat(sprintf("  Background: λ = %.1f reads/bin\n", params$background_lambda))
    cat(sprintf("  Dispersion: %.2f\n\n", params$dispersion))
  }

  cat("Technical Biases:\n")
  cat(sprintf("  GC bias strength: ±%.0f%%\n", TECH_PARAMS$gc_bias_strength * 100))
  cat(sprintf("  Low mappability: %.0f%% of genome\n",
              TECH_PARAMS$mappability_fraction_low * 100))
}


# Print Summary
print_parameter_summary()

# >>> SAVE PARAMETERS FOR REPRODUCILITY <<<

save_parameters <- function(output_dir = "results/parameters") {
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create parameter list
  all_params <- list(
    date_created = Sys.Date(),
    genome = list(size = GENOME_SIZE, bin_size = BIN_SIZE, n_bins = N_BINS),
    ctcf = CTCF_PARAMS,
    h3k27me3 = H3K27ME3_PARAMS,
    technical = TECH_PARAMS,
    analysis = ANALYSIS_PARAMS,
    random_seed = RANDOM_SEED
  )


# Saving as YAML for reproducibility
yaml::write_yaml(all_params, file.path(output_dir, "simulation_parameters.yaml"))


# Saving as RDS for exact R object preservation
saveRDS(all_params, file.path(output_dir, "simulation_parameters.rds"))

cat("\nParameters saved to:", output_dir, "\n")
}

save_parameters()




# >>> PHASE 1, TASK 3: GENERATE BASE DISTRIBUTIONS <<<


# Function to generate background read counts
generate_background <- function(n_bins, lambda, dispersion) {
  # Generate background read counts using negative binomial distribution
  # @return Vector of background read counts

  # Negative binomial in R:
  # size = 1/dispersion (smaller size = more dispersion)
  # mu = lambda (mean count)


  size_param <- 1 / dispersion


  # Generate counts
  counts <- rnbinom(n = n_bins,
                   size = size_param,
                   mu = lambda)

  return(counts)
}


# Function to create base genome structure
create_base_genome <- function(params) {
  # Create the base genome structure with background reads


  # Create bin coordinates
  genome_df <- data.frame(
    bin_id = 1:params$n_bins,
    start = seq(0, params$genome_size - params$bin_size, by = params$bin_size),
    end = seq(params$bin_size, params$genome_size, by = params$bin_size),


    # Chromosome (only simulating one)
    chr = "chr1"
  )

  # center position for each bin
  genome_df$center <- (genome_df$start + genome_df$end) / 2


  # Generate background counts
  set.seed(RANDOM_SEED)
  genome_df$background_counts <- generate_background(
    n_bins <- params$n_bins,
    lambda = params$background_lambda,
    dispersion = params$dispersion
  )


  # Initialize columns
  genome_df$true_peak <- FALSE          # mark true peak locations
  genome_df$peak_id <- NA               # store peak IDs
  genome_df$enrichment <- 1             # Fold enrichment (1 = no enrichment)


  # Add metadata
  attr(genome_df, "params") <- params
  attr(genome_df, "created") <- Sys.time()

  return(genome_df)
}


  # >>> Generate base genomes for both scenarios <<<


  # CTCF scenario
  ctcf_params <- get_scenario_params("CTCF")
  ctcf_genome <- create_base_genome(ctcf_params)

  # H3K27me3 scenario
  h3k27me3_params <- get_scenario_params("H3K27me3")
  h3k27me3_genome <- create_base_genome(h3k27me3_params)



  # >>> Examining the distributions <<<

  # Function to summarize background distribution
  summarize_background <- function(genome_df, scenario_name) {
    counts <- genome_df$background_counts

    cat(sprintf("\n>>> %s Background Distribution <<<\n", scenario_name))
    cat(sprintf("Mean counts: %.2f\n", mean(counts)))
    cat(sprintf("Variance: %.2f\n", var(counts)))
    cat(sprintf("Variance/Mean ratio: %.2f (>1 indicates overdispersion)\n",
                var(counts) / mean(counts)))
    cat(sprintf("Zeros: %d bins (%.1f%%)\n",
                sum(counts == 0),
                100 * sum(counts == 0) / length(counts)))
    cat(sprintf("Max count: %d\n", max(counts)))

    # Quantiles
    quants <- quantile(counts, probs = c(0.5, 0.5, 0.75, 0.95, 0.99))
    cat("\nQuantiles:\n")
    print(quants)
  }

  summarize_background(ctcf_genome, "CTCF")
  summarize_background(h3k27me3_genome, "H3K27me3")


  # >>> Visualize background distributions <<<

  # Visualization
  library(ggplot2)

  # Function to plot background distribution
  plot_background_distribution <- function(genome_df, scenario_name, params) {
    counts <- genome_df$background_counts

    # Create theoretical distribution for integer values only
    x_range <- 0:max(counts)
    theoretical_probs <- dnbinom(x_range,
                                 size = 1/params$dispersion,
                                 mu = params$background_lambda)
    theory_df <- data.frame(x = x_range, y = theoretical_probs)


    # Histogram
    p <- ggplot(data.frame(counts = counts), aes(x = counts)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30,
                     fill = "lightblue",
                     color = "black",
                     alpha = 0.7) +

      # overlay negative binomial
      geom_line(data = theory_df,
                aes(x = x, y = y),
                color = "red",
                linewidth = 1.5) +



      labs(
        title = sprintf("%s Background Distribution", scenario_name),
        subtitle = sprintf("λ = %.1f, dispersion = %.2f",
                           params$background_lambda,
                           params$dispersion),
        x = "Read counts per bin",
        y = "density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12)
      )

    return(p)
  }


  # Create plots
  p1 <- plot_background_distribution(ctcf_genome, "CTCF", ctcf_params)
  p2 <- plot_background_distribution(h3k27me3_genome, "H3K27me3", h3k27me3_params)


  # Display plots
  print(p1)
  print(p2)

  # Save plots
  ggsave(filename = "results/figures/background_dist_ctcf.png",
         plot = p1, width = 8, height = 6, dpi = 300)
  ggsave(filename = "results/figures/background_dist_h3k27me3.png",
         plot = p2, width = 8, height = 6, dpi = 300)


  # Function to plot background along genome
  plot_background_spatial <- function(genome_df, scenario_name,
                                      window_start = 1, window_size = 500) {
    # Plot a window of the genome
    window_end <- min(window_start + window_size - 1, nrow(genome_df))
    plot_data <- genome_df[window_start:window_end, ]

    p <- ggplot(plot_data, aes(x = bin_id, y = background_counts)) +
      geom_line(color = "darkblue", alpha = 0.7) +
      geom_point(size = 0.8, alpha = 0.5) +

      # Add smooth trend
      geom_smooth(method = "loess", span = 0.1, se = TRUE,
                  color = "red", fill = "pink", alpha = 0.3) +

      labs(
        title = sprintf("%s Background Signal (bins %d-%d)",
                        scenario_name, window_start, window_end),
        x = "Genomic bin",
        y = "Read count",
        subtitle = "Background only - no peaks added yet"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14)
      )

    return(p)
  }

  # Plot first 500 bins
  p3 <- plot_background_spatial(ctcf_genome, "CTCF")
  p4 <- plot_background_spatial(h3k27me3_genome, "H3K27me3")

  print(p3)
  print(p4)


  # Function to check spatial autocorrelation
  check_autocorrelation <- function(genome_df, max_lag = 20) {
    counts <- genome_df$background_counts

    # Calculate autocorrelation
    acf_result <- acf(counts, lag.max = max_lag, plot = FALSE)

    # Extract values (excluding lag 0)
    lags <- 1:max_lag
    correlations <- acf_result$acf[2:(max_lag + 1)]

    # Create plot
    p <- ggplot(data.frame(lag = lags, correlation = correlations),
                aes(x = lag, y = correlation)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_segment(aes(xend = lag, yend = 0), color = "darkblue") +
      geom_point(size = 3, color = "darkblue") +

      # Significance bounds (approximate)
      geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(counts)),
                 linetype = "dotted", color = "red", alpha = 0.5) +

      labs(
        title = "Autocorrelation of Background Signal",
        subtitle = "Should be near zero (independent bins)",
        x = "Lag (bins)",
        y = "Autocorrelation"
      ) +
      theme_minimal()

    return(p)
  }

  # Check both scenarios
  p5 <- check_autocorrelation(ctcf_genome) +
    ggtitle("CTCF Background Autocorrelation")
  p6 <- check_autocorrelation(h3k27me3_genome) +
    ggtitle("H3K27me3 Background Autocorrelation")

  print(p5)
  print(p6)


  # >>> Save genome objects <<<

    # Create output directory
    output_dir <- "data/simulated"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save base genomes
  saveRDS(ctcf_genome, file.path(output_dir, "ctcf_base_genome.rds"))
  saveRDS(h3k27me3_genome, file.path(output_dir, "h3k27me3_base_genome.rds"))

  cat("\n Base genomes saved to:", output_dir, "\n")
