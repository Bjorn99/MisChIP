# MisChIP


## Project Description

This project simulates ChIP-seq data with known peak locations and test different peak calling methods to understand how they work. By knowing where the "true" peaks are in simulated data, I can measure how well different approaches perform at finding them.

I wanted to understand peak calling by doing it myself. Tools are useful but you don't really understand the trade-offs involved until you understand how they work. By simulating data and implementing basic peak calling methods, I can see exactly how different approaches succeed or fail.

### Objective

Create a biologically-informed ChIP-seq simulation framework to compare multiple peak calling strategies, understand the sensitivity-specificity trade-offs, explore how different biological scenarios affect peak detection accuracy.

### What This Project Does

1. **Simulates ChIP-seq data** using negative binomial distributions to generate realistic read counts
    - Adds biological features like GC bias and variable mappability
    - Creates two scenarios: sharp peaks (CTCF) and broad domains (H3K27me3)

2. **Implements three peak calling methods** of increasing complexity:
    - Simple fold-change threshold
    - Statistical significance testing
    - Local background estimation with multiple testing correction (Benjamini-Hochberg FDR)

3. **Compares their performance** by measuring:
    - How many true peaks each method finds (sensitivity)
    - How many false peaks each method calls (specificity)
    - How these metrics change with different parameter settings
    - ROC and Precision-Recall curves to visualize performance
    - F1 scores to balance sensitivity and precision

### Learning Resources
This project draws from:

    - MIT 7.91J: Foundations of Computational and Systems Biology - for understanding ChIP-seq analysis

    - Modern Statistics for Modern Biology (Holmes & Huber) - specifically the chapter on generative models for discrete data

    - Primary literature on peak calling algorithms


### What I Want to Achieve

    - Understand the statistical foundations of peak calling (Poisson processes, hypothesis testing)

    - Learn why peak calling is challenging (noise, biases, different peak, shapes)

    - Implement core algorithms to understand how they work

    - Develop intuition for choosing the right method for different scenarios.


## Project Phases

### Phase1: Data Simulation (Ongoing)

**Goal**: Create realistic ChIP-seq data with ground truth

- Generated 10 Mb genome divided into 100 bp bins
- Added Poisson-distributed background noise
- Simulated two scenarios:

    - CTCF: Sharp, narrow peaks (~500 bp)
    - H3K27me3: Broad domains (~5000 bp)


- Added realistic technical biases:

    - GC content bias (affects sequencing efficiency)
    - Mappability issues (repetitive regions)


- Exported in standard formats (BED, wiggle)

**Key Learning**: How ChIP-seq data is structured and what makes it noisy

### Phase 2: Peak Calling Implementation

**Goal:** Code different peak calling approaches

Planned implementations:

    1. Sliding Window - Simple fixed-window approach
    
    2. MACS2 Algorithm - Dynamic local background model
    
    3. SICER - Designed for broad mark
    
    4. Simple Bayesian - Prior/posterior approach

**Key Learning:** Different statistical approaches to the same problem

### Phase 3: Evaluation Framework

**Goal:** Systematically compare methods

    - Match called peaks to true peaks
    - Calculate metrics:
    
        - Precision/Recall/F1-score
        - Spatial accuracy

    - Parameter sensitivity analysis

**Key Learning:** How to evaluate computational methods objectively


### Phase 4: Analysis and Summary

**Goal:** Synthesize findings

    - Performance comparisons
    - Method recommendations
    - Document statistical concepts learned

**Key Learning:** When to use which method and why

## Repository Structure

```mischip-seq-analysis/
├── code/
│   ├── R/
│   │   ├── analysis/          # Main analysis scripts
│   │   ├── functions/         # Reusable functions
│   │   └── visualization/     # Plotting functions
│   └── scripts/              # Setup scripts
├── data/
│   ├── simulated/            # Generated ChIP-seq data
│   ├── processed/            # Intermediate files
│   └── raw/                  # (Currently unused)
├── results/
│   ├── figures/              # All plots
│   ├── parameters/           # Configuration files
│   ├── tables/               # Performance metrics
│   └── calls/                # Peak calls from each method
├── docs/                     # Documentation and notes
└── literature/               # References
```


## Getting Started

### Prerequisites

    - R(>=4.0)
    - Required packages: tidyverse, ggplot2, yaml, scales, gridExtra

Install packages:

```Rscript code/scripts/install_packages.R```

### Running the Analysis

```# Phase 1: Generate simulated data
Rscript code/R/analysis/01_simulate_mischip.R

# Phase 2-4: Coming soon...
```

## Current Status

- [x] Phase 1: Data simulation complete
- [ ] Phase 2: Peak calling methods (next)
- [ ] Phase 3: Evaluation framework
- [ ] Phase 4: Analysis and summary


## Expected Outcomes

By completing this project, I expect to:

1. **Understand the statistics** - Why Poisson? When Negative Binomial? How do we test for enrichment?
2. **Appreciate the challenges** - Why is peak calling hard? How do biases affect results?
3. **Compare approaches** - When does MACS2 outperform sliding window? Why use SICER for broad marks?
4. **Make informed choices** - Given a new ChIP-seq dataset, know which method to use

```
**Note**: This project focuses on understanding basic peak calling principles using simulated data. We use FDR (False Discovery Rate) using **Benjamini-Hochberg** correction but not IDR (Irreproducible Discovery Rate), which is used for comparing biological replicates to assess reproducibility.

- This is a learning project - code clarity is prioritized over computational efficiency
- Each script includes detailed comments explaining the "why" behind each step
- Results folder is gitignored - run scripts to generate your own results
```