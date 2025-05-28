# MisChIP


## Project Description

This project simulates ChIP-seq data with known peak locations and test different peak calling methods to understand how they work. By knowing where the "true" peaks are in simulated data, I can measure how well different approaches perform at finding them.

### Background

After studying ChIP-seq data analysis through MIT 7.91J (Foundations of Computational and systems biology) and learning about generative models for discrete data (Modern Statistics for Modern Biology), I wanted to combine these concepts in a hands-on project. The Course taught me how ChIP-seq works, while studying discrete probability distributions showed me how to model count-based sequencing data.

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

### Learning Goals

- How discrete probability distribution model sequencing data
- Why some peaks are easier to detect than others
- How background noise affects peak calling
- The importance of multiple testing correction in genomics
- Trade-offs between finding all peak vs. avoiding false positives

**Note**: This project focuses on understanding basic peak calling principles using simulated data. We use FDR (False Discovery Rate) using **Benjamini-Hochberg** correction but not IDR (Irreproducible Discovery Rate), which is used for comparing biological replicates to assess reproducibility.

This is a learning project to understand the statistical principles behind ChIP-seq analysis through hands-on implementation and experimentation.
