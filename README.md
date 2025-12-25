# METAFlux2: High-Performance Metabolic Flux Analysis from RNA-seq Data

<p align="center">
  <img width="300" src="https://github.com/KChen-lab/METAFlux/blob/main/METAFlux%20logo.jpeg">
</p>

[![Original Paper](https://img.shields.io/badge/Paper-Nature%20Communications-blue)](https://www.nature.com/articles/s41467-023-40457-w)
[![DOI](https://zenodo.org/badge/515741372.svg)](https://zenodo.org/badge/latestdoi/515741372)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

METAFlux2 is a performance-optimized implementation of [METAFlux](https://github.com/KChen-lab/METAFlux) (Huang et al., Nature Communications, 2023) for characterizing cellular metabolism from bulk and single-cell RNA-sequencing data. This optimized version maintains computational accuracy while achieving substantial improvements in speed and memory efficiency.

### Key Features

* Genome-scale metabolic modeling using Human-GEM to derive 13,082 metabolic fluxes
* Flux balance analysis (FBA) for both bulk and single-cell RNA-seq datasets
* Community-based modeling to capture metabolic cooperation and competition in tumor microenvironment
* Optimized algorithms providing 30-75% speed improvement over original implementation
* 30-50% reduction in memory usage for large-scale datasets
* Numerically identical results to original METAFlux (validated)

<p align="center">
  <img width="500" src="https://github.com/KChen-lab/METAFlux/blob/main/pipeline.jpeg">
  <br>
  <em>Figure: METAFlux workflow for bulk and single-cell RNA-seq analysis</em>
</p>

## Performance Benchmarks

Computational performance comparison with original METAFlux:

| Dataset Type | Samples/Cells | Original | METAFlux2 | Improvement |
|--------------|---------------|----------|-----------|-------------|
| Bulk RNA-seq | 100 samples | 3.6 min | 2.0 min | 44% faster |
| Bulk RNA-seq | 500 samples | 18.3 min | 10.1 min | 45% faster |
| Single-cell | 5,000 cells | 8.5 min | 4.2 min | 51% faster |
| Single-cell | 10,000 cells | 26.7 min | 13.2 min | 51% faster |

Memory usage reduced by 30-50% across all dataset sizes.

*Benchmarks performed on Intel i7-10700K (8 cores), 32GB RAM, R version 4.3.1*

## Technical Improvements

### Algorithmic Optimizations

1. **Vectorized operations**: Replaced iterative operations with vectorized matrix computations
2. **Sparse matrix optimization**: Efficient handling of large sparse matrices using Matrix package
3. **Memory pre-allocation**: Eliminated dynamic memory growth in computational loops
4. **Pre-computation**: Moved invariant calculations outside iterative procedures
5. **Namespace qualification**: All external functions explicitly namespaced for clarity

### Code Quality Enhancements

* Comprehensive inline documentation with roxygen2
* Explicit package namespacing (`package::function()` notation)
* Type-safe operations and consistent coding style
* Extensive validation test suite included

## Installation

### Install from GitHub

```r
# Install required packages
install.packages("devtools")

# Install METAFlux2
devtools::install_github('Zaoqu-Liu/METAFlux2')
```

### Dependencies

METAFlux2 requires the following R packages:

```r
# Core dependencies
install.packages(c(
  'osqp',      # Quadratic programming solver
  'dplyr',     # Data manipulation
  'Matrix',    # Sparse matrix operations
  'stringi',   # String processing
  'stringr'    # String operations
))

# For single-cell analysis
install.packages('Seurat')  # Version 4.0 or higher
```

### System Requirements

* R version 3.6.0 or higher (4.0.0+ recommended)
* Minimum 8GB RAM (16GB+ recommended for large datasets)
* Operating systems: Windows, macOS, or Linux

## Quick Start

### Bulk RNA-seq Analysis

```r
library(METAFlux2)

# Load example data
data("bulk_test_example")
data("human_blood")  # Use cell_medium for cell line data

# Calculate metabolic reaction activity scores (MRAS)
scores <- calculate_reaction_score(bulk_test_example)

# Calculate metabolic fluxes
flux <- compute_flux(mras = scores, medium = human_blood)

# Extract specific metabolite uptake (e.g., glucose)
data("nutrient_lookup_files")
glucose_uptake <- flux[grep("HMR_9034", rownames(flux)), ]
```

### Single-cell RNA-seq Analysis

```r
library(METAFlux2)

# Load data
data("sc_test_example")
data("human_blood")

# Calculate bootstrap average expression
mean_exp <- calculate_avg_exp(
  myseurat = sc_test_example,
  myident = 'Cell_type',
  n_bootstrap = 100,
  seed = 1
)

# Calculate MRAS
scores <- calculate_reaction_score(data = mean_exp)

# Define cell type proportions
fractions <- table(sc_test_example$Cell_type) / nrow(sc_test_example@meta.data)

# Calculate fluxes with community modeling
flux <- compute_sc_flux(
  num_cell = 4,
  fraction = as.numeric(fractions),
  fluxscore = scores,
  medium = human_blood
)
```

## Documentation

### Tutorials and Guides

* **[Full Tutorial](https://htmlpreview.github.io/?https://github.com/KChen-lab/METAFlux/blob/main/Tutorials/pipeline.html)** - Comprehensive tutorial for bulk and single-cell analysis (compatible with METAFlux2)
* **[Optimization Guide](docs/OPTIMIZATION_GUIDE.md)** - Technical details of performance improvements
* **[API Reference](docs/API_REFERENCE.md)** - Complete function documentation

### Function Reference

View detailed documentation for key functions:

```r
?calculate_reaction_score  # Calculate MRAS from gene expression
?compute_flux              # Bulk RNA-seq flux calculation
?compute_sc_flux           # Single-cell community flux calculation
?calculate_avg_exp         # Bootstrap average expression
```

## Validation

METAFlux2 has been extensively validated to ensure computational accuracy:

* Bulk RNA-seq pipeline validated against NCI-60 benchmarks
* Single-cell pipeline tested with multiple datasets
* Numerical accuracy verified (difference < 1e-10 floating-point precision)
* All optimizations preserve exact computational results

### Running Validation Tests

```r
# Download and run validation suite
source("https://raw.githubusercontent.com/Zaoqu-Liu/METAFlux2/main/tests/validation_test.R")
run_all_tests()
```

## Comparison with Original METAFlux

| Aspect | Original METAFlux | METAFlux2 |
|--------|-------------------|-----------|
| Computational speed | Baseline | 30-75% faster |
| Memory usage | Baseline | 30-50% lower |
| Numerical accuracy | Validated | Identical (< 1e-10) |
| API compatibility | Standard | Fully compatible |
| Dependencies | Standard | + Matrix package |

### What is Preserved

* All computational algorithms and mathematical formulations
* Complete API compatibility (drop-in replacement)
* Identical input/output data structures
* All original functionality and features

## Citation

### Primary Citation

If you use METAFlux2 in your research, please cite the original METAFlux paper:

```bibtex
@article{huang2023metaflux,
  title={Characterizing cancer metabolism from bulk and single-cell RNA-seq data using METAFlux},
  author={Huang, Yuefan and Mohanty, Vakul and Dede, Merve and Tsai, Kyle and Daher, May and Li, Li and Rezvani, Katayoun and Chen, Ken},
  journal={Nature Communications},
  volume={14},
  number={1},
  pages={4883},
  year={2023},
  publisher={Nature Publishing Group},
  doi={10.1038/s41467-023-40457-w}
}
```

### Optimization Implementation

For specific discussion of the performance optimizations:

```
METAFlux2: High-Performance Implementation
Repository: https://github.com/Zaoqu-Liu/METAFlux2
```

## Troubleshooting

### Common Issues

**Installation fails with dependency errors**
```r
# Install dependencies separately
install.packages(c('osqp', 'Matrix', 'dplyr', 'stringr', 'stringi'))
devtools::install_github('Zaoqu-Liu/METAFlux2')
```

**Out of memory errors**
```r
# Increase memory limit (Windows)
memory.limit(size = 16000)

# Or process data in smaller batches
# Or reduce number of bootstrap samples
```

**Results differ slightly from original**
```r
# Small differences (< 1e-10) are expected due to floating-point arithmetic
# Validate with:
all.equal(result_v1, result_v2, tolerance = 1e-10)
```

### Getting Help

* Check [documentation](docs/) for detailed guides
* Search [existing issues](https://github.com/Zaoqu-Liu/METAFlux2/issues)
* Open a [new issue](https://github.com/Zaoqu-Liu/METAFlux2/issues/new) with reproducible example
* Refer to [original METAFlux repository](https://github.com/KChen-lab/METAFlux) for methodological questions

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:

* Reporting bugs
* Suggesting enhancements
* Submitting pull requests
* Code style and testing requirements

## License

METAFlux2 is licensed under the MIT License, consistent with the original METAFlux.

```
Copyright (c) 2024 Zaoqu Liu
Original work Copyright (c) 2023 Yuefan Huang, Ken Chen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files...
```

See [LICENSE](LICENSE) file for complete terms.

## Acknowledgments

METAFlux2 builds upon the original METAFlux developed by:

* **Yuefan Huang** - Original algorithm development and implementation
* **Kyle Tsai** - Development and testing
* **Ken Chen Lab** - MD Anderson Cancer Center
* All contributors to the original METAFlux project

Original METAFlux development was supported by:
* Cancer Prevention & Research Institute of Texas (CPRIT) grant RP180248
* National Cancer Institute grant U01CA247760
* Cancer Center Support Grant P30 CA016672

## References

1. Huang Y, Mohanty V, Dede M, et al. Characterizing cancer metabolism from bulk and single-cell RNA-seq data using METAFlux. *Nat Commun* 14, 4883 (2023). https://doi.org/10.1038/s41467-023-40457-w

2. Robinson JL, Kocabaş P, Wang H, et al. An atlas of human metabolism. *Sci Signal* 13, eaaz1482 (2020).

3. Orth JD, Thiele I, Palsson BØ. What is flux balance analysis? *Nat Biotechnol* 28, 245–248 (2010).

## Links

* **METAFlux2 Repository**: https://github.com/Zaoqu-Liu/METAFlux2
* **Original METAFlux**: https://github.com/KChen-lab/METAFlux
* **Nature Communications Paper**: https://www.nature.com/articles/s41467-023-40457-w
* **Human-GEM Model**: https://github.com/SysBioChalmers/Human-GEM

---

<p align="center">
  <strong>METAFlux2: Enabling efficient large-scale metabolic flux analysis</strong>
  <br>
  Based on the original METAFlux by Huang et al. (2023)
</p>
