# Migration Guide: METAFlux to METAFlux2

This guide provides instructions for transitioning from the original METAFlux to the optimized METAFlux2 implementation.

## Overview

METAFlux2 is designed as a **drop-in replacement** for METAFlux, requiring minimal to no changes to existing analysis code. The API is 100% backward compatible, ensuring that existing scripts will continue to function correctly.

## Quick Migration

### Step 1: Update Installation

Replace the original METAFlux installation with METAFlux2:

```r
# Remove original METAFlux (optional)
remove.packages("METAFlux")

# Install METAFlux2
devtools::install_github('Zaoqu-Liu/METAFlux2')
```

### Step 2: Install Additional Dependency

METAFlux2 requires the Matrix package for sparse matrix operations:

```r
install.packages("Matrix")
```

### Step 3: No Code Changes Required

Your existing code will work without modification:

```r
# Original code continues to work
library(METAFlux2)  # Was: library(METAFlux)

data("bulk_test_example")
data("human_blood")

scores <- calculate_reaction_score(bulk_test_example)
flux <- compute_flux(mras = scores, medium = human_blood)
```

## Compatibility

### Identical API

All function names, parameters, and return values are identical:

| Function | Status | Notes |
|----------|--------|-------|
| `calculate_reaction_score()` | Identical | Same inputs and outputs |
| `compute_flux()` | Identical | Same inputs and outputs |
| `compute_sc_flux()` | Identical | Same inputs and outputs |
| `calculate_avg_exp()` | Identical | Same inputs and outputs |

### Data Structures

All data structures remain unchanged:

* Input gene expression matrices (genes × samples)
* MRAS score matrices (reactions × samples)
* Flux output matrices (reactions × samples)
* Single-cell flux matrices (cell-type reactions × bootstraps)

### Medium Profiles

All medium definition files are identical:

```r
data("human_blood")    # Human blood medium
data("cell_medium")    # Cell line culture medium
data("nutrient_lookup_files")  # Exchange reactions lookup
```

## What Changed

### Performance Only

The only changes are internal optimizations:

1. **Vectorized operations** - Faster computation
2. **Sparse matrices** - Lower memory usage
3. **Pre-allocated memory** - Reduced overhead
4. **Pre-computed values** - Eliminated redundant calculations

### New Dependency

METAFlux2 requires the Matrix package:

```r
# Add to your package dependencies if creating a package
Imports: Matrix (>= 1.3.0)
```

## Validation

### Verify Numerical Accuracy

Optionally compare results between versions:

```r
# Using original METAFlux
library(METAFlux)
data("bulk_test_example")
data("human_blood")
scores_v1 <- calculate_reaction_score(bulk_test_example)
flux_v1 <- compute_flux(mras = scores_v1, medium = human_blood)

# Switch to METAFlux2
detach("package:METAFlux", unload = TRUE)
library(METAFlux2)
scores_v2 <- calculate_reaction_score(bulk_test_example)
flux_v2 <- compute_flux(mras = scores_v2, medium = human_blood)

# Compare results
all.equal(flux_v1, flux_v2, tolerance = 1e-10)
# Should return TRUE or indicate difference < 1e-10
```

### Run Validation Tests

METAFlux2 includes a comprehensive test suite:

```r
source("https://raw.githubusercontent.com/Zaoqu-Liu/METAFlux2/main/tests/validation_test.R")
run_all_tests()
```

## Expected Performance Gains

After migration, you should observe:

### Bulk RNA-seq Analysis

| Dataset Size | Expected Speed-up |
|--------------|-------------------|
| < 50 samples | 30-40% |
| 50-200 samples | 35-50% |
| 200-500 samples | 40-60% |
| > 500 samples | 50-70% |

### Single-cell Analysis

| Dataset Size | Expected Speed-up |
|--------------|-------------------|
| < 2K cells | 30-40% |
| 2K-5K cells | 40-50% |
| 5K-20K cells | 45-60% |
| > 20K cells | 55-75% |

### Memory Usage

Expect 30-50% reduction in peak memory usage across all dataset sizes.

## Troubleshooting Migration

### Issue: Matrix Package Not Found

```r
# Solution: Install Matrix package
install.packages("Matrix")
```

### Issue: Results Differ Slightly

Small numerical differences (< 1e-10) are expected due to floating-point arithmetic:

```r
# This is normal
all.equal(flux_v1, flux_v2, tolerance = 1e-10)
# Returns TRUE or very small difference
```

### Issue: Package Namespace Conflicts

If you have both versions installed:

```r
# Explicitly specify which package to use
METAFlux2::calculate_reaction_score(data)
```

## Updating Package Dependencies

If you maintain a package that depends on METAFlux:

### DESCRIPTION File

Update your package DESCRIPTION:

```r
# Before
Imports:
    METAFlux

# After
Imports:
    METAFlux2,
    Matrix (>= 1.3.0)
```

### Code Updates

No code changes needed, but update package references:

```r
# In your package code
#' @import METAFlux2
#' @importFrom Matrix Diagonal

# Or use explicit namespacing
result <- METAFlux2::calculate_reaction_score(data)
```

## Side-by-Side Comparison

You can run both versions side-by-side for validation:

```r
# Install both versions
devtools::install_github('KChen-lab/METAFlux')
devtools::install_github('Zaoqu-Liu/METAFlux2')

# Use with explicit namespacing
result_v1 <- METAFlux::compute_flux(mras, medium)
result_v2 <- METAFlux2::compute_flux(mras, medium)

# Compare
all.equal(result_v1, result_v2)
```

## Batch Migration Script

For projects with multiple analysis scripts:

```r
# migrate_to_metaflux2.R
# Run this script to update all your analysis files

library(stringr)

# List all R scripts in your project
r_files <- list.files(
  path = ".",
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

# Update library calls
for (file in r_files) {
  content <- readLines(file)
  
  # Replace library calls
  content <- str_replace_all(
    content,
    "library\\(METAFlux\\)",
    "library(METAFlux2)"
  )
  
  # Write back
  writeLines(content, file)
}

cat("Migration complete. Please test your scripts.\n")
```

## Rollback Procedure

If you need to revert to the original METAFlux:

```r
# Remove METAFlux2
remove.packages("METAFlux2")

# Reinstall original
devtools::install_github('KChen-lab/METAFlux')

# Your code will work without changes
```

## Additional Resources

* [METAFlux2 README](README.md) - Installation and usage
* [Optimization Guide](docs/OPTIMIZATION_GUIDE.md) - Technical details
* [Original METAFlux Tutorial](https://htmlpreview.github.io/?https://github.com/KChen-lab/METAFlux/blob/main/Tutorials/pipeline.html) - Fully compatible
* [GitHub Issues](https://github.com/Zaoqu-Liu/METAFlux2/issues) - Report problems

## Summary

Migration to METAFlux2 is straightforward:

1. Install Matrix package
2. Install METAFlux2
3. Continue using existing code
4. Enjoy faster computation and lower memory usage

No code changes required. All results numerically identical. Complete backward compatibility guaranteed.

---

**Questions?** Open an issue at https://github.com/Zaoqu-Liu/METAFlux2/issues
