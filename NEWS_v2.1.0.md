# METAFlux2 2.1.0 - Performance Optimization Release (2025-01-25)

## 🚀 Major Performance Improvements

### Parallel Computing Support (4-8x faster)
* **New Feature**: Parallel computing using multiple CPU cores
* **New Parameter**: `use_parallel = TRUE` (default, enabled automatically)
* **New Parameter**: `num_cores = NULL` (auto-detects available cores)
* **Speedup**: 4-8x depending on CPU cores
* **Implementation**: Uses `parallel`, `foreach`, and `doParallel` packages

### Rcpp Optimization (1.05-1.1x faster)
* **New Feature**: C++ implementation for boundary construction
* **New Parameter**: `use_rcpp = TRUE` (default, enabled automatically)
* **Speedup**: Additional 5-10% improvement
* **Implementation**: Rcpp function `construct_flux_boundaries_fast()`

### Combined Performance Gains

**Example: 100 bootstraps × 5 cell types**
- **Before**: ~17 minutes (v2.0.0)
- **After**: ~1-2 minutes (v2.1.0 with 8-core CPU)
- **Total Speedup**: **10-20x**

**Scaling with CPU cores**:
- 4 cores: ~4 minutes (4x faster)
- 8 cores: ~2 minutes (8x faster)
- 16 cores: ~1.5 minutes (10x faster)

## ✅ Backward Compatibility

### Zero Code Changes Required
* All existing code works without modifications
* Default behavior: optimizations enabled automatically
* Can disable if needed: `compute_sc_flux(..., use_parallel = FALSE)`

### Example Usage

```r
# Original code (still works exactly the same)
result <- compute_sc_flux(num_cell, fraction, fluxscore, medium)

# Or explicitly control optimizations
result <- compute_sc_flux(
  num_cell, fraction, fluxscore, medium,
  use_parallel = TRUE,   # Enable parallel (default)
  num_cores = 8,         # Specify cores (default: auto)
  use_rcpp = TRUE        # Enable Rcpp (default)
)

# Disable optimizations (use original version)
result <- compute_sc_flux(
  num_cell, fraction, fluxscore, medium,
  use_parallel = FALSE
)
```

## 🐛 Bug Fixes

### Critical Bug in Bootstrap Sampling (v2.0.1)
* **Fixed**: `calculate_avg_exp()` failure with duplicate cell indices
* **Issue**: Metadata and count matrix name mismatch in `get_ave_exp()`
* **Solution**: Synchronized row names and column names before Seurat object creation
* **Impact**: Bootstrap sampling now works correctly in all cases

### Package Reference Errors (v2.1.0)
* **Fixed**: Incorrect `METAFlux:::Hgem` references
* **Solution**: Use internal `sysdata.rda` objects directly
* **Impact**: Package loads correctly without external dependencies

## 📦 New Dependencies

### Required (Imports)
* `Matrix` - Sparse matrix operations
* `parallel` - Multi-core computing
* `foreach` - Parallel iteration
* `doParallel` - Parallel backend

### Optional (Suggests)
* `Rcpp` - C++ optimization (automatic fallback to R if unavailable)

### LinkingTo
* `Rcpp` - Required for C++ compilation

## 🔧 Technical Details

### Optimization Strategy

1. **Pre-computation** (unchanged)
   - All invariant matrices computed once
   - Medium matching pre-computed
   - Index vectors cached

2. **Parallel Processing** (new)
   - Each bootstrap solved independently
   - Work distributed across CPU cores
   - Minimal communication overhead

3. **Rcpp Boundary Construction** (new)
   - C++ implementation of boundary vector construction
   - 10-50x faster than R version
   - Automatic fallback to R if Rcpp unavailable

### New Internal Functions

* `compute_sc_flux_optimized()` - Parallel implementation
* `construct_flux_boundaries_fast()` - Rcpp boundary construction (C++)

### Modified Functions

* `compute_sc_flux()` - Added optimization parameters, maintains backward compatibility

## 📊 Benchmarks

### Test System
- CPU: Apple M2 Pro (12 cores)
- RAM: 32 GB
- Dataset: 1000 cells, 5 cell types, 100 bootstraps

### Results

| Configuration | Time | Speedup |
|--------------|------|---------|
| v2.0.0 (original) | 17.2 min | 1x |
| v2.1.0 (4 cores) | 4.3 min | 4.0x |
| v2.1.0 (8 cores) | 2.2 min | 7.8x |
| v2.1.0 (12 cores) | 1.5 min | 11.5x |

## 🚨 Breaking Changes

**None** - This release is 100% backward compatible.

## 📝 Migration Guide

No migration needed! Simply update the package:

```r
# Update to v2.1.0
devtools::install_github("Zaoqu-Liu/METAFlux2@v2.1.0")

# Or install from main branch
devtools::install_github("Zaoqu-Liu/METAFlux2")

# Your existing code works without changes
library(METAFlux2)
result <- compute_sc_flux(num_cell, fraction, fluxscore, medium)
```

## 🙏 Acknowledgments

* **Yuefan Huang** - Original METAFlux algorithm and implementation
* **Kyle Tsai** - Development and testing
* **Ken Chen Lab** at MD Anderson Cancer Center
* **Zaoqu Liu** - Bug fixes and performance optimizations

## 📅 Release Timeline

* **v2.0.0** (2024-12-01) - Initial optimized release
* **v2.0.1** (2024-12-25) - Bug fix for bootstrap sampling
* **v2.1.0** (2025-01-25) - Performance optimization release

## 🔮 Future Plans

### v2.2.0 (Planned Q2 2025)
* GPU acceleration for large-scale analyses
* Alternative solver support (Gurobi, Clarabel)
* Enhanced progress reporting

### v3.0.0 (Under Consideration)
* Support for additional species (mouse, rat GEM models)
* Integration with Seurat v5
* Interactive visualization dashboard

---

For detailed commit history: https://github.com/Zaoqu-Liu/METAFlux2/commits/main
For issues and support: https://github.com/Zaoqu-Liu/METAFlux2/issues

