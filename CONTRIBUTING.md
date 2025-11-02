# Contributing to METAFlux2

Thank you for your interest in contributing to METAFlux2! This document provides guidelines for contributing to the project.

## 🤝 Code of Conduct

By participating in this project, you agree to maintain a respectful and collaborative environment for all contributors.

## 🎯 How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When you create a bug report, include:

1. **Clear title and description**
2. **Minimal reproducible example**
3. **Expected vs actual behavior**
4. **System information** (R version, OS, package versions)
5. **Error messages** (complete stack trace)

Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md).

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When suggesting an enhancement:

1. **Use a clear and descriptive title**
2. **Provide a detailed description** of the suggested enhancement
3. **Explain why this enhancement would be useful**
4. **Provide examples** of how it would work

Use the [feature request template](.github/ISSUE_TEMPLATE/feature_request.md).

### Pull Requests

#### Before You Start

1. **Check existing PRs** to avoid duplicate work
2. **Open an issue** to discuss major changes before starting
3. **Fork the repository** and create your branch from `main`
4. **Follow coding standards** (see below)

#### Development Process

1. **Fork and clone** the repository
   ```bash
   git clone https://github.com/your-username/METAFlux2.git
   cd METAFlux2
   ```

2. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bug-fix
   ```

3. **Set up development environment**
   ```r
   library(devtools)
   install_dev_deps()
   ```

4. **Make your changes**
   - Write clear, commented code
   - Follow existing code style
   - Add tests for new functionality
   - Update documentation as needed

5. **Test your changes**
   ```r
   # Run all tests
   devtools::test()
   
   # Check package
   devtools::check()
   
   # Validate numerical accuracy
   source("tests/validation_test.R")
   run_all_tests()
   ```

6. **Commit your changes**
   ```bash
   git add .
   git commit -m "feat: add new feature" # or "fix: resolve bug"
   ```

7. **Push and create PR**
   ```bash
   git push origin feature/your-feature-name
   ```

## 📝 Coding Standards

### R Code Style

Follow the [tidyverse style guide](https://style.tidyverse.org/). Key points:

1. **Use explicit package namespacing**
   ```r
   # Good
   Matrix::Diagonal(n, 1)
   osqp::osqp(P, q, A, l, u)
   
   # Avoid
   Diagonal(n, 1)  # Without package prefix
   ```

2. **Use meaningful variable names**
   ```r
   # Good
   metabolic_scores <- calculate_reaction_score(data)
   
   # Avoid
   ms <- calculate_reaction_score(data)
   ```

3. **Comment complex logic**
   ```r
   # Calculate weighted average using cell type proportions
   weighted_flux <- flux %*% diag(proportions)
   ```

4. **Use consistent indentation** (2 spaces)

5. **Keep lines under 80 characters** when possible

### Performance Guidelines

1. **Vectorize operations** when possible
   ```r
   # Good
   result <- colSums(matrix / vector)
   
   # Avoid
   result <- apply(matrix, 2, function(x) sum(x / vector))
   ```

2. **Pre-allocate memory** for large objects
   ```r
   # Good
   result <- matrix(NA, nrow = n, ncol = m)
   for (i in 1:n) result[i, ] <- compute(i)
   
   # Avoid
   result <- NULL
   for (i in 1:n) result <- rbind(result, compute(i))
   ```

3. **Use sparse matrices** for large, sparse data
   ```r
   # Good
   sparse_mat <- Matrix::Matrix(0, nrow = n, ncol = m, sparse = TRUE)
   
   # Avoid
   dense_mat <- matrix(0, nrow = n, ncol = m)
   ```

### Documentation Standards

1. **Use roxygen2** for function documentation
   ```r
   #' Calculate metabolic reaction activity scores
   #'
   #' @param data Gene expression matrix (genes x samples)
   #' @return Matrix of metabolic reaction activity scores
   #' @export
   #' @examples
   #' data("bulk_test_example")
   #' scores <- calculate_reaction_score(bulk_test_example)
   calculate_reaction_score <- function(data) {
     # ...
   }
   ```

2. **Update NEWS.md** for user-facing changes

3. **Update README.md** if adding major features

## ✅ Testing Requirements

### For Bug Fixes

- Add a test that fails before the fix
- Ensure the test passes after the fix
- Verify existing tests still pass

### For New Features

- Add comprehensive tests covering:
  - Normal use cases
  - Edge cases
  - Error conditions
- Achieve >80% code coverage for new code
- Validate numerical accuracy against expected results

### For Performance Improvements

- Benchmark before and after
- Document performance gains
- Ensure numerical accuracy is preserved
  ```r
  # Validate results are identical
  all.equal(original_result, optimized_result, tolerance = 1e-10)
  ```

## 📊 Performance Validation

For performance-related contributions:

1. **Run benchmarks**
   ```r
   library(microbenchmark)
   
   result <- microbenchmark(
     original = original_function(...),
     optimized = optimized_function(...),
     times = 10
   )
   print(result)
   ```

2. **Profile code**
   ```r
   library(profvis)
   profvis({
     # Your code here
   })
   ```

3. **Document improvements** in PR description

## 🔬 Numerical Accuracy

**Critical**: All optimizations must preserve numerical accuracy.

```r
# Always validate
test_that("Optimized version matches original", {
  result_original <- original_function(test_data)
  result_optimized <- optimized_function(test_data)
  
  expect_equal(
    result_original, 
    result_optimized, 
    tolerance = 1e-10,
    info = "Results must be numerically identical"
  )
})
```

## 📦 Package Dependencies

### Adding New Dependencies

1. **Minimize new dependencies** - only add if truly necessary
2. **Use established packages** - prefer CRAN packages
3. **Document why needed** in PR description
4. **Update DESCRIPTION** file:
   ```
   Imports:
       newpackage (>= 1.0.0)
   ```

### Acceptable Reasons for New Dependencies

- ✅ Significant performance improvement
- ✅ Essential functionality not available in base R or existing deps
- ✅ Industry standard package widely used
- ❌ Convenience functions that can be easily implemented
- ❌ Packages with many transitive dependencies

## 🎓 Learning Resources

### Understanding the Codebase

1. Read the [original METAFlux paper](https://www.nature.com/articles/s41467-023-40457-w)
2. Review the [optimization guide](docs/OPTIMIZATION_GUIDE.md)
3. Study the [flux balance analysis fundamentals](https://www.nature.com/articles/nbt.1614)

### R Package Development

- [R Packages book](https://r-pkgs.org/) by Hadley Wickham
- [Tidyverse style guide](https://style.tidyverse.org/)
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)

## 💬 Communication

- **GitHub Issues**: For bug reports and feature requests
- **Pull Requests**: For code contributions
- **Discussions**: For questions and general discussion (if enabled)

## 🙏 Recognition

Contributors will be acknowledged in:
- README.md (for significant contributions)
- DESCRIPTION file (as contributors)
- Release notes

## 📜 License

By contributing to METAFlux2, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to METAFlux2! Your efforts help make metabolic flux analysis faster and more accessible to the research community. 🚀
