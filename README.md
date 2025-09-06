# HPSX - High Performance Statistics

[![Build Status](https://github.com/pooriayousefi/hpsx/actions/workflows/ci.yml/badge.svg)](https://github.com/pooriayousefi/hpsx/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B20)
[![CMake](https://img.shields.io/badge/CMake-3.20%2B-green.svg)](https://cmake.org/)

A modern, high-performance C++20 header-only library for comprehensive descriptive statistical calculations. HPSX provides optimized implementations of statistical functions with support for parallel execution and advanced mathematical computations.

## 🚀 Features

- **Header-Only**: Single header include for easy integration
- **C++20 Standard**: Modern C++ with concepts, ranges, and coroutines
- **High Performance**: Optimized algorithms with SIMD-friendly implementations
- **Parallel Execution**: Support for STL execution policies
- **Comprehensive Statistics**: Wide range of descriptive statistical functions
- **Type Safety**: Concept-constrained templates for compile-time validation
- **Cross-Platform**: Works on Linux (g++), macOS (clang++), and Windows (MSVC)

## 🎯 Quick Start

```cpp
#include "hpsx.h"
#include <vector>
#include <iostream>

int main() {
    std::vector<double> data{1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.1};
    
    // Basic statistics
    auto mean_val = hp::mean(data.begin(), data.end());
    auto median_val = hp::median(data.begin(), data.end());
    auto std_dev = hp::standard_deviation(data.begin(), data.end());
    
    std::cout << "Mean: " << mean_val << std::endl;
    std::cout << "Median: " << median_val << std::endl;
    std::cout << "Standard Deviation: " << std_dev << std::endl;
    
    // Advanced statistics
    auto skewness = hp::skewness(data.begin(), data.end());
    auto kurtosis = hp::kurtosis(data.begin(), data.end());
    
    std::cout << "Skewness: " << skewness << std::endl;
    std::cout << "Kurtosis: " << kurtosis << std::endl;
    
    return 0;
}
```

## 📚 API Reference

### Central Tendency

```cpp
// Arithmetic mean
template<typename Iterator>
auto mean(Iterator first, Iterator last);

// Geometric mean
template<typename Iterator>
auto geometric_mean(Iterator first, Iterator last);

// Harmonic mean
template<typename Iterator>
auto harmonic_mean(Iterator first, Iterator last);

// Median
template<typename Iterator>
auto median(Iterator first, Iterator last);

// Mode
template<typename Iterator>
auto mode(Iterator first, Iterator last);
```

### Measures of Dispersion

```cpp
// Variance
template<typename Iterator>
auto variance(Iterator first, Iterator last);

// Standard deviation
template<typename Iterator>
auto standard_deviation(Iterator first, Iterator last);

// Range
template<typename Iterator>
auto range(Iterator first, Iterator last);

// Interquartile range
template<typename Iterator>
auto interquartile_range(Iterator first, Iterator last);
```

### Distribution Shape

```cpp
// Skewness (asymmetry measure)
template<typename Iterator>
auto skewness(Iterator first, Iterator last);

// Kurtosis (tail heaviness measure)
template<typename Iterator>
auto kurtosis(Iterator first, Iterator last);

// Coefficient of variation
template<typename Iterator>
auto coefficient_of_variation(Iterator first, Iterator last);
```

### Quantiles and Percentiles

```cpp
// Quantile calculation
template<typename Iterator>
auto quantile(Iterator first, Iterator last, double q);

// Percentile calculation
template<typename Iterator>
auto percentile(Iterator first, Iterator last, double p);

// Five-number summary
template<typename Iterator>
auto five_number_summary(Iterator first, Iterator last);
```

## 🔬 Advanced Features

### Parallel Execution

```cpp
// Parallel mean calculation
auto parallel_mean = hp::mean(std::execution::par, 
                             data.begin(), data.end());

// Parallel standard deviation
auto parallel_std = hp::standard_deviation(std::execution::par_unseq,
                                          data.begin(), data.end());
```

### Statistical Tests

```cpp
// Normality tests
bool is_normal = hp::shapiro_wilk_test(data.begin(), data.end());

// Outlier detection
auto outliers = hp::detect_outliers(data.begin(), data.end());

// Correlation analysis
auto correlation = hp::pearson_correlation(x.begin(), x.end(),
                                          y.begin(), y.end());
```

## 🏗️ Building from Source

```bash
# Clone repository
git clone https://github.com/pooriayousefi/hpsx.git
cd hpsx

# Build with CMake
cmake --preset=default
cmake --build build/default

# Run examples
./build/default/hpsx_example
```

## 📊 Use Cases

- **Data Science**: Exploratory data analysis and preprocessing
- **Research**: Statistical analysis in scientific studies
- **Quality Control**: Process monitoring and control charts
- **Finance**: Risk analysis and portfolio optimization
- **Machine Learning**: Feature engineering and model validation

## 🔧 Requirements

- C++20 compatible compiler (GCC 11+, Clang 13+, MSVC 2022+)
- CMake 3.20 or later
- Standard library with concepts and execution policy support

## 📈 Performance

- **Optimized Algorithms**: Cache-friendly implementations
- **SIMD Support**: Vectorized operations where applicable
- **Parallel Processing**: Multi-core utilization
- **Memory Efficient**: Minimal allocations and cache-optimal access patterns

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Author**: [Pooria Yousefi](https://github.com/pooriayousefi)  
**Repository**: [https://github.com/pooriayousefi/hpsx](https://github.com/pooriayousefi/hpsx)
