
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tbSTATIS

![Package Version](https://img.shields.io/badge/version-0.1.0-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
<!-- badges: start --> <!-- badges: end -->
For documentation and a tutorial, checkout the [website!](https://samalatesta.github.io/tbSTATIS/index.html)

The `tbSTATIS` package implements the tuberculosis SeveriTy Assessment
Tool for Informed Stratification (TB-STATIS) to classify TB disease
severity for individuals using data collected at time of diagnosis. We
include functions to estimate the model, quantify model uncertainty, and
visualize results. We provide a detailed vignette that walks through an
example analysis using `tbSTATIS`.

For a more formal discussion of the theory behind and usage of this
method, see the following paper:

Malatesta S, Jacobson KR, Horsburgh CR, Farhat MR, Gile KJ, Kolaczyk ED,
White LF. An Integrated Data-Driven Model for Clinical Phenotyping of Tuberculosis Disease Severity.

Code to replicate all simulation results from this paper is available
[here](https://github.com/samalatesta/tbSTATISpaper).

## Installation

You can install the development version of `tbSTATIS` from
[GitHub](https://github.com/). We recommend setting the option
`build_vignettes=T` when installing so the package vignette can be
accessed in your local R environment (note: this may take a while). The
`devtools` package must be installed prior to installing `tbSTATIS`.

``` r
# install.packages("devtools")
devtools::install_github("samalatesta/tbSTATIS", build.vignette = T)
```

## Usage

To use `tbSTATIS` in your R scripts or projects, load the package using:

``` r
library(tbSTATIS)
```

For detailed information on how to use each function, please refer to
the package documentation and vignette. The vignette can be viewed
locally after package installation or the knitted html is also included
in the `\vignettes` directory.

``` r
vignette(package="tbSTATIS")
```

## Documentation

Comprehensive documentation for tbSTATIS functions is available within
R. You can access documentation using the ? operator followed by the
function name. For example:

``` r
?plot_states
```

## License

tbSTATIS is distributed under the MIT License.

## Contact

For questions, please contact:

Samantha Malatesta (<samalate@bu.edu>)
