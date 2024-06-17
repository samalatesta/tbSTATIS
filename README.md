
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TBebm

![Package Version](https://img.shields.io/badge/version-0.1.0-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
<!-- badges: start --> <!-- badges: end -->

The TBebm package implements the event-based model to classify disease
severity for TB patients using data collected at time of diagnosis. We
include functions to estimate the model, quantify model uncertainty, and
visualize results. We provide a detailed vignette that walks through an
example analysis using `TBebm`.

## Installation

You can install the development version of `TBebm` from
[GitHub](https://github.com/). We recommend setting the option
`build_vignettes=T` when installing so the package vignette can be
accessed in your local R environment. The `devtools` package must be
installed prior to installing `TBebm`.

``` r
# install.packages("devtools")
devtools::install_github("samalatesta/TBebm", build.vignette = T)
```

## Usage

To use `TBebm` in your R scripts or projects, load the package using:

``` r
library(TBebm)
```

For detailed information on how to use each function, please refer to
the package documentation and vignette. The vignette can be viewed
locally after package installation or the knitted html is also included
in the `\vignettes` directory.

``` r
vignette(package="TBebm")
```

## Documentation

Comprehensive documentation for TBebm functions is available within R.
You can access documentation using the ? operator followed by the
function name. For example:

``` r
?plot_events
```

## License

TBebm is distributed under the MIT License.

## Contact

For questions,upport, or contributions, please contact:

Samantha Malatesta (<samalate@bu.edu>)
