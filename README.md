---
title: "RTM: Regression to the Mean R Package"
output: github_document
---

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/umair-statistics/RTM)](https://github.com/umair-statistics/RTM/releases)
[![GitHub all releases](https://img.shields.io/github/downloads/umair-statistics/RTM/total)](https://github.com/umair-statistics/RTM/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

### Version 0.2.0

### Muhammad Umair and Manzoor Khan

**RTM**: A R package for computing Regression to the Mean (RTM) effect and test the treatment effect under bivariate truncated distributions

Please report any **bugs** or **suggestions** at:
<https://github.com/umair-statistics/RTM/issues>.

### Installation

You can install the most recent development version of **RTM** using the [devtools](https://github.com/r-lib/devtools) function `install_github()`.

However, you need to make sure you're set up to develop packages. This is platform specific:

* On Windows, download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
* On the Mac, make sure you have [Xcode](https://developer.apple.com/xcode/) installed.
* On Linux, make sure you have the R-dev packages installed.

You can check everything is installed correctly with the `has_devel()` function from the **devtools** package. Type the following at 
the **R** prompt:


```r
install.packages("devtools", dependencies = TRUE)    
devtools::has_devel()
```

If everything is installed correctly, the function will print some output and then return **TRUE**.

To install the **RTM** package, type the following at the **R** prompt:


```r
devtools::install_github('umair-statistics/RTM')
```
    
It is possible to install **RTM** with [GIT](https://git-scm.com/) and the **R CMD build** assuming you have GIT installed and the appropriate tools to build **R** from source.

```bash
git clone https://github.com/umair-statistics/RTM.git
R CMD build RTM
R CMD INSTALL RTM_*.tar.gz
```
## Contributing

Contributions are welcome!

- Fork the repo

- Create a feature branch

- Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
