
RTM
========

### Version 0.1.0


<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.


### Muhammad Umair, Manzoor khan

**RTM**: Estimate regression to the mean effect and test the treatment means under different bivariate truncated distributions 

Please report any **bugs** or **suggestions** at:
<https://github.com/umair1nonly/RTM/issues>.

### Installation

The stable version of the package is available for download from github at <https://github.com/umair1nonly/RTM>

You may install the most recent development version of **RTM** using the [devtools](https://github.com/r-lib/devtools) function `install_github()`.

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


