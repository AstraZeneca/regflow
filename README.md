![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red) 

<img
    src="https://github.com/osysoev/regflux/blob/main/R/RFAname.png"
    width="400">

## regflow

An R package for carrying out Regulatory Flow Analysis.


### Installation

To install **`regflow`** from this GitHub repository, 
use the function `install_github` in the 
[devtools](https://cran.r-project.org/package=devtools) package. 

```R
library( devtools )

install_github( "AstraZeneca/regflow" )   
```


### Help

Stardard help is available for the package and exported functions:
```
library( regflow )
? `regflow-package`
```


### Examples

A simple example is available in the directory `demo`: 
```
demo( package = "regflow" )
```

## Supported platforms

The package was tested on Windows and Unix platforms.

## Dependencies

The following packages need to be installed before the package can be used: biomaRt, ggplot2, ggrepel, igraph, Matrix, moments, OmnipathR, RSpectra