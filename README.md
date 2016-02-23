# Full Process Kinetics PCR Analysis
This branch contains the FPK-PCR algorithm (V06-11) originally published in ..... (link), its updates and expansions. We refer the user to the original publication for the interpretation of results and background on the procedure.  

## FPK-PCR Algortihtm   
Currently two versions of the algorithm are hosted here: the original V06-11 and V13-13 (a drastic update with model changes and expanded functionality). The most important changes since V06-11 are:
*  Addition of the baseline optimisation system, computationally intensive, but optimizes the baseline intercept so that the expected value of the residuals is as close as possible to zero (irrespective of their size). This prevents biased effiency estimates because of baseline over or under subtraction.
*  The FPK PCR model for the first phase of efficiency decline now has less degrees of freedom, this reduces the variability in the estimates.
* The Kalman filter is implemented to try to scavange as much information as possible from the groundphase.
* Kalman covariance is estimated via bootstrap 
* Resec implementation to weed out non-amplification fits. Surprisingly, the model will succeed in fitting itselve to random data (no amplification, only noise) surprisingly often. When analysing datasets that have a lot of negatives this can highly annoting, hence this threshold.
*  The start of the model fitting no longer lies at cycle(5%)+1 but now lies at cycle(7.5%), this has further reduced variability in effiency estimates
*  bootstrap functionality has been added. This is implemented via a second function (...) which is called within the FPK procedure and carries out the actual resampling. This makes for a compuationally intensive analysis, but gives you confidence intervals for the parameter estimates. 
*  debugg option added 

Below you will find a short description of the usage, arguments, and output of each version. The input requirements in terms of real time PCR data are identical, so both can be used on the same data format.

#### Version FPK-V06-11

Usage: 
```
analyse(fluo, plots = FALSE, param = FALSE, silent = TRUE, vecto = FALSE)
```

Arguments:
- `fluo` is a numerical vector containing the consecutive fluorescence measurments (one for each cycle).
- `plots`is a logical, if `TRUE` a number of informative graphical representations of the data are returned in addition to the other output. The plots are:
  * data and fitted model (actual fluorescence readings + values reconstructed by the model in green)
  * Efficiency plots: (**A**) double logarithmic transformed efficiency vs. fluorescence + bilinear model (vertical dotted lines represent the portion of the data between which the first phase was fitted). (**B**) efficiency vs. cycle + FPK model
- `param` is a Logical, if `TRUE` the function returns the individual model parameters rather than the standard output
- `silent` is a Logical, if `TRUE` no messages are returned 
- `vecto`is a Logical, if `TRUE` the data are returned as a vector rather than a list. 

Output: 
Standard output returns three values:
* `Cq`: quantification cycle, calculated at the first positive second derivative maximum of a five parameter logistic model.
* `Emax`: initial (maximal) efficiency of the reaction
* `a*i0`: initial reaction fluorescence (number of initial target copies times the fluorescence of a single amplicon), expressed in Fluorescence Units (FU)

parameter output (`param = TRUE`) returns three sets of parameters:

1 `Baseline`: parameters of the baseline (linear model)
  * `intercept` = y-intercept of the linear model
  * `slope` = slope of the linear model

2 `5PLM`: parameters of the five paramter logistic model (with sloped baseline):
 * `y0` = a numeric parameter representing the intercept of the baseline
 * `s` = a numeric parameter representing the slope of the baseline
 * `fmax` = a numeric parameter representing the horizontal asymptote on the right side (plateau)
 * `xmid` = a numeric parameter representing the input value at the center of the curve
 * `b` = a numeric scale parameter on the input axis, the 'growth rate'
 * `g` = the reciprocal of the shape parameter: affects near which asymptote maximum 'growth' occurs
```
5PLM model
Fn = fmax + s*x + (y0 - fmax) / (1 + (2^(1/g)-1) * exp(b*(x-xmid)))^g
```

3 `Bilinear`: parameters of the bilinear model:
 * `a1a` = "slope" describes together with a1b the curve of the first phase of efficiency decline
 * `a1b` = "slope" describes together with a1a the curve of the first phase of efficiency decline
 * `a2` = "slope" represents the steepness of the second phase of efficiency decline
 * `ipt` = corresponds to the horizontal position of the phase change
 * `eta` = affects the abruptness of transition between the two phases
 * `chi` = affects the y-intercept of the curve
```
Bilinear model:
Ln²(En) = chi + eta *log(exp(( a1a * (input-ipt)^2+ a1b *(input-ipt))/eta) + exp( a2 *(n-ipt)/eta))
```  
  
#### Version FPK-V13-13

Usage:
```
analyse(fluo, plots = FALSE, output = "estimates", base.line = "optimized", bootstrap = FALSE,  silent = FALSE, kalman = TRUE, debugg = FALSE, n = 1000, p = 0.05, resec.tolerance = 0.125) 
```

Arguments:
- `fluo` is a numerical vector containing the consecutive fluorescence measurments (one for each cycle).
- `plots`is a logical, if `TRUE` a number of informative graphical representations of the data are returned in addition to the other output. The plots are:
  * data and fitted model (actual fluorescence readings + values reconstructed by the model in green)
  * Efficiency plots: (**A**) double logarithmic transformed efficiency vs. fluorescence + bilinear model (vertical dotted lines represent the portion of the data between which the first phase was fitted). (**B**) efficiency vs. cycle + FPK model
- `base.line` a character vector of length 1 whose value is either "optimized", "slanted", or "flat" depending on the baseline fitting approach to be used.
- `output`  a character vector of length 1 whose value is either "estimates" (for E, ai0 & Cq), "parameters" (for model parameters as a list), or "all" (both as a concatenated vector)
- `silent` is a Logical, if `TRUE` no messages are returned 
- `kalman` is a logical, if `TRUE` the kalman filter is used to include groundphase data into the analysis
- `debugg` is a logical, if `TRUE`  intermediate values, warning messages, and diagnosic plots are returned
- `resec.tolerance` is a numerical vector of length 1. It determines a cut-off for the 5PLM residual error per cycle. This  serves as a tool in low-level presence analysis (reactions with high values are not likely to contain bona fide amplification, `NA` output is returned for reactions that exceed the threshold)
- `bootstrap` is a logical, if `TRUE` a bootstrap analysis will be performed and BCa confidence intervals will be returned for the parameter estimates. **Warning:** this is rather computationally intensive, you may want to go get a coffee or something. Better first run a few reactions at a time to see how well systems handles the punishment. Also, you may want to switch off the baseline optimization for the bootstrap (just use `base.line = "slanted"`) since its usage increases the computational load rather exponentially. 
- `n` is a numerical of length 1, the number of bootstraps that will be performed.
- `p` is a numerical of length 1, the p-value of the confidence bounds calculated.

Output
`ouptut = "estimates"` returns three values:
* `Cq`: quantification cycle, calculated at the first positive second derivative maximum of a five parameter logistic model.
* `Emax`: initial (maximal) efficiency of the reaction
* `a*i0`: initial reaction fluorescence (number of initial target copies times the fluorescence of a single amplicon), expressed in Fluorescence Units (FU)

`ouptut = "parameters"` returns four sets of parameters:

1 `Baseline`: parameters of the baseline (linear model)
  * `intercept` = y-intercept of the linear model
  * `slope` = slope of the linear model

2 `5PLM`: parameters of the five paramter logistic model (with sloped baseline):
 * `y0` = a numeric parameter representing the intercept of the baseline
 * `s` = a numeric parameter representing the slope of the baseline
 * `fmax` = a numeric parameter representing the horizontal asymptote on the right side (plateau)
 * `xmid` = a numeric parameter representing the input value at the center of the curve
 * `b` = a numeric scale parameter on the input axis, the 'growth rate'
 * `g` = the reciprocal of the shape parameter: affects near which asymptote maximum 'growth' occurs
```
5PLM model
Fn = fmax + s*x + (y0 - fmax) / (1 + (2^(1/g)-1) * exp(b*(x-xmid)))^g
```

3 `Quadratic`: parameters of ...
 * test

4 `Bilinear`: parameters of the bilinear model:
 * `a1a` = "slope" describes together with a1b the curve of the first phase of efficiency decline
 * `a1b` = "slope" describes together with a1a the curve of the first phase of efficiency decline
 * `a2` = "slope" represents the steepness of the second phase of efficiency decline
 * `ipt` = corresponds to the horizontal position of the phase change
 * `eta` = affects the abruptness of transition between the two phases
 * `chi` = affects the y-intercept of the curve
```
Bilinear model:
Ln²(En) = chi + eta *log(exp(( a1a * (input-ipt)^2+ a1b *(input-ipt))/eta) + exp( a2 *(n-ipt)/eta))
```  

`ouptut = "all"` returns all of the above parameters in a single vector (length = 20)

### Example
```
react<-c(193.6422, 209.2176,  221.9136,  233.8279,  243.1534,  245.9008,  251.7054,  253.5101,  260.7267,  258.2353,  261.5051,  265.9691,  268.3236,  273.0342,  274.9183, 281.7958,  291.9236,  312.8073,  353.2469,  425.0797,  558.4097,  801.8603, 1141.0918, 1617.8774, 2179.0478, 2796.2294, 3426.8560, 4029.0126, 4564.6015, 5065.4347, 5478.7222, 5865.2510, 6197.9947, 6502.7750, 6783.9353, 6997.8949, 7211.0589, 7422.2636, 7626.5014, 7732.2221, 7820.8930, 7925.0876, 7991.5905, 8029.5963, 8045.4035, 8045.6049, 8064.2315, 8103.9244, 8120.5715, 8104.7643, 8117.3816, 8090.2008, 8132.5172, 8078.8213, 8081.2783, 8058.4114, 8021.9082, 8020.0389, 7986.3723, 7999.9835)

analyse(react)
	
analyse(react,plots=TRUE)

```

### Disclaimer 
This function has only been tested on a limited amount of data, there may be a few bugs left.

Most problems are caused by a lack of data resulting in a suboptimal fit try running more cycles first! 60 cycles will do in most cases. Increasing primer concentration may help to increase the length of the firs phase of efficiency decline resulting in a better fit. 
