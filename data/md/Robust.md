## CRAN Task View: Robust Statistical Methods

  ----------------- ------------------------------------------
  **Maintainer:**   Martin Maechler
  **Contact:**      Martin.Maechler at R-project.org
  **Version:**      2016-08-29
  **URL:**          <https://CRAN.R-project.org/view=Robust>
  ----------------- ------------------------------------------

<div>

Robust (or \"resistant\") methods for statistics modelling have been
available in S from the very beginning in the 1980s; and then in R in
package `stats`. Examples are `median()`, `mean(*, trim =. )`, `mad()`,
`IQR()`, or also `fivenum()`, the statistic behind `boxplot()` in
package `graphics`) or `lowess()` (and `loess()`) for robust
nonparametric regression, which had been complemented by `runmed()` in
2003. Much further important functionality has been made available in
recommended (and hence present in all R versions) package
[MASS](../packages/MASS/index.html) (by Bill Venables and Brian Ripley,
see *the* book [Modern Applied Statistics with
S](http://www.stats.ox.ac.uk/pub/MASS4/) ). Most importantly, they
provide `rlm()` for robust regression and `cov.rob()` for robust
multivariate scatter and covariance.

This task view is about R add-on packages providing newer or faster,
more efficient algorithms and notably for (robustification of) new
models.

Please send suggestions for additions and extensions to the [task view
maintainer](mailto:maechler@R-project.org) .

An international group of scientists working in the field of robust
statistics has made efforts (since October 2005) to coordinate several
of the scattered developments and make the important ones available
through a set of R packages complementing each other. These should build
on a basic package with \"Essentials\", coined
[robustbase](../packages/robustbase/index.html) with (potentially many)
other packages building on top and extending the essential functionality
to particular models or applications. Further, there is the quite
comprehensive package [robust](../packages/robust/index.html), a version
of the robust library of S-PLUS, as an R package now GPLicensed thanks
to Insightful and Kjell Konis. Originally, there has been much overlap
between \'robustbase\' and \'robust\', now
[robust](../packages/robust/index.html) *depends* on
[robustbase](../packages/robustbase/index.html), the former providing
convenient routines for the casual user where the latter will contain
the underlying functionality, and provide the more advanced statistician
with a large range of options for robust modeling.

We structure the packages roughly into the following topics, and
typically will first mention functionality in packages
[robustbase](../packages/robustbase/index.html) and
[robust](../packages/robust/index.html).

-   *Regression (Linear, Generalized Linear, Nonlinear Models, incl.
    Mixed Effects)* : `lmrob()`
    ([robustbase](../packages/robustbase/index.html)) and `lmRob()`
    ([robust](../packages/robust/index.html)) where the former uses the
    latest of the fast-S algorithms and heteroscedasticity and
    autocorrelation corrected (HAC) standard errors, the latter makes
    use of the M-S algorithm of Maronna and Yohai (2000), automatically
    when there are factors among the predictors (where S-estimators (and
    hence MM-estimators) based on resampling typically badly fail). The
    `ltsReg()` and `lmrob.S()` functions are available in
    [robustbase](../packages/robustbase/index.html), but rather for
    comparison purposes. `rlm()` from
    [MASS](../packages/MASS/index.html) had been the first widely
    available implementation for robust linear models, and also one of
    the very first MM-estimation implementations.
    [robustreg](../packages/robustreg/index.html) provides very simple
    M-estimates for linear regression (in pure R). Note that Koenker\'s
    quantile regression package
    [quantreg](../packages/quantreg/index.html) contains L1 (aka LAD,
    least absolute deviations)-regression as a special case, doing so
    also for nonparametric regression via splines. Quantile regression
    (and hence L1 or LAD) for mixed effect models, is available in
    package [lqmm](../packages/lqmm/index.html), whereas an *MM-like*
    approach for robust linear **mixed effects** modeling is available
    from package [robustlmm](../packages/robustlmm/index.html). Package
    [mblm](../packages/mblm/index.html) \'s function `mblm()` fits
    median-based (Theil-Sen or Siegel\'s repeated) simple linear models.
    Package [TEEReg](../packages/TEEReg/index.html) provides trimmed
    elemental estimators for linear models. Generalized linear models
    (GLMs) are provided both via `glmrob()`
    ([robustbase](../packages/robustbase/index.html)) and `glmRob()`
    ([robust](../packages/robust/index.html)), where package
    [robustloggamma](../packages/robustloggamma/index.html) focuses on
    generalized log gamma models. Robust ordinal regression is provided
    by [rorutadis](../packages/rorutadis/index.html) (UTADIS). Robust
    Nonlinear model fitting is available through
    [robustbase](../packages/robustbase/index.html) \'s `nlrob()`.
    [multinomRob](../packages/multinomRob/index.html) fits overdispersed
    multinomial regression models for count data.
    [rgam](../packages/rgam/index.html) and
    [robustgam](../packages/robustgam/index.html) both fit robust GAMs,
    i.e., robust Generalized Additive Models.
    [drgee](../packages/drgee/index.html) fits \"Doubly Robust\"
    Generalized Estimating Equations (GEEs)
-   *Multivariate Analysis* : Here, the
    [rrcov](../packages/rrcov/index.html) package which builds (\"
    `Depends` \") on [robustbase](../packages/robustbase/index.html)
    provides nice S4 class based methods, more methods for robust
    multivariate variance-covariance estimation, and adds robust PCA
    methodology. It is extended by
    [rrcovNA](../packages/rrcovNA/index.html), providing robust
    multivariate methods for *for incomplete* or missing ( `NA`) data,
    and by [rrcovHD](../packages/rrcovHD/index.html), providing robust
    multivariate methods for *High Dimensional* data. High dimensional
    data with an emphasis on functional data are treated robustly also
    by [roahd](../packages/roahd/index.html). Here,
    [robustbase](../packages/robustbase/index.html) contains a slightly
    more flexible version, `covMcd()` than
    [robust](../packages/robust/index.html) \'s `fastmcd()`, and
    similarly for `covOGK()`. OTOH,
    [robust](../packages/robust/index.html) \'s `covRob()` has
    automatically chosen methods, notably `pairwiseQC()` for large
    dimensionality p. Package [robustX](../packages/robustX/index.html)
    for experimental, or other not yet established procedures, contains
    `BACON()` and `covNCC()`, the latter providing the neighbor variance
    estimation (NNVE) method of Wang and Raftery (2002), also available
    (slightly less optimized) in
    [covRobust](../packages/covRobust/index.html).
    [RobRSVD](../packages/RobRSVD/index.html) provides a robust
    Regularized Singular Value Decomposition.
    [mvoutlier](../packages/mvoutlier/index.html) (building on
    [robustbase](../packages/robustbase/index.html)) provides several
    methods for outlier identification in high dimensions.
    [GSE](../packages/GSE/index.html) estimates multivariate location
    and scatter in the presence of missing data.
    [FRB](../packages/FRB/index.html) performs robust inference based on
    **F** ast and **R** obust **B** ootstrap on robust estimators,
    including multivariate regression, PCA and Hotelling tests.
    [RSKC](../packages/RSKC/index.html) provides **R** obust **S** parse
    **K** -means **C** lustering.
    [robustDA](../packages/robustDA/index.html) for *robust mixture
    Discriminant Analysis* (RMDA) builds a mixture model classifier with
    noisy class labels. [robcor](../packages/robcor/index.html) computes
    robust pairwise correlations based on scale estimates, particularly
    on `FastQn()`. [covRobust](../packages/covRobust/index.html)
    provides the nearest neighbor variance estimation (NNVE) method of
    Wang and Raftery (2002). Note that robust PCA can be performed by
    using standard R\'s `princomp()`, e.g.,
    `X <- stackloss; pc.rob <- princomp(X, covmat= MASS::cov.rob(X))`
    See also the CRAN task views [Multivariate](Multivariate.html) and
    [Cluster](Cluster.html)
-   *Large Data Sets* : `BACON()` (in
    [robustX](../packages/robustX/index.html)) should be applicable for
    larger (n,p) than traditional robust covariance based outlier
    detectors. [OutlierDM](../packages/OutlierDM/index.html) detects
    outliers for replicated high-throughput data. (See also the CRAN
    task view [MachineLearning](MachineLearning.html).)
-   *Descriptive Statistics / Exploratory Data Analysis* :
    `boxplot.stats()`, etc mentioned above
-   *Time Series* :
    -   R\'s `runmed()` provides *most robust* running median filtering.
    -   Package [robfilter](../packages/robfilter/index.html) contains
        robust regression and filtering methods for univariate time
        series, typically based on repeated (weighted) median
        regressions.
    -   The [RobPer](../packages/RobPer/index.html) provides several
        methods for robust periodogram estimation, notably for
        irregularly spaced time series.
    -   Peter Ruckdeschel has started to lead an effort for a robust
        time-series package, see
        [[robust-ts]{.Rforge}](https://R-Forge.R-project.org/projects/robust-ts/)
        on R-Forge.
    -   Further, robKalman, *\"Routines for Robust Kalman Filtering
        \-\-- the ACM- and rLS-filter\"* , is being developed, see
        [[robkalman]{.Rforge}](https://R-Forge.R-project.org/projects/robkalman/)
        on R-Forge.

    Note however that these (last two items) are not yet available from
    CRAN.
-   *Econometric Models* : Econometricians tend to like HAC
    (heteroscedasticity and autocorrelation corrected) standard errors.
    For a broad class of models, these are provided by package
    [sandwich](../packages/sandwich/index.html). Note that
    `vcov(lmrob())` also uses a version of HAC standard errors for its
    robustly estimated linear models. See also the CRAN task view
    [Econometrics](Econometrics.html)
-   *Robust Methods for Bioinformatics* : There are several packages in
    the [Bioconductor project](http://www.bioconductor.org/) providing
    specialized robust methods. In addition,
    [RobLoxBioC](../packages/RobLoxBioC/index.html) provides
    infinitesimally robust estimators for preprocessing omics data.
-   *Robust Methods for Survival Analysis* : Package
    [coxrobust](../packages/coxrobust/index.html) provides robust
    estimation in the Cox model; package
    [rsig](../packages/rsig/index.html) does robust signature selection
    for survival outcomes. [OutlierDC](../packages/OutlierDC/index.html)
    detects outliers using quantile regression for censored data.
-   *Robust Methods for Surveys* : On R-forge only, package
    [[rhte]{.Rforge}](https://R-Forge.R-project.org/projects/rhte/)
    provides a robust Horvitz-Thompson estimator.
-   *Geostatistics* : Package [georob](../packages/georob/index.html)
    aims at robust geostatistical analysis of spatial data, such as
    kriging and more.
-   *Collections of several methodologies* :
    -   [WRS2](../packages/WRS2/index.html) contains robust tests for
        ANOVA and ANCOVA from Rand Wilcox\'s collection.
    -   [robeth](../packages/robeth/index.html) contains R functions
        interfacing to the extensive RobETH fortran library with many
        functions for regression, multivariate estimation and more.
-   *Other approaches to robust and resistant methodology* :
    -   The package [distr](../packages/distr/index.html) and its
        several child packages also allow to explore robust estimation
        concepts, see e.g.,
        [[distr]{.Rforge}](https://R-Forge.R-project.org/projects/distr/)
        on R-Forge.
    -   Notably, based on these, the project
        [[robast]{.Rforge}](https://R-Forge.R-project.org/projects/robast/)
        aims for the implementation of R packages for the computation of
        optimally robust estimators and tests as well as the necessary
        infrastructure (mainly S4 classes and methods) and diagnostics;
        cf. M. Kohl (2005). It includes the R packages
        [RandVar](../packages/RandVar/index.html),
        [RobAStBase](../packages/RobAStBase/index.html),
        [RobLox](../packages/RobLox/index.html),
        [RobLoxBioC](../packages/RobLoxBioC/index.html),
        [RobRex](../packages/RobRex/index.html). Further,
        [ROptEst](../packages/ROptEst/index.html), and
        [ROptRegTS](../packages/ROptRegTS/index.html).
    -   [RobustAFT](../packages/RobustAFT/index.html) computes Robust
        Accelerated Failure Time Regression for Gaussian and logWeibull
        errors.
    -   [wle](../packages/wle/index.html) *Weighted Likelihood
        Estimation* provides robustified likelihood estimation for a
        range of models, notably (generalized) regression, and time
        series (AR and fracdiff).
    -   [robumeta](../packages/robumeta/index.html) for robust variance
        meta-regression
    -   [ssmrob](../packages/ssmrob/index.html) provides robust
        estimation and inference in sample selection models.

</div>

### CRAN packages:

-   [covRobust](../packages/covRobust/index.html)
-   [coxrobust](../packages/coxrobust/index.html)
-   [distr](../packages/distr/index.html)
-   [drgee](../packages/drgee/index.html)
-   [FRB](../packages/FRB/index.html)
-   [georob](../packages/georob/index.html)
-   [GSE](../packages/GSE/index.html)
-   [lqmm](../packages/lqmm/index.html)
-   [MASS](../packages/MASS/index.html) (core)
-   [mblm](../packages/mblm/index.html)
-   [multinomRob](../packages/multinomRob/index.html)
-   [mvoutlier](../packages/mvoutlier/index.html)
-   [OutlierDC](../packages/OutlierDC/index.html)
-   [OutlierDM](../packages/OutlierDM/index.html)
-   [quantreg](../packages/quantreg/index.html)
-   [RandVar](../packages/RandVar/index.html)
-   [rgam](../packages/rgam/index.html)
-   [roahd](../packages/roahd/index.html)
-   [RobAStBase](../packages/RobAStBase/index.html)
-   [robcor](../packages/robcor/index.html)
-   [robeth](../packages/robeth/index.html)
-   [robfilter](../packages/robfilter/index.html)
-   [RobLox](../packages/RobLox/index.html)
-   [RobLoxBioC](../packages/RobLoxBioC/index.html)
-   [RobPer](../packages/RobPer/index.html)
-   [RobRex](../packages/RobRex/index.html)
-   [RobRSVD](../packages/RobRSVD/index.html)
-   [robumeta](../packages/robumeta/index.html)
-   [robust](../packages/robust/index.html) (core)
-   [RobustAFT](../packages/RobustAFT/index.html)
-   [robustbase](../packages/robustbase/index.html) (core)
-   [robustDA](../packages/robustDA/index.html)
-   [robustgam](../packages/robustgam/index.html)
-   [robustlmm](../packages/robustlmm/index.html)
-   [robustloggamma](../packages/robustloggamma/index.html)
-   [robustreg](../packages/robustreg/index.html)
-   [robustX](../packages/robustX/index.html)
-   [ROptEst](../packages/ROptEst/index.html)
-   [ROptRegTS](../packages/ROptRegTS/index.html)
-   [rorutadis](../packages/rorutadis/index.html)
-   [rrcov](../packages/rrcov/index.html) (core)
-   [rrcovHD](../packages/rrcovHD/index.html)
-   [rrcovNA](../packages/rrcovNA/index.html)
-   [rsig](../packages/rsig/index.html)
-   [RSKC](../packages/RSKC/index.html)
-   [sandwich](../packages/sandwich/index.html)
-   [ssmrob](../packages/ssmrob/index.html)
-   [TEEReg](../packages/TEEReg/index.html)
-   [wle](../packages/wle/index.html)
-   [WRS2](../packages/WRS2/index.html)

### Related links:

-   Mailing list:
-   [R Special Interest Group on Robust
    Statistics](https://stat.ethz.ch/mailman/listinfo/r-sig-robust/)
-   [Robust Statistics in R (TU Vienna)](http://cstat.tuwien.ac.at/rsr/)
-   R-Forge Project:
    [[distr]{.Rforge}](https://R-Forge.R-project.org/projects/distr/)
-   R-Forge Project:
    [[robast]{.Rforge}](https://R-Forge.R-project.org/projects/robast/)
-   R-Forge Project:
    [[robkalman]{.Rforge}](https://R-Forge.R-project.org/projects/robkalman/)
-   R-Forge Project:
    [[robust-ts]{.Rforge}](https://R-Forge.R-project.org/projects/robust-ts/)
