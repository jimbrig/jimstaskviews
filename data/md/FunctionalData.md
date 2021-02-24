## CRAN Task View: Functional Data Analysis

  ----------------- --------------------------------------------------
  **Maintainer:**   Fabian Scheipl
  **Contact:**      fabian.scheipl at stat.uni-muenchen.de
  **Version:**      2018-02-12
  **URL:**          <https://CRAN.R-project.org/view=FunctionalData>
  ----------------- --------------------------------------------------

<div>

Functional data analysis (FDA) deals with data that [\"provides
information about curves, surfaces or anything else varying over a
continuum.\"](https://en.wikipedia.org/wiki/Functional_data_analysis)
This task view catalogues available packages in this rapidly developing
field.

**General functional data analysis**

-   [fda](../packages/fda/index.html) provides functions to enable all
    aspects of functional data analysis: It includes object-types for
    functional data with corresponding functions for smoothing, plotting
    and regression models. The package includes data sets and script
    files for working examples from the book: Ramsay, J. O., Hooker,
    Giles, and Graves, Spencer (2009) \"Data Analysis with R and
    Matlab\" (Springer).
-   [fdasrvf](../packages/fdasrvf/index.html) performs alignment, PCA,
    and regression of multidimensional or unidimensional functions using
    the square-root velocity framework (Srivastava et al., 2011). This
    framework allows for elastic analysis of functional data through
    phase and amplitude separation.
-   [fdapace](../packages/fdapace/index.html) provides functional
    principal component based methods for sparsely or densely sampled
    random trajectories and time courses for functional regression and
    correlation, for longitudinal data analysis, the analysis of
    stochastic processes from samples of realized trajectories, and for
    the analysis of underlying dynamics.
-   [fda.usc](../packages/fda.usc/index.html) provides routines for
    exploratory and descriptive analysis of functional data such as
    depth measurements, outlier detection, as well as unsupervised and
    supervised classification, (univariate, nonparametric) regression
    models with a functional covariate and functional analysis of
    variance.
-   [funData](../packages/funData/index.html) provides S4 classes for
    univariate and multivariate functional and image data and utility
    functions.
-   [fds](../packages/fds/index.html) contains 19 data sets with
    functional data.
-   [rainbow](../packages/rainbow/index.html) contains functions and
    data sets for functional data display, exploratory analysis and
    outlier detection.
-   [roahd](../packages/roahd/index.html) provides methods for the
    robust analysis of univariate and multivariate functional data,
    possibly in high-dimensional cases, and hence with attention to
    computational efficiency and simplicity of use.

**Regression and classification for functional data**

-   [denseFLMM](../packages/denseFLMM/index.html) estimates functional
    linear mixed models for densely sampled data based on functional
    principal component analysis.
-   [GPFDA](../packages/GPFDA/index.html) uses functional regression as
    the mean structure and Gaussian processes as the covariance
    structure.
-   [growfunctions](../packages/growfunctions/index.html) estimates a
    collection of time-indexed functions under either of Gaussian
    process (GP) or intrinsic Gaussian Markov random field (iGMRF) prior
    formulations where a Dirichlet process mixture allows sub-groupings
    of the functions to share the same covariance or precision
    parameters. The GP and iGMRF formulations both support any number of
    additive covariance or precision terms, respectively, expressing
    either or both of multiple trend and seasonality.
-   [refund](../packages/refund/index.html) provides spline-based
    methods for roughness penalized function-on-scalar,
    scalar-on-function, and function-on-function regression as well as
    methods for functional PCA. Some of the functions are applicable to
    image data.
-   [refund.wave](../packages/refund.wave/index.html) provides methods
    for regressing scalar responses on functional or image predictors,
    via transformation to the wavelet domain and back.
-   [refund.shiny](../packages/refund.shiny/index.html) provides
    interactive plots for functional data analyses.
-   [FDboost](../packages/FDboost/index.html) implements flexible
    additive regression models and variable selection for
    scalar-on-function, function-on-scalar and function-on-function
    regression models that are fitted by a component-wise gradient
    boosting algorithm.
-   [fdaPDE](../packages/fdaPDE/index.html) contains an implementation
    of regression models with partial differential regularizations.
-   [flars](../packages/flars/index.html) implements variable selection
    for the functional linear regression with scalar response variable
    and mixed scalar/functional predictors based on the least angle
    regression approach.
-   [sparseFLMM](../packages/sparseFLMM/index.html) implements
    functional linear mixed models for irregularly or sparsely sampled
    data based on functional principal component analysis.
-   [dbstats](../packages/dbstats/index.html) provides prediction
    methods where explanatory information is coded as a matrix of
    distances between individuals. It includes distance based versions
    of `             lm           ` and `             glm,           `
    as well as nonparametric versions of both, based on local
    estimation. To apply these methods to functional data it is
    sufficient to calculate a distance matrix between the observed
    functional data.
-   [classiFunc](../packages/classiFunc/index.html) provides nearest
    neighbor and kernel-based estimation based on semimetrics for
    supervised classification of functional data.

**Clustering functional data**

-   [Funclustering](../packages/Funclustering/index.html) implements a
    model-based clustering algorithm for multivariate functional data.
-   [funFEM](../packages/funFEM/index.html) \'s algorithm (Bouveyron et
    al., 2014) allows to cluster functional data by modeling the curves
    within a common and discriminative functional subspace.
-   [funHDDC](../packages/funHDDC/index.html) provides the funHDDC
    algorithm (Bouveyron & Jacques, 2011) which allows to cluster
    functional data by modeling each group within a specific functional
    subspace.
-   [funcy](../packages/funcy/index.html) provides a unified framework
    to cluster functional data according to one of seven models. All
    models are based on the projection of the curves onto a basis.
    Method specific as well as general visualization tools are
    available.
-   [fdakma](../packages/fdakma/index.html) performs clustering and
    alignment of a multidimensional or unidimensional functional dataset
    by means of k-mean alignment.

**Registering and aligning functional data**

-   [fdasrvf](../packages/fdasrvf/index.html) performs alignment, PCA,
    and regression of multidimensional or unidimensional functions using
    the square-root velocity framework (Srivastava et al., 2011). This
    framework allows for elastic analysis of functional data through
    phase and amplitude separation.
-   [warpMix](../packages/warpMix/index.html) implements warping
    (alignment) for functional data using B-spline based mixed effects
    models.
-   [fdakma](../packages/fdakma/index.html) performs clustering and
    alignment of a multidimensional or unidimensional functional dataset
    by means of k-mean alignment.

**Time series of functional data**

-   [ftsa](../packages/ftsa/index.html) provides functions for
    visualizing, modeling, forecasting and hypothesis testing of
    functional time series.
-   [ftsspec](../packages/ftsspec/index.html) provides functions for
    estimating the spectral density operator of functional time series
    (FTS) and comparing the spectral density operator of two functional
    time series, in a way that allows detection of differences of the
    spectral density operator in frequencies and along the curve length.
-   [freqdom](../packages/freqdom/index.html) provides frequency domain
    methods for multivariate and functional time series and implements
    dynamic functional principal components and functional regression in
    the presence of temporal dependence.
-   [freqdom.fda](../packages/freqdom.fda/index.html) provides a wrapper
    for functionality of [freqdom](../packages/freqdom/index.html) for
    objects from [fda](../packages/fda/index.html)
-   [pcdpca](../packages/pcdpca/index.html) extends multivariate dynamic
    principal components to periodically correlated multivariate and
    functional time series.

**Other**

-   [ddalpha](../packages/ddalpha/index.html) implements depth-based
    classification and calculation of data depth, also for functional
    data.
-   [fpca](../packages/fpca/index.html) implements functional principal
    components for sparsely observed data in a geometric approach to MLE
    for functional principal components.
-   [fdatest](../packages/fdatest/index.html) provides an implementation
    of the Interval Testing Procedure for functional data in different
    frameworks (i.e., one or two-population frameworks, functional
    linear models) by means of different basis expansions (i.e.,
    B-spline, Fourier, and phase-amplitude Fourier).
-   [fdadensity](../packages/fdadensity/index.html) implements Petersen
    and Mueller (2016) (doi:10.1214/15-AOS1363) for the analysis of
    samples of density functions via specialized Functional Principal
    Components Analysis.
-   [geofd](../packages/geofd/index.html) provides Kriging based methods
    for predicting functional data (curves) with spatial dependence.
-   [RFgroove](../packages/RFgroove/index.html) implements variable
    selection tools for groups of variables and functional data based on
    a new grouped variable importance with random forests, implementing
    Gregorutti, B., Michel, B. and Saint Pierre, P. (2015). Grouped
    variable importance with random forests and application to multiple
    functional data analysis, *Computational Statistics and Data
    Analysis* **90** , 15-35.
-   [switchnpreg](../packages/switchnpreg/index.html) provides functions
    for estimating the parameters from the latent state process and the
    functions corresponding to the J states as proposed by De Souza and
    Heckman (2013).
-   [fdcov](../packages/fdcov/index.html) provides a variety of tools
    for the analysis of covariance operators.
-   [covsep](../packages/covsep/index.html) provides functions for
    testing if the covariance structure of 2-dimensional data is
    separable.

The Functional Data Analysis Task View is written by Fabian Scheipl,
Sonja Greven and Tore Erdmann (LMU München, Germany). Please contact
Fabian Scheipl with suggestions, additions and improvements.

</div>

### CRAN packages:

-   [classiFunc](../packages/classiFunc/index.html)
-   [covsep](../packages/covsep/index.html)
-   [dbstats](../packages/dbstats/index.html)
-   [ddalpha](../packages/ddalpha/index.html)
-   [denseFLMM](../packages/denseFLMM/index.html)
-   [fda](../packages/fda/index.html) (core)
-   [fda.usc](../packages/fda.usc/index.html) (core)
-   [fdadensity](../packages/fdadensity/index.html)
-   [fdakma](../packages/fdakma/index.html)
-   [fdapace](../packages/fdapace/index.html) (core)
-   [fdaPDE](../packages/fdaPDE/index.html)
-   [fdasrvf](../packages/fdasrvf/index.html) (core)
-   [fdatest](../packages/fdatest/index.html)
-   [FDboost](../packages/FDboost/index.html) (core)
-   [fdcov](../packages/fdcov/index.html)
-   [fds](../packages/fds/index.html)
-   [flars](../packages/flars/index.html)
-   [fpca](../packages/fpca/index.html)
-   [freqdom](../packages/freqdom/index.html)
-   [freqdom.fda](../packages/freqdom.fda/index.html)
-   [ftsa](../packages/ftsa/index.html) (core)
-   [ftsspec](../packages/ftsspec/index.html)
-   [Funclustering](../packages/Funclustering/index.html)
-   [funcy](../packages/funcy/index.html) (core)
-   [funData](../packages/funData/index.html)
-   [funFEM](../packages/funFEM/index.html)
-   [funHDDC](../packages/funHDDC/index.html)
-   [geofd](../packages/geofd/index.html)
-   [GPFDA](../packages/GPFDA/index.html)
-   [growfunctions](../packages/growfunctions/index.html)
-   [pcdpca](../packages/pcdpca/index.html)
-   [rainbow](../packages/rainbow/index.html)
-   [refund](../packages/refund/index.html) (core)
-   [refund.shiny](../packages/refund.shiny/index.html)
-   [refund.wave](../packages/refund.wave/index.html)
-   [RFgroove](../packages/RFgroove/index.html)
-   [roahd](../packages/roahd/index.html)
-   [sparseFLMM](../packages/sparseFLMM/index.html)
-   [switchnpreg](../packages/switchnpreg/index.html)
-   [warpMix](../packages/warpMix/index.html)

### Related links:

-   [Website of the canonical FDA book by Ramsay and
    Silverman](http://www.psych.mcgill.ca/misc/fda/)
-   [PACE: collection of MATLAB scripts from UC
    Davis](http://www.stat.ucdavis.edu/PACE/)
-   [WFMM: powerful software for Bayesian wavelet-based functional mixed
    models
    (C++/Matlab)](https://biostatistics.mdanderson.org/softwaredownload/SingleSoftware.aspx?Software_Id=70)
