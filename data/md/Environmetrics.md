## CRAN Task View: Analysis of Ecological and Environmental Data

  ----------------- --------------------------------------------------
  **Maintainer:**   Gavin Simpson
  **Contact:**      ucfagls at gmail.com
  **Version:**      2018-04-11
  **URL:**          <https://CRAN.R-project.org/view=Environmetrics>
  ----------------- --------------------------------------------------

<div>

#### Introduction

This Task View contains information about using R to analyse ecological
and environmental data.

The base version of R ships with a wide range of functions for use
within the field of environmetrics. This functionality is complemented
by a plethora of packages available via CRAN, which provide specialist
methods such as ordination & cluster analysis techniques. A brief
overview of the available packages is provided in this Task View,
grouped by topic or type of analysis. As a testament to the popularity
of R for the analysis of environmental and ecological data, a [special
volume](http://www.jstatsoft.org/v22/) of the *Journal of Statistical
Software* was produced in 2007.

Those useRs interested in environmetrics should consult the
[Spatial](Spatial.html) view. Complementary information is also
available in the [Multivariate](Multivariate.html),
[Phylogenetics](Phylogenetics.html), [Cluster](Cluster.html), and
[SpatioTemporal](SpatioTemporal.html) task views.

If you have any comments or suggestions for additions or improvements,
then please contact the
[maintainer](mailto:ucfagls@gmail.com?subject=Environmetrics%20Task%20View)
.

A list of available packages and functions is presented below, grouped
by analysis type.

#### General packages

These packages are general, having wide applicability to the
environmetrics field.

-   Package [EnvStats](../packages/EnvStats/index.html) is the successor
    to the S-PLUS module *EnvironmentalStats* , both by Steven Millard.
    A [user guide in the form of a
    book](https://dx.doi.org/10.1007/978-1-4614-8456-1) has recently be
    released.

#### Modelling species responses and other data

Analysing species response curves or modeling other data often involves
the fitting of standard statistical models to ecological data and
includes simple (multiple) regression, Generalised Linear Models (GLM),
extended regression (e.g. Generalised Least Squares \[GLS\]),
Generalised Additive Models (GAM), and mixed effects models, amongst
others.

-   The base installation of R provides `lm()` and `glm()` for fitting
    linear and generalised linear models, respectively.
-   Generalised least squares and linear and non-linear mixed effects
    models extend the simple regression model to account for clustering,
    heterogeneity and correlations within the sample of observations.
    Package [nlme](../packages/nlme/index.html) provides functions for
    fitting these models. The package is supported by Pinheiro & Bates
    (2000) *Mixed-effects Models in S and S-PLUS* , Springer, New York.
    An updated approach to mixed effects models, which also fits
    Generalised Linear Mixed Models (GLMM) and Generalised non-Linear
    Mixed Models (GNLMM) is provided by the
    [lme4](../packages/lme4/index.html) package, though this is
    currently beta software and does not yet allow correlations within
    the error structure.
-   Recommended package [mgcv](../packages/mgcv/index.html) fits GAMs
    and Generalised Additive Mixed Models (GAMM) with automatic
    smoothness selection via generalised cross-validation. The author of
    [mgcv](../packages/mgcv/index.html) has also written a companion
    monograph, Wood (2006) *Generalized Additive Models; An Introduction
    with R* Chapman Hall/CRC, which has an accompanying package
    [gamair](../packages/gamair/index.html).
-   Alternatively, package [gam](../packages/gam/index.html) provides an
    implementation of the S-PLUS function `gam()` that includes LOESS
    smooths.
-   Proportional odds models for ordinal responses can be fitted using
    `polr()` in the [MASS](../packages/MASS/index.html) package, of Bill
    Venables and Brian Ripley.
-   A negative binomial family for GLMs to model over-dispersion in
    count data is available in [MASS](../packages/MASS/index.html).
-   Models for overdispersed counts and proportions
    -   Package [pscl](../packages/pscl/index.html) also contains
        several functions for dealing with over-dispersed count data.
        Poisson or negative binomial distributions are provided for both
        zero-inflated and hurdle models.
    -   [aod](../packages/aod/index.html) provides a suite of functions
        to analyse overdispersed counts or proportions, plus utility
        functions to calculate e.g. AIC, AICc, Akaike weights.
-   [Detecting change points and structural changes in parametric models
    is well catered for in the
    [segmented](../packages/segmented/index.html) package and the
    [strucchange](../packages/strucchange/index.html) package
    respectively. [segmented](../packages/segmented/index.html) has
    recently been the subject of an R News article ( [R News, volume 8
    issue
    1](http://CRAN.R-project.org/doc/Rnews/Rnews_2008-1.pdf "Link to segment article in R News")
    ).]{#changepoints}

#### Tree-based models

Tree-based models are being increasingly used in ecology, particularly
for their ability to fit flexible models to complex data sets and the
simple, intuitive output of the tree structure. Ensemble methods such as
bagging, boosting and random forests are advocated for improving
predictions from tree-based models and to provide information on
uncertainty in regression models or classifiers.

Tree-structured models for regression, classification and survival
analysis, following the ideas in the CART book, are implemented in

-   recommended package [rpart](../packages/rpart/index.html)
-   [party](../packages/party/index.html) provides an implementation of
    conditional inference trees which embed tree-structured regression
    models into a well defined theory of conditional inference
    procedures

Multivariate trees are available in

-   package [party](../packages/party/index.html) can also handle
    multivariate responses.

Ensemble techniques for trees:

-   The Random Forest method of Breiman and Cutler is implemented in
    [randomForest](../packages/randomForest/index.html), providing
    classification and regression based on a forest of trees using
    random inputs
-   Package [ipred](../packages/ipred/index.html) provides functions for
    improved predictive models for classification, regression and
    survival problems.

Graphical tools for the visualization of trees are available in package
[maptree](../packages/maptree/index.html).

Packages [mda](../packages/mda/index.html) and
[earth](../packages/earth/index.html) implement Multivariate Adaptive
Regression Splines (MARS), a technique which provides a more flexible,
tree-based approach to regression than the piecewise constant functions
used in regression trees.

#### Ordination

R and add-on packages provide a wide range of ordination methods, many
of which are specialised techniques particularly suited to the analysis
of species data. The two main packages are
[ade4](../packages/ade4/index.html) and
[vegan](../packages/vegan/index.html).
[ade4](../packages/ade4/index.html) derives from the traditions of the
French school of "Analyse des Donnees" and is based on the use of the
duality diagram. [vegan](../packages/vegan/index.html) follows the
approach of Mark Hill, Cajo ter Braak and others, though the
implementation owes more to that presented in Legendre & Legendre (1988)
*Numerical Ecology, 2 ^nd^ English Edition* , Elsevier. Where the two
packages provide duplicate functionality, the user should choose
whichever framework that best suits their background.

-   Principal Components (PCA) is available via the `prcomp()` function.
    `rda()` (in package [vegan](../packages/vegan/index.html)), `pca()`
    (in package [labdsv](../packages/labdsv/index.html)) and
    `dudi.pca()` (in package [ade4](../packages/ade4/index.html)),
    provide more ecologically-orientated implementations.
-   Redundancy Analysis (RDA) is available via `rda()` in
    [vegan](../packages/vegan/index.html) and `pcaiv()` in
    [ade4](../packages/ade4/index.html).
-   Canonical Correspondence Analysis (CCA) is implemented in `cca()` in
    both [vegan](../packages/vegan/index.html) and
    [ade4](../packages/ade4/index.html).
-   Detrended Correspondence Analysis (DCA) is implemented in
    `decorana()` in [vegan](../packages/vegan/index.html).
-   Principal coordinates analysis (PCO) is implemented in `dudi.pco()`
    in [ade4](../packages/ade4/index.html), `pco()` in
    [labdsv](../packages/labdsv/index.html), `pco()` in
    [ecodist](../packages/ecodist/index.html), and `cmdscale()` in
    package [MASS](../packages/MASS/index.html).
-   Non-Metric multi-Dimensional Scaling (NMDS) is provided by
    `isoMDS()` in package [MASS](../packages/MASS/index.html) and
    `nmds()` in [ecodist](../packages/ecodist/index.html). `nmds()`, a
    wrapper function for `isoMDS()`, is also provided by package
    [labdsv](../packages/labdsv/index.html).
    [vegan](../packages/vegan/index.html) provides helper function
    `metaMDS()` for `isoMDS()`, implementing random starts of the
    algorithm and standardised scaling of the NMDS results. The approach
    adopted by [vegan](../packages/vegan/index.html) with `metaMDS()` is
    the recommended approach for ecological data.
-   Coinertia analysis is available via `coinertia()` and `mcoa()`, both
    in [ade4](../packages/ade4/index.html).
-   Co-correspondence analysis to relate two ecological species data
    matrices is available in
    [cocorresp](../packages/cocorresp/index.html).
-   Canonical Correlation Analysis (CCoA - not to be confused with CCA,
    above) is available in `cancor()` in standard package stats.
-   Procrustes rotation is available in `procrustes()` in
    [vegan](../packages/vegan/index.html) and `procuste()` in
    [ade4](../packages/ade4/index.html), with both
    [vegan](../packages/vegan/index.html) and
    [ade4](../packages/ade4/index.html) providing functions to test the
    significance of the association between ordination configurations
    (as assessed by Procrustes rotation) using permutation/randomisation
    and Monte Carlo methods.
-   Constrained Analysis of Principal Coordinates (CAP), implemented in
    `capscale()` in [vegan](../packages/vegan/index.html), fits
    constrained ordination models similar to RDA and CCA but with any
    any dissimilarity coefficient.
-   Constrained Quadratic Ordination (CQO; formerly known as Canonical
    Gaussian Ordination (CGO)) is a maximum likelihood estimation
    alternative to CCA fit by Quadratic Reduced Rank Vector GLMs.
    Constrained Additive Ordination (CAO) is a flexible alternative to
    CQO which uses Quadratic Reduced Rank Vector GAMs. These methods and
    more are provided in Thomas Yee\'s
    [VGAM](../packages/VGAM/index.html) package.
-   Fuzzy set ordination (FSO), an alternative to CCA/RDA and CAP, is
    available in package [fso](../packages/fso/index.html).
    [fso](../packages/fso/index.html) complements a recent paper on
    fuzzy sets in the journal *Ecology* by Dave Roberts (2008,
    Statistical analysis of multidimensional fuzzy set ordinations.
    *Ecology* **89(5)** , 1246-1260).
-   See also the [Multivariate](Multivariate.html) task view for
    complementary information.

#### Dissimilarity coefficients

Much ecological analysis proceeds from a matrix of dissimilarities
between samples. A large amount of effort has been expended formulating
a wide range of dissimilarity coefficients suitable for ecological data.
A selection of the more useful coefficients are available in R and
various contributed packages.

Standard functions that produce, square, symmetric matrices of pair-wise
dissimilarities include:

-   `dist()` in standard package stats
-   `daisy()` in recommended package
    [cluster](../packages/cluster/index.html)
-   `vegdist()` in [vegan](../packages/vegan/index.html)
-   `dsvdis()` in [labdsv](../packages/labdsv/index.html)
-   `Dist()` in [amap](../packages/amap/index.html)
-   `distance()` in [ecodist](../packages/ecodist/index.html)
-   a suite of functions in [ade4](../packages/ade4/index.html)
-   Package [simba](../packages/simba/index.html) provides functions for
    the calculation of similarity and multiple plot similarity measures
    with binary data (for instance presence/absence species data)

Function `distance()` in package
[analogue](../packages/analogue/index.html) can be used to calculate
dissimilarity between samples of one matrix and those of a second
matrix. The same function can be used to produce pair-wise dissimilarity
matrices, though the other functions listed above are faster.
`distance()` can also be used to generate matrices based on Gower\'s
coefficient for mixed data (mixtures of binary, ordinal/nominal and
continuous variables). Function `daisy()` in package
[cluster](../packages/cluster/index.html) provides a faster
implementation of Gower\'s coefficient for mixed-mode data than
`distance()` if a standard dissimilarity matrix is required. Function
`gowdis()` in package [FD](../packages/FD/index.html) also computes
Gower\'s coefficient and implements extensions to ordinal variables.

#### Cluster analysis

Cluster analysis aims to identify groups of samples within multivariate
data sets. A large range of approaches to this problem have been
suggested, but the main techniques are hierarchical cluster analysis,
partitioning methods, such as *k* -means, and finite mixture models or
model-based clustering. In the machine learning literature, cluster
analysis is an unsupervised learning problem.

The [Cluster](Cluster.html) task view provides a more detailed
discussion of available cluster analysis methods and appropriate R
functions and packages.

Hierarchical cluster analysis:

-   `hclust()` in standard package stats
-   Recommended package [cluster](../packages/cluster/index.html)
    provides functions for cluster analysis following the methods
    described in Kaufman and Rousseeuw (1990) *Finding Groups in data:
    an introduction to cluster analysis* , Wiley, New York
-   `hcluster()` in [amap](../packages/amap/index.html)
-   [pvclust](../packages/pvclust/index.html) is a package for assessing
    the uncertainty in hierarchical cluster analysis. It provides
    approximately unbiased *p* -values as well as bootstrap *p* -values.

Partitioning methods:

-   `kmeans()` in stats provides *k* -means clustering
-   `cmeans()` in [e1071](../packages/e1071/index.html) implements a
    fuzzy version of the *k* -means algorithm
-   Recommended package [cluster](../packages/cluster/index.html) also
    provides functions for various partitioning methodologies.

Mixture models and model-based cluster analysis:

-   [mclust](../packages/mclust/index.html) and
    [flexmix](../packages/flexmix/index.html) provide implementations of
    model-based cluster analysis.
-   [prabclus](../packages/prabclus/index.html) clusters a species
    presence-absence matrix object by calculating an MDS from the
    distances, and applying maximum likelihood Gaussian mixtures
    clustering to the MDS points. The maintainer\'s, Christian Hennig,
    web site contains several publications in ecological contexts that
    use [prabclus](../packages/prabclus/index.html), especially Hausdorf
    & Hennig (2007; [Oikos 116 (2007),
    818-828](https://dx.doi.org/10.1111/j.0030-1299.2007.15661.x) ).

#### Ecological theory

There is a growing number of packages and books that focus on the use of
R for theoretical ecological models.

-   [vegan](../packages/vegan/index.html) provides a wide range of
    functions related to ecological theory, such as diversity indices
    (including the "so-called" Hill\'s numbers \[e.g. Hill\'s N ^2^ \]
    and rarefaction), ranked abundance diagrams, Fisher\'s log series,
    Broken Stick model, Hubbell\'s abundance model, amongst others.
-   The [vegetarian](../packages/vegetarian/index.html) provides the
    diversity measures suggested by Jost ( [2006, Oikos 113(2),
    363-375](https://dx.doi.org/10.1111/j.2006.0030-1299.14714.x) ;
    [2007, Ecology 88(10),
    2427-2439](https://dx.doi.org/10.1890/06-1736.1) ).
-   [untb](../packages/untb/index.html) provides a collection of
    utilities for biodiversity data, including the simulation ecological
    drift under Hubbell\'s Unified Neutral Theory of Biodiversity, and
    the calculation of various diagnostics such as Preston curves.
-   [primer](../packages/primer/index.html) is a support software for
    Stevens ( [2009, *A Primer of Ecology with R* ,
    Springer](http://www.springer.com/life+sci/ecology/book/978-0-387-89881-0)
    ). The package provides a variety of functions for modeling
    ecological data and basic theoretical ecology, including functions
    related to demographic matrix models, metapopulation and source-sink
    models, host-parasitoid and disease models, multiple basins of
    attraction, the storage effect, neutral theory, and diversity
    partitioning.
-   Package [BiodiversityR](../packages/BiodiversityR/index.html)
    provides a GUI for biodiversity and community ecology analysis.
-   Function `betadiver()` in [vegan](../packages/vegan/index.html)
    implements all of the diversity indices reviewed in Koleff et al
    (2003; [Journal of Animal Ecology 72(3),
    367-382](https://dx.doi.org/10.1046/j.1365-2656.2003.00710.x) ).
    `betadiver()` also provides a `plot` method to produce the
    co-occurrence frequency triangle plots of the type found in Koleff
    et al (2003).
-   Function `betadisper()`, also in
    [vegan](../packages/vegan/index.html), implements Marti Anderson\'s
    distance-based test for homogeneity of multivariate dispersions
    (PERMDISP, PERMDISP2), a multivariate analogue of Levene\'s test
    (Anderson 2006; [Biometrics 62,
    245-253](https://dx.doi.org/10.1111/j.1541-0420.2005.00440.x) ).
    Anderson et al (2006; [Ecology Letters 9(6),
    683-693](https://dx.doi.org/10.1111/j.1461-0248.2006.00926.x) )
    demonstrate the use of this approach for measuring beta diversity.
-   The [FD](../packages/FD/index.html) package computes several
    measures of functional diversity indices from multiple traits.

#### Population dynamics

##### Estimating animal abundance and related parameters

This section concerns estimation of population parameters (population
size, density, survival probability, site occupancy etc.) by methods
that allow for incomplete detection. Many of these methods use data on
marked animals, variously called \'capture-recapture\',
\'mark-recapture\' or \'capture-mark-recapture\' data.

-   [Rcapture](../packages/Rcapture/index.html) fits loglinear models to
    estimate population size and survival rate from capture-recapture
    data as described by [Baillargeon and
    Rivest (2007)](http://www.jstatsoft.org/v19/i05) .
-   [secr](../packages/secr/index.html) estimates population density
    given spatially explicit capture-recapture data from traps, passive
    DNA sampling, automatic cameras, sound recorders etc. Models are
    fitted by maximum likelihood. The detection function may be
    halfnormal, exponential, cumulative gamma etc. Density surfaces may
    be fitted. Covariates of density and detection parameters are
    specified via formulae.
-   [SPACECAP](../packages/SPACECAP/index.html) provides a graphical
    interface for fitting a spatially explicit capture-recapture model
    to photographic \'capture\' data by the Bayesian method described in
    Royle et al. ( [2009, Ecology 90:
    3233-3244](https://dx.doi.org/10.1890/08-1481.1) ).
-   [DSpat](../packages/DSpat/index.html) provides analyses of
    line-transect distance sampling data in which the density surface
    and the detection function are estimated simultaneously ( [Johnson
    et al. 2009](https://dx.doi.org/10.1111/j.1541-0420.2009.01265.x) ).
-   [unmarked](../packages/unmarked/index.html) fits hierarchical models
    of occurrence and abundance to data collected on species subject to
    imperfect detection. Examples include single- and multi-season
    occupancy models, binomial mixture models, and hierarchical distance
    sampling models. The data can arise from survey methods such
    temporally replicated counts, removal sampling, double-observer
    sampling, and distance sampling. Parameters governing the state and
    observation processes can be modeled as functions of covariates.
-   Package [RMark](../packages/RMark/index.html) provides a
    formula-based R interface for the MARK package which fits a wide
    variety of capture-recapture models. See the [RMark
    website](http://www.phidot.org/software/mark/rmark/) and a [NOAA
    report](http://www.afsc.noaa.gov/Publications/ProcRpt/PR2013-01.pdf)
    (pdf) for further details.
-   Package [marked](../packages/marked/index.html) provides a framework
    for handling data and analysis for mark-recapture.
    [marked](../packages/marked/index.html) can fit Cormack-Jolly-Seber
    (CJS)and Jolly-Seber (JS) models via maximum likelihood and the CJS
    model via MCMC. Maximum likelihood estimates for the CJS model can
    be obtained using R or via a link to the Automatic Differentiation
    Model Builder software. A [description of the
    package](https://dx.doi.org/10.1111/2041-210X.12065) was published
    in Methods in Ecology and Evolution.
-   [mrds](../packages/mrds/index.html) fits detection functions to
    point and line transect distance sampling survey data (for both
    single and double observer surveys). Abundance can be estimated
    using Horvitz-Thompson-type estimators.
-   [Distance](../packages/Distance/index.html) is a simpler interface
    to [mrds](../packages/mrds/index.html) for single observer distance
    sampling surveys.
-   [dsm](../packages/dsm/index.html) fits [density surface
    models](https://dx.doi.org/10.1111/2041-210X.12105) to
    spatially-referenced distance sampling data. Count data are
    corrected using detection function models fitted using
    [mrds](../packages/mrds/index.html) or
    [Distance](../packages/Distance/index.html). Spatial models are
    constructed as in [mgcv](../packages/mgcv/index.html).

Packages [secr](../packages/secr/index.html) and
[DSpat](../packages/DSpat/index.html) can also be used to simulate data
from their respective models.

See also the [SpatioTemporal](SpatioTemporal.html) task view for
analysis of animal tracking data under *Moving objects, trajectories* .

##### Modelling population growth rates:

-   Package [popbio](../packages/popbio/index.html) can be used to
    construct and analyse age- or stage-specific matrix population
    models.

#### Environmental time series

-   Time series objects in R are created using the `ts()` function,
    though see [tseries](../packages/tseries/index.html) or
    [zoo](../packages/zoo/index.html) below for alternatives.
-   Classical time series functionality is provided by the `ar()`, and
    `arima()` functions in standard package stats for autoregressive
    (AR), moving average (MA), autoregressive moving average (ARMA) and
    integrated ARMA (ARIMA) models.
-   The [forecast](../packages/forecast/index.html) package provides
    methods and tools for displaying and analysing univariate time
    series forecasts including exponential smoothing via state space
    models and automatic ARIMA modelling
-   The [dse](../packages/dse/index.html) package provide a variety of
    more advanced estimation methods and multivariate time series
    analysis.
-   Packages [tseries](../packages/tseries/index.html) and
    [zoo](../packages/zoo/index.html) provide general handling and
    analysis of time series data.
-   Irregular time series can be handled using package
    [zoo](../packages/zoo/index.html) as well as by `irts()` in package
    [tseries](../packages/tseries/index.html).
-   [pastecs](../packages/pastecs/index.html) provides functions
    specifically tailored for the analysis of space-time ecological
    series.
-   [strucchange](../packages/strucchange/index.html) allows for
    testing, dating and monitoring of structural change in linear
    regression relationships.
-   Detecting change points in time series data \-\-- see
    [segmented](../packages/segmented/index.html) [above](#changepoints)
    .
-   The [surveillance](../packages/surveillance/index.html) package
    implements statistical methods for the modeling of and change-point
    detection in time series of counts, proportions and categorical
    data. Focus is on outbreak detection in count data time series.
-   Package [dynlm](../packages/dynlm/index.html) provides a convenient
    interface to fitting time series regressions via ordinary least
    squares
-   Package [dyn](../packages/dyn/index.html) provides a different
    approach to that of [dynlm](../packages/dynlm/index.html), which
    allows time series data to be used with any regression function
    written in the style of lm such as `lm()`, `glm()`, `loess()`,
    `rlm()` and `lqs()` from [MASS](../packages/MASS/index.html),
    `randomForest()` (package
    [randomForest](../packages/randomForest/index.html)), `rq()`
    (package [quantreg](../packages/quantreg/index.html)) amongst
    others, whilst preserving the time series information.
-   The [openair](../packages/openair/index.html) provides numerous
    tools to analyse, interpret and understand air pollution time series
    data
-   The [bReeze](../packages/bReeze/index.html) package is a collection
    of widely used methods to analyse, visualise, and interpret wind
    data. Wind resource analyses can subsequently be combined with
    characteristics of wind turbines to estimate the potential energy
    production.

Additionally, a fuller description of available packages for time series
analysis can be found in the [TimeSeries](TimeSeries.html) task view.

#### Spatial data analysis

See the [Spatial](Spatial.html) CRAN Task View for an overview of
spatial analysis in R.

#### [Extreme values]{#extremes}

[ismev](../packages/ismev/index.html) provides functions for models for
extreme value statistics and is support software for Coles (2001) *An
Introduction to Statistical Modelling of Extreme Values* , Springer, New
York. Other packages for extreme value theory include:

-   [evir](../packages/evir/index.html)
-   [evd](../packages/evd/index.html)
-   [evdbayes](../packages/evdbayes/index.html), which provides a
    Bayesian approach to extreme value theory
-   [extRemes](../packages/extRemes/index.html)
-   The [SpatialExtremes](../packages/SpatialExtremes/index.html)
    provides several approaches for modelling spatial extreme events.

#### Phylogenetics and evolution

Packages specifically tailored for the analysis of phylogenetic and
evolutionary data include:

-   [ape](../packages/ape/index.html)
-   [ouch](../packages/ouch/index.html)

The [Phylogenetics](Phylogenetics.html) task view provides more detailed
coverage of the subject area and related functions within R.

UseRs may also be interested in Paradis (2006) *Analysis of
Phylogenetics and Evolution with R* , Springer, New York, a book in the
new UseR series from Springer.

#### Soil science

Several packages are now available that implement R functions for
widely-used methods and approaches in pedology.

-   [soiltexture](../packages/soiltexture/index.html) provides functions
    for soil texture plot, classification and transformation.
-   [aqp](../packages/aqp/index.html) contains a collection of
    algorithms related to modeling of soil resources, soil
    classification, soil profile aggregation, and visualization.
-   Package [HydroMe](../packages/HydroMe/index.html) estimates the
    parameters in infiltration and water retention models by
    curve-fitting method.
-   The Soil Water project on r-forge.r-project.net provides packages
    providing soil water retention functions, soil hydraulic
    conductivity functions and pedotransfer functions to estimate their
    parameter from easily available soil properties. Two packages form
    the project:
    1.  [[soilwaterfun]{.Rforge}](https://R-Forge.R-project.org/projects/soilwaterfun/)
    2.  [[soilwaterptf]{.Rforge}](https://R-Forge.R-project.org/projects/soilwaterptf/)

#### Hydrology and Oceanography

A growing number of packages are available that implement methods
specifically related to the fields of hydrology and oceanography. Also
see the [Extreme Value](#extremes) and the [Climatology](#climatology)
sections for related packages.

-   Package [HydroMe](../packages/HydroMe/index.html) estimates the
    parameters in infiltration and water retention models by
    curve-fitting method.
-   [hydroTSM](../packages/hydroTSM/index.html) is a package for
    management, analysis, interpolation and plotting of time series used
    in hydrology and related environmental sciences.
-   [hydroGOF](../packages/hydroGOF/index.html) is a package
    implementing both statistical and graphical goodness-of-fit measures
    between observed and simulated values, mainly oriented to be used
    during the calibration, validation, and application of
    hydrological/environmental models. Related packages are
    [tiger](../packages/tiger/index.html), which allows temporally
    resolved groups of typical differences (errors) between two time
    series to be determined and visualized, and
    [qualV](../packages/qualV/index.html) which provides quantitative
    and qualitative criteria to compare models with data and to measure
    similarity of patterns
-   [hydroPSO](../packages/hydroPSO/index.html) is a model-independent
    global optimization tool for calibration of environmental and other
    real-world models that need to be executed from the system console.
    [hydroPSO](../packages/hydroPSO/index.html) implements a
    state-of-the-art PSO (SPSO-2011 and SPSO-2007 capable), with several
    fine-tuning options. The package is parallel-capable, to alleviate
    the computational burden of complex models.
-   [EcoHydRology](../packages/EcoHydRology/index.html) provides a
    flexible foundation for scientists, engineers, and policy makers to
    base teaching exercises as well as for more applied use to model
    complex eco-hydrological interactions.
-   [topmodel](../packages/topmodel/index.html) is a set of hydrological
    functions including an R implementation of the hydrological model
    TOPMODEL, which is based on the 1995 FORTRAN version by Keith Beven.
    New functionality is being developed as part of the
    [[RHydro]{.Rforge}](https://R-Forge.R-project.org/projects/rhydro/)
    package on R-Forge.
-   [dynatopmodel](../packages/dynatopmodel/index.html) is a native R
    implementation and enhancement of the Dynamic TOPMODEL, Beven and
    Freers\' (2001) extension to the semi-distributed hydrological model
    TOPMODEL (Beven and Kirkby, 1979).
-   [wasim](../packages/wasim/index.html) provides tools for data
    processing and visualisation of results of the hydrological model
    WASIM-ETH
-   Package [seacarb](../packages/seacarb/index.html) provides functions
    for calculating parameters of the seawater carbonate system.
-   Stephen Sefick\'s
    [StreamMetabolism](../packages/StreamMetabolism/index.html) package
    contains function for calculating stream metabolism characteristics,
    such as GPP, NDM, and R, from single station diurnal Oxygen curves.
-   Package [oce](../packages/oce/index.html) supports the analysis of
    Oceanographic data, including ADP measurements, CTD measurements,
    sectional data, sea-level time series, and coastline files.
-   The [nsRFA](../packages/nsRFA/index.html) package provides
    collection of statistical tools for objective (non-supervised)
    applications of the Regional Frequency Analysis methods in
    hydrology.
-   The [boussinesq](../packages/boussinesq/index.html) package is a
    collection of functions implementing the one-dimensional Boussinesq
    Equation (ground-water).
-   [rtop](../packages/rtop/index.html) is a package for geostatistical
    interpolation of data with irregular spatial support such as runoff
    related data or data from administrative units.

#### [Climatology]{#climatology}

Several packages related to the field of climatology.

-   [seas](../packages/seas/index.html) implements a number of functions
    for analysis and graphics of seasonal data.
-   [RMAWGEN](../packages/RMAWGEN/index.html) is set of S3 and S4
    functions for spatial multi-site stochastic generation of daily time
    series of temperature and precipitation making use of Vector
    Autoregressive Models.
-   [Interpol.T](../packages/Interpol.T/index.html) makes hourly
    interpolation of daily minimum and maximum temperature series for
    example when hourly time series must be downscaled from the daily
    information.

#### Palaeoecology and stratigraphic data

Several packages now provide specialist functionality for the import,
analysis, and plotting of palaeoecological data.

-   Transfer function models including weighted averaging (WA), modern
    analogue technique (MAT), Locally-weighted WA, & maximum likelihood
    (aka Gaussian logistic) regression (GLR) are provided by some or all
    of the [rioja](../packages/rioja/index.html), and
    [analogue](../packages/analogue/index.html) packages.
-   Import of common, legacy, palaeodata formats are provided by package
    [vegan](../packages/vegan/index.html) (cornell format) and
    [rioja](../packages/rioja/index.html) (cornell and Tilia format). In
    addition, [rioja](../packages/rioja/index.html) also allows for
    import of C2 model files.
-   Stratigraphic data plots can be drawn using `Stratiplot()` function
    in [analogue](../packages/analogue/index.html) and functions
    `strat.plot()` and `strat.plot.simple` in the
    [rioja](../packages/rioja/index.html) package.
-   [analogue](../packages/analogue/index.html) provides extensive
    support for developing and interpreting MAT transfer function
    models, including ROC curve analysis. Summary of stratigraphic data
    is supported via principal curves in the `prcurve()` function.
-   Constrained clustering of stratigraphic data is provided by function
    `chclust()` in the form of constrained hierarchical clustering in
    [rioja](../packages/rioja/index.html).

#### Other packages

Several other relevant contributed packages for R are available that do
not fit under nice headings.

-   [diveMove](../packages/diveMove/index.html) provides tools to
    represent, visualize, filter, analyse, and summarize time-depth
    recorder (TDR) data for research on animal diving and movement
    behaviour.
-   [latticeDensity](../packages/latticeDensity/index.html) implements
    methods for density estimation and nonparametric regression on
    irregular regions. A useful alternative to kernel density estimation
    for e.g. estimating animal densities and home ranges in regions with
    irregular boundaries or holes.
-   Andrew Robinson\'s [equivalence](../packages/equivalence/index.html)
    package provides some statistical tests and graphics for assessing
    tests of equivalence. Such tests have similarity as the alternative
    hypothesis instead of the null. The package contains functions to
    perform two one-sided t-tests (TOST) and paired t-tests of
    equivalence.
-   Thomas Petzoldt\'s [simecol](../packages/simecol/index.html) package
    provides an object oriented framework and tools to simulate
    ecological (and other) dynamic systems within R. See the [simecol
    website](http://www.simecol.de/) and a [R
    News](http://CRAN.R-project.org/doc/Rnews)
    [article](http://CRAN.R-project.org/doc/Rnews/Rnews_2003-3.pdf) on
    the package for further information.
-   Functions for circular statistics are found in
    [CircStats](../packages/CircStats/index.html) and
    [circular](../packages/circular/index.html).
-   Package [eco](../packages/eco/index.html) fits Bayesian models of
    ecological inference to 2 x 2 contingency tables.
-   Package [e1071](../packages/e1071/index.html) provides functions for
    latent class analysis, short time Fourier transform, fuzzy
    clustering, support vector machines, shortest path computation,
    bagged clustering, naive Bayes classifier, and more\...
-   Package [pgirmess](../packages/pgirmess/index.html) provides a suite
    of miscellaneous functions for data analysis in ecology.
-   [mefa](../packages/mefa/index.html) provides functions for handling
    and reporting on multivariate count data in ecology and
    biogeography.
-   Sensitivity analysis of models is provided by packages
    [sensitivity](../packages/sensitivity/index.html) and
    [fast](../packages/fast/index.html).
    [sensitivity](../packages/sensitivity/index.html) contains a
    collection of functions for factor screening and global sensitivity
    analysis of model output. [fast](../packages/fast/index.html) is an
    implementation of the Fourier Amplitude Sensitivity Test (FAST), a
    method to determine global sensitivities of a model on parameter
    changes with relatively few model runs.
-   Functions to analyze coherence, boundary clumping, and turnover
    following the pattern-based metacommunity analysis of [Leibold and
    Mikkelson (2002)](https://dx.doi.org/10.1034/j.1600-0706.2002.970210.x)
    are provided in the [metacom](../packages/metacom/index.html)
    package.
-   Growth curve estimation via noncrossing and nonparametric regression
    quantiles is implemented in package
    [quantregGrowth](../packages/quantregGrowth/index.html). A
    supporting paper is [Muggeo et
    al. (2013)](https://dx.doi.org/10.1007/s10651-012-0232-1) .
-   The [siplab](../packages/siplab/index.html) package provides an R
    platform for experimenting with spatially explicit individual-based
    vegetation models. A supporting paper is [García, O.
    (2014)](http://www.mcfns.com/index.php/Journal/article/view/6_36) .

</div>

### CRAN packages:

-   [ade4](../packages/ade4/index.html) (core)
-   [amap](../packages/amap/index.html)
-   [analogue](../packages/analogue/index.html)
-   [aod](../packages/aod/index.html)
-   [ape](../packages/ape/index.html)
-   [aqp](../packages/aqp/index.html)
-   [BiodiversityR](../packages/BiodiversityR/index.html)
-   [boussinesq](../packages/boussinesq/index.html)
-   [bReeze](../packages/bReeze/index.html)
-   [CircStats](../packages/CircStats/index.html)
-   [circular](../packages/circular/index.html)
-   [cluster](../packages/cluster/index.html) (core)
-   [cocorresp](../packages/cocorresp/index.html)
-   [Distance](../packages/Distance/index.html)
-   [diveMove](../packages/diveMove/index.html)
-   [dse](../packages/dse/index.html)
-   [dsm](../packages/dsm/index.html)
-   [DSpat](../packages/DSpat/index.html)
-   [dyn](../packages/dyn/index.html)
-   [dynatopmodel](../packages/dynatopmodel/index.html)
-   [dynlm](../packages/dynlm/index.html)
-   [e1071](../packages/e1071/index.html)
-   [earth](../packages/earth/index.html)
-   [eco](../packages/eco/index.html)
-   [ecodist](../packages/ecodist/index.html)
-   [EcoHydRology](../packages/EcoHydRology/index.html)
-   [EnvStats](../packages/EnvStats/index.html)
-   [equivalence](../packages/equivalence/index.html)
-   [evd](../packages/evd/index.html)
-   [evdbayes](../packages/evdbayes/index.html)
-   [evir](../packages/evir/index.html)
-   [extRemes](../packages/extRemes/index.html)
-   [fast](../packages/fast/index.html)
-   [FD](../packages/FD/index.html)
-   [flexmix](../packages/flexmix/index.html)
-   [forecast](../packages/forecast/index.html)
-   [fso](../packages/fso/index.html)
-   [gam](../packages/gam/index.html)
-   [gamair](../packages/gamair/index.html)
-   [hydroGOF](../packages/hydroGOF/index.html)
-   [HydroMe](../packages/HydroMe/index.html)
-   [hydroPSO](../packages/hydroPSO/index.html)
-   [hydroTSM](../packages/hydroTSM/index.html)
-   [Interpol.T](../packages/Interpol.T/index.html)
-   [ipred](../packages/ipred/index.html)
-   [ismev](../packages/ismev/index.html)
-   [labdsv](../packages/labdsv/index.html) (core)
-   [latticeDensity](../packages/latticeDensity/index.html)
-   [lme4](../packages/lme4/index.html)
-   [maptree](../packages/maptree/index.html)
-   [marked](../packages/marked/index.html)
-   [MASS](../packages/MASS/index.html) (core)
-   [mclust](../packages/mclust/index.html)
-   [mda](../packages/mda/index.html)
-   [mefa](../packages/mefa/index.html)
-   [metacom](../packages/metacom/index.html)
-   [mgcv](../packages/mgcv/index.html) (core)
-   [mrds](../packages/mrds/index.html)
-   [nlme](../packages/nlme/index.html)
-   [nsRFA](../packages/nsRFA/index.html)
-   [oce](../packages/oce/index.html)
-   [openair](../packages/openair/index.html)
-   [ouch](../packages/ouch/index.html)
-   [party](../packages/party/index.html)
-   [pastecs](../packages/pastecs/index.html)
-   [pgirmess](../packages/pgirmess/index.html)
-   [popbio](../packages/popbio/index.html)
-   [prabclus](../packages/prabclus/index.html)
-   [primer](../packages/primer/index.html)
-   [pscl](../packages/pscl/index.html)
-   [pvclust](../packages/pvclust/index.html)
-   [qualV](../packages/qualV/index.html)
-   [quantreg](../packages/quantreg/index.html)
-   [quantregGrowth](../packages/quantregGrowth/index.html)
-   [randomForest](../packages/randomForest/index.html)
-   [Rcapture](../packages/Rcapture/index.html)
-   [rioja](../packages/rioja/index.html)
-   [RMark](../packages/RMark/index.html)
-   [RMAWGEN](../packages/RMAWGEN/index.html)
-   [rpart](../packages/rpart/index.html)
-   [rtop](../packages/rtop/index.html)
-   [seacarb](../packages/seacarb/index.html)
-   [seas](../packages/seas/index.html)
-   [secr](../packages/secr/index.html)
-   [segmented](../packages/segmented/index.html)
-   [sensitivity](../packages/sensitivity/index.html)
-   [simba](../packages/simba/index.html)
-   [simecol](../packages/simecol/index.html)
-   [siplab](../packages/siplab/index.html)
-   [soiltexture](../packages/soiltexture/index.html)
-   [SPACECAP](../packages/SPACECAP/index.html)
-   [SpatialExtremes](../packages/SpatialExtremes/index.html)
-   [StreamMetabolism](../packages/StreamMetabolism/index.html)
-   [strucchange](../packages/strucchange/index.html)
-   [surveillance](../packages/surveillance/index.html)
-   [tiger](../packages/tiger/index.html)
-   [topmodel](../packages/topmodel/index.html)
-   [tseries](../packages/tseries/index.html)
-   [unmarked](../packages/unmarked/index.html)
-   [untb](../packages/untb/index.html)
-   [vegan](../packages/vegan/index.html) (core)
-   [vegetarian](../packages/vegetarian/index.html)
-   [VGAM](../packages/VGAM/index.html)
-   [wasim](../packages/wasim/index.html)
-   [zoo](../packages/zoo/index.html)

### Related links:

-   CRAN Task View: [Spatial](Spatial.html)
-   CRAN Task View: [Multivariate](Multivariate.html)
-   CRAN Task View: [Cluster](Cluster.html)
-   CRAN Task View: [Phylogenetics](Phylogenetics.html)
-   R-Forge Project:
    [[soilwaterfun]{.Rforge}](https://R-Forge.R-project.org/projects/soilwaterfun/)
-   R-Forge Project:
    [[soilwaterptf]{.Rforge}](https://R-Forge.R-project.org/projects/soilwaterptf/)
-   R-Forge Project:
    [[RHydro]{.Rforge}](https://R-Forge.R-project.org/projects/rhydro/)
-   [The vegan development site on
    R-Forge](http://vegan.r-forge.r-project.org/)
-   [Thomas Yee\'s VGAM package for
    R](http://www.stat.auckland.ac.nz/~yee/VGAM/)
-   [The cocorresp development site on
    R-Forge](http://cocorresp.r-forge.r-project.org/)
-   [The analogue development site on
    R-Forge](http://analogue.r-forge.r-project.org/)
-   [Brodgar](http://www.brodgar.com/)
-   [More information on the ade4 package can be found on the ADE4
    website](http://pbil.univ-lyon1.fr/ADE-4/)
-   [Mike Palmer\'s Ordination web site](http://ordination.okstate.edu/)
-   [Thomas Petzoldt\'s page about ecological modelling with
    R](https://wwwpub.zih.tu-dresden.de/%7Epetzoldt/)
-   [The FLR project web page for Fisheries Science
    in R.](http://www.flr-project.org/)
