## CRAN Task View: Multivariate Statistics

  ----------------- ------------------------------------------------
  **Maintainer:**   Paul Hewson
  **Contact:**      Paul.Hewson at plymouth.ac.uk
  **Version:**      2018-02-12
  **URL:**          <https://CRAN.R-project.org/view=Multivariate>
  ----------------- ------------------------------------------------

<div>

Base R contains most of the functionality for classical multivariate
analysis, somewhere. There are a large number of packages on CRAN which
extend this methodology, a brief overview is given below.
Application-specific uses of multivariate statistics are described in
relevant task views, for example whilst principal components are listed
here, ordination is covered in the [Environmetrics](Environmetrics.html)
task view. Further information on supervised classification can be found
in the [MachineLearning](MachineLearning.html) task view, and
unsupervised classification in the [Cluster](Cluster.html) task view.

The packages in this view can be roughly structured into the following
topics. If you think that some package is missing from the list, please
let me know.

**Visualising multivariate data**

-   *Graphical Procedures:* A range of base graphics (e.g. `pairs()` and
    `coplot()`) and [lattice](../packages/lattice/index.html) functions
    (e.g. `xyplot()` and `splom()`) are useful for visualising pairwise
    arrays of 2-dimensional scatterplots, clouds and 3-dimensional
    densities. `scatterplot.matrix` in the
    [car](../packages/car/index.html) provides usefully enhanced
    pairwise scatterplots. Beyond this,
    [scatterplot3d](../packages/scatterplot3d/index.html) provides 3
    dimensional scatterplots, [aplpack](../packages/aplpack/index.html)
    provides bagplots and `spin3R()`, a function for rotating 3d clouds.
    [misc3d](../packages/misc3d/index.html), dependent upon
    [rgl](../packages/rgl/index.html), provides animated functions
    within R useful for visualising densities.
    [YaleToolkit](../packages/YaleToolkit/index.html) provides a range
    of useful visualisation techniques for multivariate data. More
    specialised multivariate plots include the following: `faces()` in
    [aplpack](../packages/aplpack/index.html) provides Chernoff\'s
    faces; `parcoord()` from [MASS](../packages/MASS/index.html)
    provides parallel coordinate plots; `stars()` in graphics provides a
    choice of star, radar and cobweb plots respectively. `mstree()` in
    [ade4](../packages/ade4/index.html) and `spantree()` in
    [vegan](../packages/vegan/index.html) provide minimum spanning tree
    functionality. [calibrate](../packages/calibrate/index.html)
    supports biplot and scatterplot axis labelling.
    [geometry](../packages/geometry/index.html), which provides an
    interface to the qhull library, gives indices to the relevant points
    via `convexhulln()`. [ellipse](../packages/ellipse/index.html) draws
    ellipses for two parameters, and provides `plotcorr()`, visual
    display of a correlation matrix.
    [denpro](../packages/denpro/index.html) provides level set trees for
    multivariate visualisation. Mosaic plots are available via
    `mosaicplot()` in graphics and `mosaic()` in
    [vcd](../packages/vcd/index.html) that also contains other
    visualization techniques for multivariate categorical data.
    [gclus](../packages/gclus/index.html) provides a number of cluster
    specific graphical enhancements for scatterplots and parallel
    coordinate plots See the links for a reference to GGobi.
    [rggobi](../packages/rggobi/index.html) interfaces with GGobi.
    [xgobi](../packages/xgobi/index.html) interfaces to the XGobi and
    XGvis programs which allow linked, dynamic multivariate plots as
    well as projection pursuit. Finally,
    [iplots](../packages/iplots/index.html) allows particularly powerful
    dynamic interactive graphics, of which interactive parallel
    co-ordinate plots and mosaic plots may be of great interest.
    Seriation methods are provided by
    [seriation](../packages/seriation/index.html) which can reorder
    matrices and dendrograms.
-   *Data Preprocessing:* `summarize()` and `summary.formula()` in
    [Hmisc](../packages/Hmisc/index.html) assist with descriptive
    functions; from the same package `varclus()` offers variable
    clustering while `dataRep()` and `find.matches()` assist in
    exploring a given dataset in terms of representativeness and finding
    matches. Whilst `dist()` in base and `daisy()` in
    [cluster](../packages/cluster/index.html) provide a wide range of
    distance measures, [proxy](../packages/proxy/index.html) provides a
    framework for more distance measures, including measures between
    matrices. [simba](../packages/simba/index.html) provides functions
    for dealing with presence / absence data including similarity
    matrices and reshaping.

**Hypothesis testing**

-   [ICSNP](../packages/ICSNP/index.html) provides Hotellings T2 test as
    well as a range of non-parametric tests including location tests
    based on marginal ranks, spatial median and spatial signs
    computation, estimates of shape. Non-parametric two sample tests are
    also available from [cramer](../packages/cramer/index.html) and
    spatial sign and rank tests to investigate location, sphericity and
    independence are available in
    [SpatialNP](../packages/SpatialNP/index.html).

**Multivariate distributions**

-   *Descriptive measures:* `cov()` and `cor()` in stats will provide
    estimates of the covariance and correlation matrices respectively.
    [ICSNP](../packages/ICSNP/index.html) offers several descriptive
    measures such as `spatial.median()` which provides an estimate of
    the spatial median and further functions which provide estimates of
    scatter. Further robust methods are provided such as `cov.rob()` in
    MASS which provides robust estimates of the variance-covariance
    matrix by minimum volume ellipsoid, minimum covariance determinant
    or classical product-moment.
    [covRobust](../packages/covRobust/index.html) provides robust
    covariance estimation via nearest neighbor variance estimation.
    [robustbase](../packages/robustbase/index.html) provides robust
    covariance estimation via fast minimum covariance determinant with
    `covMCD()` and the Orthogonalized pairwise estimate of
    Gnanadesikan-Kettenring via `covOGK()`. Scalable robust methods are
    provided within [rrcov](../packages/rrcov/index.html) also using
    fast minimum covariance determinant with `covMcd()` as well as
    M-estimators with `covMest()`.
    [corpcor](../packages/corpcor/index.html) provides shrinkage
    estimation of large scale covariance and (partial) correlation
    matrices.
-   *Densities (estimation and simulation):* `mvnorm()` in MASS
    simulates from the multivariate normal distribution.
    [mvtnorm](../packages/mvtnorm/index.html) also provides simulation
    as well as probability and quantile functions for both the
    multivariate t distribution and multivariate normal distributions as
    well as density functions for the multivariate normal distribution.
    [mnormt](../packages/mnormt/index.html) provides multivariate normal
    and multivariate t density and distribution functions as well as
    random number simulation. [sn](../packages/sn/index.html) provides
    density, distribution and random number generation for the
    multivariate skew normal and skew t distribution.
    [delt](../packages/delt/index.html) provides a range of functions
    for estimating multivariate densities by CART and greedy methods.
    Comprehensive information on mixtures is given in the
    [Cluster](Cluster.html) view, some density estimates and random
    numbers are provided by `rmvnorm.mixt()` and `dmvnorm.mixt()` in
    [ks](../packages/ks/index.html), mixture fitting is also provided
    within [bayesm](../packages/bayesm/index.html). Functions to
    simulate from the Wishart distribution are provided in a number of
    places, such as `rwishart()` in
    [bayesm](../packages/bayesm/index.html) and `rwish()` in
    [MCMCpack](../packages/MCMCpack/index.html) (the latter also has a
    density function `dwish()`). `bkde2D()` from
    [KernSmooth](../packages/KernSmooth/index.html) and `kde2d()` from
    MASS provide binned and non-binned 2-dimensional kernel density
    estimation, [ks](../packages/ks/index.html) also provides
    multivariate kernel smoothing as does
    [ash](../packages/ash/index.html) and
    [GenKern](../packages/GenKern/index.html).
    [prim](../packages/prim/index.html) provides patient rule induction
    methods to attempt to find regions of high density in high
    dimensional multivariate data,
    [feature](../packages/feature/index.html) also provides methods for
    determining feature significance in multivariate data (such as in
    relation to local modes).
-   *Assessing normality:*
    [mvnormtest](../packages/mvnormtest/index.html) provides a
    multivariate extension to the Shapiro-Wilks test,
    [mvoutlier](../packages/mvoutlier/index.html) provides multivariate
    outlier detection based on robust methods.
    [ICS](../packages/ICS/index.html) provides tests for
    multi-normality. `mvnorm.etest()` in
    [energy](../packages/energy/index.html) provides an assessment of
    normality based on E statistics (energy); in the same package
    `k.sample()` assesses a number of samples for equal distributions.
    Tests for Wishart-distributed covariance matrices are given by
    `mauchly.test()` in stats.
-   *Copulas:* [copula](../packages/copula/index.html) provides routines
    for a range of (elliptical and archimedean) copulas including
    normal, t, Clayton, Frank, Gumbel,
    [fgac](../packages/fgac/index.html) provides generalised archimedian
    copula.

**Linear models**

-   From stats, `lm()` (with a matrix specified as the dependent
    variable) offers multivariate linear models, `anova.mlm()` provides
    comparison of multivariate linear models. `manova()` offers MANOVA.
    [sn](../packages/sn/index.html) provides `msn.mle()` and `mst.mle()`
    which fit multivariate skew normal and multivariate skew t models.
    [pls](../packages/pls/index.html) provides partial least squares
    regression (PLSR) and principal component regression,
    [ppls](../packages/ppls/index.html) provides penalized partial least
    squares, [dr](../packages/dr/index.html) provides dimension
    reduction regression options such as `"sir"` (sliced inverse
    regression), `"save"` (sliced average variance estimation).
    [plsgenomics](../packages/plsgenomics/index.html) provides partial
    least squares analyses for genomics.
    [relaimpo](../packages/relaimpo/index.html) provides functions to
    investigate the relative importance of regression parameters.

**Projection methods**

-   *Principal components:* these can be fitted with `prcomp()` (based
    on `svd()`, preferred) as well as `princomp()` (based on `eigen()`
    for compatibility with S-PLUS) from stats. `pc1()` in
    [Hmisc](../packages/Hmisc/index.html) provides the first principal
    component and gives coefficients for unscaled data. Additional
    support for an assessment of the scree plot can be found in
    [nFactors](../packages/nFactors/index.html), whereas
    [paran](../packages/paran/index.html) provides routines for Horn\'s
    evaluation of the number of dimensions to retain. For wide matrices,
    [gmodels](../packages/gmodels/index.html) provides `fast.prcomp()`
    and `fast.svd()`. [kernlab](../packages/kernlab/index.html) uses
    kernel methods to provide a form of non-linear principal components
    with `kpca()`. [pcaPP](../packages/pcaPP/index.html) provides robust
    principal components by means of projection pursuit.
    [amap](../packages/amap/index.html) provides further robust and
    parallelised methods such as a form of generalised and robust
    principal component analysis via `acpgen()` and `acprob()`
    respectively. Further options for principal components in an
    ecological setting are available within
    [ade4](../packages/ade4/index.html) and in a sensory setting in
    [SensoMineR](../packages/SensoMineR/index.html).
    [psy](../packages/psy/index.html) provides a variety of routines
    useful in psychometry, in this context these include `sphpca()`
    which maps onto a sphere and `fpca()` where some variables may be
    considered as dependent as well as `scree.plot()` which has the
    option of adding simulation results to help assess the observed
    data. [PTAk](../packages/PTAk/index.html) provides principal tensor
    analysis analagous to both PCA and correspondence analysis.
    [smatr](../packages/smatr/index.html) provides standardised major
    axis estimation with specific application to allometry.
-   *Canonical Correlation:* `cancor()` in stats provides canonical
    correlation. [kernlab](../packages/kernlab/index.html) uses kernel
    methods to provide robust canonical correlation with `kcca()`.
    [concor](../packages/concor/index.html) provides a number of
    concordance methods.
-   *Redundancy Analysis:* [calibrate](../packages/calibrate/index.html)
    provides `rda()` for redundancy analysis as well as further options
    for canonical correlation. [fso](../packages/fso/index.html)
    provides fuzzy set ordination, which extends ordination beyond
    methods available from linear algebra.
-   *Independent Components:* [fastICA](../packages/fastICA/index.html)
    provides fastICA algorithms to perform independent component
    analysis (ICA) and Projection Pursuit, and
    [PearsonICA](../packages/PearsonICA/index.html) uses score
    functions. [ICS](../packages/ICS/index.html) provides either an
    invariant co-ordinate system or independent components.
    [JADE](../packages/JADE/index.html) adds an interface to the JADE
    algorithm, as well as providing some diagnostics for ICA.
-   *Procrustes analysis:* `procrustes()` in
    [vegan](../packages/vegan/index.html) provides procrustes analysis,
    this package also provides functions for ordination and further
    information on that area is given in the
    [Environmetrics](Environmetrics.html) task view. Generalised
    procrustes analysis via `GPA()` is available from
    [FactoMineR](../packages/FactoMineR/index.html).

**Principal coordinates / scaling methods**

-   `cmdscale()` in stats provides classical multidimensional scaling
    (principal coordinates analysis), `sammon()` and `isoMDS()` in MASS
    offer Sammon and Kruskal\'s non-metric multidimensional scaling.
    [vegan](../packages/vegan/index.html) provides wrappers and
    post-processing for non-metric MDS. `indscal()` is provided by
    [SensoMineR](../packages/SensoMineR/index.html).

**Unsupervised classification**

-   *Cluster analysis:* A comprehensive overview of clustering methods
    available within R is provided by the [Cluster](Cluster.html) task
    view. Standard techniques include hierarchical clustering by
    `hclust()` and k-means clustering by `kmeans()` in stats. A range of
    established clustering and visualisation techniques are also
    available in [cluster](../packages/cluster/index.html), some cluster
    validation routines are available in
    [clv](../packages/clv/index.html) and the Rand index can be computed
    from `classAgreement()` in [e1071](../packages/e1071/index.html).
    Trimmed cluster analysis is available from
    [trimcluster](../packages/trimcluster/index.html), cluster ensembles
    are available from [clue](../packages/clue/index.html), methods to
    assist with choice of routines are available in
    [clusterSim](../packages/clusterSim/index.html) and hybrid
    methodology is provided by
    [hybridHclust](../packages/hybridHclust/index.html). Distance
    measures ( `edist()`) and hierarchical clustering (
    `hclust.energy()`) based on E-statistics are available in
    [energy](../packages/energy/index.html). Mahalanobis distance based
    clustering (for fixed points as well as clusterwise regression) are
    available from [fpc](../packages/fpc/index.html).
    [clustvarsel](../packages/clustvarsel/index.html) provides variable
    selection within model-based clustering. Fuzzy clustering is
    available within [cluster](../packages/cluster/index.html) as well
    as via the
    [[hopach]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/hopach.html)
    (Hierarchical Ordered Partitioning and Collapsing Hybrid) algorithm.
    [kohonen](../packages/kohonen/index.html) provides supervised and
    unsupervised SOMs for high dimensional spectra or patterns.
    [clusterGeneration](../packages/clusterGeneration/index.html) helps
    simulate clusters. The [Environmetrics](Environmetrics.html) task
    view also gives a topic-related overview of some clustering
    techniques. Model based clustering is available in
    [mclust](../packages/mclust/index.html).
-   *Tree methods:* Full details on tree methods are given in the
    [MachineLearning](MachineLearning.html) task view. Suffice to say
    here that classification trees are sometimes considered within
    multivariate methods; [rpart](../packages/rpart/index.html) is most
    used for this purpose. [party](../packages/party/index.html)
    provides recursive partitioning. Classification and regression
    training is provided by [caret](../packages/caret/index.html).
    [kknn](../packages/kknn/index.html) provides k-nearest neighbour
    methods which can be used for regression as well as classification.

**Supervised classification and discriminant analysis**

-   `lda()` and `qda()` within MASS provide linear and quadratic
    discrimination respectively. [mda](../packages/mda/index.html)
    provides mixture and flexible discriminant analysis with `mda()` and
    `fda()` as well as multivariate adaptive regression splines with
    `mars()` and adaptive spline backfitting with the `bruto()`
    function. Multivariate adaptive regression splines can also be found
    in [earth](../packages/earth/index.html).
    [rda](../packages/rda/index.html) provides classification for high
    dimensional data by means of shrunken centroids regularized
    discriminant analysis. Package [class](../packages/class/index.html)
    provides k-nearest neighbours by `knn()`,
    [knncat](../packages/knncat/index.html) provides k-nearest
    neighbours for categorical variables.
    [SensoMineR](../packages/SensoMineR/index.html) provides `FDA()` for
    factorial discriminant analysis. A number of packages provide for
    dimension reduction with the classification.
    [klaR](../packages/klaR/index.html) includes variable selection and
    robustness against multicollinearity as well as a number of
    visualisation routines. [superpc](../packages/superpc/index.html)
    provides principal components for supervised classification, whereas
    [[gpls]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/gpls.html)
    provides classification using generalised partial least squares.
    [hddplot](../packages/hddplot/index.html) provides cross-validated
    linear discriminant calculations to determine the optimum number of
    features. [ROCR](../packages/ROCR/index.html) provides a range of
    methods for assessing classifier performance. Further information on
    supervised classification can be found in the
    [MachineLearning](MachineLearning.html) task view.

**Correspondence analysis**

-   `corresp()` and `mca()` in MASS provide simple and multiple
    correspondence analysis respectively.
    [ca](../packages/ca/index.html) also provides single, multiple and
    joint correspondence analysis. `ca()` and `mca()` in
    [ade4](../packages/ade4/index.html) provide correspondence and
    multiple correspondence analysis respectively, as well as adding
    homogeneous table analysis with `hta()`. Further functionality is
    also available within [vegan](../packages/vegan/index.html)
    co-correspondence is available from
    [cocorresp](../packages/cocorresp/index.html).
    [FactoMineR](../packages/FactoMineR/index.html) provides `CA()` and
    `MCA()` which also enable simple and multiple correspondence
    analysis as well as associated graphical routines.
    [homals](../packages/homals/index.html) provides homogeneity
    analysis.

**Missing data**

-   [mitools](../packages/mitools/index.html) provides tools for
    multiple imputation, [mice](../packages/mice/index.html) provides
    multivariate imputation by chained equations
    [mvnmle](../packages/mvnmle/index.html) provides ML estimation for
    multivariate normal data with missing values,
    [mix](../packages/mix/index.html) provides multiple imputation for
    mixed categorical and continuous data.
    [pan](../packages/pan/index.html) provides multiple imputation for
    missing panel data. [VIM](../packages/VIM/index.html) provides
    methods for the visualisation as well as imputation of missing data.
    `aregImpute()` and `transcan()` from
    [Hmisc](../packages/Hmisc/index.html) provide further imputation
    methods. [monomvn](../packages/monomvn/index.html) deals with
    estimation models where the missing data pattern is monotone.

**Latent variable approaches**

-   `factanal()` in stats provides factor analysis by maximum
    likelihood, Bayesian factor analysis is provided for Gaussian,
    ordinal and mixed variables in
    [MCMCpack](../packages/MCMCpack/index.html).
    [GPArotation](../packages/GPArotation/index.html) offers GPA
    (gradient projection algorithm) factor rotation.
    [sem](../packages/sem/index.html) fits linear structural equation
    models and [ltm](../packages/ltm/index.html) provides latent trait
    models under item response theory and range of extensions to Rasch
    models can be found in [eRm](../packages/eRm/index.html).
    [FactoMineR](../packages/FactoMineR/index.html) provides a wide
    range of Factor Analysis methods, including `MFA()` and `HMFA()` for
    multiple and hierarchical multiple factor analysis as well as
    `ADFM()` for multiple factor analysis of quantitative and
    qualitative data. [tsfa](../packages/tsfa/index.html) provides
    factor analysis for time series.
    [poLCA](../packages/poLCA/index.html) provides latent class and
    latent class regression models for a variety of outcome variables.

**Modelling non-Gaussian data**

-   [MNP](../packages/MNP/index.html) provides Bayesian multinomial
    probit models, [polycor](../packages/polycor/index.html) provides
    polychoric and tetrachoric correlation matrices.
    [bayesm](../packages/bayesm/index.html) provides a range of models
    such as seemingly unrelated regression, multinomial logit/probit,
    multivariate probit and instrumental variables.
    [VGAM](../packages/VGAM/index.html) provides Vector Generalised
    Linear and Additive Models, Reduced Rank regression

**Matrix manipulations**

-   As a vector- and matrix-based language, base R ships with many
    powerful tools for doing matrix manipulations, which are
    complemented by the packages [Matrix](../packages/Matrix/index.html)
    and [SparseM](../packages/SparseM/index.html).
    [matrixcalc](../packages/matrixcalc/index.html) adds functions for
    matrix differential calculus. Some further sparse matrix
    functionality is also available from
    [spam](../packages/spam/index.html).

**Miscellaneous utilities**

-   [abind](../packages/abind/index.html) generalises `cbind()` and
    `rbind()` for arrays, `mApply()` in
    [Hmisc](../packages/Hmisc/index.html) generalises `apply()` for
    matrices and passes multiple functions. In addition to functions
    listed earlier, [sn](../packages/sn/index.html) provides operations
    such as marginalisation, affine transformations and graphics for the
    multivariate skew normal and skew t distribution.
    [mAr](../packages/mAr/index.html) provides for vector
    auto-regression. `rm.boot()` from
    [Hmisc](../packages/Hmisc/index.html) bootstraps repeated measures
    models. [psy](../packages/psy/index.html) also provides a range of
    statistics based on Cohen\'s kappa including weighted measures and
    agreement among more than 2 raters.
    [cwhmisc](../packages/cwhmisc/index.html) contains a number of
    interesting support functions which are of interest, such as
    `ellipse()`, `normalise()` and various rotation functions.
    [desirability](../packages/desirability/index.html) provides
    functions for multivariate optimisation.
    [geozoo](../packages/geozoo/index.html) provides plotting of
    geometric objects in GGobi.

</div>

### CRAN packages:

-   [abind](../packages/abind/index.html)
-   [ade4](../packages/ade4/index.html) (core)
-   [amap](../packages/amap/index.html)
-   [aplpack](../packages/aplpack/index.html)
-   [ash](../packages/ash/index.html)
-   [bayesm](../packages/bayesm/index.html)
-   [ca](../packages/ca/index.html)
-   [calibrate](../packages/calibrate/index.html)
-   [car](../packages/car/index.html)
-   [caret](../packages/caret/index.html)
-   [class](../packages/class/index.html)
-   [clue](../packages/clue/index.html)
-   [cluster](../packages/cluster/index.html) (core)
-   [clusterGeneration](../packages/clusterGeneration/index.html)
-   [clusterSim](../packages/clusterSim/index.html)
-   [clustvarsel](../packages/clustvarsel/index.html)
-   [clv](../packages/clv/index.html)
-   [cocorresp](../packages/cocorresp/index.html)
-   [concor](../packages/concor/index.html)
-   [copula](../packages/copula/index.html)
-   [corpcor](../packages/corpcor/index.html)
-   [covRobust](../packages/covRobust/index.html)
-   [cramer](../packages/cramer/index.html)
-   [cwhmisc](../packages/cwhmisc/index.html)
-   [delt](../packages/delt/index.html)
-   [denpro](../packages/denpro/index.html)
-   [desirability](../packages/desirability/index.html)
-   [dr](../packages/dr/index.html)
-   [e1071](../packages/e1071/index.html)
-   [earth](../packages/earth/index.html)
-   [ellipse](../packages/ellipse/index.html)
-   [energy](../packages/energy/index.html)
-   [eRm](../packages/eRm/index.html)
-   [FactoMineR](../packages/FactoMineR/index.html)
-   [fastICA](../packages/fastICA/index.html)
-   [feature](../packages/feature/index.html)
-   [fgac](../packages/fgac/index.html)
-   [fpc](../packages/fpc/index.html)
-   [fso](../packages/fso/index.html)
-   [gclus](../packages/gclus/index.html)
-   [GenKern](../packages/GenKern/index.html)
-   [geometry](../packages/geometry/index.html)
-   [geozoo](../packages/geozoo/index.html)
-   [gmodels](../packages/gmodels/index.html)
-   [GPArotation](../packages/GPArotation/index.html)
-   [hddplot](../packages/hddplot/index.html)
-   [Hmisc](../packages/Hmisc/index.html)
-   [homals](../packages/homals/index.html)
-   [hybridHclust](../packages/hybridHclust/index.html)
-   [ICS](../packages/ICS/index.html)
-   [ICSNP](../packages/ICSNP/index.html)
-   [iplots](../packages/iplots/index.html)
-   [JADE](../packages/JADE/index.html)
-   [kernlab](../packages/kernlab/index.html)
-   [KernSmooth](../packages/KernSmooth/index.html)
-   [kknn](../packages/kknn/index.html)
-   [klaR](../packages/klaR/index.html)
-   [knncat](../packages/knncat/index.html)
-   [kohonen](../packages/kohonen/index.html)
-   [ks](../packages/ks/index.html)
-   [lattice](../packages/lattice/index.html)
-   [ltm](../packages/ltm/index.html)
-   [mAr](../packages/mAr/index.html)
-   [MASS](../packages/MASS/index.html) (core)
-   [Matrix](../packages/Matrix/index.html)
-   [matrixcalc](../packages/matrixcalc/index.html)
-   [mclust](../packages/mclust/index.html)
-   [MCMCpack](../packages/MCMCpack/index.html)
-   [mda](../packages/mda/index.html)
-   [mice](../packages/mice/index.html)
-   [misc3d](../packages/misc3d/index.html)
-   [mitools](../packages/mitools/index.html)
-   [mix](../packages/mix/index.html)
-   [mnormt](../packages/mnormt/index.html)
-   [MNP](../packages/MNP/index.html)
-   [monomvn](../packages/monomvn/index.html)
-   [mvnmle](../packages/mvnmle/index.html)
-   [mvnormtest](../packages/mvnormtest/index.html)
-   [mvoutlier](../packages/mvoutlier/index.html)
-   [mvtnorm](../packages/mvtnorm/index.html)
-   [nFactors](../packages/nFactors/index.html)
-   [pan](../packages/pan/index.html)
-   [paran](../packages/paran/index.html)
-   [party](../packages/party/index.html)
-   [pcaPP](../packages/pcaPP/index.html)
-   [PearsonICA](../packages/PearsonICA/index.html)
-   [pls](../packages/pls/index.html)
-   [plsgenomics](../packages/plsgenomics/index.html)
-   [poLCA](../packages/poLCA/index.html)
-   [polycor](../packages/polycor/index.html)
-   [ppls](../packages/ppls/index.html)
-   [prim](../packages/prim/index.html)
-   [proxy](../packages/proxy/index.html)
-   [psy](../packages/psy/index.html)
-   [PTAk](../packages/PTAk/index.html)
-   [rda](../packages/rda/index.html)
-   [relaimpo](../packages/relaimpo/index.html)
-   [rggobi](../packages/rggobi/index.html)
-   [rgl](../packages/rgl/index.html)
-   [robustbase](../packages/robustbase/index.html)
-   [ROCR](../packages/ROCR/index.html)
-   [rpart](../packages/rpart/index.html)
-   [rrcov](../packages/rrcov/index.html)
-   [scatterplot3d](../packages/scatterplot3d/index.html)
-   [sem](../packages/sem/index.html)
-   [SensoMineR](../packages/SensoMineR/index.html)
-   [seriation](../packages/seriation/index.html)
-   [simba](../packages/simba/index.html)
-   [smatr](../packages/smatr/index.html)
-   [sn](../packages/sn/index.html)
-   [spam](../packages/spam/index.html)
-   [SparseM](../packages/SparseM/index.html)
-   [SpatialNP](../packages/SpatialNP/index.html)
-   [superpc](../packages/superpc/index.html)
-   [trimcluster](../packages/trimcluster/index.html)
-   [tsfa](../packages/tsfa/index.html)
-   [vcd](../packages/vcd/index.html)
-   [vegan](../packages/vegan/index.html) (core)
-   [VGAM](../packages/VGAM/index.html)
-   [VIM](../packages/VIM/index.html)
-   [xgobi](../packages/xgobi/index.html)
-   [YaleToolkit](../packages/YaleToolkit/index.html)

### Related links:

-   CRAN Task View: [Cluster](Cluster.html)
-   CRAN Task View: [Environmetrics](Environmetrics.html)
-   CRAN Task View: [MachineLearning](MachineLearning.html)
-   Bioconductor Package:
    [[gpls]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/gpls.html)
-   Bioconductor Package:
    [[hopach]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/hopach.html)
-   [GGobi (interactive dynamic visualisation software, available
    standalone or as an R library)](http://www.ggobi.org)
-   [Hmisc functions related to multivariate
    analysis](http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/HmiscMultivar)
-   [Psychometrics in R, Jan de
    Leeuw](http://www.cuddyvalley.org/psychoR/)
-   [qhull library](http://www.qhull.org)
