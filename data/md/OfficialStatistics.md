## CRAN Task View: Official Statistics & Survey Methodology

  ----------------- ------------------------------------------------------
  **Maintainer:**   Matthias Templ
  **Contact:**      matthias.templ at gmail.com
  **Version:**      2018-04-20
  **URL:**          <https://CRAN.R-project.org/view=OfficialStatistics>
  ----------------- ------------------------------------------------------

<div>

This CRAN task view contains a list of packages that includes methods
typically used in official statistics and survey methodology. Many
packages provide functionality for more than one of the topics listed
below. Therefore this list is not a strict categorization and packages
can be listed more than once. Certain data import/export facilities
regarding to often used statistical software tools like SPSS, SAS or
Stata are mentioned in the end of the task view.

**Complex Survey Design: Sampling and Sample Size Calculation**

-   Package [sampling](../packages/sampling/index.html) includes many
    different algorithms (Brewer, Midzuno, pps, systematic, Sampford,
    balanced (cluster or stratified) sampling via the cube method, etc.)
    for drawing survey samples and calibrating the design weights.
-   R package [surveyplanning](../packages/surveyplanning/index.html)
    includes tools for sample survey planning, including sample size
    calculation, estimation of expected precision for the estimates of
    totals, and calculation of optimal sample size allocation.
-   Package [simFrame](../packages/simFrame/index.html) includes a fast
    (compiled C-Code) version of Midzuno sampling.
-   The [pps](../packages/pps/index.html) package contains functions to
    select samples using pps sampling. Also stratified simple random
    sampling is possible as well as to compute joint inclusion
    probabilities for Sampford\'s method of pps sampling.
-   Package [stratification](../packages/stratification/index.html)
    allows univariate stratification of survey populations with a
    generalisation of the Lavallee-Hidiroglou method.
-   Package [SamplingStrata](../packages/SamplingStrata/index.html)
    offers an approach for choosing the best stratification of a
    sampling frame in a multivariate and multidomain setting, where the
    sampling sizes in each strata are determined in order to satisfy
    accuracy constraints on target estimates. To evaluate the
    distribution of target variables in different strata, information of
    the sampling frame, or data from previous rounds of the same survey,
    may be used.
-   The package
    [BalancedSampling](../packages/BalancedSampling/index.html) selects
    balanced and spatially balanced probability samples in
    multi-dimensional spaces with any prescribed inclusion
    probabilities. It also includes the local pivot method, the cube and
    local cube method and a few more methods.
-   Package [gridsample](../packages/gridsample/index.html) selects PSUs
    within user-defined strata using gridded population data, given
    desired numbers of sampled households within each PSU. The
    population densities used to create PSUs are drawn from rasters
-   Package [PracTools](../packages/PracTools/index.html) contains
    functions for sample size calculation for survey samples using
    stratified or clustered one-, two-, and three-stage sample designs
    as well as functions to compute variance components for multistage
    designs and sample sizes in two-phase designs.
-   Package
    [samplesize4surveys](../packages/samplesize4surveys/index.html)
    computes the required sample size for estimation of totals, means
    and proportions under complex sampling designs.

**Complex Survey Design: Point and Variance Estimation and Model
Fitting**

-   Package [survey](../packages/survey/index.html) works with survey
    samples. It allows to specify a complex survey design (stratified
    sampling design, cluster sampling, multi-stage sampling and pps
    sampling with or without replacement). Once the given survey design
    is specified within the function `svydesign()`, point and variance
    estimates can be computed. The resulting object can be used to
    estimate (Horvitz-Thompson-) totals, means, ratios and quantiles for
    domains or the whole survey sample, and to apply regression models.
    Variance estimation for means, totals and ratios can be done either
    by Taylor linearization or resampling (BRR, jackkife, bootstrap or
    user-defined).
-   The methods from the [survey](../packages/survey/index.html) package
    are called from package [srvyr](../packages/srvyr/index.html) using
    the dplyr syntax, i.e., piping, verbs like `group_by` and
    `summarize`, and other dplyr-inspired syntactic style when
    calculating summary statistics on survey data.
-   Package [convey](../packages/convey/index.html) extends package
    [survey](../packages/survey/index.html) \-- see the topic about
    indicators below.
-   Package [laeken](../packages/laeken/index.html) provides functions
    to estimate certain Laeken indicators (at-risk-of-poverty rate,
    quintile share ratio, relative median risk-of-poverty gap, Gini
    coefficient) including their variance for domains and strata using a
    calibrated bootstrap.
-   Package [simFrame](../packages/simFrame/index.html) allows to
    compare (user-defined) point and variance estimators in a simulation
    environment. It provides a framework for comparing different point
    and variance estimators under different survey designs as well as
    different conditions regarding missing values, representative and
    non-representative outliers.
-   The [lavaan.survey](../packages/lavaan.survey/index.html) package
    provides a wrapper function for packages
    [survey](../packages/survey/index.html) and
    [lavaan](../packages/lavaan/index.html). It can be used for fitting
    structural equation models (SEM) on samples from complex designs.
    Using the design object functionality from package
    [survey](../packages/survey/index.html), lavaan objects are re-fit
    (corrected) with the `lavaan.survey()` function of package
    [lavaan.survey](../packages/lavaan.survey/index.html). This allows
    for the incorporation of clustering, stratification, sampling
    weights, and finite population corrections into a SEM analysis.
    `lavaan.survey()` also accommodates replicate weights and multiply
    imputed datasets.
-   Package [vardpoor](../packages/vardpoor/index.html) allows to
    calculate linearisation of several nonlinear population statistics,
    variance estimation of sample surveys by the ultimate cluster
    method, variance estimation for longitudinal and cross-sectional
    measures and measures of change for any stage cluster sampling
    designs.
-   The package [rpms](../packages/rpms/index.html) fits a linear model
    to survey data in each node obtained by recursively partitioning the
    data. The algorithm accounts for one-stage of stratification and
    clustering as well as unequal probability of selection.
-   Package [svyPVpack](../packages/svyPVpack/index.html) extends
    package [survey](../packages/survey/index.html). This package deals
    with data which stem from survey designs and has been created to
    handle data from large scale assessments like PISA, PIAAC etc..

**Complex Survey Design: Calibration**

-   Package [survey](../packages/survey/index.html) allows for
    post-stratification, generalized raking/calibration, GREG estimation
    and trimming of weights.
-   The `calib()` function in package
    [sampling](../packages/sampling/index.html) allows to calibrate for
    nonresponse (with response homogeneity groups) for stratified
    samples.
-   The `calibWeights()` function in package
    [laeken](../packages/laeken/index.html) is a possible faster
    (depending on the example) implementation of parts of `calib()` from
    package [sampling](../packages/sampling/index.html).
-   The `calibSample()` function in package
    [simPop](../packages/simPop/index.html) is potential faster than the
    previous two mentioned functions, and it provides more
    user-friendlyness. `calibVars()` can be used to construct a matrix
    of binary variables for calibration. `calibPop()` is used to
    calibrate population person within household data using a simulated
    annealing approach.
-   Package [icarus](../packages/icarus/index.html) focuses on
    calibration and reweighting in survey sampling and was designed to
    provide a familiar setting in R for user of the SAS macro
    `             Calmar           `.
-   Package [reweight](../packages/reweight/index.html) allows for
    calibration of survey weights for categorical survey data so that
    the marginal distributions of certain variables fit more closely to
    those from a given population, but does not allow complex sampling
    designs.
-   The package [CalibrateSSB](../packages/CalibrateSSB/index.html)
    include a function to calculate weights and estimates for panel data
    with non-response.
-   Package [Frames2](../packages/Frames2/index.html) allows point and
    interval estimation in dual frame surveys. When two probability
    samples (one from each frame) are drawn. Information collected is
    suitably combined to get estimators of the parameter of interest.

**Editing and Visual Inspection of Microdata**

Editing tools:

-   Package [validate](../packages/validate/index.html) includes rule
    management and data validation and package
    [validatetools](../packages/validatetools/index.html) is checking
    and simplifying sets of validation rules.
-   Package [errorlocate](../packages/errorlocate/index.html) includes
    error localisation based on the principle of Fellegi and Holt. It
    supports categorical and/or numeric data and linear equalities,
    inequalities and conditional rules. The package includes a
    configurable backend for MIP-based error localization.
-   Package [editrules](../packages/editrules/index.html) convert
    readable linear (in)equalities into matrix form.
-   Package [deducorrect](../packages/deducorrect/index.html) depends on
    package [editrules](../packages/editrules/index.html) and applies
    deductive correction of simple rounding, typing and sign errors
    based on balanced edits. Values are changed so that the given
    balanced edits are fulfilled. To determine which values are changed
    the Levenstein-metric is applied.
-   The package [rspa](../packages/rspa/index.html) implements functions
    to minimally adjust numerical records so they obey (in)equation
    restrictions.
-   Package [SeleMix](../packages/SeleMix/index.html) can be used for
    selective editing for continuous scaled data. A mixture model
    (Gaussian contamination model) based on response(s) y and a depended
    set of covariates is fit to the data to quantify the impact of
    errors to the estimates.
-   Package [rrcovNA](../packages/rrcovNA/index.html) provides robust
    location and scatter estimation and robust principal component
    analysis with high breakdown point for incomplete data. It is
    therefore applicable to find representative and non-representative
    outliers.

Visual tools:

-   Package [VIM](../packages/VIM/index.html) is designed to visualize
    missing values using suitable plot methods. It can be used to
    analyse the structure of missing values in microdata using
    univariate, bivariate, multiple and multivariate plots where the
    information of missing values from specified variables are
    highlighted in selected variables. It also comes with a graphical
    user interface.
-   Package [tabplot](../packages/tabplot/index.html) provides the
    tableplot visualization method, which is used to profile or explore
    large statistical datasets. Up to a dozen of variables are shown
    column-wise as bar charts (numeric variables) or stacked bar charts
    (factors). Key aspects of the analysis with tableplots are the
    smoothness of a data distribution, the selective occurrence of
    missing values, and the distribution of correlated variables.
-   Package [treemap](../packages/treemap/index.html) provide treemaps.
    A treemap is a space-filling visualization of aggregates of data
    with hierarchical structures. Colors can be used to relate to
    highlight differences between comparable aggregates.

**Imputation**

A distinction between iterative model-based methods, k-nearest neighbor
methods and miscellaneous methods is made. However, often the criteria
for using a method depend on the scale of the data, which in official
statistics are typically a mixture of continuous, semi-continuous,
binary, categorical and count variables. In addition, measurement errors
may corrupt non-robust imputation methods. Note that only few imputation
methods can deal with mixed types of variables and only few methods
account for robustness issues.

EM-based Imputation Methods:

-   Package [mi](../packages/mi/index.html) provides iterative EM-based
    multiple Bayesian regression imputation of missing values and model
    checking of the regression models used. The regression models for
    each variable can also be user-defined. The data set may consist of
    continuous, semi-continuous, binary, categorical and/or count
    variables.
-   Package [mice](../packages/mice/index.html) provides iterative
    EM-based multiple regression imputation. The data set may consist of
    continuous, binary, categorical and/or count variables.
-   Package [mitools](../packages/mitools/index.html) provides tools to
    perform analyses and combine results from multiply-imputed datasets.
-   Package [Amelia](../packages/Amelia/index.html) provides multiple
    imputation where first bootstrap samples with the same dimensions as
    the original data are drawn, and then used for EM-based imputation.
    It is also possible to impute longitudinal data. The package in
    addition comes with a graphical user interface.
-   Package [VIM](../packages/VIM/index.html) provides EM-based multiple
    imputation (function `irmi()`) using robust estimations, which
    allows to adequately deal with data including outliers. It can
    handle data consisting of continuous, semi-continuous, binary,
    categorical and/or count variables.
-   Single imputation methods are included or called from other packages
    by the package [simputation](../packages/simputation/index.html). It
    supports regression (standard, M-estimation,
    ridge/lasso/elasticnet), hot-deck methods (powered by VIM),
    randomForest, EM-based, and iterative randomForest imputation.
-   Package [mix](../packages/mix/index.html) provides iterative
    EM-based multiple regression imputation. The data set may consist of
    continuous, binary or categorical variables, but methods for
    semi-continuous variables are missing.
-   Package [pan](../packages/pan/index.html) provides multiple
    imputation for multivariate panel or clustered data.
-   Package [norm](../packages/norm/index.html) provides EM-based
    multiple imputation for multivariate normal data.
-   Package [cat](../packages/cat/index.html) provides EM-based multiple
    imputation for multivariate categorical data.
-   Package [MImix](../packages/MImix/index.html) provides tools to
    combine results for multiply-imputed data using mixture
    approximations.
-   Package [robCompositions](../packages/robCompositions/index.html)
    provides iterative model-based imputation for compositional data
    (function `impCoda()`).
-   Package [missForest](../packages/missForest/index.html) uses the
    functionality of the randomForest to impute missing values in an
    iterative single-imputation fashion. It can deal with almost any
    kind of variables except semi-continuous ones. Even the underlying
    bootstrap approach of random forests ensures that from multiple runs
    one can get multiple imputations but the additional uncertainty of
    imputation is only considered when choosing the random forest method
    of package [mice](../packages/mice/index.html).

Nearest Neighbor Imputation Methods

-   Package [VIM](../packages/VIM/index.html) provides an implementation
    of the popular sequential and random (within a domain) hot-deck
    algorithm.
-   [VIM](../packages/VIM/index.html) also provides a fast k-nearest
    neighbor (knn) algorithm which can be used for large data sets. It
    uses a modification of the Gower Distance for numerical,
    categorical, ordered, continuous and semi-continuous variables.
-   Package [yaImpute](../packages/yaImpute/index.html) performs popular
    nearest neighbor routines for imputation of continuous variables
    where different metrics and methods can be used for determining the
    distance between observations.
-   Package [robCompositions](../packages/robCompositions/index.html)
    provides knn imputation for compositional data (function
    `impKNNa()`) using the Aitchison distance and adjustment of the
    nearest neighbor.
-   Package [rrcovNA](../packages/rrcovNA/index.html) provides an
    algorithm for (robust) sequential imputation (function `impSeq()`
    and `impSeqRob()` by minimizing the determinant of the covariance of
    the augmented data matrix. It\'s application is limited to
    continuous scaled data.
-   Package
    [[impute]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/impute.html)
    on Bioconductor impute provides knn imputation of continuous
    variables.

Copula-based Imputation Methods:

-   The S4 class package [CoImp](../packages/CoImp/index.html) imputes
    multivariate missing data by using conditional copula functions. The
    imputation procedure is semiparametric: the margins are
    non-parametrically estimated through local likelihood of low-degree
    polynomials while a range of different parametric models for the
    copula can be selected by the user. The missing values are imputed
    by drawing observations from the conditional density functions by
    means of the Hit or Miss Monte Carlo method. It works either for a
    matrix of continuous scaled variables or a matrix of discrete
    distributions.

Miscellaneous Imputation Methods:

-   Package [missMDA](../packages/missMDA/index.html) allows to impute
    incomplete continuous variables by principal component analysis
    (PCA) or categorical variables by multiple correspondence analysis
    (MCA).
-   Package [mice](../packages/mice/index.html) (function
    `mice.impute.pmm()`) and Package
    [Hmisc](../packages/Hmisc/index.html) (function `aregImpute()`)
    allow predictive mean matching imputation.
-   Package [VIM](../packages/VIM/index.html) allows to visualize the
    structure of missing values using suitable plot methods. It also
    comes with a graphical user interface.

**Statistical Disclosure Control**

Data from statistical agencies and other institutions are in its raw
form mostly confidential and data providers have to be ensure
confidentiality by both modifying the original data so that no
statistical units can be re-identified and by guaranteeing a minimum
amount of information loss.

-   Package [sdcMicro](../packages/sdcMicro/index.html) can be used for
    the generation of confidential (micro)data, i.e. for the generation
    of public- and scientific-use files. The package also comes with a
    graphical user interface.
-   Package [sdcTable](../packages/sdcTable/index.html) can be used to
    provide confidential (hierarchical) tabular data. It includes the
    HITAS and the HYPERCUBE technique and uses linear programming
    packages (Rglpk and lpSolveAPI) for solving (a large amount of)
    linear programs.
-   An interface to the package
    [sdcTable](../packages/sdcTable/index.html) is provided by package
    [easySdcTable](../packages/easySdcTable/index.html).

**Seasonal Adjustment and Forecasting**

For a more general view on time series methodology we refer to the
[TimeSeries](TimeSeries.html) task view. Only very specialized time
series packages related to complex surveys are discussed here.

-   Decomposition of time series can be done with the function
    `decompose()`, or more advanced by using the function `stl()`, both
    from the basic stats package. Decomposition is also possible with
    the `StructTS()` function, which can also be found in the stats
    package.
-   Many powerful tools can be accessed via packages
    [x12](../packages/x12/index.html) and
    [x12GUI](../packages/x12GUI/index.html) and package
    [seasonal](../packages/seasonal/index.html).
    [x12](../packages/x12/index.html) provides a wrapper function for
    the X12 binaries, which have to be installed first. It uses with a
    S4-class interface for batch processing of multiple time series.
    [x12GUI](../packages/x12GUI/index.html) provides a graphical user
    interface for the X12-Arima seasonal adjustment software. Less
    functionality but with the support of SEATS Spec is supported by
    package [seasonal](../packages/seasonal/index.html).
-   Given the large pool of individual forecasts in survey-type
    forecasting, forecast combination techniques from package
    [GeomComb](../packages/GeomComb/index.html) can be useful. It can
    also handle missing values in the time series.

**Statistical Matching and Record Linkage**

-   Package [StatMatch](../packages/StatMatch/index.html) provides
    functions to perform statistical matching between two data sources
    sharing a number of common variables. It creates a synthetic data
    set after matching of two data sources via a likelihood approach or
    via hot-deck.
-   Package [RecordLinkage](../packages/RecordLinkage/index.html)
    provides functions for linking and deduplicating data sets.
-   Package [MatchIt](../packages/MatchIt/index.html) allows nearest
    neighbor matching, exact matching, optimal matching and full
    matching amongst other matching methods. If two data sets have to be
    matched, the data must come as one data frame including a factor
    variable which includes information about the membership of each
    observation.
-   Package [stringdist](../packages/stringdist/index.html) can
    calculate various string distances based on edits
    (damerau-levenshtein, hamming, levenshtein, optimal sting
    alignment), qgrams (q-gram, cosine, jaccard distance) or heuristic
    metrics (jaro, jaro-winkler).
-   Package [XBRL](../packages/XBRL/index.html) allows the extraction of
    business financial information from XBRL Documents.

**Small Area Estimation**

-   Package [rsae](../packages/rsae/index.html) provides functions to
    estimate the parameters of the basic unit-level small area
    estimation (SAE) model (aka nested error regression model) by means
    of maximum likelihood (ML) or robust ML. On the basis of the
    estimated parameters, robust predictions of the area-specific means
    are computed (incl. MSE estimates; parametric bootstrap). The
    current version (rsae 0.4-x) does not allow for categorical
    independent variables.
-   Package [nlme](../packages/nlme/index.html) provides facilities to
    fit Gaussian linear and nonlinear mixed-effects models and
    [lme4](../packages/lme4/index.html) provides facilities to fit
    linear and generalized linear mixed-effects model, both used in
    small area estimation.
-   The [hbsae](../packages/hbsae/index.html) package provides functions
    to compute small area estimates based on a basic area or unit-level
    model. The model is fit using restricted maximum likelihood, or in a
    hierarchical Bayesian way. Auxilary information can be either counts
    resulting from categorical variables or means from continuous
    population information.
-   With package [JoSAE](../packages/JoSAE/index.html) point and
    variance estimation for the generalized regression (GREG) and a unit
    level empirical best linear unbiased prediction EBLUP estimators can
    be made at domain level. It basically provides wrapper functions to
    the [nlme](../packages/nlme/index.html) package that is used to fit
    the basic random effects models.

**Indices, Indicators, Tables and Visualisation of Indicators**

-   Package [laeken](../packages/laeken/index.html) provides functions
    to estimate popular risk-of-poverty and inequality indicators
    (at-risk-of-poverty rate, quintile share ratio, relative median
    risk-of-poverty gap, Gini coefficient). In addition, standard and
    robust methods for tail modeling of Pareto distributions are
    provided for semi-parametric estimation of indicators from
    continuous univariate distributions such as income variables.
-   Package [convey](../packages/convey/index.html) estimates variances
    on indicators of income concentration and poverty using familiar
    linearized and replication-based designs created by the
    [survey](../packages/survey/index.html) package such as the Gini
    coefficient, Atkinson index, at-risk-of-poverty threshold, and more
    than a dozen others.
-   Package [ineq](../packages/ineq/index.html) computes various
    inequality measures (Gini, Theil, entropy, among others),
    concentration measures (Herfindahl, Rosenbluth), and poverty
    measures (Watts, Sen, SST, and Foster). It also computes and draws
    empirical and theoretical Lorenz curves as well as Pen\'s parade. It
    is not designed to deal with sampling weights directly (these could
    only be emulated via `rep(x, weights)`).
-   Package [IC2](../packages/IC2/index.html) include three inequality
    indices: extended Gini, Atkinson and Generalized Entropy. It can
    deal with sampling weights and subgroup decomposition is supported.
-   Package [DHS.rates](../packages/DHS.rates/index.html) estimates key
    indicators (especially fertility rates) and their variances for the
    Demographic and Health Survey (DHS) data.
-   Functions `priceIndex()` from package
    [micEconIndex](../packages/micEconIndex/index.html) allows to
    estimate the Paasche, the Fisher and the Laspeyres price indices.
    For estimating quantities (of goods, for example), function
    `quantityIndex()` might be your friend.
-   Package [tmap](../packages/tmap/index.html) offers a layer-based way
    to make thematic maps, like choropleths and bubble maps.
-   Package [rworldmap](../packages/rworldmap/index.html) outline how to
    map country referenced data and support users in visualising their
    own data. Examples are given, e.g., maps for the world bank and UN.
    It provides also new ways to visualise maps.
-   Package [rrcov3way](../packages/rrcov3way/index.html) provides
    robust methods for multiway data analysis, applicable also for
    compositional data.
-   Package [robCompositions](../packages/robCompositions/index.html)
    methods for compositional tables including statistical tests.

**Microsimulation**

-   Using package [simPop](../packages/simPop/index.html) one can
    simulate populations from surveys based on auxiliary data with
    model-based methods or synthetic reconstruction methods.
    Hiercharical and cluster structures (such as households) can be
    considered as well as the methods takes account for samples
    collected based on complex sample designs. Calibration tools
    (iterative proportional fitting, iterative proportional updating)
    and combinatorial optimization tools (simulated annealing) are also
    available. The code is optimized for fast computations. The package
    based on a S4 class implementation. The simulated population can
    serve as basis data for microsimulation studies.
-   The [MicSim](../packages/MicSim/index.html) package includes methods
    for microsimulations. Given a initial population, mortality rates,
    divorce rates, marriage rates, education changes, etc. and their
    transition matrix can be defined and included for the simulation of
    future states of the population. The package does not contain
    compiled code but functionality to run the microsimulation in
    parallel is provided.
-   Package [sms](../packages/sms/index.html) provides facilities to
    simulate micro-data from given area-based macro-data. Simulated
    annealing is used to best satisfy the available description of an
    area. For computational issues, the calculations can be run in
    parallel mode.
-   Package [synthpop](../packages/synthpop/index.html) using regression
    tree methods to simulate synthetic data from given data. It is
    suitable to produce synthetic data when the data have no
    hierarchical and cluster information (such as households) as well as
    when the data does not collected with a complex sampling design.
-   Package [saeSim](../packages/saeSim/index.html) Tools for the
    simulation of data in the context of small area estimation.

**Additional Packages and Functionalities**

Various additional packages are available that provides certain
functionality useful in official statistics and survey methodology.

-   The [questionr](../packages/questionr/index.html) package contains a
    set of functions to make the processing and analysis of surveys
    easier. It provides interactive shiny apps and addins for data
    recoding, contingency tables, dataset metadata handling, and several
    convenience functions.

Data Import and Export:

-   Package [SAScii](../packages/SAScii/index.html) imports ASCII files
    directly into R using only a SAS input script, which is parsed and
    converted into arguments for a read.fwf call. This is useful
    whenever SAS scripts for importing data are already available.
-   The [foreign](../packages/foreign/index.html) package includes tools
    for reading data from SAS Xport (function `read.xport()`), Stata
    (function `read.dta()`), SPSS (function `read.spss()`) and various
    other formats. It provides facilities to write file to various
    formats, see function `write.foreign()`.
-   Also the package [haven](../packages/haven/index.html) imports and
    exports SAS, Stata and SPSS (function `read.spss()`) files. The
    package is more efficient for loading heavy data sets and it handles
    the labelling of variables and values in an advanced manner.
-   Also the package [Hmisc](../packages/Hmisc/index.html) provides
    tools to read data sets from SPSS (function `spss.get()`) or Stata
    (function `stata.get()`).
-   The [pxR](../packages/pxR/index.html) package provides a set of
    functions for reading and writing PC-Axis files, used by different
    statistical organizations around the globe for dissemination of
    their (multidimensional) tables.
-   With package [prevR](../packages/prevR/index.html) and it\'s
    function `import.dhs()` it is possible to directly imports data from
    the Demographic Health Survey.
-   Function `describe()` from package
    [questionr](../packages/questionr/index.html) describes the
    variables of a dataset that might include labels imported with the
    foreign or memisc packages.
-   Package [OECD](../packages/OECD/index.html) searches and extracts
    data from the OECD.
-   Package [Rilostat](../packages/Rilostat/index.html) contains tools
    to download data from the [international labour organisation
    database](http://www.ilo.org/ilostat) together with search and
    manipulation utilities. It can also import ilostat data that are
    available on their data base in SDMX format.
-   Access to Finnish open government data is provided by package
    [sorvi](../packages/sorvi/index.html)
-   Tools to download data from the Eurostat database together with
    search and manipulation utilities are included in package
    [eurostat](../packages/eurostat/index.html).
-   Package [acs](../packages/acs/index.html) downloads, manipulates,
    and presents the American Community Survey and decennial data from
    the US Census.
-   Access to data published by INEGI, Mexico\'s official statistics
    agency, is supported by package
    [inegiR](../packages/inegiR/index.html)
-   Package [cbsodataR](../packages/cbsodataR/index.html) provides
    access to Statistics Netherlands\' (CBS) open data API.
-   A wrapper for the U.S. Census Bureau APIs that returns data frames
    of Census data and metadata is implemented in package
    [censusapi](../packages/censusapi/index.html).

Misc:

-   Package [samplingbook](../packages/samplingbook/index.html) includes
    sampling procedures from the book \'Stichproben. Methoden und
    praktische Umsetzung mit R\' by Goeran Kauermann and Helmut
    Kuechenhoff (2010).
-   Package [SDaA](../packages/SDaA/index.html) is designed to reproduce
    results from Lohr, S. (1999) \'Sampling: Design and Analysis,
    Duxbury\' and includes the data sets from this book.
-   The main contributions of
    [samplingVarEst](../packages/samplingVarEst/index.html) are
    Jackknife alternatives for variance estimation of unequal
    probability with one or two stage designs.
-   Package [TeachingSampling](../packages/TeachingSampling/index.html)
    includes functionality for sampling designs and parameter estimation
    in finite populations.
-   Package [memisc](../packages/memisc/index.html) includes tools for
    the management of survey data, graphics and simulation.
-   Package [spsurvey](../packages/spsurvey/index.html) includes
    facilities for spatial survey design and analysis for equal and
    unequal probability (stratified) sampling.
-   The [FFD](../packages/FFD/index.html) package is designed to
    calculate optimal sample sizes of a population of animals living in
    herds for surveys to substantiate freedom from disease. The criteria
    of estimating the sample sizes take the herd-level clustering of
    diseases as well as imperfect diagnostic tests into account and
    select the samples based on a two-stage design. Inclusion
    probabilities are not considered in the estimation. The package
    provides a graphical user interface as well.
-   [mipfp](../packages/mipfp/index.html) provides multidimensional
    iterative proportional fitting to calibrate n-dimensional arrays
    given target marginal tables.
-   Package [MBHdesign](../packages/MBHdesign/index.html) provides
    spatially balanced designs from a set of (contiguous) potential
    sampling locations in a study region.
-   Package [quantification](../packages/quantification/index.html)
    provides different functions for quantifying qualitative survey
    data. It supports the Carlson-Parkin method, the regression
    approach, the balance approach and the conditional expectations
    method.
-   [BIFIEsurvey](../packages/BIFIEsurvey/index.html) includes tools for
    survey statistics in educational assessment including data with
    replication weights (e.g. from bootstrap).
-   [surveybootstrap](../packages/surveybootstrap/index.html) includes
    tools for using different kinds of bootstrap for estimating sampling
    variation using complex survey data.
-   Package [surveyoutliers](../packages/surveyoutliers/index.html)
    winsorize values of a variable of interest.
-   The package [univOutl](../packages/univOutl/index.html) includes
    various methods for detecting univariate outliers, e.g. the
    Hidiroglou-Berthelot method.
-   Package [extremevalues](../packages/extremevalues/index.html) is
    designed to detect univariate outliers based on modeling the bulk
    distribution.
-   Package [RRTCS](../packages/RRTCS/index.html) includes randomized
    response techniques for complex surveys.
-   Package [panelaggregation](../packages/panelaggregation/index.html)
    aggregates business tendency survey data (and other qualitative
    surveys) to time series at various aggregation levels.
-   Package [surveydata](../packages/surveydata/index.html) makes it
    easy to keep track of metadata from surveys, and to easily extract
    columns with specific questions.
-   [RcmdrPlugin.sampling](../packages/RcmdrPlugin.sampling/index.html)
    includes tools for sampling in official statistical surveys. It
    includes tools for calculating sample sizes and selecting samples
    using various sampling designs.
-   Package [mapStats](../packages/mapStats/index.html) does automated
    calculation and visualization of survey data statistics on a
    color-coded map.

</div>

### CRAN packages:

-   [acs](../packages/acs/index.html)
-   [Amelia](../packages/Amelia/index.html)
-   [BalancedSampling](../packages/BalancedSampling/index.html)
-   [BIFIEsurvey](../packages/BIFIEsurvey/index.html)
-   [CalibrateSSB](../packages/CalibrateSSB/index.html)
-   [cat](../packages/cat/index.html)
-   [cbsodataR](../packages/cbsodataR/index.html)
-   [censusapi](../packages/censusapi/index.html)
-   [CoImp](../packages/CoImp/index.html)
-   [convey](../packages/convey/index.html)
-   [deducorrect](../packages/deducorrect/index.html)
-   [DHS.rates](../packages/DHS.rates/index.html)
-   [easySdcTable](../packages/easySdcTable/index.html)
-   [editrules](../packages/editrules/index.html)
-   [errorlocate](../packages/errorlocate/index.html)
-   [eurostat](../packages/eurostat/index.html)
-   [extremevalues](../packages/extremevalues/index.html)
-   [FFD](../packages/FFD/index.html)
-   [foreign](../packages/foreign/index.html)
-   [Frames2](../packages/Frames2/index.html)
-   [GeomComb](../packages/GeomComb/index.html)
-   [gridsample](../packages/gridsample/index.html)
-   [haven](../packages/haven/index.html)
-   [hbsae](../packages/hbsae/index.html)
-   [Hmisc](../packages/Hmisc/index.html)
-   [IC2](../packages/IC2/index.html)
-   [icarus](../packages/icarus/index.html)
-   [inegiR](../packages/inegiR/index.html)
-   [ineq](../packages/ineq/index.html)
-   [JoSAE](../packages/JoSAE/index.html)
-   [laeken](../packages/laeken/index.html)
-   [lavaan](../packages/lavaan/index.html)
-   [lavaan.survey](../packages/lavaan.survey/index.html)
-   [lme4](../packages/lme4/index.html)
-   [mapStats](../packages/mapStats/index.html)
-   [MatchIt](../packages/MatchIt/index.html)
-   [MBHdesign](../packages/MBHdesign/index.html)
-   [memisc](../packages/memisc/index.html)
-   [mi](../packages/mi/index.html)
-   [mice](../packages/mice/index.html)
-   [micEconIndex](../packages/micEconIndex/index.html)
-   [MicSim](../packages/MicSim/index.html)
-   [MImix](../packages/MImix/index.html)
-   [mipfp](../packages/mipfp/index.html)
-   [missForest](../packages/missForest/index.html)
-   [missMDA](../packages/missMDA/index.html)
-   [mitools](../packages/mitools/index.html)
-   [mix](../packages/mix/index.html)
-   [nlme](../packages/nlme/index.html)
-   [norm](../packages/norm/index.html)
-   [OECD](../packages/OECD/index.html)
-   [pan](../packages/pan/index.html)
-   [panelaggregation](../packages/panelaggregation/index.html)
-   [pps](../packages/pps/index.html)
-   [PracTools](../packages/PracTools/index.html)
-   [prevR](../packages/prevR/index.html)
-   [pxR](../packages/pxR/index.html)
-   [quantification](../packages/quantification/index.html)
-   [questionr](../packages/questionr/index.html)
-   [RcmdrPlugin.sampling](../packages/RcmdrPlugin.sampling/index.html)
-   [RecordLinkage](../packages/RecordLinkage/index.html)
-   [reweight](../packages/reweight/index.html)
-   [Rilostat](../packages/Rilostat/index.html)
-   [robCompositions](../packages/robCompositions/index.html)
-   [rpms](../packages/rpms/index.html)
-   [rrcov3way](../packages/rrcov3way/index.html)
-   [rrcovNA](../packages/rrcovNA/index.html)
-   [RRTCS](../packages/RRTCS/index.html)
-   [rsae](../packages/rsae/index.html)
-   [rspa](../packages/rspa/index.html)
-   [rworldmap](../packages/rworldmap/index.html)
-   [saeSim](../packages/saeSim/index.html)
-   [samplesize4surveys](../packages/samplesize4surveys/index.html)
-   [sampling](../packages/sampling/index.html)
-   [samplingbook](../packages/samplingbook/index.html)
-   [SamplingStrata](../packages/SamplingStrata/index.html)
-   [samplingVarEst](../packages/samplingVarEst/index.html)
-   [SAScii](../packages/SAScii/index.html)
-   [SDaA](../packages/SDaA/index.html)
-   [sdcMicro](../packages/sdcMicro/index.html)
-   [sdcTable](../packages/sdcTable/index.html)
-   [seasonal](../packages/seasonal/index.html)
-   [SeleMix](../packages/SeleMix/index.html)
-   [simFrame](../packages/simFrame/index.html)
-   [simPop](../packages/simPop/index.html)
-   [simputation](../packages/simputation/index.html)
-   [sms](../packages/sms/index.html)
-   [sorvi](../packages/sorvi/index.html)
-   [spsurvey](../packages/spsurvey/index.html)
-   [srvyr](../packages/srvyr/index.html)
-   [StatMatch](../packages/StatMatch/index.html)
-   [stratification](../packages/stratification/index.html)
-   [stringdist](../packages/stringdist/index.html)
-   [survey](../packages/survey/index.html) (core)
-   [surveybootstrap](../packages/surveybootstrap/index.html)
-   [surveydata](../packages/surveydata/index.html)
-   [surveyoutliers](../packages/surveyoutliers/index.html)
-   [surveyplanning](../packages/surveyplanning/index.html)
-   [svyPVpack](../packages/svyPVpack/index.html)
-   [synthpop](../packages/synthpop/index.html)
-   [tabplot](../packages/tabplot/index.html)
-   [TeachingSampling](../packages/TeachingSampling/index.html)
-   [tmap](../packages/tmap/index.html)
-   [treemap](../packages/treemap/index.html)
-   [univOutl](../packages/univOutl/index.html)
-   [validate](../packages/validate/index.html)
-   [validatetools](../packages/validatetools/index.html)
-   [vardpoor](../packages/vardpoor/index.html)
-   [VIM](../packages/VIM/index.html)
-   [x12](../packages/x12/index.html)
-   [x12GUI](../packages/x12GUI/index.html)
-   [XBRL](../packages/XBRL/index.html)
-   [yaImpute](../packages/yaImpute/index.html)

### Related links:

-   CRAN Task View: [TimeSeries](TimeSeries.html)
-   CRAN Task View: [SocialSciences](SocialSciences.html)
-   Bioconductor Package:
    [[impute]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/impute.html)
-   [useR!2008 Tutorial: Small Area
    Estimation](http://www.R-project.org/conferences/useR-2008/tutorials/gomez.html)
