## CRAN Task View: Econometrics

  ----------------- ------------------------------------------------
  **Maintainer:**   Achim Zeileis
  **Contact:**      Achim.Zeileis at R-project.org
  **Version:**      2018-04-24
  **URL:**          <https://CRAN.R-project.org/view=Econometrics>
  ----------------- ------------------------------------------------

<div>

Base R ships with a lot of functionality useful for computational
econometrics, in particular in the stats package. This functionality is
complemented by many packages on CRAN, a brief overview is given below.
There is also a considerable overlap between the tools for econometrics
in this view and those in the task views on [Finance](Finance.html),
[SocialSciences](SocialSciences.html), and
[TimeSeries](TimeSeries.html). Furthermore, the [Finance
SIG](https://stat.ethz.ch/mailman/listinfo/R-SIG-Finance/) is a suitable
mailing list for obtaining help and discussing questions about both
computational finance and econometrics.

The packages in this view can be roughly structured into the following
topics. If you think that some package is missing from the list, please
contact the maintainer.

**Basic linear regression**

-   *Estimation and standard inference* : Ordinary least squares (OLS)
    estimation for linear models is provided by `lm()` (from stats) and
    standard tests for model comparisons are available in various
    methods such as `summary()` and `anova()`.
-   *Further inference and nested model comparisons* : Functions
    analogous to the basic `summary()` and `anova()` methods that also
    support asymptotic tests ( *z* instead of *t* tests, and Chi-squared
    instead of *F* tests) and plug-in of other covariance matrices are
    `coeftest()` and `waldtest()` in
    [lmtest](../packages/lmtest/index.html). Tests of more general
    linear hypotheses are implemented in `linearHypothesis()` and for
    nonlinear hypotheses in `deltaMethod()` in
    [car](../packages/car/index.html).
-   *Robust standard errors* : HC and HAC covariance matrices are
    available in [sandwich](../packages/sandwich/index.html) and can be
    plugged into the inference functions mentioned above.
-   *Nonnested model comparisons* : Various tests for comparing
    non-nested linear models are available in
    [lmtest](../packages/lmtest/index.html) (encompassing test, J test,
    Cox test). The Vuong test for comparing other non-nested models is
    provided by [nonnest2](../packages/nonnest2/index.html) (and
    specifically for count data regression in
    [pscl](../packages/pscl/index.html)).
-   *Diagnostic checking* : The packages
    [car](../packages/car/index.html) and
    [lmtest](../packages/lmtest/index.html) provide a large collection
    of regression diagnostics and diagnostic tests.

**Microeconometrics**

-   *Generalized linear models (GLMs)* : Many standard microeconometric
    models belong to the family of generalized linear models and can be
    fitted by `glm()` from package stats. This includes in particular
    logit and probit models for modeling choice data and Poisson models
    for count data. Effects for typical values of regressors in these
    models can be obtained and visualized using
    [effects](../packages/effects/index.html). Marginal effects tables
    for certain GLMs can be obtained using the
    [mfx](../packages/mfx/index.html) and
    [margins](../packages/margins/index.html) packages. Interactive
    visualizations of both effects and marginal effects are possible in
    [LinRegInteractive](../packages/LinRegInteractive/index.html).
-   *Binary responses* : The standard logit and probit models (among
    many others) for binary responses are GLMs that can be estimated by
    `glm()` with `family = binomial`. Bias-reduced GLMs that are robust
    to complete and quasi-complete separation are provided by
    [brglm](../packages/brglm/index.html). Discrete choice models
    estimated by simulated maximum likelihood are implemented in
    [Rchoice](../packages/Rchoice/index.html). Heteroscedastic probit
    models (and other heteroscedastic GLMs) are implemented in
    [glmx](../packages/glmx/index.html) along with parametric link
    functions and goodness-of-link tests for GLMs.
-   *Count responses* : The basic Poisson regression is a GLM that can
    be estimated by `glm()` with `family = poisson` as explained above.
    Negative binomial GLMs are available via `glm.nb()` in package
    [MASS](../packages/MASS/index.html). Another implementation of
    negative binomial models is provided by
    [aod](../packages/aod/index.html), which also contains other models
    for overdispersed data. Zero-inflated and hurdle count models are
    provided in package [pscl](../packages/pscl/index.html). A
    reimplementation by the same authors is currently under development
    in
    [[countreg]{.Rforge}](https://R-Forge.R-project.org/projects/countreg/)
    on R-Forge which also encompasses separate functions for
    zero-truncated regression, finite mixture models etc.
-   *Multinomial responses* : Multinomial models with
    individual-specific covariates only are available in `multinom()`
    from package [nnet](../packages/nnet/index.html). Implementations
    with both individual- and choice-specific variables are
    [mlogit](../packages/mlogit/index.html) and
    [mnlogit](../packages/mnlogit/index.html). Generalized multinomial
    logit models (e.g., with random effects etc.) are in
    [gmnl](../packages/gmnl/index.html). Generalized additive models
    (GAMs) for multinomial responses can be fitted with the
    [VGAM](../packages/VGAM/index.html) package. A Bayesian approach to
    multinomial probit models is provided by
    [MNP](../packages/MNP/index.html). Various Bayesian multinomial
    models (including logit and probit) are available in
    [bayesm](../packages/bayesm/index.html). Furthermore, the package
    [RSGHB](../packages/RSGHB/index.html) fits various hierarchical
    Bayesian specifications based on direct specification of the
    likelihood function.
-   *Ordered responses* : Proportional-odds regression for ordered
    responses is implemented in `polr()` from package
    [MASS](../packages/MASS/index.html). The package
    [ordinal](../packages/ordinal/index.html) provides cumulative link
    models for ordered data which encompasses proportional odds models
    but also includes more general specifications. Bayesian ordered
    probit models are provided by
    [bayesm](../packages/bayesm/index.html).
-   *Censored responses* : Basic censored regression models (e.g., tobit
    models) can be fitted by `survreg()` in
    [survival](../packages/survival/index.html), a convenience interface
    `tobit()` is in package [AER](../packages/AER/index.html). Further
    censored regression models, including models for panel data, are
    provided in [censReg](../packages/censReg/index.html). Interval
    regression models are in [intReg](../packages/intReg/index.html).
    Censored regression models with conditional heteroscedasticity are
    in [crch](../packages/crch/index.html). Furthermore, hurdle models
    for left-censored data at zero can be estimated with
    [mhurdle](../packages/mhurdle/index.html). Models for sample
    selection are available in
    [sampleSelection](../packages/sampleSelection/index.html) and
    semiparametric extensions of these are provided by
    [SemiParSampleSel](../packages/SemiParSampleSel/index.html). Package
    [matchingMarkets](../packages/matchingMarkets/index.html) corrects
    for selection bias when the sample is the result of a stable
    matching process (e.g., a group formation or college admissions
    problem).
-   *Truncated responses* : [crch](../packages/crch/index.html) for
    truncated (and potentially heteroscedastic) Gaussian, logistic, and
    t responses. Homoscedastic Gaussian responses are also available in
    [truncreg](../packages/truncreg/index.html).
-   *Fraction and proportion responses* : Fractional response models are
    in [frm](../packages/frm/index.html). Beta regression for responses
    in (0, 1) is in [betareg](../packages/betareg/index.html) and
    [gamlss](../packages/gamlss/index.html).
-   *Miscellaneous* : Further more refined tools for microeconometrics
    are provided in the [micEcon](../packages/micEcon/index.html) family
    of packages: Analysis with Cobb-Douglas, translog, and quadratic
    functions is in [micEcon](../packages/micEcon/index.html); the
    constant elasticity of scale (CES) function is in
    [micEconCES](../packages/micEconCES/index.html); the symmetric
    normalized quadratic profit (SNQP) function is in
    [micEconSNQP](../packages/micEconSNQP/index.html). The almost ideal
    demand system (AIDS) is in
    [micEconAids](../packages/micEconAids/index.html). Stochastic
    frontier analysis (SFA) is in
    [frontier](../packages/frontier/index.html) and certain special
    cases also in [sfa](../packages/sfa/index.html). Semiparametric SFA
    in is available in [semsfa](../packages/semsfa/index.html) and
    spatial SFA in [spfrontier](../packages/spfrontier/index.html) and
    [ssfa](../packages/ssfa/index.html). The package
    [bayesm](../packages/bayesm/index.html) implements a Bayesian
    approach to microeconometrics and marketing. Estimation and marginal
    effect computations for multivariate probit models can be carried
    out with [mvProbit](../packages/mvProbit/index.html). Inference for
    relative distributions is contained in package
    [reldist](../packages/reldist/index.html).

**Instrumental variables**

-   *Basic instrumental variables (IV) regression* : Two-stage least
    squares (2SLS) is provided by `ivreg()` in
    [AER](../packages/AER/index.html). Other implementations are in
    `tsls()` in package [sem](../packages/sem/index.html), in
    [ivpack](../packages/ivpack/index.html), and
    [lfe](../packages/lfe/index.html) (with particular focus on multiple
    group fixed effects).
-   *Binary responses* : An IV probit model via GLS estimation is
    available in [ivprobit](../packages/ivprobit/index.html). The
    [LARF](../packages/LARF/index.html) package estimates local average
    response functions for binary treatments and binary instruments.
-   *Panel data* : Certain basic IV models for panel data can also be
    estimated with standard 2SLS functions (see above). Dedicated IV
    panel data models are provided by
    [ivfixed](../packages/ivfixed/index.html) (fixed effects) and
    [ivpanel](../packages/ivpanel/index.html) (between and random
    effects).
-   *Miscellaneous* : [REndo](../packages/REndo/index.html) fits linear
    models with endogenous regressor using various latent instrumental
    variable approaches. [ivbma](../packages/ivbma/index.html) estimates
    Bayesian IV models with conditional Bayes factors.
    [ivlewbel](../packages/ivlewbel/index.html) implements the Lewbel
    approach based on GMM estimation of triangular systems using
    heteroscedasticity-based IVs.

**Panel data models**

-   *Panel-corrected standard errors* : A simple approach for panel data
    is to fit the pooling (or independence) model (e.g., via `lm()` or
    `glm()`) and only correct the standard errors. Different types of
    panel-corrected standard errors are available in
    [multiwayvcov](../packages/multiwayvcov/index.html),
    [clusterSEs](../packages/clusterSEs/index.html),
    [pcse](../packages/pcse/index.html),
    [clubSandwich](../packages/clubSandwich/index.html),
    [plm](../packages/plm/index.html), and
    [geepack](../packages/geepack/index.html), respectively. The latter
    two require estimation of the pooling/independence models via
    `plm()` and `geeglm()` from the respective packages (which also
    provide other types of models, see below).
-   *Linear panel models* : [plm](../packages/plm/index.html), providing
    a wide range of within, between, and random-effect methods (among
    others) along with corrected standard errors, tests, etc. Another
    implementation of several of these models is in
    [Paneldata](../packages/Paneldata/index.html). Various dynamic panel
    models are available in [plm](../packages/plm/index.html) and
    dynamic panel models with fixed effects in
    [OrthoPanels](../packages/OrthoPanels/index.html). Panel vector
    autoregressions are implemented in
    [panelvar](../packages/panelvar/index.html).
-   *Generalized estimation equations and GLMs* : GEE models for panel
    data (or longitudinal data in statistical jargon) are in
    [geepack](../packages/geepack/index.html). The
    [pglm](../packages/pglm/index.html) package provides estimation of
    GLM-like models for panel data.
-   *Mixed effects models* : Linear and nonlinear models for panel data
    (and more general multi-level data) are available in
    [lme4](../packages/lme4/index.html) and
    [nlme](../packages/nlme/index.html).
-   *Instrumental variables* : [ivfixed](../packages/ivfixed/index.html)
    and [ivpanel](../packages/ivpanel/index.html), see also above.
-   *Heterogeneous time trends* : [phtt](../packages/phtt/index.html)
    offers the possibility of analyzing panel data with large dimensions
    n and T and can be considered when the unobserved heterogeneity
    effects are time-varying.
-   *Miscellaneous* : Multiple group fixed effects are in
    [lfe](../packages/lfe/index.html). Autocorrelation and
    heteroscedasticity correction in are available in
    [wahc](../packages/wahc/index.html) and
    [panelAR](../packages/panelAR/index.html). PANIC Tests of
    nonstationarity are in [PANICr](../packages/PANICr/index.html).
    Threshold regression and unit root tests are in
    [pdR](../packages/pdR/index.html). The panel data approach method
    for program evaluation is available in
    [pampe](../packages/pampe/index.html).

**Further regression models**

-   *Nonlinear least squares modeling* : `nls()` in package stats.
-   *Quantile regression* : [quantreg](../packages/quantreg/index.html)
    (including linear, nonlinear, censored, locally polynomial and
    additive quantile regressions).
-   *Generalized method of moments (GMM) and generalized empirical
    likelihood (GEL)* : [gmm](../packages/gmm/index.html).
-   *Spatial econometric models* : The [Spatial](Spatial.html) view
    gives details about handling spatial data, along with information
    about (regression) modeling. In particular, spatial regression
    models can be fitted using [spdep](../packages/spdep/index.html) and
    [sphet](../packages/sphet/index.html) (the latter using a GMM
    approach). [splm](../packages/splm/index.html) is a package for
    spatial panel models. Spatial probit models are available in
    [spatialprobit](../packages/spatialprobit/index.html).
-   *Bayesian model averaging (BMA)* : A comprehensive toolbox for BMA
    is provided by [BMS](../packages/BMS/index.html) including flexible
    prior selection, sampling, etc. A different implementation is in
    [BMA](../packages/BMA/index.html) for linear models, generalizable
    linear models and survival models (Cox regression).
-   *Linear structural equation models* :
    [lavaan](../packages/lavaan/index.html) and
    [sem](../packages/sem/index.html). See also the
    [Psychometrics](Psychometrics.html) task view for more details.
-   *Simultaneous equation estimation* :
    [systemfit](../packages/systemfit/index.html).
-   *Nonparametric kernel methods* : [np](../packages/np/index.html).
-   *Linear and nonlinear mixed-effect models* :
    [nlme](../packages/nlme/index.html) and
    [lme4](../packages/lme4/index.html).
-   *Generalized additive models (GAMs)* :
    [mgcv](../packages/mgcv/index.html),
    [gam](../packages/gam/index.html),
    [gamlss](../packages/gamlss/index.html) and
    [VGAM](../packages/VGAM/index.html).
-   *Extreme bounds analysis* :
    [ExtremeBounds](../packages/ExtremeBounds/index.html).
-   *Miscellaneous* : The packages [VGAM](../packages/VGAM/index.html),
    [rms](../packages/rms/index.html) and
    [Hmisc](../packages/Hmisc/index.html) provide several tools for
    extended handling of (generalized) linear regression models.
    [Zelig](../packages/Zelig/index.html) is a unified easy-to-use
    interface to a wide range of regression models.

**Time series data and models**

-   The [TimeSeries](TimeSeries.html) task view provides much more
    detailed information about both basic time series infrastructure and
    time series models. Here, only the most important aspects relating
    to econometrics are briefly mentioned. Time series models for
    financial econometrics (e.g., GARCH, stochastic volatility models,
    or stochastic differential equations, etc.) are described in the
    [Finance](Finance.html) task view.
-   *Infrastructure for regularly spaced time series* : The class `"ts"`
    in package stats is R\'s standard class for regularly spaced time
    series (especially annual, quarterly, and monthly data). It can be
    coerced back and forth without loss of information to `"zooreg"`
    from package [zoo](../packages/zoo/index.html).
-   *Infrastructure for irregularly spaced time series* :
    [zoo](../packages/zoo/index.html) provides infrastructure for both
    regularly and irregularly spaced time series (the latter via the
    class `"zoo"`) where the time information can be of arbitrary class.
    This includes daily series (typically with `"Date"` time index) or
    intra-day series (e.g., with `"POSIXct"` time index). An extension
    based on [zoo](../packages/zoo/index.html) geared towards time
    series with different kinds of time index is
    [xts](../packages/xts/index.html). Further packages aimed
    particularly at finance applications are discussed in the
    [Finance](Finance.html) task view.
-   *Classical time series models* : Simple autoregressive models can be
    estimated with `ar()` and ARIMA modeling and Box-Jenkins-type
    analysis can be carried out with `arima()` (both in the stats
    package). An enhanced version of `arima()` is in
    [forecast](../packages/forecast/index.html).
-   *Linear regression models* : A convenience interface to `lm()` for
    estimating OLS and 2SLS models based on time series data is
    [dynlm](../packages/dynlm/index.html). Linear regression models with
    AR error terms via GLS is possible using `gls()` from
    [nlme](../packages/nlme/index.html).
-   *Structural time series models* : Standard models can be fitted with
    `StructTS()` in stats. Further packages are discussed in the
    [TimeSeries](TimeSeries.html) task view.
-   *Filtering and decomposition* : `decompose()` and `HoltWinters()` in
    stats. The basic function for computing filters (both rolling and
    autoregressive) is `filter()` in stats. Many extensions to these
    methods, in particular for forecasting and model selection, are
    provided in the [forecast](../packages/forecast/index.html) package.
-   *Vector autoregression* : Simple models can be fitted by `ar()` in
    stats, more elaborate models are provided in package
    [vars](../packages/vars/index.html) along with suitable diagnostics,
    visualizations etc. A Bayesian approach is available in
    [MSBVAR](../packages/MSBVAR/index.html) and panel vector
    autoregressions in [panelvar](../packages/panelvar/index.html).
-   *Unit root and cointegration tests* :
    [urca](../packages/urca/index.html),
    [tseries](../packages/tseries/index.html),
    [CADFtest](../packages/CADFtest/index.html). See also
    [pco](../packages/pco/index.html) for panel cointegration tests.
-   *Miscellaneous* :
    -   [tsDyn](../packages/tsDyn/index.html) - Threshold and smooth
        transition models.
    -   [midasr](../packages/midasr/index.html) - *MIDAS regression* and
        other econometric methods for mixed frequency time series data
        analysis.
    -   [gets](../packages/gets/index.html) - GEneral-To-Specific (GETS)
        model selection for either ARX models with log-ARCH-X errors, or
        a log-ARCH-X model of the log variance.
    -   [tsfa](../packages/tsfa/index.html) - Time series factor
        analysis.
    -   [dlsem](../packages/dlsem/index.html) - Distributed-lag linear
        structural equation models.
    -   [apt](../packages/apt/index.html) - Asymmetric price
        transmission models.

**Data sets**

-   *Textbooks and journals* : Packages
    [AER](../packages/AER/index.html),
    [Ecdat](../packages/Ecdat/index.html), and
    [wooldridge](../packages/wooldridge/index.html) contain a
    comprehensive collections of data sets from various standard
    econometric textbooks as well as several data sets from the Journal
    of Applied Econometrics and the Journal of Business & Economic
    Statistics data archives. [AER](../packages/AER/index.html) and
    [wooldridge](../packages/wooldridge/index.html) additionally provide
    extensive sets of examples reproducing analyses from the
    textbooks/papers, illustrating various econometric methods. In
    [pder](../packages/pder/index.html) a wide collection of data sets
    for panel data econometrics with R is available.
-   *Canadian monetary aggregates* :
    [CDNmoney](../packages/CDNmoney/index.html).
-   *Penn World Table* : [pwt](../packages/pwt/index.html) provides
    versions 5.6, 6.x, 7.x. Version 8.x and 9.x data are available in
    [pwt8](../packages/pwt8/index.html) and
    [pwt9](../packages/pwt9/index.html), respectively.
-   *Time series and forecasting data* : The packages
    [expsmooth](../packages/expsmooth/index.html),
    [fma](../packages/fma/index.html), and
    [Mcomp](../packages/Mcomp/index.html) are data packages with time
    series data from the books \'Forecasting with Exponential Smoothing:
    The State Space Approach\' (Hyndman, Koehler, Ord, Snyder, 2008,
    Springer) and \'Forecasting: Methods and Applications\' (Makridakis,
    Wheelwright, Hyndman, 3rd ed., 1998, Wiley) and the M-competitions,
    respectively.
-   *Empirical Research in Economics* : Package
    [erer](../packages/erer/index.html) contains functions and datasets
    for the book of \'Empirical Research in Economics: Growing up with
    R\' (Sun, forthcoming).
-   *Panel Study of Income Dynamics (PSID)* :
    [psidR](../packages/psidR/index.html) can build panel data sets from
    the Panel Study of Income Dynamics (PSID).
-   US state- and county-level panel data:
    [rUnemploymentData](../packages/rUnemploymentData/index.html).
-   World Bank data and statistics: The
    [wbstats](../packages/wbstats/index.html) package provides
    programmatic access to the World Bank API.

**Miscellaneous**

-   *Matrix manipulations* : As a vector- and matrix-based language,
    base R ships with many powerful tools for doing matrix
    manipulations, which are complemented by the packages
    [Matrix](../packages/Matrix/index.html) and
    [SparseM](../packages/SparseM/index.html).
-   *Optimization and mathematical programming* : R and many of its
    contributed packages provide many specialized functions for solving
    particular optimization problems, e.g., in regression as discussed
    above. Further functionality for solving more general optimization
    problems, e.g., likelihood maximization, is discussed in the the
    [Optimization](Optimization.html) task view.
-   *Bootstrap* : In addition to the recommended
    [boot](../packages/boot/index.html) package, there are some other
    general bootstrapping techniques available in
    [bootstrap](../packages/bootstrap/index.html) or
    [simpleboot](../packages/simpleboot/index.html) as well some
    bootstrap techniques designed for time-series data, such as the
    maximum entropy bootstrap in [meboot](../packages/meboot/index.html)
    or the `tsbootstrap()` from
    [tseries](../packages/tseries/index.html).
-   *Inequality* : For measuring inequality, concentration and poverty
    the package [ineq](../packages/ineq/index.html) provides some basic
    tools such as Lorenz curves, Pen\'s parade, the Gini coefficient and
    many more.
-   *Structural change* : R is particularly strong when dealing with
    structural changes and changepoints in parametric models, see
    [strucchange](../packages/strucchange/index.html) and
    [segmented](../packages/segmented/index.html).
-   *Exchange rate regimes* : Methods for inference about exchange rate
    regimes, in particular in a structural change setting, are provided
    by [fxregime](../packages/fxregime/index.html).
-   *Global value chains* : Tools and decompositions for global value
    chains are in [gvc](../packages/gvc/index.html) and
    [decompr](../packages/decompr/index.html).
-   *Regression discontinuity design* : A variety of methods are
    provided in the [rdd](../packages/rdd/index.html),
    [rddapp](../packages/rddapp/index.html),
    [rddtools](../packages/rddtools/index.html),
    [rdrobust](../packages/rdrobust/index.html), and
    [rdlocrand](../packages/rdlocrand/index.html) packages. The
    [rdpower](../packages/rdpower/index.html) package offers power
    calculations for regression discontinuity designs. And
    [rdmulti](../packages/rdmulti/index.html) implements analysis with
    multiple cutoffs or scores.
-   *z-Tree* : [zTree](../packages/zTree/index.html) can import data
    from the z-Tree software for developing and carrying out economic
    experiments.

</div>

### CRAN packages:

-   [AER](../packages/AER/index.html) (core)
-   [aod](../packages/aod/index.html)
-   [apt](../packages/apt/index.html)
-   [bayesm](../packages/bayesm/index.html)
-   [betareg](../packages/betareg/index.html)
-   [BMA](../packages/BMA/index.html)
-   [BMS](../packages/BMS/index.html)
-   [boot](../packages/boot/index.html)
-   [bootstrap](../packages/bootstrap/index.html)
-   [brglm](../packages/brglm/index.html)
-   [CADFtest](../packages/CADFtest/index.html)
-   [car](../packages/car/index.html) (core)
-   [CDNmoney](../packages/CDNmoney/index.html)
-   [censReg](../packages/censReg/index.html)
-   [clubSandwich](../packages/clubSandwich/index.html)
-   [clusterSEs](../packages/clusterSEs/index.html)
-   [crch](../packages/crch/index.html)
-   [decompr](../packages/decompr/index.html)
-   [dlsem](../packages/dlsem/index.html)
-   [dynlm](../packages/dynlm/index.html)
-   [Ecdat](../packages/Ecdat/index.html)
-   [effects](../packages/effects/index.html)
-   [erer](../packages/erer/index.html)
-   [expsmooth](../packages/expsmooth/index.html)
-   [ExtremeBounds](../packages/ExtremeBounds/index.html)
-   [fma](../packages/fma/index.html)
-   [forecast](../packages/forecast/index.html) (core)
-   [frm](../packages/frm/index.html)
-   [frontier](../packages/frontier/index.html)
-   [fxregime](../packages/fxregime/index.html)
-   [gam](../packages/gam/index.html)
-   [gamlss](../packages/gamlss/index.html)
-   [geepack](../packages/geepack/index.html)
-   [gets](../packages/gets/index.html)
-   [glmx](../packages/glmx/index.html)
-   [gmm](../packages/gmm/index.html)
-   [gmnl](../packages/gmnl/index.html)
-   [gvc](../packages/gvc/index.html)
-   [Hmisc](../packages/Hmisc/index.html)
-   [ineq](../packages/ineq/index.html)
-   [intReg](../packages/intReg/index.html)
-   [ivbma](../packages/ivbma/index.html)
-   [ivfixed](../packages/ivfixed/index.html)
-   [ivlewbel](../packages/ivlewbel/index.html)
-   [ivpack](../packages/ivpack/index.html)
-   [ivpanel](../packages/ivpanel/index.html)
-   [ivprobit](../packages/ivprobit/index.html)
-   [LARF](../packages/LARF/index.html)
-   [lavaan](../packages/lavaan/index.html)
-   [lfe](../packages/lfe/index.html)
-   [LinRegInteractive](../packages/LinRegInteractive/index.html)
-   [lme4](../packages/lme4/index.html)
-   [lmtest](../packages/lmtest/index.html) (core)
-   [margins](../packages/margins/index.html)
-   [MASS](../packages/MASS/index.html)
-   [matchingMarkets](../packages/matchingMarkets/index.html)
-   [Matrix](../packages/Matrix/index.html)
-   [Mcomp](../packages/Mcomp/index.html)
-   [meboot](../packages/meboot/index.html)
-   [mfx](../packages/mfx/index.html)
-   [mgcv](../packages/mgcv/index.html)
-   [mhurdle](../packages/mhurdle/index.html)
-   [micEcon](../packages/micEcon/index.html)
-   [micEconAids](../packages/micEconAids/index.html)
-   [micEconCES](../packages/micEconCES/index.html)
-   [micEconSNQP](../packages/micEconSNQP/index.html)
-   [midasr](../packages/midasr/index.html)
-   [mlogit](../packages/mlogit/index.html)
-   [mnlogit](../packages/mnlogit/index.html)
-   [MNP](../packages/MNP/index.html)
-   [MSBVAR](../packages/MSBVAR/index.html)
-   [multiwayvcov](../packages/multiwayvcov/index.html)
-   [mvProbit](../packages/mvProbit/index.html)
-   [nlme](../packages/nlme/index.html)
-   [nnet](../packages/nnet/index.html)
-   [nonnest2](../packages/nonnest2/index.html)
-   [np](../packages/np/index.html)
-   [ordinal](../packages/ordinal/index.html)
-   [OrthoPanels](../packages/OrthoPanels/index.html)
-   [pampe](../packages/pampe/index.html)
-   [panelAR](../packages/panelAR/index.html)
-   [Paneldata](../packages/Paneldata/index.html)
-   [panelvar](../packages/panelvar/index.html)
-   [PANICr](../packages/PANICr/index.html)
-   [pco](../packages/pco/index.html)
-   [pcse](../packages/pcse/index.html)
-   [pder](../packages/pder/index.html)
-   [pdR](../packages/pdR/index.html)
-   [pglm](../packages/pglm/index.html)
-   [phtt](../packages/phtt/index.html)
-   [plm](../packages/plm/index.html) (core)
-   [pscl](../packages/pscl/index.html)
-   [psidR](../packages/psidR/index.html)
-   [pwt](../packages/pwt/index.html)
-   [pwt8](../packages/pwt8/index.html)
-   [pwt9](../packages/pwt9/index.html)
-   [quantreg](../packages/quantreg/index.html)
-   [Rchoice](../packages/Rchoice/index.html)
-   [rdd](../packages/rdd/index.html)
-   [rddapp](../packages/rddapp/index.html)
-   [rddtools](../packages/rddtools/index.html)
-   [rdlocrand](../packages/rdlocrand/index.html)
-   [rdmulti](../packages/rdmulti/index.html)
-   [rdpower](../packages/rdpower/index.html)
-   [rdrobust](../packages/rdrobust/index.html)
-   [reldist](../packages/reldist/index.html)
-   [REndo](../packages/REndo/index.html)
-   [rms](../packages/rms/index.html)
-   [RSGHB](../packages/RSGHB/index.html)
-   [rUnemploymentData](../packages/rUnemploymentData/index.html)
-   [sampleSelection](../packages/sampleSelection/index.html)
-   [sandwich](../packages/sandwich/index.html) (core)
-   [segmented](../packages/segmented/index.html)
-   [sem](../packages/sem/index.html)
-   [SemiParSampleSel](../packages/SemiParSampleSel/index.html)
-   [semsfa](../packages/semsfa/index.html)
-   [sfa](../packages/sfa/index.html)
-   [simpleboot](../packages/simpleboot/index.html)
-   [SparseM](../packages/SparseM/index.html)
-   [spatialprobit](../packages/spatialprobit/index.html)
-   [spdep](../packages/spdep/index.html)
-   [spfrontier](../packages/spfrontier/index.html)
-   [sphet](../packages/sphet/index.html)
-   [splm](../packages/splm/index.html)
-   [ssfa](../packages/ssfa/index.html)
-   [strucchange](../packages/strucchange/index.html)
-   [survival](../packages/survival/index.html)
-   [systemfit](../packages/systemfit/index.html)
-   [truncreg](../packages/truncreg/index.html)
-   [tsDyn](../packages/tsDyn/index.html)
-   [tseries](../packages/tseries/index.html) (core)
-   [tsfa](../packages/tsfa/index.html)
-   [urca](../packages/urca/index.html) (core)
-   [vars](../packages/vars/index.html)
-   [VGAM](../packages/VGAM/index.html)
-   [wahc](../packages/wahc/index.html)
-   [wbstats](../packages/wbstats/index.html)
-   [wooldridge](../packages/wooldridge/index.html)
-   [xts](../packages/xts/index.html)
-   [Zelig](../packages/Zelig/index.html)
-   [zoo](../packages/zoo/index.html) (core)
-   [zTree](../packages/zTree/index.html)

### Related links:

-   CRAN Task View: [Finance](Finance.html)
-   CRAN Task View: [Optimization](Optimization.html)
-   CRAN Task View: [Psychometrics](Psychometrics.html)
-   CRAN Task View: [SocialSciences](SocialSciences.html)
-   CRAN Task View: [Spatial](Spatial.html)
-   CRAN Task View: [TimeSeries](TimeSeries.html)
-   [Mailing List: R Special Interest Group
    Finance](https://stat.ethz.ch/mailman/listinfo/R-SIG-Finance/)
-   [Journal of Statistical Software: Special Volume on \'Econometrics
    in R\' (2008)](http://www.jstatsoft.org/v27/)
-   [A Brief Guide to R for Beginners in
    Econometrics](https://mondo.su.se/access/content/user/ma@su.se/Public/)
-   [R for
    Economists](http://www.mayin.org/ajayshah/KB/R/R_for_economists.html)
-   [Book: Applied Econometrics with R (Kleiber and
    Zeileis)](http://www.springer.com/978-0-387-77316-2)
-   [Book: Using R for Introductory Econometrics
    (Heiss)](http://www.urfie.net/)
-   [Book: Hands-On Intermediate Econometrics Using R
    (Vinod)](http://www.worldscibooks.com/economics/6895.html)
