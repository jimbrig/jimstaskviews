## CRAN Task View: Bayesian Inference

  ----------------- --------------------------------------------
  **Maintainer:**   Jong Hee Park
  **Contact:**      jongheepark at snu.ac.kr
  **Version:**      2018-04-13
  **URL:**          <https://CRAN.R-project.org/view=Bayesian>
  ----------------- --------------------------------------------

<div>

Applied researchers interested in Bayesian statistics are increasingly
attracted to R because of the ease of which one can code algorithms to
sample from posterior distributions as well as the significant number of
packages contributed to the Comprehensive R Archive Network (CRAN) that
provide tools for Bayesian inference. This task view catalogs these
tools. In this task view, we divide those packages into four groups
based on the scope and focus of the packages. We first review R packages
that provide Bayesian estimation tools for a wide range of models. We
then discuss packages that address specific Bayesian models or
specialized methods in Bayesian statistics. This is followed by a
description of packages used for post-estimation analysis. Finally, we
review packages that link R to other Bayesian sampling engines such as
[JAGS](http://mcmc-jags.sourceforge.net/) ,
[OpenBUGS](http://www.openbugs.net/) ,
[WinBUGS](http://www.mrc-bsu.cam.ac.uk/software/bugs/) , and
[Stan](http://mc-stan.org/) .

**Bayesian packages for general model fitting**

-   The [arm](../packages/arm/index.html) package contains R functions
    for Bayesian inference using lm, glm, mer and polr objects.
-   [BACCO](../packages/BACCO/index.html) is an R bundle for Bayesian
    analysis of random functions. [BACCO](../packages/BACCO/index.html)
    contains three sub-packages: emulator, calibrator, and approximator,
    that perform Bayesian emulation and calibration of computer
    programs.
-   [bayesm](../packages/bayesm/index.html) provides R functions for
    Bayesian inference for various models widely used in marketing and
    micro-econometrics. The models include linear regression models,
    multinomial logit, multinomial probit, multivariate probit,
    multivariate mixture of normals (including clustering), density
    estimation using finite mixtures of normals as well as Dirichlet
    Process priors, hierarchical linear models, hierarchical multinomial
    logit, hierarchical negative binomial regression models, and linear
    instrumental variable models.
-   [bayesSurv](../packages/bayesSurv/index.html) contains R functions
    to perform Bayesian inference for survival regression models with
    flexible error and random effects distributions.
-   [DPpackage](../packages/DPpackage/index.html) contains R functions
    for Bayesian nonparametric and semiparametric models. DPpackage
    currently includes semiparametric models for density estimation, ROC
    curve analysis, interval censored data, binary regression models,
    generalized linear mixed models, and IRT type models.
-   [LaplacesDemon](../packages/LaplacesDemon/index.html) seeks to
    provide a complete Bayesian environment, including numerous MCMC
    algorithms, Laplace Approximation with multiple optimization
    algorithms, scores of examples, dozens of additional probability
    distributions, numerous MCMC diagnostics, Bayes factors, posterior
    predictive checks, a variety of plots, elicitation, parameter and
    variable importance, and numerous additional utility functions.
-   [MCMCpack](../packages/MCMCpack/index.html) provides model-specific
    Markov chain Monte Carlo (MCMC) algorithms for wide range of models
    commonly used in the social and behavioral sciences. It contains R
    functions to fit a number of regression models (linear regression,
    logit, ordinal probit, probit, Poisson regression, etc.),
    measurement models (item response theory and factor models),
    changepoint models (linear regression, binary probit, ordinal
    probit, Poisson, panel), and models for ecological inference. It
    also contains a generic Metropolis sampler that can be used to fit
    arbitrary models.
-   The [mcmc](../packages/mcmc/index.html) package consists of an R
    function for a random-walk Metropolis algorithm for a continuous
    random vector.
-   The [nimble](../packages/nimble/index.html) package provides a
    general MCMC system that allows customizable MCMC for models written
    in the BUGS/JAGS model language. Users can choose samplers and write
    new samplers. Models and samplers are automatically compiled via
    generated C++. The package also supports other methods such as
    particle filtering or whatever users write in its algorithm
    language.

**Bayesian packages for specific models or methods**

-   [abc](../packages/abc/index.html) package implements several ABC
    algorithms for performing parameter estimation and model selection.
    Cross-validation tools are also available for measuring the accuracy
    of ABC estimates, and to calculate the misclassification
    probabilities of different models.
-   [abn](../packages/abn/index.html) is a package for modelling
    multivariate data using additive Bayesian networks. It provides
    routines to help determine optimal Bayesian network models for a
    given data set, where these models are used to identify statistical
    dependencies in messy, complex data.
-   [AdMit](../packages/AdMit/index.html) provides functions to perform
    the fitting of an adapative mixture of Student-t distributions to a
    target density through its kernel function. The mixture
    approximation can be used as the importance density in importance
    sampling or as the candidate density in the Metropolis-Hastings
    algorithm.
-   The [BaBooN](../packages/BaBooN/index.html) package contains two
    variants of Bayesian Bootstrap Predictive Mean Matching to multiply
    impute missing data.
-   [bamlss](../packages/bamlss/index.html) provides an infrastructure
    for estimating probabilistic distributional regression models in a
    Bayesian framework. The distribution parameters may capture
    location, scale, shape, etc. and every parameter may depend on
    complex additive terms similar to a generalized additive model.
-   The [BAS](../packages/BAS/index.html) package implements BMA for
    regression models using g-priors and mixtures of g-priors.
    [BAS](../packages/BAS/index.html) utilizes an efficient algorithm to
    sample models without replacement.
-   The [bayesGARCH](../packages/bayesGARCH/index.html) package provides
    a function which perform the Bayesian estimation of the GARCH(1,1)
    model with Student\'s t innovations.
-   [bayesImageS](../packages/bayesImageS/index.html) is an R package
    for Bayesian image analysis using the hidden Potts model.
-   [bayesmeta](../packages/bayesmeta/index.html) is an R package to
    perform meta-analyses within the common random-effects model
    framework.
-   [BayesSummaryStatLM](../packages/BayesSummaryStatLM/index.html)
    provides two functions: one function that computes summary
    statistics of data and one function that carries out the MCMC
    posterior sampling for Bayesian linear regression models where
    summary statistics are used as input.
-   [Bayesthresh](../packages/Bayesthresh/index.html) fits a linear
    mixed model for ordinal categorical responses using Bayesian
    inference via Monte Carlo Markov Chains. Default is Nandran and Chen
    algorithm using Gaussian link function and saving just the summaries
    of the chains.
-   [BayesTree](../packages/BayesTree/index.html) implements BART
    (Bayesian Additive Regression Trees) by Chipman, George, and
    McCulloch (2006).
-   [bayesQR](../packages/bayesQR/index.html) supports Bayesian quantile
    regression using the asymmetric Laplace distribution, both
    continuous as well as binary dependent variables.
-   [BayesFactor](../packages/BayesFactor/index.html) provides a suite
    of functions for computing various Bayes factors for simple designs,
    including contingency tables, one- and two-sample designs, one-way
    designs, general ANOVA designs, and linear regression.
-   [BayesVarSel](../packages/BayesVarSel/index.html) calculate Bayes
    factors in linear models and then to provide a formal Bayesian
    answer to testing and variable selection problems.
-   [BayHaz](../packages/BayHaz/index.html) contains a suite of R
    functions for Bayesian estimation of smooth hazard rates via
    Compound Poisson Process (CPP) priors.
-   [BAYSTAR](../packages/BAYSTAR/index.html) provides functions for
    Bayesian estimation of threshold autoregressive models.
-   [bbemkr](../packages/bbemkr/index.html) implements Bayesian
    bandwidth estimation for Nadaraya-Watson type multivariate kernel
    regression with Gaussian error.
-   [BCE](../packages/BCE/index.html) contains function to estimates
    taxonomic compositions from biomarker data using a Bayesian
    approach.
-   [BCBCSF](../packages/BCBCSF/index.html) provides functions to
    predict the discrete response based on selected high dimensional
    features, such as gene expression data.
-   [bclust](../packages/bclust/index.html) builds a dendrogram with log
    posterior as a natural distance defined by the model. It is also
    capable to computing Bayesian discrimination probabilities
    equivalent to the implemented Bayesian clustering. Spike-and-Slab
    models are adopted in a way to be able to produce an importance
    measure for clustering and discriminant variables.
-   [bcp](../packages/bcp/index.html) implements a Bayesian analysis of
    changepoint problem using Barry and Hartigan product partition
    model.
-   [bisoreg](../packages/bisoreg/index.html) implements the Bayesian
    isotonic regression with Bernstein polynomials.
-   [BLR](../packages/BLR/index.html) provides R functions to fit
    parametric regression models using different types of shrinkage
    methods.
-   The [BMA](../packages/BMA/index.html) package has functions for
    Bayesian model averaging for linear models, generalized linear
    models, and survival models. The complementary package
    [ensembleBMA](../packages/ensembleBMA/index.html) uses the
    [BMA](../packages/BMA/index.html) package to create probabilistic
    forecasts of ensembles using a mixture of normal distributions.
-   [BMS](../packages/BMS/index.html) is Bayesian Model Averaging
    library for linear models with a wide choice of (customizable)
    priors. Built-in priors include coefficient priors (fixed, flexible
    and hyper-g priors), and 5 kinds of model priors.
-   [Bmix](../packages/Bmix/index.html) is a bare-bones implementation
    of sampling algorithms for a variety of Bayesian stick-breaking
    (marginally DP) mixture models, including particle learning and
    Gibbs sampling for static DP mixtures, particle learning for dynamic
    BAR stick-breaking, and DP mixture regression.
-   [bnlearn](../packages/bnlearn/index.html) is a package for Bayesian
    network structure learning (via constraint-based, score-based and
    hybrid algorithms), parameter learning (via ML and Bayesian
    estimators) and inference.
-   [BoomSpikeSlab](../packages/BoomSpikeSlab/index.html) provides
    functions to do spike and slab regression via the stochastic search
    variable selection algorithm. It handles probit, logit, poisson, and
    student T data.
-   [bqtl](../packages/bqtl/index.html) can be used to fit quantitative
    trait loci (QTL) models. This package allows Bayesian estimation of
    multi-gene models via Laplace approximations and provides tools for
    interval mapping of genetic loci. The package also contains
    graphical tools for QTL analysis.
-   [bridgesampling](../packages/bridgesampling/index.html) provides R
    functions for estimating marginal likelihoods, Bayes factors,
    posterior model probabilities, and normalizing constants in general,
    via different versions of bridge sampling (Meng and Wong, 1996).
-   [bsamGP](../packages/bsamGP/index.html) provides functions to
    perform Bayesian inference using a spectral analysis of Gaussian
    process priors. Gaussian processes are represented with a Fourier
    series based on cosine basis functions. Currently the package
    includes parametric linear models, partial linear additive models
    with/without shape restrictions, generalized linear additive models
    with/without shape restrictions, and density estimation model.
-   [bspec](../packages/bspec/index.html) performs Bayesian inference on
    the (discrete) power spectrum of time series.
-   [bspmma](../packages/bspmma/index.html) is a package for Bayesian
    semiparametric models for meta-analysis.
-   [BSquare](../packages/BSquare/index.html) models the quantile
    process as a function of predictors.
-   [bsts](../packages/bsts/index.html) is a package for time series
    regression using dynamic linear models using MCMC.
-   [BVS](../packages/BVS/index.html) is a package for Bayesian variant
    selection and Bayesian model uncertainty techniques for genetic
    association studies.
-   [catnet](../packages/catnet/index.html) is a package that handles
    discrete Bayesian network models and provides inference using the
    frequentist approach.
-   [coalescentMCMC](../packages/coalescentMCMC/index.html) provides a
    flexible framework for coalescent analyses in R.
-   [dclone](../packages/dclone/index.html) provides low level functions
    for implementing maximum likelihood estimating procedures for
    complex models using data cloning and MCMC methods.
-   [deal](../packages/deal/index.html) provides R functions for
    Bayesian network analysis; the current version of covers discrete
    and continuous variables under Gaussian network structure.
-   [deBInfer](../packages/deBInfer/index.html) provides R functions for
    Bayesian parameter inference in differential equations using MCMC
    methods.
-   [dlm](../packages/dlm/index.html) is a package for Bayesian (and
    likelihood) analysis of dynamic linear models. It includes the
    calculations of the Kalman filter and smoother, and the forward
    filtering backward sampling algorithm.
-   [EbayesThresh](../packages/EbayesThresh/index.html) implements
    Bayesian estimation for thresholding methods. Although the original
    model is developed in the context of wavelets, this package is
    useful when researchers need to take advantage of possible sparsity
    in a parameter set.
-   [eco](../packages/eco/index.html) fits Bayesian ecological inference
    models in two by two tables using MCMC methods.
-   [ebdbNet](../packages/ebdbNet/index.html) can be used to infer the
    adjacency matrix of a network from time course data using an
    empirical Bayes estimation procedure based on Dynamic Bayesian
    Networks.
-   [eigenmodel](../packages/eigenmodel/index.html) estimates the
    parameters of a model for symmetric relational data (e.g., the
    above-diagonal part of a square matrix), using a model-based
    eigenvalue decomposition and regression using MCMC methods.
-   [evdbayes](../packages/evdbayes/index.html) provides tools for
    Bayesian analysis of extreme value models.
-   [exactLoglinTest](../packages/exactLoglinTest/index.html) provides
    functions for log-linear models that compute Monte Carlo estimates
    of conditional P-values for goodness of fit tests.
-   [factorQR](../packages/factorQR/index.html) is a package to fit
    Bayesian quantile regression models that assume a factor structure
    for at least part of the design matrix.
-   [FME](../packages/FME/index.html) provides functions to help in
    fitting models to data, to perform Monte Carlo, sensitivity and
    identifiability analysis. It is intended to work with models be
    written as a set of differential equations that are solved either by
    an integration routine from deSolve, or a steady-state solver from
    rootSolve.
-   The `gbayes()` function in [Hmisc](../packages/Hmisc/index.html)
    derives the posterior (and optionally) the predictive distribution
    when both the prior and the likelihood are Gaussian, and when the
    statistic of interest comes from a two-sample problem.
-   [ggmcmc](../packages/ggmcmc/index.html) is a tool for assessing and
    diagnosing convergence of Markov Chain Monte Carlo simulations, as
    well as for graphically display results from full MCMC analysis.
-   [gRain](../packages/gRain/index.html) is a package for probability
    propagation in graphical independence networks, also known as
    Bayesian networks or probabilistic expert systems.
-   [growcurves](../packages/growcurves/index.html) is a package for
    Bayesian semi and nonparametric growth curve models that
    additionally include multiple membership random effects.
-   The [HI](../packages/HI/index.html) package has functions to
    implement a geometric approach to transdimensional MCMC methods and
    random direction multivariate Adaptive Rejection Metropolis
    Sampling.
-   The [hbsae](../packages/hbsae/index.html) package provides functions
    to compute small area estimates based on a basic area or unit-level
    model. The model is fit using restricted maximum likelihood, or in a
    hierarchical Bayesian way.
-   [iterLap](../packages/iterLap/index.html) performs an iterative
    Laplace approximation to build a global approximation of the
    posterior (using mixture distributions) and then uses importance
    sampling for simulation based inference.
-   The function `krige.bayes()` in the
    [geoR](../packages/geoR/index.html) package performs Bayesian
    analysis of geostatistical data allowing specification of different
    levels of uncertainty in the model parameters. The
    `binom.krige.bayes()` function in the
    [geoRglm](../packages/geoRglm/index.html) package implements
    Bayesian posterior simulation and spatial prediction for the
    binomial spatial model (see the [Spatial](Spatial.html) view for
    more information).
-   The [lmm](../packages/lmm/index.html) package contains R functions
    to fit linear mixed models using MCMC methods.
-   [MasterBayes](../packages/MasterBayes/index.html) is an R package
    that implements MCMC methods to integrate over uncertainty in
    pedigree configurations estimated from molecular markers and
    phenotypic data.
-   [matchingMarkets](../packages/matchingMarkets/index.html) implements
    a structural model based on a Gibbs sampler to correct for the bias
    from endogenous matching (e.g. group formation or two-sided
    matching).
-   [MCMCglmm](../packages/MCMCglmm/index.html) is package for fitting
    Generalised Linear Mixed Models using MCMC methods.
-   The `mcmcsamp()` function in [lme4](../packages/lme4/index.html)
    allows MCMC sampling for the linear mixed model and generalized
    linear mixed model.
-   The [mlogitBMA](../packages/mlogitBMA/index.html) Provides a
    modified function `bic.glm()` of the
    [BMA](../packages/BMA/index.html) package that can be applied to
    multinomial logit (MNL) data.
-   The [MNP](../packages/MNP/index.html) package fits multinomial
    probit models using MCMC methods.
-   [mombf](../packages/mombf/index.html) performs model selection based
    on non-local priors, including MOM, eMOM and iMOM priors..
-   [monomvn](../packages/monomvn/index.html) is an R package for
    estimation of multivariate normal and Student-t data of arbitrary
    dimension where the pattern of missing data is monotone.
-   [MSBVAR](../packages/MSBVAR/index.html) is an R package for
    estimating Bayesian Vector Autoregression models and Bayesian
    structural Vector Autoregression models.
-   [NetworkChange](../packages/NetworkChange/index.html) is an R
    package for change point analysis in longitudinal network data. It
    implements a hidden Markovmultilinear tensor regression model. Model
    diagnostic tools using marginal likelihoods and WAIC are provided.
-   [openEBGM](../packages/openEBGM/index.html) calculates Empirical
    Bayes Geometric Mean (EBGM) and quantile scores from the posterior
    distribution using the Gamma-Poisson Shrinker (GPS) model to find
    unusually large cell counts in large, sparse contingency tables.
-   [pacbpred](../packages/pacbpred/index.html) perform estimation and
    prediction in high-dimensional additive models, using a sparse
    PAC-Bayesian point of view and a MCMC algorithm.
-   [predmixcor](../packages/predmixcor/index.html) provides functions
    to predict the binary response based on high dimensional binary
    features modeled with Bayesian mixture models.
-   [prevalence](../packages/prevalence/index.html) provides functions
    for the estimation of true prevalence from apparent prevalence in a
    Bayesian framework. MCMC sampling is performed via JAGS/rjags.
-   [profdpm](../packages/profdpm/index.html) facilitates profile
    inference (inference at the posterior mode) for a class of product
    partition models.
-   The [pscl](../packages/pscl/index.html) package provides R functions
    to fit item-response theory models using MCMC methods and to compute
    highest density regions for the Beta distribution and the inverse
    gamma distribution.
-   The [PAWL](../packages/PAWL/index.html) package implements parallel
    adaptive Metropolis-Hastings and sequential Monte Carlo samplers for
    sampling from multimodal target distributions.
-   [PReMiuM](../packages/PReMiuM/index.html) is a package for profile
    regression, which is a Dirichlet process Bayesian clustering where
    the response is linked non-parametrically to the covariate profile.
-   [revdbayes](../packages/revdbayes/index.html) provides functions for
    the Bayesian analysis of extreme value models using direct random
    sampling from extreme value posterior distributions.
-   The [RJaCGH](../packages/RJaCGH/index.html) package implements
    Bayesian analysis of CGH microarrays using hidden Markov chain
    models. The selection of the number of states is made via their
    posterior probability computed by reversible jump Markov chain Monte
    Carlo Methods.
-   The `hitro.new()` function in
    [Runuran](../packages/Runuran/index.html) provides an MCMC sampler
    based on the Hit-and-Run algorithm in combination with the
    Ratio-of-Uniforms method.
-   [RSGHB](../packages/RSGHB/index.html) can be used to estimate models
    using a hierarchical Bayesian framework and provides flexibility in
    allowing the user to specify the likelihood function directly
    instead of assuming predetermined model structures.
-   [rstiefel](../packages/rstiefel/index.html) simulates random
    orthonormal matrices from linear and quadratic exponential family
    distributions on the Stiefel manifold using the Gibbs sampling
    method. The most general type of distribution covered is the
    matrix-variate Bingham-von Mises-Fisher distribution.
-   [RxCEcolInf](../packages/RxCEcolInf/index.html) fits the R x C
    inference model described in Greiner and Quinn (2009).
-   [SamplerCompare](../packages/SamplerCompare/index.html) provides a
    framework for running sets of MCMC samplers on sets of distributions
    with a variety of tuning parameters, along with plotting functions
    to visualize the results of those simulations.
-   [SampleSizeMeans](../packages/SampleSizeMeans/index.html) contains a
    set of R functions for calculating sample size requirements using
    three different Bayesian criteria in the context of designing an
    experiment to estimate a normal mean or the difference between two
    normal means.
-   [SampleSizeProportions](../packages/SampleSizeProportions/index.html)
    contains a set of R functions for calculating sample size
    requirements using three different Bayesian criteria in the context
    of designing an experiment to estimate the difference between two
    binomial proportions.
-   [sbgcop](../packages/sbgcop/index.html) estimates parameters of a
    Gaussian copula, treating the univariate marginal distributions as
    nuisance parameters as described in Hoff(2007). It also provides a
    semiparametric imputation procedure for missing multivariate data.
-   [SimpleTable](../packages/SimpleTable/index.html) provides a series
    of methods to conduct Bayesian inference and sensitivity analysis
    for causal effects from 2 x 2 and 2 x 2 x K tables.
-   [sna](../packages/sna/index.html), an R package for social network
    analysis, contains functions to generate posterior samples from
    Butt\'s Bayesian network accuracy model using Gibbs sampling.
-   [spBayes](../packages/spBayes/index.html) provides R functions that
    fit Gaussian spatial process models for univariate as well as
    multivariate point-referenced data using MCMC methods.
-   [spikeslab](../packages/spikeslab/index.html) provides functions for
    prediction and variable selection using spike and slab regression.
-   [spikeSlabGAM](../packages/spikeSlabGAM/index.html) implements
    Bayesian variable selection, model choice, and regularized
    estimation in (geo-)additive mixed models for Gaussian, binomial,
    and Poisson responses.
-   [spTimer](../packages/spTimer/index.html) fits, spatially predict
    and temporally forecast large amounts of space-time data using
    Bayesian Gaussian Process Models, Bayesian Auto-Regressive (AR)
    Models, and Bayesian Gaussian Predictive Processes based AR Models.
-   [stochvol](../packages/stochvol/index.html) provides efficient
    algorithms for fully Bayesian estimation of stochastic volatility
    (SV) models.
-   The [tgp](../packages/tgp/index.html) package implements Bayesian
    treed Gaussian process models: a spatial modeling and regression
    package providing fully Bayesian MCMC posterior inference for models
    ranging from the simple linear model, to nonstationary treed
    Gaussian process, and others in between.
-   [tRophicPosition](../packages/tRophicPosition/index.html) estimates
    the trophic position of a consumer relative to a baseline species.
    It implements a Bayesian approach which combines an interface to the
    JAGS MCMC library of rjags and stable isotopes.
-   [[vbmp]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/vbmp.html)
    is a package for variational Bayesian multinomial probit regression
    with Gaussian process priors. It estimates class membership
    posterior probability employing variational and sparse approximation
    to the full posterior. This software also incorporates feature
    weighting by means of Automatic Relevance Determination.
-   The `vcov.gam()` function the [mgcv](../packages/mgcv/index.html)
    package can extract a Bayesian posterior covariance matrix of the
    parameters from a fitted `gam` object.
-   [zic](../packages/zic/index.html) provides functions for an MCMC
    analysis of zero-inflated count models including stochastic search
    variable selection.

**Post-estimation tools**

-   [BayesValidate](../packages/BayesValidate/index.html) implements a
    software validation method for Bayesian softwares.
-   The [boa](../packages/boa/index.html) package provides functions for
    diagnostics, summarization, and visualization of MCMC sequences. It
    imports draws from BUGS format, or from plain matrices.
    [boa](../packages/boa/index.html) provides the Gelman and Rubin,
    Geweke, Heidelberger and Welch, and Raftery and Lewis diagnostics,
    the Brooks and Gelman multivariate shrink factors.
-   The [coda](../packages/coda/index.html) (Convergence Diagnosis and
    Output Analysis) package is a suite of functions that can be used to
    summarize, plot, and and diagnose convergence from MCMC samples.
    [coda](../packages/coda/index.html) also defines an `mcmc` object
    and related methods which are used by other packages. It can easily
    import MCMC output from WinBUGS, OpenBUGS, and JAGS, or from plain
    matrices. [coda](../packages/coda/index.html) contains the Gelman
    and Rubin, Geweke, Heidelberger and Welch, and Raftery and Lewis
    diagnostics.
-   [ramps](../packages/ramps/index.html) implements Bayesian
    geostatistical analysis of Gaussian processes using a
    reparameterized and marginalized posterior sampling algorithm.

**Packages for learning Bayesian statistics**

-   [AtelieR](../packages/AtelieR/index.html) is a GTK interface for
    teaching basic concepts in statistical inference, and doing
    elementary bayesian statistics (inference on proportions,
    multinomial counts, means and variances).
-   The [BaM](../packages/BaM/index.html) package is an R package
    associated with Jeff Gill\'s book, \"Bayesian Methods: A Social and
    Behavioral Sciences Approach, Second Edition\" (CRC Press, 2007).
-   [BayesDA](../packages/BayesDA/index.html) provides R functions and
    datasets for \"Bayesian Data Analysis, Second Edition\" (CRC Press,
    2003) by Andrew Gelman, John B. Carlin, Hal S. Stern, and Donald B.
    Rubin.
-   The [Bolstad](../packages/Bolstad/index.html) package contains a set
    of R functions and data sets for the book Introduction to Bayesian
    Statistics, by Bolstad, W.M. (2007).
-   The [LearnBayes](../packages/LearnBayes/index.html) package contains
    a collection of functions helpful in learning the basic tenets of
    Bayesian statistical inference. It contains functions for
    summarizing basic one and two parameter posterior distributions and
    predictive distributions and MCMC algorithms for summarizing
    posterior distributions defined by the user. It also contains
    functions for regression models, hierarchical models, Bayesian
    tests, and illustrations of Gibbs sampling.

**Packages that link R to other sampling engines**

-   [bayesmix](../packages/bayesmix/index.html) is an R package to fit
    Bayesian mixture models using
    [JAGS](http://mcmc-jags.sourceforge.net/) .
-   [BayesX](../packages/BayesX/index.html) provides functionality for
    exploring and visualizing estimation results obtained with the
    software package [BayesX](http://www.BayesX.org/) .
-   [Boom](../packages/Boom/index.html) provides a C++ library for
    Bayesian modeling, with an emphasis on Markov chain Monte Carlo.
-   **BRugs** provides an R interface to
    [OpenBUGS](http://www.openbugs.net/) . It works under Windows and
    Linux. **BRugs** used to be available from CRAN, now it is located
    at the [CRANextras](http://www.stats.ox.ac.uk/pub/RWin/) repository.
-   [brms](../packages/brms/index.html) implements Bayesian multilevel
    models in R using [Stan](http://mc-stan.org/) . A wide range of
    distributions and link functions are supported, allowing users to
    fit linear, robust linear, binomial, Pois- son, survival, response
    times, ordinal, quantile, zero-inflated, hurdle, and even non-linear
    models all in a multilevel context.
-   [cudaBayesreg](../packages/cudaBayesreg/index.html) provides a
    Compute Unified Device Architecture (CUDA) implementation of a
    Bayesian multilevel model for the analysis of brain fMRI data. CUDA
    is a software platform for massively parallel high-performance
    computing on NVIDIA GPUs.
-   There are two packages that can be used to interface R with
    [WinBUGS](http://www.mrc-bsu.cam.ac.uk/software/bugs/) .
    [R2WinBUGS](../packages/R2WinBUGS/index.html) provides a set of
    functions to call WinBUGS on a Windows system and a Linux system;
    [rbugs](../packages/rbugs/index.html) supports Linux systems through
    [OpenBUGS](http://www.openbugs.net/) on Linux (LinBUGS).
-   [glmmBUGS](../packages/glmmBUGS/index.html) writes BUGS model files,
    formats data, and creates starting values for generalized linear
    mixed models.
-   There are three packages that provide R interface with [Just Another
    Gibbs Sampler (JAGS)](http://mcmc-jags.sourceforge.net/) :
    [rjags](../packages/rjags/index.html),
    [R2jags](../packages/R2jags/index.html), and
    [runjags](../packages/runjags/index.html).
-   All of these BUGS engines use graphical models for model
    specification. As such, the [gR](gR.html) task view may be of
    interest.
-   [rstan](../packages/rstan/index.html) provides R functions to parse,
    compile, test, estimate, and analyze Stan models by accessing the
    header-only Stan library provided by the \`StanHeaders\' package.
    The [Stan](http://mc-stan.org/) project develops a probabilistic
    programming language that implements full Bayesian statistical
    inference via MCMC and (optionally penalized) maximum likelihood
    estimation via optimization.
-   [R2BayesX](../packages/R2BayesX/index.html) provides an R interface
    to estimate structured additive regression (STAR) models with
    \'BayesX\'.

The Bayesian Inference Task View is written by Jong Hee Park (Seoul
National University, South Korea), Andrew D. Martin (University of
Michigan, Ann Arbor, MI, USA), and Kevin M. Quinn (UC Berkeley,
Berkeley, CA, USA). Please email the [task view
maintainer](mailto:jongheepark@snu.ac.kr?subject=Bayesian%20Task%20View)
with suggestions.

</div>

### CRAN packages:

-   [abc](../packages/abc/index.html)
-   [abn](../packages/abn/index.html)
-   [AdMit](../packages/AdMit/index.html)
-   [arm](../packages/arm/index.html) (core)
-   [AtelieR](../packages/AtelieR/index.html)
-   [BaBooN](../packages/BaBooN/index.html)
-   [BACCO](../packages/BACCO/index.html) (core)
-   [BaM](../packages/BaM/index.html)
-   [bamlss](../packages/bamlss/index.html)
-   [BAS](../packages/BAS/index.html)
-   [BayesDA](../packages/BayesDA/index.html)
-   [BayesFactor](../packages/BayesFactor/index.html)
-   [bayesGARCH](../packages/bayesGARCH/index.html)
-   [bayesImageS](../packages/bayesImageS/index.html)
-   [bayesm](../packages/bayesm/index.html) (core)
-   [bayesmeta](../packages/bayesmeta/index.html)
-   [bayesmix](../packages/bayesmix/index.html)
-   [bayesQR](../packages/bayesQR/index.html)
-   [BayesSummaryStatLM](../packages/BayesSummaryStatLM/index.html)
-   [bayesSurv](../packages/bayesSurv/index.html) (core)
-   [Bayesthresh](../packages/Bayesthresh/index.html)
-   [BayesTree](../packages/BayesTree/index.html)
-   [BayesValidate](../packages/BayesValidate/index.html)
-   [BayesVarSel](../packages/BayesVarSel/index.html)
-   [BayesX](../packages/BayesX/index.html)
-   [BayHaz](../packages/BayHaz/index.html)
-   [BAYSTAR](../packages/BAYSTAR/index.html)
-   [bbemkr](../packages/bbemkr/index.html)
-   [BCBCSF](../packages/BCBCSF/index.html)
-   [BCE](../packages/BCE/index.html)
-   [bclust](../packages/bclust/index.html)
-   [bcp](../packages/bcp/index.html)
-   [bisoreg](../packages/bisoreg/index.html)
-   [BLR](../packages/BLR/index.html)
-   [BMA](../packages/BMA/index.html)
-   [Bmix](../packages/Bmix/index.html)
-   [BMS](../packages/BMS/index.html)
-   [bnlearn](../packages/bnlearn/index.html)
-   [boa](../packages/boa/index.html) (core)
-   [Bolstad](../packages/Bolstad/index.html)
-   [Boom](../packages/Boom/index.html)
-   [BoomSpikeSlab](../packages/BoomSpikeSlab/index.html)
-   [bqtl](../packages/bqtl/index.html)
-   [bridgesampling](../packages/bridgesampling/index.html)
-   [brms](../packages/brms/index.html)
-   [bsamGP](../packages/bsamGP/index.html)
-   [bspec](../packages/bspec/index.html)
-   [bspmma](../packages/bspmma/index.html)
-   [BSquare](../packages/BSquare/index.html)
-   [bsts](../packages/bsts/index.html)
-   [BVS](../packages/BVS/index.html)
-   [catnet](../packages/catnet/index.html)
-   [coalescentMCMC](../packages/coalescentMCMC/index.html)
-   [coda](../packages/coda/index.html) (core)
-   [cudaBayesreg](../packages/cudaBayesreg/index.html)
-   [dclone](../packages/dclone/index.html)
-   [deal](../packages/deal/index.html)
-   [deBInfer](../packages/deBInfer/index.html)
-   [dlm](../packages/dlm/index.html)
-   [DPpackage](../packages/DPpackage/index.html) (core)
-   [EbayesThresh](../packages/EbayesThresh/index.html)
-   [ebdbNet](../packages/ebdbNet/index.html)
-   [eco](../packages/eco/index.html)
-   [eigenmodel](../packages/eigenmodel/index.html)
-   [ensembleBMA](../packages/ensembleBMA/index.html)
-   [evdbayes](../packages/evdbayes/index.html)
-   [exactLoglinTest](../packages/exactLoglinTest/index.html)
-   [factorQR](../packages/factorQR/index.html)
-   [FME](../packages/FME/index.html)
-   [geoR](../packages/geoR/index.html)
-   [geoRglm](../packages/geoRglm/index.html)
-   [ggmcmc](../packages/ggmcmc/index.html)
-   [glmmBUGS](../packages/glmmBUGS/index.html)
-   [gRain](../packages/gRain/index.html)
-   [growcurves](../packages/growcurves/index.html)
-   [hbsae](../packages/hbsae/index.html)
-   [HI](../packages/HI/index.html)
-   [Hmisc](../packages/Hmisc/index.html)
-   [iterLap](../packages/iterLap/index.html)
-   [LaplacesDemon](../packages/LaplacesDemon/index.html)
-   [LearnBayes](../packages/LearnBayes/index.html)
-   [lme4](../packages/lme4/index.html)
-   [lmm](../packages/lmm/index.html)
-   [MasterBayes](../packages/MasterBayes/index.html)
-   [matchingMarkets](../packages/matchingMarkets/index.html)
-   [mcmc](../packages/mcmc/index.html) (core)
-   [MCMCglmm](../packages/MCMCglmm/index.html)
-   [MCMCpack](../packages/MCMCpack/index.html) (core)
-   [mgcv](../packages/mgcv/index.html)
-   [mlogitBMA](../packages/mlogitBMA/index.html)
-   [MNP](../packages/MNP/index.html)
-   [mombf](../packages/mombf/index.html)
-   [monomvn](../packages/monomvn/index.html)
-   [MSBVAR](../packages/MSBVAR/index.html)
-   [NetworkChange](../packages/NetworkChange/index.html)
-   [nimble](../packages/nimble/index.html) (core)
-   [openEBGM](../packages/openEBGM/index.html)
-   [pacbpred](../packages/pacbpred/index.html)
-   [PAWL](../packages/PAWL/index.html)
-   [predmixcor](../packages/predmixcor/index.html)
-   [PReMiuM](../packages/PReMiuM/index.html)
-   [prevalence](../packages/prevalence/index.html)
-   [profdpm](../packages/profdpm/index.html)
-   [pscl](../packages/pscl/index.html)
-   [R2BayesX](../packages/R2BayesX/index.html)
-   [R2jags](../packages/R2jags/index.html)
-   [R2WinBUGS](../packages/R2WinBUGS/index.html)
-   [ramps](../packages/ramps/index.html)
-   [rbugs](../packages/rbugs/index.html)
-   [revdbayes](../packages/revdbayes/index.html)
-   [RJaCGH](../packages/RJaCGH/index.html)
-   [rjags](../packages/rjags/index.html)
-   [RSGHB](../packages/RSGHB/index.html)
-   [RSGHB](../packages/RSGHB/index.html)
-   [rstan](../packages/rstan/index.html)
-   [rstiefel](../packages/rstiefel/index.html)
-   [runjags](../packages/runjags/index.html)
-   [Runuran](../packages/Runuran/index.html)
-   [RxCEcolInf](../packages/RxCEcolInf/index.html)
-   [SamplerCompare](../packages/SamplerCompare/index.html)
-   [SampleSizeMeans](../packages/SampleSizeMeans/index.html)
-   [SampleSizeProportions](../packages/SampleSizeProportions/index.html)
-   [sbgcop](../packages/sbgcop/index.html)
-   [SimpleTable](../packages/SimpleTable/index.html)
-   [sna](../packages/sna/index.html)
-   [spBayes](../packages/spBayes/index.html)
-   [spikeslab](../packages/spikeslab/index.html)
-   [spikeSlabGAM](../packages/spikeSlabGAM/index.html)
-   [spTimer](../packages/spTimer/index.html)
-   [stochvol](../packages/stochvol/index.html)
-   [tgp](../packages/tgp/index.html)
-   [tRophicPosition](../packages/tRophicPosition/index.html)
-   [zic](../packages/zic/index.html)

### Related links:

-   [Bayesian Statistics and Marketing
    (bayesm)](http://www.perossi.org/home/bsm-1)
-   [BayesX](http://www.BayesX.org/)
-   [Just Another Gibbs Sampler
    (JAGS)](http://mcmc-jags.sourceforge.net/)
-   [MCMCpack](http://mcmcpack.berkeley.edu/)
-   [The BUGS Project
    (WinBUGS)](http://www.mrc-bsu.cam.ac.uk/software/bugs/)
-   [OpenBUGS](http://www.openbugs.net/)
-   [BRugs in
    CRANextras](http://www.stats.ox.ac.uk/pub/RWin/src/contrib/)
-   [NIMBLE](http://r-nimble.org/)
-   [Stan](http://mc-stan.org/)
-   CRAN Task View: [gR](gR.html)
-   Bioconductor Package:
    [[vbmp]{.BioC}](https://www.Bioconductor.org/packages/release/bioc/html/vbmp.html)
-   [BOA](http://www.public-health.uiowa.edu/boa/)
