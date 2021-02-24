## CRAN Task View: Statistics for the Social Sciences

  ----------------- --------------------------------------------------
  **Maintainer:**   John Fox
  **Contact:**      jfox at mcmaster.ca
  **Version:**      2018-05-08
  **URL:**          <https://CRAN.R-project.org/view=SocialSciences>
  ----------------- --------------------------------------------------

<div>

Social scientists use a wide range of statistical methods, most of which
are not unique to the social sciences. Indeed, most statistical data
analysis in the social sciences is covered by the facilities in the base
and recommended packages, which are part of the standard R distribution.
In the package descriptions below, I identify base and recommended
packages on first mention; packages that are not specifically identified
as \"R-base\" or \"recommended\" are contributed packages.

[**Other Relevant Task Views:**]{#taskviews}

Beyond the base and contributed packages, many of the methods commonly
employed in the social sciences are covered extensively in other CRAN
task views, including the following. I will try to minimize duplicating
information present in these other task views, given here in
alphabetical order.

-   [Bayesian](Bayesian.html): Methods of Bayesian inference in a
    variety of settings of interest to social scientists, including
    mixed-effects models.
-   [Econometrics](Econometrics.html) and [Finance](Finance.html): In
    addition to methods of specific interest to economists and financial
    analysts, these task views covers a variety of commonly used
    regression models and methods, instrumental-variables estimation,
    models for panel data, and some time-series models.
-   [MetaAnalysis](MetaAnalysis.html): Methods of meta analysis for
    combining results from primary studies. If data on individuals in
    each study are available, meta analysis can be performed using
    [mixed-effects models](#mixed-models) .
-   [Multivariate](Multivariate.html): A broad, if far from exhaustive,
    catalog of methods implemented in R for analyzing multivariate data,
    from data visualization to statistical modeling, and including
    correspondence analysis for multivariate categorical data.
-   [OfficialStatistics](OfficialStatistics.html): Covers not only
    official statistics but also methods for collecting and analyzing
    data from complex sample surveys, such as the
    [survey](../packages/survey/index.html) package.
-   [Psychometrics](Psychometrics.html): Extensively covers methods of
    scale construction, including item-response theory, multidimensional
    scaling, and classical test theory, along with other topics of
    interest in the social sciences, such as structural-equation
    modeling.
-   [Spatial](Spatial.html): Methods for managing, visualizing, and
    modeling spatial data, including spatial regression analysis.
-   [SpatioTemporal](SpatioTemporal.html): Methods for representing,
    visualizing, and analyzing data with information both on time and
    location.
-   [Survival](Survival.html): Methods for survival analysis (often
    termed \"event-history analysis\" in the social sciences), beyond
    the basic and standard methods, such as for Cox regression, included
    in the recommended [survival](../packages/survival/index.html)
    package.
-   [TimeSeries](TimeSeries.html): Methods for representing,
    manipulating, visualizing, and modeling time-series data, including
    time-series regression methods.

It is noteworthy that this enumeration includes about a third of the
CRAN task views. Moreover, there are other task views of potential
interest to social scientists (such as the [Graphics](Graphics.html)
task view on statistical graphics); I suggest that you look at the [list
of all task views on CRAN](https://cran.r-project.org/web/views/) .

**Linear and Generalized Linear Models:**

Univariate and multivariate linear models are fit by the `lm` function,
generalized linear models by the `glm` function, both in the R-base
**stats** package. Beyond `summary` and `plot` methods for `lm` and
`glm` objects, there is a wide array of functions that support these
objects.

-   The generic `anova` function in the **stats** package constructs
    sequential (\"Type-I\") analysis of variance and analysis of
    deviance tables, and can also compute *F* and chisquare
    likelihood-ratio tests for nested models. (It is typical for other
    classes of statistical models in R to have `anova` methods as well,
    along with methods for other standard generics, such as `coef`, for
    returning regression coefficients; `vcov` for the coefficient
    covariance matrix; `residuals`; and `fitted` for fitted values of
    the response.) The generic `Anova` function in the
    [car](../packages/car/index.html) package (associated with Fox and
    Weisberg, *An R Companion to Applied Regression, Second Edition* ,
    Sage, 2011) constructs so-called \"Type-II\" and \"Type-III\"
    partial tests for linear, generalized linear, and many other classes
    of regression models.
-   *F* and chisquare Wald tests for a variety of hypotheses are
    available from the `coeftest` and `waldtest` functions in the
    [lmtest](../packages/lmtest/index.html) package, and the
    `linearHypothesis` function in the [car](../packages/car/index.html)
    package. All of these functions permit the use of heteroscedasticity
    and heteroscedasticity/autocorrelation-consistent covariance
    matrices, as computed, e.g., by functions in the
    [sandwich](../packages/sandwich/index.html) and
    [car](../packages/car/index.html) packages. Also see the `glh.test`
    function in the [gmodels](../packages/gmodels/index.html) package.
    Nonlinear functions of parameters can be tested via the
    `deltaMethod` function in the [car](../packages/car/index.html)
    package. The [multcomp](../packages/multcomp/index.html) package
    includes functions for multiple comparisons. The `vuong` function in
    the [pscl](../packages/pscl/index.html) package tests non-nested
    hypotheses for generalized linear and some other models. Also see
    the [rms](../packages/rms/index.html) package for tests on linear
    and generalized linear models.
-   The standard R distribution has excellent basic facilities for
    linear and generalized linear model \"diagnostics,\" including, for
    example, hat-values and deletion statistics such as studentized
    residuals and Cook\'s distances ( `hatvalues`, `rstudent`, and
    `cooks.distance`, all in the **stats** package). These are augmented
    by other packages: several functions in the
    [car](../packages/car/index.html) package, which emphasizes
    graphical methods, e.g., `crPlots` for component-plus-residual plots
    and `avPlots` for added-variable plots (among others), in addition
    to numerical diagnostics, such `vif` for (generalized)
    variance-inflation factors; the [dr](../packages/dr/index.html)
    package for dimension reduction in regression, including SIR, SAVE,
    and pHd; and the [lmtest](../packages/lmtest/index.html) package,
    which implements a variety of diagnostic tests (e.g., for
    heteroscedasticity, nonlinearity, and autocorrelation). Other
    collinearity diagnostics are in the
    [perturb](../packages/perturb/index.html) package. Diagnostics may
    also be found in the [rms](../packages/rms/index.html) package. See
    the [influence.ME](../packages/influence.ME/index.html) package for
    influential-data diagnostics for mixed-effects models.
-   Several packages contain functions that are useful for interpreting
    linear and generalized linear models that have been fit to data: The
    [qvcalc](../packages/qvcalc/index.html) packages computes \"quasi
    variances\" for factors in linear and generalized linear models (and
    more generally). The [effects](../packages/effects/index.html)
    package constructs effect displays, including, e.g., \"adjusted
    means,\" for linear, generalized linear, and many other regression
    models; diagnostic partial-residual plots are available for linear
    and generalized linear models. Similar, if somewhat less general,
    plots are available in the [visreg](../packages/visreg/index.html)
    package. The [lsmeans](../packages/lsmeans/index.html) implements
    so-called \"least-squares means\" for linear, generalized linear,
    and mixed models, and includes provisions for hypothesis tests. The
    [Zelig](../packages/Zelig/index.html) package (see under
    [\"Collections\"](#Collections) ) creates summary displays for many
    kinds of statistical models.

**Analysis of Categorical and Count Data:**

Binomial logit and probit models, as well as Poisson-regression and
loglinear models for contingency tables (including models for
\"over-dispersed\" binomial and Poisson data), can be fit with the `glm`
function in the **stats** package. For over-dispersed data, see also the
[aod](../packages/aod/index.html) package, the
[dispmod](../packages/dispmod/index.html) package, and the `glm.nb`
function in the recommended [MASS](../packages/MASS/index.html) package
(associated with Venables and Ripley, *Modern Applied Statistics in S,
Fourth Ed.* , Springer, 2002), which fits negative-binomial GLMs. The
[pscl](../packages/pscl/index.html) package includes functions for
fitting zero-inflated and hurdle regression models to count data. The
multinomial logit model is fit by the `multinom` function in the
recommended [nnet](../packages/nnet/index.html) package, and ordered
logit and probit models by the `polr` function in the
[MASS](../packages/MASS/index.html) package. Also see the
[mlogit](../packages/mlogit/index.html) for the multinomial logit model,
the [MNP](../packages/MNP/index.html) package for the multinomial probit
model, and the [multinomRob](../packages/multinomRob/index.html) package
for the analysis of overdispersed multinomial data. The
[VGAM](../packages/VGAM/index.html) package is capable of fitting a very
wide variety of fixed-effect regression models within a unified
framework, including models for ordered and unordered categorical
responses and for count data.

There are other noteworthy facilities for analyzing categorical and
count data.

-   The `table` function in the R-base **base** package and the `xtabs`
    and `ftable` functions in the **stats** package construct
    contingency tables.
-   The `chisq.test` and `fisher.test` functions in the **stats**
    package may be used to test for independence in two-way contingency
    tables.
-   The `loglm` and `loglin` functions in the
    [MASS](../packages/MASS/index.html) package fit hierarchical
    loglinear models to contingency tables, the former as a front end to
    `glm`, the latter by iterative proportional fitting.
-   See the [brglm](../packages/brglm/index.html) and
    [logistf](../packages/logistf/index.html) packages for
    bias-reduction in binomial-response GLMs (useful, e.g., in cases of
    complete separation); the [elrm](../packages/elrm/index.html)
    package, which approximates exact conditional inference in logistic
    regression; the
    [exactLoglinTest](../packages/exactLoglinTest/index.html) package
    for exact tests of loglinear models; the `clogit` function in the
    [survival](../packages/survival/index.html) package for conditional
    logistic regression; and the [vcd](../packages/vcd/index.html)
    package for graphical displays of categorical data, including mosaic
    plots.
-   The [gnm](../packages/gnm/index.html) package estimates generalized
    *nonlinear* models, and can be used, e.g., to fit certain
    specialized models to mobility tables. The
    [logmult](../packages/logmult/index.html) package provides
    convenience functions based on [gnm](../packages/gnm/index.html) to
    fit log-multiplicative (e.g., UNIDIFF) and association (e.g.,
    Goodman\'s RC) models. Also see the
    [catspec](../packages/catspec/index.html) package for estimating
    various special models for square tables.
-   As previously mentioned, the [Multivariate](Multivariate.html) task
    view covers correspondence analysis of multivariate categorical
    data.
-   See the [betareg](../packages/betareg/index.html) package for beta
    regression of data on rates and proportions, a topic closely
    associated with categorical data.

**Other Regression Models:**

It is possible to fit a very wide variety of regression models with the
facilities provided by the base and recommended packages, and an even
wider variety of models with contributed packages, in addition to those
covered extensively in [other task views](#taskviews) .

-   *Nonlinear regression:* The `nls` function in the **stats** package
    fits nonlinear models by least-squares. The
    [nlstools](../packages/nlstools/index.html) includes several
    functions for assessing models fit by `nls`.
-   [*Mixed-effects models:*]{#mixed-models} The recommended
    [nlme](../packages/nlme/index.html) package, associated with
    Pinheiro and Bates, *Mixed-Effects Models in S and S-PLUS*
    (Springer, 2000), fits linear ( `lme`) and nonlinear ( `nlme`)
    mixed-effects models, commonly used in the social sciences for
    hierarchical and longitudinal data. Generalized linear mixed-effects
    models may be fit by the `glmmPQL` function in the
    [MASS](../packages/MASS/index.html) package, or (preferably) by the
    `glmer` function in the [lme4](../packages/lme4/index.html) package.
    The [lme4](../packages/lme4/index.html) package also largely
    supersedes [nlme](../packages/nlme/index.html) for *linear* mixed
    models, via its `lmer` function. Unlike `lme`, `lmer` supports
    crossed random effects, but does not support autocorrelated or
    heteroscedastic individual-level errors. Also see the
    [lmeSplines](../packages/lmeSplines/index.html),
    [lmm](../packages/lmm/index.html), and
    [MCMCglmm](../packages/MCMCglmm/index.html) packages.
-   *Generalized estimating equations:* The
    [gee](../packages/gee/index.html) and
    [geepack](../packages/geepack/index.html) packages fit marginal
    models by generalized estimating equations; see the
    [multgee](../packages/multgee/index.html) package for GEE estimation
    of models for correlated nominal or ordinal multinomial responses.
-   *Nonparametric regression analysis:* This is one of the conspicuous
    strengths of R. The standard R distribution includes several
    functions for smoothing scatterplots, including `loess.smooth` and
    `smooth.spline`, both in the **stats** package. The `loess`
    function, also in the **stats** package, fits simple and multiple
    nonparametric-regression models by local polynomial regression.
    Generalized additive models are covered by several packages,
    including the recommended [mgcv](../packages/mgcv/index.html)
    package and the [gam](../packages/gam/index.html) package, the
    latter associated with Hastie and Tibshirani, *Generalized Additive
    Models* (Chapman and Hall, 1990); also see the
    [VGAM](../packages/VGAM/index.html) package. Some other noteworthy
    contributed packages in this area are
    [gss](../packages/gss/index.html), which fits spline regressions;
    [locfit](../packages/locfit/index.html), for local-polynomial
    regression (and also density estimation) (Loader, *Local Regression
    and Likelihood,* Springer, 1999); [sm](../packages/sm/index.html),
    for a variety of smoothing techniques, including for regression
    (Bowman and Azzalini, *Applied Smoothing Techniques for Data
    Analysis,* Oxford, 1997); [np](../packages/np/index.html), which
    implements kernel smoothing methods for mixed data types; and
    [acepack](../packages/acepack/index.html) for ACE (alternating
    conditional expectations) and AVAS (additivity and variance
    stabilization) nonparametric transformation of the response and
    explanatory variables in regression.
-   *Quantile regression:* Methods for linear, nonlinear, and
    nonparametric quantile regression are extensively provided by the
    [quantreg](../packages/quantreg/index.html) package.
-   *Regression splines:* Parametric regression splines (as opposed to
    nonparametric smoothing splines), supported by the base-R
    **splines** package, can be used by `lm`, `glm`, and other
    statistical modeling functions that employ model formulas. See the
    `bs` (B-spline) and `ns` (natural spline) functions.
-   *Very large data sets:* The [biglm](../packages/biglm/index.html)
    package can fit linear and generalized linear models to data sets
    too large to fit in memory.

**Other Statistical Methods:**

Here is a brief survey of implementations in R of other statistical
methods commonly used by social scientists.

-   *Missing Data:* Several packages implement methods for handling
    missing data by multiple imputation, including the (conspicuously
    aging) [mix](../packages/mix/index.html),
    [norm](../packages/norm/index.html), and
    [pan](../packages/pan/index.html) packages associated with Shafer,
    *Analysis of Incomplete Multivariate Data* (Chapman and Hall, 1997),
    and the newer and more actively maintained
    [Amelia](../packages/Amelia/index.html),
    [mi](../packages/mi/index.html),
    [mice](../packages/mice/index.html), and
    [mitools](../packages/mitools/index.html) packages (the latter for
    drawing inferences from multiply imputed data sets). There are also
    some facilities for missing-data imputation in the general
    [Hmisc](../packages/Hmisc/index.html) package, which is described
    below, under [\"Collections\"](#Collections) . The
    [mvnmle](../packages/mvnmle/index.html) package finds the maximum
    likelihood estimates of means and covariances assuming
    multivariate-normal data. As well, some of the structural-equation
    modeling software discussed in the
    [Psychometrics](Psychometrics.html) taskview is capable of
    maximum-likelihood estimation of regression models with missing
    data. The [VIM](../packages/VIM/index.html) package has functions
    for visualizing missing and imputed values.
-   *Bootstrapping and Other Resampling Methods:* The recommended
    package [boot](../packages/boot/index.html), associated with Davison
    and Hinkley, *Bootstrap Methods and Their Application* (Cambridge,
    1997), has excellent facilities for bootstrapping and some related
    methods. Also notable is the
    [bootstrap](../packages/bootstrap/index.html) package, associated
    with Efron and Tibshirani, *An Introduction to the Bootstrap*
    (Chapman and Hall, 1993), which has functions for bootstrapping and
    jackknifing. In addition, see the functions `Boot` and `bootCase` in
    the [car](../packages/car/index.html) package, and `nlsBoot` in the
    [nlstools](../packages/nlstools/index.html) package, along with the
    [simpleboot](../packages/simpleboot/index.html) package.
-   *Model Selection:* The `step` function in the **stats** package and
    the more broadly applicable `stepAIC` function in the
    [MASS](../packages/MASS/index.html) package perform forward,
    backward, and forward-backward stepwise selection for a variety of
    statistical models. The `regsubsets` function in the
    [leaps](../packages/leaps/index.html) package performs all-subsets
    regression. The [BMA](../packages/BMA/index.html) package performs
    Bayesian model averaging. The standard `AIC` and `BIC` functions are
    also relevant to model selection. Beyond these packages and
    functions, see the [MachineLearning](MachineLearning.html) task
    view.
-   *Social Network Analysis:* There are several packages useful for
    social network analysis, including [sna](../packages/sna/index.html)
    for sociometric analysis of networks (e.g., blockmodeling),
    [network](../packages/network/index.html) for manipulating and
    displaying network objects,
    [latentnet](../packages/latentnet/index.html) for latent position
    and cluster models for networks, [ergm](../packages/ergm/index.html)
    for exponential random graph models of networks, and the
    \"metapackage\" [statnet](../packages/statnet/index.html), all
    associated with the [statnet project](http://statnet.org/) . Also
    see the [RSiena](../packages/RSiena/index.html) and
    [PAFit](../packages/PAFit/index.html) packages for longitudinal
    social network analysis; and the
    [multiplex](../packages/multiplex/index.html) package, which
    implements algebraic procedures for the analysis of multiple social
    networks.
-   *Propensity Scores and Matching:* See the
    [Matching](../packages/Matching/index.html),
    [MatchIt](../packages/MatchIt/index.html),
    [optmatch](../packages/optmatch/index.html), and
    [PSAgraphics](../packages/PSAgraphics/index.html) packages, and the
    `matching` function in the [arm](../packages/arm/index.html) package
    (associated with Gelman and Hill, *Data Analysis Using Regression
    and Multilevel/Hierarchical Models,* Cambridge, 2007).
-   *Demographic methods* : The
    [demography](../packages/demography/index.html) package includes
    functions for constructing life tables, for analyzing mortality,
    fertility, and immigration, and for forecasting population.

[**Collections of Functions:**]{#Collections}

There are some packages that are so heterogeneous that they are
difficult to classify, yet contain functions (typically in multiple
domains) that are of interest to social scientists:

-   I have already made several references to the recommended
    [MASS](../packages/MASS/index.html) package, which is associated
    with Venables and Ripley\'s *Modern Applied Statistics With S.*
    Other recommended packages associated with this book are
    [nnet](../packages/nnet/index.html), for fitting neural networks
    (but also, as mentioned, multinomial logistic-regression models);
    [spatial](../packages/spatial/index.html) for spatial statistics;
    and [class](../packages/class/index.html), which contains functions
    for classification.
-   I\'ve also mentioned the [car](../packages/car/index.html) package,
    associated with Fox and Weisberg, *An R Companion to Applied
    Regression, Second Edition* , which has a variety of functions
    supporting regression analysis, data exploration, and data
    transformation.
-   The [Hmisc](../packages/Hmisc/index.html) and
    [rms](../packages/rms/index.html) packages (both mentioned above),
    associated with Harrell, *Regression Modeling Strategies, Second
    Edition* (Springer, 2015), provide functions for data manipulation,
    linear models, logistic-regression models, and survival analysis,
    many of them \"front ends\" to or modifications of other facilities
    in R.
-   The [Zelig](../packages/Zelig/index.html) package integrates a wide
    array of statistical models of interest to social scientists (see
    the [Zelig web site](http://gking.harvard.edu/zelig) for details).

**Acknowledgments:**

Jangman Hong contributed to the general revision of this task view, as
did other individuals who made a variety of specific suggestions.

If I have omitted something of importance not covered in one of the
other task views cited, or if a new package or function should be
mentioned here, [please let me know.](mailto:jfox@mcmaster.ca)

Compilation of this task view was partly supported by grants from the
Social Sciences and Humanities Research Council of Canada.

</div>

### CRAN packages:

-   [acepack](../packages/acepack/index.html)
-   [Amelia](../packages/Amelia/index.html)
-   [aod](../packages/aod/index.html)
-   [arm](../packages/arm/index.html)
-   [betareg](../packages/betareg/index.html)
-   [biglm](../packages/biglm/index.html)
-   [BMA](../packages/BMA/index.html)
-   [boot](../packages/boot/index.html) (core)
-   [bootstrap](../packages/bootstrap/index.html)
-   [brglm](../packages/brglm/index.html)
-   [car](../packages/car/index.html) (core)
-   [catspec](../packages/catspec/index.html)
-   [class](../packages/class/index.html)
-   [demography](../packages/demography/index.html)
-   [dispmod](../packages/dispmod/index.html)
-   [dr](../packages/dr/index.html)
-   [effects](../packages/effects/index.html) (core)
-   [elrm](../packages/elrm/index.html)
-   [ergm](../packages/ergm/index.html)
-   [exactLoglinTest](../packages/exactLoglinTest/index.html)
-   [gam](../packages/gam/index.html) (core)
-   [gee](../packages/gee/index.html)
-   [geepack](../packages/geepack/index.html)
-   [gmodels](../packages/gmodels/index.html)
-   [gnm](../packages/gnm/index.html)
-   [gss](../packages/gss/index.html)
-   [Hmisc](../packages/Hmisc/index.html) (core)
-   [influence.ME](../packages/influence.ME/index.html)
-   [latentnet](../packages/latentnet/index.html)
-   [leaps](../packages/leaps/index.html)
-   [lme4](../packages/lme4/index.html) (core)
-   [lmeSplines](../packages/lmeSplines/index.html)
-   [lmm](../packages/lmm/index.html)
-   [lmtest](../packages/lmtest/index.html) (core)
-   [locfit](../packages/locfit/index.html)
-   [logistf](../packages/logistf/index.html)
-   [logmult](../packages/logmult/index.html)
-   [lsmeans](../packages/lsmeans/index.html) (core)
-   [MASS](../packages/MASS/index.html) (core)
-   [Matching](../packages/Matching/index.html)
-   [MatchIt](../packages/MatchIt/index.html)
-   [MCMCglmm](../packages/MCMCglmm/index.html) (core)
-   [mgcv](../packages/mgcv/index.html) (core)
-   [mi](../packages/mi/index.html) (core)
-   [mice](../packages/mice/index.html) (core)
-   [mitools](../packages/mitools/index.html)
-   [mix](../packages/mix/index.html)
-   [mlogit](../packages/mlogit/index.html)
-   [MNP](../packages/MNP/index.html)
-   [multcomp](../packages/multcomp/index.html) (core)
-   [multgee](../packages/multgee/index.html)
-   [multinomRob](../packages/multinomRob/index.html)
-   [multiplex](../packages/multiplex/index.html)
-   [mvnmle](../packages/mvnmle/index.html)
-   [network](../packages/network/index.html)
-   [nlme](../packages/nlme/index.html) (core)
-   [nlstools](../packages/nlstools/index.html)
-   [nnet](../packages/nnet/index.html) (core)
-   [norm](../packages/norm/index.html)
-   [np](../packages/np/index.html)
-   [optmatch](../packages/optmatch/index.html)
-   [PAFit](../packages/PAFit/index.html)
-   [pan](../packages/pan/index.html)
-   [perturb](../packages/perturb/index.html)
-   [PSAgraphics](../packages/PSAgraphics/index.html)
-   [pscl](../packages/pscl/index.html)
-   [quantreg](../packages/quantreg/index.html) (core)
-   [qvcalc](../packages/qvcalc/index.html)
-   [rms](../packages/rms/index.html) (core)
-   [RSiena](../packages/RSiena/index.html)
-   [sandwich](../packages/sandwich/index.html) (core)
-   [simpleboot](../packages/simpleboot/index.html)
-   [sm](../packages/sm/index.html)
-   [sna](../packages/sna/index.html)
-   [spatial](../packages/spatial/index.html)
-   [statnet](../packages/statnet/index.html)
-   [survey](../packages/survey/index.html) (core)
-   [survival](../packages/survival/index.html) (core)
-   [vcd](../packages/vcd/index.html)
-   [VGAM](../packages/VGAM/index.html) (core)
-   [VIM](../packages/VIM/index.html)
-   [visreg](../packages/visreg/index.html)
-   [Zelig](../packages/Zelig/index.html)

### Related links:

-   CRAN Task View: [Bayesian](Bayesian.html)
-   CRAN Task View: [Econometrics](Econometrics.html)
-   CRAN Task View: [Finance](Finance.html)
-   CRAN Task View: [Graphics](Graphics.html)
-   CRAN Task View: [MachineLearning](MachineLearning.html)
-   CRAN Task View: [MetaAnalysis](MetaAnalysis.html)
-   CRAN Task View: [Multivariate](Multivariate.html)
-   CRAN Task View: [OfficialStatistics](OfficialStatistics.html)
-   CRAN Task View: [Psychometrics](Psychometrics.html)
-   CRAN Task View: [Spatial](Spatial.html)
-   CRAN Task View: [SpatioTemporal](SpatioTemporal.html)
-   CRAN Task View: [Survival](Survival.html)
-   CRAN Task View: [TimeSeries](TimeSeries.html)
-   [Fox and Weisberg, *An R Companion to Applied Regression, Second
    Edition*
    website](http://socserv.mcmaster.ca/jfox/Books/Companion/index.html)
-   [Harrell, *Regression Modeling Strategies, Second Edition*
    website](http://biostat.mc.vanderbilt.edu/wiki/Main/RmS)
-   [Statnet Project website](http://statnet.org/)
-   [Venables and Ripley, *Modern Applied Statistics with S, 4th ed.*
    website](http://www.stats.ox.ac.uk/pub/MASS4/)
-   [Zelig Software website](http://gking.harvard.edu/zelig)
