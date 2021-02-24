## CRAN Task View: Meta-Analysis

  ----------------- ------------------------------------------------
  **Maintainer:**   Michael Dewey
  **Contact:**      lists at dewey.myzen.co.uk
  **Version:**      2018-05-10
  **URL:**          <https://CRAN.R-project.org/view=MetaAnalysis>
  ----------------- ------------------------------------------------

<div>

This task view covers packages which include facilities for
meta-analysis of summary statistics from primary studies. The task view
does not consider the meta-analysis of individual participant data (IPD)
which can be handled by any of the standard linear modelling functions
but does include some packages which offer special facilities for IPD.

The standard meta-analysis model is a form of weighted least squares and
so any of the wide range of R packages providing weighted least squares
would in principle be able to fit the model. The advantage of using a
specialised package is that (a) it takes care of the small tweaks
necessary (b) it provides a range of ancillary functions for displaying
and investigating the model. Where the model is referred to below it is
this model which is meant.

Where summary statistics are not available a meta-analysis of
significance levels is possible. This is not completely unconnected with
the problem of adjustment for multiple comparisons but the packages
below which offer this, chiefly in the context of genetic data, also
offer additional functionality.

#### Univariate meta-analysis

*Preparing for meta-analysis*

-   The primary studies often use a range of statistics to present their
    results. Convenience functions to convert these onto a common metric
    are presented by: [compute.es](../packages/compute.es/index.html)
    which converts from various statistics to d, g, r, z and the log
    odds ratio, [MAc](../packages/MAc/index.html) which converts to
    correlation coefficients, [MAd](../packages/MAd/index.html) which
    converts to mean differences, and
    [metafor](../packages/metafor/index.html) which converts to effect
    sizes an extensive set of measures for comparative studies (such as
    binary data, person years, mean differences and ratios and so on),
    for studies of association (a wide range of correlation types), for
    non-comparative studies (proportions, incidence rates, and mean
    change). It also provides for a measure used in psychometrics
    (Cronbach\'s alpha). [esc](../packages/esc/index.html) provides a
    range of effect size calculations with partial overlap with
    [metafor](../packages/metafor/index.html) but with some extras,
    noticeably for converting test statistics, also includes a
    convenience function for collating its output for input to another
    package like [metafor](../packages/metafor/index.html) or producing
    a CSV file. [effsize](../packages/effsize/index.html) contains
    functions to compute effect sizes mean difference (Cohen\'s d and
    Hedges g), dominance matrices (Cliff\'s Delta) and stochastic
    superiority (Vargha-Delaney A).
    [psychmeta](../packages/psychmeta/index.html) provides extensive
    facilties for converting effect sizes and for correcting for a
    variety of restrictions and measurement errors.
-   [meta](../packages/meta/index.html) provides functions to read and
    work with files output by RevMan 4 and 5.
-   [metagear](../packages/metagear/index.html) provides many tools for
    the systematic review process including screening articles,
    downloading the articles, generating a PRISMA diagram, and some
    tools for effect sizes. [revtools](../packages/revtools/index.html)
    provides tools for downloading from bibliographic databases and uses
    machine learning methods to process them.
-   [metavcov](../packages/metavcov/index.html) computes the
    variance-covariance matrix for multivariate meta-analysis when
    correlations between outcomes can be provided but not between
    treatment effects, and
    [clubSandwich](../packages/clubSandwich/index.html) imputes
    variance-covariance matrix for multivariate meta-analysis
-   [metafuse](../packages/metafuse/index.html) uses a fused lasso to
    merge covariate estimates across a number of independent datasets.

*Fitting the model*

-   Four packages provide the inverse variance weighted,
    Mantel-Haenszel, and Peto methods:
    [epiR](../packages/epiR/index.html),
    [meta](../packages/meta/index.html),
    [metafor](../packages/metafor/index.html), and
    [rmeta](../packages/rmeta/index.html).
-   For binary data [metafor](../packages/metafor/index.html) provides
    the binomial-normal model.
-   For sparse binary data [exactmeta](../packages/exactmeta/index.html)
    provides an exact method which does not involve continuity
    corrections.
-   Packages which work with specific effect sizes may be more congenial
    to workers in some areas of science and include
    [MAc](../packages/MAc/index.html) and
    [metacor](../packages/metacor/index.html) which provide
    meta-analysis of correlation coefficients and
    [MAd](../packages/MAd/index.html) which provides meta-analysis of
    mean differences. [MAc](../packages/MAc/index.html) and
    [MAd](../packages/MAd/index.html) provide a range of graphics.
    [psychometric](../packages/psychometric/index.html) provides an
    extensive range of functions for the meta-analysis of psychometric
    studies.
-   [psychmeta](../packages/psychmeta/index.html) implements the
    Hunter-Schmidt method including corrections for reliability and
    range-restriction issues
-   Bayesian approaches are contained in various packages.
    [bspmma](../packages/bspmma/index.html) which provides two different
    models: a non-parametric and a semi-parametric. Graphical display of
    the results is provided. [metamisc](../packages/metamisc/index.html)
    provides a method with priors suggested by Higgins.
    [mmeta](../packages/mmeta/index.html) provides meta-analysis using
    beta-binomial prior distributions. A Bayesian approach is also
    provided by [bmeta](../packages/bmeta/index.html) which provides
    forest plots via [forestplot](../packages/forestplot/index.html) and
    diagnostic graphical output.
    [bayesmeta](../packages/bayesmeta/index.html) includes shrinkage
    estimates, posterior predictive p-values and forest plots via either
    [metafor](../packages/metafor/index.html) or
    [forestplot](../packages/forestplot/index.html). Diagnostic
    graphical output is available.
-   Some packages concentrate on providing a specialised version of the
    core meta-analysis function without providing the range of ancillary
    functions. These are: [gmeta](../packages/gmeta/index.html) which
    subsumes a very wide variety of models under the method of
    confidence distributions and also provides a graphical display,
    [metaLik](../packages/metaLik/index.html) which uses a more
    sophisticated approach to the likelihood,
    [metamisc](../packages/metamisc/index.html) which as well as the
    method of moments provides two likelihood-based methods, and
    [metatest](../packages/metatest/index.html) which provides another
    improved method of obtaining confidence intervals,
    [metaBMA](../packages/metaBMA/index.html) has a Bayesian approach
    using model averaging, a variety of priors are provided and it is
    possible for the user to define new ones.
-   [metagen](../packages/metagen/index.html) provides a range of
    methods for random effects models and also facilities for extensive
    simulation studies of the properties of those methods.
-   [metaplus](../packages/metaplus/index.html) fits random effects
    models relaxing the usual assumption that the random effects have a
    normal distribution by providing t or a mixture of normals.
-   [ratesci](../packages/ratesci/index.html) fits random effects models
    to binary data using a variety of methods for confidence intervals.
-   [RandMeta](../packages/RandMeta/index.html) estimates exact
    confidence intervals in random effects models using an efficient
    algorithm.
-   [rma.exact](../packages/rma.exact/index.html) estimates exact
    confidence intervals in random effects normal-normal models and also
    provides plots of them.
-   [clubSandwich](../packages/clubSandwich/index.html) gives
    cluster-robust variance estimates.
-   [pimeta](../packages/pimeta/index.html) implements prediction
    intervals for random effects meta-analysis.
-   [RBesT](../packages/RBesT/index.html) generates a meta-analytic
    prior for future studies. Forest plots are available as well as
    diagnostic plots. It uses Stan as its engine.

*Graphical methods*

An extensive range of graphical procedures is available.

-   Forest plots are provided in
    [forestmodel](../packages/forestmodel/index.html) (using ggplot2),
    [forestplot](../packages/forestplot/index.html),
    [meta](../packages/meta/index.html),
    [metafor](../packages/metafor/index.html),
    [psychmeta](../packages/psychmeta/index.html), and
    [rmeta](../packages/rmeta/index.html). Although the most basic plot
    can be produced by any of them they each provide their own choice of
    enhancements.
-   Funnel plots are provided in [meta](../packages/meta/index.html),
    [metafor](../packages/metafor/index.html),
    [psychometric](../packages/psychometric/index.html) and
    [rmeta](../packages/rmeta/index.html). In addition to the standard
    funnel plots an enhanced funnel plot to assess the impact of extra
    evidence is available in
    [extfunnel](../packages/extfunnel/index.html), a funnel plot for
    limit meta-analysis in [metasens](../packages/metasens/index.html),
    and [metaviz](../packages/metaviz/index.html) provides funnel plots
    in the context of visual inference.
-   Radial (Galbraith) plots are provided in
    [meta](../packages/meta/index.html) and
    [metafor](../packages/metafor/index.html).
-   L\'Abbe plots are provided in [meta](../packages/meta/index.html)
    and [metafor](../packages/metafor/index.html).
-   Baujat plots are provided in [meta](../packages/meta/index.html) and
    [metafor](../packages/metafor/index.html).
-   [metaplotr](../packages/metaplotr/index.html) provides a crosshair
    plot
-   [MetaAnalyser](../packages/MetaAnalyser/index.html) provides an
    interactive visualisation of the results of a meta-analysis.
-   [metaviz](../packages/metaviz/index.html) provides rainforestplots,
    an enhanced version of forest plots. It accepts input from
    [metafor](../packages/metafor/index.html).

*Investigating heterogeneity*

-   Confidence intervals for the heterogeneity parameter are provided in
    [metafor](../packages/metafor/index.html),
    [metagen](../packages/metagen/index.html), and
    [psychmeta](../packages/psychmeta/index.html).
-   [altmeta](../packages/altmeta/index.html) presents a variety of
    alternative methods for measuring and testing heterogeneity with a
    focus on robustness to outlying studies.
-   [hetmeta](../packages/hetmeta/index.html) calculates some extra
    measures of heterogeneity.
-   [metaforest](../packages/metaforest/index.html) investigates
    heterogeneity using random forests. Note that it has nothing to do
    with forest plots.

*Model criticism*

-   An extensive series of plots of diagnostic statistics is provided in
    [metafor](../packages/metafor/index.html).
-   [metaplus](../packages/metaplus/index.html) provides outlier
    diagnostics.
-   [psychmeta](../packages/psychmeta/index.html) provides leave-one-out
    methods.
-   [ConfoundedMeta](../packages/ConfoundedMeta/index.html) conducts a
    sensitivity analysis to estimate the proportion of studies with true
    effect sizes above a threshold.

*Investigating small study bias*

The issue of whether small studies give different results from large
studies has been addressed by visual examination of the funnel plots
mentioned above. In addition:

-   [meta](../packages/meta/index.html) and
    [metafor](../packages/metafor/index.html) provide both the
    non-parametric method suggested by Begg and Mazumdar and a range of
    regression tests modelled after the approach of Egger.
-   [xmeta](../packages/xmeta/index.html) provides a method in the
    context of multivariate meta-analysis.
-   An exploratory technique for detecting an excess of statistically
    significant studies is provided by
    [PubBias](../packages/PubBias/index.html).
-   [metamisc](../packages/metamisc/index.html) provides funnel plots
    and tests for asymmetry.

*Unobserved studies*

A recurrent issue in meta-analysis has been the problem of unobserved
studies.

-   Rosenthal\'s fail safe n is provided by
    [MAc](../packages/MAc/index.html) and
    [MAd](../packages/MAd/index.html).
    [metafor](../packages/metafor/index.html) provides it as well as two
    more recent methods by Orwin and Rosenberg.
-   Duval\'s trim and fill method is provided by
    [meta](../packages/meta/index.html) and
    [metafor](../packages/metafor/index.html).
-   [metasens](../packages/metasens/index.html) provides Copas\'s
    selection model and also the method of limit meta-analysis (a
    regression based approach for dealing with small study effects) due
    to Rücker et al.
-   [selectMeta](../packages/selectMeta/index.html) provides various
    selection models: the parametric model of Iyengar and Greenhouse,
    the non-parametric model of Dear and Begg, and proposes a new
    non-parametric method imposing a monotonicity constraint.
-   [SAMURAI](../packages/SAMURAI/index.html) performs a sensitivity
    analysis assuming the number of unobserved studies is known, perhaps
    from a trial registry, but not their outcome.
-   The [metansue](../packages/metansue/index.html) package allows the
    inclusion by multiple imputation of studies known only to have a
    non-significant result.
-   [weightr](../packages/weightr/index.html) provides facilities for
    using the weight function model of Vevea and Hedges.

*Other study designs*

-   [SCMA](../packages/SCMA/index.html) provides single case
    meta-analysis. It is part of a suite of packages dedicated to
    single-case designs.
-   [joint.Cox](../packages/joint.Cox/index.html) provides facilities
    for the meta-analysis of studies of joint time-to-event and disease
    progression.
-   [metamisc](../packages/metamisc/index.html) provides for
    meta-analysis of prognostic studies using the c statistic or the O/E
    ratio. Some plots are provided.

*Meta-analysis of significance values*

-   [metap](../packages/metap/index.html) provides some facilities for
    meta-analysis of significance values.
-   [aggregation](../packages/aggregation/index.html) provides a smaller
    subset of methods.
-   [TFisher](../packages/TFisher/index.html) provides Fisher\'s method
    using thresholding for the p-values.

Some methods are also provided in some of the genetics packages
mentioned below.

#### Multivariate meta-analysis

Standard methods outlined above assume that the effect sizes are
independent. This assumption may be violated in a number of ways: within
each primary study multiple treatments may be compared to the same
control, each primary study may report multiple endpoints, or primary
studies may be clustered for instance because they come from the same
country or the same research team. In these situations where the outcome
is multivariate:

-   [mvmeta](../packages/mvmeta/index.html) assumes the within study
    covariances are known and provides a variety of options for fitting
    random effects. [metafor](../packages/metafor/index.html) provides
    fixed effects and likelihood based random effects model fitting
    procedures. Both these packages include meta-regression,
    [metafor](../packages/metafor/index.html) also provides for
    clustered and hierarchical models.
-   [mvtmeta](../packages/mvtmeta/index.html) provides multivariate
    meta-analysis using the method of moments for random effects
    although not meta-regression,
-   [metaSEM](../packages/metaSEM/index.html) provides multivariate (and
    univariate) meta-analysis and meta-regression by embedding it in the
    structural equation framework and using OpenMx for the structural
    equation modelling. It can provide a three-level meta-analysis
    taking account of clustering and allowing for level 2 and level 3
    heterogeneity. It also provides via a two-stage approach
    meta-analysis of correlation or covariance matrices.
-   [xmeta](../packages/xmeta/index.html) provides various functions for
    multivariate meta-analysis and also for detecting publication bias.
-   [dosresmeta](../packages/dosresmeta/index.html) concentrates on the
    situation where individual studies have information on the
    dose-response relationship.
-   [robumeta](../packages/robumeta/index.html) provides robust variance
    estimation for clustered and hierarchical estimates.
-   [CIAAWconsensus](../packages/CIAAWconsensus/index.html) has a
    function for multivariate m-a in the context of atomic weights and
    estimating isotope ratios.

#### Meta-analysis of studies of diagnostic tests

A special case of multivariate meta-analysis is the case of summarising
studies of diagnostic tests. This gives rise to a bivariate, binary
meta-analysis with the within-study correlation assumed zero although
the between-study correlation is estimated. This is an active area of
research and a variety of methods are available including what is
referred to here as Reitsma\'s method, and the hierarchical summary
receiver operating characteristic (HSROC) method. In many situations
these are equivalent.

-   [mada](../packages/mada/index.html) provides various descriptive
    statistics and univariate methods (diagnostic odds ratio and Lehman
    model) as well as the bivariate method due to Reitsma. In addition
    meta-regression is provided. A range of graphical methods is also
    available.
-   [Metatron](../packages/Metatron/index.html) provides a method for
    the Reitsma model incuding the case of an imperfect reference
    standard.
-   [metamisc](../packages/metamisc/index.html) provides the method of
    Riley which estimates a common within and between correlation.
    Graphical output is also provided.
-   [bamdit](../packages/bamdit/index.html) provides Bayesian
    meta-analysis with a bivariate random effects model (using JAGS to
    implement the MCMC method). Graphical methods are provided.
-   [meta4diag](../packages/meta4diag/index.html) provides Bayesian
    inference analysis for bivariate meta-analysis of diagnostic test
    studies and an extensive range of graphical methods.
-   [CopulaREMADA](../packages/CopulaREMADA/index.html) uses a copula
    based mixed model
-   [diagmeta](../packages/diagmeta/index.html) considers the case where
    the primary studies provide analysis using multiple cut-offs.
    Graphical methods are also provided.

#### Meta-regression

Where suitable moderator variables are available they may be included
using meta-regression. All these packages are mentioned above, this just
draws that information together.

-   [metafor](../packages/metafor/index.html) provides meta-regression
    (multiple moderators are catered for). Various packages rely on
    [metafor](../packages/metafor/index.html) to provide meta-regression
    ([meta](../packages/meta/index.html),
    [MAc](../packages/MAc/index.html), and
    [MAd](../packages/MAd/index.html)) and all three of these provide
    bubble plots. [psychmeta](../packages/psychmeta/index.html) also
    uses [metafor](../packages/metafor/index.html).
-   [bmeta](../packages/bmeta/index.html),
    [metagen](../packages/metagen/index.html),
    [metaLik](../packages/metaLik/index.html),
    [metaSEM](../packages/metaSEM/index.html), and
    [metatest](../packages/metatest/index.html) also provide
    meta-regression.
-   [mvmeta](../packages/mvmeta/index.html) provides meta-regression for
    multivariate meta-analysis as do
    [metafor](../packages/metafor/index.html) and
    [metaSEM](../packages/metaSEM/index.html).
-   [metacart](../packages/metacart/index.html) integrates regression
    and classification trees into the meta-analysis framework for
    moderator selection.
-   [mada](../packages/mada/index.html) provides for the meta-regression
    of diagnostic test studies.

#### Individual participant data (IPD)

Where all studies can provide individual participant data then software
for analysis of multi-centre trials or multi-centre cohort studies
should prove adequate and is outside the scope of this task view. Other
packages which provide facilities related to IPD are:

-   [ipdmeta](../packages/ipdmeta/index.html) which uses information on
    aggregate summary statistics and a covariate of interest to assess
    whether a full IPD analysis would have more power.
-   [ecoreg](../packages/ecoreg/index.html) which is designed for
    ecological studies enables estimation of an individual level
    logistic regression from aggregate data or individual data.
-   [surrosurv](../packages/surrosurv/index.html) evaluates failure time
    surrogates in the context of IPD meta-analysis

#### Network meta-analysis

Also known as multiple treatment comparison. This is a very active area
of research and development. Note that some of the packages mentioned
above under multivariate meta-analysis can also be used for network
meta-analysis with appropriate setup.

This is provided in a Bayesian framework by
[gemtc](../packages/gemtc/index.html), which acts as a front-end to BUGS
or JAGS, and [pcnetmeta](../packages/pcnetmeta/index.html), which uses
JAGS. [nmaINLA](../packages/nmaINLA/index.html) uses integrated nested
Laplace approximations as an alternative to MCMC. It provides a number
of data-sets. [netmeta](../packages/netmeta/index.html) works in a
frequentist framework. Both
[pcnetmeta](../packages/pcnetmeta/index.html) and
[netmeta](../packages/netmeta/index.html) provide network graphs and
[netmeta](../packages/netmeta/index.html) provides a heatmap for
displaying inconsistency and heterogeneity.
[nmathresh](../packages/nmathresh/index.html) provides
decision-invariant bias adjustment thresholds and intervals the smallest
changes to the data that would result in a change of decision.

#### Genetics

There are a number of packages specialising in genetic data:
[CPBayes](../packages/CPBayes/index.html) uses a Bayesian approach to
study cross-phenotype genetic associations,
[etma](../packages/etma/index.html) proposes a new statistical method to
detect epistasis, [gap](../packages/gap/index.html) combines p-values,
[getmstatistic](../packages/getmstatistic/index.html) quantifies
systematic heterogeneity,
[MendelianRandomization](../packages/MendelianRandomization/index.html)
provides several methods for performing Mendelian randomisation analyses
with summarised data, [MetABEL](../packages/MetABEL/index.html) provides
meta-analysis of genome wide SNP association results,
[MetaIntegrator](../packages/MetaIntegrator/index.html) provides an
extensive set of functions for genetic studies,
[metaMA](../packages/metaMA/index.html) provides meta-analysis of
p-values or moderated effect sizes to find differentially expressed
genes, [MetaPath](../packages/MetaPath/index.html) performs
meta-analysis for pathway enrichment,
[MetaPCA](../packages/MetaPCA/index.html) provides meta-analysis in the
dimension reduction of genomic data,
[MetaQC](../packages/MetaQC/index.html) provides objective quality
control and inclusion/exclusion criteria for genomic meta-analysis,
[metaRNASeq](../packages/metaRNASeq/index.html) meta-analysis from
multiple RNA sequencing experiments,
[MultiMeta](../packages/MultiMeta/index.html) for meta-analysis of
multivariate GWAS results with graphics, designed to accept GEMMA
format, [MetaSKAT](../packages/MetaSKAT/index.html),
[seqMeta](../packages/seqMeta/index.html), provide meta-analysis for the
SKAT test.

#### Interfaces

[RcmdrPlugin.EZR](../packages/RcmdrPlugin.EZR/index.html) provides an
interface via the Rcmdr GUI using [meta](../packages/meta/index.html)
and [metatest](../packages/metatest/index.html) to do the heavy lifting,
[RcmdrPlugin.RMTCJags](../packages/RcmdrPlugin.RMTCJags/index.html)
provides an interface for network meta-analysis using BUGS code, and
[MAVIS](../packages/MAVIS/index.html) provides a Shiny interface using
[metafor](../packages/metafor/index.html),
[MAc](../packages/MAc/index.html), [MAd](../packages/MAd/index.html),
and [weightr](../packages/weightr/index.html).

#### Simulation

Extensive facilities for simulation are provided in
[metagen](../packages/metagen/index.html) including the ability to make
use of parallel processing.
[psychmeta](../packages/psychmeta/index.html) provides facilities for
simulation of psychometric data-sets.

#### Others

[CRTSize](../packages/CRTSize/index.html) provides meta-analysis as part
of a package primarily dedicated to the determination of sample size in
cluster randomised trials in particular by simulating adding a new study
to the meta-analysis.

[CAMAN](../packages/CAMAN/index.html) offers the possibility of using
finite semiparametric mixtures as an alternative to the random effects
model where there is heterogeneity. Covariates can be included to
provide meta-regression.

[joineRmeta](../packages/joineRmeta/index.html) provides functions for
meta-analysis of a single longitudinal and a single time-to-event
outcome from multiple studies using joint models

</div>

### CRAN packages:

-   [aggregation](../packages/aggregation/index.html)
-   [altmeta](../packages/altmeta/index.html)
-   [bamdit](../packages/bamdit/index.html)
-   [bayesmeta](../packages/bayesmeta/index.html)
-   [bmeta](../packages/bmeta/index.html)
-   [bspmma](../packages/bspmma/index.html)
-   [CAMAN](../packages/CAMAN/index.html)
-   [CIAAWconsensus](../packages/CIAAWconsensus/index.html)
-   [clubSandwich](../packages/clubSandwich/index.html)
-   [compute.es](../packages/compute.es/index.html)
-   [ConfoundedMeta](../packages/ConfoundedMeta/index.html)
-   [CopulaREMADA](../packages/CopulaREMADA/index.html)
-   [CPBayes](../packages/CPBayes/index.html)
-   [CRTSize](../packages/CRTSize/index.html)
-   [diagmeta](../packages/diagmeta/index.html)
-   [dosresmeta](../packages/dosresmeta/index.html)
-   [ecoreg](../packages/ecoreg/index.html)
-   [effsize](../packages/effsize/index.html)
-   [epiR](../packages/epiR/index.html)
-   [esc](../packages/esc/index.html)
-   [etma](../packages/etma/index.html)
-   [exactmeta](../packages/exactmeta/index.html)
-   [extfunnel](../packages/extfunnel/index.html)
-   [forestmodel](../packages/forestmodel/index.html)
-   [forestplot](../packages/forestplot/index.html)
-   [gap](../packages/gap/index.html)
-   [gemtc](../packages/gemtc/index.html)
-   [getmstatistic](../packages/getmstatistic/index.html)
-   [gmeta](../packages/gmeta/index.html)
-   [hetmeta](../packages/hetmeta/index.html)
-   [ipdmeta](../packages/ipdmeta/index.html)
-   [joineRmeta](../packages/joineRmeta/index.html)
-   [joint.Cox](../packages/joint.Cox/index.html)
-   [MAc](../packages/MAc/index.html)
-   [MAd](../packages/MAd/index.html)
-   [mada](../packages/mada/index.html)
-   [MAVIS](../packages/MAVIS/index.html)
-   [MendelianRandomization](../packages/MendelianRandomization/index.html)
-   [meta](../packages/meta/index.html) (core)
-   [meta4diag](../packages/meta4diag/index.html)
-   [MetaAnalyser](../packages/MetaAnalyser/index.html)
-   [MetABEL](../packages/MetABEL/index.html)
-   [metaBMA](../packages/metaBMA/index.html)
-   [metacart](../packages/metacart/index.html)
-   [metacor](../packages/metacor/index.html)
-   [metafor](../packages/metafor/index.html) (core)
-   [metaforest](../packages/metaforest/index.html)
-   [metafuse](../packages/metafuse/index.html)
-   [metagear](../packages/metagear/index.html)
-   [metagen](../packages/metagen/index.html)
-   [metagen](../packages/metagen/index.html)
-   [MetaIntegrator](../packages/MetaIntegrator/index.html)
-   [metaLik](../packages/metaLik/index.html)
-   [metaMA](../packages/metaMA/index.html)
-   [metamisc](../packages/metamisc/index.html)
-   [metansue](../packages/metansue/index.html)
-   [metap](../packages/metap/index.html)
-   [MetaPath](../packages/MetaPath/index.html)
-   [MetaPCA](../packages/MetaPCA/index.html)
-   [metaplotr](../packages/metaplotr/index.html)
-   [metaplus](../packages/metaplus/index.html)
-   [MetaQC](../packages/MetaQC/index.html)
-   [metaRNASeq](../packages/metaRNASeq/index.html)
-   [metaSEM](../packages/metaSEM/index.html)
-   [metasens](../packages/metasens/index.html)
-   [MetaSKAT](../packages/MetaSKAT/index.html)
-   [metatest](../packages/metatest/index.html)
-   [Metatron](../packages/Metatron/index.html)
-   [metavcov](../packages/metavcov/index.html)
-   [metaviz](../packages/metaviz/index.html)
-   [mmeta](../packages/mmeta/index.html)
-   [MultiMeta](../packages/MultiMeta/index.html)
-   [mvmeta](../packages/mvmeta/index.html)
-   [mvtmeta](../packages/mvtmeta/index.html)
-   [netmeta](../packages/netmeta/index.html)
-   [nmaINLA](../packages/nmaINLA/index.html)
-   [nmathresh](../packages/nmathresh/index.html)
-   [pcnetmeta](../packages/pcnetmeta/index.html)
-   [pimeta](../packages/pimeta/index.html)
-   [psychmeta](../packages/psychmeta/index.html)
-   [psychometric](../packages/psychometric/index.html)
-   [PubBias](../packages/PubBias/index.html)
-   [RandMeta](../packages/RandMeta/index.html)
-   [ratesci](../packages/ratesci/index.html)
-   [RBesT](../packages/RBesT/index.html)
-   [RcmdrPlugin.EZR](../packages/RcmdrPlugin.EZR/index.html)
-   [RcmdrPlugin.RMTCJags](../packages/RcmdrPlugin.RMTCJags/index.html)
-   [revtools](../packages/revtools/index.html)
-   [rma.exact](../packages/rma.exact/index.html)
-   [rmeta](../packages/rmeta/index.html)
-   [robumeta](../packages/robumeta/index.html)
-   [SAMURAI](../packages/SAMURAI/index.html)
-   [SCMA](../packages/SCMA/index.html)
-   [selectMeta](../packages/selectMeta/index.html)
-   [seqMeta](../packages/seqMeta/index.html)
-   [surrosurv](../packages/surrosurv/index.html)
-   [TFisher](../packages/TFisher/index.html)
-   [weightr](../packages/weightr/index.html)
-   [xmeta](../packages/xmeta/index.html)

### Related links:

-   CRAN Task View: [Genetics](Genetics.html)
-   CRAN Task View: [ClinicalTrials](ClinicalTrials.html)
