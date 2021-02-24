## CRAN Task View: Design of Experiments (DoE) & Analysis of Experimental Data

  ----------------- ------------------------------------------------------
  **Maintainer:**   Ulrike Groemping
  **Contact:**      groemping at bht-berlin.de
  **Version:**      2018-02-17
  **URL:**          <https://CRAN.R-project.org/view=ExperimentalDesign>
  ----------------- ------------------------------------------------------

<div>

This task view collects information on R packages for experimental
design and analysis of data from experiments. With a strong increase in
the number of relevant packages, packages that focus on analysis only
and do not make relevant contributions for design creation are no longer
added to this task view. Please feel free to suggest enhancements, and
please send information on new packages or major package updates if you
think they belong here. Contact details are given on my [Web
page](http://prof.beuth-hochschule.de/groemping/DoE) .

Experimental design is applied in many areas, and methods have been
tailored to the needs of various fields. This task view starts out with
a section on the historically earliest application area, agricultural
experimentation. Subsequently, it covers the most general packages,
continues with specific sections on industrial experimentation, computer
experiments, and experimentation in the clinical trials contexts (this
section is going to be removed eventually; experimental design packages
for clinical trials will be integrated into the clinical trials task
view), and closes with a section on various special experimental design
packages that have been developed for other specific purposes. Of
course, the division into fields is not always clear-cut, and some
packages from the more specialized sections can also be applied in
general contexts.\
You may also notice that my own experience is mainly from industrial
experimentation (in a broad sense), which may explain a somewhat biased
view on things.

### Experimental designs for agricultural and plant breeding experiments

Package [agricolae](../packages/agricolae/index.html) is by far the
most-used package from this task view (status: October 2017). It offers
extensive functionality on experimental design especially for
agricultural and plant breeding experiments, which can also be useful
for other purposes. It supports *planning* of lattice designs, factorial
designs, randomized complete block designs, completely randomized
designs, (Graeco-)Latin square designs, balanced incomplete block
designs and alpha designs. There are also various *analysis* facilities
for experimental data, e.g. treatment comparison procedures and several
non-parametric tests, but also some quite specialized possibilities for
specific types of experiments. Package
[desplot](../packages/desplot/index.html) is made for plotting the
layout of agricultural experiments. Package
[agridat](../packages/agridat/index.html) offers a large repository of
useful agricultural data sets.

### Experimental designs for general purposes

There are a few packages for creating and analyzing experimental designs
for general purposes: First of all, the standard (generalized) linear
model functions in the base package stats are of course very important
for analyzing data from designed experiments (especially functions
`lm()`, `aov()` and the methods and functions for the resulting linear
model objects). These are concisely explained in Kuhnert and Venables
(2005, p. 109 ff.); Vikneswaran (2005) points out specific usages for
experimental design (using function `contrasts()`, multiple comparison
functions and some convenience functions like `model.tables()`,
`replications()` and `plot.design()`). Lawson (2014) is a good
introductory textbook on experimental design in R, which gives many
example applications. Lalanne (2012) provides an R companion to the
well-known book by Montgomery (2005); he so far covers approximately the
first ten chapters; he does not include R\'s design generation
facilities, but mainly discusses the analysis of existing designs.
Package [GAD](../packages/GAD/index.html) handles general balanced
analysis of variance models with fixed and/or random effects and also
nested effects (the latter can only be random); they quote Underwood
(1997) for this work. The package is quite valuable, as many users have
difficulties with using the R packages for handling random or mixed
effects. Package [granova](../packages/granova/index.html) offers some
interesting non-standard graphical representations for results of
simply-structured experiments (one-way and two-way layouts, paired
data), package [ez](../packages/ez/index.html) aims at supporting
intuitive analysis and visualization of factorial experiments based on
package \"ggplot2\".

-   Package [AlgDesign](../packages/AlgDesign/index.html) creates full
    factorial designs with or without additional quantitative variables,
    creates mixture designs (i.e., designs where the levels of factors
    sum to 1=100%; lattice designs are created only) and creates D-, A-,
    or I-optimal designs exactly or approximately, possibly with
    blocking, using the Federov (1972) algorithm.
-   Package [skpr](../packages/skpr/index.html) also provides optimal
    designs (D, I, A, Alias, G, T, or E optimal); a selection of the
    optimality criteria can also be used for the stepwise creation of
    split-plot designs. The package can also assess the power of designs
    and display diagnostic plots. At the moment (October 2017), the
    algorithms used are not yet documented.
-   Package [OptimalDesign](../packages/OptimalDesign/index.html)
    likewise calculates unblocked D-, A-, or I-optimal designs (they use
    \"IV-optimal\" instead of \"I-optimal\") exactly or approximately,
    treating quantitative variables only, including mixture designs;
    this package uses different algorithms (e.g. Atkinson, Donev and
    Tobias 2007, Harman and Filova 2014), some of which rely on the
    availability of the gurobi software ( <http://www.gurobi.com/> ,
    free for academics and academic institutions) and its accompanying R
    package \"gurobi\" (not on CRAN).
-   Package [ICAOD](../packages/ICAOD/index.html) implements the
    \"Imperialist Competitive Algorithm for Optimal Designs\" for
    nonlinear models according to Masoudi, Holling and Wong (2016).
    Package [LDOD](../packages/LDOD/index.html) implements locally
    D-optimal designs for some nonlinear and generalized linear models,
    package [designGLMM](../packages/designGLMM/index.html) locally
    optimal designs for completely randomized or blocked Poisson models
    and package [PopED](../packages/PopED/index.html) provides optimal
    designs for nonlinear mixed effect models.
-   There are various further packages that deal with optimal designs of
    different types: Package [rodd](../packages/rodd/index.html)
    provides T-optimal designs, also called optimal discriminating
    designs (Dette, Melas and Shpilev 2013, Dette, Melas and Guchenko
    2014), Package [acebayes](../packages/acebayes/index.html)
    calculates optimal Bayesian designs using an approximate coordinate
    exchange algorithm, package [OBsMD](../packages/OBsMD/index.html)
    provides \"Objective Bayesian Model Discrimination in Follow-Up
    Designs\" according to Consonni and Deldossi (2015). Further optimal
    design packages for very specific purposes are listed at the end of
    this view.
-   Package [conf.design](../packages/conf.design/index.html) allows to
    create a design with certain interaction effects confounded with
    blocks (function `conf.design()`) and allows to combine existing
    designs in several ways (e.g., useful for Taguchi\'s inner and outer
    array designs in industrial experimentation).
-   Package [ibd](../packages/ibd/index.html) creates and analyses
    incomplete block designs. Packages
    [PGM2](../packages/PGM2/index.html),
    [RPPairwiseDesign](../packages/RPPairwiseDesign/index.html) and
    [CombinS](../packages/CombinS/index.html) all produce designs
    related to (resolvable) (partially) balanced incomplete block
    designs. Package [PBIBD](../packages/PBIBD/index.html) also provides
    experts with some series of partially balanced incomplete block
    designs.
-   Package [crossdes](../packages/crossdes/index.html) creates and
    analyses cross-over designs of various types (including latin
    squares, mutually orthogonal latin squares and Youden squares) that
    can for example be used in sensometrics. Package
    [Crossover](../packages/Crossover/index.html) also provides
    crossover designs; it offers designs from the literature and
    algorithmic designs, makes use of the functionality in
    [crossdes](../packages/crossdes/index.html) and in addition provides
    a GUI.
-   Package [DoE.base](../packages/DoE.base/index.html) provides full
    factorial designs with or without blocking (function `fac.design`)
    and orthogonal arrays (function `oa.design`) for main effects
    experiments (those listed by Kuhfeld 2009 up to 144 runs, plus a few
    additional ones). There is also some functionality for assessing the
    quality of orthogonal arrays, related to Groemping and Xu (2014) and
    Groemping (2017), and some analysis functionality with half-normal
    effects plots in quite general form (Groemping 2015).\
    Package [DoE.base](../packages/DoE.base/index.html) also forms the
    basis of a suite of related packages: together with
    [FrF2](../packages/FrF2/index.html) (cf. below) and
    [DoE.wrapper](../packages/DoE.wrapper/index.html), it provides the
    work horse of the GUI package
    [RcmdrPlugin.DoE](../packages/RcmdrPlugin.DoE/index.html) (beta
    version; tutorial available in Groemping 2011), which integrates
    design of experiments functionality into the R-Commander (package
    \"Rcmdr\", Fox 2005) for the benefit of those R users who cannot or
    do not want to do command line programming. The role of package
    [DoE.wrapper](../packages/DoE.wrapper/index.html) in that suite is
    to wrap functionality from other packages into the input and output
    structure of the package suite (so far for response surface designs
    with package [rsm](../packages/rsm/index.html) (cf. also below),
    design of computer experiments with packages
    [lhs](../packages/lhs/index.html) and
    [DiceDesign](../packages/DiceDesign/index.html) (cf. also below),
    and , and D-optimal designs with package
    [AlgDesign](../packages/AlgDesign/index.html) (cf. also above).
-   Package [DoE.MIParray](../packages/DoE.MIParray/index.html) creates
    optimized orthogonal arrays (or even supersaturated arrays) for
    factorial experiments. Arrays created with this package can be used
    as input to function oa.design of package
    [DoE.base](../packages/DoE.base/index.html). Note, however, that the
    package is only useful in combination with at least one of the
    commercial optimizers
    [Gurobi](http://www.gurobi.com/products/modeling-languages/r)
    (R-package gurobi delivered with the software) or
    [Mosek](https://www.mosek.com/documentation/) (R-package Rmosek
    downloadable from the vendor (an outdated version is on CRAN)).
-   Package [dae](../packages/dae/index.html) provides various utility
    functions around experimental design and R factors, e.g. a
    randomization routine that can handle various nested structures
    (according to Bailey 1981) and functions for combining several
    factors into one or dividing one factor into several factors.
    Furthermore, the package provides features for post-processing
    objects returned by the `aov()` function, e.g. extraction of Yates
    effects for 2-level experiments.
-   Package [daewr](../packages/daewr/index.html) accompanies the book
    *Design and Analysis of Experiments with R* by Lawson (2014) and
    does not only provide data sets from the book but also some
    standalone functionality that is not available elsewhere in R, e.g.
    definitive screening designs.
-   Package [OPDOE](../packages/OPDOE/index.html) accompanies the book
    *Optimal Experimental Design with R* by Rasch et al. (2011). It has
    some interesting sample size estimation functionality, but is almost
    unusable without the book (the first edition of which I would not
    recommend buying).
-   Package [blockTools](../packages/blockTools/index.html) assigns
    units to blocks in order to end up with homogeneous sets of blocks
    in case of too small block sizes; package
    [blocksdesign](../packages/blocksdesign/index.html) permits the
    creation of nested block structures.
-   There are several packages for determining sample sizes in
    experimental contexts, some of them quite general, others very
    specialized. All of these are mentioned here: packages
    [powerAnalysis](../packages/powerAnalysis/index.html),
    [powerbydesign](../packages/powerbydesign/index.html) and
    [easypower](../packages/easypower/index.html) deal with estimating
    the power, sample size and/or effect size for factorial experiments.
    Package [JMdesign](../packages/JMdesign/index.html) deals with the
    power for the special situation of jointly modeling longitudinal and
    survival data, package [PwrGSD](../packages/PwrGSD/index.html) with
    the power for group sequential designs, package
    [powerGWASinteraction](../packages/powerGWASinteraction/index.html)
    with the power for interactions in genome wide association studies,
    package [ssizeRNA](../packages/ssizeRNA/index.html) with sample size
    for RNA sequencing experiments, and package
    [ssize.fdr](../packages/ssize.fdr/index.html) for sample sizes in
    microarray experiments (requesting a certain power while limiting
    the false discovery rate).

### Experimental designs for industrial experiments

Some further packages especially handle designs for industrial
experiments that are often highly fractionated, intentionally confounded
and have few extra degrees of freedom for error.

Fractional factorial 2-level designs are particularly important in
industrial experimentation.

-   Package [FrF2](../packages/FrF2/index.html) (Groemping 2014) is the
    most comprehensive R package for their creation. It generates
    regular Fractional Factorial designs for factors with 2 levels as
    well as Plackett-Burman type screening designs. Regular fractional
    factorials default to maximum resolution minimum aberration designs
    and can be customized in various ways, supported by an incorporated
    catalogue of designs (including the designs catalogued by Chen, Sun
    and Wu 1993, and further larger designs catalogued in Block and Mee
    2005 and Xu 2009; the additional package
    [FrF2.catlg128](../packages/FrF2.catlg128/index.html) provides a
    very large complete catalogue for resolution IV 128 run designs with
    up to 23 factors for special purposes). Analysis-wise,
    [FrF2](../packages/FrF2/index.html) provides simple graphical
    analysis tools (normal and half-normal effects plots (modified from
    [BsMD](../packages/BsMD/index.html), cf. below), main effects plots
    and interaction plot matrices similar to those in Minitab software,
    and a cube plot for the combinations of three factors). It can also
    show the alias structure for regular fractional factorials of
    2-level factors, regardless whether they have been created with the
    package or not.\
    Fractional factorial 2-level plans can also be created by other R
    packages, namely [BHH2](../packages/BHH2/index.html) and
    [qualityTools](../packages/qualityTools/index.html) (but do not use
    function pbDesign from version 1.54 of that package!), or with a
    little bit more complication by packages
    [conf.design](../packages/conf.design/index.html) or
    [AlgDesign](../packages/AlgDesign/index.html). Package
    [ALTopt](../packages/ALTopt/index.html) provides optimal designs for
    accelerated life testing.
-   Package [BHH2](../packages/BHH2/index.html) accompanies the 2nd
    edition of the book by Box, Hunter and Hunter and provides various
    of its data sets. It can generate full and fractional factorial
    two-level-designs from a number of factors and a list of defining
    relations (function `ffDesMatrix()`, less comfortable than package
    FrF2). It also provides several functions for analyzing data from
    2-level factorial experiments: The function anovaPlot assesses
    effect sizes relative to residuals, and the function `lambdaPlot()`
    assesses the effect of Box-Cox transformations on statistical
    significance of effects.
-   [BsMD](../packages/BsMD/index.html) provides Bayesian charts as
    proposed by Box and Meyer (1986) as well as effects plots (normal,
    half-normal and Lenth) for assessing which effects are active in a
    fractional factorial experiment with 2-level factors; package
    [OBsMD](../packages/OBsMD/index.html) provides the functionality for
    follow- up experiments for resolving ambiguities after applying the
    Bayesian analysis of package [BsMD](../packages/BsMD/index.html).
-   Package [unrepx](../packages/unrepx/index.html) provides a battery
    of methods for the assessment of effect estimates from unreplicated
    factorial experiments, including many of the effects plots also
    present in other packages, but also further possibilities.
-   The small package [FMC](../packages/FMC/index.html) provides
    factorial designs with minimal number of level changes; the package
    does not take any measures to account for the statistical
    implications this may imply. Thus, using this package must be
    considered very risky for many experimental situations, because in
    many experiments some variability is caused by level changes. For
    such situations (and they are the rule rather than the exception),
    minimizing the level changes without taking precautions in the
    analysis will yield misleading results.
-   Package [pid](../packages/pid/index.html) accompanies an online book
    by Dunn (2010-2016) and also makes heavy use of the Box, Hunter and
    Hunter book; it provides various data sets, which are mostly from
    fractional factorial 2-level designs.

Apart from tools for planning and analysing factorial designs, R also
offers support for response surface optimization for quantitative
factors (cf. e.g. Myers and Montgomery 1995):

-   Package [rsm](../packages/rsm/index.html) supports sequential
    optimization with first order and second order response surface
    models (central composite or Box-Behnken designs), offering
    optimization approaches like steepest ascent and visualization of
    the response function for linear model objects. Also, coding for
    response surface investigations is facilitated.
-   Package [DoE.wrapper](../packages/DoE.wrapper/index.html) enhances
    design creation from package [rsm](../packages/rsm/index.html) with
    the possibilities of automatically choosing the cube portion of
    central composite designs and of augmenting an existing (fractional)
    factorial 2-level design with a star portion.
-   The small package [rsurface](../packages/rsurface/index.html)
    provides rotatable central composite designs for which the user
    specifies the minimum and maximum of the experimental variables
    instead of the corner points of the cube.
-   The small package [minimalRSD](../packages/minimalRSD/index.html)
    provides central composite and Box-Behnken designs with minimal
    number of level changes; the package does not take any measures to
    account for the statistical implications this may imply. Thus, using
    this package must be considered very risky for many experimental
    situations, because in many experiments some variability is caused
    by level changes. For such situations (and they are the rule rather
    than the exception), minimizing the level changes without taking
    precautions in the analysis will yield misleading results.
-   Package [OptimaRegion](../packages/OptimaRegion/index.html) provides
    functionality for inspecting the optimal region of a response
    surface for quadratic polynomials and thin-plate spline models and
    can compute a confidence interval for the distance between two
    optima.
-   Package [Vdgraph](../packages/Vdgraph/index.html) implements a
    variance dispersion graph (Vining 1993) for response surface designs
    created by package [rsm](../packages/rsm/index.html). Packages
    [VdgRsm](../packages/VdgRsm/index.html) and
    [vdg](../packages/vdg/index.html) provide similar functionality with
    more variety.
-   Package [qualityTools](../packages/qualityTools/index.html) can also
    create central composite designs and can visualize response
    surfaces.
-   Package [EngrExpt](../packages/EngrExpt/index.html) provides a
    collection of data sets from the book *Introductory Statistics for
    Engineering Experimentation* by Nelson, Coffin and Copeland (2003).

In some industries, mixtures of ingredients are important; these require
special designs, because the quantitative factors have a fixed total.
Mixture designs are handled by packages
[AlgDesign](../packages/AlgDesign/index.html) (function `gen.mixture`,
lattice designs), [qualityTools](../packages/qualityTools/index.html)
(function `mixDesign`, lattice designs and simplex centroid designs),
and [mixexp](../packages/mixexp/index.html) (several small functions for
simplex centroid, simplex lattice and extreme vertices designs as well
as for plotting).

Occasionally, supersaturated designs can be useful. The two small
packages [mkssd](../packages/mkssd/index.html) and
[mxkssd](../packages/mxkssd/index.html) provide fixed level and mixed
level k-circulant supersaturated designs. The aforementioned package
[DoE.MIParray](../packages/DoE.MIParray/index.html) can also provide
(small!) supersaturated arrays (by choosing resolution II), but requires
the presence of at least one of the commercial optimizers
[Gurobi](http://www.gurobi.com/products/modeling-languages/r) or
[Mosek](https://www.mosek.com/documentation/) .

### Experimental designs for computer experiments

Computer experiments with quantitative factors require special types of
experimental designs: it is often possible to include many different
levels of the factors, and replication will usually not be beneficial.
Also, the experimental region is often too large to assume that a linear
or quadratic model adequately represents the phenomenon under
investigation. Consequently, it is desirable to fill the experimental
space with points as well as possible (space-filling designs) in such a
way that each run provides additional information even if some factors
turn out to be irrelevant. The [lhs](../packages/lhs/index.html) package
provides latin hypercube designs for this purpose. Furthermore, the
package provides ways to analyse such computer experiments with emphasis
on what follow-up experiments to conduct. Another package with similar
orientation is the [DiceDesign](../packages/DiceDesign/index.html)
package, which adds further ways to construct space-filling designs and
some measures to assess the quality of designs for computer experiments.
The package [DiceKriging](../packages/DiceKriging/index.html) provides
the kriging methodology which is often used for creating meta models
from computer experiments, the package
[DiceEval](../packages/DiceEval/index.html) creates and evaluates meta
models (among others Kriging ones).

Package [MaxPro](../packages/MaxPro/index.html) provides maximum
projection designs as introduced by Joseph, Gul and Ba(2015). Package
[SLHD](../packages/SLHD/index.html) provides optimal sliced latin
hypercube designs according to Ba et al. (2015), package
[sFFLHD](../packages/sFFLHD/index.html) provides sliced full
factorial-based latin hypercube designs according to Duan et al. (2017).
Package [simrel](../packages/simrel/index.html) allows creation of
designs for computer experiments according to the Multi-level binary
replacement (MBR) strategy by Martens et al. (2010). Package
[minimaxdesign](../packages/minimaxdesign/index.html) provides minimax
designs and minimax projection designs according to Mak and Joseph
(2016).

Package [tgp](../packages/tgp/index.html) is another package dedicated
to planning and analysing computer experiments. Here, emphasis is on
Bayesian methods. The package can for example be used with various kinds
of (surrogate) models for sequential optimization, e.g. with an expected
improvement criterion for optimizing a noisy blackbox target function.
Packages [plgp](../packages/plgp/index.html) and
[dynaTree](../packages/dynaTree/index.html) enhance the functionality
offered by [tgp](../packages/tgp/index.html) with particle learning
facilities and learning for dynamic regression trees.

Package [BatchExperiments](../packages/BatchExperiments/index.html) is
also designed for computer experiments, in this case specifically for
experiments with algorithms to be run under different scenarios. The
package is described in a technical report by Bischl et al. (2012).

### Experimental designs for clinical trials

This task view only covers specific design of experiments packages
(which will eventually also be removed here); there may be some grey
areas. Please, also consult the [ClinicalTrials](ClinicalTrials.html)
task view.

-   Package [experiment](../packages/experiment/index.html) contains
    tools for clinical experiments, e.g., a randomization tool, and it
    provides a few special analysis options for clinical trials.
-   Package [ThreeArmedTrials](../packages/ThreeArmedTrials/index.html)
    provides design and analysis tools for three-armed superiority or
    non-inferiority trials. Beside the standard functionality, the
    package includes the negative Binomial response situation discussed
    in Muetze et al. (2016).
-   Package [gsDesign](../packages/gsDesign/index.html) implements group
    sequential designs, package
    [GroupSeq](../packages/GroupSeq/index.html) gives a GUI for
    probability spending in such designs, package
    [OptGS](../packages/OptGS/index.html) near-optimal balanced group
    sequential designs. Package
    [gsbDesign](../packages/gsbDesign/index.html) evaluates operating
    characteristics for group sequential Bayesian designs. Package
    [gset](../packages/gset/index.html) handles group sequential
    equivalence testing. Package
    [seqDesign](../packages/seqDesign/index.html) handles group
    sequential two-stage treatment efficacy trials with time-to-event
    endpoints.
-   Package [binseqtest](../packages/binseqtest/index.html) handles
    sequential single arm binary response trials.
-   Package [asd](../packages/asd/index.html) implements adaptive
    seamless designs (see e.g. Parsons et al. 2012).
-   Package [OptInterim](../packages/OptInterim/index.html) is for two-
    and three-stage designs for longterm binary endpoints.
-   Packages [bcrm](../packages/bcrm/index.html) offers Bayesian CRM
    designs.
-   Package [MAMS](../packages/MAMS/index.html) offers designs for
    multi-arm multi stage studies,
    [BayesMAMS](../packages/BayesMAMS/index.html) provides a Bayesian
    sample size calculations for these.
-   Package [BOIN](../packages/BOIN/index.html) provides Bayesian
    optimal interval designs, which are used in phase I clinical trials
    for finding the maximum tolerated dose.
-   The [DoseFinding](../packages/DoseFinding/index.html) package
    provides functions for the design and analysis of dose-finding
    experiments (for example pharmaceutical Phase II clinical trials);
    it combines the facilities of the \"MCPMod\" package (maintenance
    discontinued; described in Bornkamp, Pinheiro and Bretz 2009) with a
    special type of optimal designs for dose finding situations
    (MED-optimal designs, or D-optimal designs, or a mixture of both;
    cf., Dette et al. 2008).
-   Package [VNM](../packages/VNM/index.html) provides multi-objective
    optimal designs for simultaneously optimizing inference about the
    shape of the dose-response curve, ED50 and minimum effective dose
    (MED) for certain classes of logistic models.
-   Package [TEQR](../packages/TEQR/index.html) provides toxicity
    equivalence range designs (Blanchard and Longmate 2010) for phase I
    clinical trials, package
    [pipe.design](../packages/pipe.design/index.html) so-called *product
    of independent beta probabilities dose escalation* (PIPE) designs
    for phase I. Package [dfpk](../packages/dfpk/index.html) implements
    a Bayesian dose-finding design using pharmacokinetics for phase I
    trials. Package [dfcrm](../packages/dfcrm/index.html) provides
    designs for classical or TITE continual reassessment trials in
    phase I.
-   Packages [dfcomb](../packages/dfcomb/index.html) and
    [dfmta](../packages/dfmta/index.html) provide phase I/II adaptive
    dose-finding designs for combination studies or single-agent
    molecularly targeted agent, respectively.
-   Packages [ph2bayes](../packages/ph2bayes/index.html) and
    [ph2bye](../packages/ph2bye/index.html) are concerned with Bayesian
    single arm phase II trials.
-   Package [sp23design](../packages/sp23design/index.html) claims to
    offer seamless integration of phase II to III.

### Experimental designs for special purposes

Various further packages handle special situations in experimental
design:

-   Package [desirability](../packages/desirability/index.html) provides
    ways to combine several target criteria into a desirability function
    in order to simplify multi-criteria analysis; desirabilities are
    also offered as part of package
    [qualityTools](../packages/qualityTools/index.html).
-   [osDesign](../packages/osDesign/index.html) designs studies nested
    in observational studies,
    [designmatch](../packages/designmatch/index.html) can also be useful
    for this purpose.
-   [qtlDesign](../packages/qtlDesign/index.html) is for quantitative
    trait locus designs,
-   [toxtestD](../packages/toxtestD/index.html) creates optimal designs
    for binary toxicity tests,
-   [hiPOD](../packages/hiPOD/index.html) provides optimal designs for
    pooled next generation sequencing experiments,
-   [designGG](../packages/designGG/index.html) creates optimal designs
    for genetical genomics experiments (see Li et al. 2009),
-   packages [optbdmaeAT](../packages/optbdmaeAT/index.html),
    [optrcdmaeAT](../packages/optrcdmaeAT/index.html) and
    [soptdmaeA](../packages/soptdmaeA/index.html) provide optimal block
    designs, optimal row-column designs, and sequential optimal or
    near-optimal block or row-column designs for two-colour cDNA
    microarray experiments, with optimality according to an A-, MV-, D-
    or E-criterion.
-   Package [docopulae](../packages/docopulae/index.html) implements
    optimal designs for copula models according to Perrone and Mueller
    (2016),
-   [optDesignSlopeInt](../packages/optDesignSlopeInt/index.html)
    provides an optimal design for the estimation of the ratio of slope
    to intercept, and
-   Packages [edesign](../packages/edesign/index.html) and
    [MBHdesign](../packages/MBHdesign/index.html) provide spatially
    balanced designs, allowing the inclusion of prespecified (legacy)
    sites. The more elaborate package
    [geospt](../packages/geospt/index.html) allows to optimize spatial
    networks of sampling points (see e.g. Santacruz, Rubiano and Melo
    2014).
-   Package [SensoMineR](../packages/SensoMineR/index.html) contains
    special designs for sensometric studies, e.g., for the triangle
    test.
-   Package [choiceDes](../packages/choiceDes/index.html) creates choice
    designs with emphasis on discrete choice models and MaxDiff
    functionality; it is based on optimal designs. Package
    [idefix](../packages/idefix/index.html) provides D-efficient designs
    for discrete choice experiments based on the multinomial logit
    model, and individually adapted designs for the mixed multinomial
    logit model (Crabbe et al. 2014). Package
    [support.CEs](../packages/support.CEs/index.html) provides tools for
    creating stated choice designs for market research investigations,
    based on orthogonal arrays.
-   Package [odr](../packages/odr/index.html) creates optimal designs
    for cluster randomized trials under condition- and unit-specific
    cost structures.
-   Package [bioOED](../packages/bioOED/index.html) offers sensitivity
    analysis and optimal design for microbial inactivation.

### Key references for packages in this task view

-   Atkinson, A.C. and Donev, A.N. (1992). *Optimum Experimental
    Designs.* Oxford: Clarendon Press.
-   Atkinson, A.C., Donev, A.N. and Tobias, R.D. (2007). *Optimum
    Experimental Designs, with SAS.* Oxford University Press, Oxford.
-   Ba,S., Brenneman, W.A. and Myers, W.R. (2015). Optimal Sliced Latin
    Hypercube Designs. *Technometrics* **57** 479-487.
-   Bailey, R.A. (1981). A unified approach to design of experiments.
    *Journal of the Royal Statistical Society, Series A* **144**
    214-223.
-   Ball, R.D. (2005). Experimental Designs for Reliable Detection of
    Linkage Disequilibrium in Unstructured Random Population Association
    Studies. *Genetics* **170** 859-873.
-   Bischl, B., Lang, M., Mersmann, O., Rahnenfuehrer, J. and Weihs, C.
    (2012). [Computing on high performance clusters with R: Packages
    BatchJobs and
    BatchExperiments](http://www1.beuth-hochschule.de/FB_II/reports/Report-2011-004.pdf)
    . *Technical Report 1/2012* , TU Dortmund, Germany.
-   Blanchard, M.S. and Longmate, J.A. (2010). Toxicity equivalence
    range design (TEQR): A practical Phase I design. *Contemporary
    Clinical Trials* doi:10.1016/j.cct.2010.09.011.
-   Block, R. and Mee, R. (2005). Resolution IV Designs with 128 Runs.
    *Journal of Quality Technology* **37** 282-293.
-   Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). [MCPMod: An R
    Package for the Design and Analysis of Dose-Finding
    Studies](http://www.jstatsoft.org/v29/i07/paper) . *Journal of
    Statistical Software* **29** (7) 1-23.
-   Box G. E. P, Hunter, W. C. and Hunter, J. S. (2005). *Statistics for
    Experimenters* (2nd edition). New York: Wiley.
-   Box, G. E. P and R. D. Meyer (1986). An Analysis for Unreplicated
    Fractional Factorials. *Technometrics* **28** 11-18.
-   Box, G. E. P and R. D. Meyer (1993). Finding the Active Factors in
    Fractionated Screening Experiments. *Journal of Quality Technology*
    **25** 94-105.
-   Chasalow, S., Brand, R. (1995). Generation of Simplex Lattice
    Points. *Journal of the Royal Statistical Society, Series C* **44**
    534-545.
-   Chen, J., Sun, D.X. and Wu, C.F.J. (1993). A catalogue of 2-level
    and 3-level orthogonal arrays. *International Statistical Review*
    **61** 131-145.
-   Consonni, G. and Deldossi, L. (2015), Objective Bayesian model
    discrimination in follow-up experimental designs DOI
    10.1007/s11749-015-0461-3. TEST.
-   Collings, B. J. (1989). Quick Confounding. *Technometrics* **31**
    107-110.
-   Cornell, J. (2002). *Experiments with Mixtures* . Third Edition.
    Wiley.
-   Crabbe, M., Akinc, D. and Vandebroek, M. (2014). Fast algorithms to
    generate individualized designs for the mixed logit choice model.
    *Transportation Research Part B: Methodological* **60** , 1-15.
-   Daniel, C. (1959). Use of Half Normal Plots in Interpreting Two
    Level Experiments. *Technometrics* **1** 311-340.
-   Derringer, G. and Suich, R. (1980). Simultaneous Optimization of
    Several Response Variables. *Journal of Quality Technology* **12**
    214-219.
-   Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro, J. C. (2008).
    Optimal Designs for Dose Finding Studies. *Journal of the American
    Statisical Association* **103** 1225-1237.
-   Dette, H., Melas, V.B. and Shpilev, P. (2013). Robust T-optimal
    discriminating designs. *The Annals of Statistics* **41** 1693-1715.
-   Dette H., Melas V.B. and Guchenko R. (2014). Bayesian T-optimal
    discriminating designs. [ArXiv link](http://arxiv.org/abs/1412.2548)
    .
-   Duan, W., Ankenman, B.E. Sanchez, S.M. and Sanchez, P.J. (2017).
    Sliced Full Factorial-Based Latin Hypercube Designs as a Framework
    for a Batch Sequential Design Algorithm. *Technometrics* **59** ,
    11-22.
-   Dunn, K. (2010-2016). *Process Improvement Using Data* . [Online
    book.](http://learnche.org/pid)
-   Federov, V.V. (1972). *Theory of Optimal Experiments.* Academic
    Press, New York.
-   Fox, J. (2005). [The R Commander: A Basic-Statistics Graphical User
    Interface to R](http://www.jstatsoft.org/v14/i09/paper) . *Journal
    of Statistical Software* **14** (9) 1-42.
-   Gramacy, R.B. (2007). [tgp: An R Package for Bayesian Nonstationary,
    Semiparametric Nonlinear Regression and Design by Treed Gaussian
    Process Models](http://www.jstatsoft.org/v19/i09/paper) . *Journal
    of Statistical Software* **19** (9) 1-46.
-   Groemping, U. (2011). [Tutorial for designing experiments using the
    R package
    RcmdrPlugin.DoE](http://www1.beuth-hochschule.de/FB_II/reports/Report-2011-004.pdf)
    . [Reports in Mathematics, Physics and
    Chemistry](http://www1.beuth-hochschule.de/FB_II/reports/welcome.htm)
    , Department II, Beuth University of Applied Sciences Berlin.
-   Groemping, U. (2014). R Package FrF2 for Creating and Analysing
    Fractional Factorial 2-Level Designs. *Journal of Statistical
    Software* **56** (1) 1-56.
-   Groemping, U. (2015). Augmented Half Normal Effects Plots in the
    Presence of a Few Error Degrees of Freedom. *Quality and Reliability
    Engineering International* **31** , 1185-1196. DOI:
    10.1002/qre.1842.
-   Groemping, U. (2017). Frequency Tables for the Coding Invariant
    Quality Assessment of Factorial Designs. *IISE Transactions* **49**
    , 505-517.
-   Groemping, U. and Xu, H. (2014). Generalized resolution for
    orthogonal arrays. *The Annals of Statistics* **42** 918-939.
-   Harman R., Filova L. (2014): Computing efficient exact designs of
    experiments using integer quadratic programming, *Computational
    Statistics and Data Analysis* **71** 1159-1167
-   Hoaglin D., Mosteller F. and Tukey J. (eds., 1991). *Fundamentals of
    Exploratory Analysis of Variance* . Wiley, New York.
-   Jones, B. and Kenward, M.G. (1989). *Design and Analysis of
    Cross-Over Trials* . Chapman and Hall, London.
-   Johnson, M.E., Moore L.M. and Ylvisaker D. (1990). Minimax and
    maximin distance designs. *Journal of Statistical Planning and
    Inference* **26** 131-148.
-   Joseph, V. R., Gul, E., and Ba, S. (2015). Maximum Projection
    Designs for Computer Experiments. *Biometrika* **102** 371-380.
-   Kuhfeld, W. (2009). Orthogonal arrays. Website courtesy of SAS
    Institute Inc., accessed August 4th 2010. URL
    <http://support.sas.com/techsup/technote/ts723.html> .
-   Kuhnert, P. and Venables, B. (2005) *An Introduction to R: Software
    for Statistical Modelling & Computing* . URL
    [http://CRAN.R-project.org/doc/contrib/Kuhnert+Venables-R_Course_Notes.zip](../../doc/contrib/Kuhnert+Venables-R_Course_Notes.zip)
    . (PDF document (about 360 pages) of lecture notes in combination
    with the data sets and R scripts)
-   Kunert, J. (1998). Sensory Experiments as Crossover Studies. *Food
    Quality and Preference* **9** 243-253.
-   Lalanne, C. (2012). R Companion to Montgomerys Design and Analysis
    of Experiments. Manuscript, downloadable at URL
    <http://www.aliquote.org/articles/tech/dae/dae.pdf> . (The file
    accompanies the book by Montgomery 2005 (cf. below).)
-   Lawson, J. (2014). *Design and Analysis of Experiments with R.*
    Chapman and Hall/CRC, Boca Raton.
-   Lenth, R.V. (1989). Quick and Easy Analysis of Unreplicated
    Factorials. *Technometrics* **31** 469-473.
-   Lenth, R.V. (2009). [Response-Surface Methods in R, Using
    rsm](http://www.jstatsoft.org/v32/i07/paper) . *Journal of
    Statistical Software* **32** (7) 1-17.
-   Y. Li, M. Swertz, G. Vera, J. Fu, R. Breitling, and R.C. Jansen.
    [designGG](../packages/designGG/index.html): An R-package and Web
    tool for the optimal design of genetical genomics experiments. *BMC
    Bioinformatics* **10** :188
-   Mak, S., and Joseph, V.R. (2016). Minimax designs using clustering.
    *Journal of Computational and Graphical Statistics* . In revision.
-   Martens, H., Mage, I., Tondel, K., Isaeva, J., Hoy, M. and Saebo, S.
    (2010). Multi-level binary replacement (MBR) design for computer
    experiments in high-dimensional nonlinear systems, *J. Chemom.*
    **24** 748-756.
-   Masoudi, E., Holling, H. and Wong, W.-K. (2016). Application of
    imperialist competitive algorithm to find minimax and standardized
    maximin optimal designs. *Computational Statistics and Data
    Analysis* , in press. DOI: 10.1016/j.csda.2016.06.014
-   Mee, R. (2009). *A Comprehensive Guide to Factorial Two-Level
    Experimentation.* Springer, New York.
-   Montgomery, D. C. (2005, 6th ed.). *Design and Analysis of
    Experiments.* Wiley, New York.
-   Muetze,T., Munk, A. and Friede, T. (2016). Design and analysis of
    three-arm trials with negative binomially distributed endpoints.
    *Statistics in Medicine* **35** (4) 505-521.
-   Myers, R. H. and Montgomery, D. C. (1995). *Response Surface
    Methodology: Process and Product Optimization Using Designed
    Experiments.* Wiley, New York.
-   Nelson, P.R., Coffin, M. and Copeland, K.A.F. (2003). *Introductory
    Statistics for Engineering Experimentation.* Academic Press, San
    Diego.
-   Parsons N, Friede T, Todd S, Valdes Marquez E, Chataway J, Nicholas
    R, Stallard N. (2012). An R package for implementing simulations for
    seamless phase II/III clinicals trials using early outcomes for
    treatment selection. *Computational Statistics and Data Analysis*
    **56** , 1150-1160.
-   Perrone, E. and Mueller, W.G. (2016) Optimal designs for copula
    models, *Statistics* **50** (4), 917-929. DOI:
    10.1080/02331888.2015.1111892
-   Plackett, R.L. and Burman, J.P. (1946). The design of optimum
    multifactorial experiments. *Biometrika* **33** 305-325.
-   Rasch, D., Pilz, J., Verdooren, L.R. and Gebhardt, A. (2011).
    *Optimal Experimental Design with R.* Chapman and Hall/CRC.
    (caution, does not live up to its title!)
-   Rosenbaum, P. (1989). Exploratory Plots for Paired Data. *The
    American Statistician* **43** 108-109.
-   Sacks, J., Welch, W.J., Mitchell, T.J. and Wynn, H.P. (1989). Design
    and analysis of computer experiments. *Statistical Science* **4**
    409-435.
-   Santacruz, A., Rubiano, Y., Melo, C., 2014. Evolutionary
    optimization of spatial sampling networks designed for the
    monitoring of soil carbon. In: Hartemink, A., McSweeney, K. (Eds.).
    *Soil Carbon.* Series: Progress in Soil Science. (pp. 77-84).
    Springer, New York.
-   Santner T.J., Williams B.J. and Notz W.I. (2003). *The Design and
    Analysis of Computer Experiments.* Springer, New York.
-   Sen S, Satagopan JM and Churchill GA (2005). Quantitative Trait
    Locus Study Design from an Information Perspective. *Genetics*
    **170** 447-464.
-   Stein, M. (1987). Large Sample Properties of Simulations Using Latin
    Hypercube Sampling. *Technometrics* **29** 143-151.
-   Stocki, R. (2005). A Method to Improve Design Reliability Using
    Optimal Latin Hypercube Sampling. *Computer Assisted Mechanics and
    Engineering Sciences* **12** 87-105.
-   Underwood, A.J. (1997). *Experiments in Ecology: Their Logical
    Design and Interpretation Using Analysis of Variance.* Cambridge
    University Press, Cambridge.
-   Vikneswaran (2005). *An R companion to \"Experimental Design\".* URL
    [http://CRAN.R-project.org/doc/contrib/Vikneswaran-ED_companion.pdf](../../doc/contrib/Vikneswaran-ED_companion.pdf)
    . (The file accompanies the book \"Experimental Design with
    Applications in Management, Engineering and the Sciences\" by Berger
    and Maurer, 2002.)
-   Vining, G. (1993). A Computer Program for Generating Variance
    Dispersion Graphs. *Journal of Quality Technology* **25** 45-58.
    Corrigendum in the same volume, pp. 333-335.
-   Xu, H. (2009). Algorithmic Construction of Efficient Fractional
    Factorial Designs With Large Run Sizes. *Technometrics* **51**
    262-277.
-   Yin, J., Qin, R., Ezzalfani, M., Sargent, D. J., and Mandrekar, S.
    J. (2017). A Bayesian dose-finding design incorporating toxicity
    data from multiple treatment cycles. *Statistics in Medicine* **36**
    , 67-80. doi: 10.1002/sim.7134.

</div>

### CRAN packages:

-   [acebayes](../packages/acebayes/index.html)
-   [agricolae](../packages/agricolae/index.html) (core)
-   [agridat](../packages/agridat/index.html)
-   [AlgDesign](../packages/AlgDesign/index.html) (core)
-   [ALTopt](../packages/ALTopt/index.html)
-   [asd](../packages/asd/index.html)
-   [BatchExperiments](../packages/BatchExperiments/index.html)
-   [BayesMAMS](../packages/BayesMAMS/index.html)
-   [bcrm](../packages/bcrm/index.html)
-   [BHH2](../packages/BHH2/index.html)
-   [binseqtest](../packages/binseqtest/index.html)
-   [bioOED](../packages/bioOED/index.html)
-   [blocksdesign](../packages/blocksdesign/index.html)
-   [blockTools](../packages/blockTools/index.html)
-   [BOIN](../packages/BOIN/index.html)
-   [BsMD](../packages/BsMD/index.html)
-   [choiceDes](../packages/choiceDes/index.html)
-   [CombinS](../packages/CombinS/index.html)
-   [conf.design](../packages/conf.design/index.html) (core)
-   [crossdes](../packages/crossdes/index.html) (core)
-   [Crossover](../packages/Crossover/index.html)
-   [dae](../packages/dae/index.html)
-   [daewr](../packages/daewr/index.html)
-   [designGG](../packages/designGG/index.html)
-   [designGLMM](../packages/designGLMM/index.html)
-   [designmatch](../packages/designmatch/index.html)
-   [desirability](../packages/desirability/index.html)
-   [desplot](../packages/desplot/index.html)
-   [dfcomb](../packages/dfcomb/index.html)
-   [dfcrm](../packages/dfcrm/index.html)
-   [dfmta](../packages/dfmta/index.html)
-   [dfpk](../packages/dfpk/index.html)
-   [DiceDesign](../packages/DiceDesign/index.html)
-   [DiceEval](../packages/DiceEval/index.html)
-   [DiceKriging](../packages/DiceKriging/index.html)
-   [docopulae](../packages/docopulae/index.html)
-   [DoE.base](../packages/DoE.base/index.html) (core)
-   [DoE.MIParray](../packages/DoE.MIParray/index.html)
-   [DoE.wrapper](../packages/DoE.wrapper/index.html) (core)
-   [DoseFinding](../packages/DoseFinding/index.html)
-   [dynaTree](../packages/dynaTree/index.html)
-   [easypower](../packages/easypower/index.html)
-   [edesign](../packages/edesign/index.html)
-   [EngrExpt](../packages/EngrExpt/index.html)
-   [experiment](../packages/experiment/index.html)
-   [ez](../packages/ez/index.html)
-   [FMC](../packages/FMC/index.html)
-   [FrF2](../packages/FrF2/index.html) (core)
-   [FrF2.catlg128](../packages/FrF2.catlg128/index.html)
-   [GAD](../packages/GAD/index.html)
-   [geospt](../packages/geospt/index.html)
-   [granova](../packages/granova/index.html)
-   [GroupSeq](../packages/GroupSeq/index.html)
-   [gsbDesign](../packages/gsbDesign/index.html)
-   [gsDesign](../packages/gsDesign/index.html)
-   [gset](../packages/gset/index.html)
-   [hiPOD](../packages/hiPOD/index.html)
-   [ibd](../packages/ibd/index.html)
-   [ICAOD](../packages/ICAOD/index.html)
-   [idefix](../packages/idefix/index.html)
-   [JMdesign](../packages/JMdesign/index.html)
-   [LDOD](../packages/LDOD/index.html)
-   [lhs](../packages/lhs/index.html)
-   [MAMS](../packages/MAMS/index.html)
-   [MaxPro](../packages/MaxPro/index.html)
-   [MBHdesign](../packages/MBHdesign/index.html)
-   [minimalRSD](../packages/minimalRSD/index.html)
-   [minimaxdesign](../packages/minimaxdesign/index.html)
-   [mixexp](../packages/mixexp/index.html)
-   [mkssd](../packages/mkssd/index.html)
-   [mxkssd](../packages/mxkssd/index.html)
-   [OBsMD](../packages/OBsMD/index.html)
-   [odr](../packages/odr/index.html)
-   [OPDOE](../packages/OPDOE/index.html)
-   [optbdmaeAT](../packages/optbdmaeAT/index.html)
-   [optDesignSlopeInt](../packages/optDesignSlopeInt/index.html)
-   [OptGS](../packages/OptGS/index.html)
-   [OptimalDesign](../packages/OptimalDesign/index.html)
-   [OptimaRegion](../packages/OptimaRegion/index.html)
-   [OptInterim](../packages/OptInterim/index.html)
-   [optrcdmaeAT](../packages/optrcdmaeAT/index.html)
-   [osDesign](../packages/osDesign/index.html)
-   [PBIBD](../packages/PBIBD/index.html)
-   [PGM2](../packages/PGM2/index.html)
-   [ph2bayes](../packages/ph2bayes/index.html)
-   [ph2bye](../packages/ph2bye/index.html)
-   [pid](../packages/pid/index.html)
-   [pipe.design](../packages/pipe.design/index.html)
-   [plgp](../packages/plgp/index.html)
-   [PopED](../packages/PopED/index.html)
-   [powerAnalysis](../packages/powerAnalysis/index.html)
-   [powerbydesign](../packages/powerbydesign/index.html)
-   [powerGWASinteraction](../packages/powerGWASinteraction/index.html)
-   [PwrGSD](../packages/PwrGSD/index.html)
-   [qtlDesign](../packages/qtlDesign/index.html)
-   [qualityTools](../packages/qualityTools/index.html)
-   [RcmdrPlugin.DoE](../packages/RcmdrPlugin.DoE/index.html)
-   [rodd](../packages/rodd/index.html)
-   [RPPairwiseDesign](../packages/RPPairwiseDesign/index.html)
-   [rsm](../packages/rsm/index.html) (core)
-   [rsurface](../packages/rsurface/index.html)
-   [SensoMineR](../packages/SensoMineR/index.html)
-   [seqDesign](../packages/seqDesign/index.html)
-   [sFFLHD](../packages/sFFLHD/index.html)
-   [simrel](../packages/simrel/index.html)
-   [skpr](../packages/skpr/index.html) (core)
-   [SLHD](../packages/SLHD/index.html)
-   [soptdmaeA](../packages/soptdmaeA/index.html)
-   [sp23design](../packages/sp23design/index.html)
-   [ssize.fdr](../packages/ssize.fdr/index.html)
-   [ssizeRNA](../packages/ssizeRNA/index.html)
-   [support.CEs](../packages/support.CEs/index.html)
-   [TEQR](../packages/TEQR/index.html)
-   [tgp](../packages/tgp/index.html)
-   [ThreeArmedTrials](../packages/ThreeArmedTrials/index.html)
-   [toxtestD](../packages/toxtestD/index.html)
-   [unrepx](../packages/unrepx/index.html)
-   [vdg](../packages/vdg/index.html)
-   [Vdgraph](../packages/Vdgraph/index.html)
-   [VdgRsm](../packages/VdgRsm/index.html)
-   [VNM](../packages/VNM/index.html)

### Related links:

-   [Dunn, K. (2010-2016). Process Improvement Using
    Data.](http://learnche.org/pid)
-   [Kuhnert, P. and Venables, B. (2005) *An Introduction to R: Software
    for Statistical Modelling & Computing* .
    (\~4MB)](../../doc/contrib/Kuhnert+Venables-R_Course_Notes.zip)
-   [Vikneswaran (2005). An R companion to \"Experimental
    Design\".](../../doc/contrib/Vikneswaran-ED_companion.pdf)
