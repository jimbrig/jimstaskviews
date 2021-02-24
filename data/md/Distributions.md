## CRAN Task View: Probability Distributions

  ----------------- -------------------------------------------------
  **Maintainer:**   Christophe Dutang, Patrice Kiener
  **Contact:**      Christophe.Dutang at ensimag.fr
  **Version:**      2018-04-17
  **URL:**          <https://CRAN.R-project.org/view=Distributions>
  ----------------- -------------------------------------------------

<div>

For most of the classical distributions, base R provides probability
distribution functions (p), density functions (d), quantile functions
(q), and random number generation (r). Beyond this basic functionality,
many CRAN packages provide additional useful distributions. In
particular, multivariate distributions as well as copulas are available
in contributed packages.

Ultimate bibles on probability distributions are:

-   different volumes of N. L. Johnson, S. Kotz and N. Balakrishnan
    books, e.g. Continuous Univariate Distributions, Vol. 1,
-   Thesaurus of univariate discrete probability distributions by G.
    Wimmer and G. Altmann.
-   Statistical Distributions by M. Evans, N. Hastings, B. Peacock.
-   Distributional Analysis with L-moment Statistics using the R
    Environment for Statistical Computing, Asquith (2011).

The maintainer gratefully acknowledges Achim Zeileis, David Luethi,
Tobias Verbeke, Robin Hankin, Mathias Kohl, G. Jay Kerns, Kjetil
Halvorsen, William Asquith for their useful comments/suggestions. If you
think information is not accurate or not complete, please let me know.

## [Base functionnality:]{#Base}

-   Base R provides probability distribution functions `p` *foo* `()`
    density functions `d` *foo* `()`, quantile functions `q` *foo* `()`,
    and random number generation `r` *foo* `()` where *foo* indicates
    the type of distribution: beta ( *foo* = `beta`), binomial `binom`,
    Cauchy `cauchy`, chi-squared `chisq`, exponential `exp`, Fisher F
    `f`, gamma `gamma`, geometric `geom`, hypergeometric `hyper`,
    logistic `logis`, lognormal `lnorm`, negative binomial `nbinom`,
    normal `norm`, Poisson `pois`, Student t `t`, uniform `unif`,
    Weibull `weibull`. Following the same naming scheme, but somewhat
    less standard are the following distributions in base R:
    probabilities of coincidences (also known as \"birthday paradox\")
    `birthday` (only p and q), studentized range distribution `tukey`
    (only p and q), Wilcoxon signed rank distribution `signrank`,
    Wilcoxon rank sum distribution `wilcox`.
-   Probability generating function:
    [Compounding](../packages/Compounding/index.html) provides pgf for
    `xxx` distribution, inverse `xxx` distribution, first derivative of
    the `xxx` distribution, where `xxx` belongs to binomial,
    binomial-Poisson, geometric, hypergeometric, hyper-Poisson, Katti
    type H1/H2, logarithmic, logarithmic-binomial, logarithmic-Poisson,
    negative binomial, Neyman type A/B/C, Pascal-Poisson, Poisson,
    Poisson-binomial, Poisson-Lindley, Poisson-Pascal, Polya Aeppli,
    Thomas, Waring, Yule.

## [Discrete univariate distributions:]{#UnivariateDiscrete}

-   *Beta-binomial distribution* : provided in
    [VGAM](../packages/VGAM/index.html),
    [extraDistr](../packages/extraDistr/index.html). ZI/ZM beta binomial
    distributions are implemented in
    [gamlss.dist](../packages/gamlss.dist/index.html).
-   *Beta-geometric distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Binomial (including Bernoulli) distribution* : provided in
    **stats** . Zero-modified, zero-inflated, truncated versions are
    provided in [gamlss.dist](../packages/gamlss.dist/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [actuar](../packages/actuar/index.html) and in
    [VGAM](../packages/VGAM/index.html).
    [rpgm](../packages/rpgm/index.html) provides a fast random number
    generator. [LaplacesDemon](../packages/LaplacesDemon/index.html)
    provides dedicated functions for the Bernoulli distribution.\
      ---------------------- ------------- ------------- -----------------------
      *Distribution name*    *Packages*    *Functions*   *Distribution suffix*
      binomial               stats         d, p, q, r    `binom`
      zero-infl. binomial    extraDistr    d, p, q, r    `zib`
      zero-infl. binomial    VGAM          d, p, q, r    `zibinom`
      zero-infl. binomial    gamlss.dist   d, p, q, r    `ZIBI`
      zero mod. binomial     VGAM          d, p, q, r    `zabinom`
      zero mod. binomial     actuar        d, p, q, r    `zmbinom`
      zero mod. binomial     gamlss.dist   d, p, q, r    `ZABI`
      zero trunc. binomial   actuar        d, p, q, r    `ztbinom`
      trunc. binomial        extraDistr    d, p, q, r    `tbinom`
      ---------------------- ------------- ------------- -----------------------

      : Summary for Binomial-related distributions

    \
-   *Benford distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Bernoulli distribution* : provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Borel-Tanner distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Conway-Maxwell-Poisson distribution* : provided in
    [compoisson](../packages/compoisson/index.html) and
    [CompGLM](../packages/CompGLM/index.html).
-   *Cantor distribution* : fast random generator provided in
    [rpgm](../packages/rpgm/index.html).
-   *Delaporte distribution* : provided in
    [gamlss.dist](../packages/gamlss.dist/index.html) and
    [Delaporte](../packages/Delaporte/index.html).
-   *Dirac distribution* : provided in
    [distr](../packages/distr/index.html).
-   *Discrete categorical distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Discrete exponential distribution* : provided in
    [poweRlaw](../packages/poweRlaw/index.html).
-   *Discrete gamma distribution* : provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Discrete inverse Weibull distribution* :
    [DiscreteInverseWeibull](../packages/DiscreteInverseWeibull/index.html)
    provides d, p, q, r functions for the inverse Weibull as well as
    hazard rate function and moments.
-   *Discrete Laplace distribution* : The discrete Laplace distribution
    is provided in [extraDistr](../packages/extraDistr/index.html) (d,
    p, r). The skew discrete Laplace distribution has two
    parametrization (DSL and ADSL), both provided in
    [DiscreteLaplace](../packages/DiscreteLaplace/index.html) and DSL in
    [disclap](../packages/disclap/index.html).
    [LaplacesDemon](../packages/LaplacesDemon/index.html) also provides
    the DSL parametrization only.
-   *Discrete lognormal distribution* : provided in
    [poweRlaw](../packages/poweRlaw/index.html).
-   *Discrete normal distribution* : provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Discrete uniform distribution* : can be easily obtained with the
    functions `sum,cumsum,sample` and is provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Discrete Weibull distribution* : provided in
    [DiscreteWeibull](../packages/DiscreteWeibull/index.html): d, p, q,
    r, m for disc. Weib. type 1, d, p, q, r, m, h for disc. Weib. type
    3. [extraDistr](../packages/extraDistr/index.html) provides d, p, q,
    r for Type 1.
-   *Felix distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Lindley distribution* : provided in
    [gambin](../packages/gambin/index.html).
-   *Geometric distribution* : provided in **stats** . Zero-modified,
    zero-inflated, truncated versions are provided in
    [gamlss.dist](../packages/gamlss.dist/index.html),
    [actuar](../packages/actuar/index.html) and in
    [VGAM](../packages/VGAM/index.html).
    [rpgm](../packages/rpgm/index.html) provides a fast random number
    generator.
-   *Geometric (compound) Poisson distribution (also known Polya-Aeppli
    distribution)* : provided in
    [polyaAeppli](../packages/polyaAeppli/index.html).
-   *Generalized binomial distribution* : provided in
    [GenBinomApps](../packages/GenBinomApps/index.html).
-   *Generalized Hermite distribution* : provided in
    [hermite](../packages/hermite/index.html).
-   *Hypergeometric distribution* : provided in **stats** . Extented
    hypergeometric distribution can be found in
    [BiasedUrn](../packages/BiasedUrn/index.html) package, which
    provides not only p, d, q, r functions but also mean, variance, mode
    functions. Generalized hypergeometric distribution is implemented in
    [SuppDists](../packages/SuppDists/index.html). Negative
    hypergeometric distribution is provided in
    [tolerance](../packages/tolerance/index.html),
    [extraDistr](../packages/extraDistr/index.html).
-   *Lagrangian Poisson distribution* :
    [RMKdiscrete](../packages/RMKdiscrete/index.html) provides d, p, q,
    r functions for the univariate and the bivariate Lagrangian Poisson
    distribution.
-   *Lindley distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Logarithmic distribution* : This can be found in
    [extraDistr](../packages/extraDistr/index.html),
    [VGAM](../packages/VGAM/index.html),
    [actuar](../packages/actuar/index.html) and
    [gamlss.dist](../packages/gamlss.dist/index.html). Zero-modified and
    zero-truncated versions is provided in
    [actuar](../packages/actuar/index.html). A fast random generator is
    available for the logarithmic distribution is implemented in
    [Runuran](../packages/Runuran/index.html) as well as the \'density\'
    function.
-   *Poisson distribution* : provided in **stats** and in
    [poweRlaw](../packages/poweRlaw/index.html). Zero-modified,
    zero-inflated, truncated versions are provided in
    [extraDistr](../packages/extraDistr/index.html),
    [gamlss.dist](../packages/gamlss.dist/index.html),
    [actuar](../packages/actuar/index.html) and in
    [VGAM](../packages/VGAM/index.html).
    [extraDistr](../packages/extraDistr/index.html) provides the
    truncated Poisson distribution.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides the
    generalized Poisson distribution. See the mixture section such as
    the Poisson-lognormal mixture.
-   *Poisson-Lindley distribution* : provided in
    [tolerance](../packages/tolerance/index.html).
-   *Power law distribution* : provided in
    [poweRlaw](../packages/poweRlaw/index.html).
-   *Mana Clash distribution* : provided in
    [RMKdiscrete](../packages/RMKdiscrete/index.html).
-   *Negative binomial distribution* : provided in **stats** .
    Zero-modified, zero-inflated, truncated versions are provided in
    [gamlss.dist](../packages/gamlss.dist/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [actuar](../packages/actuar/index.html) and in
    [VGAM](../packages/VGAM/index.html). New parametrization of the
    negative binomial distribution is available in
    [RMKdiscrete](../packages/RMKdiscrete/index.html).
-   *Sichel distribution* : provided in
    [gamlss.dist](../packages/gamlss.dist/index.html).
-   *Skellam distribution* : provided in
    [extraDistr](../packages/extraDistr/index.html),
    [VGAM](../packages/VGAM/index.html) and
    [skellam](../packages/skellam/index.html).
-   *Waring distribution* : sampling in
    [degreenet](../packages/degreenet/index.html).
-   *Yule-Simon distribution* : provided in
    [VGAM](../packages/VGAM/index.html) and sampling in
    [degreenet](../packages/degreenet/index.html).
-   *Zeta and Haight\'s Zeta distribution* : provided in
    [VGAM](../packages/VGAM/index.html),
    [tolerance](../packages/tolerance/index.html).
-   *Zipf law* : d, p, q, r functions of the Zipf and the
    Zipf-Mandelbrot distributions are provided in
    [tolerance](../packages/tolerance/index.html),
    [VGAM](../packages/VGAM/index.html). Package
    [zipfR](../packages/zipfR/index.html) provides tools for
    distribution of word frequency, such as the Zipf distribution.

## [Discrete multivariate distributions:]{#MultivariateDiscrete}

-   *Bivariate Poisson-lognormal* : provided in
    [poilog](../packages/poilog/index.html).
-   *Dirichlet distribution* :
    [Compositional](../packages/Compositional/index.html) and
    [LaplacesDemon](../packages/LaplacesDemon/index.html) packages
    provide d, r functions as well as a fitting function for
    [Compositional](../packages/Compositional/index.html).
-   *Hyper Dirichlet distribution* : provided in
    [hyper2](../packages/hyper2/index.html) package.
-   *Multinomial distribution* : stats,
    [mc2d](../packages/mc2d/index.html),
    [extraDistr](../packages/extraDistr/index.html) packages provide d,
    r functions.
-   *Negative multinomial distribution* : A bivariate distribution with
    negative-binomial marginals is available in
    [RMKdiscrete](../packages/RMKdiscrete/index.html). The
    multiplicative multinomial distribution is implemented in
    [MM](../packages/MM/index.html).
-   *Multivariate Poisson distribution* : not yet implemented?
-   *Multivariate hypergeometric distribution* : provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Multivariate Polya distribution* : functions d, r of the Dirichlet
    Multinomial (also known as multivariate Polya) distribution are
    provided in [extraDistr](../packages/extraDistr/index.html),
    [LaplacesDemon](../packages/LaplacesDemon/index.html) and
    [Compositional](../packages/Compositional/index.html).
-   *Multivariate Ewens distribution* : not yet implemented?
-   *Truncated Stick-Breakin distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).

## [Continuous univariate distributions:]{#UnivariateContinuous}

-   *Arcsine distribution* : implemented in package
    [distr](../packages/distr/index.html).
-   *Beta distribution and its extensions* : Base R provides the d, p,
    q, r functions for this distribution (see above).
    [extraDistr](../packages/extraDistr/index.html) provides the beta
    distribution parametrized by the mean and the precision.
    [actuar](../packages/actuar/index.html) provides moments and limited
    expected values. [sadists](../packages/sadists/index.html)
    implements Gram Charlier, Edgeworth and Cornish-Fisher
    approximations for doubly non central beta distribution for
    computing d, p, q, r functions.
    [extraDistr](../packages/extraDistr/index.html) provides the
    four-parameter beta with lower and upper bounds. The generalized
    beta of the first kind (GB1) (exponentation of beta 1) is provided
    in [gamlss.dist](../packages/gamlss.dist/index.html),
    [mbbefd](../packages/mbbefd/index.html),
    [actuar](../packages/actuar/index.html). The beta prime (or beta of
    the second kind), which is the distribution of X/(1-X) when X
    follows a beta distribution of the first kind, is provided in
    [VGAM](../packages/VGAM/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [LaplacesDemon](../packages/LaplacesDemon/index.html) and
    [mc2d](../packages/mc2d/index.html). The zero and one inflated beta
    distribution can be found in
    [gamlss.dist](../packages/gamlss.dist/index.html). The generalized
    beta of the second kind (GB2) is provided in
    [gamlss.dist](../packages/gamlss.dist/index.html),
    [GB2](../packages/GB2/index.html). Several special cases of the
    generalized beta distribution are also implemented in
    [VGAM](../packages/VGAM/index.html),
    [mc2d](../packages/mc2d/index.html): Lomax, inverse Lomax, Dagum,
    Singh-Maddala, Pert distributions.
    [actuar](../packages/actuar/index.html) provides the transformed
    beta 2 distribution which includes as special cases Burr,
    loglogistic, paralogistic, generalized Pareto, Pareto, see also the
    Pareto subsection.\
      ------------------------- --------------------------------------------------- -------------------- -----------------------
      *Distribution name*       *Packages*                                          *Functions*          *Distribution suffix*
      Beta (1st kind)           stats                                               d, p, q, r           `beta`
      Beta                      [actuar](../packages/actuar/index.html)             m, mgf, lev          `beta`
      Beta                      [extraDistr](../packages/extraDistr/index.html)     d, p, q, r           `prop`
      Doubly non central beta   [sadists](../packages/sadists/index.html)           d, p, q, r           `nbeta`
      4-param beta              [extraDistr](../packages/extraDistr/index.html)     d, p, q, r           `nsbeta`
      zero-infl beta            [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r           `BEZI`
      one-infl beta             [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r           `BEOI`
      one-infl beta             [mbbefd](../packages/mbbefd/index.html)             d, p, q, r, m, ec    `oibeta`
      GB1                       [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r           `GB1`
      GB1                       [mbbefd](../packages/mbbefd/index.html)             d, p, q, r, m, ec    `gbeta`
      GB1                       [actuar](../packages/actuar/index.html)             d, p, q, r, m, lev   `genbeta`
      one-infl GB1              [mbbefd](../packages/mbbefd/index.html)             d, p, q, r, m, ec    `oigbeta`
      ------------------------- --------------------------------------------------- -------------------- -----------------------

      : Summary for Beta-related distributions

    \
    \
      --------------------- ------------------------------------------------------- -------------------- -----------------------
      *Distribution name*   *Packages*                                              *Functions*          *Distribution suffix*
      Beta (2nd kind)       [VGAM](../packages/VGAM/index.html)                     d, p, q, r           `beta`
      Beta (2nd kind)       [extraDistr](../packages/extraDistr/index.html)         d, p, q, r           `invbeta`
      Beta (2nd kind)       [LaplacesDemon](../packages/LaplacesDemon/index.html)   d, r                 `betapr`
      GB2                   [VGAM](../packages/VGAM/index.html)                     d, p, q, r           `genbetaII`
      GB2                   [gamlss.dist](../packages/gamlss.dist/index.html)       d, p, q, r           `GB2`
      GB2                   [GB2](../packages/GB2/index.html)                       d, p, q, r           `gb2`
      Trans beta 2          [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `trbeta`
      --------------------- ------------------------------------------------------- -------------------- -----------------------

      : Summary for Beta-2-related distributions

    \
-   *Benini distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Bhattacharjee (normal+uniform) distribution* : provided in package
    [extraDistr](../packages/extraDistr/index.html).
-   *Birnbaum-Saunders distribution* : provided in package
    [VGAM](../packages/VGAM/index.html) and
    [extraDistr](../packages/extraDistr/index.html).
-   *Bridge distribution* : provided in
    [bridgedist](../packages/bridgedist/index.html), as detailed in Wang
    and Louis (2003). The distribution of random intercept that allows a
    marginalized random intercept logistic regression to also be
    logistic regression.
-   *Box Cox distribution* :
    [gamlss.dist](../packages/gamlss.dist/index.html) provides the
    Box-Cox normal, the Box-Cox power exponential and the Box-Cox t
    distributions.
-   *Burr distribution* : see Pareto.
-   *Cardioid distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Cauchy distribution* : Base R provides the d, p, q, r functions for
    this distribution (see above). Other implementations are available
    in [lmomco](../packages/lmomco/index.html) and
    [sgt](../packages/sgt/index.html). The skew Cauchy distribution is
    provided in [sn](../packages/sn/index.html).
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, p,
    q, r functions for the Half-Cauchy distribution.
-   *Chen distribution* : provided in
    [reliaR](../packages/reliaR/index.html).
-   *Chi(-squared or not) distribution* : Base R provides the d, p, q, r
    functions for the chi-squared distribution, both central and
    non-central (see above). Moments, limited expected values and the
    moment generating function are provided in
    [actuar](../packages/actuar/index.html).
    [extraDistr](../packages/extraDistr/index.html) provides d, p, q, r
    functions for inverse chi-squared distribution (standard and
    scaled). Only d,r functions are available for the inverse
    chi-squared distribution in package
    [geoR](../packages/geoR/index.html) and
    [LaplacesDemon](../packages/LaplacesDemon/index.html). A fast random
    generator is available for the Chi distribution is implemented in
    [Runuran](../packages/Runuran/index.html) as well as the density
    function. The non-central Chi distribution is not yet implemented.
    The chi-bar-squared distribution is implemented in
    [emdbook](../packages/emdbook/index.html).
    [sadists](../packages/sadists/index.html) implements Gram Charlier,
    Edgeworth and Cornish-Fisher approximations for sums of non central
    chi-squared raised to powers distribution and sums of log of non
    central chi-squared for computing d, p, q, r functions.\
    \
      ---------------------------- ------------------------------------------------- ------------- -----------------------
      *Distribution name*          *Packages*                                        *Functions*   *Distribution suffix*
      Chi-squared                  stats                                             d, p, q, r    `chisq`
      Chi-squared                  [actuar](../packages/actuar/index.html)           m, mgf, lev   `chisq`
      Chi-squared                  [Runuran](../packages/Runuran/index.html)         d, r          `chisq`
      Chi-bar-squared              [emdbook](../packages/emdbook/index.html)         d, p, q, r    `chibarsq`
      Chi                          [Runuran](../packages/Runuran/index.html)         d, r          `chi`
      Inverse Chi-squared          [geoR](../packages/geoR/index.html)               d, r          `invchisq`
      Inverse Chi-squared          [extraDistr](../packages/extraDistr/index.html)   d, p, q, r    `invchisq`
      Scaled Inverse Chi-squared   [extraDistr](../packages/extraDistr/index.html)   d, p, q, r    `invchisq`
      Sum of power Chi-squared     [sadists](../packages/sadists/index.html)         d, p, q, r    `sumchisqpow`
      Sum of log Chi-squared       [sadists](../packages/sadists/index.html)         d, p, q, r    `sumlogchisq`
      ---------------------------- ------------------------------------------------- ------------- -----------------------

      : Summary for Chi-related distributions

    \
-   *Dagum distribution* : see beta.
-   *Davies distribution* : The Davies distribution is provided in
    [Davies](../packages/Davies/index.html) package.
-   *(non-central) Dunnett\'s test distribution* : provided in
    [nCDunnett](../packages/nCDunnett/index.html).
-   *Eta-mu distribution* : provided in
    [lmomco](../packages/lmomco/index.html).
    [sadists](../packages/sadists/index.html) implements Gram Charlier,
    Edgeworth and Cornish-Fisher approximations for doubly non central
    eta distribution for computing d, p, q, r functions.
-   *Exponential distribution and its extensions* : Base R provides the
    d, p, q, r functions for this distribution (see above).
    [actuar](../packages/actuar/index.html) provides additional
    functions such as the moment generating function, moments and
    limited expected values. It also has the d, p, q, r for the inverse
    exponential distribution. The shifted (or two-parameter exponential)
    and the truncated exponential distributions are implemented in
    [lmomco](../packages/lmomco/index.html) and
    [tolerance](../packages/tolerance/index.html) packages with d, p, q,
    r functions. Exponential Power distribution is also known as General
    Error Distribution: d, p, q, r functions for the power and the skew
    power exponential type 1-4 distributions are implemented in
    [gamlss.dist](../packages/gamlss.dist/index.html) and
    [lmomco](../packages/lmomco/index.html). The power exponential
    distribution is also provided in
    [normalp](../packages/normalp/index.html),
    [LaplacesDemon](../packages/LaplacesDemon/index.html) and
    [sgt](../packages/sgt/index.html). The skew power exponential is
    provided [sgt](../packages/sgt/index.html).
    [reliaR](../packages/reliaR/index.html) provides the generalized
    exponential, the inverse generalized exponential, the logistic
    exponential, the Marshall-Olkin Extended Exponential and the
    exponential extension distributions. A fast random generator is
    available for the power Exponential distribution is implemented in
    [Runuran](../packages/Runuran/index.html) as well as the density
    function. [rpgm](../packages/rpgm/index.html) provides a fast random
    number generator for the exponential distribution.\
    \
      -------------------------------------------------- --------------------------------------------------- ---------------------- ----------------------------
      *Distribution name*                                *Packages*                                          *Functions*            *Distribution suffix*
      Exponential                                        stats                                               d, p, q, r             `exp`
      Exponential                                        [actuar](../packages/actuar/index.html)             m, mgf, lev            `exp`
      Exponential                                        [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r             `EXP`
      Exponential                                        [poweRlaw](../packages/poweRlaw/index.html)         d, p, q, r             `exp`
      Exponential                                        [rpgm](../packages/rpgm/index.html)                 r                      `rpgm.rexp`
      Inverse exponential                                [actuar](../packages/actuar/index.html)             d, p, q, r, m, lev     `invexp`
      Shifted exponential                                [lmomco](../packages/lmomco/index.html)             d, p, q, r, lm, tlmr   `exp`
      Shifted exponential                                [tolerance](../packages/tolerance/index.html)       d, p, q, r             `2exp`
      Truncated exponential                              [lmomco](../packages/lmomco/index.html)             d, p, q, r, lm, tlmr   `texp`
      Truncated exponential                              [ReIns](../packages/ReIns/index.html)               d, p, q, r             `texp`
      Power exponential                                  [normalp](../packages/normalp/index.html)           d, p, q, r             `normp`
      Power exponential                                  [Runuran](../packages/Runuran/index.html)           d, r                   `exp`
      Skew power exp.                                    [lmomco](../packages/lmomco/index.html)             d, p, q, r, lm, tlmr   `aep4`
      Power and skew power exp.                          [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r             `PE, SEP`
      Generalized and inverse gen. exp.                  [reliaR](../packages/reliaR/index.html)             d, p, q, r             `gen.exp, inv.genexp`
      Logistic, Marshall-Olkin Ext. exp. and exp. ext.   [reliaR](../packages/reliaR/index.html)             d, p, q, r             `logis.exp, moee, exp.ext`
      -------------------------------------------------- --------------------------------------------------- ---------------------- ----------------------------

      : Summary for exponential-related distributions

    \
-   *Externally studentized midrange distribution* : Package
    [SMR](../packages/SMR/index.html) computes the studentized midrange
    distribution (d, p, q, r).
-   *Fisher-Snedecor (or F) distribution* : Base R provides the d, p, q,
    r functions for the F distribution, possibly with a non-central
    parameter. [sadists](../packages/sadists/index.html) implements Gram
    Charlier, Edgeworth and Cornish-Fisher approximations for doubly non
    central Fisher distribution (and product of multiple doubly non
    central Fisher distribution) for computing d, p, q, r functions.
    [flexsurv](../packages/flexsurv/index.html) provides d, p, q, r
    functions as well as hazard (h) and integrated hazard rate (i)
    functions for the generalized F distribution.
    [fpow](../packages/fpow/index.html) returns the noncentrality
    parameter of the noncentral F distribution if probability of type I
    and type II error, degrees of freedom of the numerator and the
    denominator are given.
-   *Frechet distribution* : provided in
    [VGAM](../packages/VGAM/index.html),
    [RTDE](../packages/RTDE/index.html),
    [ReIns](../packages/ReIns/index.html),
    [extraDistr](../packages/extraDistr/index.html) and
    [evd](../packages/evd/index.html). A fast random generator is
    available for the Frechet distribution is implemented in
    [Runuran](../packages/Runuran/index.html) as well as the density
    function. The truncated Frechet distribution is provided in
    [ReIns](../packages/ReIns/index.html).
-   *Friedman\'s Chi distribution* : provided in
    [SuppDists](../packages/SuppDists/index.html).
-   *Gamma distribution and its extensions* : Base R provides the d, p,
    q, r functions for this distribution (see above).
    [EnvStats](../packages/EnvStats/index.html) provides d, p, q, r
    functions of the gamma parametrized by the mean and the coefficient
    of variation. [actuar](../packages/actuar/index.html) provides d, p,
    q, r functions of the inverse, the inverse transformed and the log
    gamma distributions while [ghyp](../packages/ghyp/index.html)
    provides those functions for the variance gamma distribution.
    [extraDistr](../packages/extraDistr/index.html) and
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provide the
    inverse gamma distribution.
    [VarianceGamma](../packages/VarianceGamma/index.html) provides d, p,
    q, r functions for the variance gamma distribution as well as
    moments (skewness, kurtosis, \...).
    [VGAM](../packages/VGAM/index.html) provides d, p, q, r functions of
    the log gamma and the generalized gamma distribution. The
    generalized gamma distribution can also be found in
    [gamlss.dist](../packages/gamlss.dist/index.html).
    [reliaR](../packages/reliaR/index.html) provides the log gamma
    distribution. See Pearson III for a three-parameter gamma
    distribution with a location parameter.
    [flexsurv](../packages/flexsurv/index.html) provides d, p, q, r
    functions as well as hazard (h) and integrated hazard rate (i)
    functions for the generalized gamma distribution.
    [coga](../packages/coga/index.html) provides d, p, r functions for a
    sum of independent but not identically distributed gamma
    distributions.\
    \
      ---------------------- ------------------------------------------------------- ------------------------- -----------------------
      *Distribution name*    *Packages*                                              *Functions*               *Distribution suffix*
      Gamma                  stats                                                   d, p, q, r                `gamma`
      Gamma                  [actuar](../packages/actuar/index.html)                 m, mgf, lev               `gamma`
      Gamma                  [EnvStats](../packages/EnvStats/index.html)             d, p, q, r                `gammaAlt`
      Inverse gamma          [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev, mgf   `invgamma`
      Inverse gamma          [extraDistr](../packages/extraDistr/index.html)         d, p, q, r                `invgamma`
      Inverse gamma          [LaplacesDemon](../packages/LaplacesDemon/index.html)   d, r                      `invgamma`
      Log-gamma              [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev        `lgamma`
      Log-gamma              [VGAM](../packages/VGAM/index.html)                     d, p, q, r                `lgamma`
      Variance gamma         [ghyp](../packages/ghyp/index.html)                     d, p, q, r                `VG`
      Variance gamma         [VarianceGamma](../packages/VarianceGamma/index.html)   d, p, q, r, m             `vg`
      Generalized gamma      [flexsurv](../packages/flexsurv/index.html)             d, p, q, r, h, i          `gengamma`
      Generalized gamma      [gamlss.dist](../packages/gamlss.dist/index.html)       d, p, q, r                `GG`
      Generalized gamma      [VGAM](../packages/VGAM/index.html)                     d, p, q, r                `gengamma.stacy`
      convolution of gamma   [coga](../packages/coga/index.html)                     d, p, r                   `coga`
      ---------------------- ------------------------------------------------------- ------------------------- -----------------------

      : Summary for gamma-related distributions

    \
-   *Gaussian (or normal) distribution and its extensions* : Base R
    provides the d, p, q, r functions for this distribution (see above).
    [actuar](../packages/actuar/index.html) provides the moment
    generating function and moments. The
    [truncnorm](../packages/truncnorm/index.html) package provides d, p,
    q, r functions for the truncated gaussian distribution as well as
    functions for the first two moments.
    [mvrtn](../packages/mvrtn/index.html) provides random variates for
    left/right truncated normal distributions.
    [EnvStats](../packages/EnvStats/index.html) provides d, p, q, r
    functions for the truncated normal distribution and the
    zero-modified distribution.
    [extraDistr](../packages/extraDistr/index.html) provides the
    truncated normal. [rpgm](../packages/rpgm/index.html) provides a
    fast random number generator for the normal distribution.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, p,
    q, r functions for the Half-normal distribution.
    [lmomco](../packages/lmomco/index.html) implements the generalized
    normal distribution. The Exponentially modified Gaussian is
    available in [emg](../packages/emg/index.html),
    [gamlss.dist](../packages/gamlss.dist/index.html) and
    [retimes](../packages/retimes/index.html).
    [sn](../packages/sn/index.html) implements the skew normal
    distribution. [VGAM](../packages/VGAM/index.html) implements the
    folded and the skewed normal distribution, and
    [csn](../packages/csn/index.html) provides d, r functions for the
    closed skew normal distribution.
    [CompQuadForm](../packages/CompQuadForm/index.html) provides the
    distribution function of quadratic forms in normal variates.
    [NormalGamma](../packages/NormalGamma/index.html) provides the
    density of the sum of a gaussian and a gamma random variables.
    [NormalLaplace](../packages/NormalLaplace/index.html) provides d, p,
    q, r functions for the sum of a normal and a Laplace random
    variables, while
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    function of the sum of a normal and a Laplace random variables.\
      --------------------------------- --------------------------------------------------- --------------- -----------------------
      *Distribution name*               *Packages*                                          *Functions*     *Distribution suffix*
      Normal                            stats                                               d, p, q, r      `norm`
      Normal                            [actuar](../packages/actuar/index.html)             m, mgf          `norm`
      Truncated normal                  [truncnorm](../packages/truncnorm/index.html)       d, p, q, r, m   `truncnorm`
      Truncated normal                  [mvrtn](../packages/mvrtn/index.html)               r, m            `tn`
      Truncated normal                  [EnvStats](../packages/EnvStats/index.html)         d, p, q, r      `normTrunc`
      Truncated normal                  [extraDistr](../packages/extraDistr/index.html)     d, p, q, r      `tnorm`
      Generalized normal                [lmomco](../packages/lmomco/index.html)             d, p, q, r      `gno`
      Zero modified Gaussian            [EnvStats](../packages/EnvStats/index.html)         d, p, q, r      `zmnorm`
      Exponentially modified Gaussian   [emg](../packages/emg/index.html)                   d, p, q, r      `emg`
      Exponentially modified Gaussian   [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r      `exGAUSS`
      Exponentially modified Gaussian   [retimes](../packages/retimes/index.html)           d, p, q, r      `exgauss`
      Folded and skew normal            [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r      `SN1, SN2`
      Closed skew normal                [csn](../packages/csn/index.html)                   d, p, q, r      `csn`
      Skew normal                       [sn](../packages/sn/index.html)                     d, p, q, r      `sn`
      --------------------------------- --------------------------------------------------- --------------- -----------------------

      : Summary for Gaussian-related distributions

    \
-   *General error distribution (also known as exponential power
    distribution)* : see *exponential* item.
-   *Generalized extreme value distribution* : provided in
    [lmomco](../packages/lmomco/index.html) (d, p, q);
    [VGAM](../packages/VGAM/index.html),
    [evd](../packages/evd/index.html),
    [evir](../packages/evir/index.html),
    [FAdist](../packages/FAdist/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [EnvStats](../packages/EnvStats/index.html),
    [QRM](../packages/QRM/index.html) and
    [fExtremes](../packages/fExtremes/index.html) (d, p, q, r).
    [evdbayes](../packages/evdbayes/index.html),
    [revdbayes](../packages/revdbayes/index.html) provide d,p,q,r
    functions of the GEV distribution in a Bayesian setting.
-   *Gompertz distribution* : provided in
    [reliaR](../packages/reliaR/index.html),
    [flexsurv](../packages/flexsurv/index.html),
    [extraDistr](../packages/extraDistr/index.html).
    [flexsurv](../packages/flexsurv/index.html) also provides hazard (h)
    and integrated hazard rate (i) functions. The shifted Gompertz
    distribution is implemented in
    [extraDistr](../packages/extraDistr/index.html).
-   *Govindarajulu distribution* : provided in
    [lmomco](../packages/lmomco/index.html).
-   *Gumbel distribution* : provided in packages
    [lmomco](../packages/lmomco/index.html),
    [VGAM](../packages/VGAM/index.html),
    [gamlss.dist](../packages/gamlss.dist/index.html),
    [FAdist](../packages/FAdist/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [reliaR](../packages/reliaR/index.html),
    [QRM](../packages/QRM/index.html),
    [EnvStats](../packages/EnvStats/index.html) and
    [evd](../packages/evd/index.html).
    [actuar](../packages/actuar/index.html) provides the raw moments and
    the moment generating function (mgf) in addition to the d, p, q, r
    functions. A fast random generator is available for the Gumbel
    distribution is implemented in
    [Runuran](../packages/Runuran/index.html) as well as the density
    function. The reverse Gumbel distribution is implemented in
    [lmomco](../packages/lmomco/index.html) and
    [gamlss.dist](../packages/gamlss.dist/index.html).
-   *Huber distribution* : Huber\'s least favourable distribution
    provided in package [smoothmest](../packages/smoothmest/index.html)
    (d, r), and in [VGAM](../packages/VGAM/index.html),
    [marg](../packages/marg/index.html),
    [extraDistr](../packages/extraDistr/index.html) (d, p, q, r).
-   *(generalized) Hyperbolic distribution* :
    [fBasics](../packages/fBasics/index.html),
    [ghyp](../packages/ghyp/index.html),
    [GeneralizedHyperbolic](../packages/GeneralizedHyperbolic/index.html)
    and [HyperbolicDist](../packages/HyperbolicDist/index.html) packages
    provide d, p, q, r functions for the generalized hyperbolic
    distribution. [QRM](../packages/QRM/index.html) provides d, r
    functions for the generalized hyperbolic distribution.
    [SkewHyperbolic](../packages/SkewHyperbolic/index.html) provides the
    skewed Hyperbolic Student t-Distribution.
    [fBasics](../packages/fBasics/index.html) also implements the
    standardized generalized Hyperbolic distribution. A fast random
    generator is available for the hyperbolic distribution is
    implemented in [Runuran](../packages/Runuran/index.html) as well as
    the density function.
-   *Hyperbolic sine distribution and extension* :
    [gamlss.dist](../packages/gamlss.dist/index.html) provides the sinh
    and the asinh distributions. [ihs](../packages/ihs/index.html)
    provides the asinh distribution. Generalized Power Hyperbolic sine
    distributions are provided in
    [FatTailsR](../packages/FatTailsR/index.html).
-   *Inverse Gaussian (also known Wald) distribution* : d, p, q, and r
    functions of the inverse Gaussian are provided in
    [statmod](../packages/statmod/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [SuppDists](../packages/SuppDists/index.html) and
    [STAR](../packages/STAR/index.html).
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    functions for the inverse Gaussian distribution.
    [actuar](../packages/actuar/index.html) provides d, p, q, r, m, lev,
    mgf functions for the Inverse Gaussian distribution.
    [SuppDists](../packages/SuppDists/index.html) also provides a
    function that returns moments, skewness, kurtosis.
    [fBasics](../packages/fBasics/index.html) the normal inverse
    Gaussian and standardized normal inverse Gaussian distributions. The
    generalized inverse gaussian distribution can be found in
    [gamlss.dist](../packages/gamlss.dist/index.html),
    [QRM](../packages/QRM/index.html) and
    [HyperbolicDist](../packages/HyperbolicDist/index.html). A random
    generator is available for the (generalized) Inverse Gaussian
    distribution is implemented in
    [Runuran](../packages/Runuran/index.html) as well as the density
    function. [GIGrvg](../packages/GIGrvg/index.html) generates random
    variables from the generalized inverse Gaussian distribution.
    [frmqa](../packages/frmqa/index.html) computes p function of the
    generalized inverse Gaussian distribution.
-   *Johnson distribution* : provided in
    [SuppDists](../packages/SuppDists/index.html).
-   *K-prime distribution* : [sadists](../packages/sadists/index.html)
    implements Gram Charlier, Edgeworth and Cornish-Fisher
    approximations for K-prime distribution for computing d, p, q, r
    functions.
-   *Kappa distribution* : A 4-parameter Kappa distribution is provided
    in [lmomco](../packages/lmomco/index.html) and
    [FAdist](../packages/FAdist/index.html).
-   *Kappa-mu distribution* : provided in
    [lmomco](../packages/lmomco/index.html).
-   *Kendall\'s tau distribution* : provided in
    [SuppDists](../packages/SuppDists/index.html).
-   *Kiener distribution* : a family of distributions generalizing
    hyperbolic sine distributions (see hyperbolic sine section), d, p,
    q, r, m provided in [FatTailsR](../packages/FatTailsR/index.html).
-   *Kolmogorov distribution* : p function provided in
    [kolmim](../packages/kolmim/index.html).
-   *Kruskal Wallis distribution* : provided in
    [SuppDists](../packages/SuppDists/index.html).
-   *Kumaraswamy distribution* : provided in packages
    [VGAM](../packages/VGAM/index.html),
    [extraDistr](../packages/extraDistr/index.html) and
    [lmomco](../packages/lmomco/index.html).
-   *(Tukey) Lambda distribution and its extensions* : The generalized
    Lambda distribution (GLD) is well known for its wide range of
    shapes. The original Tukey Lambda distribution can be obtained as a
    special case of the generalized Lambda distribution. There exists
    different parametrization of GLD in the literature: RS
    (Ramberg-Schmeiser or tail-index param), FMKL
    (Freimer-Mudholkar-Kollia-Lin), FM5 (Five-parameter version of FKML
    by Gilchrist), GPD (gen. Pareto dist.) and AS (Asymmetry-steepness).
    The following packages implement such distributions (with d, p, q, r
    functions): [gld](../packages/gld/index.html) (RS, FKML, FM5, GPD),
    [Davies](../packages/Davies/index.html) (RS),
    [gb](../packages/gb/index.html) (RS),
    [lmomco](../packages/lmomco/index.html) (FMKL),
    [extraDistr](../packages/extraDistr/index.html) (original Tukey).
    [ecd](../packages/ecd/index.html) provides the elliptic lambda
    distribution and its use for financial pricing.
-   *Tukey\'s H distribution* : provided as a special case of Lambert W
    x F distribution.
-   *Lambda-prime distribution* :
    [sadists](../packages/sadists/index.html) implements Gram Charlier,
    Edgeworth and Cornish-Fisher approximations for K-prime distribution
    for computing d, p, q, r functions.
-   *Lambert W x F distribution* :
    [LambertW](../packages/LambertW/index.html) package provides d, p,
    q, r functions as well as the first 4 central moments and a qqplot.
-   *Laplace (also called double exponential distribution) and asymetric
    Laplace distribution* : provided in
    [distr](../packages/distr/index.html),
    [lmomco](../packages/lmomco/index.html),
    [LaplacesDemon](../packages/LaplacesDemon/index.html),
    [VGAM](../packages/VGAM/index.html),
    [sgt](../packages/sgt/index.html),
    [extraDistr](../packages/extraDistr/index.html) and
    [HyperbolicDist](../packages/HyperbolicDist/index.html) packages.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides the
    Laplace distribution parametrized by the precision parameter as well
    as the skew Laplace distribution. Asymetric Laplace distribution is
    implemented in [ald](../packages/ald/index.html). A fast random
    generator is available for the Laplace distribution is implemented
    in [Runuran](../packages/Runuran/index.html) as well as the density
    function. [smoothmest](../packages/smoothmest/index.html) implements
    the density and the random generator. The skew Laplace distribution
    is available in [sgt](../packages/sgt/index.html).
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides the
    log-Laplace distribution.
-   *LASSO distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Linear failure rate distribution* : provided in
    [reliaR](../packages/reliaR/index.html).
-   *Loglog distribution* : provided in
    [reliaR](../packages/reliaR/index.html)
-   *Lomax distribution* : see beta.
-   *Logistic distribution and its extensions* : Base R provides the d,
    p, q, r functions for this distribution (see above).
    [actuar](../packages/actuar/index.html) and
    [VGAM](../packages/VGAM/index.html) provide d, p, q, r functions for
    the log logistic (also called Fisk), the paralogistic and the
    inverse paralogistic distributions.
    [FAdist](../packages/FAdist/index.html) the log-logistic
    distribution with two and three parameters. The generalized logistic
    distribution (Type I, also known as skew-logistic distribution) is
    provided in [lmomco](../packages/lmomco/index.html),
    [sld](../packages/sld/index.html), [SCI](../packages/SCI/index.html)
    and [glogis](../packages/glogis/index.html). Finally,\
      ---------------------- ----------------------------------------- -------------------- -----------------------
      *Distribution name*    *Packages*                                *Functions*          *Distribution suffix*
      Logistic               stats                                     d, p, q, r           `logis`
      Logistic               [actuar](../packages/actuar/index.html)   m, mgf               `logis`
      Log logistic           [actuar](../packages/actuar/index.html)   d, p, q, r, m, lev   `llogis`
      Log logistic           [VGAM](../packages/VGAM/index.html)       d, p, q, r           `fisk`
      Log logistic           [FAdist](../packages/FAdist/index.html)   d, p, q, r           `llog, llog3`
      Paralogistic           [actuar](../packages/actuar/index.html)   d, p, q, r, m, lev   `paralogis`
      Paralogistic           [VGAM](../packages/VGAM/index.html)       d, p, q, r           `paralogistic`
      Inv. paralogistic      [actuar](../packages/actuar/index.html)   d, p, q, r, m, lev   `invparalogis`
      Inv. paralogistic      [VGAM](../packages/VGAM/index.html)       d, p, q, r           `inv.paralogistic`
      Generalized logistic   [glogis](../packages/glogis/index.html)   d, p, q, r           `glogis`
      Generalized logistic   [SCI](../packages/SCI/index.html)         d, p, q              `genlog`
      Generalized logistic   [lmomco](../packages/lmomco/index.html)   d, p, q, r           `glo`
      Generalized logistic   [sld](../packages/sld/index.html)         d, p, q, r           `sl`
      ---------------------- ----------------------------------------- -------------------- -----------------------

      : Summary for Logistic-related distributions

    \
-   *Logit-normal distribution* : provided in
    [logitnorm](../packages/logitnorm/index.html).
-   *Log-normal distribution and its extensions* : The log normal
    distribution is implemented in Base R (see above) and
    [poweRlaw](../packages/poweRlaw/index.html). The log normal
    distribution parametrized by its mean and its coefficient of
    variation is also provided in
    [EnvStats](../packages/EnvStats/index.html).
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides the
    lognormal parametrized by the precision parameter. The truncated
    lognormal distribution is provided in
    [EnvStats](../packages/EnvStats/index.html) with two possible
    parametrizations as well as in
    [ReIns](../packages/ReIns/index.html). The 3-parameter lognormal
    distribution is available in
    [lmomco](../packages/lmomco/index.html),
    [EnvStats](../packages/EnvStats/index.html) and
    [FAdist](../packages/FAdist/index.html). The package
    [loglognorm](../packages/loglognorm/index.html) implements d, p, q,
    r functions for the double lognormal distribution, as well as the
    raw moment, the expected value and the variance functions.
    [EnvStats](../packages/EnvStats/index.html) provides d, p, q, r
    functions for the zero-modified lognormal distribution with two
    possible parametrizations. [rpgm](../packages/rpgm/index.html)
    provides a fast random number generator for the lognormal
    distribution.
-   *Makeham distribution* : provided in
    [VGAM](../packages/VGAM/index.html) and
-   *Maxwell distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Minimax distribution* : provided in
    [minimax](../packages/minimax/index.html).
-   *Mittag-Leffler distribution* : d, p, q, r functions provided in
    [MittagLeffleR](../packages/MittagLeffleR/index.html).
-   *Nakagami distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Pareto distribution* : d, p, q, r functions are implemented in
    [VGAM](../packages/VGAM/index.html) for the Pareto distribution type
    IV (which includes Burr\'s distribution, Pareto type III, Pareto
    type II (also called the lomax distribution) and Pareto type I) and
    the (upper/lower) truncated Pareto distribution. In an actuarial
    context, [actuar](../packages/actuar/index.html) provides d, p, q, r
    functions as well as moments and limited expected values for the
    Pareto I and II, the inverse Pareto, the \'generalized pareto\'
    distributions, the Burr and the inverse Burr distributions, all
    special cases of the transformed beta II distribution. A fast random
    generator for the Burr and the Pareto II distribution is implemented
    in [Runuran](../packages/Runuran/index.html) as well as the density.
    [EnvStats](../packages/EnvStats/index.html) and
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, p,
    q, r functions for Pareto I distribution.
    [extremefit](../packages/extremefit/index.html) provides the Burr,
    the Pareto II, mixture of Pareto I distributions and a composite
    distribution of two Pareto I distributions.
    [lmomco](../packages/lmomco/index.html),
    [evd](../packages/evd/index.html),
    [fExtremes](../packages/fExtremes/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [QRM](../packages/QRM/index.html),
    [Renext](../packages/Renext/index.html),
    [revdbayes](../packages/revdbayes/index.html),
    [FAdist](../packages/FAdist/index.html),
    [LaplacesDemon](../packages/LaplacesDemon/index.html) and
    [evir](../packages/evir/index.html) packages implement the
    Generalized Pareto Distribution (from Extreme Value Theory), which
    is depending the shape parameter\'s value a Pareto II distribution,
    a shifted exponential distribution or a generalized beta I
    distribution.
    [ParetoPosStable](../packages/ParetoPosStable/index.html) implements
    the Pareto positive stable distribution. The extended Pareto
    distribution is implemented in [RTDE](../packages/RTDE/index.html)
    and the shifted truncated (to unit interval) Pareto is implemented
    in [mbbefd](../packages/mbbefd/index.html).
    [ReIns](../packages/ReIns/index.html) provides Burr, extended
    Pareto, generalized Pareto, Pareto 1 distributions and their
    truncated version.\
      -------------------------- ------------------------------------------------------- -------------------- -----------------------
      *Distribution name*        *Packages*                                              *Functions*          *Distribution suffix*
      Pareto I                   [VGAM](../packages/VGAM/index.html)                     d, p, q, r           `paretoI`
      Pareto I                   [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `pareto1`
      Pareto I                   [EnvStats](../packages/EnvStats/index.html)             d, p, q, r           `pareto`
      Pareto I                   [extraDistr](../packages/extraDistr/index.html)         d, p, q, r           `pareto`
      Pareto I                   [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `pareto`
      Pareto I                   [LaplacesDemon](../packages/LaplacesDemon/index.html)   d, p, q, r           `pareto`
      Trunc. Pareto I            [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `tpareto`
      Pareto II                  [VGAM](../packages/VGAM/index.html)                     d, p, q, r           `paretoII`
      Pareto II                  [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `pareto, pareto2`
      Pareto II                  [Runuran](../packages/Runuran/index.html)               d, r                 `pareto`
      Pareto II                  [extraDistr](../packages/extraDistr/index.html)         d, p, q, h           `lomax`
      Pareto II                  [extremefit](../packages/extremefit/index.html)         d, p, q, h           `pareto`
      Pareto II                  [Renext](../packages/Renext/index.html)                 d, p, q, r           `lomax`
      Pareto III                 [VGAM](../packages/VGAM/index.html)                     d, p, q, r           `paretoIII`
      Pareto IV                  [VGAM](../packages/VGAM/index.html)                     d, p, q, r           `paretoIV`
      Inverse Pareto             [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `invpareto`
      Extended Pareto            [RTDE](../packages/RTDE/index.html)                     d, p, q, r           `EPD`
      Extended Pareto            [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `epd`
      Shift. trunc. Pareto       [mbbefd](../packages/mbbefd/index.html)                 d, p, q, r, m, ec    `stpareto`
      Gen. Pareto (actuarial)    [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `genpareto`
      Gen. Pareto (EVT)          [lmomco](../packages/lmomco/index.html)                 d, p, q, r           `gpa`
      Gen. Pareto (EVT)          [evd](../packages/evd/index.html)                       d, p, q, r           `gpd`
      Gen. Pareto (EVT)          [fExtremes](../packages/fExtremes/index.html)           d, p, q, r           `gpd`
      Gen. Pareto (EVT)          [evir](../packages/evir/index.html)                     d, p, q, r           `gpd`
      Gen. Pareto (EVT)          [extraDistr](../packages/extraDistr/index.html)         d, p, q, r           `gpd`
      Gen. Pareto (EVT)          [QRM](../packages/QRM/index.html)                       d, p, q, r           `GPD`
      Gen. Pareto (EVT)          [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `gpd`
      Gen. Pareto (EVT)          [LaplacesDemon](../packages/LaplacesDemon/index.html)   d, r                 `gpd`
      Trunc. Gen. Pareto (EVT)   [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `tgpd`
      Gen. Pareto (EVT)          [revdbayes](../packages/revdbayes/index.html)           d, p, q, r           `gp`
      Gen. Pareto (EVT)          [Renext](../packages/Renext/index.html)                 d, p, q, r           `GPD`
      Burr                       [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `burr`
      Burr                       [extremefit](../packages/extremefit/index.html)         d, p, q, r           `burr`
      Burr                       [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `burr`
      Trunc. Burr                [ReIns](../packages/ReIns/index.html)                   d, p, q, r           `tburr`
      Inverse Burr               [actuar](../packages/actuar/index.html)                 d, p, q, r, m, lev   `invburr`
      -------------------------- ------------------------------------------------------- -------------------- -----------------------

      : Summary for Pareto-related distributions

    \
-   *Pearson\'s distribution* : Pearson type III available in
    [lmomco](../packages/lmomco/index.html) and
    [FAdist](../packages/FAdist/index.html). A log-Pearson type III
    distribution is also available in
    [FAdist](../packages/FAdist/index.html).
    [PearsonDS](../packages/PearsonDS/index.html) provides the d, p, q,
    r functions as well as the first four moments for the Pearson
    distributions: types I, II, III, IV, V, VI, VII.
-   *Pearson\'s Rho distribution* : provided in
    [SuppDists](../packages/SuppDists/index.html).
-   *Perks distribution* : provided in
    [VGAM](../packages/VGAM/index.html).
-   *Planck\'s distribution* : a random generator is available in
    [Runuran](../packages/Runuran/index.html).
-   *Phase-type distribution* : provided in
    [actuar](../packages/actuar/index.html), see also
    [PhaseType](../packages/PhaseType/index.html) for inference.
-   *Poisson subordinated distributions* : provided in
    [LIHNPSD](../packages/LIHNPSD/index.html) (d, p, q, r, m functions).
-   *Power distribution* : [reliaR](../packages/reliaR/index.html) and
    [poweRlaw](../packages/poweRlaw/index.html) implement the
    exponential power distribution.
-   *Proportion distribution* : this is the distribution for the
    difference between two independent beta distributions. d, p, q, r
    functions in [tolerance](../packages/tolerance/index.html).
-   *Rayleigh distribution* : provided in packages
    [VGAM](../packages/VGAM/index.html),
    [extraDistr](../packages/extraDistr/index.html) and
    [lmomco](../packages/lmomco/index.html). Generalized and logistic
    Rayleigh distributions are available in
    [reliaR](../packages/reliaR/index.html).
-   *Response time distribution* :
    [rtdists](../packages/rtdists/index.html) provides d, p, q, r
    functions for the (Ratcliff) diffusion distribution and for the
    linear ballistic accumulator (LBA) with different underlying
    drift-distributions (Normal, Gamma, Frechet, and log-normal).
-   *Rice distribution* : provided in
    [VGAM](../packages/VGAM/index.html) and
    [lmomco](../packages/lmomco/index.html).
-   *Singh-Maddala distribution* : see beta.
-   *Slash distribution* : provided in
    [lmomco](../packages/lmomco/index.html),
    [extraDistr](../packages/extraDistr/index.html) and
    [VGAM](../packages/VGAM/index.html).
-   *Spearman\'s Rho distribution* : provided in
    [SuppDists](../packages/SuppDists/index.html).
-   *Stable distribution* : d, p, q, r functions are available in
    [fBasics](../packages/fBasics/index.html) and
    [stabledist](../packages/stabledist/index.html), the functions use
    the approach of J.P. Nolan for general stable distributions.
    [MixedTS](../packages/MixedTS/index.html) provides mixed tempered
    stable distribution (d, p, q, r).
    [FMStable](../packages/FMStable/index.html) provides (d, p, q) the
    extremal or maximally skew stable and the finite moment log stable
    distributions.
-   *Student distribution and its extensions* : Base R provides the d,
    p, q, r functions for Student and non central Student distribution
    (see above). [extraDistr](../packages/extraDistr/index.html) and
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides the
    Student distribution with location and scale parameters.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, p,
    q, r functions for the Half-Student distribution.
    [sadists](../packages/sadists/index.html) implements Gram Charlier,
    Edgeworth and Cornish-Fisher approximations for doubly non central
    Student distribution for computing d, p, q, r functions. The skewed
    Student distribution is provided in
    [skewt](../packages/skewt/index.html),
    [sn](../packages/sn/index.html) and
    [gamlss.dist](../packages/gamlss.dist/index.html) packages. The
    generalized skew distribution is provided in
    [sgt](../packages/sgt/index.html). d, p, q, r functions for the
    generalized t-distribution can be found in
    [gamlss.dist](../packages/gamlss.dist/index.html).
    [fBasics](../packages/fBasics/index.html) provides d, p, q, r
    functions for the skew and the generalized hyperbolic
    t-distribution. The L-moments of the Student t (3-parameter) are
    provided in [lmomco](../packages/lmomco/index.html).
    [rpgm](../packages/rpgm/index.html) provides a fast random number
    generator for the student distribution.\
      ----------------------------- --------------------------------------------------- ------------- ---------------------------
      *Distribution name*           *Packages*                                          *Functions*   *Distribution suffix*
      Student                       stats                                               d, p, q, r    `t`
      Student with loc. and scal.   extraDistr                                          d, p, q, r    `lst`
      Student with loc. and scal.   LaplacesDemon                                       d, p, q, r    `st`
      Doubly non central St.        [sadists](../packages/sadists/index.html)           d, p, q, r    `dnt`
      Skew Student                  [skewt](../packages/skewt/index.html)               d, p, q, r    `skt`
      Skew Student                  [sn](../packages/sn/index.html)                     d, p, q, r    `st`
      Skew St. Type 1-5             [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r    `ST1, ST2, ST3, ST4, ST5`
      Gen. Student                  [gamlss.dist](../packages/gamlss.dist/index.html)   d, p, q, r    `GT`
      Gen. Hyp. Student             [fBasics](../packages/fBasics/index.html)           d, p, q, r    `ght`
      Skew Gen. Student             [sgt](../packages/sgt/index.html)                   d, p, q, r    `sgt`
      ----------------------------- --------------------------------------------------- ------------- ---------------------------

      : Summary for Student-related distributions

    \
-   *Triangle/trapezoidal distribution* : packages
    [triangle](../packages/triangle/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [mc2d](../packages/mc2d/index.html),
    [EnvStats](../packages/EnvStats/index.html) and
    [VGAM](../packages/VGAM/index.html) provide d, p, q, r functions for
    the triangle or triangular distribution, while the package
    [trapezoid](../packages/trapezoid/index.html) provides d, p, q, r
    functions for the Generalized Trapezoidal Distribution. A fast
    random generator is available for the triangle distribution is
    implemented in [Runuran](../packages/Runuran/index.html) as well as
    the density function.
-   *Tsallis or q-Exponential distribution* :
    [tsallisqexp](../packages/tsallisqexp/index.html) provides d, p, q,
    r functions for two parametrizations of the Tsallis distribution and
    also implements a left-censored version.
-   *Tweedie distribution* : the Tweedie distribution is implemented in
    package [tweedie](../packages/tweedie/index.html). Let us note that
    the Tweedie distribution is not necessarily continuous, a special
    case of it is the Poisson distribution.
-   *Uniform distribution* : d, p, q, r functions are of course provided
    in R. See section RNG for random number generation topics.
    [HI](../packages/HI/index.html) generates uniformly random points on
    a bounded convex set, in particular the unit ball.
    [KScorrect](../packages/KScorrect/index.html) provides d, p, q, r
    functions for the log-uniform distribution.
-   *Upsilon distribution* : [sadists](../packages/sadists/index.html)
    implements Gram Charlier, Edgeworth and Cornish-Fisher
    approximations for Upsilon distribution for computing d, p, q, r
    functions.
-   *Wakeby distribution* : A 5-parameter Wakeby is provided in
    [lmomco](../packages/lmomco/index.html).
-   *Weibull distribution and its extensions* : Base R provides the d,
    p, q, r functions for this distribution (see above). The inverse
    Weibull is provided in [actuar](../packages/actuar/index.html)
    package and also the moments and the limited expected value for both
    the raw and the inverse Weibull distribution.
    [FAdist](../packages/FAdist/index.html) implements the
    three-parameter Weibull distribution, while
    [reliaR](../packages/reliaR/index.html) implements the exponential
    Weibull, the flexible Weibull, the generalized power weibull, the
    Marshall-Olkin Extended Weibull and the Weibull extension
    distributions. Furthermore, [lmomco](../packages/lmomco/index.html)
    implements the Weibull distribution while
    [evd](../packages/evd/index.html) implements the reverse Weibull
    distribution. The reverse generalized extreme value distribution are
    provided in [gamlss.dist](../packages/gamlss.dist/index.html) (d, p,
    q, r) and the shifted left truncated Weibull distribution is
    provided in [Renext](../packages/Renext/index.html). The right
    truncated Weibull is provided in
    [ReIns](../packages/ReIns/index.html).

## [Continuous multivariate distributions:]{#MultivariateContinuous}

-   *Multivariate Cauchy distribution* : [sn](../packages/sn/index.html)
    provide d, p, r functions for the multivariate skew Cauchy
    distribution, while
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    functions for the multivariate Cauchy distribution parametrized
    either by sigma, by the Cholesky decomposition of sigma, by the
    precision matrix omega or by the Cholesky decomposition of omega.
-   *Multinomial Dirichlet distribution* : functions d, r are provided
    in [MCMCpack](../packages/MCMCpack/index.html),
    [mc2d](../packages/mc2d/index.html),
    [dirmult](../packages/dirmult/index.html),
    [extraDistr](../packages/extraDistr/index.html) and
    [bayesm](../packages/bayesm/index.html).
-   *Multivariate exponential distribution* : while
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    functions for the multivariate power exponential distribution
    parametrized either by sigma, or by the Cholesky decomposition of
    sigma.
-   *Multivariate Gaussian (or normal) distribution* : The multivariate
    Gaussian distribution is provided in the packages
    [mvtnorm](../packages/mvtnorm/index.html) (d, r),
    [mnormt](../packages/mnormt/index.html) (d, p, r),
    [Compositional](../packages/Compositional/index.html) (r) .
    [mvprpb](../packages/mvprpb/index.html) computes the orthant
    probability of the multivariate Gaussian distribution.
    [symmoments](../packages/symmoments/index.html) computes central and
    non-central moments of the multivariate Gaussian distribution.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    functions for the multivariate normal distribution parametrized
    either by sigma, by the Cholesky decomposition of sigma, by the
    precision matrix omega or by the Cholesky decomposition of omega.
    Futhermore, [tmvtnorm](../packages/tmvtnorm/index.html) implements
    the truncated multivariate normal distribution.
    [sparseMVN](../packages/sparseMVN/index.html) implements very fast
    algorithms to compute the density and generate random variates of a
    multivariate normal distribution for which the covariance matrix or
    precision matrix is sparse.
    [cmvnorm](../packages/cmvnorm/index.html) implements the complex
    multivariate normal distribution (d, r). Finally,
    [condMVNorm](../packages/condMVNorm/index.html) implements d, p, r
    functions for the conditional multivariate normal distribution.
    Furthermore, [sn](../packages/sn/index.html) besides providing
    facilities for their distribution functions,
    [sn](../packages/sn/index.html) allows the creation of S4 objects
    which encapsulate these distributions and provide facilities for
    plotting, summary, marginalization, conditioning, affine
    transformations of these S4 objects.
    [mnormpow](../packages/mnormpow/index.html) computes the expected
    product of the components of a multivariate Gaussian vector.
    [Compositional](../packages/Compositional/index.html) provides
    random generator for the multivariate normal distribution on the
    simplex and multivariate skew normal distribution on the simplex.
-   *Multivariate generalized hyperbolic distribution* :
    [QRM](../packages/QRM/index.html) provides d, r functions of the
    standard and the symmetric multivariate generalized hyperbolic
    distribution. [ghyp](../packages/ghyp/index.html) provides d, p, r
    functions of the standard multivariate generalized hyperbolic
    distribution.
-   *Multivariate generalized extreme value distribution* : Both
    bivariate and multivariate Extreme Value distributions as well as
    order/maxima/minima distributions are implemented in
    [evd](../packages/evd/index.html) (d, p, r).
-   *Multivariate Laplace distribution* :
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    functions for the multivariate Laplace distribution parametrized
    either by sigma, or by the Cholesky decomposition of sigma.
-   *Multivariate logistic distribution* :
    [VGAM](../packages/VGAM/index.html) package implements the bivariate
    logistic distribution.
-   *Multivariate Pareto distribution* :
    [mgpd](../packages/mgpd/index.html) provides the density for the
    multivariate generalized Pareto distribution of type II, while
    [evd](../packages/evd/index.html) provides the density for type I.
-   *Multivariate Stable distribution* : not yet implemented?
-   *Multivariate Student distribution* : The multivariate Student
    distribution is provided in the packages
    [mvtnorm](../packages/mvtnorm/index.html) (d, r),
    [mnormt](../packages/mnormt/index.html) (d, p, r),
    [Compositional](../packages/Compositional/index.html) (r),
    [QRM](../packages/QRM/index.html) (d, r). First two moments (m) and
    sampling (r) of the Truncated Multivariate t Distribution are
    provided in [TTmoment](../packages/TTmoment/index.html).
    [sn](../packages/sn/index.html) provides d, p, r functions for the
    multivariate skew t distribution.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides d, r
    functions for the multivariate Student distribution parametrized
    either by sigma, by the Cholesky decomposition of sigma, by the
    precision matrix omega or by the Cholesky decomposition of omega.

## [Mixed-type distributions:]{#MixedType}

-   *Maxwell-Boltzmann-Bose-Einstein-Fermie-Dirac (MBBEFD) distribution*
    : provided in [mbbefd](../packages/mbbefd/index.html).
-   *Mixed ordinal and normal distribution* : provided in
    [OrdNor](../packages/OrdNor/index.html).
-   *One-inflated distributions* : a generic distribution as well as
    special cases (OI-beta, OI-uniform, OI-GB1, OI-Pareto) are provided
    in [mbbefd](../packages/mbbefd/index.html). The zero and one
    inflated beta distribution can be found in
    [gamlss.dist](../packages/gamlss.dist/index.html).
-   *Zero-modified distributions* :
    [EnvStats](../packages/EnvStats/index.html) provides the
    zero-modified normal distribution and the zero-modified lognormal
    distribution.

## [Mixture of probability laws (and composite):]{#Mixture}

-   *Bernoulli-dist mixture* : d, p, q, r functions for
    Bernoulli-exponential, Bernoulli-Gamma, Bernoulli-lognormal,
    Bernoulli-Weibull distributions are provided in
    [qmap](../packages/qmap/index.html).
-   *Cauchy-polynomial quantile mixture* : d, p, q, r functions are
    provided in [Lmoments](../packages/Lmoments/index.html).
-   *composite lognormal distribution* : d, p, q, r functions are
    provided in [CompLognormal](../packages/CompLognormal/index.html).
-   *Gaussian mixture* : Functions d, r are provided in
    [mixtools](../packages/mixtools/index.html),
    [bmixture](../packages/bmixture/index.html) package when dealing
    with finite mixture models.
    [nor1mix](../packages/nor1mix/index.html),
    [extraDistr](../packages/extraDistr/index.html),
    [mclust](../packages/mclust/index.html),
    [LaplacesDemon](../packages/LaplacesDemon/index.html),
    [KScorrect](../packages/KScorrect/index.html) provides d, p, r
    functions for Gaussian mixture.
    [EnvStats](../packages/EnvStats/index.html) provides d, p, q, r
    functions for mixture of two normal distributions.
-   *Gamma Poisson* : provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Gamma mixture* : Gamma shape mixtures are implemented (d, p, r) in
    the [GSM](../packages/GSM/index.html) package,
    [bmixture](../packages/bmixture/index.html) provides d, r functions.
-   *Generic mixtures* : there is an implementation via S4-class
    UnivarMixingDistribution in package
    [distr](../packages/distr/index.html).
    [gamlss.mx](../packages/gamlss.mx/index.html) uses the
    [gamlss.dist](../packages/gamlss.dist/index.html) package.
-   *Horseshoe distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Laplace mixture distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Log normal mixture* : d, p, q, r functions are provided in
    [EnvStats](../packages/EnvStats/index.html) with two possible
    parametrizations.
-   *Normal-polynomial quantile mixture* : d, p, q, r functions are
    provided in [Lmoments](../packages/Lmoments/index.html).
-   *Pareto distribution* :
    [extremefit](../packages/extremefit/index.html) implements the
    mixture of two Pareto I distributions.
-   *Poisson Binomial distribution* :
    [poibin](../packages/poibin/index.html) implements the Poisson
    Binomial distribution.
-   *Poisson lognormal distribution* :
    [poilog](../packages/poilog/index.html) implements the Poisson
    lognormal distribution.
-   *Poisson mixture* : provided in
    [extraDistr](../packages/extraDistr/index.html).
-   *Poisson-Tweedie exponential family models* : provided in
    [poistweedie](../packages/poistweedie/index.html).
-   *Student mixture* : The [AdMit](../packages/AdMit/index.html)
    package provides d, r functions for Student mixtures in the context
    of Adaptive Mixture of Student-t distributions.
    [MitISEM](../packages/MitISEM/index.html),
    [bmixture](../packages/bmixture/index.html) package also provide d,
    r functions for mixture of Student-t distributions.
-   *von Mises Fisher (or Langevin) mixture* : The
    [movMF](../packages/movMF/index.html) package provides d, r
    functions for finite von Mises Fisher mixtures.

## [Compound, composite, discretized, exponentiated and tranformation of distributions:]{#Transform}

-   *Absolute value or half distribution* : Half-Cauchy, half normal and
    half-student are implemented both in
    [extraDistr](../packages/extraDistr/index.html) and in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Composite distribution also known as spliced distribution* :
    Composite lognormal distributions provided in
    [CompLognormal](../packages/CompLognormal/index.html). Split-normal
    (also known as the two-piece normal distribution) not yet
    implemented. Split-student provided in package
    [dng](../packages/dng/index.html).
-   *Compound distribution* : d, p, q, r, m functions are implemented by
    [Compounding](../packages/Compounding/index.html) where the parent
    distribution is any continuous distribution and the compound
    distribution is any distribution among the list: binomial,
    binomial-Poisson, geometric, hypergeometric, hyper-Poisson, Katti
    type H1/H2, logarithmic, logarithmic-binomial, logarithmic-Poisson,
    negative binomial, Neyman type A/B/C, Pascal-Poisson, Poisson,
    Poisson-binomial, Poisson-Lindley, Poisson-Pascal, Polya Aeppli,
    Thomas, Waring, Yule.
-   *Discretized distribution* :
    [distcrete](../packages/distcrete/index.html) allows discretised
    versions of continuous distribution by mapping continuous values to
    an underlying discrete grid, based on a (uniform) frequency of
    discretisation, a valid discretisation point, and an integration
    range.
-   *G transformed distribution* : implemented in
    [Newdistns](../packages/Newdistns/index.html) which includes
    Marshall Olkin G distribution, exponentiated G distribution, beta G
    distribution, gamma G distribution, Kumaraswamy G distribution,
    generalized beta G distribution, beta extended G distribution, gamma
    G distribution, gamma uniform G distribution, beta exponential G
    distribution, Weibull G distribution, log gamma G1/G2 distribution,
    exponentiated generalized G distribution, exponentiated Kumaraswamy
    G distributions, geometric exponential Poisson G distribution,
    truncated-exponential skew-symmetric G distribution, modified beta G
    distribution, and exponentiated exponential Poisson G distribution.
-   *Truncated distribution* : A generic code snippet is available [in
    the JSS](https://www.jstatsoft.org/article/view/v016c02) .
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides a
    generic function. For a given distribution, look at the
    corresponding subsection above.

## [Moments, skewness, kurtosis and etc:]{#Moments}

-   *Empirical mean, standard deviation and variance* : base R provides
    `mean()`, `sd()`, `var()` functions to compute the mean, standard
    deviation and variance, respectively.
-   *Empirical skewness* : available in
    [agricolae](../packages/agricolae/index.html),
    [e1071](../packages/e1071/index.html),
    [GLDEX](../packages/GLDEX/index.html),
    [HyperbolicDist](../packages/HyperbolicDist/index.html),
    [modeest](../packages/modeest/index.html),
    [moments](../packages/moments/index.html),
    [TSA](../packages/TSA/index.html),
    [s20x](../packages/s20x/index.html),
    [DistributionUtils](../packages/DistributionUtils/index.html),
    [EnvStats](../packages/EnvStats/index.html),
    [rpgm](../packages/rpgm/index.html) packages.
-   *Empirical kurtosis* : available in
    [agricolae](../packages/agricolae/index.html),
    [DistributionUtils](../packages/DistributionUtils/index.html),
    [e1071](../packages/e1071/index.html),
    [EnvStats](../packages/EnvStats/index.html),
    [GLDEX](../packages/GLDEX/index.html),
    [HyperbolicDist](../packages/HyperbolicDist/index.html),
    [moments](../packages/moments/index.html),
    [TSA](../packages/TSA/index.html),
    [rpgm](../packages/rpgm/index.html) packages. The raw or centered
    moments are provided in [e1071](../packages/e1071/index.html),
    [moments](../packages/moments/index.html).
-   *Empirical L-moments* : L-moments are available in
    [lmom](../packages/lmom/index.html),
    [lmomco](../packages/lmomco/index.html),
    [Lmoments](../packages/Lmoments/index.html),
    [GLDEX](../packages/GLDEX/index.html),
    [EnvStats](../packages/EnvStats/index.html), trimmed L-moments are
    available in [lmomco](../packages/lmomco/index.html), and
    [Lmoments](../packages/Lmoments/index.html), right-censored
    L-moments are available in [lmomco](../packages/lmomco/index.html),
    and cumulants in [GLDEX](../packages/GLDEX/index.html).
-   *Empirical probability weighted moments* : Probability weighted
    moments are available in
    [EnvStats](../packages/EnvStats/index.html).
-   *Mode estimation* : Package
    [modeest](../packages/modeest/index.html) provides mode estimation
    for various distributions.
-   *Order statistics* : Distribution function of the jth order
    statistic can be obtained with base R functions.
    [ORDER2PARENT](../packages/ORDER2PARENT/index.html) transforms
    distribution function of order statistics to its parent distribution
    function.
-   *Theoretical moments* : The [actuar](../packages/actuar/index.html)
    package implements raw moments, limited expected values and moment
    generating function for base R distributions.
    [HyperbolicDist](../packages/HyperbolicDist/index.html) provides the
    mean, variance, skewness, kurtosis, mode, raw and centered moments
    for the hyperbolic, the generalized hyperbolic and the generalized
    inverse Gaussian distributions.
    [GLDEX](../packages/GLDEX/index.html) also provides the mean,
    variance, skewness, kurtosis of generalized Lambda distribution.
    [mvrtn](../packages/mvrtn/index.html) provides mean, variance for
    left/right truncated normal distributions.
    [lmomco](../packages/lmomco/index.html) provides L-moments (L),
    trimmed L-moments (TL), and right-censored \[RC\] for the following
    distributions: Asymmetric Exponential Power (L), Cauchy (TL), Eta-Mu
    (L), Exponential (L), Gamma (L), Generalized Extreme Value (L),
    Generalized Lambda (L and TL), Generalized Logistic (L), Generalized
    Normal (L), Generalized Pareto (L\[RC\] and TL), Govindarajulu (L),
    Gumbel (L), Kappa (L), Kappa-Mu (L), Kumaraswamy (L), Laplace (L),
    Normal (L), 3-parameter log-Normal (L), Pearson Type III (L),
    Rayleigh (L), Reverse Gumbel (L\[RC\]), Rice/Rician (L), Slash (TL),
    3-parameter Student T (L), Truncated Exponential (L), Wakeby (L),
    and Weibull (L). Multivariate L-moments (L-comoments).

## [Random matrices:]{#Matrix}

-   *Huang-Wan distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Inverse matrix gamma distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Inverse Wishart distribution* :
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides
    inverse Wishart distribution parametrized either by Sigma or by its
    Cholesky decomposition.
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides the
    scaled inverse Wishart distribution.
-   *Marcenko-Pastur distribution* : provided in
    [RMTstat](../packages/RMTstat/index.html),
    [MCMCpack](../packages/MCMCpack/index.html) and
    [bayesm](../packages/bayesm/index.html).
-   *Matrix gamma distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Matrix normal distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Normal Inverse Wishart distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Normal Wishart distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Tracy-Widom distribution* : provided in
    [RMTstat](../packages/RMTstat/index.html),
    [MCMCpack](../packages/MCMCpack/index.html) and
    [bayesm](../packages/bayesm/index.html): supported beta values are 1
    (Gaussian Orthogonal Ensemble), 2 (Gaussian Unitary Ensemble), and 4
    (Gaussian Symplectic Ensemble).
-   *Spiked Wishart Maximum Eigenvalue Distribution* : provided in
    [RMTstat](../packages/RMTstat/index.html),
    [MCMCpack](../packages/MCMCpack/index.html) and
    [bayesm](../packages/bayesm/index.html).
-   *Wishart distributions* : Base R provides the r function for the
    Wishart distribution. the d, r functions are provided in
    [MCMCpack](../packages/MCMCpack/index.html),
    [RMTstat](../packages/RMTstat/index.html) and
    [bayesm](../packages/bayesm/index.html).
    [LaplacesDemon](../packages/LaplacesDemon/index.html) provides
    Wishart distribution parametrized either by Sigma or by its Cholesky
    decomposition.
-   *White Wishart Maximum Eigenvalue Distribution* : provided in
    [RMTstat](../packages/RMTstat/index.html),
    [MCMCpack](../packages/MCMCpack/index.html) and
    [bayesm](../packages/bayesm/index.html).
-   *Yang-Berger distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).
-   *Zellner distribution* : provided in
    [LaplacesDemon](../packages/LaplacesDemon/index.html).

## [Copulas:]{#Copulas}

-   *Unified approaches* : The packages
    [fCopulae](../packages/fCopulae/index.html),
    [copula](../packages/copula/index.html), and
    [copBasic](../packages/copBasic/index.html) provide a lot of general
    functionality for copulas. Although lacking support for many
    existing copulas themselves,
    [copBasic](../packages/copBasic/index.html) is primarily oriented
    around utility functions for the general mathematics of copulas as
    described in the well known introduction to copulas by Nelsen.
-   *Archimedean copulas* : [gumbel](../packages/gumbel/index.html) is a
    standalone package for the Gumbel copula
    [fCopulae](../packages/fCopulae/index.html) implements the 22
    Archimedean copulas of Nelsen (1998, *Introduction to Copulas* ,
    Springer-Verlag) including Gumbel, Frank, Clayton, and
    Ali-Mikhail-Haq. [VGAM](../packages/VGAM/index.html) provides
    Ali-Mikhail-Haq, Clayton, Frank, Frechet copulas.
    [copula](../packages/copula/index.html) provides Ali-Mikhail-Haq,
    Clayton, Frank, Gumbel and Joe copulas. The Frank bivariate
    distribution is available in [RTDE](../packages/RTDE/index.html).
    [CDVine](../packages/CDVine/index.html) and
    [VineCopula](../packages/VineCopula/index.html) provide Clayton,
    Gumbel, Frank, Joe, BB1, BB6, BB7 and BB8 copulas. Nested
    Archimedean copulas are available in the
    [HAC](../packages/HAC/index.html) package. Generalized Archimedean
    copulas are implemented in the [fgac](../packages/fgac/index.html)
    package. [BivarP](../packages/BivarP/index.html) provides cdf, pdf
    and survival function for Clayton, Gumbel and Frank copula.
    [QRM](../packages/QRM/index.html) provides pdf and random generator
    for Clayton, Gumbel, Frank, BB9 copula.
-   *Blomqvist copula* : provided in
    [copBasic](../packages/copBasic/index.html).
-   *Composition of copula* :
    [copBasic](../packages/copBasic/index.html) provides functions for
    composition of a single symmetric copula and composition of two
    copulas.
-   *Cubic copula* : Not yet implemented?
-   *Dirichlet copula* : Not yet implemented?
-   *Empirical copula* : provided in
    [copBasic](../packages/copBasic/index.html) and in
    [HAC](../packages/HAC/index.html).
    [GenOrd](../packages/GenOrd/index.html) provides sampling function
    for multivariate discrete random vectors with a specified
    correlation matrix.
-   *Elliptical copulas* : Gaussian, Student and Cauchy copulas are
    implemented in [fCopulae](../packages/fCopulae/index.html) for the
    bivariate cases. [copula](../packages/copula/index.html),
    [CDVine](../packages/CDVine/index.html),
    [VGAM](../packages/VGAM/index.html),
    [VineCopula](../packages/VineCopula/index.html) provide the Gaussian
    and the Student copulas. [QRM](../packages/QRM/index.html) provides
    pdf and random generator for Gaussian, Student copulas.
-   *Extreme value copulas* :
    [fCopulae](../packages/fCopulae/index.html) provides the following
    copulas Gumbel, Galambos, Husler-Reiss, Tawn, or BB5.
    [copula](../packages/copula/index.html) implements Gumbel, Galambos
    and Husler-Reiss.
-   *Eyraud-Farlie-Gumbel-Morgenstern copula* : provided in
    [VGAM](../packages/VGAM/index.html),
    [RTDE](../packages/RTDE/index.html), and
    [copula](../packages/copula/index.html).
-   *Mardia copula* : Not yet implemented?
-   *Nested copulas* : arbitrary nested versions of copulas can be
    implemented in [copula](../packages/copula/index.html).
-   *Plackett* : provided in [VGAM](../packages/VGAM/index.html),
    [copBasic](../packages/copBasic/index.html) and
    [copula](../packages/copula/index.html).
-   *Vine copulas* : Packages [CDVine](../packages/CDVine/index.html),
    [vines](../packages/vines/index.html) provide functions for C- and
    D-vine copulas and [VineCopula](../packages/VineCopula/index.html)
    for general R-vine copulas.

## [Random number generators (RNG):]{#Random}

-   *Basic functionality* : R provides several random number generators
    (RNGs). The random seed can be provided via `set.seed` and the kind
    of RNG can be specified using `RNGkind`. The default RNG is the
    Mersenne-Twister algorithm. Other generators include Wichmann-Hill,
    Marsaglia-Multicarry, Super-Duper, Knuth-TAOCP, Knuth-TAOCP-2002, as
    well as user-supplied RNGs. For normal random numbers, the following
    algorithms are available: Kinderman-Ramage, Ahrens-Dieter,
    Box-Muller, Inversion (default). In addition to the tools above,
    [setRNG](../packages/setRNG/index.html) provides an easy way to set,
    retain information about the setting, and reset the RNG.
-   *Pseudo-randomness* :
    [RDieHarder](../packages/RDieHarder/index.html) offers several dozen
    new RNGs from the GNU GSL.
    [randtoolbox](../packages/randtoolbox/index.html) provides more
    recent RNGs such as SF Mersenne-Twister and WELL, which are
    generators of Mersenne Twister type, but with improved quality
    parameters. [rngwell19937](../packages/rngwell19937/index.html)
    provides one of the WELL generators with 53 bit resolution of the
    output and allows seeding by a vector of integers of arbitrary
    length. [randaes](../packages/randaes/index.html) provides the
    deterministic part of the Fortuna cryptographic pseudorandom number
    generator (AES). [SuppDists](../packages/SuppDists/index.html)
    implements two RNGs of G. Marsaglia.
    -   Support for several independent streams:
        [rstream](../packages/rstream/index.html) focuses on multiple
        independent streams of random numbers from different sources (in
        an object oriented approach).
    -   For non-uniform generation, the
        [Runuran](../packages/Runuran/index.html) package interfaces to
        the UNU.RAN library for universal non-uniform generation as well
        as customised distributions based on polynomial interpolation of
        the inverse cumulative distribution function.
        [rust](../packages/rust/index.html) performs non-uniform random
        variate generation from unimodal (low-dimensional) multivariate
        continuous distributions, using the generalized
        ratio-of-uniforms method.
    -   [kernelboot](../packages/kernelboot/index.html) provides
        functions for random generation from univariate and multivariate
        kernel densities (in particular multivariate Gaussian kernels).
-   *Quasi-randomness* : The
    [randtoolbox](../packages/randtoolbox/index.html) provides the
    following quasi random sequences: the Sobol sequence, the Halton
    (hence Van Der Corput) sequence and the Torus sequence (also known
    as Kronecker sequence). [lhs](../packages/lhs/index.html) and
    [mc2d](../packages/mc2d/index.html) packages implement the latin
    hypercube sampling, an hybrid quasi/pseudo random method.
    [sfsmisc](../packages/sfsmisc/index.html) also provides the Halton
    sequences.
-   *True randomness* : The [random](../packages/random/index.html)
    package provides several functions that access the true random
    number service at [random.org](http://random.org) .
-   *RNG tests* : [RDieHarder](../packages/RDieHarder/index.html) offers
    numerous tests of RNGs based on a reimplementation and extension of
    Marsaglia\'s DieHarder battery.
    [randtoolbox](../packages/randtoolbox/index.html) provides basic RNG
    tests.
-   *Parallel computing* : Random-number generators for parallel
    computing are available via the
    [rlecuyer](../packages/rlecuyer/index.html) package. See the
    [HighPerformanceComputing](HighPerformanceComputing.html) task view
    for more details.

## [Miscellaneous:]{#Misc}

-   *Computation* :
    -   *Approximation of d, p, q, r functions* :
        [PDQutils](../packages/PDQutils/index.html) provides tools for
        computing the density, cumulative distribution, and quantile
        functions of a distribution when the cumulants or moments are
        given, using the classical Gram Charlier, Edgeworth and
        Cornish-Fisher approximations.
        [sadists](../packages/sadists/index.html) is a showcase for
        PDQutils, providing density, cumulative distribution, quantile,
        and random generation for the doubly non-central t, doubly
        non-central F, K-prime, Lambda-prime, Upsilon, and sum of
        (non-central) chi-squares to powers distributions.
    -   For non-uniform generation, see the
        [Runuran](../packages/Runuran/index.html) above.
    -   *Benchmark* : A set of 28 densities suitable for comparing
        nonparametric density estimators in simulation studies can be
        found in the [benchden](../packages/benchden/index.html)
        package. The densities vary greatly in degree of smoothness,
        number of modes and other properties. The package provides d,p,q
        and r functions.
-   *Non parametric models* :
    -   *Binned Empirical distributions* : The
        [HistogramTools](../packages/HistogramTools/index.html) package
        provides a number of methods for manipulating empirical data
        that has been binned into histogram form, including: (1) the
        empirical cumulative distribution function, (2) the empirical
        quantile, and (3) information loss metrics associated with
        binning.
    -   *Empirical distribution* : Base R provides functions for
        univariate analysis: (1) the empirical density (see
        density()), (2) the empirical cumulative distribution function
        (see ecdf()), (3) the empirical quantile (see quantile())
        and (4) random sampling (see sample()).
        [mded](../packages/mded/index.html) provides a function for
        measuring the difference between two independent or
        non-independent empirical distributions and returning a
        significance level of the difference.
    -   *Non Parametric distributions* :
        [spd](../packages/spd/index.html) provides the Semi Parametric
        Piecewise Distribution, while
        [fBasics](../packages/fBasics/index.html) implements spline
        smoothed distributions.
-   *Hierarchical models* : Distributions whose some parameters are no
    longer constant but random according to a particular distribution.
    [VGAM](../packages/VGAM/index.html) provides a lot of hierarchical
    models: beta/binomial, beta/geometric and beta/normal distributions.
    [bayesm](../packages/bayesm/index.html) implements: binary logit,
    linear, multivariate logit and negative binomial models. Furthermore
    [LearnBayes](../packages/LearnBayes/index.html) and
    [MCMCpack](../packages/MCMCpack/index.html) provides poisson/gamma,
    beta/binomial, normal/normal and multinomial/Dirichlet models.
-   *Distribution handling* :
    -   *Object-orientation* : General discrete and continuous
        distributions are implemented in package
        [distr](../packages/distr/index.html) respectively via S4-class
        DiscreteDistribution and AbscontDistribution providing the
        classic d, p, q and r functions.
        [distrEx](../packages/distrEx/index.html) extends available
        distributions to multivariate and conditional distributions as
        well as methods to compute useful statistics (expectation,
        variance,\...) and distances between distributions (Hellinger,
        Kolmogorov,\... distance). Finally package
        [distrMod](../packages/distrMod/index.html) provides functions
        for the computation of minimum criterion estimators (maximum
        likelihood and minimum distance estimators). See other packages
        of the distr-family
        ([distrSim](../packages/distrSim/index.html),
        [distrTEst](../packages/distrTEst/index.html),
        [distrTeach](../packages/distrTeach/index.html),
        [distrDoc](../packages/distrDoc/index.html),
        [distrEllipse](../packages/distrEllipse/index.html)).
    -   *Transformation* : Lebesgue decomposition are implemented in
        [distr](../packages/distr/index.html), as well as Convolution,
        Truncation and Huberization of distributions. Furthermore,
        [distr](../packages/distr/index.html) provides distribution of
        the maximum or minimum of two distributions. See
        Object-orientation above.
    -   *User Interface* : [AtelieR](../packages/AtelieR/index.html)
        package provides a GTK GUI for teaching basic concepts in
        statistical inference, implementing all the R base distributions
        as well as the generalized Student, the inverse Chi-square, the
        inverse gamma and the lambda-prime distributions.
-   *Transversal functions* :
    -   *Histogram, tail plots, distance estimation* :
        [DistributionUtils](../packages/DistributionUtils/index.html)
        provides log-histogram, tail plots, functions for testing
        distributions using inversion tests and the Massart inequality.
        [GMD](../packages/GMD/index.html) is a package for
        non-parametric distance measurement between two discrete
        frequency distributions.
    -   *Paremeter estimation* : [lmomco](../packages/lmomco/index.html)
        and [Lmoments](../packages/Lmoments/index.html) focus on
        univariate/multivariate (L-)moments estimation.
        [VGAM](../packages/VGAM/index.html) provides a lot of parameter
        estimation for usual and \"exotic\" distributions.
        [gaussDiff](../packages/gaussDiff/index.html) provides a
        collection difference measures for multivariate Gaussian
        probability density functions Package
        [MASS](../packages/MASS/index.html) implements the flexible
        `fitdistr` function for parameter estimations.
        [fitdistrplus](../packages/fitdistrplus/index.html) greatly
        enlarges and enhances the tools to fit any probability
        distribution. [EnvStats](../packages/EnvStats/index.html) and
        [fitteR](../packages/fitteR/index.html) also provides tools to
        fit most common distributions.
        [flexsurv](../packages/flexsurv/index.html) and
        [msm](../packages/msm/index.html) provides a quantile function
        for a generic distribution based on numerical computation based
        on a dichotomic search.

</div>

### CRAN packages:

-   [actuar](../packages/actuar/index.html) (core)
-   [AdMit](../packages/AdMit/index.html)
-   [agricolae](../packages/agricolae/index.html)
-   [ald](../packages/ald/index.html)
-   [AtelieR](../packages/AtelieR/index.html)
-   [bayesm](../packages/bayesm/index.html)
-   [benchden](../packages/benchden/index.html)
-   [BiasedUrn](../packages/BiasedUrn/index.html)
-   [BivarP](../packages/BivarP/index.html)
-   [bmixture](../packages/bmixture/index.html)
-   [bridgedist](../packages/bridgedist/index.html)
-   [CDVine](../packages/CDVine/index.html)
-   [cmvnorm](../packages/cmvnorm/index.html)
-   [coga](../packages/coga/index.html)
-   [CompGLM](../packages/CompGLM/index.html)
-   [CompLognormal](../packages/CompLognormal/index.html)
-   [compoisson](../packages/compoisson/index.html)
-   [Compositional](../packages/Compositional/index.html)
-   [Compounding](../packages/Compounding/index.html)
-   [CompQuadForm](../packages/CompQuadForm/index.html)
-   [condMVNorm](../packages/condMVNorm/index.html)
-   [copBasic](../packages/copBasic/index.html)
-   [copula](../packages/copula/index.html) (core)
-   [csn](../packages/csn/index.html)
-   [Davies](../packages/Davies/index.html)
-   [degreenet](../packages/degreenet/index.html)
-   [Delaporte](../packages/Delaporte/index.html)
-   [dirmult](../packages/dirmult/index.html)
-   [disclap](../packages/disclap/index.html)
-   [DiscreteInverseWeibull](../packages/DiscreteInverseWeibull/index.html)
-   [DiscreteLaplace](../packages/DiscreteLaplace/index.html)
-   [DiscreteWeibull](../packages/DiscreteWeibull/index.html)
-   [distcrete](../packages/distcrete/index.html)
-   [distr](../packages/distr/index.html) (core)
-   [distrDoc](../packages/distrDoc/index.html)
-   [distrEllipse](../packages/distrEllipse/index.html)
-   [distrEx](../packages/distrEx/index.html)
-   [DistributionUtils](../packages/DistributionUtils/index.html)
-   [distrMod](../packages/distrMod/index.html)
-   [distrSim](../packages/distrSim/index.html)
-   [distrTeach](../packages/distrTeach/index.html)
-   [distrTEst](../packages/distrTEst/index.html)
-   [dng](../packages/dng/index.html)
-   [e1071](../packages/e1071/index.html)
-   [ecd](../packages/ecd/index.html)
-   [emdbook](../packages/emdbook/index.html)
-   [emg](../packages/emg/index.html)
-   [EnvStats](../packages/EnvStats/index.html)
-   [evd](../packages/evd/index.html)
-   [evdbayes](../packages/evdbayes/index.html)
-   [evir](../packages/evir/index.html)
-   [extraDistr](../packages/extraDistr/index.html)
-   [extremefit](../packages/extremefit/index.html)
-   [FAdist](../packages/FAdist/index.html)
-   [FatTailsR](../packages/FatTailsR/index.html)
-   [fBasics](../packages/fBasics/index.html)
-   [fCopulae](../packages/fCopulae/index.html) (core)
-   [fExtremes](../packages/fExtremes/index.html)
-   [fgac](../packages/fgac/index.html)
-   [fitdistrplus](../packages/fitdistrplus/index.html)
-   [fitteR](../packages/fitteR/index.html)
-   [flexsurv](../packages/flexsurv/index.html)
-   [FMStable](../packages/FMStable/index.html)
-   [fpow](../packages/fpow/index.html)
-   [frmqa](../packages/frmqa/index.html)
-   [gambin](../packages/gambin/index.html)
-   [gamlss.dist](../packages/gamlss.dist/index.html) (core)
-   [gamlss.mx](../packages/gamlss.mx/index.html)
-   [gaussDiff](../packages/gaussDiff/index.html)
-   [gb](../packages/gb/index.html)
-   [GB2](../packages/GB2/index.html)
-   [GenBinomApps](../packages/GenBinomApps/index.html)
-   [GeneralizedHyperbolic](../packages/GeneralizedHyperbolic/index.html)
-   [GenOrd](../packages/GenOrd/index.html)
-   [geoR](../packages/geoR/index.html)
-   [ghyp](../packages/ghyp/index.html)
-   [GIGrvg](../packages/GIGrvg/index.html)
-   [gld](../packages/gld/index.html)
-   [GLDEX](../packages/GLDEX/index.html)
-   [glogis](../packages/glogis/index.html)
-   [GMD](../packages/GMD/index.html)
-   [GSM](../packages/GSM/index.html)
-   [gumbel](../packages/gumbel/index.html)
-   [HAC](../packages/HAC/index.html)
-   [hermite](../packages/hermite/index.html)
-   [HI](../packages/HI/index.html)
-   [HistogramTools](../packages/HistogramTools/index.html)
-   [hyper2](../packages/hyper2/index.html)
-   [HyperbolicDist](../packages/HyperbolicDist/index.html)
-   [ihs](../packages/ihs/index.html)
-   [kernelboot](../packages/kernelboot/index.html)
-   [kolmim](../packages/kolmim/index.html)
-   [KScorrect](../packages/KScorrect/index.html)
-   [LambertW](../packages/LambertW/index.html)
-   [LaplacesDemon](../packages/LaplacesDemon/index.html)
-   [LearnBayes](../packages/LearnBayes/index.html)
-   [lhs](../packages/lhs/index.html)
-   [LIHNPSD](../packages/LIHNPSD/index.html)
-   [lmom](../packages/lmom/index.html)
-   [lmomco](../packages/lmomco/index.html) (core)
-   [Lmoments](../packages/Lmoments/index.html)
-   [logitnorm](../packages/logitnorm/index.html)
-   [loglognorm](../packages/loglognorm/index.html)
-   [marg](../packages/marg/index.html)
-   [MASS](../packages/MASS/index.html)
-   [mbbefd](../packages/mbbefd/index.html)
-   [mc2d](../packages/mc2d/index.html)
-   [mclust](../packages/mclust/index.html)
-   [MCMCpack](../packages/MCMCpack/index.html)
-   [mded](../packages/mded/index.html)
-   [mgpd](../packages/mgpd/index.html)
-   [minimax](../packages/minimax/index.html)
-   [MitISEM](../packages/MitISEM/index.html)
-   [MittagLeffleR](../packages/MittagLeffleR/index.html)
-   [MixedTS](../packages/MixedTS/index.html)
-   [mixtools](../packages/mixtools/index.html)
-   [MM](../packages/MM/index.html)
-   [mnormpow](../packages/mnormpow/index.html)
-   [mnormt](../packages/mnormt/index.html) (core)
-   [modeest](../packages/modeest/index.html)
-   [moments](../packages/moments/index.html)
-   [movMF](../packages/movMF/index.html)
-   [msm](../packages/msm/index.html)
-   [mvprpb](../packages/mvprpb/index.html)
-   [mvrtn](../packages/mvrtn/index.html)
-   [mvtnorm](../packages/mvtnorm/index.html) (core)
-   [nCDunnett](../packages/nCDunnett/index.html)
-   [Newdistns](../packages/Newdistns/index.html)
-   [nor1mix](../packages/nor1mix/index.html)
-   [NormalGamma](../packages/NormalGamma/index.html)
-   [NormalLaplace](../packages/NormalLaplace/index.html)
-   [normalp](../packages/normalp/index.html)
-   [ORDER2PARENT](../packages/ORDER2PARENT/index.html)
-   [OrdNor](../packages/OrdNor/index.html)
-   [ParetoPosStable](../packages/ParetoPosStable/index.html)
-   [PDQutils](../packages/PDQutils/index.html)
-   [PearsonDS](../packages/PearsonDS/index.html) (core)
-   [PhaseType](../packages/PhaseType/index.html)
-   [poibin](../packages/poibin/index.html)
-   [poilog](../packages/poilog/index.html)
-   [poistweedie](../packages/poistweedie/index.html)
-   [polyaAeppli](../packages/polyaAeppli/index.html)
-   [poweRlaw](../packages/poweRlaw/index.html)
-   [qmap](../packages/qmap/index.html)
-   [QRM](../packages/QRM/index.html)
-   [randaes](../packages/randaes/index.html)
-   [random](../packages/random/index.html)
-   [randtoolbox](../packages/randtoolbox/index.html)
-   [RDieHarder](../packages/RDieHarder/index.html)
-   [ReIns](../packages/ReIns/index.html)
-   [reliaR](../packages/reliaR/index.html) (core)
-   [Renext](../packages/Renext/index.html)
-   [retimes](../packages/retimes/index.html)
-   [revdbayes](../packages/revdbayes/index.html)
-   [rlecuyer](../packages/rlecuyer/index.html)
-   [RMKdiscrete](../packages/RMKdiscrete/index.html)
-   [RMTstat](../packages/RMTstat/index.html)
-   [rngwell19937](../packages/rngwell19937/index.html)
-   [rpgm](../packages/rpgm/index.html)
-   [rstream](../packages/rstream/index.html)
-   [RTDE](../packages/RTDE/index.html)
-   [rtdists](../packages/rtdists/index.html)
-   [Runuran](../packages/Runuran/index.html)
-   [rust](../packages/rust/index.html)
-   [s20x](../packages/s20x/index.html)
-   [sadists](../packages/sadists/index.html)
-   [SCI](../packages/SCI/index.html)
-   [setRNG](../packages/setRNG/index.html)
-   [sfsmisc](../packages/sfsmisc/index.html)
-   [sgt](../packages/sgt/index.html)
-   [skellam](../packages/skellam/index.html)
-   [SkewHyperbolic](../packages/SkewHyperbolic/index.html)
-   [skewt](../packages/skewt/index.html)
-   [sld](../packages/sld/index.html)
-   [smoothmest](../packages/smoothmest/index.html)
-   [SMR](../packages/SMR/index.html)
-   [sn](../packages/sn/index.html)
-   [sparseMVN](../packages/sparseMVN/index.html)
-   [spd](../packages/spd/index.html)
-   [stabledist](../packages/stabledist/index.html)
-   [STAR](../packages/STAR/index.html)
-   [statmod](../packages/statmod/index.html)
-   [SuppDists](../packages/SuppDists/index.html)
-   [symmoments](../packages/symmoments/index.html)
-   [tmvtnorm](../packages/tmvtnorm/index.html)
-   [tolerance](../packages/tolerance/index.html)
-   [trapezoid](../packages/trapezoid/index.html)
-   [triangle](../packages/triangle/index.html)
-   [truncnorm](../packages/truncnorm/index.html)
-   [TSA](../packages/TSA/index.html)
-   [tsallisqexp](../packages/tsallisqexp/index.html)
-   [TTmoment](../packages/TTmoment/index.html)
-   [tweedie](../packages/tweedie/index.html)
-   [VarianceGamma](../packages/VarianceGamma/index.html)
-   [VGAM](../packages/VGAM/index.html) (core)
-   [VineCopula](../packages/VineCopula/index.html)
-   [vines](../packages/vines/index.html)
-   [zipfR](../packages/zipfR/index.html)

### Related links:

-   [Download statistics per view](http://rpackages.io/views)
-   [Advices to implement (new) distributions in
    R](http://www.rmetrics.org/Meielisalp2009/Presentations/Scott.pdf)
-   [Clickable diagram of distribution
    relationships](http://www.johndcook.com/distribution_chart.html)
-   [Diagram of discrete distribution
    relationships](http://www.stat.rice.edu/~dobelman/courses/texts/Distributions.Discrete.Kendall.jpg)
-   [Diagram of continuous distribution
    relationships](http://www.stat.rice.edu/~dobelman/courses/texts/Distributions.Chart.C&B.pdf)
-   [List and diagram of distribution
    relationship.](http://www.stat.rice.edu/~dobelman/courses/texts/leemis.distributions.2008amstat.pdf)
-   [Compendium of
    distributions.](http://www.stat.rice.edu/~dobelman/courses/DistributionCompendium.pdf)
-   [A comprehensive list of data
    types](https://en.wikipedia.org/wiki/Statistical_data_type)
-   [Journal of Statistical Software: R programs for truncated
    distributions](http://www.jstatsoft.org/v16/c02/)
-   [The \"distrXXX\"-family of
    R-packages](http://distr.r-forge.r-project.org/)
-   [Overview of vine copula
    models](http://www-m4.ma.tum.de/forschung/vine-copula-models/#c662)
