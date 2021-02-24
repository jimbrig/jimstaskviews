
# CRAN Task Views

-   [Official CRAN Task Views](https://cran.r-project.org/web/views/)
-   [Shiny
    Application](https://jimbrigapps.shinyapps.io/task-views-app/)

## Task Views

⭐ = Favorite / Recommendation 

-   [Package Development](package-development/) ⭐⭐
-   [Webservices](webservices/) ⭐⭐
-   [Finance](finance/)
-   [Mapping Tools](maptools/) ⭐
-   [High Performance Computing](hpc/)
-   [Taxonomy](taxonomy/)
-   [Databases](databases/) ⭐
-   [Model Deployment](model-deployment/)
-   [Reproducible Research](reproducible-research/) ⭐
-   [Anamaly Detection](anamoly-detection/)
-   [Security](security/)
-   [Time-Series](timeseries/)
-   [Optimization](optimization/) ⭐
-   [Graphics](graphics/)
-   [Numerical Math](numerical-math/)
-   [Open Data](opendata/) ⭐
-   [Hydrology](hydrology/)
-   [Computational environments](computational-environments/) ⭐⭐

Additionally, I included the `ctv` R package.

## Shiny App

There is also an exploratory [Shiny App](shiny_app/) that van be viewed
at <https://jimbrigapps.shinyapps.io/task-views-app/>

Quickly browse packages and their licenses using a snapshot of task
views retrieved on 2020-02-28.

*Not all packages available on CRAN are shown here, just the ones in
CRAN’s [Task Views](https://cran.r-project.org/web/views/).*

## Techinical Details and Code

Repo’s are added via [git
submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules); see
the `.gitmodules` file for details.

``` bash
git init

git submodule add git@github.com:ropensci/PackageDevelopment.git package-development
git submodule add git@github.com:ropensci/webservices.git
git submodule add git@github.com:eddelbuettel/ctv-finance.git finance
git submodule add git@github.com:ropensci/maptools.git
git submodule add git@github.com:eddelbuettel/ctv-hpc.git "high-performace-computing"
git submodule add git@github.com:ropensci/taxonomy.git
git submodule add git@github.com:terrytangyuan/ctv-databases.git "databases"
git submodule add git@github.com:terrytangyuan/ctv-model-deployment.git "model-deployments"
git submodule add git@github.com:jdblischak/reproducible-research-ctv.git "reproducible-research"
git submodule add git@github.com:pridiltal/ctv-AnomalyDetection.git "anamoly-detection"
git submodule add git@github.com:hrbrmstr/ctvsecurity.git "security"
git submodule add git@github.com:robjhyndman/ctv-TimeSeries.git "timeseries"
git submodule add git@github.com:hwborchers/ctv-optimization.git "optimization"
git submodule add git@github.com:sctyner/ctv-graphics.git "graphics"
git submodule add git@github.com:cran/ctv.git
git submodule add git@github.com:hwborchers/ctv-numericalmath.git "numerical-math"
git submodule add git@github.com:ropensci/opendata.git opendata
git submodule add git@github.com:ropensci/Hydrology.git hydrology
git submodule add git@github.com:o2r-project/ctv-computational-environments.git computational-environments
git submodule add git@github.com:r-spatial/task_views.git "spacial-and-temporal-spatial"

exit
```

In order to stay up-to-date utilize the `update_submods.sh` shell
script:

``` bash
#!/bin/sh

git submodule update --recursive
```

## Remaining CTVs not found on Github:

The remaining task views are taken from CRAN and placed into the
[data](data/) folder.

-   NOTE: I use pandoc to convert between HTML and Markdown here.

``` r
remaining <- c(
  "Distributions",
  "Bayesian",
  "Econometrics",
  "Environmetrics",
  "ExperimentalDesign",
  "ExtremeValue",
  "FunctionalData",
  "MachineLearning",
  "MetaAnalysis",
  "Multivariate",
  "NaturalLanguageProcessing",
  "OfficialStatistics",
  "Robust",
  "SocialSciences",
  "Spatial",
  "SpatioTemporal",
  "Survival",
  "WebTechnologies"
)

ctvs <- purrr::map(remaining, function(x) {
  read.ctv(system.file("ctv", paste0(x, ".ctv"), package = "ctv"))
}) %>% setNames(remaining)

purrr::walk2(
  ctvs, 
  names(ctvs), 
  ~ctv::ctv2html(x = .x, file = fs::path("data/html", paste0(.y, ".html")))
)

fs::dir_create("data/md")

purrr::walk(basename(fs::dir_ls("data/html")), function(x) {
  cmd <- paste0("pandoc -o 'data/md/", fs::path_ext_remove(x), ".md' ", x)
  cmd
  system(cmd)
})
```

## Session Info

``` r
sessioninfo::session_info()
```

    - Session info ---------------------------------------------------------------
     setting  value                       
     version  R version 4.0.3 (2020-10-10)
     os       Windows 10 x64              
     system   x86_64, mingw32             
     ui       RTerm                       
     language en-US                       
     collate  English_United States.1252  
     ctype    English_United States.1252  
     tz       America/New_York            
     date     2021-02-23                  

    - Packages -------------------------------------------------------------------
     package     * version date       lib source        
     assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.1.0)
     cli           2.3.1   2021-02-23 [1] CRAN (R 4.0.3)
     digest        0.6.27  2020-10-24 [1] CRAN (R 4.0.3)
     evaluate      0.14    2019-05-28 [1] CRAN (R 4.1.0)
     fs            1.5.0   2020-07-31 [1] CRAN (R 4.0.3)
     glue          1.4.2   2020-08-27 [1] CRAN (R 4.1.0)
     htmltools     0.5.1.1 2021-01-22 [1] CRAN (R 4.1.0)
     knitr         1.31    2021-01-27 [1] CRAN (R 4.1.0)
     magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.0.3)
     rlang         0.4.10  2020-12-30 [1] CRAN (R 4.0.3)
     rmarkdown     2.6     2020-12-14 [1] CRAN (R 4.0.3)
     sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 4.0.3)
     stringi       1.5.3   2020-09-09 [1] CRAN (R 4.0.3)
     stringr       1.4.0   2019-02-10 [1] CRAN (R 4.0.3)
     withr         2.4.1   2021-01-26 [1] CRAN (R 4.1.0)
     xfun          0.21    2021-02-10 [1] CRAN (R 4.1.0)
     yaml          2.2.1   2020-02-01 [1] CRAN (R 4.0.3)

    [1] C:/Users/Jimmy Briggs/.R/lib/4.0
    [2] C:/Program Files/R/R-4.0.3/library

## ALL Task Views Detailed Information

To get all CTVs:

``` r
library(ctv)

cdir <- system.file("ctv", package = "ctv")
ctvs <- list.files(cdir, pattern = "\\.ctv$")

rr <- sapply(ctvs,
             function(ctv) {
               cat(sprintf("%25s  ", ctv))
               R <- read.ctv(file.path(cdir, ctv))
               cat("\n")
               R
             },
             simplify = FALSE)
```

                 Bayesian.ctv  
                 ChemPhys.ctv  
           ClinicalTrials.ctv  
                  Cluster.ctv  
    DifferentialEquations.ctv  
            Distributions.ctv  
             Econometrics.ctv  
           Environmetrics.ctv  
       ExperimentalDesign.ctv  
             ExtremeValue.ctv  
                  Finance.ctv  
           FunctionalData.ctv  
                 Genetics.ctv  
                       gR.ctv  
                 Graphics.ctv  
    HighPerformanceComputing.ctv  
          MachineLearning.ctv  
           MedicalImaging.ctv  
             MetaAnalysis.ctv  
          ModelDeployment.ctv  
             Multivariate.ctv  
    NaturalLanguageProcessing.ctv  
     NumericalMathematics.ctv  
       OfficialStatistics.ctv  
             Optimization.ctv  
         Pharmacokinetics.ctv  
            Phylogenetics.ctv  
            Psychometrics.ctv  
     ReproducibleResearch.ctv  
                   Robust.ctv  
           SocialSciences.ctv  
                  Spatial.ctv  
           SpatioTemporal.ctv  
                 Survival.ctv  
               TimeSeries.ctv  
          WebTechnologies.ctv  

and for more details, including the packages within the CTVs:

``` r
for (n in names(rr)) {
  cat(n, " :\n", rep.int("=", nchar(n)), "\n", sep = '')
  knitr::knit_print(rr[[n]])
  cat("--------------------------------------------------------\n")
}
```

    Bayesian.ctv :
    ============

    CRAN Task View
    --------------
    Name:       Bayesian
    Topic:      Bayesian Inference
    Maintainer: Jong Hee Park
    Contact:    jongheepark@snu.ac.kr
    Version:    2018-04-13

    Packages:   abc, abn, AdMit, arm*, AtelieR, BaBooN, BACCO*, BaM, bamlss, BAS,
                BayesDA, BayesFactor, bayesGARCH, bayesImageS, bayesm*, bayesmeta,
                bayesmix, bayesQR, BayesSummaryStatLM, bayesSurv*, Bayesthresh,
                BayesTree, BayesValidate, BayesVarSel, BayesX, BayHaz, BAYSTAR,
                bbemkr, BCBCSF, BCE, bclust, bcp, bisoreg, BLR, BMA, Bmix, BMS,
                bnlearn, boa*, Bolstad, Boom, BoomSpikeSlab, bqtl, bridgesampling,
                brms, bsamGP, bspec, bspmma, BSquare, bsts, BVS, catnet,
                coalescentMCMC, coda*, cudaBayesreg, dclone, deal, deBInfer, dlm,
                DPpackage*, EbayesThresh, ebdbNet, eco, eigenmodel, ensembleBMA,
                evdbayes, exactLoglinTest, factorQR, FME, geoR, geoRglm, ggmcmc,
                glmmBUGS, gRain, growcurves, hbsae, HI, Hmisc, iterLap,
                LaplacesDemon, LearnBayes, lme4, lmm, MasterBayes, matchingMarkets,
                mcmc*, MCMCglmm, MCMCpack*, mgcv, mlogitBMA, MNP, mombf, monomvn,
                MSBVAR, NetworkChange, nimble*, openEBGM, pacbpred, PAWL,
                predmixcor, PReMiuM, prevalence, profdpm, pscl, R2BayesX, R2jags,
                R2WinBUGS, ramps, rbugs, revdbayes, RJaCGH, rjags, RSGHB, RSGHB,
                rstan, rstiefel, runjags, Runuran, RxCEcolInf, SamplerCompare,
                SampleSizeMeans, SampleSizeProportions, sbgcop, SimpleTable, sna,
                spBayes, spikeslab, spikeSlabGAM, spTimer, stochvol, tgp,
                tRophicPosition, zic
                (* = core package)

    --------------------------------------------------------
    ChemPhys.ctv :
    ============

    CRAN Task View
    --------------
    Name:       ChemPhys
    Topic:      Chemometrics and Computational Physics
    Maintainer: Katharine Mullen
    Contact:    katharine.mullen@stat.ucla.edu
    Version:    2018-01-24

    Packages:   ALS*, AnalyzeFMRI, AquaEnv, astro, astrochron, astrodatR, astroFns,
                astrolibR, Bchron, BioMark, bvls, celestial, CellularAutomaton,
                chemCal*, chemometrics, ChemometricsWithR, ChemoSpec, CHNOSZ,
                clustvarsel, compositions, cosmoFns, CosmoPhotoz, CRAC, dielectric,
                diffractometry, drc, drm, EEM, elasticnet, enpls, fastICA, FITSio,
                fmri, fpca, FTICRMS, homals, hyperSpec, investr, Iso*, kohonen*,
                leaps, lira, lspls, MALDIquant, minpack.lm, moonsun, nlme, nlreg,
                nnls*, OrgMassSpecR, pcaPP, Peaks, PET, planar, pls*, plspm, ppls,
                prospectr, psy, PTAk*, quantchem, rcdk, rcdklibs, represent,
                resemble, RobPer, rpubchem, sapa, SCEPtER, SCEPtERbinary, simecol,
                snapshot, solaR, som, SPADAR, speaq, spls, stellaR, stepPlr,
                subselect, TIMP, titan, titrationCurves, UPMASK, varSelRF, webchem,
                WilcoxCV
                (* = core package)

    --------------------------------------------------------
    ClinicalTrials.ctv :
    ==================

    CRAN Task View
    --------------
    Name:       ClinicalTrials
    Topic:      Clinical Trial Design, Monitoring, and Analysis
    Maintainer: Ed Zhang and Harry G. Zhang
    Contact:    Ed.Zhang.jr@gmail.com
    Version:    2018-03-26

    Packages:   adaptTest*, AGSDest, asd*, asypow, bcrm*, binomSamSize, blockrand*,
                clinfun*, clinsig, clusterPower, coin, conf.design, CRM, CRTSize*,
                dfcrm*, dfped, dfpk, DoseFinding, epibasix, ewoc, experiment*,
                FrF2, GroupSeq*, gsbDesign, gsDesign*, HH, Hmisc*,
                InformativeCensoring, ldbounds*, longpower, MChtest*, MCPMod,
                Mediana, meta, metafor, metaLik, metasens, multcomp, nppbib, PIPS*,
                PowerTOST*, pwr*, PwrGSD*, qtlDesign*, rmeta, samplesize, seqmon*,
                speff2trial*, ssanv, survival*, TEQR*, ThreeArmedTrials,
                ThreeGroups, TrialSize*
                (* = core package)

    --------------------------------------------------------
    Cluster.ctv :
    ===========

    CRAN Task View
    --------------
    Name:       Cluster
    Topic:      Cluster Analysis & Finite Mixture Models
    Maintainer: Friedrich Leisch and Bettina Gruen
    Contact:    Bettina.Gruen@jku.at
    Version:    2018-05-04

    Packages:   AdMit, ADPclust, amap, apcluster, BayesLCA, bayesm, bayesmix,
                bclust, bgmm, biclust, Bmix, bmixture, cba, cclust, CEC, CHsharp,
                clue, cluster*, clusterfly, clusterGeneration, clusterRepro,
                clusterSim, clustMixType, clustvarsel, clv, clValid, CoClust,
                compHclust, dbscan, dendextend, depmix, depmixS4, dpmixsim,
                dynamicTreeCut, e1071, edci, EMCluster, evclust, FactoClass,
                fastcluster, fclust, flashClust, flexclust*, flexCWM, flexmix*,
                fpc, FunCluster, funFEM, funHDDC, gamlss.mx, genie, GLDEX, GSM,
                HDclassif, hybridHclust, idendr0, IMIFA, isopam, kernlab, kml,
                largeVis, latentnet, LCAvarsel, lcmm, longclust, mcclust, mclust*,
                MetabolAnalyze, mixAK, mixdist, mixPHM, mixRasch, mixreg, MixSim,
                mixsmsn, mixtools, mixture, MOCCA, MoEClust, movMF, mritc, NbClust,
                nor1mix, optpart, ORIClust, pdfCluster, pmclust, poLCA, prabclus,
                prcr, PReMiuM, profdpm, protoclust, psychomix, pvclust, randomLCA,
                rebmix, rjags, Rmixmod*, RPMM, seriation, sigclust, skmeans, som,
                sparcl, tclust, teigen, treeClust, trimcluster, VarSelLCM, wle
                (* = core package)

    --------------------------------------------------------
    DifferentialEquations.ctv :
    =========================

    CRAN Task View
    --------------
    Name:       DifferentialEquations
    Topic:      Differential Equations
    Maintainer: Karline Soetaert and Thomas Petzoldt
    Contact:    karline.soetaert@nioz.nl
    Version:    2017-11-22

    Packages:   adaptivetau, bvpSolve*, cOde, CollocInfer, deSolve*, deTestSet,
                dMod, ecolMod, FME, GillespieSSA, mkin, nlmeODE, odeintr,
                PBSddesolve, PBSmodelling, phaseR, pomp, pracma, primer, QPot,
                ReacTran, rODE, rodeo, rootSolve*, rpgm, scaRabee, sde*,
                Sim.DiffProc, simecol
                (* = core package)

    --------------------------------------------------------
    Distributions.ctv :
    =================

    CRAN Task View
    --------------
    Name:       Distributions
    Topic:      Probability Distributions
    Maintainer: Christophe Dutang, Patrice Kiener
    Contact:    Christophe.Dutang@ensimag.fr
    Version:    2018-04-17

    Packages:   actuar*, AdMit, agricolae, ald, AtelieR, bayesm, benchden,
                BiasedUrn, BivarP, bmixture, bridgedist, CDVine, cmvnorm, coga,
                CompGLM, CompLognormal, compoisson, Compositional, Compounding,
                CompQuadForm, condMVNorm, copBasic, copula*, csn, Davies,
                degreenet, Delaporte, dirmult, disclap, DiscreteInverseWeibull,
                DiscreteLaplace, DiscreteWeibull, distcrete, distr*, distrDoc,
                distrEllipse, distrEx, DistributionUtils, distrMod, distrSim,
                distrTeach, distrTEst, dng, e1071, ecd, emdbook, emg, EnvStats,
                evd, evdbayes, evir, extraDistr, extremefit, FAdist, FatTailsR,
                fBasics, fCopulae*, fExtremes, fgac, fitdistrplus, fitteR,
                flexsurv, FMStable, fpow, frmqa, gambin, gamlss.dist*, gamlss.mx,
                gaussDiff, gb, GB2, GenBinomApps, GeneralizedHyperbolic, GenOrd,
                geoR, ghyp, GIGrvg, gld, GLDEX, glogis, GMD, GSM, gumbel, HAC,
                hermite, HI, HistogramTools, hyper2, HyperbolicDist, ihs,
                kernelboot, kolmim, KScorrect, LambertW, LaplacesDemon, LearnBayes,
                lhs, LIHNPSD, lmom, lmomco*, Lmoments, logitnorm, loglognorm, marg,
                MASS, mbbefd, mc2d, mclust, MCMCpack, mded, mgpd, minimax, MitISEM,
                MittagLeffleR, MixedTS, mixtools, MM, mnormpow, mnormt*, modeest,
                moments, movMF, msm, mvprpb, mvrtn, mvtnorm*, nCDunnett, Newdistns,
                nor1mix, NormalGamma, NormalLaplace, normalp, ORDER2PARENT, OrdNor,
                ParetoPosStable, PDQutils, PearsonDS*, PhaseType, poibin, poilog,
                poistweedie, polyaAeppli, poweRlaw, qmap, QRM, randaes, random,
                randtoolbox, RDieHarder, ReIns, reliaR*, Renext, retimes,
                revdbayes, rlecuyer, RMKdiscrete, RMTstat, rngwell19937, rpgm,
                rstream, RTDE, rtdists, Runuran, rust, s20x, sadists, SCI, setRNG,
                sfsmisc, sgt, skellam, SkewHyperbolic, skewt, sld, smoothmest, SMR,
                sn, sparseMVN, spd, stabledist, STAR, statmod, SuppDists,
                symmoments, tmvtnorm, tolerance, trapezoid, triangle, truncnorm,
                TSA, tsallisqexp, TTmoment, tweedie, VarianceGamma, VGAM*,
                VineCopula, vines, zipfR
                (* = core package)

    --------------------------------------------------------
    Econometrics.ctv :
    ================

    CRAN Task View
    --------------
    Name:       Econometrics
    Topic:      Econometrics
    Maintainer: Achim Zeileis
    Contact:    Achim.Zeileis@R-project.org
    Version:    2018-04-24

    Packages:   AER*, aod, apt, bayesm, betareg, BMA, BMS, boot, bootstrap, brglm,
                CADFtest, car*, CDNmoney, censReg, clubSandwich, clusterSEs, crch,
                decompr, dlsem, dynlm, Ecdat, effects, erer, expsmooth,
                ExtremeBounds, fma, forecast*, frm, frontier, fxregime, gam,
                gamlss, geepack, gets, glmx, gmm, gmnl, gvc, Hmisc, ineq, intReg,
                ivbma, ivfixed, ivlewbel, ivpack, ivpanel, ivprobit, LARF, lavaan,
                lfe, LinRegInteractive, lme4, lmtest*, margins, MASS,
                matchingMarkets, Matrix, Mcomp, meboot, mfx, mgcv, mhurdle,
                micEcon, micEconAids, micEconCES, micEconSNQP, midasr, mlogit,
                mnlogit, MNP, MSBVAR, multiwayvcov, mvProbit, nlme, nnet, nonnest2,
                np, ordinal, OrthoPanels, pampe, panelAR, Paneldata, panelvar,
                PANICr, pco, pcse, pder, pdR, pglm, phtt, plm*, pscl, psidR, pwt,
                pwt8, pwt9, quantreg, Rchoice, rdd, rddapp, rddtools, rdlocrand,
                rdmulti, rdpower, rdrobust, reldist, REndo, rms, RSGHB,
                rUnemploymentData, sampleSelection, sandwich*, segmented, sem,
                SemiParSampleSel, semsfa, sfa, simpleboot, SparseM, spatialprobit,
                spdep, spfrontier, sphet, splm, ssfa, strucchange, survival,
                systemfit, truncreg, tsDyn, tseries*, tsfa, urca*, vars, VGAM,
                wahc, wbstats, wooldridge, xts, Zelig, zoo*, zTree
                (* = core package)

    --------------------------------------------------------
    Environmetrics.ctv :
    ==================

    CRAN Task View
    --------------
    Name:       Environmetrics
    Topic:      Analysis of Ecological and Environmental Data
    Maintainer: Gavin Simpson
    Contact:    ucfagls@gmail.com
    Version:    2018-04-11

    Packages:   ade4*, amap, analogue, aod, ape, aqp, BiodiversityR, boussinesq,
                bReeze, CircStats, circular, cluster*, cocorresp, Distance,
                diveMove, dse, dsm, DSpat, dyn, dynatopmodel, dynlm, e1071, earth,
                eco, ecodist, EcoHydRology, EnvStats, equivalence, evd, evdbayes,
                evir, extRemes, fast, FD, flexmix, forecast, fso, gam, gamair,
                hydroGOF, HydroMe, hydroPSO, hydroTSM, Interpol.T, ipred, ismev,
                labdsv*, latticeDensity, lme4, maptree, marked, MASS*, mclust, mda,
                mefa, metacom, mgcv*, mrds, nlme, nsRFA, oce, openair, ouch, party,
                pastecs, pgirmess, popbio, prabclus, primer, pscl, pvclust, qualV,
                quantreg, quantregGrowth, randomForest, Rcapture, rioja, RMark,
                RMAWGEN, rpart, rtop, seacarb, seas, secr, segmented, sensitivity,
                simba, simecol, siplab, soiltexture, SPACECAP, SpatialExtremes,
                StreamMetabolism, strucchange, surveillance, tiger, topmodel,
                tseries, unmarked, untb, vegan*, vegetarian, VGAM, wasim, zoo
                (* = core package)

    --------------------------------------------------------
    ExperimentalDesign.ctv :
    ======================

    CRAN Task View
    --------------
    Name:       ExperimentalDesign
    Topic:      Design of Experiments (DoE) & Analysis of Experimental Data
    Maintainer: Ulrike Groemping
    Contact:    groemping@bht-berlin.de
    Version:    2018-02-17

    Packages:   acebayes, agricolae*, agridat, AlgDesign*, ALTopt, asd,
                BatchExperiments, BayesMAMS, bcrm, BHH2, binseqtest, bioOED,
                blocksdesign, blockTools, BOIN, BsMD, choiceDes, CombinS,
                conf.design*, crossdes*, Crossover, dae, daewr, designGG,
                designGLMM, designmatch, desirability, desplot, dfcomb, dfcrm,
                dfmta, dfpk, DiceDesign, DiceEval, DiceKriging, docopulae,
                DoE.base*, DoE.MIParray, DoE.wrapper*, DoseFinding, dynaTree,
                easypower, edesign, EngrExpt, experiment, ez, FMC, FrF2*,
                FrF2.catlg128, GAD, geospt, granova, GroupSeq, gsbDesign, gsDesign,
                gset, hiPOD, ibd, ICAOD, idefix, JMdesign, LDOD, lhs, MAMS, MaxPro,
                MBHdesign, minimalRSD, minimaxdesign, mixexp, mkssd, mxkssd, OBsMD,
                odr, OPDOE, optbdmaeAT, optDesignSlopeInt, OptGS, OptimalDesign,
                OptimaRegion, OptInterim, optrcdmaeAT, osDesign, PBIBD, PGM2,
                ph2bayes, ph2bye, pid, pipe.design, plgp, PopED, powerAnalysis,
                powerbydesign, powerGWASinteraction, PwrGSD, qtlDesign,
                qualityTools, RcmdrPlugin.DoE, rodd, RPPairwiseDesign, rsm*,
                rsurface, SensoMineR, seqDesign, sFFLHD, simrel, skpr*, SLHD,
                soptdmaeA, sp23design, ssize.fdr, ssizeRNA, support.CEs, TEQR, tgp,
                ThreeArmedTrials, toxtestD, unrepx, vdg, Vdgraph, VdgRsm, VNM
                (* = core package)

    --------------------------------------------------------
    ExtremeValue.ctv :
    ================

    CRAN Task View
    --------------
    Name:       ExtremeValue
    Topic:      Extreme Value Analysis
    Maintainer: Christophe Dutang, Kevin Jaunatre
    Contact:    Christophe.Dutang@ensimag.fr
    Version:    2017-12-26

    Packages:   copula, evd*, evdbayes, evir*, evmix, extremefit, extRemes,
                extremeStat, fExtremes, ismev, lmom, lmomco, lmomRFA, mev, POT,
                QRM, ReIns, Renext, revdbayes, RTDE, SpatialExtremes, texmex,
                threshr, VGAM
                (* = core package)

    --------------------------------------------------------
    Finance.ctv :
    ===========

    CRAN Task View
    --------------
    Name:       Finance
    Topic:      Empirical Finance
    Maintainer: Dirk Eddelbuettel
    Contact:    Dirk.Eddelbuettel@R-project.org
    Version:    2018-05-02

    Packages:   actuar, AmericanCallOpt, backtest, bayesGARCH, BCC1997,
                BenfordTests, betategarch, bizdays, BLModel, BurStFin, BurStMisc,
                CADFtest, car, ccgarch, ChainLadder, copula, CreditMetrics,
                credule, crp.CSFP, cvar, data.table, derivmkts, dlm, Dowd, dse,
                DtD, dyn, dynlm, ESG, estudy2, factorstochvol, fame, fAssets*,
                FatTailsR, fBasics*, fBonds*, fCopulae*, fExoticOptions*,
                fExtremes*, fgac, fGarch*, fImport*, financial, FinancialMath,
                FinAsym, finreportr, fmdates, fMultivar*, fNonlinear*, fOptions*,
                forecast, fPortfolio*, fracdiff, fractal, FRAPO, fRegression*,
                frmqa, fTrading*, GCPM, GetHFData, gets, GetTDData, GEVStableGarch,
                ghyp, gmm, gogarch, GUIDE, highfrequency, IBrokers, InfoTrad,
                lgarch, lifecontingencies, lmtest, longmemo, LSMonteCarlo,
                maRketSim, markovchain, MarkowitzR, matchingMarkets, MSBVAR,
                MSGARCH, mvtnorm, NetworkRiskMeasures, nlme, NMOF, obAnalytics,
                OptHedging, OptionPricing, pa, parma, pbo, PerformanceAnalytics*,
                pinbasic, portfolio, PortfolioEffectHFT, PortfolioOptim,
                portfolioSim, PortRisk, quantmod, QuantTools, ragtop, Rbitcoin,
                Rblpapi, Rcmdr, RcppQuantuccia, reinsureR, restimizeapi, Risk,
                RiskPortfolios, riskSimul, RM2006, rmgarch, RND, rpatrec, rpgm,
                RQuantLib, rugarch*, rwt, sandwich, sde, SharpeR, sharpeRratio,
                Sim.DiffProc, SmithWilsonYieldCurve, stochvol, strucchange,
                TAQMNGR, tawny, termstrc, TFX, tidyquant, timeDate*, timeSeries*,
                timsac, tis, TSdbi, tsDyn, tseries*, tseriesChaos, tsfa, TTR, tvm,
                urca*, vars, VarSwapPrice, vrtest, wavelets, waveslim, wavethresh,
                XBRL, xts*, ycinterextra, YieldCurve, Zelig, zoo*
                (* = core package)

    --------------------------------------------------------
    FunctionalData.ctv :
    ==================

    CRAN Task View
    --------------
    Name:       FunctionalData
    Topic:      Functional Data Analysis
    Maintainer: Fabian Scheipl
    Contact:    fabian.scheipl@stat.uni-muenchen.de
    Version:    2018-02-12

    Packages:   classiFunc, covsep, dbstats, ddalpha, denseFLMM, fda*, fda.usc*,
                fdadensity, fdakma, fdapace*, fdaPDE, fdasrvf*, fdatest, FDboost*,
                fdcov, fds, flars, fpca, freqdom, freqdom.fda, ftsa*, ftsspec,
                Funclustering, funcy*, funData, funFEM, funHDDC, geofd, GPFDA,
                growfunctions, pcdpca, rainbow, refund*, refund.shiny, refund.wave,
                RFgroove, roahd, sparseFLMM, switchnpreg, warpMix
                (* = core package)

    --------------------------------------------------------
    Genetics.ctv :
    ============

    CRAN Task View
    --------------
    Name:       Genetics
    Topic:      Statistical Genetics
    Maintainer: Giovanni Montana
    Contact:    g.montana@imperial.ac.uk
    Version:    2017-04-25

    Packages:   adegenet, ape, Biodem, bqtl, dlmap, gap*, GenABEL, genetics*,
                hapassoc, haplo.ccs, haplo.stats*, HardyWeinberg, hierfstat, hwde,
                ibdreg, LDheatmap, luca, ouch, pbatR, phangorn, qtl, rmetasim,
                seqinr, snp.plotter, SNPmaxsel, stepwise, tdthap, untb, wgaim
                (* = core package)

    --------------------------------------------------------
    gR.ctv :
    ======

    CRAN Task View
    --------------
    Name:       gR
    Topic:      gRaphical Models in R
    Maintainer: Soren Hojsgaard
    Contact:    sorenh@math.aau.dk
    Version:    2018-01-06

    Packages:   abn, bayesmix, BDgraph, bnlearn, bnstruct, boa, BRugs, catnet,
                coda, dclone, deal, diagram, DiagrammeR, dynamicGraph, ergm,
                FBFsearch, GeneNet, ggm, giRaph, gRain, gRbase*, gRc, gRim, huge,
                igraph, lvnet, mathgraph, MXM, ndtv, network, networkDynamic,
                parcor, pcalg, QUIC, R2OpenBUGS, R2WinBUGS, rbugs, rjags, SIN
                (* = core package)

    --------------------------------------------------------
    Graphics.ctv :
    ============

    CRAN Task View
    --------------
    Name:       Graphics
    Topic:      Graphic Displays & Dynamic Graphics & Graphic Devices & Visualization
    Maintainer: Nicholas Lewin-Koh
    Contact:    nikko@hailmail.net
    Version:    2015-01-07

    Packages:   ade4, animation, ape, aplpack, ash, biclust, Cairo, cairoDevice,
                cba, colorspace, diagram, dichromat, gclus, ggplot2*, gplots,
                gridBase, hexbin, IDPmisc, igraph, iplots, JavaGD, klaR, lattice*,
                latticeExtra, misc3d, onion, playwith, plotrix*, RColorBrewer*,
                rggobi, rgl*, RGraphics, RGtk2, RSvgDevice, RSVGTipsDevice,
                scagnostics, scatterplot3d, seriation, tkrplot, vcd*, vioplot,
                xgobi
                (* = core package)

    --------------------------------------------------------
    HighPerformanceComputing.ctv :
    ============================

    CRAN Task View
    --------------
    Name:       HighPerformanceComputing
    Topic:      High-Performance and Parallel Computing with R
    Maintainer: Dirk Eddelbuettel
    Contact:    Dirk.Eddelbuettel@R-project.org
    Version:    2018-05-07

    Packages:   aprof, batch, BatchExperiments, BatchJobs, batchtools, bayesm, bcp,
                BDgraph, biglars, biglm, bigmemory, bnlearn, caret, clustermq,
                cudaBayesreg, data.table, dclone, doFuture, doMC, doMPI, doRedis,
                doRNG, doSNOW, drake, ff, ffbase, flowr, foreach, future,
                future.BatchJobs, GAMBoost, gcbd, gpuR, GUIProfiler, h2o,
                HadoopStreaming, harvestr, HistogramTools, inline, keras, LaF,
                latentnet, lga, Matching, MonetDB.R, nws, OpenCL, orloca, parSim,
                partDSA, pbapply, pbdBASE, pbdDEMO, pbdDMAT, pbdMPI, pbdNCDF4,
                pbdPROF, pbdSLAP, peperr, permGPU, pls, pmclust, profr, proftools,
                pvclust, randomForestSRC, Rborist, Rcpp, RcppParallel, Rdsm,
                reticulate, rgenoud, Rhpc, RhpcBLASctl, RInside, rJava, rlecuyer,
                Rmpi*, RProtoBuf, rredis, rslurm, Sim.DiffProc, sitmo, snow*,
                snowfall, snowFT, speedglm, sqldf, ssgraph, STAR, tensorflow,
                tfestimators, tm, toaster, varSelRF, xgboost
                (* = core package)

    --------------------------------------------------------
    MachineLearning.ctv :
    ===================

    CRAN Task View
    --------------
    Name:       MachineLearning
    Topic:      Machine Learning & Statistical Learning
    Maintainer: Torsten Hothorn
    Contact:    Torsten.Hothorn@R-project.org
    Version:    2018-04-27

    Packages:   ahaz, arules, BART, bartMachine, BayesTree, biglasso, bigRR, bmrm,
                Boruta, bst, C50, caret, CORElearn, CoxBoost, Cubist, deepnet,
                e1071*, earth, effects, elasticnet, ElemStatLearn, evclass, evtree,
                FCNN4R, frbs, GAMBoost, gamboostLSS, gbm*, ggRandomForests, glmnet,
                glmpath, GMMBoost, gradDescent, grf, grplasso, grpreg, h2o, hda,
                hdi, hdm, ICEbox, ipred, kernlab*, klaR, lars, lasso2, LiblineaR,
                LogicForest, LogicReg, LTRCtrees, maptree, mboost*, mlr, model4you,
                MXM, ncvreg, nnet*, oem, OneR, opusminer, pamr, party, partykit,
                pdp, penalized, penalizedLDA, penalizedSVM, plotmo, quantregForest,
                randomForest*, randomForestSRC, ranger, rattle, Rborist, RcppDL,
                rda, rdetools, REEMtree, relaxo, rgenoud, RLT, Rmalschains, rminer,
                rnn, ROCR, RoughSets, rpart*, RPMM, RSNNS, RWeka, RXshrink, sda,
                SIS, spa, stabs, SuperLearner, svmpath, tensorflow, tgp, tree,
                trtf, varSelRF, vcrpart, wsrf, xgboost
                (* = core package)

    --------------------------------------------------------
    MedicalImaging.ctv :
    ==================

    CRAN Task View
    --------------
    Name:       MedicalImaging
    Topic:      Medical Image Analysis
    Maintainer: Brandon Whitcher
    Contact:    bwhitcher@gmail.com
    Version:    2018-01-24

    Packages:   adaptsmoFMRI, adimpro*, AnalyzeFMRI*, arf3DS4*, bayesImageS,
                bayesm, brainR, brainwaver, cudaBayesreg, DATforDCEMRI*, dcemriS4*,
                divest*, dpmixsim*, dti*, edfReader*, eegkit*, fmri*, fslr,
                gdimap*, mmand*, Morpho*, mritc*, neuroim*, neuRosim*, occ*,
                oro.dicom*, oro.nifti*, PET, PTAk, RNifti*, RNiftyReg*, Rvcg*,
                tractor.base*, waveslim
                (* = core package)

    --------------------------------------------------------
    MetaAnalysis.ctv :
    ================

    CRAN Task View
    --------------
    Name:       MetaAnalysis
    Topic:      Meta-Analysis
    Maintainer: Michael Dewey
    Contact:    lists@dewey.myzen.co.uk
    Version:    2018-05-10

    Packages:   aggregation, altmeta, bamdit, bayesmeta, bmeta, bspmma, CAMAN,
                CIAAWconsensus, clubSandwich, compute.es, ConfoundedMeta,
                CopulaREMADA, CPBayes, CRTSize, diagmeta, dosresmeta, ecoreg,
                effsize, epiR, esc, etma, exactmeta, extfunnel, forestmodel,
                forestplot, gap, gemtc, getmstatistic, gmeta, hetmeta, ipdmeta,
                joineRmeta, joint.Cox, MAc, MAd, mada, MAVIS,
                MendelianRandomization, meta*, meta4diag, MetaAnalyser, MetABEL,
                metaBMA, metacart, metacor, metafor*, metaforest, metafuse,
                metagear, metagen, metagen, MetaIntegrator, metaLik, metaMA,
                metamisc, metansue, metap, MetaPath, MetaPCA, metaplotr, metaplus,
                MetaQC, metaRNASeq, metaSEM, metasens, MetaSKAT, metatest,
                Metatron, metavcov, metaviz, mmeta, MultiMeta, mvmeta, mvtmeta,
                netmeta, nmaINLA, nmathresh, pcnetmeta, pimeta, psychmeta,
                psychometric, PubBias, RandMeta, ratesci, RBesT, RcmdrPlugin.EZR,
                RcmdrPlugin.RMTCJags, revtools, rma.exact, rmeta, robumeta,
                SAMURAI, SCMA, selectMeta, seqMeta, surrosurv, TFisher, weightr,
                xmeta
                (* = core package)

    --------------------------------------------------------
    ModelDeployment.ctv :
    ===================

    CRAN Task View
    --------------
    Name:       ModelDeployment
    Topic:      Model Deployment with R
    Maintainer: Yuan Tang
    Contact:    terrytangyuan@gmail.com
    Version:    2018-05-01

    Packages:   arules, arulesCBA, arulesSequences, aurelius, AzureML, dbplyr,
                domino, dplyr, FastRWeb, h2o, httpuv, ibmdbR, jug, keras, mleap,
                onnx, opencpu, plumber, pmml, pmmlTransformations, rattle,
                reticulate, RSclient, Rserve, rsparkling, sparklyr, tensorflow,
                tfestimators, tidypredict, xgboost, yhatr
    --------------------------------------------------------
    Multivariate.ctv :
    ================

    CRAN Task View
    --------------
    Name:       Multivariate
    Topic:      Multivariate Statistics
    Maintainer: Paul Hewson
    Contact:    Paul.Hewson@plymouth.ac.uk
    Version:    2018-02-12

    Packages:   abind, ade4*, amap, aplpack, ash, bayesm, ca, calibrate, car,
                caret, class, clue, cluster*, clusterGeneration, clusterSim,
                clustvarsel, clv, cocorresp, concor, copula, corpcor, covRobust,
                cramer, cwhmisc, delt, denpro, desirability, dr, e1071, earth,
                ellipse, energy, eRm, FactoMineR, fastICA, feature, fgac, fpc, fso,
                gclus, GenKern, geometry, geozoo, gmodels, GPArotation, hddplot,
                Hmisc, homals, hybridHclust, ICS, ICSNP, iplots, JADE, kernlab,
                KernSmooth, kknn, klaR, knncat, kohonen, ks, lattice, ltm, mAr,
                MASS*, Matrix, matrixcalc, mclust, MCMCpack, mda, mice, misc3d,
                mitools, mix, mnormt, MNP, monomvn, mvnmle, mvnormtest, mvoutlier,
                mvtnorm, nFactors, pan, paran, party, pcaPP, PearsonICA, pls,
                plsgenomics, poLCA, polycor, ppls, prim, proxy, psy, PTAk, rda,
                relaimpo, rggobi, rgl, robustbase, ROCR, rpart, rrcov,
                scatterplot3d, sem, SensoMineR, seriation, simba, smatr, sn, spam,
                SparseM, SpatialNP, superpc, trimcluster, tsfa, vcd, vegan*, VGAM,
                VIM, xgobi, YaleToolkit
                (* = core package)

    --------------------------------------------------------
    NaturalLanguageProcessing.ctv :
    =============================

    CRAN Task View
    --------------
    Name:       NaturalLanguageProcessing
    Topic:      Natural Language Processing
    Maintainer: Fridolin Wild, Performance Augmentation Lab (PAL, Department of Computing and Communications Technologies, Oxford Brookes University, UK
    Contact:    wild@brookes.ac.uk
    Version:    2017-11-29

    Packages:   alineR, boilerpipeR, corpora, gsubfn, gutenbergr, hunspell,
                kernlab, KoNLP, koRpus, languageR, lda, lsa, maxent, monkeylearn,
                movMF, mscstexta4r, mscsweblm4r, openNLP, ore, phonics, phonics,
                qdap, quanteda, RcmdrPlugin.temis, rel, RKEA, RTextTools, RWeka,
                skmeans, SnowballC, stm, stringdist, stringi, tau, tesseract,
                text2vec, textcat, textir, textrank, textreuse, tidytext, tm*,
                tm.plugin.alceste, tm.plugin.dc, tm.plugin.europresse,
                tm.plugin.factiva, tm.plugin.lexisnexis, tm.plugin.mail,
                tm.plugin.webmining, tokenizers, topicmodels, udpipe, wordcloud,
                wordnet, zipfR
                (* = core package)

    --------------------------------------------------------
    NumericalMathematics.ctv :
    ========================

    CRAN Task View
    --------------
    Name:       NumericalMathematics
    Topic:      Numerical Mathematics
    Maintainer: Hans W. Borchers
    Contact:    hwb@mailbox.org
    Version:    2018-04-26

    Packages:   ADPF, akima, appell, arrangements, BB, Bessel, bigIntegerAlgos,
                Brobdingnag, combinat, conicfit, contfrac, cubature, Deriv*,
                eigeninv, elliptic, expint, expm, fastGHQuad, feather, features,
                findpython, FixedPoint, fourierin, gaussquad, geigen, gmp, gsl,
                hypergeo, interp, irlba, JuliaCall, ktsolve, lamW, logOfGamma,
                magic, MASS, matlab, Matrix*, matrixcalc, MonoPoly, mpoly,
                multipol, mvQuad, nleqslv, numbers, numDeriv*, onion, optR,
                orthopolynom, Pade, partitions, permutations, polyCub, polynom*,
                PolynomF, pracma*, PRIMME, PythonInR, QZ, R.matlab, R2Cuba,
                rARPACK, Rcpp, RcppAlgos, RcppArmadillo, RcppEigen, reticulate,
                Rlinsolve, Rmpfr, rmumps, RootsExtremaInflections, rPython, Rserve,
                RSpectra, rSymPy, Ryacas, schumaker, signal, SimplicialCubature,
                SnakeCharmR, SolveLS, SparseGrid, SparseM, SphericalCubature, ssvd,
                statmod, stinepack, svd, tripack, VeryLargeIntegers, XR, XRJulia,
                XRPython, Zseq
                (* = core package)

    --------------------------------------------------------
    OfficialStatistics.ctv :
    ======================

    CRAN Task View
    --------------
    Name:       OfficialStatistics
    Topic:      Official Statistics & Survey Methodology
    Maintainer: Matthias Templ
    Contact:    matthias.templ@gmail.com
    Version:    2018-04-20

    Packages:   acs, Amelia, BalancedSampling, BIFIEsurvey, CalibrateSSB, cat,
                cbsodataR, censusapi, CoImp, convey, deducorrect, DHS.rates,
                easySdcTable, editrules, errorlocate, eurostat, extremevalues, FFD,
                foreign, Frames2, GeomComb, gridsample, haven, hbsae, Hmisc, IC2,
                icarus, inegiR, ineq, JoSAE, laeken, lavaan, lavaan.survey, lme4,
                mapStats, MatchIt, MBHdesign, memisc, mi, mice, micEconIndex,
                MicSim, MImix, mipfp, missForest, missMDA, mitools, mix, nlme,
                norm, OECD, pan, panelaggregation, pps, PracTools, prevR, pxR,
                quantification, questionr, RcmdrPlugin.sampling, RecordLinkage,
                reweight, Rilostat, robCompositions, rpms, rrcov3way, rrcovNA,
                RRTCS, rsae, rspa, rworldmap, saeSim, samplesize4surveys, sampling,
                samplingbook, SamplingStrata, samplingVarEst, SAScii, SDaA,
                sdcMicro, sdcTable, seasonal, SeleMix, simFrame, simPop,
                simputation, sms, sorvi, spsurvey, srvyr, StatMatch,
                stratification, stringdist, survey*, surveybootstrap, surveydata,
                surveyoutliers, surveyplanning, svyPVpack, synthpop, tabplot,
                TeachingSampling, tmap, treemap, univOutl, validate, validatetools,
                vardpoor, VIM, x12, x12GUI, XBRL, yaImpute
                (* = core package)

    --------------------------------------------------------
    Optimization.ctv :
    ================

    CRAN Task View
    --------------
    Name:       Optimization
    Topic:      Optimization and Mathematical Programming
    Maintainer: Stefan Theussl 
      and Hans W. Borchers
    Contact:    R-optimization@mailbox.org
    Version:    2018-05-07

    Packages:   ABCoptim, adagio, alabama*, BB, boot, bvls, cccp, cec2005benchmark,
                cec2013, CEoptim, clpAPI, CLSOCP, clue, cmaes, cmaesr, colf,
                coneproj, copulaedas, cplexAPI, crs, CVXR, dclone, DEoptim*,
                DEoptimR, desirability, dfoptim*, Dykstra, ECOSolveR, ecr, flacco,
                GA, genalg, GenSA, globalOptTests, glpkAPI, goalprog,
                GrassmannOptim, gsl, hydroPSO, igraph, irace, isotone, kernlab,
                kofnGA, lbfgs, lbfgsb3, limSolve, linprog, localsolver, LowRankQP,
                lpSolve, lpSolveAPI, lsei, ManifoldOptim, matchingMarkets,
                matchingR, maxLik, mcga, mco, metaheuristicOpt, minpack.lm, minqa,
                mize, mknapsack, mlrMBO, n1qn1, neldermead, NlcOptim, nleqslv,
                nlmrt, nloptr, nls2, nlsr, NMOF, nnls, ompr, onls, optimr,
                optimsimplex, optimx, optmatch, parma, powell, pso, psoptim, qap,
                quadprog*, quadprogXT, quantreg, rcdd, RCEIM, Rcgmin, rCMA, Rcplex,
                RcppDE, RcppNumerical, Rcsdp, Rdsdp, rgenoud, Rglpk, rLindo,
                Rmalschains, Rmosek, rneos, ROI, ROI.plugin.qpoases, rosqp, Rsolnp,
                Rsymphony, Rtnmin, Rvmmin, SACOBRA, scs, sdpt3r, smoof, sna, soma,
                subplex, tabuSearch, trust, trustOptim, TSP, ucminf*
                (* = core package)

    --------------------------------------------------------
    Pharmacokinetics.ctv :
    ====================

    CRAN Task View
    --------------
    Name:       Pharmacokinetics
    Topic:      Analysis of Pharmacokinetic Data
    Maintainer: Bill Denney
    Contact:    wdenney@humanpredictions.com
    Version:    2018-05-04

    Packages:   cpk, dfpk, mrgsolve, ncar, nmw, NonCompart, PK, PKgraph,
                PKPDmodels, pkr, PKreport, RxODE, scaRabee
    --------------------------------------------------------
    Phylogenetics.ctv :
    =================

    CRAN Task View
    --------------
    Name:       Phylogenetics
    Topic:      Phylogenetics, Especially Comparative Methods
    Maintainer: Brian O'Meara
    Contact:    omeara.brian@gmail.com
    Version:    2018-02-21

    Packages:   adephylo, adhoc, ape*, apTreeshape, BAMMtools, bayou, betapart,
                BioGeoBEARS, caper, cati, convevol, corHMM, DAMOCLES, DDD,
                dendextend, distory, diversitree, evobiR, expands, expoTree,
                geiger, geomorph, ggplot2, GUniFrac, HMPTrees, HyPhy, idendr0, ips,
                iteRates, jaatha, kdetrees, markophylo, MCMCglmm, metafor, MPSEM,
                mvMORPH, nLTT, ouch, outbreaker, OutbreakTools, OUwie, paleotree,
                paleoTS, pastis, PBD, PCPS, pegas, phangorn, phyclust, phyext2,
                phylobase, phylocanvas, phyloclim, PHYLOGR, phyloland, phylolm,
                phylotools, phyloTop, phytools, pmc, RADami, rdryad, rmetasim,
                rncl, RNeXML, rotl, rphast, Rphylip, SigTree, strap, surface,
                SYNCSA, taxize, TESS, treebase, TreePar, treeplyr, TreeSim, vegan
                (* = core package)

    --------------------------------------------------------
    Psychometrics.ctv :
    =================

    CRAN Task View
    --------------
    Name:       Psychometrics
    Topic:      Psychometric Models and Methods
    Maintainer: Patrick Mair
    Contact:    mair@fas.harvard.edu
    Version:    2018-02-28

    Packages:   ade4*, anacor*, AnalyzeFMRI, aspect, BayesFM, BayesLCA, betareg,
                BigSEM, BiplotGUI, birtr, blavaan*, bpca, BradleyTerry2, BTLLasso,
                ca*, cabootcrs, cacIRT, catR, CAvariants, CDM, cds, cIRT, classify,
                ClustVarLV, CMC, cncaGUI, cocor, cocorresp, cocron, CopyDetect,
                covLCA, ctsem, CTT*, CTTShiny, DAKS, dexter, DFIT, DIFlasso,
                difNLR, difR, DistatisR, dualScale, e1071, eba, ecodist, edstan,
                eegkit, EFAutilities, elasticnet, emIRT, equate, equateIRT, eRm*,
                esaBcv, EstCRM, ExPosition, FactoMineR, faoutlier, fastICA,
                fechner, flexmix, fourPNO, fwdmsa, GDINA, GPArotation, gSEM,
                gtheory, homals*, ica, ICC, immer, influence.SEM, irr, irtDemo,
                irtoys, IRTpp, irtProb, irtrees, IRTShiny, kcirt, kequate,
                KernSmoothIRT, kst, labdsv, latdiag, lava, lava.tobit, lavaan*,
                lavaan.survey, lba, LCAvarsel, lcda, lme4*, LNIRT, lordif, lsl,
                ltbayes, ltm*, LVMMCOR, MASS, MBESS, MCAvariants, mcIRT, MCMCglmm,
                MCMCpack, medflex, mediation, metaSEM, MIIVsem, mirt*, mirtCAT,
                missMDA, mixRasch, MLCIRTwithin, MLCM, MLDS, modelfree, mokken*,
                MplusAutomation, mpt, MPTinR, mRm, MultiLCIRT, multiplex, multiway,
                munfold, nFactors, nlme, nlsem, nsprcomp, OpenMx, optiscale,
                ordinal, pairwise, paran, pathmox, pcaPP, pcIRT, piecewiseSEM, pks,
                PLmixed, plotSEMM, plRasch, pls, plspm, poLCA, polycor, PP,
                prefmod*, profileR, pscl, psy*, psych*, psychometric, psychomix,
                psychotools, psychotree*, psyphy, PTAk, pwrRasch, qgraph,
                QuantPsyc, quickpsy, randomLCA, RaschSampler, regsem, REQS, rpf,
                rsem, sem*, semdiag, semGOF, SEMID, SEMModComp, semPlot, semPLS,
                semTools, semtree, SensoMineR, ShinyItemAnalysis, Sim.DiffProc,
                simsem, sirt, smacof*, smds, SNSequate, soc.ca, SOD,
                SparseFactorAnalysis, sparseSEM, subscore, superMDS, systemfit,
                TAM*, TestDataImputation, TestScorer, ThreeWay, TripleR, vegan*,
                VGAM, wCorr, WrightMap, xgobi, xxIRT
                (* = core package)

    --------------------------------------------------------
    ReproducibleResearch.ctv :
    ========================

    CRAN Task View
    --------------
    Name:       ReproducibleResearch
    Topic:      Reproducible Research
    Maintainer: Max Kuhn
    Contact:    mxkuhn@gmail.com
    Version:    2018-04-18

    Packages:   animation, apaStyle, apsrtable, archivist, ascii, bibtex, brew,
                checkpoint, DT, exams, formatR, formattable, highlight, highr,
                Hmisc*, htmlTable, htmltools, HTMLUtils, humanFormat, hwriter,
                kfigr, knitcitations, knitLatex, knitr*, latex2exp, lazyWeave,
                lubridate, markdown, memisc, miniCRAN, NMOF, packrat, pander,
                papeR, prettyunits, quantreg, R.cache, R.rsp, R2HTML*, R2PPT, R2wd,
                rapport, rbundler, RefManageR, ReporteRs, reporttools, resumer,
                rmarkdown, rms*, rprintf, rtf, SortableHTMLTables, sparktex,
                stargazer, suRtex, SweaveListingUtils, tables, TeachingSampling,
                texreg, tikzDevice, tth, tufterhandout, xtable*, ztable
                (* = core package)

    --------------------------------------------------------
    Robust.ctv :
    ==========

    CRAN Task View
    --------------
    Name:       Robust
    Topic:      Robust Statistical Methods
    Maintainer: Martin Maechler
    Contact:    Martin.Maechler@R-project.org
    Version:    2016-08-29

    Packages:   covRobust, coxrobust, distr, drgee, FRB, georob, GSE, lqmm, MASS*,
                mblm, multinomRob, mvoutlier, OutlierDC, OutlierDM, quantreg,
                RandVar, rgam, roahd, RobAStBase, robcor, robeth, robfilter,
                RobLox, RobLoxBioC, RobPer, RobRex, RobRSVD, robumeta, robust*,
                RobustAFT, robustbase*, robustDA, robustgam, robustlmm,
                robustloggamma, robustreg, robustX, ROptEst, ROptRegTS, rorutadis,
                rrcov*, rrcovHD, rrcovNA, rsig, RSKC, sandwich, ssmrob, TEEReg,
                wle, WRS2
                (* = core package)

    --------------------------------------------------------
    SocialSciences.ctv :
    ==================

    CRAN Task View
    --------------
    Name:       SocialSciences
    Topic:      Statistics for the Social Sciences
    Maintainer: John Fox
    Contact:    jfox@mcmaster.ca
    Version:    2018-05-08

    Packages:   acepack, Amelia, aod, arm, betareg, biglm, BMA, boot*, bootstrap,
                brglm, car*, catspec, class, demography, dispmod, dr, effects*,
                elrm, ergm, exactLoglinTest, gam*, gee, geepack, gmodels, gnm, gss,
                Hmisc*, influence.ME, latentnet, leaps, lme4*, lmeSplines, lmm,
                lmtest*, locfit, logistf, logmult, lsmeans*, MASS*, Matching,
                MatchIt, MCMCglmm*, mgcv*, mi*, mice*, mitools, mix, mlogit, MNP,
                multcomp*, multgee, multinomRob, multiplex, mvnmle, network, nlme*,
                nlstools, nnet*, norm, np, optmatch, PAFit, pan, perturb,
                PSAgraphics, pscl, quantreg*, qvcalc, rms*, RSiena, sandwich*,
                simpleboot, sm, sna, spatial, statnet, survey*, survival*, vcd,
                VGAM*, VIM, visreg, Zelig
                (* = core package)

    --------------------------------------------------------
    Spatial.ctv :
    ===========

    CRAN Task View
    --------------
    Name:       Spatial
    Topic:      Analysis of Spatial Data
    Maintainer: Roger Bivand
    Contact:    Roger.Bivand@nhh.no
    Version:    2018-04-11

    Packages:   ade4, adehabitatHR, adehabitatHS, adehabitatLT, adehabitatMA, ads,
                akima, AMOEBA, ash, aspace, automap, CARBayes, cartography,
                classInt*, cleangeo, CompRandFld, constrainedKriging, cshapes,
                dbmss, DCluster*, deldir*, diseasemapping, DSpat, ecespa,
                ExceedanceTools, fields, FieldSim, FRK, gdalUtils, gdistance, gear,
                geoaxe, geojson, geojsonio, GEOmap, geomapdata, geonames, geoR*,
                geoRglm, georob, geospacom, geosphere, geospt, geostatsp, ggmap,
                ggsn, glmmBUGS, gmt, Grid2Polygons, GriegSmith, gstat*, Guerry,
                GWmodel, gwrr, hdeco, HSAR, igraph, intamap, ipdw, landsat,
                latticeDensity, lawn, lctools, leafletR, magclass, mapdata,
                mapmisc, mapproj, maps, maptools*, mapview, marmap, MBA, McSpatial,
                micromap, ModelMap, ncdf4, ncf, ngspatial, nlme, OasisR,
                OpenStreetMap, osmar, pastecs, PBSmapping, PBSmodelling,
                plotGoogleMaps, plotKML, postGIStools, PReMiuM, ProbitSpatial,
                quickmapr, ramps, RandomFields*, rangeMapper, RArcInfo, raster*,
                rasterVis, RColorBrewer*, recmap, regress, rgbif, rgdal*, rgeos*,
                RgoogleMaps, rgrass7, RNetCDF, rpostgis, RPyGeo, RQGIS, RSAGA,
                RSurvey, rtop, rworldmap, rworldxtra, S2sls, seg, sf*, sgeostat,
                shapefiles, shp2graph, siplab, smacpod, smerc, sp*, spacetime*,
                spacom, spaMM, spanel, sparr, spatgraphs, spatial, spatial.tools,
                spatialCovariance, SpatialEpi, SpatialExtremes, SpatialPosition,
                spatialprobit, spatialsegregation, SpatialTools, spatstat*,
                spatsurv, spBayes, spBayesSurv, spcosa, spdep*, sperrorest,
                spgrass6, spgwr, sphet, spind, splancs*, splm, spm, spmoran,
                spsann, spselect, spsurvey, spTimer, SSN, starma, statebins, Stem,
                stplanr, taRifx, tgp, tmap, trip, tripack, tripEstimation,
                UScensus2000cdp, UScensus2000tract, vardiag, vec2dtransf, vegan,
                Watersheds, wkb
                (* = core package)

    --------------------------------------------------------
    SpatioTemporal.ctv :
    ==================

    CRAN Task View
    --------------
    Name:       SpatioTemporal
    Topic:      Handling and Analyzing Spatio-Temporal Data
    Maintainer: Edzer Pebesma
    Contact:    edzer.pebesma@uni-muenster.de
    Version:    2018-04-20

    Packages:   adehabitatLT*, animalTrack, argosfilter, BayesianAnimalTracker,
                BBMM, bcpa, CARBayesST, crawl, cshapes, ctmcmove, ctmm, diveMove,
                fishmove, FLightR, gapfill, GeoLight, googleVis, gstat*, lgcp,
                lme4, M3, mkde, move, moveHMM, mvtsplot, ncdf4, nlme, openair,
                pastecs, pbdNCDF4, plm, plotKML, RandomFields*, raster*, rasterVis,
                rgl, rmatio, RNetCDF, rsatscan, sf, sigloc, SimilarityMeasures,
                smam, solaR, sp*, spacetime*, spate, SpatioTemporal, spatstat,
                spBayes, sphet, splancs, splm, spTimer, stam, Stem, STMedianPolish,
                stppResid, surveillance*, trackeR, TrackReconstruction, trip*,
                tripEstimation, VTrack, wildlifeDI, xts*
                (* = core package)

    --------------------------------------------------------
    Survival.ctv :
    ============

    CRAN Task View
    --------------
    Name:       Survival
    Topic:      Survival Analysis
    Maintainer: Arthur Allignol and Aurelien Latouche
    Contact:    arthur.allignol@gmail.com
    Version:    2018-05-04

    Packages:   AdapEnetClass, addhazard, AER, aftgee, ahaz, AHR, AIM, APtools,
                asaur, asbio, aster, aster2, BaSTA, BayesPiecewiseICAR, bayesSurv,
                BayHaz, BGPhazard, Biograph, BMA, bnnSurvival, boot, bpcp,
                bshazard, bujar, casebase, censReg, CFC, clinfun, cmprsk*,
                cmprskQR, coarseDataTools, coin, compareC, compeir, compound.Cox,
                concreg, condGEE, condSURV, controlTest, CoxBoost, coxinterval,
                coxme, coxphf, coxphw, CoxRidge, coxrobust, coxsei, CPE, CPHshape,
                Cprob, CR, crrp, crrSC, crrstep, crskdiag, currentSurvival,
                Cyclops, DAAG, dblcens, discSurv, DPpackage, DStree, DTDA,
                dynamichazard, dynfrail, dynpred, dynsurv, eha*, ELYP, emplik,
                emplik2, Epi, epiR, etm, exactRankTests, FamEvent, fastcox,
                fastpseudo, FHtest, fitdistrplus, flexPM, flexrsurv, flexsurv,
                frailtyEM, frailtyHL, frailtypack, frailtySurv, gamboostMSM,
                gamlss.cens, gbm, gcerisk, gems, genSurv, glmnet, glmpath,
                globalboosttest, glrt, gof, GORCure, gss, GSSE, gte, hdnom,
                ICBayes, ICE, icenReg, ICGOR, icRSF, ICsurv, IDPSurvival, imputeYn,
                InformativeCensoring, intccr, intercure, interval, invGauss,
                ipdmeta, ipred, isoph, jackknifeKME, JM, JMbayes, joineR, joineRML,
                joint.Cox, JointModel, JPSurv, kaps, km.ci, kmconfband, kmi,
                KMsurv, landest, lava.tobit, lbiassurv, LearnBayes, LexisPlotR,
                lmec, locfit, logconcens, LogicReg, LogrankA, logspline, lpc,
                lsmeans, lss, LTRCtrees, MAMSE, maxstat, mboost, MCMCglmm,
                MCMCpack, mets, mexhaz, mfp, miCoPTCM, MicSim, MIICD, mixAK,
                mixPHM, MLEcens, MRsurv, msm, msmtools, msSurv, MST, mstate*,
                muhaz*, multcomp, multipleNCC, mvna, NADA, NestedCohort, NPHMC,
                NPMLEcmprsk, npsurv, OIsurv, OrdFacReg, OutlierDC, p3state.msm,
                paf, pamr, parfm, parfm, party, pch, pec, penalized, PenCoxFrail,
                penMSM, peperr, PermAlgo, PHeval, phmm, plac, polspline, popEpi,
                powerSurvEpi, PReMiuM, prodlim, psbcGroup, pseudo, quantreg,
                randomForestSRC, ranger, rankhazard, reda, relsurv, reReg, rhosp,
                riskRegression, risksetROC, rms*, RobustAFT, ROCt, rpart, rsig,
                rstpm2, saws, SemiCompRisks, SemiMarkov, SGL, simexaft, SimHaz,
                simMSM, simPH, SimSCRPiecewise, simsurv, smcure, SMIR,
                SmoothHazard, smoothHR, smoothSurv, SMPracticals, spatstat,
                spatsurv, spBayesSurv, SSRMST, superpc, surv2sampleComp,
                survAccuracyMeasures, survAUC, survC1, SurvCorr, survexp.fr,
                survey, Survgini, survIDINRI, survival*, survivalMPL, survivalROC,
                survJamda, SurvLong, survminer, survMisc, survPresmooth, SurvRank,
                SurvRegCensCov, survRM2, survsim, survSNP, SvyNom, TBSSurvival,
                tdROC, thregI, timereg*, timeROC, tlmec, TP.idm, TPmsm, tpr,
                TraMineR, TransModel, tranSurv, TSHRC, uniah, uniCox, VGAM,
                vitality, YPmodel
                (* = core package)

    --------------------------------------------------------
    TimeSeries.ctv :
    ==============

    CRAN Task View
    --------------
    Name:       TimeSeries
    Topic:      Time Series Analysis
    Maintainer: Rob J Hyndman
    Contact:    Rob.Hyndman@monash.edu
    Version:    2018-04-15

    Packages:   acp, AER, ArDec, arfima, astsa, autovarCore, BAYSTAR, bentcableAR,
                BETS, bfast, bigtime, BigVAR, BNPTSclust, boot, BootPR, brainwaver,
                bspec, bsts, CADFtest, carfima, carx, cents, changepoint, chron,
                cointReg, CommonTrend, costat, cts, dataseries, dCovTS, depmix,
                depmixS4, deseasonalize, dLagM, dlm, dlmodeler, dlnm, dse, dtw,
                dtwclust, dygraphs, dyn, dynlm, dynr, earlywarnings, Ecdat, ecm,
                ecp, EMD, ensembleBMA, EvalEst, events, expsmooth, factorstochvol,
                fame, fanplot, FeedbackTS, fGarch, FGN, FitAR, FitARMA, fma,
                fNonlinear, ForeCA, forecast*, ForecastComb, forecastHybrid,
                forecTheta, fpp, fpp2, fracdiff, fractal, fractalrock, freqdom,
                freqdom.fda, fts, ftsa, funtimes, GAS, gdpc, ggseas, glarma, GMDH,
                gsarima, gtop, HarmonicRegression, hht, hts, hwwntest, imputeTS,
                influxdbr, InspectChangepoint, itsmr, jmotif, KFAS, KFKSDS, kza,
                locits, lomb, LPStimeSeries, LSTS, ltsa, lubridate, mafs, MAPA,
                mAr, MAR1, mar1s, MARSS, mclcar, Mcomp, meboot, mFilter, mlVAR,
                mondate, MSBVAR, MSwM, MTS, mtsdi, MultipleBubbles, multitaper,
                mvcwt, nardl, nets, nlts, nnfor, nonlinearTseries, npst, odpc,
                opera, orderedLasso, paleoTS, partsm, pastecs, PCA4TS, pcdpca, pdc,
                pdfetch, pear, perARMA, pomp, portes, prophet, psd, PSF, ptw,
                Quandl, quantspec, rbcb, rdatamarket, RGENERATE, Rlibeemd, rmaf,
                RMAWGEN, robets, robfilter, robustarima, roll, RSEIS, Rssa, rts,
                rucrdtw, rugarch, rwt, sae2, scoringRules, SDD, sde, seas, season,
                seasonal, seasonalview, Sim.DiffProc, sleekts, smooth, sparsevar,
                spectral, spectral.methods, spTimer, stlplus, stochvol, stR,
                strucchange, stsm, stsm.class, surveillance, svars, sweep, Tcomp,
                TED, tempdisagg, tframe, thief, tibbletime, Tides, tiger, timeDate,
                TimeProjection, timesboot, timeSeries, timetk, timsac, tis, tpr,
                trend, TSA, tsbugs, TSclust, tscount, TSdbi, tsdecomp, tsdisagg2,
                TSdist, tsDyn, tseries*, tseriesChaos, tseriesEntropy, tsfa,
                tsfknn, tsibble, tsintermittent, TSMining, tsModel, tsoutliers,
                tsPI, TSrepr, TSstudio, TSTutorial, tswge, urca, uroot, VAR.etp,
                vars, VARsignR, Wats, WaveletComp, wavelets, waveslim, wavethresh,
                WeightedPortTest, wktmo, wmtsa, x12, x12GUI, x13binary, xts, yuima,
                ZIM, zoo*, ZRA
                (* = core package)

    --------------------------------------------------------
    WebTechnologies.ctv :
    ===================

    CRAN Task View
    --------------
    Name:       WebTechnologies
    Topic:      Web Technologies and Services
    Maintainer: Thomas Leeper, Scott Chamberlain, Patrick Mair, Karthik Ram, Christopher Gandrud
    Contact:    thosjleeper@gmail.com
    Version:    2018-05-10

    Packages:   abbyyR, aRxiv, aws.signature, AzureML, backblazer, bigml,
                bigrquery, boilerpipeR, boxr, captr, clarifai, curl, cymruservices,
                d3Network, datamart, dataone, datarobot, ddeploy, discgolf,
                downloader, dvn, europepmc, factualR, FastRWeb, fbRads, fiery,
                fitbitScraper, GAR, genderizeR, geoparser, ggmap, ggvis, gistr,
                git2r, gitlabr, gmailr, googleAnalyticsR, googleCloudStorageR,
                googlesheets, googleVis, graphTweets, gsheet, gtrendsR, htmltab,
                httpcache, httping, httpRequest, httpuv, httr*, imguR, instaR, jqr,
                jSonarR, jsonlite*, jsonvalidate, jug, leafletR, livechatR,
                longurl, lucr, magrittr, mailR, mime, mscstexta4r, mscsweblm4r,
                MTurkR, oai, OAIHarvester, openadds, opencage, opencpu, osi, osmar,
                osmplotr, osrm, pdftables, placement, plotGoogleMaps, plotKML,
                plumber, plusser, pubmed.mineR, pushoverr, RAdwords, randNames,
                rapport, Rbitcoin, rbitcoinchartsapi, Rblpapi, RCurl*, RDataCanvas,
                rdatacite, redcapAPI, RefManageR, repmis, request, restimizeapi,
                Rexperigen, Rfacebook, rfigshare, RForcecom, RGA, rgeolocate,
                RGoogleFit, RgoogleMaps, rio, rjson, RJSONIO, rLTP, RMixpanel,
                ROAuth, Rook, ROpenFIGI, ROpenWeatherMap, rorcid, rosetteApi,
                rplos, RPushbullet, rrefine, RSclient, rsdmx, Rserve,
                RSiteCatalyst, RSmartlyIO, RStripe, rvest, RYandexTranslate,
                RZabbix, scholar, scrapeR, searchConsoleR, sendmailR, servr,
                shiny*, shopifyr, slackr, SocialMediaMineR, soql, streamR,
                telegram, threewords, tm.plugin.webmining, transcribeR, translate,
                translateR, tumblR, tweet2r, twitteR, uaparserjs, urlshorteneR,
                urltools, V8, W3CMarkupValidator, webreadr, webshot, webutils,
                whisker, WikidataR, wikipediatrend, WikipediR, WikiSocio, WufooR,
                XML*, xml2, XML2R, yhatr, yummlyr, zendeskR
                (* = core package)

    --------------------------------------------------------
