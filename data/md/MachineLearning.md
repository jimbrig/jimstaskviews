## CRAN Task View: Machine Learning & Statistical Learning

  ----------------- ---------------------------------------------------
  **Maintainer:**   Torsten Hothorn
  **Contact:**      Torsten.Hothorn at R-project.org
  **Version:**      2018-04-27
  **URL:**          <https://CRAN.R-project.org/view=MachineLearning>
  ----------------- ---------------------------------------------------

<div>

Several add-on packages implement ideas and methods developed at the
borderline between computer science and statistics - this field of
research is usually referred to as machine learning. The packages can be
roughly structured into the following topics:

-   *Neural Networks and Deep Learning* : Single-hidden-layer neural
    network are implemented in package
    [nnet](../packages/nnet/index.html) (shipped with base R). Package
    [RSNNS](../packages/RSNNS/index.html) offers an interface to the
    Stuttgart Neural Network Simulator (SNNS). An interface to the FCNN
    library allows user-extensible artificial neural networks in package
    [FCNN4R](../packages/FCNN4R/index.html).
    [rnn](../packages/rnn/index.html) implements recurrent neural
    networks. Packages implementing deep learning flavours of neural
    networks include [deepnet](../packages/deepnet/index.html)
    (feed-forward neural network, restricted Boltzmann machine, deep
    belief network, stacked autoencoders),
    [RcppDL](../packages/RcppDL/index.html) (denoising autoencoder,
    stacked denoising autoencoder, restricted Boltzmann machine, deep
    belief network) and [h2o](../packages/h2o/index.html) (feed-forward
    neural network, deep autoencoders). An interface to
    [tensorflow](http://www.tensorflow.org) is available in
    [tensorflow](../packages/tensorflow/index.html).
-   *Recursive Partitioning* : Tree-structured models for regression,
    classification and survival analysis, following the ideas in the
    CART book, are implemented in [rpart](../packages/rpart/index.html)
    (shipped with base R) and [tree](../packages/tree/index.html).
    Package [rpart](../packages/rpart/index.html) is recommended for
    computing CART-like trees. A rich toolbox of partitioning algorithms
    is available in [Weka](http://www.cs.waikato.ac.nz/~ml/weka/) ,
    package [RWeka](../packages/RWeka/index.html) provides an interface
    to this implementation, including the J4.8-variant of C4.5 and M5.
    The [Cubist](../packages/Cubist/index.html) package fits rule-based
    models (similar to trees) with linear regression models in the
    terminal leaves, instance-based corrections and boosting. The
    [C50](../packages/C50/index.html) package can fit C5.0
    classification trees, rule-based models, and boosted versions of
    these.\
    Two recursive partitioning algorithms with unbiased variable
    selection and statistical stopping criterion are implemented in
    package [party](../packages/party/index.html) and
    [partykit](../packages/partykit/index.html). Function `ctree()` is
    based on non-parametric conditional inference procedures for testing
    independence between response and each input variable whereas
    `mob()` can be used to partition parametric models. Package
    [model4you](../packages/model4you/index.html) can be used to build
    trees based on more complex models, for example models featuring a
    treatment effect to be partitioned. Transformation trees for
    estimating discrete and continuous predictive distributions based on
    transformation models, also allowing censoring and truncation, are
    available from package [trtf](../packages/trtf/index.html).
    Extensible tools for visualizing binary trees and node distributions
    of the response are available in package
    [party](../packages/party/index.html) and
    [partykit](../packages/partykit/index.html) as well.\
    Tree-structured varying coefficient models are implemented in
    package [vcrpart](../packages/vcrpart/index.html).\
    For problems with binary input variables the package
    [LogicReg](../packages/LogicReg/index.html) implements logic
    regression. Graphical tools for the visualization of trees are
    available in package [maptree](../packages/maptree/index.html).\
    Trees for modelling longitudinal data by means of random effects is
    offered by package [REEMtree](../packages/REEMtree/index.html).
    Partitioning of mixture models is performed by
    [RPMM](../packages/RPMM/index.html).\
    Computational infrastructure for representing trees and unified
    methods for prediction and visualization is implemented in
    [partykit](../packages/partykit/index.html). This infrastructure is
    used by package [evtree](../packages/evtree/index.html) to implement
    evolutionary learning of globally optimal trees. Survival trees are
    available in various package,
    [LTRCtrees](../packages/LTRCtrees/index.html) allows for
    left-truncation and interval-censoring in addition to
    right-censoring.
-   *Random Forests* : The reference implementation of the random forest
    algorithm for regression and classification is available in package
    [randomForest](../packages/randomForest/index.html). Package
    [ipred](../packages/ipred/index.html) has bagging for regression,
    classification and survival analysis as well as bundling, a
    combination of multiple models via ensemble learning. In addition, a
    random forest variant for response variables measured at arbitrary
    scales based on conditional inference trees is implemented in
    package [party](../packages/party/index.html).
    [randomForestSRC](../packages/randomForestSRC/index.html) implements
    a unified treatment of Breiman\'s random forests for survival,
    regression and classification problems. Quantile regression forests
    [quantregForest](../packages/quantregForest/index.html) allow to
    regress quantiles of a numeric response on exploratory variables via
    a random forest approach. For binary data,
    [LogicForest](../packages/LogicForest/index.html) is a forest of
    logic regression trees (package
    [LogicReg](../packages/LogicReg/index.html). The
    [varSelRF](../packages/varSelRF/index.html) and
    [Boruta](../packages/Boruta/index.html) packages focus on variable
    selection by means for random forest algorithms. In addition,
    packages [ranger](../packages/ranger/index.html) and
    [Rborist](../packages/Rborist/index.html) offer R interfaces to fast
    C++ implementations of random forests. Reinforcement Learning Trees,
    featuring splits in variables which will be important down the tree,
    are implemented in package [RLT](../packages/RLT/index.html).
    [wsrf](../packages/wsrf/index.html) implements an alternative
    variable weighting method for variable subspace selection in place
    of the traditional random variable sampling. Random forests for
    parametric models, including forests for the estimation of
    predictive distributions, are available in packages
    [trtf](../packages/trtf/index.html) (predictive transformation
    forests, possibly under censoring and trunction),
    [model4you](../packages/model4you/index.html) (featuring a simple
    user interface for building forests based on arbitrary parametric
    models fitted in R), and [grf](../packages/grf/index.html) (an
    implementation of generalised random forests).
-   *Regularized and Shrinkage Methods* : Regression models with some
    constraint on the parameter estimates can be fitted with the
    [lasso2](../packages/lasso2/index.html) and
    [lars](../packages/lars/index.html) packages. Lasso with
    simultaneous updates for groups of parameters (groupwise lasso) is
    available in package [grplasso](../packages/grplasso/index.html);
    the [grpreg](../packages/grpreg/index.html) package implements a
    number of other group penalization models, such as group MCP and
    group SCAD. The L1 regularization path for generalized linear models
    and Cox models can be obtained from functions available in package
    [glmpath](../packages/glmpath/index.html), the entire lasso or
    elastic-net regularization path (also in
    [elasticnet](../packages/elasticnet/index.html)) for linear
    regression, logistic and multinomial regression models can be
    obtained from package [glmnet](../packages/glmnet/index.html). The
    [penalized](../packages/penalized/index.html) package provides an
    alternative implementation of lasso (L1) and ridge (L2) penalized
    regression models (both GLM and Cox models). Package
    [biglasso](../packages/biglasso/index.html) fits Gaussian and
    logistic linear models under L1 penalty when the data can\'t be
    stored in RAM. Package [RXshrink](../packages/RXshrink/index.html)
    can be used to identify and display TRACEs for a specified shrinkage
    path and to determine the appropriate extent of shrinkage.
    Semiparametric additive hazards models under lasso penalties are
    offered by package [ahaz](../packages/ahaz/index.html). A
    generalisation of the Lasso shrinkage technique for linear
    regression is called relaxed lasso and is available in package
    [relaxo](../packages/relaxo/index.html). Fisher\'s LDA projection
    with an optional LASSO penalty to produce sparse solutions is
    implemented in package
    [penalizedLDA](../packages/penalizedLDA/index.html). The shrunken
    centroids classifier and utilities for gene expression analyses are
    implemented in package [pamr](../packages/pamr/index.html). An
    implementation of multivariate adaptive regression splines is
    available in package [earth](../packages/earth/index.html). Variable
    selection through clone selection in SVMs in penalized models (SCAD
    or L1 penalties) is implemented in package
    [penalizedSVM](../packages/penalizedSVM/index.html). Various forms
    of penalized discriminant analysis are implemented in packages
    [hda](../packages/hda/index.html),
    [rda](../packages/rda/index.html), and
    [sda](../packages/sda/index.html). Package
    [LiblineaR](../packages/LiblineaR/index.html) offers an interface to
    the LIBLINEAR library. The [ncvreg](../packages/ncvreg/index.html)
    package fits linear and logistic regression models under the the
    SCAD and MCP regression penalties using a coordinate descent
    algorithm. High-throughput ridge regression (i.e., penalization with
    many predictor variables) and heteroskedastic effects models are the
    focus of the [bigRR](../packages/bigRR/index.html) package. An
    implementation of bundle methods for regularized risk minimization
    is available form package [bmrm](../packages/bmrm/index.html). The
    Lasso under non-Gaussian and heteroscedastic errors is estimated by
    [hdm](../packages/hdm/index.html), inference on low-dimensional
    components of Lasso regression and of estimated treatment effects in
    a high-dimensional setting are also contained. Package
    [SIS](../packages/SIS/index.html) implements sure independence
    screening in generalised linear and Cox models. Normal and binary
    logistic linear models under various penalties (or mixtures of
    those) can be estimated using package
    [oem](../packages/oem/index.html).
-   *Boosting and Gradient Descent* : Various forms of gradient boosting
    are implemented in package [gbm](../packages/gbm/index.html)
    (tree-based functional gradient descent boosting). Package
    [xgboost](../packages/xgboost/index.html) implements tree-based
    boosting using efficient trees as base learners for several and also
    user-defined objective functions. The Hinge-loss is optimized by the
    boosting implementation in package
    [bst](../packages/bst/index.html). Package
    [GAMBoost](../packages/GAMBoost/index.html) can be used to fit
    generalized additive models by a boosting algorithm. An extensible
    boosting framework for generalized linear, additive and
    nonparametric models is available in package
    [mboost](../packages/mboost/index.html). Likelihood-based boosting
    for Cox models is implemented in
    [CoxBoost](../packages/CoxBoost/index.html) and for mixed models in
    [GMMBoost](../packages/GMMBoost/index.html). GAMLSS models can be
    fitted using boosting by
    [gamboostLSS](../packages/gamboostLSS/index.html). An implementation
    of various learning algorithms based on Gradient Descent for dealing
    with regression tasks is available in package
    [gradDescent](../packages/gradDescent/index.html).
-   *Support Vector Machines and Kernel Methods* : The function `svm()`
    from [e1071](../packages/e1071/index.html) offers an interface to
    the LIBSVM library and package
    [kernlab](../packages/kernlab/index.html) implements a flexible
    framework for kernel learning (including SVMs, RVMs and other kernel
    learning algorithms). An interface to the SVMlight implementation
    (only for one-against-all classification) is provided in package
    [klaR](../packages/klaR/index.html). The relevant dimension in
    kernel feature spaces can be estimated using
    [rdetools](../packages/rdetools/index.html) which also offers
    procedures for model selection and prediction.
-   *Bayesian Methods* : Bayesian Additive Regression Trees (BART),
    where the final model is defined in terms of the sum over many weak
    learners (not unlike ensemble methods), are implemented in packages
    [BayesTree](../packages/BayesTree/index.html),
    [BART](../packages/BART/index.html), and
    [bartMachine](../packages/bartMachine/index.html). Bayesian
    nonstationary, semiparametric nonlinear regression and design by
    treed Gaussian processes including Bayesian CART and treed linear
    models are made available by package
    [tgp](../packages/tgp/index.html). [MXM](../packages/MXM/index.html)
    implements variable selection based on Bayesian networks.
-   *Optimization using Genetic Algorithms* : Package
    [rgenoud](../packages/rgenoud/index.html) offers optimization
    routines based on genetic algorithms. The package
    [Rmalschains](../packages/Rmalschains/index.html) implements memetic
    algorithms with local search chains, which are a special type of
    evolutionary algorithms, combining a steady state genetic algorithm
    with local search for real-valued parameter optimization.
-   *Association Rules* : Package
    [arules](../packages/arules/index.html) provides both data
    structures for efficient handling of sparse binary data as well as
    interfaces to implementations of Apriori and Eclat for mining
    frequent itemsets, maximal frequent itemsets, closed frequent
    itemsets and association rules. Package
    [opusminer](../packages/opusminer/index.html) provides an interface
    to the OPUS Miner algorithm (implemented in C++) for finding the key
    associations in transaction data efficiently, in the form of
    self-sufficient itemsets, using either leverage or lift.
-   *Fuzzy Rule-based Systems* : Package
    [frbs](../packages/frbs/index.html) implements a host of standard
    methods for learning fuzzy rule-based systems from data for
    regression and classification. Package
    [RoughSets](../packages/RoughSets/index.html) provides comprehensive
    implementations of the rough set theory (RST) and the fuzzy rough
    set theory (FRST) in a single package.
-   *Model selection and validation* : Package
    [e1071](../packages/e1071/index.html) has function `tune()` for
    hyper parameter tuning and function `errorest()`
    ([ipred](../packages/ipred/index.html)) can be used for error rate
    estimation. The cost parameter C for support vector machines can be
    chosen utilizing the functionality of package
    [svmpath](../packages/svmpath/index.html). Functions for ROC
    analysis and other visualisation techniques for comparing candidate
    classifiers are available from package
    [ROCR](../packages/ROCR/index.html). Packages
    [hdi](../packages/hdi/index.html) and
    [stabs](../packages/stabs/index.html) implement stability selection
    for a range of models, [hdi](../packages/hdi/index.html) also offers
    other inference procedures in high-dimensional models.
-   *Other procedures* : Evidential classifiers quantify the uncertainty
    about the class of a test pattern using a Dempster-Shafer mass
    function in package [evclass](../packages/evclass/index.html). The
    [OneR](../packages/OneR/index.html) (One Rule) package offers a
    classification algorithm with enhancements for sophisticated
    handling of missing values and numeric data together with extensive
    diagnostic functions. [spa](../packages/spa/index.html) combines
    feature-based and graph-based data for prediction of some response.
-   *Meta packages* : Package [caret](../packages/caret/index.html)
    provides miscellaneous functions for building predictive models,
    including parameter tuning and variable importance measures. The
    package can be used with various parallel implementations (e.g. MPI,
    NWS etc). In a similar spirit, package
    [mlr](../packages/mlr/index.html) offers a high-level interface to
    various statistical and machine learning packages. Package
    [SuperLearner](../packages/SuperLearner/index.html) implements a
    similar toolbox. The [h2o](../packages/h2o/index.html) package
    implements a general purpose machine learning platform that has
    scalable implementations of many popular algorithms such as random
    forest, GBM, GLM (with elastic net regularization), and deep
    learning (feedforward multilayer networks), among others.
-   *Elements of Statistical Learning* : Data sets, functions and
    examples from the book [The Elements of Statistical Learning: Data
    Mining, Inference, and
    Prediction](https://web.stanford.edu/~hastie/ElemStatLearn/) by
    Trevor Hastie, Robert Tibshirani and Jerome Friedman have been
    packaged and are available as
    [ElemStatLearn](../packages/ElemStatLearn/index.html).
-   *GUI* [rattle](../packages/rattle/index.html) is a graphical user
    interface for data mining in R.
-   *Visualisation (initially contributed by Brandon Greenwell)* The
    `stats::termplot()` function package can be used to plot the terms
    in a model whose predict method supports `type="terms"`. The
    [effects](../packages/effects/index.html) package provides graphical
    and tabular effect displays for models with a linear predictor
    (e.g., linear and generalized linear models). Friedman's partial
    dependence plots (PDPs), that are low dimensional graphical
    renderings of the prediction function, are implemented in a few
    packages. [gbm](../packages/gbm/index.html),
    [randomForest](../packages/randomForest/index.html) and
    [randomForestSRC](../packages/randomForestSRC/index.html) provide
    their own functions for displaying PDPs, but are limited to the
    models fit with those packages (the function `partialPlot` from
    [randomForest](../packages/randomForest/index.html) is more limited
    since it only allows for one predictor at a time). Packages
    [pdp](../packages/pdp/index.html),
    [plotmo](../packages/plotmo/index.html), and
    [ICEbox](../packages/ICEbox/index.html) are more general and allow
    for the creation of PDPs for a wide variety of machine learning
    models (e.g., random forests, support vector machines, etc.); both
    [pdp](../packages/pdp/index.html) and
    [plotmo](../packages/plotmo/index.html) support multivariate
    displays ([plotmo](../packages/plotmo/index.html) is limited to two
    predictors while [pdp](../packages/pdp/index.html) uses trellis
    graphics to display PDPs involving three predictors). By default,
    [plotmo](../packages/plotmo/index.html) fixes the background
    variables at their medians (or first level for factors) which is
    faster than constructing PDPs but incorporates less information.
    [ICEbox](../packages/ICEbox/index.html) focuses on constructing
    individual conditional expectation (ICE) curves, a refinement over
    Friedman\'s PDPs. ICE curves, as well as centered ICE curves can
    also be constructed with the `partial()` function from the
    [pdp](../packages/pdp/index.html) package.
    [ggRandomForests](../packages/ggRandomForests/index.html) provides
    ggplot2-based tools for the graphical exploration of random forest
    models (e.g., variable importance plots and PDPs) from the
    [randomForest](../packages/randomForest/index.html) and
    [randomForestSRC](../packages/randomForestSRC/index.html) packages.

[CORElearn](../packages/CORElearn/index.html) implements a rather broad
class of machine learning algorithms, such as nearest neighbors, trees,
random forests, and several feature selection methods. Similar, package
[rminer](../packages/rminer/index.html) interfaces several learning
algorithms implemented in other packages and computes several
performance measures.

</div>

### CRAN packages:

-   [ahaz](../packages/ahaz/index.html)
-   [arules](../packages/arules/index.html)
-   [BART](../packages/BART/index.html)
-   [bartMachine](../packages/bartMachine/index.html)
-   [BayesTree](../packages/BayesTree/index.html)
-   [biglasso](../packages/biglasso/index.html)
-   [bigRR](../packages/bigRR/index.html)
-   [bmrm](../packages/bmrm/index.html)
-   [Boruta](../packages/Boruta/index.html)
-   [bst](../packages/bst/index.html)
-   [C50](../packages/C50/index.html)
-   [caret](../packages/caret/index.html)
-   [CORElearn](../packages/CORElearn/index.html)
-   [CoxBoost](../packages/CoxBoost/index.html)
-   [Cubist](../packages/Cubist/index.html)
-   [deepnet](../packages/deepnet/index.html)
-   [e1071](../packages/e1071/index.html) (core)
-   [earth](../packages/earth/index.html)
-   [effects](../packages/effects/index.html)
-   [elasticnet](../packages/elasticnet/index.html)
-   [ElemStatLearn](../packages/ElemStatLearn/index.html)
-   [evclass](../packages/evclass/index.html)
-   [evtree](../packages/evtree/index.html)
-   [FCNN4R](../packages/FCNN4R/index.html)
-   [frbs](../packages/frbs/index.html)
-   [GAMBoost](../packages/GAMBoost/index.html)
-   [gamboostLSS](../packages/gamboostLSS/index.html)
-   [gbm](../packages/gbm/index.html) (core)
-   [ggRandomForests](../packages/ggRandomForests/index.html)
-   [glmnet](../packages/glmnet/index.html)
-   [glmpath](../packages/glmpath/index.html)
-   [GMMBoost](../packages/GMMBoost/index.html)
-   [gradDescent](../packages/gradDescent/index.html)
-   [grf](../packages/grf/index.html)
-   [grplasso](../packages/grplasso/index.html)
-   [grpreg](../packages/grpreg/index.html)
-   [h2o](../packages/h2o/index.html)
-   [hda](../packages/hda/index.html)
-   [hdi](../packages/hdi/index.html)
-   [hdm](../packages/hdm/index.html)
-   [ICEbox](../packages/ICEbox/index.html)
-   [ipred](../packages/ipred/index.html)
-   [kernlab](../packages/kernlab/index.html) (core)
-   [klaR](../packages/klaR/index.html)
-   [lars](../packages/lars/index.html)
-   [lasso2](../packages/lasso2/index.html)
-   [LiblineaR](../packages/LiblineaR/index.html)
-   [LogicForest](../packages/LogicForest/index.html)
-   [LogicReg](../packages/LogicReg/index.html)
-   [LTRCtrees](../packages/LTRCtrees/index.html)
-   [maptree](../packages/maptree/index.html)
-   [mboost](../packages/mboost/index.html) (core)
-   [mlr](../packages/mlr/index.html)
-   [model4you](../packages/model4you/index.html)
-   [MXM](../packages/MXM/index.html)
-   [ncvreg](../packages/ncvreg/index.html)
-   [nnet](../packages/nnet/index.html) (core)
-   [oem](../packages/oem/index.html)
-   [OneR](../packages/OneR/index.html)
-   [opusminer](../packages/opusminer/index.html)
-   [pamr](../packages/pamr/index.html)
-   [party](../packages/party/index.html)
-   [partykit](../packages/partykit/index.html)
-   [pdp](../packages/pdp/index.html)
-   [penalized](../packages/penalized/index.html)
-   [penalizedLDA](../packages/penalizedLDA/index.html)
-   [penalizedSVM](../packages/penalizedSVM/index.html)
-   [plotmo](../packages/plotmo/index.html)
-   [quantregForest](../packages/quantregForest/index.html)
-   [randomForest](../packages/randomForest/index.html) (core)
-   [randomForestSRC](../packages/randomForestSRC/index.html)
-   [ranger](../packages/ranger/index.html)
-   [rattle](../packages/rattle/index.html)
-   [Rborist](../packages/Rborist/index.html)
-   [RcppDL](../packages/RcppDL/index.html)
-   [rda](../packages/rda/index.html)
-   [rdetools](../packages/rdetools/index.html)
-   [REEMtree](../packages/REEMtree/index.html)
-   [relaxo](../packages/relaxo/index.html)
-   [rgenoud](../packages/rgenoud/index.html)
-   [RLT](../packages/RLT/index.html)
-   [Rmalschains](../packages/Rmalschains/index.html)
-   [rminer](../packages/rminer/index.html)
-   [rnn](../packages/rnn/index.html)
-   [ROCR](../packages/ROCR/index.html)
-   [RoughSets](../packages/RoughSets/index.html)
-   [rpart](../packages/rpart/index.html) (core)
-   [RPMM](../packages/RPMM/index.html)
-   [RSNNS](../packages/RSNNS/index.html)
-   [RWeka](../packages/RWeka/index.html)
-   [RXshrink](../packages/RXshrink/index.html)
-   [sda](../packages/sda/index.html)
-   [SIS](../packages/SIS/index.html)
-   [spa](../packages/spa/index.html)
-   [stabs](../packages/stabs/index.html)
-   [SuperLearner](../packages/SuperLearner/index.html)
-   [svmpath](../packages/svmpath/index.html)
-   [tensorflow](../packages/tensorflow/index.html)
-   [tgp](../packages/tgp/index.html)
-   [tree](../packages/tree/index.html)
-   [trtf](../packages/trtf/index.html)
-   [varSelRF](../packages/varSelRF/index.html)
-   [vcrpart](../packages/vcrpart/index.html)
-   [wsrf](../packages/wsrf/index.html)
-   [xgboost](../packages/xgboost/index.html)

### Related links:

-   [MLOSS: Machine Learning Open Source
    Software](http://www.MLOSS.org/)
