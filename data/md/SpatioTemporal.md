## CRAN Task View: Handling and Analyzing Spatio-Temporal Data

  ----------------- --------------------------------------------------
  **Maintainer:**   Edzer Pebesma
  **Contact:**      edzer.pebesma at uni-muenster.de
  **Version:**      2018-04-20
  **URL:**          <https://CRAN.R-project.org/view=SpatioTemporal>
  ----------------- --------------------------------------------------

<div>

This task view aims at presenting R packages that are useful for the
analysis of spatio-temporal data.

Please let the [maintainer](mailto:edzer.pebesma@uni-muenster.de) know
if something is inaccurate or missing.

The following people contributed to this task view: Roger Bivand, Achim
Zeileis, Michael Sumner, Ping Yang.

Although one could argue that all data are spatio-temporal, as they must
have been taken somewhere and at some point in time, in many cases the
spatial locations or times of observation are not registered, and
irrelevant to the purpose of the study. Here, we will address the cases
where both location *and* time of observation are registered, and
relevant for the analysis of the data. The [Spatial](Spatial.html) and
[TimeSeries](TimeSeries.html) task views shed light on spatial, and
temporal data handling and analysis, individually.

**Representing data**

-   **In long tables:** In some cases, spatio-temporal data can be held
    in tables ( `data.frame` objects), with longitude, latitude and time
    as three of the columns, or an identifier for a location or region
    and time as columns. For instance, data sets in package
    [plm](../packages/plm/index.html) for linear panel models have
    repeated observations for observational units, where these units
    often refer to spatial areas (countries, states) by an index. This
    index (a name, or number) can be matched to the spatial coordinates
    (polygons) of the corresponding area, an example of this is given by
    [Pebesma (2012, Journal of Statistical
    Software)](http://www.jstatsoft.org/v51/i07/) . As these data sets
    usually contain more than one attribute, to hold the data in a
    two-dimensional table a *long table* form is chosen, where each
    record contains the index of the observational unit, observation
    time, and all attributes.
-   **In time-wide tables:** When a single attribute is considered,
    another layout is that of the *time-wide table* , where each
    observational unit forms a record and each column an observation
    time. [googleVis](../packages/googleVis/index.html) lets you analyze
    such data in a way similar to gapminder (see links).
-   **In space-wide tables:** An example of a space-wide table is the
    Irish wind data set, obtained by `data(wind)` in package
    [gstat](../packages/gstat/index.html). It has time series as
    different columns, each column representing one location (weather
    station). The `stConstruct` function in package
    [spacetime](../packages/spacetime/index.html) accepts data in long,
    time-wide or space-wide tables.
-   **Generic classes:** Formal classes for spatio-temporal data in R
    are provided by the [spacetime](../packages/spacetime/index.html)
    package, which offers S4 classes for full space-time grids (every
    observational unit contains an observation for each observation
    time), sparse space-time grids (regular, but incomplete grids),
    irregular space-time data (each observational unit is observed at
    its own time), and has limited support for trajectory data.
    [spacetime](../packages/spacetime/index.html) classes have
    [sp](../packages/sp/index.html) and
    [xts](../packages/xts/index.html) objects as slots for the spatial
    and temporal components, and can deal with all spatial classes
    (points, lines, polygons, grids) of [sp](../packages/sp/index.html),
    regular and irregular time series, and extend the powerful methods
    (selection, aggregation, plotting coercion) from both packages.
-   **Dedicated classes:** dedicated classes are offered for:
    -   **Geostatistical data:** Package
        [SpatioTemporal](../packages/SpatioTemporal/index.html) offers
        an S3 class `STdata` which holds point observations and
        covariates that can vary in space, time, and space-time, with
        the aim of fitting and predicting a particular class of
        spatio-temporal models, described in its vignettes.
    -   **Gridded/raster data:** package
        [raster](../packages/raster/index.html) deals with sets of
        rasters (called bricks, or stacks), and a set may reflect a
        temporal sequence (use `setZ` on a brick or stack).
    -   **Lattice data:** package
        [surveillance](../packages/surveillance/index.html) provides a
        class `sts`, which holds a `SpatialPolygonsDataFrame` slot for
        the areas, and numeric slots to define a regular time series (no
        time objects, such as `POSIXct`).
    -   **Point patterns:** Package
        [stppResid](../packages/stppResid/index.html) provides a class
        `stwin` for a space-time cuboid, defining a (rectangular)
        space-time window, and class `stpp` for a spatio-temporal point
        pattern (including window). Package
        [spatstat](../packages/spatstat/index.html) provides a class
        `ppx` that deals spatial and temporal coordinate. None of the
        point pattern classes mentioned support spatial or explicit
        temporal reference systems.
    -   **Trajectory data:** Package
        [adehabitatLT](../packages/adehabitatLT/index.html) offers a
        class `ltraj` for trajectories, and methods for analyzing them;
        the packages [move](../packages/move/index.html) and
        [trip](../packages/trip/index.html) both extend
        [sp](../packages/sp/index.html) based classes for trajectories.
        A blog post on [tidy storm
        trajectories](http://r-spatial.org/r/2017/08/28/nest.html)
        points out how nested dataframes, along with geometry list
        columns of the [sf](../packages/sf/index.html) package, can be
        used to model sets of trajectories, and visualise properties at
        the set level and at the level of individual fixes.

**Analyzing data**

-   **Geostatistical data**
    -   [gstat](../packages/gstat/index.html) provides kriging, methods
        of moments variogram estimation and model fitting for a limited
        range of spatio-temporal models.
    -   [RandomFields](../packages/RandomFields/index.html) provides
        kriging, conditional simulation, and covariance functions and
        maximum likelihood function fitting for a very wide range of
        spatio-temporal covariance models.
    -   the [spTimer](../packages/spTimer/index.html) package is able to
        fit, spatially predict and temporally forecast large amounts of
        space-time data using Bayesian Gaussian Process (GP) Models,
        Bayesian Auto-Regressive (AR) Models, and Bayesian Gaussian
        Predictive Processes (GPP) based AR Models.
    -   Package [SpatioTemporal](../packages/SpatioTemporal/index.html)
        fits and predicts a particular class of spatio-temporal models,
        described in detail in its vignettes.
    -   [spBayes](../packages/spBayes/index.html) provides functions for
        fitting Bayesian dynamic space-time regression models for
        settings where space is viewed as continuous but time is taken
        to be discrete.
    -   [Stem](../packages/Stem/index.html) provides estimation of the
        parameters of a spatio-temporal model using the EM algorithm,
        estimation of the parameter standard errors using a
        spatio-temporal parametric bootstrap, spatial mapping.
    -   [spate](../packages/spate/index.html) provides spatio-temporal
        modeling of large data using a spectral SPDE approach.
    -   [pastecs](../packages/pastecs/index.html) is a package for the
        regulation, decomposition and analysis of space-time series.
    -   [STMedianPolish](../packages/STMedianPolish/index.html) analyses
        spatio-temporal data, decomposing data in n-dimensional arrays
        and using the median polish technique.
    -   R-Forge package
        [[spcopula]{.Rforge}](https://R-Forge.R-project.org/projects/spcopula/)
        provides a framework to analyze via copulas spatial and
        spatio-temporal data provided in the format of the spacetime
        package. Additionally, support for calculating different
        multivariate return periods is implemented.
    -   [solaR](../packages/solaR/index.html) is a package for computing
        solar radiation and photovoltaic systems performance.
    -   [nlme](../packages/nlme/index.html) and
        [lme4](../packages/lme4/index.html) contain functions to fit
        linear mixed models, and have facilities to model spatial and/or
        temporal effects.
-   **Point patterns**
    -   [splancs](../packages/splancs/index.html) provides methods for
        spatial and space-time point pattern analysis (khat, kernel3d,
        visualizing).
    -   [lgcp](../packages/lgcp/index.html) is a package for spatial and
        spatio-temporal modelling of point patterns using the
        log-Gaussian Cox process.
    -   [stppResid](../packages/stppResid/index.html) performs residual
        analysis on space-time point process models.
    -   [stam](../packages/stam/index.html) is an evolving package that
        target on the various methods to conduct Spatio-Temporal
        Analysis and Modelling,including Exploratory Spatio-Temporal
        Analysis and Inferred Spatio-Temporal Modelling, currently
        provides mostly kernel density estimation.
    -   [ptproc](http://www.biostat.jhsph.edu/~rpeng/software/)
        (off-CRAN) provides methods and classes for spatio-temporal
        (\"multi-dimensional\") point process.
-   **Lattice data**
    -   [surveillance](../packages/surveillance/index.html) provides
        temporal and spatio-temporal modeling and monitoring of epidemic
        phenomena.
    -   [plm](../packages/plm/index.html) fits linear panel models.
    -   [splm](../packages/splm/index.html) provides estimation and
        diagnostic testing of econometric models for spatial panel data.
    -   [sphet](../packages/sphet/index.html) fit spatial models with
        heteroskedastic innovations.
    -   [nlme](../packages/nlme/index.html) and
        [lme4](../packages/lme4/index.html) contain functions to fit
        linear mixed models, and have facilities to model spatial and/or
        temporal effects.
    -   [rsatscan](../packages/rsatscan/index.html) provides an R
        interface to the free (but non-open source) program SaTScan.
    -   [CARBayesST](../packages/CARBayesST/index.html) implements a
        class of spatio-temporal generalised linear mixed models for
        areal unit data, with inference in a Bayesian setting using
        Markov chain Monte Carlo (McMC) simulation.
    -   [gapfill](../packages/gapfill/index.html) provides tools to fill
        missing values in satellite data and to develop new gap-fill
        algorithms. The methods are tailored to data (images) observed
        at equally-spaced points in time. The package is illustrated
        with MODIS NDVI data.
-   **Moving objects, trajectories**
    -   [adehabitatLT](../packages/adehabitatLT/index.html) provides a
        collection of tools for the analysis of animal movements,
        including biased random walk simulation and home range
        estimation.
    -   [trip](../packages/trip/index.html) provides functions for
        accessing and manipulating spatial data for animal tracking.
        Filter for speed and create time spent plots from animal track
        data.
    -   [tripEstimation](../packages/tripEstimation/index.html) provides
        a Metropolis sampler and supporting functions for estimating
        animal movement from archival tags and satellite fixes. It
        further provides data handling and estimation functions for
        animal movement estimation from archival or satellite tags.
        Helper functions are included for making image summaries binned
        by time interval from MCMC simulations of point data.
    -   [diveMove](../packages/diveMove/index.html) provides utilities
        to represent, visualize, filter, analyze, and summarize
        time-depth recorder (TDR) data; miscellaneous functions for
        handling location data are also provided.
    -   [argosfilter](../packages/argosfilter/index.html) provides
        functions to filter animal satellite tracking data obtained from
        Argos. It is especially indicated for telemetry studies of
        marine animals, where Argos locations are predominantly of
        low-quality.
    -   [GeoLight](../packages/GeoLight/index.html) provides basic
        functions for global positioning based on light intensity
        measurements over time. Positioning process includes the
        determination of sun events, a discrimination of residency and
        movement periods, the calibration of period-specific data and,
        finally, the calculation of positions.
    -   [crawl](../packages/crawl/index.html): The (C)orrelated (RA)ndom
        (W)alk (L)ibrary of R functions was designed for fitting
        continuous-time correlated random walk (CTCRW) models with time
        indexed covariates. The model is fit using the Kalman-Filter on
        a state space version of the continuous-time stochastic movement
        process.
    -   [move](../packages/move/index.html) is a package for analyzing
        animal movement data; it contains functions to access movement
        data stored in [movebank](http://www.movebank.org/) as well as
        tools to visualize and statistically analyse animal movement
        data.
    -   [animalTrack](../packages/animalTrack/index.html) provides
        animal track reconstruction for high frequency 2-dimensional
        (2D) or 3-dimensional (3D) movement data. 2D and 3D animal
        tracking data can be used to reconstruct tracks through
        time/space with correction based on known positions as well as
        3D visualization of animal position and attitude.
    -   The [BBMM](../packages/BBMM/index.html) (Brownian bridge
        movement model) package provides an empirical estimate of a
        movement path using discrete location data obtained at
        relatively short time intervals. This is a continuous-time
        stochastic model of movement in which the probability of being
        in an area during the time of observation is conditioned on
        starting and ending locations. A BBMM is typically fit to animal
        location data obtained by a Global Positioning System (GPS) or
        Very High Frequency (VHF) device.
    -   The [bcpa](../packages/bcpa/index.html) package for behavioral
        change point analysis (BCPA) is a method of identifying hidden
        shifts in the underlying parameters of a time series, developed
        specifically to be applied to animal movement data which is
        irregularly sampled. The original paper on which it is based
        is: E. Gurarie, R. Andrews and K. Laidre A novel method for
        identifying behavioural changes in animal movement data (2009)
        Ecology Letters 12:5 395-408.
    -   The [smam](../packages/smam/index.html) package provides Animal
        movement models including moving-resting process with embedded
        Brownian motion, Brownian motion with measurement error.
    -   The
        [BayesianAnimalTracker](../packages/BayesianAnimalTracker/index.html)
        package provides a Bayesian melding approach to combine the GPS
        observations and Dead-Reckoned path for an accurate animal\'s
        track, or equivalently, use the GPS observations to correct the
        Dead-Reckoned path. It can take the measurement errors in the
        GPS observations into account and provide uncertainty statement
        about the corrected path. The main calculation can be done by
        the BMAnimalTrack function.
    -   The
        [TrackReconstruction](../packages/TrackReconstruction/index.html)
        package reconstructs animal tracks from magnetometer,
        accelerometer, depth and optional speed data. Designed primarily
        using data from Wildlife Computers Daily Diary tags deployed on
        northern fur seals.
    -   The [wildlifeDI](../packages/wildlifeDI/index.html) package
        provides tools for calculating a suite of indices used for
        quantifying dynamic interaction with wildlife telemetry data.
        Dynamic interaction refers to spatial-temporal associations in
        the movements of two (or more) animals.
    -   The [mkde](../packages/mkde/index.html) package provides
        functions to compute and visualize movement-based kernel density
        estimates (MKDEs) for animal utilization distributions in 2 or 3
        spatial dimensions.
    -   The [fishmove](../packages/fishmove/index.html) package provides
        functions to predict fish movement parameters based on multiple
        regression and plotting leptokurtic fish dispersal kernels (see
        Radinger and Wolter, 2013: Patterns and predictors of fish
        dispersal in rivers. Fish and Fisheries.)
    -   The [ctmcmove](../packages/ctmcmove/index.html) facilitates
        taking movement data in xyt format and pairing it with raster
        covariates within a continuous time Markov chain (CTMC)
        framework. As described in Hanks et al. (2015), this allows
        flexible modeling of movement in response to covariates (or
        covariate gradients) with model fitting possible within a
        Poisson GLM framework.
    -   The [ctmm](../packages/ctmm/index.html) provides functions for
        identifying, fitting, and applying continuous-space,
        continuous-time stochastic movement models to animal tracking
        data.
    -   The [moveHMM](../packages/moveHMM/index.html) package provides
        animal movement modelling using hidden Markov models.
        Pre-processing of tracking data, fitting HMMs to movement data,
        visualization of data and fitted model.
    -   The aim of the package [trackeR](../packages/trackeR/index.html)
        is to provide infrastructure for handling running and cycling
        data from GPS-enabled tracking devices. After extraction and
        appropriate manipulation of the training or competition
        attributes, the data are placed into session-based and
        unit-aware data objects of class trackeRdata (S3 class). The
        information in the resultant data objects can then be
        visualised, summarised, and analysed through corresponding
        flexible and extensible methods.
    -   The [VTrack](../packages/VTrack/index.html) package is designed
        to facilitate the assimilation, analysis and synthesis of animal
        location and movement data collected by the VEMCO suite of
        acoustic transmitters and receivers. As well as database and
        geographic information capabilities the principal feature of
        VTrack is the qualification and identification of ecologically
        relevant events from the acoustic detection and sensor data.
        This procedure condenses the acoustic detection database by
        orders of magnitude, greatly enhancing the synthesis of acoustic
        detection data.
    -   The
        [SimilarityMeasures](../packages/SimilarityMeasures/index.html)
        package computes four different similarity measures. The
        similarity measures included are: longest common subsequence
        (LCSS), Frechet distance, edit distance and dynamic time warping
        (DTW). Each of these similarity measures can be calculated from
        two n-dimensional trajectories, both in matrix form.
    -   The [sigloc](../packages/sigloc/index.html) package provides a
        collection of tools for estimating the location of a transmitter
        signal from radio telemetry studies using the maximum likelihood
        estimation (MLE) approach described in Lenth (1981).
    -   The [FLightR](../packages/FLightR/index.html) provides Hidden
        Markov Models for Solar Geolocation Archival Tags; it allows
        estimating positions of animal from data collected by solar
        geolocation archival tags; check the citations, which include
        [Rakhimberdiev et
        al.](http://www.movementecologyjournal.com/content/3/1/25)
    -   The [Electronic Tagging Geolocation Packages
        Repository](http://code.google.com/p/geolocation/) (off-CRAN)
        site provides a collection of statistical models (and R
        packages) to estimate position errors, movement model
        parameters, and most probable positions from tracking data.
    -   The [bsam](http://web.science.mq.edu.au/~ijonsen/code.html)
        (off-CRAN) package fits Bayesian state-space models to Argos
        satellite tracking data. Currently, models provided are DCRW
        (for location filtering), DCRWS (for location filtering and
        behavioural state estimation), and hDCRWS (a hierarchical model
        for location filtering and behavioural state estimation across
        multiple animals).
    -   The
        [[argosTrack]{.GitHub}](https://github.com/calbertsen/argosTrack/)
        (off-CRAN) package allows fitting movement models to Argos data.

**Visualization**

-   [rasterVis](../packages/rasterVis/index.html) includes a variety of
    methods that take advantage of the `z` slot of a `RasterStack` or
    `RasterBrick` object. Its
    [webpage](http://rastervis.r-forge.r-project.org/) includes several
    examples, from the hovmoller plot and horizon graph, to the density
    and histogram plots.
-   package [plotKML](../packages/plotKML/index.html) provides methods
    to convert spatio-temporal data into KML files, which can be
    displayed by external viewers, in particular Google Earth. It has a
    [gallery](http://plotkml.r-forge.r-project.org/) too.
-   package [googleVis](../packages/googleVis/index.html) provides an
    interface to show R data (tables) in the [Google Chart
    Tools](https://developers.google.com/chart/) .
    [spacetime](../packages/spacetime/index.html) has a
    [vignette](https://cran.r-project.org/web/packages/spacetime/vignettes/stgvis.html)
    demonstrating its use for spatio-temporal data.
-   Package [splancs](../packages/splancs/index.html) provides animation
    and 3D interactive plots (using [rgl](../packages/rgl/index.html))
    for displaying spatio-temporal point patterns.
-   [mvtsplot](../packages/mvtsplot/index.html) provides multivariate
    time series plots, with examples on spatio-temporal data, published
    by [Peng (2008, Journal of Statistical
    Software)](http://www.jstatsoft.org/v25/c01/) .

**Data sets**

-   Table data for fitting linear panel models are found in
    [plm](../packages/plm/index.html).
-   Package [cshapes](../packages/cshapes/index.html) contains a data
    base with country boundaries, varying over time.
-   [gstat](../packages/gstat/index.html) contains the classic Irish
    wind data.
-   [spacetime](../packages/spacetime/index.html) contains rural PM10
    air quality measurements over Germany.
-   Some parts of the Cressie and Wikle (2011) book \"Statistics for
    spatio-temporal data\" can be reproduced by `demo(CressieWikle)` in
    [spacetime](../packages/spacetime/index.html).

**Retrieving data**

Packages for retrieving data are:

-   Package [openair](../packages/openair/index.html) has tools to
    analyze, interpret and understand air pollution data, but also tools
    to download UK air quality data.
-   [ncdf4](../packages/ncdf4/index.html) and
    [RNetCDF](../packages/RNetCDF/index.html) allow reading and writing
    [netcdf](http://www.unidata.ucar.edu/software/netcdf/) files;
    [pbdNCDF4](../packages/pbdNCDF4/index.html) adds collective parallel
    read and write capability to [ncdf4](../packages/ncdf4/index.html).
-   [M3](../packages/M3/index.html) contains functions to read in and
    manipulate air quality model output from Models3-formatted files.
    This format is used by the Community Multiscale Air Quaility (CMAQ)
    model.
-   [rmatio](../packages/rmatio/index.html) is a package for reading and
    writing Matlab MAT files from R.

</div>

### CRAN packages:

-   [adehabitatLT](../packages/adehabitatLT/index.html) (core)
-   [animalTrack](../packages/animalTrack/index.html)
-   [argosfilter](../packages/argosfilter/index.html)
-   [BayesianAnimalTracker](../packages/BayesianAnimalTracker/index.html)
-   [BBMM](../packages/BBMM/index.html)
-   [bcpa](../packages/bcpa/index.html)
-   [CARBayesST](../packages/CARBayesST/index.html)
-   [crawl](../packages/crawl/index.html)
-   [cshapes](../packages/cshapes/index.html)
-   [ctmcmove](../packages/ctmcmove/index.html)
-   [ctmm](../packages/ctmm/index.html)
-   [diveMove](../packages/diveMove/index.html)
-   [fishmove](../packages/fishmove/index.html)
-   [FLightR](../packages/FLightR/index.html)
-   [gapfill](../packages/gapfill/index.html)
-   [GeoLight](../packages/GeoLight/index.html)
-   [googleVis](../packages/googleVis/index.html)
-   [gstat](../packages/gstat/index.html) (core)
-   [lgcp](../packages/lgcp/index.html)
-   [lme4](../packages/lme4/index.html)
-   [M3](../packages/M3/index.html)
-   [mkde](../packages/mkde/index.html)
-   [move](../packages/move/index.html)
-   [moveHMM](../packages/moveHMM/index.html)
-   [mvtsplot](../packages/mvtsplot/index.html)
-   [ncdf4](../packages/ncdf4/index.html)
-   [nlme](../packages/nlme/index.html)
-   [openair](../packages/openair/index.html)
-   [pastecs](../packages/pastecs/index.html)
-   [pbdNCDF4](../packages/pbdNCDF4/index.html)
-   [plm](../packages/plm/index.html)
-   [plotKML](../packages/plotKML/index.html)
-   [RandomFields](../packages/RandomFields/index.html) (core)
-   [raster](../packages/raster/index.html) (core)
-   [rasterVis](../packages/rasterVis/index.html)
-   [rgl](../packages/rgl/index.html)
-   [rmatio](../packages/rmatio/index.html)
-   [RNetCDF](../packages/RNetCDF/index.html)
-   [rsatscan](../packages/rsatscan/index.html)
-   [sf](../packages/sf/index.html)
-   [sigloc](../packages/sigloc/index.html)
-   [SimilarityMeasures](../packages/SimilarityMeasures/index.html)
-   [smam](../packages/smam/index.html)
-   [solaR](../packages/solaR/index.html)
-   [sp](../packages/sp/index.html) (core)
-   [spacetime](../packages/spacetime/index.html) (core)
-   [spate](../packages/spate/index.html)
-   [SpatioTemporal](../packages/SpatioTemporal/index.html)
-   [spatstat](../packages/spatstat/index.html)
-   [spBayes](../packages/spBayes/index.html)
-   [sphet](../packages/sphet/index.html)
-   [splancs](../packages/splancs/index.html)
-   [splm](../packages/splm/index.html)
-   [spTimer](../packages/spTimer/index.html)
-   [stam](../packages/stam/index.html)
-   [Stem](../packages/Stem/index.html)
-   [STMedianPolish](../packages/STMedianPolish/index.html)
-   [stppResid](../packages/stppResid/index.html)
-   [surveillance](../packages/surveillance/index.html) (core)
-   [trackeR](../packages/trackeR/index.html)
-   [TrackReconstruction](../packages/TrackReconstruction/index.html)
-   [trip](../packages/trip/index.html) (core)
-   [tripEstimation](../packages/tripEstimation/index.html)
-   [VTrack](../packages/VTrack/index.html)
-   [wildlifeDI](../packages/wildlifeDI/index.html)
-   [xts](../packages/xts/index.html) (core)

### Related links:

-   CRAN Task View: [Spatial](Spatial.html)
-   CRAN Task View: [TimeSeries](TimeSeries.html)
-   CRAN Task View: [Econometrics](Econometrics.html)
-   CRAN Task View: [Environmetrics](Environmetrics.html)
-   CRAN Task View: [DifferentialEquations](DifferentialEquations.html)
-   R-Forge Project:
    [[spcopula]{.Rforge}](https://R-Forge.R-project.org/projects/spcopula/)
-   [Pebesma (2012). spacetime: Spatio-Temporal Data in R. Journal of
    Statistical Software, 51(7).](http://www.jstatsoft.org/v51/i07/)
-   [Peng (2008). A Method for Visualizing Multivariate Time Series
    Data. Journal of Statistical Software, Code Snippets,
    25(1).](http://www.jstatsoft.org/v25/c01/)
-   [Gapminder world: A tool for visualization of economic and health
    indicators per country (including access to the underlying tables in
    time-wide form).](http://www.gapminder.org/world/)
