## CRAN Task View: Analysis of Spatial Data

  ----------------- -------------------------------------------
  **Maintainer:**   Roger Bivand
  **Contact:**      Roger.Bivand at nhh.no
  **Version:**      2018-04-11
  **URL:**          <https://CRAN.R-project.org/view=Spatial>
  ----------------- -------------------------------------------

<div>

Base R includes many functions that can be used for reading,
visualising, and analysing spatial data. The focus in this view is on
\"geographical\" spatial data, where observations can be identified with
geographical locations, and where additional information about these
locations may be retrieved if the location is recorded with care. Base R
functions are complemented by contributed packages, some of which are on
CRAN, and others are still in development. One active location is
[R-Forge](http://R-Forge.R-project.org/) , which lists \"Spatial Data
and Statistics\" projects in its [project
tree](http://R-Forge.R-project.org/softwaremap/trove_list.php) .
Information on R-spatial packages, especially
[sp](../packages/sp/index.html) is posted on the R-Forge rspatial
project [website](http://rspatial.R-Forge.R-project.org/) , including a
visualisation gallery. Active development of
[sp](../packages/sp/index.html) is continuing on
[[sp]{.GitHub}](https://github.com/edzer/sp/).

The contributed packages address two broad areas: moving spatial data
into and out of R, and analysing spatial data in R.

The [R-SIG-Geo](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
mailing-list is a good place to begin for obtaining help and discussing
questions about both accessing data, and analysing it. The mailing list
is a good place to search for information about relevant courses.
Further information about courses may be found under the \"Events\" tab
of [this blog](http://r-spatial.org/) .

There are a number of contributed tutorials and introductions; a recent
one is [Introduction to visualising spatial data in
R](https://cran.r-project.org/doc/contrib/intro-spatial-rl.pdf) by Robin
Lovelace and James Cheshire.

The packages in this view can be roughly structured into the following
topics. If you think that some package is missing from the list, please
let me know.

-   **Classes for spatial data** : Because many of the packages
    importing and using spatial data have had to include objects of
    storing data and functions for visualising it, an initiative is in
    progress to construct shared classes and plotting functions for
    spatial data. The [sp](../packages/sp/index.html) package is
    discussed in a note in [R
    News](http://CRAN.R-project.org/doc/Rnews/Rnews_2005-2.pdf) . A new
    package called [sf](../packages/sf/index.html) is now on CRAN, and
    is being actively developed on
    [[sfr]{.GitHub}](https://github.com/edzer/sfr/), providing Simple
    Features for R. The development of the package is being supported by
    the R Consortium. It provides simple features access for vector
    data, and as such is a modern implementation of parts of
    [sp](../packages/sp/index.html). Many other packages have become
    dependent on the [sp](../packages/sp/index.html) classes, including
    [rgdal](../packages/rgdal/index.html) and
    [maptools](../packages/maptools/index.html). The
    [rgeos](../packages/rgeos/index.html) package provides an interface
    to topology functions for [sp](../packages/sp/index.html) objects
    using [GEOS](http://trac.osgeo.org/geos/) . The
    [stplanr](../packages/stplanr/index.html) provides a
    \"SpatialLinesNetwork\" class based on objects defined in
    [sp](../packages/sp/index.html) and
    [igraph](../packages/igraph/index.html) that can be used for routing
    analysis within R. Another network package is
    [shp2graph](../packages/shp2graph/index.html). The
    [cleangeo](../packages/cleangeo/index.html) may be used to inspect
    spatial objects, facilitate handling and reporting of topology
    errors and geometry validity issues. It claims to provide a geometry
    cleaner that will fix all geometry problems, and eliminate (at least
    reduce) the likelihood of having issues when doing spatial data
    processing. The [raster](../packages/raster/index.html) package is a
    major extension of spatial data classes to virtualise access to
    large rasters, permitting large objects to be analysed, and
    extending the analytical tools available for both raster and vector
    data. Used with [rasterVis](../packages/rasterVis/index.html), it
    can also provide enhanced visualisation and interaction. The
    [spatial.tools](../packages/spatial.tools/index.html) package
    contains spatial functions meant to enhance the core functionality
    of the [raster](../packages/raster/index.html) package, including a
    parallel processing engine for use with rasters. The
    [micromap](../packages/micromap/index.html) package provides linked
    micromaps using ggplot2. The [recmap](../packages/recmap/index.html)
    package provides rectangular cartograms with rectangle sizes
    reflecting for example population; the
    [statebins](../packages/statebins/index.html) provides a simpler
    binning approach to US states. The
    [spacetime](../packages/spacetime/index.html) package extends the
    shared classes defined in [sp](../packages/sp/index.html) for
    spatio-temporal data (see [Spatio-Temporal Data in
    R](http://www.jstatsoft.org/v51/i07) ). The
    [Grid2Polygons](../packages/Grid2Polygons/index.html) converts a
    spatial object from class SpatialGridDataFrame to
    SpatialPolygonsDataFrame.

    An alternative approach to some of these issues is implemented in
    the [PBSmapping](../packages/PBSmapping/index.html) package;
    [PBSmodelling](../packages/PBSmodelling/index.html) provides
    modelling support. In addition,
    [GEOmap](../packages/GEOmap/index.html) provides mapping facilities
    directed to meet the needs of geologists, and uses the
    [geomapdata](../packages/geomapdata/index.html) package.

-   **Handling spatial data** : A number of packages have been written
    using sp classes. The [raster](../packages/raster/index.html)
    package introduces many GIS methods that now permit much to be done
    with spatial data without having to use GIS in addition to R. It may
    be complemented by [gdistance](../packages/gdistance/index.html),
    which provided calculation of distances and routes on geographic
    grids. [geosphere](../packages/geosphere/index.html) permits
    computations of distance and area to be carried out on spatial data
    in geographical coordinates. The
    [spsurvey](../packages/spsurvey/index.html) package provides a range
    of sampling functions. The [trip](../packages/trip/index.html)
    package extends sp classes to permit the accessing and manipulating
    of spatial data for animal tracking. The
    [hdeco](../packages/hdeco/index.html) package provides hierarchical
    decomposition of entropy for categorical map comparisons.
    [spcosa](../packages/spcosa/index.html) provides spatial coverage
    sampling and random sampling from compact geographical strata. The
    [magclass](../packages/magclass/index.html) offers a data class for
    increased interoperability working with spatial-temporal data
    together with corresponding functions and methods (conversions,
    basic calculations and basic data manipulation). The class
    distinguishes between spatial, temporal and other dimensions to
    facilitate the development and interoperability of tools build for
    it. Additional features are name-based addressing of data and
    internal consistency checks (e.g. checking for the right data order
    in calculations).

    The UScensus2000 suite of packages
    ([UScensus2000cdp](../packages/UScensus2000cdp/index.html),
    [UScensus2000tract](../packages/UScensus2000tract/index.html)) makes
    the use of data from the 2000 US Census more convenient. An
    important data set, Guerry\'s \"Moral Statistics of France\", has
    been made available in the [Guerry](../packages/Guerry/index.html)
    package, which provides data and maps and examples designed to
    contribute to the integration of multivariate and spatial analysis.
    The [marmap](../packages/marmap/index.html) package is designed for
    downloading, plotting and manipulating bathymetric and topographic
    data in R. [marmap](../packages/marmap/index.html) can query the
    ETOPO1 bathymetry and topography database hosted by the NOAA, use
    simple latitude-longitude-depth data in ascii format, and take
    advantage of the advanced plotting tools available in R to build
    publication-quality bathymetric maps (see the
    [PLOS](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0073051)
    paper). Modern country boundaries are provided at 2 resolutions by
    [rworldmap](../packages/rworldmap/index.html) along with functions
    to join and map tabular data referenced by country names or codes.
    Chloropleth and bubble maps are supported and general functions to
    work on user supplied maps (see [A New R package for Mapping Global
    Data](http://journal.r-project.org/archive/2011-1/RJournal_2011-1_South.pdf)
    . Higher resolution country borders are available from the linked
    package [rworldxtra](../packages/rworldxtra/index.html). Historical
    country boundaries (1946-2012) can be obtained from the
    [cshapes](../packages/cshapes/index.html) package along with
    functions for calculating distance matrices (see [Mapping and
    Measuring Country
    Shapes](http://journal.r-project.org/archive/2010-1/RJournal_2010-1_Weidmann+Skrede~Gleditsch.pdf)
    ).

    The [landsat](../packages/landsat/index.html) package with
    accompanying [JSS paper](http://www.jstatsoft.org/v43/i04) provides
    tools for exploring and developing correction tools for remote
    sensing data. [taRifx](../packages/taRifx/index.html) is a
    collection of utility and convenience functions, and some
    interesting spatial functions. The
    [gdalUtils](../packages/gdalUtils/index.html) package provides
    wrappers for the Geospatial Data Abstraction Library (GDAL)
    Utilities.

    An rOpenSci [blog
    entry](http://ropensci.org/blog/blog/2016/11/22/geospatial-suite)
    described a GeoJSON-centred approach to reading GeoJSON and WKT
    data. GeoJSON can be written and read using
    [rgdal](../packages/rgdal/index.html), and WKT by
    [rgeos](../packages/rgeos/index.html). The entry lists
    [geojson](../packages/geojson/index.html),
    [geojsonio](../packages/geojsonio/index.html),
    [geoaxe](../packages/geoaxe/index.html) and
    [lawn](../packages/lawn/index.html) among others. The
    [rgbif](../packages/rgbif/index.html) package is used to access
    Global Biodiversity Information Facility (GBIF) data). The
    [geoaxe](../packages/geoaxe/index.html) allows users to split
    \'geospatial\' objects into pieces. The
    [lawn](../packages/lawn/index.html) package is a client for
    \'Turfjs\' for \'geospatial\' analysis.

-   **Reading and writing spatial data -
    [rgdal](../packages/rgdal/index.html)** : Maps may be vector-based
    or raster-based. The [rgdal](../packages/rgdal/index.html) package
    provides bindings to [GDAL](http://www.gdal.org/) -supported raster
    formats and [OGR](http://www.gdal.org/ogr/) -supported vector
    formats. It contains functions to write raster files in supported
    formats. The package also provides
    [PROJ.4](http://trac.osgeo.org/proj/) projection support for vector
    objects ( [this site](http://spatialreference.org) provides
    searchable online PROJ.4 representations of projections). Affine and
    similarity transformations on sp objects may be made using functions
    in the [vec2dtransf](../packages/vec2dtransf/index.html) package.
    The Windows and Mac OSX CRAN binaries of
    [rgdal](../packages/rgdal/index.html) include subsets of possible
    data source drivers; if others are needed, use other conversion
    utilities, or install from source against a version of GDAL with the
    required drivers. The [rgeos](../packages/rgeos/index.html) package
    provides functions for reading and writing well-known text (WKT)
    geometry, and the [wkb](../packages/wkb/index.html) package provides
    functions for reading and writing well-known binary (WKB) geometry.

-   **Reading and writing spatial data - other packages** : There are a
    number of other packages for accessing vector data on CRAN:
    [maps](../packages/maps/index.html) (with
    [mapdata](../packages/mapdata/index.html) and
    [mapproj](../packages/mapproj/index.html)) provides access to the
    same kinds of geographical databases as S -
    [RArcInfo](../packages/RArcInfo/index.html) allows ArcInfo v.7
    binary files and \*.e00 files to be read, and
    [maptools](../packages/maptools/index.html) and
    [shapefiles](../packages/shapefiles/index.html) read and write
    ArcGIS/ArcView shapefiles; for NetCDF files,
    [ncdf4](../packages/ncdf4/index.html) or
    [RNetCDF](../packages/RNetCDF/index.html) may be used. The
    [maptools](../packages/maptools/index.html) package also provides
    helper functions for writing map polygon files to be read by
    WinBUGS, Mondrian, and the tmap command in Stata. It also provides
    interface functions between
    [PBSmapping](../packages/PBSmapping/index.html) and
    [spatstat](../packages/spatstat/index.html) and sp classes, in
    addition to [maps](../packages/maps/index.html) databases and sp
    classes. There is also an interface to GSHHS shoreline databases.
    The [gmt](../packages/gmt/index.html) package gives a simple
    interface between GMT map-making software and R.
    [geonames](../packages/geonames/index.html) is an interface to the
    [www.geonames.org](http://www.geonames.org/) service.
    [OpenStreetMap](../packages/OpenStreetMap/index.html) gives access
    to open street map raster images, and
    [osmar](../packages/osmar/index.html) provides infrastructure to
    access OpenStreetMap data from different sources, to work with the
    data in common R manner, and to convert data into available
    infrastructure provided by existing R packages.

    The [rpostgis](../packages/rpostgis/index.html) package provides
    additional functions to the \'RPostgreSQL\' package to interface R
    with a \'PostGIS\'-enabled database, as well as convenient wrappers
    to common \'PostgreSQL\' queries. The
    [postGIStools](../packages/postGIStools/index.html) package provides
    functions to convert geometry and \'hstore\' data types from
    \'PostgreSQL\' into standard R objects, as well as to simplify the
    import of R data frames (including spatial data frames) into
    \'PostgreSQL\'

    Integration with version 6.\* and of the leading open source GIS,
    GRASS, is provided in CRAN package
    [spgrass6](../packages/spgrass6/index.html), using
    [rgdal](../packages/rgdal/index.html) for exchanging data. For GRASS
    7.\*, use [rgrass7](../packages/rgrass7/index.html).
    [RPyGeo](../packages/RPyGeo/index.html) is a wrapper for Python
    access to the ArcGIS GeoProcessor, and
    [RSAGA](../packages/RSAGA/index.html) is a similar shell-based
    wrapper for SAGA commands. The [RQGIS](../packages/RQGIS/index.html)
    package establishes an interface between R and QGIS, i.e. it allows
    the user to access QGIS functionalities from the R console. It
    achieves this by using the QGIS Python API via the command line.
    Note also [this
    thread](https://stat.ethz.ch/pipermail/r-sig-geo/2016-August/024817.html)
    on an alternative R/QGIS integration.

-   **Visualisation** : For visualization, the colour palettes provided
    in the [RColorBrewer](../packages/RColorBrewer/index.html) package
    are very useful, and may be modified or extended using the
    `colorRampPalette` function provided with R. The
    [classInt](../packages/classInt/index.html) package provides
    functions for choosing class intervals for thematic cartography. The
    [tmap](../packages/tmap/index.html) package provides a modern basis
    for thematic mapping optionally using a Grammar of Graphics syntax.
    Because it has a custom grid graphics platform, it obviates the need
    to fortify geometries to use with ggplot2. The
    [mapview](../packages/mapview/index.html) package provides methods
    to view spatial objects interactively, usually on a web mapping
    base. The [quickmapr](../packages/quickmapr/index.html) package
    provides a simple method to visualize \'sp\' and \'raster\' objects,
    allows for basic zooming, panning, identifying, and labeling of
    spatial objects, and does not require that the data be in geographic
    coordinates. The [cartography](../packages/cartography/index.html)
    package allows various cartographic representations such as
    proportional symbols, choropleth, typology, flows or
    discontinuities. The [mapmisc](../packages/mapmisc/index.html)
    package is a minimal, light-weight set of tools for producing nice
    looking maps in R, with support for map projections.If the user
    wishes to place a map backdrop behind other displays, the the
    [RgoogleMaps](../packages/RgoogleMaps/index.html) package for
    accessing Google Maps(TM) may be useful.
    [ggmap](../packages/ggmap/index.html) may be used for spatial
    visualisation with Google Maps and OpenStreetMap;
    [ggsn](../packages/ggsn/index.html) provides North arrows and scales
    for such maps. The
    [plotGoogleMaps](../packages/plotGoogleMaps/index.html) package
    provides methods for the visualisation of spatial and
    spatio-temporal objects in Google Maps in a web browser.
    [plotKML](../packages/plotKML/index.html) is a package providing
    methods for the visualisation of spatial and spatio-temporal objects
    in Google Earth. A further option is
    [leafletR](../packages/leafletR/index.html), which provides basic
    web-mapping functionality to combine vector data files and online
    map tiles from different sources.

-   **Point pattern analysis** : The
    [spatial](../packages/spatial/index.html) package is a recommended
    package shipped with base R, and contains several core functions,
    including an implementation of Khat by its author, Prof. Ripley. In
    addition, [spatstat](../packages/spatstat/index.html) allows freedom
    in defining the region(s) of interest, and makes extensions to
    marked processes and spatial covariates. Its strengths are
    model-fitting and simulation, and it has a useful
    [homepage](http://www.spatstat.org/) . It is the only package that
    will enable the user to fit inhomogeneous point process models with
    interpoint interactions. The
    [spatgraphs](../packages/spatgraphs/index.html) package provides
    graphs, graph visualisation and graph based summaries to be used
    with spatial point pattern analysis. The
    [splancs](../packages/splancs/index.html) package also allows point
    data to be analysed within a polygonal region of interest, and
    covers many methods, including 2D kernel densities. The
    [smacpod](../packages/smacpod/index.html) package provides various
    statistical methods for analyzing case-control point data. The
    methods available closely follow those in chapter 6 of Applied
    Spatial Statistics for Public Health Data by Waller and Gotway
    (2004).

    [ecespa](../packages/ecespa/index.html) provides wrappers, functions
    and data for spatial point pattern analysis, used in the book on
    Spatial Ecology of the ECESPA/AEET. The functions for binning points
    on grids in [ash](../packages/ash/index.html) may also be of
    interest. The [ads](../packages/ads/index.html) package perform
    first- and second-order multi-scale analyses derived from Ripley\'s
    K-function. The [aspace](../packages/aspace/index.html) package is a
    collection of functions for estimating centrographic statistics and
    computational geometries from spatial point patterns.
    [DSpat](../packages/DSpat/index.html) contains functions for spatial
    modelling for distance sampling data and
    [spatialsegregation](../packages/spatialsegregation/index.html)
    provides segregation measures for multitype spatial point patterns.
    [GriegSmith](../packages/GriegSmith/index.html) uses the Grieg-Smith
    method on 2 dimensional spatial data. The
    [dbmss](../packages/dbmss/index.html) package allows simple
    computation of a full set of spatial statistic functions of
    distance, including classical ones (Ripley\'s K and others) and more
    recent ones used by spatial economists (Duranton and Overman\'s Kd,
    Marcon and Puech\'s M). It relies on spatstat for core calculation.
    [latticeDensity](../packages/latticeDensity/index.html) contains
    functions that compute the lattice-based density estimator of Barry
    and McIntyre, which accounts for point processes in two-dimensional
    regions with irregular boundaries and holes.

-   **Geostatistics** : The [gstat](../packages/gstat/index.html)
    package provides a wide range of functions for univariate and
    multivariate geostatistics, also for larger datasets, while
    [geoR](../packages/geoR/index.html) and
    [geoRglm](../packages/geoRglm/index.html) contain functions for
    model-based geostatistics. Variogram diagnostics may be carried out
    with [vardiag](../packages/vardiag/index.html). Automated
    interpolation using [gstat](../packages/gstat/index.html) is
    available in [automap](../packages/automap/index.html). This family
    of packages is supplemented by
    [intamap](../packages/intamap/index.html) with procedures for
    automated interpolation. A similar wide range of functions is to be
    found in the [fields](../packages/fields/index.html) package. The
    [spatial](../packages/spatial/index.html) package is shipped with
    base R, and contains several core functions. The
    [spBayes](../packages/spBayes/index.html) package fits Gaussian
    univariate and multivariate models with MCMC.
    [ramps](../packages/ramps/index.html) is a different Bayesian
    geostatistical modelling package. The
    [geospt](../packages/geospt/index.html) package contains some
    geostatistical and radial basis functions, including prediction and
    cross validation. Besides, it includes functions for the design of
    optimal spatial sampling networks based on geostatistical modelling.
    [spsann](../packages/spsann/index.html) is another package to offer
    functions to optimize sample configurations, using spatial simulated
    annealing. The [geostatsp](../packages/geostatsp/index.html) package
    offers geostatistical modelling facilities using Raster and
    SpatialPoints objects are provided. Non-Gaussian models are fit
    using INLA, and Gaussian geostatistical models use Maximum
    Likelihood Estimation. The [FRK](../packages/FRK/index.html) package
    is a tool for spatial/spatio-temporal modelling and prediction with
    large datasets. The approach, discussed in Cressie and Johannesson
    (2008), decomposes the field, and hence the covariance function,
    using a fixed set of n basis functions, where n is typically much
    smaller than the number of data points (or polygons) m.

    The [RandomFields](../packages/RandomFields/index.html) package
    provides functions for the simulation and analysis of random fields,
    and variogram model descriptions can be passed between
    [geoR](../packages/geoR/index.html),
    [gstat](../packages/gstat/index.html) and this package.
    [SpatialExtremes](../packages/SpatialExtremes/index.html) proposes
    several approaches for spatial extremes modelling using
    [RandomFields](../packages/RandomFields/index.html). In addition,
    [CompRandFld](../packages/CompRandFld/index.html),
    [constrainedKriging](../packages/constrainedKriging/index.html) and
    [geospt](../packages/geospt/index.html) provide alternative
    approaches to geostatistical modelling. The
    [spTimer](../packages/spTimer/index.html) package is able to fit,
    spatially predict and temporally forecast large amounts of
    space-time data using \[1\] Bayesian Gaussian Process (GP) Models,
    \[2\] Bayesian Auto-Regressive (AR) Models, and \[3\] Bayesian
    Gaussian Predictive Processes (GPP) based AR Models. The
    [rtop](../packages/rtop/index.html) package provides functions for
    the geostatistical interpolation of data with irregular spatial
    support such as runoff related data or data from administrative
    units. The [georob](../packages/georob/index.html) package provides
    functions for fitting linear models with spatially correlated errors
    by robust and Gaussian Restricted Maximum Likelihood and for
    computing robust and customary point and block kriging predictions,
    along with utility functions for cross-validation and for unbiased
    back-transformation of kriging predictions of log-transformed data.
    The [SpatialTools](../packages/SpatialTools/index.html) package has
    an emphasis on kriging, and provides functions for prediction and
    simulation. It is extended by
    [ExceedanceTools](../packages/ExceedanceTools/index.html), which
    provides tools for constructing confidence regions for exceedance
    regions and contour lines. The [gear](../packages/gear/index.html)
    package implements common geostatistical methods in a clean,
    straightforward, efficient manner, and is said to be a quasi reboot
    of [SpatialTools](../packages/SpatialTools/index.html). The
    [sperrorest](../packages/sperrorest/index.html) package implements
    spatial error estimation and permutation-based spatial variable
    importance using different spatial cross-validation and spatial
    block bootstrap methods. The [spm](../packages/spm/index.html)
    package provides functions for hybrid geostatistical and machine
    learning methods for spatial predictive modelling. It currently
    contains two commonly used geostatistical methods, two machine
    learning methods, four hybrid methods and two averaging methods.

    The [sgeostat](../packages/sgeostat/index.html) package is also
    available. Within the same general topical area are the
    [deldir](../packages/deldir/index.html) and
    [tripack](../packages/tripack/index.html) packages for triangulation
    and the [akima](../packages/akima/index.html) package for spline
    interpolation; the [MBA](../packages/MBA/index.html) package
    provides scattered data interpolation with multilevel B-splines. In
    addition, there are the
    [spatialCovariance](../packages/spatialCovariance/index.html)
    package, which supports the computation of spatial covariance
    matrices for data on rectangles, the
    [regress](../packages/regress/index.html) package building in part
    on [spatialCovariance](../packages/spatialCovariance/index.html),
    and the [tgp](../packages/tgp/index.html) package. The
    [Stem](../packages/Stem/index.html) package provides for the
    estimation of the parameters of a spatio-temporal model using the EM
    algorithm, and the estimation of the parameter standard errors using
    a spatio-temporal parametric bootstrap.
    [FieldSim](../packages/FieldSim/index.html) is another random fields
    simulations package. The [SSN](../packages/SSN/index.html) is for
    geostatistical modeling for data on stream networks, including
    models based on in-stream distance. Models are created using moving
    average constructions. Spatial linear models, including covariates,
    can be fit with ML or REML. Mapping and other graphical functions
    are included. The [ipdw](../packages/ipdw/index.html) provides
    functions o interpolate georeferenced point data via Inverse Path
    Distance Weighting. Useful for coastal marine applications where
    barriers in the landscape preclude interpolation with Euclidean
    distances. [RSurvey](../packages/RSurvey/index.html) may be used as
    a processing program for spatially distributed data, and is capable
    of error corrections and data visualisation.

-   **Disease mapping and areal data analysis** :
    [DCluster](../packages/DCluster/index.html) is a package for the
    detection of spatial clusters of diseases. It extends and depends on
    the [spdep](../packages/spdep/index.html) package, which provides
    basic functions for building neighbour lists and spatial weights,
    tests for spatial autocorrelation for areal data like Moran\'s I,
    and functions for fitting spatial regression models, such as SAR and
    CAR models. These models assume that the spatial dependence can be
    described by known weights. In
    [spdep](../packages/spdep/index.html), the `ME` and
    `SpatialFiltering` functions provide Moran Eigenvector model
    fitting, as do more modern functions in the
    [spmoran](../packages/spmoran/index.html) package. The
    [SpatialEpi](../packages/SpatialEpi/index.html) package provides
    implementations of cluster detection and disease mapping functions,
    including Bayesian cluster detection, and supports strata. The
    [smerc](../packages/smerc/index.html) package provides statistical
    methods for the analysis of data areal data, with a focus on cluster
    detection. The
    [diseasemapping](../packages/diseasemapping/index.html) package
    offers the formatting of population and case data, calculation of
    Standardized Incidence Ratios, and fitting the BYM model using INLA.
    Regionalisation of polygon objects is provided by
    [AMOEBA](../packages/AMOEBA/index.html): a function to calculate
    spatial clusters using the Getis-Ord local statistic. It searches
    irregular clusters (ecotopes) on a map, and by `skater` in
    [spdep](../packages/spdep/index.html). The
    [seg](../packages/seg/index.html) and
    [OasisR](../packages/OasisR/index.html) packages provide functions
    for measuring spatial segregation;
    [OasisR](../packages/OasisR/index.html) includes Monte Carlo
    simulations to test the indices. The
    [spgwr](../packages/spgwr/index.html) package contains an
    implementation of geographically weighted regression methods for
    exploring possible non-stationarity. The
    [gwrr](../packages/gwrr/index.html) package fits geographically
    weighted regression (GWR) models and has tools to diagnose and
    remediate collinearity in the GWR models. Also fits geographically
    weighted ridge regression (GWRR) and geographically weighted lasso
    (GWL) models. The [GWmodel](../packages/GWmodel/index.html) package
    contains functions for computing geographically weighted models. The
    [lctools](../packages/lctools/index.html) package provides
    researchers and educators with easy-to-learn user friendly tools for
    calculating key spatial statistics and to apply simple as well as
    advanced methods of spatial analysis in real data. These include:
    Local Pearson and Geographically Weighted Pearson Correlation
    Coefficients, Spatial Inequality Measures (Gini, Spatial Gini, LQ,
    Focal LQ), Spatial Autocorrelation (Global and Local Moran\'s I),
    several Geographically Weighted Regression techniques and other
    Spatial Analysis tools (other geographically weighted statistics).
    This package also contains functions for measuring the significance
    of each statistic calculated, mainly based on Monte Carlo
    simulations. The [sparr](../packages/sparr/index.html) package
    provides another approach to relative risks. The
    [CARBayes](../packages/CARBayes/index.html) package implements
    Bayesian hierarchical spatial areal unit models. In such models the
    spatial correlation is modelled by a set of random effects, which
    are assigned a conditional autoregressive (CAR) prior distribution.
    Examples of the models included are the BYM model as well as a
    recently developed localised spatial smoothing model. The
    [glmmBUGS](../packages/glmmBUGS/index.html) package is a helpful way
    of passing out spatial models to WinBUGS. The
    [spaMM](../packages/spaMM/index.html) package fits spatial GLMMs,
    using the Matern correlation function as the basic model for spatial
    random effects. The [PReMiuM](../packages/PReMiuM/index.html)
    package is for profile regression, which is a Dirichlet process
    Bayesian clustering model; it provides a spatial CAR term that can
    be included in the fixed effects (which are global, ie. non cluster
    specific, parameters) to account for any spatial correlation in the
    residuals. The [spacom](../packages/spacom/index.html) package
    provides tools to construct and exploit spatially weighted context
    data, and further allows combining the resulting spatially weighted
    context data with individual-level predictor and outcome variables,
    for the purposes of multilevel modelling. The
    [geospacom](../packages/geospacom/index.html) package generates
    distance matrices from shape files and represents spatially weighted
    multilevel analysis results. Spatial survival analysis is provided
    by the [spatsurv](../packages/spatsurv/index.html) - Bayesian
    inference for parametric proportional hazards spatial survival
    models - and [spBayesSurv](../packages/spBayesSurv/index.html) -
    Bayesian Modeling and Analysis of Spatially Correlated Survival
    Data - packages. The [spselect](../packages/spselect/index.html)
    package provides modelling functions based on forward stepwise
    regression, incremental forward stagewise regression, least angle
    regression (LARS), and lasso models for selecting the spatial scale
    of covariates in regression models.

-   **Spatial regression** : The choice of function for spatial
    regression will depend on the support available. If the data are
    characterised by point support and the spatial process is
    continuous, geostatistical methods may be used, or functions in the
    [nlme](../packages/nlme/index.html) package. If the support is
    areal, and the spatial process is not being treated as continuous,
    functions provided in the [spdep](../packages/spdep/index.html)
    package may be used. This package can also be seen as providing
    spatial econometrics functions, and, as noted above, provides basic
    functions for building neighbour lists and spatial weights, tests
    for spatial autocorrelation for areal data like Moran\'s I, and
    functions for fitting spatial regression models. It provides the
    full range of local indicators of spatial association, such as local
    Moran\'s I and diagnostic tools for fitted linear models, including
    Lagrange Multiplier tests. Spatial regression models that can be
    fitted using maximum likelihood include spatial lag models, spatial
    error models, and spatial Durbin models. For larger data sets,
    sparse matrix techniques can be used for maximum likelihood fits,
    while spatial two stage least squares and generalised method of
    moments estimators are an alternative. When using GMM,
    [sphet](../packages/sphet/index.html) can be used to accommodate
    both autocorrelation and heteroskedasticity. The
    [McSpatial](../packages/McSpatial/index.html) provides functions for
    locally weighted regression, semiparametric and conditionally
    parametric regression, fourier and cubic spline functions, GMM and
    linearized spatial logit and probit, k-density functions and
    counterfactuals, nonparametric quantile regression and conditional
    density functions, Machado-Mata decomposition for quantile
    regressions, spatial AR model, repeat sales models, and
    conditionally parametric logit and probit. The
    [splm](../packages/splm/index.html) package provides methods for
    fitting spatial panel data by maximum likelihood and GM. The two
    small packages [S2sls](../packages/S2sls/index.html) and
    [spanel](../packages/spanel/index.html) provide alternative
    implementations without most of the facilities of
    [splm](../packages/splm/index.html). The
    [HSAR](../packages/HSAR/index.html) package provides Hierarchical
    Spatial Autoregressive Models (HSAR), based on a Bayesian Markov
    Chain Monte Carlo (MCMC) algorithm.
    [spatialprobit](../packages/spatialprobit/index.html) make possible
    Bayesian estimation of the spatial autoregressive probit model (SAR
    probit model). The
    [ProbitSpatial](../packages/ProbitSpatial/index.html) package
    provides methods for fitting Binomial spatial probit models to
    larger data sets; spatial autoregressive (SAR) and spatial error
    (SEM) probit models are included. The
    [starma](../packages/starma/index.html) package provides functions
    to identify, estimate and diagnose a Space-Time AutoRegressive
    Moving Average (STARMA) model.

-   **Ecological analysis** : There are many packages for analysing
    ecological and environmental data. They include
    [ade4](../packages/ade4/index.html) for exploratory and Euclidean
    methods in the environmental sciences, the adehabitat family of
    packages for the analysis of habitat selection by animals
    ([adehabitatHR](../packages/adehabitatHR/index.html),
    [adehabitatHS](../packages/adehabitatHS/index.html),
    [adehabitatLT](../packages/adehabitatLT/index.html), and
    [adehabitatMA](../packages/adehabitatMA/index.html)),
    [pastecs](../packages/pastecs/index.html) for the regulation,
    decomposition and analysis of space-time series,
    [vegan](../packages/vegan/index.html) for ordination methods and
    other useful functions for community and vegetation ecologists, and
    many other functions in other contributed packages. One such is
    [tripEstimation](../packages/tripEstimation/index.html), basing on
    the classes provided by [trip](../packages/trip/index.html).
    [ncf](../packages/ncf/index.html) has entered CRAN recently, and
    provides a range of spatial nonparametric covariance functions. The
    [spind](../packages/spind/index.html) package provides functions for
    spatial methods based on generalized estimating equations (GEE) and
    wavelet-revised methods (WRM), functions for scaling by wavelet
    multiresolution regression (WMRR), conducting multi-model inference,
    and stepwise model selection.
    [rangeMapper](../packages/rangeMapper/index.html) is a package to
    manipulate species range (extent-of-occurrence) maps, mainly tools
    for easy generation of biodiversity (species richness) or
    life-history traits maps. The
    [siplab](../packages/siplab/index.html) package is a platform for
    experimenting with spatially explicit individual-based vegetation
    models. [ModelMap](../packages/ModelMap/index.html) builds on other
    packages to create models using underlying GIS data. The
    [SpatialPosition](../packages/SpatialPosition/index.html) computes
    spatial position models: Stewart potentials, Reilly catchment areas,
    Huff catchment areas. The
    [Watersheds](../packages/Watersheds/index.html) package provides
    methods for watersheds aggregation and spatial drainage network
    analysis. An off-CRAN package -
    [Rcitrus](http://www.leg.ufpr.br/Rcitrus/) - is for the spatial
    analysis of plant disease incidence. The
    [ngspatial](../packages/ngspatial/index.html) package provides tools
    for analyzing spatial data, especially non-Gaussian areal data. It
    supports the sparse spatial generalized linear mixed model of Hughes
    and Haran (2013) and the centered autologistic model of Caragea and
    Kaiser (2009). The [Environmetrics](Environmetrics.html) Task View
    contains a much more complete survey of relevant functions and
    packages.

</div>

### CRAN packages:

-   [ade4](../packages/ade4/index.html)
-   [adehabitatHR](../packages/adehabitatHR/index.html)
-   [adehabitatHS](../packages/adehabitatHS/index.html)
-   [adehabitatLT](../packages/adehabitatLT/index.html)
-   [adehabitatMA](../packages/adehabitatMA/index.html)
-   [ads](../packages/ads/index.html)
-   [akima](../packages/akima/index.html)
-   [AMOEBA](../packages/AMOEBA/index.html)
-   [ash](../packages/ash/index.html)
-   [aspace](../packages/aspace/index.html)
-   [automap](../packages/automap/index.html)
-   [CARBayes](../packages/CARBayes/index.html)
-   [cartography](../packages/cartography/index.html)
-   [classInt](../packages/classInt/index.html) (core)
-   [cleangeo](../packages/cleangeo/index.html)
-   [CompRandFld](../packages/CompRandFld/index.html)
-   [constrainedKriging](../packages/constrainedKriging/index.html)
-   [cshapes](../packages/cshapes/index.html)
-   [dbmss](../packages/dbmss/index.html)
-   [DCluster](../packages/DCluster/index.html) (core)
-   [deldir](../packages/deldir/index.html) (core)
-   [diseasemapping](../packages/diseasemapping/index.html)
-   [DSpat](../packages/DSpat/index.html)
-   [ecespa](../packages/ecespa/index.html)
-   [ExceedanceTools](../packages/ExceedanceTools/index.html)
-   [fields](../packages/fields/index.html)
-   [FieldSim](../packages/FieldSim/index.html)
-   [FRK](../packages/FRK/index.html)
-   [gdalUtils](../packages/gdalUtils/index.html)
-   [gdistance](../packages/gdistance/index.html)
-   [gear](../packages/gear/index.html)
-   [geoaxe](../packages/geoaxe/index.html)
-   [geojson](../packages/geojson/index.html)
-   [geojsonio](../packages/geojsonio/index.html)
-   [GEOmap](../packages/GEOmap/index.html)
-   [geomapdata](../packages/geomapdata/index.html)
-   [geonames](../packages/geonames/index.html)
-   [geoR](../packages/geoR/index.html) (core)
-   [geoRglm](../packages/geoRglm/index.html)
-   [georob](../packages/georob/index.html)
-   [geospacom](../packages/geospacom/index.html)
-   [geosphere](../packages/geosphere/index.html)
-   [geospt](../packages/geospt/index.html)
-   [geostatsp](../packages/geostatsp/index.html)
-   [ggmap](../packages/ggmap/index.html)
-   [ggsn](../packages/ggsn/index.html)
-   [glmmBUGS](../packages/glmmBUGS/index.html)
-   [gmt](../packages/gmt/index.html)
-   [Grid2Polygons](../packages/Grid2Polygons/index.html)
-   [GriegSmith](../packages/GriegSmith/index.html)
-   [gstat](../packages/gstat/index.html) (core)
-   [Guerry](../packages/Guerry/index.html)
-   [GWmodel](../packages/GWmodel/index.html)
-   [gwrr](../packages/gwrr/index.html)
-   [hdeco](../packages/hdeco/index.html)
-   [HSAR](../packages/HSAR/index.html)
-   [igraph](../packages/igraph/index.html)
-   [intamap](../packages/intamap/index.html)
-   [ipdw](../packages/ipdw/index.html)
-   [landsat](../packages/landsat/index.html)
-   [latticeDensity](../packages/latticeDensity/index.html)
-   [lawn](../packages/lawn/index.html)
-   [lctools](../packages/lctools/index.html)
-   [leafletR](../packages/leafletR/index.html)
-   [magclass](../packages/magclass/index.html)
-   [mapdata](../packages/mapdata/index.html)
-   [mapmisc](../packages/mapmisc/index.html)
-   [mapproj](../packages/mapproj/index.html)
-   [maps](../packages/maps/index.html)
-   [maptools](../packages/maptools/index.html) (core)
-   [mapview](../packages/mapview/index.html)
-   [marmap](../packages/marmap/index.html)
-   [MBA](../packages/MBA/index.html)
-   [McSpatial](../packages/McSpatial/index.html)
-   [micromap](../packages/micromap/index.html)
-   [ModelMap](../packages/ModelMap/index.html)
-   [ncdf4](../packages/ncdf4/index.html)
-   [ncf](../packages/ncf/index.html)
-   [ngspatial](../packages/ngspatial/index.html)
-   [nlme](../packages/nlme/index.html)
-   [OasisR](../packages/OasisR/index.html)
-   [OpenStreetMap](../packages/OpenStreetMap/index.html)
-   [osmar](../packages/osmar/index.html)
-   [pastecs](../packages/pastecs/index.html)
-   [PBSmapping](../packages/PBSmapping/index.html)
-   [PBSmodelling](../packages/PBSmodelling/index.html)
-   [plotGoogleMaps](../packages/plotGoogleMaps/index.html)
-   [plotKML](../packages/plotKML/index.html)
-   [postGIStools](../packages/postGIStools/index.html)
-   [PReMiuM](../packages/PReMiuM/index.html)
-   [ProbitSpatial](../packages/ProbitSpatial/index.html)
-   [quickmapr](../packages/quickmapr/index.html)
-   [ramps](../packages/ramps/index.html)
-   [RandomFields](../packages/RandomFields/index.html) (core)
-   [rangeMapper](../packages/rangeMapper/index.html)
-   [RArcInfo](../packages/RArcInfo/index.html)
-   [raster](../packages/raster/index.html) (core)
-   [rasterVis](../packages/rasterVis/index.html)
-   [RColorBrewer](../packages/RColorBrewer/index.html) (core)
-   [recmap](../packages/recmap/index.html)
-   [regress](../packages/regress/index.html)
-   [rgbif](../packages/rgbif/index.html)
-   [rgdal](../packages/rgdal/index.html) (core)
-   [rgeos](../packages/rgeos/index.html) (core)
-   [RgoogleMaps](../packages/RgoogleMaps/index.html)
-   [rgrass7](../packages/rgrass7/index.html)
-   [RNetCDF](../packages/RNetCDF/index.html)
-   [rpostgis](../packages/rpostgis/index.html)
-   [RPyGeo](../packages/RPyGeo/index.html)
-   [RQGIS](../packages/RQGIS/index.html)
-   [RSAGA](../packages/RSAGA/index.html)
-   [RSurvey](../packages/RSurvey/index.html)
-   [rtop](../packages/rtop/index.html)
-   [rworldmap](../packages/rworldmap/index.html)
-   [rworldxtra](../packages/rworldxtra/index.html)
-   [S2sls](../packages/S2sls/index.html)
-   [seg](../packages/seg/index.html)
-   [sf](../packages/sf/index.html) (core)
-   [sgeostat](../packages/sgeostat/index.html)
-   [shapefiles](../packages/shapefiles/index.html)
-   [shp2graph](../packages/shp2graph/index.html)
-   [siplab](../packages/siplab/index.html)
-   [smacpod](../packages/smacpod/index.html)
-   [smerc](../packages/smerc/index.html)
-   [sp](../packages/sp/index.html) (core)
-   [spacetime](../packages/spacetime/index.html) (core)
-   [spacom](../packages/spacom/index.html)
-   [spaMM](../packages/spaMM/index.html)
-   [spanel](../packages/spanel/index.html)
-   [sparr](../packages/sparr/index.html)
-   [spatgraphs](../packages/spatgraphs/index.html)
-   [spatial](../packages/spatial/index.html)
-   [spatial.tools](../packages/spatial.tools/index.html)
-   [spatialCovariance](../packages/spatialCovariance/index.html)
-   [SpatialEpi](../packages/SpatialEpi/index.html)
-   [SpatialExtremes](../packages/SpatialExtremes/index.html)
-   [SpatialPosition](../packages/SpatialPosition/index.html)
-   [spatialprobit](../packages/spatialprobit/index.html)
-   [spatialsegregation](../packages/spatialsegregation/index.html)
-   [SpatialTools](../packages/SpatialTools/index.html)
-   [spatstat](../packages/spatstat/index.html) (core)
-   [spatsurv](../packages/spatsurv/index.html)
-   [spBayes](../packages/spBayes/index.html)
-   [spBayesSurv](../packages/spBayesSurv/index.html)
-   [spcosa](../packages/spcosa/index.html)
-   [spdep](../packages/spdep/index.html) (core)
-   [sperrorest](../packages/sperrorest/index.html)
-   [spgrass6](../packages/spgrass6/index.html)
-   [spgwr](../packages/spgwr/index.html)
-   [sphet](../packages/sphet/index.html)
-   [spind](../packages/spind/index.html)
-   [splancs](../packages/splancs/index.html) (core)
-   [splm](../packages/splm/index.html)
-   [spm](../packages/spm/index.html)
-   [spmoran](../packages/spmoran/index.html)
-   [spsann](../packages/spsann/index.html)
-   [spselect](../packages/spselect/index.html)
-   [spsurvey](../packages/spsurvey/index.html)
-   [spTimer](../packages/spTimer/index.html)
-   [SSN](../packages/SSN/index.html)
-   [starma](../packages/starma/index.html)
-   [statebins](../packages/statebins/index.html)
-   [Stem](../packages/Stem/index.html)
-   [stplanr](../packages/stplanr/index.html)
-   [taRifx](../packages/taRifx/index.html)
-   [tgp](../packages/tgp/index.html)
-   [tmap](../packages/tmap/index.html)
-   [trip](../packages/trip/index.html)
-   [tripack](../packages/tripack/index.html)
-   [tripEstimation](../packages/tripEstimation/index.html)
-   [UScensus2000cdp](../packages/UScensus2000cdp/index.html)
-   [UScensus2000tract](../packages/UScensus2000tract/index.html)
-   [vardiag](../packages/vardiag/index.html)
-   [vec2dtransf](../packages/vec2dtransf/index.html)
-   [vegan](../packages/vegan/index.html)
-   [Watersheds](../packages/Watersheds/index.html)
-   [wkb](../packages/wkb/index.html)

### Related links:

-   CRAN Task View: [Environmetrics](Environmetrics.html)
-   [Rgeo: Spatial Statistics with
    R](http://geodacenter.asu.edu/projects/rsp)
-   [R-SIG-Geo mailing
    list](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
