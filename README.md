This repository contains the raw data and analysis code (in R) for the article:

Marchand, P., Boenke, M. and Green, D.M. "A stochastic movement model reproduces patterns of site fidelity and long-distance dispersal in a population of Fowler's Toads (Anaxyrus fowleri)." (in revision, Ecological Modelling)

The radiotracking data was originally collected as part of a study led by Morgan Boenke and supervised by David Green, reported in:

Boenke, M. 2011. Terrestrial habitat and ecology of Fowler’s toads (Anaxyrus fowleri). M.Sc. thesis, Department of Biology, McGill University, Montreal, Canada.

## Data files (data/ subfolder)

**radio2009.csv**, **radio2010.csv**: Radiotracking locations of Fowler's toads at Long Point, Ontario, measured over two observation seasons (summer 2009 and 2010). 

**Final Waterline.shp** (and associated files): ESRI shapefile delineating the shoreline of Lake Erie during the sampling period.

*Note*: The northing and easting values in the radiotracking files, as well as the waterline shapefile, are in UTM zone 17N coordinates with NAD83 datum.

## Code files (main folder)

**prep_day_location_data.R**: Script run prior to the main script; it combines the daytime refuge locations from the two files of radiotracking data and rotates the coordinates so that the *x*-axis is parallel to the lake shore.

**toad_move_funcs.R**: Functions called by the main script, including those implementing the simulation model, and calculating statistics from the simulation output.

**toad_move_main.R**: Main script used to perform the analyses and produce the results and graphs found in the article.


 