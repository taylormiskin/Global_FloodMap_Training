# Global FloodMapping Test case
Test case using the Automated Rating Curve (ARC) and Curve2Flood

The Automated Rating Curve (ARC) is a python-based model that generates rating-curve information (depth, velocity, top-width) along streams.  ARC is typically used for large-scale applications (e.g., the Yellowstone BAsin), but can also be applied to smaller domains.  Inputs to ARC are peak- and base-flow data for each stream reach (often from continental-scale flow models) and geospatial data (land cover, DEM, stream shapefile).  ARC uses a simple approach to simulate a bathymetry profile within the DEM (https://doi.org/10.1016/j.jhydrol.2023.129769).  Using a volume-fill approach of Manning's Equation, mulitple flow rates are simulated.  Parameters (a and b) for a power-function are then generated where depth, velocity, and top-width are calculated as a function of flow (Q).  For example, Depth = a *Q^b.

The ARC-generated rating curves can then be used to create flood maps.  Here, we use ERDC's FloodSpreader model to show how flood maps can be generated when given a flow rate.

# Conda Environment:
    Typically AutoRoute and FloodSpreader are run using Anaconda Prompts (https://www.anaconda.com/download)
    To activate the "ARC" environment open an Anaconda Command Prompt and type "conda env create -f environment_ARC.yml"

# Model and Data Sources:
    DEM data is from 
    Land Cover data is from 
    Streamlines are from GeoGLoWS - http://geoglows-v2.s3-website-us-west-2.amazonaws.com/#streams/
    Return period flow data is from GeoGLoWS - http://geoglows-v2-retrospective.s3-website-us-west-2.amazonaws.com/#return-periods/

# Python Scripts
    Curve2Flood: https://github.com/MikeFHS/curve2flood
    ARC [Automated Rating Curve]: https://github.com/MikeFHS/automated-rating-curve
    Similar Test Case created by FSH [Mike Follum]: https://sites.google.com/follumhydro.com/automated-rating-curve-arc
    Link to Github test case: https://github.com/MikeFHS/AutomatedRatingCurve_TestCase
    

# Model/Data/Files are AS-IS with No Warranty:
Any scripts, files, etc. provided in this document are for testing purposes only and are provided "as is" and without any warranty of any kind.
