__author__ = 'Dan Bramich'


##############################################
# General pipeline configuration parameters

# Location of original data sources
original_detectors_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Locations/detectors.csv'
original_links_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Locations/links.csv'
original_measurements_raw_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Measurements/loops_measurements.csv'
original_measurements_arima_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Measurements/loops_measurements_arima.csv'

# Output directory for results
output_dir = '/data/dmb20/Loop.Detector.Data.Analysis.MFD.Project.1/Results'

# Mean radius (km) of the Earth to be used for the approximate conversion of sets of longitude and latitude
# coordinates within a city to sets of x and y Cartesian coordinates (taken from the International Union of
# Geodesy and Geophysics - IUGG)
radius_earth = 6371.0
