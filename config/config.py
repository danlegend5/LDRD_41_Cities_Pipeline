__author__ = 'Dan Bramich'


##############################################
# General pipeline configuration parameters

# Location of original data sources
original_detectors_file = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Locations/detectors.csv'
original_links_file = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Locations/links.csv'
original_measurements_raw_file = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Measurements/loops_measurements.csv'
original_measurements_arima_file = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Measurements/loops_measurements_arima.csv'
original_measurements_innsbruck_file = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Measurements/loops_measurements_innsbruck_corrected.csv'
original_measurements_zurich_file1 = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Measurements/loops_measurements_zurich_corrected.csv'
original_measurements_zurich_file2 = '/data/dmb20/Loop.Detector.Original.Data.41.Cities/Loop.Detector.Measurements/loops_measurements_zurich_3min_1year.csv'

# Output directory for results
output_dir = '/data/dmb20/Loop.Detector.Reformatted.Data.41.Cities/Results'


##############################################
# Pipeline stage 1b parameters

# Minimum number of good measurements required for each loop detector
min_ngood_ld_measurements = 900
