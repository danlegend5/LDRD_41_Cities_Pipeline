__author__ = 'Dan Bramich'

# Imports
import os
from astropy.table import Table
from LDDA_MFD_Project1_Pipeline.config import config


# Given a set of loop detectors and their properties via the input parameter "ld_locations_table", this function reads in
# the associated loop detector measurements files (of data type "data_type" - acceptable values are "raw" and "ARIMA"),
# and it constructs the corresponding MFD. The function requires specification of the city "city_name" and country
# "country_name" that the loop detectors belong to.
def construct_mfd(ld_locations_table, country_name, city_name, data_type):

    # If there are no loop detectors in the input parameter "ld_locations_table", then return
    nld_locations = len(ld_locations_table)
    if nld_locations == 0: return 0

    # Determine the directory path to the loop detector measurements files based on the required data type (raw or ARIMA)
    if data_type == 'raw':
        input_dir_ld_measurements = os.path.join(config.output_dir, 's1.Loop.Detector.Measurements.Raw', country_name, city_name)
    elif data_type == 'ARIMA':
        input_dir_ld_measurements = os.path.join(config.output_dir, 's1.Loop.Detector.Measurements.ARIMA', country_name, city_name)
    else:
        return 0

    # For each loop detector
    nld_used = 0
    for i in range(nld_locations):

        # Determine the full directory path and filename of the loop detector measurements data file for the current loop
        # detector
        curr_det_id = ld_locations_table['DETECTOR_ID'][i]
        if data_type == 'raw':
            curr_ld_measurements_file = os.path.join(input_dir_ld_measurements, curr_det_id, 'measurements.raw.' + country_name + '.' + city_name + '.' + curr_det_id + '.fits')
        elif data_type == 'ARIMA':
            curr_ld_measurements_file = os.path.join(input_dir_ld_measurements, curr_det_id, 'measurements.ARIMA.' + country_name + '.' + city_name + '.' + curr_det_id + '.fits')

        # If the loop detector measurements data file for the current loop detector does not exist, then move on to the
        # next loop detector (this can be the case for "data_type = ARIMA")
        if not os.path.exists(curr_ld_measurements_file): continue

        # Read in the loop detector measurements data file for the current loop detector
        ld_measurements_table = Table.read(curr_ld_measurements_file, format = 'fits')
        nld_measurements = len(ld_measurements_table)


#### ABOVE FULLY READ AND TESTED


    #
    return 0
