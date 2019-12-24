#!/usr/bin/env python

__author__ = 'Dan Bramich'

# This script converts the filtered loop detector locations and measurements data created by stage 1a into
# an efficient two-dimensional data format for subsequent analyses. Note that the cities with data from a
# combination of loop detectors and bluetooth detectors are not processed by this script, i.e. 'groningen',
# 'melbourne', and 'utrecht'.

# Imports
import astropy.io.fits
import glob
import numpy
import os
import shutil
from astropy.table import Column
from astropy.table import Table
from astropy.time import TimeDelta
from LDRD_41_Cities_Pipeline.config import config
from LDRD_41_Cities_Pipeline.lib import general_functions

# Create the output directory for the loop detector data
output_dir_ld_data = os.path.join(config.output_dir, 's2a.Loop.Detector.Data')
print('')
print('Creating the output directory for the loop detector data: ' + output_dir_ld_data)
if os.path.exists(output_dir_ld_data): shutil.rmtree(output_dir_ld_data)
os.makedirs(output_dir_ld_data)

# Copy across the loop detector locations data files and the summary statistics text files that have already
# been produced in stage 1a
print('Copying across the loop detector locations data files and the summary statistics text files from stage 1a...')
file_list = glob.glob(os.path.join(config.output_dir, 's1a.Loop.Detector.Locations', 'detectors.LD.Flow.LD.*.fits'))
for file in file_list: shutil.copy(file, output_dir_ld_data)
shutil.copy(os.path.join(config.output_dir, 's1a.Loop.Detector.Locations', 'summary.statistics.of.detectors.per.city.txt'), output_dir_ld_data)
shutil.copy(os.path.join(config.output_dir, 's1a.Loop.Detector.Locations', 'summary.statistics.of.measurements.per.city.detector.txt'), output_dir_ld_data)

# Read in the summary statistics text file for the loop detector measurements
summary_stats_of_measurements_per_city_detector_file = os.path.join(output_dir_ld_data, 'summary.statistics.of.measurements.per.city.detector.txt')
print('')
print('Reading in the summary statistics text file: ' + summary_stats_of_measurements_per_city_detector_file)
with open(summary_stats_of_measurements_per_city_detector_file, mode = 'r') as ssfile:
    sslines = [line for line in ssfile if line.rstrip('\n')]
ss_city = []
ss_first_time_stamp = []
ss_last_time_stamp = []
ss_interval = []
for line in sslines:
    bits = line.split()
    if bits[0] == '#': continue
    ss_city.append(bits[1])
    if bits[21] != '-':
        tmp_list = [bits[5], bits[21]]
        tmp_list.sort()
        ss_first_time_stamp.append(tmp_list[0])
    else:
        ss_first_time_stamp.append(bits[5])
    if bits[22] != '-':
        tmp_list = [bits[6], bits[22]]
        tmp_list.sort()
        ss_last_time_stamp.append(tmp_list[1])
    else:
        ss_last_time_stamp.append(bits[6])
    ss_interval.append(bits[9])
ss_city = numpy.array(ss_city)
ss_first_time_stamp = numpy.array(ss_first_time_stamp)
ss_last_time_stamp = numpy.array(ss_last_time_stamp)
ss_interval = numpy.array(ss_interval)
print('No. of data lines read in: ' + str(len(ss_city)))

# Determine the list of loop detector locations data files
print('')
print('Determining the list of loop detector locations data files...')
file_list = glob.glob(os.path.join(output_dir_ld_data, 'detectors.*.fits'))
file_list.sort()
print('No. of files: ' + str(len(file_list)))

# For each loop detector locations data file
empty_str = '                                                                                '
for file in file_list:

    # Read in the loop detector locations data file
    print('')
    print('Reading in the loop detector locations data file: ' + file)
    ld_locations_table = Table.read(file, format = 'fits')
    nld_locations = len(ld_locations_table)
    print('No. of entries: ' + str(nld_locations))

    # Determine the corresponding data source, country, and city
    path, basename_with_ext = os.path.split(file)
    basename_bits = basename_with_ext.split('.')
    if len(basename_bits) == 10:
        data_source = '.'.join(basename_bits[1:7])
        country_name = basename_bits[7]
        city_name = basename_bits[8]
    else:
        data_source = '.'.join(basename_bits[1:5])
        country_name = basename_bits[5]
        city_name = basename_bits[6]

    # Determine the first and last time stamps, and the aggregation time-interval, for the current city
    # (N.B: Time stamps indicate the start of an aggregation time-interval)
    print('Defining the time grid on which to organise the loop detector measurements data...')
    selection = ss_city == city_name
    curr_first_time_stamp = numpy.sort(ss_first_time_stamp[selection])[0]
    curr_first_time_stamp = general_functions.convert_time_stamp_to_time_object(curr_first_time_stamp)
    curr_last_time_stamp = numpy.sort(ss_last_time_stamp[selection])[-1]
    curr_last_time_stamp = general_functions.convert_time_stamp_to_time_object(curr_last_time_stamp)
    curr_interval = TimeDelta(numpy.median(numpy.int32(ss_interval[selection])), format = 'sec')
    curr_half_interval = 0.5*curr_interval

    # Define the time grid on which to organise the loop detector measurements data
    ntime_bins = numpy.int32(numpy.ceil((curr_last_time_stamp + curr_half_interval - curr_first_time_stamp)/curr_interval))
    time_bin_limits = curr_first_time_stamp + (numpy.arange(0, ntime_bins + 1, 1, dtype = numpy.int32)*curr_interval)
    time_mid_points = time_bin_limits[0:ntime_bins] + curr_half_interval
    print('Earliest time bin limit: ' + str(time_bin_limits[0]))
    print('Latest time bin limit:   ' + str(time_bin_limits[-1]))
    print('Time bin interval (s):   ' + str(curr_interval))
    print('No. of time bins:        ' + str(ntime_bins))

    # For each loop detector
    flow_im = numpy.zeros((nld_locations, ntime_bins), dtype = numpy.float64) - 1.0
    if data_source != 'LD.Flow.LD.Speed': occ_im = numpy.zeros((nld_locations, ntime_bins), dtype = numpy.float64) - 1.0
    error_im = numpy.ones((nld_locations, ntime_bins), dtype = numpy.int32)
    if data_source != 'LD.Flow.LD.Occupancy': speed_im = numpy.zeros((nld_locations, ntime_bins), dtype = numpy.float64) - 1.0
    arima_flow_im = numpy.zeros((nld_locations, ntime_bins), dtype = numpy.float64) - 1.0
    if data_source != 'LD.Flow.LD.Speed': arima_occ_im = numpy.zeros((nld_locations, ntime_bins), dtype = numpy.float64) - 1.0
    arima_error_im = numpy.ones((nld_locations, ntime_bins), dtype = numpy.int32)
    if data_source != 'LD.Flow.LD.Occupancy': arima_speed_im = numpy.zeros((nld_locations, ntime_bins), dtype = numpy.float64) - 1.0
    arima_data_tag = 0
    for i in range(nld_locations):

        # Determine the full directory path and filename of the corresponding loop detector measurements
        # data file (raw)
        curr_detector_id = ld_locations_table['DETECTOR_ID'][i]
        print('Reading in the loop detector measurements data (raw and ARIMA) and incorporating them into the grid for the detector: ' + curr_detector_id)
        ld_measurements_file_raw = 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.' + curr_detector_id + '.fits'
        ld_measurements_file_raw = os.path.join(config.output_dir, 's1a.Loop.Detector.Measurements.Raw', data_source, country_name, city_name, curr_detector_id, ld_measurements_file_raw)

        # Read in the corresponding loop detector measurements data file (raw)
        ld_measurements_table_raw = Table.read(ld_measurements_file_raw, format = 'fits')
        nld_measurements_raw = len(ld_measurements_table_raw)

        # For each loop detector measurement
        curr_time_bin_sub = 0
        for j in range(nld_measurements_raw):

            # If the current loop detector measurement is flagged with "ERROR_FLAG = 1", then move on to
            # the next loop detector measurement
            if ld_measurements_table_raw['ERROR_FLAG'][j] == 1: continue

            # Determine the mid-point of the aggregation time-interval for the current loop detector
            # measurement
            date_interval_start_str = ld_measurements_table_raw['DATE'][j] + 'T' + '{:05d}'.format(ld_measurements_table_raw['INTERVAL_START'][j]).strip()
            curr_time_stamp = general_functions.convert_time_stamp_to_time_object(date_interval_start_str)
            curr_time_stamp += curr_half_interval

            # Determine which time bin the current loop detector measurement should be assigned to
            curr_time_diff = curr_time_stamp - time_bin_limits[curr_time_bin_sub + 1]
            if curr_time_diff >= 0:
                curr_time_bin_sub += numpy.int32(numpy.floor(curr_time_diff/curr_interval)) + 1

            # Check that the designated time bin has not already been filled for the current loop
            # detector
            if error_im[i, curr_time_bin_sub] == 0:
                print('ERROR - The designated time bin for the current loop detector measurement (raw) has already been filled (CITY, DETECTOR_ID, DATE, INTERVAL_START): ' + city_name + ', ' + \
                       curr_detector_id + ', ' + ld_measurements_table_raw['DATE'][j] + ', ' + str(ld_measurements_table_raw['INTERVAL_START'][j]))
                exit()

            # Assign the current loop detector measurement to the designated time bin
            flow_im[i, curr_time_bin_sub] = ld_measurements_table_raw['FLOW'][j]
            if data_source != 'LD.Flow.LD.Speed': occ_im[i, curr_time_bin_sub] = ld_measurements_table_raw['OCCUPANCY'][j]
            error_im[i, curr_time_bin_sub] = 0
            if data_source != 'LD.Flow.LD.Occupancy': speed_im[i, curr_time_bin_sub] = ld_measurements_table_raw['SPEED'][j]

        # Determine the full directory path and filename of the corresponding loop detector measurements
        # data file (ARIMA)
        ld_measurements_file_arima = 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.' + curr_detector_id + '.fits'
        ld_measurements_file_arima = os.path.join(config.output_dir, 's1a.Loop.Detector.Measurements.ARIMA', data_source, country_name, city_name, curr_detector_id, ld_measurements_file_arima)

        # If the corresponding loop detector measurements data file (ARIMA) does not exist, then move on to
        # the next loop detector
        if not os.path.exists(ld_measurements_file_arima): continue

        # Read in the corresponding loop detector measurements data file (ARIMA)
        ld_measurements_table_arima = Table.read(ld_measurements_file_arima, format = 'fits')
        nld_measurements_arima = len(ld_measurements_table_arima)
        arima_data_tag = 1

        # For each loop detector measurement
        curr_time_bin_sub = 0
        for j in range(nld_measurements_arima):

            # If the current loop detector measurement is flagged with "ERROR_FLAG = 1", then move on to
            # the next loop detector measurement
            if ld_measurements_table_arima['ERROR_FLAG'][j] == 1: continue

            # Determine the mid-point of the aggregation time-interval for the current loop detector
            # measurement
            date_interval_start_str = ld_measurements_table_arima['DATE'][j] + 'T' + '{:05d}'.format(ld_measurements_table_arima['INTERVAL_START'][j]).strip()
            curr_time_stamp = general_functions.convert_time_stamp_to_time_object(date_interval_start_str)
            curr_time_stamp += curr_half_interval

            # Determine which time bin the current loop detector measurement should be assigned to
            curr_time_diff = curr_time_stamp - time_bin_limits[curr_time_bin_sub + 1]
            if curr_time_diff >= 0:
                curr_time_bin_sub += numpy.int32(numpy.floor(curr_time_diff/curr_interval)) + 1

            # Check that the designated time bin has not already been filled for the current loop
            # detector
            if arima_error_im[i, curr_time_bin_sub] == 0:
                print('ERROR - The designated time bin for the current loop detector measurement (ARIMA) has already been filled (CITY, DETECTOR_ID, DATE, INTERVAL_START): ' + city_name + ', ' + \
                       curr_detector_id + ', ' + ld_measurements_table_arima['DATE'][j] + ', ' + str(ld_measurements_table_arima['INTERVAL_START'][j]))
                exit()

            # Assign the current loop detector measurement to the designated time bin
            arima_flow_im[i, curr_time_bin_sub] = ld_measurements_table_arima['ARIMA_FLOW'][j]
            if data_source != 'LD.Flow.LD.Speed': arima_occ_im[i, curr_time_bin_sub] = ld_measurements_table_arima['ARIMA_OCCUPANCY'][j]
            arima_error_im[i, curr_time_bin_sub] = 0
            if data_source != 'LD.Flow.LD.Occupancy': arima_speed_im[i, curr_time_bin_sub] = ld_measurements_table_arima['ARIMA_SPEED'][j]

    # Write out the two-dimensional arrays of loop detector measurements for the current city
    print('Writing out the two-dimensional arrays of loop detector measurements...')
    flow_file = os.path.join(output_dir_ld_data, 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.flow.fits')
    astropy.io.fits.writeto(flow_file, flow_im, checksum = True)
    if data_source != 'LD.Flow.LD.Speed':
        occ_file = os.path.join(output_dir_ld_data, 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.occupancy.fits')
        astropy.io.fits.writeto(occ_file, occ_im, checksum = True)
    error_file = os.path.join(output_dir_ld_data, 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.error_flag.fits')
    astropy.io.fits.writeto(error_file, error_im, checksum = True)
    if data_source != 'LD.Flow.LD.Occupancy':
        speed_file = os.path.join(output_dir_ld_data, 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.speed.fits')
        astropy.io.fits.writeto(speed_file, speed_im, checksum = True)
    if arima_data_tag == 1:
        arima_flow_file = os.path.join(output_dir_ld_data, 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.flow.fits')
        astropy.io.fits.writeto(arima_flow_file, arima_flow_im, checksum = True)
        if data_source != 'LD.Flow.LD.Speed':
            arima_occ_file = os.path.join(output_dir_ld_data, 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.occupancy.fits')
            astropy.io.fits.writeto(arima_occ_file, arima_occ_im, checksum = True)
        arima_error_file = os.path.join(output_dir_ld_data, 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.error_flag.fits')
        astropy.io.fits.writeto(arima_error_file, arima_error_im, checksum = True)
        if data_source != 'LD.Flow.LD.Occupancy':
            arima_speed_file = os.path.join(output_dir_ld_data, 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.speed.fits')
            astropy.io.fits.writeto(arima_speed_file, arima_speed_im, checksum = True)

    # Prepare a data table containing the time stamps defining the time grid
    timestamps_table_file = os.path.join(output_dir_ld_data, 'timestamps.' + data_source + '.' + country_name + '.' + city_name + '.fits')
    print('Writing out the time-stamps data table: ' + timestamps_table_file)
    timestamps_table = Table([Column(data = [empty_str]*ntime_bins, name = 'LOWER_BIN_LIMIT_TIMESTAMP'),
                              Column(data = [empty_str]*ntime_bins, name = 'UPPER_BIN_LIMIT_TIMESTAMP'),
                              Column(data = [empty_str]*ntime_bins, name = 'BIN_MIDPOINT_TIMESTAMP'),
                              Column(data = [empty_str]*ntime_bins, name = 'LOWER_BIN_LIMIT_MJD'),
                              Column(data = [empty_str]*ntime_bins, name = 'UPPER_BIN_LIMIT_MJD'),
                              Column(data = [empty_str]*ntime_bins, name = 'BIN_MIDPOINT_MJD')])
    for i in range(ntime_bins):
        timestamps_table['LOWER_BIN_LIMIT_TIMESTAMP'][i] = str(time_bin_limits[i])
        timestamps_table['UPPER_BIN_LIMIT_TIMESTAMP'][i] = str(time_bin_limits[i + 1])
        timestamps_table['BIN_MIDPOINT_TIMESTAMP'][i] = str(time_mid_points[i])
        timestamps_table['LOWER_BIN_LIMIT_MJD'][i] = '{:25.8f}'.format(time_bin_limits[i].mjd).strip()
        timestamps_table['UPPER_BIN_LIMIT_MJD'][i] = '{:25.8f}'.format(time_bin_limits[i + 1].mjd).strip()
        timestamps_table['BIN_MIDPOINT_MJD'][i] = '{:25.8f}'.format(time_mid_points[i].mjd).strip()

    # Write out the time-stamps data table
    timestamps_table.write(timestamps_table_file, format = 'fits')

# Finish
print('')
print('Finished!')
