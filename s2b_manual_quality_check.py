#!/usr/bin/env python

__author__ = 'Dan Bramich'

# This script allows the user to perform a visual quality check on the flow-occupancy empirical fundamental
# diagram for each loop detector and to reject any loop detectors for which the measurement data have
# significant/obvious problems. This script processes the filtered loop detector locations and measurements
# data created by stage 1b.

# Imports
import glob
import numpy
import os
import shutil
import matplotlib.pyplot as plt
from astropy.table import Table
from LDRD_41_Cities_Pipeline.config import config

# Create the output directory for the loop detector locations data tables
output_dir_ld_locations = os.path.join(config.output_dir, 's2b.Loop.Detector.Locations')
print('')
print('Creating the output directory for the loop detector locations data tables: ' + output_dir_ld_locations)
if os.path.exists(output_dir_ld_locations): shutil.rmtree(output_dir_ld_locations)
os.makedirs(output_dir_ld_locations)

# Create the output directory for the loop detector measurements data tables (raw)
output_dir_ld_measurements_raw = os.path.join(config.output_dir, 's2b.Loop.Detector.Measurements.Raw')
print('Creating the output directory for the loop detector measurements data tables (raw): ' + output_dir_ld_measurements_raw)
if os.path.exists(output_dir_ld_measurements_raw): shutil.rmtree(output_dir_ld_measurements_raw)
os.makedirs(output_dir_ld_measurements_raw)

# Create the output directory for the loop detector measurements data tables (ARIMA)
output_dir_ld_measurements_arima = os.path.join(config.output_dir, 's2b.Loop.Detector.Measurements.ARIMA')
print('Creating the output directory for the loop detector measurements data tables (ARIMA): ' + output_dir_ld_measurements_arima)
if os.path.exists(output_dir_ld_measurements_arima): shutil.rmtree(output_dir_ld_measurements_arima)
os.makedirs(output_dir_ld_measurements_arima)

# Initialise the text file that will be used to contain the summary statistics of the loop detectors on a per
# city basis
summary_stats_of_detectors_per_city_file = os.path.join(output_dir_ld_locations, 'summary.statistics.of.detectors.per.city.txt')
print('')
print('Initialising the text file that will be used to contain the summary statistics of the loop detectors: ' + summary_stats_of_detectors_per_city_file)
with open(summary_stats_of_detectors_per_city_file, mode = 'w') as ssfile:
    ssfile.write('# Country : City : No. Of Loop Detectors (LDs) : Min. Longitude (deg) Of LD Locations : Max. Longitude (deg) Of LD Locations : Min. Latitude (deg) Of LD Locations : ' + \
                 'Max. Latitude (deg) Of LD Locations : Min. LD Road Length (km) : Med. LD Road Length (km) : Max. LD Road Length (km) : Min. Normalised Position Of LD Along Road : ' + \
                 'Med. Normalised Position Of LD Along Road : Max. Normalised Position Of LD Along Road : No. LD Roads Classed "living_street" : No. LD Roads Classed "motorway" : ' + \
                 'No. LD Roads Classed "motorway_link" : No. LD Roads Classed "primary" : No. LD Roads Classed "primary_link" : No. LD Roads Classed "residential" : ' + \
                 'No. LD Roads Classed "secondary" : No. LD Roads Classed "secondary_link" : No. LD Roads Classed "service" : No. LD Roads Classed "tertiary" : ' + \
                 'No. LD Roads Classed "tertiary_link" : No. LD Roads Classed "trunk" : No. LD Roads Classed "trunk_link" : No. LD Roads Classed "unclassified" : ' + \
                 'Min. LD Road Speed Limit (km/h) : Med. LD Road Speed Limit (km/h) : Max. LD Road Speed Limit (km/h) : Min. No. Lanes Covered By LD : Max. No. Lanes Covered By LD : ' + \
                 'Min. Longitude (deg) Of LD Roads : Max. Longitude (deg) Of LD Roads : Min. Latitude (deg) Of LD Roads : Max. Latitude (deg) Of LD Roads\n')

# Initialise the text file that will be used to contain the summary statistics of the loop detector
# measurements (raw and ARIMA) on a per city and detector basis
summary_stats_of_measurements_per_city_detector_file = os.path.join(output_dir_ld_locations, 'summary.statistics.of.measurements.per.city.detector.txt')
print('Initialising the text file that will be used to contain the summary statistics of the loop detector measurements: ' + summary_stats_of_measurements_per_city_detector_file)
with open(summary_stats_of_measurements_per_city_detector_file, mode = 'w') as ssfile:
    ssfile.write('# Country : City : Detector ID : No. Of LD Measurements (raw) : No. Of Good LD Measurements (raw) : First Time Stamp (raw) : Last Time Stamp (raw) : No. Of Dates (raw) : ' + \
                 'Interval Lengths (s; raw) : Med. Interval Length (s; raw) : Min. Flow (veh/h; raw) : Med. Flow (veh/h; raw) : Max. Flow (veh/h; raw) : Min. Occupancy (raw) : ' + \
                 'Med. Occupancy (raw) : Max. Occupancy (raw) : Min. Average-Speed (km/h; raw) : Med. Average-Speed (km/h; raw) : Max. Average-Speed (km/h; raw) : No. Of LD Measurements (ARIMA) : ' + \
                 'No. Of Good LD Measurements (ARIMA) : First Time Stamp (ARIMA) : Last Time Stamp (ARIMA) : No. Of Dates (ARIMA) : Min. Flow (veh/h; ARIMA) : Med. Flow (veh/h; ARIMA) : ' + \
                 'Max. Flow (veh/h; ARIMA) : Min. Occupancy (ARIMA) : Med. Occupancy (ARIMA) : Max. Occupancy (ARIMA) : Min. Average-Speed (km/h; ARIMA) : Med. Average-Speed (km/h; ARIMA) : ' + \
                 'Max. Average-Speed (km/h; ARIMA)\n')

# Determine the list of loop detector locations data files
print('')
print('Determining the list of loop detector locations data files...')
file_list = glob.glob(os.path.join(config.output_dir, 's1b.Loop.Detector.Locations', 'detectors.*.fits'))
file_list.sort()
print('No. of files: ' + str(len(file_list)))

# For each loop detector locations data file
empty_str = '                                                                                '
for file in file_list:

    # Read in the loop detector locations data file
    print('')
    print('>-------------------------------------------------------------------------------------------------------------<')
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

    # For each loop detector
    selection = numpy.ones(nld_locations, dtype = bool)
    nrejected = 0
    nld_measurements_raw_all = numpy.zeros(nld_locations, dtype = numpy.int32)
    nld_measurements_raw_good = numpy.zeros(nld_locations, dtype = numpy.int32)
    first_time_stamp_raw = numpy.array([empty_str]*nld_locations)
    last_time_stamp_raw = numpy.array([empty_str]*nld_locations)
    first_time_stamp_raw[:] = '-'
    last_time_stamp_raw[:] = '-'
    ndates_raw = numpy.zeros(nld_locations, dtype = numpy.int32)
    interval_lengths_raw = numpy.array([empty_str]*nld_locations)
    interval_lengths_raw[:] = '-'
    med_interval_length_raw = numpy.zeros(nld_locations, dtype = numpy.int32)
    min_flow_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    med_flow_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    max_flow_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    min_occupancy_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    med_occupancy_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    max_occupancy_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    min_speed_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    med_speed_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    max_speed_raw = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    nld_measurements_arima_all = numpy.zeros(nld_locations, dtype = numpy.int32)
    nld_measurements_arima_good = numpy.zeros(nld_locations, dtype = numpy.int32)
    first_time_stamp_arima = numpy.array([empty_str]*nld_locations)
    last_time_stamp_arima = numpy.array([empty_str]*nld_locations)
    first_time_stamp_arima[:] = '-'
    last_time_stamp_arima[:] = '-'
    ndates_arima = numpy.zeros(nld_locations, dtype = numpy.int32) - 1
    min_flow_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    med_flow_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    max_flow_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    min_occupancy_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    med_occupancy_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    max_occupancy_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    min_speed_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    med_speed_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    max_speed_arima = numpy.zeros(nld_locations, dtype = numpy.float64) - 1.0
    for i in range(nld_locations):

        # Determine the full directory path and filename of the corresponding loop detector measurements data
        # file (raw)
        print('')
        print('Current loop detector ID: ' + ld_locations_table['DETECTOR_ID'][i])
        ld_measurements_file_raw = 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.' + ld_locations_table['DETECTOR_ID'][i] + '.fits'
        ld_measurements_file_raw = os.path.join(config.output_dir, 's1b.Loop.Detector.Measurements.Raw', data_source, country_name, city_name, ld_locations_table['DETECTOR_ID'][i], ld_measurements_file_raw)

        # Read in the corresponding loop detector measurements data file (raw)
        ld_measurements_table_raw = Table.read(ld_measurements_file_raw, format = 'fits')
        nld_measurements_raw = len(ld_measurements_table_raw)

        # Determine the subset of good loop detector measurements (raw)
        good_ld_measurements_table_raw = ld_measurements_table_raw[ld_measurements_table_raw['ERROR_FLAG'] == 0]
        ngood = len(good_ld_measurements_table_raw)
        print('No. of good measurements: ' + str(ngood))

        # Plot the flow-occupancy empirical fundamental diagram
        print('Plotting the flow-occupancy FD...')
        fig = plt.figure(dpi = 150, figsize = (12, 8))
        plt.xlim = (0.0, 1.0)
        plt.ylim = (0.0, 1.02*numpy.max(good_ld_measurements_table_raw['FLOW']))
        plt.rcParams.update({'font.size' : 5})
        plt.scatter(good_ld_measurements_table_raw['OCCUPANCY'], good_ld_measurements_table_raw['FLOW'], s = 2.0, marker = '.', color = 'black')
        plt.show()

        # Ask the user whether to accept or reject the loop detector measurements data
        decision = input('Keep? (y/n): ')
        if decision != 'y':
            selection[i] = False
            nrejected += 1
            continue

        # Compute summary statistics for the loop detector measurements (raw) for the current loop detector
        nld_measurements_raw_all[i] = nld_measurements_raw
        nld_measurements_raw_good[i] = ngood
        first_time_stamp_raw[i] = ld_measurements_table_raw['DATE'][0] + 'T' + '{:05d}'.format(ld_measurements_table_raw['INTERVAL_START'][0]).strip()
        last_time_stamp_raw[i] = ld_measurements_table_raw['DATE'][-1] + 'T' + '{:05d}'.format(ld_measurements_table_raw['INTERVAL_START'][-1]).strip()
        ndates_raw[i] = len(numpy.unique(ld_measurements_table_raw['DATE']))
        if nld_measurements_raw > 1:
            intervals_raw = ld_measurements_table_raw['INTERVAL_START'][1:nld_measurements_raw] - ld_measurements_table_raw['INTERVAL_START'][0:(nld_measurements_raw - 1)]
            intervals_raw_uniq = numpy.unique(intervals_raw)
            nintervals_raw_uniq = len(intervals_raw_uniq)
            for j in range(nintervals_raw_uniq):
                if intervals_raw_uniq[j] < 0: intervals_raw_uniq[j] += 86400
            intervals_raw_uniq = numpy.unique(intervals_raw_uniq)
            nintervals_raw_uniq = len(intervals_raw_uniq)
            for j in range(nintervals_raw_uniq):
                if j == 0:
                    intervals_raw_uniq_str = str(intervals_raw_uniq[j])
                else:
                    intervals_raw_uniq_str += ',' + str(intervals_raw_uniq[j])
                    if len(intervals_raw_uniq_str) > 80:
                        intervals_raw_uniq_str = 'TOO_MANY'
                        break
            interval_lengths_raw[i] = intervals_raw_uniq_str
            med_interval_length_raw[i] = numpy.int32(numpy.median(intervals_raw))
        else:
            interval_lengths_raw[i] = '-'
            med_interval_length_raw[i] = -1
        min_flow_raw[i] = numpy.min(good_ld_measurements_table_raw['FLOW'])
        med_flow_raw[i] = numpy.median(good_ld_measurements_table_raw['FLOW'])
        max_flow_raw[i] = numpy.max(good_ld_measurements_table_raw['FLOW'])
        min_occupancy_raw[i] = numpy.min(good_ld_measurements_table_raw['OCCUPANCY'])
        med_occupancy_raw[i] = numpy.median(good_ld_measurements_table_raw['OCCUPANCY'])
        max_occupancy_raw[i] = numpy.max(good_ld_measurements_table_raw['OCCUPANCY'])
        if data_source == 'LD.Flow.LD.Occupancy.LD.Speed':
            min_speed_raw[i] = numpy.min(good_ld_measurements_table_raw['SPEED'])
            med_speed_raw[i] = numpy.median(good_ld_measurements_table_raw['SPEED'])
            max_speed_raw[i] = numpy.max(good_ld_measurements_table_raw['SPEED'])

        # Write out the loop detector measurements (raw) data table for the current loop detector
        curr_output_dir = os.path.join(output_dir_ld_measurements_raw, data_source, country_name, city_name, ld_locations_table['DETECTOR_ID'][i])
        curr_output_file = os.path.join(curr_output_dir, 'measurements.raw.' + data_source + '.' + country_name + '.' + city_name + '.' + ld_locations_table['DETECTOR_ID'][i] + '.fits')
        os.makedirs(curr_output_dir)
        ld_measurements_table_raw.write(curr_output_file, format = 'fits')

        # Determine the full directory path and filename of the corresponding loop detector measurements data
        # file (ARIMA)
        ld_measurements_file_arima = 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.' + ld_locations_table['DETECTOR_ID'][i] + '.fits'
        ld_measurements_file_arima = os.path.join(config.output_dir, 's1b.Loop.Detector.Measurements.ARIMA', data_source, country_name, city_name, ld_locations_table['DETECTOR_ID'][i], ld_measurements_file_arima)

        # If the corresponding loop detector measurements data file (ARIMA) does not exist, then move on to the
        # next loop detector
        if not os.path.exists(ld_measurements_file_arima): continue

        # Read in the corresponding loop detector measurements data file (ARIMA)
        ld_measurements_table_arima = Table.read(ld_measurements_file_arima, format = 'fits')
        nld_measurements_arima = len(ld_measurements_table_arima)

        # Determine the subset of good loop detector measurements (ARIMA)
        good_ld_measurements_table_arima = ld_measurements_table_arima[ld_measurements_table_arima['ERROR_FLAG'] == 0]
        ngood = len(good_ld_measurements_table_arima)

        # Compute summary statistics for the loop detector measurements (ARIMA) for the current loop detector
        nld_measurements_arima_all[i] = nld_measurements_arima
        nld_measurements_arima_good[i] = ngood
        first_time_stamp_arima[i] = ld_measurements_table_arima['DATE'][0] + 'T' + '{:05d}'.format(ld_measurements_table_arima['INTERVAL_START'][0]).strip()
        last_time_stamp_arima[i] = ld_measurements_table_arima['DATE'][-1] + 'T' + '{:05d}'.format(ld_measurements_table_arima['INTERVAL_START'][-1]).strip()
        ndates_arima[i] = len(numpy.unique(ld_measurements_table_arima['DATE']))
        min_flow_arima[i] = numpy.min(good_ld_measurements_table_arima['ARIMA_FLOW'])
        med_flow_arima[i] = numpy.median(good_ld_measurements_table_arima['ARIMA_FLOW'])
        max_flow_arima[i] = numpy.max(good_ld_measurements_table_arima['ARIMA_FLOW'])
        min_occupancy_arima[i] = numpy.min(good_ld_measurements_table_arima['ARIMA_OCCUPANCY'])
        med_occupancy_arima[i] = numpy.median(good_ld_measurements_table_arima['ARIMA_OCCUPANCY'])
        max_occupancy_arima[i] = numpy.max(good_ld_measurements_table_arima['ARIMA_OCCUPANCY'])
        if data_source == 'LD.Flow.LD.Occupancy.LD.Speed':
            min_speed_arima[i] = numpy.min(good_ld_measurements_table_arima['ARIMA_SPEED'])
            med_speed_arima[i] = numpy.median(good_ld_measurements_table_arima['ARIMA_SPEED'])
            max_speed_arima[i] = numpy.max(good_ld_measurements_table_arima['ARIMA_SPEED'])

        # Write out the loop detector measurements (ARIMA) data table for the current loop detector
        curr_output_dir = os.path.join(output_dir_ld_measurements_arima, data_source, country_name, city_name, ld_locations_table['DETECTOR_ID'][i])
        curr_output_file = os.path.join(curr_output_dir, 'measurements.ARIMA.' + data_source + '.' + country_name + '.' + city_name + '.' + ld_locations_table['DETECTOR_ID'][i] + '.fits')
        os.makedirs(curr_output_dir)
        ld_measurements_table_arima.write(curr_output_file, format = 'fits')

    # Report
    print('')
    print('No. of loop detectors rejected: ' + str(nrejected))

    # Filter out the rejected loop detectors from the loop detector locations data table
    nld_locations = numpy.count_nonzero(selection)
    print('Final no. of loop detectors with acceptable data: ' + str(nld_locations))
    if nld_locations == 0:
        print('WARNING - There are no loop detectors with acceptable data for: ' + city_name)
        continue
    ld_locations_table = ld_locations_table[selection]
    nld_measurements_raw_all = nld_measurements_raw_all[selection]
    nld_measurements_raw_good = nld_measurements_raw_good[selection]
    first_time_stamp_raw = first_time_stamp_raw[selection]
    last_time_stamp_raw = last_time_stamp_raw[selection]
    ndates_raw = ndates_raw[selection]
    interval_lengths_raw = interval_lengths_raw[selection]
    med_interval_length_raw = med_interval_length_raw[selection]
    min_flow_raw = min_flow_raw[selection]
    med_flow_raw = med_flow_raw[selection]
    max_flow_raw = max_flow_raw[selection]
    min_occupancy_raw = min_occupancy_raw[selection]
    med_occupancy_raw = med_occupancy_raw[selection]
    max_occupancy_raw = max_occupancy_raw[selection]
    min_speed_raw = min_speed_raw[selection]
    med_speed_raw = med_speed_raw[selection]
    max_speed_raw = max_speed_raw[selection]
    nld_measurements_arima_all = nld_measurements_arima_all[selection]
    nld_measurements_arima_good = nld_measurements_arima_good[selection]
    first_time_stamp_arima = first_time_stamp_arima[selection]
    last_time_stamp_arima = last_time_stamp_arima[selection]
    ndates_arima = ndates_arima[selection]
    min_flow_arima = min_flow_arima[selection]
    med_flow_arima = med_flow_arima[selection]
    max_flow_arima = max_flow_arima[selection]
    min_occupancy_arima = min_occupancy_arima[selection]
    med_occupancy_arima = med_occupancy_arima[selection]
    max_occupancy_arima = max_occupancy_arima[selection]
    min_speed_arima = min_speed_arima[selection]
    med_speed_arima = med_speed_arima[selection]
    max_speed_arima = max_speed_arima[selection]

    # Write out the summary statistics for the sets of loop detector measurements (raw and ARIMA) for the current
    # city and its accepted loop detectors
    with open(summary_stats_of_measurements_per_city_detector_file, mode = 'a') as ssfile:
        for i in range(nld_locations):
            line_out = country_name + ' ' + city_name + ' ' + ld_locations_table['DETECTOR_ID'][i] + ' ' + str(nld_measurements_raw_all[i]) + ' ' + str(nld_measurements_raw_good[i]) + ' ' + \
                       first_time_stamp_raw[i] + ' ' + last_time_stamp_raw[i] + ' ' + str(ndates_raw[i]) + ' ' + interval_lengths_raw[i] + ' ' + str(med_interval_length_raw[i]) + ' ' + \
                       '{:25.2f}'.format(min_flow_raw[i]).strip() + ' ' + '{:25.2f}'.format(med_flow_raw[i]).strip() + ' ' + '{:25.2f}'.format(max_flow_raw[i]).strip() + ' ' + \
                       '{:25.4f}'.format(min_occupancy_raw[i]).strip() + ' ' + '{:25.4f}'.format(med_occupancy_raw[i]).strip() + ' ' + '{:25.4f}'.format(max_occupancy_raw[i]).strip() + ' ' + \
                       '{:25.3f}'.format(min_speed_raw[i]).strip() + ' ' + '{:25.3f}'.format(med_speed_raw[i]).strip() + ' ' + '{:25.3f}'.format(max_speed_raw[i]).strip() + ' ' + \
                       str(nld_measurements_arima_all[i]) + ' ' + str(nld_measurements_arima_good[i]) + ' ' + first_time_stamp_arima[i] + ' ' + last_time_stamp_arima[i] + ' ' + str(ndates_arima[i]) + ' ' + \
                       '{:25.2f}'.format(min_flow_arima[i]).strip() + ' ' + '{:25.2f}'.format(med_flow_arima[i]).strip() + ' ' + '{:25.2f}'.format(max_flow_arima[i]).strip() + ' ' + \
                       '{:25.4f}'.format(min_occupancy_arima[i]).strip() + ' ' + '{:25.4f}'.format(med_occupancy_arima[i]).strip() + ' ' + '{:25.4f}'.format(max_occupancy_arima[i]).strip() + ' ' + \
                       '{:25.3f}'.format(min_speed_arima[i]).strip() + ' ' + '{:25.3f}'.format(med_speed_arima[i]).strip() + ' ' + '{:25.3f}'.format(max_speed_arima[i]).strip() + '\n'
            ssfile.write(line_out)

    # Compute summary statistics for the set of accepted loop detectors for the current city
    min_longitude = '{:25.7f}'.format(numpy.min(ld_locations_table['LONGITUDE'])).strip()
    max_longitude = '{:25.7f}'.format(numpy.max(ld_locations_table['LONGITUDE'])).strip()
    min_latitude = '{:25.7f}'.format(numpy.min(ld_locations_table['LATITUDE'])).strip()
    max_latitude = '{:25.7f}'.format(numpy.max(ld_locations_table['LATITUDE'])).strip()
    min_length = '{:25.3f}'.format(numpy.min(ld_locations_table['LENGTH'])).strip()
    med_length = '{:25.3f}'.format(numpy.median(ld_locations_table['LENGTH'])).strip()
    max_length = '{:25.3f}'.format(numpy.max(ld_locations_table['LENGTH'])).strip()
    min_position = '{:25.4f}'.format(numpy.min(ld_locations_table['POSITION'])).strip()
    med_position = '{:25.4f}'.format(numpy.median(ld_locations_table['POSITION'])).strip()
    max_position = '{:25.4f}'.format(numpy.max(ld_locations_table['POSITION'])).strip()
    n_living_street = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'living_street'))
    n_motorway = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'motorway'))
    n_motorway_link = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'motorway_link'))
    n_primary = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] ==  'primary'))
    n_primary_link = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'primary_link'))
    n_residential = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'residential'))
    n_secondary = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'secondary'))
    n_secondary_link = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'secondary_link'))
    n_service = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'service'))
    n_tertiary = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'tertiary'))
    n_tertiary_link = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'tertiary_link'))
    n_trunk = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'trunk'))
    n_trunk_link = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'trunk_link'))
    n_unclassified = str(numpy.count_nonzero(ld_locations_table['ROAD_CLASS'] == 'unclassified'))
    selection = ld_locations_table['SPEED_LIMIT'] > 0.0
    if numpy.count_nonzero(selection) > 0:
        min_speed_limit = '{:25.2f}'.format(numpy.min(ld_locations_table['SPEED_LIMIT'][selection])).strip()
        med_speed_limit = '{:25.2f}'.format(numpy.median(ld_locations_table['SPEED_LIMIT'][selection])).strip()
        max_speed_limit = '{:25.2f}'.format(numpy.max(ld_locations_table['SPEED_LIMIT'][selection])).strip()
    else:
        min_speed_limit = '-1.00'
        med_speed_limit = '-1.00'
        max_speed_limit = '-1.00'
    min_nlanes = str(numpy.min(ld_locations_table['NLANES']))
    max_nlanes = str(numpy.max(ld_locations_table['NLANES']))
    selection = ld_locations_table['LINK_PTS_FLAG'] == 1
    if numpy.count_nonzero(selection) > 0:
        min_link_pts_longitude = '{:25.7f}'.format(numpy.min(ld_locations_table['LINK_PTS_LONGITUDE'][selection])).strip()
        max_link_pts_longitude = '{:25.7f}'.format(numpy.max(ld_locations_table['LINK_PTS_LONGITUDE'][selection])).strip()
        min_link_pts_latitude = '{:25.7f}'.format(numpy.min(ld_locations_table['LINK_PTS_LATITUDE'][selection])).strip()
        max_link_pts_latitude = '{:25.7f}'.format(numpy.max(ld_locations_table['LINK_PTS_LATITUDE'][selection])).strip()
    else:
        min_link_pts_longitude = '0.0000000'
        max_link_pts_longitude = '0.0000000'
        min_link_pts_latitude = '0.0000000'
        max_link_pts_latitude = '0.0000000'

    # Write out the summary statistics for the set of accepted loop detectors for the current city
    line_out = country_name + ' ' + city_name + ' ' + str(nld_locations) + ' ' + min_longitude + ' ' + max_longitude + ' ' + min_latitude + ' ' + max_latitude + ' ' + \
               min_length + ' ' + med_length + ' ' + max_length + ' ' + min_position + ' ' + med_position + ' ' + max_position + ' ' + n_living_street + ' ' + \
               n_motorway + ' ' + n_motorway_link + ' ' + n_primary + ' ' + n_primary_link + ' ' + n_residential + ' ' + n_secondary + ' ' + n_secondary_link + ' ' + \
               n_service + ' ' + n_tertiary + ' ' + n_tertiary_link + ' ' + n_trunk + ' ' + n_trunk_link + ' ' + n_unclassified + ' ' + min_speed_limit + ' ' + \
               med_speed_limit + ' ' + max_speed_limit + ' ' + min_nlanes + ' ' + max_nlanes + ' ' + min_link_pts_longitude + ' ' + max_link_pts_longitude + ' ' + \
               min_link_pts_latitude + ' ' + max_link_pts_latitude + '\n'
    with open(summary_stats_of_detectors_per_city_file, mode = 'a') as ssfile:
        ssfile.write(line_out)

    # Write out the filtered loop detector locations data table for the current city
    curr_output_file = os.path.join(output_dir_ld_locations, 'detectors.' + data_source + '.' + country_name + '.' + city_name + '.fits')
    print('Writing out the filtered loop detector locations data table: ' + curr_output_file)
    ld_locations_table.write(curr_output_file, format = 'fits')

# Finish
print('')
print('>-------------------------------------------------------------------------------------------------------------<')
print('')
print('Finished!')
