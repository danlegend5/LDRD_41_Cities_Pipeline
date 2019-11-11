#!/usr/bin/env python

__author__ = 'Dan Bramich'

# This script reads in the data on the loop detectors, the road links in which they are installed, and the
# traffic measurements that they have made, and it splits these data into FITS binary table files by city
# (and detector) while also standardising the entries. Note that no filtering of the data is performed
# since the purpose of this script is simply to better organise and standardise the data before the data
# processing and analysis begins. The data come from the publication "Understanding traffic capacity of
# urban networks" by Loder et al. (2019), except for the extra one year of data on Zurich covering 2018.

# Imports
import csv
import glob
import numpy
import os
import shutil
from astropy.table import Column
from astropy.table import Table
from LDDA_MFD_Project1_Pipeline.config import config
from LDDA_MFD_Project1_Pipeline.lib import general_functions

# Prepare a table for the loop detector locations data
print('')
print('Preparing a table for the loop detector locations data...')
nld_locations = general_functions.count_file_lines(config.original_detectors_file) - 1
empty_str = '                                                                                '
ld_locations_table = Table([Column(data = [empty_str]*nld_locations, name = 'DETECTOR_ID'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.float64), name = 'LONGITUDE'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.float64), name = 'LATITUDE'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.float64), name = 'LENGTH'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.float64), name = 'POSITION'),
                            Column(data = [empty_str]*nld_locations, name = 'ROAD_NAME'),
                            Column(data = [empty_str]*nld_locations, name = 'ROAD_CLASS'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.float64), name = 'SPEED_LIMIT'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.int32), name = 'NLANES'),
                            Column(data = numpy.zeros(nld_locations, dtype = numpy.int32), name = 'LINK_ID'),
                            Column(data = [numpy.zeros(100, dtype = numpy.float64)]*nld_locations, name = 'LINK_PTS_LONGITUDE'),
                            Column(data = [numpy.zeros(100, dtype = numpy.float64)]*nld_locations, name = 'LINK_PTS_LATITUDE'),
                            Column(data = [numpy.zeros(100, dtype = numpy.int32)]*nld_locations, name = 'LINK_PTS_FLAG'),
                            Column(data = [empty_str]*nld_locations, name = 'CITY_NAME')])

# Read in the loop detector locations data file
print('Reading in the loop detector locations data file: ' + config.original_detectors_file)
with open(config.original_detectors_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    i = 0
    for row in csv_contents:
        tmp_str = row['detid'].replace(' ', '_')
        tmp_str = tmp_str.replace('.', '_')
        tmp_str = tmp_str.replace('/', '_')
        tmp_str = tmp_str.replace('[', '_')
        tmp_str = tmp_str.replace(']', '_')
        tmp_str = tmp_str.replace('(', '_')
        tmp_str = tmp_str.replace(')', '_')
        if tmp_str[0] == '_': tmp_str = tmp_str[1:]
        if tmp_str[-1] == '_': tmp_str = tmp_str[0:(len(tmp_str) - 1)]
        ld_locations_table['DETECTOR_ID'][i] = tmp_str                               # Notes: - The ID of the loop detector.
                                                                                     #        - Spaces replaced with underscores.
                                                                                     #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
                                                                                     #        - Underscores at the beginning and end of the ID string are removed.
                                                                                     #        - The only characters that are present in the loop detector ID entries are:
                                                                                     #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        ld_locations_table['LONGITUDE'][i] = numpy.float64(row['coords.x1'])         # Notes: - Longitude (deg; EPSG:4326 or WGS84 - World Geodetic System 1984) of the loop detector.
                                                                                     #        - Longitudes are quoted in range -180.0 to 180.0 degrees.
                                                                                     #        - Minimum: -118.303722 deg
                                                                                     #        - Maximum: 145.08013040207 deg
        ld_locations_table['LATITUDE'][i] = numpy.float64(row['coords.x2'])          # Notes: - Latitude (deg; EPSG:4326 or WGS84 - World Geodetic System 1984) of the loop detector.
                                                                                     #        - Latitudes are quoted in range -90.0 to 90.0 degrees.
                                                                                     #        - Minimum: -37.8745292808731 deg
                                                                                     #        - Maximum: 54.70488 deg
        ld_locations_table['LENGTH'][i] = numpy.float64(row['length'])               # Notes: - Length (km) of the associated link road (i.e the road where the detector is located).
                                                                                     #        - Minimum: 0.007199797164467 km
                                                                                     #        - Maximum: 16.6844474508417 km
        if row['pos'] == 'NA':                                                       # Notes: - Loop location as a distance (km) along the associated link road from the downstream
            ld_locations_table['POSITION'][i] = numpy.float64(-1.0)                  #          intersection.
        else:                                                                        #        - 45 out of 1042 entries for "utrecht" are equal to 'NA'. These are replaced with the
            ld_locations_table['POSITION'][i] = numpy.float64(row['pos'])            #          the value '-1.0'.
                                                                                     #        - Minimum: 0.0 km
                                                                                     #        - Maximum: 6.67473242607813 km
        ld_locations_table['ROAD_NAME'][i] = row['road'].replace(' ', '_')           # Notes: - Name of the associated link road.
                                                                                     #        - Spaces replaced with underscores.
                                                                                     #        - 1544 entries are equal to 'NA'.
                                                                                     #        - Many strings are unintelligible.
        ld_locations_table['ROAD_CLASS'][i] = row['fclass']                          # Notes: - Classification of the associated link road.
                                                                                     #        - Set of possible entries:
                                                                                     #          'cycleway', 'footway', 'living_street', 'motorway', 'motorway_link', 'other',
                                                                                     #          'path', 'pedestrian', 'primary', 'primary_link', 'residential', 'secondary',
                                                                                     #          'secondary_link', 'service', 'tertiary', 'tertiary_link', 'trunk', 'trunk_link',
                                                                                     #          'unclassified'
        if row['limit'] == 'NA':                                                     # Notes: - Speed limit (km/h) of the associated link road as taken from OSM maps (not very
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #          accurate).
        elif row['limit'] == '50|50|30':                                             #        - 6632 entries are equal to 'NA' ("bremen", "london", "losanageles", "madrid",
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #          "melbourne", "paris", "tokyo", "utrecht"). These are replaced with the value
        elif row['limit'] == '60; 40; 60':                                           #          '-1.0'.
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #        - 3 entries corresponding to "Calle de Segovia" in "madrid" are equal to '50|50|30'.
        elif row['limit'] == '40; 60':                                               #          These indicate different speed limits in each lane. For simplicity, these are
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #          replaced with the value '-1.0'.
        else:                                                                        #        - 3 entries corresponding to "La Trobe Street" in "melbourne" are equal to
            if numpy.float64(row['limit']) == 0.0:                                   #          '60; 40; 60'. These indicate different speed limits in each lane. For simplicity,
                ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)           #          these are replaced with the value '-1.0'.
            else:                                                                    #        - 2 entries corresponding to "La Trobe Street" in "melbourne" are equal to '40; 60'.
                ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(row['limit'])   #          These indicate different speed limits in each lane. For simplicity, these are
                                                                                     #          replaced with the value '-1.0'.
                                                                                     #        - 3754 entries are equal to '0.0' ("BORDEAUX", "CONSTANCE", "GRAZ", "SANTANDER",
                                                                                     #          "augsburg", "basel", "bern", "birmingham", "bolton", "bremen", "cagliari",
                                                                                     #          "darmstadt", "duisburg", "essen", "frankfurt", "hamburg", "innsbruck", "kassel",
                                                                                     #          "luzern", "manchester", "marseille", "rotterdam", "speyer", "strasbourg",
                                                                                     #          "stuttgart", "taipeh", "tokyo", "torino", "toronto", "toulouse", "utrecht",
                                                                                     #          "vilnius", "wolfsburg", "zurich"). These are replaced with the value '-1.0'.
                                                                                     #        - Minimum: 17.6 km/h
                                                                                     #        - Maximum: 100.0 km/h
        if row['lanes'] == 'NA':                                                     # Notes: - Number of lanes covered by a loop detector. Typically 1, but in some cases this
            ld_locations_table['NLANES'][i] = numpy.int32(-1)                        #          is more than 1.
        else:                                                                        #        - 4 out of 438 entries for "BORDEAUX" are equal to 'NA'. These are replaced with
            ld_locations_table['NLANES'][i] = numpy.int32(row['lanes'])              #          the value '-1'.
                                                                                     #        - Minimum: 1
                                                                                     #        - Maximum: 9
        if row['linkid'] == 'NA':                                                    # Notes: - The ID of the associated link road as used in the file with the link road data.
            ld_locations_table['LINK_ID'][i] = numpy.int32(-1)                       #        - 605 entries are equal to 'NA' ("BORDEAUX", "CONSTANCE", "utrecht"). These are
        else:                                                                        #          replaced with the value '-1'.
            ld_locations_table['LINK_ID'][i] = numpy.int32(row['linkid'])            #        - Minimum: 0
                                                                                     #        - Maximum: 5268
        if row['citycode'] == 'BORDEAUX':                                            # Notes: - Name of the city the loop detector belongs to.
            ld_locations_table['CITY_NAME'][i] = 'bordeaux'                          #        - Four cities have strings in upper case ("BORDEAUX", "CONSTANCE", "GRAZ",
        elif row['citycode'] == 'CONSTANCE':                                         #          "SANTANDER"). These are converted to lower case.
            ld_locations_table['CITY_NAME'][i] = 'constance'                         #        - Los Angeles is spelt wrong "losanageles". This is corrected.
        elif row['citycode'] == 'GRAZ':                                              #        - 41 cities with names:
            ld_locations_table['CITY_NAME'][i] = 'graz'                              #          'augsburg', 'basel', 'bern', 'birmingham', 'bolton', 'bordeaux', 'bremen',
        elif row['citycode'] == 'SANTANDER':                                         #          'cagliari', 'constance', 'darmstadt', 'duisburg', 'essen', 'frankfurt', 'graz',
            ld_locations_table['CITY_NAME'][i] = 'santander'                         #          'groningen', 'hamburg', 'innsbruck', 'kassel', 'london', 'losangeles', 'luzern',
        elif row['citycode'] == 'losanageles':                                       #          'madrid', 'manchester', 'marseille', 'melbourne', 'munich', 'paris', 'rotterdam',
            ld_locations_table['CITY_NAME'][i] = 'losangeles'                        #          'santander', 'speyer', 'strasbourg', 'stuttgart', 'taipeh', 'tokyo', 'torino',
        else:                                                                        #          'toronto', 'toulouse', 'utrecht', 'vilnius', 'wolfsburg', 'zurich'
            ld_locations_table['CITY_NAME'][i] = row['citycode']
        i += 1
print('Read in ' + str(nld_locations) + ' rows...')

# N.B: There is one duplicated entry for "CITY_NAME = toulouse" and "DETECTOR_ID = 262" in the loop detector
# locations table. The remaining columns have different entries. Due to the ambiguity, these two entries will
# be filtered out in the next stage of the pipeline.

# Determine the set of unique city names
print('Determining the set of unique city names...')
city_names_uniq = [city for city in set(ld_locations_table['CITY_NAME'])]
city_names_uniq.sort()
ncities = len(city_names_uniq)
print('No. of unique city names: ' + str(ncities))

# Prepare a table for the loop detector links data
print('')
print('Preparing a table for the loop detector links data...')
nld_links = general_functions.count_file_lines(config.original_links_file) - 1
ld_links_table = Table([Column(data = [empty_str]*nld_links, name = 'CITY_NAME'),
                        Column(data = numpy.zeros(nld_links, dtype = numpy.int32), name = 'LINK_ID'),
                        Column(data = numpy.zeros(nld_links, dtype = numpy.int32), name = 'ORDER'),
                        Column(data = numpy.zeros(nld_links, dtype = numpy.float64), name = 'LONGITUDE'),
                        Column(data = numpy.zeros(nld_links, dtype = numpy.float64), name = 'LATITUDE')])

# Read in the loop detector links data file
print('Reading in the loop detector links data file: ' + config.original_links_file)
with open(config.original_links_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    i = 0
    for row in csv_contents:
        ld_links_table['CITY_NAME'][i] = row['city']                                 # Notes: - Name of the city the corresponding link road belongs to.
                                                                                     #        - 41 cities with the same names as above.
        ld_links_table['LINK_ID'][i] = numpy.int32(row['id'])                        # Notes: - The ID of the link road to which this entry corresponds.
                                                                                     #        - Minimum: 0
                                                                                     #        - Maximum: 5223
        ld_links_table['ORDER'][i] = numpy.int32(row['order'])                       # Notes: - The index of this point in the set of points making up the corresponding link road.
                                                                                     #        - Minimum: 1
                                                                                     #        - Maximum: 87
        ld_links_table['LONGITUDE'][i] = numpy.float64(row['long'])                  # Notes: - The longitude (deg; EPSG:4326 or WGS84 - World Geodetic System 1984) of this point.
                                                                                     #        - Longitudes are quoted in range -180.0 to 180.0 degrees.
                                                                                     #        - Minimum: -118.3055 deg
                                                                                     #        - Maximum: 144.9945 deg
        ld_links_table['LATITUDE'][i] = numpy.float64(row['lat'])                    # Notes: - The latitude (deg; EPSG:4326 or WGS84 - World Geodetic System 1984) of this point.
                                                                                     #        - Latitudes are quoted in range -90.0 to 90.0 degrees.
                                                                                     #        - Minimum: -37.82935 deg
                                                                                     #        - Maximum: 54.71058 deg
        i += 1
print('Read in ' + str(nld_links) + ' rows...')

# N.B: There are no duplicated entries for "CITY_NAME + LINK_ID + ORDER" in the loop detector links table.

# Merge the loop detector links data table into the loop detector locations data table
print('')
print('Merging the loop detector links data table into the loop detector locations data table...')
for i in range(nld_links):
    curr_city = ld_links_table['CITY_NAME'][i]
    curr_link_id = ld_links_table['LINK_ID'][i]
    selection = numpy.logical_and(ld_locations_table['CITY_NAME'] == curr_city, ld_locations_table['LINK_ID'] == curr_link_id)
    nassociated_loop_detectors = numpy.count_nonzero(selection)
    if nassociated_loop_detectors == 0: continue
    curr_order = ld_links_table['ORDER'][i] - 1
    curr_longitude = ld_links_table['LONGITUDE'][i]
    curr_latitude = ld_links_table['LATITUDE'][i]
    detector_subs = numpy.argwhere(selection).flatten()
    for j in range(len(detector_subs)):
        csub = detector_subs[j]
        if ld_locations_table['LINK_PTS_FLAG'][csub][curr_order] == 1:
            print('ERROR - Duplicated entries exist in the loop detector links data file!')
            exit()
        ld_locations_table['LINK_PTS_LONGITUDE'][csub][curr_order] = curr_longitude   # Notes: - The longitude (deg; EPSG:4326 or WGS84 - World Geodetic System 1984) of this point
                                                                                      #          on the link road associated with this loop detector.
        ld_locations_table['LINK_PTS_LATITUDE'][csub][curr_order] = curr_latitude     # Notes: - The latitude (deg; EPSG:4326 or WGS84 - World Geodetic System 1984) of this point
                                                                                      #          on the link road associated with this loop detector.
        ld_locations_table['LINK_PTS_FLAG'][csub][curr_order] = 1                     # Notes: - Flag indicating if this point is to be used in defining the link road associated
                                                                                      #          with this loop detector (0 = No; 1 = Yes).

# Create the output directory for the loop detector locations data tables
output_dir_ld_locations = os.path.join(config.output_dir, 's0.Loop.Detector.Locations')
print('')
print('Creating the output directory for the loop detector locations data tables: ' + output_dir_ld_locations)
if os.path.exists(output_dir_ld_locations): shutil.rmtree(output_dir_ld_locations)
os.makedirs(output_dir_ld_locations)

# For each city
for city_name in city_names_uniq:

    # Write out a loop detector locations data table for the current city
    print('Writing out a loop detector locations data table for: ' + city_name)
    curr_ld_locations_table = ld_locations_table[ld_locations_table['CITY_NAME'] == city_name]
    curr_ld_locations_table.remove_columns(['ROAD_NAME', 'CITY_NAME'])
    curr_ld_locations_table = curr_ld_locations_table[numpy.argsort(curr_ld_locations_table['LATITUDE'])]
    country_name = general_functions.get_country_name(city_name)
    curr_output_file = os.path.join(output_dir_ld_locations, 'detectors.' + country_name + '.' + city_name + '.fits')
    curr_ld_locations_table.write(curr_output_file, format = 'fits')

# Create the output directory for the loop detector measurements data tables (raw)
output_dir_ld_measurements_raw = os.path.join(config.output_dir, 's0.Loop.Detector.Measurements.Raw')
print('')
print('Creating the output directory for the loop detector measurements data tables (raw): ' + output_dir_ld_measurements_raw)
if os.path.exists(output_dir_ld_measurements_raw): shutil.rmtree(output_dir_ld_measurements_raw)
os.makedirs(output_dir_ld_measurements_raw)

# Read in the loop detector measurements data file (raw), and write out the entries to data files organised by
# city and detector ID
print('Reading in and splitting the loop detector measurements data file (raw): ' + config.original_measurements_raw_file)
fieldnames = ['day', 'interval', 'flow', 'occ', 'error', 'speed']
with open(config.original_measurements_raw_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        city_name = row['city']                                                      # Notes: - Name of the city that hosts the loop detector the measurement was taken with.
                                                                                     #        - 41 cities with the same names as above.
        if city_name == 'innsbruck': continue
        if city_name == 'zurich': continue
        country_name = general_functions.get_country_name(city_name)
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        curr_output_dir = os.path.join(output_dir_ld_measurements_raw, country_name, city_name, detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.raw.' + country_name + '.' + city_name + '.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)

# Read in the loop detector measurements data file for Innsbruck, and write out the entries to data files
# organised by city and detector ID
print('Reading in and splitting the loop detector measurements data file (raw): ' + config.original_measurements_innsbruck_file)
country_name = general_functions.get_country_name('innsbruck')
with open(config.original_measurements_innsbruck_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        row['day'] = row['day'][6:10] + '-' + row['day'][3:5] + '-' + row['day'][0:2]
        curr_output_dir = os.path.join(output_dir_ld_measurements_raw, country_name, 'innsbruck', detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.raw.' + country_name + '.innsbruck.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)

# Read in the loop detector measurements data files for Zurich, and write out the entries to data files
# organised by city and detector ID
print('Reading in and splitting the loop detector measurements data file (raw): ' + config.original_measurements_zurich_file1)
country_name = general_functions.get_country_name('zurich')
with open(config.original_measurements_zurich_file1, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        row['day'] = row['day'][6:10] + '-' + row['day'][3:5] + '-' + row['day'][0:2]
        row['speed'] = 'NA'
        curr_output_dir = os.path.join(output_dir_ld_measurements_raw, country_name, 'zurich', detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.raw.' + country_name + '.zurich.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)
print('Reading in and splitting the loop detector measurements data file (raw): ' + config.original_measurements_zurich_file2)
with open(config.original_measurements_zurich_file2, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        row['day'] = row['day'][6:10] + '-' + row['day'][3:5] + '-' + row['day'][0:2]
        row['speed'] = 'NA'
        curr_output_dir = os.path.join(output_dir_ld_measurements_raw, country_name, 'zurich', detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.raw.' + country_name + '.zurich.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)

# Convert the loop detector measurements data files (raw) from ".csv" files to ".fits" files, while also
# removing duplicated rows
file_list = glob.glob(os.path.join(output_dir_ld_measurements_raw, '*', '*', '*', 'measurements.raw.*.csv'))
for file in file_list:

    # Determine the output file name
    path, basename_with_ext = os.path.split(file)
    basename_bits = basename_with_ext.split('.')
    print('Writing out a loop detector measurements data table (raw) for: ' + basename_bits[3] + ',' + basename_bits[4])
    curr_output_file = os.path.join(path, basename_bits[0] + '.' + basename_bits[1] + '.' + basename_bits[2] + '.' + basename_bits[3] + '.' + basename_bits[4] + '.fits')

    # Read in the current loop detector measurements ".csv" data file
    with open(file, mode = 'r') as csv_file:
        lines = [line.rstrip('\n') for line in csv_file if line.rstrip('\n')]
    lines = lines[1:]

    # Remove any duplicated rows
    lines_uniq = [line for line in set(lines)]
    nld_measurements = len(lines_uniq)

    # Prepare a table for the loop detector measurements data
    ld_measurements_table = Table([Column(data = [empty_str]*nld_measurements, name = 'DATE'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.int32), name = 'INTERVAL_START'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'FLOW'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'OCCUPANCY'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.int32), name = 'ERROR_FLAG'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'SPEED')])

    # Ingest the loop detector measurements data into the table
    for i, line in enumerate(lines_uniq):
        line_bits = line.split(',')
        if line_bits[0] == 'NA':                                                                 # Notes: - Date on which the measurement was taken (YYYY-MM-DD local time).
            ld_measurements_table['DATE'][i] = '0000-00-00'                                      #        - Entries equal to 'NA' are replaced with the value '0000-00-00'.
        else:                                                                                    #        - Minimum: '2008-05-16'
            ld_measurements_table['DATE'][i] = line_bits[0]                                      #        - Maximum: '2018-12-31'
        if line_bits[1] == 'NA':                                                                 # Notes: - Time at the start of the measurement interval (seconds after midnight
            ld_measurements_table['INTERVAL_START'][i] = numpy.int32(-1)                         #          on the corresponding date - local time).
        else:                                                                                    #        - Entries equal to 'NA' are replaced with the value '-1'.
            if numpy.int32(line_bits[1]) > 86400:                                                #        - Minimum: 0
                ld_measurements_table['INTERVAL_START'][i] = numpy.int32(-1)                     #        - Maximum: 86400
            else:
                ld_measurements_table['INTERVAL_START'][i] = numpy.int32(line_bits[1])
        if line_bits[2] == 'NA':                                                                 # Notes: - Flow measurement (veh/hour). This is the vehicle count in the
            ld_measurements_table['FLOW'][i] = numpy.float64(-1.0)                               #          measurement interval, scaled to 1 hour.
        else:                                                                                    #        - Entries equal to 'NA' are replaced with the value '-1.0'.
            if numpy.float64(line_bits[2]) < 0.0:                                                #        - Entries with negative values are replaced with the value '-1.0'.
                ld_measurements_table['FLOW'][i] = numpy.float64(-1.0)
            else:
                ld_measurements_table['FLOW'][i] = numpy.float64(line_bits[2])
        if line_bits[3] == 'Inf':                                                                # Notes: - Fraction of time in the measurement interval that the loop detector
            ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                          #          is occupied by a vehicle.
        elif line_bits[3] == 'NA':                                                               #        - Entries equal to 'Inf' are replaced with the value '-1.0'.
            ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                          #        - Entries equal to 'NA' are replaced with the value '-1.0'.
        else:                                                                                    #        - Entries with negative values are replaced with the value '-1.0'.
            if numpy.float64(line_bits[3]) < 0.0:                                                #        - Entries greater than unity are replaced with the value '-1.0'.
                ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                      #        - Minimum: 0.0
            elif numpy.float64(line_bits[3]) > 1.0:                                              #        - Maximum: 1.0
                ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                      #        - Q: How is vehicle length accounted/corrected for?
            else:                                                                                #          A: It's not done directly. It is calibrated via the free flow branch
                ld_measurements_table['OCCUPANCY'][i] = numpy.float64(line_bits[3])              #             of the FD and the MFD.
        if line_bits[4] == 'NA':                                                                 # Notes: - Flag indicating an error with the measurement (0 = No error; 1 = Error)
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(0)                              #        - Entries equal to 'NA' are replaced with the value '0'.
        elif line_bits[4] == '2':                                                                #        - Entries equal to '2' are replaced with the value '1'.
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(1)
        else:
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(line_bits[4])
        if line_bits[5] == 'NA':                                                                 # Notes: - Average vehicle speed in the measurement interval (??UNITS??).
            ld_measurements_table['SPEED'][i] = numpy.float64(-1.0)                              #        - Some cities provide average speed measurements instead of occupancy.
        else:                                                                                    #          Otherwise, average speeds are calculated from speed = flow / density.
            if numpy.float64(line_bits[5]) < 0.0:                                                #        - Entries equal to 'NA' are replaced with the value '-1.0'.
                ld_measurements_table['SPEED'][i] = numpy.float64(-1.0)                          #        - Entries with negative values are replaced with the value '-1.0'.
            else:
                ld_measurements_table['SPEED'][i] = numpy.float64(line_bits[5])

    # Sort the loop detector measurements data table by "DATE" followed by "INTERVAL_START"
    ld_measurements_table = ld_measurements_table[numpy.argsort(ld_measurements_table, order = ['DATE', 'INTERVAL_START'])]

    # Set "ERROR_FLAG" to "1" for any duplicated entries for "DATE + INTERVAL_START" in this data table
    for i in range(nld_measurements - 1):
        curr_date_interval = ld_measurements_table['DATE'][i] + '_' + str(ld_measurements_table['INTERVAL_START'][i])
        next_date_interval = ld_measurements_table['DATE'][i + 1] + '_' + str(ld_measurements_table['INTERVAL_START'][i + 1])
        if next_date_interval == curr_date_interval:
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(1)
            ld_measurements_table['ERROR_FLAG'][i + 1] = numpy.int32(1)

    # Write out the loop detector measurements data table for the current city and detector ID
    ld_measurements_table.write(curr_output_file, format = 'fits')
    os.remove(file)

# N.B: There are some duplicated entries for "DATE + INTERVAL_START" in the loop detector measurements tables.
# These entries have had "ERROR_FLAG" set to "1".

# Create the output directory for the loop detector measurements data tables (ARIMA)
output_dir_ld_measurements_arima = os.path.join(config.output_dir, 's0.Loop.Detector.Measurements.ARIMA')
print('')
print('Creating the output directory for the loop detector measurements data tables (ARIMA): ' + output_dir_ld_measurements_arima)
if os.path.exists(output_dir_ld_measurements_arima): shutil.rmtree(output_dir_ld_measurements_arima)
os.makedirs(output_dir_ld_measurements_arima)

# Read in the loop detector measurements data file (ARIMA), and write out the entries to data files organised
# by city and detector ID
print('Reading in and splitting the loop detector measurements data file (ARIMA): ' + config.original_measurements_arima_file)
fieldnames = ['day', 'interval', 'flow', 'occ', 'error', 'speed', 'arima.flow', 'arima.occ', 'arima.speed']
with open(config.original_measurements_arima_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        city_name = row['city']                                                      # Notes: - Name of the city that hosts the loop detector the measurement was taken with.
                                                                                     #        - 41 cities with the same names as above.
        if city_name == 'innsbruck': continue
        if city_name == 'zurich': continue
        country_name = general_functions.get_country_name(city_name)
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        curr_output_dir = os.path.join(output_dir_ld_measurements_arima, country_name, city_name, detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.ARIMA.' + country_name + '.' + city_name + '.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)

# Read in the loop detector measurements data file for Innsbruck, and write out the entries to data files
# organised by city and detector ID
print('Reading in and splitting the loop detector measurements data file (ARIMA): ' + config.original_measurements_innsbruck_file)
country_name = general_functions.get_country_name('innsbruck')
with open(config.original_measurements_innsbruck_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        row['day'] = row['day'][6:10] + '-' + row['day'][3:5] + '-' + row['day'][0:2]
        row['arima.occ'] = 'NA'
        curr_output_dir = os.path.join(output_dir_ld_measurements_arima, country_name, 'innsbruck', detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.ARIMA.' + country_name + '.innsbruck.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)

# Read in the loop detector measurements data files for Zurich, and write out the entries to data files
# organised by city and detector ID
print('Reading in and splitting the loop detector measurements data file (ARIMA): ' + config.original_measurements_zurich_file1)
country_name = general_functions.get_country_name('zurich')
with open(config.original_measurements_zurich_file1, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        row['day'] = row['day'][6:10] + '-' + row['day'][3:5] + '-' + row['day'][0:2]
        row['speed'] = 'NA'
        row['arima.speed'] = 'NA'
        curr_output_dir = os.path.join(output_dir_ld_measurements_arima, country_name, 'zurich', detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.ARIMA.' + country_name + '.zurich.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)
print('Reading in and splitting the loop detector measurements data file (ARIMA): ' + config.original_measurements_zurich_file2)
with open(config.original_measurements_zurich_file2, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:
        detector_id = row['detid'].replace(' ', '_')                                 # Notes: - The ID of the loop detector the measurement was taken with.
        detector_id = detector_id.replace('.', '_')                                  #        - Spaces replaced with underscores.
        detector_id = detector_id.replace('/', '_')                                  #        - Characters '.', '/', '[', ']', '(' and ')' replaced with underscores.
        detector_id = detector_id.replace('[', '_')                                  #        - Underscores at the beginning and end of the ID string are removed.
        detector_id = detector_id.replace(']', '_')                                  #        - The only characters that are present in the loop detector ID entries are:
        detector_id = detector_id.replace('(', '_')                                  #          '_', '-', '+', '0', ..., '9', 'a', ..., 'z', 'A', ..., 'Z'
        detector_id = detector_id.replace(')', '_')
        if detector_id[0] == '_': detector_id = detector_id[1:]
        if detector_id[-1] == '_': detector_id = detector_id[0:(len(detector_id) - 1)]
        row['day'] = row['day'][6:10] + '-' + row['day'][3:5] + '-' + row['day'][0:2]
        row['speed'] = 'NA'
        row['arima.speed'] = 'NA'
        curr_output_dir = os.path.join(output_dir_ld_measurements_arima, country_name, 'zurich', detector_id)
        curr_output_file = os.path.join(curr_output_dir, 'measurements.ARIMA.' + country_name + '.zurich.' + detector_id + '.csv')
        if not os.path.exists(curr_output_file):
            os.makedirs(curr_output_dir)
            with open(curr_output_file, mode = 'w') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writeheader()
                writer.writerow(row)
        else:
            with open(curr_output_file, mode = 'a') as curr_csv_file:
                writer = csv.DictWriter(curr_csv_file, fieldnames = fieldnames, extrasaction = 'ignore')
                writer.writerow(row)

# Convert the loop detector measurements data files (ARIMA) from ".csv" files to ".fits" files, while also
# removing duplicated rows
file_list = glob.glob(os.path.join(output_dir_ld_measurements_arima, '*', '*', '*', 'measurements.ARIMA.*.csv'))
for file in file_list:

    # Determine the output file name
    path, basename_with_ext = os.path.split(file)
    basename_bits = basename_with_ext.split('.')
    print('Writing out a loop detector measurements data table (ARIMA) for: ' + basename_bits[3] + ',' + basename_bits[4])
    curr_output_file = os.path.join(path, basename_bits[0] + '.' + basename_bits[1] + '.' + basename_bits[2] + '.' + basename_bits[3] + '.' + basename_bits[4] + '.fits')

    # Read in the current loop detector measurements ".csv" data file
    with open(file, mode = 'r') as csv_file:
        lines = [line.rstrip('\n') for line in csv_file if line.rstrip('\n')]
    lines = lines[1:]

    # Remove any duplicated rows
    lines_uniq = [line for line in set(lines)]
    nld_measurements = len(lines_uniq)

    # Prepare a table for the loop detector measurements data
    ld_measurements_table = Table([Column(data = [empty_str]*nld_measurements, name = 'DATE'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.int32), name = 'INTERVAL_START'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'FLOW'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'OCCUPANCY'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.int32), name = 'ERROR_FLAG'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'SPEED'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'ARIMA_FLOW'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'ARIMA_OCCUPANCY'),
                                   Column(data = numpy.zeros(nld_measurements, dtype = numpy.float64), name = 'ARIMA_SPEED')])

    # Ingest the loop detector measurements data into the table
    for i, line in enumerate(lines_uniq):
        line_bits = line.split(',')
        if line_bits[0] == 'NA':                                                                 # Notes: - Date on which the measurement was taken (YYYY-MM-DD local time).
            ld_measurements_table['DATE'][i] = '0000-00-00'                                      #        - Entries equal to 'NA' are replaced with the value '0000-00-00'.
        else:                                                                                    #        - Minimum: '2008-05-16'
            ld_measurements_table['DATE'][i] = line_bits[0]                                      #        - Maximum: '2018-12-31'
        if line_bits[1] == 'NA':                                                                 # Notes: - Time at the start of the measurement interval (seconds after midnight
            ld_measurements_table['INTERVAL_START'][i] = numpy.int32(-1)                         #          on the corresponding date - local time).
        else:                                                                                    #        - Entries equal to 'NA' are replaced with the value '-1'.
            if numpy.int32(line_bits[1]) > 86400:                                                #        - Minimum: 0
                ld_measurements_table['INTERVAL_START'][i] = numpy.int32(-1)                     #        - Maximum: 86400
            else:
                ld_measurements_table['INTERVAL_START'][i] = numpy.int32(line_bits[1])
        if line_bits[2] == 'NA':                                                                 # Notes: - Flow measurement (veh/hour). This is the vehicle count in the
            ld_measurements_table['FLOW'][i] = numpy.float64(-1.0)                               #          measurement interval, scaled to 1 hour.
        else:                                                                                    #        - Entries equal to 'NA' are replaced with the value '-1.0'.
            if numpy.float64(line_bits[2]) < 0.0:                                                #        - Entries with negative values are replaced with the value '-1.0'.
                ld_measurements_table['FLOW'][i] = numpy.float64(-1.0)
            else:
                ld_measurements_table['FLOW'][i] = numpy.float64(line_bits[2])
        if line_bits[3] == 'Inf':                                                                # Notes: - Fraction of time in the measurement interval that the loop detector
            ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                          #          is occupied by a vehicle.
        elif line_bits[3] == 'NA':                                                               #        - Entries equal to 'Inf' are replaced with the value '-1.0'.
            ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                          #        - Entries equal to 'NA' are replaced with the value '-1.0'.
        else:                                                                                    #        - Entries with negative values are replaced with the value '-1.0'.
            if numpy.float64(line_bits[3]) < 0.0:                                                #        - Entries greater than unity are replaced with the value '-1.0'.
                ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                      #        - Minimum: 0.0
            elif numpy.float64(line_bits[3]) > 1.0:                                              #        - Maximum: 1.0
                ld_measurements_table['OCCUPANCY'][i] = numpy.float64(-1.0)                      #        - Q: How is vehicle length accounted/corrected for?
            else:                                                                                #          A: It's not done directly. It is calibrated via the free flow branch
                ld_measurements_table['OCCUPANCY'][i] = numpy.float64(line_bits[3])              #             of the FD and the MFD.
        if line_bits[4] == 'NA':                                                                 # Notes: - Flag indicating an error with the measurement (0 = No error; 1 = Error)
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(0)                              #        - Entries equal to 'NA' are replaced with the value '0'.
        elif line_bits[4] == '2':                                                                #        - Entries equal to '2' are replaced with the value '1'.
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(1)
        else:
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(line_bits[4])
        if line_bits[5] == 'NA':                                                                 # Notes: - Average vehicle speed in the measurement interval (??UNITS??).
            ld_measurements_table['SPEED'][i] = numpy.float64(-1.0)                              #        - Some cities provide average speed measurements instead of occupancy.
        else:                                                                                    #          Otherwise, average speeds are calculated from speed = flow / density.
            if numpy.float64(line_bits[5]) < 0.0:                                                #        - Entries equal to 'NA' are replaced with the value '-1.0'.
                ld_measurements_table['SPEED'][i] = numpy.float64(-1.0)                          #        - Entries with negative values are replaced with the value '-1.0'.
            else:
                ld_measurements_table['SPEED'][i] = numpy.float64(line_bits[5])
        if line_bits[6] == 'NA':                                                                 # Notes: - Flow measurement (veh/hour; ARIMA smoothed). This is the vehicle count in
            ld_measurements_table['ARIMA_FLOW'][i] = numpy.float64(-1.0)                         #          the measurement interval, scaled to 1 hour.
        else:                                                                                    #        - Entries equal to 'NA' are replaced with the value '-1.0'.
            if numpy.float64(line_bits[6]) < 0.0:                                                #        - Entries with negative values are replaced with the value '-1.0'.
                ld_measurements_table['ARIMA_FLOW'][i] = numpy.float64(-1.0)
            else:
                ld_measurements_table['ARIMA_FLOW'][i] = numpy.float64(line_bits[6])
        if line_bits[7] == 'Inf':                                                                # Notes: - Fraction of time in the measurement interval that the loop detector
            ld_measurements_table['ARIMA_OCCUPANCY'][i] = numpy.float64(-1.0)                    #          is occupied by a vehicle (ARIMA smoothed).
        elif line_bits[7] == 'NA':                                                               #        - Entries equal to 'Inf' are replaced with the value '-1.0'.
            ld_measurements_table['ARIMA_OCCUPANCY'][i] = numpy.float64(-1.0)                    #        - Entries equal to 'NA' are replaced with the value '-1.0'.
        else:                                                                                    #        - Entries with negative values are replaced with the value '-1.0'.
            if numpy.float64(line_bits[7]) < 0.0:                                                #        - Entries greater than unity are replaced with the value '-1.0'.
                ld_measurements_table['ARIMA_OCCUPANCY'][i] = numpy.float64(-1.0)                #        - Minimum: 0.0
            elif numpy.float64(line_bits[7]) > 1.0:                                              #        - Maximum: 1.0
                ld_measurements_table['ARIMA_OCCUPANCY'][i] = numpy.float64(-1.0)                #        - Q: How is vehicle length accounted/corrected for?
            else:                                                                                #          A: It's not done directly. It is calibrated via the free flow branch
                ld_measurements_table['ARIMA_OCCUPANCY'][i] = numpy.float64(line_bits[7])        #             of the FD and the MFD.
        if line_bits[8] == 'NA':                                                                 # Notes: - Average vehicle speed in the measurement interval (??UNITS??; ARIMA smoothed).
            ld_measurements_table['ARIMA_SPEED'][i] = numpy.float64(-1.0)                        #        - Some cities provide average speed measurements instead of occupancy.
        else:                                                                                    #          Otherwise, average speeds are calculated from speed = flow / density.
            if numpy.float64(line_bits[8]) < 0.0:                                                #        - Entries equal to 'NA' are replaced with the value '-1.0'.
                ld_measurements_table['ARIMA_SPEED'][i] = numpy.float64(-1.0)                    #        - Entries with negative values are replaced with the value '-1.0'.
            else:
                ld_measurements_table['ARIMA_SPEED'][i] = numpy.float64(line_bits[8])

    # Sort the loop detector measurements data table by "DATE" followed by "INTERVAL_START"
    ld_measurements_table = ld_measurements_table[numpy.argsort(ld_measurements_table, order = ['DATE', 'INTERVAL_START'])]

    # Set "ERROR_FLAG" to "1" for any duplicated entries for "DATE + INTERVAL_START" in this data table
    for i in range(nld_measurements - 1):
        curr_date_interval = ld_measurements_table['DATE'][i] + '_' + str(ld_measurements_table['INTERVAL_START'][i])
        next_date_interval = ld_measurements_table['DATE'][i + 1] + '_' + str(ld_measurements_table['INTERVAL_START'][i + 1])
        if next_date_interval == curr_date_interval:
            ld_measurements_table['ERROR_FLAG'][i] = numpy.int32(1)
            ld_measurements_table['ERROR_FLAG'][i + 1] = numpy.int32(1)

    # Write out the loop detector measurements data table for the current city and detector ID
    ld_measurements_table.write(curr_output_file, format = 'fits')
    os.remove(file)

# N.B: There are some duplicated entries for "DATE + INTERVAL_START" in the loop detector measurements tables.
# These entries have had "ERROR_FLAG" set to "1".

# Finish
print('')
print('Finished!')
