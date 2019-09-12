#!/usr/bin/env python

__author__ = 'Dan Bramich'

# This script reads in the data on the loop detectors, the road links in which they are installed, and the
# traffic measurements that they have made, and it splits these data into FITS binary table files by city
# while also standardising the entries. The data come from the publication "Understanding traffic capacity
# of urban networks" by Loder et al. (2019).

# Imports
import csv
import numpy
from astropy.table import Table
from LDDA_MFD_Project1_Pipeline.config import config
from LDDA_MFD_Project1_Pipeline.lib import general_functions

# Prepare a table for the loop detector locations data
print('')
print('Preparing a table for the loop detector locations data...')
nld_locations = general_functions.count_file_lines(config.original_detectors_file) - 1
empty_str = '                                                                                '
ld_locations_table = Table([[empty_str]*nld_locations,
                            numpy.zeros(nld_locations, dtype = numpy.float64),
                            numpy.zeros(nld_locations, dtype = numpy.float64),
                            numpy.zeros(nld_locations, dtype = numpy.float64),
                            numpy.zeros(nld_locations, dtype = numpy.float64),
                            [empty_str]*nld_locations,
                            [empty_str]*nld_locations,
                            numpy.zeros(nld_locations, dtype = numpy.float64),
                            numpy.zeros(nld_locations, dtype = numpy.int32),
                            numpy.zeros(nld_locations, dtype = numpy.int32),
                            [empty_str]*nld_locations],
                           names = ('DETECTOR_ID', 'LONGITUDE', 'LATITUDE', 'LENGTH', 'POSITION',
                                    'ROAD_NAME', 'ROAD_CLASS', 'SPEED_LIMIT', 'NLANES', 'LINK_ID', 'CITY_NAME'))

# Read in the loop detector locations data file
print('')
print('Reading in the loop detector locations data file: ' + config.original_detectors_file)
with open(config.original_detectors_file, mode = 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    i = 0
    for row in csv_contents:
        ld_locations_table['DETECTOR_ID'][i] = row['detid'].replace(' ', '_')        # Notes: - Spaces replaced with underscores
        ld_locations_table['LONGITUDE'][i] = numpy.float64(row['coords.x1'])         # Notes: - Longitudes are quoted in range -180.0 to 180.0 degrees
                                                                                     #        - Minimum: -118.303722 deg
                                                                                     #        - Maximum: 145.08013040207 deg
        ld_locations_table['LATITUDE'][i] = numpy.float64(row['coords.x2'])          # Notes: - Latitudes are quoted in range -90.0 to 90.0 degrees
                                                                                     #        - Minimum: -37.8745292808731 deg
                                                                                     #        - Maximum: 54.70488 deg
        ld_locations_table['LENGTH'][i] = numpy.float64(row['length'])               # Notes: - Minimum: 0.007199797164467 ???UNITS???
                                                                                     #        - Maximum: 16.6844474508417 ???UNITS???
        if row['pos'] == 'NA':                                                       # Notes: - 45 out of 1042 entries for "utrecht" are equal to 'NA'. These are replaced with
            ld_locations_table['POSITION'][i] = numpy.float64(-1.0)                  #          the value '-1.0'.
        else:                                                                        #        - Minimum: 0.0 ???UNITS???
            ld_locations_table['POSITION'][i] = numpy.float64(row['pos'])            #        - Maximum: 6.67473242607813 ???UNITS???
        ld_locations_table['ROAD_NAME'][i] = row['road'].replace(' ', '_')           # Notes: - Spaces replaced with underscores
                                                                                     #        - 1544 entries are equal to 'NA'
                                                                                     #        - Many strings are unintelligible
        ld_locations_table['ROAD_CLASS'][i] = row['fclass']                          # Notes: - Set of possible entries:
                                                                                     #          'cycleway', 'footway', 'living_street', 'motorway', 'motorway_link', 'other',
                                                                                     #          'path', 'pedestrian', 'primary', 'primary_link', 'residential', 'secondary',
                                                                                     #          'secondary_link', 'service', 'tertiary', 'tertiary_link', 'trunk', 'trunk_link',
                                                                                     #          'unclassified'
        if row['limit'] == 'NA':                                                     # Notes: - 6632 entries are equal to 'NA' ("bremen", "london", "losanageles", "madrid",
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #          "melbourne", "paris", "tokyo", "utrecht"). These are replaced with the value
        elif row['limit'] == '50|50|30':                                             #          '-1.0'.
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #        - 3 entries corresponding to "Calle de Segovia" in "madrid" are equal to '50|50|30'.
        elif row['limit'] == '60; 40; 60':                                           #          These indicate different speed limits in each lane. For simplicity, these are
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #          replaced with the value '-1.0'.
        elif row['limit'] == '40; 60':                                               #        - 3 entries corresponding to "La Trobe Street" in "melbourne" are equal to
            ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)               #          '60; 40; 60'. These indicate different speed limits in each lane. For simplicity,
        else:                                                                        #          these are replaced with the value '-1.0'.
            if numpy.float64(row['limit']) == 0.0:                                   #        - 2 entries corresponding to "La Trobe Street" in "melbourne" are equal to '40; 60'.
                ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(-1.0)           #          These indicate different speed limits in each lane. For simplicity, these are
            else:                                                                    #          replaced with the value '-1.0'.
                ld_locations_table['SPEED_LIMIT'][i] = numpy.float64(row['limit'])   #        - 3754 entries are equal to '0.0' ("BORDEAUX", "CONSTANCE", "GRAZ", "SANTANDER",
                                                                                     #          "augsburg", "basel", "bern", "birmingham", "bolton", "bremen", "cagliari",
                                                                                     #          "darmstadt", "duisburg", "essen", "frankfurt", "hamburg", "innsbruck", "kassel",
                                                                                     #          "luzern", "manchester", "marseille", "rotterdam", "speyer", "strasbourg",
                                                                                     #          "stuttgart", "taipeh", "tokyo", "torino", "toronto", "toulouse", "utrecht",
                                                                                     #          "vilnius", "wolfsburg", "zurich"). These are replaced with the value '-1.0'.
                                                                                     #        - Minimum: 17.6 ???UNITS???
                                                                                     #        - Maximum: 100.0 ???UNITS???
        if row['lanes'] == 'NA':                                                     # Notes: - 4 out of 438 entries for "BORDEAUX" are equal to 'NA'. These are replaced with
            ld_locations_table['NLANES'][i] = numpy.int32(-1)                        #          the value '-1'.
        else:                                                                        #        - Minimum: 1
            ld_locations_table['NLANES'][i] = numpy.int32(row['lanes'])              #        - Maximum: 9
        if row['linkid'] == 'NA':                                                    # Notes: - 605 entries are equal to 'NA' ("BORDEAUX", "CONSTANCE", "utrecht"). These are
            ld_locations_table['LINK_ID'][i] = numpy.int32(-1)                       #          replaced with the value '-1'.
        else:                                                                        #        - Minimum: 0
            ld_locations_table['LINK_ID'][i] = numpy.int32(row['linkid'])            #        - Maximum: 5268
        if row['citycode'] == 'BORDEAUX':                                            # Notes: - Four cities have strings in upper case ("BORDEAUX", "CONSTANCE", "GRAZ",
            ld_locations_table['CITY_NAME'][i] = 'bordeaux'                          #          "SANTANDER"). These are converted to lower case.
        elif row['citycode'] == 'CONSTANCE':                                         #        - Los angeles is spelt wrong "losanageles". This is corrected.
            ld_locations_table['CITY_NAME'][i] = 'constance'                         #        - 41 cities with names:
        elif row['citycode'] == 'GRAZ':                                              #          'augsburg', 'basel', 'bern', 'birmingham', 'bolton', 'bordeaux', 'bremen',
            ld_locations_table['CITY_NAME'][i] = 'graz'                              #          'cagliari', 'constance', 'darmstadt', 'duisburg', 'essen', 'frankfurt', 'graz',
        elif row['citycode'] == 'SANTANDER':                                         #          'groningen', 'hamburg', 'innsbruck', 'kassel', 'london', 'losangeles', 'luzern',
            ld_locations_table['CITY_NAME'][i] = 'santander'                         #          'madrid', 'manchester', 'marseille', 'melbourne', 'munich', 'paris', 'rotterdam',
        elif row['citycode'] == 'losanageles':                                       #          'santander', 'speyer', 'strasbourg', 'stuttgart', 'taipeh', 'tokyo', 'torino',
            ld_locations_table['CITY_NAME'][i] = 'losangeles'                        #          'toronto', 'toulouse', 'utrecht', 'vilnius', 'wolfsburg', 'zurich'
        else:
            ld_locations_table['CITY_NAME'][i] = row['citycode']
        i += 1


#### ABOVE FULLY READ AND TESTED
