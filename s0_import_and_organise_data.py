#!/usr/bin/env python

__author__ = 'Dan Bramich'

# This script reads in the data on the loop detectors, the road links in which they are installed, and the
# traffic measurements that they have made, and it splits these data into FITS binary table files by city
# while also standardising the entries. The data come from the publication "Understanding traffic capacity
# of urban networks" by Loder et al. (2019).

# Imports
import csv
#import os

# Parameters
detectors_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Locations/detectors.csv'
links_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Locations/links.csv'
measurements_file = '/data/dmb20/Loop.Detector.Data.41.Cities/Loop.Detector.Measurements/loops_measurements.csv'
output_dir = '/data/dmb20/Loop.Detector.Data.Analysis.MFD.Project.1/Results'

# Prepare a table for the loop detector data
print('')
print('Preparing a table for the loop detector data...')






# Read in the data file on the loop detectors
print('')
print('Reading in the loop detectors data file: ' + detectors_file)
with open(detectors_file, 'r') as csv_file:
    csv_contents = csv.DictReader(csv_file)
    for row in csv_contents:

#### ABOVE FULLY READ AND TESTED

        print(row)

        exit()
