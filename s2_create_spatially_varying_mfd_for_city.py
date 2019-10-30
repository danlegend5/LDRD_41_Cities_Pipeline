#!/usr/bin/env python

__author__ = 'Dan Bramich'

# For a specific city, this script constructs a spatially-varying MFD across a city area. The spatially-varying MFD is
# constructed by estimating the MFD for each region in a regular grid of overlapping circular regions. The grid spacing
# and the radius of the circular regions are user-defined.

# Imports
import argparse
import numpy
import os
from astropy.table import Table
from LDDA_MFD_Project1_Pipeline.config import config
from LDDA_MFD_Project1_Pipeline.lib import general_functions

# Set up the command line arguments
parser = argparse.ArgumentParser(description = 'This script is used to construct a spatially-varying MFD across the area '
                                               'of a specific city. The spatially-varying MFD is constructed by estimating '
                                               'the MFD for each region in a regular grid of overlapping circular regions. '
                                               'The grid spacing and the radius of the circular regions are user-defined.')
parser.add_argument('city_name', help = 'City name.')
parser.add_argument('grid_spacing', help = 'Distance (km) between the horizontal or vertical grid points.')
parser.add_argument('region_radius', help = 'Radius (km) of the circular region to be used at each grid point.')
args = parser.parse_args()

# Determine the locations of the required input data
try:
    city_name = str(args.city_name)
except:
    print('ERROR - The input parameter "city_name" could not be converted to a string...')
    exit()
country_name = general_functions.get_country_name(city_name)
if country_name == 'unknown':
    print('ERROR - The input parameter "city_name" is not recognised...')
    exit()
input_dir_ld_locations = os.path.join(config.output_dir, 's1.Loop.Detector.Locations')
input_dir_ld_measurements_raw = os.path.join(config.output_dir, 's1.Loop.Detector.Measurements.Raw')
input_dir_ld_measurements_arima = os.path.join(config.output_dir, 's1.Loop.Detector.Measurements.ARIMA')

# Check that the remaining input parameters make sense
try:
    grid_spacing = numpy.float64(args.grid_spacing)
except:
    print('ERROR - The input parameter "grid_spacing" could not be converted to a floating point number...')
    exit()
if grid_spacing <= 0.0:
    print('ERROR - The input parameter "grid_spacing" is zero or negative...')
    exit()
try:
    region_radius = numpy.float64(args.region_radius)
except:
    print('ERROR - The input parameter "region_radius" could not be converted to a floating point number...')
    exit()
if region_radius <= 0.0:
    print('ERROR - The input parameter "region_radius" is zero or negative...')
    exit()
region_radius2 = region_radius*region_radius

# Read in the loop detector locations data file
input_file_ld_locations = os.path.join(input_dir_ld_locations, 'detectors.' + country_name + '.' + city_name + '.fits')
print('')
print('Reading in the loop detector locations data file: ' + input_file_ld_locations)
ld_locations_table = Table.read(input_file_ld_locations, format = 'fits')
nld_locations = len(ld_locations_table)
print('No. of loop detectors: ' + str(nld_locations))

# Determine the limits in longitude and latitude of the roads hosting the loop detectors
selection_links = ld_locations_table['LINK_PTS_FLAG'] == 1
if numpy.count_nonzero(selection_links) == 0:
    print('ERROR - Unable to determine the limits in longitude and latitude of the roads hosting the loop detectors because information on these roads is not present in the data file...')
    exit()
min_link_pts_longitude = numpy.min(ld_locations_table['LINK_PTS_LONGITUDE'][selection_links])
max_link_pts_longitude = numpy.max(ld_locations_table['LINK_PTS_LONGITUDE'][selection_links])
min_link_pts_latitude = numpy.min(ld_locations_table['LINK_PTS_LATITUDE'][selection_links])
max_link_pts_latitude = numpy.max(ld_locations_table['LINK_PTS_LATITUDE'][selection_links])
print('Min. max. longitude of the roads hosting the loop detectors: ', min_link_pts_longitude, max_link_pts_longitude)
print('Min. max. latitude of the roads hosting the loop detectors: ', min_link_pts_latitude, max_link_pts_latitude)
for i in range(nld_locations):
    selection = ld_locations_table['LINK_PTS_FLAG'][i] == 1
    if numpy.count_nonzero(selection) <= 1:
        print('WARNING - The road hosting the following loop detector is undefined: ' + ld_locations_table['DETECTOR_ID'][i])

# Convert the loop detector locations, and the sets of coordinates that define the roads that host the loop detectors,
# from longitude and latitude to an x and y Cartesian coordinate system using the mean radius of the Earth. The
# conversion is approximate due to projection effects and the usage of the mean radius of the Earth. However, the
# differential errors in the scale of the Cartesian grid over the geographical area covered by the loop detectors and
# their host roads in any of the cities considered are less than 0.2%, and the error in the adopted radius of the Earth
# for any city is less than 0.1%.
print('')
print('Creating a grid covering the locations of the loop detectors and their host roads...')
rad_per_deg = numpy.pi/180.0
mean_ld_longitude = numpy.mean(ld_locations_table['LONGITUDE'])
mean_ld_latitude = numpy.mean(ld_locations_table['LATITUDE'])
cos_factor = numpy.cos(mean_ld_latitude*rad_per_deg)
ld_x = (config.radius_earth*rad_per_deg*cos_factor)*(ld_locations_table['LONGITUDE'] - mean_ld_longitude)
ld_y = (config.radius_earth*rad_per_deg)*(ld_locations_table['LATITUDE'] - mean_ld_latitude)
link_pts_x = (config.radius_earth*rad_per_deg*cos_factor)*(ld_locations_table['LINK_PTS_LONGITUDE'] - mean_ld_longitude)
link_pts_y = (config.radius_earth*rad_per_deg)*(ld_locations_table['LINK_PTS_LATITUDE'] - mean_ld_latitude)

# Create a grid covering the locations of the loop detectors and their host roads
min_x = numpy.minimum(numpy.min(ld_x), numpy.min(link_pts_x[selection_links]))
max_x = numpy.maximum(numpy.max(ld_x), numpy.max(link_pts_x[selection_links]))
min_y = numpy.minimum(numpy.min(ld_y), numpy.min(link_pts_y[selection_links]))
max_y = numpy.maximum(numpy.max(ld_y), numpy.max(link_pts_y[selection_links]))
range_x = max_x - min_x
range_y = max_y - min_y
mid_x = 0.5*(min_x + max_x)
mid_y = 0.5*(min_y + max_y)
nx = numpy.int32(numpy.ceil(range_x/grid_spacing))
ny = numpy.int32(numpy.ceil(range_y/grid_spacing))
lo_x = mid_x - (0.5*nx*grid_spacing)
lo_y = mid_y - (0.5*ny*grid_spacing)
nx = nx + 1
ny = ny + 1
grid_x = lo_x + (grid_spacing*numpy.arange(0, nx, 1, dtype = numpy.float64))
grid_y = lo_y + (grid_spacing*numpy.arange(0, ny, 1, dtype = numpy.float64))

# For each grid point
for i in range(nx):
    curr_grid_x = grid_x[i]
    ld_dist_x = ld_x - curr_grid_x
    ld_dist_x2 = ld_dist_x*ld_dist_x
    for j in range(ny):
        curr_grid_y = grid_y[j]
        ld_dist_y = ld_y - curr_grid_y
        ld_dist_y2 = ld_dist_y*ld_dist_y

        # Determine the set of loop detectors that lie within a circular region around the current grid point
        print('')
        print('Analysing the grid point (' + str(i) + ', ' + str(j) + ') from a grid of size (' + str(nx) + ', ' + str(ny) + ')...')
        ld_selection = (ld_dist_x2 + ld_dist_y2) <= region_radius2
        nld_selection = numpy.count_nonzero(ld_selection)
        if nld_selection == 0:
            print('There are no loop detectors within the associated circular region...')
            continue
        print('Number of loop detectors within the associated circular region: ' + str(nld_selection))


#### ABOVE FULLY READ AND TESTED


        # ANALYSE NETWORK FOR CONNECTED GROUPS OF ROADS - ONLY CONSIDER THE LARGEST CONNECTED GROUP (TOTAL ROAD LENGTH?) - REJECT LOOP DETECTORS THAT ARE NOT ON ROADS IN THIS LARGEST GROUP
        # RECORD NUMBER OF DISJOINT GROUPS OF ROADS

        # AT THIS POINT CALL FUNCTION TO CONSTRUCT MFD FOR SET OF LOOP DETECTORS...

        # FIT THE MFD AND EXTRACT CRITICAL DENSITY ETC.

        # CALCULATE UNDERLYING NETWORK PROPERTIES FOR THE REGION

        # STORE RESULTS





exit()




import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure(dpi = 150)
xlo = numpy.min(link_pts_x - 0.1)
xhi = numpy.max(link_pts_x + 0.1)
plt.xlim(xlo, xhi)
ylo = numpy.min(link_pts_y - 0.1)
yhi = numpy.max(link_pts_y + 0.1)
plt.ylim(ylo, yhi)
plt.scatter(link_pts_x, link_pts_y, s = 1, marker = '.', color = 'red')
plt.scatter(ld_x, ld_y, s = 1, marker = '.', color = 'black')
fig.savefig('plot.png')
plt.close()

