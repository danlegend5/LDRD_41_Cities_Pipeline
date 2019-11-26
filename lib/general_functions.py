__author__ = 'Dan Bramich'

# Imports
import numpy
from astropy.time import Time
from astropy.time import TimeDelta


# Convert a date string of the format YYYY-MM-DDTSSSSS, where SSSSS is the seconds after midnight, to a time
# object
def convert_time_stamp_to_time_object(time_stamp):
    time_isot = Time(time_stamp[0:11] + '00:00:00.000', format = 'isot', scale = 'tai')
    time_secs = TimeDelta(numpy.int32(time_stamp[11:16]), format = 'sec')
    return time_isot + time_secs


# Count the number of lines in a text file
def count_file_lines(filename):
    count = 0
    with open(filename, mode = 'r') as file:
        for line in file: count += 1
    return count


# Given an input city name, this function returns the name of the country to which the city belongs
def get_country_name(city_name):
    if city_name == 'augsburg': return 'germany'
    if city_name == 'basel': return 'switzerland'
    if city_name == 'bern': return 'switzerland'
    if city_name == 'birmingham': return 'uk'
    if city_name == 'bolton': return 'uk'
    if city_name == 'bordeaux': return 'france'
    if city_name == 'bremen': return 'germany'
    if city_name == 'cagliari': return 'italy'
    if city_name == 'constance': return 'germany'
    if city_name == 'darmstadt': return 'germany'
    if city_name == 'duisburg': return 'germany'
    if city_name == 'essen': return 'germany'
    if city_name == 'frankfurt': return 'germany'
    if city_name == 'graz': return 'austria'
    if city_name == 'groningen': return 'netherlands'
    if city_name == 'hamburg': return 'germany'
    if city_name == 'innsbruck': return 'austria'
    if city_name == 'kassel': return 'germany'
    if city_name == 'london': return 'uk'
    if city_name == 'losangeles': return 'usa'
    if city_name == 'luzern': return 'switzerland'
    if city_name == 'madrid': return 'spain'
    if city_name == 'manchester': return 'uk'
    if city_name == 'marseille': return 'france'
    if city_name == 'melbourne': return 'australia'
    if city_name == 'munich': return 'germany'
    if city_name == 'paris': return 'france'
    if city_name == 'rotterdam': return 'netherlands'
    if city_name == 'santander': return 'spain'
    if city_name == 'speyer': return 'germany'
    if city_name == 'strasbourg': return 'france'
    if city_name == 'stuttgart': return 'germany'
    if city_name == 'taipeh': return 'taiwan'
    if city_name == 'tokyo': return 'japan'
    if city_name == 'torino': return 'italy'
    if city_name == 'toronto': return 'canada'
    if city_name == 'toulouse': return 'france'
    if city_name == 'utrecht': return 'netherlands'
    if city_name == 'vilnius': return 'lithuania'
    if city_name == 'wolfsburg': return 'germany'
    if city_name == 'zurich': return 'switzerland'
    return 'unknown'
