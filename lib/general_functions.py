__author__ = 'Dan Bramich'


# Count the number of lines in a text file
def count_file_lines(filename):
    count = 0
    with open(filename, mode = 'r') as file:
        for line in file: count += 1
    return count
