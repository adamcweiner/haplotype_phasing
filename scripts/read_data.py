import numpy as np
#Function read_data
#input:
#	input_file: Name or path of the input file to read as numpy array
#return:
#	matrix: Numpy array to return
def read_data(input_file):
	matrix = np.genfromtxt(input_file, dtype = int, delimiter = " ")
	return matrix
