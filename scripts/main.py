import sys
import numpy as np
from read_data import read_data
from break_to_chunks import break_to_chunks
from clarks import Clark

#makes sure there is a command line argument
if len(sys.argv) != 2:
	print("Please input filename as command line argment")

input_file = sys.argv[1]

matrix = read_data(input_file) #reads data into numpy array (matrix)
chunks = break_to_chunks(matrix, 10) #breaks into chunks including overlapping chunks

# start playing with Clark's algorithm by passing in the first chunk in the list of chunks
c = Clark(chunks[0])
#c.run()

