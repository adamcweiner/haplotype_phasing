import sys
import numpy as np
from read_data import read_data
from break_to_chunks import min_chunk_size
from clarks import Clark

#makes sure there is a command line argument
if len(sys.argv) != 2:
	print("Please input filename as command line argment")

input_file = sys.argv[1]

matrix = read_data(input_file) #reads data into numpy array (matrix)
end_pos = []
next_chunk_start_pos = []
for ii in range(1,10): # example of just 10 chunks
    temp1, temp2 = min_chunk_size(matrix, next_chunk_start_pos[ii-1])
    end_pos.append(temp1)
    next_chunk_start_pos.append(temp2)
print("end_pos: ", end_pos)
print("next_chunk_start_pos: ", next_chunk_start_pos)
'''
chunks = break_to_chunks(matrix, 10) #breaks into chunks including overlapping chunks

full_hap1 = np.full((matrix.shape), None)
full_hap2 = full_hap1.copy()

# start playing with Clark's algorithm by passing in the first chunk in the list of chunks
print("len(chunks): " + str(len(chunks)))
finished_chunks = 0
unfinished_chunks = 0
for bit in chunks:
    c = Clark(bit)
    temp1, temp2, num = c.run() # "temp" haplotypes are arrays of shape (10, 50)
    # TODO: find some way to cleverly/correctly append temp1 and temp2 into either full_hap1 and full_hap2
    if num == 50:
        finished_chunks += 1
    else:
        unfinished_chunks += 1

print("finished: " + str(finished_chunks))
print("unfinished: " + str(unfinished_chunks))
'''
