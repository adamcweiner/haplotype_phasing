import sys
import numpy as np
from read_data import read_data
from break_to_chunks import smart_chunking
from clarks import Clark
from em_algorithm import em_algorithm

#makes sure there is a command line argument
if len(sys.argv) != 2:
	print("Please input filename as command line argment")

input_file = sys.argv[1]

matrix = read_data(input_file) #reads data into numpy array (matrix)
chunk_list, start_pos, end_pos = smart_chunking(matrix)
print("length of chunk list:", len(chunk_list))



finished_chunks = 0
unfinished_chunks = 0
num_haplo_phased = []
for chunk in chunk_list:
    print("current chunk shape:", chunk.shape)
    # out = em_algorithm(chunk)
    c = Clark(chunk)
    temp1, temp2, num = c.run() # "temp" haplotypes are arrays of shape (10, 50)
    # TODO: find some way to cleverly/correctly append temp1 and temp2 into either full_hap1 and full_hap2
    num_haplo_phased.append(num)
    if num == 50:
        finished_chunks += 1
    else:
        unfinished_chunks += 1

print("avg number of haplotypes phased per chunk:", sum(num_haplo_phased) / len(num_haplo_phased) )
print("finished:", finished_chunks)
print("unfinished:", unfinished_chunks)
