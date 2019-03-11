import sys
import os
import numpy as np
from read_data import read_data
from break_to_chunks import smart_chunking
from clarks import Clark
from em_algorithm import em_algorithm
from em_algorithm2 import em_algorithm2
from em_algorithm3 import em_algorithm3, one2two
from merge_chunks import merge_chunks
from time import time


start = time()

#make sure there is a command line argument
if len(sys.argv) != 2:
	print("Please input filename as command line argment")

input_file = sys.argv[1]

matrix = read_data(input_file) #reads data into numpy array (matrix)
chunk_list, start_pos, end_pos = smart_chunking(matrix, max_snps=8, end_shift=5)
print("length of chunk list:", len(chunk_list))

print("running EM without Clarks")
solved_chunks = []
for ii, chunk in enumerate(chunk_list):
    temp_chunk = one2two(chunk)
    temp_out = em_algorithm3(temp_chunk)
    print("progress:", ii / len(chunk_list))
    print("temp_chunk.shape:", temp_chunk.shape)
    solved_chunks.append(np.asarray(temp_out))


solution = merge_chunks(solved_chunks, start_pos, end_pos)
end = time()
#print(solution)
print(end - start," seconds")

filename = os.path.splitext(input_file)[0]+'.solutions.txt'
output_file = open(filename,'w')
for i in range(len(solution[0])):
	for j in range(len(solution)):
		output_file.write(str(solution[j][i]))
		if(j!=len(solution)-1):
			output_file.write(' ')
	output_file.write('\n')

'''
print("running Clarks to completion")


print("running partial clarks and then EM2")
c = Clark(chunk_list[0])
temp0, num = c.run() # "temp" haplotypes are arrays of shape (n, 2m)
print("temp0.shape:", temp0.shape)
print(temp0.dtype)
out = em_algorithm2(temp0)
print(out)

finished_chunks = 0
unfinished_chunks = 0
num_haplo_phased = []
for chunk in chunk_list:
    print("current chunk shape:", chunk.shape)
    # out = em_algorithm(chunk)
    c = Clark(chunk)
    temp1, num = c.run() # "temp" haplotypes are arrays of shape (n, 2m)
    # TODO: find some way to cleverly/correctly append temp1 and temp2 into either full_hap1 and full_hap2
    num_haplo_phased.append(num)
    if num == 50:
        finished_chunks += 1
    else:
        unfinished_chunks += 1

print("avg number of haplotypes phased per chunk:", sum(num_haplo_phased) / len(num_haplo_phased) )
print("finished:", finished_chunks)
print("unfinished:", unfinished_chunks)
'''
