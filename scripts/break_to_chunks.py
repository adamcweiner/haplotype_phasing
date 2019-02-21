# Function: break_to_chunks(df, size_of_chunk)
#input:
#	df: numpy matrix of genotypes
#	size_of_chunk: number of SNPs to include in each chunk
# return:
#	chunk_list: list of chunks including the overlap part
#example: assume for one person there is a string [1,2,3,4,5,6,7,8,9,10,11,12]
#	The list of chunks lets say of size 4 would be [1,2,3,4] [3,4,5,6] [5,6,7,8] [7,8,9,10] [9,10,11,12]
def break_to_chunks(df, size_of_chunk):
	number_of_SNPS = len(df)

	chunk_list = list()
	for i in range(0,number_of_SNPS, size_of_chunk//2):
		chunk_list.append(df[i:i+size_of_chunk])
	
	return chunk_list
		
