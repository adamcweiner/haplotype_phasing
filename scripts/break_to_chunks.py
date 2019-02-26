""" Set of functions that convert the large df of genotypes into chunks of genotypes that can be phased individually and then annealed together. """
def break_to_chunks(df, size_of_chunk):
    ''' Function: break_to_chunks(df, size_of_chunk)
input:
	df: numpy matrix of genotypes
	size_of_chunk: number of SNPs to include in each chunk
return:
	chunk_list: list of chunks including the overlap part
example: assume for one person there is a string [1,2,3,4,5,6,7,8,9,10,11,12]
	The list of chunks lets say of size 4 would be [1,2,3,4] [3,4,5,6] [5,6,7,8] [7,8,9,10] [9,10,11,12]
'''
	number_of_SNPS = len(df)

	chunk_list = list()
	for i in range(0,number_of_SNPS, size_of_chunk//2):
		chunk_list.append(df[i:i+size_of_chunk])
	
	return chunk_list

def min_chunk_size(df, start_pos):
    """ For a given df and start position, return the position at which this current chunk should end and the position at which the next chunk should begin. The current chunk should end when all 50 individuals have encountered at least one `1` in the genotype. The next chunk should begin at the lowest position (numerically) where the last `1` was seen (across all individuals). """
    ii = 0
    seen_het = np.full((df.shape[0]), False) # True if said individual has encountered a 1, False if it hasn't encountered a 1
    while np.all(seen_het) is False: # until every position is True
