""" Set of functions that convert the large df of genotypes into chunks of genotypes that can be phased individually and then annealed together. """
import numpy as np

'''def break_to_chunks(df, size_of_chunk):
    number_of_SNPS = len(df)
    chunk_list = list()
    for i in range(0,number_of_SNPS, size_of_chunk//2):
        chunk_list.append(df[i:i+size_of_chunk])

	return chunk_list
'''
def min_chunk_size(df, start_pos):
    """ For a given df and start position, return the position at which this current chunk should end and the position at which the next chunk should begin. The current chunk should end when all 50 individuals have encountered at least one `1` in the genotype. The next chunk should begin at the lowest position (numerically) where the last `1` was seen (across all individuals). """
    row_count = 0
    n_snp, m_ind = df.shape
    seen_het = np.full((m_ind), False) # True if said individual has encountered a 1, False if it hasn't encountered a 1
    while np.all(seen_het) == False: # until every position is True
        for ii in range(m_ind):
            if df[start_pos + row_count, ii] == 1: # if a 1 is encountered
                seen_het[ii] = True # assign position to True
        row_count += 1 # iterate the row count at end
    end_pos = start_pos + row_count
    return end_pos
