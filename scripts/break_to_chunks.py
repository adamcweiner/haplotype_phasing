""" Set of functions that convert the large df of genotypes into chunks of genotypes that can be phased individually and then annealed together. """
import numpy as np

'''def break_to_chunks(df, size_of_chunk):
    number_of_SNPS = len(df)
    chunk_list = list()
    for i in range(0,number_of_SNPS, size_of_chunk//2):
        chunk_list.append(df[i:i+size_of_chunk])

	return chunk_list
'''
def smart_chunking(df):
    """ Uses find_chunk_size to break the dataframe into chunks of various sizes. """
    n_snp, m_ind = df.shape

    end_pos = []
    start_pos = [0] # want to start at the 0th SNP for the first iteration
    chunk_list = [] # output list of the chunks
    ii = 0 # index for moving through next_chunk_start_pos
    while end_pos[ii] < n_snp or len(end_pos)==0: # length check is just to protect against first entry into loop
        temp_start = start_pos[ii]
        temp_end, temp_next_start = find_chunk_size(matrix, temp_start)
        end_pos.append(temp_end) # position for end of current chunk
        start_pos.append(temp_next_start) # position for start of next chunk
        chunk_list.append(df[temp_start:temp_end])
        ii += 1
     
    # TODO: handle edge cases for end of genome
     
    return chunk_list, next_chunk_start_pos
    

def find_chunk_size(df, start_pos, num_call=0):
    """ For a given df and start position, return the position at which this current chunk should end and the position at which the next chunk should begin. The current chunk should end when all 50 individuals have encountered at least one `1` in the genotype. The next chunk should begin at the lowest position (numerically) where the last `1` was seen (across all individuals). """
    row_count = 0
    n_snp, m_ind = df.shape
    seen_het = np.full((m_ind), False) # True if said individual has encountered a 1, False if it hasn't encountered a 1
    last_het = np.zeros((m_ind)) # stores the last heterozygous position seen for each individual in this chunk
    while np.all(seen_het) == False: # until every position is True
        for ii in range(m_ind):
            if df[start_pos + row_count, ii] == 1: # if a 1 is encountered
                seen_het[ii] = True # assign position to True
                last_het[ii] = start_pos + row_count # the 1 was seen in this row
        row_count += 1 # iterate the row count at end
    end_pos = start_pos + row_count
    #print(last_het)
    next_chunk_start_pos = np.min(last_het) + num_call # the next chunk needs to start at this position in order to guarantee at least one `1` in the overlap for all individuals. Shift up by num_call to account for recursive calls
    orig_start_pos = start_pos - num_call # want to compare the loop to the original start position instead of local one
    #print("num_call: ", num_call)
    #print("orig_start_pos: ", orig_start_pos)
    if next_chunk_start_pos == orig_start_pos:
        end_pos, next_chunk_start_pos = find_chunk_size(df, start_pos+1, num_call=num_call+1) # recursively call this function there is separation between the original start_pos and the start position of the next chunk
    return int(end_pos), int(next_chunk_start_pos)
