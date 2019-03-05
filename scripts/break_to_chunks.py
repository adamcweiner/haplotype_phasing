""" Set of functions that convert the large df of genotypes into chunks of genotypes that can be phased individually and then annealed together. """
import numpy as np

def break_to_chunks(df, size_of_chunk):
    number_of_SNPS = len(df)
    chunk_list = list()
    for i in range(0,number_of_SNPS, size_of_chunk//2):
        chunk_list.append(df[i:i+size_of_chunk])

    return chunk_list

def smart_chunking(df, max_snps=10):
    """ Uses find_chunk_size to break the dataframe into chunks of various sizes. """
    n_snp, m_ind = df.shape

    end_pos = []
    start_pos = [0] # want to start at the 0th SNP for the first iteration
    chunk_list = [] # output list of the chunks

    def chunk_help(it):
        """ Finds chunk positions and appends appropriate values/chunks to lists """
        temp_start = start_pos[it]
        print("*** start find_chunk_size ***")
        temp_end, temp_next_start = find_chunk_size(df, temp_start, max_snps=max_snps)
        print("*** end find_chunk_size ***")
        end_pos.append(temp_end) # position for end of current chunk
        start_pos.append(temp_next_start) # position for start of next chunk
        chunk_list.append(df[temp_start:temp_end])

    chunk_help(0) # do first pass outside of loop since end_pos is empty
    ii = 1 # index for moving through next_chunk_start_pos

    # iterate until the end_pos reaches the end of the genome
    while end_pos[ii-1] < n_snp:
        chunk_help(ii)
        ii += 1
        print(len(chunk_list))

    start_pos.pop() # get rid of the last item on the list since it's never truly found

    return chunk_list, start_pos, end_pos


def find_chunk_size(df, start_pos, num_call=0, max_snps=10):
    """ For a given df and start position, return the position at which this current chunk should end and the position at which the next chunk should begin. The current chunk should end when all 50 individuals have encountered at least one `1` in the genotype. The next chunk should begin at the lowest position (numerically) where the last `1` was seen (across all individuals). """
    print("num_call:", num_call)
    row_count = 0
    n_snp, m_ind = df.shape
    seen_het = np.full((m_ind), False) # True if said individual has encountered a 1, False if it hasn't encountered a 1
    last_het = np.zeros((m_ind)) # stores the last heterozygous position seen for each individual in this chunk
    num_ones = np.zeros((m_ind)) # stores the number of 1's seen in the chunk
    while np.all(seen_het) == False: # until every position is True
        if start_pos + row_count == n_snp: # break out of while-loop if we reached the end of the genome
            break
        elif np.max(num_ones) >= max_snps: # break out of while-loop if any of the individuals reach the upper bound on SNPs
            print("reached max SNP count for this chunk")
            break
        else:
            for ii in range(m_ind):
                if df[start_pos + row_count, ii] == 1: # if a 1 is encountered
                    seen_het[ii] = True # assign position to True
                    last_het[ii] = start_pos + row_count # the 1 was seen in this row
                    num_ones[ii] += 1 # add 1 to the total number of ones for this individual
            row_count += 1 # iterate the row count at end
    end_pos = start_pos + row_count
    next_chunk_start_pos = np.min(last_het) + num_call # the next chunk needs to start at this position in order to guarantee at least one `1` in the overlap for all individuals. Shift up by num_call to account for recursive calls
    orig_start_pos = start_pos - num_call # want to compare the loop to the original start position instead of local one
    # recursively call this function until there is separation between the original start_pos and the start position of the next chunk
    if next_chunk_start_pos == orig_start_pos:
        end_pos, next_chunk_start_pos = find_chunk_size(df, start_pos+1, num_call=num_call+1)
    print("end_pos:", end_pos)
    print("next_chunk_start_pos:", next_chunk_start_pos)
    return int(end_pos), int(next_chunk_start_pos)
