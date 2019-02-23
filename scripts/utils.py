import numpy as np

def complementary_hap(hap_seq, gen_seq):
    """ Given genotype and a haplotype sequences, find the haplotype sequence that complements the one given. """
    assert len(hap_seq) == len(gen_seq) # make sure the two arrays are of the same size
    comp_hap = np.zeros((len(hap_seq))) # fill comp_hap with just zeros to start
    for ii in range(len(hap_seq)):
        comp_hap[ii] = gen_seq[ii] - hap_seq[ii] # i.e. if g[ii] = 2 and h[ii] = 1 then c_h[ii] = 1
    return comp_hap

def num_het(gen_seq):
    """ Returns the number of heterogeneous positions (1's) in the genotype sequence. """
    count = 0
    for ii in range(len(gen_seq)):
        if gen_seq[ii] == 1: # add 1 to count if hetero site encountered
            count += 1
    return count

def find_homo_hap(gen_seq):
    """ Gives the haplotype of a homogenous genotype. All g[i]=2's turn into h[i]=1's and all g[i]=0's turn into h[i]=0's. """
    hap_seq = np.zeros((len(gen_seq))) # start with an array of 0's
    for ii in range(len(gen_seq)):
        assert gen_seq[ii] == 0 or gen_seq[ii] == 2 # make sure it's a 0 or 2
        if gen_seq[ii] == 2:
            hap_seq[ii] = 1
    return hap_seq

def valid_hap(hap_seq):
    """ Returns True if the haplotype is valid (meaning all entries are just 0's and 1's)"""
    temp = True
    for ii in range(len(hap_seq)):
        if hap_seq[ii] > 1 or hap_seq[ii] < 0: # if element is not 0 or 1
            temp = False
            break
    return temp