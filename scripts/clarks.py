""" Implement Clark's algorithm on a chunk to eliminate some of the haplotypes. """
import numpy as np

class Clark:
    """ Holds all functions relevant to Clark's algorithm. """
    def __init__(self, chunk):
        """ A single genotype `chunk` is passed into the class from the `chunks` list. If the chunks are of size 10 then each `chunk` is an array of shape (10, 50). """
        self.chunk = chunk.T # transpose this to make chunk of shape (50, 10)... this becomes easier to index
        self.N = self.chunk.shape[0] # represents the number of genotypes collected (should be 50)
        self.K = self.chunk.shape[1] # represents the length of each genotype collected (should be size_of_chunk)
        self.known_hap = set() # holds our set of known haplotypes
        
    def run_clarks(self):
        """ Calls proper functions to correctly run Clark's algorithm. """
        # add known haplotypes from homogeneous genotypes
        for n in range(self.N):
            if self.is_homogeneous(self.chunk[n]):
                temp_hap = self.find_homo_hap(self.chunk[n])
                self.known_hap.add(temp_hap) # append the haplotype to the known list


    def complementary_hap(self, haplotype, genotype):
        """ Given a genotype and a haplotype, find the haplotype that complements the one given. """
        assert len(haplotype) == len(genotype) # make sure the two arrays are of the same size
        comp_hap = np.zeros((len(haplotype))) # fill comp_hap with just zeros to start
        for ii in range(len(comp_hap)):
            comp_hap[ii] = genotype[ii] - haplotype[ii] # i.e. if g[ii] = 2 and h[ii] = 1 then c_h[ii] = 1
        return comp_hap

    def valid_hap(self, hap):
        """ Returns True if the haplotype is valid (meaning all entries are just 0's and 1's)"""
        temp = True
        for ii in range(len(hap)):
            if hap[ii] > 1 or hap < 0: # if element is not 0 or 1
                temp = False
                break
        return temp

    def is_homogeneous(self, gen):
        """ Returns True if the genotype only consists of 0's and 2's; returns False if any of the positions are heterogeneous. """
        temp = True
        for ii in range(len(gen)):
            if gen[ii] == 1 # if element is 1, break and return False
                temp = False
                break
        return temp

    def find_homo_hap(self, gen):
        """ Gives the haplotype of a homogenous genotype. All g[i]=2's turn into h[i]=1's and all g[i]=0's turn into h[i]=0's. """
        hap = np.zeros((len(gen))) # start with an array of 0's
        for ii in range(len(gen)):
            assert gen[ii] == 0 or gen[ii] == 2 # make sure it's a 0 or 2
            if gen[ii] == 2:
                hap[ii] = 1
        return hap
        
