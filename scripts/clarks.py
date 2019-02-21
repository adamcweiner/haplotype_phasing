""" Implement Clark's algorithm on a chunk to eliminate some of the haplotypes. """
import numpy as np

class Clark:
    """ Holds all functions relevant to Clark's algorithm. """
    def __init__(self, chunk):
        """ A single genotype `chunk` is passed into the class from the `chunks` list. If the chunks are of size 10 then each `chunk` is an array of shape (10, 50). """
        self.chunk = chunk

    def complementary_hap(self, haplotype, genotype):
        """ Given a genotype and a haplotype, find the haplotype that complements the one given. """
        assert len(haplotype) == len(genotype) # make sure the two arrays are of the same size
        comp_hap = np.zeros((len(haplotype))) # fill comp_hap with just zeros to start
        for ii in range(len(comp_hap)):
            comp_hap[ii] = genotype[ii] - haplotype[ii] # i.e. if g[ii] = 2 and h[ii] = 1 then c_h[ii] = 1
        return comp_hap
