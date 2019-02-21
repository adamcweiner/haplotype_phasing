""" Implement Clark's algorithm on a chunk to eliminate some of the haplotypes. """
import numpy as np
from nodes import Gen, Hap

class Clark:
    """ Holds all functions relevant to Clark's algorithm. """
    def __init__(self, chunk):
        """ A single genotype `chunk` is passed into the class from the `chunks` list. If the chunks are of size 10 then each `chunk` is an array of shape (10, 50). """
        self.chunk = chunk.T # transpose this to make chunk of shape (50, 10)... this becomes easier to index
        self.N = self.chunk.shape[0] # represents the number of genotypes collected (should be 50)
        self.K = self.chunk.shape[1] # represents the length of each genotype collected (should be size_of_chunk)
        self.g_nodes = [] # list of all our genotype nodes
        for n in range(self.N):
            self.g_nodes.append(Gen(self.chunk[n])) # create a Gen object with the proper sequence and append it to the list
        self.known_hap = [] # holds our list of known haplotype nodes

    def run(self):
        """ Calls proper functions to correctly run Clark's algorithm. """
        # add known haplotypes from homogeneous genotypes
        for n in range(self.N):
            if self.g_nodes[n].is_homogeneous():
                temp_hap = self.find_homo_hap(self.g_nodes[n])
                self.known_hap.add(temp_hap) # append the haplotype to the known list
        print(self.known_hap)



    def find_homo_hap(self, gen):
        """ Gives the haplotype of a homogenous genotype. All g[i]=2's turn into h[i]=1's and all g[i]=0's turn into h[i]=0's. """
        hap = np.zeros((len(gen))) # start with an array of 0's
        for ii in range(len(gen)):
            assert gen[ii] == 0 or gen[ii] == 2 # make sure it's a 0 or 2
            if gen[ii] == 2:
                hap[ii] = 1
        return hap
        
