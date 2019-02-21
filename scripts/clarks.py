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
            if self.g_nodes[n].n_ones == 0: # if the sequence is homogenous (all 0's and 2's)
                temp_hap_seq = self.find_homo_hap(self.g_nodes[n].seq)
                temp_ptr = self.unique_haplo(temp_hap_seq)
                if temp_ptr is None: # when temp_hap_seq is truly unique
                    temp_hap_node = Hap(temp_hap_seq) # create new node (h1)
                    g_node[n].add_hs(temp_hap_node) # add ptr to h1 and create+link h2.... TODO: ensure that this doesn't double count for homogenous genotypes
                    self.known_hap.append(temp_hap_node) # append new haplotype to the known_hap list
                # if temp_hap_seq is found as the sequence for a node in known_hap
                    # have the genotype point to this existing haplotype node
                # else
                    # create a new node according to the sequence
                    # have the g_node[n] point to this haplotype
                    # append new haplotype to the known_hap list
                temp_hap_node = Hap(temp_hap_seq) # create node
                self.known_hap.add(temp_hap) # append the haplotype node to the known list
            
        print(self.known_hap)


    def find_homo_hap(self, gen_seq):
        """ Gives the haplotype of a homogenous genotype. All g[i]=2's turn into h[i]=1's and all g[i]=0's turn into h[i]=0's. """
        hap_seq = np.zeros((len(gen_seq))) # start with an array of 0's
        for ii in range(len(gen_seq)):
            assert gen_seq[ii] == 0 or gen_seq[ii] == 2 # make sure it's a 0 or 2
            if gen_seq[ii] == 2:
                hap_seq[ii] = 1
        return hap_seq
        

    # need to create a function that checks if a haplotype sequence is unique compared to all that is found in self.known_hap
    def unique_haplo(self, hap_seq):
        """ Looks through all of self.known_hap and ensures that hap_seq is unique from all of the known haplotypes currently in the list. If the sequence already exists, it will return the address of the existing node. If the sequence is unique, it will return None. """
        temp = None
        for hap_node in self.known_hap:
            if np.array_equal(hap_node.seq, hap_seq): # returns true if two arrays are equal
                temp = hap_node # assign temp to current node and break out
                break
        return temp