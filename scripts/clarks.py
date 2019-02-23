""" Implement Clark's algorithm on a chunk to eliminate some of the haplotypes. """
import numpy as np
from utils import complementary_hap, num_het, find_homo_hap, valid_hap

class Clark:
    """ Holds all functions relevant to Clark's algorithm. """
    def __init__(self, chunk):
        """ A single genotype `chunk` is passed into the class from the `chunks` list. If the chunks are of size 10 then each `chunk` is an array of shape (10, 50). """
        chunk = chunk.T # transpose this to make chunk of shape (50, 10)... this becomes easier to index
        self.N = chunk.shape[0] # represents the number of genotypes collected (should be 50)
        self.K = chunk.shape[1] # represents the length of each genotype collected (should be size_of_chunk)
        self.chunk = np.zeros((3, self.N, self.K)) # 1st dim represents whether sequence is genotype (0), hap1 (1) or hap2 (2)
        for n in range(self.N):
            self.chunk[0, n] = chunk[n] # assign genotype sequence
        print(self.chunk)
        self.known_hap = [] # holds our list of known haplotype arrays

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

        

    # need to create a function that checks if a haplotype sequence is unique compared to all that is found in self.known_hap
    def unique_haplo(self, hap_seq):
        """ Looks through all of self.known_hap and ensures that hap_seq is unique from all of the known haplotypes currently in the list. If the sequence already exists, it will return the address of the existing node. If the sequence is unique, it will return None. """
        temp = None
        for hap_array in self.known_hap:
            if np.array_equal(hap_array, hap_seq): # returns true if two arrays are equal
                temp = hap_node # assign temp to current node and break out
                break
        return temp

    def add_hs(self, hap_seq, gen_seq):
        """ Adds the given haplotype and assigns the other haplotype haplotype assigned to the genotype. Note that h1 is a Hap object. """
        # TODO: convert to array format
        self.h1 = h1
        h2_seq = self.complementary_hap(h1.seq) # find sequence of second haplotype by finding complement
        self.h2 = Hap(h2_seq) # create Hap object given the proper sequence
        return self.h2 # return the complementary node in case it needs to be used
