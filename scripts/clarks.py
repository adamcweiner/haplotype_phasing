""" Implement Clark's algorithm on a chunk to eliminate some of the haplotypes. """
import numpy as np
from utils import complementary_hap, num_het, find_homo_hap, valid_hap, find_one_het

class Clark:
    """ Holds all functions relevant to Clark's algorithm. """
    def __init__(self, chunk):
        """ A single genotype `chunk` is passed into the class from the `chunks` list. If the chunks are of size 10 then each `chunk` is an array of shape (10, 50). """
        chunk = chunk.T # transpose this to make chunk of shape (50, 10)... this becomes easier to index
        self.N = chunk.shape[0] # represents the number of genotypes collected (should be 50)
        self.K = chunk.shape[1] # represents the length of each genotype collected (should be size_of_chunk)
        self.chunk = np.full((3, self.N, self.K), -1, dtype=int) # 1st dim represents whether sequence is genotype (0), hap1 (1) or hap2 (2)
        for n in range(self.N):
            self.chunk[0, n] = chunk[n] # assign genotype sequence
        self.known_hap = [] # holds our list of known haplotype arrays

    def run(self):
        """ Calls proper functions to correctly run Clark's algorithm. """
        num_phased = 0 # counter for the number of genotypes decoded (50 will mean the algorithm ran to completion)
        # add known haplotypes from homogeneous genotypes
        for n in range(self.N):
            temp_num_het = num_het(self.chunk[0, n]) # find number of 1's in the sequence
            if temp_num_het == 0: # if the sequence is homogenous (all 0's and 2's)
                temp_hap_seq = find_homo_hap(self.chunk[0, n]) # find the haplotype
                self.chunk[1, n] = temp_hap_seq
                self.chunk[2, n] = temp_hap_seq
                if self.unique_haplo(temp_hap_seq): # if this sequence isn't in the known list
                    self.known_hap.append(temp_hap_seq) # append the sequence to the list
                num_phased += 1
            if temp_num_het == 1: # we can assign the two haplotypes if there's only one "1" in the genotype
                temp_hap1, temp_hap2 = find_one_het(self.chunk[0, n]) # find the two haplotypes
                self.chunk[1, n] = temp_hap1
                self.chunk[2, n] = temp_hap2
                if self.unique_haplo(temp_hap1): # if this sequence isn't in the known list
                    self.known_hap.append(temp_hap1) # append the sequence to the list
                if self.unique_haplo(temp_hap2): # if this sequence isn't in the known list
                    self.known_hap.append(temp_hap2) # append the sequence to the list
                num_phased += 1


        # go through list repeatedly and see if any of the known haplotypes can be used for unphased genotypes
        enter = True
        while enter is True:
            enter = False # break out of while loop unless we phase a genotype
            for n in range(self.N):
                if self.chunk[1, n, 0] == -1: # only search if the haplotypes are still unknown
                    for hap_array in self.known_hap:
                        comp = complementary_hap(hap_array, self.chunk[0, n]) # find the complementary haplotype
                        # if the complementary haplotype is valid then use these two haplotypes         
                        if valid_hap(comp):
                            self.chunk[1, n] = hap_array # assign known haplotype
                            self.chunk[2, n] = comp # assign its complement
                            if self.unique_haplo(comp): # put the complement in the known list if it's unique
                                self.known_hap.append(comp)
                            num_phased += 1
                            enter = True # enter the while loop again since we phased a genotype
                            break # stop looking through known_hap for this given genotype

        hap_out = np.full((self.K, 2*self.N), -1, dtype=int)
        for ii in range(self.N):
            hap_out[:, 2*ii] = self.chunk[1, ii]
            hap_out[:, 2*ii+1] = self.chunk[2, ii]
        # Return the two haplotype chunks and the number of genotypes that are fully phased. Transpose the output in order to match input chunk shape to ouptut chunk shape.
        return hap_out, num_phased
        

    # need to create a function that checks if a haplotype sequence is unique compared to all that is found in self.known_hap
    def unique_haplo(self, hap_seq):
        """ Looks through all of self.known_hap and returns True if hap_seq is not found in known_hap. Returns false otherwise. """
        temp = True
        for hap_array in self.known_hap:
            if np.array_equal(hap_array, hap_seq): # returns true if two arrays are equal
                temp = False # assign temp to current node and break out
                break
        return temp
