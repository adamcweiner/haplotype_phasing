""" Create a classes for genotype and haplotype fragments with useful member variables so that we can use pointers. """

class Gen:
    """ Create a genotype object that has pointers to its two haplotypes. """
    def __init__(self, gen_seq):
        """ Pass in the given genotype and have haplotypes point to None (nullptr). """
        self.seq = gen_seq
        self.h1 = None
        self.h2 = None

    def add_hs(self, h1):
        """ Adds the given haplotype and assigns the other haplotype haplotype assigned to the genotype. Note that h1 is a Hap object. """
        assert self.valid_hap(h1.seq)
        self.h1 = h1
        h2_seq = self.complementary_hap(h1.seq) # find sequence of second haplotype by finding complement
        self.h2 = Hap(h2_seq) # create Hap object given the proper sequence

    def complementary_hap(self, hap_seq):
        """ Given genotype and a haplotype sequences, find the haplotype sequence that complements the one given. """
        assert len(hap_seq) == len(self.seq) # make sure the two arrays are of the same size
        comp_hap = np.zeros((len(hap_seq))) # fill comp_hap with just zeros to start
        for ii in range(len(hap_seq)):
            comp_hap[ii] = self.seq[ii] - hap_seq[ii] # i.e. if g[ii] = 2 and h[ii] = 1 then c_h[ii] = 1
        return comp_hap

    def is_homogeneous(self):
        """ Returns True if the genotype only consists of 0's and 2's; returns False if any of the positions are heterogeneous. """
        temp = True
        for ii in range(len(self.seq)):
            if self.seq[ii] == 1: # if element is 1, break and return False
                temp = False
                break
        return temp

class Hap:
    """ Creates a haplotype object that has a pointer to its genotype. """
    def __init__(self, hap_seq):
        """ Holds a sequence for the haplotype object. """
        assert self.valid_hap(hap_seq) # assert that the sequence is a valid haplotype
        self.seq = hap_seq

    def valid_hap(self, hap_seq):
        """ Returns True if the haplotype is valid (meaning all entries are just 0's and 1's)"""
        temp = True
        for ii in range(len(hap_seq)):
            if hap_seq[ii] > 1 or hap_seq[ii] < 0: # if element is not 0 or 1
                temp = False
                break
        return temp
