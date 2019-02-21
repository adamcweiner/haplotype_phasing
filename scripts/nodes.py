""" Create a classes for genotype and haplotype fragments with useful member variables so that we can use pointers. """

class Gen:
    """ Create a genotype object that has pointers to its two haplotypes. """
    def __init__(self, genotype):
        """ Pass in the given genotype and have haplotypes point to None (nullptr). """
        self.g = genotype
        self.h1 = None
        self.h2 = None

    def add_hs(self, h1):
        """ Create a pointer to the first haplotype assigned to the genotype. """
        assert self.valid_hap(h1)
        self.h1 = h1
        self.h2 = self.complementary_hap(h1, self.g) # assign second haplotype by finding complement

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
