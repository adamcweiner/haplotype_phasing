import numpy as np
from itertools import product

def em_algorithm2(data):
	"""
	Input: a matrix of {0,1,-1} n x 2m size, half phased data
	Output: a matrix of {0,1} n x 2m size
	"""
	n_snp, m_2ind = data.shape
	m_ind = m_2ind//2
	nequal = np.array([data[:,2*i]!=data[:,2*i+1] for i in range(m_ind)]).T
	unphased = np.array([data[:,2*i]==-1 for i in range(m_ind)]).T

	## Initialization: put all possible haplotypes in one pool
	pool = []  		# the pool of all haplotypes
	poolprob = []	# list of prob of each item in the pool
	poolptr = []	# list of pointers to the list of individuals
	indprob = []	# list of prob for each phasing solutions by hash value
	indptr = []		# list of pointers from ind phasing solutions to pool

	def pool_add(haplo, ind, hashind):
		for i in range(len(pool)):
			if np.array_equal( haplo, pool[i] ):
				poolptr[i].append((ind,hashind))
				indptr[ind][hashind].append(i)
				return
		pool.append(haplo)
		poolprob.append(1)
		poolptr.append([(ind,hashind)])
		indptr[ind][hashind].append(len(pool)-1)

	for ind in range(m_ind):
		ind_geno1 = data[:,ind*2]
		ind_geno2 = data[:,ind*2+1]
		unphased_spot = unphased[:,ind]
		no_nequal_spot = np.sum(nequal[:,ind]) <= 0
		n_ones = unphased_spot.sum()
		if no_nequal_spot:
			n_sol = 1 if n_ones<=1 else 2**(n_ones-1)
		else:
			n_sol = 1 if n_ones<=0 else 2**n_ones
		indprob.append([1/n_sol] * n_sol)
		indptr.append([[] for _ in range(n_sol)])
		seq1 = np.array([1 if i >= 1 else 0 for i in ind_geno1])
		seq2 = np.array([1 if i >= 1 else 0 for i in ind_geno2])
		lst = list(product([0,1], repeat=n_ones-1 if no_nequal_spot else n_ones))
		for idx, comb in enumerate(lst):
			if no_nequal_spot:
				comb_lst = np.array([0] + list(comb))
				compl_lst = [1]*n_ones-comb_lst
			else:
				comb_lst = np.array(list(comb))
				compl_lst = [1]*n_ones-comb_lst
			seq1[unphased_spot] = comb_lst
			seq2[unphased_spot] = compl_lst
			pool_add(seq1.copy(), ind, idx)
			pool_add(seq2.copy(), ind, idx)
	poolprob = poolprob / np.sum(poolprob)

	## Iterations
	num_iter = 0
	while num_iter < 50:
		num_iter += 1
		## Expectation
		for ipool in range(len(poolprob)):
			psum = 0
			for isol in range(len(poolptr[ipool])):
				ind, hashind = poolptr[ipool][isol]
				psum += indprob[ind][hashind]
			poolprob[ipool] = psum
		poolprob = poolprob / np.sum(poolprob)
		## Maximization
		for iind in range(len(indprob)):
			for isol in range(len(indprob[iind])):
				p1 = indptr[iind][isol][0]
				if len(indptr[iind][isol])>1:
					p2 = indptr[iind][isol][1]
					indprob[iind][isol] *= poolprob[p1]*poolprob[p2]
				else:
					indprob[iind][isol] *= poolprob[p1]*poolprob[p1]
			indprob[iind] = indprob[iind]/sum(indprob[iind])
		## Test convergence
		indconv = [np.any((np.array(prob)>0.99)==True) for prob in indprob]
		if all(indconv):
			print('converge at', num_iter)
			break
	
	## Output phased genotypes
	res = []
	for ind in range(m_ind):
		pos = np.argmax(indprob[ind])
		chrom1 = pool[indptr[ind][pos][0]]
		chrom2 = pool[indptr[ind][pos][1]]
		res.append(np.vstack((chrom1.copy(), chrom2.copy())).T)		
	return np.hstack(res)
