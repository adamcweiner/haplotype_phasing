import numpy as np
from itertools import product

def em_algorithm(data):
	"""
	Input: a matrix of {0,1,2} n x m size, n SNPs and m individuals
	Output: a matrix of {0,1} n x 2m size
	"""
	## Initialization: put all possible haplotypes in one pool
	n_snp, m_ind = data.shape
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
		ind_geno = data[:,ind]
		ones_idx = (ind_geno==1).reshape(-1)
		n_ones = ones_idx.sum()
		n_sol = 1 if n_ones<=1 else 2**(n_ones-1)
		indprob.append([1/n_sol] * n_sol)
		indptr.append([[] for _ in range(n_sol)])
		seq = np.array([1 if i >= 1 else 0 for i in ind_geno])
		lst = list(product([0,1], repeat=n_ones))
		for idx, comb in enumerate(lst):
			comb_lst = np.array(list(comb))
			seq[ones_idx] = comb_lst
			hashind = min(idx, 2**n_ones-1-idx)
			pool_add(seq.copy(), ind, hashind)
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
		ind_geno = data[:,ind]
		ones_idx = (ind_geno==1).reshape(-1)
		n_ones = ones_idx.sum()
		seq = np.array([1 if i >= 1 else 0 for i in ind_geno])
		lst = list(product([0,1], repeat=n_ones))
		seq[ones_idx] = lst[pos]
		chrom1 = seq.copy()
		seq[ones_idx] = lst[2**n_ones-1-pos]
		chrom2 = seq.copy()
		res.append(np.vstack((chrom1, chrom2)).T)
	return np.hstack(res)
