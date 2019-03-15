import numpy as np
from itertools import product

def em_algorithm3(data):
	""" 
	A more efficient EM Algorithm
	Input: n x 2m of {1,0,-1} half-phased (use the tool one2two() to convert {0,1,2} to this format)
	Output: n x 2m of {1,0} phased
	"""

	n_snp, m_2ind = data.shape
	m_ind = m_2ind//2
	
	pool = [] 		# pool of haplotypes, shape: ? x n, initial content placeholder
	poolprob = []	# list of prob of each item in the pool
	poolptr = []	# list of pointers to the list of individuals
	indprob = [[] for _ in range(m_ind)]	# list of prob for each phasing solutions by hash value
	indptr = [[] for _ in range(m_ind)]		# list of pointers from ind phasing solutions to pool

	est_unphased = (data == -1)[:,::2]
	est_diff = data[:,::2] != data[:,1::2]
	ind_unph_sum = np.sum(est_unphased, axis=0)
	ind_diff_sum = np.sum(est_diff, axis=0)
	order = np.argsort(ind_unph_sum)

	def compatible(haplo, chrom):
		if np.any(haplo[chrom==1]!=1):
			return False
		if np.any(haplo[chrom==0]!=0):
			return False
		return True

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

	for idx in order[range(m_ind)]:
		unph_sum = ind_unph_sum[idx]
		diff_sum = ind_diff_sum[idx]
		unphased_spot = est_unphased[:,idx]
		chrom1 = data[:,idx*2]
		chrom2 = data[:,idx*2+1]

		if unph_sum <= 0:
			indprob[idx] = [1]
			indptr[idx] = [[]]
			pool_add(chrom1,idx,0)
			if diff_sum>0:
				pool_add(chrom2,idx,0)
			continue

		compat_chrom1 = [compatible(entry, chrom1) for entry in pool]
		compat_chrom2 = [compatible(entry, chrom2) for entry in pool]
		if np.sum(compat_chrom1)<=0 and np.sum(compat_chrom2)<=0:
			# find no compatible, add all possibilities to the pool
			n_solg2 = unph_sum-1 if diff_sum <= 0 else unph_sum
			indprob[idx] = [1/2**n_solg2] * 2**n_solg2
			indptr[idx] = [[] for _ in range(2**n_solg2)]
			fill = list(product([0,1], repeat=n_solg2))
			for hashind, comb in enumerate(fill):
				comb_lst = np.array([0]+list(comb)) if diff_sum <= 0 else np.array(list(comb))
				compl_lst = [1]*unph_sum-comb_lst
				chrom1[unphased_spot] = comb_lst
				chrom2[unphased_spot] = compl_lst
				pool_add(chrom1.copy(), idx, hashind)
				pool_add(chrom2.copy(), idx, hashind)
		else:
			# find compatible, make sure their complements are also in the pool
			fill1 = np.array([pool[i][unphased_spot] for i in range(len(pool)) if compat_chrom1[i]])
			fill1 = fill1.reshape(-1,unph_sum)
			fill1c = [1]*unph_sum - fill1
			fill2 = np.array([pool[i][unphased_spot] for i in range(len(pool)) if compat_chrom2[i]])
			fill2 = fill2.reshape(-1,unph_sum)
			fill2c = [1]*unph_sum - fill2
			compat_pool = np.block([[fill1,fill1c],[fill2c,fill2]])
			compat_pool = np.unique(compat_pool, axis=0)
			n_sol = compat_pool.shape[0]
			if diff_sum <= 0:
				n_sol = n_sol//2
				compat_pool = compat_pool[:n_sol,:]	
			indprob[idx] = [1/n_sol] * n_sol
			indptr[idx] = [[] for _ in range(n_sol)]
			for i in range(n_sol):
				comb_lst = compat_pool[i,0:unph_sum]
				compl_lst = compat_pool[i,unph_sum:]
				chrom1[unphased_spot] = comb_lst
				chrom2[unphased_spot] = compl_lst
				pool_add(chrom1.copy(), idx, i)
				pool_add(chrom2.copy(), idx, i)
	print('Pool size ', len(pool))

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
				if len(indptr[iind][isol]) >1:
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
		chrom2 = pool[indptr[ind][pos][1]] if len(indptr[ind][pos])>1  else pool[indptr[ind][pos][0]]
		res.append(np.vstack((chrom1.copy(), chrom2.copy())).T)
	return np.hstack(res)


def one2two(data):
	n_snp, m_ind = data.shape
	two = (data==2)*1 + (data==1)*-1
	out = [two[:,i//2] for i in range(m_ind*2)]
	return np.vstack(out).T
