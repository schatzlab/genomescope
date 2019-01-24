#python project.py GENOME_LEN DUP_RATE HET_RATE KMER_LEN PLOIDY
from collections import Counter
from collections import defaultdict
from itertools import product
from itertools import permutations
import math
import networkx as nx
import numpy as np
import os.path
import pickle
import scipy.special
import sys

def get_repeatpartitiongroups(p):
	repeatpartitiongroups = defaultdict(list)
	for partition_length_left in range(1, p+1):
		for partition_length_right in range(1, p+1):
			for mapping_length in range(min(partition_length_left, partition_length_right)+1):
				size = partition_length_left + partition_length_right - mapping_length
				repeatpartitiongroups[size].append(np.array([partition_length_left, mapping_length, partition_length_right]))
	return [value for (key,value) in sorted(repeatpartitiongroups.items(), reverse = True)]

def get_repeatpartitiongraph(repeatpartitiongroups):
	
	def is_parent_of(repeatpartition1, repeatpartition2):
		diff = repeatpartition1 - repeatpartition2
		if (any(np.array_equal(diff, x) for x in [[1,0,0], [0,-1,0], [0,0,1]])) or (np.array_equal(diff, [1,1,1]) and (repeatpartition1[1] > 1)):
			return True
		else:
			return False
	
	G = nx.DiGraph()
	for i in range(len(repeatpartitiongroups)-1):
		for repeatpartition1 in repeatpartitiongroups[i]:
			for repeatpartition2 in repeatpartitiongroups[i+1]:
				if is_parent_of(repeatpartition1, repeatpartition2):
					G.add_edge(tuple(repeatpartition1), tuple(repeatpartition2))
	return G

def get_descendants(repeatpartitiongroups, graph):
	descendants = {}
	for repeatpartitiongroup in repeatpartitiongroups:
		for repeatpartition in repeatpartitiongroup:
			descendants[tuple(repeatpartition)] = nx.descendants(graph, tuple(repeatpartition))
	return descendants

def dd():
	return defaultdict(int)

def get_mobius(repeatpartitiongroups, graph, descendants):
	
	def get_cardinality(repeatpartition1, repeatpartition2):
		
		def get_size(repeatpartition):
			size = repeatpartition[0] + repeatpartition[2] - repeatpartition[1]
			return size
		
		def get_k_setpartitions(collection, k):
			if k == 1:
				return [[collection]]
			if k == 0 or len(collection) < k:
				return []
			
			def visit(n, a):
				ps = [[] for i in xrange(k)]
				for j in xrange(n):
					ps[a[j + 1]].append(collection[j])
				return ps
			
			def f(mu, nu, sigma, n, a):
				if mu == 2:
					yield visit(n, a)
				else:
					for v in f(mu - 1, nu - 1, (mu + sigma) % 2, n, a):
						yield v
				if nu == mu + 1:
					a[mu] = mu - 1
					yield visit(n, a)
					while a[nu] > 0:
						a[nu] = a[nu] - 1
						yield visit(n, a)
				elif nu > mu + 1:
					if (mu + sigma) % 2 == 1:
						a[nu - 1] = mu - 1
					else:
						a[mu] = mu - 1
					if (a[nu] + sigma) % 2 == 1:
						for v in b(mu, nu - 1, 0, n, a):
							yield v
					else:
						for v in f(mu, nu - 1, 0, n, a):
							yield v
					while a[nu] > 0:
						a[nu] = a[nu] - 1
						if (a[nu] + sigma) % 2 == 1:
							for v in b(mu, nu - 1, 0, n, a):
								yield v
						else:
							for v in f(mu, nu - 1, 0, n, a):
								yield v
			
			def b(mu, nu, sigma, n, a):
				if nu == mu + 1:
					while a[nu] < mu - 1:
						yield visit(n, a)
						a[nu] = a[nu] + 1
					yield visit(n, a)
					a[mu] = 0
				elif nu > mu + 1:
					if (a[nu] + sigma) % 2 == 1:
						for v in f(mu, nu - 1, 0, n, a):
							yield v
					else:
						for v in b(mu, nu - 1, 0, n, a):
							yield v
					while a[nu] < mu - 1:
						a[nu] = a[nu] + 1
						if (a[nu] + sigma) % 2 == 1:
							for v in f(mu, nu - 1, 0, n, a):
								yield v
						else:
							for v in b(mu, nu - 1, 0, n, a):
								yield v
					if (mu + sigma) % 2 == 1:
						a[nu - 1] = 0
					else:
						a[mu] = 0
				if mu == 2:
					yield visit(n, a)
				else:
					for v in b(mu - 1, nu - 1, (mu + sigma) % 2, n, a):
						yield v
			
			n = len(collection)
			a = [0] * (n + 1)
			for j in xrange(1, k + 1):
				a[n - k + j] = j - 1
			return f(k, n, 0, n, a)
		
		def valid_setpartition(repeatpartition1, setpartition, repeatpartition2):
			partition_length_left1, mapping_length1, partition_length_right1 = repeatpartition1
			partition_length_left2, mapping_length2, partition_length_right2 = repeatpartition2
			
			left = range(partition_length_left1 - mapping_length1)
			across = range(partition_length_left1 - mapping_length1, partition_length_left1)
			right = range(partition_length_left1, partition_length_left1 + partition_length_right1 - mapping_length1)
			
			numleftparts = 0
			numacrossparts = 0
			numrightparts = 0
			
			for part in setpartition:
				numinleft = len(set(part).intersection(left))
				numinacross = len(set(part).intersection(across))
				numinright = len(set(part).intersection(right))
				if (numinleft > 0 and numinright > 0) or (numinacross > 0):
					numleftparts += 1
					numacrossparts += 1
					numrightparts += 1
				elif numinleft > 0:
					numleftparts += 1
				elif numinright > 0:
					numrightparts += 1
			
			if partition_length_left2 == numleftparts and mapping_length2 == numacrossparts and partition_length_right2 == numrightparts:
				return True
			else:
				return False
		
		return len([setpartition for setpartition in get_k_setpartitions(range(get_size(repeatpartition1)), get_size(repeatpartition2)) if valid_setpartition(repeatpartition1, setpartition, repeatpartition2)])
	
	filename = 'mobius'+str(p)
	if os.path.exists(filename):
		mobius = pickle.load(open(filename, 'rb'))
		return mobius
	mobius = defaultdict(dd)
	source = tuple(repeatpartitiongroups[0][0])
	for repeatpartition1 in nx.dfs_postorder_nodes(graph, source):
		mobius[repeatpartition1][repeatpartition1] = 1
		for repeatpartition2 in descendants[repeatpartition1]:
			for repeatpartition3 in mobius[repeatpartition2]:
				mobius[repeatpartition1][repeatpartition3] -= get_cardinality(repeatpartition1, repeatpartition2) * mobius[repeatpartition2][repeatpartition3]
	pickle.dump(mobius, open(filename,'wb'))
	return mobius

def get_parameters(r, d, k, p, mobius):
	
	def Theta(n):
		a = [0 for i in range(n + 1)]
		k = 1
		y = n - 1
		while k != 0:
		    x = a[k - 1] + 1
		    k -= 1
		    while 2 * x <= y:
		        a[k] = x
		        y -= x
		        k += 1
		    l = k + 1
		    while x <= y:
		        a[k] = x
		        a[l] = y
		        yield a[:k + 2]
		        x += 1
		        y -= 1
		    a[k] = x + y
		    y = x + y - 1
		    yield a[:k + 1]
	
	def Chi(partition1, partition2):
		mappings = []
		maxr = min(len(partition1), len(partition2))
		for r in range(maxr + 1):
			permutations1 = set(permutations(partition1, r))
			permutations2 = set(permutations(partition2, r))
			cartesian = product(permutations1, permutations2)
			mappings.extend(list(set([tuple(sorted(zip(*item))) for item in cartesian])))
		return mappings
	
	def psi(partition1, mapping, partition2):
		combined_partition = []
		remaining_partitions = partition1 + partition2
		for arrow in mapping:
			combined_partition.extend([sum(arrow)])
			for partition in arrow:
				remaining_partitions.remove(partition)
		combined_partition.extend(remaining_partitions)
		return combined_partition
	
	def U(partition):
		return np.array([Counter(partition)[i] for i in range(1, 2*sum(partition)+1)])
	
	def D(partition):
		return np.array([Counter(partition)[i] for i in range(1, sum(partition)+1)])
	
	def c(partition):
		n = sum(partition)
		num = math.factorial(n)
		counts = Counter(partition)
		den = np.prod([(math.factorial(x)**counts[x]) * math.factorial(counts[x]) for x in set(partition)])
		return num/den
	
	def c_prime(partition_left, mapping, partition_right):
		
		def fallingfactorial(x, n):
			return np.prod([x-k for k in range(n)])
		
		def getcounters(partition_left, mapping, partition_right):
			leftcounter = []
			leftmapcounter = []
			rightmapcounter = []
			rightcounter = []
			for arrow in mapping:
				leftmapcounter.append(arrow[0])
				rightmapcounter.append(arrow[1])
			leftcounter = Counter(partition_left)
			leftmapcounter = Counter(leftmapcounter)
			mapcounter = Counter(mapping)
			rightmapcounter = Counter(rightmapcounter)
			rightcounter = Counter(partition_right)
			return (leftcounter, leftmapcounter, mapcounter, rightmapcounter, rightcounter)
		
		leftcounter, leftmapcounter, mapcounter, rightmapcounter, rightcounter = getcounters(partition_left, mapping, partition_right)
		numwaysleftend = np.prod([fallingfactorial(leftcounter[part], leftmapcounter[part]) for part in leftmapcounter])
		numwaysrightend = np.prod([fallingfactorial(rightcounter[part], rightmapcounter[part]) for part in rightmapcounter])
		maprepeats = np.prod([math.factorial(value) for value in mapcounter.values()])
		return c(partition_left) * c(partition_right) * numwaysleftend * numwaysrightend / maprepeats
	
	def h(r, k, p, partition):
		
		def stirling1(n, k):
			if n == 0 and k == 0:
				return 1
			if (n > 0 and k == 0) or (n == 0 and k > 0):
				return 0
			else:
				return -(n-1)*stirling1(n-1,k)+stirling1(n-1,k-1)
		
		def f(r, k, p, partition_length):
			het_het_num = 4**(partition_length) - 4
			het_het_denom = 4**(p) - 4
			het_het = float(het_het_num)/het_het_denom
			assert het_het <= 1
			return (r*het_het + (1-r))**k
		
		partition_length = len(partition)
		return c(partition)*sum([stirling1(partition_length, partition_length_i) * f(r, k, p, partition_length_i) for partition_length_i in range(partition_length, 0, -1)])
	
	def h_prime(r, k, p, partition_left, mapping, partition_right, mobius):
		
		def f_prime(r, k, p, partition_length_left, mapping_length, partition_length_right):
			a = partition_length_left - mapping_length
			b = partition_length_right - mapping_length
			partition_prime_length = partition_length_left + partition_length_right - mapping_length
			if mapping_length == 0:
				het_het_num = 4**(partition_prime_length)-4*(4**a+4**b)+16
			else:
				het_het_num = 4**(partition_prime_length)-4*(4**a+4**b)+4
			het_het_denom = (4**(p)-4)*(4**(p)-4)
			het_het = float(het_het_num)/het_het_denom
			assert het_het <= 1
			if mapping_length == 0:
				het_homo_num_left = 4**a - 4
				het_homo_num_right = 4**b - 4
			else:
				het_homo_num_left = 4**a - 1
				het_homo_num_right = 4**b - 1
			het_homo_num = het_homo_num_left + het_homo_num_right
			het_homo_denom = 4**(p)-4
			het_homo = float(het_homo_num)/(het_homo_denom)
			assert het_homo_num_left <= het_homo_denom
			assert het_homo_num_right <= het_homo_denom
			return (r**2*(het_het) + r*(1-r)*(het_homo) + (1-r)**2)**k
		
		partition_length_left = len(partition_left)
		mapping_length = len(mapping)
		partition_length_right = len(partition_right)
		return c_prime(partition_left, mapping, partition_right)*sum([coeff_i * f_prime(r, k, p, partition_length_left_i, mapping_length_i, partition_length_right_i) for ((partition_length_left_i, mapping_length_i, partition_length_right_i), coeff_i) in mobius[(partition_length_left, mapping_length, partition_length_right)].items()])
	
	return (1-d)*sum([(U(partition) * h(r, k, p, partition)) for partition in Theta(p)]) + d*sum([sum([sum([(D(psi(partition1, mapping, partition2)) * h_prime(r, k, p, partition1, mapping, partition2, mobius)) for mapping in Chi(partition1, partition2)]) for partition2 in Theta(p)]) for partition1 in Theta(p)])

if __name__ == '__main__':
	G = int(sys.argv[1])
	d = float(sys.argv[2])
	r = float(sys.argv[3])
	k = int(sys.argv[4])
	p = int(sys.argv[5])
	
	repeatpartitiongroups = get_repeatpartitiongroups(p)
	graph = get_repeatpartitiongraph(repeatpartitiongroups)
	descendants = get_descendants(repeatpartitiongroups, graph)
	mobius = get_mobius(repeatpartitiongroups, graph, descendants)
	parameters = get_parameters(r, d, k, p, mobius)
	
	print('Model:')
	for i,x in enumerate(G*parameters):
		print(str(i+1)+' '+str(int(x)))
