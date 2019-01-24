#def chromosomeMutator(chromosome,heterozygosity):
#	import random
#	A="CGT"
#	C="AGT"
#	G="ACT"
#	T="ACG"
#	mutlist=list(chromosome)
#	for i in range(0,len(chromosome)-1):
#		chance=random.random()
#		if chance <= heterozygosity:
#			currentbp=mutlist[i]
#			num = random.randint(0,2)
#			if currentbp == "A":
#				newbp=A[num]
#			elif currentbp == "C":
#				newbp=C[num]
#			elif currentbp == "G":
#				newbp=G[num]
#			else:
#				newbp=T[num]
#			mutlist[i]=newbp
#	mutatedchromosome="".join(mutlist)
#	return mutatedchromosome
def chromosomeMutator(chromosome,heterozygosity,p):
	import random
	mutate = {"A": "ACGT", "C": "CAGT", "G": "GACT", "T": "TACG"}
	mutatedchromosomes = {}
	for j in range(p):
		mutatedchromosomes[j] = list(chromosome)
	for i in range(len(chromosome)):
		chance=random.random()
		if chance <= heterozygosity:
			nums = []
			while len(set(nums)) <= 1: #while all the same
				nums = [random.randint(0,3) for _ in range(p)]
			for j in range(p):
				mutatedchromosomes[j][i] = mutate[mutatedchromosomes[j][i]][nums[j]]
	for j in range(p):
		mutatedchromosomes[j]="".join(mutatedchromosomes[j])
	return mutatedchromosomes
