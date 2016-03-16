def chromosomeMutator(chromosome,heterozygosity):
	import random
	A="CGT"
	C="AGT"
	G="ACT"
	T="ACG"
	mutlist=list(chromosome)
	for i in range(0,len(chromosome)-1):
		chance=random.random()
		if chance <= heterozygosity:
			currentbp=mutlist[i]
			num = random.randint(0,2)
			if currentbp == "A":
				newbp=A[num]
			elif currentbp == "C":
				newbp=C[num]
			elif currentbp == "G":
				newbp=G[num]
			else:
				newbp=T[num]
			mutlist[i]=newbp
	mutatedchromosome="".join(mutlist)
	return mutatedchromosome
