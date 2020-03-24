#Given a reference sequence, introduce a given rate of heterozygosity and create Illumina reads based on the wanted coverage level, readlength, and read duplication level.

#!/usr/bin/python
import os
import os.path
import gzip
import sys
#modify to your need:
#sys.path.append("/seq/schatz/fritz/kmers/simul/")
from chromosomeMutator import *
import math
import random
mutate = {"A": "CGT", "C": "GTA", "G": "TAC", "T": "ACG"}

num_params = {1:0, 2:1, 3:2, 4:4, 5:6, 6:10} #full model
cl=1 #compression level
k=int(sys.argv[1])
p=int(sys.argv[2])
coverage = int(sys.argv[3])
readlength= int(sys.argv[4])
br = int(sys.argv[5])
r = [float(sys.argv[i]) for i in range(6, 6+num_params[p])]
genomefile = sys.argv[6+num_params[p]]
error = float(sys.argv[7+num_params[p]])
topology_model = sys.argv[8+num_params[p]]
d = float(sys.argv[9+num_params[p]])
top = float(sys.argv[10+num_params[p]])
limit=300

genomename = genomefile
fileprefix = genomename+"_k"+str(k)+"_p"+str(p)+"_cov"+str(coverage)+"_rl"+str(readlength)+"_br"+str(br)+"_het"+'_'.join(['%.6f' % x for x in r])+"_err"+'%.3f' % error +"_d"+ '%.3f' % d + "_top" + str(top)
illuminafile=fileprefix+"_" + topology_model + "_reads.fa.gz"
if os.path.isfile(illuminafile)!=True:
	f = open(genomefile)
	genomestring = "".join(f.readlines()[1:]) #assumes has name as first line, and only one sequence
	#print(genomestrings[0])
	genomestring = genomestring.replace("\n","")
	#genomestring = genomestring.replace("a","A")
	genomestring = genomestring.upper()
	f.close()
	length = len(genomestring)
	genomestring = genomestring + genomestring[:int(length*d) + 1] #append sequence to make it repetitive
	length = len(genomestring)
	print(length)
	#for i in range(1,p):
	#	genomestrings[i]=chromosomeMutator(genomestrings[0],r)
	genomestrings = chromosomeMutator(genomestring,r,p, topology_model, top) #genomestrings is a dictionary from i to i-th (homologous set of) chromosome(s) (i is 1-indexed)
	illumina_reads=gzip.open(illuminafile, 'wb', compresslevel=cl)
	x=1
	iterate = p*length * coverage / readlength #num_reads
	print("Creating Illumina Reads...")
	while(iterate>0):
		start = random.randint(0,length-readlength)
		end = start + readlength
		selectset = random.randint(1,p)
		seqstring=genomestrings[selectset][start:end]
		#biasnum=numpy.random.poisson(br)
		#biasnum=math.ceil(biasnum)
		biasnum = 1 #added to get rid of bias
		#Simulate error, each base in the sequence string has a possibility of error.
		#For each base in the read, generate a number, if it's less than the error alter that base
		read=seqstring
		for i in range(0,int(biasnum)):
			seqstring=read
			if error!=0:
				A="CGT"
				C="GTA"
				G="TAC"
				T="ACG"
				for j in range(readlength):
					chance=random.random()
					if chance <= error:
						currentbp=seqstring[j]
						num = random.randint(0,2)
						if currentbp == "A":
							newbp=A[num]
						elif currentbp == "C":
							newbp=C[num]
						elif currentbp == "G":
							newbp=G[num]
						else:
							newbp=T[num]
						seqstring=seqstring[:j]+newbp+seqstring[j+1:]
		##if biasnum!=0:
		##	for i in range(0,int(biasnum)):
			illumina_reads.write((">"+str(x)+"_"+str(i)+"\n").encode()) #dedented 1
			illumina_reads.write(("%s\n" % seqstring).encode()) #dedented 1
		iterate-=biasnum
		x+=1
	illumina_reads.close()
else:
	print("Reads Already Exist")

