#Given a reference sequence, introduce a  given rate of heterozygosity and create Illumina reads based on the wanted coverage level, readlength, and read duplication level.

#!/usr/bin/python
import os
import os.path
import gzip
import sys
#modify to your need:
#sys.path.append("/seq/schatz/fritz/kmers/simul/")
from chromosomeMutator import *
import math
import numpy
import random
#Command Line Parameters
# parameteranalysis.py k p coverage readlength br genomefile
cl=1 #compression level
k=int(sys.argv[1])
p=int(sys.argv[2])
coverage = int(sys.argv[3])
readlength= int(sys.argv[4])
br = int(sys.argv[5])
r = float(sys.argv[6])
genomefile = sys.argv[7]
error = float(sys.argv[8])
limit=300
#q = genomefile.rfinqd("/")
#genomename = genomefile[q+1:]
#q = genomename.rfind(".")
#genomename = genomename[:q]
genomename = genomefile
fileprefix = genomename+"_p"+str(p)+"_het"+str(r)+"_br"+str(br)+"_rl"+str(readlength)+"_cov"+str(coverage)+"_err"+str(error)
illuminafile=fileprefix+"_reads.fa.gz"
if os.path.isfile(illuminafile)!=True:
	f = open(genomefile)
	genomestring = "".join(f.readlines()[1:])
	#print(genomestrings[0])
	genomestring = genomestring.replace("\n","")
	#genomestring = genomestring.replace("a","A")
	genomestring = genomestring.upper()
	f.close()
	length = len(genomestring)
	print(length)
	#for i in range(1,p):
	#	genomestrings[i]=chromosomeMutator(genomestrings[0],r)
	genomestrings = chromosomeMutator(genomestring,r,p)
	illumina_reads=gzip.open(illuminafile, 'wb', compresslevel=cl)
	x=1
	iterate = length * coverage / readlength
	print "Creating Illumina Reads..."
	while(iterate>0):
		start = random.randint(1,length-readlength)
		end = start + readlength
		selectset = random.randint(0,p-1)
		seqstring=genomestrings[selectset][start:end]
		biasnum=numpy.random.poisson(br)
		biasnum=math.ceil(biasnum)
		#Simulate error, each base in the sequence string has a possibility of error.
		#For each base in the read, generate a number, if it's less than the error alter that base
		read=seqstring
		for i in range(0,int(biasnum)):
			seqstring=read
			if error!=0:
				A="CGT"
				C="AGT"
				G="ACT"
				T="ACG"
				for j in range(0,readlength):
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
				illumina_reads.write(">"+str(x)+"_"+str(i)+"\n")
				illumina_reads.write("%s\n" % seqstring)
		iterate-=biasnum
		x+=1
	illumina_reads.close()
else:
	print "Reads Already Exist"
