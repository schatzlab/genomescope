#Given a reference sequence, introduce a  given rate of heterozygosity and create Illumina reads based on the wanted coverage level, readlength, and read duplication level.

#!/usr/bin/python
import os
import os.path
import gzip
import sys
#modify to your need:
sys.path.append("/seq/schatz/fritz/kmers/simul/")
from chromosomeMutator import *
import math
import numpy
import random
#Command Line Parameters
# parameteranalysis.py k coverage readlength br genomefile
cl=1 #compression level
k=int(sys.argv[1])
coverage = int(sys.argv[2])
readlength= int(sys.argv[3])
br = int(sys.argv[4])
r = float(sys.argv[5])
genomefile = sys.argv[6]
error = float(sys.argv[7])
limit=300
#q = genomefile.rfinqd("/")
#genomename = genomefile[q+1:]
#q = genomename.rfind(".")
#genomename = genomename[:q]
genomename = genomefile
fileprefix = genomename+"_het"+str(r)+"_br"+str(br)+"_rl"+str(readlength)+"_cov"+str(coverage)+"_err"+str(error)
illuminafile=fileprefix+"_reads.fa.gz"
if os.path.isfile(illuminafile)!=True:
	f = open(genomefile)
	genomestring = "".join(f.readlines()[1:])
	genomestring = genomestring.replace("\n","")
	genomestring = genomestring.replace("a","A")
	genomestring = genomestring.upper()
	f.close()
	length = len(genomestring)
	mutatedstring=chromosomeMutator(genomestring,r)
	illumina_reads=gzip.open(illuminafile, 'wb', compresslevel=cl)
	x=1
	iterate = length * coverage / readlength
	print "Creating Illumina Reads..."
	while(iterate>0):
		start = random.randint(1,length-readlength)
		end = start + readlength
		selectset = random.random()
		if(selectset<.5):
			seqstring=genomestring[start:end]
		if(selectset>=.5):
			seqstring=mutatedstring[start:end]
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
				for i in range(0,readlength):
					chance=random.random()
					if chance <= error:
						currentbp=seqstring[i]
						num = random.randint(0,2)
						if currentbp == "A":
							newbp=A[num]
						elif currentbp == "C":
							newbp=C[num]
						elif currentbp == "G":
							newbp=G[num]
						else:
							newbp=T[num]
						seqstring=seqstring[:i]+newbp+seqstring[i+1:]
		##if biasnum!=0:
		##	for i in range(0,int(biasnum)):
				illumina_reads.write(">"+str(x)+"_"+str(i)+"\n")
				illumina_reads.write("%s\n" % seqstring)
		iterate-=biasnum
		x=x+1
	illumina_reads.close()
else:
	print "Reads Already Exist"

