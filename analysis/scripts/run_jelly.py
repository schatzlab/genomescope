#!/usr/bin/python2
import os
import os.path
import sys
illuminafile = sys.argv[1]
k=int(sys.argv[2])
limit=1000000

histfilename = illuminafile+sys.argv[2]+".hist"
outputname = illuminafile+sys.argv[2]+"_output"
mergename = illuminafile+sys.argv[2]+"_merge"

os.system("jellyfish count -m %s -o %s -c 3 -C -s 10000000000 -t 15 %s" %(k, outputname,illuminafile))
        #More than one output file, merge them
if os.path.isfile(outputname+"_1"):
	print "More than one output file detected, merging with jellyfish..."
	os.system("jellyfish merge -o %s %s_*"%(mergename,outputname))	
	print "Creating merged histogram file with jellyfish..."
	os.system("jellyfish histo -h %s -o %s %s"%(limit,histfilename,mergename))
        #Otherwise, make the histogram with only the first output file                  
else:
	print "Only one output file detected, creating histogram file with jellyfish..."        
	os.system("jellyfish histo -h %s -o %s %s"%(limit,histfilename,outputname))

