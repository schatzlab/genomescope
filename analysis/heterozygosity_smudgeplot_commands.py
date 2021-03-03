import itertools

num_het = {1:0, 2:1, 3:2, 4:4, 5:6, 6:10} #full model
t = 0 #topology = 0 (full model)
genome_filename = "random_genome_010Mbp.fa"
k = 21
coverage = 40
read_length = 100
br = 0
err = 0.01
model = "Partition"

COMMANDS_file = open("HETEROZYGOSITY_smudgeplot_COMMANDS", 'w')
for i in range(1,7):
  for p in [2,3,4,5,6]:
    if i in [1, 5]:
      COMMANDS_file.write("parallel -j 17 < HETEROZYGOSITY_smudgeplot_COMMANDS"+str(i)+"_p"+str(p)+"\n")
    else:
      COMMANDS_file.write("bash HETEROZYGOSITY_smudgeplot_COMMANDS"+str(i)+"_p"+str(p)+"\n")

COMMANDS_file.close()

p_to_rs = {2:[[0.005*x] for x in range(1,51)], 3:[[0.005*x+0.004, 0.001] for x in range(50)], 4:[[0, 0.005*x+0.003, 0.001, 0.001] for x in range(50)] + [[0.005*x+0.003, 0, 0.001, 0.001] for x in range(50)], 5:[[0.005*x+0.002, 0, 0.001, 0, 0.001, 0.001] for x in range(50)], 6:[[0.005*x+0.001, 0, 0, 0.001, 0, 0, 0.001, 0, 0.001, 0.001] for x in range(50)]}

for p in [2, 3, 4, 5, 6]:
  files = [open("HETEROZYGOSITY_smudgeplot_COMMANDS" + str(i) + "_p" + str(p),'w') for i in range(1, 7)]
  lines = []
  for d in [0.1]: #10% repetitiveness
    rs = p_to_rs[p] #heterozygosity from 0.5% to 25% in 0.5% increments
    for r in rs:
      lines.append("python parameteranalysis.py " + str(k) + " " + str(p) + " " + str(coverage) + " " + str(read_length) + " " + str(br) + " " + ' '.join(['%.6f' % x for x in r]) + " " + genome_filename + " " + str(err) + " " + model + " " + '%.3f' % d + " " + str(t) + "\n")
      readfilenameprefix = genome_filename + "_k" + str(k) + "_p" + str(p) + "_cov" + str(coverage) + "_rl" + str(read_length) + "_br" + str(br) + "_het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(float(t)) + "_" + model + "_reads"
      lines.append("kmc -k" + str(k) + " -m256 -ci1 -cs10000 -fm " + readfilenameprefix + ".fa.gz " + readfilenameprefix + " tmp/\n")
      lines.append("kmc_tools transform " + readfilenameprefix + " histogram " + readfilenameprefix + ".histo\n")
      lines.append("kmc_tools transform " + readfilenameprefix + " -ci10 -cx300 reduce " + readfilenameprefix + "_L10_U300\n")
      lines.append("smudge_pairs " + readfilenameprefix + "_L10_U300 " + readfilenameprefix + "_L10_U300_coverages.tsv\n")
      lines.append("smudgeplot.py plot " + readfilenameprefix + "_L10_U300_coverages.tsv -o " + readfilenameprefix + "_L10_U300_coverages\n")
      for i in range(6):
        files[i].write(lines[i])
      lines = []

  for i in range(6):
    files[i].close()
