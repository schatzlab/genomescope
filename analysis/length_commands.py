import itertools

num_het = {1:0, 2:1, 3:2, 4:4, 5:6, 6:10} #full model
t = 0 #topology = 0 (full model)
k = 21
coverage = 15
read_length = 100
br = 0
err = 0.01
model = "Partition"
d = 0.10

COMMANDS_file = open("LENGTH_COMMANDS", 'w')
for i in range(1,8):
  for p in [3,4,5,6]:
    COMMANDS_file.write("parallel -j 4 < LENGTH_COMMANDS"+str(i)+"_p"+str(p)+"\n")

for p in [3, 4, 5, 6]:
  files = [open("LENGTH_COMMANDS" + str(i) + "_p" + str(p),'w') for i in range(1, 8)]
  lines = []
  for genome_filename in ["random_genome_001Mbp.fa", "random_genome_010Mbp.fa", "random_genome_100Mbp.fa", "random_genome_001Gbp.fa"]:
    rs = [[0.001/(num_het[p]) for _ in range(num_het[p])]] #heterozygosity 2.0%, equally divided among each het form
    for r in rs:
      lines.append("python parameteranalysis.py " + str(k) + " " + str(p) + " " + str(coverage) + " " + str(read_length) + " " + str(br) + " " + ' '.join(['%.6f' % x for x in r]) + " " + genome_filename + " " + str(err) + " " + model + " " + '%.3f' % d + " " + str(t) + "\n")
      readfilenameprefix = genome_filename + "_k" + str(k) + "_p" + str(p) + "_cov" + str(coverage) + "_rl" + str(read_length) + "_br" + str(br) + "_het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(float(t)) + "_" + model + "_reads"
      lines.append("jellyfish count -C -m21 -s 1000000000 -t15 <(zcat " + readfilenameprefix + ".fa.gz) -o " + readfilenameprefix + ".jf" + "\n")
      lines.append("jellyfish histo -t15 " + readfilenameprefix + ".jf > " + readfilenameprefix + ".histo" + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 0" + " -o LENGTH_OUTPUT_p" + str(p) + "_transform_exp0 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 1" + " -o LENGTH_OUTPUT_p" + str(p) + "_transform_exp1 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 2" + " -o LENGTH_OUTPUT_p" + str(p) + "_transform_exp2 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 3" + " -o LENGTH_OUTPUT_p" + str(p) + "_transform_exp3 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      for i in range(7):
        files[i].write(lines[i])

  for i in range(7):
    files[i].close()
