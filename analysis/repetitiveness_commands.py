import itertools

num_het = {1:0, 2:1, 3:2, 4:4, 5:6, 6:10} #full model
t = 0 #topology = 0 (full model)
genome_filename = "random_genome.fa"
k = 21
coverage = 15
read_length = 100
br = 0
err = 0.01
model = "Partition"

COMMANDS_file = open("REPETITIVENESS_COMMANDS", 'w')
for i in range(1,8):
  for p in [3,4,5,6]:
    COMMANDS_file.write("parallel -j 17 < REPETITIVENESS_COMMANDS"+str(i)+"_p"+str(p)+"\n")

COMMANDS_file.close()

for p in [3, 4, 5, 6]:
  files = [open("REPETITIVENESS_COMMANDS" + str(i) + "_p" + str(p),'w') for i in range(1, 8)]
  lines = []
  for d in [0.01*x for x in range(51)]: #repetitiveness from 0% to 50% in 1.0% increments
    rs = [[0.02/(num_het[p]) for _ in range(num_het[p])]] #heterozygosity 2.0%, equally divided among each het form
    for r in rs:
      lines.append("python parameteranalysis.py " + str(k) + " " + str(p) + " " + str(coverage) + " " + str(read_length) + " " + str(br) + " " + ' '.join(['%.6f' % x for x in r]) + " " + genome_filename + " " + str(err) + " " + model + " " + '%.3f' % d + " " + str(t) + "\n")
      readfilenameprefix = genome_filename + "_k" + str(k) + "_p" + str(p) + "_cov" + str(coverage) + "_rl" + str(read_length) + "_br" + str(br) + "_het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(float(t)) + "_" + model + "_reads"
      lines.append("jellyfish count -C -m21 -s 1000000000 -t15 <(zcat " + readfilenameprefix + ".fa.gz) -o " + readfilenameprefix + ".jf" + "\n")
      lines.append("jellyfish histo -t15 " + readfilenameprefix + ".jf > " + readfilenameprefix + ".histo" + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 0" + " -o REPETITIVENESS_OUTPUT_p" + str(p) + "_transform_exp0 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 1" + " -o REPETITIVENESS_OUTPUT_p" + str(p) + "_transform_exp1 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 2" + " -o REPETITIVENESS_OUTPUT_p" + str(p) + "_transform_exp2 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      lines.append("genomescope.R -i " + readfilenameprefix + ".histo -k" + str(k) + " -p" + str(p) + " -l10 -t 0 --transform_exp 3" + " -o REPETITIVENESS_OUTPUT_p" + str(p) + "_transform_exp3 -n " + "het" + '_'.join(['%.6f' % x for x in r]) + "_err" + '%.3f' % err + "_d" + '%.3f' % d + "_top" + str(t) + " --testing --true_params=" + '%.3f' % d + "," + ",".join(['%.6f' % x for x in r]) + "," + str(t) + "\n")
      for i in range(7):
        files[i].write(lines[i])

  for i in range(7):
    files[i].close()
