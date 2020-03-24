#python make_random_genome.py length output

#requires Python 3.6 or higher
import random
import sys

length = int(sys.argv[1])
output = sys.argv[2]

f = open(output, 'w')
num_lines_before_last = length//50
num_chars_last_line = length - num_lines_before_last*50

for i in range(num_lines_before_last):
  f.write("".join(random.choices("actg", k=50))+"\n")

f.write("".join(random.choices("actg", k=num_chars_last_line)))
f.close()
