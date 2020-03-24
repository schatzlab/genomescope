def chromosomeMutator(chromosome,heterozygosities,p, topology_model, top): #"Partition", p=1 to p=6
  import random
  mutate = {"A": "ACGTAC", "C": "CGTACG", "G": "GTACGT", "T": "TACGTA"}
  mutatedchromosomes = {}
  length = len(chromosome)
  if topology_model == "Partition":
    if p == 1:
      mutatedchromosomes[1] = chromosome
    if p == 2:
      mutatedchromosomes[1] = list(chromosome)
      mutatedchromosomes[2] = list(chromosome)
      r = [1-heterozygosities[0]] + heterozygosities
      for i in range(length):
        chance = random.random()
        if chance <= r[0]:
          nums = [0, 0]
        else:
          nums = [0, 1]
        mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
        mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
      mutatedchromosomes[1] = "".join(mutatedchromosomes[1])
      mutatedchromosomes[2] = "".join(mutatedchromosomes[2])
    if p == 3:
      mutatedchromosomes[1] = list(chromosome)
      mutatedchromosomes[2] = list(chromosome)
      mutatedchromosomes[3] = list(chromosome)
      r = [1-heterozygosities[0]-heterozygosities[1]] + heterozygosities
      for i in range(length):
        chance = random.random()
        if chance <= r[0]:
          nums = [0, 0, 0] #how to mutate
        elif chance <= r[0] + r[1]:
          nums = [0, 0, 1]
        else:
          nums = [0, 1, 2]
        mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
        mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
        mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
      mutatedchromosomes[1] = "".join(mutatedchromosomes[1])
      mutatedchromosomes[2] = "".join(mutatedchromosomes[2])
      mutatedchromosomes[3] = "".join(mutatedchromosomes[3])
    if p == 4:
      mutatedchromosomes[1] = list(chromosome)
      mutatedchromosomes[2] = list(chromosome)
      mutatedchromosomes[3] = list(chromosome)
      mutatedchromosomes[4] = list(chromosome)
      if top != 0:
        r = [1 - heterozygosities[0] - heterozygosities[1] - heterozygosities[2]] + heterozygosities
        for i in range(length):
          chance = random.random()
          if chance <= r[0]:
            nums = [0, 0, 0, 0]
          elif chance <= r[0] + r[1]:
            if top == 1:
              nums = [0, 0, 0, 1]
            if top == 2:
              nums = [0, 0, 1, 1]
          elif chance <= r[0] + r[1] + r[2]:
            nums = [0, 0, 1, 2]
          else:
            nums = [0, 1, 2, 3]
          mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
          mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
          mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
          mutatedchromosomes[4][i] = mutate[mutatedchromosomes[4][i]][nums[3]]
      else:
        r = [1 - heterozygosities[0] - heterozygosities[1] - heterozygosities[2] - heterozygosities[3]] + heterozygosities
        for i in range(length):
          chance = random.random()
          if chance <= r[0]:
            nums = [0, 0, 0, 0]
          elif chance <= r[0] + r[1]:
            nums = [0, 0, 0, 1]
          elif chance <= r[0] + r[1] + r[2]:
            nums = [0, 0, 1, 1]
          elif chance <= r[0] + r[1] + r[2] + r[3]:
            nums = [0, 0, 1, 2]
          else:
            nums = [0, 1, 2, 3]
          mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
          mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
          mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
          mutatedchromosomes[4][i] = mutate[mutatedchromosomes[4][i]][nums[3]]
      mutatedchromosomes[1] = "".join(mutatedchromosomes[1])
      mutatedchromosomes[2] = "".join(mutatedchromosomes[2])
      mutatedchromosomes[3] = "".join(mutatedchromosomes[3])
      mutatedchromosomes[4] = "".join(mutatedchromosomes[4])
    if p == 5:
      mutatedchromosomes[1] = list(chromosome)
      mutatedchromosomes[2] = list(chromosome)
      mutatedchromosomes[3] = list(chromosome)
      mutatedchromosomes[4] = list(chromosome)
      mutatedchromosomes[5] = list(chromosome)
      if top!=0:
        r = [1 - heterozygosities[0] - heterozygosities[1] - heterozygosities[2] - heterozygosities[3]] + heterozygosities
        for i in range(length):
          chance = random.random()
          if chance <= sum(r[:1]):
            nums = [0, 0, 0, 0, 0]
          elif chance <= sum(r[:2]):
            if top in [1, 2]:
              nums = [0, 0, 0, 0, 1]
            if top in [3, 4, 5]:
              nums = [0, 0, 0, 1, 1]
          elif chance <= sum(r[:3]):
            if top in [1,3]:
              nums = [0, 0, 0, 1, 2]
            if top in [2]:
              nums = [0, 0, 1, 1, 2]
            if top in [4, 5]:
              nums = [0, 0, 1, 2, 2]
          elif chance <= sum(r[:4]):
            if top in [1, 2, 3, 4]:
              nums = [0, 0, 1, 2, 3]
            if top in [5]:
              nums = [0, 1, 2, 3, 3]
          else:
            nums = [0, 1, 2, 3, 0]
          mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
          mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
          mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
          mutatedchromosomes[4][i] = mutate[mutatedchromosomes[4][i]][nums[3]]
          mutatedchromosomes[5][i] = mutate[mutatedchromosomes[5][i]][nums[4]]
      else:
        r = [1 - sum(heterozygosities)] + heterozygosities
        for i in range(length):
          chance = random.random()
          if chance <= sum(r[:1]):
            nums = [0, 0, 0, 0, 0]
          elif chance <= sum(r[:2]):
            nums = [0, 0, 0, 0, 1]
          elif chance <= sum(r[:3]):
            nums = [0, 0, 0, 1, 1]
          elif chance <= sum(r[:4]):
            nums = [0, 0, 0, 1, 2]
          elif chance <= sum(r[:5]):
            nums = [0, 0, 1, 1, 2]
          elif chance <= sum(r[:6]):
            nums = [0, 0, 1, 2, 3] #nums = [0, 0, 1, 2, 2]
          #elif chance <= sum(r[:7]):
          #  nums = [0, 0, 1, 2, 3]
          #elif chance <= sum(r[:8]):
          #  nums = [0, 1, 2, 3, 3]
          else:
            nums = [0, 1, 2, 3, 4]
          mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
          mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
          mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
          mutatedchromosomes[4][i] = mutate[mutatedchromosomes[4][i]][nums[3]]
          mutatedchromosomes[5][i] = mutate[mutatedchromosomes[5][i]][nums[4]]
      mutatedchromosomes[1] = "".join(mutatedchromosomes[1])
      mutatedchromosomes[2] = "".join(mutatedchromosomes[2])
      mutatedchromosomes[3] = "".join(mutatedchromosomes[3])
      mutatedchromosomes[4] = "".join(mutatedchromosomes[4])
      mutatedchromosomes[5] = "".join(mutatedchromosomes[5])
    if p == 6:
      mutatedchromosomes[1] = list(chromosome)
      mutatedchromosomes[2] = list(chromosome)
      mutatedchromosomes[3] = list(chromosome)
      mutatedchromosomes[4] = list(chromosome)
      mutatedchromosomes[5] = list(chromosome)
      mutatedchromosomes[6] = list(chromosome)
      if top!=0:
        r = [1 - heterozygosities[0] - heterozygosities[1] - heterozygosities[2] - heterozygosities[3] - heterozygosities[4]] + heterozygosities
        for i in range(length):
          chance = random.random()
          if chance <= sum(r[:1]):
            nums = [0, 0, 0, 0, 0, 0]
          elif chance <= sum(r[:2]):
            if top in [1, 2, 3, 4, 5]:
              nums = [0, 0, 0, 0, 0, 1]
            if top in [6, 7, 8, 9, 10, 11, 12, 13]:
              nums = [0, 0, 0, 0, 1, 1]
            if top in [14, 15, 16]:
              nums = [0, 0, 0, 1, 1, 1]
          elif chance <= sum(r[:3]):
            if top in [1, 2, 6, 7]:
              nums = [0, 0, 0, 0, 1, 2]
            if top in [3, 4, 5, 14, 15, 16]:
              nums = [0, 0, 0, 1, 1, 2]
            if top in [8, 9, 10]:
              nums = [0, 0, 0, 1, 2, 2]
            if top in [11, 12, 13]:
              nums = [0, 0, 1, 1, 2, 2]
          elif chance <= sum(r[:4]):
            if top in [1, 3, 6, 8, 14]:
              nums = [0, 0, 0, 1, 2, 3]
            if top in [2, 7, 11]:
              nums = [0, 0, 1, 1, 2, 3]
            if top in [4, 5, 15, 16]:
              nums = [0, 0, 1, 2, 2, 3]
            if top in [9, 10, 12, 13]:
              nums = [0, 0, 1, 2, 3, 3]
          elif chance <= sum(r[:5]):
            if top in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15]:
              nums = [0, 0, 1, 2, 3, 0]
            if top in [5, 16]:
              nums = [0, 1, 2, 3, 3, 0]
            if top in [10, 13]:
              nums = [0, 1, 2, 3, 0, 0]
          else:
            nums = [0, 1, 2, 3, 0, 1]
          mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
          mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
          mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
          mutatedchromosomes[4][i] = mutate[mutatedchromosomes[4][i]][nums[3]]
          mutatedchromosomes[5][i] = mutate[mutatedchromosomes[5][i]][nums[4]]
          mutatedchromosomes[6][i] = mutate[mutatedchromosomes[6][i]][nums[5]]
      else:
        r = [1 - sum(heterozygosities)] + heterozygosities
        for i in range(length):
          chance = random.random()
          if chance <= sum(r[:1]):
            nums = [0, 0, 0, 0, 0, 0]
          elif chance <= sum(r[:2]):
            nums = [0, 0, 0, 0, 0, 1]
          elif chance <= sum(r[:3]):
            nums = [0, 0, 0, 0, 1, 1]
          elif chance <= sum(r[:4]):
            nums = [0, 0, 0, 1, 1, 1]
          elif chance <= sum(r[:5]):
            nums = [0, 0, 0, 0, 1, 2]
          elif chance <= sum(r[:6]):
            nums = [0, 0, 0, 1, 1, 2]
          elif chance <= sum(r[:7]):
            nums = [0, 0, 1, 1, 2, 2] #nums = [0, 0, 0, 1, 2, 2]
          elif chance <= sum(r[:8]):
            nums = [0, 0, 0, 1, 2, 3] #nums = [0, 0, 1, 1, 2, 2]
          elif chance <= sum(r[:9]):
            nums = [0, 0, 1, 1, 2, 3] #nums = [0, 0, 0, 1, 2, 3]
          elif chance <= sum(r[:10]):
            nums = [0, 0, 1, 2, 3, 4] #nums = [0, 0, 1, 1, 2, 3]
          #elif chance <= sum(r[:11]):
          #  nums = [0, 0, 1, 2, 2, 3]
          #elif chance <= sum(r[:12]):
          #  nums = [0, 0, 1, 2, 3, 3]
          #elif chance <= sum(r[:13]):
          #  nums = [0, 0, 1, 2, 3, 4]
          #elif chance <= sum(r[:14]):
          #  nums = [0, 1, 2, 3, 3, 4]
          #elif chance <= sum(r[:15]):
          #  nums = [0, 1, 2, 3, 4, 4]
          else:
            nums = [0, 1, 2, 3, 4, 5]
          mutatedchromosomes[1][i] = mutate[mutatedchromosomes[1][i]][nums[0]]
          mutatedchromosomes[2][i] = mutate[mutatedchromosomes[2][i]][nums[1]]
          mutatedchromosomes[3][i] = mutate[mutatedchromosomes[3][i]][nums[2]]
          mutatedchromosomes[4][i] = mutate[mutatedchromosomes[4][i]][nums[3]]
          mutatedchromosomes[5][i] = mutate[mutatedchromosomes[5][i]][nums[4]]
          mutatedchromosomes[6][i] = mutate[mutatedchromosomes[6][i]][nums[5]]
      mutatedchromosomes[1] = "".join(mutatedchromosomes[1])
      mutatedchromosomes[2] = "".join(mutatedchromosomes[2])
      mutatedchromosomes[3] = "".join(mutatedchromosomes[3])
      mutatedchromosomes[4] = "".join(mutatedchromosomes[4])
      mutatedchromosomes[5] = "".join(mutatedchromosomes[5])
      mutatedchromosomes[6] = "".join(mutatedchromosomes[6])  
  return mutatedchromosomes
