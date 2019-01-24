#!/usr/bin/env Rscript

## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
##
## This is the automated script for computing genome characteristics
## from a kmer histogram file, k-mer size, and readlength

## Number of rounds before giving up
NUM_ROUNDS=4

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Max rounds on NLS
MAX_ITERATIONS=200

## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
SCORE_CLOSE = 0.20 #originally 0.20

## Overrule heterozygosity if there is a large difference in het rate
SCORE_HET_FOLD_DIFFERENCE = 10 #originally 10

## Print out VERBOSEging messages (0/1)
VERBOSE = 0

## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
options(warn=-1)

## Colors for plots
COLOR_BGCOLOR  = "light grey"
COLOR_HIST     = "#56B4E9"
COLOR_2pPEAK    = "black"
COLOR_pPEAK    = "#F0E442"
COLOR_ERRORS   = "#D55E00"
COLOR_KMERPEAK = "black"
COLOR_RESIDUAL = "purple"
COLOR_COVTHRES = "red"

library('minpack.lm')
library('argparse')

## Given mean +/- stderr, report min and max value within 2 SE
###############################################################################

min_max <- function(table){
  ##return (c( abs(table[1]) - 2*abs(table[2]) , abs(table[1])+ 2*abs(table[2])))
  return (c(max(0,table[1] - 2*table[2]), table[1]+ 2*table[2]))
}
min_max1 <- function(table){
  return (c(max(0,table[1] - 2*table[2]), min(1, table[1]+ 2*table[2])))
}

## Use nls to fit 2p peak model
###############################################################################

nls_peak<-function(x, y, k, p, estKmercov, estLength, max_iterations){
  model = NULL

  if (VERBOSE) { cat("trying nls_peak standard algorithm\n") }

  predict1 = function(k, d, kmercov, bias, x)
  {
    r0 = 1
    t0 = r0**k
    s0 = t0
    alpha_1 = (1-d)*(s0)
    alpha_2 = d*(s0)
    alpha_1 * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size=kmercov*2 / bias, mu = kmercov*2)
  }

#AB -> AA
  predict2 = function(r1, k, d, kmercov, bias, x)
  {
    r0 = 1-r1 #aa
    if (r0 < 0) {return(0)}
    t0 = r0**k #AA
    s0 = t0 #AA
    s1 = 1-t0 #AB
    alpha_1 = (1-d)*(2*s1) + d*(2*s0*s1 + 2*s1**2)
    alpha_2 = (1-d)*(s0) + d*(s1**2)
    alpha_3 = d*(2*s0*s1)
    alpha_4 = d*(s0**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
  }

#ABC -> AAB -> AAA
  predict3 = function(r1, r2, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2 #aaa
    if (r0 < 0) {return(0)}
    t0 = r0**k #AAA
    t1 = (r0+r1)**k #AAA + AAB
    s0 = t0 #AAA
    s1 = t1-t0 #AAB
    s2 = 1-t1 #ABC
    alpha_1 = (1-d)*(s1+3*s2) + d*(2*s0*s1 + 4*s0*s2 + 2*s1**2 + 6*s1*s2 + 4*s2**2)
    alpha_2 = (1-d)*(s1) + d*(s2**2)
    alpha_3 = (1-d)*(s0) + d*(2*s1*s2)
    alpha_4 = d*(2*s0*s2 + s1**2)
    alpha_5 = d*(2*s0*s1)
    alpha_6 = d*(s0**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

#ABCD -> AABC -> (AAAB, AABB) -> AAAA
  predict4 = function(r1, r2, r3, r4, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2-r3-r4 #aaaa
    #if (r0 < 0) {return(0)}
    t0 = r0**k #AAAA
    t1 = (r0+r1)**k #AAAA + AAAB
    t2 = (r0+r2)**k #AAAA + AABB
    t3 = (r0+r1+r2+r3)**k #AAAA + AAAB + AABB + AABC
    s0 = t0 #AAAA
    s1 = t1-t0 #AAAB
    s2 = t2-t0 #AABB
    s3 = t3-t2-t1+t0 #AABC
    s4 = 1-t3 #ABCD
    alpha_1 = (1-d)*(s1+2*s3+4*s4) + d*(2*s0*s1 + 4*s0*s3 + 6*s0*s4 + 2*s1**2 + 2*s1*s2 + 6*s1*s3 + 8*s1*s4 + 4*s2*s3 + 6*s2*s4 + 4*s3**2 + 10*s3*s4 + 6*s4**2)
    alpha_2 = (1-d)*(2*s2+s3) + d*(2*s0*s2 + 2*s1*s2 + 2*s2**2 + 2*s2*s3 + 2*s2*s4 + s4**2)
    alpha_3 = (1-d)*(s1) + d*(2*s2*s4 + 2*s3*s4)
    alpha_4 = (1-d)*(s0) + d*(2*s1*s4 + s2**2 + 2*s2*s3 + s3**2)
    alpha_5 = d*(2*s0*s4 + 2*s1*s2 + 2*s1*s3)
    alpha_6 = d*(2*s0*s2 + 2*s0*s3 + s1**2)
    alpha_7 = d*(2*s0*s1)
    alpha_8 = d*(s0**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
  }

#ACBD -> (ABAC, AABC, ACBB)-> (AAAB, ABAA, AABB) -> AAAA
  predict4 = function(r1, r2, r3, k, d, kmercov, bias, x)
  {
    #Am --(r1)--> Bm
    #Am --(r2)--> Ap
    #Bm --(r3)--> Bp
    if (r1 < r2 | r2 < r3) {return(0)}
    raaaa = (1-r1)*(1-r2)*(1-r3) #aaaa
    raaab = (1-r1)*(1-r2)*(r3)
    rabaa = (1-r1)*(r2)*(1-r3)
    raabb = r1*(1-r2)*(1-r3)
    rabac = (1-r1)*r2*r3
    raabc = r1*(1-r2)*r3
    racbb = r1*r2*(1-r3)
    racbd = r1*r2*r3
    tAAAA = (raaaa)**k
    tAAAB = (raaab + raaaa)**k
    tABAA = (rabaa + raaaa)**k
    tAABB = (raabb + raaaa)**k
    tABAC = (rabac + raaab + rabaa + raaaa)**k
    tAABC = (raabc + raaab + raabb + raaaa)**k
    tACBB = (racbb + rabaa + raabb + raaaa)**k
    tACBD = (racbd + rabac + raabc + racbb + raaab + rabaa + raabb + raaaa)**k
    sAAAA = tAAAA
    sAAAB = tAAAB - tAAAA
    sABAA = tABAA - tAAAA
    sAABB = tAABB - tAAAA
    sABAC = tABAC - tAAAB - tABAA + tAAAA
    sAABC = tAABC - tAAAB - tAABB + tAAAA
    sACBB = tACBB - tABAA - tAABB + tAAAA
    sACBD = tACBD - tABAC - tAABC - tACBB + tAAAB + tABAA + tAABB - tAAAA
    #if (r0 < 0) {return(0)}
    #t0 = r0**k #AAAA
    #t1 = (r0+r1)**k #AAAA + AAAB
    #t2 = (r0+r2)**k #AAAA + AABB
    #t3 = (r0+r1+r2+r3)**k #AAAA + AAAB + AABB + AABC
    #s0 = t0 #AAAA
    #s1 = t1-t0 #AAAB
    #s2 = t2-t0 #AABB
    #s3 = t3-t2-t1+t0 #AABC
    #s4 = 1-t3 #ABCD
    alpha_1 = (1-d)*(sAAAB + sABAA + 2*sABAC + 2*sAABC + 2*sACBB + 4*sACBD) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 2*sAAAB*sAABB + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABB*sAABC + 4*sAABC**2 + 2*sAAAA*sABAA + 4*sAAAB*sABAA + 2*sAABB*sABAA + 6*sAABC*sABAA + 2*sABAA**2 + 4*sAAAA*sABAC + 6*sAAAB*sABAC + 4*sAABB*sABAC + 8*sAABC*sABAC + 6*sABAA*sABAC + 4*sABAC**2 + 2*sAAAA*sACBB + 4*sAAAB*sACBB + 2*sAABB*sACBB + 6*sAABC*sACBB + 4*sABAA*sACBB + 6*sABAC*sACBB + 2*sACBB**2 + 6*sAAAA*sACBD + 8*sAAAB*sACBD + 6*sAABB*sACBD + 10*sAABC*sACBD + 8*sABAA*sACBD + 10*sABAC*sACBD + 8*sACBB*sACBD + 6*sACBD**2)
    alpha_2 = (1-d)*(2*sAABB + sABAC + sAABC + sACBB) + d*(2*sAAAA*sAABB + 2*sAAAB*sAABB + 2*sAABB**2 + 2*sAABB*sAABC + 2*sAABB*sABAA + 2*sAABB*sABAC + 2*sAAAA*sACBB + 2*sAAAB*sACBB + 4*sAABB*sACBB + 2*sAABC*sACBB + 2*sABAA*sACBB + 2*sABAC*sACBB + 3*sACBB**2 + 2*sAABB*sACBD + 4*sACBB*sACBD + sACBD**2)
    alpha_3 = (1-d)*(sAAAB + sABAA) + d*(2*sAABB*sACBB + 2*sAABC*sACBB + 2*sABAC*sACBB + 2*sAABB*sACBD + 2*sAABC*sACBD + 2*sABAC*sACBD)
    alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2 + 2*sAABB*sABAC + 2*sAABC*sABAC + sABAC**2 + 2*sAAAB*sACBB + 2*sABAA*sACBB + 2*sAAAB*sACBD + 2*sABAA*sACBD)
    alpha_5 = d*(2*sAAAB*sAABB + 2*sAAAB*sAABC + 2*sAABB*sABAA + 2*sAABC*sABAA + 2*sAAAB*sABAC + 2*sABAA*sABAC + 2*sAAAA*sACBB + 2*sAAAA*sACBD)
    alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABB + 2*sAAAA*sAABC + 2*sAAAB*sABAA + sABAA**2 + 2*sAAAA*sABAC)
    alpha_7 = d*(2*sAAAA*sAAAB + 2*sAAAA*sABAA)
    alpha_8 = d*(sAAAA**2)
#    alpha_1 = (1-d)*(s1+2*s3+4*s4) + d*(2*s0*s1 + 4*s0*s3 + 6*s0*s4 + 2*s1**2 + 2*s1*s2 + 6*s1*s3 + 8*s1*s4 + 4*s2*s3 + 6*s2*s4 + 4*s3**2 + 10*s3*s4 + 6*s4**2)
#    alpha_2 = (1-d)*(2*s2+s3) + d*(2*s0*s2 + 2*s1*s2 + 2*s2**2 + 2*s2*s3 + 2*s2*s4 + s4**2)
#    alpha_3 = (1-d)*(s1) + d*(2*s2*s4 + 2*s3*s4)
#    alpha_4 = (1-d)*(s0) + d*(2*s1*s4 + s2**2 + 2*s2*s3 + s3**2)
#    alpha_5 = d*(2*s0*s4 + 2*s1*s2 + 2*s1*s3)
#    alpha_6 = d*(2*s0*s2 + 2*s0*s3 + s1**2)
#    alpha_7 = d*(2*s0*s1)
#    alpha_8 = d*(s0**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
  }

##ABCD -> AABC -> AAAB -> AAAA
#  predict4 = function(r1, r3, r4, k, d, kmercov, bias, x)
#  {
#    r2 = 0
#    r0 = 1-r1-r2-r3-r4 #aaaa
#    t0 = r0**k #AAAA
#    t1 = (r0+r1)**k #AAAA + AAAB
#    t2 = (r0+r2)**k #AAAA + AABB
#    t3 = (r0+r1+r2+r3)**k #AAAA + AAAB + AABB + AABC
#    s0 = t0 #AAAA
#    s1 = t1-t0 #AAAB
#    s2 = t2-t0 #AABB
#    s3 = t3-t2-t1+t0 #AABC
#    s4 = 1-t3 #ABCD
#    alpha_1 = (1-d)*(s1+2*s3+4*s4) + d*(2*s0*s1 + 4*s0*s3 + 6*s0*s4 + 2*s1**2 + 2*s1*s2 + 6*s1*s3 + 8*s1*s4 + 4*s2*s3 + 6*s2*s4 + 4*s3**2 + 10*s3*s4 + 6*s4**2)
#    alpha_2 = (1-d)*(2*s2+s3) + d*(2*s0*s2 + 2*s1*s2 + 2*s2**2 + 2*s2*s3 + 2*s2*s4 + s4**2)
#    alpha_3 = (1-d)*(s1) + d*(2*s2*s4 + 2*s3*s4)
#    alpha_4 = (1-d)*(s0) + d*(2*s1*s4 + s2**2 + 2*s2*s3 + s3**2)
#    alpha_5 = d*(2*s0*s4 + 2*s1*s2 + 2*s1*s3)
#    alpha_6 = d*(2*s0*s2 + 2*s0*s3 + s1**2)
#    alpha_7 = d*(2*s0*s1)
#    alpha_8 = d*(s0**2)
#    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
#    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
#    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
#    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
#    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
#    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
#    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
#    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
#  }

  predict5 = function(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2-r3-r4-r5-r6
    if (r0 < 0) {return(0)}
    t0 = (r0)**k
    t1 = (r0+r1)**k
    t2 = (r0+r2)**k
    t3 = (r0+r1+r2+r3)**k
    t4 = (r0+r1+r2+r4)**k
    t5 = (r0+r1+r2+r3+r4+r5)**k
    s0 = t0
    s1 = t1-t0
    s2 = t2-t0
    s3 = t3-t1-t2+t0
    s4 = t4-t1-t2+t0
    s5 = t5-t3-t4+t1+t2-t0
    s6 = 1-t5
    alpha_1  = (1-d)*(s1 + 2*s3 + s4 + 3*s5 + 5*s6) + d*(2*s0*s1 + 2*s1**2 + 2*s1*s2 + 4*s0*s3 + 6*s1*s3 + 4*s2*s3 + 4*s3**2 + 2*s0*s4 + 4*s1*s4 + 2*s2*s4 + 6*s3*s4 + 2*s4**2 + 6*s0*s5 + 8*s1*s5 + 6*s2*s5 + 10*s3*s5 + 8*s4*s5 + 6*s5**2 + 8*s0*s6 + 10*s1*s6 + 8*s2*s6 + 12*s3*s6 + 10*s4*s6 + 14*s5*s6 + 8*s6**2)
    alpha_2  = (1-d)*(s2 + 2*s4 + s5) + d*(2*s0*s2 + 2*s1*s2 + 2*s2**2 + 2*s2*s3 + 2*s0*s4 + 2*s1*s4 + 4*s2*s4 + 2*s3*s4 + 2*s4**2 + 2*s2*s5 + 2*s4*s5 + 2*s2*s6 + 2*s4*s6 + s6**2)
    alpha_3  = (1-d)*(s2 + s3) + d*(2*s4*s6 + 2*s5*s6)
    alpha_4  = (1-d)*(s1) + d*(s4**2 + 2*s4*s5 + s5**2 + 2*s2*s6 + 2*s3*s6)
    alpha_5  = (1-d)*(s0) + d*(2*s2*s4 + 2*s3*s4 + 2*s2*s5 + 2*s3*s5 + 2*s1*s6)
    alpha_6  = d*(s2**2 + 2*s2*s3 + s3**2 + 2*s1*s4 + 2*s1*s5 + 2*s0*s6)
    alpha_7  = d*(2*s1*s2 + 2*s1*s3 + 2*s0*s4 + 2*s0*s5)
    alpha_8  = d*(s1**2 + 2*s0*s2 + 2*s0*s3)
    alpha_9  = d*(2*s0*s1)
    alpha_10 = d*(s0**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)
  }

  predict6 = function(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x)
  {
    r0  = 1-r1-r2-r3-r4-r5-r6-r7-r8-r9-r10
    if (r0 < 0) {return(0)}
    t0  = (r0)**k
    t1  = (r0+r1)**k
    t2  = (r0+r2)**k
    t3  = (r0+r3)**k
    t4  = (r0+r1+r2+r4)**k
    t5  = (r0+r1+r2+r3+r5)**k
    t6  = (r0+r2+r6)**k
    t7  = (r0+r1+r2+r3+r4+r5+r7)**k
    t8  = (r0+r1+r2+r3+r4+r5+r6+r8)**k
    t9  = (r0+r1+r2+r3+r4+r5+r6+r7+r8+r9)**k
    s0  = t0
    s1  = t1-t0
    s2  = t2-t0
    s3  = t3-t0
    s4  = t4-t1-t2+t0
    s5  = t5-t1-t2-t3+2*t0
    s6  = t6-t2
    s7  = t7-t4-t5+t1+t2-t0
    s8  = t8-t4-t5-t6+t1+2*t2-t0
    s9 = t9-t7-t8+t4+t5-t2-t1+t0
    s10 = 1-t9
    alpha_1  = (1-d)*(4*s9 + 6*s10 + s1 + 2*s4 + s5 + 3*s7 + 2*s8) + d*(8*s0*s9 + 8*s9**2 + 10*s0*s10 + 18*s9*s10 + 10*s10**2 + 2*s0*s1 + 10*s9*s1 + 12*s10*s1 + 2*s1**2 + 8*s9*s2 + 10*s10*s2 + 2*s1*s2 + 8*s9*s3 + 10*s10*s3 + 2*s1*s3 + 4*s0*s4 + 12*s9*s4 + 14*s10*s4 + 6*s1*s4 + 4*s2*s4 + 4*s3*s4 + 4*s4**2 + 2*s0*s5 + 10*s9*s5 + 12*s10*s5 + 4*s1*s5 + 2*s2*s5 + 2*s3*s5 + 6*s4*s5 + 2*s5**2 + 8*s9*s6 + 10*s10*s6 + 2*s1*s6 + 4*s4*s6 + 2*s5*s6 + 6*s0*s7 + 14*s9*s7 + 16*s10*s7 + 8*s1*s7 + 6*s2*s7 + 6*s3*s7 + 10*s4*s7 + 8*s5*s7 + 6*s6*s7 + 6*s7**2 + 4*s0*s8 + 12*s9*s8 + 14*s10*s8 + 6*s1*s8 + 4*s2*s8 + 4*s3*s8 + 8*s4*s8 + 6*s5*s8 + 4*s6*s8 + 10*s7*s8 + 4*s8**2)
    alpha_2  = (1-d)*(s9 + s2 + s5 + 3*s6 + 2*s8) + d*(s10**2 + 2*s0*s2 + 2*s9*s2 + 2*s10*s2 + 2*s1*s2 + 2*s2**2 + 2*s2*s3 + 2*s2*s4 + 2*s0*s5 + 2*s9*s5 + 2*s10*s5 + 2*s1*s5 + 4*s2*s5 + 2*s3*s5 + 2*s4*s5 + 2*s5**2 + 4*s0*s6 + 4*s9*s6 + 4*s10*s6 + 4*s1*s6 + 6*s2*s6 + 4*s3*s6 + 4*s4*s6 + 6*s5*s6 + 4*s6**2 + 2*s2*s7 + 2*s5*s7 + 4*s6*s7 + 2*s0*s8 + 2*s9*s8 + 2*s10*s8 + 2*s1*s8 + 4*s2*s8 + 2*s3*s8 + 2*s4*s8 + 4*s5*s8 + 6*s6*s8 + 2*s7*s8 + 2*s8**2)
    alpha_3  = (1-d)*(2*s3 + s5 + s7) + d*(2*s9*s10 + 2*s0*s3 + 2*s9*s3 + 2*s10*s3 + 2*s1*s3 + 2*s2*s3 + 2*s3**2 + 2*s3*s4 + 2*s3*s5 + 2*s10*s6 + 2*s3*s6 + 2*s3*s7 + 2*s10*s8 + 2*s3*s8)
    alpha_4  = (1-d)*(s2 + s4) + d*(s9**2 + 2*s10*s3 + 2*s10*s5 + 2*s9*s6 + s6**2 + 2*s10*s7 + 2*s9*s8 + 2*s6*s8 + s8**2)
    alpha_5  = (1-d)*(s1) + d*(2*s10*s2 + 2*s9*s3 + 2*s10*s4 + 2*s9*s5 + 2*s3*s6 + 2*s5*s6 + 2*s9*s7 + 2*s6*s7 + 2*s3*s8 + 2*s5*s8 + 2*s7*s8)
    alpha_6  = (1-d)*(s0) + d*(2*s10*s1 + 2*s9*s2 + s3**2 + 2*s9*s4 + 2*s3*s5 + s5**2 + 2*s2*s6 + 2*s4*s6 + 2*s3*s7 + 2*s5*s7 + s7**2 + 2*s2*s8 + 2*s4*s8)
    alpha_7  = d*(2*s0*s10 + 2*s9*s1 + 2*s2*s3 + 2*s3*s4 + 2*s2*s5 + 2*s4*s5 + 2*s1*s6 + 2*s2*s7 + 2*s4*s7 + 2*s1*s8)
    alpha_8  = d*(2*s0*s9 + s2**2 + 2*s1*s3 + 2*s2*s4 + s4**2 + 2*s1*s5 + 2*s0*s6 + 2*s1*s7 + 2*s0*s8)
    alpha_9  = d*(2*s1*s2 + 2*s0*s3 + 2*s1*s4 + 2*s0*s5 + 2*s0*s7)
    alpha_10 = d*(s1**2 + 2*s0*s2 + 2*s0*s4)
    alpha_11 = d*(2*s0*s1)
    alpha_12 = d*(s0**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

#try p=1 model
  if (p == 1) {
    try(model <- 
      nlsLM(y ~ length*predict1(k, d, kmercov, bias, x),
      start = list(d=0.001, kmercov=estKmercov, bias = 0.5, length = estLength),
      lower = c(0, 0, 0, 0),
      upper = c(1, Inf, Inf, Inf),
      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) { cat("retrying nls_peak with port algorithm\n") }

      try(model <- 
      nlsLM(y ~ length*predict1(k, d, kmercov, bias, x),
      start = list(d=0, kmercov=estKmercov, bias = 0.5, length=estLength),
      lower = c(0, 0, 0, 0),
      upper = c(1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 1
      model$het   = c(0,0)
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet  = 0
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
    }
  }

#try p=2 model
  if (p == 2) {
    try(model <- 
      nlsLM(y ~ length*predict2(r1, k, d, kmercov, bias, x),
      start = list(d=0.001, r1=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/2),
      lower = c(0, 0, 0, 0, 0),
      upper = c(1, 1, Inf, Inf, Inf),
      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) { cat("retrying nls_peak with port algorithm\n") }

      try(model <- 
      nlsLM(y ~ length*predict2(r1, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, kmercov=estKmercov, bias = 0.5, length=estLength/2),
      lower = c(0, 0, 0, 0, 0),
      upper = c(1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 2
      #print('HIIIIIII')
      #print(model_sum$coefficients['r1',])
      model$het1  = min_max1(model_sum$coefficients['r1',])
      #print(model$het1)
      model$het   = model$het1
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet  = model$ahet1
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
    }
  }

#try p=3 model
  if (p==3) {
    try(model <- 
    nlsLM(y ~ length*predict3(r1, r2, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/3),
    lower = c(0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) { cat("retrying nls_peak with port algorithm\n") }

      try(model <- 
      nlsLM(y ~ length*predict3(r1, r2, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, kmercov=estKmercov, bias = 0.5, length=estLength/3),
      lower = c(0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 3
      model$het1  = min_max1(model_sum$coefficients['r1',])
      model$het2  = min_max1(model_sum$coefficients['r2',])
      model$het   = model$het1 + model$het2
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet2 = model_sum$coefficients['r2',][[1]]
      model$ahet  = model$ahet1 + model$ahet2
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
      #if (is.null(model) || model_3$m$deviance() *10 < model$m$deviance() || (model_3$m$deviance() < model$m$deviance() && !(model_3$het > 2*model$het) && !(model_3$d > 2*model$d))){
      #  model <- model_3
      #  model$p = 3
      #}
    }
  }

#try p=4 model
  if (p==4) {
    try(model <- 
    nlsLM(y ~ length*predict4(r1, r2, r3, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2 = 0.001, r3=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/4),
    lower = c(0, 0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) { cat("retrying nls_peak with port algorithm\n") }

      try(model <- 
      nlsLM(y ~ length*predict4(r1, r2, r3, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, r3=0, kmercov=estKmercov, bias = 0.5, length=estLength/4),
      lower = c(0, 0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	  }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 4
      model$het1  = min_max1(model_sum$coefficients['r1',])
      model$het2  = min_max1(model_sum$coefficients['r2',])
      #model$het2  = c(0,0)
      model$het3  = min_max1(model_sum$coefficients['r3',])
      #model$het4  = min_max1(model_sum$coefficients['r4',])
      model$het4  = c(0,0)
      model$het   = model$het1 + model$het2 + model$het3 + model$het4
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet2 = model_sum$coefficients['r2',][[1]]
      #model$ahet2 = 0
      model$ahet3 = model_sum$coefficients['r3',][[1]]
      #model$ahet4 = model_sum$coefficients['r4',][[1]]
      model$ahet4 = 0
      model$ahet  = model$ahet1 + model$ahet2 + model$ahet3 + model$ahet4
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]

      #if (is.null(model) || model_4$m$deviance() *10 < model$m$deviance() || (model_4$m$deviance() < model$m$deviance() && !(model_4$het > 2*model$het) && !(model_4$d > 2*model$d))){
      #  model <- model_4
      #  model$p = 4
      #}
    }
  }

#try p=5 model
  if (p==5) {
    try(model <- 
    nlsLM(y ~ length*predict5(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2=0.001, r3=0.001, r4=0.001, r5=0.001, r6=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/5),
    lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) { cat("retrying nls_peak with port algorithm\n") }

      try(model <- 
      nlsLM(y ~ length*predict5(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, r3=0, r4=0, r5=0, r6=0, kmercov=estKmercov, bias = 0.5, length=estLength/5),
      lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	  }

    if (!is.null(model))
    {
      model_sum   = summary(model)
      model$p     = 5
      model$het1  = min_max1(model_sum$coefficients['r1',])
      model$het2  = min_max1(model_sum$coefficients['r2',])
      model$het3  = min_max1(model_sum$coefficients['r3',])
      model$het4  = min_max1(model_sum$coefficients['r4',])
      model$het5  = min_max1(model_sum$coefficients['r5',])
      model$het6  = min_max1(model_sum$coefficients['r6',])
      model$het   = model$het1 + model$het2 + model$het3 + model$het4 + model$het5 + model$het6
      model$homo  = 1-model$het
      model$dups  = min_max(model_sum$coefficients['bias',])
      model$kcov  = min_max(model_sum$coefficients['kmercov',])
      model$mlen  = min_max(model_sum$coefficients['length',])
      model$md    = min_max1(model_sum$coefficients['d',])
      model$ahet1 = model_sum$coefficients['r1',][[1]]
      model$ahet2 = model_sum$coefficients['r2',][[1]]
      model$ahet3 = model_sum$coefficients['r3',][[1]]
      model$ahet4 = model_sum$coefficients['r4',][[1]]
      model$ahet5 = model_sum$coefficients['r5',][[1]]
      model$ahet6 = model_sum$coefficients['r6',][[1]]
      model$ahet  = model$ahet1 + model$ahet2 + model$ahet3 + model$ahet4 + model$ahet5 + model$ahet6
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
      #if (is.null(model) || model_5$m$deviance() *10 < model$m$deviance() || (model_5$m$deviance() < model$m$deviance() && !(model_5$het > 2*model$het) && !(model_5$d > 2*model$d))){
      #  model <- model_5
      #  model$p = 5
      #}
    }
  }

#try p=6 model
  if (p==6) {
    try(model <- 
    nlsLM(y ~ length*predict6(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x),
    start = list(d=0.001, r1=0.001, r2=0.001, r3=0.001, r4=0.001, r5=0.001, r6=0.001, r7=0.001, r8=0.001, r9=0.001, r10=0.001, kmercov=estKmercov, bias = 0.5, length = estLength/6),
    lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    upper = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
    control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model) == "try-error" || is.null(model)){

      if (VERBOSE) { cat("retrying nls_peak with port algorithm\n") }

      try(model <- 
      nlsLM(y ~ length*predict6(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x),
      start = list(d=0, r1=0, r2=0, r3=0, r4=0, r5=0, r6=0, r7=0, r8=0, r9=0, r10=0, kmercov=estKmercov, bias = 0.5, length=estLength/6),
      lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      upper = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, Inf, Inf, Inf),
      algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	  }

    if (!is.null(model))
    {
      model_sum    = summary(model)
      model$p      = 6
      model$het1   = min_max1(model_sum$coefficients['r1',])
      model$het2   = min_max1(model_sum$coefficients['r2',])
      model$het3   = min_max1(model_sum$coefficients['r3',])
      model$het4   = min_max1(model_sum$coefficients['r4',])
      model$het5   = min_max1(model_sum$coefficients['r5',])
      model$het6   = min_max1(model_sum$coefficients['r6',])
      model$het7   = min_max1(model_sum$coefficients['r7',])
      model$het8   = min_max1(model_sum$coefficients['r8',])
      model$het9   = min_max1(model_sum$coefficients['r9',])
      model$het10  = min_max1(model_sum$coefficients['r10',])
      model$het    = model$het1 + model$het2 + model$het3 + model$het4 + model$het5 + model$het6 + model$het7 + model$het8 + model$het9 + model$het10
      model$homo   = 1-model$het
      model$dups   = min_max(model_sum$coefficients['bias',])
      model$kcov   = min_max(model_sum$coefficients['kmercov',])
      model$mlen   = min_max(model_sum$coefficients['length',])
      model$md     = min_max1(model_sum$coefficients['d',])
      model$ahet1  = model_sum$coefficients['r1',][[1]]
      model$ahet2  = model_sum$coefficients['r2',][[1]]
      model$ahet3  = model_sum$coefficients['r3',][[1]]
      model$ahet4  = model_sum$coefficients['r4',][[1]]
      model$ahet5  = model_sum$coefficients['r5',][[1]]
      model$ahet6  = model_sum$coefficients['r6',][[1]]
      model$ahet7  = model_sum$coefficients['r7',][[1]]
      model$ahet8  = model_sum$coefficients['r8',][[1]]
      model$ahet9  = model_sum$coefficients['r9',][[1]]
      model$ahet10 = model_sum$coefficients['r10',][[1]]
      model$ahet   = model$ahet1 + model$ahet2 + model$ahet3 + model$ahet4 + model$ahet5 + model$ahet6 + model$ahet7 + model$ahet8 + model$ahet9 + model$ahet10
      model$ahomo = 1-model$ahet
      model$adups = model_sum$coefficients['bias',][[1]]
      model$akcov = model_sum$coefficients['kmercov',][[1]]
      model$amlen = model_sum$coefficients['length',][[1]]
      model$amd   = model_sum$coefficients['d',][[1]]
      #if (is.null(model) || model_6$m$deviance() *10 < model$m$deviance() || (model_6$m$deviance() < model$m$deviance() && !(model_6$het > 2*model$het) && !(model_6$d > 2*model$d))){
      #  model <- model_6
      #  model$p = 6
      #}
    }
  }

  print(model)

  return(model)
}

## score model by number and percent of residual errors after excluding sequencing errors
#########################################################################################

score_model<-function(kmer_hist_orig, nls, round, foldername){
  x = kmer_hist_orig[[1]]
  y = kmer_hist_orig[[2]]

  pred=predict(nls, newdata=data.frame(x))
  model_sum=summary(nls)
  p=nls$p
  kcovfloor = max(1, floor(min_max(model_sum$coefficients['kmercov',])[[1]]))

  ## Compute error rate, by counting kmers unexplained by model through first peak
  ## truncate errors as soon as it goes to zero, dont allow it to go back up
  error_xcutoff = kcovfloor
  error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1) #this may not work if gaps in simulated data
  if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

  error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]

  first_zero = -1

  for (i in 1:error_xcutoff_ind)
  {
    if (first_zero == -1)
    {
      if (error_kmers[i] < 1.0)
      {
        first_zero = i
        if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
      }
    }
    else
    {
      error_kmers[i] = 0
    }
  }

  if (first_zero == -1)
  {
    first_zero = error_xcutoff_ind
  }

  ## The fit is residual sum of square error, excluding sequencing errors
  model_fit_all    = c(sum(as.numeric(y[first_zero:length(y)]     - pred[first_zero:length(y)])     ** 2), first_zero, x[length(y)])
  model_fit_full   = c(sum(as.numeric(y[first_zero:(min(length(y),(2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y), (2*p+1)*kcovfloor))]) ** 2), first_zero, (min(length(y), (2*p+1)*kcovfloor)))
  model_fit_unique = c(sum(as.numeric(y[first_zero:((p+1)*kcovfloor)] - pred[first_zero:((p+1)*kcovfloor)]) ** 2), first_zero, ((p+1)*kcovfloor))

  ## The score is the percentage of unexplained kmers, excluding sequencing errors
  model_fit_allscore    = c(1-sum(abs(as.numeric(y[first_zero:length(y)]     - pred[first_zero:length(y)])))     / sum(as.numeric(y[first_zero:length(y)])),     first_zero, x[length(y)])
  model_fit_fullscore   = c(1-sum(abs(as.numeric(y[first_zero:(min(length(y), (2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y), (2*p+1)*kcovfloor))]))) / sum(as.numeric(y[first_zero:(min(length(y), (2*p+1)*kcovfloor))])), first_zero, (min(length(y), (2*p+1)*kcovfloor)))
  model_fit_uniquescore = c(1-sum(abs(as.numeric(y[first_zero:((p+1)*kcovfloor)] - pred[first_zero:((p+1)*kcovfloor)]))) / sum(as.numeric(y[first_zero:((p+1)*kcovfloor)])), first_zero, ((p+1)*kcovfloor))

  fit = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                   full = model_fit_full,     fullscore = model_fit_fullscore,
                   unique = model_fit_unique, uniquescore = model_fit_uniquescore)

  return (fit)
}

## Pick between the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
###############################################################################

eval_model<-function(kmer_hist_orig, nls0, nls1, p, round, foldername){
  nls0score = -1
  nls1score = -1

  ## Evaluate the score the nls0
  if (!is.null(nls0))
  {
    nls0score = score_model(kmer_hist_orig, nls0, round+0.1, foldername)

    if(VERBOSE){ cat(paste("nls0score$all:\t", nls0score$all[[1]], "\n"))}

    if (VERBOSE)
    {
      mdir = paste(foldername, "/round", round, ".1", sep="")
      dir.create(mdir, showWarnings=FALSE)
      report_results(kmer_prof_orig,kmer_prof_orig, k, p, (list(nls0, nls0score)) , mdir)
    }
  }
  else
  {
    if (VERBOSE) { cat("nls0score failed to converge\n") }
  }


  ## Evaluate the score of nls1
  if (!is.null(nls1))
  {
    nls1score = score_model(kmer_hist_orig, nls1, round+0.2, foldername)

    if(VERBOSE){ cat(paste("nls1score$all:\t", nls1score$all[[1]], "\n"))}

    if (VERBOSE)
    {
      mdir = paste(foldername, "/round", round, ".2", sep="")
      dir.create(mdir, showWarnings=FALSE)
      report_results(kmer_prof_orig, kmer_prof_orig, k, p, (list(nls1, nls1score)) , mdir)
    }
  }
  else
  {
    if (VERBOSE) { cat("nls1score failed to converge\n") }
  }

  ## Return the better of the scores
  if (!is.null(nls0))
  {
    if (!is.null(nls1))
    {
      pdiff = abs(nls0score$all[[1]] - nls1score$all[[1]]) / max(nls0score$all[[1]], nls1score$all[[1]])

      if (pdiff < SCORE_CLOSE)
      {
        het0 = nls0$ahet
        het1 = nls1$ahet

        if (het1 * SCORE_HET_FOLD_DIFFERENCE < het0)
        {
          if (VERBOSE) { cat(paste("returning nls0, similar score, higher het\n")) }
          return (list(nls0, nls0score)) #originally nls0
        }
        else if (het0 * SCORE_HET_FOLD_DIFFERENCE < het1)
        {
          if (VERBOSE) { cat(paste("returning nls1, similar score, higher het\n")) }
          return (list(nls1, nls1score)) #originally nls1
        }
      }

      if (nls0score$all[[1]] < nls1score$all[[1]])
      {
        if (VERBOSE) { cat(paste("returning nls0, better score\n")) }
        return (list(nls0, nls0score))
      }
      else
      {
        if (VERBOSE) { cat(paste("returning nls1, better score\n")) }
        return (list(nls1, nls1score))
      }
    }
    else
    {
      if (VERBOSE) { cat(paste("returning nls0, nls1 fail\n")) }
      return (list(nls0, nls0score))
    }
  }

  if (VERBOSE) { cat(paste("returning nls1 by default\n")) }
  return (list(nls1, nls1score))
}

## Wrapper function to try fitting 2p peak model with 2 forms
###############################################################################

estimate_Genome_peak2<-function(kmer_hist_orig, x, y, k, p, estKmercov, round, foldername){
  ## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
  #numofReads   = sum(as.numeric(x)*as.numeric(y))/(readlength-k+1)
  numofKmers = sum(as.numeric(x)*as.numeric(y))
  if (estKmercov==-1){
    estKmercov1  = x[which(y==max(y))][1]
  } else {
    estKmercov1 = estKmercov
  }
  #estCoverage1 = estKmercov1*readlength/(readlength-k+1)
  #estLength1   = numofReads*readlength/estCoverage1
  estLength1 = numofKmers/estKmercov1

  if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov1, "\n")) }

  nls0 = nls_peak(x, y, k, p, estKmercov1, estLength1, MAX_ITERATIONS)

  if (VERBOSE) { print(summary(nls0)) }

  if (estKmercov==-1){
    ploidies = c(1,2,3,4,5,6)
    for (i in ploidies){
      if (i<=p){

        ## next we i-th the max kmercoverage (in the case of ploidy >= i)
        estKmercov2  = estKmercov1 / i
        #estCoverage2 = estKmercov2*readlength/(readlength-k+1)
        #estLength2   = numofReads*readlength/estCoverage2
        estLength2 = numofKmers/estKmercov2

        if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov2, "\n")) }

        nls1 = nls_peak(x, y, k, p, estKmercov2, estLength2, MAX_ITERATIONS)

        if (VERBOSE) { print(summary(nls1)) }

        if (i<p){
          nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername)[[1]]
        }

      }

    }
    return(eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername))
  }

  return(eval_model(kmer_hist_orig, nls0, nls0, p, round, foldername))
}

## Format numbers
###############################################################################
bp_format<-function(num) {paste(formatC(round(num),format="f",big.mark=",", digits=0), "bp",sep=" ")}

percentage_format<-function(num) {paste(signif(num,6)*100,"%",sep="")}

X_format<-function(num) {paste(signif(num,4),"X",sep="")}

## Report results and make plots
###############################################################################

report_results<-function(kmer_hist,kmer_hist_orig, k, p, container, foldername){

  predict1_unique = function(k, d, kmercov, bias, x)
  {
    r0 = 1
    t0 = r0**k
    s0 = t0
    alpha_1_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)
  }

#AB -> AA
  predict2_unique = function(r1, k, d, kmercov, bias, x)
  {
    r0 = 1-r1
    t0 = r0**k
    s0 = t0
    s1 = 1-t0
    alpha_1_unique = (1-d)*(2*s1)
    alpha_2_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)
  }

#ABC -> AAB -> AAA
  predict3_unique = function(r1, r2, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2
    t0 = r0**k
    t1 = (r0+r1)**k
    s0 = t0
    s1 = t1-t0
    s2 = 1-t1
    alpha_1_unique = (1-d)*(s1+3*s2)
    alpha_2_unique = (1-d)*(s1)
    alpha_3_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)
  }

#ABCD -> AABC -> (AAAB, AABB) -> AAAA
  predict4_unique = function(r1, r2, r3, r4, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2-r3-r4
    t0 = r0**k
    t1 = (r0+r1)**k
    t2 = (r0+r2)**k
    t3 = (r0+r1+r2+r3)**k
    s0 = t0
    s1 = t1-t0
    s2 = t2-t0
    s3 = t3-t2-t1+t0
    s4 = 1-t3
    alpha_1_unique = (1-d)*(s1+2*s3+4*s4)
    alpha_2_unique = (1-d)*(2*s2+s3)
    alpha_3_unique = (1-d)*(s1)
    alpha_4_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
  }

##ABCD -> AABC -> AAAB -> AAAA
#  predict4_unique = function(r1, r3, r4, k, d, kmercov, bias, x)
#  {
#    r2 = 0
#    r0 = 1-r1-r2-r3-r4
#    t0 = r0**k
#    t1 = (r0+r1)**k
#    t2 = (r0+r2)**k
#    t3 = (r0+r1+r2+r3)**k
#    s0 = t0
#    s1 = t1-t0
#    s2 = t2-t0
#    s3 = t3-t2-t1+t0
#    s4 = 1-t3
#    alpha_1_unique = (1-d)*(s1+2*s3+4*s4)
#    alpha_2_unique = (1-d)*(2*s2+s3)
#    alpha_3_unique = (1-d)*(s1)
#    alpha_4_unique = (1-d)*(s0)
#    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
#    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
#    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
#    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
#  }

  predict5_unique = function(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2-r3-r4-r5-r6
    t0 = (r0)**k
    t1 = (r0+r1)**k
    t2 = (r0+r2)**k
    t3 = (r0+r1+r2+r3)**k
    t4 = (r0+r1+r2+r4)**k
    t5 = (r0+r1+r2+r3+r4+r5)**k
    s0 = t0
    s1 = t1-t0
    s2 = t2-t0
    s3 = t3-t1-t2+t0
    s4 = t4-t1-t2+t0
    s5 = t5-t3-t4+t1+t2-t0
    s6 = 1-t5
    alpha_1_unique  = (1-d)*(s1 + 2*s3 + s4 + 3*s5 + 5*s6) + d*(2*s0*s1 + 2*s1**2 + 2*s1*s2 + 4*s0*s3 + 6*s1*s3 + 4*s2*s3 + 4*s3**2 + 2*s0*s4 + 4*s1*s4 + 2*s2*s4 + 6*s3*s4 + 2*s4**2 + 6*s0*s5 + 8*s1*s5 + 6*s2*s5 + 10*s3*s5 + 8*s4*s5 + 6*s5**2 + 8*s0*s6 + 10*s1*s6 + 8*s2*s6 + 12*s3*s6 + 10*s4*s6 + 14*s5*s6 + 8*s6**2)
    alpha_2_unique  = (1-d)*(s2 + 2*s4 + s5) + d*(2*s0*s2 + 2*s1*s2 + 2*s2**2 + 2*s2*s3 + 2*s0*s4 + 2*s1*s4 + 4*s2*s4 + 2*s3*s4 + 2*s4**2 + 2*s2*s5 + 2*s4*s5 + 2*s2*s6 + 2*s4*s6 + s6**2)
    alpha_3_unique  = (1-d)*(s2 + s3) + d*(2*s4*s6 + 2*s5*s6)
    alpha_4_unique  = (1-d)*(s1) + d*(s4**2 + 2*s4*s5 + s5**2 + 2*s2*s6 + 2*s3*s6)
    alpha_5_unique  = (1-d)*(s0) + d*(2*s2*s4 + 2*s3*s4 + 2*s2*s5 + 2*s3*s5 + 2*s1*s6)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  predict6_unique = function(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x)
  {
    r0  = 1-r1-r2-r3-r4-r5-r6-r7-r8-r9-r10
    t0  = (r0)**k
    t1  = (r0+r1)**k
    t2  = (r0+r2)**k
    t3  = (r0+r3)**k
    t4  = (r0+r1+r2+r4)**k
    t5  = (r0+r1+r2+r3+r5)**k
    t6  = (r0+r2+r6)**k
    t7  = (r0+r1+r2+r3+r4+r5+r7)**k
    t8  = (r0+r1+r2+r3+r4+r5+r6+r8)**k
    t9  = (r0+r1+r2+r3+r4+r5+r6+r7+r8+r9)**k
    s0  = t0
    s1  = t1-t0
    s2  = t2-t0
    s3  = t3-t0
    s4  = t4-t1-t2+t0
    s5  = t5-t1-t2-t3+2*t0
    s6  = t6-t2
    s7  = t7-t4-t5+t1+t2-t0
    s8  = t8-t4-t5-t6+t1+2*t2-t0
    s9 = t9-t7-t8+t4+t5-t2-t1+t0
    s10 = 1-t9
    alpha_1_unique  = (1-d)*(4*s9 + 6*s10 + s1 + 2*s4 + s5 + 3*s7 + 2*s8)
    alpha_2_unique  = (1-d)*(s9 + s2 + s5 + 3*s6 + 2*s8)
    alpha_3_unique  = (1-d)*(2*s3 + s5 + s7)
    alpha_4_unique  = (1-d)*(s2 + s4)
    alpha_5_unique  = (1-d)*(s1)
    alpha_6_unique  = (1-d)*(s0)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  x=kmer_hist_orig[[1]]
  y=kmer_hist_orig[[2]]
  model = container[[1]]

  #automatically zoom into the relevant regions of the plot, ignore first 15 positions
  xmax=length(x)
  start=which(y == min(y[1:TYPICAL_ERROR]))
  zoomx=x[start:(xmax-1)]
  zoomy=y[start:(xmax-1)]

  ## allow for a little space above max value past the noise
  y_limit = max(zoomy[start:length(zoomy)])*1.1

  x_limit = which(y == max(y[start:length(zoomx)])) * 3

  if (min(zoomy) > zoomy[1]){
    x_limit=max(which(zoomy<zoomy[1])[2],600)
  }

  if (!is.null(model))
  {
     model_sum=summary(model)
     kcov = min_max(model_sum$coefficients['kmercov',])[1]
     x_limit = max(kcov*(2*p+1.1), x_limit)
  }

  ## Uncomment this to enforce a specific number
  # x_limit=150

  ## Features to report
  het=c(-1,-1)
  homo=c(-1,-1)
  het1=c(-1,-1)
  het2=c(-1,-1)
  het3=c(-1,-1)
  het4=c(-1,-1)
  het5=c(-1,-1)
  het6=c(-1,-1)
  het7=c(-1,-1)
  het8=c(-1,-1)
  het9=c(-1,-1)
  het10=c(-1,-1)
  total_len=c(-1,-1)
  repeat_len=c(-1,-1)
  unique_len=c(-1,-1)
  dups=c(-1,-1)
  error_rate=c(-1,-1)
  model_status="fail"

  model_fit_unique      = c(0,0,0)
  model_fit_full        = c(0,0,0)
  model_fit_all         = c(0,0,0)
  model_fit_allscore    = c(0,0,0)
  model_fit_fullscore   = c(0,0,0)
  model_fit_uniquescore = c(0,0,0)

  plot_size=2000
  font_size=1.2
  resolution=300

  ## Plot the distribution, and hopefully with the model fit
  png(paste(foldername, "/", args$name_prefix, "_plot.png", sep=""),
  width=plot_size, height=plot_size, res=resolution)
  par(mar = c(5.1,4.1,5.1,2.1))
  plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
  xlab="Coverage", ylab="Frequency", ylim=c(0,y_limit), xlim=c(0,x_limit),
  cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
  rect(0, 0, x_limit*1.1 , y_limit*1.1, col=COLOR_BGCOLOR)
  points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
#  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
#    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
#  }
  box(col="black")

  ## Make a second plot in log space over entire range
  png(paste(foldername, "/", args$name_prefix, "_plot.log.png", sep=""),
  width=plot_size, height=plot_size, res=resolution)
  par(mar = c(5.1,4.1,5.1,2.1))
  plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
  xlab="Coverage", ylab="Frequency", log="xy",
  cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  rect(1e-10, 1e-10, max(kmer_hist_orig[[1]])*10 , max(kmer_hist_orig[[2]])*10, col=COLOR_BGCOLOR)
  points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
  }
  box(col="black")


  if(!is.null(model))
  {
    x=kmer_hist[[1]]
    y=kmer_hist[[2]]

    ## The model converged!
    pred=predict(model, newdata=data.frame(x))

    ## Compute the genome characteristics
    model_sum=summary(model)
    print(model_sum)

    ## save the model to a file
    capture.output(model_sum, file=paste(foldername,"/", args$name_prefix, "_model.txt", sep=""))

    ## Identify key values
    het   = model$het
    ahet  = model$ahet
    homo  = model$homo
    ahomo = model$ahomo
    if (p >= 2)
    {
      het1  = model$het1
      ahet1 = model$ahet1
    }
    if (p >= 3)
    {
      het2  = model$het2
      ahet2 = model$ahet2
    }
    if (p >= 4)
    {
      het3  = model$het3
      het4  = model$het4
      ahet3 = model$ahet3
      ahet4 = model$ahet4
    }
    if (p >= 5)
    {
      het5  = model$het5
      het6  = model$het6
      ahet5 = model$ahet5
      ahet6 = model$ahet6
    }
    if (p >= 6)
    {
      het7   = model$het7
      het8   = model$het8
      het9   = model$het9
      het10  = model$het10
      ahet7  = model$ahet7
      ahet8  = model$ahet8
      ahet9  = model$ahet9
      ahet10 = model$ahet10
    }

    dups = model$dups
    kcov = model$kcov
    mlen = model$mlen
    md   = model$md

    adups = model$adups
    akcov = model$akcov
    amlen = model$amlen
    amd   = model$amd

    ## Compute error rate, by counting kmers unexplained by model through first peak
    ## truncate errors as soon as it goes to zero, dont allow it to go back up
    error_xcutoff = max(1, floor(kcov[1]))
    error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
    if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

    error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]

    first_zero = -1

    for (i in 1:error_xcutoff_ind)
    {
      if (first_zero == -1)
      {
        if (error_kmers[i] < 1.0)
        {
          first_zero = i
          if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
        }
      }
      else
      {
        error_kmers[i] = 0
      }
    }

    if (first_zero == -1)
    {
      first_zero = error_xcutoff_ind
    }

    ## Rather than "0", set to be some very small number so log-log plot looks okay
    error_kmers = pmax(error_kmers, 1e-10)

    total_error_kmers = sum(as.numeric(error_kmers) * as.numeric(x[1:error_xcutoff_ind]))
    total_kmers = sum(as.numeric(x)*as.numeric(y))

    error_rate = 1-(1-(total_error_kmers/total_kmers))**(1/k)
    error_rate = c(error_rate, error_rate)

    total_len = (total_kmers-total_error_kmers)/(p*kcov)

    ## find kmers that fit the p peak model (no repeats)
    if (p==1)
    {
      unique_hist = amlen*predict1_unique(k, amd, akcov, adups, x)
    }
    if (p==2)
    {
      unique_hist = amlen*predict2_unique(ahet1, k, amd, akcov, adups, x)
    }
    if (p==3)
    {
      unique_hist = amlen*predict3_unique(ahet1, ahet2, k, amd, akcov, adups, x)
    }
    if (p==4)
    {
      unique_hist = amlen*predict4_unique(ahet1, ahet2, ahet3, ahet4, k, amd, akcov, adups, x)
    }
    if (p==5)
    {
      unique_hist = amlen*predict5_unique(ahet1, ahet2, ahet3, ahet4, ahet5, ahet6, k, amd, akcov, adups, x)
    }
    if (p==6)
    {
      unique_hist = amlen*predict6_unique(ahet1, ahet2, ahet3, ahet4, ahet5, ahet6, ahet7, ahet8, ahet9, ahet10, k, amd, akcov, adups, x)
    }

    unique_kmers = sum(as.numeric(x)*as.numeric(unique_hist))
    repeat_kmers = total_kmers - unique_kmers - total_error_kmers

    repeat_len=repeat_kmers/(p*kcov)
    unique_len=unique_kmers/(p*kcov)

    score = container[[2]]

    model_fit_allscore    = score$allscore
    model_fit_fullscore   = score$fullscore
    model_fit_uniquescore = score$uniquescore

    model_fit_all    = score$all
    model_fit_full   = score$full
    model_fit_unique = score$unique

    residual = y - pred

    if (p==1){hetline = paste0("A:",     format(100*ahomo, digits=3), "%")}
    if (p==2){hetline = paste0("AA:",     format(100*ahomo, digits=3), "% ",
                               "AB:",     format(100*ahet1, digits=3), "%")}
    if (p==3){hetline = paste0("AAA:",    format(100*ahomo, digits=3), "% ",
                               "AAB:",    format(100*ahet1, digits=3), "% ",
                               "ABC:",    format(100*ahet2, digits=3), "%")}
    if (p==4){hetline = paste0("AAAA:",   format(100*ahomo, digits=3), "% ",
                               "AAAB:",   format(100*ahet1, digits=3),"% ",
                               "AABB:",   format(100*ahet2, digits=3), "% ",
                               "AABC:",   format(100*ahet3, digits=3), "% ",
                               "ABCD:",   format(100*ahet4, digits=3), "%")}
    if (p==5){hetline = paste0("AAAAA:",  format(100*ahomo, digits=3), "% ",
                               "AAAAB:",  format(100*ahet1, digits=3), "% ",
                               "AAABB:",  format(100*ahet2, digits=3), "% ",'\n',
                               "AAABC:",  format(100*ahet3, digits=3), "% ",
                               "AABBC:",  format(100*ahet4, digits=3), "% ",
                               "AABCD:",  format(100*ahet5, digits=3), "% ",
                               "ABCDE:",  format(100*ahet6, digits=3), "%")}
    if (p==6){hetline = paste0("AAAAAA:", format(100*ahomo, digits=3), "% ",
                               "AAAAAB:", format(100*ahet1, digits=3), "% ",
                               "AAAABB:", format(100*ahet2, digits=3), "% ",
                               "AAABBB:", format(100*ahet3, digits=3), "% ",
                               "AAAABC:", format(100*ahet4, digits=3), "% ",'\n',
                               "AAABBC:", format(100*ahet5, digits=3), "% ",
                               "AABBCC:", format(100*ahet6, digits=3), "% ",
                               "AAABCD:", format(100*ahet7, digits=3), "% ",
                               "AABBCD:", format(100*ahet8, digits=3), "% ",
                               "AABCDE:", format(100*ahet9, digits=3), "% ",
                               "ABCDEF:", format(100*ahet10, digits=3), "%")}
    print(hetline)
    ## Finish Log plot
    title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
    "bp",
    " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
    "% ", "\n",
    hetline, "\n",
    " kcov:", format(akcov, digits=3),
    " err:",   format(100*error_rate[1], digits=3),
    "% ",
    " dup:",  format(adups, digits=3),
    " ",
    " k:",   format(k, digits=3),
    " p:",   format(p, digits=3),
    sep=""),
    cex.main=.85)

    ## Mark the modes of the peaks
    abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

    ## Draw just the unique portion of the model
    lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
    lines(x, pred, col=COLOR_2pPEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

    if (VERBOSE) { lines(x, residual, col=COLOR_RESIDUAL, lwd=3) }

    ## Add legend
    if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
    {
      legend(exp(.65 * log(max(x))), 1.0 * max(y),
      legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
      lty=c("solid", "solid", "solid", "solid", "dashed"),
      lwd=c(3,3,3,3,3),
      col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
      bg="white")
    }
    else
    {
      legend("topright",
      ##legend(exp(.65 * log(max(x))), 1.0 * max(y),
      legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
      lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
      lwd=c(3,3,3,3,2,3),
      col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
      bg="white")
    }

    dev.set(dev.next())

    ## Finish Linear Plot
    title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
    "bp",
    " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
    "% ", "\n",
    hetline, "\n",
    " kcov:", format(akcov, digits=3),
    " err:",   format(100*error_rate[1], digits=3),
    "% ",
    " dup:",  format(adups, digits=3),
    " ",
    " k:",   format(k, digits=3),
    " p:",   format(p, digits=3),
    sep=""),
    cex.main=.85)

    ## Mark the modes of the peaks
    abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

    ## Draw just the unique portion of the model
    lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
    lines(x, pred, col=COLOR_2pPEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

    if (VERBOSE) { lines(x, residual, col=COLOR_RESIDUAL, lwd=3) }

    ## Add legend
    legend(.65 * x_limit, 1.0 * y_limit,
    legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
    lty=c("solid", "solid", "solid", "solid", "dashed"),
    lwd=c(3,3,3,3,2),
    col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
    bg="white")

    model_status="done"

    cat(paste("Model converged het:", format(ahet, digits=3),
    " kcov:", format(akcov, digits=3),
    " err:", format(error_rate[1], digits=3),
    " model fit:", format(adups, digits=3),
    " len:", round(total_len[1]), "\n", sep=""))
  }
  else
  {
    title("\nFailed to converge")
    dev.set(dev.next())
    title("\nFailed to converge")
    cat("Failed to converge.", file=paste(foldername,"/", args$name_prefix, "_model.txt", sep=""))
    cat("Failed to converge\n")
  }

  dev.off()
  dev.off()

  ## Write key values to summary file
  summaryFile <- paste(foldername,"/", args$name_prefix, "_summary.txt",sep="")

  format_column_1 = "%-30s"
  format_column_2 = "%-18s"
  format_column_3 = "%-18s"

  cat(paste("GenomeScope version 2.0", sep=""), file=summaryFile, sep="\n")
  cat(paste("input file = ", args$input, sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste("k = ", k,sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste("p = ", p,sep=""), file=summaryFile, sep="\n", append=TRUE)
  if (args$lambda!=-1) {
    cat(paste("initial kmercov estimate = ", args$lambda, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  if (args$max_kmercov!=-1) {
    cat(paste("max_kmercov = ", args$max_kmercov, sep=""), file=summaryFile, sep="\n", append=TRUE)
  }
  cat(paste("\n",sprintf(format_column_1,"property"),               sprintf(format_column_2,"min"),                              sprintf(format_column_3,"max"), sep=""),                                     file=summaryFile, sep="\n", append=TRUE)
  if (p==1)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (A)"),            sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==2)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (AA)"),           sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (AB)"),         sprintf(format_column_2,percentage_format(het1[1])),         sprintf(format_column_3,percentage_format(het1[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==3)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (AAA)"),          sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not AAA)"),    sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAB"),                       sprintf(format_column_2,percentage_format(het1[1])),         sprintf(format_column_3,percentage_format(het1[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"ABC"),                       sprintf(format_column_2,percentage_format(het2[1])),         sprintf(format_column_3,percentage_format(het2[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==4)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (AAAA)"),         sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not AAAA)"),   sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAAB"),                      sprintf(format_column_2,percentage_format(het1[1])),         sprintf(format_column_3,percentage_format(het1[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABB"),                      sprintf(format_column_2,percentage_format(het2[1])),         sprintf(format_column_3,percentage_format(het2[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABC"),                      sprintf(format_column_2,percentage_format(het3[1])),         sprintf(format_column_3,percentage_format(het3[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"ABCD"),                      sprintf(format_column_2,percentage_format(het4[1])),         sprintf(format_column_3,percentage_format(het4[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==5)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (AAAAA)"),        sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_2,"Heterozygous (not AAAAA)"),  sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAAAB"),                     sprintf(format_column_2,percentage_format(het1[1])),         sprintf(format_column_3,percentage_format(het1[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAABB"),                     sprintf(format_column_2,percentage_format(het2[1])),         sprintf(format_column_3,percentage_format(het2[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAABC"),                     sprintf(format_column_2,percentage_format(het3[1])),         sprintf(format_column_3,percentage_format(het3[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABBC"),                     sprintf(format_column_2,percentage_format(het4[1])),         sprintf(format_column_3,percentage_format(het4[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABCD"),                     sprintf(format_column_2,percentage_format(het5[1])),         sprintf(format_column_3,percentage_format(het5[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"ABCDE"),                     sprintf(format_column_2,percentage_format(het6[1])),         sprintf(format_column_3,percentage_format(het6[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
  }
  if (p==6)
  {
    cat(paste(sprintf(format_column_1,"Homozygous (AAAAAA)"),       sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"Heterozygous (not AAAAAA)"), sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])), sep=""),                 file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAAAAB"),                    sprintf(format_column_2,percentage_format(het1[1])),         sprintf(format_column_3,percentage_format(het1[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAAABB"),                    sprintf(format_column_2,percentage_format(het2[1])),         sprintf(format_column_3,percentage_format(het2[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAABBB"),                    sprintf(format_column_2,percentage_format(het3[1])),         sprintf(format_column_3,percentage_format(het3[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAAABC"),                    sprintf(format_column_2,percentage_format(het4[1])),         sprintf(format_column_3,percentage_format(het4[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAABBC"),                    sprintf(format_column_2,percentage_format(het5[1])),         sprintf(format_column_3,percentage_format(het5[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABBCC"),                    sprintf(format_column_2,percentage_format(het6[1])),         sprintf(format_column_3,percentage_format(het6[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AAABCD"),                    sprintf(format_column_2,percentage_format(het7[1])),         sprintf(format_column_3,percentage_format(het7[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABBCD"),                    sprintf(format_column_2,percentage_format(het8[1])),         sprintf(format_column_3,percentage_format(het8[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"AABCDE"),                    sprintf(format_column_2,percentage_format(het9[1])),         sprintf(format_column_3,percentage_format(het9[2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
    cat(paste(sprintf(format_column_1,"ABCDEF"),                    sprintf(format_column_2,percentage_format(het10[1])),        sprintf(format_column_3,percentage_format(het10[2])), sep=""),               file=summaryFile, sep="\n", append=TRUE)
  }
  cat(paste(sprintf(format_column_1,"Genome Haploid Length"), sprintf(format_column_2,bp_format(total_len[2])),                  sprintf(format_column_3,bp_format(total_len[1])), sep=""),                   file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Repeat Length"),  sprintf(format_column_2,bp_format(repeat_len[2])),                 sprintf(format_column_3,bp_format(repeat_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Unique Length"),  sprintf(format_column_2,bp_format(unique_len[2])),                 sprintf(format_column_3,bp_format(unique_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Model Fit "),            sprintf(format_column_2,percentage_format(model_fit_allscore[1])), sprintf(format_column_3,percentage_format(model_fit_fullscore[1])), sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Read Error Rate"),       sprintf(format_column_2,percentage_format(error_rate[1])),         sprintf(format_column_3,percentage_format(error_rate[2])), sep=""),          file=summaryFile, sep="\n", append=TRUE)
  if (VERBOSE)
  {
    cat(paste("\nPercent Kmers Modeled (All Kmers) = ",  percentage_format(model_fit_allscore[1]),    " [", model_fit_allscore[2],    ", ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Percent Kmers Modeled (Full Model) = ",   percentage_format(model_fit_fullscore[1]),   " [", model_fit_fullscore[2],   ", ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Percent Kmers Modeled (Unique Kmers) = ", percentage_format(model_fit_uniquescore[1]), " [", model_fit_uniquescore[2], ", ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

    cat(paste("\nModel RSSE (All Kmers) = ",  model_fit_all[1],    " [", model_fit_all[2],    ", ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model RSSE (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   ", ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model RSSE (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], ", ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)	
  }
  ## Finalize the progress
  progressFilename=paste(foldername,"/", args$name_prefix, "_progress.txt",sep="")
  cat(model_status, file=progressFilename, sep="\n", append=TRUE)
}

## Main program starts here
###############################################################################
parser <- ArgumentParser()
parser$add_argument("-v", "--version", action="store_true", default=FALSE, help="print the version and exit")
parser$add_argument("-i", "--input", help = "input histogram file")
parser$add_argument("-o", "--output", help = "output directory name")
parser$add_argument("-k", "--kmer_size", type = "integer", default = 21, help = "kmer size used to calculate kmer spectra [default 21]")
parser$add_argument("-p", "--ploidy", type = "integer", default = 2, help = "ploidy for model to use [default 2]")
parser$add_argument("-l", "--lambda", "--kcov", "--kmercov", type = "integer", default=-1, help = "optional initial kmercov estimate for model to use")
parser$add_argument("-m", "--max_kmercov", type = "integer", default=-1, help = "optional maximum kmer coverage threshold (kmers with coverage greater than max_kmercov are ignored by the model)")
parser$add_argument("-n", "--name_prefix", default = "OUTPUT", help = "name prefix for output files [default OUTPUT]")
parser$add_argument("--verbose", action="store_true", default=FALSE, help = "optional flag to print messages during execution")

args <- parser$parse_args()
version_message <- "GenomeScope 2.0\n"

if (args$version) {
  cat(version_message)
  quit()
}

#args<-commandArgs(TRUE)

#if (length(args) < 5) {
#  cat("USAGE: genomescope.R histogram_file k-mer_length read_length ploidy output_dir [lambda] [kmer_max] [verbose]\n")
#} else {
if (is.null(args$input) | is.null(args$output)) {
  cat("USAGE: genomescope.R -i input_histogram_file -k kmer_length -p ploidy -o output_dir\n")
  cat("OPTIONAL PARAMETERS: -l lambda -m max_kmercov -n 'name_prefix' --verbose\n")
  cat("HELP: genomescope.R --help\n")
} else {

  ## Load the arguments from the user
  histfile   <- args$input #args[[1]]
  k          <- args$kmer_size #as.numeric(args[[2]])
  #readlength <- as.numeric(args[[3]])
  p          <- args$ploidy #as.numeric(args[[4]])
  foldername <- args$output #args[[5]]
  estKmercov <- args$lambda
  maxCovGenomeLen <- args$max_kmercov

  #estKmercov = -1
  #if ((length(args) >= 6))
  #{
  #  estKmercov = as.numeric(args[[6]])
  #}
  #maxCovGenomeLen = -1

  #if ((length(args) >= 7))
  #{
  #  maxCovGenomeLen = as.numeric(args[[7]])
  #}

  #if ((length(args) == 8) && (as.numeric(args[[8]] == 1))) { VERBOSE = 1 }

  ## values for testing
  #histfile <- "~/build/genomescope/simulation/simulation_results/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"
  #k <- 21
  #readlength <- 100
  #foldername <- "~/build/genomescope/simulation/simulation_analysis/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"

  #if (k > readlength) { stop("K cannot be greater than readlength") }

  #cat(paste("GenomeScope analyzing ", histfile, " k=", k, " p=", p, " readlen=", readlength, " outdir=", foldername, "\n", sep=""))
  cat(paste("GenomeScope analyzing ", histfile, " k=", k, " p=", p, " outdir=", foldername, "\n", sep=""))

  dir.create(foldername, showWarnings=FALSE)

  kmer_prof <- read.csv(file=histfile,sep="", header=FALSE)

  minkmerx = 1;
  if (kmer_prof[1,1] == 0)
  {
    if (VERBOSE) { cat("Histogram starts with zero, reseting minkmerx\n");  }
    minkmerx = 2;
  }

  kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
  kmer_prof_orig <- kmer_prof

  ## Initialize the status
  progressFilename <- paste(foldername,"/", args$name_prefix, "_progress.txt",sep="")
  cat("starting", file=progressFilename, sep="\n")

  ## try to find the local minimum between errors and the first (heterozygous) peak
  start <- tail(which(kmer_prof[1:TYPICAL_ERROR,2]==min(kmer_prof[1:TYPICAL_ERROR,2])),n=1)

  maxCovIndex = -1

  ## Figure out which kmers to exclude, if any
  if(maxCovGenomeLen == -1)
  {
    maxCovIndex <- length(kmer_prof[,1])
  }
  else
  {
    ## Figure out the index we should use for this coverage length
    x <- kmer_prof[,1]
    maxCovIndex <- length(x[x<=maxCovGenomeLen])
  }

  if (VERBOSE) { cat(paste("using maxCovGenomeLen:", maxCovGenomeLen, " with index:", maxCovIndex, "trying 2p peak model... \n")) }

  # terminate after NUM_ROUND iterations, store best result so far in container
  round <- 0
  best_container <- list(NULL,0)

  while(round < NUM_ROUNDS)
  {
    cat(paste("round", round, "trimming to", start, "trying 2p peak model... "), file=progressFilename, sep="", append=TRUE)
    if (VERBOSE) { cat(paste("round", round, "trimming to", start, "trying 2p peak model... \n")) }

    ## Reset the input trimming off low frequency error kmers
    kmer_prof=kmer_prof_orig[1:maxCovIndex,]
    x <- kmer_prof[start:maxCovIndex,1]
    y <- kmer_prof[start:maxCovIndex,2]

    model_peaks <- estimate_Genome_peak2(kmer_prof, x, y, k, p, estKmercov, round, foldername)

    if (!is.null(model_peaks[[1]]))
    {
      cat(paste("converged. score: ", model_peaks[[2]]$all[[1]]), file=progressFilename, sep="\n", append=TRUE)

      if (VERBOSE)
      {
        mdir = paste(foldername, "/round", round, sep="")
        dir.create(mdir, showWarnings=FALSE)
        report_results(kmer_prof,kmer_prof_orig, k, p, model_peaks, mdir)
      }
    }
    else
    {
     cat(paste("unconverged"), file=progressFilename, sep="\n", append=TRUE)
    }

    #check if this result is better than previous
    if (!is.null(model_peaks[[1]]))
    {
      if (is.null(best_container[[1]]))
      {
        if (VERBOSE) { cat(paste("no previous best, updating best")) }
        best_container = model_peaks
      }
      else
      {
        pdiff = abs(model_peaks[[2]]$all[[1]] - best_container[[2]]$all[[1]]) / max(model_peaks[[2]]$all[[1]], best_container[[2]]$all[[1]])

        if (pdiff < SCORE_CLOSE)
        {
          hetm = model_peaks[[1]]$ahet
          hetb = best_container[[1]]$ahet

          if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm)
          {
            if (VERBOSE) { cat(paste("model has significantly higher heterozygosity but similar score, overruling")) }
            best_container = model_peaks #originally best_container = model_peaks
          }
          else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb)
          {
            if (VERBOSE) { cat(paste("previous best has significantly higher heterozygosity and similar score, keeping")) }
            #best_container = model_peaks #originally blank
          }
          else if (model_peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
          {
            if (VERBOSE) { cat(paste("score is marginally better but het rate is not extremely different, updating")) }
            best_container = model_peaks
          }
        }
        else if (model_peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
        {
          if (VERBOSE) { cat(paste("score is significantly better, updating")) }
          best_container = model_peaks
        }
      }
    }

    ## Ignore a larger number of kmers as errors
    start <- start + START_SHIFT
    round <- round + 1
  }
  ## Report the results, note using the original full profile
  report_results(kmer_prof,kmer_prof_orig, k, p, best_container, foldername)
  if (!is.null(best_container[[1]])){
    print('model score')
    print(best_container[[2]]$all[[1]])
    print(best_container[[1]]$m$deviance())
  }
}
