#' Produce model estimated (p=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict1_1 = function(k, d, kmercov, bias, x)
{
  r0 = 1
  t0 = r0**k
  s0 = t0
  alpha_1 = (1-d)*(s0)
  alpha_2 = d*(s0)
  alpha_1 * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size=kmercov*2 / bias, mu = kmercov*2)
}

#' Produce model estimated (p=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict1_unique = function(k, d, kmercov, bias, x)
{
  r0 = 1
  t0 = r0**k
  s0 = t0
  alpha_1_unique = (1-d)*(s0)
  alpha_1_unique * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)
}

#AB -> AA
#' Produce model estimated (p=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1 A numeric corresponding to the nucleotide heterozygosity ab.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict2_1 = function(r1, k, d, kmercov, bias, x)
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

#AB -> AA
#' Produce model estimated (p=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1 A numeric corresponding to the nucleotide heterozygosity ab.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict2_1_unique = function(r1, k, d, kmercov, bias, x)
{
  r0 = 1-r1
  if (r0 < 0) {return(0)}
  t0 = r0**k
  s0 = t0
  s1 = 1-t0
  alpha_1_unique = (1-d)*(2*s1)
  alpha_2_unique = (1-d)*(s0)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)
}

#ABC -> AAB -> AAA
#' Produce model estimated (p=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1,r2 Numerics corresponding to the nucleotide heterozygosities aab and abc respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict3_1 = function(r1, r2, k, d, kmercov, bias, x)
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

#ABC -> AAB -> AAA
#' Produce model estimated (p=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1,r2 Numerics corresponding to the nucleotide heterozygosities aab and abc respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict3_1_unique = function(r1, r2, k, d, kmercov, bias, x)
{
  r0 = 1-r1-r2
  if (r0 < 0) {return(0)}
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



#AAAA -> AAAB -> AABC -> ABCD
#' Produce model estimated (p=4, topology=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaab,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabc, and abcd respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict4_1 = function(raaab, raabc, rabcd, k, d, kmercov, bias, x)
{
  raaaa = 1-raaab-raabc-rabcd
  if (raaaa < 0) {return(0)}
  tAAAA = raaaa**k
  tAAAB = (raaaa+raaab)**k
  tAABC = (raaaa+raaab+raabc)**k
  sAAAA = tAAAA
  sAAAB = tAAAB-tAAAA
  sAABC = tAABC-tAAAB
  sABCD = 1-tAABC
  alpha_1 = (1-d)*(sAAAB + 2*sAABC + 4*sABCD) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABC**2 + 6*sAAAA*sABCD + 8*sAAAB*sABCD + 10*sAABC*sABCD + 6*sABCD**2)
  alpha_2 = (1-d)*(sAABC) + d*(sABCD**2)
  alpha_3 = (1-d)*(sAAAB) + d*(2*sAABC*sABCD)
  alpha_4 = (1-d)*(sAAAA) + d*(sAABC**2 + 2*sAAAB*sABCD)
  alpha_5 = d*(2*sAAAB*sAABC + 2*sAAAA*sABCD)
  alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABC)
  alpha_7 = d*(2*sAAAA*sAAAB)
  alpha_8 = d*(sAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
}

#AAAA -> AAAB -> AABC -> ABCD
#' Produce model estimated (p=4, topology=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaab,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabc, and abcd respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict4_1_unique = function(raaab, raabc, rabcd, k, d, kmercov, bias, x)
{
  raaaa = 1-raaab-raabc-rabcd
  if (raaaa < 0) {return(0)}
  tAAAA = raaaa**k
  tAAAB = (raaaa+raaab)**k
  tAABC = (raaaa+raaab+raabc)**k
  sAAAA = tAAAA
  sAAAB = tAAAB-tAAAA
  sAABC = tAABC-tAAAB
  sABCD = 1-tAABC
  alpha_1_unique = (1-d)*(sAAAB + 2*sAABC + 4*sABCD)
  alpha_2_unique = (1-d)*(sAABC)
  alpha_3_unique = (1-d)*(sAAAB)
  alpha_4_unique = (1-d)*(sAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
}

#AAAA -> AABB -> AABC -> ABCD
#' Produce model estimated (p=4, topology=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aabb, aabc, and abcd respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict4_2 = function(raabb, raabc, rabcd, k, d, kmercov, bias, x)
{
  raaaa = 1-raabb-raabc-rabcd
  if (raaaa < 0) {return(0)}
  tAAAA = raaaa**k
  tAABB = (raaaa+raabb)**k
  tAABC = (raaaa+raabb+raabc)**k
  sAAAA = tAAAA
  sAABB = tAABB-tAAAA
  sAABC = tAABC-tAABB
  sABCD = 1-tAABC
  alpha_1 = (1-d)*(2*sAABC + 4*sABCD) + d*(4*sAAAA*sAABC + 4*sAABB*sAABC + 4*sAABC**2 + 6*sAAAA*sABCD + 6*sAABB*sABCD + 10*sAABC*sABCD + 6*sABCD**2)
  alpha_2 = (1-d)*(2*sAABB + sAABC) + d*(2*sAAAA*sAABB + 2*sAABB**2 + 2*sAABB*sAABC + 2*sAABB*sABCD + sABCD**2)
  alpha_3 = (1-d)*(0) + d*(2*sAABB*sABCD + 2*sAABC*sABCD)
  alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2)
  alpha_5 = d*(2*sAAAA*sABCD)
  alpha_6 = d*(2*sAAAA*sAABB + 2*sAAAA*sAABC)
  alpha_7 = d*(0)
  alpha_8 = d*(sAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
}

#AAAA -> AABB -> AABC -> ABCD
#' Produce model estimated (p=4, topology=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aabb, aabc, and abcd respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict4_2_unique = function(raabb, raabc, rabcd, k, d, kmercov, bias, x)
{
  raaaa = 1-raabb-raabc-rabcd
  if (raaaa < 0) {return(0)}
  tAAAA = raaaa**k
  tAABB = (raaaa+raabb)**k
  tAABC = (raaaa+raabb+raabc)**k
  sAAAA = tAAAA
  sAABB = tAABB-tAAAA
  sAABC = tAABC-tAABB
  sABCD = 1-tAABC
  alpha_1_unique = (1-d)*(2*sAABC + 4*sABCD)
  alpha_2_unique = (1-d)*(2*sAABB + sAABC)
  alpha_3_unique = (1-d)*(0)
  alpha_4_unique = (1-d)*(sAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
}

#' Produce model estimated (p=5) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1,r2,r3,r4,r5,r6 Numerics corresponding to the nucleotide heterozygosities aaaab, aaabb, aaabc, aabbc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
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

#' Produce model estimated (p=5, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1,r2,r3,r4,r5,r6 Numerics corresponding to the nucleotide heterozygosities aaaab, aaabb, aaabc, aabbc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_unique = function(r1, r2, r3, r4, r5, r6, k, d, kmercov, bias, x)
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

#' Produce model estimated (p=6) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1,r2,r3,r4,r5,r6,r7,r8,r9,r10 Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabb, aaabbb, aaaabc, aaabbc, aabbcc, aaabcd, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
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

#' Produce model estimated (p=6, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param r1,r2,r3,r4,r5,r6,r7,r8,r9,r10 Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabb, aaabbb, aaaabc, aaabbc, aabbcc, aaabcd, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_unique = function(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, k, d, kmercov, bias, x)
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
