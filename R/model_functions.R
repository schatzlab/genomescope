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
predict3_0 = function(r1, r2, k, d, kmercov, bias, x)
{
  predict3_1(r1, r2, k, d, kmercov, bias, x)
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

predict3_0_unique = function(r1, r2, k, d, kmercov, bias, x)
{
  predict3_1_unique(r1, r2, k, d, kmercov, bias, x)
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

#AAAA -> (AAAB, AABB) -> AABC -> ABCD
#' Produce model estimated (p=4, full model) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaab,raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabb, aabc, and abcd respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict4_0 = function(raaab, raabb, raabc, rabcd, k, d, kmercov, bias, x)
{
  raaaa = 1-raaab-raabb-raabc-rabcd
  if (raaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAA = raaaa
    sAAAB = raaab
    sAABB = raabb
    sAABC = raabc
    sABCD = rabcd
  } else {
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raaab+raabb+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAAAB-tAABB+tAAAA
    sABCD = 1-tAABC
  }
  alpha_1 = (1-d)*(sAAAB + 2*sAABC + 4*sABCD) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 2*sAAAB*sAABB + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABB*sAABC + 4*sAABC**2 + 6*sAAAA*sABCD + 8*sAAAB*sABCD + 6*sAABB*sABCD + 10*sAABC*sABCD + 6*sABCD**2)
  alpha_2 = (1-d)*(2*sAABB + sAABC) + d*(2*sAAAA*sAABB + 2*sAAAB*sAABB + 2*sAABB**2 + 2*sAABB*sAABC + 2*sAABB*sABCD + sABCD**2)
  alpha_3 = (1-d)*(sAAAB) + d*(2*sAABB*sABCD + 2*sAABC*sABCD)
  alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2 + 2*sAAAB*sABCD)
  alpha_5 = d*(2*sAAAB*sAABB + 2*sAAAB*sAABC + 2*sAAAA*sABCD)
  alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABB + 2*sAAAA*sAABC)
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

#AAAA -> (AAAB, AABB) -> AABC -> ABCD
#' Produce model estimated (p=4, full model, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaab,raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabb, aabc, and abcd respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict4_0_unique = function(raaab, raabb, raabc, rabcd, k, d, kmercov, bias, x)
{
  raaaa = 1-raaab-raabb-raabc-rabcd
  if (raaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAA = raaaa
    sAAAB = raaab
    sAABB = raabb
    sAABC = raabc
    sABCD = rabcd
  } else {
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raaab+raabb+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAAAB-tAABB+tAAAA
    sABCD = 1-tAABC
  }
  alpha_1_unique = (1-d)*(sAAAB + 2*sAABC + 4*sABCD)
  alpha_2_unique = (1-d)*(2*sAABB + sAABC)
  alpha_3_unique = (1-d)*(sAAAB)
  alpha_4_unique = (1-d)*(sAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
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
  if (KMER_RATES) {
    sAAAA = raaaa
    sAAAB = raaab
    sAABC = raabc
    sABCD = rabcd
  } else {
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABC = (raaaa+raaab+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABC = tAABC-tAAAB
    sABCD = 1-tAABC
  }
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
  if (KMER_RATES) {
    sAAAA = raaaa
    sAAAB = raaab
    sAABC = raabc
    sABCD = rabcd
  } else {
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABC = (raaaa+raaab+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABC = tAABC-tAAAB
    sABCD = 1-tAABC
  }
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
  if (KMER_RATES) {
    sAAAA = raaaa
    sAABB = raabb
    sAABC = raabc
    sABCD = rabcd
  } else {
    tAAAA = raaaa**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raabb+raabc)**k
    sAAAA = tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAABB
    sABCD = 1-tAABC
  }
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
  if (KMER_RATES) {
    sAAAA = raaaa
    sAABB = raabb
    sAABC = raabc
    sABCD = rabcd
  } else {
    tAAAA = raaaa**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raabb+raabc)**k
    sAAAA = tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAABB
    sABCD = 1-tAABC
  }
  alpha_1_unique = (1-d)*(2*sAABC + 4*sABCD)
  alpha_2_unique = (1-d)*(2*sAABB + sAABC)
  alpha_3_unique = (1-d)*(0)
  alpha_4_unique = (1-d)*(sAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
}

#' Produce model estimated (p=5, full model) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaab,raaabb,raaabc,raabbc,raabcc,raabcd,rabcdd,rabcde Numerics corresponding to the nucleotide heterozygosities.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_0 = function(raaaab, raaabb, raaabc, raabbc, raabcc, raabcd, rabcdd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaaab-raaabb-raaabc-raabbc-raabcc-raabcd-rabcdd-rabcde
  if (raaaaa < 0) {return(0)}
  tAAAAA = (raaaaa)**k
  tAAAAB = (raaaaa+raaaab)**k
  tAAABB = (raaaaa+raaabb)**k
  tAAABC = (raaaaa+raaaab+raaabb+raaabc)**k
  tAABBC = (raaaaa+raaaab+raabbc)**k
  tAABCC = (raaaaa+raaabb+raabcc)**k
  tAABCD = (raaaaa+raaaab+raaabb+raaabc+raabbc+raabcc+raabcd)**k
  tABCDD = (raaaaa+raaabb+raabcc+rabcdd)**k
  sAAAAA = tAAAAA
  sAAAAB = tAAAAB-tAAAAA
  sAAABB = tAAABB-tAAAAA
  sAAABC = tAAABC-tAAAAB-tAAABB+tAAAAA
  sAABBC = tAABBC-tAAAAB
  sAABCC = tAABCC-tAAABB
  sAABCD = tAABCD-tAAABC-tAABBC-tAABCC+tAAAAB+tAAABB
  sABCDD = tABCDD-tAABCC
  sABCDE = 1-tAABCD-tABCDD+tAABCC
  alpha_1  = (1-d)*(sAAAAB + 2*sAAABC + sAABBC + sAABCC + 3*sAABCD + 3*sABCDD + 5*sABCDE) + d*(2*sAAAAA*sAAAAB + 2*sAAAAB**2 + 2*sAAAAB*sAAABB + 4*sAAAAA*sAAABC + 6*sAAAAB*sAAABC + 4*sAAABB*sAAABC + 4*sAAABC**2 + 2*sAAAAA*sAABBC + 4*sAAAAB*sAABBC + 2*sAAABB*sAABBC + 6*sAAABC*sAABBC + 2*sAABBC**2 + 2*sAAAAA*sAABCC + 4*sAAAAB*sAABCC + 2*sAAABB*sAABCC + 6*sAAABC*sAABCC + 4*sAABBC*sAABCC + 2*sAABCC**2 + 6*sAAAAA*sAABCD + 8*sAAAAB*sAABCD + 6*sAAABB*sAABCD + 10*sAAABC*sAABCD + 8*sAABBC*sAABCD + 8*sAABCC*sAABCD + 6*sAABCD**2 + 4*sAAAAA*sABCDD + 6*sAAAAB*sABCDD + 4*sAAABB*sABCDD + 8*sAAABC*sABCDD + 6*sAABBC*sABCDD + 6*sAABCC*sABCDD + 10*sAABCD*sABCDD + 4*sABCDD**2 + 8*sAAAAA*sABCDE + 10*sAAAAB*sABCDE + 8*sAAABB*sABCDE + 12*sAAABC*sABCDE + 10*sAABBC*sABCDE + 10*sAABCC*sABCDE + 14*sAABCD*sABCDE + 12*sABCDD*sABCDE + 8*sABCDE**2)
  alpha_2  = (1-d)*(sAAABB + 2*sAABBC + 2*sAABCC + sAABCD + sABCDD) + d*(2*sAAAAA*sAAABB + 2*sAAAAB*sAAABB + 2*sAAABB**2 + 2*sAAABB*sAAABC + 2*sAAAAA*sAABBC + 2*sAAAAB*sAABBC + 4*sAAABB*sAABBC + 2*sAAABC*sAABBC + 2*sAABBC**2 + 2*sAAAAA*sAABCC + 2*sAAAAB*sAABCC + 4*sAAABB*sAABCC + 2*sAAABC*sAABCC + 4*sAABBC*sAABCC + 2*sAABCC**2 + 2*sAAABB*sAABCD + 2*sAABBC*sAABCD + 2*sAABCC*sAABCD + 2*sAAAAA*sABCDD + 2*sAAAAB*sABCDD + 4*sAAABB*sABCDD + 2*sAAABC*sABCDD + 4*sAABBC*sABCDD + 4*sAABCC*sABCDD + 2*sAABCD*sABCDD + 3*sABCDD**2 + 2*sAAABB*sABCDE + 2*sAABBC*sABCDE + 2*sAABCC*sABCDE + 4*sABCDD*sABCDE + sABCDE**2)
  alpha_3  = (1-d)*(sAAABB + sAAABC) + d*(2*sAABBC*sABCDD + 2*sAABCC*sABCDD + 2*sAABCD*sABCDD + 2*sAABBC*sABCDE + 2*sAABCC*sABCDE + 2*sAABCD*sABCDE)
  alpha_4  = (1-d)*(sAAAAB) + d*(sAABBC**2 + 2*sAABBC*sAABCC + sAABCC**2 + 2*sAABBC*sAABCD + 2*sAABCC*sAABCD + sAABCD**2 + 2*sAAABB*sABCDD + 2*sAAABC*sABCDD + 2*sAAABB*sABCDE + 2*sAAABC*sABCDE)
  alpha_5  = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABBC + 2*sAAABC*sAABBC + 2*sAAABB*sAABCC + 2*sAAABC*sAABCC + 2*sAAABB*sAABCD + 2*sAAABC*sAABCD + 2*sAAAAB*sABCDD + 2*sAAAAB*sABCDE)
  alpha_6  = d*(sAAABB**2 + 2*sAAABB*sAAABC + sAAABC**2 + 2*sAAAAB*sAABBC + 2*sAAAAB*sAABCC + 2*sAAAAB*sAABCD + 2*sAAAAA*sABCDD + 2*sAAAAA*sABCDE)
  alpha_7  = d*(2*sAAAAB*sAAABB + 2*sAAAAB*sAAABC + 2*sAAAAA*sAABBC + 2*sAAAAA*sAABCC + 2*sAAAAA*sAABCD)
  alpha_8  = d*(sAAAAB**2 + 2*sAAAAA*sAAABB + 2*sAAAAA*sAAABC)
  alpha_9  = d*(2*sAAAAA*sAAAAB)
  alpha_10 = d*(sAAAAA**2)
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

#' Produce model estimated (p=5, full model, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaab,raaabb,raaabc,raabbc,raabcc,raabcd,rabcdd,rabcde Numerics corresponding to the nucleotide heterozygosities.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_0_unique = function(raaaab, raaabb, raaabc, raabbc, raabcc, raabcd, rabcdd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaaab-raaabb-raaabc-raabbc-raabcc-raabcd-rabcdd-rabcde
  if (raaaaa < 0) {return(0)}
  tAAAAA = (raaaaa)**k
  tAAAAB = (raaaaa+raaaab)**k
  tAAABB = (raaaaa+raaabb)**k
  tAAABC = (raaaaa+raaaab+raaabb+raaabc)**k
  tAABBC = (raaaaa+raaaab+raabbc)**k
  tAABCC = (raaaaa+raaabb+raabcc)**k
  tAABCD = (raaaaa+raaaab+raaabb+raaabc+raabbc+raabcc+raabcd)**k
  tABCDD = (raaaaa+raaabb+raabcc+rabcdd)**k
  sAAAAA = tAAAAA
  sAAAAB = tAAAAB-tAAAAA
  sAAABB = tAAABB-tAAAAA
  sAAABC = tAAABC-tAAAAB-tAAABB+tAAAAA
  sAABBC = tAABBC-tAAAAB
  sAABCC = tAABCC-tAAABB
  sAABCD = tAABCD-tAAABC-tAABBC-tAABCC+tAAAAB+tAAABB
  sABCDD = tABCDD-tAABCC
  sABCDE = 1-tAABCD-tABCDD+tAABCC
  alpha_1_unique  = (1-d)*(sAAAAB + 2*sAAABC + sAABBC + sAABCC + 3*sAABCD + 3*sABCDD + 5*sABCDE)
  alpha_2_unique  = (1-d)*(sAAABB + 2*sAABBC + 2*sAABCC + sAABCD + sABCDD)
  alpha_3_unique  = (1-d)*(sAAABB + sAAABC)
  alpha_4_unique  = (1-d)*(sAAAAB)
  alpha_5_unique  = (1-d)*(sAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
}

#' Produce model estimated (p=5, topology=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaab,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aaabc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_1 = function(raaaab, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaaab-raaabc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAAAB = raaaab
    sAAABC = raaabc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAAABC = (raaaaa+raaaab+raaabc)**k
    tAABCD = (raaaaa+raaaab+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAAABC = tAAABC-tAAAAB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
  }
  alpha_1 = (1-d)*(sAAAAB + 2*sAAABC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAAAAB + 2*sAAAAB**2 + 4*sAAAAA*sAAABC + 6*sAAAAB*sAAABC + 4*sAAABC**2 + 6*sAAAAA*sAABCD + 8*sAAAAB*sAABCD + 10*sAAABC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 10*sAAAAB*sABCDE + 12*sAAABC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
  alpha_2 = (1-d)*(sAABCD) + d*(sABCDE**2)
  alpha_3 = (1-d)*(sAAABC) + d*(2*sAABCD*sABCDE)
  alpha_4 = (1-d)*(sAAAAB) + d*(sAABCD**2 + 2*sAAABC*sABCDE)
  alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABC*sAABCD + 2*sAAAAB*sABCDE)
  alpha_6 = d*(sAAABC**2 + 2*sAAAAB*sAABCD + 2*sAAAAA*sABCDE)
  alpha_7 = d*(2*sAAAAB*sAAABC + 2*sAAAAA*sAABCD)
  alpha_8 = d*(sAAAAB**2 + 2*sAAAAA*sAAABC)
  alpha_9 = d*(2*sAAAAA*sAAAAB)
  alpha_10= d*(sAAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
  alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
  alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
}

#' Produce model estimated (p=5, topology=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaab,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aaabc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_1_unique = function(raaaab, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaaab-raaabc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAAAB = raaaab
    sAAABC = raaabc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAAABC = (raaaaa+raaaab+raaabc)**k
    tAABCD = (raaaaa+raaaab+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAAABC = tAAABC-tAAAAB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
  }
  alpha_1_unique = (1-d)*(sAAAAB + 2*sAAABC + 3*sAABCD + 5*sABCDE)
  alpha_2_unique = (1-d)*(sAABCD)
  alpha_3_unique = (1-d)*(sAAABC)
  alpha_4_unique = (1-d)*(sAAAAB)
  alpha_5_unique = (1-d)*(sAAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
}

#' Produce model estimated (p=5, topology=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaab,raabbc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aabbc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_2 = function(raaaab, raabbc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaaab-raabbc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAAAB = raaaab
    sAABBC = raabbc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAABBC = (raaaaa+raaaab+raabbc)**k
    tAABCD = (raaaaa+raaaab+raabbc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAABBC = tAABBC-tAAAAB
    sAABCD = tAABCD-tAABBC
    sABCDE = 1-tAABCD
  }
  alpha_1 = (1-d)*(sAAAAB + sAABBC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAAAAB + 2*sAAAAB**2 + 2*sAAAAA*sAABBC + 4*sAAAAB*sAABBC + 2*sAABBC**2 + 6*sAAAAA*sAABCD + 8*sAAAAB*sAABCD + 8*sAABBC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 10*sAAAAB*sABCDE + 10*sAABBC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
  alpha_2 = (1-d)*(2*sAABBC + sAABCD) + d*(2*sAAAAA*sAABBC + 2*sAAAAB*sAABBC + 2*sAABBC**2 + 2*sAABBC*sAABCD + 2*sAABBC*sABCDE + sABCDE**2)
  alpha_3 = (1-d)*(0) + d*(2*sAABBC*sABCDE + 2*sAABCD*sABCDE)
  alpha_4 = (1-d)*(sAAAAB) + d*(sAABBC**2 + 2*sAABBC*sAABCD + sAABCD**2)
  alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAAAB*sABCDE)
  alpha_6 = d*(2*sAAAAB*sAABBC + 2*sAAAAB*sAABCD + 2*sAAAAA*sABCDE)
  alpha_7 = d*(2*sAAAAA*sAABBC + 2*sAAAAA*sAABCD)
  alpha_8 = d*(sAAAAB**2)
  alpha_9 = d*(2*sAAAAA*sAAAAB)
  alpha_10= d*(sAAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
  alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
  alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
}

#' Produce model estimated (p=5, topology=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaab,raabbc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aabbc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_2_unique = function(raaaab, raabbc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaaab-raabbc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAAAB = raaaab
    sAABBC = raabbc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAABBC = (raaaaa+raaaab+raabbc)**k
    tAABCD = (raaaaa+raaaab+raabbc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAABBC = tAABBC-tAAAAB
    sAABCD = tAABCD-tAABBC
    sABCDE = 1-tAABCD
  }
  alpha_1_unique = (1-d)*(sAAAAB + sAABBC + 3*sAABCD + 5*sABCDE)
  alpha_2_unique = (1-d)*(2*sAABBC + sAABCD)
  alpha_3_unique = (1-d)*(0)
  alpha_4_unique = (1-d)*(sAAAAB)
  alpha_5_unique = (1-d)*(sAAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
}

#' Produce model estimated (p=5, topology=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabb,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aaabc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_3 = function(raaabb, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaabb-raaabc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAABB = raaabb
    sAAABC = raaabc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAAABC = (raaaaa+raaabb+raaabc)**k
    tAABCD = (raaaaa+raaabb+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAAABC = tAAABC-tAAABB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
  }
  alpha_1 = (1-d)*(2*sAAABC + 3*sAABCD + 5*sABCDE) + d*(4*sAAAAA*sAAABC + 4*sAAABB*sAAABC + 4*sAAABC**2 + 6*sAAAAA*sAABCD + 6*sAAABB*sAABCD + 10*sAAABC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 8*sAAABB*sABCDE + 12*sAAABC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
  alpha_2 = (1-d)*(sAAABB + sAABCD) + d*(2*sAAAAA*sAAABB + 2*sAAABB**2 + 2*sAAABB*sAAABC + 2*sAAABB*sAABCD + 2*sAAABB*sABCDE + sABCDE**2)
  alpha_3 = (1-d)*(sAAABB + sAAABC) + d*(2*sAABCD*sABCDE)
  alpha_4 = (1-d)*(0) + d*(sAABCD**2 + 2*sAAABB*sABCDE + 2*sAAABC*sABCDE)
  alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABCD + 2*sAAABC*sAABCD)
  alpha_6 = d*(sAAABB**2 + 2*sAAABB*sAAABC + sAAABC**2 + 2*sAAAAA*sABCDE)
  alpha_7 = d*(2*sAAAAA*sAABCD)
  alpha_8 = d*(2*sAAAAA*sAAABB + 2*sAAAAA*sAAABC)
  alpha_9 = d*(0)
  alpha_10= d*(sAAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
  alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
  alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
}

#' Produce model estimated (p=5, topology=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabb,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aaabc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_3_unique = function(raaabb, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaabb-raaabc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAABB = raaabb
    sAAABC = raaabc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAAABC = (raaaaa+raaabb+raaabc)**k
    tAABCD = (raaaaa+raaabb+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAAABC = tAAABC-tAAABB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
  }
  alpha_1_unique = (1-d)*(2*sAAABC + 3*sAABCD + 5*sABCDE)
  alpha_2_unique = (1-d)*(sAAABB + sAABCD)
  alpha_3_unique = (1-d)*(sAAABB + sAAABC)
  alpha_4_unique = (1-d)*(0)
  alpha_5_unique = (1-d)*(sAAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
}

#' Produce model estimated (p=5, topology=4) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabb,raabcc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_4 = function(raaabb, raabcc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaabb-raabcc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAABB = raaabb
    sAABCC = raabcc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tAABCD = (raaaaa+raaabb+raabcc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sAABCD = tAABCD-tAABCC
    sABCDE = 1-tAABCD
  }
  alpha_1 = (1-d)*(sAABCC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAABCC + 2*sAAABB*sAABCC + 2*sAABCC**2 + 6*sAAAAA*sAABCD + 6*sAAABB*sAABCD + 8*sAABCC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 8*sAAABB*sABCDE + 10*sAABCC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
  alpha_2 = (1-d)*(sAAABB + 2*sAABCC + sAABCD) + d*(2*sAAAAA*sAAABB + 2*sAAABB**2 + 2*sAAAAA*sAABCC + 4*sAAABB*sAABCC + 2*sAABCC**2 + 2*sAAABB*sAABCD + 2*sAABCC*sAABCD + 2*sAAABB*sABCDE + 2*sAABCC*sABCDE + sABCDE**2)
  alpha_3 = (1-d)*(sAAABB) + d*(2*sAABCC*sABCDE + 2*sAABCD*sABCDE)
  alpha_4 = (1-d)*(0) + d*(sAABCC**2 + 2*sAABCC*sAABCD + sAABCD**2 + 2*sAAABB*sABCDE)
  alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABCC + 2*sAAABB*sAABCD)
  alpha_6 = d*(sAAABB**2 + 2*sAAAAA*sABCDE)
  alpha_7 = d*(2*sAAAAA*sAABCC + 2*sAAAAA*sAABCD)
  alpha_8 = d*(2*sAAAAA*sAAABB)
  alpha_9 = d*(0)
  alpha_10= d*(sAAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
  alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
  alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
}

#' Produce model estimated (p=5, topology=4, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabb,raabcc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, aabcd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_4_unique = function(raaabb, raabcc, raabcd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaabb-raabcc-raabcd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAABB = raaabb
    sAABCC = raabcc
    sAABCD = raabcd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tAABCD = (raaaaa+raaabb+raabcc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sAABCD = tAABCD-tAABCC
    sABCDE = 1-tAABCD
  }
  alpha_1_unique = (1-d)*(sAABCC + 3*sAABCD + 5*sABCDE)
  alpha_2_unique = (1-d)*(sAAABB + 2*sAABCC + sAABCD)
  alpha_3_unique = (1-d)*(sAAABB)
  alpha_4_unique = (1-d)*(0)
  alpha_5_unique = (1-d)*(sAAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
}

#' Produce model estimated (p=5, topology=5) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabb,raabcc,rabcdd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, abcdd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_5 = function(raaabb, raabcc, rabcdd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaabb-raabcc-rabcdd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAABB = raaabb
    sAABCC = raabcc
    sABCDD = rabcdd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tABCDD = (raaaaa+raaabb+raabcc+rabcdd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sABCDD = tABCDD-tAABCC
    sABCDE = 1-tABCDD
  }
  alpha_1 = (1-d)*(sAABCC + 3*sABCDD + 5*sABCDE) + d*(2*sAAAAA*sAABCC + 2*sAAABB*sAABCC + 2*sAABCC**2 + 4*sAAAAA*sABCDD + 4*sAAABB*sABCDD + 6*sAABCC*sABCDD + 4*sABCDD**2 + 8*sAAAAA*sABCDE + 8*sAAABB*sABCDE + 10*sAABCC*sABCDE + 12*sABCDD*sABCDE + 8*sABCDE**2)
  alpha_2 = (1-d)*(sAAABB + 2*sAABCC + sABCDD) + d*(2*sAAAAA*sAAABB + 2*sAAABB**2 + 2*sAAAAA*sAABCC + 4*sAAABB*sAABCC + 2*sAABCC**2 + 2*sAAAAA*sABCDD + 4*sAAABB*sABCDD + 4*sAABCC*sABCDD + 3*sABCDD**2 + 2*sAAABB*sABCDE + 2*sAABCC*sABCDE + 4*sABCDD*sABCDE + sABCDE**2)
  alpha_3 = (1-d)*(sAAABB) + d*(2*sAABCC*sABCDD + 2*sAABCC*sABCDE)
  alpha_4 = (1-d)*(0) + d*(sAABCC**2 + 2*sAAABB*sABCDD + 2*sAAABB*sABCDE)
  alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABCC)
  alpha_6 = d*(sAAABB**2 + 2*sAAAAA*sABCDD + 2*sAAAAA*sABCDE)
  alpha_7 = d*(2*sAAAAA*sAABCC)
  alpha_8 = d*(2*sAAAAA*sAAABB)
  alpha_9 = d*(0)
  alpha_10= d*(sAAAAA**2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
  alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
  alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
}

#' Produce model estimated (p=5, topology=5, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabb,raabcc,rabcdd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, abcdd, and abcde respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict5_5_unique = function(raaabb, raabcc, rabcdd, rabcde, k, d, kmercov, bias, x)
{
  raaaaa = 1-raaabb-raabcc-rabcdd-rabcde
  if (raaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAA = raaaaa
    sAAABB = raaabb
    sAABCC = raabcc
    sABCDD = rabcdd
    sABCDE = rabcde
  } else {
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tABCDD = (raaaaa+raaabb+raabcc+rabcdd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sABCDD = tABCDD-tAABCC
    sABCDE = 1-tABCDD
  }
  alpha_1_unique = (1-d)*(sAABCC + 3*sABCDD + 5*sABCDE)
  alpha_2_unique = (1-d)*(sAAABB + 2*sAABCC + sABCDD)
  alpha_3_unique = (1-d)*(sAAABB)
  alpha_4_unique = (1-d)*(0)
  alpha_5_unique = (1-d)*(sAAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
}

#' Produce model estimated (p=6, full model) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaaabb,raaabbb,raaaabc,raaabbc,raaabcc,raabbcc,raaabcd,raabbcd,raabccd,raabcdd,raabcde,rabcdde,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_0 = function(raaaaab, raaaabb, raaabbb, raaaabc, raaabbc, raaabcc, raabbcc, raaabcd, raabbcd, raabccd, raabcdd, raabcde, rabcdde, rabcdee, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaaabb-raaabbb-raaaabc-raaabbc-raaabcc-raabbcc-raaabcd-raabbcd-raabccd-raabcdd-raabcde-rabcdde-rabcdee-rabcdef
  if (raaaaaa < 0) {return(0)}
  taaaaaa  = (raaaaaa)**k
  taaaaab  = (raaaaaa+raaaaab)**k
  taaaabb  = (raaaaaa+raaaabb)**k
  taaabbb  = (raaaaaa+raaabbb)**k
  taaaabc  = (raaaaaa+raaaaab+raaaabb+raaaabc)**k
  taaabbc  = (raaaaaa+raaaaab+raaabbb+raaabbc)**k
  taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
  taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
  taaabcd  = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raaabcc+raaabcd)**k
  taabbcd  = (raaaaaa+raaaaab+raaaabb+raaaabc+raabbcc+raabbcd)**k
  taabccd  = (raaaaaa+raaaaab+raaabbb+raaabbc+raabccd)**k
  taabcdd  = (raaaaaa+raaaabb+raaabcc+raabbcc+raabcdd)**k
  taabcde  = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raaabcc+raabbcc+raaabcd+raabbcd+raabccd+raabcdd+raabcde)**k
  tabcdde  = (raaaaaa+raaaaab+raaabbb+raaabbc+raabccd+rabcdde)**k
  tabcdee  = (raaaaaa+raaaabb+raaabcc+raabbcc+raabcdd+rabcdee)**k
  sAAAAAA  = taaaaaa
  sAAAAAB  = taaaaab-taaaaaa
  sAAAABB  = taaaabb-taaaaaa
  sAAABBB  = taaabbb-taaaaaa
  sAAAABC  = taaaabc-taaaaab-taaaabb+taaaaaa
  sAAABBC  = taaabbc-taaaaab-taaabbb+taaaaaa
  sAAABCC  = taaabcc-taaaabb
  sAABBCC  = taabbcc-taaaabb
  sAAABCD  = taaabcd-taaaabc-taaabbc-taaabcc+taaaaab+taaaabb
  sAABBCD  = taabbcd-taaaabc-taabbcc+taaaabb
  sAABCCD  = taabccd-taaabbc
  sAABCDD  = taabcdd-taaabcc-taabbcc+taaaabb
  sAABCDE  = taabcde-taaabcd-taabbcd-taabccd-taabcdd+taaaabc+taaabbc+taaabcc+taabbcc-taaaabb
  sABCDDE  = tabcdde-taabccd
  sABCDEE  = tabcdee-taabcdd
  sABCDEF  = 1-taabcde-tabcdde-tabcdee+taabccd+taabcdd
  alpha_1  = (1-d)*(sAAAAAB + 2*sAAAABC + sAAABBC + sAAABCC + 3*sAAABCD + 2*sAABBCD + 2*sAABCCD + 2*sAABCDD + 4*sAABCDE + 4*sABCDDE + 4*sABCDEE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAB*sAAAABB + 4*sAAAAAA*sAAAABC + 6*sAAAAAB*sAAAABC + 4*sAAAABB*sAAAABC + 4*sAAAABC**2 + 2*sAAAAAB*sAAABBB + 4*sAAAABC*sAAABBB + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAAABB*sAAABBC + 6*sAAAABC*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAAABCC + 4*sAAAAAB*sAAABCC + 2*sAAAABB*sAAABCC + 6*sAAAABC*sAAABCC + 2*sAAABBB*sAAABCC + 4*sAAABBC*sAAABCC + 2*sAAABCC**2 + 6*sAAAAAA*sAAABCD + 8*sAAAAAB*sAAABCD + 6*sAAAABB*sAAABCD + 10*sAAAABC*sAAABCD + 6*sAAABBB*sAAABCD + 8*sAAABBC*sAAABCD + 8*sAAABCC*sAAABCD + 6*sAAABCD**2 + 2*sAAAAAB*sAABBCC + 4*sAAAABC*sAABBCC + 2*sAAABBC*sAABBCC + 2*sAAABCC*sAABBCC + 6*sAAABCD*sAABBCC + 4*sAAAAAA*sAABBCD + 6*sAAAAAB*sAABBCD + 4*sAAAABB*sAABBCD + 8*sAAAABC*sAABBCD + 4*sAAABBB*sAABBCD + 6*sAAABBC*sAABBCD + 6*sAAABCC*sAABBCD + 10*sAAABCD*sAABBCD + 4*sAABBCC*sAABBCD + 4*sAABBCD**2 + 4*sAAAAAA*sAABCCD + 6*sAAAAAB*sAABCCD + 4*sAAAABB*sAABCCD + 8*sAAAABC*sAABCCD + 4*sAAABBB*sAABCCD + 6*sAAABBC*sAABCCD + 6*sAAABCC*sAABCCD + 10*sAAABCD*sAABCCD + 4*sAABBCC*sAABCCD + 8*sAABBCD*sAABCCD + 4*sAABCCD**2 + 4*sAAAAAA*sAABCDD + 6*sAAAAAB*sAABCDD + 4*sAAAABB*sAABCDD + 8*sAAAABC*sAABCDD + 4*sAAABBB*sAABCDD + 6*sAAABBC*sAABCDD + 6*sAAABCC*sAABCDD + 10*sAAABCD*sAABCDD + 4*sAABBCC*sAABCDD + 8*sAABBCD*sAABCDD + 8*sAABCCD*sAABCDD + 4*sAABCDD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 8*sAAAABB*sAABCDE + 12*sAAAABC*sAABCDE + 8*sAAABBB*sAABCDE + 10*sAAABBC*sAABCDE + 10*sAAABCC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABBCC*sAABCDE + 12*sAABBCD*sAABCDE + 12*sAABCCD*sAABCDE + 12*sAABCDD*sAABCDE + 8*sAABCDE**2 + 6*sAAAAAA*sABCDDE + 8*sAAAAAB*sABCDDE + 6*sAAAABB*sABCDDE + 10*sAAAABC*sABCDDE + 6*sAAABBB*sABCDDE + 8*sAAABBC*sABCDDE + 8*sAAABCC*sABCDDE + 12*sAAABCD*sABCDDE + 6*sAABBCC*sABCDDE + 10*sAABBCD*sABCDDE + 10*sAABCCD*sABCDDE + 10*sAABCDD*sABCDDE + 14*sAABCDE*sABCDDE + 6*sABCDDE**2 + 6*sAAAAAA*sABCDEE + 8*sAAAAAB*sABCDEE + 6*sAAAABB*sABCDEE + 10*sAAAABC*sABCDEE + 6*sAAABBB*sABCDEE + 8*sAAABBC*sABCDEE + 8*sAAABCC*sABCDEE + 12*sAAABCD*sABCDEE + 6*sAABBCC*sABCDEE + 10*sAABBCD*sABCDEE + 10*sAABCCD*sABCDEE + 10*sAABCDD*sABCDEE + 14*sAABCDE*sABCDEE + 12*sABCDDE*sABCDEE + 6*sABCDEE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 10*sAAAABB*sABCDEF + 14*sAAAABC*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 12*sAAABCC*sABCDEF + 16*sAAABCD*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABBCD*sABCDEF + 14*sAABCCD*sABCDEF + 14*sAABCDD*sABCDEF + 18*sAABCDE*sABCDEF + 16*sABCDDE*sABCDEF + 16*sABCDEE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + sAAABBC + sAAABCC + 3*sAABBCC + 2*sAABBCD + 2*sAABCCD + 2*sAABCDD + sAABCDE + sABCDDE + sABCDEE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAAAB*sAAAABB + 2*sAAAABB**2 + 2*sAAAABB*sAAAABC + 2*sAAAABB*sAAABBB + 2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 4*sAAAABB*sAAABBC + 2*sAAAABC*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAAABCC + 2*sAAAAAB*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAAABC*sAAABCC + 2*sAAABBB*sAAABCC + 4*sAAABBC*sAAABCC + 2*sAAABCC**2 + 2*sAAAABB*sAAABCD + 2*sAAABBC*sAAABCD + 2*sAAABCC*sAAABCD + 4*sAAAAAA*sAABBCC + 4*sAAAAAB*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAAAABC*sAABBCC + 4*sAAABBB*sAABBCC + 6*sAAABBC*sAABBCC + 6*sAAABCC*sAABBCC + 4*sAAABCD*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABBCD + 2*sAAAAAB*sAABBCD + 4*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAAABBB*sAABBCD + 4*sAAABBC*sAABBCD + 4*sAAABCC*sAABBCD + 2*sAAABCD*sAABBCD + 6*sAABBCC*sAABBCD + 2*sAABBCD**2 + 2*sAAAAAA*sAABCCD + 2*sAAAAAB*sAABCCD + 4*sAAAABB*sAABCCD + 2*sAAAABC*sAABCCD + 2*sAAABBB*sAABCCD + 4*sAAABBC*sAABCCD + 4*sAAABCC*sAABCCD + 2*sAAABCD*sAABCCD + 6*sAABBCC*sAABCCD + 4*sAABBCD*sAABCCD + 2*sAABCCD**2 + 2*sAAAAAA*sAABCDD + 2*sAAAAAB*sAABCDD + 4*sAAAABB*sAABCDD + 2*sAAAABC*sAABCDD + 2*sAAABBB*sAABCDD + 4*sAAABBC*sAABCDD + 4*sAAABCC*sAABCDD + 2*sAAABCD*sAABCDD + 6*sAABBCC*sAABCDD + 4*sAABBCD*sAABCDD + 4*sAABCCD*sAABCDD + 2*sAABCDD**2 + 2*sAAAABB*sAABCDE + 2*sAAABBC*sAABCDE + 2*sAAABCC*sAABCDE + 4*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAABCCD*sAABCDE + 2*sAABCDD*sAABCDE + 2*sAAAAAA*sABCDDE + 2*sAAAAAB*sABCDDE + 4*sAAAABB*sABCDDE + 2*sAAAABC*sABCDDE + 2*sAAABBB*sABCDDE + 4*sAAABBC*sABCDDE + 4*sAAABCC*sABCDDE + 2*sAAABCD*sABCDDE + 6*sAABBCC*sABCDDE + 4*sAABBCD*sABCDDE + 4*sAABCCD*sABCDDE + 4*sAABCDD*sABCDDE + 2*sAABCDE*sABCDDE + 3*sABCDDE**2 + 2*sAAAAAA*sABCDEE + 2*sAAAAAB*sABCDEE + 4*sAAAABB*sABCDEE + 2*sAAAABC*sABCDEE + 2*sAAABBB*sABCDEE + 4*sAAABBC*sABCDEE + 4*sAAABCC*sABCDEE + 2*sAAABCD*sABCDEE + 6*sAABBCC*sABCDEE + 4*sAABBCD*sABCDEE + 4*sAABCCD*sABCDEE + 4*sAABCDD*sABCDEE + 2*sAABCDE*sABCDEE + 6*sABCDDE*sABCDEE + 3*sABCDEE**2 + 2*sAAAABB*sABCDEF + 2*sAAABBC*sABCDEF + 2*sAAABCC*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + 2*sAABCCD*sABCDEF + 2*sAABCDD*sABCDEF + 4*sABCDDE*sABCDEF + 4*sABCDEE*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(2*sAAABBB + sAAABBC + sAAABCC + sAAABCD) + d*(2*sAAAAAA*sAAABBB + 2*sAAAAAB*sAAABBB + 2*sAAAABB*sAAABBB + 2*sAAAABC*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAAABCC + 2*sAAABBB*sAAABCD + 2*sAAABBB*sAABBCC + 2*sAAABBB*sAABBCD + 2*sAAABBB*sAABCCD + 2*sAAABBB*sAABCDD + 2*sAAABBB*sAABCDE + 2*sAAABBB*sABCDDE + 2*sAABBCC*sABCDDE + 2*sAABBCD*sABCDDE + 2*sAABCCD*sABCDDE + 2*sAABCDD*sABCDDE + 2*sAABCDE*sABCDDE + 2*sAAABBB*sABCDEE + 2*sAABBCC*sABCDEE + 2*sAABBCD*sABCDEE + 2*sAABCCD*sABCDEE + 2*sAABCDD*sABCDEE + 2*sAABCDE*sABCDEE + 2*sAAABBB*sABCDEF + 2*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + 2*sAABCCD*sABCDEF + 2*sAABCDD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB + sAAAABC) + d*(sAABBCC**2 + 2*sAABBCC*sAABBCD + sAABBCD**2 + 2*sAABBCC*sAABCCD + 2*sAABBCD*sAABCCD + sAABCCD**2 + 2*sAABBCC*sAABCDD + 2*sAABBCD*sAABCDD + 2*sAABCCD*sAABCDD + sAABCDD**2 + 2*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAABCCD*sAABCDE + 2*sAABCDD*sAABCDE + sAABCDE**2 + 2*sAAABBB*sABCDDE + 2*sAAABBC*sABCDDE + 2*sAAABCC*sABCDDE + 2*sAAABCD*sABCDDE + 2*sAAABBB*sABCDEE + 2*sAAABBC*sABCDEE + 2*sAAABCC*sABCDEE + 2*sAAABCD*sABCDEE + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF + 2*sAAABCC*sABCDEF + 2*sAAABCD*sABCDEF)
  alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBB*sAABBCC + 2*sAAABBC*sAABBCC + 2*sAAABCC*sAABBCC + 2*sAAABCD*sAABBCC + 2*sAAABBB*sAABBCD + 2*sAAABBC*sAABBCD + 2*sAAABCC*sAABBCD + 2*sAAABCD*sAABBCD + 2*sAAABBB*sAABCCD + 2*sAAABBC*sAABCCD + 2*sAAABCC*sAABCCD + 2*sAAABCD*sAABCCD + 2*sAAABBB*sAABCDD + 2*sAAABBC*sAABCDD + 2*sAAABCC*sAABCDD + 2*sAAABCD*sAABCDD + 2*sAAABBB*sAABCDE + 2*sAAABBC*sAABCDE + 2*sAAABCC*sAABCDE + 2*sAAABCD*sAABCDE + 2*sAAAABB*sABCDDE + 2*sAAAABC*sABCDDE + 2*sAAAABB*sABCDEE + 2*sAAAABC*sABCDEE + 2*sAAAABB*sABCDEF + 2*sAAAABC*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2 + 2*sAAABBB*sAAABCC + 2*sAAABBC*sAAABCC + sAAABCC**2 + 2*sAAABBB*sAAABCD + 2*sAAABBC*sAAABCD + 2*sAAABCC*sAAABCD + sAAABCD**2 + 2*sAAAABB*sAABBCC + 2*sAAAABC*sAABBCC + 2*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAAAABB*sAABCCD + 2*sAAAABC*sAABCCD + 2*sAAAABB*sAABCDD + 2*sAAAABC*sAABCDD + 2*sAAAABB*sAABCDE + 2*sAAAABC*sAABCDE + 2*sAAAAAB*sABCDDE + 2*sAAAAAB*sABCDEE + 2*sAAAAAB*sABCDEF)
  alpha_7  = d*(2*sAAAABB*sAAABBB + 2*sAAAABC*sAAABBB + 2*sAAAABB*sAAABBC + 2*sAAAABC*sAAABBC + 2*sAAAABB*sAAABCC + 2*sAAAABC*sAAABCC + 2*sAAAABB*sAAABCD + 2*sAAAABC*sAAABCD + 2*sAAAAAB*sAABBCC + 2*sAAAAAB*sAABBCD + 2*sAAAAAB*sAABCCD + 2*sAAAAAB*sAABCDD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDDE + 2*sAAAAAA*sABCDEE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAABB*sAAAABC + sAAAABC**2 + 2*sAAAAAB*sAAABBB + 2*sAAAAAB*sAAABBC + 2*sAAAAAB*sAAABCC + 2*sAAAAAB*sAAABCD + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCCD + 2*sAAAAAA*sAABCDD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAB*sAAAABB + 2*sAAAAAB*sAAAABC + 2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC + 2*sAAAAAA*sAAABCC + 2*sAAAAAA*sAAABCD)
  alpha_10 = d*(sAAAAAB**2 + 2*sAAAAAA*sAAAABB + 2*sAAAAAA*sAAAABC)
  alpha_11 = d*(2*sAAAAAA*sAAAAAB)
  alpha_12 = d*(sAAAAAA**2)
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


#' Produce model estimated (p=6, full model, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaaabb,raaabbb,raaaabc,raaabbc,raaabcc,raabbcc,raaabcd,raabbcd,raabccd,raabcdd,raabcde,rabcdde,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_0_unique = function(raaaaab, raaaabb, raaabbb, raaaabc, raaabbc, raaabcc, raabbcc, raaabcd, raabbcd, raabccd, raabcdd, raabcde, rabcdde, rabcdee, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaaabb-raaabbb-raaaabc-raaabbc-raaabcc-raabbcc-raaabcd-raabbcd-raabccd-raabcdd-raabcde-rabcdde-rabcdee-rabcdef
  if (raaaaaa < 0) {return(0)}
  taaaaaa  = (raaaaaa)**k
  taaaaab  = (raaaaaa+raaaaab)**k
  taaaabb  = (raaaaaa+raaaabb)**k
  taaabbb  = (raaaaaa+raaabbb)**k
  taaaabc  = (raaaaaa+raaaaab+raaaabb+raaaabc)**k
  taaabbc  = (raaaaaa+raaaaab+raaabbb+raaabbc)**k
  taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
  taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
  taaabcd  = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raaabcc+raaabcd)**k
  taabbcd  = (raaaaaa+raaaaab+raaaabb+raaaabc+raabbcc+raabbcd)**k
  taabccd  = (raaaaaa+raaaaab+raaabbb+raaabbc+raabccd)**k
  taabcdd  = (raaaaaa+raaaabb+raaabcc+raabbcc+raabcdd)**k
  taabcde  = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raaabcc+raabbcc+raaabcd+raabbcd+raabccd+raabcdd+raabcde)**k
  tabcdde  = (raaaaaa+raaaaab+raaabbb+raaabbc+raabccd+rabcdde)**k
  tabcdee  = (raaaaaa+raaaabb+raaabcc+raabbcc+raabcdd+rabcdee)**k
  sAAAAAA  = taaaaaa
  sAAAAAB  = taaaaab-taaaaaa
  sAAAABB  = taaaabb-taaaaaa
  sAAABBB  = taaabbb-taaaaaa
  sAAAABC  = taaaabc-taaaaab-taaaabb+taaaaaa
  sAAABBC  = taaabbc-taaaaab-taaabbb+taaaaaa
  sAAABCC  = taaabcc-taaaabb
  sAABBCC  = taabbcc-taaaabb
  sAAABCD  = taaabcd-taaaabc-taaabbc-taaabcc+taaaaab+taaaabb
  sAABBCD  = taabbcd-taaaabc-taabbcc+taaaabb
  sAABCCD  = taabccd-taaabbc
  sAABCDD  = taabcdd-taaabcc-taabbcc+taaaabb
  sAABCDE  = taabcde-taaabcd-taabbcd-taabccd-taabcdd+taaaabc+taaabbc+taaabcc+taabbcc-taaaabb
  sABCDDE  = tabcdde-taabccd
  sABCDEE  = tabcdee-taabcdd
  sABCDEF  = 1-taabcde-tabcdde-tabcdee+taabccd+taabcdd
  alpha_1_unique  = (1-d)*(sAAAAAB + 2*sAAAABC + sAAABBC + sAAABCC + 3*sAAABCD + 2*sAABBCD + 2*sAABCCD + 2*sAABCDD + 4*sAABCDE + 4*sABCDDE + 4*sABCDEE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + sAAABBC + sAAABCC + 3*sAABBCC + 2*sAABBCD + 2*sAABCCD + 2*sAABCDD + sAABCDE + sABCDDE + sABCDEE)
  alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC + sAAABCC + sAAABCD)
  alpha_4_unique  = (1-d)*(sAAAABB + sAAAABC)
  alpha_5_unique  = (1-d)*(sAAAAAB)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_1 = function(raaaaab, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaaabc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAAABC  = raaaabc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taaabcd  = (raaaaaa+raaaaab+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAAAAB + 2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 4*sAAAAAA*sAAAABC + 6*sAAAAAB*sAAAABC + 4*sAAAABC**2 + 6*sAAAAAA*sAAABCD + 8*sAAAAAB*sAAABCD + 10*sAAAABC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 12*sAAAABC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 14*sAAAABC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAABCDE) + d*(sABCDEF**2)
  alpha_3  = (1-d)*(sAAABCD) + d*(2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABC) + d*(sAABCDE**2 + 2*sAAABCD*sABCDEF)
  alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABCD*sAABCDE + 2*sAAAABC*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCD**2 + 2*sAAAABC*sAABCDE + 2*sAAAAAB*sABCDEF)
  alpha_7  = d*(2*sAAAABC*sAAABCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABC**2 + 2*sAAAAAB*sAAABCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAB*sAAAABC + 2*sAAAAAA*sAAABCD)
  alpha_10 = d*(sAAAAAB**2 + 2*sAAAAAA*sAAAABC)
  alpha_11 = d*(2*sAAAAAA*sAAAAAB)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_1_unique = function(raaaaab, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaaabc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAAABC  = raaaabc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taaabcd  = (raaaaaa+raaaaab+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
}
  alpha_1_unique  = (1-d)*(sAAAAAB + 2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAABCDE)
  alpha_3_unique  = (1-d)*(sAAABCD)
  alpha_4_unique  = (1-d)*(sAAAABC)
  alpha_5_unique  = (1-d)*(sAAAAAB)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_2 = function(raaaaab, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaaabc-raabbcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAAABC  = raaaabc
    sAABBCD  = raabbcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taabbcd  = (raaaaaa+raaaaab+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAAAAB + 2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 4*sAAAAAA*sAAAABC + 6*sAAAAAB*sAAAABC + 4*sAAAABC**2 + 4*sAAAAAA*sAABBCD + 6*sAAAAAB*sAABBCD + 8*sAAAABC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 12*sAAAABC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 14*sAAAABC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAABBCD + 2*sAAAAAB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAABBCD**2 + 2*sAABBCD*sAABCDE + 2*sAABBCD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(0) + d*(2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABC) + d*(sAABBCD**2 + 2*sAABBCD*sAABCDE + sAABCDE**2)
  alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAAABC*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABC*sAABBCD + 2*sAAAABC*sAABCDE + 2*sAAAAAB*sABCDEF)
  alpha_7  = d*(2*sAAAAAB*sAABBCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABC**2 + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAB*sAAAABC)
  alpha_10 = d*(sAAAAAB**2 + 2*sAAAAAA*sAAAABC)
  alpha_11 = d*(2*sAAAAAA*sAAAAAB)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_2_unique = function(raaaaab, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaaabc-raabbcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAAABC  = raaaabc
    sAABBCD  = raabbcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taabbcd  = (raaaaaa+raaaaab+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAAAAB + 2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(2*sAABBCD + sAABCDE)
  alpha_3_unique  = (1-d)*(0)
  alpha_4_unique  = (1-d)*(sAAAABC)
  alpha_5_unique  = (1-d)*(sAAAAAB)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_3 = function(raaaaab, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaabbc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAABBC  = raaabbc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taaabcd  = (raaaaaa+raaaaab+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAAAAB + sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 6*sAAAAAA*sAAABCD + 8*sAAAAAB*sAAABCD + 8*sAAABBC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 10*sAAABBC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 12*sAAABBC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAABBC + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 2*sAAABBC*sAAABCD + 2*sAAABBC*sAABCDE + 2*sAAABBC*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABBC + sAAABCD) + d*(2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(0) + d*(sAABCDE**2 + 2*sAAABBC*sABCDEF + 2*sAAABCD*sABCDEF)
  alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBC*sAABCDE + 2*sAAABCD*sAABCDE)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBC**2 + 2*sAAABBC*sAAABCD + sAAABCD**2 + 2*sAAAAAB*sABCDEF)
  alpha_7  = d*(2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(2*sAAAAAB*sAAABBC + 2*sAAAAAB*sAAABCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABBC + 2*sAAAAAA*sAAABCD)
  alpha_10 = d*(sAAAAAB**2)
  alpha_11 = d*(2*sAAAAAA*sAAAAAB)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_3_unique = function(raaaaab, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaabbc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAABBC  = raaabbc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taaabcd  = (raaaaaa+raaaaab+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAAAAB + sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAABBC + sAABCDE)
  alpha_3_unique  = (1-d)*(sAAABBC + sAAABCD)
  alpha_4_unique  = (1-d)*(0)
  alpha_5_unique  = (1-d)*(sAAAAAB)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=4) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_4 = function(raaaaab, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaabbc-raabccd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 6*sAAAAAB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 10*sAAABBC*sAABCDE + 12*sAABCCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAAAAB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAABBC*sAABCDE + 2*sAABCCD*sAABCDE + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABBC) + d*(2*sAABCCD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAABCCD*sAABCDE + sAABCDE**2 + 2*sAAABBC*sABCDEF)
  alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBC*sAABCCD + 2*sAAABBC*sAABCDE)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBC**2 + 2*sAAAAAB*sABCDEF)
  alpha_7  = d*(2*sAAAAAB*sAABCCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(2*sAAAAAB*sAAABBC + 2*sAAAAAA*sAABCCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABBC)
  alpha_10 = d*(sAAAAAB**2)
  alpha_11 = d*(2*sAAAAAA*sAAAAAB)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=4, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_4_unique = function(raaaaab, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaabbc-raabccd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE)
  alpha_3_unique  = (1-d)*(sAAABBC)
  alpha_4_unique  = (1-d)*(0)
  alpha_5_unique  = (1-d)*(sAAAAAB)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=5) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, abcdde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_5 = function(raaaaab, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaabbc-raabccd-rabcdde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sABCDDE  = rabcdde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaaaab+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
  }
  alpha_1  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 6*sAAAAAB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 6*sAAAAAA*sABCDDE + 8*sAAAAAB*sABCDDE + 8*sAAABBC*sABCDDE + 10*sAABCCD*sABCDDE + 6*sABCDDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 16*sABCDDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAAAAB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAAAAA*sABCDDE + 2*sAAAAAB*sABCDDE + 4*sAAABBC*sABCDDE + 4*sAABCCD*sABCDDE + 3*sABCDDE**2 + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + 4*sABCDDE*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABBC) + d*(2*sAABCCD*sABCDDE + 2*sAABCCD*sABCDEF)
  alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAAABBC*sABCDDE + 2*sAAABBC*sABCDEF)
  alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBC*sAABCCD)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBC**2 + 2*sAAAAAB*sABCDDE + 2*sAAAAAB*sABCDEF)
  alpha_7  = d*(2*sAAAAAB*sAABCCD + 2*sAAAAAA*sABCDDE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(2*sAAAAAB*sAAABBC + 2*sAAAAAA*sAABCCD)
  alpha_9  = d*(2*sAAAAAA*sAAABBC)
  alpha_10 = d*(sAAAAAB**2)
  alpha_11 = d*(2*sAAAAAA*sAAAAAB)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=5, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaaab,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, abcdde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_5_unique = function(raaaaab, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaaab-raaabbc-raabccd-rabcdde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAAAB  = raaaaab
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sABCDDE  = rabcdde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaaaab+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
  }
  alpha_1_unique  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE)
  alpha_3_unique  = (1-d)*(sAAABBC)
  alpha_4_unique  = (1-d)*(0)
  alpha_5_unique  = (1-d)*(sAAAAAB)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=6) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_6 = function(raaaabb, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaaabc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAAABC  = raaaabc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taaabcd  = (raaaaaa+raaaabb+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAAAABC + 4*sAAAABB*sAAAABC + 4*sAAAABC**2 + 6*sAAAAAA*sAAABCD + 6*sAAAABB*sAAABCD + 10*sAAAABC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 12*sAAAABC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 14*sAAAABC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAABB*sAAAABC + 2*sAAAABB*sAAABCD + 2*sAAAABB*sAABCDE + 2*sAAAABB*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABCD) + d*(2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB + sAAAABC) + d*(sAABCDE**2 + 2*sAAABCD*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABCD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAAABC*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCD**2 + 2*sAAAABB*sAABCDE + 2*sAAAABC*sAABCDE)
  alpha_7  = d*(2*sAAAABB*sAAABCD + 2*sAAAABC*sAAABCD + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAABB*sAAAABC + sAAAABC**2 + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABCD)
  alpha_10 = d*(2*sAAAAAA*sAAAABB + 2*sAAAAAA*sAAAABC)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=6, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_6_unique = function(raaaabb, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaaabc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAAABC  = raaaabc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taaabcd  = (raaaaaa+raaaabb+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + sAABCDE)
  alpha_3_unique  = (1-d)*(sAAABCD)
  alpha_4_unique  = (1-d)*(sAAAABB + sAAAABC)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=7) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_7 = function(raaaabb, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaaabc-raabbcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAAABC  = raaaabc
    sAABBCD  = raabbcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taabbcd  = (raaaaaa+raaaabb+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAAAABC + 4*sAAAABB*sAAAABC + 4*sAAAABC**2 + 4*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 8*sAAAABC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 12*sAAAABC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 14*sAAAABC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + 2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAABB*sAAAABC + 2*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAABBCD**2 + 2*sAAAABB*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAABBCD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(0) + d*(2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB + sAAAABC) + d*(sAABBCD**2 + 2*sAABBCD*sAABCDE + sAABCDE**2)
  alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEF + 2*sAAAABC*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAAAABB*sAABCDE + 2*sAAAABC*sAABCDE)
  alpha_7  = d*(2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAABB*sAAAABC + sAAAABC**2 + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(0)
  alpha_10 = d*(2*sAAAAAA*sAAAABB + 2*sAAAAAA*sAAAABC)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=7, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_7_unique = function(raaaabb, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaaabc-raabbcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAAABC  = raaaabc
    sAABBCD  = raabbcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taabbcd  = (raaaaaa+raaaabb+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + 2*sAABBCD + sAABCDE)
  alpha_3_unique  = (1-d)*(0)
  alpha_4_unique  = (1-d)*(sAAAABB + sAAAABC)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=8) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaabcc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_8 = function(raaaabb, raaabcc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaabcc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAABCC  = raaabcc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taaabcd  = (raaaaaa+raaaabb+raaabcc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAAABCD  = taaabcd-taaabcc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAABCC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABCC + 2*sAAAABB*sAAABCC + 2*sAAABCC**2 + 6*sAAAAAA*sAAABCD + 6*sAAAABB*sAAABCD + 8*sAAABCC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 10*sAAABCC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 12*sAAABCC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + sAAABCC + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAAAA*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAABCC**2 + 2*sAAAABB*sAAABCD + 2*sAAABCC*sAAABCD + 2*sAAAABB*sAABCDE + 2*sAAABCC*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAABCC*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABCC + sAAABCD) + d*(2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB) + d*(sAABCDE**2 + 2*sAAABCC*sABCDEF + 2*sAAABCD*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABCC*sAABCDE + 2*sAAABCD*sAABCDE + 2*sAAAABB*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCC**2 + 2*sAAABCC*sAAABCD + sAAABCD**2 + 2*sAAAABB*sAABCDE)
  alpha_7  = d*(2*sAAAABB*sAAABCC + 2*sAAAABB*sAAABCD + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABCC + 2*sAAAAAA*sAAABCD)
  alpha_10 = d*(2*sAAAAAA*sAAAABB)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=8, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaabcc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_8_unique = function(raaaabb, raaabcc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaabcc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAABCC  = raaabcc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taaabcd  = (raaaaaa+raaaabb+raaabcc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAAABCD  = taaabcd-taaabcc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAABCC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + sAAABCC + sAABCDE)
  alpha_3_unique  = (1-d)*(sAAABCC + sAAABCD)
  alpha_4_unique  = (1-d)*(sAAAABB)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=9) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaabcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_9 = function(raaaabb, raaabcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaabcc-raabcdd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAABCC  = raaabcc
    sAABCDD  = raabcdd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABCC + 2*sAAAABB*sAAABCC + 2*sAAABCC**2 + 4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAAABCC*sAABCDD + 4*sAABCDD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 10*sAAABCC*sAABCDE + 12*sAABCDD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 12*sAAABCC*sABCDEF + 14*sAABCDD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAAAA*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAABCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAAABCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAABB*sAABCDE + 2*sAAABCC*sAABCDE + 2*sAABCDD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAABCC*sABCDEF + 2*sAABCDD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABCC) + d*(2*sAABCDD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB) + d*(sAABCDD**2 + 2*sAABCDD*sAABCDE + sAABCDE**2 + 2*sAAABCC*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABCC*sAABCDD + 2*sAAABCC*sAABCDE + 2*sAAAABB*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCC**2 + 2*sAAAABB*sAABCDD + 2*sAAAABB*sAABCDE)
  alpha_7  = d*(2*sAAAABB*sAAABCC + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABCDD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABCC)
  alpha_10 = d*(2*sAAAAAA*sAAAABB)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=9, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaabcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_9_unique = function(raaaabb, raaabcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaabcc-raabcdd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAABCC  = raaabcc
    sAABCDD  = raabcdd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sAABCDE)
  alpha_3_unique  = (1-d)*(sAAABCC)
  alpha_4_unique  = (1-d)*(sAAAABB)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=10) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaabcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, abcdee, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_10 = function(raaaabb, raaabcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaabcc-raabcdd-rabcdee-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAABCC  = raaabcc
    sAABCDD  = raabcdd
    sABCDEE  = rabcdee
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raaabcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
  }
  alpha_1  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sABCDEE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABCC + 2*sAAAABB*sAAABCC + 2*sAAABCC**2 + 4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAAABCC*sAABCDD + 4*sAABCDD**2 + 6*sAAAAAA*sABCDEE + 6*sAAAABB*sABCDEE + 8*sAAABCC*sABCDEE + 10*sAABCDD*sABCDEE + 6*sABCDEE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 12*sAAABCC*sABCDEF + 14*sAABCDD*sABCDEF + 16*sABCDEE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sABCDEE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAAAA*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAABCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAAABCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAAAA*sABCDEE + 4*sAAAABB*sABCDEE + 4*sAAABCC*sABCDEE + 4*sAABCDD*sABCDEE + 3*sABCDEE**2 + 2*sAAAABB*sABCDEF + 2*sAAABCC*sABCDEF + 2*sAABCDD*sABCDEF + 4*sABCDEE*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(sAAABCC) + d*(2*sAABCDD*sABCDEE + 2*sAABCDD*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB) + d*(sAABCDD**2 + 2*sAAABCC*sABCDEE + 2*sAAABCC*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABCC*sAABCDD + 2*sAAAABB*sABCDEE + 2*sAAAABB*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCC**2 + 2*sAAAABB*sAABCDD)
  alpha_7  = d*(2*sAAAABB*sAAABCC + 2*sAAAAAA*sABCDEE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABCDD)
  alpha_9  = d*(2*sAAAAAA*sAAABCC)
  alpha_10 = d*(2*sAAAAAA*sAAAABB)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=10, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raaabcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, abcdee, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_10_unique = function(raaaabb, raaabcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raaabcc-raabcdd-rabcdee-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAAABCC  = raaabcc
    sAABCDD  = raabcdd
    sABCDEE  = rabcdee
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raaabcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
  }
  alpha_1_unique  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sABCDEE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sABCDEE)
  alpha_3_unique  = (1-d)*(sAAABCC)
  alpha_4_unique  = (1-d)*(sAAAABB)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=11) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raabbcc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_11 = function(raaaabb, raabbcc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raabbcc-raabbcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAABBCC  = raabbcc
    sAABBCD  = raabbcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabbcd  = (raaaaaa+raaaabb+raabbcc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABBCD  = taabbcd-taabbcc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 4*sAABBCC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 8*sAABBCC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 4*sAAAAAA*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 6*sAABBCC*sAABBCD + 2*sAABBCD**2 + 2*sAAAABB*sAABCDE + 4*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAAAABB*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(0) + d*(2*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB) + d*(sAABBCC**2 + 2*sAABBCC*sAABBCD + sAABBCD**2 + 2*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + sAABCDE**2)
  alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCC + 2*sAAAABB*sAABBCD + 2*sAAAABB*sAABCDE)
  alpha_7  = d*(2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(0)
  alpha_10 = d*(2*sAAAAAA*sAAAABB)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=11, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raabbcc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabbcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_11_unique = function(raaaabb, raabbcc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raabbcc-raabbcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAABBCC  = raabbcc
    sAABBCD  = raabbcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabbcd  = (raaaaaa+raaaabb+raabbcc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABBCD  = taabbcd-taabbcc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABBCD + sAABCDE)
  alpha_3_unique  = (1-d)*(0)
  alpha_4_unique  = (1-d)*(sAAAABB)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=12) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raabbcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_12 = function(raaaabb, raabbcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raabbcc-raabcdd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAABBCC  = raabbcc
    sAABCDD  = raabcdd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(2*sAABCDD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAABBCC*sAABCDD + 4*sAABCDD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 8*sAABBCC*sAABCDE + 12*sAABCDD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABCDD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 4*sAAAAAA*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAABBCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAABB*sAABCDE + 4*sAABBCC*sAABCDE + 2*sAABCDD*sAABCDE + 2*sAAAABB*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(0) + d*(2*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB) + d*(sAABBCC**2 + 2*sAABBCC*sAABCDD + sAABCDD**2 + 2*sAABBCC*sAABCDE + 2*sAABCDD*sAABCDE + sAABCDE**2)
  alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCC + 2*sAAAABB*sAABCDD + 2*sAAAABB*sAABCDE)
  alpha_7  = d*(2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABCDD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(0)
  alpha_10 = d*(2*sAAAAAA*sAAAABB)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=12, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raabbcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_12_unique = function(raaaabb, raabbcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raabbcc-raabcdd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAABBCC  = raabbcc
    sAABCDD  = raabcdd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(2*sAABCDD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sAABCDE)
  alpha_3_unique  = (1-d)*(0)
  alpha_4_unique  = (1-d)*(sAAAABB)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=13) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raabbcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, abcdee, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_13 = function(raaaabb, raabbcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raabbcc-raabcdd-rabcdee-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAABBCC  = raabbcc
    sAABCDD  = raabcdd
    sABCDEE  = rabcdee
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raabbcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
  }
  alpha_1  = (1-d)*(2*sAABCDD + 4*sABCDEE + 6*sABCDEF) + d*(4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAABBCC*sAABCDD + 4*sAABCDD**2 + 6*sAAAAAA*sABCDEE + 6*sAAAABB*sABCDEE + 6*sAABBCC*sABCDEE + 10*sAABCDD*sABCDEE + 6*sABCDEE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABCDD*sABCDEF + 16*sABCDEE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sABCDEE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 4*sAAAAAA*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAABBCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAAAA*sABCDEE + 4*sAAAABB*sABCDEE + 6*sAABBCC*sABCDEE + 4*sAABCDD*sABCDEE + 3*sABCDEE**2 + 2*sAAAABB*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF + 4*sABCDEE*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(0) + d*(2*sAABBCC*sABCDEE + 2*sAABCDD*sABCDEE + 2*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF)
  alpha_4  = (1-d)*(sAAAABB) + d*(sAABBCC**2 + 2*sAABBCC*sAABCDD + sAABCDD**2)
  alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEE + 2*sAAAABB*sABCDEF)
  alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCC + 2*sAAAABB*sAABCDD)
  alpha_7  = d*(2*sAAAAAA*sABCDEE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABCDD)
  alpha_9  = d*(0)
  alpha_10 = d*(2*sAAAAAA*sAAAABB)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=13, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaaabb,raabbcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, abcdee, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_13_unique = function(raaaabb, raabbcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaaabb-raabbcc-raabcdd-rabcdee-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAAABB  = raaaabb
    sAABBCC  = raabbcc
    sAABCDD  = raabcdd
    sABCDEE  = rabcdee
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raabbcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
  }
  alpha_1_unique  = (1-d)*(2*sAABCDD + 4*sABCDEE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sABCDEE)
  alpha_3_unique  = (1-d)*(0)
  alpha_4_unique  = (1-d)*(sAAAABB)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=14) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabbb,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_14 = function(raaabbb, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaabbb-raaabbc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAABBB  = raaabbb
    sAAABBC  = raaabbc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taaabcd  = (raaaaaa+raaabbb+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 6*sAAAAAA*sAAABCD + 6*sAAABBB*sAAABCD + 8*sAAABBC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAABBB*sAABCDE + 10*sAAABBC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAABBC + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAABBC*sAAABCD + 2*sAAABBC*sAABCDE + 2*sAAABBC*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(2*sAAABBB + sAAABBC + sAAABCD) + d*(2*sAAAAAA*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAAABCD + 2*sAAABBB*sAABCDE + 2*sAAABBB*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(0) + d*(sAABCDE**2 + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF + 2*sAAABCD*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABBB*sAABCDE + 2*sAAABBC*sAABCDE + 2*sAAABCD*sAABCDE)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2 + 2*sAAABBB*sAAABCD + 2*sAAABBC*sAAABCD + sAAABCD**2)
  alpha_7  = d*(2*sAAAAAA*sABCDEF)
  alpha_8  = d*(2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC + 2*sAAAAAA*sAAABCD)
  alpha_10 = d*(0)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=14, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabbb,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aaabcd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_14_unique = function(raaabbb, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaabbb-raaabbc-raaabcd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAABBB  = raaabbb
    sAAABBC  = raaabbc
    sAAABCD  = raaabcd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taaabcd  = (raaaaaa+raaabbb+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAABBC + sAABCDE)
  alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC + sAAABCD)
  alpha_4_unique  = (1-d)*(0)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=15) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabbb,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_15 = function(raaabbb, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaabbb-raaabbc-raabccd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAABBB  = raaabbb
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
  }
  alpha_1  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 4*sAAABBB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAABBB*sAABCDE + 10*sAAABBC*sAABCDE + 12*sAABCCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAABBB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAABBC*sAABCDE + 2*sAABCCD*sAABCDE + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(2*sAAABBB + sAAABBC) + d*(2*sAAAAAA*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAABCCD + 2*sAAABBB*sAABCDE + 2*sAAABBB*sABCDEF + 2*sAABCCD*sABCDEF + 2*sAABCDE*sABCDEF)
  alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAABCCD*sAABCDE + sAABCDE**2 + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABBB*sAABCCD + 2*sAAABBC*sAABCCD + 2*sAAABBB*sAABCDE + 2*sAAABBC*sAABCDE)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2)
  alpha_7  = d*(2*sAAAAAA*sABCDEF)
  alpha_8  = d*(2*sAAAAAA*sAABCCD + 2*sAAAAAA*sAABCDE)
  alpha_9  = d*(2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC)
  alpha_10 = d*(0)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=15, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabbb,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, aabcde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_15_unique = function(raaabbb, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaabbb-raaabbc-raabccd-raabcde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAABBB  = raaabbb
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sAABCDE  = raabcde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
  }
  alpha_1_unique  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE)
  alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC)
  alpha_4_unique  = (1-d)*(0)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}

#' Produce model estimated (p=6, topology=16) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabbb,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, abcdde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_16 = function(raaabbb, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaabbb-raaabbc-raabccd-rabcdde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAABBB  = raaabbb
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sABCDDE  = rabcdde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaabbb+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
  }
  alpha_1  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 4*sAAABBB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 6*sAAAAAA*sABCDDE + 6*sAAABBB*sABCDDE + 8*sAAABBC*sABCDDE + 10*sAABCCD*sABCDDE + 6*sABCDDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 16*sABCDDE*sABCDEF + 10*sABCDEF**2)
  alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAABBB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAAAAA*sABCDDE + 2*sAAABBB*sABCDDE + 4*sAAABBC*sABCDDE + 4*sAABCCD*sABCDDE + 3*sABCDDE**2 + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + 4*sABCDDE*sABCDEF + sABCDEF**2)
  alpha_3  = (1-d)*(2*sAAABBB + sAAABBC) + d*(2*sAAAAAA*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAABCCD + 2*sAAABBB*sABCDDE + 2*sAABCCD*sABCDDE + 2*sAAABBB*sABCDEF + 2*sAABCCD*sABCDEF)
  alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAAABBB*sABCDDE + 2*sAAABBC*sABCDDE + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF)
  alpha_5  = (1-d)*(0) + d*(2*sAAABBB*sAABCCD + 2*sAAABBC*sAABCCD)
  alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2)
  alpha_7  = d*(2*sAAAAAA*sABCDDE + 2*sAAAAAA*sABCDEF)
  alpha_8  = d*(2*sAAAAAA*sAABCCD)
  alpha_9  = d*(2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC)
  alpha_10 = d*(0)
  alpha_11 = d*(0)
  alpha_12 = d*(sAAAAAA**2)
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

#' Produce model estimated (p=6, topology=16, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
#'
#' @param raaabbb,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, abcdde, and abcdef respectively.
#' @param k An integer corresponding to the kmer length.
#' @param d A numeric corresponding to the repetitiveness.
#' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
#' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
#' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
#' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
#' @export
predict6_16_unique = function(raaabbb, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
{
  raaaaaa  = 1-raaabbb-raaabbc-raabccd-rabcdde-rabcdef
  if (raaaaaa < 0) {return(0)}
  if (KMER_RATES) {
    sAAAAAA  = raaaaaa
    sAAABBB  = raaabbb
    sAAABBC  = raaabbc
    sAABCCD  = raabccd
    sABCDDE  = rabcdde
    sABCDEF  = rabcdef
  } else {
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaabbb+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
  }
  alpha_1_unique  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF)
  alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE)
  alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC)
  alpha_4_unique  = (1-d)*(0)
  alpha_5_unique  = (1-d)*(0)
  alpha_6_unique  = (1-d)*(sAAAAAA)
  alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
}
