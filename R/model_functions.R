#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
predict1_unique = function(k, d, kmercov, bias, x)
{
  r0 = 1
  t0 = r0**k
  s0 = t0
  alpha_1_unique = (1-d)*(s0)
  alpha_1_unique * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)
}

#AB -> AA
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#AB -> AA
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#ABC -> AAB -> AAA
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#ABCD -> AABC -> (AAAB, AABB) -> AAAA
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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


#ACBD -> (ABAC, AABC, ACBB)-> (AAAB, ABAA, AABB) -> AAAA
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
predict4 = function(r1, r2, r3, k, d, kmercov, bias, x)
{
  #Am --(r1)--> Bm
  #Am --(r2)--> Ap
  #Bm --(r3)--> Bp
  if (r1 < r2 | r2 < r3) {return(0)}
  #if (r2 < r3) {return(0)}
  #if (r1 < 0 | r2 < 0 | r3 < 0 | d < 0) {return(0)}
  #if (r1 > 1 | r2 > 1 | r3 > 1 | d > 1) {return(0)}
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
  alpha_1 = (1-d)*(sAAAB + sABAA + 2*sABAC + 2*sAABC + 2*sACBB + 4*sACBD) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 2*sAAAB*sAABB + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABB*sAABC + 4*sAABC**2 + 2*sAAAA*sABAA + 4*sAAAB*sABAA + 2*sAABB*sABAA + 6*sAABC*sABAA + 2*sABAA**2 + 4*sAAAA*sABAC + 6*sAAAB*sABAC + 4*sAABB*sABAC + 8*sAABC*sABAC + 6*sABAA*sABAC + 4*sABAC**2 + 2*sAAAA*sACBB + 4*sAAAB*sACBB + 2*sAABB*sACBB + 6*sAABC*sACBB + 4*sABAA*sACBB + 6*sABAC*sACBB + 2*sACBB**2 + 6*sAAAA*sACBD + 8*sAAAB*sACBD + 6*sAABB*sACBD + 10*sAABC*sACBD + 8*sABAA*sACBD + 10*sABAC*sACBD + 8*sACBB*sACBD + 6*sACBD**2)
  alpha_2 = (1-d)*(2*sAABB + sABAC + sAABC + sACBB) + d*(2*sAAAA*sAABB + 2*sAAAB*sAABB + 2*sAABB**2 + 2*sAABB*sAABC + 2*sAABB*sABAA + 2*sAABB*sABAC + 2*sAAAA*sACBB + 2*sAAAB*sACBB + 4*sAABB*sACBB + 2*sAABC*sACBB + 2*sABAA*sACBB + 2*sABAC*sACBB + 3*sACBB**2 + 2*sAABB*sACBD + 4*sACBB*sACBD + sACBD**2)
  alpha_3 = (1-d)*(sAAAB + sABAA) + d*(2*sAABB*sACBB + 2*sAABC*sACBB + 2*sABAC*sACBB + 2*sAABB*sACBD + 2*sAABC*sACBD + 2*sABAC*sACBD)
  alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2 + 2*sAABB*sABAC + 2*sAABC*sABAC + sABAC**2 + 2*sAAAB*sACBB + 2*sABAA*sACBB + 2*sAAAB*sACBD + 2*sABAA*sACBD)
  alpha_5 = d*(2*sAAAB*sAABB + 2*sAAAB*sAABC + 2*sAABB*sABAA + 2*sAABC*sABAA + 2*sAAAB*sABAC + 2*sABAA*sABAC + 2*sAAAA*sACBB + 2*sAAAA*sACBD)
  alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABB + 2*sAAAA*sAABC + 2*sAAAB*sABAA + sABAA**2 + 2*sAAAA*sABAC)
  alpha_7 = d*(2*sAAAA*sAAAB + 2*sAAAA*sABAA)
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

#ACBD -> (ABAC, AABC, ACBB)-> (AAAB, ABAA, AABB) -> AAAA
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
predict4_unique = function(r1, r2, r3, k, d, kmercov, bias, x)
{
  #Am --(r1)--> Bm
  #Am --(r2)--> Ap
  #Bm --(r3)--> Bp
  if (r1 < r2 | r2 < r3) {return(0)}
  #if (r2 < r3) {return(0)}
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
  alpha_1_unique = (1-d)*(sAAAB + sABAA + 2*sABAC + 2*sAABC + 2*sACBB + 4*sACBD)
  alpha_2_unique = (1-d)*(2*sAABB + sABAC + sAABC + sACBB)
  alpha_3_unique = (1-d)*(sAAAB + sABAA)
  alpha_4_unique = (1-d)*(sAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
}

##ACBD -> (ABAC, AABC, ACBB)-> (AAAB, ABAA, AABB) -> AAAA
#  predict4 = function(sAAAA, sAAAB, sAABB, sAABC, k, d, kmercov, bias, x)
#  {
    #Am --(r1)--> Bm
    #Am --(r2)--> Ap
    #Bm --(r3)--> Bp
    #if (r1 < r2 | r2 < r3) {return(0)}
    #if (r2 < r3) {return(0)}
    #if (r1 < 0 | r2 < 0 | r3 < 0 | d < 0) {return(0)}
    #if (r1 > 1 | r2 > 1 | r3 > 1 | d > 1) {return(0)}
    #raaaa = (1-r1)*(1-r2)*(1-r3) #aaaa
    #raaab = (1-r1)*(1-r2)*(r3)
    #rabaa = (1-r1)*(r2)*(1-r3)
    #raabb = r1*(1-r2)*(1-r3)
    #rabac = (1-r1)*r2*r3
    #raabc = r1*(1-r2)*r3
    #racbb = r1*r2*(1-r3)
    #racbd = r1*r2*r3
    #tAAAA = (raaaa)**k
    #tAAAB = (raaab + raaaa)**k
    #tABAA = (rabaa + raaaa)**k
    #tAABB = (raabb + raaaa)**k
    #tABAC = (rabac + raaab + rabaa + raaaa)**k
    #tAABC = (raabc + raaab + raabb + raaaa)**k
    #tACBB = (racbb + rabaa + raabb + raaaa)**k
    #tACBD = (racbd + rabac + raabc + racbb + raaab + rabaa + raabb + raaaa)**k
    #sAAAA = tAAAA
    #sAAAB = tAAAB - tAAAA
    #sABAA = tABAA - tAAAA
    #sAABB = tAABB - tAAAA
    #sABAC = tABAC - tAAAB - tABAA + tAAAA
    #sAABC = tAABC - tAAAB - tAABB + tAAAA
    #sACBB = tACBB - tABAA - tAABB + tAAAA
    #sACBD = tACBD - tABAC - tAABC - tACBB + tAAAB + tABAA + tAABB - tAAAA
#    alpha_1 = (1-d)*(sAAAB + 2*sAABC) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 2*sAAAB*sAABB + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABB*sAABC + 4*sAABC**2)
#    alpha_2 = (1-d)*(2*sAABB + sAABC) + d*(2*sAAAA*sAABB + 2*sAAAB*sAABB + 2*sAABB**2 + 2*sAABB*sAABC)
#    alpha_3 = (1-d)*(sAAAB) + d*(0)
#    alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2)
#    alpha_5 = d*(2*sAAAB*sAABB + 2*sAAAB*sAABC)
#    alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABB + 2*sAAAA*sAABC)
#    alpha_7 = d*(2*sAAAA*sAAAB)
#    alpha_8 = d*(sAAAA**2)
#    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
#    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
#    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
#    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
#    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
#    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
#    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
#    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
#  }

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

#Am ->(r1)-> Ap ->(r2)-> Bm ->(r3)-> Bp
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
predict4 = function(r1, r2, r3, k, d, kmercov, bias, x)
{
  #enforce r1 <= r2 <= r3
  if (r1 > r2 || r2 > r3) {return(0)}
  #enforce r1 <= r3 <= r2
  #if (r1 > r3 || r3 > r2) {return(0)}
  raaaa = (1-r1)*(1-r2)*(1-r3)
  raaab = (1-r1)*(1-r2)*r3
  raabb = (1-r1)*r2*(1-r3)
  rabbb = r1*(1-r2)*(1-r3)
  raabc = (1-r1)*r2*r3
  rabbc = r1*(1-r2)*r3
  rabcc = r1*r2*(1-r3)
  rabcd = r1*r2*r3
  tAAAA = raaaa**k
  tAAAB = (raaaa+raaab)**k
  tAABB = (raaaa+raabb)**k
  tABBB = (raaaa+rabbb)**k
  tAABC = (raaaa+raaab+raabb+raabc)**k
  tABBC = (raaaa+raaab+rabbb+rabbc)**k
  tABCC = (raaaa+raabb+rabbb+rabcc)**k
  sAAAA = tAAAA
  sAAAB = tAAAB-tAAAA
  sAABB = tAABB-tAAAA
  sABBB = tABBB-tAAAA
  sAABC = tAABC-tAAAB-tAABB+tAAAA
  sABBC = tABBC-tAAAB-tABBB+tAAAA
  sABCC = tABCC-tAABB-tABBB+tAAAA
  sABCD = 1-tAABC-tABBC-tABCC+tAAAB+tAABB+tABBB-tAAAA
  alpha_1 = (1-d)*(sAAAB + 2*sAABC + sABBB + 2*sABBC + 2*sABCC + 4*sABCD) + d*(2*sAAAA*sAAAB + 2*sAAAB^2 + 2*sAAAB*sAABB + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABB*sAABC + 4*sAABC^2 + 2*sAAAB*sABBB + 4*sAABC*sABBB + 2*sAAAA*sABBC + 4*sAAAB*sABBC + 2*sAABB*sABBC + 6*sAABC*sABBC + 2*sABBB*sABBC + 2*sABBC^2 + 2*sAAAA*sABCC + 4*sAAAB*sABCC + 2*sAABB*sABCC + 6*sAABC*sABCC + 2*sABBB*sABCC + 4*sABBC*sABCC + 2*sABCC^2 + 6*sAAAA*sABCD + 8*sAAAB*sABCD + 6*sAABB*sABCD + 10*sAABC*sABCD + 6*sABBB*sABCD + 8*sABBC*sABCD + 8*sABCC*sABCD + 6*sABCD^2)
  alpha_2 = (1-d)*(2*sAABB + sAABC + sABBC + sABCC) + d*(2*sAAAA*sAABB + 2*sAAAB*sAABB + 2*sAABB^2 + 2*sAABB*sAABC + 2*sAABB*sABBB + sABBB^2 + 2*sAAAA*sABBC + 2*sAAAB*sABBC + 4*sAABB*sABBC + 2*sAABC*sABBC + 4*sABBB*sABBC + 3*sABBC^2 + 2*sAAAA*sABCC + 2*sAAAB*sABCC + 4*sAABB*sABCC + 2*sAABC*sABCC + 4*sABBB*sABCC + 6*sABBC*sABCC + 3*sABCC^2 + 2*sAABB*sABCD + 2*sABBB*sABCD + 4*sABBC*sABCD + 4*sABCC*sABCD + sABCD^2)
  alpha_3 = (1-d)*(sAAAB + sABBB) + d*(2*sAAAA*sABBB + 2*sAAAB*sABBB + 4*sAABB*sABBB + 4*sAABC*sABBB + 2*sABBB^2 + 2*sAABB*sABBC + 2*sAABC*sABBC + 2*sABBB*sABBC + 2*sAABB*sABCC + 2*sAABC*sABCC + 2*sABBB*sABCC + 2*sAABB*sABCD + 2*sAABC*sABCD + 2*sABBB*sABCD)
  alpha_4 = (1-d)*(sAAAA) + d*(sAABB^2 + 2*sAABB*sAABC + sAABC^2 + 2*sAAAB*sABBB + 2*sAAAB*sABBC + 2*sAAAB*sABCC + 2*sAAAB*sABCD)
  alpha_5 = d*(2*sAAAB*sAABB + 2*sAAAB*sAABC + 2*sAAAA*sABBB + 2*sAAAA*sABBC + 2*sAAAA*sABCC + 2*sAAAA*sABCD)
  alpha_6 = d*(sAAAB^2 + 2*sAAAA*sAABB + 2*sAAAA*sAABC)
  alpha_7 = d*(2*sAAAA*sAAAB)
  alpha_8 = d*(sAAAA^2)
  alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
  alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
  alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
  alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
  alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
}

#Am ->(r1)-> Ap ->(r2)-> Bm ->(r3)-> Bp
#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
predict4_unique = function(r1, r2, r3, k, d, kmercov, bias, x)
{
  #enforce r1 <= r2 <= r3
  if (r1 > r2 || r2 > r3) {return(0)}
  #enforce r1 <= r3 <= r2
  #if (r1 > r3 || r3 > r2) {return(0)}
  raaaa = (1-r1)*(1-r2)*(1-r3)
  raaab = (1-r1)*(1-r2)*r3
  raabb = (1-r1)*r2*(1-r3)
  rabbb = r1*(1-r2)*(1-r3)
  raabc = (1-r1)*r2*r3
  rabbc = r1*(1-r2)*r3
  rabcc = r1*r2*(1-r3)
  rabcd = r1*r2*r3
  tAAAA = raaaa**k
  tAAAB = (raaaa+raaab)**k
  tAABB = (raaaa+raabb)**k
  tABBB = (raaaa+rabbb)**k
  tAABC = (raaaa+raaab+raabb+raabc)**k
  tABBC = (raaaa+raaab+rabbb+rabbc)**k
  tABCC = (raaaa+raabb+rabbb+rabcc)**k
  sAAAA = tAAAA
  sAAAB = tAAAB-tAAAA
  sAABB = tAABB-tAAAA
  sABBB = tABBB-tAAAA
  sAABC = tAABC-tAAAB-tAABB+tAAAA
  sABBC = tABBC-tAAAB-tABBB+tAAAA
  sABCC = tABCC-tAABB-tABBB+tAAAA
  sABCD = 1-tAABC-tABBC-tABCC+tAAAB+tAABB+tABBB-tAAAA
  alpha_1_unique = (1-d)*(sAAAB + 2*sAABC + sABBB + 2*sABBC + 2*sABCC + 4*sABCD)
  alpha_2_unique = (1-d)*(2*sAABB + sAABC + sABBC + sABCC)
  alpha_3_unique = (1-d)*(sAAAB + sABBB)
  alpha_4_unique = (1-d)*(sAAAA)
  alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
  alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
  alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
  alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
}

#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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

#' something
#'
#' @param something
#' @return something
#' @seealso something
#' @export
#' @examples
#' blah blah blah
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
