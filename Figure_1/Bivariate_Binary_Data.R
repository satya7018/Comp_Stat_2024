library(xlsx)

library(bindata)
 # correlation parameter
rho <- 0.2
# number of sample size change here for different n
n<-400

ncolumns <- 16
r<-10000
ns <- matrix(NA, r, ncolumns)
true<-c(0.2, 0.3, 0.5, -0.5)
et1_AA<-true[1]+true[2]+true[3]
p1_AA<-exp(et1_AA)/(1+exp(et1_AA))

et1_AB<-true[1]+true[2]+true[3]
p1_AB<-exp(et1_AB)/(1+exp(et1_AB))

et1_BA<-true[1]+true[2]-true[3]
p1_BA<-exp(et1_BA)/(1+exp(et1_BA))

et1_BB<-true[1]+true[2]-true[3]
p1_BB<-exp(et1_BB)/(1+exp(et1_BB))

et2_AA<-true[1]-true[2]+true[3]+true[4]
p2_AA<-exp(et2_AA)/(1+exp(et2_AA))

et2_AB<-true[1]-true[2]-true[3] +true[4]
p2_AB<-exp(et2_AB)/(1+exp(et2_AB))

et2_BA<-true[1]-true[2]+true[3]-true[4]
p2_BA<-exp(et2_BA)/(1+exp(et2_BA))

et2_BB<-true[1]-true[2]-true[3]-true[4]
p2_BB<-exp(et2_BB)/(1+exp(et2_BB))


  
for (t in 1:r) {
  ## Construct a binary correlation matrix
x1<-0
x2<-0
x3<-0
x4<-0

m <- matrix(c(1,rho,rho,1), ncol=2) 


n11_AA<-0
n01_AA<-0
n10_AA<-0
n00_AA<-0

n11_AB<-0
n01_AB<-0
n10_AB<-0
n00_AB<-0

n11_BA<-0
n01_BA<-0
n10_BA<-0
n00_BA<-0

n11_BB<-0
n01_BB<-0
n10_BB<-0
n00_BB<-0
c<-n/4
x1 <- rmvbin(n/4, margprob = c(p1_AA, p2_AA), bincorr = m) 
x2 <- rmvbin(n/4, margprob = c(p1_AB, p2_AB), bincorr = m)
x3 <- rmvbin(n/4, margprob = c(p1_BA, p2_BA), bincorr = m)
x4 <- rmvbin(n/4, margprob = c(p1_BB, p2_BB), bincorr = m)
#AA
for (i in 1:c) {
  if (x1[i,1]==1 && x1[i,2]==1)
    n11_AA<-n11_AA+1
}
for (i in 1:c) {
  if (x1[i,1]==0 && x1[i,2]==1)
    n01_AA =n01_AA+1
}
for (i in 1:c) {
  if (x1[i,1]==1 && x1[i,2]==0)
    n10_AA =n10_AA+1
}
for (i in 1:c) {
  if (x1[i,1]==0 && x1[i,2]==0)
    n00_AA =n00_AA+1
}
#AB
for (i in 1:c) {
  if (x2[i,1]==1 && x2[i,2]==1)
    n11_AB<-n11_AB+1
}
for (i in 1:c) {
  if (x2[i,1]==0 && x2[i,2]==1)
    n01_AB =n01_AB+1
}
for (i in 1:c) {
  if (x2[i,1]==1 && x2[i,2]==0)
    n10_AB =n10_AB+1
}
for (i in 1:c) {
  if (x2[i,1]==0 && x2[i,2]==0)
    n00_AB =n00_AB+1
}

#BA
for (i in 1:c) {
  if (x3[i,1]==1 && x3[i,2]==1)
    n11_BA<-n11_BA+1
}
for (i in 1:c) {
  if (x3[i,1]==0 && x3[i,2]==1)
    n01_BA =n01_BA+1
}
for (i in 1:c) {
  if (x3[i,1]==1 && x3[i,2]==0)
    n10_BA =n10_BA+1
}
for (i in 1:c) {
  if (x3[i,1]==0 && x3[i,2]==0)
    n00_BA =n00_BA+1
}

#BB
for (i in 1:c) {
  if (x4[i,1]==1 && x4[i,2]==1)
    n11_BB<-n11_BB+1
}
for (i in 1:c) {
  if (x4[i,1]==0 && x4[i,2]==1)
    n01_BB =n01_BB+1
}
for (i in 1:c) {
  if (x4[i,1]==1 && x4[i,2]==0)
    n10_BB =n10_BB+1
}
for (i in 1:c) {
  if (x4[i,1]==0 && x4[i,2]==0)
    n00_BB =n00_BB+1
}


ns[t , ]<-c(n11_AA, n01_AA, n10_AA, n00_AA, 
      n11_AB, n01_AB, n10_AB, n00_AB, 
      n11_BA, n01_BA, n10_BA, n00_BA, 
      n11_BB, n01_BB, n10_BB, n00_BB)
}
#Define the path of directory appropriately

write.xlsx(ns, '/Users/satya/Desktop/COST_Reproducible_Codes/Figure_1/Likelihood_data.xlsx')







