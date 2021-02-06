

############### Function to compute the kernel matrix
K_matrix = function(X, h)
{
  
  return(exp(-as.matrix(dist(X))^2/(2*h^2)))
  
}


sim.fcn = function(X, Y, Z, B, dof = 5, tp = tp, tn = tn, ngrid, cores)
{
  m = nrow(Y)
  
  n = ncol(Y)
  if(length(X) != n) return("Error: X and Y dimensions don't match")

############### Additional Samples (Generating Xb and an additional B samples for ordering)
BZ <- bs(Z, df=dof, intercept = FALSE)
Lm <- lm(X ~ BZ)
res <- resid(Lm)
index <- sapply(rep(n,2*B),function(x) sample(x))
Xb  <- matrix(res[index], n, 2*B) + Lm$fitted.values 




#############Compute T1.cor, T2.cor, T1.hsic, T2.hsic
T1.cor <- T2.cor <- T1.hsic <- T2.hsic <- vector()
T1b.cor <- T2b.cor <- T1b.hsic <- T2b.hsic <- matrix(NA, nrow = m , ncol = 2*B) 

resX <- scale(resid(lm(X ~ BZ)))
resXb <- scale(resid(lm(Xb ~ BZ)))  
H = diag(rep(1,n)) - 1/n * matrix(rep(1, n^2), nrow = n, ncol = n)
hx = median(as.matrix(dist(X)))
KX = H %*% K_matrix(X, hx) %*% H
hz = median(as.matrix(dist(Z)))
hxz = median(as.matrix(dist(cbind(X,Z))))
KZ = H%*%K_matrix(Z, hz) %*% H
R = 0.001 * solve(KZ + 0.001*diag(rep(1,n)))
XZ = cbind(X,Z)
KXZ = R %*% H%*% K_matrix(XZ, hxz)%*%H %*% R


output = sapply(mclapply(1:m, function(j){

    y = Y[j,]
  resY <- scale(resid(lm(y ~ BZ))) 
  T1.cor[j] <- abs(mean(resY*resX))
  T2.cor[j] <- abs(cor(X,y))
  T1b.cor[j,] <- abs(colMeans(as.vector(resY) * resXb))
  T2b.cor[j,] <- abs(cor(Xb,y))
  
  hy = median(as.matrix(Y[j,]))
  KY = H %*% K_matrix(Y[j,], hy) %*% H
  KYZ = R %*% KY %*% R
  
  T1.hsic[j] <- (1/n)*sum(KXZ*KYZ)
  T2.hsic[j] <- (1/n)*sum(KX* KY)
  # mid = Sys.time()
  for(b in 1:(2*B))
  {
    hxzb = median(as.matrix(dist(cbind(Xb[,b],Z))))  ########OIII. CHECK THIS BIT BEFORE COPYING
    XZb = cbind(Xb[,b],Z)
    KXZb =  R %*% H%*% K_matrix(XZb, hxzb)%*%H %*% R
    KXb = H %*%K_matrix(Xb[,b], hx) %*% H
    T1b.hsic[j,b] <- (1/n)*sum(KXZb* KYZ)
    T2b.hsic[j,b] <- (1/n)*sum(KXb*KY)
  }
  
  T1b.hsic <-c(T1b.hsic[j,], T1.hsic[j])
  T2b.hsic <-c(T2b.hsic[j,], T2.hsic[j])
  
  T1b.cor <- c(T1b.cor[j,], T1.cor[j])
  T2b.cor <- c(T2b.cor[j,], T2.cor[j])
  return(c(T1b.hsic, T2b.hsic, T1b.cor, T2b.cor))
  
}, mc.cores = 1), 'c') 

B1 = 2*B
T1b.hsic = t(output[1:B1,])
T2b.hsic = t(output[(B1+2): (2*B1 + 1),])
T1b.cor = t(output[(2*B1+3):(3*B1+2),])
T2b.cor = t(output[(3*B1 + 4): (4*B1 + 3),])

T1.hsic = output[B1+1,]
T2.hsic = output[2*B1 + 2,]
T1.cor = output[3*B1 + 3,]
T2.cor = output[4*B1 + 4,]
  
  
 
 


T1b.hsic <-cbind(T1b.hsic, T1.hsic)
T2b.hsic <-cbind(T2b.hsic, T2.hsic)

T1b.cor <- cbind(T1b.cor, T1.cor)
T2b.cor <- cbind(T2b.cor, T2.cor)


###############1. ORDERING ALGO WITH DATA POINTS
############### COMPUTE F#
F_matrix = matrix(nrow = m, ncol = 2)
for(i in 1:m)
{
  F_matrix[i,1] = mean((T1b.hsic[,51:100] > T1.hsic[i]) * (T2b.hsic[,51:100] > T2.hsic[i]))
  F_matrix[i,2] = mean((T1b.cor[,51:100]> T1.cor[i]) * (T2b.cor[,51:100]> T2.cor[i]) )
    
}


###############SORT INDICES ACCORDING TO F#
ind.hsic = m - rank(F_matrix[,1], ties.method = "random") + 1
ind.cor = m - rank(F_matrix[,2], ties.method = "random") + 1

############### SELECT INDICES
k = 1
l = 1
fdp_init_hsic = 2
fdp_init_cor = 2
while(fdp_init_hsic > 0.05 || fdp_init_cor > 0.05)
{
  if(fdp_init_hsic > 0.05)
  {
    t1 = T1.hsic[ind.hsic == k]
    t2 = T2.hsic[ind.hsic == k]
    num = mean((T1b.hsic[,1:50]>t1) * (T2b.hsic[,1:50]>t2))
    den = max(c(1/m,mean((T1.hsic > t1)*(T2.hsic > t2))))
    fdp_init_hsic = num/den  ##put the fdp calculation in here. It will make code faster
    k = k + 1
  }
  
  if(fdp_init_cor > 0.05)
  {
    t1 = T1.cor[ind.cor == l]
    t2 = T2.cor[ind.cor== l]
    fdp_init_cor = mean((T1b.cor[,1:50]>t1) * (T2b.cor[,1:50]>t2))/max(c(1/m,mean((T1.cor > t1)*(T2.cor > t2))))
    
    l = l+1
  }
  
  
}

index.hsic = k 
index.cor = l 


#########Compute fdr and power
t1_star_hsic = T1.hsic[ind.hsic == index.hsic]
t2_star_hsic = T2.hsic[ind.hsic == index.hsic]

t1_star_cor =  T1.cor[ind.cor == index.cor]
t2_star_cor = T2.cor[ind.cor == index.cor]

rej.hsic = (T1.hsic > t1_star_hsic)* (T2.hsic > t2_star_hsic)
rej.cor = (T1.cor > t1_star_cor) * (T2.cor > t2_star_cor)

pow.hsic <- sum(tp*rej.hsic)/sum(tp)
fdr.hsic <- sum(tn*rej.hsic)/max(1,sum(rej.hsic)) 

pow.cor <- sum(tp*rej.cor)/sum(tp)
fdr.cor <- sum(tn*rej.cor)/max(1,sum(rej.cor)) 



##################2. ORDERING ALGO WITH EXTERNAL GRID

##########Compute F#
f_comp = function(t1c, t2c, method = "hsic")
{
  if(method == "hsic")
  {
    f = mean((T1b.hsic[,51:100] > t1) * (T2b.hsic[,51:100] > t2)) 
  }
    
  if(method == "cor")
  {
    f = mean((T1b.cor[,51:100] > t1) * (T2b.cor[,51:100] > t2)) 
  }
  
    
 
  return(f)
}


#DEFINING GRIDS FOR COR
t1c = seq(min(abs(T1.cor)), max(T1.cor), len = ngrid)
t2c = seq(min(abs(T2.cor)), max(T2.cor), len = ngrid)

#DEFINING GRIDS FOR HSIC

t1h = seq(min(abs(T1.hsic)), max(T1.hsic), len = ngrid)
t2h = seq(min(abs(T2.hsic)), max(T2.hsic), len = ngrid)


##Figure this out
registerDoMC(cores = 4)
rr <- foreach(t1 = rep(t1h, each = length(t2h)), t2 = rep(t2h, length(t1h)), .combine = cbind) %dopar% f_comp(t1, t2, "hsic")   
FDP <- unname(unlist(rr[1, ]))
NP <- unname(unlist(rr[2, ]))


F_matrix.h = matrix(nrow = ngrid, ncol = ngrid)
F_matrix.c = matrix(nrow = ngrid, ncol = ngrid)
for(i in 1:200)
{

for (j in 1:200)
{
  F_matrix[i,1] = mean((T1b.hsic[,51:100] > T1.hsic[i]) * (T2b.hsic[,51:100] > T2.hsic[i]))
  F_matrix[i,2] = mean((T1b.cor[,51:100]> T1.cor[i]) * (T2b.cor[,51:100]> T2.cor[i]) )
  
}

}




return(c(fdr.cor, fdr.hsic, pow.cor, pow.hsic))
}
