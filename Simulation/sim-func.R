#####################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances,
# which are used to support numerical simulation studies.
#####################################################################################


#####################################################################################
# Functions for generation of simulated data
# This section includes four functions: 
# generate.data()  tridiag.cor()  Power.law.network()  nearest.neighbor.network()
#####################################################################################
generate.data = function(N,Mu0.list,Theta0.list,Sigma0.list){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: generate.data
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the simulated data.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R packages: mvtnorm
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ N: K0 * 1 vector, the sample sizes of subgroups.
  ## @ Mu0.list: a list including K0 mean vectors (p * 1).
  ## @ Theta0.list: a list including K0 precision matrices (p * p).
  ## @ Sigma0.list: a list including K0 correlation matrices (p * p).
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "whole.data" including:
  ## @ L0: n * 1 vector, the subgroup labels to which each sample belongs.
  ## @ Mu0: K0 * p matrix, K0 mean vectors.
  ## @ Theta0: K0 * p * p array, K0 precision matrices.
  ## @ data: n * p matrix, the design matrix.
  ## @ n_all: int, the total sample size.
  ## @ K0: int, the true number of subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  K0 = length(Mu0.list)
  p = length(Mu0.list[[1]])
  Mu0=matrix(0,K0,p);L0=NULL;Theta0=array(0, dim = c(p, p, K0));data=NULL
  for (k in 1:K0) {
    Mu0[k,] <- Mu0.list[[k]]
    L0 <- c(L0,rep(k,N[k]))
  }
  for (k in 1:K0) {
    Theta0[,,k] <- as.matrix(Theta0.list[[k]])
  }
  for (k in 1:K0) {
    data <- rbind(data,mvrnorm(N[k],Mu0[k,],Sigma0.list[[k]]))
  }
  n_all = dim(data)[1]
  whole.data=list()
  whole.data$L0=L0
  whole.data$Mu0=Mu0
  whole.data$Theta0=Theta0
  whole.data$data=data
  whole.data$n_all=n_all
  whole.data$K0=K0
  return(whole.data)
}

generate.data.XZcor = function(n.XZ,N1,Mu0.list.X,Theta0.list.X,Sigma0.list.X,N2,Mu0.list.Z,Theta0.list.Z,Sigma0.list.Z,c=0.05){
  whole.data.X <- generate.data(N1,Mu0.list.X,Theta0.list.X,Sigma0.list.X)
  whole.data.Z <- generate.data(N2,Mu0.list.Z,Theta0.list.Z,Sigma0.list.Z)
  K2 <- length(N2)
  q <- length(Mu0.list.X[[1]])
  p <- length(Mu0.list.Z[[1]])
  data=NULL
  for (k in 1:K2) {
    k1 <- n.XZ[k,1]
    k2 <- n.XZ[k,2]
    Mu.XZ <- c(Mu0.list.X[[k1]],Mu0.list.Z[[k2]])
    Sigma.XZ <- as.matrix(bdiag(Sigma0.list.X[[k1]],Sigma0.list.Z[[k2]]))
    Sigma.XZ[1:q, (q+1):(q+p)] <- c*min(abs(Sigma.XZ[Sigma.XZ!=0]))
    Sigma.XZ[(q+1):(q+p), 1:q] <- c*min(abs(Sigma.XZ[Sigma.XZ!=0]))
    data <- rbind(data, mvrnorm(N2[k], Mu.XZ, Sigma.XZ))
  }
  data.X <- data[,1:q]
  data.Z <- data[,q+1:p]
  whole.data.X$data <- data.X
  whole.data.Z$data <- data.Z
  
  return(list(whole.data.X=whole.data.X, whole.data.Z=whole.data.Z))
}

tridiag.cor = function(p, s=4, rho=0.3){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: tridiag.cor
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the tri-diagonal precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ rho: The value of non-zero elements on non-diagonal elements.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A: The precision matrix.
  ## ---------------------------------------------------------------------------------------------------------------
  
  m=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  sig.G = matrix(0,m,m)
  for (j in 1:m) {
    ncol0 = c(j-1,j,j+1)
    sig = c(rho, 1, rho)
    ncol = ncol0[which(ncol0>0 & ncol0<=m)]
    sig = sig[which(ncol0>0 & ncol0<=m)]
    sig.G[j,ncol] = sig
  }
  submatrix=list()
  for (ss in 1:s) {
    submatrix[[ss]] = sig.G
  }
  A=submatrix[[1]]
  if(s > 1){
    for (ss in 2:s) {
      A=bdiag(A,submatrix[[ss]])
    }
  }
  return(A)
}

Power.law.network = function(p,s=10,umin=0.1,umax=0.4,I.list){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: Power.law.network
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the s-block power-law precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ umin: The lower bound of non-zero elements on non-diagonal elements.
  ## @ umax: The upper bound of non-zero elements on non-diagonal elements.
  ## @ I2: The replacement blocks for the precision matrix of the second subgroup.
  ## @ I3: The replacement blocks for the precision matrix of the third subgroup.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A list including The precision matrices of three subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=list()
  for (ss in 1:s) {
    g = sample_pa(pp, m=2,directed = F)
    Eg = as.data.frame(as_edgelist(g))
    subi = diag(1,pp)
    for (qq in 1:dim(Eg)[1]) {
      i=Eg[qq,1];j=Eg[qq,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    for (i in 1:pp) {
      subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.5
      # subi[i,i] = max(sum(abs(subi[i,setdiff(1:pp,i)]))+0.1,1)
    }
    submatrix[[ss]]=subi
  }
  A=submatrix[[1]]
  for (ss in 2:s) {
    A=bdiag(A,submatrix[[ss]])
  }
  A = as.matrix(A)
  
  matrix.list <- I.list
  A.list <- I.list
  matrix.list[[1]] <- submatrix
  A.list[[1]] <- A
  n.I <- length(I.list)
  for (ni in 2:n.I) {
    Ii <- I.list[[ni]]
    submatrixi=submatrix
    for (ii in 1:length(Ii)) {
      iii=Ii[ii]
      g = ba.game(pp, m=2,directed = F)
      Eg = as.data.frame(get.edgelist(g))
      subi = diag(1,pp)
      for (q in 1:dim(Eg)[1]) {
        i=Eg[q,1];j=Eg[q,2]
        ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
        subi[i,j]=ij;subi[j,i]=ij
      }
      for (i in 1:pp) {
        subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.5
        # subi[i,i] = max(sum(abs(subi[i,setdiff(1:pp,i)]))+0.1,1)
      }
      submatrixi[[iii]] = subi
    }
    matrix.list[[ni]] <- submatrixi
    
    A=submatrixi[[1]]
    for (ss in 2:s) {
      A=bdiag(A,submatrixi[[ss]])
    }
    A = as.matrix(A)
    A.list[[ni]] <- A
  }
  
  return(A.list)
}

nearest.neighbor.network = function(p,s=10,near=3,umin=0.1,umax=0.4,I.list){
  
  pp=p/s
  near = min(near,pp)
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=list()
  for (ss in 1:s) {
    subi = diag(1,pp,pp)
    point.matrix = matrix(runif(2*pp,0,1),pp,2)
    num.nearest = matrix(0,pp,near)
    for (j in 1:pp) {
      distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
      corj = setdiff(1:pp,j)[order(distancej)[1:near]]
      num.nearest[j,] = corj
    }
    for (oi in 1:dim(num.nearest)[1]) {
      for (oj in 1:dim(num.nearest)[2]) {
        i=oi;j=num.nearest[oi,oj]
        ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
        subi[i,j]=ij;subi[j,i]=ij
      }
    }
    for (i in 1:pp) {
      subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.5
    }
    submatrix[[ss]]=subi
  }
  A=submatrix[[1]]
  for (ss in 2:s) {
    A=bdiag(A,submatrix[[ss]])
  }
  A = as.matrix(A)
  
  matrix.list <- I.list
  A.list <- I.list
  matrix.list[[1]] <- submatrix
  A.list[[1]] <- A
  n.I <- length(I.list)
  for (ni in 2:n.I) {
    Ii <- I.list[[ni]]
    submatrixi=submatrix
    if(length(Ii)>1 | Ii[1]!=0){
      for (ii in 1:length(Ii)) {
        iii=Ii[ii]
        subi = diag(1,pp,pp)
        point.matrix = matrix(runif(2*pp,0,1),pp,2)
        num.nearest = matrix(0,pp,near)
        for (j in 1:pp) {
          distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
          corj = setdiff(1:pp,j)[order(distancej)[1:near]]
          num.nearest[j,] = corj
        }
        for (oi in 1:dim(num.nearest)[1]) {
          for (oj in 1:dim(num.nearest)[2]) {
            i=oi;j=num.nearest[oi,oj]
            ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
            subi[i,j]=ij;subi[j,i]=ij
          }
        }
        for (i in 1:pp) {
          subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.5
        }
        submatrixi[[iii]] = subi
      }
    }
    matrix.list[[ni]] <- submatrixi
    
    A=submatrixi[[1]]
    for (ss in 2:s) {
      A=bdiag(A,submatrixi[[ss]])
    }
    A = as.matrix(A)
    A.list[[ni]] <- A
  }
  
  
  return(A.list)
}


################################################################################
# Functions for evaluating performances of proposed methods
################################################################################
Esti.error = function(data, mu_hat, Theta_hat, Mu0, Theta0, K0, L.mat, L0, prob, memb=0){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Esti.error
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Gauging performance of the proposed and alternative approaches
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: 
  ##            R functions: f.den.vec() 
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: ( (q+1)*(p+1)-1 ) * s matrix, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ mu_hat: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups.
  ## @ Theta_hat: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ other input parameters: Similar to the previous definition.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## The vector including:
  ## @ K: The estimated number of subgroups.
  ## @ CE: The sub-grouping error
  ## @ CME: The mean squared error (MSE) for the mean vectors.
  ## @ PME: The mean squared error (MSE) for the precision matrices.
  ## @ TPR/FPR: The true and false positive rates for the off-diagonal elements of the precision matrices.
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = dim(data)[1]
  if(K_hat == K0){
    num = rep(0,K0)
    numk = NULL
    for (k in 1:K0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    mu_hat = mu_hat[num,]
    Theta_hat = Theta_hat[,,num]
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0)^2))/K0
    PME = sqrt(sum((Theta_hat - Theta0)^2))/K0
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0!=0) == 2,3,sum) - p) / (apply(Theta0!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K0
    FPk = apply((Theta_hat!=0) + (Theta0==0) == 2,3,sum) / apply(Theta0==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K0
    
    #  CE
    if(length(memb)>1)
    {
      member = memb
      } else if(length(L.mat) == 1){
      f.mat = matrix(0,n_all,K0)
      L.mat = matrix(0,n_all,K0)
      for(k.ind in 1:K_hat) {
        f.mat[,k.ind]=f.den.vec( data, as.numeric(mu_hat[k.ind,]), Theta_hat[,,k.ind] )
      }
      for(k.ind in 1:K0) {
        for(i in 1:n_all) {
          L.mat[i,k.ind] = prob[k.ind] * f.mat[i,k.ind] / prob %*% f.mat[i,]
        }
      }
        member = apply(L.mat,1,which.max)
      } else {
        member = apply(L.mat,1,which.max)
    }
    
    
    
    aa = L0
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  } else{
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      mu_hatk = mu_hat[k,]
      errork.mu = apply((t(Mu0) - mu_hatk)^2,2,sum)
      Theta_hatk = Theta_hat[,,k]
      errork.Theta = apply((Theta0 - rep(Theta_hatk,K0))^2,3,sum)
      errork = errork.mu + errork.Theta
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    Mu0.re = Mu0[num,]
    Theta0.re = Theta0[,,num]
    if(K_hat == 1){
      Mu0.re = as.matrix(t(Mu0.re))
      Theta0.re = as.array(Theta0.re)
      dim(Theta0.re) <- c(p,p,K_hat)
    }
    
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0.re)^2))/K_hat
    PME = sqrt(sum((Theta_hat - Theta0.re)^2))/K_hat
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0.re!=0) == 2,3,sum) - p) / (apply(Theta0.re!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K_hat
    FPk = apply((Theta_hat!=0) + (Theta0.re==0) == 2,3,sum) / apply(Theta0.re==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K_hat
    
    #  CE
    num = rep(0,K0)
    numk = NULL
    L.hat = rep(0,n_all)
    for (k in 1:K0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
      L.hat[which(L0 == k)] = numk
    }
    if(K_hat < K0){L.hat=L0}
    
    #  CE
    if(length(memb)>1){
      member = memb
    } else if(length(L.mat) == 1){
      f.mat = matrix(0,n_all,K_hat)
      L.mat = matrix(0,n_all,K_hat)
      for(k.ind in 1:K_hat) {
        f.mat[,k.ind]=f.den.vec( data, as.numeric(mu_hat[k.ind,]), Theta_hat[,,k.ind] )                  
      }
      for(k.ind in 1:K_hat) {
        for(i in 1:n_all) {
          L.mat[i,k.ind] = prob[k.ind] * f.mat[i,k.ind] / prob %*% f.mat[i,]
        }
      }
      
      member = apply(L.mat,1,which.max)
    } else {
      member = apply(L.mat,1,which.max)
    }
    
    
    aa = L.hat
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  }
  # index = as.data.frame(t(c(K_hat,CE,CME/sqrt(p),PME/p,TPR,FPR)))
  index = as.data.frame(t(c(K_hat,CE,CME,PME,TPR,FPR)))
  names(index) = c("K","CE","CME","PME","TPR","FPR")
  
  return(index)
}

sum.index = function(index.rep, K0=3, p=10){
  o = dim(index.rep)[1]
  sum.K = as.numeric(index.rep$K)
  sum.CE = as.numeric(index.rep$CE)
  sum.CME = as.numeric(index.rep$CME)
  sum.PME = as.numeric(index.rep$PME)
  sum.TPR = as.numeric(index.rep$TPR)
  sum.FPR = as.numeric(index.rep$FPR)
  
  per = sum(sum.K==K0)/o
  sum.K = round(c(mean(sum.K),sd(sum.K)),2)
  sum.CE = round(c(mean(sum.CE),sd(sum.CE)),3)
  sum.CME = round(c(mean(sum.CME),sd(sum.CME)),3)
  sum.PME = round(c(mean(sum.PME),sd(sum.PME)),3)
  sum.TPR = round(c(mean(sum.TPR),sd(sum.TPR)),3)
  sum.FPR = round(c(mean(sum.FPR),sd(sum.FPR)),3)
  sum.K = paste0(sum.K[1],"(",sum.K[2],")")
  sum.CE = paste0(sum.CE[1],"(",sum.CE[2],")")
  sum.CME = paste0(sum.CME[1],"(",sum.CME[2],")")
  sum.PME = paste0(sum.PME[1],"(",sum.PME[2],")")
  sum.TPR = paste0(sum.TPR[1],"(",sum.TPR[2],")")
  sum.FPR = paste0(sum.FPR[1],"(",sum.FPR[2],")")
  
  sum = as.data.frame(t(c(per,sum.K,sum.CE,sum.CME,sum.PME,sum.TPR,sum.FPR)))
  names(sum) = c("per","K","CE","CME","PME","TPR","FPR")
  return(sum)
  
}

sum.index.XZ = function(index.rep.XZ, K1=3, K2=5){
  ncol = dim(index.rep.XZ)[2]/2
  index.X = sum.index(index.rep.XZ[,1:ncol], K1)
  index.Z = sum.index(index.rep.XZ[,ncol+1:ncol], K2)
  sum = cbind(index.X, index.Z)
  return(sum)
}

sum.index.list = function(index.rep.list, K1=3, K2=5, row.names){
  n_row = length(index.rep.list)
  n_col = dim(sum.index.XZ(index.rep.list[[1]], K1=K1, K2=K2))[2]
  sum = as.data.frame(matrix(0,n_row,n_col))
  names(sum) = names(sum.index.XZ(index.rep.list[[1]], K1=K1, K2=K2))
  
  for (i in 1:n_row) {
    sum[i,] = sum.index.XZ(index.rep.list[[i]], K1=K1, K2=K2)
  }
  sum.kmean = sum[1:2,]
  sum.kmean[,c(5:7,12:14)] <- ""
  sum0 = rbind(sum.kmean,sum)
  
  row.names(sum0) = row.names
  return(sum0)
  
}

n1 <- 300               # The sample size of each subgroup
q <- 10                 # The dimension of the precision matrix X
K1 <- 3                 # The true number of subgroups
N1 <- c(100,200,200)    # The sample sizes of K0 subgroups
n2 <- 100               # The sample size of each sub-subgroup
p <- 100                # The dimension of the precision matrix Z
K2 <- 5                 # The true number of sub-subgroups
N2 <- rep(n2,K2)        # The sample sizes of K0 subgroups
K <- 8                  # The given upper bound of K0.
mu.list <- c(2,3)
network.type.list <- c("tridiagonal","nearest-neighbor","power-law")
n.XZ <- as.matrix(cbind(c(1,2,2,3,3),c(1,2,3,4,5)))

row.names <- c("kmeans.givenK","kmeans.selectK",
               "kmeans.givenK.JGL","kmeans.selectK.JGL",
               "pro","pro.refit.JGL","oracle")

# ------------- The propose method --------------
lambda <- genelambda.obo4(nlambda1=1,lambda1_max=0.2,lambda1_min=0.1,
                          nlambda2=5,lambda2_max=1.5,lambda2_min=0.1,
                          nlambda3=5,lambda3_max=4,lambda3_min=1,
                          nlambda4=2,lambda4_max=0.05,lambda4_min=0.01)
lambda.refit.X <- genelambda.obo4(nlambda2=10,lambda2_max=0.5,lambda2_min=0.2)
lambda.refit.Z <- genelambda.obo4(nlambda2=10,lambda2_max=0.15,lambda2_min=0.03)


#### alternatives #####
lambda.km.X <- genelambda.obo4(nlambda2=10,lambda2_max=0.6,lambda2_min=0.3)
lambda.km.Z <- genelambda.obo4(nlambda2=10,lambda2_max=0.15,lambda2_min=0.05)


