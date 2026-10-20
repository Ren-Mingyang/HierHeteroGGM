
####################### alternative methods ############
kmeans_Hier_selectK <- function(data.X, data.Z, K_min=2,K_max=6, k.index= "gap")
  {
  n = dim(data.X)[1]
  q = dim(data.X)[2]
  p = dim(data.Z)[2]
  nb_clust = NbClust(data.X, distance = "euclidean",min.nc=K_min, max.nc=K_max, 
                     method = "kmeans",index = k.index, alphaBeale=0.1)
  member.X = nb_clust$Best.partition
  K.X = length(unique(member.X))
  Theta.X = array(0, dim = c(q, q, K.X))
  mu.hat.X = matrix(0, K.X, q)
  prob.X = rep(0,K.X)
  cluster.data.X = vector("list", length = K.X)
  member = as.data.frame(cbind(1:n, member.X, 0))
  names(member) = c("num","memb.X","memb.Z")
  member.Z = NULL
  S.K = 0
  for(k in 1:K.X) {
    cluster.data.X[[k]] = data.X[member.X==k,]
    mu.hat.X[k,] = apply(cluster.data.X[[k]],2,mean)
    prob.X[k] = sum(member.X == k)/length(member.X)
    
    nb_clust.Zk = NbClust(data.Z[member.X==k,], distance = "euclidean",min.nc=K_min, max.nc=K_max, 
                       method = "kmeans",index = k.index, alphaBeale=0.1)
    member.Zk = nb_clust.Zk$Best.partition
    K.Zk = length(unique(member.Zk))
    member.Zk = member.Zk + S.K
    member$memb.Z[member.X==k] = member.Zk
    S.K = K.Zk + S.K
  }
  member.list = list(member.table=member,member.X=member$memb.X,member.Z=member$memb.Z)
  return(member.list)
}

kmeans_Hier_givenK <- function(data.X, data.Z, K.X.k=2, K.Z.k=c(2,2))
  {
  n = dim(data.X)[1]
  q = dim(data.X)[2]
  p = dim(data.Z)[2]
  nb_clust = kmeans(data.X, K.X.k)
  member.X <- nb_clust$cluster
  member = as.data.frame(cbind(1:n, member.X, 0))
  names(member) = c("num","memb.X","memb.Z")
  member.Z = NULL
  S.K = 0
  for(k in 1:K.X.k) {
    nb_clust.Zk = kmeans(data.Z[member.X==k,], K.Z.k[k])
    member.Zk = nb_clust.Zk$cluster
    K.Zk = K.Z.k[k]
    member.Zk = member.Zk + S.K
    member$memb.Z[member.X==k] = member.Zk
    S.K = K.Zk + S.K
  }
  member.list = list(member.table=member,member.X=member$memb.X,member.Z=member$memb.Z)
  return(member.list)
}

kmeans.Hier.JGL.selectK <- function(data.X, data.Z, lambda.X, lambda.Z, K_min=2,K_max=6, k.index= "gap")
  {
  member.list = kmeans_Hier_selectK(data.X, data.Z, K_min=K_min,K_max=K_max, k.index= k.index)
  member.X = member.list$member.X
  member.Z = member.list$member.Z
  result.X = tuning.lambda.JGL(lambda.X, data.X, member.X)
  result.Z = tuning.lambda.JGL(lambda.Z, data.Z, member.Z)
  return(list(result.X=result.X, result.Z=result.Z))
}

kmeans.Hier.JGL.givenK <- function(data.X, data.Z, lambda.X, lambda.Z, K.X.k=2, K.Z.k=c(2,2))
  {
  member.list = kmeans_Hier_givenK(data.X, data.Z, K.X.k=K.X.k, K.Z.k=K.Z.k)
  member.X = member.list$member.X
  member.Z = member.list$member.Z
  result.X = tuning.lambda.JGL(lambda.X, data.X, member.X)
  result.Z = tuning.lambda.JGL(lambda.Z, data.Z, member.Z)
  return(list(result.X=result.X, result.Z=result.Z))
}

oracle.JGL <- function(whole.data.X, whole.data.Z, lambda.X, lambda.Z)
  {
  member.X=whole.data.X$L0
  data.X=whole.data.X$data
  member.Z=whole.data.Z$L0
  data.Z=whole.data.Z$data
  result.X = tuning.lambda.JGL(lambda.X, data.X, member.X)
  result.Z = tuning.lambda.JGL(lambda.Z, data.Z, member.Z)
  return(list(result.X=result.X, result.Z=result.Z))
}

Esti.error.JGL <- function(data, mu_hat, Theta_hat, Mu0, Theta0, K0, L0, member)
  {
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = dim(data)[1]
  if(K_hat == K0){
    num = rep(0,K0)
    numk = NULL
    L.hat = rep(0,n_all)
    for (k in 1:K0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
      L.hat[which(L0 == k)] = numk
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
  } else{
    # print("Failure to identify K0")
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      mu_hatk = mu_hat[k,]
      errork.mu = apply((t(Mu0) - mu_hatk)^2,2,sum)
      Theta_hatk = Theta_hat[,,k]
      errork.Theta = apply((Theta0 - rep(Theta_hatk,K0))^2,3,sum)
      errork = errork.mu + errork.Theta
      numk = which(errork.Theta == max(errork.Theta))[1]
      # numk = which(errork == min(errork))[1]
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
    if(K_hat == 2){
      m12=member[L0==1]
      if(sum(m12==1)<sum(m12==2)){
        n.L.hat = c(which(L.hat==1),which(L.hat==2))
        L.hat[n.L.hat]=c(L0[L0==2],L0[L0==1])
      }
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

Esti.error.JGL.all <- function(whole.data.X, whole.data.Z, result.XZ)
  {
  result.X = result.XZ$result.X
  result.Z = result.XZ$result.Z
  
  # data.X
  result = result.X
  whole.data = whole.data.X
  
  L0=whole.data$L0
  Mu0=whole.data$Mu0
  Theta0=whole.data$Theta0
  data=whole.data$data
  n_all=whole.data$n_all
  K0=whole.data$K0
  
  lambda = result$lambda
  L = length(lambda)
  nameindex = c("lambda","K","CE","CME","PME","TPR","FPR")
  index = as.data.frame(matrix(0,L,length(nameindex)))
  names(index) = nameindex
  
  Mu_hat.list = result$Mu_hat.list
  Theta_hat.list = result$Theta_hat.list
  prob.list = result$prob.list
  member.list = result$member.list
  for (l in 1:L) {
    lam2 = lambda[l]
    index[l,1] = round(lam2,2)
    index[l,-1] = Esti.error.JGL(data, Mu_hat.list[[l]], Theta_hat.list[[l]], 
                                 Mu0, Theta0, K0, L0, member.list[[l]])
  }
  Opt_num = result$Opt_num
  Opt_index = as.data.frame(index[Opt_num,])
  Opt.kmeans = Opt_index
  Opt.kmeans[5:7] = 0
  result.index.X = list(Opt.kmeans=Opt.kmeans, Opt_index=Opt_index,index.all=index)
  
  
  # data.Z
  result = result.Z
  whole.data = whole.data.Z
  
  L0=whole.data$L0
  Mu0=whole.data$Mu0
  Theta0=whole.data$Theta0
  data=whole.data$data
  n_all=whole.data$n_all
  K0=whole.data$K0
  
  lambda = result$lambda
  L = length(lambda)
  nameindex = c("lambda","K","CE","CME","PME","TPR","FPR")
  index = as.data.frame(matrix(0,L,length(nameindex)))
  names(index) = nameindex
  
  Mu_hat.list = result$Mu_hat.list
  Theta_hat.list = result$Theta_hat.list
  prob.list = result$prob.list
  member.list = result$member.list
  for (l in 1:L) {
    lam2 = lambda[l]
    index[l,1] = round(lam2,2)
    index[l,-1] = Esti.error.JGL(data, Mu_hat.list[[l]], Theta_hat.list[[l]], 
                                 Mu0, Theta0, K0, L0, member.list[[l]])
  }
  Opt_num = result$Opt_num
  Opt_index = as.data.frame(index[Opt_num,])
  Opt.kmeans = Opt_index
  Opt.kmeans[5:7] = 0
  result.index.Z = list(Opt.kmeans=Opt.kmeans, Opt_index=Opt_index,index.all=index)
  
  return(result.XZ.all = list(result.index.X = result.index.X, 
                              result.index.Z=result.index.Z,
                              X.kmeans = result.index.X$Opt.kmeans[,-1],
                              X.JGL = result.index.X$Opt_index[,-1],
                              Z.kmeans = result.index.Z$Opt.kmeans[,-1],
                              Z.JGL = result.index.Z$Opt_index[,-1]))
}
