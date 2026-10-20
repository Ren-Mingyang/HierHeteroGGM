##################################################################################
# This document includes main functions of proposed methods, 
# which are used to support numerical simulation studies and real data analysis
##################################################################################



############################# Functions for main algorithms ############################
HHGGM = function(data.X, data.Z, K, lambda1 = 0.5, lambda2 = 0.2, lambda3 = 2, lambda4 = 2, a = 3, rho = 1, 
                eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5, initialization=T, initialize, 
                average=F, asymmetric=T, local_appro=T, cutoff.X=5, cutoff.Z=15, hier.enf=F){

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: HHGGM
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            The key function of Gaussian graphical model-based heterogeneity analysis via penalized fusion:
  ##            identifying K_0 and reconstructing the network structure.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: Update_Theta_Hier()    AMA_XI_Hier()   mcp_d()   S_soft()    MCP_soft()    f.den.vec()    Symmetrize()    cut_diff_ama()
  ##            R packages: NbClust
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ lambda1: a float value, the tuning parameter controlling the sparse of the mean parameter. 
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ niter: int, Maximum number of cycles of the EM algorithm, the default setting is 20.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ initialization: the logical variable, whether to calculate the initial value, the default setting is T,
  ##                                     if initialization = F, the initial value uses initialize.
  ## @ initialize: a given initial value used if initialization = F.
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among 
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list HHGGM_res including:
  ## @ mu: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups (K0_hat is the number of identified subgroups).
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ niter: int, the actual number of cycles of the algorithm.
  ## @ diff_Xi: a float value, the L2-norm difference between subgroups.
  ## @ prob0: K0_hat * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ L.mat0: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## @ Theta0: p * p * K array, the estimated original precision matrices of given K subgroups.
  ## @ Xi0: p * p * K array, the estimated original dual variables of the precision matrices of given K subgroups in ADMM.
  ## @ mu0: K * p matrix, the estimated original mean vectors of K subgroups.
  ## @ group: a list, the sequence number partition of original K subgroups compressed to K0_hat.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ fit.error: a float value, the value of the loss function (without penalty function) corresponding the choice of given tuning parameters.
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  n <- as.integer(dim(data.X)[1])
  q <- as.integer(dim(data.X)[2])
  p <- as.integer(dim(data.Z)[2])
  K_c <- combn(K,2)
  if (initialization) {
    set.seed(1)
    out.initial = initialize_fuc(cbind(data.X,data.Z),K)
    prob = out.initial$prob
    mu = out.initial$Mu
    Theta = out.initial$Theta
    memb = out.initial$memb
    L.mat = matrix(0,n,K)
    for(jj in 1:n) L.mat[jj, memb[jj]]=1
    
    mu.X <- mu[,1:q]; mu.Z <- mu[,q+(1:p)]
    Theta.X <- Theta[1:q,1:q,]; Theta.Z <- Theta[q+(1:p),q+(1:p),]
    L.mat.X <- L.mat; L.mat.Z <- L.mat
    prob.X <- prob; prob.Z <- prob
  } else {
    Theta = initialize$Theta
    mu = initialize$Mu
    prob = initialize$prob
    L.mat = initialize$L.mat
    memb = apply(L.mat,1,which.max)
    
    mu.X <- mu[,1:q]; mu.Z <- mu[,q+(1:p)]
    Theta.X <- Theta[1:q,1:q,]; Theta.Z <- Theta[q+(1:p),q+(1:p),]
    L.mat.X <- L.mat; L.mat.Z <- L.mat
    prob.X <- prob; prob.Z <- prob
  }
  
  # EM algorithm 
  f.mat.X = matrix(0,n,K)
  L.mat.X.old = matrix(0,n,K)
  mu.X.old = matrix(0, K, q)
  Theta.X.old = array(0, dim = c(q, q, K))
  f.mat.Z = matrix(0,n,K)
  L.mat.Z.old = matrix(0,n,K)
  mu.Z.old = matrix(0, K, p)
  Theta.Z.old = array(0, dim = c(p, p, K))
  
  
  t = 0
  diff_mu = 10
  diff_theta = 10
  while( diff_theta >= eps  && t < niter )
  {  
    prob.X.old = prob.X
    prob.Z.old = prob.Z
    mu.X.old = mu.X
    mu.Z.old = mu.Z
    Theta.X.old = Theta.X
    Theta.Z.old = Theta.Z
    L.mat.X.old = L.mat.X
    L.mat.Z.old = L.mat.Z
    
    # calculate pdfs
    for(k.ind in 1:K) {
      f.mat.X[,k.ind]=f.den.vec( data.X, as.numeric(mu.X.old[k.ind,]), Theta.X.old[,,k.ind] ) 
      f.mat.Z[,k.ind]=f.den.vec( data.Z, as.numeric(mu.Z.old[k.ind,]), Theta.Z.old[,,k.ind] ) 
    }
    
    # update L and pi
    for(k.ind in 1:K) {
      for(i in 1:n) {
        L.mat.X[i,k.ind] = prob.X.old[k.ind] * f.mat.X[i,k.ind] / prob.X.old %*% f.mat.X[i,]
        L.mat.Z[i,k.ind] = prob.Z.old[k.ind] * f.mat.Z[i,k.ind] / prob.Z.old %*% f.mat.Z[i,]
      }
      prob.X[k.ind] = mean(L.mat.X[,k.ind])
      prob.Z[k.ind] = mean(L.mat.Z[,k.ind])
    }
    nK.X = apply(L.mat.X,2,sum) + 0.00001
    nK.Z = apply(L.mat.Z,2,sum) + 0.00001
    
    
    Theta_kk_diff.X <- rep(0,dim(K_c)[2])
    Theta_kk_diff.Z <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {
      Theta_kk_diff.X[l] <- sum((Theta.X.old[,,K_c[1,l]] - Theta.X.old[,,K_c[2,l]])^2)
      Theta_kk_diff.Z[l] <- sum((Theta.Z.old[,,K_c[1,l]] - Theta.Z.old[,,K_c[2,l]])^2)
    }
    
    # update mean vectors                       
    for(k.ind in 1:K) {  
      for(j in 1:q){   
        tmp.X = t(t(data.X) - as.numeric(mu.X[k.ind,])) %*% Theta.X.old[,j,k.ind] + mu.X[k.ind,j] * Theta.X.old[j,j,k.ind]
        hj.X = t(L.mat.X[,k.ind])%*%tmp.X
        tau_k.X = sqrt(apply((t(mu.X[k.ind,] - t(mu.X[-k.ind,])))^2,1,sum) + Theta_kk_diff.X[ceiling(which(K_c == k.ind)/2)])
        tau_k.XZ = sqrt(apply((t(mu.X[k.ind,] - t(mu.X[-k.ind,])))^2,1,sum) + Theta_kk_diff.X[ceiling(which(K_c == k.ind)/2)]
                        + apply((t(mu.Z[k.ind,] - t(mu.Z[-k.ind,])))^2,1,sum) + Theta_kk_diff.Z[ceiling(which(K_c == k.ind)/2)])
        v_k_tilde.X = sum( (mcp_d(tau_k.X, lambda4, a) / (tau_k.X+0.00001) + mcp_d(tau_k.XZ, lambda3, a) / (tau_k.XZ+0.00001) ) * mu.X[-k.ind,j])
        v_k_hat.X = sum( mcp_d(tau_k.X, lambda4, a) / (tau_k.X+0.00001) + mcp_d(tau_k.XZ, lambda3, a) / (tau_k.XZ+0.00001) )
        mcp_lambda1.X = mcp_d(mu.X[k.ind,j], lambda1, a)
        if(local_appro){
          mu.X[k.ind,j] = (hj.X + n*v_k_tilde.X) / (n*v_k_hat.X + nK.X[k.ind] * Theta.X.old[j,j,k.ind] + n*mcp_lambda1.X/(abs(mu.X[k.ind,j])+0.00001))
        } 
      }
      for(j in 1:p){   
        tmp.Z = t(t(data.Z) - as.numeric(mu.Z[k.ind,])) %*% Theta.Z.old[,j,k.ind] + mu.Z[k.ind,j] * Theta.Z.old[j,j,k.ind]
        hj.Z = t(L.mat.Z[,k.ind])%*%tmp.Z
        tau_k.XZ = sqrt(apply((t(mu.X[k.ind,] - t(mu.X[-k.ind,])))^2,1,sum) + Theta_kk_diff.X[ceiling(which(K_c == k.ind)/2)]
                        + apply((t(mu.Z[k.ind,] - t(mu.Z[-k.ind,])))^2,1,sum) + Theta_kk_diff.Z[ceiling(which(K_c == k.ind)/2)])
        v_k_tilde.Z = sum(mcp_d(tau_k.XZ, lambda3, a) / (tau_k.XZ+0.00001) * mu.Z[-k.ind,j])
        v_k_hat.Z = sum( mcp_d(tau_k.XZ, lambda3, a) / (tau_k.XZ+0.00001) )
        mcp_lambda1.Z = mcp_d(mu.Z[k.ind,j], lambda1, a)
        if(local_appro){
          mu.Z[k.ind,j] = (hj.Z + n*v_k_tilde.Z) / (n*v_k_hat.Z + nK.Z[k.ind] * Theta.Z.old[j,j,k.ind] + n*mcp_lambda1.Z/(abs(mu.Z[k.ind,j])+0.00001))
        } 

      }
      
    }
    mu.X[abs(mu.X) < 1e-3] <- 0
    mu.Z[abs(mu.Z) < 1e-3] <- 0
    

    # update precision matrices
    mu_kk_diff.X <- rep(0,dim(K_c)[2])
    mu_kk_diff.Z <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {
      mu_kk_diff.X[l] <- sum((mu.X[K_c[1,l],] - mu.X[K_c[2,l],])^2)
      mu_kk_diff.Z[l] <- sum((mu.Z[K_c[1,l],] - mu.Z[K_c[2,l],])^2)
    }

    # define the pseudo sample covariance matrices \tilde{S}
    S.X = array(0, dim = c(q, q, K))
    S.Z = array(0, dim = c(p, p, K))
    for (k.ind in 1:K) {
      L_ikx = sqrt(L.mat.X[,k.ind])*t(t(data.X) - mu.X[k.ind,])
      S.X[,,k.ind] = t(L_ikx) %*% L_ikx / nK.X[k.ind]
      L_ikz = sqrt(L.mat.Z[,k.ind])*t(t(data.Z) - mu.Z[k.ind,])
      S.Z[,,k.ind] = t(L_ikz) %*% L_ikz / nK.Z[k.ind]
    }
    
    Theta_out = Update_Theta_Hier(S.X,S.Z,nK.X,nK.Z,lambda2,lambda3,lambda4, mu_kk_diff.X,mu_kk_diff.Z, K_c, tol=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA)
    Theta.X = Theta_out$Theta.X
    Xi.X = Theta_out$Xi.X
    Theta.Z = Theta_out$Theta.Z
    Xi.Z = Theta_out$Xi.Z
    V_kk.X = Theta_out$V_kk.X
    V_kk.Z = Theta_out$V_kk.Z
    
    t = t + 1
    diff_mu = norm(mu.X.old-mu.X,type="2")/(norm(mu.X,type="2")+0.001) + norm(mu.Z.old-mu.Z,type="2")/(norm(mu.Z,type="2")+0.001)
    diff_theta = norm(Theta.X.old-Theta.X,type="2")/(norm(Theta.X,type="2")+0.001) + norm(Theta.Z.old-Theta.Z,type="2")/(norm(Theta.Z,type="2")+0.001)
  }
  
  group_final.X=cut_diff_ama(mu_kk_diff.X ,K_c,K,cutoff=cutoff.X)
  group_final.Z=cut_diff_ama(mu_kk_diff.Z ,K_c,K,cutoff=cutoff.Z)
  if(hier.enf){
    group_final.Z0 = list()
    group_final.XZ = list()
    K_0.X = length(group_final.X)
    K_0.Z = length(group_final.Z)
    v=1
    for (k1 in 1:K_0.X) {
      gX.K = group_final.X[[k1]]
      vecZ = NULL
      vv=1
      for (k2 in 1:K_0.Z) {
        gZ.K = group_final.Z[[k2]]
        gXZ.K = gZ.K[!is.na(match(gZ.K,gX.K))]
        if(length(gXZ.K)>0){
          group_final.Z0[[v]] = gXZ.K
          vecZ[vv] = v
          vv=vv+1
          v=v+1
        }
      }
      group_final.XZ[[k1]] = vecZ
    }
    group_final.Z = group_final.Z0
  }else{ group_final.XZ = NULL }
  group_final.X
  group_final.Z
  group_final.XZ
  mu.X
  mu.Z
  
  K_0.X = length(group_final.X)
  mu_final.X = matrix(0, K_0.X, q)
  Theta_final.X = array(0, dim = c(q, q, K_0.X))
  Xi_final.X = array(0, dim = c(q, q, K_0.X))
  prob0.X = rep(0,K_0.X)
  L.mat0.X = matrix(0,n,K_0.X)
  
  K_0.Z = length(group_final.Z)
  mu_final.Z = matrix(0, K_0.Z, p)
  Theta_final.Z = array(0, dim = c(p, p, K_0.Z))
  Xi_final.Z = array(0, dim = c(p, p, K_0.Z))
  prob0.Z = rep(0,K_0.Z)
  L.mat0.Z = matrix(0,n,K_0.Z)
  
  for (l in 1:K_0.X) {
    gg.X = group_final.X[[l]]
    prob0.X[l] = sum(prob.X[gg.X])
    L.mat0.X[,l] = apply(as.matrix(L.mat.X[,gg.X]),1,sum)
    if(length(gg.X) > 1){
      mu_final.X[l,] = apply(mu.X[gg.X,],2,mean)
      
      if(!average){
        Theta_final.X[,,l] = Theta.X[,,gg.X[which.max(nK.X[gg.X])]] 
        Xi_final.X[,,l] = Xi.X[,,gg.X[which.max(nK.X[gg.X])]]       
      } else {
        Theta_final.X[,,l] = Theta.X[,,gg.X[1]]/length(gg.X)
        Xi_final.X[,,l] = Xi.X[,,gg.X[1]]/length(gg.X)
        for (gi in gg.X[-1]) {
          Theta_final.X[,,l] = Theta_final.X[,,l] + Theta.X[,,gi]/length(gg.X)
          Xi_final.X[,,l] = Xi_final.X[,,l] + Xi.X[,,gi]/length(gg.X)
        }
      }
      
    }else{ mu_final.X[l,] = mu.X[gg.X,]; Theta_final.X[,,l] = Theta.X[,,gg.X]; Xi_final.X[,,l] = Xi.X[,,gg.X] }
  }
  
  if(asymmetric){
    for(k in 1:K_0.X) {
      Theta_final.X[,,k] = Symmetrize(Theta_final.X[,,k])
      Xi_final.X[,,k] = Symmetrize(Xi_final.X[,,k])
    }
  }
  
  for (l in 1:K_0.Z) {
    gg.Z = group_final.Z[[l]]
    prob0.Z[l] = sum(prob.Z[gg.Z])
    L.mat0.Z[,l] = apply(as.matrix(L.mat.Z[,gg.Z]),1,sum)
    if(length(gg.Z) > 1){
      mu_final.Z[l,] = apply(mu.Z[gg.Z,],2,mean)
      
      if(!average){
        Theta_final.Z[,,l] = Theta.Z[,,gg.Z[which.max(nK.Z[gg.Z])]] 
        Xi_final.Z[,,l] = Xi.Z[,,gg.Z[which.max(nK.Z[gg.Z])]]       
      } else {
        Theta_final.Z[,,l] = Theta.Z[,,gg.Z[1]]/length(gg.Z)
        Xi_final.Z[,,l] = Xi.Z[,,gg.Z[1]]/length(gg.Z)
        for (gi in gg.Z[-1]) {
          Theta_final.Z[,,l] = Theta_final.Z[,,l] + Theta.Z[,,gi]/length(gg.Z)
          Xi_final.Z[,,l] = Xi_final.Z[,,l] + Xi.Z[,,gi]/length(gg.Z)
        }
      }
      
    }else{ mu_final.Z[l,] = mu.Z[gg.Z,]; Theta_final.Z[,,l] = Theta.Z[,,gg.Z]; Xi_final.Z[,,l] = Xi.Z[,,gg.Z]}
  }
  
  if(asymmetric){
    for(k in 1:K_0.Z) {
      Theta_final.Z[,,k] = Symmetrize(Theta_final.Z[,,k])
      Xi_final.Z[,,k] = Symmetrize(Xi_final.Z[,,k])
    }
  }
  
  Theta_final.X[abs(Theta_final.X) < 1e-3] <- 0
  member.X = apply(L.mat0.X,1,which.max)
  member.Z = apply(L.mat0.Z,1,which.max)
  BIC.res.X = BIC(data.X, mu_final.X, Xi_final.X, L.mat0.X)
  BIC.res.Z = BIC(data.Z, mu_final.Z, Xi_final.Z, L.mat0.Z)
  
  if(hier.enf){
    memb.XZ = hier.enf.HHGGM(data.X, data.Z, group_final.XZ, mu_final.X, mu_final.Z, 
                             Theta_final.X, Theta_final.Z, prob0.X, prob0.Z)
    
    member.X = memb.XZ$member.X
    member.Z = memb.XZ$member.Z
    L.mat0.X = memb.XZ$L.mat0.X
    L.mat0.Z = memb.XZ$L.mat0.Z
    prob0.X = memb.XZ$prob0.X
    prob0.X = memb.XZ$prob0.X
  }
  
  
  HHGGM_res <- list();
  HHGGM_res$niter <- t; 
  HHGGM_res$mu.X <-  mu_final.X; HHGGM_res$Theta.X <- Theta_final.X
  HHGGM_res$Xi.X <- Xi_final.X; HHGGM_res$diff_Xi.X <- V_kk.X
  HHGGM_res$prob0.X <- prob0.X; HHGGM_res$L.mat0.X <- L.mat0.X;
  HHGGM_res$Theta0.X <- Theta.X; HHGGM_res$member.X <- member.X;
  HHGGM_res$Xi0.X <- Xi.X; HHGGM_res$mu0.X <- mu.X;  HHGGM_res$group.X <- group_final.X
  HHGGM_res$bic.X <- BIC.res.X$bic;  HHGGM_res$fit.error.X <- BIC.res.X$fit.error;  HHGGM_res$df.X <- BIC.res.X$df
  
  HHGGM_res$mu.Z <-  mu_final.Z; HHGGM_res$Theta.Z <- Theta_final.Z
  HHGGM_res$Xi.Z <- Xi_final.Z; HHGGM_res$diff_Xi.Z <- V_kk.Z
  HHGGM_res$prob0.Z <- prob0.Z; HHGGM_res$L.mat0.Z <- L.mat0.Z;
  HHGGM_res$Theta0.Z <- Theta.Z; HHGGM_res$member.Z <- member.Z;
  HHGGM_res$Xi0.Z <- Xi.Z; HHGGM_res$mu0.Z <- mu.Z;  HHGGM_res$group.Z <- group_final.Z
  HHGGM_res$bic.Z <- BIC.res.Z$bic;  HHGGM_res$fit.error.Z <- BIC.res.Z$fit.error;  HHGGM_res$df.Z <- BIC.res.Z$df
  HHGGM_res$group <- group_final.XZ
  return(HHGGM_res)
}


Update_Theta_Hier = function(S.X, S.Z, nK.X, nK.Z, lambda2, lambda3, lambda4, 
                        mu_kk_diff.X, mu_kk_diff.Z, K_c, a = 3, rho=1, 
                        maxiter=10, maxiter.AMA=5, tol=1e-2, rho.increment=1){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: Update_Theta_Hier
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Updating the precision matrices using the ADMM algorithm after updating the mean vectors in the EM algorithm. 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: AMA_XI_Hier()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ S: p * p * K array, the estimated pseudo sample covariance matrices of given K subgroups in the iteration of the EM algorithm.
  ## @ nK: K * 1 vector, the estimated sample sizes of K subgroups in the iteration of the EM algorithm.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ mu_kk_diff: L2-norm differences in the mean vector between different subgroups.
  ## @ K_c: All combinations of natural numbers within K.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ tol: a float value, algorithm termination threshold.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list Theta_out including:
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ V_kk: L2-norm differences in the precision matrices between different subgroups.
  ## @ iters: int, the actual number of cycles of the algorithm.
  ## @ diff: a float value, the L1-norm difference between Theta and Xi.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  q = dim(S.X)[1]
  p = dim(S.Z)[1]
  K = dim(S.X)[3]
  nkk.X = rep(1,K)
  nkk.Z = rep(1,K)
  # nkk.X = nK.X
  # nkk.Z = nK.Z
  # initialize Theta:
  Theta.X = array(0, dim = c(q, q, K))
  Theta.Z = array(0, dim = c(p, p, K))
  for (k in 1:K) { Theta.X[,,k] = diag(1,q); Theta.Z[,,k] = diag(1,p)}
  # initialize Xi:
  Xi.X = array(0, dim = c(q, q, K))
  Xi.Z = array(0, dim = c(p, p, K))
  # initialize Phi:
  Phi.X = array(0, dim = c(q, q, K))
  Phi.Z = array(0, dim = c(p, p, K))
  
  iter=0
  diff_value = 10
  diff_value_Xi = 10
  while( iter<maxiter && !(diff_value < tol || diff_value_Xi < tol) )
  {
    Theta.prev.X = Theta.X
    Theta.prev.Z = Theta.Z
    Xi.prev.X = Xi.X
    Xi.prev.Z = Xi.Z
    Theta.X = Theta.prev.X
    Theta.Z = Theta.prev.Z
    # update Theta:
    for(k.ind in 1:K){
      Sa.X=S.X[,,k.ind] - rho*Xi.X[,,k.ind]/nkk.X[k.ind] + rho*Phi.X[,,k.ind]/nkk.X[k.ind]
      edecomp.X = eigen(Sa.X)
      D.X = edecomp.X$values
      if(is.complex(D.X)){Sa.X=Symmetrize(Sa.X);edecomp.X = eigen(Sa.X);D.X = edecomp.X$values}
      V.X = edecomp.X$vectors
      D2.X = nkk.X[k.ind]/(2*rho) * ( -D.X + sqrt(D.X^2 + 4*(rho/nkk.X[k.ind]) ) )
      Theta.X[,,k.ind] = V.X %*% diag(D2.X) %*% t(V.X)
      
      Sa.Z=S.Z[,,k.ind] - rho*Xi.Z[,,k.ind]/nkk.Z[k.ind] + rho*Phi.Z[,,k.ind]/nkk.Z[k.ind]
      edecomp.Z = eigen(Sa.Z)
      D.Z = edecomp.Z$values
      if(is.complex(D.Z)){Sa.Z=Symmetrize(Sa.Z);edecomp.Z = eigen(Sa.Z);D.Z = edecomp.Z$values}
      V.Z = edecomp.Z$vectors
      D2.Z = nkk.Z[k.ind]/(2*rho) * ( -D.Z + sqrt(D.Z^2 + 4*(rho/nkk.Z[k.ind]) ) )
      Theta.Z[,,k.ind] = V.Z %*% diag(D2.Z) %*% t(V.Z)
    }
    # update Xi:
    # define M matrices:
    M.X = array(0, dim = c(q, q, K))
    M.Z = array(0, dim = c(p, p, K))
    for(k in 1:K){ M.X[,,k] = Theta.X[,,k] + Phi.X[,,k]; M.Z[,,k] = Theta.Z[,,k] + Phi.Z[,,k] }
    
    Xi_out_list = AMA_XI_Hier(M.X,M.Z,K_c,lambda2,lambda3,lambda4, 
                         mu_kk_diff.X, mu_kk_diff.Z,maxiter=maxiter.AMA)
    Xi.X = Xi_out_list$Xi.X
    Xi.X[abs(Xi.X) < 1e-3] <- 0
    Xi.Z = Xi_out_list$Xi.Z
    Xi.Z[abs(Xi.Z) < 1e-3] <- 0
    V.X = Xi_out_list$V.X
    V.Z = Xi_out_list$V.Z
    V_kk.X = round(apply(V.X^2,3,sum),4)
    V_kk.Z = round(apply(V.Z^2,3,sum),4)
    
    # update the dual variable Phi:
    Phi.X = Phi.X + Theta.X - Xi.X
    Phi.Z = Phi.Z + Theta.Z - Xi.Z
    iter = iter+1
    diff_value = sum(abs(Theta.X - Theta.prev.X)) / sum(abs(Theta.prev.X)) + sum(abs(Theta.Z - Theta.prev.Z)) / sum(abs(Theta.prev.Z))
    diff_value_Xi = sum(abs(Xi.X - Xi.prev.X)) / (sum(abs(Xi.prev.X))+0.001) + sum(abs(Xi.Z - Xi.prev.Z)) / (sum(abs(Xi.prev.Z))+0.001)
    
    # increment rho by a constant factor:
    rho = rho*rho.increment
  }
  
  
  diff = sum(abs(Theta.X-Xi.X)) + sum(abs(Theta.Z-Xi.Z))
  Theta_out = list(Theta.X=Theta.X,Xi.X=Xi.X,Theta.Z=Theta.Z,Xi.Z=Xi.Z,
                   V_kk.X=V_kk.X,V_kk.Z=V_kk.Z,diff=diff,iters=iter)
  return(Theta_out)
}


AMA_XI_Hier = function(M.X, M.Z, K_c, lambda2, lambda3, lambda4, 
                  mu_kk_diff.X, mu_kk_diff.Z, a = 3, 
                  kappa=1, maxiter=5, tol=1e-2){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: AMA_XI_Hier
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Solving (A.11) in Supplementary Materials using S-AMA algorithm, which is the key step of updating the precision matrices. 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: S_soft()   MCP_soft()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ B: The output of the previous step in ADMM algorithm (see (A.11) in Supplementary Materials).
  ## @ K_c: All combinations of natural numbers within K.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ mu_kk_diff: L2-norm differences in the mean vector between different subgroups.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ kappa: a float value, the penalty parameter in S-AMA algorithm, the default setting is 1.
  ## @ tol: a float value, algorithm termination threshold.
  ## @ maxiter: int, Maximum number of cycles of the algorithm.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list Xi_out including:
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ V: p * p * K0_hat array, the estimated dual variables of Xis in S-AMA.
  ## @ Delta: Lagrange multiplier in in S-AMA.
  ## @ iters: int, the actual number of cycles of the algorithm.
  ## @ diff: a float value, the convergence criteria value.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  q = dim(M.X)[1]
  p = dim(M.Z)[1]
  K = dim(M.X)[3]
  num_Kc = dim(K_c)[2]
  # initialize Xi:
  Xi.X = array(0, dim = c(q, q, K))
  Xi.Z = array(0, dim = c(p, p, K))
  for (k in 1:K) { Xi.X[,,k] =  diag(1,q); Xi.Z[,,k] =  diag(1,p)}
  diag_index.X = Xi.X
  diag_index.Z = Xi.Z
  
  # initialize V:
  V.X = array(0, dim = c(q, q, num_Kc))
  V.Z = array(0, dim = c(p, p, num_Kc))
  # initialize Delta:
  Delta.X = array(0, dim = c(q, q, num_Kc))
  Delta.Z = array(0, dim = c(p, p, num_Kc))
  
  U.X = array(0, dim = c(q, q, K))
  U.Z = array(0, dim = c(p, p, K))
  Omega.X = array(0, dim = c(q, q, num_Kc))
  Omega.Z = array(0, dim = c(p, p, num_Kc))
  e_k12 = matrix(0,K,num_Kc)
  for (i in 1:K) {e_k12[i,which(K_c[1,] == i)] = 1;e_k12[i,which(K_c[2,] == i)] = -1}
  
  # iterations
  iter=0
  diff_value = 10
  # lambda4 = lambda4/10
  while( iter<maxiter && diff_value > tol )
  {
    V.prev.X = V.X
    Xi.prev.X = Xi.X
    Delta.prev.X = Delta.X
    V.prev.Z = V.Z
    Xi.prev.Z = Xi.Z
    Delta.prev.Z = Delta.Z
    
    for (i in 1:q) {
      for (j in 1:q) {
        U.X[i,j,] = M.X[i,j,] + apply(t(Delta.X[i,j,] * t(e_k12)),1,sum)
      }
    }
    Xi.X = S_soft(U.X,lambda2)/(1-1/a) * (abs(U.X) <= a*lambda2 & diag_index.X==0) + U.X * (abs(U.X) > a*lambda2 | diag_index.X==1)
    
    for (i in 1:p) {
      for (j in 1:p) {
        U.Z[i,j,] = M.Z[i,j,] + apply(t(Delta.Z[i,j,] * t(e_k12)),1,sum)
      }
    }
    Xi.Z = S_soft(U.Z,lambda2)/(1-1/a) * (abs(U.Z) <= a*lambda2 & diag_index.Z==0) + U.Z * (abs(U.Z) > a*lambda2 | diag_index.Z==1)
    
    for (l in 1:num_Kc) {
      Omega.X[,,l] = Xi.X[,,K_c[1,l]] - Xi.X[,,K_c[2,l]] - Delta.X[,,l]/kappa
      Omega.Z[,,l] = Xi.Z[,,K_c[1,l]] - Xi.Z[,,K_c[2,l]] - Delta.Z[,,l]/kappa
      Omega.r1 = sqrt(mu_kk_diff.X[l] + mu_kk_diff.Z[l] + sum(Omega.X[,,l]^2) + sum(Omega.Z[,,l]^2))
      Omega.r2 = sqrt(mu_kk_diff.X[l] + mu_kk_diff.Z[l] + sum(Omega.X[,,l]^2))
      soft1 = MCP_soft0(Omega.r1,lambda3,kappa=kappa)
      soft2 = MCP_soft0(Omega.r2,lambda4,kappa=kappa)
      if(Omega.r1 > a*lambda3 & Omega.r2 > a*lambda4){
        V.Z[,,l] = Omega.Z[,,l]
        V.X[,,l] = Omega.X[,,l]
      } else if(Omega.r1 <= a*lambda3 & soft1*Omega.r2 > a*lambda4){
        V.Z[,,l] = soft1*Omega.Z[,,l]
        V.X[,,l] = soft1*Omega.X[,,l]
      } else if(sum(Omega.Z[,,l]^2) + (soft2*Omega.r2)^2 > (a*lambda3)^2 & Omega.r1 <= a*lambda4){
        V.Z[,,l] = Omega.Z[,,l]
        V.X[,,l] = soft2*Omega.X[,,l]
      } else{
        V.Z[,,l] = Omega.Z[,,l]/(1+mcp_d(Omega.r1,lambda3)/kappa/(Omega.r1+10^(-7)))
        V.X[,,l] = Omega.X[,,l]/(1+mcp_d(Omega.r1,lambda3)/kappa/(Omega.r1+10^(-7))+mcp_d(Omega.r2,lambda4)/kappa/(Omega.r2+10^(-7)))
      }
    
      Delta.X[,,l] = Delta.X[,,l] + kappa * ( V.X[,,l] - Xi.X[,,K_c[1,l]] + Xi.X[,,K_c[2,l]] )* (sum(V.X[,,l]^2) > 0 )
      Delta.Z[,,l] = Delta.Z[,,l] + kappa * ( V.Z[,,l] - Xi.Z[,,K_c[1,l]] + Xi.Z[,,K_c[2,l]] )* (sum(V.Z[,,l]^2) > 0 )

    }
    
    iter = iter+1
    diff_value = sum(abs(Xi.X - Xi.prev.X)^2) / sum(abs(Xi.prev.X)^2) + sum(abs(Xi.Z - Xi.prev.Z)^2) / sum(abs(Xi.prev.Z)^2)
  }
  Xi_out = list(Xi.X=Xi.X,V.X=V.X,Delta.X=Delta.X,Xi.Z=Xi.Z,V.Z=V.Z,Delta.Z=Delta.Z,diff=diff_value,iters=iter)
  return(Xi_out)
}



############################## Functions for tuning lambdas ############################
BIC = function(data, mu_hat, Theta_hat, L.mat, adjust.bic=T){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: BIC
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the adaptive BIC-type criterion.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: f.den.vec()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu_hat: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups.
  ## @ Theta_hat: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ L.mat: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list P including:
  ## @ fit.error: a float value, the value of the loss function (without penalty function).
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  n = nrow(data)
  K = nrow(mu_hat)
  
  # fitting error
  pi_vec = apply(L.mat, 2, sum)/n
  fit.error_mat = matrix(0, n, K)
  for(k in 1:K) {
    fit.error_mat[,k] = pi_vec[k] * f.den.vec( data, as.numeric(mu_hat[k,]), Theta_hat[,,k] )                  
  }
  fit0 = apply(fit.error_mat, 1, sum)
  fit.error = sum(log( fit0 + min(fit0[fit0>0]) ))
  fit.error = - 2*fit.error
  
  # degrees of freedom
  for(i in 1:K){
    Theta_hat[upper.tri(Theta_hat[, , i], diag = T)] = 0
  }
  
  df =  log(n) * length(which(mu_hat != 0)) + 2 * length(which(Theta_hat != 0))
  if(adjust.bic){
    df =  log(n) * length(which(mu_hat != 0)) + log(n) * length(which(Theta_hat != 0))
  }
  bic = fit.error + df
  P = list()
  P$fit.error = fit.error
  P$df = df
  P$bic = bic
  return(P)
}

genelambda.obo4 = function(nlambda1=10,lambda1_max=1,lambda1_min=0.05,
                           nlambda2=10,lambda2_max=1,lambda2_min=0.01,
                           nlambda3=10,lambda3_max=5,lambda3_min=0.5,
                           nlambda4=5,lambda4_max=1,lambda4_min=0.1){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: genelambda.obo
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating a sequence of the tuning parameters (lambda1, lambda2, and lambda3).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ nlambda1, nlambda2, and nlambda3: The numbers of lambda 1 2 3.
  ## @ lambda1_min, lambda2_min, and lambda3_min: The minimum values of lambda 1 2 3.
  ## @ lambda1_max, lambda2_max, and lambda3_max: The maximum values of lambda 1 2 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ lambda: a sequence of the tuning parameters (lambda1, lambda2, and lambda3).   
  ## -----------------------------------------------------------------------------------------------------------------
  
  lambda1 = exp(seq(log(lambda1_max),log(lambda1_min),len= nlambda1))
  lambda2 =exp(seq(log(lambda2_max),log(lambda2_min),len= nlambda2))
  lambda3 =exp(seq(log(lambda3_max),log(lambda3_min),len= nlambda3))
  lambda4 =exp(seq(log(lambda4_max),log(lambda4_min),len= nlambda4))
  lambda = list(lambda1=lambda1,lambda2=lambda2,lambda3=lambda3,lambda4=lambda4)
  return(lambda)
}

tuning.lambda.HHGGM = function(data.X, data.Z, lambda, K, initial.selection="K-means", 
                              initialize, average=F, asymmetric=T, eps = 5e-2, maxiter=10, maxiter.AMA=5,
                              local_appro=T, trace = T, cutoff.X=5, cutoff.Z=20, 
                              hier.enf=F, merge.max=F, hier.enf.refitting=F){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: tuning.lambda.HHGGM
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Searching and selecting the optional tuning parameters under the adaptive BIC-type criterion
  ##            using the proposed method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: HHGGM.refit()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ lambda: a list, the sequences of the tuning parameters (lambda1, lambda2, and lambda3).
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ initial.selection: the different initial values from two clustering methods, which can be selected from c("K-means","dbscan").
  ## @ initialize: A given initial values, which should be given when initial.selection is not in c("K-means","dbscan").
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among 
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## @ trace: the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "result" including:
  ## @ Opt_lambda: the selected optional tuning parameters.
  ## @ Mu_hat.list: the estimated mean vectors of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ Theta_hat.list: the estimated precision matrices of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ prob.list: the estimated mixture probabilities of subgroups corresponding all choices of given tuning parameters.
  ## @ member.list: subgroup labels to which each sample belongs corresponding all choices of given tuning parameters.
  ## @ L.mat.list: the estimated probability that each sample belongs to each subgroup corresponding all choices of given tuning parameters.
  ## @ Opt_aBIC: the optional BIC value.
  ## @ Opt_num: the position of the optimal parameter for all given parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  q <- as.integer(dim(data.X)[2])
  p <- as.integer(dim(data.Z)[2])
  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  lambda3 = lambda$lambda3
  lambda4 = lambda$lambda4
  L1 = length(lambda1)
  L2 = length(lambda2)
  L3 = length(lambda3)
  L4 = length(lambda4)
  L = L1+L2+L3+L4

  aBIC.X = rep(0,L)
  aBIC.Z = rep(0,L)
  fit.error.X.v = rep(0,L)
  df.X.v = rep(0,L)
  fit.error.Z.v = rep(0,L)
  df.Z.v = rep(0,L)
  n_all = dim(data.X)[1]
  # initialize
  if(initial.selection=="K-means"){
    set.seed(1)
    out.initial = initialize_fuc(cbind(data.X,data.Z),K)
    memb = out.initial$memb
    L.mat = matrix(0,n_all,K)
    for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
    out.initial$L.mat = L.mat
  } else if(initial.selection=="dbscan"){
    out.initial = initialize_fuc.dbscan(cbind(data.X,data.Z),K)
    memb = out.initial$memb
    L.mat = matrix(0,n_all,K)
    for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
    out.initial$L.mat = L.mat
  } else {out.initial = initialize}
  if(hier.enf.refitting){cutoff.X=1}
  if(L == 4){
    l=1
    aBIC.X = rep(0,l)
    aBIC.Z = rep(0,l)
    fit.error.X.v = rep(0,l)
    df.X.v = rep(0,l)
    fit.error.Z.v = rep(0,l)
    df.Z.v = rep(0,l)
    Mu_hat.list.X = list()
    Theta_hat.list.X = list()
    prob.list.X = list()
    L.mat.list.X = list()
    member.list.X = list()
    
    Mu_hat.list.Z = list()
    Theta_hat.list.Z = list()
    prob.list.Z = list()
    L.mat.list.Z = list()
    member.list.Z = list()
    
    group.list = list()
    
    lam1 = lambda1;lam2 = lambda2;lam3 = lambda3;lam4 = lambda4;
    PP = HHGGM(data.X, data.Z, K, lam1, lam2, lam3, lam4, initialization=F, initialize=out.initial, average=average, 
               asymmetric=asymmetric, local_appro=local_appro, cutoff.X=cutoff.X, cutoff.Z=cutoff.Z, hier.enf=hier.enf)
    
    mu_hat.X=PP$mu.X;Theta_hat.X=PP$Xi.X;L.mat.X = PP$L.mat0.X;group.X = PP$group.X;prob.X = PP$prob0.X;aBIC.X[l] = PP$bic.X; member.X = PP$member.X
    Mu_hat.list.X[[l]]=mu_hat.X; Theta_hat.list.X[[l]]=Theta_hat.X; prob.list.X[[l]]=prob.X; L.mat.list.X[[l]] = L.mat.X; member.list.X[[l]]=member.X
    mu_hat.Z=PP$mu.Z;Theta_hat.Z=PP$Xi.Z;L.mat.Z = PP$L.mat0.Z;group.Z = PP$group.Z;prob.Z = PP$prob0.Z;aBIC.Z[l] = PP$bic.Z; member.Z = PP$member.Z
    Mu_hat.list.Z[[l]]=mu_hat.Z; Theta_hat.list.Z[[l]]=Theta_hat.Z; prob.list.Z[[l]]=prob.Z; L.mat.list.Z[[l]] = L.mat.Z; member.list.Z[[l]]=member.Z
    aBIC = aBIC.X + aBIC.Z
    
    fit.error.X.v[l] = PP$fit.error.X;  df.X.v[l] = PP$df.X
    fit.error.Z.v[l] = PP$fit.error.Z;  df.Z.v[l] = PP$df.Z
    group.list[[l]] = PP$group
    
    } else {
    Mu_hat.list.X = list()
    Theta_hat.list.X = list()
    prob.list.X = list()
    L.mat.list.X = list()
    member.list.X = list()
    
    Mu_hat.list.Z = list()
    Theta_hat.list.Z = list()
    prob.list.Z = list()
    L.mat.list.Z = list()
    member.list.Z = list()
    group.list = list()
    # search lam3
    lam1 = median(lambda1);lam2 = median(lambda2);lam4 = 0
    for (l in 1:L3) {
      lam3 = lambda3[l]
      
      PP = HHGGM(data.X, data.Z, K, lam1, lam2, lam3, lam4, initialization=F, initialize=out.initial, average=average, 
                 asymmetric=asymmetric, local_appro=local_appro, cutoff.X=cutoff.X, cutoff.Z=cutoff.Z, hier.enf=hier.enf)
      mu_hat.X=PP$mu.X;Theta_hat.X=PP$Xi.X;L.mat.X = PP$L.mat0.X;group.X = PP$group.X;prob.X = PP$prob0.X;aBIC.X[l] = PP$bic.X; member.X = PP$member.X
      Mu_hat.list.X[[l]]=mu_hat.X; Theta_hat.list.X[[l]]=Theta_hat.X; prob.list.X[[l]]=prob.X; L.mat.list.X[[l]] = L.mat.X; member.list.X[[l]]=member.X
      mu_hat.Z=PP$mu.Z;Theta_hat.Z=PP$Xi.Z;L.mat.Z = PP$L.mat0.Z;group.Z = PP$group.Z;prob.Z = PP$prob0.Z;aBIC.Z[l] = PP$bic.Z; member.Z = PP$member.Z
      Mu_hat.list.Z[[l]]=mu_hat.Z; Theta_hat.list.Z[[l]]=Theta_hat.Z; prob.list.Z[[l]]=prob.Z; L.mat.list.Z[[l]] = L.mat.Z; member.list.Z[[l]]=member.Z
      fit.error.X.v[l] = PP$fit.error.X;  df.X.v[l] = PP$df.X
      fit.error.Z.v[l] = PP$fit.error.Z;  df.Z.v[l] = PP$df.Z
      group.list[[l]] = PP$group
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 lam4 ="),c(round(lam1,2),round(lam2,2),round(lam3,2),round(lam4,2)),":"),
                paste("K1=",as.numeric(dim(Theta_hat.X)[3]), "K2=",as.numeric(dim(Theta_hat.Z)[3]))))
      }
    }
    aBIC.X[is.na(aBIC.X)] = 10^10
    aBIC.X[aBIC.X==0] = abs(aBIC.X[aBIC.X!=0][1])*10
    aBIC.Z[is.na(aBIC.Z)] = 10^10
    aBIC.Z[aBIC.Z==0] = abs(aBIC.Z[aBIC.Z!=0][1])*10
    aBIC = aBIC.X + aBIC.Z
    n_lam3 = which(aBIC[1:L3] == min(aBIC[1:L3]))[1];lam3 = lambda3[n_lam3]
    
    # search lam2
    for (l2 in 1:L2) {
      lam2 = lambda2[l2];l = L3+l2
      
      PP = HHGGM(data.X, data.Z, K, lam1, lam2, lam3, lam4, initialization=F, initialize=out.initial, average=average, 
                 asymmetric=asymmetric, local_appro=local_appro, cutoff.X=cutoff.X, cutoff.Z=cutoff.Z, hier.enf=hier.enf)
      mu_hat.X=PP$mu.X;Theta_hat.X=PP$Xi.X;L.mat.X = PP$L.mat0.X;group.X = PP$group.X;prob.X = PP$prob0.X;aBIC.X[l] = PP$bic.X; member.X = PP$member.X
      Mu_hat.list.X[[l]]=mu_hat.X; Theta_hat.list.X[[l]]=Theta_hat.X; prob.list.X[[l]]=prob.X; L.mat.list.X[[l]] = L.mat.X; member.list.X[[l]]=member.X
      mu_hat.Z=PP$mu.Z;Theta_hat.Z=PP$Xi.Z;L.mat.Z = PP$L.mat0.Z;group.Z = PP$group.Z;prob.Z = PP$prob0.Z;aBIC.Z[l] = PP$bic.Z; member.Z = PP$member.Z
      Mu_hat.list.Z[[l]]=mu_hat.Z; Theta_hat.list.Z[[l]]=Theta_hat.Z; prob.list.Z[[l]]=prob.Z; L.mat.list.Z[[l]] = L.mat.Z; member.list.Z[[l]]=member.Z
      
      fit.error.X.v[l] = PP$fit.error.X;  df.X.v[l] = PP$df.X
      fit.error.Z.v[l] = PP$fit.error.Z;  df.Z.v[l] = PP$df.Z
      group.list[[l]] = PP$group
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 lam4 ="),c(round(lam1,2),round(lam2,2),round(lam3,2),round(lam4,2)),":"),
                paste("K1=",as.numeric(dim(Theta_hat.X)[3]), "K2=",as.numeric(dim(Theta_hat.Z)[3]))))
      }
    }
    aBIC.X[is.na(aBIC.X)] = 10^10
    aBIC.X[aBIC.X==0] = abs(aBIC.X[aBIC.X!=0][1])*10
    aBIC.Z[is.na(aBIC.Z)] = 10^10
    aBIC.Z[aBIC.Z==0] = abs(aBIC.Z[aBIC.Z!=0][1])*10
    aBIC = aBIC.X + aBIC.Z
    n_lam2 = which(aBIC[(L3+1):(L3+L2)] == min(aBIC[(L3+1):(L3+L2)]))[1];lam2 = lambda2[n_lam2]
    # search lam1
    for (l1 in 1:L1) {
      lam1 = lambda1[l1];l = L3+L2+l1
      
      PP = HHGGM(data.X, data.Z, K, lam1, lam2, lam3, lam4, initialization=F, initialize=out.initial, average=average, 
                 asymmetric=asymmetric, local_appro=local_appro, cutoff.X=cutoff.X, cutoff.Z=cutoff.Z, hier.enf=hier.enf)
      mu_hat.X=PP$mu.X;Theta_hat.X=PP$Xi.X;L.mat.X = PP$L.mat0.X;group.X = PP$group.X;prob.X = PP$prob0.X;aBIC.X[l] = PP$bic.X; member.X = PP$member.X
      Mu_hat.list.X[[l]]=mu_hat.X; Theta_hat.list.X[[l]]=Theta_hat.X; prob.list.X[[l]]=prob.X; L.mat.list.X[[l]] = L.mat.X; member.list.X[[l]]=member.X
      mu_hat.Z=PP$mu.Z;Theta_hat.Z=PP$Xi.Z;L.mat.Z = PP$L.mat0.Z;group.Z = PP$group.Z;prob.Z = PP$prob0.Z;aBIC.Z[l] = PP$bic.Z; member.Z = PP$member.Z
      Mu_hat.list.Z[[l]]=mu_hat.Z; Theta_hat.list.Z[[l]]=Theta_hat.Z; prob.list.Z[[l]]=prob.Z; L.mat.list.Z[[l]] = L.mat.Z; member.list.Z[[l]]=member.Z
      
      fit.error.X.v[l] = PP$fit.error.X;  df.X.v[l] = PP$df.X
      fit.error.Z.v[l] = PP$fit.error.Z;  df.Z.v[l] = PP$df.Z
      group.list[[l]] = PP$group
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 lam4 ="),c(round(lam1,2),round(lam2,2),round(lam3,2),round(lam4,2)),":"),
                paste("K1=",as.numeric(dim(Theta_hat.X)[3]), "K2=",as.numeric(dim(Theta_hat.Z)[3]))))
      }
    }
    aBIC.X[is.na(aBIC.X)] = 10^10
    aBIC.X[aBIC.X==0] = abs(aBIC.X[aBIC.X!=0][1])*10
    aBIC.Z[is.na(aBIC.Z)] = 10^10
    aBIC.Z[aBIC.Z==0] = abs(aBIC.Z[aBIC.Z!=0][1])*10
    aBIC = aBIC.X + aBIC.Z
    n_lam1 = which(aBIC[(L3+L2+1):(L3+L2+L1)] == min(aBIC[(L3+L2+1):(L3+L2+L1)]))[1];lam1 = lambda1[n_lam1]
    
    # search lam4
    for (l4 in 1:L4) {
      lam4 = lambda4[l4];l = L3+L2+L1+l4
      
      PP = HHGGM(data.X, data.Z, K, lam1, lam2, lam3, lam4, initialization=F, initialize=out.initial, average=average, 
                 asymmetric=asymmetric, local_appro=local_appro, cutoff.X=cutoff.X, cutoff.Z=cutoff.Z, hier.enf=hier.enf)
      mu_hat.X=PP$mu.X;Theta_hat.X=PP$Xi.X;L.mat.X = PP$L.mat0.X;group.X = PP$group.X;prob.X = PP$prob0.X;aBIC.X[l] = PP$bic.X; member.X = PP$member.X
      Mu_hat.list.X[[l]]=mu_hat.X; Theta_hat.list.X[[l]]=Theta_hat.X; prob.list.X[[l]]=prob.X; L.mat.list.X[[l]] = L.mat.X; member.list.X[[l]]=member.X
      mu_hat.Z=PP$mu.Z;Theta_hat.Z=PP$Xi.Z;L.mat.Z = PP$L.mat0.Z;group.Z = PP$group.Z;prob.Z = PP$prob0.Z;aBIC.Z[l] = PP$bic.Z; member.Z = PP$member.Z
      Mu_hat.list.Z[[l]]=mu_hat.Z; Theta_hat.list.Z[[l]]=Theta_hat.Z; prob.list.Z[[l]]=prob.Z; L.mat.list.Z[[l]] = L.mat.Z; member.list.Z[[l]]=member.Z
      
      fit.error.X.v[l] = PP$fit.error.X;  df.X.v[l] = PP$df.X
      fit.error.Z.v[l] = PP$fit.error.Z;  df.Z.v[l] = PP$df.Z
      group.list[[l]] = PP$group
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 lam4 ="),c(round(lam1,2),round(lam2,2),round(lam3,2),round(lam4,2)),":"),
                paste("K1=",as.numeric(dim(Theta_hat.X)[3]), "K2=",as.numeric(dim(Theta_hat.Z)[3]))))
      }
    }
    aBIC.X[is.na(aBIC.X)] = 10^10
    aBIC.X[aBIC.X==0] = abs(aBIC.X[aBIC.X!=0][1])*10
    aBIC.Z[is.na(aBIC.Z)] = 10^10
    aBIC.Z[aBIC.Z==0] = abs(aBIC.Z[aBIC.Z!=0][1])*10
    aBIC = aBIC.X + aBIC.Z
    n_lam4 = which(aBIC[(L3+L2+L1+1):(L3+L2+L1+L4)] == min(aBIC[(L3+L2+L1+1):(L3+L2+L1+L4)]))[1];lam4 = lambda4[n_lam4]
    
    }
  
  
  if(sum(is.na(aBIC)) == length(aBIC)){
    aBIC = rep(10^10,length(aBIC))
  }
  
  fit.error.X.v[is.na(fit.error.X.v)] = 10^10
  fit.error.Z.v[is.na(fit.error.Z.v)] = 10^10
  df.X.v[is.na(df.X.v)] = 10^10
  df.Z.v[is.na(df.Z.v)] = 10^10
  scale.aBIC = fit.error.X.v/q + fit.error.Z.v/p + df.X.v/(q+q^2) + df.Z.v/(p+p^2)
  if(sum(is.na(scale.aBIC)) == length(scale.aBIC)){
    scale.aBIC = aBIC
  } else {
    scale.aBIC[is.na(scale.aBIC)] = 10^10
  }
  
  K.list <- rep(1,length(Theta_hat.list.X))
  for (l in 1:length(Theta_hat.list.X)) {
    K.list[l] <- as.numeric(dim(Theta_hat.list.X[[l]])[3])
  }
  aBIC[which(K.list == 1)] <- 10^10
  scale.aBIC[which(K.list == 1)] <- 10^10
  
  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)
  Opt_lambda = c(lam1,lam2,lam3,lam4)
  Opt_K1 = as.numeric(dim(Theta_hat.list.X[[n_lam]])[3])
  Opt_K2 = as.numeric(dim(Theta_hat.list.Z[[n_lam]])[3])
  
  n_lam.s = which(scale.aBIC == min(scale.aBIC))[1]
  Opt_aBIC.s = min(scale.aBIC)
  Opt_K1.s = as.numeric(dim(Theta_hat.list.X[[n_lam.s]])[3])
  Opt_K2.s = as.numeric(dim(Theta_hat.list.Z[[n_lam.s]])[3])
  memb.mat = NULL
  
  result = list(Opt_lambda=Opt_lambda, Opt_aBIC=Opt_aBIC, Opt_num=n_lam, Opt_K1=Opt_K1, Opt_K2=Opt_K2,
                Opt_aBIC.s=Opt_aBIC.s, Opt_num.s=n_lam.s, Opt_K1.s=Opt_K1.s, Opt_K2.s=Opt_K2.s,
                Mu_hat.list.X=Mu_hat.list.X, Theta_hat.list.X=Theta_hat.list.X, prob.list.X=prob.list.X, 
                member.list.X=member.list.X, L.mat.list.X=L.mat.list.X, BIC.X=aBIC.X, fit.error.X=fit.error.X.v, df.X=df.X.v,
                Mu_hat.list.Z=Mu_hat.list.Z, Theta_hat.list.Z=Theta_hat.list.Z, prob.list.Z=prob.list.Z, 
                member.list.Z=member.list.Z, L.mat.list.Z=L.mat.list.Z, BIC.Z=aBIC.Z, fit.error.Z=fit.error.Z.v, df.Z=df.Z.v,
                memb.mat=memb.mat, group.list = group.list)
  return(result)
}

Hier.Hetero.GGM = function(data.X, data.Z, lambda, lambda.refit.X, lambda.refit.Z, K, initial.selection="K-means", initialize, 
                           eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5, initialization=T, average=F, asymmetric=T, 
                           local_appro=T, trace = F, cutoff.X=5, cutoff.Z=20, refit.method = "JGL", 
                           BIC.scale=F, hier.enf=F, merge.max=F, refitting=T, hier.enf.refitting=F){
  
  set.seed(1)
  PP = tuning.lambda.HHGGM(data.X, data.Z, lambda, K, initial.selection=initial.selection, initialize, 
                           average=average, asymmetric=asymmetric, eps = eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA,
                           local_appro=local_appro, trace = trace, cutoff.X=cutoff.X, cutoff.Z=cutoff.Z, 
                           hier.enf=hier.enf, merge.max=merge.max, hier.enf.refitting=hier.enf.refitting)
  K_hat.Z = PP$Opt_K2
  K_hat.X = PP$Opt_K1
  
  if(K_hat.Z == 1){
    print("only one cluster")
    return(PP)
  } else if(refitting){
    if(refit.method == "GGMPF"){
      PP.refit.GGM <- HHGGM.refit(PP, data.X, data.Z, lambda.refit.X, lambda.refit.Z, initial.selection=initial.selection, initialize, BIC.scale=BIC.scale)
      PP.refit.JGL=NULL
    } else if(refit.method == "JGL"){
      PP.refit.JGL <- JGL.refit(PP, data.X, data.Z, lambda.refit.X, lambda.refit.Z, BIC.scale=BIC.scale, hier.enf=hier.enf)
      PP.refit.GGM=NULL
    } else if(refit.method == "BOTH"){
      PP.refit.GGM <- HHGGM.refit(PP, data.X, data.Z, lambda.refit.X, lambda.refit.Z, initial.selection=initial.selection, initialize, BIC.scale=BIC.scale)
      PP.refit.JGL <- JGL.refit(PP, data.X, data.Z, lambda.refit.X, lambda.refit.Z, BIC.scale=BIC.scale)
    } 
    hier.enf.res = NULL
    if(refit.method == "JGL" & hier.enf.refitting){
      hier.enf.res <- hier.enf.refit(PP.refit.JGL, data.X, data.Z)
    }
    return(list(pro=PP, hier.enf.res=hier.enf.res,
                PP.refit.GGM=PP.refit.GGM, PP.refit.JGL=PP.refit.JGL))
  } else {return(PP)}
    
}


HHGGM.refit = function(PP, data.X, data.Z, lambda.refit.X, lambda.refit.Z, initial.selection="K-means", initialize, BIC.scale=T)
{
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: HHGGM.refit
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Refitting when K0_hat is identified.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: HHGGM()    Update_Theta_Hier()    AMA_XI_Hier()
  ##            R packages: NbClust
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ lambda1: a float value, the tuning parameter controlling the sparse of the mean parameter. 
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ niter: int, Maximum number of cycles of the algorithm, the default setting is 20.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ initialization: the logical variable, whether to calculate the initial value, the default setting is T,
  ##                                     if initialization = F, the initial value uses initialize.
  ## @ initialize: A given initial value used if initialization = F.
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among 
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list PP.refit including:
  ## @ mu: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups (K0_hat is the number of identified subgroups).
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ niter: int, the actual number of cycles of the algorithm.
  ## @ diff_Xi: a float value, the L2-norm difference between subgroups.
  ## @ prob0: K0_hat * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ L.mat0: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## @ Theta0: p * p * K array, the estimated original precision matrices of given K subgroups.
  ## @ Xi0: p * p * K array, the estimated original dual variables of the precision matrices of given K subgroups in ADMM.
  ## @ mu0: K * p matrix, the estimated original mean vectors of K subgroups.
  ## @ group: a list, the sequence number partition of original K subgroups compressed to K0_hat.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ fit.error: a float value, the value of the loss function (without penalty function) corresponding the choice of given tuning parameters.
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  if(BIC.scale){
    K_hat.Z = PP$Opt_K2.s
    K_hat.X = PP$Opt_K1.s
  } else {
    K_hat.Z = PP$Opt_K2
    K_hat.X = PP$Opt_K1
  }
  
  
    # X
    lambda1 = lambda.refit.X$lambda1
    lambda2 = lambda.refit.X$lambda2
    L1 = length(lambda1)
    L2 = length(lambda2)
    L = L1+L2
    
    aBIC = rep(0,L)
    data.refit = data.X
    n_all = dim(data.refit)[1]
    K = K_hat.X
    if(initial.selection=="K-means"){
      set.seed(1)
      out.initial = initialize_fuc(data.refit,K)
      memb = out.initial$memb
      L.mat = matrix(0,n_all,K)
      for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
      out.initial$L.mat = L.mat
    } else if(initial.selection=="dbscan"){
      out.initial = initialize_fuc.dbscan(data.refit,K)
      memb = out.initial$memb
      L.mat = matrix(0,n_all,K)
      for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
      out.initial$L.mat = L.mat
    } else {out.initial = initialize}
    
    # search lam2
    for (l2 in 1:L2) {
      lam2 = lambda2[l2];l = l2;lam1 = 0
      PP.refit = FGGM(data.refit, K, lam1, lam2, 0, initialization=F, initialize=out.initial)
      aBIC[l2] = PP.refit$bic
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam2 = which(aBIC[1:L2] == min(aBIC[1:L2]))[1];lam2 = lambda2[n_lam2]
    # search lam1
    for (l1 in 1:L1) {
      lam1 = lambda1[l1];l = L2+l1
      PP.refit = FGGM(data.refit, K, lam1, lam2, 0, initialization=F, initialize=out.initial)
      aBIC[l] = PP.refit$bic
    }
    n_lam1 = which(aBIC[(L2+1):(L2+L1)] == min(aBIC[(L2+1):(L2+L1)]))[1];lam1 = lambda1[n_lam1]
    
    PP.refit.X = FGGM(data.refit, K, lam1, lam2, 0)
    opt_Mu_hat.X <- PP.refit.X$mu
    opt_Theta_hat.X <- PP.refit.X$Xi
    opt_L.mat.X <- PP.refit.X$L.mat0
    opt_prob.X <- PP.refit.X$prob0
    
    
    # Z
    lambda1 = lambda.refit.Z$lambda1
    lambda2 = lambda.refit.Z$lambda2
    L1 = length(lambda1)
    L2 = length(lambda2)
    L = L1+L2
    
    aBIC = rep(0,L)
    data.refit = data.Z
    n_all = dim(data.refit)[1]
    K = K_hat.Z
    if(initial.selection=="K-means"){
      set.seed(1)
      out.initial = initialize_fuc(data.refit,K)
      memb = out.initial$memb
      L.mat = matrix(0,n_all,K)
      for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
      out.initial$L.mat = L.mat
    } else if(initial.selection=="dbscan"){
      out.initial = initialize_fuc.dbscan(data.refit,K)
      memb = out.initial$memb
      L.mat = matrix(0,n_all,K)
      for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
      out.initial$L.mat = L.mat
    } else {out.initial = initialize}
    
    # search lam2
    for (l2 in 1:L2) {
      lam2 = lambda2[l2];l = l2;lam1 = 0
      PP.refit = FGGM(data.refit, K, lam1, lam2, 0, initialization=F, initialize=out.initial)
      aBIC[l2] = PP.refit$bic
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam2 = which(aBIC[1:L2] == min(aBIC[1:L2]))[1];lam2 = lambda2[n_lam2]
    # search lam1
    for (l1 in 1:L1) {
      lam1 = lambda1[l1];l = L2+l1
      PP.refit = FGGM(data.refit, K, lam1, lam2, 0, initialization=F, initialize=out.initial)
      aBIC[l] = PP.refit$bic
    }
    n_lam1 = which(aBIC[(L2+1):(L2+L1)] == min(aBIC[(L2+1):(L2+L1)]))[1];lam1 = lambda1[n_lam1]
    
    PP.refit.Z = FGGM(data.refit, K, lam1, lam2, 0)
    opt_Mu_hat.Z <- PP.refit.Z$mu
    opt_Theta_hat.Z <- PP.refit.Z$Xi
    opt_L.mat.Z <- PP.refit.Z$L.mat0
    opt_prob.Z <- PP.refit.Z$prob0
    
    PP.refit <- list(results=PP,results.refit.X=PP.refit.X,results.refit.Z=PP.refit.Z, K1=K_hat.X, K2=K_hat.Z)
    return(PP.refit)
}

JGL.refit = function(PP, data.X, data.Z, lambda.X, lambda.Z, BIC.scale=F, hier.enf=F)
{
  if(BIC.scale){
    K_hat.Z = PP$Opt_K2.s
    K_hat.X = PP$Opt_K1.s
    member.X = apply(PP$L.mat.list.X[[PP$Opt_num.s]],1,which.max)
    member.Z = apply(PP$L.mat.list.Z[[PP$Opt_num.s]],1,which.max)
    if(hier.enf){ 
      member.X = PP$member.list.X
      member.Z = PP$member.list.Z
    }
  } else {
    K_hat.Z = PP$Opt_K2
    K_hat.X = PP$Opt_K1
    member.X = apply(PP$L.mat.list.X[[PP$Opt_num]],1,which.max)
    member.Z = apply(PP$L.mat.list.Z[[PP$Opt_num]],1,which.max)
    if(hier.enf){ 
      member.X = PP$member.list.X[[PP$Opt_num]]
      member.Z = PP$member.list.Z[[PP$Opt_num]]
      }
  }
  result.X = tuning.lambda.JGL(lambda.X, data.X, member.X)
  result.Z = tuning.lambda.JGL(lambda.Z, data.Z, member.Z)
  PP.refit <- list(results=PP, results.refit.X=result.X, results.refit.Z=result.Z, 
                   K1=K_hat.X, K2=K_hat.Z)
  return(PP.refit)
}

memb.JGL <- function(data.X, member.X, lambda2)
{
  n = dim(data.X)[1]
  q = dim(data.X)[2]
  K.X = length(unique(member.X))
  t.X.memb = as.data.frame(table(member.X))
  t.X.memb[,1] = as.numeric(as.character(t.X.memb[,1]))
  mu.hat.X = matrix(0, K.X, q)
  prob.X = rep(0,K.X)
  cluster.data.X = vector("list", length = K.X)
  for(k in 1:K.X) {
    if(sum(member.X==t.X.memb[k,1]) > 1){
      cluster.data.X[[k]] = data.X[member.X==t.X.memb[k,1],]
    } else {
      cluster.data.X[[k]] = data.X[unique(c(1,which(member.X==t.X.memb[k,1]),n)),]
    }
    mu.hat.X[k,] = apply(cluster.data.X[[k]],2,mean)
    prob.X[k] = sum(member.X == t.X.memb[k,1])/length(member.X)
  }
  
  Theta.X = array(0, dim = c(q, q, K.X))
  out.jgl.X <- JGL(Y=cluster.data.X, penalty="group",lambda2,0,return.whole.theta=TRUE)
  for(k in 1:K.X) {
    Theta.X[,,k] = out.jgl.X$theta[[k]]
  }
  
  bic.X = BIC.JGL(data.X, mu.hat.X, Theta.X, prob.X)
  JGL_res <- list();JGL_res$mu <- mu.hat.X;JGL_res$Theta <- Theta.X
  JGL_res$prob <- prob.X;JGL_res$bic <- bic.X$bic; JGL_res$member <- member.X
  return(JGL_res)
  
}

tuning.lambda.JGL <- function(lambda, data, member)
{
  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  lambda3 = lambda$lambda3
  L1 = 0
  L2 = length(lambda2)
  L3 = 0
  L = L1+L2+L3
  
  aBIC = rep(0,L)
  
  K = length(unique(member))
  p = dim(data)[2]
  Theta = array(0, dim = c(p, p, K))
  mu.hat = matrix(0, K, p)
  prob = rep(0,K)
  cluster.data = vector("list", length = K)
  for(k in 1:K) {
    cluster.data[[k]] = data[member==k,]
    mu.hat[k,] = apply(cluster.data[[k]],2,mean)
    prob[k] = sum(member == k)/length(member)
  }
  
  
  Mu_hat.list = list()
  Theta_hat.list = list()
  prob.list = list()
  L.mat.list = list()
  member.list = list()
  # search lam2
  for (l2 in 1:L2) {
    lam2 = lambda2[l2];l = l2
    SCAN = memb.JGL(data, lambda2=lam2, member)
    mu_hat <- SCAN$mu; Theta_hat <- SCAN$Theta; prob=SCAN$prob
    aBIC[l] = SCAN$bic; member <- SCAN$member
    Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; 
    prob.list[[l]]=prob; member.list[[l]] <- member
  }
  aBIC[is.na(aBIC)] = abs(aBIC[!is.na(aBIC)][1])*10
  
  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)
  
  result = list(Mu_hat.list=Mu_hat.list,Theta_hat.list=Theta_hat.list,
                prob.list=prob.list, member.list = member.list,
                Opt_aBIC=Opt_aBIC,BIC=aBIC,Opt_num=n_lam, lambda=lambda2)
  return(result)
}

BIC.JGL <- function(data, mu_hat, Theta_hat, prob)
{
  n = nrow(data)
  K = nrow(mu_hat)
  
  # fitting error
  pi_vec = prob
  fit.error_mat = matrix(0, n, K)
  for(k in 1:K) {
    fit.error_mat[,k] = pi_vec[k] * f.den.vec( data, as.numeric(mu_hat[k,]), Theta_hat[,,k] )                  
  }
  fit0 = apply(fit.error_mat, 1, sum)
  fit.error = sum(log( fit0 + min(fit0[fit0>0]) ))
  fit.error = - 2*fit.error
  
  
  # degrees of freedom
  for(i in 1:K){
    Theta_hat[upper.tri(Theta_hat[, , i], diag = T)] = 0
  }
  
  df =  log(n) * length(which(mu_hat != 0)) + 2 * length(which(Theta_hat != 0))
  bic = fit.error + df
  P = list()
  P$fit.error = fit.error
  P$df = df
  P$bic = bic
  return(P)
}

hier.enf.refit <- function(PP.refit.JGL, data.X, data.Z){
  res0 <- PP.refit.JGL
  results.X <- res0$results.refit.X
  results.Z <- res0$results.refit.Z
  o.x <- results.X$Opt_num
  o.z <- results.Z$Opt_num
  
  mu_final.X = results.X$Mu_hat.list[[o.x]]
  Theta_final.X = results.X$Theta_hat.list[[o.x]]
  prob.X = results.X$prob.list[[o.x]]
  
  mu_final.Z = results.Z$Mu_hat.list[[o.z]]
  Theta_final.Z = results.Z$Theta_hat.list[[o.z]]
  prob.Z = results.Z$prob.list[[o.z]]
  
  Opt_K1 = as.numeric(dim(Theta_final.X)[3])
  Opt_K2 = as.numeric(dim(Theta_final.Z)[3])
  
  n = dim(data.X)[1]
  f.mat.X = matrix(0,n,Opt_K1)
  f.mat.Z = matrix(0,n,Opt_K2)
  L.mat.X = matrix(0,n,Opt_K1)
  L.mat.Z = matrix(0,n,Opt_K2)
  # calculate pdfs
  for(k.ind in 1:Opt_K1) {
    f.mat.X[,k.ind]=f.den.vec( data.X, as.numeric(mu_final.X[k.ind,]), Theta_final.X[,,k.ind] )
  }
  for(k.ind in 1:Opt_K2) {
    f.mat.Z[,k.ind]=f.den.vec( data.Z, as.numeric(mu_final.Z[k.ind,]), Theta_final.Z[,,k.ind] ) 
  }
  # update L and pi
  for(i in 1:n) {
    for(k.ind in 1:Opt_K1) {
      L.mat.X[i,k.ind] = prob.X[k.ind] * f.mat.X[i,k.ind] / prob.X %*% f.mat.X[i,]
    }
    for(k.ind in 1:Opt_K2) {
      L.mat.Z[i,k.ind] = prob.Z[k.ind] * f.mat.Z[i,k.ind] / prob.Z %*% f.mat.Z[i,]
    }
  }
  
  member.Z = apply(L.mat.Z,1,which.max)
  
  prob.XZ = matrix(0, Opt_K1, Opt_K2)
  memb.XZ = prob.XZ
  for (k.ind in 1:Opt_K2) {
    prob.XZ[,k.ind] = apply(L.mat.X[member.Z==k.ind,],2,mean)
    memb.XZ[,k.ind] = as.numeric(prob.XZ[,k.ind] == max(prob.XZ[,k.ind]))
  }
  
  q <- as.integer(dim(data.X)[2])
  memb.XZ.vec = apply(memb.XZ,2,which.max)
  t.memb.XZ = as.data.frame(table(memb.XZ.vec))
  t.memb.XZ[,1] = as.numeric(as.character(t.memb.XZ[,1]))
  K1 = dim(t.memb.XZ)[1]
  mu_final.X0 = matrix(0, K1, q)
  Theta_final.X0 = array(0, dim = c(q, q, K1))
  for (k1 in 1:K1) {
    mu_final.X0[k1,] = mu_final.X[t.memb.XZ[k1,1],]
    Theta_final.X0[,,k1] = Theta_final.X[,,t.memb.XZ[k1,1]]
  }
  member.X = rep(0,n)
  memb.XZ0 = as.matrix(memb.XZ[apply(memb.XZ,1,function(x) sum(x)>0),])
  if(min(dim(memb.XZ0)) == 1){memb.XZ0 = t(memb.XZ0)}
  for (k in 1:Opt_K2) {
    member.X[member.Z == k] = which(memb.XZ0[,k] == 1)
  }
  prob.X = rep(0,K1)
  for (k1 in 1:K1) {
    prob.X[k1] = sum(member.X == k1) / n
  }
  
  results.Z.opt = list(Theta_hat=Theta_final.Z,Mu_hat=mu_final.Z,
                       member=member.Z,prob=prob.Z)
  results.X.opt = list(Theta_hat=Theta_final.X0,Mu_hat=mu_final.X0,
                       member=member.X,prob=prob.X)
  return(list(Opt_K1=K1, Opt_K2=Opt_K2,
              results.X.opt=results.X.opt, results.Z.opt=results.Z.opt))
}

hier.enf.HHGGM <- function(data.X, data.Z, group_final.XZ, mu_final.X, mu_final.Z, 
                          Theta_final.X, Theta_final.Z, prob0.X, prob0.Z){
  prob.X = prob0.X
  prob.Z = prob0.Z
  
  Opt_K1 = as.numeric(dim(Theta_final.X)[3])
  Opt_K2 = as.numeric(dim(Theta_final.Z)[3])
  
  n = dim(data.X)[1]
  f.mat.X = matrix(0,n,Opt_K1)
  f.mat.Z = matrix(0,n,Opt_K2)
  L.mat.X = matrix(0,n,Opt_K1)
  L.mat.Z = matrix(0,n,Opt_K2)
  for(k.ind in 1:Opt_K1) {
    f.mat.X[,k.ind]=f.den.vec( data.X, as.numeric(mu_final.X[k.ind,]), Theta_final.X[,,k.ind] )
  }
  for(i in 1:n) {
    for(k.ind in 1:Opt_K1) {
      L.mat.X[i,k.ind] = prob.X[k.ind] * f.mat.X[i,k.ind] / prob.X %*% f.mat.X[i,]
    }
  }
  member.X = apply(L.mat.X,1,which.max)
  
  
  for (k1 in 1:Opt_K1) {
    n.X.K = which(member.X == k1)
    data.Z.k = data.Z[n.X.K,]
    Opt_K2.k1 = length(group_final.XZ[[k1]])
    
    for(k.ind0 in 1:Opt_K2.k1) {
      k.ind = group_final.XZ[[k1]][k.ind0]
      f.mat.Z[n.X.K,k.ind]=f.den.vec( data.Z.k, as.numeric(mu_final.Z[k.ind,]), Theta_final.Z[,,k.ind] ) 
    }
    for(i in n.X.K) {
      for(k.ind0 in 1:Opt_K2.k1) {
        k.ind = group_final.XZ[[k1]][k.ind0]
        L.mat.Z[i,k.ind] = prob.Z[k.ind] * f.mat.Z[i,k.ind] / prob.Z %*% f.mat.Z[i,]
      }
    }
  }
  
  member.Z = apply(L.mat.Z,1,which.max)
  
  prob0.X = rep(0,Opt_K1)
  for (k1 in 1:Opt_K1) {
    prob0.X[k1] = sum(member.X == k1) / n
  }
  prob0.Z = rep(0,Opt_K2)
  for (k2 in 1:Opt_K2) {
    prob0.Z[k2] = sum(member.Z == k2) / n
  }
  
  return(list(member.Z=member.Z, prob0.Z=prob0.Z, 
              member.X=member.X, prob0.X=prob0.X,
              L.mat0.X=L.mat.X, L.mat0.Z=L.mat.Z))
}

############################# Some fundamental supporting functions ############################
mcp_d = function(x,lambda,a=3){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: mcp_d
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the derivative of the MCP
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ x: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ rho: the derivative of the MCP.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  if(lambda!=0){
    rho <- lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}

S_soft = function(z,lambda){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: S_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Define the soft threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable.
  ## @ lambda: a float value, the tuning parameter.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the soft threshold operator.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  return((abs(z) - lambda)*(abs(z) - lambda > 0)*sign(z))
}

MCP_soft = function(z,lambda,a=3){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: MCP_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Define the MCP threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the MCP threshold operator.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  return( S_soft(z,lambda)/(1-1/a) * (abs(z) - a*lambda <= 0) + z * (abs(z) - a*lambda > 0) )
}

MCP_soft0 = function(z,lambda,a=3,kappa=1){
  return( (1-lambda/kappa/z)/(1-1/a/kappa) * (1-lambda/kappa/z > 0) )
}

f.den.vec = function(data, mu, Theta){ 
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: f.den.vec
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            calculate the density function values at each sample point.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu1: p * 1 vector, the mean vector.
  ## @ Omega1: p * p matrix, the precision matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ fdensity: The density function values at each sample point.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = length(mu)
  fden = as.numeric( (2*pi)^(-p/2) * (det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))) )
  return(fden)
}

Symmetrize = function(X){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Symmetrize
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Symmetrize the precision matrices using the symmetrization strategy of Cai et al. (2016).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ X: p * p matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ X: a symmetrical matrix.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = dim(X)[1]
  for(i in 1:p){
    for(j in i:p){
      if(X[i,j] < X[j, i]){
        X[j, i] = X[i, j]
      }else{
        X[i, j] = X[j, i]
      }
    }
  }
  return(X)
}

cut_diff_ama = function(V_kk,K_c,K,cutoff=0.01, trace=F){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: cut_diff_ama
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Truncation of small differences in parameter estimates between subgroups.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ V_kk: L2-norm differences in the precision matrices between different subgroups.
  ## @ K_c: All combinations of natural numbers within K.
  ## @ K: int, a selected upper bound of K_0.
  ## @ cutoff: a float value, a given cut-off value.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ K_group_final: the partition of the original K subgroups.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  V_kk_num = which(V_kk < cutoff)
  K_group_final = list()
  if(length(V_kk_num) > 0){
    K_group = list()
    for (j in 1:length(V_kk_num)) {
      K_group[[j]] = K_c[,V_kk_num[j]]
    }
    outnum = setdiff(1:K,Reduce(union,K_group))
    if(length(outnum) > 0){
      for (j in 1:length(outnum)) {
        K_group[[length(V_kk_num)+j]] = outnum[j]
      }
      K_group[[length(V_kk_num)+j+1]] = K
    } else{
      K_group[[length(V_kk_num)+1]] = K
    }
    
    kk = 1
    repeat{
      repeat{
        K_group_old=K_group
        k_del = NULL
        for (kkk in setdiff(1:length(K_group),1) ) {
          if(length(Reduce(intersect,list(K_group[[1]],K_group_old[[kkk]]))) > 0){
            K_group[[1]] = sort(unique(c(K_group[[1]],K_group_old[[kkk]])))
            k_del = c(k_del,kkk)
          }
        }
        if(length(k_del) > 0){
          for (j in sort(k_del,decreasing = T)) {
            K_group[[j]] = NULL
          }
        }
        if(length(K_group_old) == length(K_group)){break}
      }
      K_group_final[[kk]] = K_group[[1]]
      if(kk==1 && length(K_group) == 1){if(trace){print("Warning: Only one cluster!")};break}
      if(length(K_group) == 2){K_group_final[[kk+1]] = K_group[[2]];break}
      if(kk>1 && length(K_group) == 1){break}
      K_group[[1]] = NULL
      kk = kk+1
    }
  }else {
    for (k in 1:K) {
      K_group_final[[k]] = k
    }
  }
  
  return(K_group_final)
}

memb.cluster = function(memb){
  n = length(memb)
  capgfl.matrix = matrix(0,n,n)
  for(i in 1:(n-1)){
    for (j in (i+1):(n)) {
      capgfl.matrix[i,j] <- as.numeric(memb[i] == memb[j])
    }
  }
  capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
  return(capgfl.matrix2)
}
############################# Functions for generating the initial values ############################
initialize_fuc = function(data, K){  
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: initialize_fuc
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the initial values using K-means.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  Mu <- matrix(0, K, p)
  kmeans.clust <- kmeans(data, K)
  memb <- kmeans.clust$cluster
  prob <- kmeans.clust$size/n
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[memb == k, , drop = FALSE]) )
    S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    if(qr(S[,,k])$rank < dim(S[,,k])[2]){
      aa <- diag(S[,,k])
      S[,,k]  <- S[,,k]/2
      diag(S[,,k]) <- aa
    }
    Theta[,,k] <- solve(S[,,k])
  }
  
  int.res <- list()
  int.res$prob <-  prob
  int.res$Mu <-  Mu
  int.res$Theta <- Theta
  int.res$S <- S
  int.res$memb <- memb
  return(int.res)
}

initialize_fuc.dbscan = function(data, K){  
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: initialize_fuc.dbscan
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the initial values using hierarchical clustering.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: 
  ##            R packages: NbClust
  ## -----------------------------------------------------------------------------------------------------------------
  
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  
  # initialization via dbscan clustering## 
  Mu <- matrix(0, K, p)
  hc <- hclust(dist(data,method = "euclidean"),method = "ward.D2")
  memb <- cutree(hc,k=K)
  
  prob <- rep(0,K)
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[memb == k, , drop = FALSE]) )
    S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    Theta[,,k] <- solve(S[,,k] + diag(1,p))
    prob[k] <- sum(memb == k)/n
  }
  
  P <- list()
  P$prob <-  prob
  P$Mu <-  Mu
  P$Theta <- Theta
  P$S <- S
  P$memb <- memb
  return(P)
}
