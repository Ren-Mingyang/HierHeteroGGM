if(network.type == "tridiagonal"){
  ################ The true parameters X ################
  nonnum <- 2
  mu00 <- c(rep(0,2*nonnum),rep(0,q-2*nonnum))
  mu01 <- c(rep(mue+0.5,2*nonnum),rep(0,q-2*nonnum))
  mu02 <- c(rep(-mue-0.5,2*nonnum),rep(0,q-2*nonnum))
  Theta00 <- tridiag.cor(q, s=5, rho=0)
  Theta01 <- tridiag.cor(q, s=5, rho=-0.2)
  Theta02 <- tridiag.cor(q, s=5, rho=0.2)
  sigma00 <- solve(Theta00)
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  Mu0.list.X <- list(mu00,mu01,mu02)
  Sigma0.list.X <- list(sigma00,sigma01,sigma02)
  Theta0.list.X <- list(Theta00,Theta01,Theta02)
  
  ################ The true parameters Z ################
  nonnum <- 4
  mu00 <- c(rep(0,2*nonnum),rep(0,p-2*nonnum))
  mu01 <- c(rep(mue,nonnum),rep(-mue,nonnum),rep(0,p-2*nonnum))
  mu02 <- c(rep(-mue,nonnum),rep(mue,nonnum),rep(0,p-2*nonnum))
  mu03 <- c(rep(mue,2*nonnum),rep(0,p-2*nonnum))
  mu04 <- c(rep(-mue,2*nonnum),rep(0,p-2*nonnum))
  Theta00 <- tridiag.cor(p, s=10, rho=0)
  Theta01 <- tridiag.cor(p, s=10, rho=0.2)
  Theta02 <- tridiag.cor(p, s=10, rho=-0.3)
  Theta03 <- tridiag.cor(p, s=10, rho=-0.2)
  Theta04 <- tridiag.cor(p, s=10, rho=0.3)
  sigma00 <- solve(Theta00)
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  sigma03 <- solve(Theta03)
  sigma04 <- solve(Theta04)
  Mu0.list.Z <- list(mu00,mu01,mu02,mu03,mu04)
  Sigma0.list.Z <- list(sigma00,sigma01,sigma02,sigma03,sigma04)
  Theta0.list.Z <- list(Theta00,Theta01,Theta02,Theta03,Theta04)

}


if(network.type == "nearest-neighbor"){
  ################ The true parameters X ################
  nonnum <- 2
  mu00 <- c(rep(0,2*nonnum),rep(0,q-2*nonnum))
  mu01 <- c(rep(mue+0.5,2*nonnum),rep(0,q-2*nonnum))
  mu02 <- c(rep(-mue-0.5,2*nonnum),rep(0,q-2*nonnum))
  I.list <- list(I1=c(0),I2=c(1),I3=c(2))
  set.seed(1)
  A.list <- nearest.neighbor.network(q,s=5,near=3,I.list=I.list)
  Theta00 <- A.list[[1]]
  Theta01 <- A.list[[2]]
  Theta02 <- A.list[[3]]
  sigma00 <- solve(Theta00)
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  Mu0.list.X <- list(mu00,mu01,mu02)
  Sigma0.list.X <- list(sigma00,sigma01,sigma02)
  Theta0.list.X <- list(Theta00,Theta01,Theta02)
  
  ################ The true parameters Z ################
  nonnum <- 4
  mu00 <- c(rep(0,2*nonnum),rep(0,p-2*nonnum))
  mu01 <- c(rep(mue,nonnum),rep(-mue,nonnum),rep(0,p-2*nonnum))
  mu02 <- c(rep(-mue,nonnum),rep(mue,nonnum),rep(0,p-2*nonnum))
  mu03 <- c(rep(mue,2*nonnum),rep(0,p-2*nonnum))
  mu04 <- c(rep(-mue,2*nonnum),rep(0,p-2*nonnum))
  I.list <- list(I1=c(0),I2=c(1),I3=c(2),I4=c(3),I5=c(4))
  set.seed(1)
  A.list <- nearest.neighbor.network(p,s=10,near=3,I.list=I.list)
  Theta00 <- A.list[[1]]
  Theta01 <- A.list[[2]]
  Theta02 <- A.list[[3]]
  Theta03 <- A.list[[4]]
  Theta04 <- A.list[[5]]
  # deii <- matrix(0,length(A.list),length(A.list))
  # for (ii in 1:length(A.list)) {
  #   for (jj in 1:length(A.list)){
  #     deii[ii,jj] <- sum((A.list[[ii]] - A.list[[jj]])^2)
  #   }
  # }
  sigma00 <- solve(Theta00)
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  sigma03 <- solve(Theta03)
  sigma04 <- solve(Theta04)
  Mu0.list.Z <- list(mu00,mu01,mu02,mu03,mu04)
  Sigma0.list.Z <- list(sigma00,sigma01,sigma02,sigma03,sigma04)
  Theta0.list.Z <- list(Theta00,Theta01,Theta02,Theta03,Theta04)
  
}


if(network.type == "power-law"){
  ################ The true parameters X ################
  nonnum <- 2
  mu00 <- c(rep(0,2*nonnum),rep(0,q-2*nonnum))
  mu01 <- c(rep(mue+0.5,2*nonnum),rep(0,q-2*nonnum))
  mu02 <- c(rep(-mue-0.5,2*nonnum),rep(0,q-2*nonnum))
  I.list <- list(I1=c(0),I2=c(1),I3=c(2))
  set.seed(1)
  A.list <- Power.law.network(q,s=5,umin=0.1,umax=0.4,I.list)
  Theta00 <- A.list[[1]]
  Theta01 <- A.list[[2]]
  Theta02 <- A.list[[3]]
  # deii <- matrix(0,length(A.list),length(A.list))
  # for (ii in 1:length(A.list)) {
  #   for (jj in 1:length(A.list)){
  #     deii[ii,jj] <- sum((A.list[[ii]] - A.list[[jj]])^2)
  #   }
  # }
  sigma00 <- solve(Theta00)
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  Mu0.list.X <- list(mu00,mu01,mu02)
  Sigma0.list.X <- list(sigma00,sigma01,sigma02)
  Theta0.list.X <- list(Theta00,Theta01,Theta02)
  
  ################ The true parameters Z ################
  nonnum <- 4
  mu00 <- c(rep(0,2*nonnum),rep(0,p-2*nonnum))
  mu01 <- c(rep(mue,nonnum),rep(-mue,nonnum),rep(0,p-2*nonnum))
  mu02 <- c(rep(-mue,nonnum),rep(mue,nonnum),rep(0,p-2*nonnum))
  mu03 <- c(rep(mue,2*nonnum),rep(0,p-2*nonnum))
  mu04 <- c(rep(-mue,2*nonnum),rep(0,p-2*nonnum))
  I.list <- list(I1=c(0),I2=c(1),I3=c(2),I4=c(3),I5=c(4))
  set.seed(1)
  A.list <- Power.law.network(p,s=20,umin=0.1,umax=0.4,I.list)
  Theta00 <- A.list[[1]]
  Theta01 <- A.list[[2]]
  Theta02 <- A.list[[3]]
  Theta03 <- A.list[[4]]
  Theta04 <- A.list[[5]]
  # deii <- matrix(0,length(A.list),length(A.list))
  # for (ii in 1:length(A.list)) {
  #   for (jj in 1:length(A.list)){
  #     deii[ii,jj] <- sum((A.list[[ii]] - A.list[[jj]])^2)
  #   }
  # }
  sigma00 <- solve(Theta00)
  sigma01 <- solve(Theta01)
  sigma02 <- solve(Theta02)
  sigma03 <- solve(Theta03)
  sigma04 <- solve(Theta04)
  Mu0.list.Z <- list(mu00,mu01,mu02,mu03,mu04)
  Sigma0.list.Z <- list(sigma00,sigma01,sigma02,sigma03,sigma04)
  Theta0.list.Z <- list(Theta00,Theta01,Theta02,Theta03,Theta04)
  
}




