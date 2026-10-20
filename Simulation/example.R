################################################################################################
#Codes for conducting simulation studies
################################################################################################
rm(list = ls(all = TRUE))
ls()
###############################
setwd("./")
library(MASS)
library(NbClust)
library(JGL)
library(Matrix)
library(igraph)
library(HeteroGGM)
source("../main_function/function.R")
source("../main_function/alternatives.R")
source("sim-func.R")


network.type = network.type.list[3] # 3:power-law; 2:nearest-neighbor; 1:tridiagonal
mue = mu.list[1]                    #1: mue=2; 2: mue=3
source("simulation-parameter.R")


################ The implementation and evaluation ################
################ Generating simulated data ################
whole.data = generate.data.XZcor(n.XZ,N1,Mu0.list.X,Theta0.list.X,Sigma0.list.X,N2,Mu0.list.Z,Theta0.list.Z,Sigma0.list.Z)
whole.data.X = whole.data$whole.data.X
whole.data.Z = whole.data$whole.data.Z
data.X = whole.data.X$data
data.Z = whole.data.Z$data

# kmeans.Hier.JGL.givenK
result.kmeans.Hier.JGL.givenK = kmeans.Hier.JGL.givenK(data.X, data.Z, lambda.km.X, lambda.km.Z, K.X.k=3, K.Z.k=c(1,2,2))
index.kmeans.Hier.JGL.givenK = Esti.error.JGL.all(whole.data.X, whole.data.Z, result.kmeans.Hier.JGL.givenK)
index.givenk.JGL = cbind(index.kmeans.Hier.JGL.givenK$X.JGL,index.kmeans.Hier.JGL.givenK$Z.JGL)

# kmeans.Hier.JGL.selectK
result.kmeans.Hier.JGL.selectK = kmeans.Hier.JGL.selectK(data.X, data.Z, lambda.km.X, lambda.km.Z, K_min=2,K_max=6, k.index= "gap")
index.kmeans.Hier.JGL.selectK = Esti.error.JGL.all(whole.data.X, whole.data.Z, result.kmeans.Hier.JGL.selectK)
index.selectk.JGL = cbind(index.kmeans.Hier.JGL.selectK$X.JGL,index.kmeans.Hier.JGL.selectK$Z.JGL)

# oracle.JGL
result.oracle.JGL = oracle.JGL(whole.data.X, whole.data.Z, lambda.km.X, lambda.km.Z)
index.oracle.JGL = Esti.error.JGL.all(whole.data.X, whole.data.Z, result.oracle.JGL)
index.oracle = cbind(index.oracle.JGL$X.JGL,index.oracle.JGL$Z.JGL)

# proposed
res00 = Hier.Hetero.GGM(data.X, data.Z, lambda, lambda.refit.X, lambda.refit.Z, K)
res0 = res00$PP.refit.JGL
results.X = res0$results.refit.X
results.Z = res0$results.refit.Z
pro.index.refit.X = Esti.error(data.X, results.X$Mu_hat.list[[results.X$Opt_num]], results.X$Theta_hat.list[[results.X$Opt_num]], whole.data.X$Mu0, whole.data.X$Theta0, K1, L.mat=c(0,0), whole.data.X$L0, results.X$prob0, results.X$member.list[[results.X$Opt_num]])
pro.index.refit.Z = Esti.error(data.Z, results.Z$Mu_hat.list[[results.Z$Opt_num]], results.Z$Theta_hat.list[[results.Z$Opt_num]], whole.data.Z$Mu0, whole.data.Z$Theta0, K2, L.mat=c(0,0), whole.data.Z$L0, results.Z$prob0, results.Z$member.list[[results.Z$Opt_num]])
index.pro = cbind(pro.index.refit.X, pro.index.refit.Z)


index.givenk.JGL
index.selectk.JGL
index.pro
index.oracle


