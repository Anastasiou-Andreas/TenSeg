#===============================================================
# Author: Andreas Anastasiou
# Project: Tensor time series change-point detection in cryptocurrency network data

# Description: The code below provides the TenSeg routine as well as some main
# functions that are necessary in in order to be able to run the simulation study
# in the paper.

#===============================================================
rm(list = ls())   # Clear environment

library(rTensor)
library(ccid)
library(Matrix)
library(Rcpp)
library(RcppEigen)
sourceCpp("helper_funcs.cpp")
#==============================
# Some utility functions
#==============================
kronecker_sum_sparse <- function(A,B){
  m = NCOL(A)
  n = NCOL(B)
  return(as(kronecker(Diagonal(n),A) + kronecker(B,Diagonal(m)), "dgCMatrix"))
}

kway_kronecker_sum_sparse <- function(matrices_list){
  K <- length(matrices_list)
  if (K == 1) {
    return(matrices_list[[1]])
  } else{
    A <- kronecker_sum_sparse(matrices_list[[1]],matrices_list[[2]])
    if (K == 2) {
      return(A)
    } else {
      for (k in 3:K){
        A <- kronecker_sum_sparse(A,matrices_list[[k]])
      }
      return(A)
    }
  }
}

# Generate sylvester matrix (standardized) normal variable with relationship
# AX + XB =  C
gen_sylvester_data <- function(N, dims, Omega_sqrt, ac_coeff = 0){
  Omega_sqrt <- as(Omega_sqrt, "dgCMatrix")
  
  if (N == 1){
    Z = t(stdmvrnormArma(N,prod(dims)))
    X = mvrnormEigen(Z,Omega_sqrt)
    X = array(X[,1],dim = c(dims,1))
    X = as.tensor(X)
  } else{
    Z = t(stdmvrnormArma(N,prod(dims)))
    X = mvrnormEigen(Z,Omega_sqrt)
    X = apply(X,1,scale)
    if (ac_coeff !=0){
      for (i in 1:dim(X)[2]){
        for (j in 2:dim(X)[1]){
          X[j,i] = ac_coeff * X[j-1,i] + X[j,i]
        }
      }
    }
    X = foreach(i=1:nrow(X),.combine=function(x,y) abind(x,y,rev.along=1)) %do% array(X[i,],dim=c(dims,1))
    X = as.tensor(X)
  }
  return(X)
}

# Generate sparse inverse covariance matrix that has auto-regressive structure
# and its corresponding covariance matrix
AR_model <- function(m, rho) {
  cov = matrix(0, nrow = m, ncol = m)
  inv = matrix(0, nrow = m, ncol = m)
  
  cov = rho^abs(row(cov) - col(cov))
  
  diag(inv) = 1+rho^2
  inv[1,1] = 1
  inv[m,m] = 1
  diag(inv[-1,]) = -rho
  diag(inv[,-1]) = -rho
  inv = inv / (1-rho^2)
  
  return(list(cov=cov, inv=inv))
}

# Generate sparse inverse covariance matrix that has star-block structure and
# its corresponding covariance matrix
SB_model <- function(m, rho, num_subgraph=4, tol=1e-5){
  size_subgraph = floor(m/num_subgraph)
  unequal_graph = (m != size_subgraph * num_subgraph)
  
  list_subgraph = list()
  list_invsubgraph = list()
  for (i in 1:num_subgraph){
    if(i==num_subgraph && unequal_graph){
      size_subgraph = m - size_subgraph*(i-1)
    }
    central_node = sample(1:size_subgraph,1)
    subgraph = matrix(rho^2, nrow=size_subgraph, ncol = size_subgraph)
    subgraph[,central_node] = rho
    subgraph[central_node,] = rho
    diag(subgraph)=1
    inv_subgraph = list(solve(subgraph))
    subgraph = list(subgraph)
    list_subgraph = c(list_subgraph, subgraph)
    list_invsubgraph = c(list_invsubgraph, inv_subgraph)
  }
  
  B = as.matrix(bdiag(list_subgraph))
  B[abs(B)<tol] = 0
  B_inv = as.matrix(bdiag(list_invsubgraph))
  B_inv[abs(B_inv)<tol] = 0
  
  return(list(cov=B, inv=B_inv))
}

# Generate sparse inverse covariance matrix that has Erdos-Renyi model
# and its corresponding covariance matrix
ER_model <- function(m, dens, wmin = 0.1, wmax = 0.3, corr=TRUE, tol=1e-5, tau_B=NA){
  B_inv = 0.25*diag(m)
  row_sample = sample(1:m, size = dens)
  col_sample = sample(1:m, size = dens)
  rand_weight = runif(dens, min = wmin, max = wmax)
  
  for (i in 1:dens) {
    B_inv[row_sample[i], col_sample[i]] = B_inv[row_sample[i], col_sample[i]] - rand_weight[i]
    B_inv[col_sample[i], row_sample[i]] = B_inv[row_sample[i], col_sample[i]]
    B_inv[row_sample[i], row_sample[i]] = B_inv[row_sample[i], row_sample[i]] + rand_weight[i]
    B_inv[col_sample[i], col_sample[i]] = B_inv[col_sample[i], col_sample[i]] + rand_weight[i]
  }
  
  B = solve(B_inv)
  
  
  if(corr){
    diag_B = 1/sqrt(diag(B))
    B = diag(diag_B) %*% B %*% diag(diag_B)
    B_inv = solve(B)
  }
  
  B[abs(B)<tol] = 0
  B_inv[abs(B_inv)<tol] = 0
  return(list(cov=B, inv=B_inv))
}

#============================================================================
# Our main TenSeg function when directly used for a Tensor as an input. This is
# not the case in the simulation study where we first created and decomposed
# the tensors using the functions tensor_creation_decomposition_XX, where XX
# could take the value AR, ER, SB depending on the type of changes we introduce
# in the tensor. Then, all the relevant functions suitable for the decomposed data
# are run on the decomposed outcome. This saves an enormous amount of computational
# time in the simulations.
#============================================================================

TenSeg <- function(tensor = NULL, decomp = "cp", par_decom = NULL, points = 10, thr = 2.1){
if (decomp == "cp"){
  decom <- rTensor:::cp(tensor, num_components = par_decom)
}
if (decomp == "hosvd"){
  decom <- rTensor:::hosvd(tensor, ranks = c(z.tnsr@modes[-4],par_decom))
}
  z <- ccid:::detect.ic(decom$U[[4]], approach = "infinity", pointsgen = points, th_max = thr)
  return(list("Change-point locations" = z$changepoints, "Number of Change-points" = z$no.of.cpts))
}


# The following function builds tensors according to values given in the input arguments
# with precision matrices that undergo changes in the AR structure. We then employ
# a decomposition method chosen by the user and we return the decomposed data, as well as
# the percent of the Frobenius norm explained by the approximation.

tensor_creation_decomposition_AR <- function(dims = c(20,20,20), N = 100, size = 200, dens = NULL, decomp = "cp", ac = 0, par_decom = NULL){
  K <- length(dims)
  cpt.no <- length(N)
  X <- list()
  if (length(dens) != (cpt.no+1)){stop("density and change-point number do not agree")}
  if(cpt.no == 0){
    Psis <- list()
    for (k in 1:length(dims)){
      Psis <- append(Psis,list(AR_model(dims[k],dens[1])$inv))
    }
    Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
    X[[1]] <- gen_sylvester_data(size,dims,Omega_sqrt, ac_coeff = ac)
    z = X[[1]]@data
  }
  else{
    Psis <- list()
    for (k in 1:length(dims)){
      Psis <- append(Psis,list(AR_model(dims[k],dens[1])$inv))
    }
    Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
    X[[1]] <- gen_sylvester_data(N[1],dims,Omega_sqrt, ac_coeff = ac)
    if(cpt.no == 1){
      Psis <- list()
      for (k in 1:length(dims)){
        Psis <- append(Psis,list(AR_model(dims[k],dens[2])$inv))
      }
      Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
      X[[2]] <- gen_sylvester_data(size - N[1],dims,Omega_sqrt, ac_coeff = ac)
      z = abind(X[[1]]@data, X[[2]]@data)
    }
    else{
      for (i in 2:(cpt.no)){
        Psis <- list()
        for (k in 1:length(dims)){
          Psis <- append(Psis,list(AR_model(dims[k],dens[i])$inv))
        }
        Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
        X[[i]] <- gen_sylvester_data(N[i] - N[i-1],dims,Omega_sqrt, ac_coeff = ac)
      }
      Psis <- list()
      for (k in 1:length(dims)){
        Psis <- append(Psis,list(AR_model(dims[k],dens[cpt.no+1])$inv))
      }
      Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
      X[[cpt.no+1]] <- gen_sylvester_data(size - N[cpt.no],dims,Omega_sqrt, ac_coeff = ac)
      z <- X[[1]]@data
      for(i in 2:(cpt.no+1)){
        z <- abind(z, X[[i]]@data)
      }
    }
  }
  z.tnsr = as.tensor(z)
  if (decomp == "cp"){
    cpz <- cp(z.tnsr, num_components = par_decom)
    return(list(data = cpz$U[[4]], norm = cpz$norm_percent))
  }
  if (decomp == "hosvd"){
    hosvd1 <- hosvd(z.tnsr, ranks = c(z.tnsr@modes[-4],par_decom))
    return(list(data = hosvd1$U[[4]], norm = hosvd1$fnorm_resid))
  }
  if (decomp == "tucker"){
    library(tensorr)
    tucker1 <- tucker(z.tnsr, ranks = c(1,10,10,dim(z)[4]))
    return(list(data = tucker1$U[[4]], norm = tucker1$norm_percent))
  }
}

# The following function builds tensors according to values given in the input arguments
# with precision matrices that undergo changes in the ER structure. We then employ
# a decomposition method chosen by the user and we return the decomposed data, as well as
# the percent of the Frobenius norm explained by the approximation.

tensor_creation_decomposition_ER <- function(dims = c(20,20,20), N = 100, size = 200,dens = 10, wmin = NULL, wmax = NULL, ac = 0, decomp = "cp", par_decom = NULL){
  K <- length(dims)
  cpt.no <- length(N)
  X <- list()
  if (length(wmin) != length(wmax)){stop("The lengths of the parameter vectors do not agree.")}
  if (length(wmin) != (cpt.no+1)){stop("The lengths of the parameter vectors do not agree.")}
  if(cpt.no == 0){
    Psis <- list()
    for (k in 1:length(dims)){
      Psis <- append(Psis,list(ER_model(dims[k], dens, wmin, wmax, corr=TRUE)$inv))
    }
    Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
    X[[1]] <- gen_sylvester_data(size,dims,Omega_sqrt, ac_coeff = ac)
    z = X[[1]]@data
  }
  else{
    Psis <- list()
    for (k in 1:length(dims)){
      Psis <- append(Psis,list(ER_model(dims[k], dens, wmin[1], wmax[1], corr=TRUE)$inv))
    }
    Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
    X[[1]] <- gen_sylvester_data(N[1],dims,Omega_sqrt, ac_coeff = ac)
    if(cpt.no == 1){
      Psis <- list()
      for (k in 1:length(dims)){
        Psis <- append(Psis,list(ER_model(dims[k], dens, wmin[2], wmax[2], corr=TRUE)$inv))
      }
      Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
      X[[2]] <- gen_sylvester_data(size - N[1],dims,Omega_sqrt, ac_coeff = ac)
      z = abind(X[[1]]@data, X[[2]]@data)
    }
    else{
      for (i in 2:(cpt.no)){
        Psis <- list()
        for (k in 1:length(dims)){
          Psis <- append(Psis,list(ER_model(dims[k], dens, wmin[i], wmax[i], corr=TRUE)$inv))
        }
        Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
        X[[i]] <- gen_sylvester_data(N[i] - N[i-1],dims,Omega_sqrt, ac_coeff = ac)
      }
      Psis <- list()
      for (k in 1:length(dims)){
        Psis <- append(Psis,list(ER_model(dims[k], dens, wmin[cpt.no + 1], wmax[cpt.no + 1], corr=TRUE)$inv))
      }
      Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
      X[[cpt.no+1]] <- gen_sylvester_data(size - N[cpt.no],dims,Omega_sqrt, ac_coeff = ac)
      z <- X[[1]]@data
      for(i in 2:(cpt.no+1)){
        z <- abind(z, X[[i]]@data)
      }
    }
  }
  z.tnsr = as.tensor(z)
  if (decomp == "cp"){
    cpz <- cp(z.tnsr, num_components = par_decom)
    return(list(data = cpz$U[[4]], norm = cpz$norm_percent))
  }
  if (decomp == "hosvd"){
    hosvd1 <- hosvd(z.tnsr, ranks = c(z.tnsr@modes[-4],par_decom))
    return(list(data = hosvd1$U[[4]], norm = hosvd1$fnorm_resid))
  }
  if (decomp == "tucker"){
    library(tensorr)
    tucker1 <- tucker(z.tnsr, ranks = c(1,10,10,dim(z)[4])) # Have to specify the ranks
    return(list(data = tucker1$U[[4]], norm = tucker1$norm_percent))
  }
}

# The following function builds tensors according to values given in the input arguments
# with precision matrices that undergo changes in the SB structure. We then employ
# a decomposition method chosen by the user and we return the decomposed data, as well as
# the percent of the Frobenius norm explained by the approximation.

tensor_creation_decomposition_SB <- function(dims = c(20,20,20), N = 100, size = 200, rho= NULL, num_subgraph = NULL, ac = 0, decomp = "cp", par_decom = NULL){
  K <- length(dims)
  cpt.no <- length(N)
  X <- list()
  if (length(rho) != length(num_subgraph)){stop("The lengths of the parameter vectors do not agree.")}
  if (length(rho) != (cpt.no+1)){stop("The lengths of the parameter vectors do not agree.")}
  if(cpt.no == 0){
    Psis <- list()
    for (k in 1:length(dims)){
      Psis <- append(Psis,list(SB_model(dims[k], rho, num_subgraph)$inv))
    }
    Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
    X[[1]] <- gen_sylvester_data(size,dims,Omega_sqrt, ac_coeff = ac)
    z = X[[1]]@data
  }
  else{
    Psis <- list()
    for (k in 1:length(dims)){
      Psis <- append(Psis,list(SB_model(dims[k], rho[1], num_subgraph[1])$inv))
    }
    Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
    X[[1]] <- gen_sylvester_data(N[1],dims,Omega_sqrt, ac_coeff = ac)
    if(cpt.no == 1){
      Psis <- list()
      for (k in 1:length(dims)){
        Psis <- append(Psis,list(SB_model(dims[k], rho[2], num_subgraph[2])$inv))
      }
      Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
      X[[2]] <- gen_sylvester_data(size - N[1],dims,Omega_sqrt, ac_coeff = ac)
      z = abind(X[[1]]@data, X[[2]]@data)
    }
    else{
      for (i in 2:(cpt.no)){
        Psis <- list()
        for (k in 1:length(dims)){
          Psis <- append(Psis,list(SB_model(dims[k], rho[i], num_subgraph[i])$inv))
        }
        Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
        X[[i]] <- gen_sylvester_data(N[i] - N[i-1],dims,Omega_sqrt, ac_coeff = ac)
      }
      Psis <- list()
      for (k in 1:length(dims)){
        Psis <- append(Psis,list(SB_model(dims[k], rho[cpt.no + 1], num_subgraph[cpt.no + 1])$inv))
      }
      Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
      X[[cpt.no+1]] <- gen_sylvester_data(size - N[cpt.no],dims,Omega_sqrt, ac_coeff = ac)
      z <- X[[1]]@data
      for(i in 2:(cpt.no+1)){
        z <- abind(z, X[[i]]@data)
      }
    }
  }
  z.tnsr = as.tensor(z)
  if (decomp == "cp"){
    cpz <- cp(z.tnsr, num_components = par_decom)
    return(list(data = cpz$U[[4]], norm = cpz$norm_percent))
  }
  if (decomp == "hosvd"){
    hosvd1 <- hosvd(z.tnsr, ranks = c(z.tnsr@modes[-4],par_decom))
    return(list(data = hosvd1$U[[4]], norm = hosvd1$fnorm_resid))
  }
  if (decomp == "tucker"){
    library(tensorr)
    tucker1 <- tucker(z.tnsr, ranks = c(1,10,10,dim(z)[4])) # Have to specify the ranks
    return(list(data = tucker1$U[[4]], norm = tucker1$norm_percent))
  }
}
