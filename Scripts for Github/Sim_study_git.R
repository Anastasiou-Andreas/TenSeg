#===============================================================
# Author: Andreas Anastasiou
# Project: Change-point detection on tensors

# Description: The code below carries out the simulation study as explained in
# the paper. First, the user needs to run all the functions in the main.R script

#===============================================================
rm(list = ls())   # Clear environment

library(rTensor)
library(ccid)
library(hdbinseg)
library(cpt.cov)

source(main.R)


#=============================================================
# Functions necessary for the Barnett & Onnela (2016) method.
#=============================================================
func_2_segmentation <- function(dataset, offset = 10,alpha_value){
  n <- nrow(dataset)
  tt = 0
  
  res_CPD = res_pvalue = final_result <- c()
  
  data_work = dataset
  
  res_CPD = 0
  while (res_CPD < 5){
    cp.null.out = CPdetection(as.matrix(data_work),offset = offset,
                              display = FALSE)
    
    res_CPD = cp.null.out$CPind
    res_pvalue = cp.null.out$pvalue
  }
  
  out2 <- data.frame (CPD_vector = c(0, res_CPD, nrow(data_work)),
                      pvalues = c(0, res_pvalue, 0))
  
  newps = res_pvalue
  
  while (any(newps < alpha_value)){
    cutpoints <- out2$CPD_vector
    N <- length (cutpoints)
    newps <- c()
    
    for (i in 1:(N-1))
    {
      try({
        
        if ((out2$pvalues[i]) < alpha_value & (out2$pvalues[i+1]) < alpha_value){
          
          new_data <- data_work[(cutpoints[i]+1):cutpoints[i+1], ]
          
          new_data.boot <- CPdetection(as.matrix(new_data), display = FALSE, offset = offset)
          
          
          new_ch <- new_data.boot$CPind + cutpoints[i]
          
          
          newp <- new_data.boot$pvalue
          
          newps <- c (newps, newp)
          
          out2 <- rbind (out2, c(new_ch, newp))
          
        }
      })
      
    }
    out2 <- out2[order(out2$CPD_vector), ]
    
  }
  
  out <- out2[out2$CPD_vector != 0 & out2$CPD_vector != n &
                out2$pvalues<alpha_value, ]
  
  final_result <- rbind(final_result, out)
  
  return(final_result)
  
}

#=====================================================================
# Below you can find all the functions that will be used for the main
# simulation study. Each function below corresponds to a method used in the
# study.
#=====================================================================

# Barnett & Onnela (2016)
rev.sim.small_BO <- function(x) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  
  print("Barnett")
  z <- func_2_segmentation(x, alpha_value = 0.05)
  
  
  BO <- new("cpt.est")
  if(length(z$CPD_vector) == 0){BO@cpt = 0
  BO@nocpt = 0}
  else {BO@cpt= z$CPD_vector
  BO@nocpt <- length(BO@cpt)}
  BO@time <- system.time(  z <- func_2_segmentation(x, alpha_value = 0.05))[[3]]
  
  
  list(BO = BO)
}

# Ryan & Killick (2023)
rev.sim.small_Ratio <- function(x) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  
  print("Ratio")
  result.fisher = bin.seg(x, c(0,nrow(x)), matrix.dist.test.stat, threshold=qnorm(1-.05/(500^2)), minseglen=10, c())
  z = result.fisher %>% bin.seg.to.cpt(qnorm(1-.05/(500^2)))
  Ratio <- new("cpt.est")
  if(length(z$cpts) == 0){Ratio@cpt = 0
  Ratio@nocpt = 0}
  else {Ratio@cpt= z$cpts
  Ratio@nocpt <- length(Ratio@cpt)}
  Ratio@time <- system.time(z <- bin.seg(x, c(0,nrow(x)), matrix.dist.test.stat, threshold=qnorm(1-.05/(500^2)), minseglen=10, c()) %>% bin.seg.to.cpt(qnorm(1-.05/(500^2))))[[3]]
  list(Ratio = Ratio)
}

# Galeano & Pena (2007)
rev.sim.small_GP <- function(x) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  
  print("GP")
  result.GP = bin.seg(x, c(0,nrow(x)), cpt.cov:::galeano.cusum.stat, threshold=qnorm(1-.05/(500^2)), minseglen=10, c())
  z = result.GP %>% bin.seg.to.cpt(log(nrow(x)))
  GP <- new("cpt.est")
  if(length(z$cpts) == 0){GP@cpt = 0
  GP@nocpt = 0}
  else {GP@cpt= z$cpts
  GP@nocpt <- length(GP@cpt)}
  GP@time <- system.time(z <- bin.seg(x, c(0,nrow(x)), cpt.cov:::galeano.cusum.stat, threshold=qnorm(1-.05/(500^2)), minseglen=10, c()) %>% bin.seg.to.cpt(log(nrow(x))))[[3]]
  list(GP = GP)
}

# Wang et. al (2021)
rev.sim.small_WYR <- function(x) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  
  print("WYR")
  WYR.thresh = wang.threshold(x)
  result.WYR = bin.seg(x, c(0,nrow(x)), wang.stat, threshold=wang.thresh, minseglen=10, c())
  z = result.WYR %>% bin.seg.to.cpt(wang.thresh)
  WYR <- new("cpt.est")
  if(length(z$cpts) == 0){WYR@cpt = 0
  WYR@nocpt = 0}
  else {WYR@cpt= z$cpts
  WYR@nocpt <- length(WYR@cpt)}
  WYR@time <- system.time(z <- bin.seg(x, c(0,nrow(x)), wang.stat, threshold=wang.threshold(x), minseglen=10, c()) %>% bin.seg.to.cpt(wang.thresh))[[3]]
  list(WYR = WYR)
}

# Cho & Fryzlewicz (2015)
rev.sim.small_sbs <- function(x) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  print("SBS")
  z <- sbs.alg(t(x),cp.type = 2)
  cpt.z = z$ecp
  
  sbs <- new("cpt.est")
  if(length(cpt.z) == 0){sbs@cpt = 0
  sbs@nocpt = 0}
  else{sbs@cpt <- cpt.z
  sbs@nocpt <- length(cpt.z)}
  sbs@time <- system.time(sbs.alg(t(x),cp.type = 2))[[3]]
  
  list(sbs=sbs)
}

# The TenSeg method
rev.sim.small_TenSeg <- function(x, points = 10, thr_max_ic = 2.1) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  print("ID_Linf_IC")
  z <- ccid::detect.ic(x, approach = "infinity", pointsgen = points, th_max = thr_max_ic, alpha_gen = 0.1)
  cptLinf_IC = z$changepoints
  IDLinf_IC <- new("cpt.est")
  if(is.na(cptLinf_IC)){IDLinf_IC@cpt = 0
  IDLinf_IC@nocpt = 0}
  else {IDLinf_IC@cpt= cptLinf_IC
  IDLinf_IC@nocpt <- length(cptLinf_IC)}
  IDLinf_IC@time <- system.time(ccid::detect.ic(x, approach = "infinity", pointsgen = points, th_max = thr_max_ic, alpha_gen = 0.1))[[3]]
  
  list(IDLinf_IC = IDLinf_IC)
}

#===========================================================================
# Code for the simulation study when the input is an already decomposed tensor
# and the precision matrices could undergo changes in either AR, ER, or SB
# structure.
#===========================================================================
simulation_study_competitors <- function(decomposed_tensor = NULL,decomp_sim = "cp", true.cpt=NULL, m = 100, seed = NULL) {
  setClass("est.eval", representation(cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  sbs <- new("est.eval")
  BO <- new("est.eval")
  Ratio <- new("est.eval")
  GP <- new("est.eval")
  WYR <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(c(diff(true.cpt), dim(decomposed_tensor[[1]]$data)[1]))
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- decomposed_tensor[[i]]$data
    ## SBS
    est <- rev.sim.small_sbs(x)
    
    sbs@dnc[i] <- est$sbs@nocpt - no.of.cpt
    sbs@cpt[[i]] <- est$sbs@cpt
    sbs@diff <- abs(matrix(est$sbs@cpt,nrow=no.of.cpt,ncol=length(est$sbs@cpt),byr=TRUE)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$sbs@cpt),byr=F))
    sbs@dh[i] <- max(apply(sbs@diff,1,min),apply(sbs@diff,2,min))/ns
    sbs@time[i] <- est$sbs@time
    
    # BO
    est <- rev.sim.small_BO(x)
    
    BO@dnc[i] <- est$BO@nocpt - no.of.cpt
    BO@cpt[[i]] <- est$BO@cpt
    BO@diff <- abs(matrix(est$BO@cpt,nrow=no.of.cpt,ncol=length(est$BO@cpt),byr=TRUE)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$BO@cpt),byr=F))
    BO@dh[i] <- max(apply(BO@diff,1,min),apply(BO@diff,2,min))/ns
    BO@time[i] <- est$BO@time
    
    # Ratio
    est <- rev.sim.small_Ratio(x)
    
    Ratio@dnc[i] <- est$Ratio@nocpt - no.of.cpt
    Ratio@cpt[[i]] <- est$Ratio@cpt
    Ratio@diff <- abs(matrix(est$Ratio@cpt,nrow=no.of.cpt,ncol=length(est$Ratio@cpt),byr=TRUE)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Ratio@cpt),byr=F))
    Ratio@dh[i] <- max(apply(Ratio@diff,1,min),apply(Ratio@diff,2,min))/ns
    Ratio@time[i] <- est$Ratio@time
    
    # GP
    est <- rev.sim.small_GP(x)
    
    GP@dnc[i] <- est$GP@nocpt - no.of.cpt
    GP@cpt[[i]] <- est$GP@cpt
    GP@diff <- abs(matrix(est$GP@cpt,nrow=no.of.cpt,ncol=length(est$GP@cpt),byr=TRUE)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$GP@cpt),byr=F))
    GP@dh[i] <- max(apply(GP@diff,1,min),apply(GP@diff,2,min))/ns
    GP@time[i] <- est$GP@time
    
    # WYR
    est <- rev.sim.small_WYR(x)
    
    WYR@dnc[i] <- est$WYR@nocpt - no.of.cpt
    WYR@cpt[[i]] <- est$WYR@cpt
    WYR@diff <- abs(matrix(est$WYR@cpt,nrow=no.of.cpt,ncol=length(est$WYR@cpt),byr=TRUE)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$WYR@cpt),byr=F))
    WYR@dh[i] <- max(apply(WYR@diff,1,min),apply(WYR@diff,2,min))/ns
    WYR@time[i] <- est$WYR@time
    
    gc()
  }
  list(sbs = sbs, BO = BO, Ratio = Ratio, GP = GP, WYR = WYR)
}


simulation_study_TenSeg <- function(decomposed_tensor = NULL, decomp_sim = "cp",points_sim = 10, true.cpt=NULL, m = 100, seed = NULL, size_sim) {
  setClass("est.eval", representation(cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  IDLinf_IC <- new("est.eval")
  
  no.of.cpt <- length(true.cpt)
  ns <- max(c(diff(true.cpt), dim(decomposed_tensor[[1]]$data)[1]))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- decomposed_tensor[[i]]$data
    
    ## ID
    est <- rev.sim.small_TenSeg(x, points = points_sim)
    
    IDLinf_IC@dnc[i] <- est$IDLinf_IC@nocpt - no.of.cpt
    IDLinf_IC@cpt[[i]] <- est$IDLinf_IC@cpt
    IDLinf_IC@diff <- abs(matrix(est$IDLinf_IC@cpt,nrow=no.of.cpt,ncol=length(est$IDLinf_IC@cpt),byr=TRUE)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$IDLinf_IC@cpt),byr=F))
    IDLinf_IC@dh[i] <- max(apply(IDLinf_IC@diff,1,min),apply(IDLinf_IC@diff,2,min))/ns
    IDLinf_IC@time[i] <- est$IDLinf_IC@time
    
    gc()
  }
  list(IDLinf_IC = IDLinf_IC)
}

#==========================================================================
# Creation of all the decomposed data used in the simulation study. We create
# all the combinations explained in the paper, where the tensor has dimensionality
# 20 x 20 x 20 x T with various different numbers and locations for the change-points.
# as explained below.
# 
# In order to run the code you need relevant functions from the main.R script
#==========================================================================

# AR CHANGES
#
# The scenario where either the CP or HOSVD decompositions are employed and
# the precision matrices are generated from an AR1($\rho$) random graph model.
# We define a seed for reproducibility purposes.
#
# No change-points and 20 components are kept after the decomposition 
seed.temp <- 1
Tensors_0_cpt_20_components <- list()
Tensors_0_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_0_cpt_20_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, par_decom = 20)
  Tensors_0_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, decomp = "hosvd", par_decom = 20)
}

# No change-points and 10 components are kept after the decomposition 
Tensors_0_cpt_10_components <- list()
Tensors_0_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_0_cpt_10_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, par_decom = 10)
  Tensors_0_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, decomp = "hosvd", par_decom = 10)
}

# No change-points and 5 components are kept after the decomposition 
Tensors_0_cpt_5_components <- list()
Tensors_0_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_0_cpt_5_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, par_decom = 5)
  Tensors_0_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, decomp = "hosvd", par_decom = 5)
}

# One change-point and 20 components are kept after the decomposition
Tensors_1_cpt_20_components <- list()
Tensors_1_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_1_cpt_20_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), par_decom = 20)
  Tensors_1_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), decomp = "hosvd", par_decom = 20)
}

# One change-point and 10 components are kept after the decomposition
Tensors_1_cpt_10_components <- list()
Tensors_1_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_1_cpt_10_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), par_decom = 10)
  Tensors_1_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), decomp = "hosvd", par_decom = 10)
}

# One change-point and 5 components are kept after the decomposition
Tensors_1_cpt_5_components <- list()
Tensors_1_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_1_cpt_5_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), par_decom = 5)
  Tensors_1_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), decomp = "hosvd", par_decom = 5)
}

# Four change-points and 20 components are kept after the decomposition
Tensors_4_cpt_20_components <- list()
Tensors_4_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_4_cpt_20_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), par_decom = 20)
  Tensors_4_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), decomp = "hosvd", par_decom = 20)
}

# Four change-points and 10 components are kept after the decomposition
Tensors_4_cpt_10_components <- list()
Tensors_4_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_4_cpt_10_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), par_decom = 10)
  Tensors_4_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), decomp = "hosvd", par_decom = 10)
}

# Four change-points and 5 components are kept after the decomposition
Tensors_4_cpt_5_components <- list()
Tensors_4_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_4_cpt_5_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), par_decom = 5)
  Tensors_4_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), decomp = "hosvd", par_decom = 5)
}

# Ten change-points and 20 components are kept after the decomposition
Tensors_10_cpt_20_components <- list()
Tensors_10_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_10_cpt_20_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), par_decom = 20)
  Tensors_10_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), decomp = "hosvd", par_decom = 20)
}

# Ten change-points and 10 components are kept after the decomposition
Tensors_10_cpt_10_components <- list()
Tensors_10_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_10_cpt_10_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), par_decom = 10)
  Tensors_10_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), decomp = "hosvd", par_decom = 10)  
}

# Ten change-points and 5 components are kept after the decomposition
Tensors_10_cpt_5_components <- list()
Tensors_10_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  Tensors_10_cpt_5_components[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), par_decom = 5)
  Tensors_10_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), decomp = "hosvd", par_decom = 5)
}


#==============================================================================
# After the creation of the decomposed data, we now run the simulation study
# functions in order to investigate the performance of all the methods in all
# scenarios tested in the paper related to AR changes.
#==============================================================================

# Results for Tenseg related to changes in the AR structure when the decomposition
# used is the CP one
SIM0_cp_TenSeg_20_components <- simulation_study_AR_TenSeg(tensor = Tensors_0_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SIM1_cp_TenSeg_20_components <- simulation_study_AR_TenSeg(tensor = Tensors_1_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SIM4_cp_TenSeg_20_components <- simulation_study_AR_TenSeg(tensor = Tensors_4_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM10_cp_TenSeg_20_components <- simulation_study_AR_TenSeg(tensor = Tensors_10_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SIM0_cp_TenSeg_10_components <- simulation_study_AR_TenSeg(tensor = Tensors_0_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SIM1_cp_TenSeg_10_components <- simulation_study_AR_TenSeg(tensor = Tensors_1_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SIM4_cp_TenSeg_10_components <- simulation_study_AR_TenSeg(tensor = Tensors_4_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM10_cp_TenSeg_10_components <- simulation_study_AR_TenSeg(tensor = Tensors_10_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), size = 660, m = 100, seed = seed.temp)

SIM0_cp_TenSeg_5_components <- simulation_study_AR_TenSeg(tensor = Tensors_0_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SIM1_cp_TenSeg_5_components <- simulation_study_AR_TenSeg(tensor = Tensors_1_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SIM4_cp_TenSeg_5_components <- simulation_study_AR_TenSeg(tensor = Tensors_4_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM10_cp_TenSeg_5_components <- simulation_study_AR_TenSeg(tensor = Tensors_10_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), size = 660, m = 100, seed = seed.temp)

# Results for the competitors related to changes in the AR structure when the decomposition
# used is the CP one
SIM0_cp_comp_20_components <- simulation_study_AR_competitors(tensor = Tensors_0_cpt_20_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
SIM1_cp_comp_20_components <- simulation_study_AR_competitors(tensor = Tensors_1_cpt_20_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
SIM4_cp_comp_20_components <- simulation_study_AR_competitors(tensor = Tensors_4_cpt_20_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM10_cp_comp_20_components <- simulation_study_AR_competitors(tensor = Tensors_10_cpt_20_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SIM0_cp_comp_10_components <- simulation_study_AR_competitors(tensor = Tensors_0_cpt_10_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
SIM1_cp_comp_10_components <- simulation_study_AR_competitors(tensor = Tensors_1_cpt_10_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
SIM4_cp_comp_10_components <- simulation_study_AR_competitors(tensor = Tensors_4_cpt_10_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM10_cp_comp_10_components <- simulation_study_AR_competitors(tensor = Tensors_10_cpt_10_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SIM0_cp_comp_5_components <- simulation_study_AR_competitors(tensor = Tensors_0_cpt_5_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
SIM1_cp_comp_5_components <- simulation_study_AR_competitors(tensor = Tensors_1_cpt_5_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
SIM4_cp_comp_5_components <- simulation_study_AR_competitors(tensor = Tensors_4_cpt_5_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM10_cp_comp_5_components <- simulation_study_AR_competitors(tensor = Tensors_10_cpt_5_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

# Results for Tenseg related to changes in the AR structure when the decomposition
# used is the HOSVD one
SIM0_cp_TenSeg_hosvd_20 <- simulation_study_AR_TenSeg(tensor = Tensors_0_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SIM0_cp_TenSeg_hosvd_10 <- simulation_study_AR_TenSeg(tensor = Tensors_0_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SIM0_cp_TenSeg_hosvd_5 <- simulation_study_AR_TenSeg(tensor = Tensors_0_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)

SIM1_cp_TenSeg_hosvd_20 <- simulation_study_AR_TenSeg(tensor = Tensors_1_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SIM1_cp_TenSeg_hosvd_10 <- simulation_study_AR_TenSeg(tensor = Tensors_1_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SIM1_cp_TenSeg_hosvd_5 <- simulation_study_AR_TenSeg(tensor = Tensors_1_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)

SIM4_cp_TenSeg_hosvd_20 <- simulation_study_AR_TenSeg(tensor = Tensors_4_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM4_cp_TenSeg_hosvd_10 <- simulation_study_AR_TenSeg(tensor = Tensors_4_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SIM4_cp_TenSeg_hosvd_5 <- simulation_study_AR_TenSeg(tensor = Tensors_4_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)

SIM10_cp_TenSeg_hosvd_20 <- simulation_study_AR_TenSeg(tensor = Tensors_10_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
SIM10_cp_TenSeg_hosvd_10 <- simulation_study_AR_TenSeg(tensor = Tensors_10_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
SIM10_cp_TenSeg_hosvd_5 <- simulation_study_AR_TenSeg(tensor = Tensors_10_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

#=================
# ER CHANGES
#
# The scenario where either the CP or HOSVD decompositions are employed and
# the precision matrices are generated from an Erdos-Renyi random graph (ER).
# We define a seed for reproducibility purposes.
#
# No change-points and 20 components are kept after the decomposition 

seed.temp <- 1
ER_Tensors_0_cpt_20_components <- list()
ER_Tensors_0_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_0_cpt_20_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "cp", par_decom = 20)
  ER_Tensors_0_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "hosvd", par_decom = 20)
}

# No change-points and 10 components are kept after the decomposition
ER_Tensors_0_cpt_10_components <- list()
ER_Tensors_0_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_0_cpt_10_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "cp", par_decom = 10)
  ER_Tensors_0_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "hosvd", par_decom = 10)
}

# No change-points and 5 components are kept after the decomposition
ER_Tensors_0_cpt_5_components <- list()
ER_Tensors_0_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_0_cpt_5_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "cp", par_decom = 5)
  ER_Tensors_0_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "hosvd", par_decom = 5)
}

# One change-point and 20 components are kept after the decomposition
ER_Tensors_1_cpt_20_components <- list()
ER_Tensors_1_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_1_cpt_20_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "cp", par_decom = 20)
  ER_Tensors_1_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "hosvd", par_decom = 20)
}

# One change-point and 10 components are kept after the decomposition
ER_Tensors_1_cpt_10_components <- list()
ER_Tensors_1_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_1_cpt_10_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "cp", par_decom = 10)
  ER_Tensors_1_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "hosvd", par_decom = 10)
}

# One change-point and 5 components are kept after the decomposition
ER_Tensors_1_cpt_5_components <- list()
ER_Tensors_1_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_1_cpt_5_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "cp", par_decom = 5)
  ER_Tensors_1_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "hosvd", par_decom = 5)
}

# Four change-points and 20 components are kept after the decomposition
ER_Tensors_4_cpt_20_components <- list()
ER_Tensors_4_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_4_cpt_20_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9), decomp = "cp", par_decom = 20)
  ER_Tensors_4_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9), decomp = "hosvd", par_decom = 20)
}

# Four change-points and 10 components are kept after the decomposition
ER_Tensors_4_cpt_10_components <- list()
ER_Tensors_4_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_4_cpt_10_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9), decomp = "cp", par_decom = 10)
  ER_Tensors_4_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9), decomp = "hosvd", par_decom = 10)
}

# Four change-points and 5 components are kept after the decomposition
ER_Tensors_4_cpt_5_components <- list()
ER_Tensors_4_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_4_cpt_5_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9), decomp = "cp", par_decom = 5)
  ER_Tensors_4_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9), decomp = "hosvd", par_decom = 5)
}

# Ten change-points and 20 components are kept after the decomposition
ER_Tensors_10_cpt_20_components <- list()
ER_Tensors_10_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_10_cpt_20_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9), decomp = "cp", par_decom = 20)
  ER_Tensors_10_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9), decomp = "hosvd", par_decom = 20)
}

# Ten change-points and 10 components are kept after the decomposition
ER_Tensors_10_cpt_10_components <- list()
ER_Tensors_10_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_10_cpt_10_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9), decomp = "cp", par_decom = 10)
  ER_Tensors_10_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9), decomp = "hosvd", par_decom = 10)
}

# Ten change-points and 5 components are kept after the decomposition
ER_Tensors_10_cpt_5_components <- list()
ER_Tensors_10_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  ER_Tensors_10_cpt_5_components[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9), decomp = "cp", par_decom = 5)
  ER_Tensors_10_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8,0.05,0.8), wmax = c(0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9), decomp = "hosvd", par_decom = 5)
}



#==============================================================================
# After the creation of the decomposed data, we now run the simulation study
# functions in order to investigate the performance of all the methods in all
# scenarios tested in the paper related to ER changes.
#==============================================================================

# Results for Tenseg related to changes in the ER structure when the decomposition
# used is the CP one

ER_SIM0_cp_TenSeg_20_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_0_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM1_cp_TenSeg_20_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_1_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM4_cp_TenSeg_20_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_4_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM10_cp_TenSeg_20_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_10_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

ER_SIM0_cp_TenSeg_10_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_0_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM1_cp_TenSeg_10_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_1_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM4_cp_TenSeg_10_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_4_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM10_cp_TenSeg_10_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_10_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

ER_SIM0_cp_TenSeg_5_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_0_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM1_cp_TenSeg_5_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_1_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM4_cp_TenSeg_5_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_4_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM10_cp_TenSeg_5_components <- simulation_study_ER_TenSeg(tensor = ER_Tensors_10_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)


# Results for the competitors related to changes in the ER structure when the decomposition
# used is the CP one
ER_SIM0_cp_comp_20_components <- simulation_study_ER_competitors(tensor = ER_Tensors_0_cpt_20_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM1_cp_comp_20_components <- simulation_study_ER_competitors(tensor = ER_Tensors_1_cpt_20_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM4_cp_comp_20_components <- simulation_study_ER_competitors(tensor = ER_Tensors_4_cpt_20_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM10_cp_comp_20_components <- simulation_study_ER_competitors(tensor = ER_Tensors_10_cpt_20_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

ER_SIM0_cp_comp_10_components <- simulation_study_ER_competitors(tensor = ER_Tensors_0_cpt_10_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM1_cp_comp_10_components <- simulation_study_ER_competitors(tensor = ER_Tensors_1_cpt_10_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM4_cp_comp_10_components <- simulation_study_ER_competitors(tensor = ER_Tensors_4_cpt_10_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM10_cp_comp_10_components <- simulation_study_ER_competitors(tensor = ER_Tensors_10_cpt_10_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

ER_SIM0_cp_comp_5_components <- simulation_study_ER_competitors(tensor = ER_Tensors_0_cpt_5_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM1_cp_comp_5_components <- simulation_study_ER_competitors(tensor = ER_Tensors_1_cpt_5_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM4_cp_comp_5_components <- simulation_study_ER_competitors(tensor = ER_Tensors_4_cpt_5_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM10_cp_comp_5_components <- simulation_study_ER_competitors(tensor = ER_Tensors_10_cpt_5_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)


# Results for Tenseg related to changes in the ER structure when the decomposition
# used is the HOSVD one

ER_SIM0_cp_TenSeg_hosvd_20 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_0_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM0_cp_TenSeg_hosvd_10 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_0_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
ER_SIM0_cp_TenSeg_hosvd_5 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_0_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)

ER_SIM1_cp_TenSeg_hosvd_20 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_1_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM1_cp_TenSeg_hosvd_10 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_1_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
ER_SIM1_cp_TenSeg_hosvd_5 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_1_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)

ER_SIM4_cp_TenSeg_hosvd_20 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_4_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM4_cp_TenSeg_hosvd_10 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_4_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
ER_SIM4_cp_TenSeg_hosvd_5 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_4_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)

ER_SIM10_cp_TenSeg_hosvd_20 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_10_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
ER_SIM10_cp_TenSeg_hosvd_10 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_10_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
ER_SIM10_cp_TenSeg_hosvd_5 <- simulation_study_ER_TenSeg(tensor = ER_Tensors_10_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)


#=================
# SB CHANGES
#
# The scenario where either the CP or HOSVD decompositions are employed and
# the precision matrices are generated from a Star-Block (SB) random graph model.
# We define a seed for reproducibility purposes.
#
# No change-points and 20 components are kept after the decomposition 
seed.temp <- 1
SB_Tensors_0_cpt_20_components <- list()
SB_Tensors_0_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_0_cpt_20_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "cp", par_decom = 20)
  SB_Tensors_0_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "hosvd", par_decom = 20)
}

# No change-points and 10 components are kept after the decomposition
SB_Tensors_0_cpt_10_components <- list()
SB_Tensors_0_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_0_cpt_10_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "cp", par_decom = 10)
  SB_Tensors_0_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "hosvd", par_decom = 10)
}

# No change-points and 5 components are kept after the decomposition
SB_Tensors_0_cpt_5_components <- list()
SB_Tensors_0_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_0_cpt_5_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "cp", par_decom = 5)
  SB_Tensors_0_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "hosvd", par_decom = 5)
}

# One change-point and 20 components are kept after the decomposition
SB_Tensors_1_cpt_20_components <- list()
SB_Tensors_1_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_1_cpt_20_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "cp", par_decom = 20)
  SB_Tensors_1_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "hosvd", par_decom = 20)
}

# One change-point and 10 components are kept after the decomposition
SB_Tensors_1_cpt_10_components <- list()
SB_Tensors_1_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_1_cpt_10_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "cp", par_decom = 10)
  SB_Tensors_1_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "hosvd", par_decom = 10)
}

# One change-point and 5 components are kept after the decomposition
SB_Tensors_1_cpt_5_components <- list()
SB_Tensors_1_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_1_cpt_5_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "cp", par_decom = 5)
  SB_Tensors_1_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "hosvd", par_decom = 5)
}

# Four change-points and 20 components are kept after the decomposition
SB_Tensors_4_cpt_20_components <- list()
SB_Tensors_4_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_4_cpt_20_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "cp", par_decom = 20)
  SB_Tensors_4_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "hosvd", par_decom = 20)
}

# Four change-points and 10 components are kept after the decomposition
SB_Tensors_4_cpt_10_components <- list()
SB_Tensors_4_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_4_cpt_10_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "cp", par_decom = 10)
  SB_Tensors_4_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "hosvd", par_decom = 10)
}

# Four change-points and 5 components are kept after the decomposition
SB_Tensors_4_cpt_5_components <- list()
SB_Tensors_4_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_4_cpt_5_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "cp", par_decom = 5)
  SB_Tensors_4_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "hosvd", par_decom = 5)
}

# Ten change-points and 20 components are kept after the decomposition
SB_Tensors_10_cpt_20_components <- list()
SB_Tensors_10_cpt_hosvd_20 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_10_cpt_20_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "cp", par_decom = 20)
  SB_Tensors_10_cpt_hosvd_20[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "hosvd", par_decom = 20)
}

# Ten change-points and 10 components are kept after the decomposition
SB_Tensors_10_cpt_10_components <- list()
SB_Tensors_10_cpt_hosvd_10 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_10_cpt_10_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "cp", par_decom = 10)
  SB_Tensors_10_cpt_hosvd_10[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "hosvd", par_decom = 10)
}

# Ten change-points and 5 components are kept after the decomposition
SB_Tensors_10_cpt_5_components <- list()
SB_Tensors_10_cpt_hosvd_5 <- list()
set.seed(seed.temp)
for (i in 1:100){
  print(i)
  SB_Tensors_10_cpt_5_components[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "cp", par_decom = 5)
  SB_Tensors_10_cpt_hosvd_5[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "hosvd", par_decom = 5)
}


#==============================================================================
# After the creation of the decomposed data, we now run the simulation study
# functions in order to investigate the performance of all the methods in all
# scenarios tested in the paper related to SB changes.
#==============================================================================

# Results for Tenseg related to changes in the SB structure when the decomposition
# used is the CP one
SB_SIM0_cp_TenSeg_20_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_0_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM1_cp_TenSeg_20_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_1_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM4_cp_TenSeg_20_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_4_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM10_cp_TenSeg_20_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_10_cpt_20_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SB_SIM0_cp_TenSeg_10_components <- sim+ulation_study_SB_TenSeg(tensor = SB_Tensors_0_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM1_cp_TenSeg_10_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_1_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM4_cp_TenSeg_10_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_4_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM10_cp_TenSeg_10_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_10_cpt_10_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SB_SIM0_cp_TenSeg_5_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_0_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM1_cp_TenSeg_5_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_1_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM4_cp_TenSeg_5_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_4_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM10_cp_TenSeg_5_components <- simulation_study_SB_TenSeg(tensor = SB_Tensors_10_cpt_5_components, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)


# Results for the competitors related to changes in the SB structure when the decomposition
# used is the CP one
SB_SIM0_cp_comp_20_components <- simulation_study_SB_competitors(tensor = SB_Tensors_0_cpt_20_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM1_cp_comp_20_components <- simulation_study_SB_competitors(tensor = SB_Tensors_1_cpt_20_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM4_cp_comp_20_components <- simulation_study_SB_competitors(tensor = SB_Tensors_4_cpt_20_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM10_cp_comp_20_components <- simulation_study_SB_competitors(tensor = SB_Tensors_10_cpt_20_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SB_SIM0_cp_comp_10_components <- simulation_study_SB_competitors(tensor = SB_Tensors_0_cpt_10_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM1_cp_comp_10_components <- simulation_study_SB_competitors(tensor = SB_Tensors_1_cpt_10_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM4_cp_comp_10_components <- simulation_study_SB_competitors(tensor = SB_Tensors_4_cpt_10_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM10_cp_comp_10_components <- simulation_study_SB_competitors(tensor = SB_Tensors_10_cpt_10_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)

SB_SIM0_cp_comp_5_components <- simulation_study_SB_competitors(tensor = SB_Tensors_0_cpt_5_components, decomp_sim = "cp", true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM1_cp_comp_5_components <- simulation_study_SB_competitors(tensor = SB_Tensors_1_cpt_5_components, decomp_sim = "cp", true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM4_cp_comp_5_components <- simulation_study_SB_competitors(tensor = SB_Tensors_4_cpt_5_components, decomp_sim = "cp", true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM10_cp_comp_5_components <- simulation_study_SB_competitors(tensor = SB_Tensors_10_cpt_5_components, decomp_sim = "cp", true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)


# Results for Tenseg related to changes in the SB structure when the decomposition
# used is the HOSVD one
SB_SIM0_cp_TenSeg_hosvd_20 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_0_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM0_cp_TenSeg_hosvd_10 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_0_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)
SB_SIM0_cp_TenSeg_hosvd_5 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_0_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = integer(0), m = 100, seed = seed.temp)

SB_SIM1_cp_TenSeg_hosvd_20 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_1_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM1_cp_TenSeg_hosvd_10 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_1_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)
SB_SIM1_cp_TenSeg_hosvd_5 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_1_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = 100, m = 100, seed = seed.temp)

SB_SIM4_cp_TenSeg_hosvd_20 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_4_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM4_cp_TenSeg_hosvd_10 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_4_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)
SB_SIM4_cp_TenSeg_hosvd_5 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_4_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(100,150,200,250), m = 100, seed = seed.temp)

SB_SIM10_cp_TenSeg_hosvd_20 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_10_cpt_hosvd_20, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
SB_SIM10_cp_TenSeg_hosvd_10 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_10_cpt_hosvd_10, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
SB_SIM10_cp_TenSeg_hosvd_5 <- simulation_study_SB_TenSeg(tensor = SB_Tensors_10_cpt_hosvd_5, decomp_sim = "hosvd", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 100, seed = seed.temp)
