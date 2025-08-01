#===============================================================
# Author: Andreas Anastasiou
# Project: Tensor time series change-point detection in cryptocurrency network data

# Description: The code below carries out the simulation study related to the
# proposed TenSeg method when the tensor data exhibit temporal dependence.
# The addition of temporal dependence makes the identification of changes far
# more challenging. More specifically, we work under the scenario where the noise
# in the simulated data is generated from an AR(1) process with a correlation
# coefficient set at 0.7. We also present results for high serial correlation (coefficient
# equal to 0.9). First, the user needs to run all the functions in the main.R script
#===============================================================



#==========================================================================
# Creation of all the decomposed data used in the simulation study. We create
# all the combinations explained in the paper, where the tensor has dimensionality
# 20 x 20 x 20 x T with various different numbers and locations for the change-points.
# as explained below.
#==========================================================================

# AR CHANGES
#
# The scenario where the CP decomposition is employed and the precision
# matrices are generated from an AR1($\rho$) random graph model. As
# mentioned, the noise is from an AR(1) process instead of white noise.
# We define a seed for reproducibility purposes.
#
# No change-points and 20 components are kept after the decomposition 
seed.temp <- 1
Tensors_0_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_0_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, ac = 0.7, par_decom = 20)
}

# No change-points and 10 components are kept after the decomposition 
Tensors_0_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_0_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, ac = 0.7,par_decom = 10)
}

# No change-points and 5 components are kept after the decomposition 
Tensors_0_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_0_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = integer(0), size = 200, dens = 0.2, ac = 0.7,par_decom = 5)
}

# One change-point and 20 components are kept after the decomposition
Tensors_1_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_1_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), ac = 0.7,par_decom = 20)
}

# One change-point and 10 components are kept after the decomposition
Tensors_1_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_1_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), ac = 0.7, par_decom = 10)
}

# One change-point and 5 components are kept after the decomposition
Tensors_1_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_1_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = 100, size = 200, dens = c(0.2,0.8), ac = 0.7, par_decom = 5)
}

# Four change-points and 20 components are kept after the decomposition
Tensors_4_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_4_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), ac = 0.7, size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), par_decom = 20)
}

# Four change-points and 10 components are kept after the decomposition
Tensors_4_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_4_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), ac = 0.7,size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), par_decom = 10)
}

# Four change-points and 5 components are kept after the decomposition
Tensors_4_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_4_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(100,150,200,250), ac = 0.7,size = 300, dens = c(0.2,0.8,0.2,0.8,0.2), par_decom = 5)
}

# Ten change-points and 20 components are kept after the decomposition
Tensors_10_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_10_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), ac = 0.7,size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), par_decom = 20)
}

# Ten change-points and 10 components are kept after the decomposition
Tensors_10_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_10_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), ac = 0.7,size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), par_decom = 10)
}

# Ten change-points and 5 components are kept after the decomposition
Tensors_10_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  Tensors_10_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_AR(dims = c(20,20,20), N = c(60,120,180,240,300,360,420,480,540,600), ac = 0.7,size = 660, dens = c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), par_decom = 5)
}


## FOR AR changes
SIM0_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_0_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
SIM1_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_1_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
SIM4_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_4_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
SIM10_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_10_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)

SIM0_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_0_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
SIM1_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_1_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
SIM4_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_4_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
SIM10_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_10_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), size = 660, m = 500, seed = seed.temp)

SIM0_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_0_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
SIM1_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_1_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
SIM4_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_4_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
SIM10_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = Tensors_10_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), size = 660, m = 500, seed = seed.temp)

# ER CHANGES
#
# The scenario where the CP decomposition is employed and the precision
# matrices are generated from an Erdos-Renyi (ER) random graph model. As
# mentioned, the noise is from an AR(1) process instead of white noise.
# We define a seed for reproducibility purposes.
#
# No change-points and 20 components are kept after the decomposition 
ER_Tensors_0_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_0_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "cp", par_decom = 20, ac = 0.7)
}

# No change-points and 10 components are kept after the decomposition 
ER_Tensors_0_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_0_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "cp", par_decom = 10, ac = 0.7)
}

# No change-points and 5 components are kept after the decomposition 
ER_Tensors_0_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_0_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = integer(0), size = 200, dens = 20, wmin = 0.7, wmax = 0.9, decomp = "cp", par_decom = 5, ac = 0.7)
}

# One change-point and 20 components are kept after the decomposition 
ER_Tensors_1_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_1_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "cp", par_decom = 20, ac = 0.7)
}

# One change-point and 10 components are kept after the decomposition 
ER_Tensors_1_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_1_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "cp", par_decom = 10, ac = 0.7)
}

# One change-point and 5 components are kept after the decomposition 
ER_Tensors_1_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_1_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = 100, size = 200, dens = 20, wmin = c(0.7, 0.1), wmax = c(0.9, 0.2), decomp = "cp", par_decom = 5, ac = 0.7)
}

# Four change-points and 20 components are kept after the decomposition
ER_Tensors_4_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_4_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.7,0.1,0.7,0.1,0.7), wmax = c(0.9,0.2,0.9,0.2,0.9), decomp = "cp", par_decom = 20, ac = 0.7)
}

# Four change-points and 10 components are kept after the decomposition
ER_Tensors_4_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_4_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.7,0.1,0.7,0.1,0.7), wmax = c(0.9,0.2,0.9,0.2,0.9), decomp = "cp", par_decom = 10, ac =0.7)
}

# Four change-points and 5 components are kept after the decomposition
ER_Tensors_4_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_4_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, dens = 20, wmin = c(0.7,0.1,0.7,0.1,0.7), wmax = c(0.9,0.2,0.9,0.2,0.9), decomp = "cp", par_decom = 5, ac=0.7)
}

# Ten change-points and 20 components are kept after the decomposition
ER_Tensors_10_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_10_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.7,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.7), wmax = c(0.9,0.2,0.9,0.2,0.9,0.2,0.9,0.2,0.9,0.2,0.9), decomp = "cp", par_decom = 20, ac=0.7)
}

# Ten change-points and 10 components are kept after the decomposition
ER_Tensors_10_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_10_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.7,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.7), wmax = c(0.9,0.2,0.9,0.2,0.9,0.2,0.9,0.2,0.9,0.2,0.9), decomp = "cp", par_decom = 10, ac=0.7)
}

# Ten change-points and 5 components are kept after the decomposition
ER_Tensors_10_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  ER_Tensors_10_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_ER(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, dens = 20, wmin = c(0.7,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.7), wmax = c(0.9,0.2,0.9,0.2,0.9,0.2,0.9,0.2,0.9,0.2,0.9), decomp = "cp", par_decom = 5, ac=0.7)
}

ER_SIM0_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_0_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
ER_SIM1_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_1_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
ER_SIM4_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_4_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
ER_SIM10_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_10_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)

ER_SIM0_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_0_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
ER_SIM1_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_1_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
ER_SIM4_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_4_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
ER_SIM10_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_10_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)

ER_SIM0_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_0_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
ER_SIM1_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_1_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
ER_SIM4_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_4_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
ER_SIM10_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = ER_Tensors_10_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)

# SB CHANGES
#
# The scenario where the CP decomposition is employed and the precision
# matrices are generated from a Star-Block (SB) random graph model. As
# mentioned, the noise is from an AR(1) process instead of white noise.
# We define a seed for reproducibility purposes.
#
# No change-points and 20 components are kept after the decomposition 
SB_Tensors_0_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_0_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "cp", par_decom = 20, ac = 0.7)
}

# No change-points and 10 components are kept after the decomposition
SB_Tensors_0_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_0_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "cp", par_decom = 10, ac = 0.7)
}

# No change-points and 5 components are kept after the decomposition
SB_Tensors_0_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_0_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = integer(0), size = 200, rho = 0.8, num_subgraph = 4, decomp = "cp", par_decom = 5, ac=0.7)
}

# One change-point and 20 components are kept after the decomposition
SB_Tensors_1_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_1_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "cp", par_decom = 20, ac=0.7)
}

# One change-point and 10 components are kept after the decomposition
SB_Tensors_1_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_1_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "cp", par_decom = 10, ac=0.7)
}

# One change-point and 5 components are kept after the decomposition
SB_Tensors_1_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_1_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = 100, size = 200, rho = c(0.8,0.2), num_subgraph = c(4,2), decomp = "cp", par_decom = 5, ac=0.7)
}

# Four change-points and 20 components are kept after the decomposition
SB_Tensors_4_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_4_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "cp", par_decom = 20, ac=0.7)
}

# Four change-points and 10 components are kept after the decomposition
SB_Tensors_4_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_4_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "cp", par_decom = 10, ac=0.7)
}

# Four change-points and 5 components are kept after the decomposition
SB_Tensors_4_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_4_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(100, 150, 200, 250), size = 300, rho = c(0.8, 0.2, 0.8, 0.2, 0.8), num_subgraph = c(4,2,4,2,4), decomp = "cp", par_decom = 5, ac=0.7)
}

# Ten change-points and 20 components are kept after the decomposition
SB_Tensors_10_cpt_20_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_10_cpt_20_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "cp", par_decom = 20, ac=0.7)
}

# Ten change-points and 10 components are kept after the decomposition
SB_Tensors_10_cpt_10_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_10_cpt_10_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "cp", par_decom = 10, ac=0.7)
}

# Ten change-points and 5 components are kept after the decomposition
SB_Tensors_10_cpt_5_components_ar_0.7 <- list()
set.seed(seed.temp)
for (i in 1:500){
  print(i)
  SB_Tensors_10_cpt_5_components_ar_0.7[[i]] <- tensor_creation_decomposition_SB(dims = c(20,20,20), N = c(60, 120, 180, 240,300,360,420,480,540,600), size = 660, rho = c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8), num_subgraph = c(4,2,4,2,4,2,4,2,4,2,4), decomp = "cp", par_decom = 5, ac=0.7)
}

SB_SIM0_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_0_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
SB_SIM1_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_1_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
SB_SIM4_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_4_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
SB_SIM10_cp_TenSeg_20_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_10_cpt_20_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)

SB_SIM0_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_0_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
SB_SIM1_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_1_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
SB_SIM4_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_4_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
SB_SIM10_cp_TenSeg_10_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_10_cpt_10_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)

SB_SIM0_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_0_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = integer(0), m = 500, seed = seed.temp)
SB_SIM1_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_1_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = 100, m = 500, seed = seed.temp)
SB_SIM4_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_4_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(100,150,200,250), m = 500, seed = seed.temp)
SB_SIM10_cp_TenSeg_5_components_ar_0.7 <- simulation_study_TenSeg(decomposed_tensor = SB_Tensors_10_cpt_5_components_ar_0.7, decomp_sim = "cp", points_sim = 10, true.cpt = c(60,120,180,240,300,360,420,480,540,600), m = 500, seed = seed.temp)
