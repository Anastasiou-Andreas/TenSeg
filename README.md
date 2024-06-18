# Scope

This repository contains R code that can be used to replicate the results in the paper: ``_Tensor time series change-point detection in cryptocurrency network data_''

There are three parts in this repository:

1. The main function that implements the TenSeg algorithm. In the same file, the functions creating tensors when the precision matrices undertake possible changes in either the AR1($\rho$), or the Star-Block (SB), or the Erdos-Renyi structure.
2. The functions for the comparative simulation study.
3. The functions for the results of the simulation study related to only TenSeg when serial correlation was added to the data.

The scripts for the first part are organised in the folder _Main functions_, while the scripts for the last two parts are organised in the folder _Simulations_.
