# SACut: an R package for Stochastic Approximation Cut Algorithm (SACut)
[![](https://travis-ci.com/MathBilibili/Stochastic-approximation-cut-algorithm.svg?branch=master)](https://travis-ci.com/MathBilibili/Stochastic-approximation-cut-algorithm)

This is an R package to conduct the Stochastic Approximation Cut Algorithm. It also contains code that replicates the results in the paper.
<img align="right" width="250" height="250" src="https://user-images.githubusercontent.com/24710640/81212775-495b2880-8fcd-11ea-9319-52ac4fd15f4f.png">

## Authors
Yang Liu and Robert Goudie

MRC Biostatistics Unit, University of Cambridge

## Introduction of SACut
The potential effect of partial misspecification of Bayesian modelling is a concern. Recent studies have proposed the idea of modularized models, and the cut model is proposed to prevent feedback of information from the suspect module. This leads to the cut posterior distribution which normally does not have a closed-form. Previous studies have proposed algorithms to sample from this distribution, but these algorithms have unclear theoretical convergence properties. To address this convergence problem, the novel Stochastic Approximation Cut algorithm (SACut) is proposed as an alternative. The algorithm is divided into two parallel chains. The main chain targets an approximation of the cut distribution; the auxiliary chain forms a proposal distribution used in the main chain. We prove convergence of the samples drawn by the proposed algorithm and present the exact limit. Although SACut is biased, since the main chain does not target the exact cut distribution, we prove this bias can be reduced geometrically by increasing a user-chosen tuning parameter. In addition, parallel processing can be easily adopted for SACut, unlike existing algorithms. This greatly reduces computation time.
