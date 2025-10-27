##  Particle Filter Estimation with Two-step RJMCMC Resampling for Sequential Hierarchical Bayesian Model 

### Overview
This repository provides the implementation code for our novel  with **Two-Stage Resampling Particle Filter** for estimating **Sequential Hierarchical Bayesian Model (SHBM)**, addressing key challenges in Data Assimilation (DA) involving uncertain dynamic model structures and varying parameter dimensions.

### Key Features
- **Sequential Hierarchical Bayesian Model (SHBM)**: A flexible framework integrating state space modeling with hierarchical parameter structures to handle model uncertainty
- **Two-Stage Resampling Scheme**:
  - **Stage 1**: Bootstrap particle filter for efficient preliminary resampling
  - **Stage 2**: RJMCMC-based resampling for particle diversification in trans-dimensional spaces
- **Simultaneous Capabilities**:
  - Automatic model structure identification
  - Parameter estimation in variable-dimensional spaces
  - State assimilation and prediction

### Problem Solved
Traditional DA methods often suffer from:
- **Filter divergence** due to incorrect model assumptions
- **Particle impoverishment** in high-dimensional parameter spaces
- **Inability to handle** varying parameter dimensions across different model structures

Our method effectively overcomes these limitations through the two-stage resampling approach, maintaining particle diversity while preventing filter divergence.

### Applications
- Advection equation models
- Lorenz 96 systems
- Spatiotemporal process studies
- Systems with structural uncertainties and parameter dimension variations

### Reference
For theoretical foundations, methodological details, and experimental results, please refer to our accompanying publication.

## Installation

You can install the development version of SIStree from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC")
```




## 📚 File Description

### 🎯 Core Algorithm Files
- **RJ2S_basic_function_and_setting.R** - Provides basic functions for PF with two-stage RJMCMC resampling method, including proposal functions for four transition types and the main DA process function `DA_Process(Enstat, EnPara, Mthd)`


### 📊 Analysis Scripts
- **Loop_R2S_prior_compare2_alpha_and_c.R** - Sensitivity analysis for parameter alpha, comparing MSEs and MSPEs of the proposed method under different alpha values through independent repeated simulations
- **Loop_RU2S_prior_compare_3_c.R** - Sensitivity analysis for parameter c, comparing MSEs and MSPEs of the proposed method under different c values through independent repeated simulations
- **Loop_R2S_Non_para_method_compare.R** - Comparison of three parameter distribution estimation methods, evaluating Acceptance rate, ESS, MSEs, and MSPEs through independent repeated simulations
- **Loop_R2S_Non_para_bandwidth_compare.R** - Bandwidth analysis, comparing Acceptance rate, ESS, MSEs, and MSPEs for different bandwidths in parameter estimation method through independent repeated simulations

### 📈 Plotting Scripts
- **Plot_RJ2s_prior_compare_3_c.R** - Generates plots for parameter c sensitivity analysis based on results from "Loop RU2S prior compare 3_c.R", corresponding to the right subplot of Figure 2 in the associated paper
- **Plot_R2s_prior_compare_2_alpha.R** - Generates plots for parameter alpha sensitivity analysis based on results from "Loop R2S prior_compare2 alpha and c.R", corresponding to the left subplot of Figure 2 in the associated paper

### 🔧 Utility Functions
- **Func_R2S_test_acceptance_and_ess.R** - Contains functions to calculate monitoring metrics such as acceptance rate and ESS for each of the four transition types at each time step
- **Func_R2S_Para_KDE.R** - Kernel density estimation functions used for non-parametric estimation
- **Func_R2S_add_perturbation.R** - Contains functions to add perturbations to duplicated particles after the first resampling step
- **Func_RU2S_Para_estimation_proposal.R** - Parameter distribution estimation methods and corresponding DA main function `DA_Process_EstP`
- **Func_RU2S_Non_para_estimation_proposal.R** - Non-parametric distribution function estimation methods and corresponding DA main function `DA_Process_NonP`

## Usage

We show how to use the code with examples.

## 🚀 Quick Start

### 1. Load the necessary scripts
```r
R_path <- ''

source(paste0(R_path, 'RJ2S_basic function and setting.R'))
source(paste0(R_path, 'Func_RJ2S  add_perturbation.R'))
source(paste0(R_path, 'Func_RJ2S  Para KDE.R'))
source(paste0(R_path, 'Func_RJ2S  test acceptance and ess.R'))
source(paste0(R_path, 'Func_RJ2S Para estimation proposal.R'))
source(paste0(R_path, 'Func_RJ2S Non_para estimation proposal.R'))

```

### 2. 
