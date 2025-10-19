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

## Usage

We show how to use the SIStree package with examples.
