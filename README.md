# Promotion Time Cure Model for Interval-Censored Survival Data

This repository contains R implementations of the Promotion Time Cure Model (PTCM) for interval-censored survival data, incorporating neural networks and conventional regression approaches for incidence modeling.

The framework uses an Expectation-Maximization (EM) algorithm for parameter estimation, with a Weibull proportional hazards model for the latency component.

## Methods

The following incidence modeling approaches are implemented:

- Logistic regression
- Natural spline regression
- Neural networks (NN)

The implementations estimate uncured probabilities and survival quantities for the overall population and susceptible subpopulation.

## Repository Contents

### 1. `NN_IC_PCM_Simulation.R`

Implements a simulation study for interval-censored promotion time cure models.

The script includes:

- Generation of interval-censored survival observations.
- Three incidence probability settings, including linear and nonlinear relationships.
- Weibull proportional hazards latency modeling.
- EM-based estimation.
- Training and testing data partitioning.
- Evaluation of estimated uncured probabilities and survival functions.
- Calculation of bias, absolute bias, mean squared error (MSE), and AUC.

The implemented methods can be selected through the `methods` variable. The current script is configured to run logistic regression.

### 2. `NN_IC_PCM_Data_With_Plots.R`

Implements a real-data application using a smoking cessation dataset.

The script includes:

- Data preparation and covariate standardization.
- Logistic regression, spline regression, and neural network incidence modeling.
- Weibull proportional hazards latency modeling.
- Stratified 10-fold cross-validation.
- Estimation of uncured probabilities and survival functions.
- Incidence probability surface visualization.
- Optional bootstrap-based inference.

## R Packages

The scripts use packages including:

`survival`, `neuralnet`, `splines`, `MASS`, `caTools`, `dplyr`, `pROC`, `e1071`, `xtable`, `writexl`, `readxl`, `PRROC`, `rpart`, `randomForest`, and `xgboost`.

Install the required packages before running the scripts.

## Running the Code

1. Download or clone the repository.
2. Open the relevant `.R` script in R or RStudio.
3. Install the required R packages.
4. Configure the simulation parameters or specify the path to the smoking cessation dataset.
5. Run the selected analysis.

**Note:** The real-data script requires an external smoking cessation dataset that is not included in this repository.

## Research Status

This repository contains research code for the development and evaluation of promotion time cure models under interval censoring. The implementations may be updated as the research progresses.

## Author

**Alimatu Agongo**

Department of Mathematics  
University of Texas at Arlington
