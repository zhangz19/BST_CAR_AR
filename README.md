# BST_CAR_AR
Bayesian Spatio-temporal model with separable CAR-AR covariance structure. 

* The code is developed for: Hamaamin, Y. A., Nejadhashemi, A. P., Zhang, Z., Giri, S. and Woznicki, S. A. (2016). Bayesian Regression and Neuro Fuzzy Methods Reliability Assessment for Estimating Streamflow. _Water_. 8(7), 287.

* The model is: Y = X*beta + Zu + epsilon, with u: NT by 1 spatio-temporal random effects with covariance being the kronecker product of A: AR(1) model correlation for time series; and D: Conditional Auto-regressive(CAR) model for spatial data Both assume Markovian structure and have closed-form, sparse inverse. 

* This function also evaluates DIC4 for mixed-effects model for comparison. This function can be used for functional regression using spike-and-slab prior for wavelet. It contains some useful util functions for AR and CAR model such as matrix inverse, determinant and Cholesky inverse by taking full computational advantage of AR and CAR covariance structure. 

* Please contact the authors if there are any questions or implementation issues: Zhen Zhang, zhangz19@galton.uchicago.edu. Date coded: 2013-6-17
