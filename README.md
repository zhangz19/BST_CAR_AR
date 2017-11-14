# BST_CAR_AR
Bayesian Spatio-temporal model with separable CAR-AR covariance structure. 

* The code is developed for: Hamaamin, Y. A., Nejadhashemi, A. P., Zhang, Z., Giri, S. and Woznicki, S. A. (2016). Bayesian Regression and Neuro Fuzzy Methods Reliability Assessment for Estimating Streamflow. _Water_. 8(7), 287.

* The input: 
  - Y: NT*1 column vector with T blocks: each block Yt has observation for N locations, i.e. Y=[Y1',...YT']';
In other woard, reshape(Y, [N,T]) should give you an N-by-T matrix with each row being the time series (1-T) for each location (1-N). You can plot(Y') to plot the N time series.
  - Likewise, X is NT*p with p covariates orgnized in the same way as Y.
  - W is an N-by*N neighborhood matrix which is binary: W(i,j)=1 if location i is adjacent with location j (sharing borderline), and 0 otherwise. If your data is non-spatial, set W=eye(N) and turn nmodel=1 or 3 with non-spatial component.
  - load the example data for demo ("load yourData.mat"). 

* The model is: Y = X*beta + Zu + epsilon, with u: NT by 1 spatio-temporal random effects with covariance being the kronecker product of A: AR(1) model correlation for time series; and D: Conditional Auto-regressive(CAR) model for spatial data Both assume Markovian structure and have closed-form, sparse inverse. 

* This function also evaluates DIC4 for mixed-effects model for comparison. This function can be used for functional regression using spike-and-slab prior for wavelet. It contains some useful util functions for AR and CAR model such as matrix inverse, determinant and Cholesky inverse by taking full computational advantage of AR and CAR covariance structure. 

* Please contact the authors if there are any questions or implementation issues: Zhen Zhang, zhangquake@gmail.com. Date coded: 2013-6-17
