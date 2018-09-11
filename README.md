# DTRlearn2

## R Package of Statistical Learning Methods for Optimizing Dynamic Treatment Regimes

We provide a comprehensive software to estimate general K-stage DTRs from SMARTs with Q-learning and a variety of outcome-weighted learning methods. Penalizations are allowed for variable selection and model regularization. With the outcome-weighted learning scheme, different loss functions - SVM hinge loss, SVM ramp loss, binomial deviance loss, and L2 loss - are adopted to solve the weighted classification problem at each stage; augmentation in the outcomes is allowed to improve efficiency. The estimated DTR can be easily applied to a new sample for individualized treatment recommendations or DTR evaluation.

Author: Yuan Chen, Ying Liu, Donglin Zeng, Yuanjia Wang 

Maintainer: Yuan Chen <irene.yuan.chen@gmail.com>, <yc3281@columbia.edu>


### Installation in R
First install the "devtools" package:

install.packages("devtools")

Then install the "DTRlearn2" package from github:

library(devtools)

install_github("ychen178/DTRlearn2")
