# DTRlearn2

## Statistical Learning Methods for Optimizing Dynamic Treatment Regimes

We provide a comprehensive software to estimate general K-stage DTRs from SMARTs with Q-learning and a variety of outcome-weighted learning methods. Penalizations are allowed for variable selection and model regularization. With the outcome-weighted learning scheme, different loss functions - SVM hinge loss, SVM ramp loss, binomial deviance loss, and L2 loss - are adopted to solve the weighted classification problem at each stage; augmentation in the outcomes is allowed to improve efficiency. The estimated DTR can be easily applied to a new sample for individualized treatment recommendations or DTR evaluation.

- Author: **Yuan Chen<sup>a</sup>, Ying Liu<sup>b</sup>, Donglin Zeng<sup>c</sup>, Yuanjia Wang<sup>a,d</sup>** 
- Maintainer: **Yuan Chen** (<irene.yuan.chen@gmail.com>, <yc3281@cumc.columbia.edu>)
- Affiliations: 
  + 1. Department of Biostatistics, Mailman School of Public Health, Columbia University, New York**
  + 2. Department of Psychiatry, Columbia University Irving Medical Center, New York**
  + 3. Department of Biostatistics, University of North Carolina, Chapel Hill, North Carolina**
  + 4. Department of Psychiatry, Columbia University, New York

#### Reference
Chen, Yuan, Ying Liu, Donglin Zeng, and Yuanjia Wang. "Statistical Learning Methods for Optimizing Dynamic Treatment Regimes in Subgroup Identification." In Design and Analysis of Subgroups with Biopharmaceutical Applications, pp. 271-297. Springer, Cham, 2020.




### Installation in R
First install the "devtools" package:

install.packages("devtools")

Then install the "DTRlearn2" package from github:

library(devtools)

install_github("ychen178/DTRlearn2")
