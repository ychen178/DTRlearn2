
########### Simulation for the SMARTER paper (modification from previous)
# modifications: 1) A1 and A2 to depend on variables
# 2) add gold-standard: Study 2 as a SMART to learn the DTR
# 3) to have less variability for the error term in the outcome
# 4) to add a gold standard (Study 2 as SMART to learn the whole rule)
#### some notes:
#  don't add R1 to R2, in which case Naive 1 would be very case since R2 is mainly determined by R1
#  in sensitivity analysis, don't just scale the decision function, since the recommneded treatments wont change

library(R.utils)
sourceDirectory("/Users/YChen/Dropbox/MINE/00_CU/Projects/04_Software Paper/02_code/01_R/DTRlearn2/R")
library(glmnet)
library(Matrix)
library(kernlab)
library(MASS)


#######  define function to run Naive1, Naive2, our method + gold standard
# add gold-standard (deem Study 2 as a SMART)
runSim4 = function(simfun1,simfun2, X_t_s1,X_t_s2,dat_t, Vdat, nstart,nend, n_sample1,n_sample2, c, one,two,three,four)  {
  for (i in nstart : nend ) {

    print(c("i", i))
    t = i
    set.seed(2*t-1)
    dat_tr1 = simfun1 (n_sample1) # data for learn stage 1
    set.seed(2*t)
    dat_tr2 = simfun2 (n_sample2) # data for learn stage 2

    X_tr1_s1 = scale(dat_tr1$X)
    X_tr1_s2 = scale( cbind(dat_tr1$X, dat_tr1$A[[1]], dat_tr1$R1, dat_tr1$X * dat_tr1$A[[1]], dat_tr1$R1 * dat_tr1$A[[1]]) ) # used for predict the Q function
    X_tr2_s1 = scale(dat_tr2$X)
    X_tr2_s2 = scale( cbind(dat_tr2$X, dat_tr2$A[[1]], dat_tr2$R1, dat_tr2$X * dat_tr2$A[[1]], dat_tr2$R1 * dat_tr2$A[[1]]) )

    # use generated probability instead of empirical proportion
    p_tr1_s1 = dat_tr1$pi[[1]]
    p_tr2_s2 = dat_tr2$pi[[2]]

    owl_s2 = owl(H=X_tr2_s2, AA=dat_tr2$A[[2]], RR=dat_tr2$R[[2]], n=n_sample2, K=1, pi=p_tr2_s2, loss='hinge', augment=F, c=c, m=3)
    naive1_s1 = owl(H=X_tr1_s1, AA=dat_tr1$A[[1]], RR=dat_tr1$R[[1]], n=n_sample1, K=1, pi=p_tr1_s1, loss='hinge', augment=F, c=c, m=3)
    naive2_s1 = owl(H=X_tr1_s1, AA=dat_tr1$A[[1]], RR=dat_tr1$R[[1]] + dat_tr1$R[[2]], n=n_sample1, K=1, pi=p_tr1_s1, loss='hinge', augment=F, c=c, m=3)

    ql_s2 = ql(H=X_tr2_s2, AA=dat_tr2$A[[2]], RR=dat_tr2$R[[2]], K=1, pi=p_tr2_s2, m=3)
    # since dat1 all drop out, the outcome is the Q function (since no intermediate outcome)
    Q_tr1 = predict(ql_s2, X_tr1_s2, K=1)$Q[[1]]
    owl_s1 = owl(H=X_tr1_s1, AA=dat_tr1$A[[1]], RR=dat_tr1$R[[1]] + Q_tr1, n=n_sample1, K=1, pi=p_tr1_s1, loss='hinge', augment=F, c=c, m=3)
    # use Study2 as a SMART to learn the rule
    gold = owl(H=list(X_tr2_s1, X_tr2_s2), AA=dat_tr2$A, RR=dat_tr2$R, n=n_sample2, K=2, pi=dat_tr2$pi, loss='hinge', augment=T, c=c, m=3)

    # predict treatments at stage 1 & 2
    pred_a2 = predict(owl_s2, X_t_s2, K=1)$treatment[[1]]
    naive1_a1 = predict(naive1_s1, X_t_s1, K=1)$treatment[[1]]
    naive2_a1 = predict(naive2_s1, X_t_s1, K=1)$treatment[[1]]
    pred_a1 = predict(owl_s1, X_t_s1, K=1)$treatment[[1]]
    pred_gold = predict(gold, list(X_t_s1, X_t_s2), K=2)$treatment

    # calculate the value function
    naive1_v = mean( (pred_a2==dat_t$A[[2]] & naive1_a1==dat_t$A[[1]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    naive2_v = mean( (pred_a2==dat_t$A[[2]] & naive2_a1==dat_t$A[[1]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    owl_v = mean( (pred_a2==dat_t$A[[2]] & pred_a1==dat_t$A[[1]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    gold_v = mean( (pred_gold[[1]]==dat_t$A[[1]] & pred_gold[[2]]==dat_t$A[[2]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )

    Vdat[i, one] = naive1_v
    Vdat[i, two] = naive2_v
    Vdat[i, three] = owl_v
    Vdat[i, four] = gold_v
    print(c(naive1_v, naive2_v, owl_v, gold_v))
  }
  return(Vdat)
}


#########################################################################################
#  use the new settings  Sim1.2 & Sim2.2 for simulating two studies
#########################################################################################


#-------------------------  "Sim1" (modified Sim12)   ------------------------------#
Sim1 = function(n) {
  mu = rep(0, 10)
  sigma1 = ifelse(diag(5), 1, 0.2)
  sigma2 = diag(5)
  sigma = rbind( cbind(sigma1, matrix(0, 5, 5)), cbind(matrix(0, 5, 5), sigma2) )
  X = mvrnorm(n=n, mu=mu, Sigma=sigma)

  pA1 = 1 / (1 + exp(-0.5 * X[,1] + X[,2]))  # added
  A1 = 2 * rbinom(n=n, size=1, prob=pA1) - 1
  decfun1 = X[,1] + X[,2] + X[,2]^2
  R1 = decfun1 * A1 + rnorm(n, 0, 0.2)  # was sd=1
  y1 = 2 * (decfun1 > 0) - 1  # optimal

  pA2 = 1 / (1 + exp(0.1 * R1))  # added, assignment depends on R1
  A2 = 2 * rbinom(n, 1, pA2) - 1
  decfun2 = R1 + 2*X[,3] + X[,4]^2 + 1.5*A1*X[,4] - 0.2  # R1 + 2 X3 + X4^2 + 1.5*A1 - 0.2
  R2 = 3*X[,3] + 2*A1*X[,5]  + decfun2 * A2 + rnorm(n, 0, 0.2)  # was sd=1
  y2 = 2 * (decfun2 > 0) - 1  # optimal

  list(X=X, A=list(A1,A2), R=list(R1,R2), optA=list(y1,y2), pi=list(pA1,pA2) )
}

##  create a big test set for "Sim1" --  optimal:
{
set.seed(1)
n_test = 20000
dat1_t = Sim1(n_test)

  # check the optimal    # 5.352199
  mean( (dat1_t$R[[1]] + dat1_t$R[[2]]) * (dat1_t$A[[1]]==dat1_t$optA[[1]]) * (dat1_t$A[[2]]==dat1_t$optA[[2]]) /
          (dat1_t$pi[[1]] * dat1_t$pi[[2]]) )   # note used the known randomization probability

X1_t_s1 = scale(dat1_t$X)
X1_t_s2 = scale( cbind(dat1_t$X, dat1_t$A[[1]], dat1_t$R1, dat1_t$X * dat1_t$A[[1]], dat1_t$R1 * dat1_t$A[[1]]) )

  summary(dat1_t$R[[1]])
  summary(dat1_t$R[[2]])
  par(mfrow=c(1,2))
  hist(dat1_t$R[[1]])
  hist(dat1_t$R[[2]])

}
#---------------- end of setting up "Sim1" ------------#
## run Sim12
V_Sim1 = matrix(NA, nsim, 16)
clinear = 2^(-5:5)

print("Sim1")
V_Sim1 = runSim4(Sim1,Sim1,  X1_t_s1,X1_t_s2,dat1_t, V_Sim1, nstart=201,nend=500,  n_sample1=100,n_sample2=100, c=clinear, 1,2,3,4)
V_Sim1 = runSim4(Sim1,Sim1,  X1_t_s1,X1_t_s2,dat1_t, V_Sim1, nstart=201,nend=500,  n_sample1=100,n_sample2=200, c=clinear, 5,6,7,8)
V_Sim1 = runSim4(Sim1,Sim1,  X1_t_s1,X1_t_s2,dat1_t, V_Sim1, nstart=201,nend=500,  n_sample1=200,n_sample2=100, c=clinear, 9,10,11,12)
V_Sim1 = runSim4(Sim1,Sim1,  X1_t_s1,X1_t_s2,dat1_t, V_Sim1, nstart=201,nend=500,  n_sample1=200,n_sample2=200, c=clinear, 13,14,15,16)

    apply(V_Sim1[1:500,], 2, mean, na.rm=T)
    apply(V_Sim1[1:500,], 2, sd, na.rm=T)
    apply(V_Sim1[1:500,], 2, median, na.rm=T)






#-----------------------------  "Sim2"  -------------------------------#
# add latent variables U1 and U2 to "Sim1"

Sim2 = function(n) {
  mu = rep(0, 10)
  sigma1 = ifelse(diag(5), 1, 0.2)
  sigma2 = diag(5)
  sigma = rbind( cbind(sigma1, matrix(0, 5, 5)), cbind(matrix(0, 5, 5), sigma2) )
  X = mvrnorm(n=n, mu=mu, Sigma=sigma)
  # add latent vars
  U1 = rnorm(n, 1, 1)
  U2 = rnorm(n, 0.5, 1)

  pA1 = 1 / (1 + exp(-0.5 * X[,1] + X[,2]))  # added
  A1 = 2 * rbinom(n=n, size=1, prob=pA1) - 1
  decfun1 = X[,1] + X[,2] + X[,2]^2 + U1
  R1 = decfun1 * A1 + rnorm(n, 0, 0.2)  # was sd=1
  y1 = 2 * (decfun1 > 0) - 1  # optimal

  pA2 = 1 / (1 + exp(0.1 * R1))  # added, assignment depends on R1
  A2 = 2 * rbinom(n, 1, pA2) - 1
  decfun2 = R1 + 2*X[,3] + X[,4]^2 + 1.5*A1*X[,4] - 0.2 + U2  # R1 + 2*X[,3] + X[,4]^2 + 1.5*A1 - 0.2
  R2 = 3*X[,3] + 2*A1*X[,5]  + decfun2 * A2 + rnorm(n, 0, 0.2)
  y2 = 2 * (decfun2 > 0) - 1  # optimal

  list(X=X, A=list(A1,A2), R=list(R1,R2), optA=list(y1,y2), pi=list(pA1,pA2) )
}
##  create a big test set
{
  set.seed(1)
  n_test = 20000
  dat2_t = Sim2(n_test)

  # check the optimal value function   6.786551
  mean( (dat2_t$R[[1]] + dat2_t$R[[2]]) * (dat2_t$A[[1]]==dat2_t$optA[[1]]) * (dat2_t$A[[2]]==dat2_t$optA[[2]])
        / (dat2_t$pi[[1]] * dat2_t$pi[[2]]) )

  X2_t_s1 = scale(dat2_t$X)
  X2_t_s2 = scale( cbind(dat2_t$X, dat2_t$A[[1]], dat2_t$R1, dat2_t$X * dat2_t$A[[1]], dat2_t$R1 * dat2_t$A[[1]]) )

  summary(dat2_t$R[[1]])
  summary(dat2_t$R[[2]])
  par(mfrow=c(1,2))
  hist(dat2_t$R[[1]])
  hist(dat2_t$R[[2]])
}


## Run "Sim2"
V_Sim2 = matrix(NA, nsim, 16)
V_Sim2 = runSim4(Sim2,Sim2,  X2_t_s1,X2_t_s2,dat2_t,  V_Sim2, nstart=201,nend=500, n_sample1=100,n_sample2=100,  c=clinear, 1,2,3,4)
V_Sim2 = runSim4(Sim2,Sim2,  X2_t_s1,X2_t_s2,dat2_t,  V_Sim2, nstart=201,nend=500, n_sample1=100,n_sample2=200,  c=clinear, 5,6,7,8)
V_Sim2 = runSim4(Sim2,Sim2,  X2_t_s1,X2_t_s2,dat2_t,  V_Sim2, nstart=201,nend=500, n_sample1=200,n_sample2=100,  c=clinear, 9,10,11,12)
V_Sim2 = runSim4(Sim2,Sim2,  X2_t_s1,X2_t_s2,dat2_t,  V_Sim2, nstart=201,nend=500, n_sample1=200,n_sample2=200,  c=clinear, 13,14,15,16)

    apply(V_Sim2[1:200,], 2, mean, na.rm=T)
    apply(V_Sim2[1:200,], 2, sd, na.rm=T)
    apply(V_Sim2[1:50,], 2, median, na.rm=T)



#-------------------------   Sensitivity Analyses   ------------------------------#
# The outcome functions differ in study 2 by a factor alpha (one term times alpha)
# test set from Study 1 (same test set as in "Sim1")
Sim1a = function(n, alpha) {
  mu = rep(0, 10)
  sigma1 = ifelse(diag(5), 1, 0.2)
  sigma2 = diag(5)
  sigma = rbind( cbind(sigma1, matrix(0, 5, 5)), cbind(matrix(0, 5, 5), sigma2) )
  X = mvrnorm(n=n, mu=mu, Sigma=sigma)

  pA1 = 1 / (1 + exp(-0.5 * X[,1] + X[,2]))  # added
  A1 = 2 * rbinom(n=n, size=1, prob=pA1) - 1
  decfun1 = X[,1] + (1+alpha) * (X[,2] + X[,2]^2)
  R1 = decfun1 * A1 + rnorm(n, 0, 0.2)  # was sd=1
  y1 = 2 * (decfun1 > 0) - 1  # optimal

  pA2 = 1 / (1 + exp(0.1 * R1))  # added, assignment depends on R1
  A2 = 2 * rbinom(n, 1, pA2) - 1
  decfun2 = (1+alpha) * (R1 + 2*X[,3] + X[,4]^2) + 1.5*A1*X[,4] - 0.2
  R2 = 3*X[,3] + 2*A1*X[,5]  + decfun2 * A2 + rnorm(n, 0, 0.2)  # was sd=1
  y2 = 2 * (decfun2 > 0) - 1  # optimal

  list(X=X, A=list(A1,A2), R=list(R1,R2), optA=list(y1,y2), pi=list(pA1,pA2) )
}


####### Define a function to run simulations with input alpha & predict on two test sets (one same with Study 1 and one same with Study 2)
# X_t_s1_a,X_t_s2_a,dat_t_a  is the test set with alpha
runSim5 = function(simfun1,alpha1, simfun2,alpha2,  X_t_s1,X_t_s2,dat_t,  X_t_s1_a,X_t_s2_a,dat_t_a,  Vdat, nstart,nend, n_sample1,n_sample2, c, one,two,three,four,  five,six,seven,eight)  {
  for (i in nstart : nend ) {

    print(c("i", i))
    t = i
    set.seed(2*t-1)
    dat_tr1 = simfun1 (n_sample1, alpha1) # data for learn stage 1
    set.seed(2*t)
    dat_tr2 = simfun2 (n_sample2, alpha2) # data for learn stage 2

    X_tr1_s1 = scale(dat_tr1$X)
    X_tr1_s2 = scale( cbind(dat_tr1$X, dat_tr1$A[[1]], dat_tr1$R1, dat_tr1$X * dat_tr1$A[[1]], dat_tr1$R1 * dat_tr1$A[[1]]) ) # used for predict the Q function
    X_tr2_s1 = scale(dat_tr2$X)
    X_tr2_s2 = scale( cbind(dat_tr2$X, dat_tr2$A[[1]], dat_tr2$R1, dat_tr2$X * dat_tr2$A[[1]], dat_tr2$R1 * dat_tr2$A[[1]]) )

    # use generated probability instead of empirical proportion
    p_tr1_s1 = dat_tr1$pi[[1]]
    p_tr2_s2 = dat_tr2$pi[[2]]

    owl_s2 = owl(H=X_tr2_s2, AA=dat_tr2$A[[2]], RR=dat_tr2$R[[2]], n=n_sample2, K=1, pi=p_tr2_s2, loss='hinge', augment=F, c=c, m=3)
    naive1_s1 = owl(H=X_tr1_s1, AA=dat_tr1$A[[1]], RR=dat_tr1$R[[1]], n=n_sample1, K=1, pi=p_tr1_s1, loss='hinge', augment=F, c=c, m=3)
    naive2_s1 = owl(H=X_tr1_s1, AA=dat_tr1$A[[1]], RR=dat_tr1$R[[1]] + dat_tr1$R[[2]], n=n_sample1, K=1, pi=p_tr1_s1, loss='hinge', augment=F, c=c, m=3)

    ql_s2 = ql(H=X_tr2_s2, AA=dat_tr2$A[[2]], RR=dat_tr2$R[[2]], K=1, pi=p_tr2_s2, m=3)
    # since dat1 all drop out, the outcome is the Q function (since no intermediate outcome)
    Q_tr1 = predict(ql_s2, X_tr1_s2, K=1)$Q[[1]]
    owl_s1 = owl(H=X_tr1_s1, AA=dat_tr1$A[[1]], RR=dat_tr1$R[[1]] + Q_tr1, n=n_sample1, K=1, pi=p_tr1_s1, loss='hinge', augment=F, c=c, m=3)
    # use Study2 as a SMART to learn the rule
    gold = owl(H=list(X_tr2_s1, X_tr2_s2), AA=dat_tr2$A, RR=dat_tr2$R, n=n_sample2, K=2, pi=dat_tr2$pi, loss='hinge', augment=T, c=c, m=3)

    # predict treatments at stage 1 & 2 on the first test set
    pred_a2 = predict(owl_s2, X_t_s2, K=1)$treatment[[1]]
    naive1_a1 = predict(naive1_s1, X_t_s1, K=1)$treatment[[1]]
    naive2_a1 = predict(naive2_s1, X_t_s1, K=1)$treatment[[1]]
    pred_a1 = predict(owl_s1, X_t_s1, K=1)$treatment[[1]]
    pred_gold = predict(gold, list(X_t_s1, X_t_s2), K=2)$treatment
    # predict treatments at stage 1 & 2 on the second test set
    pred_a2_a = predict(owl_s2, X_t_s2_a, K=1)$treatment[[1]]
    naive1_a1_a = predict(naive1_s1, X_t_s1_a, K=1)$treatment[[1]]
    naive2_a1_a = predict(naive2_s1, X_t_s1_a, K=1)$treatment[[1]]
    pred_a1_a = predict(owl_s1, X_t_s1_a, K=1)$treatment[[1]]
    pred_gold_a = predict(gold, list(X_t_s1_a, X_t_s2), K=2)$treatment

    ## calculate the value function
    # for prediction on the first test set
    naive1_v = mean( (pred_a2==dat_t$A[[2]] & naive1_a1==dat_t$A[[1]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    naive2_v = mean( (pred_a2==dat_t$A[[2]] & naive2_a1==dat_t$A[[1]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    owl_v = mean( (pred_a2==dat_t$A[[2]] & pred_a1==dat_t$A[[1]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    gold_v = mean( (pred_gold[[1]]==dat_t$A[[1]] & pred_gold[[2]]==dat_t$A[[2]]) * (dat_t$R[[1]] + dat_t$R[[2]]) / (dat_t$pi[[1]] * dat_t$pi[[2]]) )
    # for prediction on the second test set
    naive1_v_a = mean( (pred_a2_a==dat_t_a$A[[2]] & naive1_a1_a==dat_t_a$A[[1]]) * (dat_t_a$R[[1]] + dat_t_a$R[[2]]) / (dat_t_a$pi[[1]] * dat_t_a$pi[[2]]) )
    naive2_v_a = mean( (pred_a2_a==dat_t_a$A[[2]] & naive2_a1_a==dat_t_a$A[[1]]) * (dat_t_a$R[[1]] + dat_t_a$R[[2]]) / (dat_t_a$pi[[1]] * dat_t_a$pi[[2]]) )
    owl_v_a = mean( (pred_a2_a==dat_t_a$A[[2]] & pred_a1_a==dat_t_a$A[[1]]) * (dat_t_a$R[[1]] + dat_t_a$R[[2]]) / (dat_t_a$pi[[1]] * dat_t_a$pi[[2]]) )
    gold_v_a = mean( (pred_gold_a[[1]]==dat_t_a$A[[1]] & pred_gold_a[[2]]==dat_t_a$A[[2]]) * (dat_t_a$R[[1]] + dat_t_a$R[[2]]) / (dat_t_a$pi[[1]] * dat_t_a$pi[[2]]) )

    Vdat[i, one] = naive1_v
    Vdat[i, two] = naive2_v
    Vdat[i, three] = owl_v
    Vdat[i, four] = gold_v

    Vdat[i, five] = naive1_v_a
    Vdat[i, six] = naive2_v_a
    Vdat[i, seven] = owl_v_a
    Vdat[i, eight] = gold_v_a

    print(c(naive1_v, naive2_v, owl_v, gold_v, "  ", naive1_v_a,naive2_v_a,owl_v_a,gold_v_a))
  }
  return(Vdat)
}



### create test sets
# case 1: Test set same with Study 1 (use dat1_t)
# Case 2: test set same with Study 2 (depends o妻子的浪漫旅行n alpha)
n_test = 20000

# a==0.5 test set
set.seed(1)
dat_a0.5_t = Sim1a(n_test, 0.5)
X1_a0.5_t_s1 = scale(dat_a0.5_t$X)
X1_a0.5_t_s2 = scale( cbind(dat_a0.5_t$X, dat_a0.5_t$A[[1]], dat_a0.5_t$R1, dat_a0.5_t$X * dat_a0.5_t$A[[1]], dat_a0.5_t$R1 * dat_a0.5_t$A[[1]]) )
    # check the optimal    # 7.903557
    mean( (dat_a0.5_t$R[[1]] + dat_a0.5_t$R[[2]]) * (dat_a0.5_t$A[[1]]==dat_a0.5_t$optA[[1]]) * (dat_a0.5_t$A[[2]]==dat_a0.5_t$optA[[2]]) /
            (dat_a0.5_t$pi[[1]] * dat_a0.5_t$pi[[2]]) )

# a==1 test set
set.seed(1)
dat_a1_t = Sim1a(n_test, 1)
X1_a1_t_s1 = scale(dat_a1_t$X)
X1_a1_t_s2 = scale( cbind(dat_a1_t$X, dat_a1_t$A[[1]], dat_a1_t$R1, dat_a1_t$X * dat_a1_t$A[[1]], dat_a1_t$R1 * dat_a1_t$A[[1]]) )
    # check the optimal    # 11.11476
    mean( (dat_a1_t$R[[1]] + dat_a1_t$R[[2]]) * (dat_a1_t$A[[1]]==dat_a1_t$optA[[1]]) * (dat_a1_t$A[[2]]==dat_a1_t$optA[[2]]) /
            (dat_a1_t$pi[[1]] * dat_a1_t$pi[[2]]) )

###### Run sensitivity analyses
# a==0.5
V_sen_a0.5 = matrix(NA, nsim, 36)
clinear = c(2^(-5: 5))

V_sen_a0.5 = runSim5(Sim1a,0, Sim1a,0.5,  X1_t_s1,X1_t_s2,dat1_t,  X1_a0.5_t_s1,X1_a0.5_t_s2,dat_a0.5_t,  V_sen_a0.5, nstart=101,nend=500,  n_sample1=100,n_sample2=100,  c=clinear,  1,2,3,4,  21,22,23,24)
V_sen_a0.5 = runSim5(Sim1a,0, Sim1a,0.5,  X1_t_s1,X1_t_s2,dat1_t,  X1_a0.5_t_s1,X1_a0.5_t_s2,dat_a0.5_t,  V_sen_a0.5, nstart=1,nend=500,  n_sample1=200,n_sample2=100,  c=clinear,  5,6,7,8,  25,26,27,28)
V_sen_a0.5 = runSim5(Sim1a,0, Sim1a,0.5,  X1_t_s1,X1_t_s2,dat1_t,  X1_a0.5_t_s1,X1_a0.5_t_s2,dat_a0.5_t,  V_sen_a0.5, nstart=1,nend=500,  n_sample1=100,n_sample2=200,  c=clinear,  9,10,11,12,  29,30,31,32)
V_sen_a0.5 = runSim5(Sim1a,0, Sim1a,0.5,  X1_t_s1,X1_t_s2,dat1_t,  X1_a0.5_t_s1,X1_a0.5_t_s2,dat_a0.5_t,  V_sen_a0.5, nstart=101,nend=500,  n_sample1=200,n_sample2=200,  c=clinear,  13,14,15,16,  33,34,35,36)

    apply(V_sen_a0.5[1:500,], 2, mean, na.rm=T)
    apply(V_sen_a0.5, 2, sd, na.rm=T)
    apply(V_sen_a0.5, 2, median, na.rm=T)

  # try to have bigger N1 and N2
    V_sen_a0.5_more = matrix(NA, nsim, 36)




#### run when a==1
V_sen_a1 = matrix(NA, nsim, 36)
clinear = c(2^(-5 : 5))

V_sen_a1 = runSim5(Sim1a,0, Sim1a,1,  X1_t_s1,X1_t_s2,dat1_t,  X1_a1_t_s1,X1_a1_t_s2,dat_a1_t,  V_sen_a1, nstart=31,nend=500,  n_sample1=100,n_sample2=100, c=clinear,  1,2,3,4,  21,22,23,24)
V_sen_a1 = runSim5(Sim1a,0, Sim1a,1,  X1_t_s1,X1_t_s2,dat1_t,  X1_a1_t_s1,X1_a1_t_s2,dat_a1_t,  V_sen_a1, nstart=1,nend=500,  n_sample1=200,n_sample2=100, c=clinear,  5,6,7,8,  25,26,27,28)
V_sen_a1 = runSim5(Sim1a,0, Sim1a,1,  X1_t_s1,X1_t_s2,dat1_t,  X1_a1_t_s1,X1_a1_t_s2,dat_a1_t,  V_sen_a1, nstart=1,nend=500,  n_sample1=100,n_sample2=200, c=clinear,  9,10,11,12,  29,30,31,32)
V_sen_a1 = runSim5(Sim1a,0, Sim1a,1,  X1_t_s1,X1_t_s2,dat1_t,  X1_a1_t_s1,X1_a1_t_s2,dat_a1_t,  V_sen_a1, nstart=101,nend=500,  n_sample1=200,n_sample2=200, c=clinear,  13,14,15,16,   33,34,35,36)

    apply(V_sen_a1, 2, mean, na.rm=T)
    apply(V_sen_a1[1:60,], 2, sd, na.rm=T)
    apply(V_sen_a1, 2, median, na.rm=T)






### save the R objects
setwd("/Users/YChen/Dropbox/MINE/00_CU/Projects/05_SMARTer/04_results")
save(V_Sim1, V_Sim2, V_sen_a0.5, V_sen_a1, file="SimRes_sep4.RData")




#########################################################################################
#   Plot  simulation results
#########################################################################################

# load the simulation results
load("SimRes.RData")

######### plot Setting 1
###  create dataframes for ggplots (grouped boxplot)
V_Sim1_plot <- data.frame(
  x = c(V_Sim1[,1:16]),  # convert by col
  y = rep(c("N1=100\nN2=100", "N1=100\nN2=200", "N1=200\nN2=100", "N1=200\nN2=200"), each=500*4),
  z=rep(rep(c("1", "2", "3", "4"), each=500), 4),
  stringsAsFactors = F)
V_Sim1_plot$y <- factor(V_Sim1_plot$y,levels = c("N1=100\nN2=100", "N1=200\nN2=100", "N1=100\nN2=200", "N1=200\nN2=200"))  # reorder y

library(ggplot2)
dev.off()
ggplot(V_Sim1_plot, aes(y, x, fill=factor(z))) +
  geom_hline(yintercept=5.352199, linetype="dashed", color="red") +
  geom_boxplot() +
  labs(title="Setting 1", x="", y="Value function", fill="Method") +
  scale_fill_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))


######### plot Setting 2
V_Sim2_plot <- data.frame(
  x = c(V_Sim2[,1:16]),
  y = rep(c("N1=100\nN2=100", "N1=100\nN2=200", "N1=200\nN2=100", "N1=200\nN2=200"), each=500*4),
  z=rep(rep(c("1", "2", "3", "4"), each=500), 4),
  stringsAsFactors = T)
V_Sim2_plot$y <- factor(V_Sim2_plot$y,levels = c("N1=100\nN2=100", "N1=200\nN2=100", "N1=100\nN2=200", "N1=200\nN2=200"))

ggplot(V_Sim2_plot, aes(y, x, fill=factor(z))) +
  geom_hline(yintercept=6.786551, linetype="dashed", color="red") +
  geom_boxplot() +
  labs(title="Setting 2", x="", y="Value function", fill="Method") +
  scale_fill_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))



##### Sensitivity (4*2 plots)
###  Plot 1 series (test set same as Study 1)
V_sen_plot1 = data.frame(
  valfun = c(apply(V_Sim1[,1:16],2,mean, na.rm=T), apply(V_sen_a0.5[,1:16],2,mean, na.rm=T), apply(V_sen_a1[,1:16],2,mean, na.rm=T)),
  alpha = rep(c(0, 0.5, 1), each=16),
  group = rep(rep(c("N1=100\nN2=100", "N1=200\nN2=100", "N1=100\nN2=200", "N1=200\nN2=200"), 3), each=4),
  method = rep(c("1", "2", "3", "4"), 12))
V_sen_plot1$group <- factor(V_sen_plot1$group,levels = c("N1=100\nN2=100", "N1=200\nN2=100", "N1=100\nN2=200", "N1=200\nN2=200"))

# sensitivity plot 11 (N1=100, N2=100)
sen_plot11 = ggplot(data=V_sen_plot1[V_sen_plot1$group=="N1=100\nN2=100",],
               aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + coord_cartesian( ylim = c(4,5.2)) +
  geom_hline(yintercept=5.352199, linetype="dashed", color="red") +
  labs(title="N1=100, N2=100", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot11

# sensitivity plot 12 (N1=200, N2=100)
sen_plot12 = ggplot(data=V_sen_plot1[V_sen_plot1$group=="N1=200\nN2=100",],
               aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + coord_cartesian( ylim = c(4,5.2)) +
  geom_hline(yintercept=5.352199, linetype="dashed", color="red") +
  labs(title="N1=200, N2=100", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot12

# sensitivity plot 13 (N1=100, N2=200)
sen_plot13 = ggplot(data=V_sen_plot1[V_sen_plot1$group=="N1=100\nN2=200",],
               aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + coord_cartesian( ylim = c(4,5.2)) +
  geom_hline(yintercept=5.352199, linetype="dashed", color="red") +
  labs(title="N1=100, N2=200", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot13

# sensitivity plot 14 (N1=200, N2=200)
sen_plot14 = ggplot(data=V_sen_plot1[V_sen_plot1$group=="N1=200\nN2=200",],
               aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + coord_cartesian( ylim = c(4,5.2)) +
  geom_hline(yintercept=5.352199, linetype="dashed", color="red") +
  labs(title="N1=200, N2=200", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot14

#install.packages("gridExtra")
library(gridExtra)
sen_plot1 = grid.arrange(sen_plot11, sen_plot12, sen_plot13, sen_plot14, nrow = 2)


###  Plot 2 series (test set same as Study 2)
V_sen_plot2 = data.frame(
  valfun = c(apply(V_Sim1[,1:16],2,mean, na.rm=T), apply(V_sen_a0.5[,21:36],2,mean, na.rm=T), apply(V_sen_a1[,21:36],2,mean, na.rm=T)),
  alpha = rep(c(0, 0.5, 1), each=16),
  group = rep(rep(c("N1=100\nN2=100", "N1=200\nN2=100", "N1=100\nN2=200", "N1=200\nN2=200"), 3), each=4),
  method = rep(c("1", "2", "3", "4"), 12))
V_sen_plot2$group <- factor(V_sen_plot2$group,levels = c("N1=100\nN2=100", "N1=200\nN2=100", "N1=100\nN2=200", "N1=200\nN2=200"))

# sensitivity plot 11 (N1=100, N2=100)
sen_plot21 = ggplot(data=V_sen_plot2[V_sen_plot2$group=="N1=100\nN2=100",],
                    aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + # coord_cartesian( ylim = c(4.5,5.8)) +
  labs(title="N1=100, N2=100", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot21

# sensitivity plot 12 (N1=200, N2=100)
sen_plot22 = ggplot(data=V_sen_plot2[V_sen_plot2$group=="N1=200\nN2=100",],
                    aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + #coord_cartesian(ylim = c(4.5,5.8)) +
  labs(title="N1=200, N2=100", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot22

# sensitivity plot 13 (N1=100, N2=200)
sen_plot23 = ggplot(data=V_sen_plot2[V_sen_plot2$group=="N1=100\nN2=200",],
                    aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + #coord_cartesian(ylim = c(4.5,5.8)) +
  labs(title="N1=100, N2=200", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot23

# sensitivity plot 14 (N1=200, N2=200)
sen_plot24 = ggplot(data=V_sen_plot2[V_sen_plot2$group=="N1=200\nN2=200",],
                    aes(x=alpha, y=valfun, group=method, color=method)) +
  geom_line() + #coord_cartesian(ylim = c(4.5,5.8)) +
  labs(title="N1=200, N2=200", x="alpha", y="Value function", fill="Method") +
  scale_color_discrete(labels=c("Naive 1", "Naive 2", "Integrated", "SMART"))
sen_plot24

# combine
sen_plot2 = grid.arrange(sen_plot21, sen_plot22, sen_plot23, sen_plot24, nrow = 2)



# # boxplot (not used -- hard to group and create grouped label)
# boxplot(V_Sim22[,1:12], at = c(2,3,4, 6,7,8, 9,10,11, 13,14,15))
# boxplot(V_Sim12[,1:12], las = 2, at = c(2,3,4, 5.5,6.5,7.5, 9,10,11, 12.5,13.5,14.5),
#         col=c("lightsteelblue1", "cornflowerblue","darkolivegreen3",  "lightsteelblue1", "cornflowerblue","darkolivegreen3",
#               "lightsteelblue1", "cornflowerblue","darkolivegreen3",  "lightsteelblue1", "cornflowerblue","darkolivegreen3"))


# spagetti plot for setting 5


















