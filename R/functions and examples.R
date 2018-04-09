############## Simulation of Allele specific Expression from Normal Distribution ##############
#' Function of simulation allele specific expression from normal distribution
#' @param n, number of samples, default value is 1000;
#' @param maf, minor allele frequency, default at 0.1;
#' @param mean.u, mean of U, default at 1;
#' @param sd.u, standard deviation of U, default at 1;
#' @param mean.v, mean of V, default at 1;
#' @param sd.v, standard deviation of V, default at 1;
#' @param mean.w, mean of W, default at 1;
#' @param sd.w, standard deviation of W, default at 1;
#' @param sd.e, standard deviation of measurement error, default at sqrt(2)/4;
#' @param rho, correlation between U and V, default at 0.5;
#' @param beta, causal effect size, default at 1;
#'
#' @import MASS
#' @details
#' Data are generated with model
#' \code{Z1} and \code{Z2} are genotypes, take values on 0 or 1.
#' \code{Z1 = rbinom(n, 1, maf)}
#' \code{Z1 + Z2 = 1} are heterozygous indiviudals, we can observe allele specific
#' expression \code{X1} and \code{X2}; \code{Z1 + Z2 = 0} or \code{2} are homozygous
#' individuals, we can only observe total expression \code{X}. For heterozygous
#' individuals,
#' \code{X1 = U + Z1*W1 + e1; T2 = U + Z2*W2 + e2; X = X1 + X2;}
#' and for homozygous individuals, we can observe
#' \code{X = 2U + Z1*W1 + Z2*W2 + e;}
#' where \code{W1} and \code{W2} are independent effect of genotypes.
#' \code{W1 = rnorm(n, mean.w, sd.w); W2 = rnorm(n, mean.w, sd.w);}
#' Outcome \code{Y = V + beta* (2U + Z1*W1 + Z2*W2)}, where \code{beta} is the causal
#' effect and \code{V} are confounder effect,
#' \code{(U, V) = mvrnorm(n, mu = c(mean.u, mean.v), Sigma = matrix(c(sd.u^2, rho*sd.u*sd.v, rho*sd.u*sd.v, sd.v^2),nrow = 2))}
#'
#' @return a list of allele specific expression: Y, X, X1, X2, Z1, Z2
#'
#' @examples
#' library(MASS)
#' Sim.out = sim_allelicMR()
#' str(Sim.out)
#' @export
sim_allelicMR = function(n = 1000, maf = 0.1, mean.u = 1, sd.u = 1, mean.v = 1, sd.v = 1, mean.w = 1, sd.w = 1, sd.e = sqrt(2)/4, rho = 0.5, beta = 1){
  theta = c(sigma.e = sd.e^2, sigma.v = sd.v^2, sigma.u = sd.u^2, sigma.w = sd.w^2, rho = rho, mu.u = mean.u, mu.v = mean.v, mu.w = mean.w);
  Z = rbinom(n, 2, maf)
  Z1 = rep(0,n); Z2 = rep(0,n)
  G1 = which(Z == 0)
  G2 = which(Z == 1)
  G3 = which(Z == 2)
  Z1[G2] = 1; Z1[G3] = 1; Z2[G3] = 1;
  N1 = length(G1); N2 = length(G2); N3 = length(G3)
  W1 = rnorm(n, mean.w, sd.w);
  W2 = rnorm(n, mean.w, sd.w);
  UV = mvrnorm(n, mu = c(mean.u, mean.v), Sigma = matrix(c(sd.u^2, rho*sd.u*sd.v, rho*sd.u*sd.v, sd.v^2),nrow = 2));
  U = UV[,1];
  V = UV[,2];
  R = 2*U + W1*Z1 + W2*Z2; R1 = U + W1*Z1; R2 = U + W2*Z2;
  X = R; X1 = R1; X2 = R2;
  X1[G2] = R1[G2] + rnorm(N2, 0, sd.e)
  X2[G2] = R2[G2] + rnorm(N2, 0, sd.e)
  X[G1] = R[G1] + rnorm(N1, 0, sd.e)
  X[G3] = R[G3] + rnorm(N3, 0, sd.e)
  X[G2] = X1[G2] + X2[G2]
  Y = beta*R + V;
  return(list(Y = Y, X = X, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2))
}

#' Function of Allele-Specific Mendelian Expression Estimation
#'
#' Estimate allele-specific Mendelian randomization parameters from expression
#'
#' @param exprs.list, List of expression: X1, X2 are allele specific expression level,
#' X is total expression level, Y is the outcome, Z1 and Z2 are genotypes.
#' @return Estimated parameters, including causal effect.
#'
#' @seealso
#' \code{\link{feasible.hin}},  \code{\link{VarCov}},
#' \code{\link{Likelihood.beta}},
#'
#' @importFrom alabama numDeriv MASS
#'
#' @examples
#' Sim.out = sim_allelicMR()
#' Est.AMR = AllelicMR_est(Sim.out$Y, Sim.out$X, Sim.out$X1, Sim.out$X2, Sim.out$Z1, Sim.out$Z2)
#' Est.AMR
AllelicMR_est = function(exprs.list){
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2; Y = exprs.list$Y;
  X1 = exprs.list$X1; X2 = exprs.list$X2;

  Xz = lm(X~Z)$fitted
  Yz = lm(Y~Z)$fitted
  beta.tsls = lm(Yz ~ Xz)$coef[2] # Two-stage least squares
  theta0 = feasible.hin(VarCov(beta.tsls, X, X1, X2, Y, Z1, Z2))
  OP = constrOptim.nl(c(theta0,beta.tsls), Likelihood.beta, gr = NULL, hin = hin.beta,
                      X = X, X1 = X1, X2 = X2, Y = Y, Z1 = Z1, Z2 = Z2,
                      control.outer = list(mu0 = 100000, itmax = 100000))
  beta.mle = OP$par[10]
  Par = OP$par[1:9]
  return(list(Par = Par, beta.asmr = beta.mle))
}


#' Function of feasible starting point
#'
#' Transform possible parameters to feasible area
#'
#'@param theta, possible parameters for likelihood estimation
#'
#'@return adjusted feasible parameters.
#'
feasible.hin = function(theta){
  sigma.e = theta[1];
  sigma.v = theta[2];
  sigma.uw = theta[3];
  sigma.uw0 = theta[4];
  sigma.uw1 = theta[5];
  rho.uv = theta[6];

  mu.u = theta[7]; mu.v = theta[8]; mu.w = theta[9];
  if(sigma.e<0){
    sigma.e = 10^{-3}
  }
  if(sigma.v<0){
    sigma.v = 10^{-3}
  }
  if(min(sigma.uw0, sigma.uw1)<0){
    if(sigma.uw0<0){
      sigma.uw0=10^{-3}
    }
    if(sigma.uw1<0){
      sigma.uw1=10^{-3}
    }
    sigma.uw = 10^{-4}
  }
  if(min(sigma.uw0, sigma.uw1)<sigma.uw){
    sigma.uw = min(sigma.uw0, sigma.uw1)/2
  }
  if(min(sigma.uw0, sigma.uw1) - 2*sigma.uw + sqrt(min(sigma.uw0, sigma.uw1) - sigma.uw) - abs(rho.uv/sqrt(sigma.v))<0){
    if(min(sigma.uw0, sigma.uw1) - 2*sigma.uw + sqrt(min(sigma.uw0, sigma.uw1) - sigma.uw)<0){
      sigma.uw = min(sigma.uw0, sigma.uw1)/2
    }
    rho.uv = abs(sqrt(sigma.v)*(min(sigma.uw0, sigma.uw1) - 2*sigma.uw + sqrt(min(sigma.uw0, sigma.uw1) - sigma.uw)))/2
  }
  theta.adj = c(sigma.e, sigma.v, sigma.uw, sigma.uw0, sigma.uw1, rho.uv, mu.u, mu.v, mu.w);
  return(theta.adj)
}

#' Function of feasible beta
#'
#' @param theta
hin.beta = function(theta, Sim.out){
  X1 = Sim.out$X1; X2 = Sim.out$X2; Y = Sim.out$Y;
  Z1 = Sim.out$Z1; Z2 = Sim.out$Z2;

  sigma.e = theta[1];
  sigma.v = theta[2];
  sigma.uw = theta[3];
  sigma.uw0 = theta[4];
  sigma.uw1 = theta[5];
  rho.uv = theta[6];

  mu.u = theta[7]; mu.v = theta[8]; mu.w = theta[9];
  Beta = theta[10];

  h = rep(NA, 7);
  h[1] = sigma.uw0
  h[2] = sigma.uw1
  h[3] = sigma.uw1 - sigma.uw
  h[4] = sigma.uw0 - sigma.uw
  if(min(sigma.uw1, sigma.uw0)-sigma.uw<0){
    h[5]=-1
  }else{
    if(sigma.v<0){
      h[5] = -1
    }else{
      h[5] = min(sigma.uw0, sigma.uw1) - 2*sigma.uw + sqrt(min(sigma.uw0, sigma.uw1) - sigma.uw) -
        abs(rho.uv/sqrt(sigma.v))
    }
  }
  h[6] = sigma.e
  h[7] = sigma.v
  return(h)
}


#' Function for calculating log-likelihood
#'
#' @param theta, all parameters;
#' @param exprs.list, List of outcome Y, allele specific expression X1, X2, and total expression X, and genotype Z1, Z3;
#'
#' @return
#'
Likelihood.beta = function(theta, exprs.list){
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2; Y = exprs.list$Y;
  X1 = exprs.list$X1; X2 = exprs.list$X2;

  sigma.e = theta[1];
  sigma.v = theta[2];
  sigma.uw = theta[3];
  sigma.uw0 = theta[4];
  sigma.uw1 = theta[5];
  rho.uv = theta[6];

  mu.u = theta[7]; mu.v = theta[8]; mu.w = theta[9];
  Beta = theta[10];

  Z = Z1+Z2;
  G1 = which(Z == 0)
  G2 = which(Z == 1)
  G3 = which(Z == 2)

  N1 = length(G1)
  N2 = length(G2)
  N3 = length(G3)

  Sigma1 = matrix(c(sigma.v + Beta^2*sigma.e, 2*rho.uv-Beta*sigma.e,
                    2*rho.uv-Beta*sigma.e, 2*sigma.uw0 + sigma.e), nrow = 2)
  if(det(Sigma1)>0){
    R1 = Y[G1]- Beta*X[G1]-mu.v
    R2 = X[G1]-2*mu.u
    iSigma1 = matrix(c(Sigma1[2,2], -Sigma1[1,2],
                       -Sigma1[2,1], Sigma1[1,1]), nrow = 2)/det(Sigma1)
    L1 = is.zero(N1)*(  N1* log(det(Sigma1))+ sum(R1^2)*iSigma1[1,1] + sum(R2^2)*iSigma1[2,2] + 2*iSigma1[1,2]*sum(R1*R2)   )
  }else{
    L1 = 10^8
  }
  Sigma2 = matrix(c(sigma.v+ 2*Beta^2*sigma.e, 2*rho.uv - 2*Beta*sigma.e, 0,
                    2*rho.uv- 2*Beta*sigma.e,  sigma.uw0 + sigma.uw1 + 2*sigma.e, sigma.uw1-sigma.uw0,
                    0, sigma.uw1-sigma.uw0, sigma.uw1 + sigma.uw0 - 4*sigma.uw + 2*sigma.e), nrow = 3)
  if(det(Sigma2)>0){
    iSigma2 = matrix(c(det(Sigma2[2:3, 2:3]), -det(Sigma2[2:3,c(1,3)]),
                       det(Sigma2[2:3, 1:2]), -det(Sigma2[c(1,3), 2:3]),
                       det(Sigma2[c(1,3), c(1,3)]), -det(Sigma2[c(1,3), 1:2]),
                       det(Sigma2[1:2, 2:3]), -det(Sigma2[1:2,c(1,3)]),
                       det(Sigma2[1:2, 1:2]) ) , nrow = 3)/det(Sigma2)
    R1 = Y[G2] - Beta*X[G2] - mu.v
    R2 = X[G2] - 2*mu.u - mu.w
    R3 = X1[G2] - X2[G2] - mu.w
    L2 = is.zero(N2)*(N2*log(det(Sigma2))+ iSigma2[1,1]*sum(R1^2) + iSigma2[2,2]*sum(R2^2) +
                        iSigma2[3,3]*sum(R3^2) + 2*iSigma2[1,2]*sum(R1*R2) + 2*iSigma2[1,3]*sum(R1*R3) +
                        2*iSigma2[2,3]*sum(R2*R3))
  }else{
    L2 = 10^8
  }

  Sigma3 = matrix(c(sigma.v+Beta^2*sigma.e, 2*rho.uv-Beta*sigma.e,
                    2*rho.uv-Beta*sigma.e, 2*sigma.uw1 + sigma.e), nrow = 2)
  if(det(Sigma3)>0){
    iSigma3 = matrix(c(Sigma3[2,2], - Sigma3[2,1],
                       -Sigma3[1,2], Sigma3[1,1]), nrow = 2)/det(Sigma3)
    R1 = Y[G3] - Beta*X[G3] - mu.v
    R2 = X[G3] - 2*mu.u - 2*mu.w
    L3 = is.zero(N3)*(N3*log(det(Sigma3)) + sum(R1^2)*iSigma3[1,1] + sum(R2^2)*iSigma3[2,2] +
                        2*iSigma3[1,2]*sum(R1*R2) )
  }else{
    L3 = 10^8
  }

  L = L1 + L2 + L3
  return(L)
}

#' Function for calculating log-likelihood, stage 1
#'
#' This function calculate the log-likelihood of first stage
#' \code{X = 2*U + Z1*W1 + Z2*W2 + e} for homozygous individuals and
#' allele specific expression \code{X1 = U + Z1*W1 + e1;}
#' \code{X2 = Z2*W2 + e2} for heterozygous individuals.
#'
#' @param theta, first stage parameters, including sigma.e, sigma.u, sigma.w, mu.u and mu.w.
#' @param exprs.list, expression list, but Y is not used in this function;
#'
#'
#' @return Log-likelihoodof the first stage.
#'
Likelihood.1 = function(theta, exprs.list){
  X = exprs.list$X; X1 = exprs.list$X1; X2 = exprs.list$X2;
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2;

  sigma.e = theta[1];
  sigma.u = theta[2];
  sigma.w = theta[3];
  mu.u = theta[4]; mu.w = theta[5];

  Z = Z1+Z2;
  G1 = which(Z == 0)
  G2 = which(Z == 1)
  G3 = which(Z == 2)

  N1 = length(G1)
  N2 = length(G2)
  N3 = length(G3)

  Sigma1 = 4*sigma.u + sigma.e
  if(is.na(Sigma1)|| Sigma1<=0){
    L1 = 10^8
  }else{
    R = X[G1]-2*mu.u;
    L1 = is.zero(N1)*(  N1*log(Sigma1)+ sum(R^2)/Sigma1  )
  }

  Sigma2 = matrix(c(sigma.u+sigma.w+sigma.e, sigma.u, sigma.u, sigma.u+sigma.e), nrow = 2);
  if((det(Sigma2)<=0)||is.na(det(Sigma2))){
    L2 = 10^8
  }else{
    iSigma2 = matrix(c(Sigma2[2,2], -Sigma2[1,2],
                       -Sigma2[2,1], Sigma2[1,1]), nrow = 2)/det(Sigma2)
    R1 = X1[G2] - mu.u - mu.w
    R2 = X2[G2] - mu.u
    L2 = is.zero(N2)*(  N2* log(det(Sigma2)) + sum(R1^2)*iSigma2[1,1] + sum(R2^2)*iSigma2[2,2] + 2*iSigma2[1,2]*sum(R1*R2)   )
  }

  Sigma3 = 4*sigma.u + 2*sigma.w + sigma.e
  if(is.na(Sigma3)|| Sigma3<=0){
    L1 = 10^8
  }else{
    R = X[G3] - 2*mu.u - 2*mu.w;
    L3 = is.zero(N3)*(  N3*log(Sigma3)+ sum(R^2)/Sigma3  )
  }
  L = L1 + L2 + L3
  return(L)
}

#' Function for calculating log-likelihood of second stage
#' This function take estimated parameters from first stage, then estimate the
#' second stage parameters sigma.v, rho, mu.v and beta by maxizing full log-likelihood.
#' @param theta, parameters for second stage;
#' @param theta.x, estimated parameters from first stage;
#' @param exprs.list, List of expressions: total expression X, allele specific expression X1 and X2,
#' genotype Z1, Z2, and outcome Y;
#'
#' @return Log-likelihood
Likelihood.2 = function(theta, theta.x, exprs.list){
  X = exprs.list$X; X1 = exprs.list$X1; X2 = exprs.list$X2;
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2;
  Y = exprs.list$Y;

  sigma.e = theta.x[1];
  sigma.u = theta.x[2];
  sigma.w = theta.x[3];
  mu.u = theta.x[4]; mu.w = theta.x[5];

  sigma.v = theta[1];
  rho = theta[2];
  mu.v = theta[3];
  Beta = theta[4];

  Z = Z1+Z2;
  G1 = which(Z == 0)
  G2 = which(Z == 1)
  G3 = which(Z == 2)

  N1 = length(G1)
  N2 = length(G2)
  N3 = length(G3)

  Sigma1 = matrix(c(sigma.v + Beta^2*sigma.e, 2*rho*sqrt(sigma.u*sigma.v) - Beta*sigma.e,
                    2*rho*sqrt(sigma.u*sigma.v) - Beta*sigma.e, 4*sigma.u + sigma.e), nrow = 2)
  if((det(Sigma1)<=0)||is.na(det(Sigma1))){
    L1 = 10^8
  }else{
    R1 = Y[G1]- Beta*X[G1]-mu.v
    R2 = X[G1]-2*mu.u
    iSigma1 = matrix(c(Sigma1[2,2], -Sigma1[1,2],
                       -Sigma1[2,1], Sigma1[1,1]), nrow = 2)/det(Sigma1)
    L1 = is.zero(N1)*(  N1* log(det(Sigma1))+ sum(R1^2)*iSigma1[1,1] + sum(R2^2)*iSigma1[2,2] + 2*iSigma1[1,2]*sum(R1*R2)   )
  }
  Sigma2 = matrix(c(sigma.v+ 2*Beta^2*sigma.e, 2*rho*sqrt(sigma.u)*sqrt(sigma.v) - 2*Beta*sigma.e, 0,
                    2*rho*sqrt(sigma.u)*sqrt(sigma.v)- 2*Beta*sigma.e,  4*sigma.u + sigma.w + 2*sigma.e, sigma.w,
                    0, sigma.w, sigma.w + 2*sigma.e), nrow = 3)
  if((det(Sigma2)<=0)||is.na(det(Sigma2))){
    L2 = 10^8
  }else{
    iSigma2 = matrix(c(det(Sigma2[2:3, 2:3]), -det(Sigma2[2:3,c(1,3)]),
                       det(Sigma2[2:3, 1:2]), -det(Sigma2[c(1,3), 2:3]),
                       det(Sigma2[c(1,3), c(1,3)]), -det(Sigma2[c(1,3), 1:2]),
                       det(Sigma2[1:2, 2:3]), -det(Sigma2[1:2,c(1,3)]),
                       det(Sigma2[1:2, 1:2]) ) , nrow = 3)/det(Sigma2)
    R1 = Y[G2] - Beta*X[G2] - mu.v
    R2 = X[G2] - 2*mu.u - mu.w
    R3 = X1[G2] - X2[G2] - mu.w
    L2 = is.zero(N2)*(   N2*log(det(Sigma2))+ iSigma2[1,1]*sum(R1^2) + iSigma2[2,2]*sum(R2^2) + iSigma2[3,3]*sum(R3^2) + 2*iSigma2[1,2]*sum(R1*R2) + 2*iSigma2[1,3]*sum(R1*R3) + 2*iSigma2[2,3]*sum(R2*R3)   )

  }

  Sigma3 = matrix(c(sigma.v + Beta^2*sigma.e, 2*rho*sqrt(sigma.u)*sqrt(sigma.v) - Beta*sigma.e,
                    2*rho*sqrt(sigma.u)*sqrt(sigma.v) - Beta*sigma.e, 4*sigma.u + 2*sigma.w + sigma.e), nrow = 2)
  if(det(Sigma3)<=0 || is.na(det(Sigma3))){
    L3 = 10^8
  }else{
    iSigma3 = matrix(c(Sigma3[2,2], - Sigma3[2,1],
                       -Sigma3[1,2], Sigma3[1,1]), nrow = 2)/det(Sigma3)
    R1 = Y[G3] - Beta*X[G3] - mu.v
    R2 = X[G3] - 2*mu.u - 2*mu.w
    L3 = is.zero(N3)*( N3*log(det(Sigma3)) + sum(R1^2)*iSigma3[1,1] + sum(R2^2)*iSigma3[2,2] + 2*iSigma3[1,2]*sum(R1*R2) )
  }

  L = L1 + L2 + L3
  return(L)
}


#' Function for calculating possible starting point for all parameters
#'
#' Calculate possible starting point for all parameters
#'
#' @param Beta, possible estimated beta;
#' @param exprs.list, List of outcome Y, allele specific expression X1, X2, and total expression X, and genotype Z1, Z3;
#'
#' @return starting points theta;
#'
#' @examples
#' Sim.out = sim_allelic()
#' Xz = lm(X~Z)$fitted
#' Yz = lm(Y~Z)$fitted
#' beta.tsls = lm(Yz ~ Xz)$coef[2] # Two-stage least squares
#' theta = VarCov(beta.tsls, Sim.out)
VarCov = function(Beta, exprs.list){
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2; Y = exprs.list$Y;
  X1 = exprs.list$X1; X2 = exprs.list$X2;

  Z = Z1+Z2; # Genotype
  G1 = which(Z == 0);  # Number of individuals in each genotypes
  G2 = which(Z == 1);
  G3 = which(Z == 2);

  S1 = cov(cbind(Y[G1] - Beta *X[G1], X[G1]));
  S2 = cov(cbind(Y[G2] - Beta *X[G2], X[G2], X1[G2] - X2[G2]));
  S3 = cov(cbind(Y[G3] - Beta *X[G3], X[G3]));
  M1 = colMeans(cbind(Y[G1] - Beta *X[G1], X[G1]));
  M2 = colMeans(cbind(Y[G2] - Beta *X[G2], X[G2], X1[G2] - X2[G2]));
  M3 = colMeans(cbind(Y[G3] - Beta *X[G3], X[G3]));
  MS = list(M1, S1, M2, S2, M3, S3);
  MY = c(rep(MS[[1]], N1), rep(MS[[3]], N2), rep(MS[[5]], N3))
  M =cbind(MY, matrix( c(rep( c(1,0,0,0,2,0), N1), rep( c(1,0,0,0,2,1,0,0,1), N2), rep( c(1,0,0,0,2,2), N3)),
                       ncol =  3 , byrow = TRUE ))
  colnames(M) = c('Y', 'mu.v', 'mu.u', 'mu.w')
  M = as.data.frame(M)
  mu = lm(Y ~ -1 + mu.u + mu.v + mu.w, data = M)$coef
  SY = c(rep(MS[[2]][-2], N1), rep(MS[[4]][c(-3,-4,-7,-8)], N2), rep(MS[[6]][-2], N3))
  S = cbind(SY, matrix(c( rep( c(Beta^2,1,0,0,0,0, -Beta,0,0,0,0,2, 1,0,0,2,0,0), N1),
                          rep( c(2*Beta^2,1,0,0,0,0, -Beta,0,0,0,0,2, 2,0,0,1,1,0, 0,0,0,-1,1,0, 2,0,-4,1,1,0), N2),
                          rep(c(Beta^2,1,0,0,0,0, -Beta,0,0,0,0,2, 1,0,0,0,2,0), N3) ), ncol = 6, byrow = T))
  colnames(S) = c('Y', 'sigma.e', 'sigma.v', 'sigma.uw', 'sigma.uw0', 'sigma.uw1', 'rho.uv')
  S = as.data.frame(S)
  sigma = lm(Y ~ sigma.e + sigma.v + sigma.uw + sigma.uw0 + sigma.uw1 + rho.uv, data = S)$coef[-1]
  theta = c(sigma, mu)
  return(theta)
}

#' Function for calculating possible starting point for first stage parameters
#'
#' Calculate possible starting point for first stage parameters
#'
#' @param Beta, possible estimated beta;
#' @param exprs.list, List of outcome Y, allele specific expression X1, X2, and total expression X, and genotype Z1, Z3;
#'
#' @return starting points theta;
#'
#' @examples
#' Sim.out = sim_allelic()
#' theta = VarCov.1(Sim.out)
VarCov.1 = function(exprs.list){
  X = exprs.list$X; X1 = exprs.list$X1; X2 = exprs.list$X2;
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2;

  Z = Z1+Z2;
  G1 = which(Z == 0);
  G2 = which(Z == 1);
  G3 = which(Z == 2);

  S1 = var(X[G1]);
  S2 = cov(cbind(X1[G2], X2[G2]));
  S3 = var(X[G3]);
  M1 = mean(X[G1]);
  M2 = colMeans(cbind(X1[G2], X2[G2]));
  M3 = mean(X[G3]);
  MS = list(M1, S1, M2, S2, M3, S3);
  MY = c(rep(MS[[1]], N1), rep(MS[[3]], N2), rep(MS[[5]], N3))
  M =cbind(MY, matrix( c(rep( c(2, 0), N1), rep( c(1,1,0,1), N2), rep( c(2,2), N3)) , ncol =  2 , byrow = TRUE ))
  colnames(M) = c('Y', 'mu.u', 'mu.w')
  M = as.data.frame(M)
  mu = lm(Y ~ -1 + mu.u + mu.w, data = M)$coef
  SY = c(rep(MS[[2]], N1), rep(MS[[4]][-2], N2), rep(MS[[6]], N3))
  S = cbind(SY, matrix(c( rep( c(1, 4, 0), N1), rep( c(1,1,1, 0,1,0, 1,1,0), N2), rep(c(1,4,2), N3) ), ncol = 3, byrow = T))
  colnames(S) = c('Y', 'sigma.e', 'sigma.u', 'sigma.w')
  S = as.data.frame(S)
  sigma = lm(Y ~ sigma.e + sigma.u + sigma.w, data = S)$coef[-1]
  sigma[1:3][sigma[1:3]<0] = 0.001

  theta = c(sigma, mu)
  return(theta)
}

#' Function of iszero
#'
#' is.zero
#'
#' @param x
#' @return as.numeric(x!=0)
is.zero = function(x){
  if(x==0)
    a = 0
  if(x!=0)
    a = 1
  return(a)
}


#' Function of expit
#'
#' expit
#'
#' @param x numeric
#' @return expit(x)
#'
expit = function(x){
  a = exp(x)/(1+exp(x))
  return(a)
}


#' Function of two-stage maximize likelihood estimation
#'
#' Two-stage MLE
#'
#' @param exprs.list List of outcome Y, allele specific expression X1, X2, and total expression X, and genotype Z1, Z3;
#' @return causal effect and other parameters
#' @import MASS
tsMLE_est = function(exprs.list){
  Z1 = exprs.list$Z1; Z2 = exprs.list$Z2; Y = exprs.list$Y;
  X1 = exprs.list$X1; X2 = exprs.list$X2;
  Xz = lm(X~Z)$fitted
  Yz = lm(Y~Z)$fitted
  beta.tsls = lm(Yz ~ Xz)$coef[2] # Two-stage least squares

  # 2 stage likelihood
  ui1 = matrix(rep(0, 15), nrow = 3)
  ui1[1,1] = 1; ui1[2,2] = 1; ui1[3,3] = 1;
  ci1 = c(0, 0, 0)

  ui2 = matrix(rep(0, 12), nrow = 3)
  ui2[1,1] = 1; ui2[2,2] = 1; ui2[3,2] = -1;
  ci2 = c(0, -1, -1)

  thetax0 = VarCov.1(X, X1, X2, Z1, Z2)
  OP1 = constrOptim(thetax0, Likelihood.1, gr = NULL, ui = ui1, ci = ci1, X = X,
                    X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, control = list(maxit = 100000))
  thetax = OP1$par
  OP2 = constrOptim(c(theta0[c(2,5,7)], beta.tsls), Likelihood.2, gr = NULL, ui = ui2, ci = ci2, theta.x = thetax,
                    X =X, X1 = X1, X2 = X2, Y = Y, Z1 = Z1, Z2 = Z2, control = list(maxit = 100000))
  beta.tsmle = OP2$par[4];
  par = c(OP1$par[1], OP2$par[1], OP1$par[2:3], OP2$par[2], OP1$par[4], OP2$par[3], OP1$par[5]);
  return(Par = par, beta.tsmle = beta.tsmle)
}

