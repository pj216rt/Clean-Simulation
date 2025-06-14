library(tidyverse)
library(splines) #for working with b splines
library(mvtnorm) #multivariate normal
library(spam)
library(MASS)
library(MCMCpack) #for the inverse wishart distribution
require(gridExtra)
library(bayesplot)
library(abind)

##univariate functional portion
beta.gen.fun.uni <- function(grid.size, knots, parameters, deg){
  x <- seq(0, grid.size, length.out = grid.size)
  
  #b spline basis matrix
  bs_basis <- bs(x, degree = deg)  
  
  beta_list <- vector("list", 3)
  
  #going to loop over each dimension
  for (dim in 1:1) {
    random_coeffs_matrix <- matrix(rnorm(ncol(bs_basis)*parameters), 
                                   nrow = ncol(bs_basis), ncol = parameters)
    
    #compute the beta matrix for dimension dim
    beta_matrix <- bs_basis %*% random_coeffs_matrix
    
    #store this in the list
    beta_list[[dim]] <- t(beta_matrix)
  }
  
  #combine all three matrices column wise
  coefs <- do.call(cbind, beta_list)
  
  return(coefs)
}

#generate smooth beta functions for the trivariate functional data
beta.gen.fun.tri <- function(grid.size, knots, parameters, deg){
  x <- seq(0, grid.size, length.out = grid.size)
  
  #b spline basis matrix
  bs_basis <- bs(x, degree = deg)  
  
  beta_list <- vector("list", 3)
  
  #going to loop over each dimension
  for (dim in 1:3) {
    random_coeffs_matrix <- matrix(rnorm(ncol(bs_basis)*parameters), 
                                   nrow = ncol(bs_basis), ncol = parameters)
    
    #compute the beta matrix for dimension dim
    beta_matrix <- bs_basis %*% random_coeffs_matrix
    
    #store this in the list
    beta_list[[dim]] <- t(beta_matrix)
  }
  
  #combine all three matrices column wise
  coefs <- do.call(cbind, beta_list)
  
  return(coefs)
}

#function to plot and save simulated cofficient plots
generated.coef.plot <- function(generated.beta.values){
  ##function accepta generate beta values
  
  #get the number of grid time points in each of the 3 dimensions, prepare for 
  #plotting them by getting dimension labels ready
  cols_per_dim <- ncol(generated.beta.values) /3
  dim_labels <- c("X", "Y", "Z")
  
  df <- as.data.frame(generated.beta.values) %>%
    mutate(param = row_number()) %>%
    pivot_longer(
      cols = -param,
      names_to = "colname",
      values_to = "value"
    ) %>%
    mutate(
      col_index = as.integer(gsub("V", "", colname)),
      dim = dim_labels[ceiling(col_index / cols_per_dim)],
      time = (col_index - 1) %% cols_per_dim + 1
    )
  
  coef.plot <- ggplot(df, aes(x = time, y = value, group = param, color = factor(param))) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    facet_wrap(~dim, nrow = 1, strip.position = "top") +
    labs(
      x = "Time",
      y = "Coefficient Value",
      title = "Functional Regression Coefficient Trajectories",
      subtitle = "Colored by Parameter",
      color = "Parameter"
    ) +
    ##the theme minimal base sets all fonbt size to 14pts
    ##the strip text controls the facet labels
    ##element_text for the plot_title centers and makes bold the title
    ##also centers the subtitle, but does not bold it
    ##also tweaks the spacing between the facets
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 13),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.spacing = unit(1, "lines")
    ) +
    scale_color_viridis_d(option = "D")
  
  return(coef.plot)
}


#generate univariate data
gen.1d.data <- function(num.subj = 50, num.visits = 5,
                        beta.vals, var.x = 1, var.z=1){
  #accepts num.subj = the number of subjects
  #num.visits = the number of repeated visits/repeated measurements
  #beta.vals = the matrix of true coefficient curves
  #var.x = the variance of the covariates
  #var.z = the variance of the random effects
  
  #error message
  #make sure you supply beta vals
  if (missing(beta.vals)) {
    message("Beta values not supplied.")
  }
  
  #get the grid size and number of parameters from the beta coefficient matrix
  grid.size <- ncol(beta.vals)
  n.params <- nrow(beta.vals)
  
  #compute total observations
  tot.obs <- num.subj*num.visits
  
  #generate the X design matrix
  #each subject receives their own set of covariates
  X <- matrix(rnorm(n=(num.subj*n.params), sd = sqrt(var.x)), 
              nrow = num.subj, ncol = n.params) %>%
    as.data.frame() %>%
    mutate(subj = row_number())
  
  X$subj = as.factor(X$subj)
  
  #X1 replicates each subject's covariates num.visits times to put it in long format
  X1 <- X[rep(seq_len(nrow(X)), each = num.visits), ]
  
  #create design matrix dropping intercept and the column subj
  X.des <- model.matrix(~. -1 -subj, data = X1)
  
  #generate fixed effects by multiplying design matrix by beta values
  fixef = as.matrix(X.des) %*% as.matrix(beta.vals)
  
  #random effects
  
  #create random effect design matrix using subject IDs
  Z.des = model.matrix( ~ 0 + subj + (-1):subj, data = X1)
  
  #creates a matrix where each row is a subject, and each column corresponds to
  #a point on the grid.  Basically each row is a set of different random deviations
  #over the grid
  subj.ranef <- matrix(rnorm(n=(num.subj*grid.size), sd=sqrt(var.z)), 
                       nrow = num.subj, ncol = grid.size)
  ranef <- Z.des %*% subj.ranef
  
  #level 1 residuals, adding Gaussian noise for each observation
  eps <- matrix(rnorm(n = tot.obs * grid.size), nrow = tot.obs, ncol = grid.size)
  
  #add random effects and epsilon errors to the observed data
  #Yij.true = a true latent trajectory??
  Yij.true <- fixef + ranef
  
  #Yij.obs = add level 1 error
  Yij.obs <- fixef + ranef + eps
  
  #return "observed" data and other stuff in a named list
  output.list <- list()
  output.list$obs <- Yij.obs
  output.list$true <- Yij.true
  output.list$x_design <- X.des
  output.list$z_design <- Z.des
  output.list$raw.data <- X1
  return(output.list)
}

test1 <- gen.1d.data(beta.vals = test)


#generate 3 dimensional data, assuming a balanced simulation design
gen.3d.data <- function(num.subj = 50, num.visits = 5,
                        beta.vals, var.x = 1, var.z=1){
  #accepts num.subj = the number of subjects
  #num.visits = the number of repeated visits/repeated measurements
  #beta.vals = the matrix of true coefficient curves
  #var.x = the variance of the covariates
  #var.z = the variance of the random effects
  
  #error message
  #make sure you supply beta vals
  if (missing(beta.vals)) {
    message("Beta values not supplied.")
  }
  
  #get the grid size and number of parameters from the beta coefficient matrix
  grid.size.3d <- ncol(beta.vals)
  n.params <-nrow(beta.vals)
  
  #compute total observations
  tot.obs <- num.subj*num.visits
  
  #generate the X design matrix
  #each subject receives their own set of covariates
  X <- matrix(rnorm(n=(num.subj*n.params), sd = sqrt(var.x)), 
              nrow = num.subj, ncol = n.params) %>%
    as.data.frame() %>%
    mutate(subj = row_number())
  
  X$subj = as.factor(X$subj)
  
  #X1 replicates each subject's covariates num.visits times to put it in long format
  X1 <- X[rep(seq_len(nrow(X)), each = num.visits), ]
  
  #create design matrix dropping intercept and the column subj
  X.des <- model.matrix(~. -1 -subj, data = X1)
  
  #generate fixed effects by multiplying design matrix by beta values
  fixef = as.matrix(X.des) %*% as.matrix(beta.vals)
  
  #random effects
  
  #create random effect design matrix using subject IDs
  Z.des = model.matrix( ~ 0 + subj + (-1):subj, data = X1)
  
  #creates a matrix where each row is a subject, and each column corresponds to
  #a point on the grid.  Basically each row is a set of different random deviations
  #over the grid
  subj.ranef <- matrix(rnorm(n=(num.subj*grid.size.3d), sd=sqrt(var.z)), 
                       nrow = num.subj, ncol = grid.size.3d)
  ranef <- Z.des %*% subj.ranef
  
  #level 1 residuals, adding Gaussian noise for each observation
  eps <- matrix(rnorm(n = tot.obs * grid.size.3d), nrow = tot.obs, ncol = grid.size.3d)
  
  #add random effects and epsilon errors to the observed data
  #Yij.true = a true latent trajectory??
  Yij.true <- fixef + ranef
  
  #Yij.obs = add level 1 error
  Yij.obs <- fixef + ranef + eps
  
  #return "observed" data and other stuff in a named list
  output.list <- list()
  output.list$obs <- Yij.obs
  output.list$true <- Yij.true
  output.list$x_design <- X.des
  output.list$z_design <- Z.des
  output.list$raw.data <- X1
  return(output.list)
  
}

#goldsmiths uni variate Gibbs sampler
goldsmiths.one.d.gibbs <- function(Y, fixed_design.mat, random_design.mat, 
                                   Kt, N.iter = 1000, N.burn = 200, alpha = .1){
  ##This Gibbs sampler accepts a Y matrix composed of concatenated observations
  ##for one dimension, a fixed effect design matrix, a random effect design
  ##matrix, the number of knots to be fit for each dimension, as well as overrideable 
  ##defaults for the number of iterations and number of burn in iterations.  One can
  ##also change the alpha term.
  
  set.seed(1)
  
  ## fixed and random effect design matrices
  # W.des = model.matrix( fixef.form, data = data)
  # W.des <- W.des[rep(seq_len(nrow(W.des)), length(id)), ]
  # Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  # W.des = as.spam(W.des)
  # print(dim(W.des))
  # Z.des = as.spam(Z.des)
  ## fixed and random effect design matrices
  #W.des = model.matrix( fixef.form, data = data)
  #Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  W.des = as.spam(fixed_design.mat)
  Z.des = as.spam(random_design.mat)
  
  I = dim(Z.des)[2]
  D = dim(Y)[2]
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  
  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))
  
  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = W.des[firstobs,]
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  IIP = kronecker(diag(1, I, I), P.mat)
  WIk = kronecker(Wi, diag(1, Kt, Kt))
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP
  
  
  ## initial estimation and hyperparameter choice
  vec.bz = solve(kronecker(t(Z.des)%*% Z.des, t(Theta) %*% Theta)) %*% t(kronecker(Z.des, Theta)) %*% Y.vec
  bz = matrix(vec.bz, nrow = Kt, ncol = I)
  
  w.temp = kronecker(t(Wi), diag(1, Kt, Kt))
  vec.bw = solve(tWIW) %*% tWI %*% vec.bz
  bw = matrix(vec.bw, Kt, p)
  
  Yhat = as.matrix(Z.des %*% t(bz) %*% t(Theta))
  varhat = var(as.vector(Y - Yhat))
  
  Psi = diag(varhat*IJ, D, D)
  v = IJ
  inv.sig = solve(Psi/v)
  
  Az = I*Kt / 2
  Bz = sum(diag((t(bz) - Wi %*% t(bw)) %*% P.mat %*% t(t(bz) - Wi %*% t(bw))))
  
  Aw = Kt / 2
  Bw = sapply(1:p, function(u) max(1, sum(diag( t(bw[,u]) %*% P.mat %*% (bw[,u])))))
  
  
  ## matrices to store within-iteration estimates 
  BW = array(NA, c(Kt, p, N.iter))
  BW[,,1] = bw
  BZ = array(NA, c(Kt, I, N.iter))
  BZ[,,1] = bz
  INV.SIG = array(NA, c(D, D, N.iter))
  INV.SIG[,,1] = inv.sig 
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW[1,] = lambda.bw = Aw/Bw
  LAMBDA.BZ = rep(NA, N.iter)
  LAMBDA.BZ[1] = lambda.ranef = Az/Bz
  
  y.post = array(NA, dim = c(IJ, D, (N.iter - N.burn)))
  
  cat("Beginning Sampler \n")
  pb <- txtProgressBar(min = 0, max = N.iter, initial = 0, style = 3)
  for(i in 1:N.iter){
    setTxtProgressBar(pb,i)
    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################
    
    for(subj in 1:length(unique(SUBJ))){
      
      t.designmat.Z = t(kronecker(rep(1, Ji[subj]), Theta))
      
      sigma = solve(t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% t(t.designmat.Z) + 
                      (lambda.ranef) * P.mat)
      mu = sigma %*% (t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),]))) +
                        ((lambda.ranef) * P.mat) %*% bw %*% t(Wi[subj,]))
      
      bz[,subj] = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = 1)
    }
    ranef.cur = Z.des %*% t(bz) %*% t(Theta)
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    
    sigma = solve( kronecker(diag(lambda.bw), P.mat) + ((lambda.ranef) * tWIW) )
    mu = sigma %*% ((lambda.ranef) * tWI %*% as.vector(bz))
    
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = p)
    
    beta.cur = t(bw) %*% t(Theta)
    
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    resid.cur = Y - ranef.cur
    inv.sig = solve(riwish(v + IJ, Psi + t(resid.cur) %*% resid.cur))
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      b.post = Bw[term] + 1/2 * bw[,term] %*% P.mat %*% bw[,term]
      lambda.bw[term] = rgamma(1, a.post, b.post)
    }
    
    ## lambda for random effects
    a.post = Az + I*Kt/2
    b.post = Bz + .5 * sum(sapply(1:I, function(u) (t(bz[,u]) - Wi[u,] %*% t(bw)) %*% P.mat %*% t(t(bz[,u]) - Wi[u,] %*% t(bw)) ))
    lambda.ranef = rgamma(1, a.post, b.post)
    
    ###############################################################
    ## save this iteration's parameters
    ###############################################################
    
    BW[,,i] = as.matrix(bw)
    BZ[,,i] = as.matrix(bz)
    
    INV.SIG[,,i] = inv.sig
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.BZ[i] = lambda.ranef
    
    if(i > N.burn){
      y.post[,,i - N.burn] = ranef.cur
    }
    
  }
  close(pb)
  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################
  
  ## main effects
  beta.post = array(NA, dim = c(p, D, (N.iter - N.burn)))
  for(n in 1:(N.iter - N.burn)){
    beta.post[,,n] = t(BW[,, n + N.burn]) %*% t(Theta)
  }
  
  beta.pm = apply(beta.post, c(1,2), mean)
  beta.LB = apply(beta.post, c(1,2), quantile, c(.025))
  beta.UB = apply(beta.post, c(1,2), quantile, c(.975))
  
  ## random effects
  b.pm = matrix(NA, nrow = I, ncol = D)
  for(i in 1:I){
    b.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
    for(n in 1:(N.iter - N.burn)){
      b.post[n,] = BZ[,i, n + N.burn] %*% t(Theta)
    }
    b.pm[i,] = apply(b.post, 2, mean)
  }
  
  ## covariance matrix
  sig.pm = solve(apply(INV.SIG, c(1,2), mean))
  
  ## export fitted values
  ranef.pm = Z.des %*% b.pm
  Yhat = apply(y.post, c(1,2), mean)
  
  #extracing posterior samples
  #the fixed and random effect spline coefficients
  BW.samples = BW[,, (N.burn + 1):N.iter]
  BZ.samples = BZ[,, (N.burn + 1):N.iter]
  
  #the lambda samples
  lambda_bw_samples <- LAMBDA.BW[(N.burn+1):N.iter, ] #(N.iter-N.burn)*p
  lambda_bz_samples <- LAMBDA.BZ[(N.burn+1):N.iter]  #(N.iter-N.burn)
  inv_sig_samples   <- INV.SIG[,, (N.burn+1):N.iter]  #D*D*(N.iter-N.burn)
  
  #convert precisions, to variances
  sigma2_w_samples <- 1 / lambda_bw_samples #dim should be (800*p)
  sigma2_z_samples <- 1 / lambda_bz_samples
  
  #convert precision matrix to covariance matrix
  #Sigma_samples <- apply(inv_sig_samples, 3, solve)
  #abind is a generalization of rbind and cbind.  Accepts a sequence of vectors, matrices,
  #and arrays and returns an array of the same dimension
  #Sigma_samples <- abind::abind(lapply((N.burn+1):N.iter, function(i) solve(inv_sig_samples[,,i])), along = 3)

  S <- dim(inv_sig_samples)[3]
  Sigma_samples <- abind::abind(
    lapply(1:S, function(i) solve(inv_sig_samples[,,i])),
    along = 3
  )
  
  #return things
  ret = list(beta.pm, beta.LB, beta.UB, ranef.pm, sig.pm, Yhat, BW.samples, BZ.samples, Sigma_samples,
             sigma2_z_samples, sigma2_w_samples, P.mat, v, Psi, Az, Bz, Aw, Bw, Wi)
  names(ret) = c("beta.pm", "beta.LB", "beta.UB", "ranef.pm", "sig.pm", "Yhat", "BW_samples", "BZ_samples",
                 "Sigma_samples", "sigma_Z_samples", "sigma_W_samples", "Penalty_Matrix", "inv.wish.df",
                 "inv.wish.scale", "SigmaZ.a", "SigmaZ.b", "SigmaW.a", "SigmaW.b", "fixed_effect_design")
  
  ret
}


test2 <- goldsmiths.one.d.gibbs(Y = test1$obs, fixed_design.mat = test1$x_design,
                                random_design.mat = test1$z_design, Kt=5)

#function compute log prior of the functional model
univariate.log.prior <- function(Sigma.samples, Sigma.prior.df, Sigma.prior.scale, sigm2.samps, 
                                 sigmaZ.hyper1, sigmaZ.hyper2, sigmaw2.samps, sigmaW.hypers1, 
                                 sigmaW.hypers2, comp.penalty.mat, BW_s, BZ_s, Wi){
  #get parameters for everything
  S <- dim(Sigma.samples)[3]
  Kt <- nrow(BW_s)
  p <- ncol(BW_s)
  I <- ncol(BZ_s)
  print(dim(BZ_s))
  
  #init log prior values
  log_priors <- numeric(S)
  
  #extract the sth sample
  for (s in 1:S) {
    #sigma_w^2 portion
    lambda_bw_s <- sigmaw2.samps[s,]
    lp_sigma_w2 <- sum(log(mapply(MCMCpack::dinvgamma, lambda_bw_s, 
                                  shape = sigmaW.hypers1, scale = sigmaW.hypers2)))
    
    #sigma_z^2 portion
    lambda_bz_s <- sigm2.samps[s]
    lp_sigma_z2 <- log(MCMCpack::dinvgamma(lambda_bz_s, shape = sigmaZ.hyper1, scale = sigmaZ.hyper2))
    
    #the Sigma portion
    Sigma_s  <- Sigma.samples[,,s]
    lp_Sigma <- log(MCMCpack::diwish(Sigma_s, v = Sigma.prior.df, S = Sigma.prior.scale))
    
    #BW portion
    lp_BW <- 0
    for (k in 1:p) {
      bw_k <- BW_s[,k,s]  
      cov_k <- lambda_bw_s[k]*solve(comp.penalty.mat)
      lp_BW <- lp_BW + dmvnorm(bw_k, mean = rep(0, Kt), sigma = cov_k, log = TRUE)
    }
    
    #BZ portion
    lp_BZ <- 0
    for (i in 1:I) {
      mu_i <- BW_s[,,s] %*% t(Wi[i, ])  #mean vector of length Kt
      cov_i <- lambda_bz_s * solve(comp.penalty.mat)
      bz_i <- BZ_s[,i,s]
      lp_BZ <- lp_BZ + dmvnorm(bz_i, mean = as.vector(mu_i), sigma = cov_i, log = TRUE)
    }
    
    #add everything together
    log_priors[s] <- lp_sigma_w2 + lp_sigma_z2 + lp_Sigma + lp_BW + lp_BZ
  }

  return(log_priors)
}

#testin the log prior function
test3 <- univariate.log.prior(Sigma.samples = test2$Sigma_samples, Sigma.prior.df = test2$inv.wish.df,
                              Sigma.prior.scale = test2$inv.wish.scale, sigm2.samps = test2$sigma_Z_samples, 
                              sigmaZ.hyper1 = test2$SigmaZ.a, 
                              sigmaZ.hyper2 = test2$SigmaZ.b, sigmaw2.samps = test2$sigma_W_samples,
                              sigmaW.hypers1 = test2$SigmaW.a, sigmaW.hypers2 = test2$SigmaW.b,
                              comp.penalty.mat = test2$Penalty_Matrix, BW_s = test2$BW_samples,
                              BZ_s = test2$BZ_samples, Wi = test2$fixed_effect_design)


plot(test3)

#function for the log likelihood of the model
univaraite.log.likelihood <- function(){
  
}

#Goldsmith's modified Gibbs sampler for 3 dimensions
goldsmiths.three.d.gibbs <- function(Y, fixed_design.mat, random_design.mat, 
                                Kt, N.iter = 1000, N.burn = 200, alpha = .1){
  ##This Gibbs sampler accepts a Y matrix composed of concatenated observations
  ##for the 3 dimensions, a fixed effect design matrix, a random effect design
  ##matrix, the number of knots to be fit for each dimension, as well as overrideable 
  ##defaults for the number of iterations and number of burn in iterations.  One can
  ##also change the alpha term.
  
  
  set.seed(1)
  
  ## fixed and random effect design matrices
  # W.des = model.matrix( fixef.form, data = data)
  # W.des <- W.des[rep(seq_len(nrow(W.des)), length(id)), ]
  # Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  # W.des = as.spam(W.des)
  # print(dim(W.des))
  # Z.des = as.spam(Z.des)
  ## fixed and random effect design matrices
  #W.des = model.matrix( fixef.form, data = data)
  #Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  W.des = as.spam(fixed_design.mat)
  Z.des = as.spam(random_design.mat)
  
  I = dim(Z.des)[2]
  D = dim(Y)[2]/3 #divide by three here
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  Gamma = kronecker(diag(1,3), Theta)
  
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  # not doing anything to the penalty matrix here
  #need to wait until we declare the variance terms first
  
  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))
  
  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  
  Wi = W.des[firstobs, ]
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  IIP = kronecker(kronecker(diag(1, I, I), P.mat), diag(1,3)) #kronecker this by a 3*3 matrix
  WIk = kronecker(Wi, diag(1, 3*Kt, 3*Kt)) #double the knots
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP
  
  # initial estimation and hyperparameter choice
  vec.bz = solve(kronecker(t(Z.des)%*% Z.des, t(Gamma) %*% Gamma)) %*% t(kronecker(Z.des, Gamma)) %*% Y.vec
  bz = matrix(vec.bz, nrow = 3*Kt, ncol = I) #double the knots
  
  w.temp = kronecker(t(Wi), diag(1, Kt, Kt))
  vec.bw = solve(tWIW) %*% tWI %*% vec.bz
  bw = matrix(vec.bw, nrow = 3*Kt, ncol = p) #need to double the knots
  
  Yhat = as.matrix(Z.des %*% t(bz) %*% t(Gamma))
  varhat = var(as.vector(Y - Yhat))
  
  Psi = diag(varhat*IJ, 3*D, 3*D)
  v = IJ
  inv.sig = solve(Psi/v)
  
  Az = I*Kt / 2  #set Az equal to Az1 equal to Az2
  Az1 = I*Kt / 2 
  Az2 = I*Kt / 2
  Bz = sum(diag((t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ])) %*% P.mat %*% t(t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ]))))
  Bz1 = sum(diag((t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ]))))
  Bz2 = sum(diag((t(bz[((2*Kt)+1):(3*Kt), ]) - Wi %*% t(bw[((2*Kt)+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[((2*Kt)+1):(3*Kt), ]) - Wi %*% t(bw[((2*Kt)+1):(3*Kt), ]))))
  
  Aw = Kt / 2  #set Aw equal to Aw1
  Aw1 = Kt / 2
  Aw2 = Kt / 2
  Bw = sapply(1:p, function(u) max(1, sum(diag( t(bw[1:Kt,u]) %*% P.mat %*% (bw[1:Kt,u])))))
  Bw1 = sapply(1:p, function(u) max(1, sum(diag( t(bw[(Kt+1):(2*Kt),u]) %*% P.mat %*% (bw[(Kt+1):(2*Kt),u])))))
  Bw2 = sapply(1:p, function(u) max(1, sum(diag( t(bw[((2*Kt)+1):(3*Kt),u]) %*% P.mat %*% (bw[((2*Kt)+1):(3*Kt),u])))))
  
  ## matrices to store within-iteration estimates
  BW = array(NA, c(3*Kt, p, N.iter)) #change t0 2*Kt
  BW[,,1] = bw
  BZ = array(NA, c(3*Kt, I, N.iter)) #change to 2*Kt
  BZ[,,1] = bz
  INV.SIG = array(NA, c(3*D, 3*D, N.iter)) #change to 2D
  INV.SIG[,,1] = inv.sig
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW[1,] = lambda.bw = Aw/Bw
  LAMBDA.BW.1 = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW.1[1,] = lambda.bw.1 = Aw1/Bw1
  LAMBDA.BW.2 = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW.2[1,] = lambda.bw.2 = Aw2/Bw2
  LAMBDA.BZ = rep(NA, N.iter)
  LAMBDA.BZ[1] = lambda.ranef = Az/Bz
  LAMBDA.BZ.1 = rep(NA, N.iter)
  LAMBDA.BZ.1[1] = lambda.ranef.1 = Az1/Bz1
  LAMBDA.BZ.2 = rep(NA, N.iter)
  LAMBDA.BZ.2[1] = lambda.ranef.2 = Az2/Bz2
  
  y.post = array(NA, dim = c(IJ, 3*D, (N.iter - N.burn))) #change dimensions to 2*D
  
  cat("Beginning Sampler \n")
  pb <- txtProgressBar(min = 0, max = N.iter, initial = 0, style = 3)
  for(i in 1:N.iter){
    setTxtProgressBar(pb,i)
    # if(i %% 100 == 0){ #adding a print statement to tell us where we are
    #   print(i)
    # }
    
    #stick these outside of the subj loop?
    combined <- diag(c(lambda.ranef, lambda.ranef.1, lambda.ranef.2))
    two.P = kronecker(P.mat, combined)  #constructing P \otimes two variances
    
    #sigma_w_k.pre.kron <- c(LAMBDA.BW[i, ], LAMBDA.BW.1[i, ])
    sigma_w_k.pre.kron <- c(lambda.bw, lambda.bw.1, lambda.bw.2)
    #print(sigma_w_k.pre.kron)
    sigma_w_k.pre.kron.1 <- diag(sigma_w_k.pre.kron, nrow = length(sigma_w_k.pre.kron))
    sigma_w_k.post.kron <- kronecker(sigma_w_k.pre.kron.1, P.mat)
    #print(dim(sigma_w_k.post.kron))
    
    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################
    for(subj in 1:length(unique(SUBJ))){
      
      t.designmat.Z = t(kronecker(rep(1, Ji[subj]), Gamma))   #change to Gamma
      
      
      #print(dim(two.P))
      sigma = solve(t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% t(t.designmat.Z) +
                      two.P)
      mu = sigma %*% (t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),]))) +
                        (two.P) %*% bw %*% t(Wi[subj,]))
      
      bz[,subj] = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = 3*Kt, ncol = 1)
    }
    ranef.cur = Z.des %*% t(bz) %*% t(Gamma)
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    #error in evaluating the argument 'a' in selecting a method for function 
    #'solve': NAs in argument 5 and 'NAOK = FALSE' (dotCall64)
    
    #maybe here?
    #IIP.1 <- kronecker(kronecker(diag(1, I, I), P.mat), combined)
    IIP.1 <- kronecker(diag(1, I, I), two.P)
    tWIW.1 <- t(WIk) %*% IIP.1 %*% WIk
    
    tWI.1 = t(WIk) %*% IIP.1
    
    sigma = solve(tWIW.1 + sigma_w_k.post.kron)
    #print(determinant(sigma))
    mu = sigma %*% (tWI.1 %*% as.vector(bz))
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = 3*Kt, ncol = p)
    # 
    beta.cur = t(bw) %*% t(Gamma)
    
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    resid.cur = Y - ranef.cur
    inv.sig = solve(riwish(v + IJ, Psi + t(resid.cur) %*% resid.cur))
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      a.post.1 = Aw1 + Kt/2
      a.post.2 = Aw2 + Kt/2
      b.post = Bw[term] + 1/2 * bw[(1:Kt),term] %*% P.mat %*% bw[(1:Kt),term]
      b.post.1 = Bw1[term] + 1/2 * bw[(Kt+1):(2*Kt),term] %*% P.mat %*% bw[(Kt+1):(2*Kt),term]
      b.post.2 = Bw1[term] + 1/2 * bw[(2*Kt+1):(3*Kt),term] %*% P.mat %*% bw[(2*Kt+1):(3*Kt),term]
      lambda.bw[term] = rgamma(1, a.post, b.post)
      lambda.bw.1[term] = rgamma(1, a.post.1, b.post.1)
      lambda.bw.2[term] = rgamma(1, a.post.2, b.post.2)
    }
    
    ## lambda for random effects
    a.post = Az + I*Kt/2
    a.post.1 = Az1 + I*Kt/2
    a.post.2 = Az2 + I*Kt/2
    b.post = Bz + .5 * sum(sapply(1:I, function(u) (t(bz[(1:Kt),u]) - Wi[u,] %*% t(bw[(1:Kt), ])) %*% P.mat %*% t(t(bz[(1:Kt),u]) - Wi[u,] %*% t(bw[(1:Kt), ])) ))
    b.post.1 = Bz1 + .5 * sum(sapply(1:I, function(u) (t(bz[(Kt+1):(2*Kt),u]) - Wi[u,] %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt),u]) - Wi[u,] %*% t(bw[(Kt+1):(2*Kt), ])) ))
    b.post.2 = Bz1 + .5 * sum(sapply(1:I, function(u) (t(bz[(2*Kt+1):(3*Kt),u]) - Wi[u,] %*% t(bw[(2*Kt+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[(2*Kt+1):(3*Kt),u]) - Wi[u,] %*% t(bw[(2*Kt+1):(3*Kt), ])) ))
    lambda.ranef = rgamma(1, a.post, b.post)
    lambda.ranef.1 = rgamma(1, a.post.1, b.post.1)
    lambda.ranef.2 = rgamma(1, a.post.2, b.post.2)
    
    ###############################################################
    ## save this iteration's parameters
    ###############################################################
    
    BW[,,i] = as.matrix(bw)
    BZ[,,i] = as.matrix(bz)
    
    INV.SIG[,,i] = inv.sig
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.BW.1[i,] = lambda.bw.1
    LAMBDA.BW.2[i,] = lambda.bw.2
    LAMBDA.BZ[i] = lambda.ranef
    LAMBDA.BZ.1[i] = lambda.ranef.1
    LAMBDA.BZ.2[i] = lambda.ranef.2
    
    if(i > N.burn){
      y.post[,,i - N.burn] = ranef.cur
    }
  }
  close(pb)
  
  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################
  
  #Save the relevant posterior summaries
  
  ## main effects
  beta.post = array(NA, dim = c(p, 3*D, (N.iter - N.burn)))
  for(n in 1:(N.iter - N.burn)){
    beta.post[,,n] = t(BW[,, n + N.burn]) %*% t(Gamma)
  }
  
  beta.pm = apply(beta.post, c(1,2), mean)
  beta.LB = apply(beta.post, c(1,2), quantile, c(.025))
  beta.UB = apply(beta.post, c(1,2), quantile, c(.975))
  
  ## random effects
  b.pm = matrix(NA, nrow = I, ncol = 3*D)
  for(i in 1:I){
    b.post = matrix(NA, nrow = (N.iter - N.burn), ncol = 3*D)
    for(n in 1:(N.iter - N.burn)){
      b.post[n,] = BZ[,i, n + N.burn] %*% t(Gamma)
    }
    b.pm[i,] = apply(b.post, 2, mean)
  }
  
  ## covariance matrix
  sig.pm = solve(apply(INV.SIG, c(1,2), mean))
  
  ## export fitted values
  ranef.pm = Z.des %*% b.pm
  Yhat = apply(y.post, c(1,2), mean)
  
  ####
  ###mcmc diagnostics
  #b_w matrix
  b.w.mcmc <- data.frame(matrix(BW, nrow=dim(BW)[3], byrow=TRUE))
  b.w.mcmc <- mcmc(tail(b.w.mcmc, n = (N.iter - N.burn)))
  print(effectiveSize(b.w.mcmc)[which.min(effectiveSize(b.w.mcmc))])
  
  #b_z matrix
  b.z.mcmc <- data.frame(matrix(BZ, nrow=dim(BZ)[3], byrow=TRUE))
  b.z.mcmc <- mcmc(tail(b.z.mcmc, n = (N.iter - N.burn)))
  print(effectiveSize(b.z.mcmc)[which.min(effectiveSize(b.z.mcmc))])
  
  #lambda bz terms
  lambda.bzs <- mcmc(tail(cbind(LAMBDA.BZ, LAMBDA.BZ.1, LAMBDA.BZ.2), n = (N.iter - N.burn)))
  print(effectiveSize(lambda.bzs))
  
  #lambda bw terms
  lambda.bws <- mcmc(tail(cbind(LAMBDA.BW, LAMBDA.BW.1, LAMBDA.BW.2), n = (N.iter - N.burn)))
  print(effectiveSize(lambda.bws))
  
  #returning what we want to return
  ret = list(beta.pm, beta.LB, beta.UB, ranef.pm, sig.pm, Yhat, BW, BZ)
  names(ret) = c("beta.pm", "beta.LB", "beta.UB", "ranef.pm", "sig.pm", "Yhat", "B_W", "B_Z")
  
  ret
}


#plot the posterior output of the Gibbs samplers
#getting the posterior mean and bounds
clean.gibbs.output.univariate <- function(posterior.ouptut){
  #posterior mean
  post.means <- posterior.ouptut$beta.pm %>% as.data.frame() %>%
    mutate(parameters = row_number()) %>% 
    pivot_longer(!parameters, names_to = "time", values_to = "value")
  post.means$time <- gsub('V', '', post.means$time)
  post.means$time <- as.numeric(post.means$time)
  
  #posterior lower bound
  post.lower <- posterior.ouptut$beta.LB %>% as.data.frame() %>%
    mutate(parameters = row_number()) %>% 
    pivot_longer(!parameters, names_to = "time", values_to = "LB")
  post.lower$time <- gsub('V', '', post.lower$time)
  post.lower$time <- as.numeric(post.lower$time)
  
  #posterior upper bound
  post.upper <- posterior.ouptut$beta.UB %>% as.data.frame() %>%
    mutate(parameters = row_number()) %>% 
    pivot_longer(!parameters, names_to = "time", values_to = "UB")
  post.upper$time <- gsub('V', '', post.upper$time)
  post.upper$time <- as.numeric(post.upper$time)
  
  output <- post.means
  output$LB <- post.lower$LB
  output$UB <- post.upper$UB
  
  return(output)
}

#actually plotting the regression coefficients with CIs
plot.univariate.gibbs <- function(cleaned.posterior){
  plot1 <- ggplot(cleaned.posterior, aes(x = time, y = value, group = parameters)) +
    geom_ribbon(aes(ymin = LB, ymax = UB), fill = "grey80", alpha = 0.5) +
    geom_line(aes(color = factor(parameters))) +
    labs(x = "Time", y = "Value", color = "Parameter", 
         title = "Posterior Mean of Functional Regression Coefficients with 95% Credible Intervals") +
    theme_minimal() 
  
  return(plot1)
}


plot1 <- clean.gibbs.output.univariate(posterior.ouptut = test2)
plot2 <- plot.univariate.gibbs(cleaned.posterior = plot1)
plot2
