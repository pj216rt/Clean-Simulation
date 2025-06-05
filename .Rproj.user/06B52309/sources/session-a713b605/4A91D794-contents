library(tidyverse)
library(splines) #for working with b splines
library(mvtnorm) #multivariate normal

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
  X <- matrix(rnorm(n=(num.subj*n.params), sd = sqrt(var.x)), 
              nrow = num.subj, ncol = n.params) %>%
    as.data.frame() %>%
    mutate(subj = row_number())
  
  X$subj = as.factor(X$subj)
  X1 <- X[rep(seq_len(nrow(X)), each = num.visits), ]
  X.des <- model.matrix(~. -1 -subj, data = X1)
  
  #generate fixed effects
  fixef = as.matrix(X.des) %*% as.matrix(beta.vals)
  
  #random effects
  Z.des = model.matrix( ~ 0 + subj + (-1):subj, data = X1)
  subj.ranef <- matrix(rnorm(n=(num.subj*grid.size.3d), sd=sqrt(var.z)), 
                       nrow = num.subj, ncol = grid.size.3d)
  ranef <- Z.des %*% subj.ranef
  
  #level 1 residuals
  eps <- rnorm(n = tot.obs)
  
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
