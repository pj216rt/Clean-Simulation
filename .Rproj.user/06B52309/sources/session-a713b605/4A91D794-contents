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

