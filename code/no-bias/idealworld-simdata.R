### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This function generates data, (co)variance matrices and mean structures for population and sample

simdata <- function(lambda1, lambda2, lambda3, alpha, nobs, n, q, k) {
  
  # Fixed values for all simulations
  phi <- c(1)                                                             # latent variable variances
  nu.1 <- c(0,0,0)  													                              # Intercepts for all indicators are 0 in group 1
  nu.2 <- c(0,0,0)												                                # Intercept 3 varies over condition
  
  # Matrices
  alpha.1 <- matrix(0, nrow = n, ncol = n)               					          # latent mean group 1 is 0
  alpha.2 <- matrix(alpha, nrow = n, ncol = n)               				        # latent mean group 2 varies over condition
  
  lambda <- matrix(c(lambda1,lambda2,lambda3), nrow = q, ncol = n)          # lambdas (factor loadings)
  errorvar <- 1 - lambda[1:q] * lambda[1:q]								                  # error variance 
  theta <- matrix(c(errorvar[1], 0, 0, 
                    0, errorvar[2], 0, 
                    0, 0, errorvar[3]), nrow = q, ncol = q)                 # matrix of error variances
  phi <- matrix(phi, nrow = n, ncol = n)									                  # matrix of latent variable variance
  
  # Simulating data
  eta.1 <- mvrnorm(nobs, mu = alpha.1, Sigma = phi) 							          # factor scores group 1 and 2
  eta.2 <- mvrnorm(nobs, mu = alpha.2, Sigma = phi)
  
  epsilon.1 <- mvrnorm(nobs, rep(0, q), Sigma = theta)						          # residual scores group 1 and 2
  epsilon.2 <- mvrnorm(nobs, rep(0, q), Sigma = theta) 
  
  df.1 <- as.data.frame(eta.1 %*% t(lambda) + epsilon.1) 					          # measurement equations group 1 and 2
  
  # add the intercepts (zero in this case so not strictly necessary)
  df.1[,1] <- df.1[,1] + nu.1[1]
  df.1[,2] <- df.1[,2] + nu.1[2]
  df.1[,3] <- df.1[,3] + nu.1[3]
  
  df.2 <- as.data.frame(eta.2 %*% t(lambda) + epsilon.2)
  
  df.2[,1] <- df.2[,1] + nu.2[1]
  df.2[,2] <- df.2[,2] + nu.2[2] 
  df.2[,3] <- df.2[,3] + nu.2[3]  
  
  ### If you write out the complete measurement equation like this:
  ### df.2 <- as.data.frame(nu.2 + eta.2 %*% t(lambda) + epsilon.2)
  ### the biased intercept (the third one) will not add to only the 
  ### third column, but subsequently to column 1, column 2, and column 3.
  ### as such, make dataframe in two steps, adding the intercept vector last
  
  df.1[,4] <- rep(1,nobs); df.2[,4] <- rep(2,nobs)							            # assigning group values per group
  df <- rbind(df.1,df.2)													                          # make one dataframe
  colnames(df) <- c("V1","V2","V3","G")										                  # assign column names
  
  # Population covariance and mean structure
  sigma.pop.1 <- lambda %*% phi %*% t(lambda) + theta 					          	# population covariance group 1 and 2
  sigma.pop.2 <- lambda %*% phi %*% t(lambda) + theta
  
  mu.pop.1 <- nu.1 + lambda %*% alpha.1									                    # population mean structure group 1 and 2
  mu.pop.2 <- nu.2 + lambda %*% alpha.2
  
  # Sample covariance and mean structure
  sigma.sample.1 <- cov(df.1)[1:q,1:q]								                  		# sample covariance group 1 and 2
  sigma.sample.2 <- cov(df.2)[1:q,1:q]
  
  mu.sample.1 <- as.matrix(colMeans(df.1)[1:q])					              			# sample mean structure group 1 and 2
  mu.sample.2 <- as.matrix(colMeans(df.2)[1:q])
  
  colnames <- c("V1","V2","V3")										                      		# column names for all covariances and meanstructures
  colnames(sigma.pop.1) <- colnames(sigma.pop.2) <- colnames(sigma.sample.1) <- colnames(sigma.sample.2) <- colnames
  
  data <- list(df,sigma.pop.1,sigma.pop.2,mu.pop.1,mu.pop.2,sigma.sample.1,sigma.sample.2,mu.sample.1,mu.sample.2)        # all data in a list
  names(data) <- c("df","sigma.pop.1","sigma.pop.2","mu.pop.1","mu.pop.2","sigma.sample.1","sigma.sample.2","mu.sample.1","mu.sample.2")
  
  return(data)
}
