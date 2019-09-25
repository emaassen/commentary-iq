
#iters <- 20

# Function to save warnings and errors ------------------------------------
myTryCatch <- function(expr) {
  warn <- err <- "none"                                                     # if no errors occur, let warnings and errors be "none"
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e                                                             # if an error occurs, save it in err
    }), warning=function(w) {
      warn <<- w                                                            # if a warning occurs, save it in warn
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)                                
}

# Simulate data function --------------------------------------------------
simdata <- function(lambda1, lambda2, lambda3, nus, alpha, nobs, n, q, k) {
  
  # Fixed values for all simulations
  phi <- c(1,1)                                                             # latent variable variances
  nu.1 <- c(0,0,0)  													                              # Intercepts for all indicators are 0 in group 1
  nu.2 <- nus           										                                # Intercept 3 varies over condition
  
  # Matrices
  alpha.1 <- matrix(0, nrow = n, ncol = n)               					          # latent mean group 1 is 0
  alpha.2 <- matrix(alpha, nrow = n, ncol = n)               				        # latent mean group 2 varies over condition
  
  lambda <- matrix(c(lambda1,lambda2,lambda3), nrow = q, ncol = n)          # lambdas (factor loadings)
  errorvar <- 1-lambda[1:q]*lambda[1:q]										                  # error variance 
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
  df.2[,3] <- df.2[,3] + nu.2[3]  # add the biased intercept
  
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

# CFI function ------------------------------------------------------------
cfi <- function(data, nobs, pow.lv) {
  
  # Sample variance covariance matrix and sample mean structure
  sigma.sample <- list(data$"sigma.sample.1",data$"sigma.sample.2")
  mu.sample <- list(data$"mu.sample.1",data$"mu.sample.2")
  
  # Independence model
  indepmodel <- '!no covariances between manifest variables
  V1 ~~ 0*V2
  V1 ~~ 0*V3
  V2 ~~ 0*V3
  
  !no covariances between latent variables
  F1 ~~ 0*F2
  F1 ~~ 0*F3
  F2 ~~ 0*F3
  
  !loadings fixed to 1
  F1 =~ 1*V1
  F2 =~ 1*V2
  F3 =~ 1*V3
  
  !factor variances fixed to 1
  F1 ~~ 1*F1
  F2 ~~ 1*F2
  F3 ~~ 1*F3
  
  !estimate residual variances
  V1 ~~ NA*V1
  V2 ~~ NA*V2
  V3 ~~ NA*V3
  
  !fix latent means to 0
  F1 ~ c(0,0)*1
  F2 ~ c(0,0)*1
  F3 ~ c(0,0)*1'
  
  # Fit independence model on sample variance covariance matrix and sample mean structure
  mod.indep <- myTryCatch(cfa(indepmodel, sample.cov = sigma.sample, sample.mean = mu.sample, sample.nobs = c(nobs, nobs), group.equal = c("intercepts", "residuals")))
  
  # CFI estimate of difference independence model with model 1 (= model with estimated mean difference between group 1 (=0) and 2 (=alpha))
  cfi <- tryCatch({															
    
    chi.in <- fitmeasures(mod.indep$`value`)["chisq"]					            	# chisquare independence model fit
    df.in <- fitmeasures(mod.indep$`value`)["df"]						              	# df independence model fit
    chi.full <- fitmeasures(unlist(pow.lv)$mod.1.value)["chisq"] 	      		# chisquare model 1 fit
    df.full <- fitmeasures(unlist(pow.lv)$mod.1.value)["df"]		        		# df model 1 fit
    
    cfi <- as.numeric(1 - (chi.full-df.full) / (chi.in-df.in))		        		# CFI formula
    
    cfi <- sapply(cfi, function(x) {x[x < 0] <- 0; x}) 					          	# if cfi is < 0, change to 0
    cfi <- sapply(cfi, function(x) {x[x > 1] <- 1; x}) 				          		# if cfi is > 1, change to 1
    
    cfi <- list(cfi)						                  	# print CFI and warning and error messages
    
  }, error = function(err) {
    
    cfi <- list(NA)					                  		# if error occurs because of no convergence so no fitmeasures can be estimated, let CFI = NA
    
  })
  
  return(cfi)
  
}

# Power functions ---------------------------------------------------------
power.cfa.th <- function(model2, data, nobs) {
  
  # To test this function step by step:
  #data <- datalist
  
  ### Theoretical power
  
  # Population variance covariance structure and population mean vector
  sigma.pop <- list(data$"sigma.pop.1",data$"sigma.pop.2")
  mu.pop <- list(data$"mu.pop.1",data$"mu.pop.2")
  
  mod.th.2 <- tryCatch({
    
    # Theoretical power
    
    # Fit model2, assuming no intercept difference, on population sigma and population mu
    mod.th.2 <- cfa(model2, sample.cov = sigma.pop, sample.mean = mu.pop, sample.nobs = c(nobs, nobs), group.equal = c("loadings", "intercepts"))
  
    # Chisquare estimate is NCP, calculate theoretical power
    ncp <- fitmeasures(mod.th.2)["chisq"]
    power.cfa.th <- 1 - pchisq(qchisq(.95, df = 1), df = 1, ncp = ncp)
    
    mod.th.2 <- list(ncp,power.cfa.th)
    
  }, error = function(err) {
    
    mod.th.2 <- list(NA,NA)
    
  })
  
  return(mod.th.2)
}

    
power.cfa.cal <- function(model1, model2, data, nobs) {
  
  # To test this function step by step:
  #data <- datalist
  
  # CALCULATED POWER
  # Sample variance covariance structure and sample mean vector
  sigma.sample <- list(data$"sigma.sample.1",data$"sigma.sample.2")
  mu.sample <- list(data$"mu.sample.1",data$"mu.sample.2")
    
  # Fit model1 (estimating intercept) and model2 (restricting intercept)
  mod.cal.1 <- myTryCatch(cfa(model1, sample.cov = sigma.sample, sample.mean = mu.sample, sample.nobs = c(nobs, nobs), group.equal = c("loadings")))
  mod.cal.2 <- myTryCatch(cfa(model2, sample.cov = sigma.sample, sample.mean = mu.sample, sample.nobs = c(nobs, nobs), group.equal = c("loadings")))

  # Save warning and error messages of both models
  mod.1.warn <- as.character(mod.cal.1$warning[1])
  mod.2.warn <- as.character(mod.cal.2$warning[1])
  mod.1.err <- as.character(mod.cal.1$error[1])
  mod.2.err <- as.character(mod.cal.2$error[1])
  
  # Save no of iterations of both models
  lv.iter.m1 <- lavTech(mod.cal.1$value, what = "iterations") 
  lv.iter.m2 <- lavTech(mod.cal.2$value, what = "iterations")  
  
  #parameterEstimates(mod.cal.1$value)
  #parameterEstimates(mod.cal.2$value)
  
  fitmeasures.m1 <- tryCatch ({
    
    # Save chisquare values and RMSEA values from model 1
    lv.chi.1 <- fitMeasures(mod.cal.1[[1]])["chisq"]    	          				# chi square model 1
    lv.rmsea.1 <- fitMeasures(mod.cal.1[[1]])["rmsea"]  		          			# RMSEA model 1
    
    fitmeasures.m1 <- list(lv.chi.1, lv.rmsea.1)
    
    }, error = function(err) { 
    
    lv.chi.1 <- NA           												                        # if non-convergence, estimates = NA	
    lv.rmsea.1 <- NA
    
    fitmeasures.m1 <- list(lv.chi.1, lv.rmsea.1)
    
    })
  
  fitmeasures.m2 <- tryCatch ({
    
    # Save chisquare values and RMSEA values from model 2
    lv.chi.2 <- fitMeasures(mod.cal.2[[1]])["chisq"]    	          				# chi square model 1
    lv.rmsea.2 <- fitMeasures(mod.cal.2[[1]])["rmsea"]  		          			# RMSEA model 1
    
    fitmeasures.m2 <- list(lv.chi.2,lv.rmsea.2)
    
  }, error = function(err) { 
    
    lv.chi.2 <- NA           												                        # if non-convergence, estimates = NA	
    lv.rmsea.2 <- NA
    
    fitmeasures.m2 <- list(lv.chi.2,lv.rmsea.2)
    
  })
  
  # Tally significant differences (p < .05) found in likelihood ratio test: pwr to detect intercept difference
  power.cfa.cal <- tryCatch ({
    power.cfa.cal <- lavTestLRT(mod.cal.1[[1]],mod.cal.2[[1]])$`Pr(>Chisq)`[2] < .05
  }, error = function(err) { power.cfa.cal <- 0 })						            	# if error occurs due to non convergence of model 1 or 2, power = NA
  
  ### Save results in list
  power.cfa.list <- list(power.cfa.cal, mod.cal.1, mod.cal.2, as.numeric(fitmeasures.m1[[1]]), as.numeric(fitmeasures.m2[[1]]), as.numeric(fitmeasures.m1[[2]]), as.numeric(fitmeasures.m2[[2]]), mod.1.warn, mod.2.warn, mod.1.err, mod.2.err,lv.iter.m1,lv.iter.m2)
  names(power.cfa.list) <- c("pow.cal","mod.1","mod.2","lv.chi.1","lv.chi.2","lv.rmsea.1","lv.rmsea.2","mod.1.warn","mod.2.warn","mod.1.err","mod.2.err","mod.1.iter","mod.2.iter")
  
  return(power.cfa.list)
  
}

power.sanova <- function(data, df, nobs, n, q, k) {
  
  ### Theoretical power
  
  # Calculating theoretical effect size
  lq.g1 <- nobs - 1
  lq.g2 <- nobs - 1
  var.p.1 <- sum(data$"sigma.pop.1") 									                  	# variance sumscore group 1 and 2
  var.p.2 <- sum(data$"sigma.pop.2")
  pooled.sd <- lq.g1 * var.p.1 + lq.g2 * var.p.2 				            			# pooled SD calclulation
  pooled.sd <- sqrt(pooled.sd/(lq.g1 + lq.g2)) 
  
  sd.pop <- matrix(c(pooled.sd),nrow=1, ncol=1) 					            		# population standard deviation matrix
  mu.diff.pop <- sum(data$"mu.pop.2")-sum(data$"mu.pop.1")                # population mean difference sumscore
  smd.th <- solve(sd.pop) %*% mu.diff.pop       					            		# standardized mean difference calculation (Xiaofeng Steven Liu formula 6.6)
  
  # Theoretical power
  f <- smd.th * 0.5  													                          	# Cohen p.276: 2 groups equal N means f = cohen's d * 0.5
  df1 <- k - 1
  df2 <- (nobs - 1) * k
  ncp <- f^2 * (df1 + df2 + 1)  								                    			# Cohen p.414
  
  power.sanova.th <- 1 - pf(qf(.95, df1, df2), df1, df2, ncp) 				
  
  ### Calculated power
  
  # Calculated effect size
  lq.g1 <- nobs - 1
  lq.g2 <- nobs - 1
  var.s.1 <- sum(data$"sigma.sample.1") 						                			# variance sumscore group 1 and 2
  var.s.2 <- sum(data$"sigma.sample.2")
  pooled.sd <- lq.g1 * var.s.1 + lq.g2 * var.s.2 			            				# pooled SD calculation
  pooled.sd <- sqrt(pooled.sd/(lq.g1 + lq.g2)) 
  
  sd.sample <- matrix(c(pooled.sd),nrow=n, ncol=n)				          			# sample standard deviation matrix
  mu.diff.sample <- sum(data$"mu.sample.2") - sum(data$"mu.sample.1") 		# sample mean difference sumscore 
  smd.cal <- solve(sd.sample) %*% mu.diff.sample  						          	# standardized mean difference calculation (Xiaofeng Steven Liu formula 6.6)
  
  # Calculated coverage
  smd.se.cal <- sqrt(((nobs + nobs) / (nobs * nobs) + (smd.cal^2 / (2*(nobs + nobs - 2)))) * ((nobs + nobs) / (nobs + nobs - 2))) # standard error standardized mean difference
  cill.sanova <- as.numeric(smd.cal - (1.96 * smd.se.cal))				      	# lower limit CI sample
  ciul.sanova <- as.numeric(smd.cal + (1.96 * smd.se.cal))				      	# upper limit CI sample
  cov.sanova <- cill.sanova < smd.th & smd.th < ciul.sanova				      	# does population SMD fall within sample CI?
  
  # Calculated power
  df[,"sum"] <- apply(df[,1:q],1,function(x) {sum(x)})				        		# sumscore per respondent
  
  anovatable <- anova(lm(df[,"sum"] ~ G, data=df))							          # anova of sumscore
  power.sanova.cal <- anovatable$`Pr(>F)`[1] < .05							          # tally how often significance difference between groups on sumscore
  
  ### Save results in list
  power.sanova.list <- list(smd.th,power.sanova.th,smd.cal,smd.se.cal,cill.sanova,ciul.sanova,cov.sanova,power.sanova.cal)
  names(power.sanova.list) <- c("smd.th","pow.th","smd.cal","smd.se.cal","cill","ciul","cov","pow.cal")
  
  return(power.sanova.list)
} 

power.anova <- function(data, df, nobs, q, k) {
  
  ### Theoretical power
  
  # Theoretical effect size
  sd.pop <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=q, ncol=q) 				        	# standard deviations indicator variables
  mu.diff.pop <- (data$"mu.pop.2")-(data$"mu.pop.1")                      # mean difference indicator variables 
  smd.th <- solve(sd.pop) %*% mu.diff.pop                      		    		# standardized mean difference calculation (Xiaofeng Steven Liu formula 6.6)
  
  # Theoretical power
  f <- smd.th * 0.5  													                          	# Cohen p.276: 2 groups equal N means f = cohen's d * 0.5
  df1 <- k - 1
  df2 <- (nobs - 1) * k
  ncp <- f^2 * (df1 + df2 + 1)  										                    	# Cohen p.414
  
  power.anova.th <- 1 - pf(qf(.95, df1, df2), df1, df2, ncp) 
  
  ### Calculated power
  power.anova.cal <- c()												                        	# empty vector
  
  # Calculated effect size
  lq.g1 <- nobs - 1
  lq.g2 <- nobs - 1
  var.1 <- data$"sigma.sample.1" 										                    	# variance group 1 and 2
  var.2 <- data$"sigma.sample.2"
  pooled.sd <- lq.g1 * var.1 + lq.g2 * var.2 				              				# pooled variance
  pooled.sd <- pooled.sd/(lq.g1 + lq.g2)									
  sd.pool <- sqrt(diag(pooled.sd))									                  		# pooled SD calculation
  
  sd.sample <- matrix(c(sd.pool[1],0,0,0,sd.pool[2],0,0,0,sd.pool[3]),nrow=q, ncol=q)	# sample standard deviation matrix
  mu.diff.sample <- (data$"mu.sample.2")-(data$"mu.sample.1") 					  # sample mean difference
  
  smd.cal <- solve(sd.sample) %*% mu.diff.sample  					          		# standardized mean difference calculation (Xiaofeng Steven Liu formula 6.6)
  
  # Calculated coverage
  smd.se.cal <- sqrt(((nobs + nobs) / (nobs * nobs) + (smd.cal^2 / (2 * (nobs + nobs - 2)))) * ((nobs + nobs) / (nobs + nobs - 2)))																	   # standard error standardized mean difference
  cill.anova <- as.numeric(smd.cal - (1.96 * smd.se.cal))					        # lower limit CI sample
  ciul.anova <- as.numeric(smd.cal + (1.96 * smd.se.cal))					        # lower limit CI sample
  cov.anova <- cill.anova < smd.th & smd.th < ciul.anova 					        # does population SMD fall within sample CI?
  
  for (q in 1:q) {
    
    # Calculated power
    anovatable <- anova(lm(df[,q] ~ G, data=df))							            # anova of indicators
    power.anova.cal <- c(power.anova.cal,anovatable$`Pr(>F)`[1] < .05)		# tally how often significance difference between groups on indicators
    
    # Calculating Cohen's d with pooled standard deviation from raw data:
    
    # q.g1 <- df[df$G == '1',][q]
    # q.g2 <- df[df$G == '2',][q]
    # lq.g1 <- nrow(q.g1) - 1
    # lq.g2 <- nrow(q.g2) - 1
    # pooled.sd <- lq.g1 * var(q.g1) + lq.g2 * var(q.g2)
    # pooled.sd <- sqrt(pooled.sd/(lq.g1 + lq.g2))
    # meandif  <- sapply(q.g2,mean) - sapply(q.g1,mean)
    # smd.cal.q <- meandif/pooled.sd
    # smd.cal  <- c(smd.cal,smd.cal.q)
    
  }
  
  ### Save results in list
  power.anova.list <- list(smd.th,power.anova.th,smd.cal,smd.se.cal,cill.anova,ciul.anova,cov.anova,power.anova.cal) 
  names(power.anova.list) <- c("smd.th","pow.th","smd.cal","smd.se.cal","cill","ciul","cov","pow.cal")
  
  return(power.anova.list)
} 

power.manova <- function(data, df, nobs, q, k) {
  
  ### Theoretical power
  grandmean <- (data$"mu.pop.1" + data$"mu.pop.2") / k   				        	# grand mean for each variable over all groups
  
  # between group SSCP matrix
  sscp.b <- nobs * (((data$"mu.pop.1" - grandmean) %*% t(data$"mu.pop.1" - grandmean)) + ((data$"mu.pop.2" - grandmean) %*% t(data$"mu.pop.2" - grandmean)))
  
  # # calculating each value of the sscp.b matrix by hand
  # # grand mean
  # grandmean <- (data$"mu.pop.1" + data$"mu.pop.2") / k   				       	# grand mean for each variable over all groups
  #
  # # sum of squares
  # ss.x1 <- (nobs * (data$"mu.pop.1"[1] - grandmean[1])^2) + (nobs * (data$"mu.pop.2"[1] - grandmean[1])^2) # sum of squares x1, x2, x3
  # ss.x2 <- (nobs * (data$"mu.pop.1"[2] - grandmean[2])^2) + (nobs * (data$"mu.pop.2"[2] - grandmean[2])^2)
  # ss.x3 <- (nobs * (data$"mu.pop.1"[3] - grandmean[3])^2) + (nobs * (data$"mu.pop.2"[3] - grandmean[3])^2)
  #
  # # cross products
  # cp.x1x2 <- ((data$"mu.pop.1"[1] - grandmean[1]) *  (data$"mu.pop.1"[2] - grandmean[2]) * nobs) + ((data$"mu.pop.2"[1] - grandmean[1]) *  (data$"mu.pop.2"[2] - grandmean[2]) * nobs)
  # cp.x2x3 <- ((data$"mu.pop.1"[2] - grandmean[2]) *  (data$"mu.pop.1"[3] - grandmean[3]) * nobs) + ((data$"mu.pop.2"[2] - grandmean[2]) *  (data$"mu.pop.2"[3] - grandmean[3]) * nobs)
  # cp.x1x3 <- ((data$"mu.pop.1"[1] - grandmean[1]) *  (data$"mu.pop.1"[3] - grandmean[3]) * nobs) + ((data$"mu.pop.2"[1] - grandmean[1]) *  (data$"mu.pop.2"[3] - grandmean[3]) * nobs)
  #
  # # between groups SSCP matrix
  # sscp.b <- matrix(c(ss.x1, cp.x1x2, cp.x1x3, cp.x1x2, ss.x2, cp.x2x3, cp.x1x3, cp.x2x3, ss.x3), ncol=q, nrow=q)
  
  # the covariance matrix is the created from the SSCP matrix by dividing each element in the matrix by N - 1.
  # so from covariance matrix to SSCP is multiplying each element in the matrix by N - 1
  # since the population covariances are the same for group 1 and 2, only the one for group 1 will be used
  # do not do this when the population covariance matrices are different!!!
  
  sscp.t <- data$"sigma.pop.1" * ((nobs*k) - 1)		  					          	# total SSCP matrix
  sscp.e <- sscp.t - sscp.b                           					        	# within group SSCP matrix
  
  wilks.th <- det(sscp.e) / det(sscp.t)               					        	# population wilks lambda   
  
  ky <- q                                             					        	# number of variables in the Y set;     Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 490
  kx <- k - 1                                         					        	# number of groups - 1;                 Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 490
  s <- sqrt((ky^2 * kx^2 - 4) / (ky^2 + kx^2 - 5))    					        	# Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471 formula 10.1.9
  
  f2.th <- (wilks.th ^ (-1 / s)) - 1                  					        	# effect size;                          Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 473 formula 10.2.2
  
  # another way to get to effect size using sscp.e, sscp.b and s
  # install.packages("pracma");library(pracma)
  # f2 <- (nthroot(det(sscp.e+sscp.b),s) - nthroot(det(sscp.e),s)) / nthroot(det(sscp.e),s) # Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 473 formula 10.2.2
  
  kc <- 0                                           						          # number of variables in the set that is being partialled;             Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471
  ka <- 0                                           						          # number of variables in the set that is being partialled;             Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471
  kg <- 0                                           						          # number of variables in the set used for model 2 error reduction;     Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471
  # kc, ka, and kg are zero when the set does not exist for the type of association or error model in question
  
  m <- (nobs * 2) - max(kc, (ka + kg)) - ((ky + kx + 3) / 2)              # Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471 formula 10.1.8
  
  df1 <- ky * kx                                    						          # numerator df;                         Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471 formula 10.1.6
  df2 <- m * s + 1 - (df1 / 2)                      						          # denominator df;                       Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471 formula 10.1.7
  
  ncp <- f2.th * (df1 + df2 + 1)                      						        # noncentrality parameter noncentral F; Cohen 1988 - Statistical Power Analysis for the Social Sciences p. 471 formula 10.3.1
  
  pillai.th <- tr(sscp.b %*% (solve(sscp.b + sscp.e))) 						        # Save Pillai for possible recalculations (eg in Gpower)
  
  power.manova.th  <- 1 - pf(qf(.95, df1, df2), df1, df2, ncp)
  
  ### Calculated power
  
  # Calculated power
  pow.man <- summary(manova(cbind(V1, V2, V3) ~ G, data = df), test="Wilks")	# MANOVA group 1 and 2
  power.manova.cal <- pow.man$stats[1,6] < .05									          # tally how often significance difference between groups 
  
  # Calculated effect size
  wilks.cal <- pow.man$stats[1,2]												                  # Calculated Wilks lambda
  f2.cal <- (wilks.cal ^ (-1 / s)) - 1										              	# calculated effect size f2
  
  # Save Pillai for possible recalculations (eg in Gpower)
  pillai.cal <- summary(manova(cbind(V1, V2, V3) ~ G, data = df), test="Pillai")$stats[1,2]
  
  # Calculated coverage
  f2.se.cal <- as.numeric(sqrt(((nobs + nobs) / (nobs * nobs) + (f2.cal^2 / (2 * (nobs + nobs - 2)))) * ((nobs + nobs) / (nobs + nobs - 2)))) # standard error standardized mean difference
  cill.manova <- as.numeric(f2.cal - (1.96 * f2.se.cal))				      		# lower limit CI sample
  ciul.manova <- as.numeric(f2.cal + (1.96 * f2.se.cal))				      		# upper limit CI sample
  cov.manova <- cill.manova < f2.th & f2.th < ciul.manova				      		# does population effect size f2 fall within sample CI?
  
  ### Save results in list
  power.manova.list <- list(power.manova.th,f2.th,f2.cal,f2.se.cal,cill.manova,ciul.manova,cov.manova,power.manova.cal,wilks.th,wilks.cal,pillai.th,pillai.cal)
  names(power.manova.list) <- c("pow.th","f2.th","f2.cal","f2.se.cal","cill","ciul","cov","pow.cal","wilks.th","wilks.cal","pillai.th","pillai.cal")
  
  return(power.manova.list)
  
}


# Simulation function -----------------------------------------------------
do_sim = function(pos, cond, iters) {
  
  # to run this function step by step, run these two lines:
  #cond <- conditions
  #pos <- 1
  
  # Fixed values for all simulations
  n <- 1            									 	                                 	# no LVs       
  q <- 3            									 					                         	# no indicators    
  k <- 2            									 					                         	# no groups
  
  # Select condition by changing pos and selecting conditions from cond
  lambda1 <- cond$lambdas[pos][[1]][1]
  lambda2 <- cond$lambdas[pos][[1]][2]
  lambda3 <- cond$lambdas[pos][[1]][3]
  nu3 <- cond$nu3[pos]
  alpha <- cond$alpha[pos]
  nobs <- cond$nobs[pos]
  
  set.seed(as.numeric(rownames(conditions)[pos]))																	
  
  # Empty vectors & matrices ------------------------------------------------
  lv.est.m1 <- lv.se.m1 <- lv.cill <- lv.ciul <- lv.cov <- lv.rmsea.m1 <- lv.ncp.th <- lv.pow.th <- lv.pow.cal <- c()
  lv.l1.est.m1 <- lv.l1.se.m1 <- lv.l2.est.m1 <- lv.l2.se.m1 <- lv.l3.est.m1 <- lv.l3.se.m1 <- c()
  lv.r1.est.m1 <- lv.r1.se.m1 <- lv.r2.est.m1 <- lv.r2.se.m1 <- lv.r3.est.m1 <- lv.r3.se.m1 <- c()
  lv.r4.est.m1 <- lv.r4.se.m1 <- lv.r5.est.m1 <- lv.r5.se.m1 <- lv.r6.est.m1 <- lv.r6.se.m1 <- c()
  lv.rlv2.est.m1 <- lv.rlv2.se.m1 <- c() 
  lv.i1.est.m1 <- lv.i1.se.m1 <- lv.i2.est.m1 <- lv.i2.se.m1 <- lv.i3.est.m1 <- lv.i3.se.m1 <- lv.i4.est.m1 <- lv.i4.se.m1 <- c()
  lv.chi.m1 <- c()
  lv.est.m2 <- lv.se.m2 <- c()
  lv.l1.est.m2 <- lv.l1.se.m2 <- lv.l2.est.m2 <- lv.l2.se.m2 <- lv.l3.est.m2 <- lv.l3.se.m2 <- c()
  lv.r1.est.m2 <- lv.r1.se.m2 <- lv.r2.est.m2 <- lv.r2.se.m2 <- lv.r3.est.m2 <- lv.r3.se.m2 <- c()
  lv.r4.est.m2 <- lv.r4.se.m2 <- lv.r5.est.m2 <- lv.r5.se.m2 <- lv.r6.est.m2 <- lv.r6.se.m2 <- c()
  lv.rlv2.est.m2 <- lv.rlv2.se.m2 <- c()
  lv.i1.est.m2 <- lv.i1.se.m2 <- lv.i2.est.m2 <- lv.i2.se.m2 <- lv.i3.est.m2 <- lv.i3.se.m2 <- c()
  lv.chi.m2 <- lv.rmsea.m2 <- c()
  lv.iter.m1 <- lv.iter.m2 <- c()
  lv.cfi.est <- c()
  sanova.est <- sanova.se <- sanova.cill <- sanova.ciul <- sanova.cov <- sanova.est.th <- sanova.pow.th <- sanova.pow.cal <- c()
  anova.1.est <- anova.2.est <- anova.3.est <- c()
  anova.1.se <- anova.2.se <- anova.3.se <- c()
  anova.1.cill <- anova.2.cill <- anova.3.cill <- c()
  anova.1.ciul <- anova.2.ciul <- anova.3.ciul <- c()
  anova.1.cov <- anova.2.cov <- anova.3.cov <- c() 
  anova.1.est.th <- anova.2.est.th <- anova.3.est.th <- c()
  anova.1.pow.th <- anova.2.pow.th <- anova.3.pow.th <- c() 
  anova.1.pow.cal <- anova.2.pow.cal <- anova.3.pow.cal <- c()
  manova.est <- manova.se <- manova.cill <- manova.ciul <- manova.cov <- manova.est.th <- manova.pow.th <- manova.pow.cal <- c()
  manova.wilks.th <- manova.wilks.cal <- manova.pillai.th <- manova.pillai.cal <- c()  
  count.conv <- count.lv <- count.ov <- c()
  
  for (i in 1:iters) {
    
    datalist <- simdata(lambda1, lambda2, lambda3, nu3, alpha, nobs, n, q, k)		  	# list with all data, covar matrices and mean structures
    df <- datalist[["df"]]											                    		# raw data
    
    # MGCFA -------------------------------------------------------------
    
    # Make sure model 1 = model with LV mean difference estimated, model 2 = no LV mean difference
    
    # Model where one latent mean difference is estimated and intercept 3 is estimated 
    model1 <- ' !Estimate all loadings
    F =~ NA * V1 + NA * V2 + NA * V3
    
    !Fix one LV variance, estimate one
    F ~~ c(1,NA) * F
    
    !Fix one LV mean to zero, estimate one
    F ~ c(0,NA) * 1 
    
    !Restrict first two intercepts to be equal across groups, estimate 3rd 
    V1 ~ c(i1,i1) * 1
    V2 ~ c(i2,i2) * 1
    V3 ~ c(NA,NA) * 1'
    
    # Covariances = 18 
    # Parameters = 3 loadings, 1 LV variance, 1 LV mean, *6 residuals, 4 intercepts = 12 (df = 6)
    
    # Model where one latent mean difference is estimated and intercepts are restricted across groups
    model2 <- ' !Estimate all loadings
    F =~ NA * V1 + NA * V2 + NA * V3
    
    !Fix one LV variance, estimate one
    F ~~ c(1,NA) * F
    
    !Fix one LV mean to zero, estimate one
    F ~ c(0,NA) * 1
    
    !Restrict first two intercepts to be equal across groups, estimate 3rd 
    V1 ~ c(i1,i1) * 1
    V2 ~ c(i2,i2) * 1
    V3 ~ c(i3,i3) * 1'
    
    # Parameters = 3 loadings, 1 LV variance, 1 LV mean, *3 residuals, 3 intercepts = 11 (df = 3)
    
    # Calculate NCP, theoretical power, based on model2
    pow.lv.th <- power.cfa.th(model2 = model2, data = datalist, nobs = nobs)
    
    # Calculate NCP, theoretical power, calculated power, and give model1 and model2 fit, chi square and rmsea estimates and possible warning and error messages
    pow.lv.cal <- power.cfa.cal(model1 = model1, model2 = model2, data = datalist, nobs = nobs)
    
    # Calculate CFI and save possible warning and error messages
    cfi.res <- cfi(data = datalist, nobs = nobs, pow.lv = pow.lv.cal)
    
    # LV estimates ------------------------------------------------------------
    # Model 1 estimates (with mean difference and intercept estimated)
    lv.est.m1 <- c(lv.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"est"])  	
    lv.se.m1 <- c(lv.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"se"])     
    lv.cill <- c(lv.cill,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.lower"])    
    lv.ciul <- c(lv.ciul,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.upper"]) 
    lv.cov <- c(lv.cov,as.numeric(parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.lower"] < alpha & parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.upper"] > alpha))						 # mean estimate coverage
    lv.cfi.est <- c(lv.cfi.est,cfi.res[[1]])  											                    # cfi estimate
    lv.rmsea.m1 <- c(lv.rmsea.m1,unlist(pow.lv.cal)$"lv.rmsea.1") 							            # rmsea model 1
    lv.ncp.th <- c(lv.ncp.th,unlist(pow.lv.th)[[1]]) 									                  	# noncentrality parameter used to calculate theoretical power
    lv.pow.th <- c(lv.pow.th,unlist(pow.lv.th)[[2]]) 
    lv.pow.cal <- c(lv.pow.cal,as.numeric(unlist(pow.lv.cal)$"pow.cal")) 				          	# calculated power / type 1 error
    
    lv.i4.est.m1 <- c(lv.i4.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[19,"est"]) # this is the biased intercept       
    lv.i4.se.m1 <- c(lv.i4.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[19,"se"])
    
    # Additional estimates
    lv.l1.est.m1 <- c(lv.l1.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[1,"est"])  # factor loading estimates + standard errors
    lv.l1.se.m1 <- c(lv.l1.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[1,"se"])     
    lv.l2.est.m1 <- c(lv.l2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[2,"est"])  
    lv.l2.se.m1 <- c(lv.l2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[2,"se"])     
    lv.l3.est.m1 <- c(lv.l3.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[3,"est"])  
    lv.l3.se.m1 <- c(lv.l3.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[3,"se"])     
    lv.r1.est.m1 <- c(lv.r1.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[9,"est"])   # residual estimates + standard errors
    lv.r1.se.m1 <- c(lv.r1.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[9,"se"])     
    lv.r2.est.m1 <- c(lv.r2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[10,"est"])  
    lv.r2.se.m1 <- c(lv.r2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[10,"se"])     
    lv.r3.est.m1 <- c(lv.r3.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[11,"est"])  
    lv.r3.se.m1 <- c(lv.r3.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[11,"se"])     
    lv.r4.est.m1 <- c(lv.r4.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[20,"est"])
    lv.r4.se.m1 <- c(lv.r4.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[20,"se"])
    lv.r5.est.m1 <- c(lv.r5.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[21,"est"])
    lv.r5.se.m1 <-c(lv.r5.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[21,"se"])
    lv.r6.est.m1 <- c(lv.r6.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[22,"est"])
    lv.r6.se.m1 <- c(lv.r6.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[22,"se"])
    lv.rlv2.est.m1 <- c(lv.rlv2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[15,"est"])      # group 2 latent variable variance estimate + standard errors (variance LV group 1 = restricted to 1)
    lv.rlv2.se.m1 <-c(lv.rlv2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[15,"se"])       
    lv.i1.est.m1 <- c(lv.i1.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[6,"est"])           # intercept estimates + standard errors
    lv.i1.se.m1 <- c(lv.i1.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[6,"se"])            
    lv.i2.est.m1 <- c(lv.i2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[7,"est"])        
    lv.i2.se.m1 <- c(lv.i2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[7,"se"])
    lv.i3.est.m1 <- c(lv.i3.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[8,"est"])        
    lv.i3.se.m1 <- c(lv.i3.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[8,"se"])
    lv.chi.m1 <- c(lv.chi.m1,unlist(pow.lv.cal)$"lv.chi.1")									              	      # chisquare model 1
    lv.iter.m1 <- c(lv.iter.m1,unlist(pow.lv.cal)$"mod.1.iter")                                   # iterations model 1
    
    #  LV restricted estimates ------------------------------------------------
    # Model 2 estimates (no intercept differences)
    lv.est.m2 <- c(lv.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[16,"est"])	      # mean estimate group 2 
    lv.se.m2 <- c(lv.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[16,"se"])           # mean estimate se
    lv.l1.est.m2 <- c(lv.l1.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[12,"est"])   # factor loading estimates + standard errors
    lv.l1.se.m2 <- c(lv.l1.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[12,"se"])
    lv.l2.est.m2 <- c(lv.l2.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[13,"est"])
    lv.l2.se.m2 <- c(lv.l2.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[13,"se"])
    lv.l3.est.m2 <- c(lv.l3.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[14,"est"])
    lv.l3.se.m2 <- c(lv.l3.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[14,"se"])
    lv.r1.est.m2 <- c(lv.r1.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[4,"est"])   # residual estimates + standard errors
    lv.r1.se.m2 <- c(lv.r1.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[4,"se"])
    lv.r2.est.m2 <- c(lv.r2.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[5,"est"])
    lv.r2.se.m2 <- c(lv.r2.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[5,"se"])
    lv.r3.est.m2 <- c(lv.r3.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[6,"est"])
    lv.r3.se.m2 <- c(lv.r3.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[6,"se"])
    lv.r4.est.m2 <- c(lv.r4.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[20,"est"])
    lv.r4.se.m2 <- c(lv.r4.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[20,"se"])
    lv.r5.est.m2 <- c(lv.r5.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[21,"est"])
    lv.r5.se.m2 <-c(lv.r5.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[21,"se"])
    lv.r6.est.m2 <- c(lv.r6.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[22,"est"])
    lv.r6.se.m2 <- c(lv.r6.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[22,"se"])
    lv.rlv2.est.m2 <- c(lv.rlv2.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[15,"est"])  # group 2 latent variable variance estimate + standard errors (variance LV group 1 = restricted to 1)
    lv.rlv2.se.m2 <- c(lv.rlv2.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[15,"se"])
    lv.i1.est.m2 <- c(lv.i1.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[20,"est"])   # intercept estimates + standard errors
    lv.i1.se.m2 <- c(lv.i1.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[20,"se"])
    lv.i2.est.m2 <- c(lv.i2.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[21,"est"])
    lv.i2.se.m2 <- c(lv.i2.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[21,"se"])
    lv.i3.est.m2 <- c(lv.i3.est.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[22,"est"])
    lv.i3.se.m2 <- c(lv.i3.se.m2,parameterEstimates(unlist(pow.lv.cal)$mod.2.value)[22,"se"])
    lv.chi.m2 <- c(lv.chi.m2,unlist(pow.lv.cal)$"lv.chi.2")                                       # chisquare model 2
    lv.rmsea.m2 <- c(lv.rmsea.m2,unlist(pow.lv.cal)$"lv.rmsea.2")                                 # rmsea model 2
    lv.iter.m2 <- c(lv.iter.m2,unlist(pow.lv.cal)$"mod.2.iter")                                   # iterations model 2
    
    # Nonconvergence and warning counter
    count.conv <- c(count.conv, as.numeric(lavTech(unlist(pow.lv.cal)$mod.1.value, what = "converged")))
    
    if (unlist(pow.lv.cal)$"mod.1.warn" == "lavaan WARNING: some estimated lv variances are negative") {
        count.lv <- c(count.lv,1)
        count.ov <- c(count.ov,0)
      } else {
        if (unlist(pow.lv.cal)$"mod.1.warn" == "lavaan WARNING: some estimated ov variances are negative") {
          count.lv <- c(count.lv,0)
          count.ov <- c(count.ov,1)
        } else {
          count.lv <- c(count.lv,0)
          count.ov <- c(count.ov,0)
        }
        }
    
    # ANOVA sum ----------------------------------------------------
    
    # Calculate theoretical effect size, theoretical power, actual effect size + se, coverage, calculated power
    pow.sanova <- power.sanova(data = datalist, df = df, nobs = nobs, n = n, q = q, k = k)
    
    # Save results in appropriate vectors
    sanova.est.th <- c(sanova.est.th,as.numeric(unlist(pow.sanova)["smd.th"]))                # sum anova population effect size
    sanova.est <- c(sanova.est,unlist(pow.sanova)[["smd.cal"]])                               # sum anova sample effect size
    sanova.se <- c(sanova.se,unlist(pow.sanova)[["smd.se.cal"]])                              # sum anova sample effect size standard error
    sanova.cill <- c(sanova.cill,unlist(pow.sanova)[["cill"]])                                # sum anova effect size CI lower limit
    sanova.ciul <- c(sanova.ciul,unlist(pow.sanova)[["ciul"]])                                # sum anova effect size CI upper limit
    sanova.cov <- c(sanova.cov,as.numeric(unlist(pow.sanova)[["cov"]]))                       # sum anova coverage
    sanova.pow.th <- c(sanova.pow.th,unlist(pow.sanova)[["pow.th"]])                          # sum anova theoretical power 
    sanova.pow.cal <- c(sanova.pow.cal,as.numeric(unlist(pow.sanova)[["pow.cal"]]))           # sum anova calculated power 
    
    # ANOVA separate ----------------------------------------------------
    
    # Calculate theoretical effect size, theoretical power, actual effect size + se, coverage, calculated power
    pow.anova <- power.anova(data = datalist, df = df, nobs = nobs, q = q, k = k)
    
    # Save results in appropriate vectors 
    anova.1.est.th <- c(anova.1.est.th,unlist(pow.anova)[["smd.th1"]])                        # indicator 1 population effect size
    anova.1.est <- c(anova.1.est,unlist(pow.anova)[["smd.cal1"]])                             # indicator 1 sample effect size
    anova.1.se <- c(anova.1.se,unlist(pow.anova)[["smd.se.cal1"]])                            # indicator 1 sample effect size standard error
    anova.1.cill <- c(anova.1.cill,unlist(pow.anova)[["cill1"]])                              # indicator 1 effect size CI lower limit
    anova.1.ciul <- c(anova.1.ciul,unlist(pow.anova)[["ciul1"]])                              # indicator 1 effect size CI upper limit
    anova.1.cov <- c(anova.1.cov,as.numeric(unlist(pow.anova)[["cov1"]]))                     # indicator 1 coverage
    anova.1.pow.th <- c(anova.1.pow.th,unlist(pow.anova)[["pow.th1"]])                        # indicator 1 theoretical power
    anova.1.pow.cal <- c(anova.1.pow.cal,as.numeric(unlist(pow.anova)[["pow.cal1"]]))         # indicator 1 calculated power
    
    anova.2.est.th <- c(anova.2.est.th,unlist(pow.anova)[["smd.th2"]])                        # indicator 2 population effect size
    anova.2.est <- c(anova.2.est,unlist(pow.anova)[["smd.cal2"]])                             # indicator 2 sample effect size
    anova.2.se <- c(anova.2.se,unlist(pow.anova)[["smd.se.cal2"]])                            # indicator 2 sample effect size standard error
    anova.2.cill <- c(anova.2.cill,unlist(pow.anova)[["cill2"]])                              # indicator 2 effect size CI lower limit
    anova.2.ciul <- c(anova.2.ciul,unlist(pow.anova)[["ciul2"]])                              # indicator 2 effect size CI upper limit
    anova.2.cov <- c(anova.2.cov,as.numeric(unlist(pow.anova)[["cov2"]]))                     # indicator 2 coverage
    anova.2.pow.th <- c(anova.2.pow.th,unlist(pow.anova)[["pow.th2"]])                        # indicator 2 theoretical power
    anova.2.pow.cal <- c(anova.2.pow.cal,as.numeric(unlist(pow.anova)[["pow.cal2"]]))         # indicator 2 calculated power
    
    anova.3.est.th <- c(anova.3.est.th,unlist(pow.anova)[["smd.th3"]])                        # indicator 3 population effect size
    anova.3.est <- c(anova.3.est,unlist(pow.anova)[["smd.cal3"]])                             # indicator 3 sample effect size
    anova.3.se <- c(anova.3.se,unlist(pow.anova)[["smd.se.cal3"]])                            # indicator 3 sample effect size standard error
    anova.3.cill <- c(anova.3.cill,unlist(pow.anova)[["cill3"]])                              # indicator 3 effect size CI lower limit
    anova.3.ciul <- c(anova.3.ciul,unlist(pow.anova)[["ciul3"]])                              # indicator 3 effect size CI upper limit
    anova.3.cov <- c(anova.3.cov,as.numeric(unlist(pow.anova)[["cov3"]]))                     # indicator 3 coverage
    anova.3.pow.th <- c(anova.3.pow.th,unlist(pow.anova)[["pow.th3"]])                        # indicator 3 theoretical power
    anova.3.pow.cal <- c(anova.3.pow.cal,as.numeric(unlist(pow.anova)[["pow.cal3"]]))         # indicator 3 calculated power
    
    # MANOVA ------------------------------------------------------------
    # Calculate theoretical effect size, theoretical power, actual effect size + se, coverage, calculated power
    pow.manova <- power.manova(data = datalist, df = df, nobs = nobs, q = q, k = k)
    
    # Save results in appropriate vectors
    manova.est.th <- c(manova.est.th,unlist(pow.manova)[["f2.th"]])                           # manova population effect size f2
    manova.est <- c(manova.est,unlist(pow.manova)[["f2.cal"]])                                # manova sample effect size f2
    manova.se <- c(manova.se,unlist(pow.manova)[["f2.se.cal"]])                               # manova sample effect size f2 standard error
    manova.cill <- c(manova.cill,unlist(pow.manova)[["cill"]])                                # manova sample effect size f2 CI lower limit
    manova.ciul <- c(manova.ciul,unlist(pow.manova)[["ciul"]])                                # manova sample effect size f2 CI upper limit
    manova.cov <- c(manova.cov,as.numeric(unlist(pow.manova)[["cov"]]))                       # manova coverage
    manova.pow.th <- c(manova.pow.th,unlist(pow.manova)[["pow.th"]])                          # manova theoretical power
    manova.pow.cal <- c(manova.pow.cal,as.numeric(unlist(pow.manova)[["pow.cal"]]))           # manova calculated power
    manova.wilks.th <- c(manova.wilks.th,unlist(pow.manova)[["wilks.th"]])                    # manova theoretical wilks lambda
    manova.wilks.cal <- c(manova.wilks.cal,unlist(pow.manova)[["wilks.cal"]])                 # manova calculated wilks lambda
    manova.pillai.th <- c(manova.pillai.th,unlist(pow.manova)[["pillai.th"]])                 # manova theoretcal pillais trace
    manova.pillai.cal <- c(manova.pillai.cal,unlist(pow.manova)[["pillai.cal"]])              # manova calculated pillais trace
    
    
  }
  
  # Save results ------------------------------------------------------------
  results.cond <- cbind.data.frame(rep(as.numeric(rownames(conditions)[pos]),iters), rep(lambda1,iters), rep(lambda2,iters), rep(lambda3,iters), rep(nu3,iters), rep(alpha,iters), rep(nobs,iters), 
                                   lv.est.m1, lv.se.m1, lv.cill, lv.ciul, lv.cov, lv.cfi.est, lv.rmsea.m1, lv.ncp.th, lv.pow.th, lv.pow.cal,
                                   anova.1.est, anova.1.se, anova.1.cill, anova.1.ciul, anova.1.cov, anova.1.est.th, anova.1.pow.th, anova.1.pow.cal,
                                   anova.2.est, anova.2.se, anova.2.cill, anova.2.ciul, anova.2.cov, anova.2.est.th, anova.2.pow.th, anova.2.pow.cal,
                                   anova.3.est, anova.3.se, anova.3.cill, anova.3.ciul, anova.3.cov, anova.3.est.th, anova.3.pow.th, anova.3.pow.cal,
                                   sanova.est, sanova.se, sanova.cill, sanova.ciul, sanova.cov, sanova.est.th, sanova.pow.th, sanova.pow.cal,
                                   manova.est, manova.se, manova.cill, manova.ciul, manova.cov, manova.est.th, manova.pow.th, manova.pow.cal,
                                   manova.wilks.th, manova.wilks.cal, manova.pillai.th, manova.pillai.cal,
                                   lv.l1.est.m1, lv.l1.se.m1, lv.l2.est.m1, lv.l2.se.m1, lv.l3.est.m1, lv.l3.se.m1,
                                   lv.r1.est.m1, lv.r1.se.m1, lv.r2.est.m1, lv.r2.se.m1, lv.r3.est.m1, lv.r3.se.m1, 
                                   lv.r4.est.m1, lv.r4.se.m1, lv.r5.est.m1, lv.r5.se.m1, lv.r6.est.m1, lv.r6.se.m1, 
                                   lv.rlv2.est.m1, lv.rlv2.se.m1,
                                   lv.i1.est.m1, lv.i1.se.m1, lv.i2.est.m1, lv.i2.se.m1, lv.i3.est.m1, lv.i3.se.m1, lv.i4.est.m1, lv.i4.se.m1,
                                   lv.chi.m1, 
                                   lv.est.m2, lv.se.m2,
                                   lv.l1.est.m2, lv.l1.se.m2, lv.l2.est.m2, lv.l2.se.m2, lv.l3.est.m2, lv.l3.se.m2,
                                   lv.r1.est.m2, lv.r1.se.m2, lv.r2.est.m2, lv.r2.se.m2, lv.r3.est.m2, lv.r3.se.m2, 
                                   lv.r4.est.m2, lv.r4.se.m2, lv.r5.est.m2, lv.r5.se.m2, lv.r6.est.m2, lv.r6.se.m2, 
                                   lv.rlv2.est.m2, lv.rlv2.se.m2,
                                   lv.i1.est.m2, lv.i1.se.m2, lv.i2.est.m2, lv.i2.se.m2, lv.i3.est.m2, lv.i3.se.m2, lv.chi.m2, lv.rmsea.m2, 
                                   lv.iter.m1, lv.iter.m2,
                                   count.conv, count.lv, count.ov)
  
  # write individual iterations to file
  colnames(results.cond)[1] <- "cond"; colnames(results.cond)[2] <- "lambda1"; colnames(results.cond)[3] <- "lambda2"; colnames(results.cond)[4] <- "lambda3"; colnames(results.cond)[5] <- "nu3"; colnames(results.cond)[6] <- "alpha2"; colnames(results.cond)[7] <- "n"
  results.cond <- as.data.frame(results.cond)
  #write.table(results.cond, file = paste("../intercept-results/results-individual/result_inter_", results.cond[1,"cond"], "_", results.cond[1,"lambda1"], "_", results.cond[1,"lambda2"], "_", results.cond[1,"lambda3"], "_", results.cond[1,"nu3"], "_", results.cond[1,"alpha2"], "_", results.cond[1,"n"], ".dat", sep = ""), row.names=F, col.names=T)
  write.table(results.cond, file = paste("C:/Users/s421506/tiu/research/mgcfa-intercepts/intercept-results/results-individual/result_inter_", results.cond[1,"cond"], "_", results.cond[1,"lambda1"], "_", results.cond[1,"lambda2"], "_", results.cond[1,"lambda3"], "_", results.cond[1,"nu3"], "_", results.cond[1,"alpha2"], "_", results.cond[1,"n"], ".dat", sep = ""), row.names=F, col.names=T)
  
  # assign NA to string variables first, otherwise mean cannot be calculated
  results.cond <- sapply(results.cond, as.numeric)
  
  # take mean of all iterations per condition 
  results.summarized <- as.numeric(apply(results.cond,2,mean,na.rm=T))

  return(results.summarized)
  
}

# objects that should be used for all clusters
clusterExport(cl, list("packages","conditions","myTryCatch","simdata","cfi","power.cfa.th","power.cfa.cal","power.sanova","power.anova","power.manova","results.summarized"))

# load library packages on all clusters
clusterEvalQ(cl,sapply(packages,library,character.only=T))                  

# run simulation over all clusters
results.summarized <- clusterApply(cl, 1:nrow(conditions), do_sim, cond = conditions, iters = iters) 

# shut down the nodes
stopCluster(cl) 

# transform results.summarized (now in list form) to data frame
results.summarized <- as.data.frame(do.call(rbind, results.summarized))

# assign colnames to dataframe
clmns <- c("cond", "lambda1", "lambda2", "lambda3", "nu3", "alpha2", "n", "lv.est.m1", "lv.se.m1", "lv.cill", "lv.ciul", "lv.cov", 
           "lv.cfi.est", "lv.rmsea.m1", "lv.ncp.th", "lv.pow.th", "lv.pow.cal", "anova.1.est", "anova.1.se", "anova.1.cill" , "anova.1.ciul", 
           "anova.1.cov", "anova.1.est.th", "anova.1.pow.th", "anova.1.pow.cal", "anova.2.est", "anova.2.se", "anova.2.cill" , "anova.2.ciul", 
           "anova.2.cov", "anova.2.est.th", "anova.2.pow.th", "anova.2.pow.cal", "anova.3.est", "anova.3.se", "anova.3.cill", "anova.3.ciul", 
           "anova.3.cov", "anova.3.est.th", "anova.3.pow.th", "anova.3.pow.cal", "sanova.est", "sanova.se" , "sanova.cill", "sanova.ciul", 
           "sanova.cov", "sanova.est.th" , "sanova.pow.th", "sanova.pow.cal", "manova.est", "manova.se", "manova.cill", "manova.ciul", 
           "manova.cov" , "manova.est.th" , "manova.pow.th", "manova.pow.cal", "manova.wilks.th", "manova.wilks.cal", "manova.pillai.th", 
           "manova.pillai.cal" , "lv.l1.est.m1" , "lv.l1.se.m1", "lv.l2.est.m1", "lv.l2.se.m1", "lv.l3.est.m1" , "lv.l3.se.m1", 
           "lv.r1.est.m1", "lv.r1.se.m1", "lv.r2.est.m1" , "lv.r2.se.m1", "lv.r3.est.m1" , "lv.r3.se.m1",  
           "lv.r4.est.m1" , "lv.r4.se.m1", "lv.r5.est.m1" , "lv.r5.se.m1","lv.r6.est.m1" , "lv.r6.se.m1",
           "lv.rlv2.est.m1", "lv.rlv2.se.m1" , 
           "lv.i1.est.m1", "lv.i1.se.m1", "lv.i2.est.m1", "lv.i2.se.m1" , "lv.i3.est.m1", "lv.i3.se.m1", 
           "lv.chi.m1" , 
           "lv.est.m2", "lv.se.m2", 
           "lv.l1.est.m2" , "lv.l1.se.m2", "lv.l2.est.m2", "lv.l2.se.m2" , "lv.l3.est.m2", "lv.l3.se.m2", 
           "lv.r1.est.m2", "lv.r1.se.m2", "lv.r2.est.m2", "lv.r2.se.m2", "lv.r3.est.m2" , "lv.r3.se.m2", 
           "lv.rlv2.est.m2", "lv.rlv2.se.m2", 
           "lv.i1.est.m2", "lv.i1.se.m2", "lv.i2.est.m2", "lv.i2.se.m2", "lv.i3.est.m2", "lv.i3.se.m2", "lv.chi.m2", "lv.rmsea.m2", 
           "lv.iter.m1", "lv.iter.m2", 
           "count.conv", "count.lv", "count.ov")

colnames(results.summarized) <- clmns

# order results.summarized by condition 
results.summarized <- results.summarized[order(results.summarized[,"cond"]),]

# write summarized results to file
write.table(results.summarized, file = "C:/Users/s421506/tiu/research/mgcfa-intercepts/intercept-results/results-individual/intercept_summarized_results.dat", row.names=F, col.names=T)

# if results.summarized is not saved correctly ----------------------------
# we need to load all individual files and summarize those.
# rm(list = ls())
# options(scipen=999)
# x = c("data.table")
# sapply(x,library,character.only=T)
# results.summarized <- c()
# 
# # Load all individual result files ----------------------------------------
# filelist = list.files(pattern = "result_")
# datalist = lapply(filelist, fread)
# 
# for (i in seq_along(datalist)){
#   
#   # assign NA to string variables first, otherwise mean cannot be calculated
#   datalist[[i]] <- sapply(datalist[[i]], as.numeric)
#   
#   # take mean of all iterations per condition
#   results.summarized.new <- apply(datalist[[i]],2,mean,na.rm=T)
#   
#   # assign it to summarized results
#   results.summarized <- rbind(results.summarized,results.summarized.new)
#   
# }
# 
# # assign column names to results.summarized and sort by cond
# colnames(results.summarized) <- colnames(datalist[[1]]) 
# results.summarized <- results.summarized[order(results.summarized[,1]),]
# 
# # write summarized results to file
# write.table(results.summarized, file = "idealworld_summarized_results1.dat", row.names=F, col.names=T)
