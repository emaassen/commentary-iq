### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This function generates theoretical and sampled power estimates for MANOVA,
### as well as an effect size, standard error, CIs, and coverage. 

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

