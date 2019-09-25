### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This function generates theoretical and sampled power estimates for sum ANOVAs,
### as well as an effect size, standard error, CIs, and coverage. 

power.sanova <- function(data, df, nobs, n, q, k) {
  
  ### Theoretical power
  
  # Calculating theoretical effect size
  lq.g1 <- nobs - 1
  lq.g2 <- nobs - 1
  var.p.1 <- sum(data$"sigma.pop.1") 									                  	# variance sumscore group 1 and 2
  var.p.2 <- sum(data$"sigma.pop.2")
  pooled.sd <- lq.g1 * var.p.1 + lq.g2 * var.p.2 				            			# pooled SD calclulation
  pooled.sd <- sqrt((pooled.sd/(lq.g1 + lq.g2))) 
  
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
  pooled.sd <- sqrt((pooled.sd/(lq.g1 + lq.g2)))
  
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