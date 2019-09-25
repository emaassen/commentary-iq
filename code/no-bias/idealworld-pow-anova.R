### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This function generates theoretical and sampled power estimates for all three 
### ANOVAs seperately, as well as an effect size, standard error, CIs, and coverage. 

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
  var.1 <- data$"sigma.sample.1"										                      # variance group 1 and 2
  var.2 <- data$"sigma.sample.2"
  pooled.sd <- lq.g1 * var.1 + lq.g2 * var.2 				              				# pooled variance
  pooled.sd <- pooled.sd/(lq.g1 + lq.g2)								                  # pooled SD
  sd.pool <- sqrt(diag(pooled.sd))
  
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