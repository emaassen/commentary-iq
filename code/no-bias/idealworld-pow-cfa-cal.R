### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This function generates calculated power estimates for the CFA model,
### as well as an effect size, standard error, CIs, coverage, Chisquare and RMSEA values. 


power.cfa.cal <- function(model1, model2, data, nobs) {
  
  # To test this function step by step:
  #data <- datalist
  
  # CALCULATED POWER
  # Sample variance covariance structure and sample mean vector
  sigma.sample <- list(data$"sigma.sample.1",data$"sigma.sample.2")
  mu.sample <- list(data$"mu.sample.1",data$"mu.sample.2")
  
  # Fit model1 (estimating intercept) and model2 (restricting intercept)
  mod.cal.1 <- myTryCatch(cfa(model1, sample.cov = sigma.sample, sample.mean = mu.sample, sample.nobs = c(nobs, nobs), group.equal = c("loadings", "intercepts")))
  mod.cal.2 <- myTryCatch(cfa(model2, sample.cov = sigma.sample, sample.mean = mu.sample, sample.nobs = c(nobs, nobs), group.equal = c("loadings", "intercepts")))
  
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