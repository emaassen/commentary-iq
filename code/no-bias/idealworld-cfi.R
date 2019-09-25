### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### First, an independence model is generated and fitted to the sampled data
### The correct model (model 1) is fitted to the sampled data
### Estimates used from both models are used to calculate the CFI estimate

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