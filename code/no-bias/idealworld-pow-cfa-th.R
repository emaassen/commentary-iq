### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This function generates theoretical power estimates for the CFA model,
### as well as an effect size, standard error, CIs, coverage, Chisquare and RMSEA values. 

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