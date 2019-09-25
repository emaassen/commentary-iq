### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
rm(list = ls())                                               # remove all objects from current workspace
options(scipen=999, warn=1)                                   # no scientific notation + print all warnings
packages <- c("parallel","MASS","lavaan","psych")             # all packages needed to run simulation
#apply(packages,install.packages(packages),character.only=T)  # install all packages                      
sapply(packages,library,character.only=T)                     # load all packages
nworkers <- detectCores() -1                                  # number of cores computer - needed for parallel processing
cl <- makePSOCKcluster(nworkers)                              # create cluster based on number of cores -1 - needed for parallel processing

source("idealworld-cfi.R")
source("idealworld-mytrycatch.R")
source("idealworld-pow-anova.R")
source("idealworld-pow-cfa-cal.R")
source("idealworld-pow-cfa-th.R")
source("idealworld-pow-manova.R")
source("idealworld-pow-sanova.R")
source("idealworld-simdata.R")

iters <- 1000                                                               # no of iterations
results.summarized <- c()                                                   # empty dataframe to store all summarized results per condition

lambdas <- list(c(0.1,0.1,0.1), c(0.1,0.1,0.3), c(0.1,0.1,0.5),             # list with all variations of lambda (factor loadings) per 3 indicators
                c(0.1,0.1,0.7), c(0.1,0.3,0.3), c(0.1,0.3,0.5),
                c(0.1,0.3,0.7), c(0.1,0.5,0.5), c(0.1,0.5,0.7),
                c(0.1,0.7,0.7), c(0.2,0.2,0.2), c(0.3,0.3,0.3),
                c(0.3,0.3,0.5), c(0.3,0.3,0.7), c(0.3,0.5,0.5),
                c(0.3,0.5,0.7), c(0.3,0.7,0.7), c(0.4,0.4,0.4),
                c(0.5,0.5,0.5), c(0.5,0.5,0.7), c(0.5,0.7,0.7), 
                c(0.6,0.6,0.6), c(0.7,0.7,0.7), c(0.8,0.8,0.8), 
                c(0.9,0.9,0.9), c(0.8,0.7,0.6), c(0.7,0.6,0.5))

alpha <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)                                             # all alphas (latent means for group 2)

nobs <- c(25, 50, 75, 100, 150, 200, 250, 300, 400, 500)                             # all sample sizes

phi <- c(1)

nus <- list(c(0,0,0))

cond <- expand.grid(nobs=nobs, lambdas = lambdas, nus = nus, alpha = alpha, phi = phi) 

# it's possible to reproduce only one condition or subsets of conditions, eg only no 430, and to vary in the no of iterations, eg:
#cond <- cond[c(741:743),]
#iters <- 5

# Simulation function -----------------------------------------------------
do_sim = function(pos, cond, iters) {
  
  # to run this function step by step, run this line:
  #pos <- 1
  
  # Fixed values for all simulations
  n <- 1            									 	                                 	# no LVs       
  q <- 3            									 					                         	# no indicators    
  k <- 2            									 					                         	# no groups
  
  # Select condition by changing pos and selecting conditions from cond
  lambda1 <- cond$lambdas[pos][[1]][1]
  lambda2 <- cond$lambdas[pos][[1]][2]
  lambda3 <- cond$lambdas[pos][[1]][3]
  alpha <- cond$alpha[pos]
  nobs <- cond$nobs[pos]
  
  set.seed(as.numeric(rownames(cond)[pos]))															
  
  # Empty vectors & matrices ------------------------------------------------
  lv.est.m1 <- lv.se.m1 <- lv.cill <- lv.ciul <- lv.cov <- lv.rmsea.m1 <- lv.ncp.th <- lv.pow.th <- lv.pow.cal <- c()
  lv.l1.est.m1 <- lv.l1.se.m1 <- lv.l2.est.m1 <- lv.l2.se.m1 <- lv.l3.est.m1 <- lv.l3.se.m1 <- c()
  lv.r1.est.m1 <- lv.r1.se.m1 <- lv.r2.est.m1 <- lv.r2.se.m1 <- lv.r3.est.m1 <- lv.r3.se.m1 <- c()
  lv.r4.est.m1 <- lv.r4.se.m1 <- lv.r5.est.m1 <- lv.r5.se.m1 <- lv.r6.est.m1 <- lv.r6.se.m1 <- c()
  lv.rlv2.est.m1 <- lv.rlv2.se.m1 <- c() 
  lv.i1.est.m1 <- lv.i1.se.m1 <- lv.i2.est.m1 <- lv.i2.se.m1 <- lv.i3.est.m1 <- lv.i3.se.m1 <- c()
  lv.chi.m1 <- c()
  lv.iter.m1 <- c()
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
  count.conv <- warn.m1 <- warn.m2 <- err.m1 <- err.m2 <- c()
  
  for (i in 1:iters) {
    
    datalist <- simdata(lambda1,lambda2,lambda3,alpha,nobs,n,q,k)		  	# list with all data, covar matrices and mean structures
    df <- datalist[["df"]]											                    		# raw data
    
    # MGCFA -------------------------------------------------------------
    
    # Make sure model 1 = model with LV mean difference estimated, model 2 = no LV mean difference
    
    # Model where one latent mean difference is estimated
    model1 <- ' !Estimate all loadings
    F =~ NA*V1 + NA*V2 + NA*V3
    
    !Fix one LV variance, estimate one
    F ~~ c(1,NA) * F
    
    !Fix one LV mean to zero, estimate one
    F ~ c(0,NA) * 1 '
    
    # Covariances = 18
    # Parameters = 3 loadings, 3 intercepts, 3 residuals, 1 LV residual variance, 1 LV mean = 11 (df = 7)
    
    # Model where there is no latent mean difference assumed between groups
    model2 <- ' !Estimate all loadings
    F =~ NA*V1 + NA*V2 + NA*V3
    
    !Fix one LV variance, estimate one
    F ~~ c(1,NA) * F
    
    !Fix both LV means to zero, estimate none
    F ~ c(0,0) * 1'
    
    # Covariances = 18
    # Parameters = 3 loadings, 3 intercepts, 3 residuals, 1 LV residual variance = 10 (df = 8)
    
    
    # If you want to estimate both latent residual variances and fix the first loading to 1, replace model1 and 2 by:
    
    #model1 <- ' F =~ V1 + V2 + V3
    #F ~ c(0,NA) * 1'
    
    # Covariances = 18
    # Parameters = 2 loadings, 3 intercepts, 3 residuals, 2 LV residual variances, 1 LV mean = 11 (df = 7)
    
    #model2 <- ' F =~ V1 + V2 + V3
    #F ~ c(0,0) * 1'
    
    # Covariances = 18
    # Parameters = 2 loadings, 3 intercepts, 3 residuals, 2 LV residual variances = 10 (df = 8)
    
    # Calculate NCP, theoretical power, calculated power, and give model1 and model2 fit, chi square and rmsea estimates and possible warning and error messages
    pow.lv.th <- power.cfa.th(model2 = model2, data = datalist, nobs = nobs)
    
    pow.lv.cal <- power.cfa.cal(model1 = model1, model2 = model2, data = datalist, nobs = nobs)
    
    # Calculate CFI and save possible warning and error messages
    cfi.res <- cfi(data = datalist, nobs = nobs, pow.lv = pow.lv.cal)
    
    # LV estimates ------------------------------------------------------------
    # Model 1 estimates (with mean difference estimated)
    lv.est.m1 <- c(lv.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"est"]) 	# mean estimate group 2 (mean estimate group 1 restricted to 0)
    lv.se.m1 <- c(lv.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"se"])  	  # mean estimate se
    lv.cill <- c(lv.cill,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.lower"]) # mean estimate CI 
    lv.ciul <- c(lv.ciul,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.upper"])
    lv.cov <- c(lv.cov,as.numeric(parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.lower"] < alpha & parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[16,"ci.upper"] > alpha))							 # mean estimate coverage
    lv.cfi.est <- c(lv.cfi.est,cfi.res[[1]])  											                    # cfi estimate
    lv.rmsea.m1 <- c(lv.rmsea.m1,unlist(pow.lv.cal)$"lv.rmsea.1") 							            # rmsea model 1
    lv.ncp.th <- c(lv.ncp.th,unlist(pow.lv.th)[[1]]) 									                  	# noncentrality parameter used to calculate theoretical power
    lv.pow.th <- c(lv.pow.th,unlist(pow.lv.th)[[2]]) 								                	# theoretical power estimate
    lv.pow.cal <- c(lv.pow.cal,as.numeric(unlist(pow.lv.cal)$"pow.cal")) 				          	# calculated power / type 1 error
    
    # Additional estimates
    lv.l1.est.m1 <- c(lv.l1.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[1,"est"]) # factor loading estimates + standard errors
    lv.l1.se.m1 <- c(lv.l1.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[1,"se"])
    lv.l2.est.m1 <- c(lv.l2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[2,"est"])
    lv.l2.se.m1 <- c(lv.l2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[2,"se"])
    lv.l3.est.m1 <- c(lv.l3.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[3,"est"])
    lv.l3.se.m1 <- c(lv.l3.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[3,"se"])
    lv.r1.est.m1 <- c(lv.r1.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[6,"est"]) # residual estimates + standard errors
    lv.r1.se.m1 <- c(lv.r1.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[6,"se"])
    lv.r2.est.m1 <- c(lv.r2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[7,"est"])
    lv.r2.se.m1 <- c(lv.r2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[7,"se"])
    lv.r3.est.m1 <- c(lv.r3.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[8,"est"])
    lv.r3.se.m1 <- c(lv.r3.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[8,"se"])
    lv.r4.est.m1 <- c(lv.r4.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[17,"est"])
    lv.r4.se.m1 <- c(lv.r4.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[17,"se"])
    lv.r5.est.m1 <- c(lv.r5.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[18,"est"])
    lv.r5.se.m1 <-c(lv.r5.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[18,"se"])
    lv.r6.est.m1 <- c(lv.r6.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[19,"est"])
    lv.r6.se.m1 <- c(lv.r6.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[19,"se"])
    lv.rlv2.est.m1 <- c(lv.rlv2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[15,"est"]) # group 2 latent variable variance estimate + standard errors (variance LV group 1 = restricted to 1)
    lv.rlv2.se.m1 <- c(lv.rlv2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[15,"se"])
    lv.i1.est.m1 <- c(lv.i1.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[9,"est"]) # intercept estimates + standard errors
    lv.i1.se.m1 <- c(lv.i1.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[9,"se"])
    lv.i2.est.m1 <- c(lv.i2.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[10,"est"])
    lv.i2.se.m1 <- c(lv.i2.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[10,"se"])
    lv.i3.est.m1 <- c(lv.i3.est.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[11,"est"])
    lv.i3.se.m1 <- c(lv.i3.se.m1,parameterEstimates(unlist(pow.lv.cal)$mod.1.value)[11,"se"])
    lv.chi.m1 <- c(lv.chi.m1,unlist(pow.lv.cal)$"lv.chi.1")									              	      # chisquare model 1
    lv.iter.m1 <- c(lv.iter.m1,unlist(pow.lv.cal)$"mod.1.iter")                                   # iterations model 1
    
    # Convergence counter + Save warnings and error messages from model 1 and 2
    count.conv <- c(count.conv, as.numeric(lavTech(unlist(pow.lv.cal)$mod.1.value, what = "converged")))
    warn.m1 <- c(warn.m1, unlist(pow.lv.cal)$"mod.1.warn")
    err.m1 <- c(err.m1, unlist(pow.lv.cal)$"mod.1.err")
    warn.m2 <- c(warn.m2, unlist(pow.lv.cal)$"mod.2.warn")
    err.m2 <- c(err.m2, unlist(pow.lv.cal)$"mod.2.err")
    
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
  results.cond <- cbind.data.frame(rep(as.numeric(rownames(cond)[pos]),iters), rep(lambda1,iters), rep(lambda2,iters), rep(lambda3,iters), rep(alpha,iters), rep(nobs,iters),
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
                                   lv.i1.est.m1, lv.i1.se.m1, lv.i2.est.m1, lv.i2.se.m1, lv.i3.est.m1, lv.i3.se.m1, lv.chi.m1, 
                                   lv.iter.m1, 
                                   count.conv, warn.m1, err.m1, warn.m2, err.m2)
  
  # write individual iterations to file
  colnames(results.cond)[1] <- "cond"; colnames(results.cond)[2] <- "lambda1"; colnames(results.cond)[3] <- "lambda2"; colnames(results.cond)[4] <- "lambda3"; colnames(results.cond)[5] <- "alpha2"; colnames(results.cond)[6] <- "n"
  results.cond <- as.data.frame(results.cond)
  
  if(i == 1){
    write.table(results.cond, file = paste("../results/res_ideal_", results.cond[1,"cond"], "_", results.cond[1,"lambda1"], "_", results.cond[1,"lambda2"], "_", results.cond[1,"lambda3"], "_", results.cond[1,"alpha2"], "_", results.cond[1,"n"], ".dat", sep = ""), row.names=F, col.names=T)
  }
  
  if(i > 1){
    write.table(results.cond, file = paste("../results/res_ideal_", results.cond[1,"cond"], "_", results.cond[1,"lambda1"], "_", results.cond[1,"lambda2"], "_", results.cond[1,"lambda3"], "_", results.cond[1,"alpha2"], "_", results.cond[1,"n"], ".dat", sep = ""), row.names=F, col.names=T, append=T)
  }
  
  # assign NA to string variables first, otherwise mean cannot be calculated
  #results.cond <- sapply(results.cond, as.numeric)
  
  # remove the last four columns (i.e., string variables of warnings and errors) because mean estimate is not accurate
  results.cond <- results.cond[,1:(length(results.cond)-4)]
  
  # take mean of all iterations per condition 
  results.summarized <- as.numeric(apply(results.cond,2,mean,na.rm=T))

  return(results.summarized)
  
}

# objects that should be used for all clusters
clusterExport(cl, list("packages","cond","myTryCatch","simdata","cfi","power.cfa.th","power.cfa.cal","power.sanova","power.anova","power.manova","results.summarized"))

# load library packages on all clusters
clusterEvalQ(cl,sapply(packages,library,character.only=T))                  

# run simulation over all clusters
results.summarized <- clusterApply(cl, 1:nrow(cond), do_sim, cond = cond, iters = iters) 

# shut down the nodes
stopCluster(cl) 

# transform results.summarized (now in list form) to data frame
results.summarized <- as.data.frame(do.call(rbind, results.summarized))

# assign colnames to dataframe
clmns <- c("cond", "lambda1", "lambda2", "lambda3", "alpha2", "n", 
           "lv.est.m1", "lv.se.m1", "lv.cill", "lv.ciul", "lv.cov", "lv.cfi.est", "lv.rmsea.m1", "lv.ncp.th", "lv.pow.th", "lv.pow.cal", 
           "anova.1.est", "anova.1.se", "anova.1.cill", "anova.1.ciul", "anova.1.cov", "anova.1.est.th", "anova.1.pow.th", "anova.1.pow.cal", 
           "anova.2.est", "anova.2.se", "anova.2.cill", "anova.2.ciul", "anova.2.cov", "anova.2.est.th", "anova.2.pow.th", "anova.2.pow.cal",
           "anova.3.est", "anova.3.se", "anova.3.cill", "anova.3.ciul",  "anova.3.cov", "anova.3.est.th", "anova.3.pow.th", "anova.3.pow.cal", 
           "sanova.est", "sanova.se" , "sanova.cill", "sanova.ciul", "sanova.cov", "sanova.est.th", "sanova.pow.th", "sanova.pow.cal", 
           "manova.est", "manova.se", "manova.cill", "manova.ciul", "manova.cov", "manova.est.th", "manova.pow.th", "manova.pow.cal", 
           "manova.wilks.th", "manova.wilks.cal", "manova.pillai.th", "manova.pillai.cal", 
           "lv.l1.est.m1", "lv.l1.se.m1", "lv.l2.est.m1", "lv.l2.se.m1", "lv.l3.est.m1", "lv.l3.se.m1", 
           "lv.r1.est.m1",  "lv.r1.se.m1", "lv.r2.est.m1", "lv.r2.se.m1", "lv.r3.est.m1", "lv.r3.se.m1", 
           "lv.r4.est.m1",  "lv.r4.se.m1", "lv.r5.est.m1", "lv.r5.se.m1", "lv.r6.est.m1", "lv.r6.se.m1",
           "lv.rlv2.est.m1", "lv.rlv2.se.m1", 
           "lv.i1.est.m1",  "lv.i1.se.m1", "lv.i2.est.m1", "lv.i2.se.m1", "lv.i3.est.m1", "lv.i3.se.m1", "lv.chi.m1", 
           "lv.iter.m1", "count.conv")

colnames(results.summarized) <- clmns

# order results.summarized by condition 
results.summarized <- results.summarized[order(results.summarized[,"cond"]),]

# write summarized results to file
write.table(results.summarized, file = "../results/idealworld_summarized_results.dat", row.names=F, col.names=T)

# write session information to file
writeLines(capture.output(sessionInfo()), "../results/idealworld-sessioninfo.txt")

# if results.summarized is not saved correctly ----------------------------
# we need to load all individual files and summarize those.
# rm(list = ls())
# options(scipen=999)
# x = c("data.table")
# sapply(x,library,character.only=T)
# results.summarized <- c()
# 
# # Load all individual result files 
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
