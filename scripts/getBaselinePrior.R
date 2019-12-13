########## get hyperprior value for baseline treatment in logit scale by seperate mcmc simulation 
##  baseline treatment: non active control, coded as arm 1 in the dataset 
##  Input:
##   Datasets needed:
##    MTCdata: MTCdata file
##    Options needed:
##    dataType: the type of your MTCdata, available options are “Arm” and “Contrast”
##    save_res: logical, whether to save the output into your local directory
##    Values needed: defaults rarely change
##    hyperSDInUnif: The prior distribution for heterogeneity parameter sigma is Unif(0, hyperSDInUnif)
##    hyperVarInNormMean: The prior distribution for the log odds is N(0, 1/hyperVarInNormMean)
##    parameters for running MCMC: i.e. n.adapt, n.iter, n.chains, n.thin
##    n.burnin_prop: The proportion you want to do burning for a chain.
##  Return:
##    smry.RData object contains the summary statistics from mcmc simulation 
getBaselinePrior <- function(MTCdata, dataType = "Arm", hyperSDInUnif = 5, hyperVarInNormMean = 1e-4, n.adapt = 5000,
                             n.iter = 10000, n.chains = 3, n.burnin_prop = 0.5, n.thin = 1, save_res = T){
  n.burnin <- floor(n.iter * n.burnin_prop)
  n.thin <-  max(1,floor((n.iter-n.burnin)/100000))
  
  ## jags initial value 
  init.jags <- list(sd.m = 0.2)
  
  ## parameters to be monitored in jags 
  parms <- c("R","R.new","m","mu.new","sd.m","tau.m")
  
  set.seed(20190912)
  if(dataType == "Arm"){
    non_active <- MTCdata[MTCdata$Arm1 == 1, ]
    r <- unlist(non_active[,"Number.of.Event.in.arm.1"])
    n <- unlist(round(non_active[,"Total.number.in.arm.1"],0))
    ns <- dim(non_active)[1]
    ## jags data 
    data.jags <- list(hyperSDInUnif = hyperSDInUnif, hyperVarInNormMean = hyperVarInNormMean,
                      r = r, n = n, ns = ns)
    ## run jags 
    jags.m <- jags.model(file = "./scripts/jagsCodes/getBaselinePrior.bug",
                         data = data.jags, inits = init.jags,
                         n.chains = n.chains, n.adapt = n.adapt)
  }
  
  if(dataType == "Contrast"){
    non_active <- as.numeric(unlist(MTCdata[,"PLA.lo"]))
    mu <- non_active[!is.na(non_active)]
    ns <- length(mu)
    ## jags data 
    data.jags <- list(hyperSDInUnif = hyperSDInUnif, hyperVarInNormMean = hyperVarInNormMean,
                      mu=mu, ns = ns)
    ## run jags 
    jags.m <- jags.model(file = "./scripts/jagsCodes/getBaselinePrior_cb.bug",
                         data = data.jags, inits = init.jags,
                         n.chains = n.chains, n.adapt = n.adapt)
  }
  update(jags.m, n.iter = n.burnin)
  jags.out <- coda.samples(model=jags.m,variable.names=parms,n.iter=n.iter,thin=n.thin)
  smry <- summary(jags.out)
  baselinePriorSmry <- cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")])
  if(save_res == T){
    ## save summary output smry as a RData object 
    save(baselinePriorSmry, file = "./scripts/RDataObjects/baselinePriorSummary.RData")
  }
  return(baselinePriorSmry)
}
