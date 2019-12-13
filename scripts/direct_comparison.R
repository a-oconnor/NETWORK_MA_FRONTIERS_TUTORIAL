### direct comparison on each trial 
### Input: 
###   pairwise: pairwise comparison generated from pairwise_comp()
###   baselinePriorSmry: baseline prior summary generated from getBaselinePrior()
###   dataType: the type of your MTCdata, available options are “Arm” and “Contrast”
###   hyperSDInUnif: The prior distribution for heterogeneity parameter sigma is Unif(0, hyperSDInUnif)
###   hyperVarInNormMean: The prior distribution for the log odds is N(0, 1/hyperVarInNormMean)
###   parameters for running MCMC: i.e. n.adapt, n.iter, n.chains, n.thin
###   n.burnin_prop: The proportion you want to do burning for a chain.
###   save_res: logical, whether to save the output into your local directory
### Return: 
###   a dataframe containing direct comparisons summary
###   If save_res is True, also return a csv file contain direct comparison mean & sd: direct_comparison.csv

direct_comparison <- function(pairwise, baselinePriorSmry, dataType = "Arm", hyperSDInUnif = 2, 
                              hyperVarInNormMean = 1e-4, n.adapt = 5000, n.iter = 10000, 
                              n.chains = 3, n.burnin_prop = 0.5, n.thin = 1, save_res = T){
  if( !(dataType %in% c("Arm", "Contrast") ) ){
    return(cat("Please input the correct data type, available options are \"Arm\" and \"Contrast\""))
  }
  n.burnin <- floor(n.iter * n.burnin_prop)
  n.thin <-  max(1,floor((n.iter-n.burnin)/100000))
  compar <- unique(pairwise$Comparison)
  lc <- length(compar)
  ddir <- matrix(0, lc, 2)
  meanmA <- baselinePriorSmry[3,1] #mean of m in the baselinePriorSmry object
  precmA <- baselinePriorSmry[6,1] #mean of tau in the baselinePriorSmry object
  set.seed(20190710)
  auxt<-NULL
  if(dataType=="Arm"){
    
    for (i in 1:lc){
      idx <- which(pairwise$Comparison==compar[i])
      lidx <- length(idx)
      r <- matrix(1,lidx,2)
      n <- matrix(1,lidx,2)
      t <- matrix(1,lidx,2)
      r[,1] <- pairwise[idx,"Number.of.Event.in.arm.1"]
      r[,2] <- pairwise[idx,"Number.of.Event.in.arm.2"]
      n[,1] <- pairwise[idx,"Total.number.in.arm.1"]
      n[,2] <- pairwise[idx,"Total.number.in.arm.2"]
      t[,1] <- pairwise[idx,"Arm1"]
      t[,2] <- pairwise[idx,"Arm2"]
      t <- t(apply(t,1,order))
      na <- pairwise[idx,"rms"]
      NT <- 2
      NS <- lidx  
      
      if(NS==1){
        ## only one trial for this comparison 
        n <- n[1,]
        t <- t[1,]
        r <- r[1,]
        
        ## jags data 
        data.jags <- list(NT = NT, r = r, n = n, t = t, na = na, meanmA = meanmA, 
                          precmA = precmA, hyperVarInNormMean = hyperVarInNormMean,
                          hyperSDInUnif = hyperSDInUnif)
        
        ## jags initial value 
        init.jags <- list(sd = 0.2, d = c(NA,0))
        
        ## parameters to be monitored in jags 
        parms <- c("rk","T","best","d","lor","or","sd" )
        
        ## run jags  
        jags.m <- jags.model(file = "./scripts/jagsCodes/directComparison/MTC_direct1.bug",
                             data = data.jags,
                             inits = init.jags,
                             n.chains = n.chains, n.adapt = n.adapt)
        update(jags.m, n.iter = n.burnin)
        jags.out <- coda.samples(model=jags.m,
                                 variable.names=parms,
                                 n.iter=n.iter,thin=n.thin)
        ## gelman diagnostics 
        ## gelman.diag(jags.out)
        model_smry <- summary(jags.out)
      }
      else{
        ## jags data 
        data.jags <- list(NT = NT, NS = NS, r = r, n = n, t = t, na = na, meanmA = meanmA, 
                          precmA = precmA, hyperVarInNormMean = hyperVarInNormMean,
                          hyperSDInUnif = hyperSDInUnif)
        
        ## jags initial value 
        init.jags <- list(sd = 0.2, d = c(NA,0))
        
        ## parameters to be monitored in jags 
        parms <- c("rk","T","best","d","lor","or","sd" )
        
        ## run jags 
        jags.m <- jags.model(file = "./scripts/jagsCodes/directComparison/MTC_direct.bug",
                             data = data.jags,
                             inits = init.jags,
                             n.chains = n.chains, n.adapt = n.adapt)
        update(jags.m, n.iter = n.burnin)
        jags.out <- coda.samples(model=jags.m,
                                 variable.names=parms,
                                 n.iter=n.iter,thin=n.thin)
        ## gelman diagnostics 
        ## gelman.diag(jags.out)
        model_smry <- summary(jags.out)
      }
      
      idx <- which(rownames(model_smry$statistics) == "d[2]")
      ddir[i,1:2] <- model_smry$statistics[idx, c("Mean","SD")]
      auxt[[i]] <- jags.out
      names(auxt)[[i]] = compar[i]
    }
  }
  
  
  if(dataType=="Contrast"){
    
    for (i in 1:lc){
      idx <- which(pairwise$Comparison==compar[i])
      lidx <- length(idx)
      y <- matrix(1,lidx,2)
      se <- matrix(1,lidx,2)
      t <- matrix(1,lidx,2)
      y[,2] <- pairwise[idx,"lor.2"]
      se[,2] <- pairwise[idx,"se.2"]
      t[,1] <- pairwise[idx,"Arm1"]
      t[,2] <- pairwise[idx,"Arm2"]
      t <- t(apply(t,1,order))
      na <- pairwise[idx,"rms"]
      NT <- 2
      NS <- lidx  
      
      if(NS==1){
        ## only one trial for this comparison 
        y <- y[1,]
        t <- t[1,]
        se <- se[1,]
        
        ## jags data 
        data.jags <- list(NT = NT, y = y, se = se, t = t, meanmA = meanmA, 
                          precmA = precmA, hyperVarInNormMean = hyperVarInNormMean,hyperSDInUnif = hyperSDInUnif)
        
        ## jags initial value 
        init.jags <- list(sd = 0.2, d = c(NA,0))
        
        ## parameters to be monitored in jags 
        parms <- c("rk","T","best","d","lor","or","sd" )
        
        ## run jags  
        jags.m <- jags.model(file = "./scripts/jagsCodes/directComparison/MTC_direct1_cb.bug",
                             data = data.jags,
                             inits = init.jags,
                             n.chains = n.chains, n.adapt = n.adapt)
        update(jags.m, n.iter = n.burnin)
        jags.out <- coda.samples(model=jags.m,
                                 variable.names=parms,
                                 n.iter=n.iter,thin=n.thin)
        ## gelman diagnostics 
        ## gelman.diag(jags.out)
        model_smry <- summary(jags.out)
      }
      else{
        ## jags data 
        data.jags <- list(NT = NT, NS = NS, y = y, se = se, t = t,  meanmA = meanmA, 
                          precmA = precmA, hyperVarInNormMean = hyperVarInNormMean,hyperSDInUnif = hyperSDInUnif)
        
        ## jags initial value 
        init.jags <- list(sd = 0.2, d = c(NA,0))
        
        ## parameters to be monitored in jags 
        parms <- c("rk","T","best","d","lor","or","sd" )
        
        ## run jags 
        jags.m <- jags.model(file = "./scripts/jagsCodes/directComparison/MTC_direct_cb.bug",
                             data = data.jags,
                             inits = init.jags,
                             n.chains = n.chains, n.adapt = n.adapt)
        update(jags.m, n.iter = n.burnin)
        jags.out <- coda.samples(model=jags.m,
                                 variable.names=parms,
                                 n.iter=n.iter,thin=n.thin)
        ## gelman diagnostics 
        ## gelman.diag(jags.out)
        model_smry <- summary(jags.out)
      }
      
      idx <- which(rownames(model_smry$statistics) == "d[2]")
      ddir[i,1:2] <- model_smry$statistics[idx, c("Mean","SD")]
      auxt[[i]] <- jags.out
      names(auxt)[[i]] = compar[i]
    }
  }
  res_smry <- NULL
  res_smry$compar <- compar
  res_smry$Mean <- ddir[,1]
  res_smry$SD <- ddir[,2]
  if(save_res == T){
    ## save jags summary output
    save(model_smry, file = "./scripts/RDataObjects/directComparisonSummary.RData")
    ## jags output 
    save(jags.out, file = "./scripts/RDataObjects/directComparison.RData")
    write.table(res_smry ,file="./data/direct_comparison.csv",sep=",",
                row.names=FALSE, quote=FALSE, col.names=TRUE)
  }
  return(as.data.frame(res_smry))
}


