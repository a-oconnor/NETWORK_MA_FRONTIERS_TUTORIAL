#################### run general model:random effects model
## Input: 
##     Datasets needed:
##     MTCdata: MTCdata file
##     baselinePriorSmry: baseline prior summary generated from getBaselinePrior()
##     Options needed:
##     good_event: Indicator, 1 if the events are good, 0 if the events are bad, would influence the rank
##     dataType: the type of your MTCdata, available options are “Arm” and “Contrast”
##     Values needed: defaults rarely change
##     hyperSDInUnif: The prior distribution for heterogeneity parameter sigma is Unif(0, hyperSDInUnif)
##     hyperVarInNormMean: The prior distribution for the log odds is N(0, 1/hyperVarInNormMean)
##     parameters for running MCMC: i.e. n.adapt, n.iter, n.chains, n.thin
##     n.burnin_prop: The proportion you want to do burning for a chain.
## Return:
##     mcmc simulation and its summary in RData object
##     model deviance information: generalModeldev_result.RData
##     A dataframe with general model summary
################################################################################
general_model <- function(MTCdata, baselinePriorSmry, good_event = 0, dataType = "Arm",
                         hyperSDInUnif = 5, hyperVarInNormMean = 1e-4, n.adapt = 5000,
                         n.iter = 10000, n.chains = 3, n.burnin_prop = 0.5, n.thin = 1){
  if( !(dataType %in% c("Arm", "Contrast") ) ){
    return(cat("Please input the correct data type, available options are \"Arm\" and \"Contrast\""))
  }
  n.burnin <- floor(n.iter * n.burnin_prop)
  n.thin <-  max(1,floor((n.iter-n.burnin)/100000))
  meanmA <- baselinePriorSmry[3,1] #mean of m in the baselinePriorSmry dataframe
  precmA <- baselinePriorSmry[6,1] #mean of tau in the baselinePriorSmry dataframe
  
  ## seed 
  set.seed(20190914)
  
  nc <- max(MTCdata$Number.of.arms) ## maximum no.treatment in each trial over all trails 
  MTCdata <- as.data.frame(MTCdata)
  MTCdata <- MTCdata[order(MTCdata$Number.of.arms),]
  nr <- dim(MTCdata)[1]  ## no.of studies  
  na <- MTCdata[,"Number.of.arms"] ## no. of arms on each trial 

  
  ## parameters to be monitored in jags 
  parms <- c("rk","T","best","d","lor","or","sd","resdev","RRa","RR","totresdev")
  
  if(dataType == "Arm"){
    r <- matrix(1, nr, nc)
    n <- matrix(1, nr, nc)
    t <- matrix(1, nr, nc)
    
    for(i in 1:nc){
      rColName <- paste("Number.of.Event.in.arm.", i, sep = "")
      r[,i] <- round(as.numeric(MTCdata[,rColName]))
      nColName <- paste("Total.number.in.arm.", i, sep = "")
      n[,i] <- round(as.numeric(MTCdata[,nColName]))
      tColName <- paste("Arm", i, sep = "")
      t[,i] <- round(as.numeric(MTCdata[,tColName]))
    }
    
    NS <- nr ## no. of trials
    NT <- length(unique(as.vector(t[!is.na(t)]))) ## no. of treatment after removing NA
    
    ## jags initial value
    init.jags <- list(sd = 0.2, d = c(NA, rep(0,(NT-1))))

    ## jags data
    data.jags <- list(hyperSDInUnif = hyperSDInUnif, hyperVarInNormMean = hyperVarInNormMean,
                      meanmA = meanmA, precmA = precmA,
                      r = r, n = n, t = t,
                      NS = NS, NT = NT, na = na, good_event = good_event
    )
    
    
    ## run jags 
    jags.m <- jags.model(file = "./scripts/jagsCodes/generalModel/MTC_generalModel.bug",
                         data = data.jags, inits = init.jags,
                         n.chains = n.chains, n.adapt = n.adapt)
  }
  
  if(dataType == "Contrast"){
    y <- matrix(1,nr,nc)
    se <- matrix(1,nr,nc)
    t <- matrix(1, nr, nc)
    V <- as.numeric(MTCdata$V)
    NS <- rep(0,4)
    NS[1] <- nr
    
    for(i in 1:nc){
      tColName <- paste("Arm", i, sep = "")
      t[,i] <- round(as.numeric(MTCdata[,tColName]))
      if(i > 1){
        yColName <- paste("lor.",i,sep = "")
        y[,i] <- as.numeric(MTCdata[,yColName])
        seColName <- paste("se.", i, sep = "")
        se[,i] <- as.numeric(MTCdata[,seColName])
        NS[i] <- sum(MTCdata$Number.of.arms==i) ## NS[1] is the total number of studies, NS[i] (i>1) is the total number of i-arm studies
      }
    }
    NT <- length(unique(as.vector(t[!is.na(t)]))) ## no. of treatment after removing NA
    
    ## jags initial value
    init.jags <- list(sd = 0.2, d = c(NA, rep(0,(NT-1))))
    
    ## jags data
    data.jags <- list(hyperSDInUnif = hyperSDInUnif, hyperVarInNormMean = hyperVarInNormMean,
                      meanmA = meanmA, precmA = precmA,
                      y = y, se = se, t = t, V = V,
                      NS = NS, NT = NT, na = na, good_event = good_event)
    
    ## run jags 
    jags.m <- jags.model(file = "./scripts/jagsCodes/generalModel/MTC_generalModel_cb.bug",
                         data = data.jags, inits = init.jags,
                         n.chains = n.chains, n.adapt = n.adapt)
    
  }
  update(jags.m, n.iter = n.burnin)
  jags.out <- coda.samples(model=jags.m,variable.names=parms,n.iter=n.iter,thin=n.thin)
  general_model_smry <- summary(jags.out)
  general_model_smry <- cbind(general_model_smry$statistics[,c("Mean","SD")],
                              general_model_smry$quantiles[,c("2.5%","50%","97.5%")])
  
  ## save summary output smry as a RData object and csv file
  save(general_model_smry, file = "./scripts/RDataObjects/generalModelSummary.RData")
  write.csv(as.data.frame(round(general_model_smry,3)), file = "./data/general_model_smry.csv", row.names = F)
  ## save jags output 
  save(jags.out, file = "./scripts/RDataObjects/generalModelJagasOutput.RData")
  
  ## save deviance result 
  dev_result <- dic.samples(jags.m, n.iter = n.burnin, type = 'pD')
  save(dev_result, file = "./scripts/RDataObjects/generalModeldev_result.RData")
  
  ## return general model summary as a dataframe and jags output
  return(list(general_model_smry = as.data.frame(general_model_smry), jags.out = jags.out))
}
