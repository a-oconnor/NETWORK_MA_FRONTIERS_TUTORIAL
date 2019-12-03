### Input:
###   MTCdata: MTCdata file
###   map_txt: a dataframe which maps each treatment to a number
###   arm_name: select the treatment name type to output, available options are “abbr” or “full”
###   save_res: logical, whether to save the output into your local directory
### Return:
###   string for mtm network node names: i.e. treatment name (associated no. of trials)
###   Also save csv file treatmentFrequencySummary.csv if save_res is true
mtmNetWorkNodeName <- function(MTCdata, map_txt, arm_name = "abbr", save_res = F){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  if(arm_name == "abbr"){
  treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(abbr) %>% dplyr::pull()
  }
  if(arm_name == "full"){
    treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(name) %>% dplyr::pull()
  }
  ## no. of treatment 
  NT <- length(treatNms)
  ## treatment frequency - the col numbers are the cols with the numerical txt indictors 
  maxArm <- max(MTCdata$Number.of.arms) ## get the maximum number of arms in a study
  ntFreq <- NULL
  for(i in 1:NT){
    ntFreq[i] <- (MTCdata[,c(paste0("Arm", 1:maxArm))] == i) %>% sum(na.rm = T)
  }
  ## combine treatment name with frequency 
  treat_new <- ""
  for(i in 1:NT){
    treat_new[i] <- paste(treatNms[i], "(", ntFreq[i], ")", sep = "")
  }
  
  ## write treatment frequency out 
  num_trea <- cbind(treatNms,ntFreq)
  if(save_res == T){
    write.csv(num_trea, file = "./data/treatmentFrequencySummary.csv", row.names = FALSE)
  }
  numTrials <- data.frame(matrix(nrow = 1, ncol = maxArm-1))
  for(i in 2:maxArm){
    numTrials[1,(i-1)] <- sum(MTCdata$Number.of.arms == i)
  }
  colnames(numTrials) <- paste0("num", 2:maxArm, "ArmsTrials")
  
  return(list(treat_new = treat_new, 
              treat_freq = ntFreq,
              totalTrialArms  = sum(ntFreq),
              numTrials = numTrials))
}
