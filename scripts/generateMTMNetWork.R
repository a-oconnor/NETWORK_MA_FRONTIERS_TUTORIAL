## This network generate pair-wise for network plot 
## For each study it consists of the treatment in the first and the second arm of a
## comparison t1 and t2, respectively and the corresponding randomized sample size, n1 and n2(if 
## data type is Arm-level):
## If study contains more than two arms, includes all possible pair-wise comparison 
## Input:
##    MTCdata: MTCdata file
##    dataType: the type of your MTCdata, available options are “Arm” and “Contrast”
## Return: 
##    dataframe contains t1, t2 (Contrast data type) or t1, t2, n1, n2 (Arm data type) for network plotting 
generateMTMNetWork <- function(MTCdata, dataType = "Arm"){
  if( !(dataType %in% c("Arm", "Contrast") ) ){
    return(cat("Please input the correct data type, available options are \"Arm\" and \"Contrast\""))
  }
  maxArm <- max(MTCdata$Number.of.arms)
  if(dataType == "Arm"){
    ans <- data.frame(matrix(ncol = 4, nrow = 0))
    names(ans) <- c("t1", "t2", "n1", "n2")
    for(i in 2:maxArm){
      tmp_all <- MTCdata[MTCdata$Number.of.arms == i, ]
      for(j in 1:(i-1)){
        for(k in (j+1):i){
          tmp <- tmp_all[, c(paste0("Arm", j), paste0("Arm", k), paste0("Total.number.in.arm.", j), 
                             paste0("Total.number.in.arm.", k))]
          names(tmp) <- c("t1", "t2", "n1", "n2")
          ans <- rbind(ans, tmp)
        }
      }
    }
    ## order the data 
    ans <- ans[order(ans$n1), ]
  }
  if(dataType == "Contrast"){
    ans <- data.frame(matrix(ncol = 2, nrow = 0))
    names(ans) <- c("t1", "t2")
    for(i in 2:maxArm){
      tmp_all <- MTCdata[MTCdata$Number.of.arms == i, ]
      for(j in 1:(i-1)){
        for(k in (j+1):i){
          tmp <- tmp_all[, c(paste0("Arm", j), paste0("Arm", k))]
          names(tmp) <- c("t1", "t2")
          ans <- rbind(ans, tmp)
        }
      }
    }
  }
  rownames(ans) <- NULL
  return(ans)
}
    
    
    

  

