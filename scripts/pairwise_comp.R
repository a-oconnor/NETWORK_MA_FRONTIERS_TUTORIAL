## Input:
##     Datasets needed:
##     MTCdata: MTCdata file
##     Options needed:
##     dataType: the type of your MTCdata, available options are “Arm” and “Contrast”
##     save_res: logical, whether to save the output into your local directory
## Return:
##     A dataframe containing the pariwise comparison
##     Also save a .csv file pairwise_comparison.csv if save_res is true

pairwise_comp <- function(MTCdata, dataType, save_res = T){
  if( !(dataType %in% c("Arm", "Contrast") ) ){
    return(cat("Please input the correct data type, available options are \"Arm\" and \"Contrast\""))
  }
  maxArm <- max(MTCdata$Number.of.arms)
  if(dataType == "Arm"){
    ans <- data.frame(matrix(ncol = 8, nrow = 0))
    names(ans) <- c("Number.of.Event.in.arm.1", "Number.of.Event.in.arm.2", "Total.number.in.arm.1", "Total.number.in.arm.2","Arm.1","Arm.2","Arm1","Arm2")
    for(i in 2:maxArm){
      tmp_all <- MTCdata[MTCdata$Number.of.arms == i, ]
      for(j in 1:(i-1)){
        for(k in (j+1):i){
          tmp <- tmp_all[, c(paste0("Number.of.Event.in.arm.",j), paste0("Number.of.Event.in.arm.",k), 
                             paste0("Total.number.in.arm.", j), paste0("Total.number.in.arm.", k), 
                         paste0("Arm.", j), paste0("Arm.", k), paste0("Arm", j), paste0("Arm", k))]
          names(tmp) <- c("Number.of.Event.in.arm.1", "Number.of.Event.in.arm.2", 
                          "Total.number.in.arm.1", "Total.number.in.arm.2","Arm.1","Arm.2","Arm1","Arm2")
          ans <- rbind(ans, tmp)
        }
      }
    }
  }
  if(dataType == "Contrast"){
    ans <- data.frame(matrix(ncol = 8, nrow = 0))
    names(ans) <- c("lor.2","se.2","Arm.1","Arm.2","Arm1","Arm2")
    for(i in 2:maxArm){
      tmp_all <- MTCdata[MTCdata$Number.of.arms == i, ]
      for(j in 1:(i-1)){
        for(k in (j+1):i){
          if(j > 1){
            tmp <- data.frame(matrix(ncol = 6, nrow = nrow(tmp_all)))
            names(tmp) <- c("lor.2","se.2","Arm.1","Arm.2","Arm1","Arm2")
            tmp$lor.2 <- (tmp_all[, paste0("lor.",k)] - tmp_all[, paste0("lor.",j)])  %>% t() %>% as.vector()
            tmp$se.2 <- ifelse(class(tmp_all)[1] == "tbl_df", sqrt(tmp_all[, paste0("se.",k)]^2 + tmp_all[, paste0("se.",j)]^2 - 2 * tmp_all[, "V"]) %>% pull, sqrt(tmp_all[, paste0("se.",k)]^2 + tmp_all[, paste0("se.",j)]^2 - 2 * tmp_all[, "V"]))
            tmp$Arm.1 <- ifelse(class(tmp_all)[1] == "tbl_df",tmp_all[, paste0("Arm.",j)] %>% pull,tmp_all[, paste0("Arm.",j)])
            tmp$Arm.2 <- ifelse(class(tmp_all)[1] == "tbl_df",tmp_all[, paste0("Arm.",k)] %>% pull,tmp_all[, paste0("Arm.",k)])
            tmp$Arm1 <- ifelse(class(tmp_all)[1] == "tbl_df", tmp_all[, paste0("Arm",j)] %>% pull,tmp_all[, paste0("Arm",j)])
            tmp$Arm2 <-  ifelse(class(tmp_all)[1] == "tbl_df", tmp_all[, paste0("Arm",k)] %>% pull, tmp_all[, paste0("Arm",k)])
          }else{
            tmp <- tmp_all[, c(paste0("lor.",k), paste0("se.",k), paste0("Arm.", j), 
                               paste0("Arm.", k), paste0("Arm", j), paste0("Arm", k))]
            names(tmp) <- c("lor.2","se.2","Arm.1","Arm.2","Arm1","Arm2")
          }
          ans <- rbind(ans, tmp)
        }
      }
    }
  }
  ans$rms <- rep(2,nrow(ans))
  ans$Comparison <- paste(ans$Arm1,"vs",ans$Arm2,sep = "")
  if(dataType == "Arm"){cols <- c(1:4, 7:8)}
  if(dataType == "Contrast"){cols <- c(1:2, 5:6)}
  ans[,cols] <- apply(ans[,cols], 2, function(x) as.numeric(as.character(x)))
  ans <- ans[order(ans$Arm1), ]
  ans <- as.data.frame(ans)
  if(save_res == T){
    write.csv(ans, file = "./data/pairwise_comparison.csv",row.names = FALSE)
  }
  return(ans)
}

