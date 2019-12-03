## network geometry summary: PIE index, c_score, etc
## Input:
##   Datasets needed:
##   MTCdata: MTCdata file
##   Options:
##   save_res: logical, whether to save the output into your local directory
## Return:
##   PIE index
##   network geometry dataframe containing study and treatment
##   network co-occurance matrix
##   c_score


network_geometry <- function(MTCdata, save_res = T){
  maxArm <- max(MTCdata$Number.of.arms)
  network_df <- data.frame(study = rep(MTCdata$study, maxArm))
  network_df$treatment <- MTCdata[,paste0("Arm.",1:maxArm)] %>% as.matrix %>% as.vector()
  network_df <- network_df[!is.na(network_df$treatment),]
  PIE_index <- function(dat,Trtname){
    S <- nrow(dat)
    N <- length(Trtname)
    pVecs <- rep(0,N)
    for(i in 1:N){
      trt <- Trtname[i]
      pVecs[i] <- (sum(dat$treatment==trt)/S)^2
    }
    PIE <- (S/(S-1))*(1-sum(pVecs))
    PIE_max <- S*(N-1)/(N*(S-1))
    return(list(PIE=PIE,PIE_max=PIE_max))
  }
  Trtname <- unique(network_df$treatment)
  PIE <- PIE_index(network_df,Trtname = Trtname)
  network_df <- network_df[order(network_df$study),]
  network_cooc_mat <- matrix(0,nrow = length(unique(network_df$treatment)),ncol = length(unique(network_df$study)))
  for(i in 1:length(unique(network_df$study))){
    index <- match(network_df[network_df$study==unique(network_df$study)[i],"treatment"],sort(unique(network_df$treatment)))
    network_cooc_mat[index,i] <- 1
  }
  if(save_res == T){
  write.csv(network_df,'./data/Network_geometry.csv',row.names = F)
  write.csv(network_cooc_mat,file = "./data/network_cooc_mat.csv",row.names = F)
  }
  cooc <- cooc_null_model(network_cooc_mat,suppressProg = T)
  summary(cooc)
  return(list(PIE = PIE, network_df = network_df, network_cooc_mat = network_cooc_mat, 
              c_score = c_score(network_cooc_mat)))
}

