### generate a list containing two dataframes. One is the probability matrix that a treatment is better than
### another based on posterior samples of treatment risk (probability). The other one summarize the probability that each
### treatment can be the best or the worst.
### Input:
###   Jags output (not summary)
###   map_txt
###   good_event: Indicator, 1 if the events are good, 0 if the events are bad.
###   arm_name: select the treatment name type to output, available options are “abbr” or “full”
###   save_res: logical, whether to save the output into your local directory
### Return:
###   A list with two dataframes.One is the probability matrix that a treatment is better than
###   another based on posterior samples of treatment risk (probability). The other one summarize the probability that each
###   treatment can be the best or the worst.
###   Will save two csv files if save_res = T

prob_rank_table <- function(jags.out, map_txt, good_event = 0, arm_name = "abbr", save_res = T){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  db <- jags.out
  T_idx <- grep("T\\[", colnames(db[[1]]))
  T_db <- NULL
  for(i in 1:length(db)){
    T_db <- rbind(T_db,db[[i]][,T_idx])
  }
  NT <- length(T_idx)
  T_matrix <-  matrix(0,nrow = NT, ncol = NT)
  if(arm_name == "abbr"){
  treatNmsAbrev <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(abbr) %>% pull
  }
  if(arm_name == "full"){
    treatNmsAbrev <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(name) %>% pull
  }
  for(i in 1:(NT)){
    for(j in i:NT){
      if(j > i){
        T_matrix[i,j] <- as.character(round(mean(T_db[,i]<T_db[,j]),3))
        T_matrix[j,i] <- as.character(round(1 - as.numeric(T_matrix[i,j]),3))
      }else if(j == i){
        T_matrix[i,j] <- treatNmsAbrev[i]
      }
    }
  }
  if(good_event == 1){
    T_matrix <- t(T_matrix)
  }
  rk_idx <- grep("rk\\[", colnames(db[[1]]))
  rk_db <- NULL
  for(i in 1:length(db)){
    rk_db <- rbind(rk_db,db[[i]][,rk_idx])
  }
  
  rk_smry <- data.frame(Name = character(NT), Rank_1_prob = numeric(NT), Rank_least_prob = numeric(NT))
  rk_smry$Name <- treatNmsAbrev
  for(i in 1:NT){
    rk_smry$Rank_1_prob[i] <- round(mean(rk_db[,i]==1),3)
    rk_smry$Rank_least_prob[i] <- round(mean(rk_db[,i]==NT),3)
  }
  if(save_res == T){
  write.table(T_matrix, file = "./data/prob_comp.csv",
              quote = FALSE,
              sep = ",",
              row.names = FALSE,col.names = FALSE)
  
  write.table(rk_smry, file = "./data/rank_smry.csv",
              quote = FALSE,
              sep = ",",
              row.names = FALSE,col.names = TRUE)
  }
  return(list(prob_comp = T_matrix, rank_smry = rk_smry))
}

