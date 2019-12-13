### generate mean ranking matrix
### Input: 
###   map_txt:a dataframe which maps each treatment to a number
###   model_smry: general model summary generated from general_model(). direct comparison generated 
###   arm_name: select the treatment name type to output, available options are “abbr” or “full”
###   save_res: logical, whether to save the output into your local directory
### Return: 
###   mean ranking matrix, will save as a .csv file if save_res = T 

mean_rank_matrix <- function(map_txt, model_smry, arm_name = "abbr", save_res = T){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  db <- model_smry
  if(arm_name == "full"){
  treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(name)
  }
  if(arm_name == "abbr"){
    treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(abbr)
  }
  NT <- nrow(treatNms)
  ord <- seq(1,NT)
  rk_idx <- grep("rk\\[", rownames(db))
  rk_db  <- as.data.frame(db[rk_idx,])
  rk_db <- rk_db[order(-rk_db$Mean),]
  tmp_trt <- gsub("rk\\[","", rownames(rk_db))
  tmp_trt <- gsub("\\]$", "", tmp_trt)
  rk_db <- rk_db %>% mutate(treatOder = as.numeric(tmp_trt))
  rk_db <- rk_db %>% mutate(trtNm = sapply(treatOder, function(s){as.character(treatNms[s,1])}))
  rk_db <- rk_db[,c(7,1,2,3,4,5)]
  rk_db$Mean <- round(rk_db$Mean, 2)
  rk_db$SD <- round(rk_db$SD, 2)
  if(save_res == T){
    write.csv(rk_db, file = "./data/mean_ranking_matrix.csv", quote = FALSE,
              row.names = FALSE)
  }
  return(rk_db)
}

