## generate results of MTC comparison of effect measure for all possible treatments
## Input:
##    map_txt: a dataframe which maps each treatment to a number
##    model_smry: general model summary generated from general_model(). direct comparison generated 
##    measure: the effect measure users want to compare, options are "RR" for risk ratio, "OR" for odds ratio, "LOR" for log odds ratio.
##    stat: the statistics users want to use for comparison, options are "mean" and "median"
##    num_digits: the number of decimal places in the output 
##    save_res: logical, whether to save the output into your local directory
## Return:
##    a matrix with upper triangle to be mean or median, lower triangle to be 95% credible interval
##    If save_res = T, an csv file: MTC_comparison_risk_ratio_all_treat.csv will be save in ./data folder

MTC_comparison <- function(map_txt, model_smry, measure = "RR", stat = "mean", num_digits = 3, save_res = T){
  if(!(measure %in% c("RR", "OR", "LOR"))){
    return(cat("Please input a correct measure, available options are \"RR\", \"OR\" and \"LOR\"", quote = F))
  }
  if(!(stat %in% c("mean", "median"))){
    return(cat("Please input a correct stat, available options are \"mean\" and \"median\" ", quote = F))
  }
  treatNMsAbrev <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(abbr) %>% unique()
  NT <- nrow(treatNMsAbrev)
  ord <- seq(1,NT)
  treatOrder <- data.frame(ord = ord, abrev = treatNMsAbrev, stringsAsFactors = FALSE)
  if(measure == "RR"){idx <- grep("^RR\\[", rownames(model_smry))}
  if(measure == "OR"){idx <- grep("^or\\[", rownames(model_smry))}
  if(measure == "LOR"){idx <- grep("^lor\\[", rownames(model_smry))}
  db <- model_smry[idx, ]
  comp_matrix <- matrix(nrow = NT, ncol = NT)
  name_db <- gsub("[a-zA-Z]", "", rownames(db))
  for(i in 1:(NT)){
    for(j in i:NT){
      if(j == i){
        comp_matrix[i,j] <- as.character(treatNMsAbrev[i,1])
      }else if(j > i){
        compNm <- paste("[",j,",",i,"]", sep = "")
        current_idx <- name_db == compNm
        if(stat == "mean"){comp_stat <- round(db[current_idx, 1],num_digits)}
        if(stat == "median"){comp_stat <- round(db[current_idx, 4],num_digits)}
        comp_matrix[i,j] <-  comp_stat
        comp_lb <- round(db[current_idx,3], num_digits) ## lower bound
        comp_ub <- round(db[current_idx,5], num_digits) ## upper bound 
        comp_matrix[j,i] <- paste("(", comp_lb, "_", comp_ub,")", sep = "")
      }
    }
  }
  if(save_res == T){
    write.table(comp_matrix, file = paste0("./data/MTC_comparison_all_", measure, "_", stat,".csv"),
                quote = FALSE,
                sep = ",",
                row.names = FALSE,col.names = FALSE)
  }
  return(comp_matrix)
}


