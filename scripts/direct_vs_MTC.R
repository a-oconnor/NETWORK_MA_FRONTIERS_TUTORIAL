### Generate table: direct comparison, MTC, and only indirect comparison 
### Input: 
###   map_txt: a dataframe which maps each treatment to a number
###   model_smry: general model summary generated from general_model(). direct comparison generated 
###   direct_compr: direct comparison generated from direct_comparison()
###   arm_name: select the treatment name type to output, available options are “abbr” or “full”
###   save_res: logical, whether to save the output into your local directory
### Return: 
###   a dataframe containing the comparison of direct comparison, MTC, and only indirect comparison 
### Also save csv file ResultOfIndicretAnddirectComparison.csv if save_res is true

direct_vs_MTC <- function(map_txt, model_smry, direct_compr, arm_name = "abbr", save_res = T){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  d_MTC_idx <- grep("d\\[", rownames(model_smry))
  d_MTC_idx <- model_smry[d_MTC_idx, c("Mean", "SD")]
  d_MTC_idx <- as.data.frame(d_MTC_idx)
  d_MTC_idx$Nm <- rownames(d_MTC_idx)
  if(arm_name == "abbr"){
  treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(abbr) %>% pull
  }
  if(arm_name == "full"){
    treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(name) %>% pull
  }
  ord <- as.character(1:length(treatNms))
  treatOrder <- data.frame(treatNms = treatNms, ord = ord, stringsAsFactors = FALSE)
  convertToTreatNms <- function(s){
    ## input string in direct comparison
    ## convert into treatname, 
    s_split <- strsplit(s, "vs")[[1]]
    pre <- s_split[1]
    after <- s_split[2]
    preNm <- treatOrder %>% filter(ord == pre) %>%.[,"treatNms"]
    afterNm <- treatOrder %>% filter(ord == after) %>%.[,"treatNms"]
    nm <- paste(preNm, "vs", afterNm, sep = " ")
  }
  calMTCCompr <- function(s){
    ## calculate MTC comparison mean and associated sd, o1 vs active -- > d[active], otherwise 
    s_split <- strsplit(s, "vs")[[1]]
    pre <- s_split[1]
    aft <- s_split[2]
    if(pre == "01"){
      nm <- ifelse(grepl("0",aft),gsub("^0","",aft),aft) ## match same format as output from general model: 
      idxInMTC <- paste("d[",nm,"]", sep = "")
      d_MTC_mean <- d_MTC_idx %>%filter(Nm == idxInMTC) %>% .[, "Mean"] ## extract mean 
      d_MTC_sd   <- d_MTC_idx %>% filter(Nm == idxInMTC) %>% .[, "SD"] ## extract sd 
      return(list(d_MTC = d_MTC_mean, sd_MTC = d_MTC_sd))
    }else{
      preNm <- ifelse(grepl("0",pre),gsub("^0","",pre),pre)
      aftNm <- ifelse(grepl("0",aft),gsub("^0","",aft),aft)
      preInMTC <- paste("d[",preNm,"]", sep = "")
      aftInMTC <- paste("d[",aftNm,"]", sep = "")
      d_MTC_mean <- d_MTC_idx %>% filter(Nm == aftInMTC) %>% .[, "Mean"] - 
        d_MTC_idx %>% filter(Nm == preInMTC) %>% .[, "Mean"]
      sd_pre <- d_MTC_idx %>% filter(Nm == preInMTC) %>% .[, "SD"]
      sd_aft <- d_MTC_idx %>% filter(Nm == aftInMTC) %>% .[, "SD"]
      d_MTC_sd <- sqrt((sd_pre)^2 + (sd_aft)^2 -  sd_pre * sd_aft)
      return(list(d_MTC = d_MTC_mean, sd_MTC = d_MTC_sd))
    }
  }
  comr_str <- as.character(direct_compr$compar)
  ## string to treat names 
  comparInTreatNm <- unlist(lapply(comr_str, convertToTreatNms))
  ## get MTC mean and MTC sd
  d_MTC <- NULL
  sd_MTC <- NULL
  for(i in 1:length(comr_str)){
    tmp <- calMTCCompr(comr_str[i])
    d_MTC[i] <- tmp$d_MTC
    sd_MTC[i] <- tmp$sd_MTC
  }
  
  df <- data.frame(Comparison = comparInTreatNm, d_dir = direct_compr$Mean, sd_dir = direct_compr$SD,
                   d_MTC = d_MTC, sd_MTC = sd_MTC)
  
  ## only indirect comparison 
  df <- df %>% mutate(d_rest = (d_MTC/((sd_MTC)^2) - d_dir/((sd_dir)^2))/(1.0/((sd_MTC)^2) - 1.0/((sd_dir)^2))
  )
  df <- df %>% mutate(sd_rest = sqrt(1.0/(1.0/(sd_MTC^2) - 1.0/(sd_dir^2))))
  df <- df %>% mutate(w = d_dir - d_rest)
  
  df <- df %>% mutate(sd_w = sqrt(sd_dir^2 + sd_rest^2))
  
  df <- df %>% mutate(z_score = w/sd_w)
  df <- df %>% mutate(p_value = 2 * pnorm(-abs(z_score)))
  
  ## write result into csv file 
  df[,2:11] <- round(df[,2:11],2)
  if(save_res == T){
  write.csv(df, file = "./data/ResultOfIndicretAnddirectComparison.csv", quote = FALSE,
            row.names = FALSE)
  }
  return(df)
}

