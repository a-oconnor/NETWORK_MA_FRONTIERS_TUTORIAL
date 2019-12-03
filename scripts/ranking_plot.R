### generate the ranking plot
### Input:
###   MTCdata: MTCdata file
###   map_txt: a dataframe which maps each treatment to a number
###   model_smry: general model summary generated from general_model(). direct comparison generated 
###   arm_name: select the treatment name type to output, available options are “abbr” or “full”
###   cex: the font size
### Return:
###   A ranking plot

ranking_plot <- function(MTCdata, map_txt, model_smry, arm_name = "abbr", cex = NULL){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  networkNode <-  mtmNetWorkNodeName(MTCdata = MTCdata, map_txt = map_txt, arm_name = arm_name)
  
  idx <- grep("rk\\[", rownames(model_smry))
  nt <- length(idx)
  m.r <- model_smry[idx,"Mean"]
  sd.r <- model_smry[idx,"SD"]
  prec.r <- sd.r^(-2)
  nam.r <- names(model_smry[idx,"Mean"])
  cilb <- model_smry[idx,"2.5%"]
  ciub <- model_smry[idx,"97.5%"]
  
  ## treatment name 
  trea_n <- networkNode$treat_new
  
  
  pf <- data.frame(m.r,sd.r,prec.r,cilb,ciub,trea_n)
  #Order the previous data frame to made the forest plot in order
  pf.o <- pf[order(pf[,1]),]
  if(nt < 5){
    at_seg <- 1:nt
  }else{
    at_seg <- round(quantile(1:nt),0)
  }
  forest.nds(pf.o[,1],pf.o[,2],main='',refline=NA,efac=1,
             alim=c(0,(nt+1)), at = at_seg,
             steps=10, xlim=c(0,(nt+1)),
             ci.lb=as.numeric(pf.o[,4]),ci.ub=as.numeric(pf.o[,5]),
             xlab='Ranking',
             slab=pf.o[,6],digits = 2, cex = cex)
}
