## plot the network. You can set the node text size and node size and add plot title 
## Input:
##    MTCdata: MTCdata file
##    map_txt: a dataframe which maps each treatment to a number
##    dataType: the type of your MTCdata, available options are “Arm” and “Contrast”
##    for other arguments, please see mtm.networkplot.fun()
## Return:
##    the network plot

network_plot <- function(MTCdata, map_txt,  dataType = "Arm", 
                         percomparison = T, nodetextsize=.8, nodesize =0.8,  graphtitle='', 
                         vertex.col = "red",...){
  if( !(dataType %in% c("Arm", "Contrast") ) ){
    return(cat("Please input the correct data type, available options are \"Arm\" and \"Contrast\""))
  }
  pairwise_comp <- generateMTMNetWork(MTCdata, dataType = dataType)
  networkNode <-  mtmNetWorkNodeName(MTCdata, map_txt = map_txt)
  treat_new <- networkNode$treat_new
  if(dataType=="Arm"){
    total <- apply(cbind(pairwise_comp$n1, pairwise_comp$n2), 1, sum)
    mtm.networkplot.fun(pairwise_comp$t1, pairwise_comp$t2, percomparison = percomparison, 
                        nameoftreatments = treat_new, 
                        VAR2 = total, nodetextsize=nodetextsize, nodesize = nodesize, 
                        graphtitle = graphtitle, vertex.col = vertex.col)
  }
  if(dataType=="Contrast"){
    mtm.networkplot.fun(pairwise_comp$t1, pairwise_comp$t2, percomparison = percomparison, 
                        nameoftreatments = treat_new, 
                        nodetextsize=nodetextsize, nodesize = nodesize, 
                        graphtitle = graphtitle, vertex.col = vertex.col)
  }
}
