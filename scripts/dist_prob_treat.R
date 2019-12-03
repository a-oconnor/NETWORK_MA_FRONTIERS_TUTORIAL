### generate plot for the distribution of probability of retreatmeng for each treatment in the MTC  
### Input:
###   jags.out: jags output (the posterior samples, not the summary)
###   map_txt: a dataframe which maps each treatment to a number
###   arm_name: select the treatment name type to output, available options are “abbr” or “full”
###   layout: the layout of your figure, i.e c(1,2), c(2,2)
###   num_treat_per_plot: the number of treatment shown in each plot
###   treat_interested: the group of treatments for which you want to plot probability distribution. The default is "all", you can also input a vector of treatment numbers to specify the treatments you want to plot.
###   legend_pos: the position (coordinates) where you want to add your legend
###   other ploting parameters: i.e. ylim, main, ylab, xlab, etc
### Return:
###  A plot for the distribution of probability of retreatmeng for each treatment in the MTC 

dist_prob_treat <- function(jags.out, map_txt, arm_name = "full", layout, num_treat_per_plot, treat_interested = "all",ylim = c(0,5),main="", ylab = "Density", xlab = "Probability of the event", legend_pos = c(0.5,4.5),cex = 0.8, cex.axis = 1){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  if(num_treat_per_plot >= 8){
    return(cat("The number of treatment per plot is too large, please pick a smaller one"))
  }
  db <- NULL
  for(i in length(jags.out)){
    db <- rbind(db, jags.out[[i]])
  }
  ## trt Idx 
  trt_idx <- grep("T\\[", colnames(db))
  ## extract treatment database
  db_trt <- db[,trt_idx]
  if(treat_interested == "all"){
  d <- seq(1,length(trt_idx))
  }else{
    d <- treat_interested
  }
  for(i in d){
    assign(paste('d', i, sep=''), density(as.vector(db_trt[, i]), from = 0, to =1))
  }
  if(layout[1]*layout[2]*(num_treat_per_plot-1) < (length(d) - 1)){
    return(cat("Some treatments cannot show up, please change the settings, either \"layout\" or \"num_treat_per_plot\" "))
  }
  s = seq(0,500,25)
  par(mfrow=layout,mar=c(4, 4, 1, 1) + 0.1)
  pch_map <- c(24, 22, 19, 3, 12, 10, 15)
  d_tmp <- c(d[-1], rep(0, layout[1]*layout[2] * (num_treat_per_plot - 1) - length(d[-1])))
  treat_mat <- matrix(d_tmp, nrow = layout[1]*layout[2], ncol = num_treat_per_plot - 1, byrow = T)
  for(i in 1:(layout[1]*layout[2])){
      d_onerow <- treat_mat[i,]
      if(d_onerow[1] == 0){break}
      plot(d1, xlim = c(0,1), ylim = ylim, main = main, ylab = ylab, xlab = xlab, cex = cex, type = "l",
           cex.axis = cex.axis)
      points(d1$x[s],d1$y[s], pch=pch_map[1])
      for(j in d_onerow){
        if(j == 0){break}
        lines(get(paste0("d",j)), main = "", xlab = "", ylab = "")
        points(get(paste0("d",j))$x[s], get(paste0("d",j))$y[s], pch = pch_map[which(treat_mat[i,] == j) + 1] )
      }
      if(arm_name == "full"){
      legend(legend_pos[1], legend_pos[2], 
             legend = map_txt$name[order(map_txt$number)][c(1,d_onerow[d_onerow!=0])], 
             pch = pch_map[c(1:(length(d_onerow[d_onerow!=0])+1) )],bty='n',cex=cex)
      }
      if(arm_name == "abbr"){
        legend(legend_pos[1], legend_pos[2], 
               legend = map_txt$abbr[order(map_txt$number)][c(1,d_onerow[d_onerow!=0])], 
               pch = pch_map[c(1:(length(d_onerow[d_onerow!=0])+1) )],bty='n',cex=cex)
      }
      if(any(d_onerow == 0)){break}
  }
  par(mfrow = c(1,1))
}





