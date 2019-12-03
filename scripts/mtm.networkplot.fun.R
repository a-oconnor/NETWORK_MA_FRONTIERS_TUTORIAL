## reference: 
## code was download from http://www.mtm.uoi.gr/index.php/how-to-do-an-mtm/10-how-to-do-an-mtm/15-graphicaldescriptionofanetwork
#####################################################
#########FUNCTION TO RUN THE NETWORK PLOT############
#####################################################
mtm.networkplot.fun <- function(t1,t2,percomparison,nameoftreatments,VAR1,VAR2,graphtitle,thickness,nodetextsize,nodesize,
                                vertex.col = "red",...) {
  # The purpose of this function is the multiple-treatments network representation 
  # in a graphical manner using nodes and edges. Every node will have different
  # size and every edge will have different thickness with respect to the totalrandomized 
  # sample size of comparisons. 
  
  # IMPORTANT NOTE: In order to use this function  "network" library is required
  # WARNING: When the user supply as the first two input arguments the two treatments
  # then it is mandatory to set the third Boolean argument as T. 
  # In the alternative case, where the first three input arguments corresponds to the ID
  # of the treatment, the treatments and the third argument should be put as F.
  
  # t1 refers to the first treatment where,
  # t2 refers to the second treatment and
  #______________________________________________________________________
  
  options(warn=-1)
  if (missing(t1)){stop("Need to specify the first argument.")}
  
  if (missing(t2)){stop("Need to specify the second argument.")}
  
  #_______________________________________________________________________
  library(network)# this function needs package "network" in order to work
  #_______________________________________________________________________
  
  
  if (percomparison==T){#if the first two arguments are the treatments and the 
    #second two arguments are the corresponding randomized sample sizes
    t1<-t1
    t2<-t2
  }
  #in this part we transform the two arguments: ID and Treatment
  #in the same as previous form
  else{
    ID<-t1
    Treat<-t2    
    treat <- as.double(as.factor(Treat))
    id <- ID
    nt <- length(unique(treat))
    ns <- length(unique(ID))
    na <- table(match(id, unique(id)))
    #_____Check whether the treatments are 1 to nt
    if(max(sort(unique(treat)) - c(1:nt)) > 0) {
      treat <- as.factor(treat)
      levels(treat) <- c(1:nt)
      treat <- as.double(treat)
    }
    TT <- matrix(NA, nrow = ns, ncol = max(na))
    
    for(i in 1:ns) {
      TT[i, 1:na[i]] <- treat[id == unique(id)[i]]
    }
    u <- TT
    new.id <- 1:length(table(id))  
    TT <- (apply(TT, 1, sort))
    na <- na
    u <- c()
    torepeat <- na[na > 2]
    # if there are more than 3 arms
    if(sum(na > 2))
    {
      for(i in 1:(sum(na > 2))) 
      {
        u <- rbind(u, t(combn(unlist(TT[na > 2][i]),2)))
      }
    }
    
    TT<- rbind(matrix(unlist(TT[na == 2]), ncol = 2, byrow = T), u)
    t1<-TT[,1]
    t2<-TT[,2]    
  }# end of the else that transforms (ID,Treatment) to (t1,t2)
  #_____________________________________________________________________
  
  numoftreatments=max(cbind(t1,t2));
  
  # Default adjacency-style matrix initialization
  mat_treat<-matrix(0,numoftreatments,numoftreatments);
  
  if (missing(nameoftreatments))
  {
    nam="";
    for(i in 1:numoftreatments) 
    { 
      nam[i] <- paste("treat",i,sep=".")
      assign(nam[i],1:i)
    }
    colnames(mat_treat)<-nam;
    rownames(mat_treat)<-nam;
  }
  else
  {
    colnames(mat_treat)<-nameoftreatments;
    rownames(mat_treat)<-nameoftreatments;
  }
  
  # Adjusted adjacency-style matrix initialization (for thickness of links)
  thickness_new<-matrix(0,numoftreatments,numoftreatments);
  
  if (missing(nameoftreatments))
  {
    nam="";
    for(i in 1:numoftreatments) 
    { 
      nam[i] <- paste("treat",i,sep=".")
      assign(nam[i],1:i)
    }
    colnames(thickness_new)<-nam;
    rownames(thickness_new)<-nam;
  }
  else
  {
    colnames(thickness_new)<-nameoftreatments;
    rownames(thickness_new)<-nameoftreatments;
  }
  #_____________________________________________________________________
  # based on the frequency of treatments we construct the node's thickness
  # divisor has fixed value equal to 5 for aesthetical reasons
  if (missing(nodesize)){
    divisor=5} else {
      divisor=5/nodesize
    }
  
  if (missing(VAR1)){
    nodethickness<-table(c(t1,t2))/divisor }else{
      nodethickness<-VAR1
    }
  
  if (missing(thickness)){
    thickness=10
  }
  
  #combined treatments (for default link thickness)
  tr<-cbind(t1,t2)
  
  for (i in 1:length(t1))
  {
    mat_treat[tr[i,1],tr[i,2]]<-mat_treat[tr[i,1],tr[i,2]]+1
  }
  
  mat_treatb<-t(mat_treat)
  mat_treat<-mat_treat+mat_treatb
  
  if(percomparison==T){
    if(missing(VAR2)){
      thickness=10} else {
        tr2<-tr #combined treatments (for adjusted link thickness) when dataset is formed by study
        n<-VAR2 
        for (i in 1:length(t1))
        {
          thickness_new[tr2[i,1],tr2[i,2]]<-thickness_new[tr2[i,1],tr2[i,2]]+n[i]
        }
        
        thickness_newb<-t(thickness_new)
        thickness_new<-thickness_new+thickness_newb
      }
  }
  
  if(percomparison==F){
    if(missing(VAR2)){
      thickness=10} else {
        tr2<-tr #combined treatments (for adjusted link thickness) when dataset is formed by arm
        Total<-VAR2
        N <- matrix(NA, nrow = ns, ncol = max(na))    
        for(i in 1:ns) {
          N[i, 1:na[i]] <- Total[id == unique(id)[i]]
        }
        v <- N	
        N <- (apply(N, 1, sort))
        v <- c()
        torepeat <- na[na > 2]
        # if there are more than 3 arms
        if(sum(na > 2))
        {
          for(i in 1:(sum(na > 2))) 
          {
            v <- rbind(v, t(combn(unlist(N[na > 2][i]),2)))
          }
        }
        N<- rbind(matrix(unlist(N[na == 2]), ncol = 2, byrow = T), v)  
        n1<-N[,1]
        n2<-N[,2]
        VAR2<-apply(cbind(n1,n2),1,sum)
        n<-VAR2 
        for (i in 1:length(t1))
        {
          thickness_new[tr2[i,1],tr2[i,2]]<-thickness_new[tr2[i,1],tr2[i,2]]+n[i]
        }
        
        thickness_newb<-t(thickness_new)
        thickness_new<-thickness_new+thickness_newb
      }
  }
  
  # Based on package "network" we construct the netdata
  if (missing(VAR2)){
    netdata<-network(mat_treat,directed=FALSE,cex=.5)
    mat_treat_new=mat_treat*thickness/max(mat_treat)} else { # network plotting (default for thickness of link)
      netdata<-network(thickness_new,directed=FALSE)
      mat_treat_new=thickness_new*thickness/max(thickness_new) # network plotting (adjusted for thickness of link)
    }
  
  if (missing(nodetextsize)){
    par(cex=0.8)}else
    {
      par(cex=nodetextsize)
    }
  
  if (missing(graphtitle)){
    par(mai=c(0,0,0,0))
    par(cex=0.8)
    pltitle<-""
  }else{
    par(mai=c(0.2,0.2,0.2,0.2))
    par(cex=0.8)
    pltitle=graphtitle
  }
  
  #plot of the network
  plot(netdata,mode="circle",displaylabels=TRUE,vertex.cex=nodethickness,boxed.label=FALSE,edge.lwd=mat_treat_new,cex=.4,
       vertex.col = vertex.col)
  
  #_____________________________________________________________________
  
  title(main = list(pltitle, col="blue", font=2))
}
