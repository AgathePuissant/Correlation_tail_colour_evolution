library(RRphylo)
library(ggtree)
library(RColorBrewer)


plot_shifts <- function(RR, subtree){
  search.shift(RR, status.type="clade",filename="ss automatic")->SSauto
  
  # ov<-overfitRR(RR=RR,y=meanphen_M, swap.args =list(si=0.2,si2=0.2),
            # shift.args = list(node=rownames(SSauto$single.clades)),
            # nsim=10)
  
  p<-ggtree(subtree,layout="circular")
  
  clades <-SSauto$single.clades
  
  posshifts = rbind(as.numeric(rownames(clades[clades$rate.difference>0,])),clades[clades$rate.difference>0,1])
  pshifts = data.frame(t(posshifts))
  colnames(pshifts) <- c("Node","Value")
  negshifts = rbind(as.numeric(rownames(clades[clades$rate.difference<0,])),(-1)*clades[clades$rate.difference<0,1])
  nshifts = data.frame(t(negshifts))
  colnames(nshifts)<-c("Node","Value")
  
  for (i in c(1:dim(pshifts)[1])){
    p <- p + geom_point2(aes_(subset=(p$data$node %in% pshifts[i,]$Node)),size=pshifts[i,]$Value*10, shape=21, fill='blue',alpha = 0.5)
  }
  for (i in c(1:dim(nshifts)[1])){
    p <- p + geom_point2(aes_(subset=(p$data$node %in% nshifts[i,]$Node)),size=nshifts[i,]$Value*10, shape=21, fill='red',alpha = 0.5) 
  }
    
  p <- p + scale_size(limits = c(0,0.1), range = c(1, 7))+
    theme(legend.position='none')
  p
  return(p)
}

overfitRR_shifts <- function(RR, subtree, meanphen){
  search.shift(RR, status.type="clade",filename="ss automatic")->SSauto
  search.shift(RR, status.type="clade",filename="ss automatic", node = rownames(SSauto$single.clades))->SS
  print(SS$single.clades)
  
  ov<-overfitRR(RR=RR,y=meanphen, swap.args =list(si=0.2,si2=0.2),
  shift.args = list(node=rownames(SS$single.clades)),
  nsim=10)
  print(ov)
  
  return(ov)
}

plot_rates <-function(RR,subtree, log_tr = F){
  rates <- as.data.frame(RR$rates)
  colnames(rates)<-c("V1")
  if (log_tr==T){
    rates$V1 = log(rates$V1+1+abs(min(rates$V1)))
    rates$V1[is.na(rates$V1)] <- min(rates$V1[is.finite(rates$V1)])
  }
  
  
  p <- ggtree(subtree, layout='circular') + 
    geom_tree(aes(color=rates$V1), continuous = 'colour', size=1) + 
    # scale_color_distiller(palette = "YlOrRd", trans = "reverse")+
    scale_color_gradientn(colours=c("lightblue", 'cyan', 'green', 'orange', 'red'),limits=c(min(rates$V1),max(rates$V1))) +
    theme(legend.position = c(.05, .85))+
    geom_tiplab()
  
  p
  return(p)
}

