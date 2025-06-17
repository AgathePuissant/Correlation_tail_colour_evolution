sex="M"
view="D"

#---------------------Load libraries--------------------
source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")
source("./code/basis_functions/f_RR.R")

library(geomorph)
library(caper)
library(RRphylo)
library(dispRity)
library(stringr)
library(picante)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)
library(viridis)
library(RColorBrewer)

#-------------RRphylo loop (very long to run)-------------------
tour=1

R2 = rep(NA, tour)

Rates_colour = list()
Rates_shape = list()

for (t in c(1:tour)){
  print(t)
  
  #-----------------------Load and clean data--------------------
  
  list_get_phenotype = get_phenotype_plusshape(sex,
                                               view, 
                                               mode = 'random', 
                                               level = "sp", 
                                               path_coords = "./data/pca_embeddings_wotailsDV.csv", 
                                               path_data_photos = "./data/data_photos.csv", 
                                               fileshape="./data/pcaH_onlytails.RData")
  meanphen <- list_get_phenotype[[1]]
  data <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  meanphen_shape <- list_get_phenotype[[9]]
  
  rm(list=c("list_get_phenotype"))
  
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly=F, tree_path = "./data/Papilionidae_MCC_clean.tre")
  
  subtree <- list_match[[1]]
  meanphen <- list_match[[2]]
  
  subtree_ori = subtree
  subtree_shape = match.phylo.data(subtree, meanphen_shape)$phy
  subtree_shape = match.phylo.data(subtree_shape, meanphen)$phy
  
  subtree = match.phylo.data(subtree_shape, meanphen)$phy

  meanphen_shape = match.phylo.data(subtree, meanphen_shape)$data
  meanphen = match.phylo.data(subtree_shape, meanphen)$data
  
  #-----------------------RRphylo estimates--------------------
  
  RR_shape <- RRphylo(tree = subtree_shape, y=meanphen_shape)
  RR_colour <- RRphylo(tree = subtree_shape, y=meanphen)
  
  
  estimates_shape <- c(RR_shape$aces[,1],meanphen_shape[,1])
  names(estimates_shape) <- c(names(RR_shape$aces[,1]),rownames(meanphen_shape))
  RR_shape_cov <- RRphylo(tree=subtree_shape,y=meanphen_shape,cov=estimates_shape)
  
  estimates_colour <- c(RR_colour$aces[,1],meanphen[,1])
  names(estimates_colour) <- c(names(RR_colour$aces[,1]),rownames(meanphen))
  RR_colour_cov <- RRphylo(tree=subtree_shape,y=meanphen,cov=estimates_colour)
  
  #for sp level
  subtree_shape <- RR_shape_cov$tree
  
  
  fitrates <- lm(RR_shape_cov$rates~RR_colour_cov$rates)
  print(summary(fitrates))
  
  Rates_colour[[t]] <- RR_colour_cov$rates
  Rates_shape[[t]] <- RR_shape_cov$rates
  
  R2[t] <- summary(fitrates)$adj.r.squared
  
}


#-----------------Load pre-run results--------------
load("./results/RR_MD_100t.RData")
#-------------------Plot results-------------------------

ggplot()+geom_density(aes(x=R2), fill="lightblue", alpha = 0.5)+theme_classic()
summary(R2)


tomatchtree = rownames(Rates_shape[[1]])[str_detect(rownames(Rates_shape[[1]]),"_")]
names(tomatchtree) = tomatchtree
subtree = read.tree("./data/Papilionidae_MCC_clean.tre")
subtree_shape=match.phylo.data(subtree, tomatchtree)$phy

# Extract the edge information
edges_df <- as.data.frame(subtree_shape$edge)
edges_df <- cbind(edges_df,subtree_shape$edge.length)
colnames(edges_df) <- c("parent", "node","length")
# Get tip labels
tip_labels <- subtree_shape$tip.label
# Map node numbers to tip labels
node_to_tip <- data.frame(node = 1:length(tip_labels), tip_label = tip_labels)
# Merge with edges to associate nodes with tip labels
edges_df <- edges_df %>%
  left_join(node_to_tip, by = c("node" = "node"))
edges_df$newnode = edges_df$node
edges_df[!is.na(edges_df$tip_label),]$newnode = edges_df[!is.na(edges_df$tip_label),]$tip_label


df_rates_colour = as.data.frame(Rates_colour)
df_rates_colour = apply(df_rates_colour,1,FUN = mean)
df_rates_colour_log = log(1+df_rates_colour) #log scale for better visualization
df_rates_colour_log= as.data.frame(df_rates_colour_log)
df_rates_colour_log$edge_number = rownames(df_rates_colour_log)
colnames(df_rates_colour_log) <- c("colour_rate","newnode")
# Combine edge information with the evolution rates dataframe
df_rates_colour_log <- df_rates_colour_log %>%
  left_join(edges_df, by = "newnode")


df_rates_shape = as.data.frame(Rates_shape)
df_rates_shape = apply(df_rates_shape,1,FUN = mean)
df_rates_shape_log = log(1+df_rates_shape)
df_rates_shape_log= as.data.frame(df_rates_shape_log)
df_rates_shape_log$edge_number = rownames(df_rates_shape_log)
colnames(df_rates_shape_log) <- c("shape_rate","newnode")
# Combine edge information with the evolution rates dataframe
df_rates_shape_log <- df_rates_shape_log %>%
  left_join(edges_df, by = "newnode")


p1 <- ggtree(subtree_shape) %<+% df_rates_shape_log + 
  geom_tree(aes(color=shape_rate), size=1)+
  scale_color_viridis()+
  labs(color = "Rates of shape evolution" )


p2 <- ggtree(subtree_shape) %<+% df_rates_colour_log + 
  geom_tree(aes(color=colour_rate), size=1)


d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1+      
  ggnewscale::new_scale_color() + 
  geom_tree(data=d2,aes(color=colour_rate), size=1)+      
  scale_color_viridis(option="magma")+
  labs(color = "Rates of colour pattern evolution" )


pp



#Variation

df_rates_colour_var = as.data.frame(Rates_colour)
df_rates_colour_var = apply(df_rates_colour_var,1,FUN = sd)
df_rates_colour_var_log = log(1+df_rates_colour_var)
df_rates_colour_var_log= as.data.frame(df_rates_colour_var_log)
df_rates_colour_var_log$edge_number = rownames(df_rates_colour_var_log)
colnames(df_rates_colour_var_log) <- c("colour_rate","newnode")
# Combine edge information with the evolution rates dataframe
df_rates_colour_var_log <- df_rates_colour_var_log %>%
  left_join(edges_df, by = "newnode")

df_rates_shape_var = as.data.frame(Rates_shape)
df_rates_shape_var = apply(df_rates_shape_var,1,FUN = sd)
df_rates_shape_var_log = log(1+df_rates_shape_var) #df_rates_shape_var #
df_rates_shape_var_log= as.data.frame(df_rates_shape_var_log)
df_rates_shape_var_log$edge_number = rownames(df_rates_shape_var_log)

colnames(df_rates_shape_var_log) <- c("shape_rate","newnode")
# Combine edge information with the evolution rates dataframe
df_rates_shape_var_log <- df_rates_shape_var_log %>%
  left_join(edges_df, by = "newnode")


p1 <- ggtree(subtree_shape) %<+% df_rates_shape_var_log + 
  geom_tree(aes(color=shape_rate), size=1)+
  scale_color_viridis(option="turbo")+
  labs(color = "Standard deviation of rates \n of shape evolution" )


p2 <- ggtree(subtree_shape) %<+% df_rates_colour_var_log + 
  geom_tree(aes(color=colour_rate), size=1)


d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1+      
  ggnewscale::new_scale_color() + 
  geom_tree(data=d2,aes(color=colour_rate), size=1)+
  scale_color_viridis(option="turbo")+
  labs(color = "Standard deviation of rates \n of colour pattern evolution" )

pp


nodename = df_rates_colour_var_log$parent
node_newname = df_rates_colour_var_log$newnode
edgelength = df_rates_colour_log$length



df_plot = data.frame(xshape = df_rates_shape, 
                     ycolour = df_rates_colour,
                     xshape_var = df_rates_shape_var,
                     ycolour_var = df_rates_colour_var,
                     nodename=nodename,
                     node_newnames = node_newname,
                     length=edgelength)
df_plot$node_newname[str_detect(df_plot$node_newname,"[0-9]")]=""
df_plot$genus = gsub("_[a-z]*$","",df_plot$node_newname)


p1 <- ggplot(df_plot,aes(x=xshape,y=ycolour, label=nodename))+
  theme_classic()+
  geom_errorbar(aes(ymin=ycolour-ycolour_var,
                    ymax=ycolour+ycolour_var),
                width=.3,
                alpha=0.2)+
  geom_errorbarh(aes(xmin=xshape-xshape_var,
                     xmax=xshape+xshape_var),
                 height=.3,
                 alpha=0.2)+
  geom_point(fill = "grey",colour="black", size=3, shape=21)+ #`Node age`+
  ylab("Rates of colour pattern evolution")+
  xlab("Rates of tail shape evolution")
p1

