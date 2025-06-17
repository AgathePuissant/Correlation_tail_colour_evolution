#---------------------Load libraries--------------------------------------
source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")

hush=function(code){
  sink("NUL")
  tmp = code
  sink()
  return(tmp)
}

library(ggplot2)
library(plyr)
library(motmot)
library(geomorph)
library(progressr)
library(future)
library(future.apply)
library(parallel)
library(motmot)
library(geomorph)
library(readxl)
#--------------Load data---------------------

listsex=c("M","F")
view=c("D","V")
nloop = 1000

Digitalisation <- read_excel("./data/Digitalisation.xlsx")
Digitalisation$tipsgenre = paste0(Digitalisation$Genre,"_",Digitalisation$Espèce)

# Set up future for parallel execution
plan(multisession, workers = 4)

# Progress handler
handlers(global = TRUE)

# Data storage for results
results <- list()

#-------------Intra-specific re-sampling------------------------
# Parallelized loop with progressr
with_progress({
  p <- progressr::progressor(steps = nloop * 4)
  
  results <- future_lapply(seq_len(nloop), function(i) {
    result_list <- list()
    compteur <- 1
    
    for (j in seq_along(listsex)) {
      for (k in seq_along(view)) {
        
        list_get_phenotype <- hush(get_phenotype_plusshape(
          sex = listsex[j],
          view = view[k], 
          mode = 'random', 
          level = "sp", 
          path_coords = "./data/pca_embeddings_wotailsDV.csv", 
          path_data_photos = "./data/data_photos.csv", 
          fileshape = "./data/pcaH_onlytails.RData"))
        
        meanphen <- list_get_phenotype[[1]]
        data <- list_get_phenotype[[2]]
        meanphen_shape <- list_get_phenotype[[9]]
        lamb = list_get_phenotype[[10]]
        lam = lamb[1]
        select_genre <- unique(data$Genre)
        
        if (sum(length(unique(data[data$Genre %in% select_genre, ]$genresp))) > 2) {
          list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly = F, 
                                   tree_path = "./data/Papilionidae_MCC_clean.tre", genre = select_genre)
          
          saved_subtree <- list_match[[1]]
          meanphen <- list_match[[2]]
          meanphen_shape <- meanphen_shape[match(rownames(meanphen), rownames(meanphen_shape)), ]
          meanphen <- meanphen[rownames(meanphen) != "Ornithoptera_meridionalis" & rownames(meanphen) != "Ornithoptera_paradisea", ]
          meanphen_shape <- meanphen_shape[rownames(meanphen), ]
          
          subtree <- hush(match.phylo.data(saved_subtree, meanphen)$phy)
          subtree <- transformPhylo(subtree, model = "lambda", lambda = lam)
          
          IT <- phylo.integration(print.progress = F, A = as.matrix(meanphen), A2 = as.matrix(meanphen_shape), phy = subtree)
          
          if (IT$P.value < 0.05) {
            result_list[[compteur]] <- data.frame(R = as.numeric(IT$r.pls), sex = listsex[j], view = view[k], tailedonly = F)
            compteur <- compteur + 1
          }
          
          # Tailed only condition
          Tail_sub2 <- Digitalisation[Digitalisation$SEXE == listsex[j], c(2,13, 14), , F]
          Tail_sub2 <- as.data.frame(Tail_sub2)
          colnames(Tail_sub2)<-c("id","Tail","tipsgenre")
          Tail_sub2 = merge(Tail_sub2, data, by=c("id"))
          Tail_sub2 = Tail_sub2[,c(1,2,3)]
          Tail_sub2 = as.data.frame(Tail_sub2 %>% group_by(tipsgenre.x) %>% dplyr::summarise(Tail = unique(Tail)))
          rownames(Tail_sub2) <- Tail_sub2$tipsgenre.x
          Tail_sub2 <- Tail_sub2[(Tail_sub2$Tail == "T"),,F]
          lam = lamb[2]
          meanphen <- meanphen[rownames(meanphen) %in% rownames(Tail_sub2), ]
          meanphen_shape <- meanphen_shape[rownames(meanphen_shape) %in% rownames(Tail_sub2), ]
          
          
          if (dim(meanphen)[1] > 1) {
            subtree <- hush(match.phylo.data(saved_subtree, meanphen)$phy)
            subtree <- transformPhylo(subtree, model = "lambda", lambda = lam)
            
            IT2 <- phylo.integration(print.progress = F, A = as.matrix(meanphen), A2 = as.matrix(meanphen_shape), phy = subtree)
            
            if (IT2$P.value < 0.05) {
              result_list[[compteur]] <- data.frame(R = as.numeric(IT2$r.pls), sex = listsex[j], view = view[k], tailedonly = T)
              compteur <- compteur + 1
            }
          }
        }
        
        p()  # Update progress
      }
    }
    return(do.call(rbind, result_list))
  })
})

# Combine the results into a single dataframe
tableau_stockage <- do.call(rbind, results)

saveRDS(tableau_stockage,"./results/pPLS_table.RData")
#--------------Plot results-------------------------------

subset_tab = tableau_stockage
subset_tab = subset_tab[complete.cases(subset_tab),]
p_meds <- ddply(subset_tab, .(sex,tailedonly,view), summarise, med = round(median(R),2))

ggplot(data=subset_tab, aes(x=tailedonly, y = R, fill=sex))+
  geom_violin(alpha=0.2)+ 
  geom_boxplot(width=0.1)+ 
  facet_grid(sex~view)+
  theme_classic()+ 
  geom_text(data = p_meds, aes(x = tailedonly, y = med, label = med, color=sex), 
            size = 5, vjust = -2)+
  ylim(c(0,1))



tableau_stockage %>%
  group_by(sex,view,tailedonly) %>%
  dplyr::summarize(Mean = mean(R, na.rm=TRUE))

