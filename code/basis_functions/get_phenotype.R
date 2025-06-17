library(dplyr)
library(StatMatch)
library(ade4)
library(stringr)
library(readr)
library(sampling)

harmonize <- function(path_data_photos= "data_photos.csv", path_coords = "pca_embeddings.csv", level="ssp", fileshape){
  data<-read.csv(file=path_coords, sep =';', header = T)
  data_loc<-read_delim(path_data_photos, delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
 
  oldnames <- colnames(data)
  
  data<-merge(x = data, y = data_loc[,c("id","Genre","SEXE","sp","loc","ssp","Subgenus","form")], by = "id") #,"lat","lon","altitude"
  
  colnames(data) <- c(oldnames,"Genre","SEXE","sp","loc","ssp","Subgenus","form") #,"lat","lon","altitude"
  
  data$ssp[data$ssp=='-'] <- data$sp[data$ssp=='-']
  data$Subgenus[data$Subgenus=='-'] <- data$Genre[data$Subgenus=='-']
  data$form[is.na(data$form)] <- "NA"#data$ssp[is.na(data$form)] #MODIF
  
  data$tipsgenre <- paste(data$Genre,data$sp,sep="_")
  if (level=="ssp"){
    data$genresp <- paste(data$Genre,data$sp,data$ssp,sep="_")
  }
  else if (level=="sp"){
    data$genresp <- paste(data$Genre,data$sp,sep="_")
  }
  else if (level=="form"){
    data$genresp <- paste(data$Genre,data$sp,data$form,sep="_")#MODIF #data$genresp <- paste(data$Genre,data$sp,data$ssp,data$form,sep="_")
  }
  
  data<-data %>% 
    mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
  
  data$sp<-as.factor(data$sp)
  data$Genre<-as.factor(data$Genre)
  data$tipsgenre <- as.factor(data$tipsgenre)
  data$genresp <- as.factor(data$genresp)
  data$ssp <- as.factor(data$ssp)
  data$loc<-as.factor(data$loc)
  data$Subgenus<-as.factor(data$Subgenus)
  
  genresp_keep <- table(data[,c("tipsgenre","SEXE","view")])
  
  
  
  data_shape = load(file=fileshape)
  
  data_shape = data.frame(pcaH$x)
  
  rownames(data_shape) = gsub(".nef","", rownames(data_shape))
  rownames(data_shape) = gsub("\\.","", rownames(data_shape))
  rownames(data_shape) = gsub("D|V","", rownames(data_shape))
  rownames(data_shape) = gsub("EL0","EL", rownames(data_shape))
  data_shape = cbind(rownames(data_shape), data_shape)
  colnames(data_shape) <- c("id", colnames(data_shape)[2:length(colnames(data_shape))])
  
  oldnames <- colnames(data_shape)
  data_shape<-merge(x = data_shape, y = data_loc[,c("id","Genre","SEXE","sp","loc","ssp","Subgenus","form")], by = "id")
  
  
  
  colnames(data_shape) <- c(oldnames,"Genre","SEXE","sp","loc","ssp","Subgenus","form") #,"lat","lon","altitude"
  
  data_shape$ssp[data_shape$ssp=='-'] <- data_shape$sp[data_shape$ssp=='-']
  data_shape$Subgenus[data_shape$Subgenus=='-'] <- data_shape$Genre[data_shape$Subgenus=='-']
  data_shape$form[is.na(data_shape$form)] <- "NA" #MODIF #data_shape$ssp[is.na(data_shape$form)]
  
  data_shape$tipsgenre <- paste(data_shape$Genre,data_shape$sp,sep="_")
  if (level=="ssp"){
    data_shape$genresp <- paste(data_shape$Genre,data_shape$sp,data_shape$ssp,sep="_")
  }
  else if (level=="sp"){
    data_shape$genresp <- data_shape$tipsgenre
  }
  else if (level=="form"){
    data_shape$genresp <- paste(data_shape$Genre,data_shape$sp,data_shape$form,sep="_")#MODIF #data_shape$genresp <- paste(data_shape$Genre,data_shape$sp,data_shape$ssp,data_shape$form,sep="_")
  }
  
  
  data_shape = data_shape[data_shape$id %in% data$id,]
  data_shape = data_shape[!is.na(rownames(data_shape)),]
  
  
  return(list(data,data_shape))
}

get_phenotype_plusshape <- function(sex, view, wing=0, side=0, level = "ssp", mode = 'mean', path_data_photos = "data_photos.csv", path_coords = "pca_embeddings.csv", indices = 0, mode_indices="id", fileshape = "./data/pcaH.RData"){
  
  data_all<-harmonize(path_data_photos = path_data_photos, path_coords = path_coords, level=level, fileshape = fileshape)
  data <- data_all[[1]]
  data_shape<-data_all[[2]]
  
  if(mode_indices=="ssp"){
    data <- data %>%                                        # Create ID by group
      group_by(Genre, sp, ssp) %>%
      dplyr::mutate(ID = cur_group_id())
  }
  else {
    data <- data %>%                                        # Create ID by group
      group_by(id) %>%
      dplyr::mutate(ID = cur_group_id())
  }
  
  data_sex <- data[ which(data$SEXE %in% sex & data$view %in% view), ]
  
  if (wing != 0){
    data_sex <- data_sex[ which(data_sex$wing %in% wing), ]
    }
  if (side != 0){
    data_sex <- data_sex[ which(data_sex$side %in% side), ]
    }
  
  data_sex = data_sex[order(data_sex$genresp),]
  
  data_shape = data_shape[data_shape$id %in% data_sex$id,]
  
  if (mode == 'mean'){
    meanphen <- aggregate(dplyr::select( data.frame(data_sex), contains('coord')), list(data_sex$genresp), FUN=mean)
    meanphen2 <- meanphen[,-1]
    rownames(meanphen2) <- meanphen[,1]
    meanphen <- meanphen2
    ids = NaN
    
    meanphen_shape <- aggregate(dplyr::select( data.frame(data_shape), contains('Comp')), list(data_shape$genresp), FUN=mean)
    meanphen_shape2 <- meanphen_shape[,-1]
    rownames(meanphen_shape2) <- meanphen_shape[,1]
    meanphen_shape <- meanphen_shape2
    ids = NaN
    
    meanphen_shape = meanphen_shape[rownames(meanphen),,drop=F]
  }
  else if (mode == 'random'){

    data_shape$genresp = as.factor(data_shape$genresp)
    x <- sampling::strata(data = data_shape, stratanames = "genresp", size = rep(1,nlevels(data_shape$genresp)),method="srswor")
    ids = data_shape[x$ID_unit,]$id
    meanphen_shape = data.frame(data_shape[x$ID_unit,])
    meanphen_shape = meanphen_shape[!is.na(meanphen_shape$genresp),]
    rownames(meanphen_shape) <- meanphen_shape$genresp
    data_shape = data_shape[x$ID_unit,]
    
    meanphen = data.frame(data_sex[match(ids,data_sex$id),])
    meanphen = meanphen[!is.na(meanphen$genresp),]
    rownames(meanphen) <- meanphen$genresp
    data_sex = data.frame(data_sex[match(ids,data_sex$id),])
    
    
    meanphen = meanphen[,grep("coord", colnames(meanphen))]
    meanphen_shape = meanphen_shape[,grep("Comp", colnames(meanphen_shape))]
    
  }
  else if (mode == 'indices'){
    
    data_sex2 = data.frame(data_sex[data_sex$ID %in% indices,])
    x <- sampling::strata(data = data_sex2, stratanames = "ID", size = rep(1,nlevels(data_sex2$genresp)),method="srswor")
    meanphen = data.frame(data_sex2[x$ID_unit,])
    rownames(meanphen) <- meanphen$genresp
    meanphen = dplyr::select(meanphen, contains('coord'))
    ids = NaN
  }
  
  listloc <- 0
  
  genre_data<-unique(data_sex$Genre)
  sp_data <- unique(data_sex$sp)
  tipsgenre_data <- unique(data_sex$tipsgenre)
  genresp_data <- unique(data_sex$genresp)
  
  lambdas = matrix(c(0.8852552, 0.8741705, 0.8577764, 0.8404956, 0.8164, 0.7898786, 0.8559746, 0.8882898), ncol=4, nrow=2)
  colnames(lambdas) = c("MD","MV","FD","FV")
  rownames(lambdas) = c("all", "tailed")
  lam = c(lambdas[1,paste0(sex,view)],lambdas[2,paste0(sex,view)])
  
  return (list(meanphen, data_sex, genre_data, sp_data, genresp_data, tipsgenre_data, listloc, ids, meanphen_shape, lam))
  
}



get_whole_data <- function(level = "form",
                           mode = 'mean', 
                           path_data_photos = "data_photos.csv", 
                           path_coords = "pca_embeddings.csv",
                           indices = 0,
                           mode_indices="id",
                           fileshape = "./data/pcaH.RData"){
  
  list_trees = list()
  list_tip_trees = c(0,0,0,0)
  #Males
  
  list_get_phenotype = get_phenotype_plusshape(c("M"),c("D"), mode = mode, level = level, path_coords = path_coords, path_data_photos = path_data_photos, fileshape=fileshape)
  meanphen <- list_get_phenotype[[1]]
  data <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  meanphen_shape <- list_get_phenotype[[9]]
  rm(list=c("list_get_phenotype"))
  
  select_genre = unique(data$Genre)
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly=T, tree_path = "./data/Papilionidae_MCC_clean.tre", 
                           genre=select_genre)#c("Graphium","Iphiclides","Protesilaus","Protographium","Mimoides","Parides"))
  
  subtree <- list_match[[1]]
  
  list_trees[[1]] = subtree
  list_tip_trees[1] = length(subtree$tip.label)
  
  meanphen <- list_match[[2]]
  
  meanphen_MD = meanphen
  meanphen_shape_MD = meanphen_shape
  
  list_get_phenotype = get_phenotype_plusshape(c("M"),c("V"), mode = mode, level = level, path_coords = path_coords, path_data_photos = path_data_photos, fileshape=fileshape)
  meanphen <- list_get_phenotype[[1]]
  data <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  meanphen_shape <- list_get_phenotype[[9]]
  rm(list=c("list_get_phenotype"))
  
  select_genre = unique(data$Genre)
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly=T, tree_path = "./data/Papilionidae_MCC_clean.tre", 
                           genre=select_genre)#c("Graphium","Iphiclides","Protesilaus","Protographium","Mimoides","Parides"))
  
  subtree <- list_match[[1]]
  
  list_trees[[2]] = subtree
  list_tip_trees[2] = length(subtree$tip.label)
  
  meanphen <- list_match[[2]]
  
  meanphen_MV = meanphen
  meanphen_shape_MV = meanphen_shape
  
  
  #Females
  
  list_get_phenotype = get_phenotype_plusshape(c("F"),c("D"), mode = mode, level = level, path_coords = path_coords, path_data_photos = path_data_photos, fileshape=fileshape)
  meanphen <- list_get_phenotype[[1]]
  data <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  meanphen_shape <- list_get_phenotype[[9]]
  rm(list=c("list_get_phenotype"))
  
  select_genre = unique(data$Genre)
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly=T, tree_path = "./data/Papilionidae_MCC_clean.tre", 
                           genre=select_genre)#c("Graphium","Iphiclides","Protesilaus","Protographium","Mimoides","Parides"))
  
  subtree <- list_match[[1]]
  
  list_trees[[3]] = subtree
  list_tip_trees[3] = length(subtree$tip.label)
  
  meanphen <- list_match[[2]]
  
  meanphen_FD = meanphen
  meanphen_shape_FD = meanphen_shape
  
  list_get_phenotype = get_phenotype_plusshape(c("F"),c("V"), mode = mode, level = level, path_coords = path_coords, path_data_photos = path_data_photos, fileshape=fileshape)
  meanphen <- list_get_phenotype[[1]]
  data <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  meanphen_shape <- list_get_phenotype[[9]]
  rm(list=c("list_get_phenotype"))
  
  select_genre = unique(data$Genre)
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly=T, tree_path = "./data/Papilionidae_MCC_clean.tre", 
                           genre=select_genre)#c("Graphium","Iphiclides","Protesilaus","Protographium","Mimoides","Parides"))
  
  subtree <- list_match[[1]]
  
  list_trees[[4]] = subtree
  list_tip_trees[4] = length(subtree$tip.label)
  
  meanphen <- list_match[[2]]
  
  
  meanphen_FV = meanphen
  meanphen_shape_FV = meanphen_shape
  
  
  rownames(meanphen_MD) = sapply(rownames(meanphen_MD), function(x) paste0(x,"MD"))
  rownames(meanphen_FD) = sapply(rownames(meanphen_FD), function(x) paste0(x,"FD"))
  rownames(meanphen_MV) = sapply(rownames(meanphen_MV), function(x) paste0(x,"MV"))
  rownames(meanphen_FV) = sapply(rownames(meanphen_FV), function(x) paste0(x,"FV"))
  
  rownames(meanphen_shape_MD) = sapply(rownames(meanphen_shape_MD), function(x) paste0(x,"MD"))
  rownames(meanphen_shape_FD) = sapply(rownames(meanphen_shape_FD), function(x) paste0(x,"FD"))
  rownames(meanphen_shape_MV) = sapply(rownames(meanphen_shape_MV), function(x) paste0(x,"MV"))
  rownames(meanphen_shape_FV) = sapply(rownames(meanphen_shape_FV), function(x) paste0(x,"FV"))
  
  meanphen_all = rbind(meanphen_MD, meanphen_MV, meanphen_FD, meanphen_FV)
  meanphen_shape_all = rbind(meanphen_shape_MD, meanphen_shape_MV, meanphen_shape_FD, meanphen_shape_FV)
  
  rownames(meanphen_all) <- sapply(rownames(meanphen_all),function(x) iconv(x, from = 'latin1', to = 'ASCII//TRANSLIT'))
  rownames(meanphen_shape_all) <- sapply(rownames(meanphen_shape_all),function(x) iconv(x, from = 'latin1', to = 'ASCII//TRANSLIT'))
  
  
  return(list(meanphen_all, meanphen_shape_all))
}