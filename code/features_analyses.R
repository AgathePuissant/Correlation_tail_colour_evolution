#-------------------- Load libraries ---------------------
library(mvMORPH)
library(broom)
library(dplyr)
library(writexl)
library(grid)
library(sjPlot)
library(readr)
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(geomorph)
library(GEOmap)
library(sp)
library(pracma)
library(readr)
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(caper)
library(picante)
library(ape)
library(phylolm)
library(nlme)
library(rr2)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_patches_granularity.R")

# Function to remove duplicated dimnames along a specified dimension
remove_duplicated_dimnames <- function(arr, dim_index) {
  dim_names <- dimnames(arr)[[dim_index]]
  
  # Find unique names (first occurrence only)
  unique_indices <- !duplicated(dim_names)
  
  # Filter array to keep only unique slices along the specified dimension
  if (dim_index == 1) {
    arr_filtered <- arr[unique_indices, , , drop = FALSE]
  } else if (dim_index == 2) {
    arr_filtered <- arr[, unique_indices, , drop = FALSE]
  } else if (dim_index == 3) {
    arr_filtered <- arr[, , unique_indices, drop = FALSE]
  }
  
  # Update dimnames to reflect the filtered array
  new_dimnames <- dimnames(arr)
  new_dimnames[[dim_index]] <- dim_names[unique_indices]
  dimnames(arr_filtered) <- new_dimnames
  
  return(arr_filtered)
}

# Function to extract coefficients from a phylolm object
extract_phylolm_info <- function(model, sex, view, onlytailed, gran_or_pat) {
  summ <- summary(model)
  coefs <- as.data.frame(summ$coefficients)
  coefs$term <- rownames(coefs)
  coefs$sex <- sex
  coefs$view <- view
  coefs$onlytailed <- onlytailed
  coefs$gran_or_pat <- gran_or_pat
  coefs <- coefs[, c("term", "Estimate", "StdErr", "t.value", "p.value", "sex", "view", "onlytailed", "gran_or_pat")]
  return(coefs)
}

#--------------Define variables---------------------
model_results <- list()
to_test_sex = c("M","F")
to_test_view = c("D","V")
to_test_onlytailed = c(T,F) #True : only tailed dataset, False : All species dataset
to_test_granularity = c(T,F) #True : test for granularity, False : test for spots measurements

height_scale = 2848 / 500 #Size of original picture vs size of image used for patch detection

#----------------Analyses loop-------------------
for (s in to_test_sex){
  
  sex=s
  
  for (v in to_test_view){
    
    view=v
    
    
    #---------- Get mean shape per species and calculate lambda --------------------
    
    nlm = 51
    gpa = loadRData("./data/coord_gpa.RData")
    
    data_photos <- read_delim("./data/data_photos_old.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
    
    alterne=two.d.array(gpa)
    colnames(alterne)=paste0("LM",c(1:ncol(alterne)))
    alterne=as.data.frame(alterne)
    alterne$id = rownames(alterne)
    alterne = merge(alterne,data_photos,by=c("id"),all.x=T)
    
    alterne$tipsgenre = paste(alterne$Genre,alterne$sp,alterne$form, sep = "_")
    # Convert character columns to numeric
    alterne <- alterne %>%
      dplyr::mutate(across(starts_with("LM"), as.numeric))
    alterne <- alterne %>%
      dplyr::group_by(tipsgenre, SEXE) %>%
      dplyr::summarise(across(starts_with("LM"), mean, na.rm = TRUE))
    alterne <- as.data.frame(alterne)
    alterne <- alterne[alterne$SEXE==sex,,F]
    alterne = alterne[!is.na(alterne$tipsgenre),,F]
    rownames(alterne)<-alterne$tipsgenre
    
    alterne=alterne[complete.cases(alterne),]
    alterne$tipsgenre <- NULL
    alterne$SEXE <- NULL
    alterne = arrayspecs(alterne, p=nlm, k=2)
    
    Sp_info = dimnames(alterne)[[3]]
    
    #-------------------- Get patches data ---------------------
    patches_df <- get_patches_data(sex=sex, view=view, path_patches = "./data/patches_df.csv")
    patches_df$x = height_scale*patches_df$x
    patches_df$y = 2848 - height_scale*patches_df$y
    
    patches_df$distance_centroid = NA
    patches_df$box = NA
    patches_df$AR_window = NA
    #-------------------- Get raw landmarks data ---------------------
    list_LM = clean_landmark_data() 
    files_HW = list_LM[[1]]
    HW = list_LM[[2]]
    files_FW = list_LM[[3]]
    FW = list_LM[[4]]
    
    #------------Create wing cell------------
    
    Aprime = c(19,1,2,18)
    A = c(18,2,3,17)
    B = c(16,3,4,15)
    C = c(15,4,5,14)
    D = c(14,5,6,13)
    E = c(13,6,7,12)
    Fe = c(12,7,8,11)
    Fprime = c(11,8,9,10)
    AF = c(17,16,15,14,13,12,11,10)
    
    G = c(17,1,2,18)
    H = c(18,2,3,18)
    I = c(17,3,4,16)
    J = c(16,4,5,15)
    K = c(15,5,6,14)
    L = c(14,6,7,13)
    M = c(13,7,8,12)
    N = c(12,8,9,11)
    O = c(11,9,10,11)
    
    boxAprime = HW[Aprime,,]
    boxA = HW[A,,]
    boxB = HW[B,,]
    boxC = HW[C,,]
    boxD = HW[D,,]
    boxE = HW[E,,]
    boxF = HW[Fe,,]
    boxFprime = HW[Fprime,,]
    boxAF = HW[AF,,]
    
    boxG = FW[G,,]
    boxH = FW[H,,]
    boxI = FW[I,,]
    boxJ = FW[J,,]
    boxK = FW[K,,]
    boxL = FW[L,,]
    boxM = FW[M,,]
    boxN = FW[N,,]
    boxO = FW[O,,]
    
    list_boxes_Post = list(boxAprime,boxA,boxB,boxC,boxD,boxE,boxF,boxFprime,boxAF)
    
    list_names_box_Post = c("boxAprime","boxA","boxB","boxC","boxD","boxE","boxF","boxFprime","boxAF")
    
    #------------- Compute distance to centroid ------------------------------------
    
    
    for (i in c(1:length(patches_df$id))){
      
        individu = match(patches_df[i,]$id,gsub("-","",gsub("_H[L|R].tps","",str_split(files_HW,"/", simplify=T)[,10])))
        
        for (b in c(1:length(list_boxes_Post))){
          list_distance_centroid = deviation_barycentre(box=list_boxes_Post[[b]], point = c(patches_df[i,]$x,patches_df[i,]$y), ind=individu, patches_df[i,]$side)
          
          distance_centroid = list_distance_centroid[[1]]
          AR_window = list_distance_centroid[[2]]
          
          if(!is.null(distance_centroid)){
            patches_df[i,]$distance_centroid = distance_centroid
            patches_df[i,]$AR_window = AR_window
            patches_df[i,]$box = list_names_box_Post[b]
          }
        }

    }
    
    patches_df$distance_centroid_std = patches_df$distance_centroid/patches_df$AR_window
    
    
    #---------Granularity or patches ?---------
    for (gop in to_test_granularity){
      gran_or_pat=gop #true if granularity else patches
      #---------- Get granularity analysis data --------------------------------------
      granularity_df <- get_granularity_data(sex,view)
      #---------------------- Mean measures per species ------------------------------
      
      granularity_df_sp = granularity_df %>% group_by(tipsgenre) %>% dplyr::summarise(mean_maxFreq = mean(maxFreq, na.rm=TRUE),
                                                                                      mean_maxPower = mean(maxPower, na.rm=TRUE),
                                                                                      mean_propPower = mean(propPower, na.rm=TRUE),
                                                                                      mean_lumMean = mean(lumMean, na.rm=T),
                                                                                      mean_sumPower = mean(sumPower, na.rm=T),
                                                                                      tail_status = tail(names(sort(table(clust))), 1),
                                                                                      binary_tail_status = tail(names(sort(table(Tail))), 1),
                                                                                      genus = tail(names(sort(table(Genre))), 1),
                                                                                      A_wing = mean(A_wing,na.rm=T))
      
      granularity_df_sp <- as.data.frame(granularity_df_sp)
      rownames(granularity_df_sp)<-granularity_df_sp$tipsgenre
      
      patches_df_sp = patches_df %>% group_by(tipsgenre) %>% dplyr::summarise(mean_dc_std = mean(distance_centroid_std, na.rm=T),
                                                                              mean_dcontours_std = mean(distance_contours_std, na.rm=T),
                                                                              mean_dcu2_std = mean(distance_CU2_std, na.rm=T),
                                                                              mean_ma_std = mean(major_axis_std, na.rm=T))
      patches_df_sp = as.data.frame(patches_df_sp)
      rownames(patches_df_sp)<-patches_df_sp$genresp
      granularity_df_sp$mean_maxFreq = granularity_df_sp$mean_maxFreq/granularity_df_sp$A_wing
      
      if(gran_or_pat==T){
        merged_df_sp<- granularity_df_sp
      }else{
        merged_df_sp <- merge(patches_df_sp,granularity_df_sp, by="tipsgenre")
        merged_df_sp$mean_dc_std=merged_df_sp$mean_dc_std/merged_df_sp$A_wing
        merged_df_sp$mean_dcontours_std=merged_df_sp$mean_dcontours_std/merged_df_sp$A_wing
        merged_df_sp$mean_dcu2_std=merged_df_sp$mean_dcu2_std/merged_df_sp$A_wing
      }
      rownames(merged_df_sp) <- merged_df_sp$tipsgenre
      #------------ Set restrictions -------------------------------------------------
      
      for (tailed in to_test_onlytailed){
        onlytailed=tailed
        #------------ Match after restrictions -----------------------------------------
        shape_data = alterne
        
        matched_df_sp = merged_df_sp[merged_df_sp$tipsgenre %in% dimnames(shape_data)[[3]],]
        
        matched_df_sp = matched_df_sp[complete.cases(matched_df_sp),]
        
        
        if (onlytailed==T){
          matched_df_sp= matched_df_sp[matched_df_sp$binary_tail_status %in% c("T"),,F]
        }
        
        
        matched_df_sp$genresp = gsub("_NA|_[a-z]*$","",matched_df_sp$tipsgenre)
        # Removing rows with duplicated values
        matched_df_sp <- matched_df_sp %>%
          group_by(genresp) %>%
          filter(n() == 1) %>%
          ungroup()
        matched_df_sp$tipsgenre = matched_df_sp$genresp
        matched_df_sp = as.data.frame(matched_df_sp)
        rownames(matched_df_sp) <- matched_df_sp$tipsgenre
        
        list_match_tree <- match_tree(add_poly=F,meanphen_match = matched_df_sp, data_match = matched_df_sp)
        subtree <- list_match_tree[[1]]
        
        
        
        dimnames(shape_data)[[3]]<-gsub("_NA|_[a-z]*$","",dimnames(shape_data)[[3]])
        shape_data <- remove_duplicated_dimnames(shape_data, 3)
        
        # Use which() to get the indices of the names in the third dimension
        indices <- which(dimnames(shape_data)[[3]] %in% matched_df_sp$tipsgenre)
        
        # Reorder the indices to match the order in name_vector
        indices <- match(matched_df_sp$tipsgenre, dimnames(shape_data)[[3]])
        
        # Extract the corresponding slices from the array
        shape_data <- shape_data[, , indices]
        
        dim(shape_data)
        Sp_info = dimnames(shape_data)[[3]]
        
        #------------------ Phylogenetic regression ------------------------------------
        
        gdf <- geomorph.data.frame(shape_data = shape_data,
                                   species = Sp_info,
                                   A_wing = log(matched_df_sp$A_wing),
                                   genus = as.factor(matched_df_sp$genus),
                                   bin_tail = as.factor(matched_df_sp$binary_tail_status)) # geomorph data frame
        
        
        allometry <- procD.lm(shape_data ~ A_wing, data = gdf,
                              verbose=F)
        summary(allometry)
        
        
        if(is.null(allometry$aov.table$`Pr(>F)`)){
          allometry$aov.table$`Pr(>F)` = 1
        }
        
        PCA <- gm.prcomp(shape_data)
        PC <- PCA$x[,1]
        
        
        preds <- shape.predictor(shape_data, x= PC, Intercept = FALSE,
                                 pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
        
        
        par(mfrow=c(1,2))
        plotRefToTarget(mshape(shape_data), preds$pred1, method="TPS")
        plotRefToTarget(mshape(shape_data), preds$pred2, method="TPS")
        par(mfrow=c(1,1))
        
        meanshape = mshape(shape_data)
        
        
        df_comp = matched_df_sp
        
        # Scale predictors and response
        df_comp <- df_comp %>% dplyr::mutate(across(where(is.numeric), scale))
        
        df_comp$PC = PC[rownames(df_comp)]
        
        
        #Reorder PC so that positive values are elongated tails
        if(interlmkdist(preds$pred1, c(14,48)) < interlmkdist(preds$pred2, c(14,48))){
          df_comp$PC = df_comp$PC*-1
          saved = preds$pred1
          preds$pred1 = preds$pred2
          preds$pred2 = saved
        }
        
        
        comp_data = comparative.data(phy=subtree,data=df_comp, names.col=tipsgenre)
        
        
        #----------------Model comparison---------------
        
        
        # Define your response variable and predictor variables
        
        if (gop==F){
          predictors <- c(
            "mean_dcontours_std",
            "mean_dcu2_std",
            "mean_dc_std"
          )  
        }else{
          predictors <- c(
            "mean_maxFreq",
            "mean_maxPower",
            "mean_propPower"
          )  
        }
        response <- "PC"
        
        spnames = comp_data$data$genresp
        
        # Create the formula for main effects and two-way interactions only
        if(allometry$aov.table$`Pr(>F)`[1]<0.05){
          main_effects <- paste(c("A_wing", predictors), collapse = " + ") # "Csize",  
        }else{
          main_effects <- paste(predictors, collapse = " + ")
        }
        
        if (length(predictors)>1){
          two_way_interactions <- combn(predictors, 2, function(x) paste(x, collapse = ":"))
          two_way_interactions_formula <- paste(two_way_interactions, collapse = " + ")
          # Combine the main effects and two-way interactions into the formula
          new_formula <- as.formula(paste(response, "~", main_effects, "+", two_way_interactions_formula))
          
        }else{
          new_formula <- as.formula(paste(response, "~", main_effects))
        }
       
        pat.pgls = phylostep(new_formula,
                             data=comp_data$data,
                             model="lambda",
                             phy=subtree,
                             direction="both")
        summary(pat.pgls)
        
        new_formula2 = update(formula(pat.pgls$formula), . ~ .)
        pat.pgls = phylolm(new_formula2,data=comp_data$data,model="lambda",phy=subtree)
        summary(pat.pgls)
        
        
        # Extract and store
        model_name <- paste(sex, view, onlytailed,gran_or_pat, sep = "_")
        model_results[[model_name]] <- extract_phylolm_info(pat.pgls, sex, view, onlytailed,gran_or_pat)
        
      }
      
    }
  }
  
}

# Combine into a single data frame and export
all_results <- do.call(rbind, model_results)
write_xlsx(list("phylolm_results" = all_results), path = "./results/phylolm_regression_results.xlsx")
