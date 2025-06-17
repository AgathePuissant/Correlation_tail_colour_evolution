loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

get_patches_data <- function(path_patches = "./data/patches_df.csv", 
                             path_data_photos = "./data/data_photos_old.csv", 
                             path_Species = "./data/Digitalisation.xlsx", 
                             path_cluster_shape = "./data/", 
                             sex=0, 
                             view=0, 
                             wing=0,
                             tail_clust = "clust"){
  
  patches <- read.csv(path_patches, sep=";")
  data_photos = data_photos <- read_csv2(path_data_photos)
  Species <- read_excel(path_Species)
  clust <- loadRData(paste0(path_cluster_shape,"cluster_shape_",sex,"V.RData"))
  
  
  data_photos = data_photos[,colnames(data_photos) %in% c("id","Genre","sp","ssp","form","SEXE")]
  data_photos$ssp[data_photos$ssp=='-'] <- data_photos$sp[data_photos$ssp=='-']
  data_photos$form[is.na(data_photos$form)] <- "NA"#data_photos$ssp[is.na(data_photos$form)]
  data_photos$tipsgenre <- paste(data_photos$Genre,data_photos$sp,data_photos$form,sep="_")
  data_photos$genresp <- paste(data_photos$Genre,data_photos$sp,sep="_")
  
  Species = Species[,colnames(Species) %in% c("CODE","Tail")]
  colnames(Species) <- c("id","Tail")
  Species$id = gsub("-","",Species$id)
  Species_info = Species
  Species_info$tipsgenre = data_photos[match(Species_info$id,data_photos$id),]$tipsgenre
  data_photos$tipsgenre <- NULL
  #put all data frames into list
  patches <- list(patches, data_photos, Species_info)
  
  #merge all data frames in list
  patches = patches %>% reduce(inner_join, by='id')
  patches = patches[!(patches$SEXE=="Ind" | is.na(patches$SEXE)),]
  
  
  patches$clust <- clust[paste0(patches$tipsgenre,sex,view),1]
  
  patches$major_axis_std = patches$major_axis
  patches$distance_contours_std = patches$distance_contours
  patches$distance_CU2_std = patches$distance_CU2
  patches$A_wing=NULL
  
  if (sex!=0){
    patches = patches[patches$SEXE==sex,]
  }
  if (view!=0){
    patches = patches[patches$view==view,]
  }
  if (wing!=0){
    patches = patches[patches$wing==gsub("_[a-z]*$","",wing),]
  }
  
  return(patches)
}



clean_landmark_data <- function(path_LM_HW = "C:/Users/agaca/Mon Drive/Ariane/Morphométrie/Data brutes/Hindwings/Brutes",
                                path_LM_FW = "C:/Users/agaca/Mon Drive/Ariane/Morphométrie/Data brutes/Forewings/Brutes",
                                nlines_HW = 175,
                                nlines_FW = 141){
  files_HW = list.files(path_LM_HW, full.names=T)
  files_HW = files_HW[grep(".tps",files_HW)]
  
  # Define the minimum number of lines/landmarks
  min_lines <- nlines_HW  # Adjust this number as needed
  
  # Initialize an empty list to store filtered files
  filtered_files <- list()
  
  # Loop through each tps file
  for (file in files_HW) {
    # Read the lines from the file
    lines <- readLines(file)
    # Count the number of lines/landmarks
    num_lines <- length(lines)
    
    remove_id = grep("ID|SCALE",lines)
    if (length(remove_id)>0){
      lines = lines[-remove_id]
    }
    lines = gsub(",",".",lines)
    
    
    
    # Check if the number of lines is greater than or equal to the minimum required
    if (num_lines >= min_lines) {
      
      writeLines(lines, file)
      # If it meets the criteria, add it to the filtered list
      filtered_files <- c(filtered_files, file)
    }
  }
  
  files_HW = unlist(filtered_files)
  
  HW = readmulti.tps(files_HW,readcurves=T,negNA=T)
  
  # dim(HW)
  
  # Check if all elements along the third dimension are NA
  all_na <- apply(HW, MARGIN = 3, FUN = function(x) all(is.na(x)))

  # Filter out elements along the third dimension that are not all NA
  HW <- HW[, , !all_na]
  files_HW <- files_HW[!all_na]
  # dim(HW)
  
  
  
  
  files_FW = list.files(path_LM_FW, full.names=T)
  files_FW = files_FW[grep(".tps",files_FW)]
  
  
  # Define the minimum number of lines/landmarks
  min_lines <- nlines_FW  # Adjust this number as needed
  
  # Initialize an empty list to store filtered files
  filtered_files <- list()
  
  # Loop through each tps file
  for (file in files_FW) {
    # Read the lines from the file
    lines <- readLines(file)
    # Count the number of lines/landmarks
    num_lines <- length(lines)
    
    remove_id = grep("ID|SCALE",lines)
    if (length(remove_id)>0){
      lines = lines[-remove_id]
    }
    lines = gsub(",",".",lines)
    
    
    
    
    # Check if the number of lines is greater than or equal to the minimum required
    if (num_lines >= min_lines) {
      writeLines(lines, file)
      # If it meets the criteria, add it to the filtered list
      filtered_files <- c(filtered_files, file)
    }
  }
  
  files_FW = unlist(filtered_files)
  
  FW = readmulti.tps(files_FW,readcurves=T)
  
  # dim(FW)
  
  # Check if all elements along the third dimension are NA
  all_na <- apply(FW, MARGIN = 3, FUN = function(x) all(is.na(x)))
  
  # Filter out elements along the third dimension that are not all NA
  FW <- FW[, , !all_na]
  # dim(FW)
  files_FW <- files_FW[!all_na]
  
  return(list(files_HW,HW,files_FW,FW))
}


# Function to calculate aspect ratio of a polygon
calculate_aspect_ratio <- function(coords) {
  # Extract x and y coordinates
  x <- coords$x
  y <- coords$y
  
  # Calculate bounding box
  bbox <- c(min(x), min(y), max(x), max(y))
  
  # Calculate width and height of bounding box
  width <- bbox[3] - bbox[1]
  height <- bbox[4] - bbox[2]
  
  # Calculate aspect ratio
  aspect_ratio <- max(width,height)/min(width,height)
  
  return(aspect_ratio)
}

deviation_barycentre = function(box, point, ind, LorR){
  
  H = list()
  
  H$x = box[,1,ind]
  
  H$y = box[,2,ind]
  # H
  
  
  
  p=list(x=point[1],y=point[2])
  
  # plot(c(0,4288), c(0,2848), type='n')
  # polygon(H, col=NULL, border='blue')
  # points(p)
  
  if (point.in.polygon(p$x, p$y, H$x, H$y)){
    
    # plot(c(0,4288), c(0,2848), type='n')
    # polygon(H, col=NULL, border='blue')
    # points(p)
    AR = calculate_aspect_ratio(H)
    # print(AR)
    
    centroid_x <- mean(H$x)
    centroid_y <- mean(H$y)
    
    distance_centroid = sqrt(sum((point - c(centroid_x, centroid_y))^2))  #/ AR
    
    return(list(distance_centroid, AR))
  }
}


get_granularity_data <- function(sex,
                                 view,
                                 wing_side = "Post",
                                 path_granularity_results = "./data/pattern_granularity_results.csv",
                                 path_data_photos = "./data/data_photos_old.csv",
                                 path_species = "./data/Digitalisation.xlsx"){
  
  data_photos = data_photos <- read_csv2(path_data_photos)
  Species <- read_excel(path_species)
  clust <- loadRData(paste0("./data/cluster_shape_",sex,"V.RData"))
  butterfly_areas <- read.csv("./data/butterfly_areas.csv", sep=";")
  
  pattern_granularity_results <- read.csv(path_granularity_results, sep=";")
  pattern_granularity_results$view = sapply(pattern_granularity_results$File,function(x) str_sub(x, -1))
  pattern_granularity_results$id <- gsub("D|V","",pattern_granularity_results$File)
  
  data_photos = data_photos[,colnames(data_photos) %in% c("id","Genre","sp","ssp","form","SEXE")]
  data_photos$ssp[data_photos$ssp=='-'] <- data_photos$sp[data_photos$ssp=='-']
  data_photos$form[is.na(data_photos$form)] <- "NA" #data_photos$ssp[is.na(data_photos$form)]
  data_photos$tipsgenre <- paste(data_photos$Genre,data_photos$sp,data_photos$form,sep="_")#,data_photos$ssp,data_photos$form,sep="_")
  data_photos$genresp <- paste(data_photos$Genre,data_photos$sp,sep="_")#,data_photos$ssp,data_photos$form,sep="_")
  
  Species = Species[,colnames(Species) %in% c("CODE","Tail")]
  colnames(Species) <- c("id","Tail")
  Species$id = gsub("-","",Species$id)
  Species_info = Species
  Species_info$tipsgenre = data_photos[match(Species_info$id,data_photos$id),]$tipsgenre
  data_photos$tipsgenre <- NULL
  #put all data frames into list
  df_list <- list(pattern_granularity_results, data_photos, Species_info)
  
  #merge all data frames in list
  df_list = df_list %>% reduce(inner_join, by='id')
  df_list = df_list[df_list$SEXE==sex,]
  df_list = df_list[df_list$view==view,]
  
  df_list$clust <- clust[df_list$genresp,1]
  
  butterfly_areas = butterfly_areas[butterfly_areas$view==view & butterfly_areas$wing==wing_side,,F]
  butterfly_areas = butterfly_areas %>% dplyr::group_by(id) %>% dplyr::summarise(A_wing=mean(A_wing,na.rm=T))
    
  df_list$A_wing = butterfly_areas[match(df_list$id,butterfly_areas$id),]$A_wing
  
  df_list$maxPower = as.numeric(df_list$maxPower)
  df_list$propPower = as.numeric(df_list$propPower)
  df_list <- df_list[!is.na(df_list$maxFreq),]
  
  df_list <- df_list[!is.na(df_list$tipsgenre),]
  
  return(df_list)
}
