sex="M"
view="D"
tailed_only=T

#-----------Load libraries and functions--------------------
library(ggplot2)
library(dplyr)
library(ggimage)
library(geomorph)
library(motmot)
library(dplyr)
library(stringr)
library(mvMORPH)
library(deeptime)
library(tidyverse)
library(readxl)
library(ggrepel)
library(viridis)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")

#-----------Load data------------

list_get_phenotype = get_phenotype_plusshape(sex,
                                             view, 
                                             mode = 'mean', 
                                             level = "form", 
                                             path_coords = "./data/pca_embeddings_wotailsDV.csv", 
                                             path_data_photos = "./data/data_photos.csv",
                                             fileshape="./data/pcaH_onlytails.RData") #_4
meanphen <- list_get_phenotype[[1]]
data <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
meanphen_shape <- list_get_phenotype[[9]]
lamb = list_get_phenotype[[10]]
if (tailed_only==T){
  lam= lamb[2]
}else{
  lam=lamb[1]
}
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data, add_poly=T, tree_path = "./data/Papilionidae_MCC_clean.tre")

subtree <- list_match[[1]]
meanphen <- list_match[[2]]

meanphen_shape <- meanphen_shape[match(rownames(meanphen),rownames(meanphen_shape)),]
meanphen_shape <- meanphen_shape[complete.cases(meanphen_shape), ]
meanphen <- meanphen[match(rownames(meanphen_shape),rownames(meanphen)),]
meanphen <- meanphen[complete.cases(meanphen), ]
meanphen <- meanphen[rownames(meanphen_shape),]
meanphen_shape <- meanphen_shape[rownames(meanphen),]

sp_biomes <- read.csv("./data/sp_biomes.csv", sep=";")
sp_biomes$N.A=NULL
sp_biomes = sp_biomes[rowSums(sp_biomes==0)!=14,,F]
biomes = as.data.frame(apply(sp_biomes,1,function(x) names(x)[which.max(x)]))
colnames(biomes)<-("biomes")

hostplant <- read.csv("./data/hostplant.csv", sep=";")
hostplant = hostplant %>%
  pivot_longer(cols = c(2:13)) %>%
  filter(value == 1)
hostplant$Species <- str_replace(hostplant$Species,"[0-9]","")
hostplant<-hostplant[!duplicated(hostplant$Species),]



Digitalisation <- read_excel("./data/Digitalisation.xlsx")
Digitalisation$tipsgenre = paste0(Digitalisation$Genre,"_",Digitalisation$Espèce,"_",Digitalisation$Forme)
tail_class <- Digitalisation[Digitalisation$SEXE==sex,c(13,14),,F]
tail_class<-tail_class[!duplicated(tail_class),,F]
tail_class <- as.data.frame(tail_class)
rownames(tail_class) <- tail_class$tipsgenre
tail_class <- tail_class[rownames(meanphen),]$Tail 
names(tail_class) <- rownames(meanphen) 
tail_class <- tail_class[complete.cases(tail_class)]
tail_class <- tail_class[tail_class=="T"]

if (tailed_only==T){
  
  meanphen <- meanphen[names(tail_class),]
  meanphen_shape <- meanphen_shape[names(tail_class),]
}

#---------------------run PLS-----------------------

PLS <- two.b.pls(A1 = as.matrix(meanphen), A2 = as.matrix(meanphen_shape))
plot(PLS)
summary(PLS)

#---------------------Plot PLS------------------------
nax = 1

coord_pls_color <- PLS$XScores[,nax]
coord_pls_shape <- PLS$YScores[,nax]

data_pls <- cbind(coord_pls_color,coord_pls_shape)
data_pls <- as.data.frame(data_pls)
data_pls$spnames<- (rownames(data_pls))
data_pls$biomes <- biomes[match(gsub("_NA|_[a-z]*$","",data_pls$spnames),rownames(biomes)),]
data_pls$hp <- hostplant[match(gsub("_NA|_[a-z]*$","",data_pls$spnames),hostplant$Species),]$name

data_pls$genresp = gsub("_NA|_[a-z]*$","",rownames(data_pls))
# Removing rows with duplicated values
data_pls <- data_pls %>%
  group_by(genresp) %>%
  filter(n() == 1) %>%
  ungroup()
data_pls = as.data.frame(data_pls)
rownames(data_pls) <- data_pls$genresp
data_pls$genresp <- NULL

list_match_tree <- match_tree(add_poly=F,meanphen_match = data_pls, data_match = data)
subtree <- list_match_tree[[1]]

if (sex=="M"){
  load("./data/tail_length_males.RData")
  tail_length = meanphen_Tail
}else{
  load("./data/tail_length_females.RData")
  tail_length=meanphen_Tail_female
}
data_pls$tail_length = tail_length[match(gsub("_NA|_[a-z]*$","",data_pls$spnames),rownames(tail_length)),1]
data_pls$tailed = tail_class[match(data_pls$spnames,names(tail_class))]
data_pls$tailed = is.na(data_pls$tailed)

gg <- ggplot(data_pls) +
  geom_phylomorpho(subtree, aes(x = coord_pls_color, y = coord_pls_shape))+
  geom_point(aes(x=coord_pls_color,y=coord_pls_shape, color=tail_length), size=5)+
  theme_test()+
  xlab("PLS coordinates for colour pattern")+
  ylab("PLS coordinates for tail shape")+
  scale_colour_viridis_c(option="magma")

gg



#-------Test for ecological variables---------------------


# Get the eigenvalues
eigenvalues <- PLS$svd$d
explained_variance <- eigenvalues^2 / sum(eigenvalues^2)
cumulative_variance <- cumsum(explained_variance)
threshold <- 0.8
num_axes <- which(cumulative_variance >= threshold)[1]

#Clean data for PCM
data_pls$genresp = rownames(data_pls)
# Removing rows with duplicated values
data_pls <- data_pls %>%
  group_by(genresp) %>%
  filter(n() == 1) %>%
  ungroup()
data_pls = as.data.frame(data_pls)
rownames(data_pls) <- data_pls$genresp
data_pls$genresp <- NULL
multi_pls = as.data.frame(cbind(PLS$XScores[,c(1:num_axes)],PLS$YScores[,c(1:num_axes)]))
multi_pls$genresp = gsub("_NA|_[a-z]*$","",rownames(multi_pls))
# Removing rows with duplicated values
multi_pls <- multi_pls %>%
  group_by(genresp) %>%
  filter(n() == 1) %>%
  ungroup()
multi_pls = as.data.frame(multi_pls)
rownames(multi_pls) <- multi_pls$genresp
multi_pls$genresp <- NULL
list_match_tree <- match_tree(add_poly=F,meanphen_match = multi_pls, data_match = data)
subtree <- list_match_tree[[1]]
multi_pls = as.matrix(multi_pls)
data_pls = data_pls[rownames(multi_pls),,F]


# hostplant##########################
data_pls = data_pls[complete.cases(data_pls$hp),,F]
multi_pls = multi_pls[rownames(data_pls),]
sum(rownames(multi_pls) %in% subtree$tip.label)
sum(rownames(data_pls) %in% subtree$tip.label)
subtree = match.phylo.data(subtree, data_pls)$phy

gg <- ggplot(data_pls) +
  geom_phylomorpho(subtree, aes(x = coord_pls_color, y = coord_pls_shape))+
  geom_point(aes(x=coord_pls_color,y=coord_pls_shape, color=hp), size=5)+
  theme_test()+
  xlab("PLS coordinates for colour pattern")+
  ylab("PLS coordinates for tail shape")
gg

#not run
# mvmorph_data = list(multi_pls = multi_pls, hp=as.factor(data_pls$hp))
# subtree = match.phylo.data(subtree,data_pls)$phy
# resglslambda = mvgls(data = mvmorph_data, multi_pls~hp, tree=subtree, model="lambda")
# res_manova = manova.gls(resglslambda)
# res_manova

#biomes######################

biome_palette <- c(
  "Deserts...Xeric.Shrublands" = "#E69F00",  # Vibrant orange for deserts
  "Mediterranean.Forests..Woodlands...Scrub" = "#56B4E9",  # Bright sky blue
  "Montane.Grasslands...Shrublands" = "#009E73",  # Bold teal-green for highlands
  "Temperate.Broadleaf...Mixed.Forests" = "#F0E442",  # Bright yellow
  "Temperate.Conifer.Forests" = "#0072B2",  # Strong deep blue
  "Temperate.Grasslands..Savannas...Shrublands" = "#D55E00",  # Bold burnt orange
  "Tropical...Subtropical.Coniferous.Forests" = "#CC79A7",  # Distinct magenta
  "Tropical...Subtropical.Dry.Broadleaf.Forests" = "#8E44AD",  # Strong purple
  "Tropical...Subtropical.Grasslands..Savannas...Shrublands" = "#F781BF",  # Pink for warm tropical savannas
  "Tropical...Subtropical.Moist.Broadleaf.Forests" = "#009688"  # Bright cyan-turquoise
)


gg <- ggplot(data_pls) +
  geom_phylomorpho(subtree, aes(x = coord_pls_color, y = coord_pls_shape))+
  geom_point(aes(x=coord_pls_color,y=coord_pls_shape, color=biomes), size=5)+
  theme_test()+
  xlab("PLS coordinates for colour pattern")+
  ylab("PLS coordinates for tail shape")+
  coord_fixed(ratio=50)

gg

sum(rownames(multi_pls) == rownames(data_pls))/length(rownames(data_pls))
data_pls = data_pls[!is.na(data_pls$biomes),,F]
subtree = match.phylo.data(subtree,data_pls)$phy
multi_pls = multi_pls[rownames(data_pls),]

#not run
# mvmorph_data = list(multi_pls = multi_pls,biomes = data_pls$biomes)
# resglslambda = mvgls(data = mvmorph_data, multi_pls~biomes, tree=subtree, model="lambda")
# res_manova = manova.gls(resglslambda)
# res_manova



#---------------Plot with images----------------------------

# Select a random sample of 1 from each group in 'genresp'
randomsample <- data %>%
  group_by(genresp) %>%
  sample_n(1)

# Generate file paths
randomsample <- randomsample %>%
  mutate(id = paste0('thumbnails/', id, view,'.png'))

# Update row names to match 'genresp'
rownames(randomsample) <- randomsample$genresp

# Map the file paths to the 'data_pls' dataframe
data_pls$id <- randomsample[match(data_pls$spnames, rownames(randomsample)), "id"]$id

# Remove rows from 'data_pls' where 'id' is NA (no corresponding image)
data_pls <- data_pls %>%
  filter(!is.na(id))

# Offset factor to move point
tail_x_offset <- 1.2
tail_y_offset <- 0

data_pls <- data_pls %>%
  mutate(
    tail_marker_x = coord_pls_color + tail_x_offset,
    tail_marker_y = coord_pls_shape - tail_y_offset
  )

ggplot(data_pls) +
  # Tail marker as a point
  # adjust point size as needed
  geom_point(aes(x = coord_pls_color, y = coord_pls_shape, color = tail_length), size=3) +
  # Image of butterfly wings without tails
  geom_image(aes(x = tail_marker_x, y = tail_marker_y, image = id),
             by = "height", size = 0.08, asp = 1, alpha=0.5) +
  
  scale_color_viridis_c(name = "Tail length", option = "plasma") +
  coord_fixed(ratio = 1) +
  theme_test() +
  theme(aspect.ratio = 0.5)
