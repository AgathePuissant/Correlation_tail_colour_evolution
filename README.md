[![DOI](https://zenodo.org/badge/590114695.svg)](https://doi.org/10.5281/zenodo.7665032)

detect_spots.py: Custom Python code for contrasting spot detection.

# In the folder code:
R code for recreating the analyses and plots. 
- features_analyses.R: code for the phylogenetic regression of spots measurements and granularity measurements on tail shape.
- PLS_analyses.R: code to perform PLS between tail shape and colour pattern ML embeddings and plot the result with thumbnail images.
- pPLS_analyses.R: code to perform phylogenetic PLS with intra-specific resampling to take account of phenotypic forms within species.
- RR_analyses.R: code to perform RRphylo evolutionary rates estimation with intra-specific resampling to take account of phenotypic forms.

# In the data folder:
- butterfly_areas.csv: wing area with tails removed, in pixels, for each specimen image.
- cluster_shape_FV: class of tail shapes for each species (females).
- cluster_shape_MV: class of tail shapes for each species (males)
- coord_gpa: coordinates of LM after Procrustes alignment for each specimen.
- data_photos and data_photos_old: csv containing information about the photographs (species, genus etc.)
- Digitalisation.xlsx: file containing information about the photographs (binary tail class)
- pca_embeddings_wotailsDV.csv: ML embeddings coordinates on the first 24 axis of the PCA.
- pcaH_onlytails.RData: Tail shape coordinates on the PCA axis.
- sp_biomes.csv: species main biomes.
- tail_length_females.RData: tail length of species (females).
- tail_length_males.RData: tail length of species (males).

# Results:
- phylolm_regression_results.csv: csv containing the results for all the features analyses.
- pPLS_table.RData: table containing the results of phylogenetic PLS for all the intra-specific resampling.
- RR_FD_100t.RData: results of the RRPhylo estimates for all the intra-specific resampling (M: Males, F: Females, D: Dorsal, V: Ventral)

# Thumbnails:
-Images of specimens with tails removed, resized and flipped for easier plotting of the results.

