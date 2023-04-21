rm(list = ls())
gc()

if(Sys.info()["user"] == "junyinglim"){
  main.dir = "OneDrive - National University of Singapore/PEEBL FYP Projects/Avian Functional Diversity Loss/"
} else {
  main.dir = "C:/Users/leewq/Documents/FYP Working Folder"  
}

raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

library(letsR)
library(terra)
library(sf)
library(ade4)
library(glue)
library(ggplot2)
library(ggrepel)
library(broom)
library(tidyverse)




set.seed(500)


## Creating original site x species matrix======================================
# read in tdwg level 3 shapefile
tdwg_l1 = vect(file.path(raw.dir, "SpatialData/TDWG/level1/level1.shp"))

# read in ranges of all frugivorous birds
# frug_birds is already filtered from BOTW.gdb for birds that are frugivorous, 
# presence and seasonality are 1,2,3 while origin is 1,2 and all spatial features with
# same species names are dissolved into the same feature
# all species are extant, with data deficient species included
frug_birds = vect(file.path(data.dir,"Spatial/Bird Range 3/all_frug_range_filtered_dissolved.shp"))


# create grid cells with 1 degree x 1 degree dimension
world_grid = st_make_grid(as(as(tdwg_l1, "Spatial"), "sf"),
                          1,
                          crs = st_crs(tdwg_l1),
                          what = "polygons",
                          square = TRUE)
world_grid_vect = vect(world_grid)

# create grid id
world_grid_vect$grid_id = seq(1,nrow(world_grid_vect),1)

# create original site x species matrix
grid_occ_overlap = relate(x = world_grid_vect,
                          y = frug_birds,
                          relation = "intersects") #contains all created grids

# replace column and grid names of site x species matrix
bird_names = frug_birds$sci_name
colnames(grid_occ_overlap) = bird_names
grid_ids = world_grid_vect$grid_id
rownames(grid_occ_overlap) = grid_ids

# subset to grids that have at least 1 species
grid_occ_overlap_filled = grid_occ_overlap[as.logical(rowSums(grid_occ_overlap)),]
world_grid_vect_filled = world_grid_vect[as.logical(rowSums(grid_occ_overlap))]

# converting site x species to a matrix with 1s and 0s
site_species_orig = as.matrix(grid_occ_overlap_filled) * 1

# sub-setting to species which appear in at least one grid
site_species_orig = site_species_orig[,as.logical(colSums(site_species_orig))]
all_traits = values(frug_birds)
all_traits[,c("HL","WT")] = sapply(all_traits[,c("HL","WT")], as.numeric)


## Simulating defaunation ======================================================

# create function that produces new site x species matrix based on simulation scenario
defaunate = function(occ, bird_data, threat1, threat2 = NULL){
  # Arguments
  #   occ: site x species matrix
  #   bird_data: bird traits data frame
  #   threat: threat type ("HL", "WT", "BI", "P", "CC")
  # Returns
  #   New site x species matrix based on simulation
  
  threat_types = c("HL", "WT", "BI", "P", "CC")
  
  # checking for invalid entry of threat1
  if (!(threat1 %in% threat_types)){
    print("Invalid threat1 type")
    return (F)
  }
  if (is.null(threat2)){
    # change presence/absence of species that have threat > 5 to 0 in matrix
    birds_defaun = bird_data$sci_name[bird_data[,threat1] > 5]
    for (i in 1:length(birds_defaun)){
      occ[,birds_defaun[i]] = 0
    }
    return (occ)
  }else if (!(threat2 %in% threat_types)){
    # if threat2 is entered, check for invalid entry of threat2
    print("Invalid threat2 type")
    return (F)
  }else { # if threat 2 is entered and is valid
    birds_defaun = bird_data$sci_name[bird_data[,threat1] > 5 | bird_data[,threat2] > 5]
    for (i in 1:length(birds_defaun)){
      occ[,birds_defaun[i]] = 0
    }
    return (occ)
  }
  
}

# simulate habitat loss and wildlife trade
site_species_HL = defaunate(occ = site_species_orig, bird_data = all_traits, 
                            threat1 = "HL")
site_species_WT = defaunate(occ = site_species_orig, bird_data = all_traits, 
                            threat1 = "WT")
site_species_HL_WT = defaunate(occ = site_species_orig, bird_data = all_traits, 
                            threat1 = "HL",  threat2 = "WT")



# species present in all grids in all scenarios
all_grid_species = list()
for (i in 1:nrow(site_species_orig)){
  orig_sp = colnames(site_species_orig)[as.logical(site_species_orig[i,])]
  HL_sp = colnames(site_species_HL)[as.logical(site_species_HL[i,])]
  WT_sp = colnames(site_species_WT)[as.logical(site_species_WT[i,])]
  HL_WT_sp = colnames(site_species_HL_WT)[as.logical(site_species_HL_WT[i,])]
  
  HL_sp_lost = orig_sp[!(orig_sp %in% HL_sp)]
  WT_sp_lost = orig_sp[!(orig_sp %in% WT_sp)]
  HL_WT_sp_lost = orig_sp[!(orig_sp %in% HL_WT_sp)]
  
  grid_species = list(orig_sp = orig_sp,HL_sp = HL_sp,WT_sp = WT_sp,HL_WT_sp = HL_WT_sp,
                      HL_sp_lost = HL_sp_lost, WT_sp_lost = WT_sp_lost, HL_WT_sp_lost = HL_WT_sp_lost)
  
  all_grid_species[[i]] = grid_species
}
names(all_grid_species) = rownames(site_species_orig)



# magnitude scores for all species in all cells
all_mag = list()
for (i in 1:nrow(site_species_orig)){
  cell_sp = colnames(site_species_orig)[as.logical(site_species_orig[i,])]
  cell_mag = all_traits[all_traits$sci_name %in% cell_sp,c("sci_name","HL","WT")]
  all_mag[[i]] = cell_mag
}
names(all_mag) = rownames(site_species_orig)

# average magnitudes for all cells
mean_mag = lapply(all_mag, FUN = function(x){
  cell_mean = apply(x[,c("HL","WT")], MARGIN = 2, FUN = mean)
}) %>%
  do.call("rbind",.) %>%
  as.data.frame(.)

mean_mag$grid_id = names(all_mag)
colnames(mean_mag) = c('mean_HL','mean_WT','grid_id')


## Preparing traits ============================================================

# obtaining species x traits matrix
rownames(all_traits) = all_traits$sci_name
all_traits$Beak.Lengt = as.numeric(all_traits$Beak.Lengt) # length of culmen
all_traits$Mass = as.numeric(all_traits$Mass)
all_traits$Beak.Width = as.numeric(all_traits$Beak.Width)
all_traits$Kipps.Dist = as.numeric(all_traits$Kipps.Dist)
all_traits$Wing.Lengt = as.numeric(all_traits$Wing.Lengt)
all_traits$Kipps.Index = all_traits$Kipps.Dist / all_traits$Wing.Lengt
fd_traits = all_traits[,c("Mass","Beak.Lengt","Beak.Width","Kipps.Index")]

# subset species x traits matrix to original species
fd_traits_orig = fd_traits[match(colnames(site_species_orig),
                                             rownames(fd_traits)),]



## PCoA ========================================================================
tol = 1e-07

## creating distance matrix 'x.dist'
# scaling traits
x.s <- apply(fd_traits_orig, 2, scale, center = TRUE, scale = TRUE)
x.dist <- dist(x.s)

x.rn <- row.names(fd_traits_orig)
attr(x.dist, "Labels" ) <- x.rn

# ordination
n <- attr(x.dist, "Size")

# create a matrix of distances with transformations
A <- matrix(0, ncol = n, nrow = n)

A[row(A) > col(A)] = -0.5 * x.dist^2 

A = A + t(A) # t() transposes the dataframe

## ordination
# gower's double-centering
G = bicenter.wt(A)
e = eigen(G, symmetric = TRUE) 

vectors <- e$vectors
eig <- e$values

# check if 'x.dist' is Euclidean or not
w0 <- eig[n] / eig[1]
if (w0 > -tol) r <- sum(eig > (eig[1] * tol)) else r <- length(eig)

# PCoA axes; coordinates of all species in morphospace
vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), r)
dimnames(vectors) <- list(colnames(site_species_orig), NULL)


## Functional Dispersion =======================================================

# loading fdFunc script for funcdiv function
source(file.path(main.dir,"src/fdFunc1.R"))

# FD calculations
res_orig = funcdiv(vectors = vectors, eig = eig, a = site_species_orig)
res_HL = funcdiv(vectors = vectors, eig = eig, a = site_species_HL)
res_WT = funcdiv(vectors = vectors, eig = eig, a = site_species_WT)
res_HL_WT = funcdiv(vectors = vectors, eig = eig, a = site_species_HL_WT)


orig_sp = sum(as.logical(colSums(site_species_orig)))
HL_sp = sum(as.logical(colSums(site_species_HL)))
WT_sp = sum(as.logical(colSums(site_species_WT)))
HL_WT_sp = sum(as.logical(colSums(site_species_HL_WT)))

# total number of species globally for each scenario
global_sr = c(orig_sp,HL_sp,WT_sp,HL_WT_sp)
# species richness at cell level for each scenario
cell_sr = as.data.frame(cbind(sr_orig = rowSums(site_species_orig), 
                              sr_HL = rowSums(site_species_HL), 
                              sr_WT = rowSums(site_species_WT), 
                              sr_HL_WT = rowSums(site_species_HL_WT)))


## Null models =================================================================

# loading nullModel script for null_mod function
source(file.path(main.dir,"src/nullModel.R"))

null_HL = null_mod(site_species_orig, cell_sr$sr_HL,vectors, eig)
null_HL_mat = do.call("rbind",null_HL)

null_WT = null_mod(site_species_orig, cell_sr$sr_WT,vectors, eig)
null_WT_mat = do.call("rbind",null_WT)

null_HL_WT = null_mod(site_species_orig, cell_sr$sr_HL_WT,vectors, eig)
null_HL_WT_mat = do.call("rbind",null_HL_WT)

z_HL = (res_HL$FDis - rowMeans(null_HL_mat))/apply(null_HL_mat, MARGIN = 1, FUN = sd)
z_WT = (res_WT$FDis - rowMeans(null_WT_mat))/apply(null_WT_mat, MARGIN = 1, FUN = sd)
z_HL_WT = (res_HL_WT$FDis - rowMeans(null_HL_WT_mat))/apply(null_HL_WT_mat, MARGIN = 1, FUN = sd)


# collating results
z_scores = cbind(z_HL, z_WT, z_HL_WT)
colnames(z_scores) = c("z_HL", "z_WT", "z_HL_WT")
f_disp = cbind(res_orig$FDis, res_HL$FDis, res_WT$FDis, res_HL_WT$FDis, 
               rowMeans(null_HL_mat), rowMeans(null_WT_mat), 
               rowMeans(null_HL_WT_mat))
colnames(f_disp) = c("fdisp_orig", "fdisp_HL", "fdisp_WT", "fdisp_HL_WT", 
                     "fdisp_HL_null", "fdisp_WT_null", "fdisp_HL_WT_null")
f_disp_z = as.data.frame(cbind(f_disp, z_scores))


# checking if sequence of values in f_disp_z is same as grid_ids in world_grid_vect_filled
any(values(world_grid_vect_filled)$grid_id != as.numeric(rownames(f_disp_z)))


# creating data frame with all required values
data_final = values(world_grid_vect_filled)
data_final$fortify_id = rownames(data_final)
data_final = cbind(data_final, f_disp_z, cell_sr[match(data_final$grid_id,
                                                       rownames(cell_sr)),])



## using cmdscale function for pcoa ============================================
pcoa = cmdscale(x.dist, k = 4, eig = T, add = T)

positions = as.data.frame(pcoa$points)
colnames(positions) = c("pcoa1", "pcoa2", "pcoa3", "pcoa4")

# adding taxonomy data
bird_tax = all_traits[,c("sci_name", "Family1", "Order1")]
bird_tax = bird_tax[match(bird_tax$sci_name,rownames(positions)),]
positions = cbind(positions, Family = bird_tax$Family1, Order = bird_tax$Order1)


# calculating percentage variance explained by each pcoa axis
percent_explained = (100 * pcoa$eig / sum(pcoa$eig))
rounded_pe = round(percent_explained[1:4],2)

labs = c(glue("PCoA 1 ({rounded_pe[1]}%)"),
         glue("PCoA 2 ({rounded_pe[2]}%)"))

# plot pcoa
positions_tibble = as_tibble(positions, rownames = "species")
pcoa_graph = ggplot(data = positions_tibble, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(color = Order)) +
  scale_colour_manual(values = c('#6B6100','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                          '#e31a1c','#fdbf6f','#ff7f00','#007E71','#6a3d9a',
                                          '#FFFF37','#b15928','#2A00FF','#DE0093')) +
                                            labs(x = labs[1], y = labs[2]) +
  geom_text_repel(aes(label=species)) + 
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11))
pcoa_graph


# plot percentage explained
pe_graph = tibble(pe = cumsum(rounded_pe),
                  axis = 1:length(rounded_pe)) %>%
  ggplot(aes(x=axis,y=pe)) +
  geom_line(size = 1) +
  geom_point()+
  coord_cartesian(xlim = c(1,4), ylim = c(0,100)) +
  labs(x = "PCoA Axis", y = "Cumulative Percentage Explained") + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))
pe_graph




## Comparing average traits of morphologically unique species vs average =======
HL_defaun_sp = colnames(site_species_HL)[!as.logical(colSums(site_species_HL))]
WT_defaun_sp = colnames(site_species_WT)[!as.logical(colSums(site_species_WT))]
HL_WT_defaun_sp = colnames(site_species_HL_WT)[!as.logical(colSums(site_species_HL_WT))]

fd_traits_orig
fd_traits_HL = fd_traits_orig[rownames(fd_traits_orig)%in%HL_defaun_sp,]
fd_traits_WT = fd_traits_orig[rownames(fd_traits_orig)%in%WT_defaun_sp,]
fd_traits_HL_WT = fd_traits_orig[rownames(fd_traits_orig)%in%HL_WT_defaun_sp,]

fd_traits_HL$Type = rep('HL', nrow(fd_traits_HL))
fd_traits_WT$Type = rep('WT', nrow(fd_traits_WT))
fd_traits_orig$Type = rep('All', nrow(fd_traits_orig))
fd_traits_all = rbind(fd_traits_orig,fd_traits_HL,fd_traits_WT)


# mass
library(plyr)
mu_mass = ddply(fd_traits_all, "Type", summarise, grp.mean=mean(log(Mass)))

dens_mass = ggplot(fd_traits_all, aes(x = log(Mass), color = Type, fill = Type)) +
  geom_density(alpha = 0.2) +
  geom_vline(data = mu_mass, aes(xintercept=grp.mean, color=Type),
           linetype="dashed")

# beak length
mu_bl = ddply(fd_traits_all, "Type", summarise, grp.mean=mean(log(Beak.Lengt)))

dens_bl = ggplot(fd_traits_all, aes(x = log(Beak.Lengt), color = Type, fill = Type)) +
  geom_density(alpha = 0.2) +
  geom_vline(data = mu_bl, aes(xintercept=grp.mean, color=Type),
             linetype="dashed")

# beak width
mu_bw = ddply(fd_traits_all, "Type", summarise, grp.mean=mean(log(Beak.Width)))

dens_bw = ggplot(fd_traits_all, aes(x = log(Beak.Width), color = Type, fill = Type)) +
  geom_density(alpha = 0.2) +
  geom_vline(data = mu_bw, aes(xintercept=grp.mean, color=Type),
             linetype="dashed")

# kipp's index
mu_kipps = ddply(fd_traits_all, "Type", summarise, grp.mean=mean(Kipps.Index))

dens_kipps = ggplot(fd_traits_all, aes(x = Kipps.Index, color = Type, fill = Type)) +
  geom_density(alpha = 0.2) +
  geom_vline(data = mu_kipps, aes(xintercept=grp.mean, color=Type),
             linetype="dashed")


## Morphological uniqueness ~ IUCN threat category and HL/WT====================

# creating a function that takes in the vectors of all frugivores in morphospace
# and returns their Euclidean distance to centroid
morpho_uniq = function(vectors){
  centroid = apply(vectors, 2, mean)
  dist.pos = sweep(vectors, 2, centroid) 
  dist.pos = rowSums(dist.pos^2)
  zij = sqrt(abs(dist.pos))
  return(zij)
}

all_frug_uniq = data.frame(uniq = morpho_uniq(vectors))
IUCN_threat = all_traits$status[match(rownames(all_frug_uniq),all_traits$sci_name)]
HL_mag = all_traits$HL[match(rownames(all_frug_uniq),all_traits$sci_name)]
WT_mag = all_traits$WT[match(rownames(all_frug_uniq),all_traits$sci_name)]
all_frug_uniq = cbind(all_frug_uniq, IUCN_threat, HL, WT)
all_frug_uniq$HL_bin = ifelse(all_frug_uniq$HL > 5, TRUE, FALSE)
all_frug_uniq$WT_bin = ifelse(all_frug_uniq$WT > 5, TRUE, FALSE)
all_frug_uniq$IUCN_threat = factor(all_frug_uniq$IUCN_threat)

HL = all_frug_uniq[all_frug_uniq$HL_bin == T,]
WT = all_frug_uniq[all_frug_uniq$WT_bin == T,]

# all_frug_uniq_melt = melt(all_frug_uniq,)

ggplot(data = all_frug_uniq, aes(x=IUCN_threat, y=uniq)) +
  geom_boxplot()

ggplot(data = all_frug_uniq, aes(y=uniq)) +
  geom_boxplot()
ggplot(data = HL, aes(y=uniq)) +
  geom_boxplot()
ggplot(data = WT, aes(y=uniq)) +
  geom_boxplot()

fd_traits[match(colnames(site_species_orig),
               rownames(fd_traits)),]



## saving files ================================================================
## saving world_grid_vect_filled
# saveRDS(world_grid_vect_filled, file = file.path(data.dir,'Defaunation/world_grid_vect_filled_1degree.rds'))

## saving site_species_orig
# saveRDS(site_species_orig, file = file.path(data.dir,'Defaunation/site_species_orig_1degree.rds'))

## saving all_traits
# saveRDS(all_traits, file = file.path(data.dir,'Defaunation/all_traits_1degree.rds'))

## saving mean HL/WT magnitudes per cell
# saveRDS(mean_mag, file = file.path(data.dir,'Defaunation/mean_mag_1degree.rds'))

## saving out results
# saveRDS(data_final, file = file.path(results.dir,'Defaunation/data_final_1degree.rds'))
# write.csv(data_final, file = file.path(results.dir,'Defaunation/data_final_1degree.csv'))

## saving world_grid_vect_fill with final data
# world_grid_vect_filled_data = world_grid_vect_filled
# values(world_grid_vect_filled_data) = left_join(values(world_grid_vect_filled_data), data_final, by = c("grid_id"="grid_id"))
# writeVector(x = world_grid_vect_filled_data, filename = file.path(results.dir,'grids_1degree'), filetype = "ESRI Shapefile")

# # fortify grid cells and joining final data and mean_magnitudes, and saving as RDS
# mean_mag$grid_id = as.numeric(mean_mag$grid_id)
# grid_fortified = fortify(as(world_grid_vect_filled, "Spatial")) %>%
#   left_join(., data_final, by = c("id"="fortify_id")) %>%
#   left_join(.,mean_mag, by = c("grid_id"="grid_id"))
# saveRDS(grid_fortified, file = file.path(results.dir,'Defaunation/grid_fortified.rds'))

# # reading in dissolved tdwg level 1, fortify, and save as RDS
# tdwg_l1_dissolved = vect(file.path(raw.dir, "SpatialData/TDWG/level1/level1_dissolved.shp"))
# tdwg_l1_fortified = fortify(as(tdwg_l1_dissolved, "Spatial"))
# saveRDS(tdwg_l1_fortified, file = file.path(data.dir,'Defaunation/tdwg_l1_fortified.rds'))

## save pcoa graphs in RDS
# saveRDS(pcoa_graph, file = file.path(results.dir,'Defaunation/pcoa_graph.rds'))
# saveRDS(pe_graph, file = file.path(results.dir,'Defaunation/pe_graph.rds'))

## save density plots for morphological traits in RDS
# saveRDS(dens_mass, file = file.path(results.dir,'Defaunation/dens_mass.rds'))
# saveRDS(dens_bl, file = file.path(results.dir,'Defaunation/dens_bl.rds'))
# saveRDS(dens_bw, file = file.path(results.dir,'Defaunation/dens_bw.rds'))
# saveRDS(dens_kipps, file = file.path(results.dir,'Defaunation/dens_kipps.rds'))
