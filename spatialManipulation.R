rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Desktop/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

library(letsR)
library(rgdal)
library(terra)
library(sf)
library(ggplot2)
library(picante)
library(broom)
library(tidyverse)

set.seed(500)

# # loading workspace
# load(file.path(data.dir,"fd.RData"))

## Creating original site x species matrix======================================

# read in tdwg level 3 shapefile
tdwg_l3 = vect(file.path(raw.dir, "SpatialData/TDWG/level3/level3.shp"))

# frug_birds is already filtered from BOTW.gdb for birds that are frugivorous, 
# presence, origin, and seasonality are all 1,2,3, and all spatial features with
# same species names are dissolved into the same feature
# all species are extant, with data deficient species included
frug_birds = vect(file.path(data.dir,"Spatial/Bird Range/all_frug_range_dissolved.shp"))

# remove species that do not intersect at all with tdwg_l3
overlap_check = terra::relate(x = frug_birds, y = tdwg_l3, relation = "intersects")

frug_birds_overlap = frug_birds[as.logical(rowSums(overlap_check))]

# create grid cells with 1 degree x 1 degree dimension
world_grid = st_make_grid(as(as(tdwg_l3, "Spatial"), "sf"),
                          1,
                          crs = st_crs(tdwg_l3),
                          what = "polygons",
                          square = TRUE) 

# limit grid to the spatial extent of tdwg_l3
world_grid_vect = vect(world_grid)
world_grid_vect_subset = world_grid_vect[as.logical(rowSums(relate(world_grid_vect, tdwg_l3, relation = "intersects")))]
# create grid id
world_grid_vect_subset$grid_id = seq(1,nrow(world_grid_vect_subset),1)

# create original site x species matrix
grid_occ_overlap = relate(x = world_grid_vect_subset,
                          y = frug_birds_overlap,
                          relation = "intersects") #contains all created grids

# replace column and grid names of site x species matrix
bird_names = frug_birds_overlap$sci_name
colnames(grid_occ_overlap) = bird_names
grid_ids = world_grid_vect_subset$grid_id
rownames(grid_occ_overlap) = grid_ids


# subset to grids that have at least 1 species
grid_occ_overlap_filled = grid_occ_overlap[as.logical(rowSums(grid_occ_overlap)),]
world_grid_vect_filled = world_grid_vect_subset[as.logical(rowSums(grid_occ_overlap))]

# converting site x species to a matrix with 1s and 0s
site_species_orig = as.matrix(grid_occ_overlap_filled) * 1

# sub-setting to species which appear in at least one grid
site_species_orig = site_species_orig[,as.logical(colSums(site_species_orig))]

## do i even need to subset to grids with at least 1 species?

  
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
  }else { # if threat 2 is entered as is valid
    birds_defaun = bird_data$sci_name[bird_data[,threat1] > 5 | bird_data[,threat2] > 5]
    for (i in 1:length(birds_defaun)){
      occ[,birds_defaun[i]] = 0
    }
    return (occ)
  }
  
}

# simulate habitat loss and wildlife trade
site_species_HL = defaunate(occ = site_species_orig, bird_data = values(frug_birds_overlap), 
                            threat1 = "HL")
site_species_WT = defaunate(occ = site_species_orig, bird_data = values(frug_birds_overlap), 
                            threat1 = "WT")
site_species_HL_WT = defaunate(occ = site_species_orig, bird_data = values(frug_birds_overlap), 
                            threat1 = "HL",  threat2 = "WT")


## FD Analysis =================================================================

# obtaining species x traits matrix
all_traits = values(frug_birds_overlap)
rownames(all_traits) = all_traits$sci_name
all_traits$Beak.Lengt = as.numeric(all_traits$Beak.Lengt) # length of culmen
all_traits$Mass = as.numeric(all_traits$Mass)
all_traits$Beak.Width = as.numeric(all_traits$Beak.Width)
all_traits$Kipps.Dist = as.numeric(all_traits$Kipps.Dist)
all_traits$Wing.Lengt = as.numeric(all_traits$Wing.Lengt)
all_traits$Kipps.Index = all_traits$Kipps.Dist / all_traits$Wing.Lengt
fd_traits = all_traits[,c("Mass","Beak.Lengt","Beak.Width","Kipps.Index")]

# subset species x traits matrix to original species
fd_traits_orig <- fd_traits[match(colnames(site_species_orig),
                                             rownames(fd_traits)),]


# subset site x species matrices to those plots with at least 3 species?
# currently, all scenario site x species matrices have same number of sites and species
# species x traits using original species


# loading fdFunc script for funcdiv function
source(file.path(main.dir,"src/fdFunc.R"))

# FD calculations
res_orig = funcdiv(x = fd_traits_orig, a = site_species_orig, original = TRUE)
res_HL = funcdiv(x = fd_traits_orig, a = site_species_HL)
res_WT = funcdiv(x = fd_traits_orig, a = site_species_WT)
res_HL_WT = funcdiv(x = fd_traits_orig, a = site_species_HL_WT)


## Null models =================================================================
# original assemblage
null_orig = list()
for(i in 1:100){
  null_orig[[i]] = randomizeMatrix(site_species_orig,null.model = "trialswap",
                                   iterations = 1000000)
}

null_res_orig = lapply(null_orig, FUN = funcdiv, x = fd_traits_orig, original = TRUE)

null_res_orig_fdis = lapply(null_res_orig, FUN = function(x) {x$FDis})

null_res_orig_fdis_mat = do.call("rbind",null_res_orig_fdis)

z_scores_orig = (res_orig$FDis - colMeans(null_res_orig_fdis_mat)) / apply(null_res_orig_fdis_mat, MARGIN = 2, FUN = sd)

# habitat loss defaunation assemblage
null_HL = list()
for(i in 1:100){
  null_HL[[i]] = randomizeMatrix(site_species_HL,null.model = "trialswap",
                                   iterations = 1000000)
}

null_res_HL = lapply(null_HL, FUN = funcdiv, x = fd_traits_orig)

null_res_HL_fdis = lapply(null_res_HL, FUN = function(x) {x$FDis})

null_res_HL_fdis_mat = do.call("rbind",null_res_HL_fdis)

z_scores_HL = (res_HL$FDis - colMeans(null_res_HL_fdis_mat)) / apply(null_res_HL_fdis_mat, MARGIN = 2, FUN = sd)

# wildlife trade defaunation assemblage
null_WT = list()
for(i in 1:100){
  null_WT[[i]] = randomizeMatrix(site_species_WT,null.model = "trialswap",
                                   iterations = 1000000)
}

null_res_WT = lapply(null_WT, FUN = funcdiv, x = fd_traits_orig)

null_res_WT_fdis = lapply(null_res_WT, FUN = function(x) {x$FDis})

null_res_WT_fdis_mat = do.call("rbind",null_res_WT_fdis)

z_scores_WT = (res_WT$FDis - colMeans(null_res_WT_fdis_mat)) / apply(null_res_WT_fdis_mat, MARGIN = 2, FUN = sd)


# habitat loss + wildlife trade defaunation assemblage
null_HL_WT = list()
for(i in 1:100){
  null_HL_WT[[i]] = randomizeMatrix(site_species_HL_WT,null.model = "trialswap",
                                 iterations = 1000000)
}

# change 100 to 999
# richness null model
# mclapply --> run things faster by parallelization

null_res_HL_WT = lapply(null_HL_WT, FUN = funcdiv, x = fd_traits_orig)

null_res_HL_WT_fdis = lapply(null_res_HL_WT, FUN = function(x) {x$FDis})

null_res_HL_WT_fdis_mat = do.call("rbind",null_res_HL_WT_fdis)

z_scores_HL_WT = (res_HL_WT$FDis - colMeans(null_res_HL_WT_fdis_mat)) / apply(null_res_HL_WT_fdis_mat, MARGIN = 2, FUN = sd)

# hist(z_scores_orig[1:50])
# hist(z_scores_orig[51:100])

z_scores = cbind(z_scores_orig, z_scores_HL, z_scores_WT, z_scores_HL_WT)
colnames(z_scores) = c("z_orig", "z_HL", "z_WT", "z_HL_WT")
f_disp = cbind(res_orig$FDis, res_HL$FDis, res_WT$FDis, res_HL_WT$FDis, 
               colMeans(null_res_orig_fdis_mat), colMeans(null_res_HL_fdis_mat), 
               colMeans(null_res_WT_fdis_mat), colMeans(null_res_HL_WT_fdis_mat))
colnames(f_disp) = c("fdisp_orig", "fdisp_HL", "fdisp_WT", "fdisp_HL_WT", 
                     "fdisp_orig_null", "fdisp_HL_null", "fdisp_WT_null", "fdisp_HL_WT_null")

f_disp_z = as.data.frame(cbind(f_disp, z_scores))

# checking if sequence of values in f_disp_z is same as grid_ids in world_grid_vect_filled
any(values(world_grid_vect_filled)$grid_id != as.numeric(rownames(f_disp_z)))
values(world_grid_vect_filled) = cbind(values(world_grid_vect_filled), f_disp_z)
values(world_grid_vect_filled)$fortify_id = rownames(values(world_grid_vect_filled))

# save.image(file = file.path(data.dir, "fd.RData"))



orig_sp = sum(as.logical(colSums(site_species_orig)))
HL_sp = sum(as.logical(colSums(site_species_HL)))
WT_sp = sum(as.logical(colSums(site_species_WT)))
HL_WT_sp = sum(as.logical(colSums(site_species_HL_WT)))

# total number of species globally for each scenario
global_sr = c(orig_sp,HL_sp,WT_sp,HL_WT_sp)
# species richness at cell level for each scenario
cell_sr = as.data.frame(cbind(sr_orig = res_orig$uniq.sp, sr_HL = res_HL$uniq.sp, 
                              sr_WT = res_WT$uniq.sp, sr_HL_WT = res_HL_WT$uniq.sp))
nrow(cell_sr)
nrow(cell_sr[cell_sr$sr_orig != cell_sr$sr_HL_WT,])

values(world_grid_vect_filled) = cbind(values(world_grid_vect_filled), 
                                       cell_sr[match(values(world_grid_vect_filled)$grid_id,
                                                     rownames(cell_sr)),])


## Plotting Maps ============================================
# fortify tdwg and grids and merge the data with the fortified data frame
tdwg_l3_fortified = fortify(as(tdwg_l3, "Spatial"))
grid_fortified = fortify(as(world_grid_vect_filled, "Spatial")) %>%
  left_join(., values(world_grid_vect_filled), by = c("id"="fortify_id"))

# original species richness
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_orig))

# HL species richness
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_HL))

# WT species richness
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_WT))

# HL_WT species richness
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_HL_WT))

# z_orig
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = ))

# ggplot() +
#   geom_polygon(data = fortify(as(tdwg_l3, "Spatial")), aes(x = long, y = lat, group = group), color = "black", fill = "white") +
#   geom_polygon(data = fortify(as(world_grid_vect_filled, "Spatial")), aes(x = long, y = lat, group = group), color = "red", fill = "red")
# 
# 
# ggplot() +
#   geom_polygon(data = fortify(as(tdwg_l3, "Spatial")), aes(x = long, y = lat, group = group), color = "black", fill = "white") +
#   geom_polygon(data = fortify(as(world_grid_vect_filled, "Spatial")), aes(x = long, y = lat, group = group), color = "red", fill = "red")

