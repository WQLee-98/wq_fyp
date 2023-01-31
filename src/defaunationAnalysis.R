rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Documents/FYP Working Folder"
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
library(viridis)


set.seed(500)

## Creating original site x species matrix======================================

# read in tdwg level 3 shapefile
tdwg_l3 = vect(file.path(raw.dir, "SpatialData/TDWG/level3/level3.shp"))

# read in ranges of all frugivorous birds
# frug_birds is already filtered from BOTW.gdb for birds that are frugivorous, 
# presence, origin, and seasonality are all 1,2,3, and all spatial features with
# same species names are dissolved into the same feature
# all species are extant, with data deficient species included
frug_birds = vect(file.path(data.dir,"Spatial/Bird Range 3/all_frug_range_filtered_dissolved.shp"))


# create grid cells with 5 degree x 5 degree dimension
world_grid = st_make_grid(as(as(tdwg_l3, "Spatial"), "sf"),
                          5,
                          crs = st_crs(tdwg_l3),
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

# saving site_species_orig
# saveRDS(site_species_orig, file = file.path(data.dir,'Defaunation/site_species_orig.rds'))

# saving all_traits
all_traits = values(frug_birds)
# saveRDS(all_traits, file = file.path(data.dir,'Defaunation/all_traits.rds'))


## read these in
# site_species_orig = read_rds(file.path(data.dir,'Defaunation/site_species_orig.rds'))
# all_traits = read_rds(file.path(data.dir,'Defaunation/all_traits.rds'))


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
site_species_HL = defaunate(occ = site_species_orig, bird_data = all_traits, 
                            threat1 = "HL")
site_species_WT = defaunate(occ = site_species_orig, bird_data = all_traits, 
                            threat1 = "WT")
site_species_HL_WT = defaunate(occ = site_species_orig, bird_data = all_traits, 
                            threat1 = "HL",  threat2 = "WT")



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
fd_traits_orig <- fd_traits[match(colnames(site_species_orig),
                                             rownames(fd_traits)),]


# subset site x species matrices to those plots with at least 3 species?
# currently, all scenario site x species matrices have same number of sites and species
# species x traits using original species



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

# create a matrix of distances with transformations (not sure why)
A <- matrix(0, ncol = n, nrow = n)

A[row(A) > col(A)] = -0.5 * x.dist^2 #why?

A = A + t(A) # t() transposes the dataframe

## ordination
# gower's double-centering
G = bicenter.wt(A)
e = eigen(G, symmetric = TRUE) # i assume this does the ordination

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

# saving out results
# saveRDS(data_final, file = file.path(results.dir,'Defaunation/data_final.rds'))
# write.csv(data_final, file = file.path(results.dir,'Defaunation/data_final.csv'))

## Plotting Maps ===============================================================
# fortify tdwg and grids and merge the data with the fortified data frame
tdwg_l3_fortified = fortify(as(tdwg_l3, "Spatial"))
grid_fortified = fortify(as(world_grid_vect_filled, "Spatial")) %>%
  left_join(., data_final, by = c("id"="fortify_id"))


# in progress!
ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_orig)) +
  scale_fill_viridis(option = 'mako') +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  theme_void()


# original species richness
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_orig)) +
  theme_void()

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

# z_HL
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_HL))

# z_WT
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_WT))

# z_HL_WT
ggplot() +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_HL_WT))
