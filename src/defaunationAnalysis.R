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
library(ade4)
library(glue)
library(ggplot2)
library(ggrepel)
library(picante)
library(broom)
library(tidyverse)
library(RColorBrewer)


set.seed(500)



## to create maps, read these in and proceed to 'Plotting maps' section
tdwg_l3 = vect(file.path(raw.dir, "SpatialData/TDWG/level3/level3.shp"))
site_species_orig = read_rds(file.path(data.dir,'Defaunation/site_species_orig.rds'))
all_traits = read_rds(file.path(data.dir,'Defaunation/all_traits.rds'))
world_grid_vect_filled = read_rds(file.path(data.dir, 'Defaunation/world_grid_vect_filled.rds'))
world_grid_vect_filled = vect(world_grid_vect_filled)
data_final = read_rds(file.path(results.dir,'Defaunation/data_final.rds'))



## Creating original site x species matrix======================================
# read in tdwg level 3 shapefile
tdwg_l3 = vect(file.path(raw.dir, "SpatialData/TDWG/level3/level3.shp"))

# read in ranges of all frugivorous birds
# frug_birds is already filtered from BOTW.gdb for birds that are frugivorous, 
# presence and seasonality are 1,2,3 while origin is 1,2 and all spatial features with
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

# saving world_grid_vect_filled
# saveRDS(world_grid_vect_filled, file = file.path(data.dir,'Defaunation/world_grid_vect_filled.rds'))

# saving site_species_orig
# saveRDS(site_species_orig, file = file.path(data.dir,'Defaunation/site_species_orig.rds'))

# saving all_traits
all_traits = values(frug_birds)
# saveRDS(all_traits, file = file.path(data.dir,'Defaunation/all_traits.rds'))





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



################################################################################
# using cmdscale function for pcoa
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

# save graphs
# ggsave(filename = "pcoa_graph_defaun.png", plot = pcoa_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "pe_graph_defaun.png", plot = pe_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
################################################################################


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


## species richness

# creating breaks in species richness data for graph plotting
grid_fortified$sr_orig_breaks = cut(grid_fortified$sr_orig, breaks = c(-Inf,5,10,25,50,100,150,200,Inf),
                                    labels = c("≤ 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","151 - 200","> 200"))
grid_fortified$sr_HL_breaks = cut(grid_fortified$sr_HL, breaks = c(-Inf,5,10,25,50,100,150,200,Inf),
                                  labels = c("≤ 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","151 - 200","> 200"))
grid_fortified$sr_WT_breaks = cut(grid_fortified$sr_WT, breaks = c(-Inf,5,10,25,50,100,150,200,Inf),
                                  labels = c("≤ 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","151 - 200","> 200"))
grid_fortified$sr_HL_WT_breaks = cut(grid_fortified$sr_HL_WT, breaks = c(-Inf,5,10,25,50,100,150,200,Inf),
                                     labels = c("≤ 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","151 - 200","> 200"))

brewer.pal(8, "BuGn")
sr_colour = c("#F7FCFD", "#E5F5F9", "#CCECE6", "#99D8C9", "#66C2A4",
                            "#41AE76", "#238B45", "#005824")

# original species richness
orig_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_orig_breaks), color = "grey") +
  scale_fill_manual(values = sr_colour, name = "Species Richness", 
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Pre-Defaunation") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
  coord_map(projection = "mollweide") # mollweide projection has issues

orig_sr_map


# HL species richness
HL_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_HL_breaks), color = "grey") +
  scale_fill_manual(values = sr_colour, name = "Species Richness",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HL_sr_map


# WT species richness
WT_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_WT_breaks), color = "grey") +
  scale_fill_manual(values = sr_colour, name = "Species Richness",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

WT_sr_map


# HL_WT species richness
HL_WT_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_HL_WT_breaks), color = "grey") +
  scale_fill_manual(values = sr_colour, name = "Species Richness",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss + Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HL_WT_sr_map



## species lost
grid_fortified$HL_sp_lost = grid_fortified$sr_orig - grid_fortified$sr_HL
grid_fortified$WT_sp_lost = grid_fortified$sr_orig - grid_fortified$sr_WT
grid_fortified$HL_WT_sp_lost = grid_fortified$sr_orig - grid_fortified$sr_HL_WT

grid_fortified$HL_sp_lost_breaks = cut(grid_fortified$HL_sp_lost, breaks = c(-Inf,0,5,10,15,20,Inf),
                                       labels = c("0","1 - 5","6 - 10","10 - 15","15 - 20","> 20"))
grid_fortified$WT_sp_lost_breaks = cut(grid_fortified$WT_sp_lost, breaks = c(-Inf,0,5,10,15,20,Inf),
                                       labels = c("0","1 - 5","6 - 10","10 - 15","15 - 20","> 20"))
grid_fortified$HL_WT_sp_lost_breaks = cut(grid_fortified$HL_WT_sp_lost, breaks = c(-Inf,0,5,10,15,20,Inf),
                                       labels = c("0","1 - 5","6 - 10","10 - 15","15 - 20","> 20"))
brewer.pal(6, "Reds")
sp_lost_colour = c("#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15")

# HL
HL_sp_lost_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = HL_sp_lost_breaks), color = "grey") +
  scale_fill_manual(values = sp_lost_colour, name = "Number of Species",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Species Lost (Habitat Loss)") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HL_sp_lost_map


# WT
WT_sp_lost_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = WT_sp_lost_breaks), color = "grey") +
  scale_fill_manual(values = sp_lost_colour, name = "Number of Species Lost",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Species Loss (Wildlife Trade)") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

WT_sp_lost_map

# HL + WT
HL_WT_sp_lost_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = HL_WT_sp_lost_breaks), color = "grey") +
  scale_fill_manual(values = sp_lost_colour, name = "Number of Species Lost",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Species Loss (Habitat Loss + Wildlife Trade)") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HL_WT_sp_lost_map



## functional dispersion
grid_fortified$fdisp_orig_breaks = cut(grid_fortified$fdisp_orig, breaks = c(-Inf,0.5,1,1.5,2,2.5,3,Inf),
                                       labels = c("0 - 0.5","0.6 - 1","1.1 - 1.5","1.6 - 2","2.1 - 2.5","2.6 - 3","> 3"))

orig_fdisp_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = fdisp_orig_breaks), color = "grey") +
  scale_fill_brewer(palette = "Blues", name = "Functional Dispersion",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Pre-Defaunation Functional Disperson") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
  
orig_fdisp_map

# boxplots for functional dispersion values and z-scores
fdisp_boxplot_data = grid_fortified[,c('fdisp_orig','fdisp_HL', 'fdisp_WT','fdisp_HL_WT')]
fdisp_boxplot_data_melt = melt(fdisp_boxplot_data, na.rm = T)
fdisp_boxplot = ggplot(data = fdisp_boxplot_data_melt, aes(x=variable, y=value)) +
  geom_boxplot() +
  labs(title = "Functional Disperson", x = 'Scenario', y = 'Functional Dispersion') +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("Pre-Defaunation","HL", "WT","HL+WT"))
fdisp_boxplot


z_boxplot_data = grid_fortified[,c('z_HL','z_WT','z_HL_WT')]
z_boxplot_data_melt = melt(z_boxplot_data, na.rm = T)
z_boxplot = ggplot(data = z_boxplot_data_melt, aes(x=variable, y=value)) +
  geom_boxplot()+
  labs(title = "Standardised Effect Size", x = 'Scenario', y = 'SES') +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("HL", "WT","HL+WT"))
z_boxplot


## z-scores
# creating breaks for z-scores
grid_fortified$z_HL_breaks = cut(grid_fortified$z_HL, breaks = c(-Inf,-4,-3,-2,-1,0,1,2,Inf),
                                 labels = c("≤ -4","≤ -3","≤ -2","≤ -1","≤ 0","≤ 1","≤ 2","> 2"))
grid_fortified$z_WT_breaks = cut(grid_fortified$z_WT, breaks = c(-Inf,-4,-3,-2,-1,0,1,2,Inf),
                                 labels = c("≤ -4","≤ -3","≤ -2","≤ -1","≤ 0","≤ 1","≤ 2","> 2"))
grid_fortified$z_HL_WT_breaks = cut(grid_fortified$z_HL_WT, breaks = c(-Inf,-4,-3,-2,-1,0,1,2,Inf),
                                    labels = c("≤ -4","≤ -3","≤ -2","≤ -1","≤ 0","≤ 1","≤ 2","> 2"))

z_colour = c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63')

# z_HL
HL_z_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_HL_breaks), color = "grey") +
  scale_fill_manual(values = z_colour, name = "SES", na.value = 'white',
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HL_z_map


# z_WT
WT_z_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_WT_breaks), color = "grey") +
  scale_fill_manual(values = z_colour, name = "SES", na.value = 'white',
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

WT_z_map


# z_HL_WT
HL_WT_z_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_HL_WT_breaks), color = "grey") +
  scale_fill_manual(values = z_colour, name = "SES", na.value = 'white',
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l3_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss + Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HL_WT_z_map


# experimenting with cowplot
# library(cowplot)
# x <- plot_grid(HL_WT_z_map, HL_WT_z_map, nrow = 2, ncol = 1, labels = "auto")
# ggsave(x, filename = "nhffjyj.pdf", height = 6, width = 3, path = figures.dir)


# saving maps
# ggsave(filename = "sr_orig_defaun.png", plot = orig_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sr_HL_defaun.png", plot = HL_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sr_WT_defaun.png", plot = WT_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sr_HL_WT_defaun.png", plot = HL_WT_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "sp_lost_HL_defaun.png", plot = HL_sp_lost_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sp_lost_WT_defaun.png", plot = WT_sp_lost_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sp_lost_HL_WT_defaun.png", plot = HL_WT_sp_lost_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "fdisp_orig_defaun.png", plot = orig_fdisp_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_HL_defaun.png", plot = HL_z_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_WT_defaun.png", plot = WT_z_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_HL_WT_defaun.png", plot = HL_WT_z_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")

# ggsave(filename = "fdisp_boxplot.png", plot = fdisp_boxplot, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_boxplot.png", plot = z_boxplot, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
