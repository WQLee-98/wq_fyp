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
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(cowplot)


# reading in spatial files and results
tdwg_l1_fortified = readRDS(file.path(data.dir,'Defaunation/tdwg_l1_fortified.rds'))
grid_fortified = readRDS(file.path(results.dir,'Defaunation/grid_fortified.rds'))
all_traits = readRDS(file.path(data.dir,'Defaunation/all_traits_1degree.rds'))


## pcoa graph ==================================================================
pcoa_graph = readRDS(file.path(results.dir,'Defaunation/pcoa_graph.rds'))
pe_graph = readRDS(file.path(results.dir,'Defaunation/pe_graph.rds'))
pcoa_graph
pe_graph

## morphology density plots ====================================================
dens_mass = readRDS(file.path(results.dir,'Defaunation/dens_mass.rds'))
dens_bl = readRDS(file.path(results.dir,'Defaunation/dens_bl.rds'))
dens_bw = readRDS(file.path(results.dir,'Defaunation/dens_bw.rds'))
dens_kipps = readRDS(file.path(results.dir,'Defaunation/dens_kipps.rds'))

dens_mass
dens_bl
dens_bw
dens_kipps

## species richness ============================================================

# creating breaks in species richness data for graph plotting
grid_fortified$sr_orig_breaks = cut(grid_fortified$sr_orig, breaks = c(-Inf,1,5,10,25,50,100,150,Inf),
                                    labels = c("1","2 - 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","> 150"))
grid_fortified$sr_HL_breaks = cut(grid_fortified$sr_HL, breaks = c(-Inf,1,5,10,25,50,100,150,Inf),
                                  labels = c("1","2 - 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","> 150"))
grid_fortified$sr_WT_breaks = cut(grid_fortified$sr_WT, breaks = c(-Inf,1,5,10,25,50,100,150,Inf),
                                  labels = c("1","2 - 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","> 150"))
grid_fortified$sr_HL_WT_breaks = cut(grid_fortified$sr_HL_WT, breaks = c(-Inf,1,5,10,25,50,100,150,Inf),
                                     labels = c("1","2 - 5","6 - 10","11 - 25","26 - 50","51 - 100","101 - 150","> 150"))

brewer.pal(9, "BuGn")
sr_colour = c("#E5F5F9", "#CCECE6", "#99D8C9", "#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B")
                       
# original species richness
orig_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_orig_breaks)) +
  scale_fill_manual(values = sr_colour, name = "Species Richness", 
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Pre-Defaunation") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

orig_sr_map


# HL species richness
HL_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_HL_breaks)) +
  scale_fill_manual(values = sr_colour, name = "Species Richness",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_sr_map


# WT species richness
WT_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_WT_breaks)) +
  scale_fill_manual(values = sr_colour, name = "Species Richness",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

WT_sr_map


# HL_WT species richness
HL_WT_sr_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = sr_HL_WT_breaks)) +
  scale_fill_manual(values = sr_colour, name = "Species Richness",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss & Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_WT_sr_map



## species lost ================================================================
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
sp_lost_colour = c("grey90", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15")

# HL
HL_sp_lost_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = HL_sp_lost_breaks)) +
  scale_fill_manual(values = sp_lost_colour, name = "Number of Species with Threat Magnitude > 5",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), 
        plot.title = element_text(hjust = 0.5))  + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_sp_lost_map


# WT
WT_sp_lost_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = WT_sp_lost_breaks)) +
  scale_fill_manual(values = sp_lost_colour, name = "Number of Species with Threat Magnitude > 5",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), 
        plot.title = element_text(hjust = 0.5))  + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

WT_sp_lost_map

# HL|WT
HL_WT_sp_lost_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = HL_WT_sp_lost_breaks)) +
  scale_fill_manual(values = sp_lost_colour, name = "Number of Species with Threat Magnitude > 5",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss & Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), 
        plot.title = element_text(hjust = 0.5))  + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_WT_sp_lost_map



## functional dispersion =======================================================
grid_fortified$fdisp_orig_breaks = cut(grid_fortified$fdisp_orig, breaks = c(-Inf,0,0.5,1,1.5,2,2.5,3,Inf),
                                       labels = c("0","0.1 - 0.5","0.6 - 1","1.1 - 1.5","1.6 - 2","2.1 - 2.5","2.6 - 3","> 3"))
grid_fortified$fdisp_HL_breaks = cut(grid_fortified$fdisp_HL, breaks = c(-Inf,0,0.5,1,1.5,2,2.5,3,Inf),
                                     labels = c("0","0.1 - 0.5","0.6 - 1","1.1 - 1.5","1.6 - 2","2.1 - 2.5","2.6 - 3","> 3"))
grid_fortified$fdisp_WT_breaks = cut(grid_fortified$fdisp_WT, breaks = c(-Inf,0,0.5,1,1.5,2,2.5,3,Inf),
                                     labels = c("0","0.1 - 0.5","0.6 - 1","1.1 - 1.5","1.6 - 2","2.1 - 2.5","2.6 - 3","> 3"))
grid_fortified$fdisp_HL_WT_breaks = cut(grid_fortified$fdisp_HL_WT, breaks = c(-Inf,0,0.5,1,1.5,2,2.5,3,Inf),
                                        labels = c("0","0.1 - 0.5","0.6 - 1","1.1 - 1.5","1.6 - 2","2.1 - 2.5","2.6 - 3","> 3"))
brewer.pal(9,"Blues")
fd_colour = c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")

# original fdisp
orig_fdisp_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = fdisp_orig_breaks)) +
  scale_fill_manual(values = fd_colour, name = "Functional Dispersion", na.value = "grey80",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Pre-Defaunation") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

orig_fdisp_map


# HL fdisp
HL_fdisp_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = fdisp_HL_breaks)) +
  scale_fill_manual(values = fd_colour, name = "Functional Dispersion", na.value = "grey80",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_fdisp_map


# WT fdisp
WT_fdisp_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = fdisp_WT_breaks)) +
  scale_fill_manual(values = fd_colour, name = "Functional Dispersion", na.value = "grey80",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

WT_fdisp_map

# HL & WT fdisp
HL_WT_fdisp_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = fdisp_HL_WT_breaks)) +
  scale_fill_manual(values = fd_colour, name = "Functional Dispersion", na.value = "grey80",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3.5,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss & Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_WT_fdisp_map



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
  scale_x_discrete(labels = c("Pre-Defaunation","HL", "WT","HL&WT"))
fdisp_boxplot


z_boxplot_data = grid_fortified[,c('z_HL','z_WT','z_HL_WT')]
z_boxplot_data_melt = melt(z_boxplot_data, na.rm = T)
z_boxplot = ggplot(data = z_boxplot_data_melt, aes(x=variable, y=value)) +
  geom_boxplot()+
  labs(title = "Standard Effect Size", x = 'Scenario', y = 'SES') +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("HL", "WT","HL&WT"))
z_boxplot


## z-scores=====================================================================
# creating breaks for z-scores
grid_fortified$z_HL_breaks = cut(grid_fortified$z_HL, breaks = c(-Inf,-4,-3,-2,-1,0,1,2,Inf),
                                 labels = c("≤ -4","≤ -3","≤ -2","≤ -1","≤ 0","≤ 1","≤ 2","> 2"))
grid_fortified$z_WT_breaks = cut(grid_fortified$z_WT, breaks = c(-Inf,-4,-3,-2,-1,0,1,2,Inf),
                                 labels = c("≤ -4","≤ -3","≤ -2","≤ -1","≤ 0","≤ 1","≤ 2","> 2"))
grid_fortified$z_HL_WT_breaks = cut(grid_fortified$z_HL_WT, breaks = c(-Inf,-4,-3,-2,-1,0,1,2,Inf),
                                    labels = c("≤ -4","≤ -3","≤ -2","≤ -1","≤ 0","≤ 1","≤ 2","> 2"))
brewer.pal(10, "RdBu")
z_colour = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#92C5DE", "#2166AC", "#053061")

# z_HL
HL_z_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_HL_breaks)) +
  scale_fill_manual(values = z_colour, name = "Functional Dispersion Lost", na.value = 'grey80',
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_z_map


# z_WT
WT_z_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_WT_breaks)) +
  scale_fill_manual(values = z_colour, name = "Functional Dispersion Lost", na.value = 'grey80',
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

WT_z_map


# z_HL_WT
HL_WT_z_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = z_HL_WT_breaks)) +
  scale_fill_manual(values = z_colour, name = "Functional Dispersion Lost", na.value = 'grey80',
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss & Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

HL_WT_z_map



## mean threat magnitudes ======================================================
# creating breaks for mean threat magnitudes
grid_fortified$mean_HL_breaks = cut(grid_fortified$mean_HL, breaks = c(-Inf,0,1,2,3,4,5,6,Inf),
                                       labels = c("0","≤ 1","≤ 2","≤ 3","≤ 4","≤ 5","≤ 6","≤ 7"))
grid_fortified$mean_WT_breaks = cut(grid_fortified$mean_WT, breaks = c(-Inf,0,1,2,3,4,5,6,Inf),
                                       labels = c("0","≤ 1","≤ 2","≤ 3","≤ 4","≤ 5","≤ 6","≤ 7"))

brewer.pal(9,"Purples")
mean_threat_colour = c("#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D")

# mean_HL
mean_HL_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = mean_HL_breaks)) +
  scale_fill_manual(values = mean_threat_colour, name = "Mean Threat Magnitude",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Habitat Loss") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

mean_HL_map

mean_WT_map = ggplot() +
  geom_polygon(data = grid_fortified, aes(x = long, y = lat, group = group, fill = mean_WT_breaks)) +
  scale_fill_manual(values = mean_threat_colour, name = "Mean Threat Magnitude",
                    guide = guide_legend(direction = 'horizontal',
                                         title.position = 'top',
                                         title.hjust = 0.5,
                                         label.position = 'bottom',
                                         label.hjust = 0.5,
                                         nrow = 1,
                                         keywidth = 3,
                                         keyheight = 1)) +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  labs(title = "Wildlife Trade") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

mean_WT_map




## creating figures ============================================================
# species lost for all scenarios [WIP]
sp_lost_legend = get_legend(HL_WT_sp_lost_map + theme(legend.position = "bottom")
)

sp_lost_maps = plot_grid(HL_sp_lost_map + theme(legend.position = "none"), 
                         WT_sp_lost_map + theme(legend.position = "none"), 
                         HL_WT_sp_lost_map + theme(legend.position = "none"),
                         sp_lost_legend, 
                         ncol = 1, labels = c("A","B","C",""),
                         rel_heights = c(1,1,1,0.1),
                         align = "v")
  

## saving maps =================================================================
save_plot(file.path(figures.dir, "sp_lost_maps_test.png"),sp_lost_maps,
          base_height = 8, base_width = 4,bg="white")


# ggsave(filename = "pcoa_graph_defaun.png", plot = pcoa_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "pe_graph_defaun.png", plot = pe_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")


# ggsave(filename = "sr_orig_defaun_1degree.png", plot = orig_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sr_HL_defaun_1degree.png", plot = HL_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sr_WT_defaun_1degree.png", plot = WT_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sr_HL_WT_defaun_1degree.png", plot = HL_WT_sr_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "sp_lost_HL_defaun_1degree.png", plot = HL_sp_lost_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sp_lost_WT_defaun_1degree.png", plot = WT_sp_lost_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "sp_lost_HL_WT_defaun_1degree.png", plot = HL_WT_sp_lost_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "fdisp_orig_defaun_1degree.png", plot = orig_fdisp_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "fdisp_HL_defaun_1degree.png", plot = HL_fdisp_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "fdisp_WT_defaun_1degree.png", plot = WT_fdisp_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "fdisp_HL_WT_defaun_1degree.png", plot = HL_WT_fdisp_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "z_HL_defaun_1degree.png", plot = HL_z_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_WT_defaun_1degree.png", plot = WT_z_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_HL_WT_defaun_1degree.png", plot = HL_WT_z_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "fdisp_boxplot_1degree.png", plot = fdisp_boxplot, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "z_boxplot_1degree.png", plot = z_boxplot, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# 
# ggsave(filename = "mean_HL_map_1degree.png", plot = mean_HL_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
# ggsave(filename = "mean_WT_map_1degree.png", plot = mean_WT_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
