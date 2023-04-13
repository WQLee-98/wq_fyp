rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Documents/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

library(dplyr)
library(reshape2)
library(purrr)
library(nlme)
library(MuMIn)
library(bipartite)
library(glue)
library(ggplot2)
library(ggrepel)
library(letsR)
library(terra)
library(sf)

## Data preparation ============================================================
# read Fricke's network data (without unidentified taxa)
network = read.csv(file.path(raw.dir, "Fricke/network_long.csv"), header = T)
# filter birds
network_bird = network[network$animal.class == 'Aves',]
length(unique(network_bird$animal.id)) # 1268 bird species

birds_network = unique(network_bird$animal.id)

# read AVONET trait data
AVONET_traits = read.csv(file.path(raw.dir, "AVONET1_Birdlife.csv"), header = T)
birds_AVONET = unique(AVONET_traits$Species1) 

sum(birds_network %in% birds_AVONET) # only 1125 species out of 1268 species in network data are in the list of AVONET frugivores

# subset network data to bird species that are in AVONET
network_AVONET = network_bird[network_bird$animal.id %in% birds_AVONET,]


# subset to interactions only (value == 1)
network_int = network_AVONET[network_AVONET$value == 1,] # 258 networks, 16177 interactions, 1086 bird species and 1911 plant species
network_no_int = network_AVONET[network_AVONET$value == 0,]
length(unique(network_no_int$animal.id)[which(!(unique(network_no_int$animal.id) %in% unique(network_int$animal.id)))]) # 39 recorded bird species with no interactions


# data exploration
prop.table(table(network_int$animal.native.status))*100 # 15697 native bird interactions vs 480 introduced bird interactions
length(unique(network_int[network_int$animal.native.status == "native",]$animal.id)) # 1070 native species
length(unique(network_int[network_int$animal.native.status == "introduced",]$animal.id)) # 31 introduced species


# all interactions split by network into a list
network_int_list = split(network_int, f = network_int$net.id)



# creating list of unique bird species in each network
network_unique_bird_list = lapply(network_int_list, FUN = function(x){
  bird_data = x[,c('animal.id', 'animal.native.status','latitude', 'longitude', 'net.id', 'locality.id', 'locality', 'region')]
  bird_data_no_dup = bird_data[!duplicated(bird_data),]
  return(bird_data_no_dup)
})

# creating list of unique plant species in each network
network_unique_plant_list = lapply(network_int_list, FUN = function(x){
  plant_data = x[,c('plant.id', 'plant.native.status','latitude', 'longitude', 'net.id', 'locality.id', 'locality', 'region')]
  plant_data_no_dup = plant_data[!duplicated(plant_data),]
  return(plant_data_no_dup)
})

# creating dataframe with percentage of native animal species and number of animal species in each network
network_unique_bird_info = lapply(network_unique_bird_list, FUN = function(x){
  net.id = x[1,'net.id']
  native_bird = prop.table(table(x$animal.native.status))*100
  native_bird = native_bird['native']
  nb.bird.sp = nrow(x)
  return(cbind(net.id,native_bird,nb.bird.sp))
}) %>%
  do.call("rbind",.)

# creating dataframe with percentage of native plant species and number of plant species in each network
network_unique_plant_info = lapply(network_unique_plant_list, FUN = function(x){
  native_plant = prop.table(table(x$plant.native.status))*100
  native_plant = native_plant['native']
  nb.plant.sp = nrow(x)
  return(cbind(native_plant,nb.plant.sp))
}) %>%
  do.call("rbind",.)

network_species_info = as.data.frame(cbind(network_unique_bird_info,network_unique_plant_info))
rownames(network_species_info) = 1:nrow(network_species_info)
network_species_info[,2:5] = sapply(network_species_info[,2:5],as.numeric)

# networks with >3 species and predominantly native bird and plant species
selected = network_species_info[network_species_info$native_bird > 50 & network_species_info$nb.bird.sp > 3
                                & network_species_info$native_plant > 50 & network_species_info$nb.plant.sp > 3,]
selected = selected[!is.na(selected$net.id),]




# subsetting to networks with >3 species and predominantly native bird and plant species
network_int_final = network_int[network_int$net.id %in% selected$net.id,] # 15263 interactions, 210 networks, 1041 bird species, 1824 plant species
network_int_list_final = network_int_list[selected$net.id]
network_unique_bird_list_final = network_unique_bird_list[selected$net.id]
network_unique_plant_list_final = network_unique_plant_list[selected$net.id]


## Functional Diversity ========================================================
# creating site x species matrix
site_species = acast(network_int_final, formula = net.id~animal.id, fill = 0, 
                     fun.aggregate = length, value.var = "animal.id")
site_species = ifelse(site_species > 0, 1, 0)

# creating species x traits matrix
rownames(AVONET_traits) = AVONET_traits$Species1
AVONET_traits$Beak.Length_Culmen = as.numeric(AVONET_traits$Beak.Length_Culmen) # length of culmen
AVONET_traits$Mass = as.numeric(AVONET_traits$Mass)
AVONET_traits$Beak.Width = as.numeric(AVONET_traits$Beak.Width)
AVONET_traits$Kipps.Distance = as.numeric(AVONET_traits$Kipps.Distance)
AVONET_traits$Wing.Length = as.numeric(AVONET_traits$Wing.Length)
AVONET_traits$Kipps.Index = AVONET_traits$Kipps.Distance / AVONET_traits$Wing.Length

# subset species to those in network_int_final
network_bird_traits = AVONET_traits[AVONET_traits$Species1 %in% unique(network_int_final$animal.id),]
network_bird_traits = network_bird_traits[,c("Mass","Beak.Length_Culmen",
                                             "Beak.Width","Kipps.Index")]
network_bird_traits = network_bird_traits[match(colnames(site_species),
                                                rownames(network_bird_traits)),]


## functional dispersion
source(file.path(main.dir,"src/fdFunc2.R"))

# pcoa on global scale
fd_res_global = funcdiv_global(x = network_bird_traits, a = site_species, original = TRUE)
dist.cent.global = fd_res_global$dist.cent

# pcoa on network scale
fd_res_network = funcdiv_network(x = network_bird_traits, a = site_species, original = TRUE)
dist.cent.network = fd_res_network$dist.cent

# create new list with animal.id, net.id, and dist.cent
network_dist = list()
for (i in 1:length(network_unique_bird_list_final)){
  network_dist[[i]] = network_unique_bird_list_final[[i]]
  network_dist[[i]]$dist.cent.global = dist.cent.global[[i]][match(network_unique_bird_list_final[[i]]$animal.id, names(dist.cent.global[[i]]))]
  network_dist[[i]]$dist.cent.network = dist.cent.network[[i]][match(network_unique_bird_list_final[[i]]$animal.id, names(dist.cent.network[[i]]))]
}


## Seed Dispersal Reliance ======================================================
# create bird x plant matrix for each network
bird_plant = lapply(network_int_list_final, FUN = acast,formula = animal.id~plant.id, 
                    fill = 0, fun.aggregate = length, value.var = "plant.id") %>%
  lapply(., FUN = function(x){
    as.matrix((x>0)+0)
  })

# create function that calculates the Seed Dispersal Reliance of each species in each network
int_uniq = function(x){
  # Arguments
  #   x: bird x plant matrix for a network
  unique_int = c()
  for (i in 1:nrow(x)){
    # average number of dispersers per plant for the plants that species i interacts with
    sp_uniq = 1/(sum(colSums(x)*x[i,]))/sum(x[i,])
    unique_int = c(unique_int,sp_uniq)
  }
  names(unique_int) = rownames(x)
  return(unique_int)
}

bird_uniq = lapply(bird_plant, FUN = int_uniq)

# merging bird_uniq with network_dist
bird_data = list()
for (i in 1:length(network_dist)){
  bird_data[[i]] = network_dist[[i]]
  bird_data[[i]]$int.uniq = bird_uniq[[i]][match(bird_data[[i]]$animal.id, names(bird_uniq[[i]]))]
}

# normalising data
for (i in 1:length(bird_data)){
  bird_data[[i]]$dist.cent.global.norm = scale(bird_data[[i]]$dist.cent.global)
  bird_data[[i]]$dist.cent.network.norm = scale(bird_data[[i]]$dist.cent.network)
  bird_data[[i]]$int.uniq.norm = scale(bird_data[[i]]$int.uniq)
}


# merge all dataframes into one
bird_data_final = do.call("rbind",bird_data) # 3705 data points


range(bird_data_final$dist.cent.global)
range(bird_data_final$dist.cent.network)
range(bird_data_final$int.uniq)
range(bird_data_final$dist.cent.global.norm)
range(bird_data_final$dist.cent.network.norm)
range(bird_data_final$int.uniq.norm)



## Model fitting (linear mixed effects model) ==================================
# making net.id a factor
bird_data_final$net.id = factor(bird_data_final$net.id)

# saving bird_data_final
# saveRDS(bird_data_final, file = file.path(results.dir,'Uniqueness/bird_data_final.rds'))
# bird_data_final = readRDS(file.path(results.dir,'Uniqueness/bird_data_final.rds'))


plot(int.uniq~dist.cent.global, data = bird_data_final)
plot(int.uniq~dist.cent.network, data = bird_data_final)

plot(int.uniq.norm~dist.cent.network.norm, data = bird_data_final)
plot(int.uniq.norm~dist.cent.global.norm, data = bird_data_final)


### network pcoa =========================
## random effects model ===
# random intercept model
network_rand_int_model = lme(int.uniq.norm ~ 1, random = ~1|net.id, data = bird_data_final)
summary(network_rand_int_model)
r.squaredGLMM(network_rand_int_model)
# marginal: variance explained by fixed factors
# conditional: variance explained by fixed and random factors (entire model)


plot(network_rand_int_model) # check heteroscedasticity in whole model
plot(network_rand_int_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
plot(network_rand_int_model, net.id~resid(., scaled=TRUE))
plot(network_rand_int_model,resid(.,scaled=TRUE)~dist.cent.network.norm|net.id, abline = 0)
qqnorm(network_rand_int_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
qqnorm(network_rand_int_model,~ranef(.)) # check normality of random effects
plot(ranef(network_rand_int_model)) # random effects



plot(int.uniq.norm~dist.cent.network.norm, data = bird_data_final)
network_rand_int_model_coef = network_rand_int_model$coefficients$random$net.id
for (i in 1:length(network_rand_int_model_coef)){
  abline(a = network_rand_int_model_coef[i], b = 0, col = "blue")
}
abline(a = network_rand_int_model$coefficients$fixed, b = 0, col = "red")



# random slope model
network_rand_slope_model = lme(int.uniq.norm ~ 1, random = ~1 + dist.cent.network.norm|net.id, data = bird_data_final)
summary(network_rand_slope_model)
r.squaredGLMM(network_rand_slope_model)


plot(network_rand_slope_model) # check heteroscedasticity in whole model
plot(network_rand_slope_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
plot(network_rand_slope_model, net.id~resid(., scaled=TRUE))
plot(network_rand_slope_model,resid(.,scaled=TRUE)~dist.cent.network.norm|net.id, abline = 0)
qqnorm(network_rand_slope_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
qqnorm(network_rand_slope_model,~ranef(.)) # check normality of random effects
plot(ranef(network_rand_slope_model)) # random effects


plot(int.uniq.norm~dist.cent.network.norm, data = bird_data_final)
network_rand_slope_model_coef = network_rand_slope_model$coefficients$random$net.id
for (i in 1:nrow(network_rand_slope_model_coef)){
  abline(coef = network_rand_slope_model_coef[i,], col = "blue")
}
abline(a = network_rand_slope_model$coefficients$fixed, b = 0, col = "red")


## fixed + random effects model===
# random intercept model
network_fixed_rand_int_model = lme(int.uniq.norm~dist.cent.network.norm, random = ~1|net.id, data = bird_data_final)
summary(network_fixed_rand_int_model)
r.squaredGLMM(network_fixed_rand_int_model)


plot(network_fixed_rand_int_model) # check heteroscedasticity in whole model
plot(network_fixed_rand_int_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
plot(network_fixed_rand_int_model, net.id~resid(., scaled=TRUE))
plot(network_fixed_rand_int_model,resid(.,scaled=TRUE)~dist.cent.network.norm|net.id, abline = 0)
qqnorm(network_fixed_rand_int_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
qqnorm(network_fixed_rand_int_model,~ranef(.)) # check normality of random effects
plot(ranef(network_fixed_rand_int_model)) # random effects



plot(int.uniq.norm~dist.cent.network.norm, data = bird_data_final)
network_fixed_rand_int_model_coef = network_fixed_rand_int_model$coefficients$random$net.id
for (i in 1:length(network_fixed_rand_int_model_coef)){
  abline(a = network_fixed_rand_int_model_coef[i], b = network_fixed_rand_int_model$coefficients$fixed[2], col = "blue")
}
abline(coef = network_fixed_rand_int_model$coefficients$fixed, col = "red")


# random slope model
network_fixed_rand_slope_model = lme(int.uniq.norm~dist.cent.network.norm, random = ~1 + dist.cent.network.norm|net.id, data = bird_data_final)
summary(network_fixed_rand_slope_model)
r.squaredGLMM(network_fixed_rand_slope_model)

# test <- coefficients(network_fixed_rand_slope_model)
# str(network_fixed_rand_slope_model)
# 
# ggplot() + geom_point(aes(y = int.uniq.norm, x = dist.cent.network.norm), 
#                       data = bird_data_final) +
#   geom_abline(aes(intercept = 0.2187702, slope = 0.15), col = "red") +
#   geom_abline(aes(intercept = `(Intercept)`, slope = dist.cent.network.norm),
#               data =test )

plot(network_fixed_rand_slope_model) # check heteroscedasticity in whole model
plot(network_fixed_rand_slope_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
plot(network_fixed_rand_slope_model, net.id~resid(., scaled=TRUE))
plot(network_fixed_rand_slope_model,resid(.,scaled=TRUE)~dist.cent.network.norm|net.id, abline = 0)
qqnorm(network_fixed_rand_slope_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
qqnorm(network_fixed_rand_slope_model,~ranef(.)) # check normality of random effects
plot(ranef(network_fixed_rand_slope_model)) # random effects

# exploring random effect intercept and slope
library(lme4)
library(lattice)
refit_network_fixed_rand_slope_model = lmer(int.uniq.norm~dist.cent.network.norm + 
                                              (dist.cent.network.norm|net.id), data = bird_data_final)
dotplot(ranef(refit_network_fixed_rand_slope_model, condVar = T))



plot(int.uniq.norm~dist.cent.network.norm, data = bird_data_final)
network_fixed_rand_slope_model_coef = network_fixed_rand_slope_model$coefficients$random$net.id
for (i in 1:nrow(network_fixed_rand_slope_model_coef)){
  abline(coef = network_fixed_rand_slope_model_coef[i,], col = "blue")
}
abline(coef = network_fixed_rand_slope_model$coefficients$fixed, col = "red")



model_sel = model.sel(network_rand_int_model,network_rand_slope_model,
                       network_fixed_rand_int_model,network_fixed_rand_slope_model, rank = "AIC")


# model simplification
network_fixed_rand_slope_model_REML = lme(int.uniq.norm~dist.cent.network.norm, 
                                     random = ~1 + dist.cent.network.norm|net.id, 
                                     method = "REML",
                                     data = bird_data_final)
network_fixed_rand_int_model_REML = lme(int.uniq.norm~dist.cent.network.norm, 
                                        random = ~1|net.id, 
                                        method = "REML",
                                        data = bird_data_final)
anova(network_fixed_rand_slope_model_REML,network_fixed_rand_int_model_REML)
# maximal model has lower AIC


network_fixed_rand_slope_model_ML = lme(int.uniq.norm~dist.cent.network.norm, 
                                          random = ~1 + dist.cent.network.norm|net.id, 
                                          method = "ML",
                                          data = bird_data_final)
network_rand_slope_model_ML = lme(int.uniq.norm ~ 1, 
                                  random = ~1 + dist.cent.network.norm|net.id, 
                                  method = "ML",
                                  data = bird_data_final)

anova(network_fixed_rand_slope_model_ML, network_rand_slope_model_ML)
# maximal model has lower AIC



# # correlation between latitude and slope
# res_ex = bird_data_final[,c('net.id','latitude')]
# res_ex = res_ex[!duplicated(res_ex),]
# res_ex$intercept = as.data.frame(network_fixed_rand_slope_model_coef)[match(rownames(network_fixed_rand_slope_model_coef),res_ex$net.id),][,1]
# res_ex$slope = as.data.frame(network_fixed_rand_slope_model_coef)[match(rownames(network_fixed_rand_slope_model_coef),res_ex$net.id),][,2]
# lat_model = lm(slope~latitude, data = res_ex)
# summary(lat_model)
# 
# nestedness = lapply(bird_plant, FUN = nested, method = 'NODF2') %>%
#   do.call("rbind",.)
# 
# res_ex$nestedness = as.data.frame(nestedness)[match(rownames(nestedness),res_ex$net.id),]
# nest_model = lm(slope~nestedness, data = res_ex)
# summary(nest_model)
# 
# lat_nest_model = lm(slope~latitude*nestedness, data = res_ex)
# summary(lat_nest_model)








## Plotting graphs =============================================================
all_networks_scatter = ggplot(data = bird_data_final, aes(x = dist.cent.network.norm, y = int.uniq.norm)) +
  geom_point() +
  labs(x = "Morphological Uniqueness", y = "Seed Dispersal Reliance") +
  geom_abline(aes(slope = network_fixed_rand_slope_model$coefficients$fixed[2], 
                  intercept = network_fixed_rand_slope_model$coefficients$fixed[1]), 
              col = 'red', size = 1.5, show.legend = F) +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14))

all_networks_scatter


# ggsave(filename = "all_networks_scatter.png", plot = all_networks_scatter, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")






## PCoA global
# using cmdscale function for pcoa
pcoa_global = cmdscale(fd_res_global$dist.mat, k = 4, eig = T, add = T)

global_positions = as.data.frame(pcoa_global$points)
colnames(global_positions) = c("pcoa1", "pcoa2", "pcoa3", "pcoa4")

# adding taxonomy data
bird_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
bird_tax = bird_tax[bird_tax$Species1 %in% rownames(global_positions),]
bird_tax = bird_tax[match(rownames(global_positions),bird_tax$Species1),]
global_positions = cbind(global_positions, Family = bird_tax$Family1, Order = bird_tax$Order1)

# calculating percentage variance explained by each pcoa axis
percent_explained = (100 * pcoa_global$eig / sum(pcoa_global$eig))
rounded_pe = round(percent_explained[1:4],2)

labs_global = c(glue("PCoA 1 ({rounded_pe[1]}%)"),
                glue("PCoA 2 ({rounded_pe[2]}%)"))

# plot pcoa
global_positions_tibble = as_tibble(global_positions, rownames = "species")
pcoa_global_graph = ggplot(data = global_positions_tibble, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(color = Order)) +
  geom_text_repel(aes(label=species), data = subset(global_positions_tibble, pcoa1<(-7.5))) +
  labs(x = labs_global[1], y = labs_global[2]) + 
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11))
pcoa_global_graph

# ggsave(filename = "pcoa_global.png", plot = pcoa_global_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")


# scale_colour_manual(values = c('#6B6100','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#                                         '#e31a1c','#fdbf6f','#ff7f00','#007E71','#6a3d9a',
#                                         '#FFFF37','#b15928','#2A00FF','#DE0093')) +
#                                           labs(x = labs[1], y = labs[2]) +


# plot percentage explained
pe_global_graph = tibble(pe = cumsum(rounded_pe),
                  axis = 1:length(rounded_pe)) %>%
  ggplot(aes(x=axis,y=pe)) +
  geom_line(size = 1) +
  geom_point()+
  coord_cartesian(xlim = c(1,4), ylim = c(0,100)) +
  labs(x = "PCoA Axis", y = "Cumulative Percentage Explained") + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))
pe_global_graph



slope_coef_ordered = as.data.frame(network_fixed_rand_slope_model_coef) %>%
  arrange(desc(dist.cent.network.norm))


## PCoA network 
# finding max slope
network_fixed_rand_slope_model_coef[which.max(network_fixed_rand_slope_model_coef[,2]),] # Dehling 2017 San Pedro 1 
fd_res_network$FDis["Dehling 2017 San Pedro 1"]

# max slope
pcoa_dehling = cmdscale(fd_res_network$dist.mat$`Dehling 2017 San Pedro 1`, k = 4, eig = T, add = T)
dehling_positions = as.data.frame(pcoa_dehling$points)
colnames(dehling_positions) = c("pcoa1", "pcoa2", "pcoa3", "pcoa4")

# adding taxonomy data
dehling_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
dehling_tax = dehling_tax[dehling_tax$Species1 %in% rownames(dehling_positions),]
dehling_tax = dehling_tax[match(rownames(dehling_positions),dehling_tax$Species1),]
dehling_positions = cbind(dehling_positions, Family = dehling_tax$Family1, Order = dehling_tax$Order1)

# calculating percentage variance explained by each pcoa axis
dehling_percent_explained = (100 * pcoa_dehling$eig / sum(pcoa_dehling$eig))
dehling_rounded_pe = round(dehling_percent_explained[1:4],2)


labs_dehling = c(glue("PCoA 1 ({dehling_rounded_pe[1]}%)"),
                 glue("PCoA 2 ({dehling_rounded_pe[2]}%)"))

# plot pcoa
dehling_positions_tibble = as_tibble(dehling_positions, rownames = "species")
dehling_centroid = apply(pcoa_dehling$points,2,FUN = mean)
pcoa_dehling_graph = ggplot(data = dehling_positions_tibble, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(color = Order)) +
  geom_text_repel(aes(label=ifelse(pcoa1 < (-2) | pcoa2 < (-2), as.character(species),''))) +
  labs(x = labs_dehling[1], y = labs_dehling[2], title = 'Dehling 2017 San Pedro 1 Network') +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        title = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 20))
pcoa_dehling_graph

# ggsave(filename = "pcoa_dehling.png", plot = pcoa_dehling_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")

dehling_scatter_data = subset(bird_data_final, net.id == 'Dehling 2017 San Pedro 1')
dehling_scatter_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
dehling_scatter_tax = dehling_scatter_tax[dehling_scatter_tax$Species1 %in% dehling_scatter_data$animal.id,]
dehling_scatter_tax = dehling_scatter_tax[match(dehling_scatter_data$animal.id,dehling_scatter_tax$Species1),]
dehling_scatter_data = cbind(dehling_scatter_data, Family = dehling_scatter_tax$Family1, Order = dehling_scatter_tax$Order1)


dehling_scatter = ggplot() +
  geom_abline(aes(slope = network_fixed_rand_slope_model$coefficients$random$net.id['Dehling 2017 San Pedro 1',][2], 
                  intercept = network_fixed_rand_slope_model$coefficients$random$net.id['Dehling 2017 San Pedro 1',][1]), 
              color = "red", size = 1.5, show.legend = F) +
  geom_point(data=dehling_scatter_data, aes(x = dist.cent.network.norm, y = int.uniq.norm, color = Order)) +
  labs(x = "Morphological Uniqueness", y = "Seed Dispersal Reliance", title = 'Dehling 2017 San Pedro 1 Network') +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),        
        title = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 20))

dehling_scatter

# ggsave(filename = "dehling_scatter.png", plot = dehling_scatter, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")


# finding min slope
network_fixed_rand_slope_model_coef[which.min(network_fixed_rand_slope_model_coef[,2]),] # Sethi 2012, also has max intercept
fd_res_network$FDis["Sethi 2012"]

# min slope
pcoa_sethi = cmdscale(fd_res_network$dist.mat$`Sethi 2012`, k = 4, eig = T, add = T)
sethi_positions = as.data.frame(pcoa_sethi$points)
colnames(sethi_positions) = c("pcoa1", "pcoa2", "pcoa3", "pcoa4")

# adding taxonomy data
sethi_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
sethi_tax = sethi_tax[sethi_tax$Species1 %in% rownames(sethi_positions),]
sethi_tax = sethi_tax[match(rownames(sethi_positions),sethi_tax$Species1),]
sethi_positions = cbind(sethi_positions, Family = sethi_tax$Family1, Order = sethi_tax$Order1)

# calculating percentage variance explained by each pcoa axis
sethi_percent_explained = (100 * pcoa_sethi$eig / sum(pcoa_sethi$eig))
sethi_rounded_pe = round(sethi_percent_explained[1:4],2)


labs_sethi = c(glue("PCoA 1 ({sethi_rounded_pe[1]}%)"),
                 glue("PCoA 2 ({sethi_rounded_pe[2]}%)"))

# plot pcoa
sethi_positions_tibble = as_tibble(sethi_positions, rownames = "species")
sethi_centroid = apply(pcoa_sethi$points,2,FUN = mean)
pcoa_sethi_graph = ggplot(data = sethi_positions_tibble, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(color = Order)) +
  geom_text_repel(aes(label=ifelse(pcoa1 > 1 | pcoa2 > 1 | pcoa2 < (-1), as.character(species),''))) +
  labs(x = labs_sethi[1], y = labs_sethi[2], title = "Sethi 2012 Network") +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        title = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 20))
pcoa_sethi_graph

# ggsave(filename = "pcoa_sethi.png", plot = pcoa_sethi_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")


sethi_scatter_data = subset(bird_data_final, net.id == 'Sethi 2012')
sethi_scatter_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
sethi_scatter_tax = sethi_scatter_tax[sethi_scatter_tax$Species1 %in% sethi_scatter_data$animal.id,]
sethi_scatter_tax = sethi_scatter_tax[match(sethi_scatter_data$animal.id,sethi_scatter_tax$Species1),]
sethi_scatter_data = cbind(sethi_scatter_data, Family = sethi_scatter_tax$Family1, Order = sethi_scatter_tax$Order1)



sethi_scatter = ggplot(data=sethi_scatter_data, aes(x = dist.cent.network.norm, y = int.uniq.norm, color = Order)) +
  labs(x = "Morphological Uniqueness", y = "Seed Dispersal Reliance", title = "Sethi 2012 Network") +
  geom_abline(aes(slope = network_fixed_rand_slope_model$coefficients$random$net.id['Sethi 2012',][2], 
                  intercept = network_fixed_rand_slope_model$coefficients$random$net.id['Sethi 2012',][1]),
              color = 'red', size = 1.5, show.legend = F) +
  geom_point() +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        title = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 20))

sethi_scatter

# ggsave(filename = "sethi_scatter.png", plot = sethi_scatter, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")




# Quitian 2018 2000Natural
fd_res_network$FDis["Quitian 2018 2000Natural"]

# max slope
pcoa_Quitian = cmdscale(fd_res_network$dist.mat$`Quitian 2018 2000Natural`, k = 4, eig = T, add = T)
Quitian_positions = as.data.frame(pcoa_Quitian$points)
colnames(Quitian_positions) = c("pcoa1", "pcoa2", "pcoa3", "pcoa4")

# adding taxonomy data
Quitian_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
Quitian_tax = Quitian_tax[Quitian_tax$Species1 %in% rownames(Quitian_positions),]
Quitian_tax = Quitian_tax[match(rownames(Quitian_positions),Quitian_tax$Species1),]
Quitian_positions = cbind(Quitian_positions, Family = Quitian_tax$Family1, Order = Quitian_tax$Order1)

# calculating percentage variance explained by each pcoa axis
Quitian_percent_explained = (100 * pcoa_Quitian$eig / sum(pcoa_Quitian$eig))
Quitian_rounded_pe = round(Quitian_percent_explained[1:4],2)


labs_Quitian = c(glue("PCoA 1 ({Quitian_rounded_pe[1]}%)"),
                 glue("PCoA 2 ({Quitian_rounded_pe[2]}%)"))

# plot pcoa
Quitian_positions_tibble = as_tibble(Quitian_positions, rownames = "species")
Quitian_centroid = apply(pcoa_Quitian$points,2,FUN = mean)
pcoa_Quitian_graph = ggplot(data = Quitian_positions_tibble, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(color = Order)) +
  geom_text_repel(aes(label=ifelse(pcoa1 < (-2) | pcoa2 < (-2), as.character(species),''))) +
  labs(x = labs_Quitian[1], y = labs_Quitian[2], title = 'Quitian 2018 2000Natural Network') +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        title = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 20))
pcoa_Quitian_graph

# ggsave(filename = "pcoa_Quitian.png", plot = pcoa_Quitian_graph, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")

Quitian_scatter_data = subset(bird_data_final, net.id == 'Quitian 2018 2000Natural')
Quitian_scatter_tax = AVONET_traits[,c("Species1", "Family1", "Order1")]
Quitian_scatter_tax = Quitian_scatter_tax[Quitian_scatter_tax$Species1 %in% Quitian_scatter_data$animal.id,]
Quitian_scatter_tax = Quitian_scatter_tax[match(Quitian_scatter_data$animal.id,Quitian_scatter_tax$Species1),]
Quitian_scatter_data = cbind(Quitian_scatter_data, Family = Quitian_scatter_tax$Family1, Order = Quitian_scatter_tax$Order1)


Quitian_scatter = ggplot() +
  geom_abline(aes(slope = network_fixed_rand_slope_model$coefficients$random$net.id['Quitian 2018 2000Natural',][2], 
                  intercept = network_fixed_rand_slope_model$coefficients$random$net.id['Quitian 2018 2000Natural',][1]), 
              color = "red", size = 1.5, show.legend = F) +
  geom_point(data=Quitian_scatter_data, aes(x = dist.cent.network.norm, y = int.uniq.norm, color = Order)) +
  labs(x = "Morphological Uniqueness", y = "Seed Dispersal Reliance", title = 'Quitian 2018 2000Natural Network') +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size = 14),        
        title = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 20))

Quitian_scatter

# ggsave(filename = "Quitian_scatter.png", plot = Quitian_scatter, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")



# locations of networks
locations = lapply(network_unique_bird_list_final, FUN = function(x){
  x[1,c("net.id","latitude","longitude")]
}) %>%
  do.call("rbind",.)

tdwg_l1 = vect(file.path(raw.dir, "SpatialData/TDWG/level1/level1_dissolved.shp"))
tdwg_l1_fortified = fortify(as(tdwg_l1, "Spatial"))

network_location_map = ggplot() +
  geom_polygon(data = tdwg_l1_fortified, aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  geom_point(data = locations, aes(x = longitude, y = latitude), shape = 21, fill = '#33a02c', color = 'black', size = 2.5, stroke = 0.8) +
  labs(title = "Network Locations") +
  theme_void() +
  theme(title = element_text(face = 'bold'), legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  coord_map(projection = "mollweide", xlim = c(-180,180), ylim = c(-55,90))

network_location_map


# ggsave(filename = "network_location_map.png", plot = network_location_map, path = figures.dir, width = 25,
#        height = 15, units = "cm", dpi = "retina")
















## Experimenting ===============================================================
# ### global pcoa =========================
# ## random effects model ===
# # random intercept model
# global_rand_int_model = lme(int.uniq.norm ~ 1, random = ~1|net.id, data = bird_data_final)
# summary(global_rand_int_model)
# r.squaredGLMM(global_rand_int_model)
# # marginal: variance explained by fixed factors
# # conditional: variance explained by fixed and random factors (entire model)
# 
# 
# plot(global_rand_int_model) # check heteroscedasticity in whole model
# plot(global_rand_int_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
# plot(global_rand_int_model, net.id~resid(., scaled=TRUE))
# plot(global_rand_int_model,resid(.,scaled=TRUE)~dist.cent.global.norm|net.id, abline = 0)
# qqnorm(global_rand_int_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
# qqnorm(global_rand_int_model,~ranef(.)) # check normality of random effects
# plot(ranef(global_rand_int_model)) # random effects
# 
# 
# 
# plot(int.uniq.norm~dist.cent.global.norm, data = bird_data_final)
# global_rand_int_model_coef = global_rand_int_model$coefficients$random$net.id
# for (i in 1:length(global_rand_int_model_coef)){
#   abline(a = global_rand_int_model_coef[i], b = 0, col = "blue")
# }
# abline(a = global_rand_int_model$coefficients$fixed, b = 0, col = "red")
# 
# 
# 
# # random slope model
# global_rand_slope_model = lme(int.uniq.norm ~ 1, random = ~1 + dist.cent.global.norm|net.id, data = bird_data_final)
# summary(global_rand_slope_model)
# r.squaredGLMM(global_rand_slope_model)
# 
# 
# plot(global_rand_slope_model) # check heteroscedasticity in whole model
# plot(global_rand_slope_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
# plot(global_rand_slope_model, net.id~resid(., scaled=TRUE))
# plot(global_rand_slope_model,resid(.,scaled=TRUE)~dist.cent.global.norm|net.id, abline = 0)
# qqnorm(global_rand_slope_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
# qqnorm(global_rand_slope_model,~ranef(.)) # check normality of random effects
# plot(ranef(global_rand_slope_model)) # random effects
# 
# 
# plot(int.uniq.norm~dist.cent.global.norm, data = bird_data_final)
# global_rand_slope_model_coef = global_rand_slope_model$coefficients$random$net.id
# for (i in 1:nrow(global_rand_slope_model_coef)){
#   abline(coef = global_rand_slope_model_coef[i,], col = "blue")
# }
# abline(a = global_rand_slope_model$coefficients$fixed, b = 0, col = "red")
# 
# 
# ## fixed + random effects model===
# # random intercept model
# global_fixed_rand_int_model = lme(int.uniq.norm~dist.cent.global.norm, random = ~1|net.id, data = bird_data_final)
# summary(global_fixed_rand_int_model)
# r.squaredGLMM(global_fixed_rand_int_model)
# 
# 
# plot(global_fixed_rand_int_model) # check heteroscedasticity in whole model
# plot(global_fixed_rand_int_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
# plot(global_fixed_rand_int_model, net.id~resid(., scaled=TRUE))
# plot(global_fixed_rand_int_model,resid(.,scaled=TRUE)~dist.cent.global.norm|net.id, abline = 0)
# qqnorm(global_fixed_rand_int_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
# qqnorm(global_fixed_rand_int_model,~ranef(.)) # check normality of random effects
# plot(ranef(global_fixed_rand_int_model)) # random effects
# 
# 
# 
# plot(int.uniq.norm~dist.cent.global.norm, data = bird_data_final)
# global_fixed_rand_int_model_coef = global_fixed_rand_int_model$coefficients$random$net.id
# for (i in 1:length(global_fixed_rand_int_model_coef)){
#   abline(a = global_fixed_rand_int_model_coef[i], b = global_fixed_rand_int_model$coefficients$fixed[2], col = "blue")
# }
# abline(coef = global_fixed_rand_int_model$coefficients$fixed, col = "red")
# 
# 
# # random slope model
# global_fixed_rand_slope_model = lme(int.uniq.norm~dist.cent.global.norm, random = ~1 + dist.cent.global.norm|net.id, data = bird_data_final)
# summary(global_fixed_rand_slope_model)
# r.squaredGLMM(global_fixed_rand_slope_model)
# 
# 
# 
# plot(global_fixed_rand_slope_model) # check heteroscedasticity in whole model
# plot(global_fixed_rand_slope_model,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
# plot(global_fixed_rand_slope_model, net.id~resid(., scaled=TRUE))
# plot(global_fixed_rand_slope_model,resid(.,scaled=TRUE)~dist.cent.global.norm|net.id, abline = 0)
# qqnorm(global_fixed_rand_slope_model, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
# qqnorm(global_fixed_rand_slope_model,~ranef(.)) # check normality of random effects
# plot(ranef(global_fixed_rand_slope_model)) # random effects
# 
# 
# 
# plot(int.uniq.norm~dist.cent.global.norm, data = bird_data_final)
# global_fixed_rand_slope_model_coef = global_fixed_rand_slope_model$coefficients$random$net.id
# for (i in 1:nrow(global_fixed_rand_slope_model_coef)){
#   abline(coef = global_fixed_rand_slope_model_coef[i,], col = "blue")
# }
# abline(coef = global_fixed_rand_slope_model$coefficients$fixed, col = "red")
# 
# 
# 
# 
# # correlation between latitude and slope
# res_ex = bird_data_final[,c('net.id','latitude')]
# res_ex = res_ex[!duplicated(res_ex),]
# res_ex$intercept = as.data.frame(global_fixed_rand_slope_model_coef)[match(rownames(global_fixed_rand_slope_model_coef),res_ex$net.id),][,1]
# res_ex$slope = as.data.frame(global_fixed_rand_slope_model_coef)[match(rownames(global_fixed_rand_slope_model_coef),res_ex$net.id),][,2]
# lat_model = lm(slope~latitude, data = res_ex)
# summary(lat_model)
# 
# nestedness = lapply(bird_plant, FUN = nested, method = 'NODF2') %>%
#   do.call("rbind",.)
# 
# res_ex$nestedness = as.data.frame(nestedness)[match(rownames(nestedness),res_ex$net.id),]
# nest_model = lm(slope~nestedness, data = res_ex)
# summary(nest_model)
# 
# lat_nest_model = lm(slope~latitude*nestedness, data = res_ex)
# summary(lat_nest_model)
# 

## exploring gls================================================================
# gls_model = gls(int.uniq.norm~dist.cent.network.norm, data = bird_data_final)
# vf1 = varIdent(form = ~1|net.id)
# gls_model_weighted = gls(int.uniq.norm~dist.cent.network.norm + net.id, data = bird_data_final, weights = vf1)
# anova(gls_model,gls_model_weighted,network_fixed_rand_slope_model)
