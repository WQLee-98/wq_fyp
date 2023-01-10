rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Desktop/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

library(dplyr)
library(reshape2)
library(purrr)
library(nlme)
library(lme4)
library(MuMIn)

## Data preparation ============================================================
# read Fricke's network data (without unidentified taxa)
network = read.csv(file.path(raw.dir, "Fricke/network_long.csv"), header = T)
# filter birds
network_bird = network[network$animal.class == 'Aves',]
length(unique(network_bird$animal.id)) # 1268 bird species

birds_network = unique(network$animal.id)

# read AVONET frugivores trait data
AVONET_frug = readRDS(file.path(data.dir, "AVONET_frug_traits.rds"))
birds_AVONET = unique(AVONET_frug$Species1)

sum(birds_network %in% birds_AVONET) # only 327 species out of 1268 species in network data are in the list of AVONET frugivores

# subset network data to frugivorous species that are in AVONET
network_frug = network_bird[network_bird$animal.id %in% birds_AVONET,] # 159 networks, 327 birds species, 1956 plant species

# subset to interactions only (value == 1)
network_int = network_frug[network_frug$value == 1,] # 158 networks, 5443 interactions, 313 bird species and 1324 plant species
network_no_int = network_frug[network_frug$value == 0,]
unique(network_no_int$animal.id)[which(!(unique(network_no_int$animal.id) %in% unique(network_int$animal.id)))]
# 14 frugivorous birds species recorded in original network but do not interact with any plants

# all interactions split by network into a list
network_int_list = split(network_int, f = network_int$net.id)

# check number of natives so that data is sensible

# creating list of unique species in each network
network_unique_frug_list = lapply(network_int_list, FUN = function(x){
  animal.id = unique(x$animal.id)
  net.id = rep(x[1,'net.id'],length(animal.id))
  locality.id = rep(x[1,'locality.id'],length(animal.id))
  locality = rep(x[1,'locality'],length(animal.id))
  region = rep(x[1,'region'],length(animal.id))
  return(as.data.frame(cbind(animal.id,net.id,locality.id,locality,region)))
})


## Functional Diversity ========================================================
# creating site x species matrix
site_species = acast(network_int, formula = net.id~animal.id, fill = 0, 
                     fun.aggregate = length, value.var = "animal.id")
site_species = ifelse(site_species > 0, 1, 0)

# creating species x traits matrix
rownames(AVONET_frug) = AVONET_frug$Species1
AVONET_frug$Beak.Length_Culmen = as.numeric(AVONET_frug$Beak.Length_Culmen) # length of culmen
AVONET_frug$Mass = as.numeric(AVONET_frug$Mass)
AVONET_frug$Beak.Width = as.numeric(AVONET_frug$Beak.Width)
AVONET_frug$Kipps.Distance = as.numeric(AVONET_frug$Kipps.Distance)
AVONET_frug$Wing.Length = as.numeric(AVONET_frug$Wing.Length)
AVONET_frug$Kipps.Index = AVONET_frug$Kipps.Distance / AVONET_frug$Wing.Length

# subset species to those in network_int
network_frug_traits = AVONET_frug[AVONET_frug$Species1 %in% unique(network_int$animal.id),]
network_frug_traits = network_frug_traits[,c("Mass","Beak.Length_Culmen",
                                             "Beak.Width","Kipps.Index")]
network_frug_traits = network_frug_traits[match(colnames(site_species),
                                                rownames(network_frug_traits)),]


# functional dispersion
source(file.path(main.dir,"src/fdFunc.R"))
fd_res = funcdiv(x = network_frug_traits, a = site_species, original = TRUE)

dist.cent = fd_res$dist.cent

# create new list with animal.id, net.id, and dist.cent
network_dist = list()
for (i in 1:length(network_unique_frug_list)){
  if(nrow(network_unique_frug_list[[i]])>1){
    network_dist[[i]] = network_unique_frug_list[[i]][match(names(dist.cent[[i]]), network_unique_frug_list[[i]]$animal.id),]
    network_dist[[i]]$dist.cent = dist.cent[[i]]
  }
  else{
    network_dist[[i]] = network_unique_frug_list[[i]]
    network_dist[[i]]$dist.cent = 0
  }
}

## some networks have 1 frugivorous bird only so dist to centre is 0


## Interaction uniqueness ======================================================
# create bird x plant matrix for each network
bird_plant = lapply(network_int_list, FUN = acast,formula = animal.id~plant.id, 
                    fill = 0, fun.aggregate = length, value.var = "plant.id")

# create function that calculates the interaction uniqueness of each species in each network
int_uniq = function(x){
  # Arguments
  #   x: bird x plant matrix for a network
  unique_int = c()
  for (i in 1:nrow(x)){
    # average number of dispersers per plant for the plants that species i interacts with
    sp_uniq = sum(colSums(x)*x[i,])/sum(x[i,])
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
  bird_data[[i]]$bird.uniq = bird_uniq[[i]][match(bird_data[[i]]$animal.id, names(bird_uniq[[i]]))]
}

# merge all dataframes into one
bird_data_final = reduce(bird_data, full_join)


## Building a linear mixed effects model =======================================
# making net.id a factor
bird_data_final$net.id = factor(bird_data_final$net.id)

# model building
# random intercept model
model1 = lme(bird.uniq~dist.cent, random = ~1|net.id, data = bird_data_final)
# random intercept and slope model
model2 = lme(bird.uniq~dist.cent, random = ~1 + dist.cent|net.id, data = bird_data_final)

summary(model1)
summary(model2)

plot(bird.uniq~dist.cent, data = bird_data_final)

plot(model1) # check heteroscedasticity in whole model
plot(model1,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
plot(model1, net.id~resid(., scaled=TRUE))
plot(model1,resid(.,scaled=TRUE)~dist.cent|net.id, abline = 0)
qqnorm(model1, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
qqnorm(model1,~ranef(.)) # check normality of random effects
plot(ranef(model1)) # random effects

r.squaredGLMM(model1)

# calculate range of distance to centroid of each network
# scale both y and x axis within each network

plot(model2) # check heteroscedasticity in whole model
plot(model2,resid(.,scaled=TRUE)~fitted(.)|net.id,abline=0) # check heteroscedasticity in each network
plot(model2, net.id~resid(., scaled=TRUE))
plot(model2,resid(.,scaled=TRUE)~dist.cent|net.id, abline = 0)
qqnorm(model2, ~resid(., type = "n")|net.id,abline=c(0,1)) # check normality of residuals nested in each network
qqnorm(model2,~ranef(.)) # check normality of random effects
plot(ranef(model2)) # random effects

r.squaredGLMM(model2)



# # remove data points with dist.cent == 0
# bird_data_test = subset(bird_data_final, dist.cent != 0)
# model3 = lme(bird.uniq~dist.cent, random = ~1|net.id, data = bird_data_test)
# # random intercept and slope model
# model4 = lme(bird.uniq~dist.cent, random = ~1 + dist.cent|net.id, data = bird_data_test)
# 
# summary(model3)
# summary(model4)

# plot(bird.uniq~dist.cent, data = bird_data_test)
# plot(model3)
# plot(model4)

# ANCOVA model
model5 = lm(bird.uniq~dist.cent*net.id, data = bird_data_final)
summary(model5)
AIC(model5,model2)

