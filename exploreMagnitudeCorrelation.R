rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Desktop/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

library(ggplot2)
library(plyr)
library(dplyr)

frug_data = readRDS(file.path(data.dir,"frug_masterlist.rds"))
head(frug_data)
iucn_red_list = data.frame(RL_categ = c("Least Concern", "Near Threatened", "Vulnerable", 
                                         "Endangered", "Critically Endangered", 
                                         "Extinct in the Wild", "Data Deficient"),
                            ExtinctionRisk = c(0, 0, 1, 1, 1, NA, NA))

frug_data_final = frug_data %>%
  # Merge with extinction risk scores
  left_join(iucn_red_list, by= c("status" = "RL_categ")) %>%
  # Remove species that are either extinct or data deficient
  subset(!is.na(ExtinctionRisk)) %>%
  # Remove fully migratory species
  # subset(Migration %in% c(1,2)) %>%
  # Only keep the following columns (has to come AFTER the previous step for obvious reasons!)
  dplyr::select(Species1, Family1, Order1, Habitat, Primary.Lifestyle,
         Mass, Centroid.Latitude, status, ExtinctionRisk, HL, WT, BI, P, CC) %>%
  # Only include rows where the following columns are complete.
  dplyr::filter(complete.cases(.[c("ExtinctionRisk", "Mass", "Centroid.Latitude", 
                                   "Primary.Lifestyle")]))

habitat_loss = frug_data_final[frug_data_final$HL >0,]
wildlife_trade = frug_data_final[frug_data_final$WT >0,]

all_HL_mod <- glm(ExtinctionRisk ~ HL, data = frug_data_final, family = "binomial")
all_WT_mod <- glm(ExtinctionRisk ~ WT, data = frug_data_final, family = "binomial")
HL_mod <- glm(ExtinctionRisk ~ HL, data = habitat_loss, family = "binomial")
WT_mod <- glm(ExtinctionRisk ~ WT, data = wildlife_trade, family = "binomial")
summary(all_HL_mod)
summary(all_WT_mod)
summary(HL_mod)
summary(WT_mod)
