rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Desktop/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")


## Exploring relationship between threat level of frugivores vs other 
## feeding guilds ==============================================================
# read in IUCN threat status of all birds
all_bird = read.csv(file.path(raw.dir,"IUCN Red List/all_bird_IUCN_2022/assessments.csv"), header = T)
all_bird_IUCN = all_bird[c('scientificName','redlistCategory')]

AVONET = read.csv(file.path(raw.dir,"AVONET1_Birdlife.csv"),header = T) 

AVONET_IUCN = merge(AVONET, all_bird_IUCN,by.x = "Species1", by.y = "scientificName")
# species with no matching IUCN threat categories are removed

# removing extinct and data deficient species
AVONET_IUCN_extant_noDD = subset(AVONET_IUCN, redlistCategory != 'Data Deficient' 
                                 & redlistCategory !='Extinct' & redlistCategory 
                                 != 'Extinct in the Wild')

# making redlistCategory and Trophic.Niche factors
AVONET_IUCN_extant_noDD$redlistCategory = factor(AVONET_IUCN_extant_noDD$redlistCategory)
levels(AVONET_IUCN_extant_noDD$redlistCategory) = c("Critically Endangered",'Endangered',
                                        'Vulnerable','Near Threatened','Least Concern')
AVONET_IUCN_extant_noDD$Trophic.Niche = factor(AVONET_IUCN_extant_noDD$Trophic.Niche)

# creating tables of counts and propotions of redlist category and trophic niche
tab = table(AVONET_IUCN_extant_noDD$redlistCategory,AVONET_IUCN_extant_noDD$Trophic.Niche)
proptab = prop.table(tab,2)*100
barplot(tab, xlab = "Trophic Niche", legend = rownames(tab))
barplot(proptab, xlab = "Trophic Niche", legend = rownames(proptab))
  


## Exploring relationship between magnitude of threats vs threat category 
## of frugivores ===============================================================
frug_data = readRDS(file.path(data.dir,"frug_masterlist.rds"))
head(frug_data)

frug_data_extant_noDD = subset(frug_data, status != 'Data Deficient' & status 
                               != 'Extinct in the Wild')

frug_data_extant_noDD$HL = factor(frug_data_extant_noDD$HL)
frug_data_extant_noDD$WT = factor(frug_data_extant_noDD$WT)
frug_data_extant_noDD$BI = factor(frug_data_extant_noDD$BI)
frug_data_extant_noDD$P = factor(frug_data_extant_noDD$P)
frug_data_extant_noDD$CC = factor(frug_data_extant_noDD$CC)
frug_data_extant_noDD$status = factor(frug_data_extant_noDD$status, 
                                      levels = c('Critically Endangered','Endangered',
                                                 'Vulnerable','Near Threatened','Least Concern'))

# maybe take magnitudes as numeric rather than ordinal?

# levels(frug_data_extant_noDD$HL) = c('0','1','2','3','4','5','6','7','8','9')
# levels(frug_data_extant_noDD$WT) = c('0','1','2','3','4','5','6','7','8','9')
# levels(frug_data_extant_noDD$BI) = c('0','1','2','3','4','5','6','7','8','9')
# levels(frug_data_extant_noDD$P) = c('0','1','2','3','4','5','6','7','8','9')
# levels(frug_data_extant_noDD$CC) = c('0','1','2','3','4','5','6','7','8','9')

tab_HL = table(frug_data_extant_noDD$status, frug_data_extant_noDD$HL)
proptab_HL = prop.table(tab_HL,2)*100
barplot(tab_HL, xlab = "Threat Magnitude", legend = rownames(tab_HL), main = 'Habitat Loss')
barplot(proptab_HL, xlab = "Threat Magnitude", legend = rownames(proptab_HL), main = 'Habitat Loss')

tab_WT = table(frug_data_extant_noDD$status, frug_data_extant_noDD$WT)
proptab_WT = prop.table(tab_WT,2)*100
barplot(tab_WT, xlab = "Threat Magnitude", legend = rownames(tab_WT), main = 'Wildlife Trade')
barplot(proptab_WT, xlab = "Threat Magnitude", legend = rownames(proptab_WT), main = 'Wildlife Trade')
