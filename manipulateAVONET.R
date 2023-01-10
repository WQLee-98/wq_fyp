rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Desktop/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

# read AVONET database
AVONET = read.csv(file.path(raw.dir,"AVONET1_Birdlife.csv"),header = T) 

# filter out frugivores and clean
all_frug_trait = AVONET[AVONET$Trophic.Niche == "Frugivore",]
all_frug_trait = all_frug_trait[!is.na(all_frug_trait$Sequence),] 

#save data into RDS file
saveRDS(all_frug_trait,file.path(data.dir, "AVONET_frug_traits.rds"))
