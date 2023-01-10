rm(list = ls())
gc()

main.dir = "C:/Users/leewq/Desktop/FYP Working Folder"
raw.dir = file.path(main.dir,"raw")
data.dir = file.path(main.dir,"data")
results.dir = file.path(main.dir,"results")
figures.dir = file.path(main.dir,"figures")

library(tidyverse)


## Retrieving threats from IUCN database =======================================
# reading AVONET database of frugivores only
frug_traits = readRDS(file.path(data.dir, "AVONET_frug_traits.rds"))
frug_names = frug_traits$Species1

# create function that retrieves the list of threats of each frugivorous 
# species and compiles them into a list
retrieve = function(names){
  # Arguments
  #   names: vector of frugivorous species names
  # Returns
  #   List of list of threats for each frugivore
  
  # create empty list for storing queries
  threats = list()
  
  # query every species in 'names' in Red List API for their threats and 
  # append them to 'threats'
  for(species in names){
    species_threat = rredlist::rl_threats(name = species, key = "d04f4de4b444fcff51b39fd77f7c8e88a09fa94f3cb123cdc8504a8d732bec50")
    threats = append(threats, list(species_threat))
  }
  
  return (threats)
}

#frug_threats = retrieve(frug_names) # do not uncomment
#saveRDS(frug_threats, file.path(data.dir, "frug_threats.rds")) # do not uncomment

frug_threats = readRDS(file.path(data.dir, "frug_threats.rds"))

## data exploration ============================================================
# magnitude score
mag = lapply(frug_threats, FUN = function(x){return(x$result$score)})
all_mag = c()
for(i in 1:length(mag)){
  if (!is.null(mag[[i]])){
    all_mag = c(all_mag,mag[[i]])
  }
}
uniq_mag = unique(all_mag)

# past impact. currently no threat
mag_PI = which(lapply(frug_threats, FUN = function(x){any(x$result$score == "Past Impact")})==TRUE)
# unknown score. will calculate minimum score
mag_unknown = which(lapply(frug_threats, FUN = function(x){any(x$result$score == "Unknown")})==TRUE)
# NA score. no clue as to why they are recorded as NA. will assume no threat
mag_NA = which(lapply(frug_threats, FUN = function(x){any(is.na(x$result$score))})==TRUE)


# scope score
scope = lapply(frug_threats, FUN = function(x){return(x$result$scope)})
all_scope = c()
for(i in 1:length(scope)){
  if (!is.null(scope[[i]])){
    all_scope = c(all_scope,scope[[i]])
  }
}
uniq_scope = unique(all_scope)

scope_unknown = which(lapply(frug_threats, FUN = function(x){any(x$result$scope == "Unknown")})==TRUE)
scope_NA = which(lapply(frug_threats, FUN = function(x){any(is.na(x$result$scope))})==TRUE)

# severity score
severity = lapply(frug_threats, FUN = function(x){return(x$result$severity)})
all_severity = c()
for(i in 1:length(severity)){
  if (!is.null(severity[[i]])){
    all_severity = c(all_severity,severity[[i]])
  }
}
uniq_severity = unique(all_severity)

severity_unknown = which(lapply(frug_threats, FUN = function(x){any(x$result$severity == "Unknown")})==TRUE)
severity_NA = which(lapply(frug_threats, FUN = function(x){any(is.na(x$result$severity))})==TRUE)

# timing score
timing = lapply(frug_threats, FUN = function(x){return(x$result$timing)})
all_timing = c()
for(i in 1:length(timing)){
  if (!is.null(timing[[i]])){
    all_timing = c(all_timing,timing[[i]])
  }
}
uniq_timing = unique(all_timing)

timing_unknown = which(lapply(frug_threats, FUN = function(x){any(x$result$timing == "Unknown")})==TRUE)

## magnitude extraction ========================================================
# create function that extracts the maximum threat magnitude for each of the 
# threat categories (Habitat Loss, Wildlife Trade, Biological Invasion, Pollution
# Climate Change)
extract = function(species){
  # Arguments
  #   species: list of vector of species name (1st element) and data frame of threats (2nd element)
  # Returns
  #   data frame with highest magnitudes of all threats
  
  species_threats = data.frame(species = species[[1]], HL = 0, WT = 0, BI = 0, 
                               P = 0, CC = 0)
  
  # if the species has no threats, return data frame with all threat magnitude 0
  if (is.null(nrow(species[[2]]))){
    return(species_threats)
  }else{
    # storing IUCN threat codes into key/value pairs
    keys = c('1.1'='HL','1.2'='HL','1.3'='HL','2.1.1'='HL','2.1.2'='HL','2.1.3'='HL',
             '2.1.4'='HL','2.2.1'='HL','2.2.2'='HL','2.2.3'='HL','2.3.1'='HL','2.3.2'='HL',
             '2.3.3'='HL','2.3.4'='HL','2.4.1'='HL','2.4.2'='HL','2.4.3'='HL','3.1'='HL',
             '3.2'='HL','3.3'='HL','4.1'='HL','4.2'='HL','4.3'='HL','5.3.3'='HL',
             '5.3.4'='HL','5.3.5'='HL','6.1'='HL','6.2'='HL','6.3'='HL','7.1.1'='HL',
             '7.1.2'='HL','7.1.3'='HL','7.2.1'='HL','7.2.2'='HL','7.2.3'='HL',
             '7.2.4'='HL','7.2.5'='HL','7.2.6'='HL','7.2.7'='HL','7.2.8'='HL',
             '7.2.9'='HL','7.2.10'='HL','7.2.11'='HL','7.3'='HL', 
             '5.1.1'='WT','5.1.2'='WT','5.1.3'='WT','5.1.4'='WT','5.2.1'='WT',
             '5.2.2'='WT','5.2.3'='WT','5.2.4'='WT','5.3.1'='WT','5.3.2'='WT',
             '5.4.1'='WT','5.4.2'='WT','5.4.3'='WT','5.4.4'='WT','5.4.5'='WT',
             '5.4.6'='WT',
             '8.1.1'='BI', '8.1.2'='BI','8.2.1'='BI','8.2.2'='BI','8.3'='BI',
             '8.4.1'='BI','8.4.2'='BI','8.5'='BI',
             '9.1.1'='P','9.1.2'='P','9.1.3'='P','9.2.1'='P','9.2.2'='P','9.2.3'='P',
             '9.3.1'='P','9.3.2'='P','9.3.3'='P','9.3.4'='P','9.4'='P','9.5.1'='P',
             '9.5.2'='P','9.5.3'='P','9.5.4'='P','9.6.1','9.6.2'='P','9.6.3'='P',
             '9.6.4'='P',
             '11.1'='CC','11.2'='CC','11.3'='CC','11.4'='CC','11.5'='CC')
    
    # key/value pair for magnitude,timing,scope, and severity scores
    mag_score = c("No/Negligible Impact: 2" = 2,"Low Impact: 3" = 3,
               "Low Impact: 4" = 4,"Low Impact: 5" = 5,"Medium Impact: 6" = 6,
               "Medium Impact: 7" = 7,"High Impact: 8" = 8,"High Impact: 9" = 9, 
               "Past Impact" = 0)
    timing_score = c("Future" = 1, "Ongoing" = 3, "Unknown" = 3, 
                     "Past, Unlikely to Return" = 1, "Past, Likely to Return" = 1)
    scope_score = c("Minority (<50%)" = 1, "Majority (50-90%)" = 2, 
                    "Whole (>90%)" =3, "Unknown" = 1)
    severity_score = c("No decline" = 0,"Negligible declines" = 0, 
                       "Slow, Significant Declines" = 1, "Causing/Could cause fluctuations" = 1, 
                       "Rapid Declines" = 2, "Very Rapid Declines" = 3, "Unknown" = 0)
    
    # checking if threat code is in the vectors above, and extracting highest magnitude
    for (i in 1:nrow(species[[2]])){
      if (species[[2]][i,'code'] %in% names(keys)){
        prev_mag = species_threats[1,keys[species[[2]][i,'code']]]
        if(species[[2]][i,'score'] %in% names(mag_score)){
          curr_mag = mag_score[species[[2]][i,'score']]
        }else { # for 'Unknown' or NA magnitudes
          curr_mag = timing_score[species[[2]][i,'timing']] + scope_score[species[[2]][i,'scope']] 
          + severity_score[species[[2]][i,'severity']]
        }
        species_threats[1,keys[species[[2]][i,'code']]] = max(c(prev_mag,curr_mag))
      }
    }
    return (species_threats)
  }
}

# use 'extract' and then merge all lists into one data frame
frug_mag = lapply(frug_threats, FUN = extract) %>%
  reduce(full_join)

## Adding IUCN Statuses ========================================================
all_bird = read.csv(file.path(raw.dir,"IUCN Red List/all_bird_IUCN_2022/assessments.csv"), header = T)
frug_mag_stat = data.frame(frug_mag)

for (i in 1:nrow(frug_mag)){
  if (frug_mag_stat[i,'species'] %in% all_bird$scientificName){
    frug_mag_stat[i,'status'] = all_bird[all_bird$scientificName == frug_mag_stat[i,'species'],
                                         'redlistCategory']
  }else {
    frug_mag_stat[i,'status'] = "Unknown"
  }
}

# [Manual cleaning] changing Parophasma galinieri status
frug_mag_stat[frug_mag_stat$species == "Parophasma galinieri", "status"] = "Least Concern"

frug_clean = merge(frug_traits,frug_mag_stat,by.x = "Species1", by.y = "species", all.x = T)

## Saving files ================================================================
# saveRDS(frug_clean,file.path(data.dir,"frug_masterlist.rds"))
# write.csv(frug_clean,file.path(data.dir,"frug_masterlist.csv"))
