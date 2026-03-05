#set working directory 
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024")

#install package to connect to BOLD API
#remotes::install_github("ropensci/bold")

#load bold package  
library(bold)

#mine BOLD for fish in Alberta and Columbia 
fish <- bold_specimens(taxon='Actinopterygii', geo=c('Alberta', 'British Columbia'), format='tsv')

#save all mined data to csv file 
write.csv(fish, "BOLD_Actinopterygii.csv", row.names=FALSE, quote=FALSE)

#list of unique fish species 
fish_spp <- unique(fish$species_name)

#save just list of unique species to csv file 
write.csv(fish_spp, "BOLD_spp_Dec132024.csv", row.names=FALSE, sep = '\t', quote=FALSE)

#mine BOLD for Arthropods in Alberta and Columbia 
insects <- bold_specimens(taxon='Arthropoda', geo=c('Alberta'), format='tsv')
insects_bc <- bold_specimens(taxon='Arthropoda', geo=c('British Columbia'), format='tsv')

#list of unique fish species 
insect_ab_uniq <- unique(insects$species_name)
insect_bc_uniq <- unique(insects_bc$species_name)
#save just list of unique species to csv file 
write.csv(insect_ab_uniq, "BOLD_spp_arth_AB_Dec142024.csv", row.names=FALSE, quote=FALSE)
write.csv(insect_bc_uniq, "BOLD_spp_arth_BC_Dec142024.csv", row.names=FALSE, quote=FALSE)


all_spp_arth <- c(insect_ab_uniq, insect_bc_uniq)
all_spp_uniq_arth <- unique(all_spp_arth)
write.csv(all_spp_uniq_arth, "BOLD_spp_arthALL_Dec132024.csv", row.names=FALSE, quote=FALSE)

####################################################################

#install package 
remotes::install_github("ropensci/rfishbase")

#install packages 
library("rfishbase")
library("dplyr")

###if command is not working, run this 
options("duckdbfs_use_nightly" = FALSE)

s <- fb_tbl("countrysub") %>% filter(C_Code %in% "124") %>% filter(CSub_Code %in% c("CA-AB", "CA-BC")) %>% select(SpecCode, CSub_Code, CurrentPresence)
spp <- left_join(fb_tbl("species"),  
                 s,  
                 relationship = "many-to-many")  
species <- spp %>% filter(CurrentPresence %in% c("present", "possible"))%>% filter(Fresh %in% c("1")) %>% mutate(sci_name = paste(Genus, Species))


#list of unique fish species 
fishbase_species <- unique(species$sci_name)

#save just lisit of unique species to csv file 
write.csv(unique(species$sci_name), "Fishbase_spp_freshwater_Jan272025.csv", row.names=FALSE, quote=FALSE)

####################################################################

all_spp <- c(fish_spp, fishbase_species)

all_spp_uniq <- unique(all_spp)

write.csv(all_spp_uniq, "BOLD_Fishbase_spp_freshwater_Jan272025.csv", row.names=FALSE, quote=FALSE)





