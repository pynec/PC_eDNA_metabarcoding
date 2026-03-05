##############################################################################
#Script to analyze MetaWorks output for metabarcoding data 
#Written by Cassandre Pyne and Yoamel Milian-Garcia, 2024
##############################################################################

#load packages 
library(tidyverse)
library(data.table)
library(ggplot2)
library(foreach)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(openxlsx)

#set working directory 
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024")

alberta_fishes <- read.csv("ALBERTA FISHES_SPECIES LIST.csv")
colnames(alberta_fishes)[2] <- "Species"
alberta_fishes$Species <-  gsub("\xca", "_", alberta_fishes$Species)
alberta_fishes$Species <- gsub("_", " ", alberta_fishes$Species)

fishes <- read.csv("BOLD_Fishbase_spp_freshwater_Jan272025.csv")
colnames(fishes)[1] <- "Species"

arth <- read.csv("BOLD_spp_arthALL_Dec132024.csv")
colnames(arth)[1] <- "Species"

AlbFilt <- function(df){
  x <- df
  x$Species <- lapply(x$Species, str_extract, pattern = "(\\w+\\s+\\w+)")
  x$Species <- as.character(x$Species)
  y <- x[x$Species %in% alberta_fishes$Species, ]
  return(y)
}

Filt <- function(df){
  x <- df
  x$Species <- lapply(x$Species, str_extract, pattern = "(\\w+\\s+\\w+)")
  x$Species <- as.character(x$Species)
  y <- x[x$Species %in% fishes$Species, ]
  return(y)
}

RegFilt_mb <- function(df){
  x <- df
  x$species <- lapply(x$species, str_extract, pattern = "(\\w+\\s+\\w+)")
  x$species <- as.character(x$species)
  
  y <- x[x$species %in% fishes$Species, ]
  return(y)
}

FiltArth <- function(df){
  x <- df
  x$Species <- lapply(x$Species, str_extract, pattern = "(\\w+\\s+\\w+)")
  x$Species <- as.character(x$Species)
  y <- x[x$Species %in% arth$Species, ]
  return(y)
}

#create new reads file for metabar
resultstoreads <- function(results){
  results$SampleName <- gsub("_S1.*", "", results$SampleName)
  results$SampleName <- gsub("_S5.*", "", results$SampleName)
  results$SampleName <- gsub("_S3.*", "", results$SampleName)
  results$SampleName <- gsub("_S2.*", "", results$SampleName)
  results$SampleName <- gsub("_S8.*", "", results$SampleName)
  results$SampleName <- gsub("_S4.*", "", results$SampleName)
  results$SampleName <- gsub("_S6.*", "", results$SampleName)
  results$SampleName <- gsub("_S7.*", "", results$SampleName)
  results$SampleName <- gsub("_S9.*", "", results$SampleName)
  
  id <- unique(results$SampleName)
  esv <- results$GlobalESV
  df <- data.frame(id=id)
  
  ##automate the creation of the reads file for metabaR 
  ##takes the output of MetaWorks and extracts esv information and the associated taxonomy 
  for(i in 1:length(unique(esv))){
    vec = rep(0, length(id))
    esv <- unique(esv)
    esv_of_interest <- grepl(esv[i], results$GlobalESV)
    index <- which(esv_of_interest)
    for(j in 1:length(index)){
      id_val <- results$SampleName[index[j]]
      id_index <- which(id == id_val)
      esv_size <- results$ESVsize[index[j]]
      vec[id_index] <- esv_size
    }
    df <- cbind(df, vec)
    names(df)[names(df) == "vec"] <- esv[i]
    
  }
  return(df)
}

##############################################################################
half12s <- read.csv("PARKS_CANADA_MiFish_12S_RAW_HALF_AND_C_RESULTS.csv", sep = ",")

half12s$Species <- gsub("_", " ", half12s$Species)

##filter by bootstrap
bs <- half12s[which(half12s$sBP >= 0.97), ,drop = FALSE]
common_bs_half12s <- Filt(bs)

#only the species confirmed with MiFish 
data_hbs <- half12s[which(half12s$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | half12s$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common_bs_half12s)
common_bs <- common_bs[!grepl("Sebastes pinniger", common_bs$Species),]

common_bs$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_bs$Species) 
common_bs$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_bs$Species) 
common_bs$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_bs$Species) 
common_bs$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_bs$Species) 
common_bs$Species <- gsub("Coregonus artedi", "Coregonus spp.", common_bs$Species) 

df <- resultstoreads(common_bs)
write.csv(df, "metabar_reads_12Smifish_half_filter_fish.csv")

##############################################################################
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024")

data <- read.csv("PARKS_CANADA_VertCOI_HALF_RAW_DATA_RESULTS.csv", sep = ",")
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]
common <- Filt(hbs)

data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)
common_bs <- common_bs[!grepl("Sebastes pinniger", common_bs$Species),]

common_bs$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_bs$Species) 
common_bs$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_bs$Species) 
common_bs$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_bs$Species) 
common_bs$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_bs$Species) 
common_bs$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_bs$Species) 

df <- resultstoreads(common_bs)
write.csv(df, "metabar_reads_vertCOI_half_filter_fish.csv")

##############################################################################
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024")

data <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_MiFish_RESULTS.csv", sep = ",")
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]
common <- Filt(hbs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)

common_graph_filtered$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered$Species) 

df <- resultstoreads(common_bs)
write.csv(df, "metabar_reads_12Smifish_fish.csv")

##############################################################################
##############################################################################
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024/12Smifish_whole_fish")

#load in metabar files 
data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_12Smifish_half_filter_fish.txt",
                                file_motus = "metabar_motus_12Smifish_half_fish.txt",
                                file_pcrs = "metabar_pcrs_12Smifish_half.txt",
                                file_samples = "metabar_samples_12Smifish_half.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_vertCOI_half_filter_fish.txt",
                                file_motus = "metabar_motus_vertCOI_half_fish.txt",
                                file_pcrs = "metabar_pcrs_vertCOI_half_fish.txt",
                                file_samples = "metabar_samples_vertCOI_half.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_12Smifish_fish.txt",
                                file_motus = "metabar_motus_12Smifish_fish.txt",
                                file_pcrs = "metabar_pcrs_12Smifish_fish.txt",
                                file_samples = "metabar_samples_12Smifish.txt")

#rarefaction file 
data.raref = hill_rarefaction(data, nboot = 20, nsteps = 10)

mifish_half_fish <- gghill_rarefaction(data.raref) 
vert_half <- gghill_rarefaction(data.raref)
mifish_whole <- gghill_rarefaction(data.raref)

mifish_half_fish <- mifish_half_fish + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
vert_half <- vert_half + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
mifish_whole <- mifish_whole + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#include more descriptive labels 
labels = c(D0 = "Species Richness", D1 = "Shannon Index", D2 = "Simpson Index", coverage = "Coverage")

rarefaction_12shalf <- mifish_half_fish  +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")
rarefaction_vertCOIhalf <- vert_half +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")
rarefaction_12Smifish <- mifish_whole +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")


x <- ggarrange(rarefaction_12shalf, rarefaction_vertCOIhalf, rarefaction_12Smifish, ncol=3, nrow=1, labels=c("A", "B", "C"))

#ggsave(x,
#      filename = "plots.pdf",
#       device = "pdf",
#       height = 7, width = 15, units = "in")

