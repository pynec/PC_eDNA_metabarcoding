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
library(cowplot)
library(VennDiagram)
library(ggVennDiagram)

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

##############################
#BARPLOTS
##############################
##load in data and clean sample names 
#data <- read.csv("PARKS_CANADA_12Svert_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
#data <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_12Svert_RESULTS.csv", sep = ",")
#data <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_C_COI_INVERT_RESULTS.csv", sep = ",")
#data <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_MiFish_RESULTS.csv", sep = ",")
#data <- read.csv("PARKS_CANADA_InvertCOI_HALF_AND_C_RAW_DATA_RESULTS.csv", sep = ",")
#data <- read.csv("PARKS_CANADA_MiFish_12S_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
#data <- read.csv("PARKS_CANADA_VertCOI_HALF_RAW_DATA_RESULTS.csv", sep = ",")

data$SampleName <- gsub("_InvertCOI.*", "", data$SampleName)
data$SampleName <- gsub("_VertCOI", "", data$SampleName)

##############################
#12S vert half 
##############################
data <- read.csv("PARKS_CANADA_12Svert_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
data$SampleName <- gsub("_S5.*", "", data$SampleName)
data$SampleName <- gsub("_S3.*", "", data$SampleName)
data$SampleName <- gsub("_S2.*", "", data$SampleName)
data$SampleName <- gsub("_S8.*", "", data$SampleName)
data$SampleName <- gsub("_S4.*", "", data$SampleName)
data$SampleName <- gsub("_S6.*", "", data$SampleName)
data$SampleName <- gsub("_S7.*", "", data$SampleName)
data$SampleName <- gsub("_S9.*", "", data$SampleName)
data$SampleName <- gsub("_S1.*", "", data$SampleName)
data$SampleName <- gsub("12S_", "", data$SampleName)
data$SampleName <- gsub("_12S", "", data$SampleName)
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]

#filter by species in geographic area 
common <- Filt(hbs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)

common_bs <- common_bs[!grepl("NC", common_bs$SampleName),]
common_bs <- common_bs[!grepl("TW", common_bs$SampleName),]

common_graph_filtered <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_graph_filtered$Species) 
                           
numSpecies <- length(unique(common_graph_filtered$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample 
common_filtered <- ggplot(data = common_graph_filtered) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
common_filtered + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

##############################
#12S vert all 
##############################
data <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_12Svert_RESULTS.csv", sep = ",")
data$SampleName <- gsub("_S5.*", "", data$SampleName)
data$SampleName <- gsub("_S3.*", "", data$SampleName)
data$SampleName <- gsub("_S2.*", "", data$SampleName)
data$SampleName <- gsub("_S8.*", "", data$SampleName)
data$SampleName <- gsub("_S4.*", "", data$SampleName)
data$SampleName <- gsub("_S6.*", "", data$SampleName)
data$SampleName <- gsub("_S7.*", "", data$SampleName)
data$SampleName <- gsub("_S9.*", "", data$SampleName)
data$SampleName <- gsub("_S1.*", "", data$SampleName)
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]

#filter by species in geographic area 
common <- Filt(hbs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)

common_bs <- common_bs[!grepl("NC", common_bs$SampleName),]
common_bs <- common_bs[!grepl("TW", common_bs$SampleName),]

common_graph_filtered <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered$Species) 


numSpecies <- length(unique(common_graph_filtered$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample 
common_filtered <- ggplot(data = common_graph_filtered) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
common_filtered + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

##############################
#12S mifish all 
##############################
data <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_MiFish_RESULTS.csv", sep = ",")
data$SampleName <- gsub("_S5.*", "", data$SampleName)
data$SampleName <- gsub("_S3.*", "", data$SampleName)
data$SampleName <- gsub("_S2.*", "", data$SampleName)
data$SampleName <- gsub("_S8.*", "", data$SampleName)
data$SampleName <- gsub("_S4.*", "", data$SampleName)
data$SampleName <- gsub("_S6.*", "", data$SampleName)
data$SampleName <- gsub("_S7.*", "", data$SampleName)
data$SampleName <- gsub("_S9.*", "", data$SampleName)
data$SampleName <- gsub("_S1.*", "", data$SampleName)
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]

#filter by species in geographic area 
common <- Filt(hbs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)

common_bs <- common_bs[!grepl("NC", common_bs$SampleName),]
common_bs <- common_bs[!grepl("TW", common_bs$SampleName),]

common_graph_filtered <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered$Species) 


numSpecies <- length(unique(common_graph_filtered$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample 
common_filtered <- ggplot(data = common_graph_filtered) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
common_filtered + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

##############################
#12S mifish half  
##############################
data <- read.csv("PARKS_CANADA_MiFish_12S_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
data$SampleName <- gsub("_S5.*", "", data$SampleName)
data$SampleName <- gsub("_S3.*", "", data$SampleName)
data$SampleName <- gsub("_S2.*", "", data$SampleName)
data$SampleName <- gsub("_S8.*", "", data$SampleName)
data$SampleName <- gsub("_S4.*", "", data$SampleName)
data$SampleName <- gsub("_S6.*", "", data$SampleName)
data$SampleName <- gsub("_S7.*", "", data$SampleName)
data$SampleName <- gsub("_S9.*", "", data$SampleName)
data$SampleName <- gsub("_S1.*", "", data$SampleName)
data$SampleName <- gsub("12S_", "", data$SampleName)
data$SampleName <- gsub("_12S", "", data$SampleName)
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]

#filter by species in geographic area 
common <- Filt(hbs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)

common_bs <- common_bs[!grepl("NC", common_bs$SampleName),]
common_bs <- common_bs[!grepl("TW", common_bs$SampleName),]

common_graph_filtered <- common_bs %>%
  group_by(SampleName) %>%
  count(Species) 

common_graph_filtered$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus artedi", "Coregonus spp.", common_graph_filtered$Species) 

numSpecies <- length(unique(common_graph_filtered$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

all_levels <- sort(unique(rbind(common_graph_filtered, common_graph_filtered_v)$Species))
common_graph_filtered$Species <- factor(common_graph_filtered$Species, levels = all_levels)

#set up colors for common legend in combined plot 
mycolors <- c(
  "Catostomus catostomus" = "#8DD3C7",
  "Catostomus commersonii" = "#D5EFBA",
  "Coregonus spp." = "#EDECBD",
  "Esox lucius" = "#C3C0D6",
  "Hiodon alosoides" = "#DF9AA1",
  "Lota lota" = "#E48883",
  "Notropis atherinoides" = "#96A8C1",
  "Notropis hudsonius" = "#B8B29F",
  "Oncorhynchus spp." = "#F6B762",
  "Perca flavescens" = "#C7D267",
  "Percopsis omiscomaycus" = "#CDD796",
  "Prosopium williamsoni" = "#FCCDE5",
  "Rhinichthys cataractae" = "#7570b3",
  "Sander vitreus" = "#A6A1B4"
)

common_graph_filtered$Species  <- factor(common_graph_filtered$Species,
                                         levels = names(mycolors))

common_graph_filtered_v$Species <- factor(common_graph_filtered_v$Species,
                                          levels = names(mycolors))


##plot common spp (filtered) by sample with specific color/species pairs (combined legend)
common_filtered <- ggplot(data = common_graph_filtered) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(
    values = mycolors,
    breaks = names(mycolors),
    drop = FALSE) + theme(
    legend.text = element_text(face = "italic")) + 
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()


##############################
#vertCOI
##############################
data <- read.csv("PARKS_CANADA_VertCOI_HALF_RAW_DATA_RESULTS.csv", sep = ",")
data$SampleName <- gsub("_S5.*", "", data$SampleName)
data$SampleName <- gsub("_S3.*", "", data$SampleName)
data$SampleName <- gsub("_S2.*", "", data$SampleName)
data$SampleName <- gsub("_S8.*", "", data$SampleName)
data$SampleName <- gsub("_S4.*", "", data$SampleName)
data$SampleName <- gsub("_S6.*", "", data$SampleName)
data$SampleName <- gsub("_S7.*", "", data$SampleName)
data$SampleName <- gsub("_S9.*", "", data$SampleName)
data$SampleName <- gsub("_S1.*", "", data$SampleName)
data$SampleName <- gsub("VertCOI_", "", data$SampleName)
data$SampleName <- gsub("_VertCOI", "", data$SampleName)
data$Species <- gsub("_", " ", data$Species)

##filter by bootstrap
hbs <- data[which(data$sBP >= 0.97), ,drop=FALSE]

#filter by species in geographic area 
common <- Filt(hbs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common)

common_bs <- common_bs[!grepl("NC", common_bs$SampleName),]
common_bs <- common_bs[!grepl("TW", common_bs$SampleName),]
common_bs <- common_bs[!grepl("Sebastes pinniger", common_bs$Species),]

common_graph_filtered_v <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered_v$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered_v$Species) 
common_graph_filtered_v$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_v$Species) 
common_graph_filtered_v$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered_v$Species) 
common_graph_filtered_v$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_graph_filtered_v$Species) 
common_graph_filtered_v$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_graph_filtered_v$Species) 

common_graph_filtered_v$Species <- factor(common_graph_filtered_v$Species, levels = all_levels)

numSpecies <- length(unique(common_graph_filtered_v$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample with species/color pairs (combined legend)
common_filtered_v <- ggplot(data = common_graph_filtered_v) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(
    values = mycolors,
    breaks = names(mycolors),
    drop = FALSE) + theme(legend.text = element_text(face = "italic")) + 
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

#create new dataframe with species from both 12S and COI 
combined_data <- rbind(common_graph_filtered, common_graph_filtered_v)

#empty plot with combined legend 
legend_plot <- ggplot(combined_data,
                     aes(x = SampleName, fill = Species, y = n)) + geom_bar(stat = "identity") +
  scale_fill_manual(
    values = mycolors,
    breaks = names(mycolors),
    drop = FALSE) + theme(legend.text = element_text(face = "italic"))

#extract just the combined legend
legend <- get_legend(legend_plot)

#plot combined plot (12S and COI) with combined legend 
ggarrange(
  ggarrange(common_filtered + theme(legend.position="none"), common_filtered_v + theme(legend.position="none"), ncol = 2,
    labels = c("A", "B"),
    label.x = c(0, 0),   # force both labels to left
    label.y = 1),
  legend,
  ncol = 2,
  widths = c(4, 1))

##############################
#VENN DIAGRAMS
##############################
half <- read.csv("PARKS_CANADA_MiFish_12S_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
half$SampleName <- gsub("_S1.*", "", half$SampleName)
half$SampleName <- gsub("_S5.*", "", half$SampleName)
half$SampleName <- gsub("_S3.*", "", half$SampleName)
half$SampleName <- gsub("_S2.*", "", half$SampleName)
half$SampleName <- gsub("_S8.*", "", half$SampleName)
half$SampleName <- gsub("_S4.*", "", half$SampleName)
half$SampleName <- gsub("_S6.*", "", half$SampleName)
half$SampleName <- gsub("_S7.*", "", half$SampleName)
half$SampleName <- gsub("_S9.*", "", half$SampleName)
half$SampleName <- gsub("12S_", "", half$SampleName)
half$SampleName <- gsub("_12S", "", half$SampleName)

half$Species <- gsub("_", " ", half$Species)

##filter by bootstrap
bs <- half[which(half$sBP >= 0.97), ,drop = FALSE]
common_bs_half <- Filt(bs)

#only the species confirmed with MiFish 
data_hbs <- half[which(half$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | half$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common_bs_half)

common_bs <- common_bs[!grepl("Sebastes pinniger", common_bs$Species),]

common_graph_filtered <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_graph_filtered$Species) 
common_graph_filtered$Species <- gsub("Coregonus artedi", "Coregonus spp.", common_graph_filtered$Species) 

#########################################
whole <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_MiFish_RESULTS.csv", sep = ",")
whole$SampleName <- gsub("_S1.*", "", whole$SampleName)
whole$SampleName <- gsub("_S5.*", "", whole$SampleName)
whole$SampleName <- gsub("_S3.*", "", whole$SampleName)
whole$SampleName <- gsub("_S2.*", "", whole$SampleName)
whole$SampleName <- gsub("_S8.*", "", whole$SampleName)
whole$SampleName <- gsub("_S4.*", "", whole$SampleName)
whole$SampleName <- gsub("_S6.*", "", whole$SampleName)
whole$SampleName <- gsub("_S7.*", "", whole$SampleName)
whole$SampleName <- gsub("_S9.*", "", whole$SampleName)

whole$Species <- gsub("_", " ", whole$Species)

##filter by bootstrap
bs <- whole[which(whole$sBP >= 0.97), ,drop = FALSE]
common_bs_whole <- Filt(bs)

#only the species confirmed with MiFish 
data_hbs <- whole[which(whole$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | whole$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common_bs_whole)

common_bs <- common_bs[!grepl("Sebastes pinniger", common_bs$Species),]

common_graph_filtered_whole <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered_whole$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered_whole$Species) 
common_graph_filtered_whole$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_whole$Species) 
common_graph_filtered_whole$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered_whole$Species) 
common_graph_filtered_whole$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_graph_filtered_whole$Species) 
common_graph_filtered_whole$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_graph_filtered_whole$Species) 
common_graph_filtered_whole$Species <- gsub("Coregonus artedi", "Coregonus spp.", common_graph_filtered_whole$Species) 

#########################################
COI <- read.csv("PARKS_CANADA_VertCOI_HALF_RAW_DATA_RESULTS.csv", sep = ",")
COI$SampleName <- gsub("_S1.*", "", COI$SampleName)
COI$SampleName <- gsub("_S5.*", "", COI$SampleName)
COI$SampleName <- gsub("_S3.*", "", COI$SampleName)
COI$SampleName <- gsub("_S2.*", "", COI$SampleName)
COI$SampleName <- gsub("_S8.*", "", COI$SampleName)
COI$SampleName <- gsub("_S4.*", "", COI$SampleName)
COI$SampleName <- gsub("_S6.*", "", COI$SampleName)
COI$SampleName <- gsub("_S7.*", "", COI$SampleName)
COI$SampleName <- gsub("_S9.*", "", COI$SampleName)
COI$SampleName <- gsub("_VertCOI", "", COI$SampleName)

COI$Species <- gsub("_", " ", COI$Species)

##filter by bootstrap
bs <- COI[which(COI$sBP >= 0.97), ,drop = FALSE]
common_bs_COI <- Filt(bs)

#only the species confirmed with MiFish 
data_hbs <- COI[which(COI$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | COI$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common_bs_COI)

common_bs <- common_bs[!grepl("Sebastes pinniger", common_bs$Species),]

common_graph_filtered_COI <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered_COI$Species <- gsub("Oncorhynchus tshawytscha", "Oncorhynchus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Oncorhynchus gorbuscha", "Oncorhynchus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Oncorhynchus mykiss", "Oncorhynchus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus artedi", "Coregonus spp.", common_graph_filtered_COI$Species) 

######################################################################

#unique species for 12S half vs whole vs COI
h <- unique(common_graph_filtered$Species)
w <- unique(common_graph_filtered_whole$Species)
cc <- unique(common_graph_filtered_COI$Species)
c <- cc[cc != "Sebastes pinniger"] 

#v <- venn.diagram(list(SET1 = h, SET2 = w, SET3 = c), fill = c("red", "green", "grey"), alpha = c(0.5,0.5,0.5), cat.cex = 1.5, filename=NULL, category.names = c("12S_half", "12S_whole", "COI_half"))
#v <- venn.diagram(list(SET1 = h, SET2 = w, SET3 = c), fill = c("red", "green", "grey"), alpha = c(0.5,0.5,0.5), cat.cex = 1.5, filename=NULL, category.names = c("12S_half", "12S_whole", "COI_half"))

#lists including number of species per unique species (h vs w vs cc)
one <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve")
two <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "thirteen")
three <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "fourteen", "fifteen")

#venn diagram with matching colours to word diagram 
v <- venn.diagram(list(SET1 = one, SET2 = two, SET3 = three), fill = c("#70AD47", "#FFC000", "#4472C4"), alpha = c(0.5,0.5,0.5), cat.cex = 1.5, filename=NULL, category.names = c("12S half", "12S whole", "COI half"))

#plot 
grid.newpage()
grid.draw(v)

##############################
#VENN DIAGRAMS: invert COI
##############################

whole <- read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_C_COI_INVERT_RESULTS.csv", sep = ",")
half <- read.csv("PARKS_CANADA_InvertCOI_HALF_AND_C_RAW_DATA_RESULTS.csv", sep = ",")

whole$SampleName <- gsub("_S1.*", "", whole$SampleName)
whole$SampleName <- gsub("_S5.*", "", whole$SampleName)
whole$SampleName <- gsub("_S3.*", "", whole$SampleName)
whole$SampleName <- gsub("_S2.*", "", whole$SampleName)
whole$SampleName <- gsub("_S8.*", "", whole$SampleName)
whole$SampleName <- gsub("_S4.*", "", whole$SampleName)
whole$SampleName <- gsub("_S6.*", "", whole$SampleName)
whole$SampleName <- gsub("_S7.*", "", whole$SampleName)
whole$SampleName <- gsub("_S9.*", "", whole$SampleName)
whole$SampleName <- gsub("_InvertCOI.*", "", whole$SampleName)

half$SampleName <- gsub("_S1.*", "", half$SampleName)
half$SampleName <- gsub("_S5.*", "", half$SampleName)
half$SampleName <- gsub("_S3.*", "", half$SampleName)
half$SampleName <- gsub("_S2.*", "", half$SampleName)
half$SampleName <- gsub("_S8.*", "", half$SampleName)
half$SampleName <- gsub("_S4.*", "", half$SampleName)
half$SampleName <- gsub("_S6.*", "", half$SampleName)
half$SampleName <- gsub("_S7.*", "", half$SampleName)
half$SampleName <- gsub("_S9.*", "", half$SampleName)
half$SampleName <- gsub("_InvertCOI.*", "", half$SampleName)

half$Species <- gsub("_", " ", half$Species)
whole$Species <- gsub("_", " ", whole$Species)

##filter by bootstrap
bs <- whole[which(whole$sBP >= 0.97), ,drop = FALSE]
common_bs_whole <- FiltArth(bs)
#common_bs_whole <- Filt(bs)

bs <- half[which(half$sBP >= 0.97), ,drop = FALSE]
common_bs_half <- FiltArth(bs)
#common_bs_half <- Filt(bs)

#unique half vs whole species 
h <- unique(common_bs_half$Species)
w <- unique(common_bs_whole$Species)

v <- venn.diagram(list(SET1 = h, SET2 = w), fill = c("red", "green"), cat.dist = c(0.055, 0.055), cat.pos = c(-27, 27), alpha = c(0.5,0.5), cat.cex = 1.5, filename=NULL, category.names = c("invertCOI_half", "invertCOI_whole"))

#plot 
grid.newpage()
grid.draw(v)

##############################################################################
#STOMACH CONTENTS
##############################################################################
data <- read.csv("RESULTS_12S_FASTQ_STOMAC_CONTENT_2025.csv", sep = ",")
#data <- read.csv("RESULTS_VertCOI_FASTQ_STOMACH_CONTENT_2025.csv", sep = ",")
#data <- read.csv("RESULTS_InvertCOI_STOMACH_CONTENT_2025.csv", sep = ",")

data$SampleName <- gsub("_S1.*", "", data$SampleName)
data$SampleName <- gsub("_S5.*", "", data$SampleName)
data$SampleName <- gsub("_S3.*", "", data$SampleName)
data$SampleName <- gsub("_S2.*", "", data$SampleName)
data$SampleName <- gsub("_S8.*", "", data$SampleName)
data$SampleName <- gsub("_S4.*", "", data$SampleName)
data$SampleName <- gsub("_S6.*", "", data$SampleName)
data$SampleName <- gsub("_S7.*", "", data$SampleName)
data$SampleName <- gsub("_S9.*", "", data$SampleName)
data$SampleName <- gsub("_InvertCOI.*", "", data$SampleName)
data$SampleName <- gsub("_VertCOI.*", "", data$SampleName)
data$SampleName <- gsub("12S_", "", data$SampleName)
data$SampleName <- gsub("_12S.*", "", data$SampleName)
data$Species <- gsub("_", " ", data$Species)

################################################
#12S
################################################
bs <- data[which(data$sBP >= 0.97), ,drop = FALSE]
common_bs_12S <- Filt(bs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
common_bs <- rbind(data_hbs, common_bs_12S)

common_graph_filtered_12S <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered_12S <- common_graph_filtered_12S[!grepl("NC", common_graph_filtered_12S$SampleName),]
common_graph_filtered_12S$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_12S$Species) 

numSpecies <- length(unique(common_graph_filtered_12S$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample 
common_filtered <- ggplot(data = common_graph_filtered_12S) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
c <- common_filtered + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

#sample per location 
qs <- c("48_SC_I", "42_SC_I", "44_SC_I", "37_SC_I", "40_SC_I", "50_SC_I", "43_SC_I", "1_SC_I", "34_SC_I", "45_SC_I", "47_SC_I", "46_SC_I", "38_SC_I", "41_SC_I", "49_SC_I", "2_SC_I", "35_SC_I", "U_SC_I")
pr <- c("107_SC_I", "105_SC_I", "104_SC_I", "39_SC_I", "114_SC_I")
ja <- c("181_SC_I", "174_SC_I", "180_SC_I", "176_SC_I", "178_SC_I", "166_SC_I", "171_SC_I", "177_SC_I", "179_SC_I", "167_SC_I")

common_bs$location <- with(common_bs, ifelse(SampleName %in% qs, 'qs',
                              ifelse(SampleName %in% pr, 'pr', 'ja')))

common_graph_filtered_12S_location <- common_bs %>%
  group_by(location) %>%
  count(Species)

common_graph_filtered_12S_location$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_12S_location$Species) 
common_graph_filtered_12S_location$location <- gsub("ja", "Jackfish Area", common_graph_filtered_12S_location$location)
common_graph_filtered_12S_location$location <- gsub("qs", "Quatre Fourches", common_graph_filtered_12S_location$location)
common_graph_filtered_12S_location$location <- gsub("pr", "Peace River", common_graph_filtered_12S_location$location)

##plot common spp (filtered) by sample location  
common_filtered_location <- ggplot(data = common_graph_filtered_12S_location) +
  geom_bar(
    mapping = aes(x = location, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Location", y = "Proportion")+
  scale_fill_manual(values = mycolors) 
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) 
f <- common_filtered_location + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip() 

ggarrange(c, f,
          ncol = 2, nrow = 1, # Arrange in 2 columns, 1 row
          common.legend = TRUE, # Use a single common legend
          legend = "right", labels=c("A", "B")     # Position the legend at the bottom
)

################################################
#vertCOI
################################################ 
bs <- data[which(data$sBP >= 0.97), ,drop = FALSE]
#common_bs_COI <- Filt(bs)

#only the species confirmed with MiFish 
data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
#common_bs <- rbind(data_hbs, common_bs_COI)
common_bs <- rbind(data_hbs, bs)

common_graph_filtered_COI <- common_bs %>%
  group_by(SampleName) %>%
  count(Species)

common_graph_filtered_COI <- common_graph_filtered_COI[!grepl("NC", common_graph_filtered_COI$SampleName),]
common_graph_filtered_COI$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus clupeaformis", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus maraena", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus brienzii", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus migratorius", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_graph_filtered_COI$Species) 
common_graph_filtered_COI$Species <- gsub("Coregonus zuerichensis", "Coregonus spp.", common_graph_filtered_COI$Species) 

numSpecies <- length(unique(common_graph_filtered_COI$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample 
common_filtered <- ggplot(data = common_graph_filtered_COI) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
cf <- common_filtered + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

##plot common spp (filtered) by sample location  
qs <- c("48_SC_I", "42_SC_I", "44_SC_I", "37_SC_I", "40_SC_I", "50_SC_I", "43_SC_I", "1_SC_I", "34_SC_I", "45_SC_I", "47_SC_I", "46_SC_I", "38_SC_I", "41_SC_I", "49_SC_I", "2_SC_I", "35_SC_I", "U_SC_I")
pr <- c("107_SC_I", "105_SC_I", "104_SC_I", "39_SC_I", "114_SC_I")
ja <- c("181_SC_I", "174_SC_I", "180_SC_I", "176_SC_I", "178_SC_I", "166_SC_I", "171_SC_I", "177_SC_I", "179_SC_I", "167_SC_I")

common_bs$location <- with(common_bs, ifelse(SampleName %in% qs, 'qs',
                                             ifelse(SampleName %in% pr, 'pr', 'ja')))

common_graph_filtered_vertCOI_location <- common_bs %>%
  group_by(location) %>%
  count(Species)

common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus clupeaformis", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus maraena", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus brienzii", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus migratorius", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus pollan", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$Species <- gsub("Coregonus zuerichensis", "Coregonus spp.", common_graph_filtered_vertCOI_location$Species) 
common_graph_filtered_vertCOI_location$location <- gsub("ja", "Jackfish Area", common_graph_filtered_vertCOI_location$location)
common_graph_filtered_vertCOI_location$location <- gsub("qs", "Quatre Fourches", common_graph_filtered_vertCOI_location$location)
common_graph_filtered_vertCOI_location$location <- gsub("pr", "Peace River", common_graph_filtered_vertCOI_location$location)

numSpecies <- length(unique(common_graph_filtered_vertCOI_location$Species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

common_filtered_location <- ggplot(data = common_graph_filtered_vertCOI_location) +
  geom_bar(
    mapping = aes(x = location, fill = Species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Location", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
cfl <- common_filtered_location + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

ggarrange(cf, cfl,
          ncol = 2, nrow = 1, # Arrange in 2 columns, 1 row
          common.legend = TRUE, # Use a single common legend
          legend = "right", labels=c("A", "B")     # Position the legend at the bottom
)

################################################
#invertCOI
################################################ 
bs <- data[which(data$sBP >= 0.97), ,drop = FALSE]

common_graph_filtered_invert <- bs %>%
  group_by(SampleName) %>%
  count(Family)

numSpecies <- length(unique(common_graph_filtered_invert$Family))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

##plot common spp (filtered) by sample 
common_filtered <- ggplot(data = common_graph_filtered_invert) +
  geom_bar(
    mapping = aes(x = SampleName, fill = Family, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Family") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
common_filtered + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

##plot common spp (filtered) by sample location  
qs <- c("48_SC_I", "42_SC_I", "44_SC_I", "37_SC_I", "40_SC_I", "50_SC_I", "43_SC_I", "1_SC_I", "34_SC_I", "45_SC_I", "47_SC_I", "46_SC_I", "38_SC_I", "41_SC_I", "49_SC_I", "2_SC_I", "35_SC_I", "U_SC_I")
pr <- c("107_SC_I", "105_SC_I", "104_SC_I", "39_SC_I", "114_SC_I")
ja <- c("181_SC_I", "174_SC_I", "180_SC_I", "176_SC_I", "178_SC_I", "166_SC_I", "171_SC_I", "177_SC_I", "179_SC_I", "167_SC_I")

bs$location <- with(bs, ifelse(SampleName %in% qs, 'qs',
                                             ifelse(SampleName %in% pr, 'pr', 'ja')))

common_graph_filtered_invertCOI_location <- bs %>%
  group_by(location) %>%
  count(Family)

common_graph_filtered_invertCOI_location$location <- gsub("ja", "Jackfish Area", common_graph_filtered_invertCOI_location$location)
common_graph_filtered_invertCOI_location$location <- gsub("qs", "Quatre Fourches", common_graph_filtered_invertCOI_location$location)
common_graph_filtered_invertCOI_location$location <- gsub("pr", "Peace River", common_graph_filtered_invertCOI_location$location)

common_filtered_location <- ggplot(data = common_graph_filtered_invertCOI_location) +
  geom_bar(
    mapping = aes(x = location, fill = Family, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Location", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Family") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
common_filtered_location + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

##############################################################################
#STOMACH CONTENTS mBRAVE
##############################################################################
df<-readr::read_tsv("STOMACH_CONTENT_mBRAVE_30-05-2025.tsv")
  
vertCOI <- df[grepl("VertCOI", df$sampleId),]
vertCOI$sampleId <- gsub("-VertCOI.*", "", vertCOI$sampleId)

vertCOI <- vertCOI %>%
  filter(species != "")

#unique(vertCOI$species)
#common_bs_vert <- RegFilt_mb(vertCOI)

##plot common spp (filtered) by sample location  
qs <- c("48-SC-I", "42-SC-I", "44-SC-I", "37-SC-I", "40-SC-I", "50-SC-I", "43-SC-I", "1-SC-I", "34-SC-I", 
        "45-SC-I", "47-SC-I", "46-SC-I", "38-SC-I", "41-SC-I", "49-SC-I", "2-SC-I", "35-SC-I", "U-SC-I")
pr <- c("107-SC-I", "105-SC-I", "104-SC-I", "39-SC-I", "114-SC-I")
ja <- c("181-SC-I", "174-SC-I", "180-SC-I", "176-SC-I", "178-SC-I", "166-SC-I", "171-SC-I", "177-SC-I", "179-SC-I", "167-SC-I")

vertCOI$location <- with(vertCOI, ifelse(sampleId %in% qs, 'qs',
                               ifelse(sampleId %in% pr, 'pr', 'ja')))

common_graph_filtered_vert <- vertCOI %>%
  group_by(location) %>%
  count(species)

#####################################################

common_graph_filtered_vert <- common_bs_vert %>%
  group_by(sampleId) %>%
  count(species)

common_graph_filtered_vert$species <- gsub("Coregonus laurettae, Coregonus sp. 177_ryapushka, Coregonus ladogae, Coregonus autumnalis, Coregonus nasus, Coregonus lucinensis, Coregonus peled, Coregonus pollan, Coregonus sp. 193_peled, Coregonus maraena, Coregonus sp. 197_peled, Coregonus sp. 179_ryapushka, Coregonus albula, Coregonus sp., Coregonus sp. 196_peled, Coregonus sp. 173_ryapushka, Coregonus sp. 198_peled, Coregonus hoyi, Coregonus sp. 187_ryapushka, Coregonus sardinella, Coregonus zenithicus, Coregonus nigripinnis, Coregonus artedi, Coregonus sp. 184_ryapushka, Coregonus fontanae, Coregonus kiyi", "Coregonus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Coregonus sp. 170_sig, Coregonus cf. suidteri, Coregonus muksun, Coregonus migratorius, Coregonus clupeaformis, Coregonus pravdinellus, Coregonus autumnalis, Coregonus nasus, Coregonus alpinus, Coregonus nobilis, Coregonus palaea, Coregonus sp. 186_muksun, Coregonus lavaretus, Coregonus peled, Coregonus bavaricus, Coregonus arenicolus, Coregonus pidschian, Coregonus sp. 208_muksun, Coregonus sp. 206_muksun, Coregonus sp. 190_sig, Coregonus sp. 180_sig, Coregonus zugensis, Coregonus candidus, Coregonus maraena, Coregonus holsatus, Coregonus oxyrinchus, Coregonus heglingus, Coregonus atterensis, Coregonus maraenoides, Coregonus sp., Coregonus zuerichensis, Coregonus sp. TAXChiem, Coregonus macrophthalmus, Coregonus sp. 189_peled, Coregonus albula, Coregonus sp. 195_muksun, Coregonus sp. 181_sig, Coregonus duplex, Coregonus sp. 188_sig, Coregonus fatioi, Coregonus albellus, Coregonus widegreni, Coregonus sp. 199_muksun, Coregonus sp. cluncaformis, Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Coregonus ussuriensis", "Coregonus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Coregonus tugun", "Coregonus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Hiodon alosoides T2537", "Hiodon alosoides", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Esox lucius, Esox cisalpinus", "Esox spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Esox aquitanicus", "Esox spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Lynceus sp. 1 NA", "Lynceus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Ichthyborus sp. NM-2010", "Ichthyborus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Mus musculus, Mus musculus C57BL/6J, Mus domesticus, Mus macedonicus", "Mus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Mus musculus", "Mus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Mus spretus", "Mus spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Canis anthus, Canis lycaon, Canis familiaris, Canis lupus, Canis aureus", "Canis spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Homo sapiens, Homo neanderthalensis, Homo heidelbergensis, Homo denisova", "Homo spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Homo sapiens", "Homo spp.", common_graph_filtered_vert$species) 
common_graph_filtered_vert$species <- gsub("Catostomus catostomus, Catostomus commersonii", "Catostomus spp.", common_graph_filtered_vert$species) 

common_graph_filtered_vert$location <- gsub("ja", "Jackfish Area", common_graph_filtered_vert$location)
common_graph_filtered_vert$location <- gsub("qs", "Quatre Fourches", common_graph_filtered_vert$location)
common_graph_filtered_vert$location <- gsub("pr", "Peace River", common_graph_filtered_vert$location)


numSpecies <- length(unique(common_graph_filtered_vert$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

ggplot(data = common_graph_filtered_vert) +
  geom_bar(
    mapping = aes(x = sampleId, fill = species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Species") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

common_filtered_location <- ggplot(data = common_graph_filtered_vert) +
  geom_bar(
    mapping = aes(x = location, fill = species, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "location", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "species") + scale_x_discrete(guide = guide_axis(n.dodge = 1))
common_filtered_location + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()


##############################################################################

invertCOI <- df[grepl("InvertCOI", df$sampleId),]
invertCOI$sampleId <- gsub("-InvertCOI.*", "", invertCOI$sampleId)

invertCOI <- invertCOI %>%
  filter(species != "")

qs <- c("48-SC-I", "42-SC-I", "44-SC-I", "37-SC-I", "40-SC-I", "50-SC-I", "43-SC-I", "1-SC-I", "34-SC-I", 
        "45-SC-I", "47-SC-I", "46-SC-I", "38-SC-I", "41-SC-I", "49-SC-I", "2-SC-I", "35-SC-I", "U-SC-I")
pr <- c("107-SC-I", "105-SC-I", "104-SC-I", "39-SC-I", "114-SC-I")
ja <- c("181-SC-I", "174-SC-I", "180-SC-I", "176-SC-I", "178-SC-I", "166-SC-I", "171-SC-I", "177-SC-I", "179-SC-I", "167-SC-I")

invertCOI$location <- with(invertCOI, ifelse(sampleId %in% qs, 'qs',
                                         ifelse(sampleId %in% pr, 'pr', 'ja')))

common_graph_filtered_vert <- invertCOI %>%
  group_by(location) %>%
  count(family)

common_graph_filtered_invert <- invertCOI %>%
  group_by(sampleId) %>%
  count(family)

common_graph_filtered_vert$location <- gsub("ja", "Jackfish Area", common_graph_filtered_vert$location)
common_graph_filtered_vert$location <- gsub("qs", "Quatre Fourches", common_graph_filtered_vert$location)
common_graph_filtered_vert$location <- gsub("pr", "Peace River", common_graph_filtered_vert$location)

numSpecies <- length(unique(common_graph_filtered_invert$family))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(numSpecies)

sample <- ggplot(data = common_graph_filtered_invert) +
  geom_bar(
    mapping = aes(x = sampleId, fill = family, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Sample", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Family") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

location <- ggplot(data = common_graph_filtered_vert) +
  geom_bar(
    mapping = aes(x = location, fill = family, y = n),
    position = "fill", stat='identity')+ 
  labs(x = "Location", y = "Proportion")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "Family") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) + theme(axis.text = element_text(size = 8)) + theme(legend.text = element_text(face="italic")) + coord_flip()

ggarrange(sample, location,
          ncol = 2, nrow = 1, # Arrange in 2 columns, 1 row
          common.legend = TRUE, # Use a single common legend
          legend = "right", labels=c("A", "B")     # Position the legend at the bottom
)

####################################
#to get info for table 4
####################################

#data_hbs <- data[which(data$Species %in% c("Rhinichthys cataractae", "Couesius plumbeus") | data$Genus %in% c("Oncorhynchus", "Coregonus")),]
#data_hbs <- data[which(data$Species %in% c("Catostomus commersonii", "Percopsis omiscomaycus", "Notropis hudsonius", "Catostomus catostomus")),]
#e <- common[common$SampleName %in% c("BEL_1", "BEL_2", "BEL_3"),]
#e <- data_hbs[data_hbs$SampleName %in% c("BEL_1", "BEL_2", "BEL_3"),]
#e <- e[c("SampleName", "Species", "ESVsize", "sBP")]
#t <- e %>%
#  group_by(Species)
#aggregate(ESVsize~Species, e, FUN=mean)
