##load in library 
library(metabaR)

#set working directory 
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024")

results <-read.csv("PARKS_CANADA_12Svert_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
results <-read.csv("PARKS_CANADA_InvertCOI_HALF_AND_C_RAW_DATA_RESULTS.csv", sep = ",")
results <-read.csv("PARKS_CANADA_VertCOI_HALF_RAW_DATA_RESULTS.csv", sep = ",")
results <-read.csv("PARKS_CANADA_MiFish_12S_RAW_HALF_AND_C_RESULTS.csv", sep = ",")
results <-read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_C_COI_INVERT_RESULTS.csv", sep = ",")
results <-read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_12Svert_RESULTS.csv", sep = ",")
results <-read.csv("PARKS_CANADA_FASTQ_ALL_REPLICATES_MiFish_RESULTS.csv", sep = ",")

results <-read.csv("RESULTS_12S_FASTQ_STOMAC_CONTENT_2025.csv", sep = ",")
results <-read.csv("RESULTS_VertCOI_FASTQ_STOMACH_CONTENT_2025.csv", sep = ",")
results <-read.csv("RESULTS_InvertCOI_STOMACH_CONTENT_2025.csv", sep = ",")

##clean sample name 
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
write.table(df, "metabar_reads_invertCOI_SC.txt", row.names=FALSE, sep = '\t', quote=FALSE)

#######################################################
#2020
#######################################################

# Load requested package for plotting
library(ggplot2)
library(reshape2)
#setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024/12Svert_half_filter")
#setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024/12Svert_whole")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_invertCOI_half.txt",
                                file_motus = "metabar_motus_invertCOI_half.txt",
                                file_pcrs = "metabar_pcrs_invertCOI_half.txt",
                                file_samples = "metabar_samples_invertCOI_half.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_vertCOI_half.txt",
                                file_motus = "metabar_motus_vertCOI_half.txt",
                                file_pcrs = "metabar_pcrs_vertCOI_half.txt",
                                file_samples = "metabar_samples_vertCOI_half.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_12sMifish_half.txt",
                                file_motus = "metabar_motus_12Smifish_half.txt",
                                file_pcrs = "metabar_pcrs_12Smifish_half.txt",
                                file_samples = "metabar_samples_12Smifish_half.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_invert_COI.txt",
                                file_motus = "metabar_motus_invertCOI.txt",
                                file_pcrs = "metabar_pcrs_invertCOI.txt",
                                file_samples = "metabar_samples_invertCOI.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_12Svert.txt",
                                file_motus = "metabar_motus_12Svert.txt",
                                file_pcrs = "metabar_pcrs_12Svert.txt",
                                file_samples = "metabar_samples_12Svert.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_12Svert_half.txt",
                                file_motus = "metabar_motus_12Svert_half.txt",
                                file_pcrs = "metabar_pcrs_12Svert_half.txt",
                                file_samples = "metabar_samples_12Svert_half.txt")

data <- tabfiles_to_metabarlist(file_reads = "metabar_reads_12mifish.txt",
                                file_motus = "metabar_motus_12Smifish.txt",
                                file_pcrs = "metabar_pcrs_12Smifish.txt",
                                file_samples = "metabar_samples_12Smifish.txt")


# Compute the number of reads per pcr
data$pcrs$nb_reads <- rowSums(data$reads)
# Compute the number of motus per pcr
data$pcrs$nb_motus <- rowSums(data$reads>0)

check1 <- melt(data$pcrs[,c("control_type", "nb_reads", "nb_motus")])
labels = c(nb_reads = "Number of Reads", nb_motus = "Number of Motus")
check1[is.na(check1)] <- "sample"

#plot initial quality 
reads_16S <- ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y", labeller = labeller(variable = labels)) + 
  theme(axis.text.x = element_text(angle=45, h=1))

check1$control_type <- gsub("positive", "field negative", check1$control_type)

ggplot(data$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")


#plot of pcr plate (and associated number of reads)
ggpcrplate(data, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")

#list of all tag/indices used in the experiment 
tag.list <- as.character(unique(data$pcrs$tag_rev))
ggpcrtag(data,
         legend_title = "# of reads per PCR",
         FUN = function(m) {rowSums(m$reads)},
         taglist = tag.list)


summary_metabarlist(data)

data.raref = hill_rarefaction(data, nboot = 20, nsteps = 10)

##output rarefaction information 
write.csv(data.raref$hill_table, "COI_rarefaction.csv")

mifish_half <- gghill_rarefaction(data.raref) 
invert_half <- gghill_rarefaction(data.raref)
vert12S_half <- gghill_rarefaction(data.raref)
vert_half <- gghill_rarefaction(data.raref)

mifish_whole <- gghill_rarefaction(data.raref)
invert_whole <- gghill_rarefaction(data.raref)
vert12S_whole <- gghill_rarefaction(data.raref)

#set up plots 
mifish_half <- mifish_half + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
invert_half <- invert_half + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
vert12S_half <- vert12S_half + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
vert_half <- vert_half + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mifish_whole <- mifish_whole + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
invert_whole <- invert_whole + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
vert12S_whole <- vert12S_whole + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

data.raref = hill_rarefaction(data, nboot = 20, nsteps = 10)

#more descriptive labels
labels = c(D0 = "Species Richness", D1 = "Shannon Index", D2 = "Simpson Index", coverage = "Coverage")

rarefaction_12shalf <- mifish_half  +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")
rarefaction_invertCOIhalf <- invert_half +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")
rarefaction_vertCOIhalf <- vert_half +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")
rarefaction_12Sverthalf <- vert12S_half +
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
rarefaction_invertCOI <- invert_whole +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")
rarefaction_12svert <- vert12S_whole +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), strip.text = element_text(size=7), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")

x <- ggarrange(rarefaction_12shalf, rarefaction_invertCOIhalf, rarefaction_12Sverthalf, rarefaction_vertCOIhalf, rarefaction_12Smifish, rarefaction_12svert, rarefaction_invertCOI, ncol=2, nrow=3, labels=c("A", "B", "C", "D", "E", "F", "G"))
ggarrange(rarefaction_12Smifish, rarefaction_invertCOI, rarefaction_12svert, ncol=3, nrow=1, labels=c("A", "B", "C"))

#ggsave(x, 
       #filename = "plots.pdf",
      # device = "pdf",
      # height = 13, width = 9.5, units = "in")

# Identifying extraction contaminants
data <- contaslayer(data, control_types = "extraction", output_col = "not_an_extraction_conta")

#######################################################
#Stomach contents
#######################################################

# Load requested package for plotting
library(ggplot2)
library(reshape2)
setwd("/Users/cassandrepyne/Documents/HANNER_LAB/PARKS_CANADA_2024/SC_invertCOI")

data <- tabfiles_to_metabarlist(file_motus = "metabar_motus_invertCOI_SC.txt",
                                file_pcrs = "metabar_pcrs_invertCOI_SC.txt",
                                file_samples = "metabar_samples_invertCOI_SC.txt",
                                file_reads = "metabar_reads_invertCOI_SC.txt")

data <- tabfiles_to_metabarlist(file_motus = "metabar_motus_12S_SC.txt",
                                file_pcrs = "metabar_pcrs_12S_SC.txt",
                                file_samples = "metabar_samples_12S_SC.txt",
                                file_reads = "metabar_reads_12S_SC.txt")

data <- tabfiles_to_metabarlist(file_motus = "metabar_motus_vertCOI_SC.txt",
                                file_pcrs = "metabar_pcrs_vertCOI_SC.txt",
                                file_samples = "metabar_samples_vertCOI_SC.txt",
                                file_reads = "metabar_reads_vertCOI_SC.txt")

# Compute the number of reads per pcr
data$pcrs$nb_reads <- rowSums(data$reads)
# Compute the number of motus per pcr
data$pcrs$nb_motus <- rowSums(data$reads>0)

check1 <- melt(data$pcrs[,c("control_type", "nb_reads", "nb_motus")])

labels = c(nb_reads = "Number of Reads", nb_motus = "Number of Motus")
check1[is.na(check1)] <- "sample"

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y", labeller = labeller(variable = labels)) + 
  theme(axis.text.x = element_text(angle=45, h=1))


#plot of pcr plate (and associated number of reads)
ggpcrplate(data, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")

# Here the list of all tag/indices used in the experiment 
tag.list <- as.character(unique(data$pcrs$tag_rev))
ggpcrtag(data,
         legend_title = "# of reads per PCR",
         FUN = function(m) {rowSums(m$reads)},
         taglist = tag.list)


summary_metabarlist(data)

data.raref = hill_rarefaction(data, nboot = 20, nsteps = 10)
head(data.raref$hill_table)

#gghill_rarefaction(data.raref) 

labels = c(D0 = "Species Richness", D1 = "Shannon Index", D2 = "Simpson Index", coverage = "Coverage")

rarefaction_12S <- gghill_rarefaction(data.raref) +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") + theme(axis.text=element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")

# Plot the unweighted distribution of MOTUs similarity scores 

ggplot(data$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 1e3, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all MOTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
data$pcrs$seqdepth_ok <- ifelse(data$pcrs$nb_reads < 1e3, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(data$pcrs$seqdepth_ok[data$pcrs$type=="sample"]) 
nrow(data$pcrs[data$pcrs$type=="sample",])



