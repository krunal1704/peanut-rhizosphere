library(dada2); packageVersion("dada2")
set.seed(100)

#saving plots. width and height variables storing mm length of A4 landscape page to be used throughout
width = 297
height = 210

path <- "E:/Ph.D/krunal/metagenomic data/1. Anand/Krunal/metagenomic sequences/krunal" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])

#trimming was done to remove 2 nucleotide in R1 and 11 in R2 reads. Further, primer sequences were removed by triming 17 in R1 and 21 in R2 reads.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE, maxEE=c(3,3), truncLen = c(249,240), trimLeft = c(17,21), compress=TRUE, multithread=FALSE)

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

#learn error rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#denoise data
dadaFs <- dada(filtFs, err=errF, multithread = FALSE)
dadaRs <- dada(filtRs, err=errR, multithread = FALSE)

#merge pairs and extract sequences with merged length in range 400-430
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:430]
table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

saveRDS(seqtab2, "seqtab-run1.rds")

#remove bimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)

#track reads through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(x = track, file = "reads-stats.txt", sep = "\t", quote = FALSE)

saveRDS(seqtab.nochim, "16s.rds")
save.image("16s.RData")

# Assign taxonomy using GTDB120 database
taxGTDB <- assignTaxonomy(seqtab.nochim, tryRC = TRUE, refFasta = "GTDB_bac120_arc122_ssu_r202_fullTaxo.fa", minBoot = 80)
saveRDS(taxGTDB, "taxaGTDB.rds")

#ADD SPECIES information using same database
speciesGTDB <- addSpecies (taxGTDB, tryRC = TRUE, refFasta = "GTDB_bac120_arc122_ssu_r202_Species.fa", n = 10000)
saveRDS(speciesGTDB, "speciesGTDB.rds")

write.table(taxGTDB, file = "taxaGTDB120.tsv", sep = "\t")
write.table(speciesGTDB, file = "speciesGTDB120.tsv", sep = "\t")

#summaries assignments for each taxonomic levels and plot them
GTDB120 <- apply(speciesGTDB, 2, function(x) length(which(!is.na(x))))
read.counts <- as.data.frame(rbind(GTDB120))
read.counts$Database <- row.names(read.counts)
read.counts <- reshape2::melt(data = read.counts, id = "Database")
dbplot <- ggplot(read.counts,aes(x=Database,y=value,fill=variable)) + geom_bar(stat="identity",position = "identity", alpha = 0.8, width = 0.2) + scale_y_log10() + labs(fill = "Taxonomy assignment level", y = "Reads assigned") + annotation_logticks(sides = "lr") + theme(legend.position = "bottom") + theme_classic() + scale_fill_brewer(palette = "Dark2")
ggpubr::ggsave(filename = "database-assignment.pdf", plot = dbplot, device = "pdf", width = height, height = width, units = "mm")

###metadata comparison####
#read metadata information file containing information of groups and environmental parameters. Also compare all parameters across groups.
sample.data <- read.delim("sample-data.txt")
rownames(sample.data) <- sample.data[,1]


#load required packages
library(phyloseq)
library(microbiome)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(data.table)
library(randomcoloR)
library(tidyr)
library(scales)
library(vegan)
library(ggConvexHull)
library(grid)
library(RColorBrewer)
library(ggnewscale)
library(plyr)
library(dplyr)
library(pairwiseAdonis)
library(DESeq2)
library(ggrepel)
library(corrplot)
library(ape)
library(gridExtra)
library(stringr)
library(UpSetR)
library(viridis)
library(ComplexHeatmap)

#Find out significance of all physico-chemical parameters
all.kruskal <- ggpubr::compare_means(formula = c(pH, EC, OC, P2O5, K2O, S, Mn, Cu, Fe, Zn) ~ Farm, data = sample.data, method = "kruskal.test", p.adjust.method = "BH")
write.table(x = all.kruskal, file = "Type-comparison-metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#

#Physico-chemical parameter represent in plot with significance per variable
#create data frame
farmparameter <- sample.data[,c(1,6,8:17)] %>% group_by(Farm)
variable = rep(farmparameter$variable, each= 1) #change variable
Farm = rep(c(farmparameter$Name))
value = farmparameter$variable #change variable
data=data.frame(variable, Farm, value)
#plot bar-plot
pH <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "pH") + ylab("pH")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 2.6e-06", colour = "black", y.position = 8.5, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
EC <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "EC") + ylab(expression(paste("EC m", Omega,"/cm")))+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 1.5e-05", colour = "black", y.position = 1.2, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
OC <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "OC") + ylab("OC (%)")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 3e-05", colour = "black", y.position = 1, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
P2O5 <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(expression(paste("P"[2],O[5]))) + ylab(expression(paste("P"[2],O[5], (Kg/Ac))))+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 0.11", colour = "black", y.position = 12, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
K2O <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(expression(paste("K"[2],O))) + ylab(expression(paste("K"[2],O (Kg/Ac))))+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 2.0e-07", colour = "black", y.position = 600, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
S <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "S") + ylab("S (ppm)")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 0.75", colour = "black", y.position = 11, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
Zn <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "Zn") + ylab("Zn (ppm)")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 2.5e-06", colour = "black", y.position = 2, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
Fe <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "Fe") + ylab("Fe (ppm)")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 1.8e-07", colour = "black", y.position = 6, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
Mn <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "Mn") + ylab("Mn (ppm)")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 1.3e-06", colour = "black", y.position = 9, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")
Cu <- ggplot(data, aes(x=Farm, colour = sample.data$Farm, y=value, fill=variable)) + geom_bar(stat="identity", fill="white", position=position_dodge()) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ggtitle(label = "Cu") + ylab("Cu (ppm)")+ theme(plot.title = element_text(size = 15, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold")) + guides(fill="none") + geom_bracket(data = data, label = "Kruskal-Wallis, p - 2.6e-06", colour = "black", y.position = 3, xmin = "F03-1", xmax = "F03-1", label.size = 6) + scale_color_manual(values = color.beta) + labs(colour = "Farm")

pc <- ggarrange(pH, EC, OC, P2O5, K2O, S, Zn, Fe, Mn, Cu, nrow = 5, ncol = 2)
ggsave(filename = "Physico-chemical propeties of soil.pdf", plot = pc, device = "pdf", width = height * 3.5, height = width * 2.5, units = "mm")


#PCA plot of soil property
soil.pca <- prcomp(sample.data[, c(8:17)], center = TRUE, scale. = TRUE)
soil.pca.plot <- autoplot(soil.pca, data = sample.data, colour= 'Farm', frame = TRUE, label = TRUE, label.size = 5)
#soil.pca.plot <- autoplot(soil.pca, data = sample.data, colour= 'District', frame = TRUE, label = TRUE, label.size = 5, frame.type = 'norm')
ggsave(filename = "soil.pca.pdf", plot = soil.pca.plot, device = "pdf", width = width*1.25, height = height*1.25, units = "mm")


#create phyloseq object with metadata information and taxonomy from GTDB
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sample.data), tax_table(speciesGTDB))

#renaming sequence as ASV names to something more comfortable. "ASV" as prefix followed by numbers.
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
saveRDS(object = ps, file = "ps.RDS")

#to check distribution of reads per sample and counts per ASVs. This will be used to determine further filters.
sort(sample_sums(ps))
table(taxa_sums(ps))

#removing samples with read counts less than 10,000 (except F10-1,4 and F04-2) and more than 1,70,000. This will be the final object analysed throughout.
ps1 <- prune_samples(sample_sums(ps) > 9600, ps)
ps1
ps1.1 <- prune_samples(sample_sums(ps1) < 170000, ps1)
ps1.1
saveRDS(object = ps1.1, file = "ps1.RDS")

#removing ASVs with total count less than 21
ps2 <- prune_taxa(taxa_sums(ps1.1) > 20, ps1.1)
ps2
saveRDS(object = ps2, file = "ps2.RDS")

#plotting rarefaction plot
rarefaction.plot <- ranacapa::ggrare(ps2, step = 500, color = "Farm", label = "Name", se = FALSE )
ggsave(filename = "rarefactionplot.pdf", plot = rarefaction.plot, device = "pdf", width = width*1.25, height = height*1.25, units = "mm")

ps2.rel <- microbiome::transform(ps2, "compositional")

#########ALPHA DIVERSITY
#create a metadata data frame to add diversity values.
metadata <- data.frame(sample_data(ps2))
#calculating diversity
div <- microbiome::alpha(x = ps2, index = "all")
metadata$ShannonDiversity = div$diversity_shannon
metadata$ObservedASV = div$observed

compare_means(formula = ObservedASV ~ Farm, data = metadata, method = "wilcox.test", p.adjust.method = "BH")
compare_means(formula = ShannonDiversity ~ Farm, data = metadata, method = "wilcox.test", p.adjust.method = "BH")

color.beta <- rev(c(brewer.pal(n = 12, name = "Paired"), "gray28", "blue3", "deeppink3"))
color.beta[5] <- "yellow3"

ggplot(data = metadata, mapping = aes(x = Farm, y = ObservedASV, color = District)) + geom_boxplot() + stat_compare_means(method = "kruskal.test", label.y = 2000, label.x = 1)
ggplot(data = metadata, mapping = aes(x = Farm, y = ShannonDiversity, color = District)) + geom_boxplot() + stat_compare_means(method = "kruskal.test", label.y = 7, label.x = 1)
alphaObserved <- ggplot(data = metadata, mapping = aes(x = Farm, y = ObservedASV, color = District)) + geom_boxplot() + stat_compare_means(method = "kruskal.test", label.y = 2000, label.x = 1) + scale_color_manual(values = color.beta)
alphaShannon <- ggplot(data = metadata, mapping = aes(x = Farm, y = ShannonDiversity, color = District)) + geom_boxplot() + stat_compare_means(method = "kruskal.test", label.y = 7, label.x = 1) + scale_color_manual(values = color.beta)
palpha <- ggarrange(alphaObserved, alphaShannon, ncol = 1, align = "hv")
palpha
ggsave(filename = "alpha-plot.pdf", plot = palpha, device = "pdf", width = width*1.25, height = height*1.25, units = "mm")


########plotting taxonomy
#####PHYLUM
# define the levels to glom 
ps2.Phylum.rel <- tax_glom(physeq = ps2.rel, taxrank = "Phylum", NArm = FALSE)
#adding phylum names and changing names of NA to Unclassified taxa
taxa_names(ps2.Phylum.rel) <- tax_table(ps2.Phylum.rel)[, 2]
taxa_names(ps2.Phylum.rel)[is.na(taxa_names(ps2.Phylum.rel))] <- "Unclassified Phylum"
#create dataframe
ps2.Phylum.rel.df <- data.table(psmelt(ps2.Phylum.rel))
#group df by taxa and calculate mean rel. abundance
ps2.Phylum.rel.df[, mean := mean(Abundance, na.rm = TRUE), by = "OTU"]
#calculate the phyla-wise group significance
phylum.type <- ggpubr::compare_means(formula = Abundance ~ Farm, data = ps2.Phylum.rel.df[ps2.Phylum.rel.df$mean > 0.0001,], method = "kruskal.test", p.adjust.method = "BH", group.by = "OTU")
phylum.type.wilcox <- ggpubr::compare_means(formula = Abundance ~ Farm, data = ps2.Phylum.rel.df, method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")

write.table(x = phylum.type, file = "phylum-type.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Change name of taxa less than 1%
ps2.Phylum.rel.df[(mean <= 0.01), OTU := "Less Abundant Phyla"]
# Creating df with summarized lesser abundant taxa abundance
ps2.Phylum.rel.df <- ps2.Phylum.rel.df[, sum(Abundance), by = list(OTU,Sample,Description,Farm)]
colnames(ps2.Phylum.rel.df)[5] <- "Abundance"
#plot box plot
Rhizo.Phylum.box <- ggplot(ps2.Phylum.rel.df, aes(x=OTU, y=Abundance*100, fill = factor(Farm))) + geom_point(position=position_dodge(width=0.75), show.legend = FALSE, size=1) + geom_boxplot(alpha = 0.7) + scale_y_log10() + theme(legend.position = "left") + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Phylum") + xlab("") + theme(axis.text.x = element_text(size = 10), plot.margin=unit(c(1,0.25,0.25,0.25), "cm")) + labs(fill = "Farm")
ggsave(filename = "Rhizo-phylum.jpeg", plot = Rhizo.Phylum.box, device = "jpeg", width = width*1.25, height = height*1, units = "mm")
#get color codes, store in object and remember to not run that code again or different colors everytime
colcodes.phylum <- distinctColorPalette(length(unique(ps2.Phylum.rel.df$OTU))+13)
#plotting stacked bar chart
Phylum.bar <- ggplot(data=ps2.Phylum.rel.df, aes(x=Sample, y=Abundance*100, fill=OTU)) + geom_bar(position="stack", stat="identity", color = "black", size = 0.05, width = 1)  + xlab("Samples") + ylab("Relative abundance of Phyla")   + scale_fill_manual(values = colcodes.phylum) + labs(tag = "A", fill = "Phylum") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + coord_cartesian(expand = FALSE ) + guides(fill=guide_legend(ncol=1)) + scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))
ggsave(filename = "Phylum-bar-plot.pdf", plot = Phylum.bar, device = "pdf", width = width, height = height, units = "mm")
##

#####Genus
# # define the levels to glom 
ps2.Genus.rel <- tax_glom(physeq = ps2.rel, taxrank = "Genus", NArm = FALSE)
#melt to dataframe and convert factor to character
ps2.Genus.rel.df <- data.table(psmelt(ps2.Genus.rel))
ps2.Genus.rel.df$Kingdom <- as.character(ps2.Genus.rel.df$Kingdom)
ps2.Genus.rel.df$Phylum <- as.character(ps2.Genus.rel.df$Phylum)
ps2.Genus.rel.df$Class <- as.character(ps2.Genus.rel.df$Class)
ps2.Genus.rel.df$Order <- as.character(ps2.Genus.rel.df$Order)
ps2.Genus.rel.df$Family <- as.character(ps2.Genus.rel.df$Family)
ps2.Genus.rel.df$Genus <- as.character(ps2.Genus.rel.df$Genus)
#NA entries in Genera were renamed to "highest annotated level"_X
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Phylum),]$Phylum <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Phylum),]$Kingdom, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Class),]$Class <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Class),]$Phylum, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Order),]$Order <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Order),]$Class, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Family),]$Family <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Family),]$Order, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Genus),]$Genus <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Genus),]$Family, "_X")

genus.type <- ggpubr::compare_means(formula = Abundance ~ Farm, data = ps2.Genus.rel.df, method = "kruskal.test", p.adjust.method = "BH", group.by = "Genus")
genus.type.wilcox <- ggpubr::compare_means(formula = Abundance ~ Farm, data = ps2.Genus.rel.df, method = "wilcox.test", p.adjust.method = "BH", group.by = "Genus")

write.table(x = genus.type, file = "genus-type.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = genus.type.wilcox, file = "genus-type-wilcoxon.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


ps2Genus.df = ps2.Genus.rel.df
##group df by taxa and calculate mean rel. abundance and maximum abundance per genus
ps2.Genus.rel.df[, mean := ((mean(Abundance, na.rm = TRUE))), by = "Genus"]
#ps2.Genus.rel.df[, maximum := max(Abundance[Abundance > 0], na.rm = TRUE), by = "Genus"]

#To get a better picture, mean > 0.01(1%) was considered.
#Those not meeting these criteria were renamed as Lesser abundant genera
ps2.Genus.rel.df[(mean < 0.01), Genus := " Less Abundant Genera (n=645)"]
#Creating df with summarized lesser abundant taxa abundance
ps2.Genus.rel.df <- ps2.Genus.rel.df[, sum(Abundance), by = list(Genus,Sample,Description,Farm)]
colnames(ps2.Genus.rel.df)[5] <- "Abundance"
# #plot box plot
Rhizo.Genus.box <- ggplot(ps2.Genus.rel.df[Genus != " Less Abundant Genera (n=645)"], aes(x=Genus, y=Abundance*100, fill = factor(Farm))) + geom_point(position=position_dodge(width=0.75), show.legend = FALSE, size = 1) + geom_boxplot(alpha = 0.7) + scale_y_log10() + theme(legend.position = "left", axis.text.x = element_text(size = 7, face = "bold")) + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Genus") + xlab("") + theme(plot.margin=unit(c(1,0.25,0.25,0.5), "cm")) + labs(fill = "Farm")
ggsave(filename = "Rhizo.Genus.box.jpeg", plot = Rhizo.Genus.box, device = "jpeg", width = 325, height = 210, units = "mm")
# get color codes, store in object and remember to not run that code again or different colors everytime
colcodes.genus <- distinctColorPalette(length(unique(ps2.Genus.rel.df$Genus))+7)
#plotting stacked bar chart
Genus.bar <- ggplot(data=ps2.Genus.rel.df, aes(x=Sample, y=Abundance*100, fill=Genus)) + geom_bar(position="stack", stat="identity", color = "black", size=0.05, width = 1)  + xlab("Samples") + ylab("Relative abundance of Genera")   + scale_fill_manual(values = colcodes.genus) + labs(tag = "B", fill = "Genus") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + coord_cartesian(expand = FALSE, ylim = c(0,55) ) + guides(fill=guide_legend(ncol=1)) + scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))
Genus.bar
ggsave(filename = "Genus-bar-plot.pdf", plot = Genus.bar, device = "pdf",dpi = 300, width = width, height = height, units = "mm")

taxo <- ggarrange(Phylum.bar, Genus.bar, ncol = 1, align = "hv")
ggsave(filename = "taxo.pdf", plot = taxo, device = "pdf", width = 250*1.25, height = 210*1.25, units = "mm")


#SPECIES & ASV
#get a species level ps object and identify top 15 species to plot separately
ps2.Species.rel <- tax_glom(physeq = ps2.rel, taxrank = "Species", NArm = FALSE)
Species.name <- names(sort(taxa_sums(subset_taxa(physeq = ps2.Species.rel, !is.na(Species))), decreasing = TRUE)[1:15])
#separate top 15 species in a new ps object and remove them from original ps object. Both of this will be plotted.
ps2.Species <- subset_taxa(physeq = ps2.Species.rel, rownames(tax_table(ps2.Species.rel)) %in% Species.name)
`%notin%` <- Negate(`%in%`)
ps2.Species.rel <- subset_taxa(physeq = ps2.Species.rel, rownames(tax_table(ps2.Species.rel)) %notin% Species.name)
#rename species name to include genus name as well
ps2.Species@tax_table@.Data[,7] = paste(ps2.Species@tax_table@.Data[,6], ps2.Species@tax_table@.Data[,7], sep = " ")

#plot the heatmap
heatmap.species <- plot_heatmap(physeq = ps2.Species, sample.order = "Name", taxa.label = "Species")
heatmap.species <- heatmap.species + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10), axis.text.y = element_text(size = 7)) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), labels[i]))}) + aes(fill = Abundance * 100) + theme(plot.margin = unit(c(1,5,5,5), "pt")) + labs(fill = "Abundance") + scale_fill_gradient(limits=c(0.0005,75) , trans=log_trans(4), na.value = "black", breaks = c(0.001,0.01,0.1,1,10,50)) + theme(legend.position="none")
heatmap.ASV <- plot_heatmap(physeq = ps2.Species.rel, sample.order = "Name", max.label = 1000)
heatmap.ASV <- heatmap.ASV  + aes(fill = Abundance * 100) + labs(fill = "Abundance") + scale_fill_gradient(limits=c(0.0005,75), trans=log_trans(4), na.value = "black", breaks = c(0.001,0.01,0.1,1,10,50)) + theme(legend.position="left") + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + theme(plot.margin = unit(c(5,5,1,25), "pt"))
heatmapPlot <- ggarrange(heatmap.ASV, heatmap.species, nrow = 2, ncol = 1, heights = c(4,1), align = "hv" )
ggsave(filename = "heatmap-ASV-Species.pdf", plot = heatmapPlot, device = "pdf", width = height*2, height = width*2, units = "mm")


########BETA DIVERSITY

#calculating bray-curtis distance and NMDS ordination
ps2.bray.dist = phyloseq::distance(ps2.rel, method="bray")
ps2.bray.ordination = ordinate(ps2.rel, method="NMDS", distance=ps2.bray.dist)
#fitting all the environmental variables(physical properties and nutrient composition) with calculated NMDS and identifying significantly associated
ef.ps2 <- envfit(ps2.bray.ordination, metadata[,8:17], permu = 999)
ef.ps2.df <- as.data.frame(scores(ef.ps2, display = "vectors"))
ef.ps2.df <- cbind(ef.ps2.df, EnvironFactor = rownames(ef.ps2.df), pvalue = ef.ps2$vectors$pvals, color = "black")
ef.ps2.df$color <- as.character(ef.ps2.df$color)
ef.ps2.df[ef.ps2.df$pvalue < 0.05,]$color <- "red"

#perform permanova on calculated distance across Type of samples
ps2.adonis <- vegan::adonis2(ps2.bray.dist ~ Farm, data=metadata)
#taken from microbiomeseq package (http://www.github.com/umerijaz/microbiomeSeq)
ps2.adn_pvalue <- ps2.adonis[[1]][["Pr(>F)"]][1]
ps2.adn_rsquared <- round(ps2.adonis[[1]][["R2"]][1],3)
#use the bquote function to format adonis results to be annotated on the ordination plot.
ps2.signi_label <- paste(cut(ps2.adn_pvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ".")))
ps2.adn_res_format <- bquote(atop(atop("PERMANOVA", R^2==~.(ps2.adn_rsquared)), atop("p-value="~.(ps2.adn_pvalue)~.(ps2.signi_label), phantom())))
#extract NMDS stress for plotting
ps2.stress.label <- paste("NMDS Stress = ", round(ps2.bray.ordination$stress, 3))
#plot ordination. Shapes are numbers 0-9 for respective collections.
color.beta <- rev(c(brewer.pal(n = 12, name = "Paired"), "gray28", "blue3", "deeppink3"))
color.beta[5] <- "yellow3"
#plot ordination
betaplot.ps2 <- plot_ordination(ps2.rel, ps2.bray.ordination, color = "Farm", shape = "District") + theme(aspect.ratio = 1) + geom_point(size = 3) +  geom_convexhull(inherit.aes = TRUE, alpha = 0.1, show.legend = FALSE ) + ggtitle("NMDS plot using Bray-Curtis distance for relative abundance of all rhizosphere samples") + labs(color = "Farm") + scale_color_manual(values = color.beta)
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.ps2 <- betaplot.ps2 + annotation_custom(grob = textGrob(label = ps2.stress.label, hjust = 0, gp = gpar(fontsize = 15, fontface = "bold")), xmin = -0.5, xmax = -0.5, ymin = 0.3, ymax = 0.3)
#add PERMANOVA results on the graph....adjust the position based on the actual graph. 
betaplot.ps2 <- betaplot.ps2 + annotation_custom(grob = textGrob(label = ps2.adn_res_format, hjust = 0, gp = gpar(fontsize = 15, fontface = "bold")), xmin = 0.15, xmax = 0.15, ymin = 0.3, ymax = 0.3)
#add the results of environment fit. Plot arrows representing each variable. Red color indicates significant association
betaplot.ps2 <- betaplot.ps2 + geom_segment(data = ef.ps2.df,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),arrow = arrow(length = unit(0.25, "cm")), colour = ef.ps2.df$color, inherit.aes = FALSE) +  geom_text(data = ef.ps2.df, aes(x = NMDS1*1.05, y = NMDS2*1.05, label = EnvironFactor), size = 3, inherit.aes = FALSE)
betaplot.ps2

ggsave(filename = "beta-plot-all.pdf", plot = betaplot.ps2, device = "pdf", width = width, height = height, units = "mm")
#plotting stressplot for calculated NMDS ordinates
stressplot(ps2.bray.ordination)
#ordisurf
ps2.ordisurf.df=data.frame(x=ps2.bray.ordination$point[,1],y=ps2.bray.ordination$point[,2],Farm=metadata[,"Farm"])
ps2.ordisurf.plot <- list()
variables <- colnames(sample.data)[8:17]
for (i in 1:length(variables))
{
  var <- variables[i]
  ordi<- vegan::ordisurf(ps2.bray.ordination,metadata[,eval(var)],plot = FALSE, bs="ds")
  ordi.grid <- ordi$grid #extracts the ordisurf object
  #str(ordi.grid) #it's a list though - cannot be plotted as is
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  ordi.mite.na <- data.frame(na.omit(ordi.mite))
  ordiplot <- ggplot() + stat_contour(data = ordi.mite.na, mapping = aes(x = x, y = y, z = z, colour = ..level..),positon="identity", size = 1) + scale_colour_gradientn(colours = viridis(n = 100, option = "H", begin = 0.05, end = 0.8))  + ggplot2::labs(colour = paste(var)) + xlab("NMDS1") + ylab("NMDS2") + theme(aspect.ratio=1) + new_scale_color() + geom_point(data=ps2.ordisurf.df, aes(x,y, color = Farm),size=2, inherit.aes = FALSE) + scale_color_manual(values=color.beta, guide=FALSE) + geom_convexhull(inherit.aes = FALSE, data = ps2.ordisurf.df, mapping = aes(x = x, y = y, fill = Farm), alpha=0.25, show.legend = FALSE) + scale_fill_manual(values=color.beta, guide=FALSE) + ggtitle(label = paste("Ordisurf for envirnomental parameter ",var, sep = "")) + guides(colour = guide_legend(order = 1, nrow = 5, byrow = T))  
  ps2.ordisurf.plot[[eval(var)]] <- ordiplot
}
ps2.ordisurf.plot.all <- ggarrange(nrow = 4, ncol = 3, plotlist = ps2.ordisurf.plot)
ggsave(filename = "beta-plot-ordisurf.pdf", plot = ps2.ordisurf.plot.all, device = "pdf", width = height * 2.2, height = width * 2.5, units = "mm")


#pairwise adonis for bray-curtis distance
pairwise.adonis.bray <- pairwise.adonis(x = ps2.bray.dist, factors = metadata$Farm)
#modify dataframe to be upper matrix
pairwise.adonis.bray$pairs <- as.character(pairwise.adonis.bray$pairs)
pairwise.adonis.bray <- separate(data = pairwise.adonis.bray[,c(1,6)], col = pairs, into = c("g1", "g2"), sep = " vs ")
pairwise.adonis.bray <- spread(data = pairwise.adonis.bray, g2, p.value)
rownames(pairwise.adonis.bray) <- pairwise.adonis.bray[,1]
pairwise.adonis.bray <- pairwise.adonis.bray[,-1]
write.table(x = pairwise.adonis.bray, file = "Rhizo.pairwise.adonis.bray.txt", quote = FALSE, sep = "\t", na = "", row.names = TRUE, col.names = TRUE)
##


#####core microbiome
#core microbiome for genera
#genus level glommed phyloseq object already exists. However, need to modify the taxa names, so will make a new phyloseq object
ps2Rhizogenus <- ps2.Genus.rel
Rhizo.taxa.genus <- tax_table(ps2Rhizogenus)
Rhizo.taxa.genus <- as.data.frame(Rhizo.taxa.genus@.Data)
colname <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
Rhizo.taxa.genus[colname] <- sapply(Rhizo.taxa.genus[colname],as.character)
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Phylum),]$Phylum <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Phylum),]$Kingdom, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Class),]$Class <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Class),]$Phylum, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Order),]$Order <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Order),]$Class, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Family),]$Family <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Family),]$Order, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Genus),]$Genus <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Genus),]$Family, "_X")
taxa_names(ps2Rhizogenus) <- Rhizo.taxa.genus[, 6]
Rhizo.genus.core <- plot_core(x = ps2Rhizogenus, min.prevalence = 0.4, prevalences = seq(.05, 1, .05), detections = round(10^seq(log10(0.001), log10(1), length = 50), digits = 4), plot.type = 'heatmap', colours = rev(brewer.pal(5, "Spectral"))) +  labs(x = "Detection Threshold", y = "Genus")  + ggtitle("Genera distribution at minimum prevalence 0.4 across rhizosphere samples") + theme(plot.title = element_text (size = 20), axis.text.x = element_text(size = 10))
ggsave(filename = "Core-microbiome-rhizosphere-genus.pdf", plot = Rhizo.genus.core, device = "pdf", width = height *2, height = width*2, units = "mm")
#Extract the core microbiome genera for correlations
rhizo.core <- core(ps2Rhizogenus, detection = 0.001, prevalence = .4)

#correlation among all core genera
Rhizo.ps2Genus.df <- psmelt(rhizo.core)
Rhizo.genus.corr <- Rhizo.ps2Genus.df[,c("Sample", "OTU", "Abundance")]
Rhizo.genus.corr <- spread(data = Rhizo.genus.corr, key = "OTU", value = "Abundance")
rownames(Rhizo.genus.corr) <- Rhizo.genus.corr$Sample
Rhizo.genus.corr <- Rhizo.genus.corr[,-1]
res <- cor(Rhizo.genus.corr)
res2<- Hmisc::rcorr(as.matrix(Rhizo.genus.corr))
pdf(file = "Rhizosphere-genus-0.001-correlation.pdf", width = height*25, height = width*15)
corrplot(res2$r, p.mat = res2$P, sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 0.5, order = "hclust", method = "square", addrect = 7)
dev.off()

############______############


#for cd approch
path2 <- "E:/Ph.D/krunal/metagenomic data/1. Anand/Krunal/culture depended sequences/culture all media data"
list.files(path2)

fnFs2 <- sort(list.files(path2, pattern="_R1_001.fastq", full.names = TRUE))
fnRs2 <- sort(list.files(path2, pattern="_R2_001.fastq", full.names = TRUE))
sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)

filtFs2 <- file.path(path2, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path2, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
names(filtFs2) <- sample.names2
names(filtRs2) <- sample.names2

plotQualityProfile(fnFs2[1:10])
plotQualityProfile(fnRs2[1:10])

out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, truncLen=c(249,245),trimLeft = c(17,21),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
head(out2)

plotQualityProfile(filtFs2[1:10])
plotQualityProfile(filtRs2[1:10])

errF2 <- learnErrors(filtFs2, multithread = FALSE)
errR2 <- learnErrors(filtRs2, multithread = FALSE)

plotErrors(errF2, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)

dadaFs2 <- dada(filtFs2, err=errF2, multithread = FALSE)
dadaRs2 <- dada(filtRs2, err=errR2, multithread = FALSE)

dadaFs2[[1]]

mergers2 <- mergePairs(dadaFs2, filtFs2, dadaRs2, filtRs2, verbose=TRUE)
seqtab.2 <- makeSequenceTable(mergers2)
table(nchar(getSequences(seqtab.2)))

seqtab2.2 <- seqtab.2[,nchar(colnames(seqtab.2)) %in% 400:430]
table(nchar(getSequences(seqtab2.2)))
hist(nchar(getSequences(seqtab2.2)), main = "Histogram of merged read length distribution", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

saveRDS(seqtab2.2, "seqtab-run2.rds")

# for comparitive analysis both runs, CI (seqtab-run1) and CD (seqtab-run2) were merge
seqtab2 <- readRDS("seqtab-run1.rds")
seqtab2.2 <- readRDS("seqtab-run2.rds")
st.all <- mergeSequenceTables(seqtab2, seqtab2.2)


seqtab.all <- removeBimeraDenovo(st.all, method="consensus", multithread=FALSE)
sum(seqtab.all)/(sum(st.all))
rowSums(seqtab.all)
saveRDS(seqtab.all, "seqtab.all.rds")

taxGTDB2 <- assignTaxonomy(seqtab.all, tryRC = TRUE, refFasta = "GTDB_bac120_arc122_ssu_r202_fullTaxo.fa", minBoot = 80)
saveRDS(taxGTDB2, "taxaall.rds")

speciesGTDB2 <- addSpecies (taxGTDB2, tryRC = TRUE, refFasta = "GTDB_bac120_arc122_ssu_r202_Species.fa", n = 10000)
saveRDS(speciesGTDB2, "speciesall.rds")
write.table(speciesGTDB2, file = "speciesall.tsv", sep = "\t")


metadata2 <- read.delim("metadata.tsv")
rownames(metadata2) <- metadata2[,1]


ps3 <- phyloseq(otu_table(seqtab.all, taxa_are_rows=FALSE), tax_table(speciesGTDB2), sample_data(metadata2))

dna <- Biostrings::DNAStringSet(taxa_names(ps3))
names(dna) <- taxa_names(ps3)
ps3 <- merge_phyloseq(ps3, dna)
taxa_names(ps3) <- paste0("ASV", seq(ntaxa(ps3)))
ps3
saveRDS(object = ps3, file = "ps3.RDS")

sort(sample_sums(ps3))
table(taxa_sums(ps3))

#removing samples with read counts less than 10,000 (except F10-1,4 and F04-2). This will be the final object analysed throughout.but sample F10-1 retain 
ps4 <- prune_samples(sample_sums(ps3) > 9500, ps3)
ps4
saveRDS(object = ps4, file = "ps4.RDS")

#removing ASVs with toal count less than 30
ps5 <- prune_taxa(taxa_sums(ps4) > 30, ps4)
ps5
saveRDS(object = ps5, file = "ps5.RDS")


#plotting rarefaction plot
rarefaction.plot.all <- ranacapa::ggrare(ps5, step = 500, color = "Type",label = "Name", se = FALSE )
ggsave(filename = "rarefactionplotall.pdf", plot = rarefaction.plot.all, device = "pdf", width = width*1.25, height = height*1.25, units = "mm")


#######Alpha diversity
#plot
alpha <- plot_richness(ps5, measures = c("Observed", "Shannon"), x = "Type", color = "Farm") + xlab("Approach")
#save
ggsave(filename = "alpha-plot.pdf", plot = alpha, device = "pdf", width = width, height = height, units = "mm")

#create a metadata data frame to add diversity values.
metadata2 <- data.frame(sample_data(ps5))

#extract values for statistical comparison
div <- microbiome::alpha(x = ps5, index = "all")
metadata2$ShannonDiversity = div$diversity_shannon
metadata2$ObservedASV = div$observed
##Comparisons
#comparison between type of approach
compare_means(formula = ObservedASV ~ Type, data = metadata2, method = "wilcox.test", p.adjust.method = "BH")
compare_means(formula = ShannonDiversity ~ Type, data = metadata2, method = "wilcox.test", p.adjust.method = "BH")
#comparison between collection
compare_means(formula = ObservedASV ~ Description, data = metadata2, method = "wilcox.test", p.adjust.method = "BH")
compare_means(formula = ShannonDiversity ~ Description, data = metadata2, method = "wilcox.test", p.adjust.method = "BH")
###

###Taxonomy exploration
#summarizing phyla level assignments
all <- as.data.frame(table(tax_table(ps5)[,2], useNA = "ifany"), stringsAsFactors = FALSE)
#summarizing phyla from CD samples only; include samples from CD only, remove taxa with sum(taxa) =0, make phylum level count and summarize
cd <- as.data.frame(table(tax_table(prune_taxa(taxa_sums(subset_samples(physeq = ps5, Type == "CD")) > 0, subset_samples(physeq = ps5, Type == "CD")))[,2], useNA = "ifany"), stringsAsFactors = FALSE)
#summarizing and plotting
colnames(all) <- c("Phylum", "All samples")
all$Phylum[is.na(all$Phylum)] <- "Unclassified phylum"
colnames(cd) <- c("Phylum", "CD samples")
ASV.count <- gather(data = merge(x = all, y = cd, all.x = TRUE), key = "Samples", value = "Frequency", -Phylum)
count.plot <- ggplot(data=ASV.count, aes(x = Phylum, y = Frequency*10, fill = Samples)) +   geom_bar(stat="identity", position=position_dodge()) + scale_y_log10(labels=function(x)x/10) + annotation_logticks(sides = "l") + theme(axis.text.x = element_text(angle = 45, size = 11, hjust = 0.8)) + theme(legend.position = c(0.5, 0.9), legend.background = element_rect(fill = NA), legend.direction = "horizontal") + geom_text(aes(label=Frequency), vjust=0, color="blue", position = position_dodge(0.9), size=3.8) + scale_fill_manual(values=c('#999999','#E69F00')) + labs(tag = "A", x = "Phyla", y = "Count of unique ASVs")
#count.plot <- ggplot(data=ASV.count, aes(x = Phylum, y = Frequency*10, fill = Samples)) +   geom_bar(stat="identity", position=position_dodge()) + scale_y_log10(labels=function(x)x/10) + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + theme(legend.position = c(0.5, 0.9), legend.background = element_rect(fill = NA), legend.direction = "horizontal") + geom_text(aes(label=Frequency), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5) + scale_fill_manual(values=c('#999999','#E69F00')) + labs(x = "Phyla", y = "Count of unique ASVs")
#count.plot <- ggplot(data=ASV.count, aes(x = Phylum, y = Frequency*10, fill = Samples)) +   geom_bar(stat="identity", position=position_dodge()) + scale_y_log10(labels=function(x)x/10) + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 3 == 1, '', ifelse(i %% 3 == 2,'\n','\n\n' )), labels[i]))}) + theme(legend.position = c(0.5, 0.9), legend.background = element_rect(fill = NA), legend.direction = "horizontal") + geom_text(aes(label=Frequency), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5) + scale_fill_manual(values=c('#999999','#E69F00')) + labs(x = "Phyla", y = "Count of unique ASVs")
ggsave(filename = "ASV-count.pdf", plot = count.plot, device = "pdf", width = height*1.6, height = width, units = "mm")


###### Plot phyla level distribution
#Agglomerate taxonomy to phylum level including not assigned ASVs on relative abundances of ps5
ps5.Phylum.rel <- tax_glom(physeq = microbiome::transform(ps5, "compositional"), taxrank = "Phylum", NArm = FALSE)
#adding phylum names and changing names of NA to Unclassified phylum
taxa_names(ps5.Phylum.rel) <- tax_table(ps5.Phylum.rel)[, 2]
taxa_names(ps5.Phylum.rel)[is.na(taxa_names(ps5.Phylum.rel))] <- "Unclassified Phylum"
##plot box plot for top 15 phyla
TopPhyla <- names(sort(taxa_sums(ps5.Phylum.rel), TRUE)[1:15])
psdf <- data.table(psmelt(prune_taxa(TopPhyla, ps5.Phylum.rel)))
colcodes.phylum <- distinctColorPalette(length(unique(psdf$OTU))+5) #commented to avoid selecting new colors everytime
Phylum.box <- ggplot(psdf, aes(x=OTU, y=Abundance*100, fill = Type))  + geom_boxplot(position = position_dodge(preserve = "single"), size = 0.2) + scale_fill_manual(values = colcodes.phylum) + scale_y_log10() + theme(legend.position = "right") + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Phylum") + xlab("Phyla") + theme(plot.margin=unit(c(1,0.25,0.25,0.75), "cm"), axis.text.x = element_text(size=10)) + labs(tag = "B", fill = "Media") + geom_vline(xintercept = seq(1.5, 14.5, by = 1))
Phylum.box
ggsave(filename = "Phylum-box-plot.pdf", plot = Phylum.box, device = "pdf", width = height, height = width/2, units = "mm")

count.phyla <- ggarrange(count.plot, Phylum.box, ncol = 1, align = "hv")
ggsave(filename = "count-phyla.jpeg", plot = count.phyla, device = "jpeg", width = 250*1.25, height = 210*1.25, units = "mm")

######## Converting to binary(presence-absence)
ps5.binary <- transform_sample_counts(ps5, function(x) ifelse(x>0,1,0))

#PCoA on bray-curtis distance
binary.bray <- distance(physeq = ps5.binary, method = "bray")
binary.ordination <- ordinate(ps5.binary, "PCoA", binary.bray)
ordinationplot <- phyloseq::plot_ordination(physeq = ps5.binary, ordination = binary.ordination, type = "samples", color = "Farm", shape = "Type") + geom_point(size=3) + theme(aspect.ratio=0.8)
ggsave(filename = "binary-ordination-plot.pdf", plot = ordinationplot, device = "pdf", width = height, height =width, units = "mm")
#checking PERMANOVA for significant differences
vegan::adonis2(formula = binary.bray ~ Type, data = metadata2, method = "bray")
vegan::adonis2(formula = binary.bray ~ Description, data = metadata2, method = "bray")

#prepare data frame for comparison. Add group wise values
taxa_table <- as.data.frame(tax_table(ps5.binary))
binary.table <- merge(taxa_table, (t(ps5.binary@otu_table@.Data)), all = TRUE, by = "row.names")
#adding/summarizing data group-wise
binary.table = binary.table %>% mutate(CI = select(., 'F01-1':'F15-5') %>% rowSums(na.rm = TRUE))
binary.table = binary.table %>% mutate(CD = select(., 'F01':'F15') %>% rowSums(na.rm = TRUE))

CD <- as.character(binary.table[binary.table$CD > 0, "Row.names"])
CI <- as.character(binary.table[binary.table$CI > 0, "Row.names"])

#plot upset plot
upset(fromList(list(CD = CD, CI = CI)))
#save the plot manually

#export CD and CI exclusive ASV table to plot Krona externally and compare taxonomy
binary.table.CD <- binary.table[binary.table$CI == 0,]
write.table(x = binary.table.CD, file = "binary.table.CD.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
binary.table.CI <- binary.table[binary.table$CI > 0,]
write.table(x = binary.table.CI, file = "binary.table.CI.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

#### Heatmap
#taxa ordering to keep CD exclusive ASVs on top
names.taxa <- taxa_names(ps5.binary)
names.taxa <- names.taxa[!names.taxa %in% binary.table.CD$Row.names]
taxa.name.order <- c(binary.table.CD$Row.names, names.taxa)
#plot heatmap
ASV.pa.heatmap <- plot_heatmap(ps5.binary, method = NULL, trans = NULL, sample.order = metadata$Name, taxa.order = rev(taxa.name.order), max.label = 6000, low = "white", na.value = "white", high = "#000033") + xlab("Samples") + ylab("ASVs")
ggsave(filename = "ASV-pa-heatmap.pdf", plot = ASV.pa.heatmap, device = "pdf", width = width, height = height, units = "mm")


########## Genus level
#agglomerate to genus level
ps5.Genus <- tax_glom(physeq = ps5, taxrank = "Genus", NArm = FALSE)
#extract taxonomy file
genus.taxonomy <- as.data.frame(tax_table(ps5.Genus)[,1:6], stringsAsFactors = FALSE)
#Many of the taxa are assigned to higher levels than genus. "_X" is added to the end of such taxa. Repitition of suffix means assignment at higher level. For example, if taxa is assigned till family level its genus name will have single suffix; if assigned till class level its genus name will have 3 suffixes.
genus.taxonomy[is.na(genus.taxonomy$Phylum),]$Phylum <- paste0(genus.taxonomy[is.na(genus.taxonomy$Phylum),]$Kingdom, "_X")
genus.taxonomy[is.na(genus.taxonomy$Class),]$Class <- paste0(genus.taxonomy[is.na(genus.taxonomy$Class),]$Phylum, "_X")
genus.taxonomy[is.na(genus.taxonomy$Order),]$Order <- paste0(genus.taxonomy[is.na(genus.taxonomy$Order),]$Class, "_X")
genus.taxonomy[is.na(genus.taxonomy$Family),]$Family <- paste0(genus.taxonomy[is.na(genus.taxonomy$Family),]$Order, "_X")
genus.taxonomy[is.na(genus.taxonomy$Genus),]$Genus <- paste0(genus.taxonomy[is.na(genus.taxonomy$Genus),]$Family, "_X")
#Convert phyloseq object to presence-absence
ps5.Genus.binary <- transform_sample_counts(ps5.Genus, function(abund) 1*(abund>0))
#adding modified taxonomy information to phyloseq object
taxa_names(ps5.Genus.binary) <- genus.taxonomy[,6]
tax_table(ps5.Genus.binary)[,6] <- genus.taxonomy[,6]

#extracting taxa information with presence-absence information
genus.binary.table <- merge(ps5.Genus.binary@tax_table@.Data, (t(ps5.Genus.binary@otu_table@.Data)), all = TRUE, by = "row.names")
#adding/summarizing group-wise counts
genus.binary.table = genus.binary.table %>% mutate(CI = select(., 'F01-1':'F15-5') %>% rowSums(na.rm = TRUE))
genus.binary.table = genus.binary.table %>% mutate(CD = select(., 'F01':'F15') %>% rowSums(na.rm = TRUE))

genus.CD <- as.character(genus.binary.table[genus.binary.table$CD > 0, "Row.names"])
genus.CI <- as.character(genus.binary.table[genus.binary.table$CI > 0, "Row.names"])

upset(fromList(list(genus.CD = genus.CD, genus.CI = genus.CI)))
#save plot manually


#extract CD exclusive genera for observation and checking media-wise detection
genus.binary.table.CD <- genus.binary.table[genus.binary.table$CI == 0,]

################THE END#############
