### LIBRARY DIRECTORY ##
#.libPaths(c("C:\\Users\\501596\\Dropbox\\PhD Hull\\PhD docs\\Thesis\\R_stats\\rlib", 
#            .libPaths("C:\\Users\\501596\\Dropbox\\PhD Hull\\PhD docs\\Thesis\\R_stats\\rlib")))
.libPaths(c("Dropbox\\R_stats\\rlib", .libPaths("Dropbox\\R_stats\\rlib")))

## SETTING UP THE DIRECTORY ##
#setwd("C:/Users/501596/Dropbox/PhD Hull/PhD docs/Thesis/C.floridanus/") 	#directory from Uni PC
setwd("Dropbox/PhD Hull/PhD docs/Thesis/C.floridanus/") 			#directory from pc

dir()
rm(list=ls())
ls()

################################
##  CHECKING AND LOADING THE	##
##    REQUIRED PACKAGES		    ##
################################
pack.list = c("aod","bipartite","dplyr","EcoSimR","ggplot2","lattice","lme4",
              "lmtest","MASS","msa","plyr","reshape2","stringi","vegan","zoo")
new.packages = pack.list[!(pack.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pack.list, require, character.only=T)

sessionInfo()

#### CHECKING NUMBER OF SAMPLES OUT OF METABEAT ####

crang.raw = read.csv("R/source/COI_17092018_merged-only_nonchimera_cl97cov5_blast_min97_ali0.8_ref-by-taxonomy-readcounts.blast.tsv",
                     header=T, sep="\t", row.names=1)
head(crang.raw)

taxonomy.frame = crang.raw$taxomomy
crang.new=crang.raw[,!(colnames(crang.raw)%in% 'taxomomy')]

crang.new=t(crang.new)
crang.new=data.frame(crang.new)

#### FILTER THRESHOLD

rowSums(crang.new)
colSums(crang.new)

# We apply an initial filter of 10 reads 

crang.new[][crang.new[,] <= 10] <- 0

# Contamination in Positive samples.

pos = grep("POS", rownames(crang.new), ignore.case=T)
crang.new[pos,]

# We then use the contamination filter

f = crang.new[pos,"unassigned"]/rowSums(crang.new[pos,])
crang.ratio = crang.new

crang.ratio = crang.ratio[,]/rowSums(crang.ratio)
crang.ratio[is.na(crang.ratio)] <- 0

crang.ratio[,][crang.ratio[]<=f[1]] <- 0

#### Investigating the results

head(crang.ratio)

# We remove the unassigned reads (which belongs to other taxa aoutside of Crangonyx), 
# and the positive taxa (Osmia bicornis)

crang.ratio=crang.ratio[,!c(colnames(crang.ratio)%in%"Osmia_bicornis")]
crang.ratio=crang.ratio[,!c(colnames(crang.ratio)%in%"unassigned")]

head(crang.ratio)

#### Saving ratio file into csv
write.csv(crang.ratio, "R/metabarcoding_Crangonyx_ratio.csv")

# Exploring the data

rownames(crang.ratio)

rb.reads = grep("RB", rownames(crang.new))
ch.reads = grep("CH", rownames(crang.new))

rb.ratio = grep("RB", rownames(crang.ratio))
ch.ratio = grep("CH", rownames(crang.ratio))

rowSums(crang.new[ch.reads,])
rowSums(crang.new[rb.reads,])

crang.ch.reads = subset(crang.new[ch.reads,])
crang.ch.ratio = subset(crang.ratio[ch.ratio,])

crang.ch.reads = crang.ch.reads[,!c(colnames(crang.ch.reads)%in%"Osmia_bicornis")]
crang.ch.reads = crang.ch.reads[,!c(colnames(crang.ch.reads)%in%"unassigned")]

## Setting up the presence/absence data frame

crang.pres = as.data.frame(crang.ch.ratio)
crang.pres[,][crang.pres[] >= 0.9] <- 1
crang.pres[,][crang.pres[] <= 0.9] <- 0

colSums(crang.pres)

## Retaining only samples where Cp or Cf were detected
## Plus setting up melted frames for pres/abs, N reads and ratio

crang.pres=subset(crang.pres, subset = c(rowSums(crang.pres) > 0))

crang.pres$sample_ID = gsub(".nc.blast$", "", rownames(crang.pres))
pres.melt=melt(crang.pres)
head(pres.melt)

crang.ch.reads$sample_ID = gsub(".nc.blast$", "", rownames(crang.ch.reads))
melt.reads = melt(crang.ch.reads)

crang.ch.ratio$sample_ID = gsub(".nc.blast$", "", rownames(crang.ch.ratio))
ratio.melt = melt(crang.ch.ratio)

## WRITE OUTPUT ###

write.csv(crang.pres, "R/CH_Crangonyx_presence.csv")

# colorblind-friendly color palette
cbcolor=c("#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999")

ggplot(ratio.melt, aes(sample_ID, value)) +
  geom_bar(stat="identity", position="stack", aes(fill=variable)) +
  scale_fill_manual(values = cbcolor, name = "") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=0.3,vjust=0.3)) + 
  labs(x = "Sample ID", y = "DNA reads") +
  theme(legend.position = "top")
  

#### PRINTING ALIGNMENT OF SANGER SEQUENCES ####

getwd()

seq=readDNAMultipleAlignment("c.flor-c.pse_full_SANGER.fasta")

msaCheckNames(seq, ".")

print(seq, show="alignment", showNames=T)
rownames(seq)

?msaPrettyPrint

msaPrettyPrint(seq, output="tex", #shadingMode="diverse", showNames="left", 
#               showConsensus="none", consensusColors = "ColdHot",
#               consensusThreshold = 40,
#               showLogo="top", showLegend=T, logoColors="rasmol", 
               paperWidth = 15, paperHeight = 13,
               askForOverwrite=FALSE)

### ALIGNMENT OF SEQUENCES FOR THE RAXML TREE ###

dir()

## Presenting the alignment from the Sanger sequences
sanger_seq=readDNAMultipleAlignment("c.flor-c.pse_full_SANGER.fasta")

rownames(sanger_seq)

## Modifying alignment names so that they can be reordered ###
#rownames(tree_seq) = gsub("^.{4,9}[0-9]{1,2}_", "", rownames(tree_seq), ignore.case=T)

summary(sanger_seq)

print(sanger_seq, show="alignment", showNames=T, order="align")

#ordered = order(rownames(tree_seq))
#test_order = msa(tree_seq, order="aligned")
?msaPrettyPrint

msaPrettyPrint(sanger_seq, output="tex",# shadingMode = "identical", showNames="left",
               #shadingColors = "greens", 
               #showConsensus="none", consensusThreshold = 40,
               #showLogo="top", showLegend=T, logoColors="rasmol", 
               paperWidth = 17, paperHeight = 12,
               askForOverwrite=FALSE)

## Presenting the alignment from the tree
tree_seq=readDNAMultipleAlignment("MT-CO1@ALL@gt90_trimmed_aln_forR.fasta")

rownames(tree_seq)

## Modifying alignment names so that they can be reordered ###
rownames(tree_seq) = gsub("^.{6,9}[0-9]{1,2}_", "", rownames(tree_seq), ignore.case=T)
rownames(tree_seq)

summary(tree_seq)

print(tree_seq, show="alignment", showNames=T, order="align")

#ordered = order(rownames(tree_seq))
#test_order = msa(tree_seq, order="aligned")

msaPrettyPrint(tree_seq, output="pdf", shadingMode="identical", showNames="left", 
               showConsensus="none", consensusColors = "ColdHot",
               #consensusThreshold = 40,
               showLogo="top", showLegend=T, logoColors="rasmol", 
               paperWidth = 17, paperHeight = 12,
               askForOverwrite=FALSE)
