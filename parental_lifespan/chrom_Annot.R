######################################
#########SCRIPT TO DO CHROMOMAP#######
######################################


library(RIdeogram)
library(chromoMap)

#We load both the full dataset and the classified Agonists/Antagonists

paired_dataset <- read.csv("full_dataset_paired.csv", header = FALSE, sep = "\t")
summary(paired_dataset)
View(paired_dataset)

annotation_domains <- read.csv("annotation_domains.csv", header = FALSE, sep = "\t")


#We have to reformat what we have
chromoMap("hum_chr_coord.tsv","annotation_domains.csv",
          data_based_color_map = T,
          data_type = "categorical",
          segment_annotation = TRUE,
          data_colors = list(c("red","green","yellow","blue","grey","orange","pink","purple","brown","darkgreen")),
          legend = T, top_margin = 30,
          ch_gap = 10,
          left_margin = 40,
          lg_y = 180, lg_x = 100,
          chr_color = "black")

#With RIdeogram

require(RIdeogram)



pleios_domain <- read.table("annotation_domains.csv", sep = "\t", header = F, stringsAsFactors = F)
pleios_domain <- pleios_domain[,c("V5","V5","V2","V3","V4")]
names(pleios_domain) <- c("Type", "Shape", "Chr", "Start", "End")
pleios_domain$Chr <- gsub("[a-z]+", "\\1", pleios_domain$Chr)
pleios_domain$Shape <- ifelse(pleios_domain$Type == "cardio", "square", "circle")
data(human_karyotype, package="RIdeogram")
pleios_domain$color <- ifelse(pleios_domain$Type == "cardio", "6a3d9a", "ff7f00")
ideogram(karyotype = human_karyotype, label = pleios_domain, label_type = "marker")
convertSVG("chromosome.svg", device = "png")

Cstack_info()
ulimit -s
##########################NOTGOOD

#With karyoplot

library(karyoploteR)

pleios_domain <- read.table("annotation_domains.csv", sep = "\t", header = F, stringsAsFactors = F)
pleios_domain <- pleios_domain[,c("V5","V5","V2","V3","V4")]
names(pleios_domain) <- c("Type", "Shape", "Chr", "Start", "End")
cardio_domains <- pleios_domain[which(pleios_domain$Type == "cardio"),]
skeletal_domains <- pleios_domain[which(pleios_domain$Type == "skeletal"),]
infectious_domains <- pleios_domain[which(pleios_domain$Type == "infectious"),]
immune_domains <- pleios_domain[which(pleios_domain$Type == "immune"),]
metabolic_domains <- pleios_domain[which(pleios_domain$Type == "metabolic"),]
endocrine_domains <- pleios_domain[which(pleios_domain$Type == "endocrine"),]
nervous_domains <- pleios_domain[which(pleios_domain$Type == "nervous"),]
kp <- plotKaryotype(genome="hg19")
cardio <- toGRanges(data.frame(chr=cardio_domains$Chr, start=cardio_domains$Start,
                                          end=cardio_domains$Start+1000))
kpPlotRegions(kp, cardio, col="red")
kpPlotMarkers(kp, chr=cardio_domains$Chr, x=cardio_domains$Start, labels = cardio_domains$Type, 
              label.color = "red", ignore.chromosome.ends	= TRUE, text.orientation = "horizontal")
kpPlotDensity(kp, data=cardio, col="red")


skeletal <- toGRanges(data.frame(chr=skeletal_domains$Chr, start=skeletal_domains$Start,
                               end=skeletal_domains$Start+1))
kpPlotRegions(kp, skeletal, col="blue")
infectious <- toGRanges(data.frame(chr=infectious_domains$Chr, start=infectious_domains$Start,
                               end=infectious_domains$Start+1))
kpPlotRegions(kp, infectious, col="yellow")
immune <- toGRanges(data.frame(chr=immune_domains$Chr, start=immune_domains$Start,
                               end=immune_domains$Start+1))
kpPlotRegions(kp, immune, col="green")
cardio <- toGRanges(data.frame(chr=cardio_domains$Chr, start=cardio_domains$Start,
                               end=cardio_domains$Start+1))
kpPlotRegions(kp, cardio, col="purple")
cardio <- toGRanges(data.frame(chr=cardio_domains$Chr, start=cardio_domains$Start,
                               end=cardio_domains$Start+1))
kpPlotRegions(kp, cardio, col="brown")
cardio <- toGRanges(data.frame(chr=cardio_domains$Chr, start=cardio_domains$Start,
                               end=cardio_domains$Start+1))
kpPlotRegions(kp, cardio, col="orange")
cardio <- toGRanges(data.frame(chr=cardio_domains$Chr, start=cardio_domains$Start,
                               end=cardio_domains$Start+1))
kpPlotRegions(kp, cardio, col="pink")
cardio <- toGRanges(data.frame(chr=cardio_domains$Chr, start=cardio_domains$Start,
                               end=cardio_domains$Start+1))
kpPlotRegions(kp, cardio, col="darkgreen")



###########################################################################
#PERMUTATITION STATISTICAL EXCESS TEST OF GRANGE ANALYSIS##################
###########################################################################


#INSTALLING NEEDED PACKAGES
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("regioneR")
BiocManager::install("GenomicRanges")


library(GenomicRanges)
library(regioneR)
#REFORMATTING PLEIOTROPIES DATASET FOR GENOMIC RANGES

regions <- as.data.frame(cbind(as.character(pleios_domain$Chr), 
                               as.double(as.character(pleios_domain$Start)), 
                               as.double(as.character(pleios_domain$End))))
colnames(regions) <- c("chr", "start", "end")

swapped_start <- ifelse(as.numeric(as.character(regions$start)) > as.numeric(as.character(regions$end)), 
                        as.numeric(as.character(regions$end)),
                        as.numeric(as.character(regions$start)))

swapped_end <- ifelse(as.numeric(as.character(regions$start)) > as.numeric(as.character(regions$end)), 
                      as.numeric(as.character(regions$start)),
                      as.numeric(as.character(regions$end)))

new_regions <- as.data.frame(cbind(as.character(pleios_domain$Chr), 
                                   as.double(as.character(swapped_start)), 
                                   as.double(as.character(swapped_end))))
colnames(new_regions) <- c("chr", "start", "end")




#CREATION OF INTERESTING RANGES OF THE GENOME FOR LIFESPAN IN CATALOG

gr <- GRanges(
  seqnames = Rle(as.character(new_regions$chr), rep(1, length(new_regions$chr))),
  ranges = IRanges(start = as.numeric(as.character(new_regions$start)), 
                   end = as.numeric(as.character(new_regions$end))),
  strand = rep("+", length(new_regions$chr)))

new12_regions <- new_regions[which(new_regions$chr == "chr12"),]

gr_12 <- GRanges(
  seqnames = Rle(as.character(new12_regions$chr), rep(1, length(new12_regions$chr))),
  ranges = IRanges(start = as.numeric(as.character(new12_regions $start)), 
                   end = as.numeric(as.character(new12_regions $end))),
  strand = rep("+", length(new12_regions$chr)))


#USING HUMAN GENOME AS RANDOM CONTROL FOR PERMUTATION

human.genome <- getGenomeAndMask("hg19", mask=NA)$genome
random_gr <- createRandomRegions(nregions=1000000, length.mean=10000, length.sd=20000, genome=human.genome,
                         non.overlapping=FALSE)

#Now the RESAMPLING permnutation

randomizeRegions(gr, genome="hg19")

#custom evaluation

ownoverlap <- function(A, genome) {
  A <- toGRanges(A)
  chr_vector <- as.character(seqnames(A))
  #unit_out <- table(startsWith(chr_vector, "chr4"))[1]
  unit_out <- table(chr_vector == "chr22")[1]
  return(unit_out)
}


chr_ve
pt <- permTest(A=gr, randomize.function=randomizeRegions, genome=human.genome,
               evaluate.function=ownoverlap, ntimes = 1000)
pt
summary(pt)
plot(pt, plotType="Tailed")

B <- createRandomRegions(nregions=10000, length.mean=10000, length.sd=20000, genome=human.genome,
                         non.overlapping=FALSE)
A <- B[sample(20)]

+
#Permutation of subset
ov <- overlapPermTest(A=gr_12, B=gr, ntimes=1000, genome=human.genome, non.overlapping=FALSE)
ov
summary(ov)
plot(ov)




BiocManager::install("rGREAT")
=======
##########################3333
#############################
#####SNP GENOME-WIDE ENRICHMENT SOFTWARE########3
################################33

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(rGREAT)
library(VariantAnnotation)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")



GRANG_Obj <- makeGRangesFromDataFrame(new_regions)
loc <- locateVariants(GRANG_Obj , TxDb.Hsapiens.UCSC.hg19.knownGene,
                      AllVariants())
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

select(txdb, keys = loc$GENEID, columns="TXNAME", keytype="GENEID")


job = submitGreatJob(GRANG_Obj)
tb = getEnrichmentTables(job)


#After enrichment, let's gonna retrieve some plots

res = plotRegionGeneAssociationGraphs(job)
plotRegionGeneAssociationGraphs(job, type = 1)
res = plotRegionGeneAssociationGraphs(job, ontology = "GO Molecular Function")









library(gprofiler2)
dataset_gene_symbols <- gsnpense(c(as.character(paired_dataset$V1), as.character(paired_dataset$V8)), filter_na = FALSE)
hist(table(as.vector(as.character(dataset_gene_symbols$gene_names))))
sum(table(as.vector(as.character(dataset_gene_symbols$gene_names))))
=======
BiocManager::install("traseR")
