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
