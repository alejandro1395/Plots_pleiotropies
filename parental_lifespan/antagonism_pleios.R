############################
############################
#READ THE TABLE WITH PAIRED#
###########################

library(ggplot2)
library(gridExtra)
library("plotrix")
library("RColorBrewer")
library(tidyverse)

paired_dataset <- read.csv("freqs_Agon_Antagon.tsv", header = TRUE, sep = "\t")
summary(paired_dataset)
sum(colSums(paired_dataset[,c("Agonistic","Antagonistic")]))

#main barplot
paired_dataset.sort <-paired_dataset[order(rowSums(paired_dataset[,c("Antagonistic", "Agonistic")]),decreasing=T),]
data_matrix <- t(as.matrix(paired_dataset.sort[,c("Agonistic", "Antagonistic")]))
barplot(data_matrix, names.arg = paired_dataset.sort$Domains, beside = TRUE, 
        legend.text = c("Agonistic", "Antagonistic"), col = c("lavender", "gold"))

####
#ALL DISEASES PLOTS
#######
Agonistic <- read.csv("Agonistic", header = FALSE, sep = "\t")
Antagonistic <-read.csv("Antagonistic", header = FALSE, sep = "\t")
summary(Agonistic)
summary(Antagonistic)

#read onsets
Onsets <- read.csv("onsets.txt", header = FALSE, sep = "\t")
summary(Onsets)


all_Onsets<- c(as.numeric(as.character(Onsets$V2)), as.numeric(as.character(Onsets$V4)))
all_Onsets <- all_Onsets[ !is.na(all_Onsets) ]
all_diseases <- c(as.character(Onsets$V1), as.character(Onsets$V3))
elements_2_remove <- c("")
all_diseases <- all_diseases[!(all_diseases %in% elements_2_remove)]
diseases_onsets <- as.data.frame(cbind(all_Onsets, all_diseases))




#agonistic filters

agon_onsets <- c()
agon_diseases <- c(as.character(Agonistic$V2), as.character(Agonistic$V10))
elements_2_remove <- c("")
agon_diseases <- agon_diseases[!(agon_diseases %in% elements_2_remove)]
for (i in 1:length(agon_diseases)){
subset <- diseases_onsets[which(diseases_onsets$all_diseases == agon_diseases[i]),]
agon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

agon_total <- as.data.frame(cbind(agon_onsets, agon_diseases))
table_agon <- table(agon_total$agon_diseases)

#antagonistic filters

antagon_onsets <- c()
antagon_diseases <- c(as.character(Antagonistic$V2), as.character(Antagonistic$V10))
elements_2_remove <- c("")
antagon_diseases <- antagon_diseases[!(antagon_diseases %in% elements_2_remove)]
for (i in 1:length(antagon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == antagon_diseases[i]),]
  antagon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

antagon_total <- as.data.frame(cbind(antagon_onsets, antagon_diseases))
table_antagon <- table(antagon_total$antagon_diseases)

#All together

agon_total$type <- "Agonist"
colnames(agon_total) <- c("onset", "disease", "type")
antagon_total$type <- "Antagonist"
colnames(antagon_total) <- c("onset", "disease", "type")

total <- rbind(agon_total, antagon_total)

#PYRAMID
sorted_total <- total[order(as.numeric(as.character(total$onset))),]
diseaselabels<-unique(sorted_total$disease)

agon_plot_counts <- c()
antagon_plot_counts <- c()
for (i in 1:length(diseaselabels)){
  if(diseaselabels[i] %in%  names(table_antagon)){
    antagon_plot_counts[i] <- table_antagon[names(table_antagon) == diseaselabels[i]]
  }
  else{
    antagon_plot_counts[i] <- 0
    names(antagon_plot_counts[i]) <- diseaselabels[i]
  }
  result2 <- try(table_agon[names(table_agon) == diseaselabels[i]])
  if(diseaselabels[i] %in% names(table_agon)){
    agon_plot_counts[i] <- table_agon[names(table_agon) == diseaselabels[i]]
  }
  else{
    agon_plot_counts[i] <- 0
    names(agon_plot_counts[i]) <- diseaselabels[i]
  }}

final_matrix <- as.data.frame(cbind(as.character(diseaselabels), agon_plot_counts, antagon_plot_counts))




pyramid.plot(as.numeric(as.character(final_matrix$agon_plot_counts)),
             as.numeric(as.character(final_matrix$antagon_plot_counts)),
             labels=(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)
             
pyramid.plot(rev(as.numeric(as.character(final_matrix$agon_plot_counts))),
            rev(as.numeric(as.character(final_matrix$antagon_plot_counts))),
            labels=rev(as.character(final_matrix$V1)), lxcol=c("lavender"),
                    rxcol=c("gold"),laxlab=c(0,5,10,20),
                    raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)


##########################################3
#QUEDA HACER PLOT DE ANTAGONISTA EXCESO CON EDAD (CHI-CUADRADO entre agon early agon late antagon early antagon late)

onsets_dimension_total <- rbind(antagon_total, agon_total)
pvalues <- c()

for (i in 10:60){
antagon_late_subset <- onsets_dimension_total[which(as.numeric(as.character(onsets_dimension_total$onset))>i &
                                               as.character(onsets_dimension_total$type)=="Antagonist"),]
antagon_late_count <- length(antagon_late_subset$onset)
antagon_early_subset <- onsets_dimension_total[which(as.numeric(as.character(onsets_dimension_total$onset))<i &
                                                      as.character(onsets_dimension_total$type)=="Antagonist"),]
antagon_early_count <- length(antagon_early_subset$onset)
agon_late_subset <- onsets_dimension_total[which(as.numeric(as.character(onsets_dimension_total$onset))>i &
                                                      as.character(onsets_dimension_total$type)=="Agonist"),]
agon_late_count <- length(agon_late_subset$onset)
agon_early_subset <- onsets_dimension_total[which(as.numeric(as.character(onsets_dimension_total$onset))<i &
                                                      as.character(onsets_dimension_total$type)=="Agonist"),]
agon_early_count <- length(agon_early_subset$onset)
pl <- matrix(c(agon_early_count, agon_late_count, antagon_early_count, antagon_late_count), ncol = 2)
pv <- as.double(chisq.test(pl)["p.value"])
pvalues<- append(pvalues, pv)
}

#plot it

Significance_variation <- -log10(pvalues)
length(Significance_variation)
Ages_varx <- seq(10, 60)
plot(x=Ages_varx, y=Significance_variation)
abline(h = -log10(0.05), col="red", lwd=3, lty=2)
abline(v=c(34,40), col=c("blue", "blue"), lty=c(1,1), lwd=c(3, 3))







##################################################
################################################
#################################################3
#INVERTED PLOT###################################
##########################################3

paired_dataset <- read.csv("freqs_Agon_Antagon_inverted.tsv", header = TRUE, sep = "\t")
summary(paired_dataset)
sum(colSums(paired_dataset[,c("Agonistic","Antagonistic")]))

#main barplot
paired_dataset.sort <-paired_dataset[order(rowSums(paired_dataset[,c("Antagonistic", "Agonistic")]),decreasing=T),]
data_matrix <- t(as.matrix(paired_dataset.sort[,c("Agonistic", "Antagonistic")]))
barplot(data_matrix, names.arg = paired_dataset.sort$Domains, beside = TRUE, 
        legend.text = c("Agonistic", "Antagonistic"), col = c("lavender", "gold"))






##################################################
################################################
#################################################3
#CARDIO NORMAL PYRAMID PLOT###################################
##########################################3


Agonistic <- read.csv("Agonistic_cardio", header = FALSE, sep = "\t")
Antagonistic <-read.csv("Antagonistic_cardio", header = FALSE, sep = "\t")
summary(Agonistic)
summary(Antagonistic)

#read onsets
Onsets <- read.csv("onsets.txt", header = FALSE, sep = "\t")
summary(Onsets)


all_Onsets<- c(as.numeric(as.character(Onsets$V2)), as.numeric(as.character(Onsets$V4)))
all_Onsets <- all_Onsets[ !is.na(all_Onsets) ]
all_diseases <- c(as.character(Onsets$V1), as.character(Onsets$V3))
elements_2_remove <- c("")
all_diseases <- all_diseases[!(all_diseases %in% elements_2_remove)]
diseases_onsets <- as.data.frame(cbind(all_Onsets, all_diseases))




#agonistic filters

agon_onsets <- c()
agon_diseases <- c(as.character(Agonistic$V2), as.character(Agonistic$V10))
elements_2_remove <- c("")
agon_diseases <- agon_diseases[!(agon_diseases %in% elements_2_remove)]
for (i in 1:length(agon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == agon_diseases[i]),]
  agon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

agon_total <- as.data.frame(cbind(agon_onsets, agon_diseases))
table_agon <- table(agon_total$agon_diseases)

#antagonistic filters

antagon_onsets <- c()
antagon_diseases <- c(as.character(Antagonistic$V2), as.character(Antagonistic$V10))
elements_2_remove <- c("")
antagon_diseases <- antagon_diseases[!(antagon_diseases %in% elements_2_remove)]
for (i in 1:length(antagon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == antagon_diseases[i]),]
  antagon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

antagon_total <- as.data.frame(cbind(antagon_onsets, antagon_diseases))
table_antagon <- table(antagon_total$antagon_diseases)

#All together

agon_total$type <- "Agonist"
colnames(agon_total) <- c("onset", "disease", "type")
antagon_total$type <- "Antagonist"
colnames(antagon_total) <- c("onset", "disease", "type")

total <- rbind(agon_total, antagon_total)

#PYRAMID
sorted_total <- total[order(as.numeric(as.character(total$onset))),]
diseaselabels<-unique(sorted_total$disease)

agon_plot_counts <- c()
antagon_plot_counts <- c()
for (i in 1:length(diseaselabels)){
  if(diseaselabels[i] %in%  names(table_antagon)){
    antagon_plot_counts[i] <- table_antagon[names(table_antagon) == diseaselabels[i]]
  }
  else{
    antagon_plot_counts[i] <- 0
    names(antagon_plot_counts[i]) <- diseaselabels[i]
  }
  result2 <- try(table_agon[names(table_agon) == diseaselabels[i]])
  if(diseaselabels[i] %in% names(table_agon)){
    agon_plot_counts[i] <- table_agon[names(table_agon) == diseaselabels[i]]
  }
  else{
    agon_plot_counts[i] <- 0
    names(agon_plot_counts[i]) <- diseaselabels[i]
  }}

final_matrix <- as.data.frame(cbind(as.character(diseaselabels), agon_plot_counts, antagon_plot_counts))




pyramid.plot(as.numeric(as.character(final_matrix$agon_plot_counts)),
             as.numeric(as.character(final_matrix$antagon_plot_counts)),
             labels=(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)

pyramid.plot(rev(as.numeric(as.character(final_matrix$agon_plot_counts))),
             rev(as.numeric(as.character(final_matrix$antagon_plot_counts))),
             labels=rev(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)






##################################################
################################################
#################################################3
#Cardio disease INVERTED PLOT###################################
##########################################3

Agonistic <- read.csv("Agonistic_cardio_inverted", header = FALSE, sep = "\t")
Antagonistic <-read.csv("Antagonistic_cardio_inverted", header = FALSE, sep = "\t")
summary(Agonistic)
summary(Antagonistic)

#read onsets
Onsets <- read.csv("onsets.txt", header = FALSE, sep = "\t")
summary(Onsets)


all_Onsets<- c(as.numeric(as.character(Onsets$V2)), as.numeric(as.character(Onsets$V4)))
all_Onsets <- all_Onsets[ !is.na(all_Onsets) ]
all_diseases <- c(as.character(Onsets$V1), as.character(Onsets$V3))
elements_2_remove <- c("")
all_diseases <- all_diseases[!(all_diseases %in% elements_2_remove)]
diseases_onsets <- as.data.frame(cbind(all_Onsets, all_diseases))




#agonistic filters

agon_onsets <- c()
agon_diseases <- c(as.character(Agonistic$V2), as.character(Agonistic$V10))
elements_2_remove <- c("")
agon_diseases <- agon_diseases[!(agon_diseases %in% elements_2_remove)]
for (i in 1:length(agon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == agon_diseases[i]),]
  agon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

agon_total <- as.data.frame(cbind(agon_onsets, agon_diseases))
table_agon <- table(agon_total$agon_diseases)

#antagonistic filters

antagon_onsets <- c()
antagon_diseases <- c(as.character(Antagonistic$V2), as.character(Antagonistic$V10))
elements_2_remove <- c("")
antagon_diseases <- antagon_diseases[!(antagon_diseases %in% elements_2_remove)]
for (i in 1:length(antagon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == antagon_diseases[i]),]
  antagon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

antagon_total <- as.data.frame(cbind(antagon_onsets, antagon_diseases))
table_antagon <- table(antagon_total$antagon_diseases)

#All together

agon_total$type <- "Agonist"
colnames(agon_total) <- c("onset", "disease", "type")
antagon_total$type <- "Antagonist"
colnames(antagon_total) <- c("onset", "disease", "type")

total <- rbind(agon_total, antagon_total)

#PYRAMID
sorted_total <- total[order(as.numeric(as.character(total$onset))),]
diseaselabels<-unique(sorted_total$disease)

agon_plot_counts <- c()
antagon_plot_counts <- c()
for (i in 1:length(diseaselabels)){
  if(diseaselabels[i] %in%  names(table_antagon)){
    antagon_plot_counts[i] <- table_antagon[names(table_antagon) == diseaselabels[i]]
  }
  else{
    antagon_plot_counts[i] <- 0
    names(antagon_plot_counts[i]) <- diseaselabels[i]
  }
  result2 <- try(table_agon[names(table_agon) == diseaselabels[i]])
  if(diseaselabels[i] %in% names(table_agon)){
    agon_plot_counts[i] <- table_agon[names(table_agon) == diseaselabels[i]]
  }
  else{
    agon_plot_counts[i] <- 0
    names(agon_plot_counts[i]) <- diseaselabels[i]
  }}

final_matrix <- as.data.frame(cbind(as.character(diseaselabels), agon_plot_counts, antagon_plot_counts))




pyramid.plot(as.numeric(as.character(final_matrix$agon_plot_counts)),
             as.numeric(as.character(final_matrix$antagon_plot_counts)),
             labels=(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)

pyramid.plot(rev(as.numeric(as.character(final_matrix$agon_plot_counts))),
             rev(as.numeric(as.character(final_matrix$antagon_plot_counts))),
             labels=rev(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)











############################
############################
#CURATED PLOT OF DOMAINS#
###########################

library(ggplot2)
library(gridExtra)
library("plotrix")
library("RColorBrewer")
library(tidyverse)

paired_dataset <- read.csv("freqs_Agon_Antagon_curated.tsv", header = TRUE, sep = "\t")
summary(paired_dataset)
sum(colSums(paired_dataset[,c("Agonistic","Antagonistic")]))

#main barplot
paired_dataset.sort <-paired_dataset[order(rowSums(paired_dataset[,c("Antagonistic", "Agonistic")]),decreasing=T),]
data_matrix <- t(as.matrix(paired_dataset.sort[,c("Agonistic", "Antagonistic")]))
barplot(data_matrix, names.arg = paired_dataset.sort$Domains, beside = TRUE, 
        legend.text = c("Agonistic", "Antagonistic"), col = c("lavender", "gold"))



####
#ALL DISEASES PLOTS
#######
Agonistic <- read.csv("Agonistic_curated", header = FALSE, sep = "\t")
Antagonistic <-read.csv("Antagonistic_curated", header = FALSE, sep = "\t")
summary(Agonistic)
summary(Antagonistic)

#read onsets
Onsets <- read.csv("onsets.txt", header = FALSE, sep = "\t")
summary(Onsets)


all_Onsets<- c(as.numeric(as.character(Onsets$V2)), as.numeric(as.character(Onsets$V4)))
all_Onsets <- all_Onsets[ !is.na(all_Onsets) ]
all_diseases <- c(as.character(Onsets$V1), as.character(Onsets$V3))
elements_2_remove <- c("")
all_diseases <- all_diseases[!(all_diseases %in% elements_2_remove)]
diseases_onsets <- as.data.frame(cbind(all_Onsets, all_diseases))




#agonistic filters

agon_onsets <- c()
agon_diseases <- c(as.character(Agonistic$V2), as.character(Agonistic$V10))
elements_2_remove <- c("")
agon_diseases <- agon_diseases[!(agon_diseases %in% elements_2_remove)]
for (i in 1:length(agon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == agon_diseases[i]),]
  agon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

agon_total <- as.data.frame(cbind(agon_onsets, agon_diseases))
table_agon <- table(agon_total$agon_diseases)

#antagonistic filters

antagon_onsets <- c()
antagon_diseases <- c(as.character(Antagonistic$V2), as.character(Antagonistic$V10))
elements_2_remove <- c("")
antagon_diseases <- antagon_diseases[!(antagon_diseases %in% elements_2_remove)]
for (i in 1:length(antagon_diseases)){
  subset <- diseases_onsets[which(diseases_onsets$all_diseases == antagon_diseases[i]),]
  antagon_onsets[i] <- as.numeric(as.character(subset$all_Onsets[1]))
}

antagon_total <- as.data.frame(cbind(antagon_onsets, antagon_diseases))
table_antagon <- table(antagon_total$antagon_diseases)

#All together

agon_total$type <- "Agonist"
colnames(agon_total) <- c("onset", "disease", "type")
antagon_total$type <- "Antagonist"
colnames(antagon_total) <- c("onset", "disease", "type")

total <- rbind(agon_total, antagon_total)

#PYRAMID
sorted_total <- total[order(as.numeric(as.character(total$onset))),]
diseaselabels<-unique(sorted_total$disease)

agon_plot_counts <- c()
antagon_plot_counts <- c()
for (i in 1:length(diseaselabels)){
  if(diseaselabels[i] %in%  names(table_antagon)){
    antagon_plot_counts[i] <- table_antagon[names(table_antagon) == diseaselabels[i]]
  }
  else{
    antagon_plot_counts[i] <- 0
    names(antagon_plot_counts[i]) <- diseaselabels[i]
  }
  result2 <- try(table_agon[names(table_agon) == diseaselabels[i]])
  if(diseaselabels[i] %in% names(table_agon)){
    agon_plot_counts[i] <- table_agon[names(table_agon) == diseaselabels[i]]
  }
  else{
    agon_plot_counts[i] <- 0
    names(agon_plot_counts[i]) <- diseaselabels[i]
  }}

final_matrix <- as.data.frame(cbind(as.character(diseaselabels), agon_plot_counts, antagon_plot_counts))




pyramid.plot(as.numeric(as.character(final_matrix$agon_plot_counts)),
             as.numeric(as.character(final_matrix$antagon_plot_counts)),
             labels=(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)

pyramid.plot(rev(as.numeric(as.character(final_matrix$agon_plot_counts))),
             rev(as.numeric(as.character(final_matrix$antagon_plot_counts))),
             labels=rev(as.character(final_matrix$V1)), lxcol=c("lavender"),
             rxcol=c("gold"),laxlab=c(0,5,10,20),
             raxlab=c(0,5,10,20),top.labels=c("Agonist","Disease","Antagonist"),gap=20, space=0.5, labelcex=1)
