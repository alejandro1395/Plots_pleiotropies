##########################################
##########################################
###RScript to make disease DOMAINS study#################
###########################################

#IMPORT libraries

library(plyr)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(ggplot2)
library(qqman)
library(qualityTools)
library(Haplin)
library(gridBase)
library(tidyverse)
library(car)
library(wesanderson)
library(reshape2)
library("gplots")
library(Hmisc)

###############
###############
#PHENOTIPIC DB#
###############

Freqs_table <- read.csv("freqs_total.tsv", header = TRUE, sep = "\t")
View(Freqs_table)


#adjust FDR
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))


#ggplot version
histogram_domains <- as.data.frame(cbind(names(combinations),combinations, corrected_pvalues))
ggplot(data=histogram_domains, 
       aes(x=reorder(V1, -as.numeric(as.character(combinations))), 
           y=as.numeric(as.character(combinations)), 
           label = ifelse(as.numeric(as.character(corrected_pvalues)) < 0.05, "*", ""))) + 
  geom_bar(stat = "identity", fill="darkred") + xlab("Disease domains") +
  ylab("Pleiotropies count") +
  geom_text(vjust = 0) +
  theme_light()















#Test of goodness of fit

scalar1 <- function(x) {x / sqrt(sum(x^2))}
vect_freqs_norm <- Freqs_table$Exp_prob/sum(Freqs_table$Exp_prob)
Freqs_table$Norm_exp_prob <- vect_freqs_norm

res <- chisq.test(Freqs_table$Total_observed, p = Freqs_table$Norm_exp_prob)
res$expected

Table_pairs <- data.frame(Combinations=Freqs_table$Combined_diseases,
                          Observed=Freqs_table$Total_observed, 
                          Expected=as.numeric((Freqs_table$Norm_exp_prob)*sum(Freqs_table$Total_observed)), 
                          stringsAsFactors=FALSE)
options(scipen=999)
View(Table_pairs)


#PLOT DATAFRAME
library(gridExtra)
library(grid)

par(mfrow=c(2,4))


#THEME
myt <- ttheme_default(base_size = 8,
                      # Use hjust and x to left justify the text
                      # Alternate the row fill colours
                      core = list(fg_params=list(hjust = 1, x=1),
                                  bg_params=list(fill=c("white", "grey"))),
                      
                      # Change column header to white text and red background
                      colhead = list(fg_params=list(col="white"),
                                     bg_params=list(fill="black"))
)

grid.newpage()
grid.draw(tableGrob(format(Table_pairs[1:15, ], big.mark=","), theme=myt  ,rows=NULL))
grid.draw(tableGrob(format(Table_pairs[16:30, ], big.mark=","), theme=myt  ,rows=NULL))
grid.draw(tableGrob(format(Table_pairs[31:45, ], big.mark=","), theme=myt  ,rows=NULL))
grid.draw(tableGrob(format(Table_pairs[46:60, ], big.mark=","), theme=myt  ,rows=NULL))
grid.draw(tableGrob(format(Table_pairs[61:75, ], big.mark=","), theme=myt  ,rows=NULL))
grid.draw(tableGrob(format(Table_pairs[76:90, ], big.mark=","), theme=myt  ,rows=NULL))
grid.draw(tableGrob(format(Table_pairs[91:105, ], big.mark=","), theme=myt  ,rows=NULL))

#Post-Hoc test

#Test for Normality (Shapiro Test)

shapiro.test(Freqs_table$Total_observed)
Freqs_table <- Freqs_table[order(Freqs_table$Total_observed, decreasing = TRUE),]


pvalues <- c()
for (i in 1:nrow(Freqs_table)){
  observed_posthoc <- c(Freqs_table$Total_observed[i], sum(Freqs_table$Total_observed[-i]))
  expected_posthoc <- c(Freqs_table$Norm_exp_prob[i], sum(Freqs_table$Norm_exp_prob[-i]))
  pv <- as.double(chisq.test(observed_posthoc, p=expected_posthoc)["p.value"])
  pvalues<- append(pvalues, pv)
  print(Freqs_table$Combined_diseases[i])
  print(chisq.test(observed_posthoc, p=expected_posthoc))
}

#adjust FDR
names(pvalues) <- as.character(Freqs_table$Combined_diseases)
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))
Freqs_table$Corrected_pvalues <- corrected_pvalues
melted_Freqs <- as.data.frame(cbind(as.character(Freqs_table$Combined_diseases), Freqs_table$Total_observed, 
                                    Freqs_table$Norm_exp_prob*sum(Freqs_table$Total_observed), Freqs_table$Corrected_pvalues))
colnames(melted_Freqs) <- c("Combined_diseases", "Total_Observed", "Expected_count", "Corrected_pvalues")
melted_Freqs <- melted_Freqs[which(as.numeric(as.character(melted_Freqs$Total_Observed))>0), ]
melted_Freqs <- melt(melted_Freqs, id=c("Combined_diseases", "Corrected_pvalues"))
colnames(melted_Freqs) <- c("Combined_diseases", "Corrected_pvalues", "Type", "Count")
melted_Freqs$Pos <- rep(as.numeric(as.character(melted_Freqs$Count[as.character(melted_Freqs$Type) == "Total_Observed"])), 2)

##ggplot bars with stars

ggplot(data=melted_Freqs,
       aes(x=reorder(Combined_diseases, -as.numeric(as.character(Count))), 
           y=as.numeric(as.character(Count)), 
           #label = ifelse(as.numeric(as.character(Corrected_pvalues)) < 0.05, "*", ""),
           fill = as.factor(Type))) + 
  geom_bar(stat = "identity", position = "dodge") + xlab("Disease domain combinations") +
  ylab("Pleiotropies count") +
  #geom_text(vjust = 0) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Count type") +
  scale_fill_manual(values=c("darkred", "red")) +
  geom_text(vjust = 0, aes(y = Pos + 4,
                           label = ifelse(as.numeric(as.character(melted_Freqs$Corrected_pvalues)) <= 0.0005, "***", 
                                          ifelse(as.numeric(as.character(melted_Freqs$Corrected_pvalues)) <= 0.05, "*", ""))))



sorted_corrected_pvalues <- sort(corrected_pvalues)
pleios_df <- as.data.frame(cbind(names(sorted_corrected_pvalues), as.numeric(sorted_corrected_pvalues)))
colnames(pleios_df) <- c("Combined_diseases", "Corrected_pval")
View(pleios_df) 

#histogram
par(mfrow=c(1,1))
ggplot(data=pleios_df[1:25,], aes(x= reorder(Combined_diseases, as.numeric(as.character(Corrected_pval))), 
                                  y=round(as.numeric(as.character(Corrected_pval)), 3), fill=Combined_diseases), 
       fill=Combined_diseases) + scale_fill_hue(l=10, c=45) +
  geom_bar(stat="identity", width=0.2) + guides(fill=FALSE) +
  theme_minimal() + theme(aspect.ratio = 0.9/1, 
                          axis.text.x = element_text(angle = 40, hjust = 1)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") + ylim(0, 1) + ylab("Corrected_p_values") +
  xlab("Combined diseases")


#heatmap


pleios_df[c('domainA','domainB')] <- colsplit(pleios_df$Combined_diseases,'_',c('domainA','domainB'))
pleios_matrix <- as.data.frame(cbind(pleios_df$domainA, pleios_df$domainB, 
                                     as.numeric(as.character(pleios_df$Corrected_pval))), stringsAsFactors = FALSE)
subset_repeated <- as.data.frame(cbind(pleios_df$domainB, pleios_df$domainA, as.numeric(as.character(pleios_df$Corrected_pval))))
subset_repeated <- subset_repeated[which(as.character(subset_repeated$V1) != as.character(subset_repeated$V2)),]
pleios_matrix_df <- rbind(pleios_matrix, subset_repeated)
pleios_matrix_df <- pleios_matrix_df[order(pleios_matrix_df$V3, decreasing = TRUE),]
pleios_matrix <- acast(pleios_matrix_df, 
                      V1~V2, value.var="V3")
pleios_matrix[is.na(pleios_matrix)] = 1
mode(pleios_matrix) = "numeric"












############################################
#######FISHER PRESENCE/ABSENCE APPROACH#####
############################################

#We import both datasets
paired <- read.csv("dataset_paired.csv",header=TRUE,sep="\t")
x<- read.csv("SinglePleiosLoc.tsv",header=TRUE,sep="\t")
single <- x[,c("skin", "skeletal", "respiratory", "reproductive", "cancer", "immune", "infectious", "kidney", "cardio", "metabolic",
               "headneck", "nervous", "eye", "endocrine")]
column_levels <- c("skin","skeletal","respiratory","reproductive","cancer","immune","infectious","kidney","cardio","metabolic",
                   "headneck","nervous","eye","endocrine")

summary(paired)
summary(single)


#LOOP AND CONCATENATE TO DO FISHER test

result = NULL
for(i in 1:(ncol(single))){
  for(j in i:ncol(single)){
    if (all(c(column_levels[i], column_levels[j]) == c("skin", "skin")) | 
        all(c(column_levels[i], column_levels[j]) == c("skin", "skeletal"))) {
      next
    }
    else if (column_levels[i] == column_levels[j]){
      paired_columnA <- ifelse(paired$DomainA == column_levels[i], 1, 0)
      paired_columnB <- ifelse(paired$DomainB == column_levels[j], 1, 0)}
    else{
      paired_columnA <- ifelse(paired$DomainA == column_levels[i] | paired$DomainB == column_levels[i], 1, 0)
      paired_columnB <- ifelse(paired$DomainA == column_levels[j] | paired$DomainB == column_levels[j], 1, 0)}
    print(column_levels[i])
    print(column_levels[j])
    print(summary(paired_columnA))
    print(summary(single[,i]))
    columnA <- c(single[,i], paired_columnA)
    columnB <- c(single[,j], paired_columnB)
    ft = fisher.test(columnA,columnB)
    result = rbind(result,c(column_levels[i],column_levels[j],ft$p.value,ft$estimate,ft$conf.int))
  }
}

result = as.data.frame(result)
names(result) = c("group1","group2","fisher.p","or","lci.or","uci.or")
result$domains <- paste(result$group1, result$group2, sep="_")
result$fisher.corrected.p <- p.adjust(as.numeric(as.character(result$fisher.p)), method = "bonferroni", 
                                      n = length(result$fisher.p))
result

#histogram ODDS RATIO
par(mfrow=c(1,1))
ggplot(data=result[round(as.numeric(as.character(result$or))) > 1,], 
       aes(x= reorder(domains, as.numeric(as.character(or))), 
           y=round(as.numeric(as.character(or)), 3), fill=domains), 
       fill=Combined_diseases) + scale_fill_hue(l=10, c=45) +
  geom_bar(stat="identity", width=0.2) + guides(fill=FALSE) +
  theme_minimal() + theme(aspect.ratio = 0.9/1, 
                          axis.text.x = element_text(angle = 40, hjust = 1)) + ylab("ODDs ratio") +
  xlab("Combined_domains")

#histogram P-VALUES
ggplot(data=result[as.numeric(as.character(result$fisher.corrected.p)) < 0.05,], 
       aes(x= reorder(domains, as.numeric(as.character(fisher.corrected.p))), 
           y=as.numeric(as.character(fisher.corrected.p)), fill=domains), 
       fill=Combined_diseases) + scale_fill_hue(l=10, c=45) +
  geom_bar(stat="identity", width=0.2) + guides(fill=FALSE) +
  theme_minimal() + theme(aspect.ratio = 0.9/1, 
                          axis.text.x = element_text(angle = 40, hjust = 1)) + ylab("corrected p-value") +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
  xlab("Combined_domains")











######################################
######################################
######################################
####comorbidity data from FAES########
######################################



Comorbidities_table <- read.csv("FAERS_net.csv", header = TRUE, sep = "\t")
View(Comorbidities_table)

#Validation heatmap

subset_comorb <- as.data.frame(cbind(as.character(Comorbidities_table$dis1_SOC), 
                                     as.character(Comorbidities_table$dis2_SOC)))
colnames(subset_comorb) <- c("domainA", "domainB")
View(subset_comorb)

freqs_comorb <- count(subset_comorb, vars = c("domainA", "domainB"))
filtered_freqs_comorb <- freqs_comorb[!(freqs_comorb$domainA=="null" | freqs_comorb$domainB=="null"),]
View(filtered_freqs_comorb)

filtered_freqs_comorb <- filtered_freqs_comorb[order(filtered_freqs_comorb$freq, decreasing = TRUE),]
comorb_matrix <- acast(filtered_freqs_comorb, 
                       domainA~domainB, value.var="freq")
comorb_matrix[is.na(comorb_matrix)] = 1
mode(comorb_matrix) = "numeric"
sort

my_palette <- colorRampPalette(c("lightyellow", "yellow","orange", "darkred"))
heatmap.2(-comorb_matrix, Rowv=FALSE, trace="none", col = heat.colors(5), scale='none', symm =T,
          cexRow=0.7, cexCol=0.7, margins = c(17,17), srtRow = 0, srtCol = 45, offsetRow = 0, offsetCol = 0)




