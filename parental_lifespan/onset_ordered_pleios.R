############################
############################
#READ THE TABLE WITH PAIRED#
###########################

library(ggplot2)
library(gridExtra)
library(plotrix)
library(reshape2)

paired_dataset <- read.csv("full_dataset_paired.csv", header = FALSE, sep = "\t")
summary(paired_dataset)
View(paired_dataset)

#pleiotropies below onset 40 and above onset 40

all_Onsets<- c(as.numeric(as.character(paired_dataset$V5)), as.numeric(as.character(paired_dataset$V12)))
all_Onsets <- all_Onsets[ !is.na(all_Onsets) ]
all_domains<- c(as.character(paired_dataset$V3), as.character(paired_dataset$V10))
elements_2_remove <- c("")
all_domains <- all_domains[!(all_domains %in% elements_2_remove)]

#LOOP for thresholds between 40 and 50
par(mfrow = c(5, 4), mar=c(3,3,3,3))


for (i in 40:59){
x <- c(paste("Diseases < ",i,"y"), paste("Diseases >= ",i,"y"))
y <- c(length(all_Onsets[all_Onsets<i]), length(all_Onsets[all_Onsets>i]))
names(y) <- x
barplot(y,
        xlab="Groups", col=c("lightgreen","darkred"))
}

#boxplot
qqnorm(table(all_Onsets, all_domains))

par(mfrow = c(1,1))
boxplot(all_Onsets~all_domains, col="gold")

#histogram withb ggplot of two categories

data_onsets <- as.data.frame(cbind(all_domains, all_Onsets))

plot_list = list()
for (i in 40:59){
data_onsets$all_filtered_onsets <- ifelse(as.numeric(as.character(data_onsets$all_Onsets)) < i, 
                                          paste("Disease < ",i,"y"), paste("Disease >=",i,"y"))
plot <- ggplot() +
  geom_bar(data=data_onsets
                 , aes(x=all_filtered_onsets, fill=all_filtered_onsets)) +
  scale_fill_brewer(palette = 5) +
  facet_wrap(~all_domains) + 
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot_list[[i-39]] = plot
}

do.call(grid.arrange, plot_list)

#All diseases ordered by Onset (sorted by Onset)
all_Onsets<- c(as.numeric(as.character(paired_dataset$V5)), as.numeric(as.character(paired_dataset$V12)))
all_Onsets <- all_Onsets[ !is.na(all_Onsets) ]
all_diseases <- c(as.character(paired_dataset$V2), as.character(paired_dataset$V9))
elements_2_remove <- c("")
all_diseases <- all_diseases[!(all_diseases %in% elements_2_remove)]


diseases_onsets <- as.data.frame(cbind(all_Onsets, all_diseases))
diseases_onsets <- diseases_onsets[order(all_Onsets),]
total_diseases <- table(diseases_onsets$all_diseases)

expected_diseases <- rep(sum(total_diseases)/length(total_diseases), length(total_diseases))

pvalues <- c()
observed_posthoc<- c()
expected_posthoc <- c()
for (i in 1:nrow(total_diseases)){
  observed_posthoc <- c(total_diseases[i], sum(total_diseases[-i]))
  expected_posthoc <- c(expected_diseases[i], sum(expected_diseases[-i]))
  pv <- as.double(chisq.test(observed_posthoc, p=expected_posthoc/sum(total_diseases))["p.value"])
  pvalues<- append(pvalues, pv)
}

#adjust FDR
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))
diseases_onsets$corrected_pvalues <- corrected_pvalues[diseases_onsets$all_diseases]
counts <- tapply(diseases_onsets$all_diseases,diseases_onsets$all_diseases,length)
diseases_onsets$Freq <- counts[diseases_onsets$all_diseases]
unique_diseases <- unique(diseases_onsets)

ggplot(data=unique_diseases, aes(x=reorder(all_diseases, -as.numeric(as.character(all_Onsets))),
                                 y=as.numeric(as.character(Freq)),
                                 fill = "gold", 
                                 label =  ifelse(as.numeric(as.character(corrected_pvalues)) < 0.05, "*", ""))) +
  geom_bar(stat="identity") + xlab("Decreasing Onset diseases") +
  scale_colour_gradient2() +
  coord_flip()+
  geom_text(vjust = 0.5, hjust = 0.0001) +
  theme_classic() + theme(legend.position = "none") + ylab("Disease Counts") + xlab("Ordered diseases by increasing Onser")


