############################
############################
#READ THE TABLE WITH PAIRED#
###########################

library(ggplot2)
library(ggsignif)

single_dataset <- read.csv("SinglePleiosLoc.tsv", header = TRUE, sep = "\t")
summary(single_dataset)

#SUM TOTAL OF DOMAIN COLUMNS

all_domains <- sort(colSums(single_dataset[,c(6:19)]),decreasing=TRUE)
elements_2_remove <- c("")
all_domains = all_domains[!(all_domains %in% elements_2_remove)]
barplot(all_domains, col="darkgreen")

#Post-Hoc test for obtaining significant asterisks for ggplot

#In random panorama:

expected_domains <- rep(sum(all_domains)/length(all_domains), length(all_domains))

pvalues <- c()
observed_posthoc<- c()
expected_posthoc <- c()
for (i in 1:length(all_domains)){
  observed_posthoc <- c(all_domains[i], sum(all_domains[-i]))
  expected_posthoc <- c(expected_domains[i], sum(expected_domains[-i]))
  pv <- as.double(chisq.test(observed_posthoc, p=expected_posthoc/sum(all_domains))["p.value"])
  pvalues<- append(pvalues, pv)
}

#adjust FDR
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))


#ggplot version
histogram_domains <- as.data.frame(cbind(names(all_domains),all_domains, corrected_pvalues))
ggplot(data=histogram_domains, 
       aes(x=reorder(V1, -as.numeric(as.character(all_domains))), 
           y=as.numeric(as.character(all_domains)), 
           label = ifelse(as.numeric(as.character(corrected_pvalues)) < 0.05, "*", ""))) + 
  geom_bar(stat = "identity", fill="darkgreen") + xlab("Disease domains") +
  ylab("Pleiotropies count") +
  geom_text(vjust = 0) +
  theme_light()
