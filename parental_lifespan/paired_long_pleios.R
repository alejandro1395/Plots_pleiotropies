############################
############################
#READ THE TABLE WITH PAIRED#
###########################

library(ggplot2)
library(ggsignif)

paired_dataset <- read.csv("full_dataset_paired.csv", header = FALSE, sep = "\t")
summary(paired_dataset)

#concatenate domains on both sides

all_domains <- c(as.character(paired_dataset$V3), as.character(paired_dataset$V10))
elements_2_remove <- c("")
all_domains = all_domains[!(all_domains %in% elements_2_remove)]
tab <- table(all_domains)
barplot(sort(table(all_domains), decreasing = TRUE), col="darkgreen")


#Post-Hoc test for obtaining significant asterisks for ggplot

#In random panorama:

ordered_domains <- sort(table(all_domains), decreasing=TRUE)
expected_domains <- rep(sum(ordered_domains)/length(ordered_domains), length(ordered_domains))

pvalues <- c()
observed_posthoc<- c()
expected_posthoc <- c()
for (i in 1:nrow(ordered_domains)){
  observed_posthoc <- c(ordered_domains[i], sum(ordered_domains[-i]))
  expected_posthoc <- c(expected_domains[i], sum(expected_domains[-i]))
  pv <- as.double(chisq.test(observed_posthoc, p=expected_posthoc/sum(ordered_domains))["p.value"])
  pvalues<- append(pvalues, pv)
}

#adjust FDR
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))


#ggplot version
histogram_domains <- as.data.frame(cbind(names(ordered_domains),ordered_domains, corrected_pvalues))
ggplot(data=histogram_domains, 
       aes(x=reorder(V1, -as.numeric(as.character(ordered_domains))), 
           y=as.numeric(as.character(ordered_domains)), 
           label = ifelse(as.numeric(as.character(corrected_pvalues)) < 0.05, "*", ""))) + 
         geom_bar(stat = "identity", fill="darkgreen") + xlab("Disease domains") +
  ylab("Pleiotropies count") +
  geom_text(vjust = 0) +
  theme_light()
  





