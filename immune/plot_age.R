##########################################
##########################################
###RScript to analize localization of pleios#
###########################################

#INPUT FILES

Ages_varx <- seq(10, 60)
length(Ages_varx)
pvalues <- c()

for (i in 10:60){
  folder <- paste("Age_threeshold_", i, sep="")
  path <- paste(folder, "pleios_count", sep="/")
  file <- read.csv(path, header = FALSE, sep = "\t")
  colnames(file) <- c("Agonistic", "Antagonistic")
  rownames(file) <- c("Early_Early|Late_Late", "Early_Late")
  pv <- as.double(chisq.test(file)["p.value"])
  pvalues<- append(pvalues, pv)
}

pvalues

length(pvalues)
Sign_vary <- -log10(pvalues)
length(Sign_vary)
plot(x=Ages_varx, y=Sign_vary) +
abline(h=-log10(0.05), col="red", lty=2)
title(main= "Immune pleiotropies excess vs age ") 
rect(xleft = 35, ybottom = 0, xright = 58, ytop = 5, density = NULL, angle = 45,
     col= rgb(0,0,1.0,alpha=0.08), border = NULL)

