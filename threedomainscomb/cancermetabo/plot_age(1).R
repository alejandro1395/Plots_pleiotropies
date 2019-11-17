###########################################
###########################################
###RScript to analize location of pleios###
###########################################

#INPUT FILES

Ages_varx <- seq(10, 60)
length(Ages_varx)
pvalues <- c()

#I use fisher test because it is designed for small values in cells, small Ns
for (i in 10:60){
  folder <- paste("Age_threeshold_", i, sep="")
  path <- paste(folder, "pleios_count", sep="/")
  file <- read.csv(path, header = FALSE, sep = "\t")
  probs <- file$V2
  observed <- file$V1
  total <- sum(file$V1)
  expected <- probs*total
  table <- cbind(observed, expected)
  #colnames(file) <- c("Agonistic", "Antagonistic")
  #rownames(file) <- c("Early_Early|Late_Late", "Early_Late")
  pv <- as.double(chisq.test(table)["p.value"])
  pvalues<- append(pvalues, pv)
}

pvalues
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))


length(pvalues)
Sign_vary <- -log10(corrected_pvalues)
length(Sign_vary)
plot(x=Ages_varx, y=Sign_vary) +
abline(h=-log10(0.05), col="red", lty=2)
title(main= "CancerMetabolic combined AP Early-Late vs age ") 
#rect(xleft = 33, ybottom = 0, xright = 60, ytop = 8, density = NULL, angle = 45,
    # col= rgb(0,0,1.0,alpha=0.08), border = NULL)

#Next_Step
smoothingSpline = smooth.spline(Ages_varx, Sign_vary, spar=0.55)
lines(smoothingSpline)

