##########################################
##########################################
###RScript to analize location of pleios#
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
plot(x=Ages_varx, y=Sign_vary)
