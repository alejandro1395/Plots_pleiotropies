##########################################
##########################################
###RScript to analize location of pleios#
###########################################

#INPUT FILES

chr_coord <- read.csv("annotate_file.tsv", header = FALSE, sep = "\t")
View(chr_coord)

hum_chr_coord <- read.csv("hum_chr_coord.tsv", header = FALSE, sep = "\t")
View(hum_chr_coord)

library(chromoMap)
chromoMap("hum_chr_coord.tsv","annotate_file.tsv",
          data_based_color_map = T,
          data_type = "numeric")

require(chromoMap)
chromoMap("hum_chr_coord.tsv","annotate_file.tsv",
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("blue","yellow")),
          legend = TRUE)
