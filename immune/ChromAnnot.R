##########################################
##########################################
###RScript to analize location of pleios#
###########################################

chr_coord <- read.csv("annotate_file.tsv", header = FALSE, sep = "\t")
View(chr_coord)

library(chromoMap)
chromoMap("hum_chr_coord.tsv","annotate_file.tsv",
          data_based_color_map = T,
          data_type = "numeric")

library(chromoMap)
chromoMap("hum_chr_coord.tsv","annotate_file.tsv", 
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("blue","yellow")),
          legend = TRUE,
          lg_x = 50)
          
          
          
          
          canvas_width = 1000, canvas_height = 2000,
          top_margin = 100,
          left_margin = 150,
          lg_x = 0,
          lg_y = 0)
