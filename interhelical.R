#install.packages("ggplot2")
#install.packages("svglite")
###----------------------------------------------------------###
###Script to generate scatteplots for interhelical distances-###
###----------------------------------------------------------###

###__author__ = "A.J. Preto"---------------------------------###
###__email__ = "martinsgomes.jose@gmail.com"-----------------###
###__group__ = "Data-Driven Molecular Design"----------------###
###__group_leader__ = "Irina S. Moreira"---------------------###
###__project__ = "GPCRs"-------------------------------------###

setwd("C:/Users/marti/OneDrive/Desktop/silverio")
change_header <- function(input_cols){
  col_names <- c()
  for (name in input_cols){
    col_names <- c(col_names, substr(name,1,nchar(name)-3))
  }
  return(col_names)
}
TM_distance_calc <- function(input_table){
  table <- read.csv(input_table, sep = ";", row.names = 1)
  colnames(table) <- change_header(colnames(table))
  T3_T6 <- c("3.50","X6.30")
  T3_T7 <- c("3.50","X7.53")
  distance_T3_T6 <- table[T3_T6[1],T3_T6[2]]
  distance_T3_T7 <- table[T3_T7[1],T3_T7[2]]
  return(c(distance_T3_T6, distance_T3_T7))
}
iterate_files <- function(input_folder){
  helical_distances <- matrix(ncol = 4, nrow = 0)
  colnames(helical_distances) <- c("receptor","partner","TM3-TM6","TM3-TM7")
  for(files in list.files(results)){
    if (startsWith(files,"weinstein_intra_chain") == TRUE){
      complex_name <- strsplit(files,"_")
      if(length(complex_name[[1]]) == 5){
        only_complex <- paste(complex_name[[1]][4],complex_name[[1]][5],sep = "_")
        }
      if(length(complex_name[[1]]) == 4){
        only_complex <- complex_name[[1]][4]
        }
      final_complex <- substr(only_complex, 1, nchar(only_complex) - 4)
      receptor <- strsplit(final_complex, "-")[[1]][1]
      partner <- strsplit(final_complex, "-")[[1]][2]
      path_name <- paste(results, files,sep = "/")
      TM_distances <- TM_distance_calc(path_name)
      row <- c(receptor, partner, TM_distances)
      helical_distances <- rbind(helical_distances, row)}
  }
  return(helical_distances)
}

results <- "processed_results"
colors_vec <- c(rep("#1a1aff", 14),rep("#ff9933", 10),rep("#33ff33", 14),rep("#cc6666", 14))
receptors_list <- c("DOR","KOR","MOR","NOP")
color_index <- as.character(c("#1a1aff","#ff9933","#33ff33","#cc6666"))
results_df <- data.frame(iterate_files(results))
results_df$partner <- gsub("ARR", "Arr", results_df$partner)
results_df
inter_dis <- ggplot(data = results_df, aes(x=as.numeric(as.character(TM3.TM6)), y=as.numeric(as.character(TM3.TM7)), label = partner, color = receptor)) +
  geom_point(size = 8, alpha = 1/4, stat = "identity") + scale_color_manual(breaks =  results_df$receptor, values = unique(as.character(color_index))) + 
  geom_text(hjust=0, vjust=0, show.legend = FALSE,color = colors_vec) +
  scale_x_continuous(limits = c(13, 18.5)) +
  theme_minimal() +
  xlab("TM3-TM6 distance") + ylab("TM3-TM7 distance")
write.csv(results_df, file = "summary/inter_dist_gprot")
ggsave(file="images/interhelical_gprot.svg", plot=inter_dis, width=10, height=8)