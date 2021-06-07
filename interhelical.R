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

library(ggplot2)
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
iterate_files <- function(input_folder, Gprot = TRUE){
  helical_distances <- matrix(ncol = 4, nrow = 0)
  colnames(helical_distances) <- c("receptor","partner","TM3-TM6","TM3-TM7")
  for(files in list.files(input_folder)){
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
      if ((substr(partner,1,1) == "G") & (Gprot == TRUE)){
        helical_distances <- rbind(helical_distances, row)
      }
      else if((substr(partner,1,1) != "G") & (Gprot == FALSE)){
        helical_distances <- rbind(helical_distances, row)
      }
      }
  }
  return(helical_distances)
}

results <- "processed_results"
receptors_list <- c("DOR","KOR","MOR","NOP")
color_index <- as.character(c("#1a1aff","#ff9933","#33ff33","#cc6666"))

colors_vec_1 <- c(rep("#1a1aff", 14),rep("#ff9933", 10),rep("#33ff33", 14),rep("#cc6666", 14))
results_df_gprot <- data.frame(iterate_files(results, Gprot = TRUE))
inter_dist_gprot <- ggplot(data = results_df_gprot, aes(x=as.numeric(as.character(TM3.TM6)), y=as.numeric(as.character(TM3.TM7)), label = partner, color = receptor)) +
  geom_point(size = 8, alpha = 1/4, stat = "identity") + scale_color_manual(breaks =  results_df_gprot$receptor, values = colors_vec_1) + 
  geom_text(hjust = 0, vjust=0, show.legend = FALSE,color = colors_vec_1) +
  scale_x_continuous(limits = c(13, 18.5)) +
  theme_minimal() +
  xlab("TM3-TM6 distance") + ylab("TM3-TM7 distance")
write.csv(results_df_gprot, file = "summary/inter_dist_gprot")
ggsave(file="images/interhelical_gprot.svg", plot=inter_dist_gprot, width=10, height=8)

colors_vec_2 <- c(rep("#1a1aff", 4),rep("#ff9933", 4),rep("#33ff33", 4),rep("#cc6666", 4))
results_df_arr <- data.frame(iterate_files(results, Gprot = FALSE))
results_df_arr$partner <- gsub("ARR", "Arr", results_df_arr$partner)
inter_dist_arr <- ggplot(data = results_df_arr, aes(x=as.numeric(as.character(TM3.TM6)), y=as.numeric(as.character(TM3.TM7)), label = partner, color = receptor)) +
  geom_point(size = 8, alpha = 1/4, stat = "identity") + scale_color_manual(breaks =  results_df_arr$receptor, values = colors_vec_2) + 
  geom_text(hjust=0, vjust=0, show.legend = FALSE,color = colors_vec_2) +
  scale_x_continuous(limits = c(13, 18.5)) +
  theme_minimal() +
  xlab("TM3-TM6 distance") + ylab("TM3-TM7 distance")
write.csv(results_df_arr, file = "summary/inter_dist_arr")
ggsave(file="images/interhelical_arr.svg", plot = inter_dist_arr, width=10, height=8)