###----------------------------------------------------------###
###Script to generate circular graphs from distance matrices-###
###----------------------------------------------------------###

#install.packages("circlize")
library(circlize)
library(stringr)
###__author__ = "A.J. Preto"---------------------------------###
###__email__ = "martinsgomes.jose@gmail.com"-----------------###
###__group__ = "Data-Driven Molecular Design"----------------###
###__group_leader__ = "Irina S. Moreira"---------------------###
###__project__ = "GPCRs"-------------------------------------###

circle_graph_builder <- function(distance,file_name,len,TM1,ICL1,TM2,TM3,ICL2,TM4,TM5,ICL3,TM6,TM7,H8,partner_info,partner_name){
  
  ###---Loading the table---###
  df = read.table(file_name, sep=';',header = TRUE, stringsAsFactors = FALSE, row.names=1, check.names = FALSE)
  df <-data.matrix(df)
  cut_off_list <- list()
  
  ###---Turn table values into binaries depending on distance cutoff value---###
  df <- ifelse(df<distance,1,df)
  df <- ifelse(df>distance,0,df)
  
  ###---Assign colors by cutoff value of distance---#
  colorizer <- function(TM1,ICL1,TM2,TM3,ICL2,TM4,TM5,ICL3,TM6,TM7,H8){
    
    color_vec <- c()
    color_vec <- c(color_vec,replicate(TM1[1] - 1,'cornsilk'))
    color_vec <- c(color_vec,replicate(TM1[2]-TM1[1],'#FFFF00'))
    color_vec <- c(color_vec,replicate(ICL1[2]-ICL1[1],'#C65911'))
    color_vec <- c(color_vec,replicate(TM2[2]-TM2[1],'#FFC000'))
    color_vec <- c(color_vec,replicate(TM3[1]-TM2[2],'cornsilk'))
    color_vec <- c(color_vec,replicate(TM3[2]-TM3[1],'#92D050'))
    color_vec <- c(color_vec,replicate(ICL2[2]-ICL2[1],'#375623'))
    color_vec <- c(color_vec,replicate(TM4[2]-TM4[1],'#00B050'))
    color_vec <- c(color_vec,replicate(TM5[1]-TM4[2],'cornsilk'))
    color_vec <- c(color_vec,replicate(TM5[2]-TM5[1],'#00B0F0'))
    color_vec <- c(color_vec,replicate(ICL3[2]-ICL3[1],'#333F4F'))
    color_vec <- c(color_vec,replicate(TM6[2]-TM6[1],'#0070C0'))
    color_vec <- c(color_vec,replicate(TM7[1]-TM6[2],'cornsilk'))
    color_vec <- c(color_vec,replicate(TM7[2]-TM7[1],'#7030A0'))
    color_vec <- c(color_vec,replicate(H8[2]-H8[1],'#C00000'))   
    return(color_vec)
  }
  
  colorizer_partner <- function(partner_dict, partner_name, interaction_vec, partner_color_dict){
    color_vec <- c()
    for (entry in interaction_vec){
      check_residue <- FALSE
      for (substructure in names(partner_dict)){
        for (residue in partner_dict[[substructure]]){
          if (as.numeric(residue) == entry){
            color_vec <- c(color_vec, partner_color_dict[[substructure]])
            check_residue <- TRUE
          }
        }
      }
      if (check_residue == FALSE){color_vec <- c(color_vec, 'cornsilk')}
    }
    
    return(color_vec)
  }
  
  ###---Make below cutoff values white/non visible, first color all links, then subtract the ones out of the threshold---#
  col_df = rand_color(length(df), transparency = 0.5)
  
  col_df[df == 0] = "#00000000"
  dim(col_df) = dim(df)
  
  ###---GPCRs section minimum distance count---###
  register <- c()
  count <- 0
  for (row in row.names(df)){
    count = count + 1
    tick_tack = 0
    for (instance in df[count,]){
      if(instance == 1){
        tick_tack = tick_tack + 1
        break}
    }
    if (tick_tack == 0){
      register <- c(register,count)
    }
  }
  
  ###---Other protein section minimum distance count---###
  register_2 <- c()
  count_2 <- 0
  for (col in colnames(df)){
    count_2 = count_2 + 1
    tick_tack_2 = 0
    for (instance in df[,count_2]){
      if(instance == 1){
        tick_tack_2 = tick_tack_2 + 1
        
        break}
    }
    if (tick_tack_2 == 1){
      register_2 <- c(register_2,count_2)
    }
  }
  ###---Colorize using previous counts---###
  original_colors <- colorizer(TM1,ICL1,TM2,TM3,ICL2,TM4,TM5,ICL3,TM6,TM7,H8)
  original_colors <- original_colors[-register]
  arrestin_color_dict <- c("Triple element" = "#ffff00", "Polar core" = "#ffc000", "Finger loop" = "#70ad47",
                           "Middle loop" = "#00b0f0", "C-loop"="#0070c0", "Lariat loop" = "#ff0000", "Not defined domain" = "cornsilk")
  gprotein_color_dict <- c("HN"="#aee39a","hns1"="#2d595a","S1"="#24ffcd","s1h1"="#589d90",
                           "H1"="#dee0ff","h1ha"="#666699","HA"="#da73f8","hahb"="#ad2270",
                           "HB"="#8575db","hbhc"="#21247b","HC"="#72e5ef","hchd"="#d65f41",
                           "HD"="#852405","hdhe"="#c4937b","HE"="#b1e63b","hehf"="#6c9f30",
                           "HF"="#65f112","hfs2"="#c65911","S2"="#4795e0","s2s3"="#542cd6",
                           "S3"="#e8c346","s3h2"="#eb1241","H2"="#894167","h2s4"="#e313ee",
                           "S4"="#f38ab6","s4h3"="#315d05","H3"="#775a18","h3s5"="#c00000",
                           "S5"="#ff0000","s5hg"="#ffc000","HG"="#ffff00","hgh4"="#7030a0",
                           "H4"="#0070c0","h4s6"="#00b0f0","S6"="#00b050","s6h5"="#b7b7b7",
                           "H5"="#757171")
  if (partner_name %in% ARR_list){usable_partner_colors <- arrestin_color_dict}
  else {usable_partner_colors <- gprotein_color_dict}
  non_GPCR_sector <- colorizer_partner(partner_info, partner_name, register_2, usable_partner_colors)
  ###---Final color vector is built with color for each existing section---###
  colors <- c(original_colors,non_GPCR_sector)
  ###----------------------###
  ###---Plot the Diagram---###
  ###----------------------###
  
  ###---Open the file for saving---###
  file_name <- paste(target_dir,file_name,'.svg')
  new_file_name <- gsub(".csv","",file_name) 
  svg(new_file_name, width=11, height = 9)
  ###---Plot the Diagram using the color vector for sector filling---#
  image <- chordDiagram(df,col = col_df ,annotationTrack = "grid", preAllocateTracks = 1,
                        grid.col = colors)
  #, big.gap = 10, small.gap = 0.2)
  ###---Residue labels---###
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    circos.axis(h = "top", labels.cex = 0.1, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  
  ###---Color labels---###
  legend("bottomleft", cex = 0.8, pch = 0.5, legend = c("TM1","TM2","TM3","TM4","TM5","TM6","TM7","ICL1","ICL2","ICL3","H8","Other"),
         fill=c('#FFFF00','#FFC000','#92D050','#00B050','#00B0F0','#0070C0','#7030A0','#C65911','#375623','#333F4F','#C00000',"cornsilk"))
  #if (partner_name %in% ARR_list){legend("topright", cex = 0.8, pch = 0.5, legend = c(names(usable_partner_colors), "Not defined domain"),
  #       fill=c(usable_partner_colors, 'cornsilk'))}
  legend("topright", cex = 0.8, pch = 0.5, legend = names(usable_partner_colors),
         fill=usable_partner_colors)
  dev.off()}

cutoff_table_builder <- function(distance,file_name){
  
  ###---Loading the table---###
  df = read.table(file_name, sep=';',header = TRUE, stringsAsFactors = FALSE, row.names=1)
  df <-data.matrix(df)
  
  #---Turn table values into binaries depending on distance cutoff value---###
  df <- ifelse(df<distance,1,df)
  df <- ifelse(df>distance,0,df)
  register <- c()
  new_file_name <- gsub("'","",file_name)
  new_file_name <- gsub(".csv","",new_file_name)
  register <- c(register, new_file_name)
  count <- 0
  for (row in row.names(df)){
    count = count + 1
    tick_tack = 0
    for (instance in df[count,]){
      if(instance == 1){
        tick_tack = tick_tack + 1
        position = paste(count,row[1])
        register <- c(register,position) 
        break}    
    }
  }
  return(register)}

correct_factor <- function(input_factor){
  print(input_factor)
  print(as.numeric(as.character(input_factor)))
  ###---Correct the factors into numeric from a pseudo dictionary made of a named list---###
  return(as.numeric(as.character(input_factor)))}

run_graph <- function(input_target, input_csv_file, receptors_list, sub_list, distance, partner_table, partner_name, receptor_entry_split = 4){
  
  ###---Load the weinstein numbering template---###
  template <- read.table(input_csv_file, sep = ";")
  partner_template <- read.table(partner_table, sep = ";")
  current_receptor <- strsplit(strsplit(input_target,split = "-")[[1]][1],split = "_")[[1]][receptor_entry_split]
  pseudo_dictionary_partner <- vector(mode = "list", length = 0)
  ###---Order the partner information vectors---###
  for (partner_row in 1:nrow(partner_template)){
    current_row <- partner_template[partner_row,]
    proper_row <- current_row[!is.na(current_row)]
    if (proper_row[1] == partner_name){
      pseudo_dictionary_partner[proper_row[2]] <- list(proper_row[3:length(proper_row)])}
  }
  ###---Generate an object similar to a Python dictionary---###
  pseudo_dictionary <- vector(mode = "list", length = 11)
  names(pseudo_dictionary) <- sub_list
  struct_count <- 1
  ###---Fill the dictionary with the substructure ranges from the template---###
  for (row in 1:nrow(template)){for(receptor in receptors_list){for(substruct in sub_list){
    if ((struct_count == 2) | (struct_count == 5) | (struct_count == 8)){
      loop_end <- template[row + 1,3]
      pseudo_dictionary[[struct_count]] <- list(end, loop_end)
      struct_count <- struct_count + 1}
    if((template[row,1] == receptor) & (template[row,2] == substruct) & (template[row,1]==current_receptor)){
      start <- template[row,3]
      end <- template[row,4]
      pseudo_dictionary[[struct_count]] <- list(start, end)
      struct_count <- struct_count + 1}
  }}}
  ###---Update length and call the graphic deployer function---###
  corrected_length <- correct_factor(pseudo_dictionary$H8[[2]]) - 1
  
  circle_graph_builder(distance, input_target, corrected_length,
                       c(correct_factor(pseudo_dictionary$TM1[[1]]),correct_factor(pseudo_dictionary$TM1[[2]])),
                       c(correct_factor(pseudo_dictionary$ICL1[[1]]),correct_factor(pseudo_dictionary$ICL1[[2]])),
                       c(correct_factor(pseudo_dictionary$TM2[[1]]),correct_factor(pseudo_dictionary$TM2[[2]])), 
                       c(correct_factor(pseudo_dictionary$TM3[[1]]),correct_factor(pseudo_dictionary$TM3[[2]])),
                       c(correct_factor(pseudo_dictionary$ICL2[[1]]),correct_factor(pseudo_dictionary$ICL2[[2]])),
                       c(correct_factor(pseudo_dictionary$TM4[[1]]),correct_factor(pseudo_dictionary$TM4[[2]])),
                       c(correct_factor(pseudo_dictionary$TM5[[1]]),correct_factor(pseudo_dictionary$TM5[[2]])),
                       c(correct_factor(pseudo_dictionary$ICL3[[1]]),correct_factor(pseudo_dictionary$ICL3[[2]])),
                       c(correct_factor(pseudo_dictionary$TM6[[1]]),correct_factor(pseudo_dictionary$TM6[[2]])),
                       c(correct_factor(pseudo_dictionary$TM7[[1]]),correct_factor(pseudo_dictionary$TM7[[2]])),
                       c(correct_factor(pseudo_dictionary$H8[[1]]),correct_factor(pseudo_dictionary$H8[[2]])),
                       pseudo_dictionary_partner, partner_name)}

###Variable Initialization-----------------------------------###
home = 'C:/Users/marti/OneDrive/Desktop/silverio/processed_results'
setwd(home)
target_dir = 'C:/Users/marti/OneDrive/Desktop/silverio/images/' 
weinstein_template <- "C:/Users/marti/OneDrive/Desktop/silverio/templates/weinstein_numbering_opioids.csv"
weinstein_template_ARR <- "C:/Users/marti/OneDrive/Desktop/silverio/templates/arrestins_final_opioids.csv"
weinstein_template_gprot <- "C:/Users/marti/OneDrive/Desktop/silverio/templates/g_proteins_final_opioids.csv"
DR_list <- c("DOR", "KOR", "MOR", "NOP")
ARR_list <- c("Arr2", "Arr3","Arr2_6U1N","Arr3_6U1N","Arr2_6PWC","Arr3_6PWC")
only_loops <- c("ICL1","ICL2","ICL3")
TM_list <- c("TM1","ICL1","TM2","TM3","ICL2","TM4","TM5","ICL3","TM6","TM7","H8")
distance_value = 8

###Run for all .csv "weinstein" files in the folder--------------###
###Do not forget to erase the output_interacts file--------------###
for(file in list.files(getwd(), pattern = glob2rx('weinstein_inter_chain_*.csv'))){
  interact <- cutoff_table_builder(distance_value,file)
  write.table(matrix(interact,nrow = 1),'output_interacts.csv', append = TRUE, sep =';',col.names = FALSE,row.names = FALSE)
  partner_name <- strsplit(strsplit(file,"-")[[1]][2],".csv")[[1]][1]
  if(partner_name %in% ARR_list){partner_file <- weinstein_template_ARR
  run_graph(file,weinstein_template, DR_list, TM_list, distance_value, partner_file, partner_name, receptor_entry_split = 6)}
  else{partner_file <- weinstein_template_gprot
  run_graph(file,weinstein_template, DR_list, TM_list, distance_value, partner_file, partner_name, receptor_entry_split = 6)}
}