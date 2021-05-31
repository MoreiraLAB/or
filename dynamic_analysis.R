###----------------------------------------------------------###
###Script to generate scatteplots for interhelical distances-###
###----------------------------------------------------------###

###__author__ = "J.G. Almeida & C.A. Barreto"----------------###
###__email__ = "jose.​gcp.​almeida@gmail.​com"------------------###
###__group__ = "Data-Driven Molecular Design"----------------###
###__group_leader__ = "Irina S. Moreira"---------------------###
###__project__ = "GPCRs"-------------------------------------###
# Importing the necessary libraries --------------------------------------
library(bio3d)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(cowplot)
library(svglite)

# Defining a few functions/constants --------------------------------------
RUN_MODELS <- T

#tm_idx_list object should be a list containing  a vector "start" 
#containing the start of each TM and a vector "end" 
#containing the end of each TM 
tm_idx_list <- data.frame(
  "DOR_start" = c(6,49,84,128,173,216,260,290),
  "DOR_end" = c(44,78,119,154,210,254,288,302),
  "MOR_start" = c(6,42,77,121,166,209,252,282),
  "MOR_end" = c(37,71,112,147,203,247,280,293),
  "KOR_start" = c(6,41,76,120,167,211,254,283),
  "KOR_end" = c(36,70,111,146,205,249,282,295),
  "NOP_start" = c(6,42,77,121,166,209,252,282),
  "NOP_end" = c(37,71,112,147,203,247,280,293)
)

motifs <- c("TM1","TM2","TM3","TM4","TM5","TM6","TM7","H8")

relevant_complexes <- list(
  DOR = c("Gi1", "Gi2", "Gi3", "Go", "Gob", "Gz",
          "Gq_6DDF", "Gq_6OIJ", "G11_6DDF", "G11_6OIJ", "G14_6DDF", "G14_6OIJ", "G15_6DDF", "G15_6OIJ",
          "ARR2_6PWC", "ARR3_6PWC", "ARR2_6U1N", "ARR3_6U1N"),
  MOR = c("Gi1", "Gi2", "Gi3", "Go", "Gob", "Gz",
          "Gslo", "Gssh",
          "Gq_6DDF", "Gq_6OIJ", "G11_6DDF", "G11_6OIJ", "G15_6DDF", "G15_6OIJ",
          "ARR2_6PWC", "ARR3_6PWC", "ARR2_6U1N", "ARR3_6U1N"),
  KOR = c("Gi1", "Gi2", "Gi3", "Go", "Gob", "Gz",
          "Gslo", "Gssh",
          "G12_6DDF", "G12_6OIJ",
          "ARR2_6PWC", "ARR3_6PWC", "ARR2_6U1N", "ARR3_6U1N"),
  NOP = c("Gi1", "Gi2", "Gi3", "Go", "Gob", "Gz",
          "Gq_6DDF", "Gq_6OIJ", "G11_6DDF", "G11_6OIJ", "G14_6DDF", "G14_6OIJ",
          "G12_6DDF", "G12_6OIJ",
          "ARR2_6PWC", "ARR3_6PWC", "ARR2_6U1N", "ARR3_6U1N")
)

dxr_name_list <- c("DOR","MOR","KOR","NOP")

gprot_name_list <- c("Gi1","Gi2","Gi3","Go","Gob","Gz",
                     "Gslo","Gssh",
                     "Gq_6DDF","Gq_6OIJ","G11_6DDF","G11_6OIJ","G14_6DDF","G14_6OIJ","G15_6DDF","G15_6OIJ",
                     "G12_6DDF","G12_6OIJ")

arrestin_name_list <- c("ARR2_6PWC","ARR3_6PWC","ARR2_6U1N","ARR3_6U1N")

ORDER_FINAL <- c(
  arrestin_name_list,
  gprot_name_list
)

#Function to calculate BC values per TM
get_tm_bc <- function(complex_list,template_string,dxr_string,partner_string){
  dxr_monomer_name <- grep(pattern = template_string,names(dxr_monomers),value = T) %>%
    grep(pattern = dxr_string,value = T)
  dxr_complex_name <- grep(pattern = sprintf('%s\\.',partner_string),
                           names(complex_list),value = T) %>%
    grep(pattern = gsub("-",".*",dxr_string),value = T)
  
  modes_monomer <- nma(dxr_monomers[[dxr_monomer_name]])
  modes_complex <- nma(complex_list[[dxr_complex_name]])
  
  cov_monomer <- cov.nma(modes_monomer)
  cov_complex <- cov.nma(modes_complex)
  
  tm_dxr <- colnames(tm_idx_list) %>% 
    grep(pattern = str_split(dxr_string,pattern = '-')[[1]][1])
  bc_list <- list()
  for (i in 1:nrow(tm_idx_list)) {
    a <- tm_idx_list[i,tm_dxr[1]]
    b <- tm_idx_list[i,tm_dxr[2]]
    bc_list[[i]] <- bhattacharyya(a = cov_monomer[a:b,a:b],cov_complex[a:b,a:b])  
  }
  
  output_list <- list(
    nma_complex = modes_complex, # modes for the complex
    nma_monomer = modes_monomer, # modes for the monomer
    cov_monomer = cov_monomer, # nma covariance for the monomer
    cov_complex = cov_complex, # nma covariance for the complex
    bc_list = unlist(bc_list) # BC for each TM
  )
  return(output_list)
}

filter_relevant_group <- function(df) {
  rel_dxr <- df$DxR[1]
  (df$Partner %in% relevant_complexes[[rel_dxr]]) %>%
    return
}

# Loading all PDB files ---------------------------------------------------
all_files <- list.files(FINAL_MODEL_PATH,pattern = "*pdb",full.names = T)

dxr_monomers <- list()
g_prot_complexes <- list()
arrestin_complexes <- list()
for (pdb_file in all_files) {
  pdb_file_base <- basename(pdb_file)
  if (grepl("[A-Z]O[P-R]",pdb_file_base) == TRUE) {
    if (grepl("G",pdb_file_base)) {
      g_prot_complexes[[pdb_file_base]] <- read.pdb(pdb_file)
    } else if (grepl("ARR",pdb_file_base)) {
      arrestin_complexes[[pdb_file_base]] <- read.pdb(pdb_file)
    } else dxr_monomers[[pdb_file_base]] <- read.pdb(pdb_file)
  }
}

# Computing BC for arrestins ----------------------------------------------

if (RUN_MODELS == TRUE) {
  arrestin_bc_list <- list()
  for (arrestin_string in arrestin_name_list) {
    for (dxr_string in dxr_name_list) {
      if (grepl("6U1N",arrestin_string)) {
        template_string <- "6U1N"
      } else {
        template_string <- "6PWC"
      }
      
      out_name <- paste(dxr_string,arrestin_string,sep='-')
      arrestin_bc_list[[out_name]] <- get_tm_bc(complex_list=arrestin_complexes,
                                                template_string=template_string,
                                                dxr_string=dxr_string,
                                                partner_string=arrestin_string)
    }
  }
  arrestin_bc_scores <- lapply(arrestin_bc_list, function(x) x$bc_list) %>% 
    do.call(what = rbind)
  
  arrestin_fluc_fold_change <- lapply(names(arrestin_bc_list), function(name) {
    x <- arrestin_bc_list[[name]]
    name_split <- strsplit(name,'-') %>% unlist()
    dxr_name <- name_split[[1]]
    se <- grep(dxr_name,colnames(tm_idx_list))
    start <- tm_idx_list[,se[1]]
    end <- tm_idx_list[,se[2]]
    fluc_ff <- c()
    for (i in 1:8) {
      sub_start <- start[i]
      sub_end <- end[i]
      cmp <- x$nma_complex$fluctuations[sub_start:sub_end]
      mon <- x$nma_monomer$fluctuations[sub_start:sub_end]
      fluc_ff <- c(fluc_ff,(cmp / mon) %>% mean)
    }
    return(fluc_ff)
  }) %>%
    do.call(what = rbind) %>%
    as.data.frame()
  colnames(arrestin_fluc_fold_change) <- paste0("Fluc_",motifs)
  arrestin_fluc_fold_change$Complex <- names(arrestin_bc_list)
  arrestin_fluc_fold_change$DXR <- sapply(arrestin_fluc_fold_change$Complex,
                                          function(x) unlist(strsplit(x,'-'))[[1]])
  arrestin_fluc_fold_change$Partner <- sapply(arrestin_fluc_fold_change$Complex,
                                              function(x) unlist(strsplit(x,'-'))[[2]])
  
  save(arrestin_bc_scores,arrestin_fluc_fold_change,file = "arrestin-dynamics.RData")
} else {
  load("arrestin-dynamics.RData")
}

# Computing BC for G-proteins ---------------------------------------------

if (RUN_MODELS == TRUE) {
  gprot_bc_list <- list()
  for (dxr_string in dxr_name_list) {
    template_string <- names(dxr_monomers) %>%
      grep(pattern = dxr_string,value = T) %>%
      grep(pattern = "6DDF",value = T) %>%
      gsub(pattern = ".pdb",replacement = "") %>%
      gsub(pattern = dxr_string,replacement = "") %>%
      gsub(pattern = '_',replacement = "")
    sub_gprot <- names(g_prot_complexes) %>% 
      grep(pattern = dxr_string,value = T) %>% 
      gsub(pattern = ".pdb",replacement = "") %>%
      gsub(pattern = dxr_string,replacement = "") %>%
      gsub(pattern = '-',replacement = "")
      
    for (gprot_name in sub_gprot) {
      out_name <- paste(dxr_string,gprot_name,sep='-')
      gprot_bc_list[[out_name]] <- get_tm_bc(g_prot_complexes,template_string,dxr_string,gprot_name)
    }
  }
  
  gprot_bc_scores <- lapply(gprot_bc_list, function(x) x$bc_list) %>% 
    do.call(what = rbind)
  gprot_fluc_fold_change <- lapply(names(gprot_bc_list), function(name) {
    x <- gprot_bc_list[[name]]
    name_split <- strsplit(name,'-') %>% unlist()
    dxr_name <- name_split[[1]]
    se <- grep(dxr_name,colnames(tm_idx_list))
    start <- tm_idx_list[,se[1]]
    end <- tm_idx_list[,se[2]]
    fluc_ff <- c()
    for (i in 1:8) {
      sub_start <- start[i]
      sub_end <- end[i]
      cmp <- x$nma_complex$fluctuations[sub_start:sub_end]
      mon <- x$nma_monomer$fluctuations[sub_start:sub_end]
      fluc_ff <- c(fluc_ff,(cmp / mon) %>% mean)
    }
    return(fluc_ff)
  }) %>%
    do.call(what = rbind) %>%
    as.data.frame()
  colnames(gprot_fluc_fold_change) <- paste0("Fluc_",motifs)
  gprot_fluc_fold_change$Complex <- names(gprot_bc_list)
  gprot_fluc_fold_change$DXR <- sapply(gprot_fluc_fold_change$Complex,
                                       function(x) unlist(strsplit(x,'-'))[[1]])
  gprot_fluc_fold_change$Partner <- sapply(gprot_fluc_fold_change$Complex,
                                           function(x) unlist(strsplit(x,'-'))[[2]])
  
  save(gprot_bc_scores,gprot_fluc_fold_change, file = "gprot-dynamics.RData")
} else {
  load("gprot-dynamics.RData")
}

# Combining and analysing everything --------------------------------------
THEME_MINIMAL_BASE_SIZE <- 5.7

all_bc_scores <- rbind(arrestin_bc_scores,gprot_bc_scores) %>% 
  data.frame
colnames(all_bc_scores) <- c("TM1","TM2","TM3","TM4","TM5","TM6","TM7","HX8")
dxr_partner <- rownames(all_bc_scores) %>%
  sapply(function(x) unlist(strsplit(x,'-')),simplify = F) %>%
  lapply(function(x) {
    x <- c(x[1],x[length(x)])
  }) %>% 
  do.call(what=rbind) 
colnames(dxr_partner) <- c("Receptor","Partner")
all_bc_scores <- cbind(all_bc_scores,dxr_partner) 
all_bc_scores <- all_bc_scores[all_bc_scores %>%
  apply(1, function(x) {
    x["Partner"] %in% unlist(relevant_complexes[[x["Receptor"]]]) %>% 
      return
  }),]

all_fluc_fold_changes <- rbind(arrestin_fluc_fold_change,gprot_fluc_fold_change)
colnames(all_fluc_fold_changes) <- c("Fluc_TM1","Fluc_TM2","Fluc_TM3","Fluc_TM4",
                                     "Fluc_TM5","Fluc_TM6","Fluc_TM7","Fluc_H8",
                                     "Complex","Receptor","Partner")

all_bc_scores_long <- all_bc_scores %>%
  gather(key = "Motif",value = "Value",TM1,TM2,TM3,TM4,TM5,TM6,TM7,HX8) %>% 
  mutate(Motif = ifelse(Motif == 'HX8','H8',Motif)) %>% 
  transmute(Receptor = factor(Receptor),
            Partner = factor(Partner,
                             levels = c(rev(gprot_name_list),
                                        rev(arrestin_name_list))),
            Motif = factor(Motif,
                           levels = motifs),
            Value = Value) 
  
flexibility_heatmap <- all_bc_scores_long %>%
  ggplot(aes(x = Receptor,y = Partner,fill = Value)) + 
  geom_tile() + 
  geom_label(aes(label = round(Value,2)),
             size = 1.2,
             fill = "white",
             label.r = unit(0,"lines"),
             label.size = 0,
             label.padding = unit(0.3,"mm"),
             alpha = 0.7) + 
  facet_wrap(~ Motif,ncol = 4) + 
  theme_minimal(base_size = THEME_MINIMAL_BASE_SIZE) + 
  theme(axis.ticks = element_blank(),panel.grid = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.1,"cm")) +
  scale_fill_viridis_c(name = NULL) +
  xlab("") +
  ylab("") +
  ggtitle("Flexibility change")

flexibility_heatmap_no_label <- all_bc_scores_long %>%
  ggplot(aes(x = Receptor,y = Partner,fill = Value)) + 
  geom_tile() + 
  facet_wrap(~ Motif,ncol = 4) + 
  theme_minimal(base_size = THEME_MINIMAL_BASE_SIZE) + 
  theme(axis.ticks = element_blank(),panel.grid = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.1,"cm")) +
  scale_fill_viridis_c(name = NULL) +
  xlab("") +
  ylab("") +
  ggtitle("Flexibility change")

flexibility_heatmap_no_label + 
  ggsave(filename = "heatmap-dynamics-no-label.svg",width = 5,height = 3.5,dpi=600) + 
  ggsave(filename = "heatmap-dynamics-no-label.pdf",width = 5,height = 3.5,dpi=600) +
  ggsave(filename = "heatmap-dynamics-no-label.png",width = 5,height = 3.5,dpi=600)


pc <- cmdscale(dist(all_bc_scores[,c(2,3,4,5,6,7)]))
dim_red <- cbind(
  data.frame(pc),
  data.frame(all_bc_scores$Receptor,all_bc_scores$Partner))
colnames(dim_red) <- c("x","y","Receptor","Partner")
dim_red <- dim_red %>% 
  merge(y = all_fluc_fold_changes,by = c("Receptor","Partner"))
dim_red$Average_fold_change <- rowMeans(dim_red[,grep("Fluc",colnames(dim_red))])

dim_red_scatter <- ggplot(data = dim_red,aes(x = x,y = y,colour = Receptor)) +
  geom_point(alpha = 0.8,aes(size = Average_fold_change * 0.3)) + 
  geom_text_repel(color = 'black',size = 1.5,aes(label = paste(Receptor,Partner,sep = '-')),
                  min.segment.length = 0.1,
                  segment.alpha = 0.4,
                  segment.size = 0.3) + 
  xlab("PrCoord1") + 
  ylab("PrCoord2") +
  scale_color_manual(breaks = c("DOR","MOR","KOR","NOP"),
                     values = c("#1A1AFF","#FF9933","#33FF33","#CC6666"),
                                name = NULL) +
  theme_minimal(base_size = THEME_MINIMAL_BASE_SIZE) +
  scale_size_continuous(guide = F,range = c(0.5,2))

dim_red_scatter %>% ggsave(filename = "pca-dynamics.svg",width = 7,height = 5.5,
                           dpi = 600) 
dim_red_scatter %>% ggsave(filename = "pca-dynamics.pdf",width = 7,height = 5.5,
                           dpi = 600) 

all_fluc_fold_changes_long <- all_fluc_fold_changes %>%
  gather(value = "Value",key = "Motif",Fluc_TM1,Fluc_TM2,Fluc_TM3,
         Fluc_TM4,Fluc_TM5,Fluc_TM6,Fluc_TM7,Fluc_H8) %>% 
  transmute(Receptor = factor(Receptor),
            Partner = factor(Partner,
                             levels = c(rev(gprot_name_list),
                                        rev(arrestin_name_list))),
            Motif = gsub("Fluc_","",Motif) %>% factor(levels = motifs),
            Value = Value) 

fluc_heatmap <- ggplot(data = all_fluc_fold_changes_long,aes(x = Receptor,y = Partner,fill = Value)) +
  geom_tile() + 
  geom_label(aes(label = round(Value,2)),
             size = 1.2,
             fill = "white",
             label.r = unit(0,"lines"),
             label.size = 0,
             label.padding = unit(0.3,"mm"),
             alpha = 0.7) + 
  facet_wrap(~ Motif,ncol = 4) + 
  theme_minimal(base_size = THEME_MINIMAL_BASE_SIZE) + 
  theme(axis.ticks = element_blank(),panel.grid = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.1,"cm")) +
  scale_fill_viridis_c(name = NULL) +
  xlab("") +
  ylab("") +
  ggtitle("Average fluctuation fold change") 

fluc_heatmap_no_label <- ggplot(data = all_fluc_fold_changes_long,aes(x = Receptor,y = Partner,fill = Value)) +
  geom_tile() + 
  facet_wrap(~ Motif,ncol = 4) + 
  theme_minimal(base_size = THEME_MINIMAL_BASE_SIZE) + 
  theme(axis.ticks = element_blank(),panel.grid = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.1,"cm")) +
  scale_fill_viridis_c(name = NULL) +
  xlab("") +
  ylab("") +
  ggtitle("Average fluctuation fold change")

fluc_heatmap_no_label + 
  ggsave(filename = "heatmap-fluc-no-label.pdf",width = 5,height = 3.5,dpi=600) + 
  ggsave(filename = "heatmap-fluc-no-label.svg",width = 5,height = 3.5,dpi=600) +
  ggsave(filename = "heatmap-fluc-no-label.png",width = 5,height = 3.5,dpi=600)
fluc_heatmap +
  ggsave(filename = "heatmap-fluc.svg",width = 5,height = 3.5,dpi=600) +
  ggsave(filename = "heatmap-fluc.pdf",width = 5,height = 3.5,dpi=600) +
  ggsave(filename = "heatmap-fluc.png",width = 5,height = 3.5,dpi=600)


final_plot <- plot_grid(
  plot_grid(fluc_heatmap,
            flexibility_heatmap,
            nrow = 1,
            labels = c("A","B"),
            label_size = 10,
            hjust = 0.0),
  plot_grid(ggplot() + geom_blank() + theme_minimal(),
            dim_red_scatter,
            ggplot() + geom_blank() + theme_minimal(),
            rel_widths = c(0.2,0.6,0.2),
            labels = c("","C",""),
            label_size = 10,
            hjust = 0.1,
            nrow = 1),
  labels = NULL,ncol = 1,
  rel_heights = c(1,1))
final_plot %>% ggsave(filename = "final-dynamics.svg",height = 6.5,width = 8,dpi = 600)
final_plot %>% ggsave(filename = "final-dynamics.pdf",height = 6.5,width = 8,dpi = 600)
final_plot %>% ggsave(filename = "final-dynamics.png",height = 6.5,width = 8,dpi = 600)

final_plot_no_labels <- plot_grid(
  plot_grid(fluc_heatmap_no_label,
            flexibility_heatmap_no_label,
            nrow = 1,
            labels = c("A","B"),
            label_size = 10,
            hjust = 0.0),
  plot_grid(ggplot() + geom_blank() + theme_minimal(),
            dim_red_scatter,
            ggplot() + geom_blank() + theme_minimal(),
            rel_widths = c(0.2,0.6,0.2),
            labels = c("","C",""),
            label_size = 10,
            hjust = 0.1,
            nrow = 1),
  labels = NULL,ncol = 1,
  rel_heights = c(1,1))
final_plot_no_labels %>% ggsave(filename = "images/final-dynamics-no-label.svg",height = 6.5,width = 8,dpi = 600)
final_plot_no_labels %>% ggsave(filename = "images/final-dynamics-no-label.pdf",height = 6.5,width = 8,dpi = 600)
final_plot_no_labels %>% ggsave(filename = "images/final-dynamics-no-label.png",height = 6.5,width = 8,dpi = 600)
