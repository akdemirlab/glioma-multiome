library(dplyr)
genome_cn_ascat <- read.table("/Users/tjmitchell1/catalyst_final/p12_updated_svs.txt", 
header = TRUE,
sep = "\t", 
stringsAsFactors = FALSE)
nrow(genome_cn_ascat)

genome_cn_ascat$start[genome_cn_ascat$start < 1]<-1


###Filter out sex chromosomes###
genome_cn_ascat <- genome_cn_ascat[genome_cn_ascat$chromosome != 'chrX',]
nrow(genome_cn_ascat)

genome_cn_ascat <- genome_cn_ascat[genome_cn_ascat$chromosome != 'chrY',]



##### filter out cells that are mostly diploid or WGD for plotting ####

genome_cn_ascat_find_filters<-genome_cn_ascat

genome_cn_ascat_find_filters$position <- paste0(genome_cn_ascat_find_filters$chrom,':' ,genome_cn_ascat_find_filters$start, '-' ,genome_cn_ascat_find_filters$end )  # No space between names

colnames(genome_cn_ascat)
total_cn<-genome_cn_ascat_find_filters%>%
     dplyr::select(name,position, total_copy_number) %>%
     tidyr::pivot_wider(names_from = position, values_from = total_copy_number)
total_cn<-as.data.frame(total_cn)
dim(total_cn)

rownames(total_cn)<-total_cn[,1]
total_cn <- total_cn[, -1]

upper_bound<-0.75
lowerbound<-0.25
total_cn$number_of_cn2<-(rowSums(total_cn==2))
total_cn$cn2ratio<-total_cn$number_of_cn2/ncol(total_cn)
cells_close_to_diploid<-row.names(total_cn[total_cn$cn2ratio<upper_bound,])
cells_too_aneuploid<-row.names(total_cn[total_cn$cn2ratio>lowerbound,])
length(cells_close_to_diploid)
length(cells_too_aneuploid)

cells_to_keep<-intersect(cells_close_to_diploid,cells_too_aneuploid)
nrow(total_cn)


### make matrix rougly proportional to segment length ###

library(GenomicRanges)
GRList<-list()
temp<-genome_cn_ascat
for (i in unique(genome_cn_ascat$name)) {
  CN_table <- temp[which(temp$name == i),]
  CN_table<-GenomicRanges::makeGRangesFromDataFrame(CN_table,start.field="start",end.field="end",keep.extra.columns=T)
  GRList[[i]]<-CN_table
}

win<-10000000
combined_matrix <- list()
chr<-c()
row_split <- c()
my_group <- c()
row_chr_map <- list()
row_start_map <- list()
for (i in names(GRList)) {
  data <- as.data.frame(GRList[[i]])
  # Calculate number of bins (floor, not round)
  n.times <- ceiling((as.numeric(data$end) - as.numeric(data$start))/win)
  new_matrix <- data$total_copy_number[rep(seq_len(nrow(data)), n.times)]
  combined_matrix[[i]] <- new_matrix
  # Chromosomes (only once)
  if (i == names(GRList)[1]) {
    chr <- data$seqnames
    chr <- chr[rep(seq_len(length(chr)), n.times)]
  }
}

chr_level = paste0('chr',c(1:22))
chr_update = factor(chr, levels = chr_level)
combined_matrix <- do.call(rbind, combined_matrix)


chromosomes_numbas <- gsub("^chr", "", chr)

chromosomes_plotting_new <- factor(chromosomes_numbas, levels = unique(chromosomes_numbas))


ncol(combined_matrix)
combined_matrix_matrix<-as.matrix(combined_matrix)

ncol(combined_matrix_matrix)
combined_matrix_matrix[combined_matrix_matrix >4] <- 5
combined_matrix_matrix[combined_matrix_matrix < 0] <- 0

library(ComplexHeatmap)



col <- c("0" = "#290AD8FF", 
"1" = "#3FA0FFFF", 
"2" = "white", 
"3" = "#F1E201FF",
"4" = "#FF9933FF",
"5" = "#A50021FF")



combined_matrix_matrix<-subset(combined_matrix_matrix, `%in%`(rownames(combined_matrix_matrix), cells_to_keep))

extracted_TIMEPOINT <- gsub(".*SC521_(.*)_W.*", "\\1", rownames(combined_matrix_matrix))
combined_matrix_matrix2<-combined_matrix_matrix
combined_matrix_matrix2$timepoint<-extracted_TIMEPOINT

extracted_TIMEPOINT[extracted_TIMEPOINT=='4A']<-"4"
extracted_TIMEPOINT[extracted_TIMEPOINT=='4D']<-"4"
extracted_TIMEPOINT[extracted_TIMEPOINT=='5A']<-"5"
extracted_TIMEPOINT[extracted_TIMEPOINT=='5D']<-"5"




row_anno3 <- rowAnnotation(
  Timepoint = extracted_TIMEPOINT, row_title = NULL, 
  col = list(
    Timepoint = c("3" = "#FED789FF", "4" = "#72874EFF",  "5"="#476F84FF")),
  annotation_legend_param = list(
     Timepoint = list(
      title = "Timepoint",            
      labels = c("P-12p", "P-12r1", "P-12r2"), 
      at = c("3", "4", "5")           
    )))



pdf('/Users/tjmitchell1/catalyst_final/p12/p12_updated_SVs_050125_main_fig.pdf', width=12, height =8)
Heatmap(combined_matrix_matrix,
cluster_rows = TRUE,
cluster_columns = FALSE, 
show_row_names = FALSE, 
show_column_names = FALSE, 
name = " Total \n\ Copy Number",
column_split = chromosomes_plotting_new, 
col = col,
column_title_gp = gpar(fontsize = 10), 
column_gap = unit(0.008, "npc"),
border=TRUE,
width = unit(20, "cm"),  # Adjust the width
height = unit(14, "cm"), 
left_annotation = row_anno3,
row_split = extracted_TIMEPOINT,  # Group rows based on cluster IDs
use_raster=FALSE, heatmap_legend_param = list(at =c("0","1", "2", "3", "4","5"), labels = c("0","1", "2", "3", "4","> 4"))) 
dev.off()



