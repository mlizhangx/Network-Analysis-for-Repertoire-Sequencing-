# test on mac

library(igraph)
library(reshape)
library(Matrix)
library(dplyr)
library(stringr)

setwd("./output/")
set.seed(9999)


#### match raw data with distance matrix ####
  
  sample_id <- "TRB-Pt-10-2"
  print(sample_id)
  dir.create(paste0("./", sample_id,"_hamming"), showWarnings = FALSE)
  input_data <- read.table(file = paste0("../data/TRB-Pt-10-2-500ng-15-04-2020-gDNA_S48.clones.txt"), sep = '\t', header = TRUE)
  
  dim(input_data) 
  length(unique(input_data$aaSeqCDR3)) 
  
  sub_input <- input_data[,c("nSeqCDR3", "aaSeqCDR3", "cloneCount",	"cloneFraction")]
  
  agg_count_on_aa <- aggregate(sub_input[,c( "cloneCount",	"cloneFraction")],
            by = list(sub_input$aaSeqCDR3),
            FUN = sum)
  colnames(agg_count_on_aa) <- c("aaSeqCDR3", "cloneCount", "cloneFraction")
  
  
  unique_aa_data <- as.data.frame(table(sub_input$aaSeqCDR3))
  colnames(unique_aa_data)<- c("aaSeqCDR3","unique_nucleotide_count")
  dim(unique_aa_data)

  sub_input_unique_count <- merge(agg_count_on_aa, unique_aa_data, by = "aaSeqCDR3")
  dim(sub_input_unique_count)
  
  input_raw_data <- sub_input_unique_count
  dim(input_raw_data)
  input_raw_data_unique <- input_raw_data[!duplicated(input_raw_data$aaSeqCDR3), ] # count if not AA count but one nucleotide count
  dim(input_raw_data_unique) # 3776  
  input_raw_data <- input_raw_data_unique
  
  start_time <- Sys.time()
  # need remove the seq withtou AA info, and AA length < 3
  input_raw_data$AA_length <- nchar(as.character(input_raw_data$aaSeqCDR3))
  sub_input_raw <- input_raw_data[ input_raw_data$aaSeqCDR3 != ""  & input_raw_data$AA_length > 2 ,]
  input_AA <- as.character(sub_input_raw$aaSeqCDR3)
  print(length(input_AA))
  if(length(input_AA) == 1) {next}
  
  
  write.table(input_AA, paste0("./",  sample_id, "_hamming/","string_info.txt"), quote = FALSE,  col.names=FALSE,  row.names=FALSE, sep="\t")
  setwd(paste0("./", sample_id, "_hamming/"))
  # source python code
  system("python3 ../../script/distance_matrix_function_v9_hamming_distance.py")
  
  # rename the hamming distance file by sample_id and Vgene
  file.rename("hamming_distance_martrix_degree_lt0.mtx", paste0("distance_matrix_betaCDR3_Europe_covid_data_", sample_id ,".mtx"))
  # test_read_mtx <- readMM("distance_matrix_betaCDR3_Europe_covid_data_Pt-1-1.mtx")
  # dim(test_read_mtx)
  select_col_id <- read.csv("select_col_id.csv", header = FALSE)
  file.remove("select_col_id.csv")
  file.remove("string_info.txt")
  index_cluster_gt1 <- as.logical(select_col_id$V1)
  meta_data <- sub_input_raw[index_cluster_gt1,]
  print(dim(meta_data))
  write.csv(meta_data, file =  paste0("./", "raw_data_betaCDR3_Europe_covid_data_", sample_id, ".csv"))
  
  end_time <- Sys.time()
  print(end_time - start_time)
  setwd("../../")
  cat(" \n")

