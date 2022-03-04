library(ggplot2)
library(ggseqlogo)
library(RecordLinkage)
library(igraph)
library(reshape)
library(RColorBrewer)
library(Matrix)
library(sparklyr)
library(dplyr)
library(stringr)
library(ggraph)

setwd("./output/")


  sample_id <- "TRB-Pt-10-2"
  print(sample_id)
  
  filename_hamming <- paste0("./",  sample_id,"_hamming","/distance_matrix_betaCDR3_Europe_covid_data_", sample_id,".mtx")
  
  if(file.exists(filename_hamming) == 0) next
  input_mtx <- readMM(filename_hamming)
  print(dim(input_mtx))
  input_meta <- paste0("./", sample_id,"_hamming", "/raw_data_betaCDR3_Europe_covid_data_", sample_id,".csv")
  meta_lv <- read.csv(input_meta)
  dim(meta_lv)
  
  
  #### network analysis ####
  input_raw_data <- meta_lv
  input_distance_matrix <- input_mtx
  
  dim(input_distance_matrix)
  dim(input_raw_data)

  net <- graph_from_adjacency_matrix(input_distance_matrix)
  net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T))
  
  #### subset network by degree > 0 ####
  deg <- degree(net)
  # print(table(deg))
  # print(sample_id)
  # deg > 0
  dim(input_distance_matrix)
  #   sub_select_matrix <- input_distance_matrix[deg > 0,deg > 0]
  sub_select_matrix <- input_distance_matrix
  # dim(sub_select_matrix)
  
  #   sub_input_raw_data <- input_raw_data_sel[deg > 0,]
  sub_input_raw_data <- input_raw_data
  dim(sub_input_raw_data)
  
  # convert the matrix into igraph format
  net <- graph_from_adjacency_matrix(sub_select_matrix, weighted=TRUE)
  net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T)) 
  
  # add cluster_id in raw data  
  cfg <- cluster_fast_greedy(as.undirected(net))
  sub_input_raw_data$cluster_id <- cfg$membership
  sub_input_raw_data$deg <- deg
  
  dir.create(paste0("./",sample_id,"_hamming"), showWarnings = FALSE)
  pdf(paste0("./", sample_id,"_hamming", "/Europe_covid_data_",  sample_id, ".pdf"))
  ### best pick layout
  set.seed(9999)
  l <- layout_components(net)
  
  # match figure a on expand vs non-expand
  print(ggraph(net,layout = l)+
          geom_edge_link0(width=0.1,colour="grey",alpha =1 )+
          geom_node_point(aes(color = 'red'),size=1) +
          # geom_node_point(aes(color = cell_type,size = cloneCount)) + scale_size(range = c(0.1,log(max(cloneCount))/2.5)) +
          theme_graph(base_family="sans") +
          ggtitle(paste0('Beta chain, Patient: ', sample_id)) )
  
  
  # CloneFraction <- as.numeric(sub_input_raw_data$cloneFraction)
  # print(ggraph(net,layout = l)+
  #         geom_edge_link0(width=0.1,colour="grey",alpha =1 )+
  #         geom_node_point(aes(color = 'red',size=CloneFraction)) +
  #         # geom_node_point(aes(color = cell_type,size = cloneCount)) + scale_size(range = c(0.1,log(max(cloneCount))/2.5)) +
  #         theme_graph(base_family="sans") +
  #         ggtitle(paste0('Beta chain, Patient: ', sample_id)) )
  
  CloneCount <- as.numeric(sub_input_raw_data$cloneCount)
  Degree <- as.character(sub_input_raw_data$deg)
  print(ggraph(net,layout = l)+
          geom_edge_link0(width=0.1,colour="grey",alpha =1 )+
          # geom_node_point(aes(color = 'red', size = CloneCount)) +
          geom_node_point(aes(color = Degree,size = CloneCount)) + scale_size(range = c(0.1,log(max(CloneCount))/2.5)) +
          theme_graph(base_family="sans") +
          ggtitle(paste0('Beta chain, Patient: ', sample_id)) )
  
  # add cluster id in plot
  membership_label <- as.factor(sub_input_raw_data$cluster_id)
  membership_label[duplicated(membership_label)] <- ""
  plot_cluster_id <- ggraph(net,layout = l)+
                             geom_edge_link0(colour="grey")+
                             geom_node_point(aes(color = Degree,size = CloneCount)) + scale_size(range = c(0.1,log(max(CloneCount))/2.5)) +
                             theme_graph(base_family="sans") +
                             # theme(legend.position = "none") +
                             geom_node_text(aes(label=membership_label),size =2.5)+
                             ggtitle(paste0('Beta chain, Patient: ', sample_id)) 
  print(plot_cluster_id)
  
  dev.off()
  
  
  #### add network statistics ####
  # input data: sub_input_raw_data
  print(dim(sub_input_raw_data))
  
  # sub_input_raw_data$betaCDR3_length <- apply(sub_input_raw_data[,2,drop=FALSE],2,nchar)[,1]
  sub_input_raw_data$betaCDR3_length <- apply(sub_input_raw_data[,2,drop=FALSE],2,nchar)[,1]
  
  sub_input_raw_data$deg <- degree(net)
  sub_input_raw_data$transitivity <- transitivity(net, type="local")
  # sub_input_raw_data$closeness <- closeness(net, mode="all", weights=NA)
  # sub_input_raw_data$centr_clo_res <- centr_clo(net, mode="all", normalized=T)$res
  sub_input_raw_data$eigen_centrality <- eigen_centrality(net, directed=T, weights=NA)$vector
  sub_input_raw_data$centr_eigen <- centr_eigen(net, directed=T, normalized=T)$vector
  sub_input_raw_data$betweenness <- betweenness(net, directed=T, weights=NA)
  sub_input_raw_data$centr_betw <- centr_betw(net, directed=T, normalized=T)$res
  sub_input_raw_data$authority_score <- authority_score(net, weights=NA)$vector
  sub_input_raw_data$coreness <- coreness(net, mode="all")
  sub_input_raw_data$page_rank <- page_rank(net)$vector
  
  # cfg <- cluster_fast_greedy(as.undirected(net))
  sub_input_raw_data$membership_stat <- cfg$membership
  
  write.csv(sub_input_raw_data, file = paste0("./", sample_id,"_hamming", "/Europe_covid_data_cell_meta_data_w_stat_", sample_id, ".csv"), row.names=FALSE)
  
  
  
  
  
  #### create membership table for each patient ####
  dim(sub_input_raw_data)
  table(sub_input_raw_data$membership_stat)
  membership_table <- as.data.frame(table(sub_input_raw_data$membership_stat))
  colnames(membership_table) <- c("membership","node_count")
  
  
  
  total_membership <- length(membership_table$membership)
  # membership_table$motif_top_deg_alpha <- ""
  
  membership_table$motif_top_deg_beta <- ""
  membership_table$max_deg_within_cluster <- ""
  membership_table$deg_mean <- 0
  membership_table$motif_top_count_beta <- ""
  membership_table$max_count_within_cluster <- ""
  membership_table$cluster_total_count <- 0
  # membership_table$motif_freq_50_alpha <- ""
  membership_table$motif_freq_50_beta <- ""
  
  membership_table$betaCDR3_length <- 0
  # membership_table$alphaCDR3_length <- 0
  
  
  # membership_table$Count_PRE_INFUSION <- 0
  # membership_table$Count_DOSE_2 <- 0
  # social network properites
  # membership_table$deg_avg <- 0
  membership_table$diam_length <- 0
  membership_table$assortativity <- 0
  membership_table$transitivity <- 0
  membership_table$edge_density <- 0
  membership_table$centr_degree <- 0
  membership_table$centr_clo <- 0
  membership_table$eigen_centrality <- 0
  membership_table$centr_eigen <- 0
  
  
  for(membership_id in 1:total_membership) {
    # print(paste0(membership_id, " out of ", total_membership))
    # membership_id <- 1
    membership_table[membership_table$membership == membership_id,]$deg_mean  <- round(mean(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$deg),2)
    membership_table[membership_table$membership == membership_id,]$betaCDR3_length  <- round(mean(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$betaCDR3_length),2)
    # membership_table[membership_table$membership == membership_id,]$alphaCDR3_length  <- round(mean(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$alphaCDR3_length),2)
    
    membership_table[membership_table$membership == membership_id,]$cluster_total_count  <- sum(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$cloneCount)
    # membership_table[membership_table$membership == membership_id,]$Count_PRE_INFUSION  <- sum(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$Count_PRE_INFUSION)
    # membership_table[membership_table$membership == membership_id,]$Count_DOSE_2  <- sum(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$Count_DOSE_2)
    max_deg <- max(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$deg)
    # membership_table[membership_table$membership == membership_id,]$motif_top_deg_alpha  <- as.character(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$deg_mean == max_deg,]$alphaCDR3[1])
    membership_table[membership_table$membership == membership_id,]$max_deg_within_cluster <- max_deg
    
    membership_table[membership_table$membership == membership_id,]$motif_top_deg_beta  <- as.character(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$deg == max_deg,]$aaSeqCDR3[1])
    
    max_count <- max(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$cloneCount)
    membership_table[membership_table$membership == membership_id,]$max_count_within_cluster <- max_count
    membership_table[membership_table$membership == membership_id,]$motif_top_count_beta  <- as.character(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$cloneCount == max_count,]$aaSeqCDR3[1])
    
    #### get the representative motif ####
    # take long time to run, comment out
    # betaCDR3_length_round <- round(membership_table[membership_table$membership == membership_id,]$betaCDR3_length)
    # for(i in 1:betaCDR3_length_round) {
    #   # print(i)
    #   string_list <- sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$betaCDR3_length == betaCDR3_length_round,]$aminoAcid
    #   freq_table <- as.data.frame(table(substring(string_list, i,i)))
    #   select_letter <- ifelse(dim(freq_table[freq_table$Freq > length(string_list)/2,])[1] == 1, as.character(freq_table[freq_table$Freq > length(string_list)/2,]$Var1), "*")
    #   # print(select_letter)
    #   select_letter_index <- paste0("char_", i)
    #   assign(select_letter_index, select_letter)
    # }
    # membership_table[membership_table$membership == membership_id,]$motif_freq_50_beta  <- paste(noquote(mget(mixedsort(ls(pattern= "char_")))), collapse = "")
    # rm(list = ls(pattern= "char_"))
    
    
    
    ### network properties ###
    input_matrix <- sub_select_matrix
    input_data_membership <- sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]
    input_data_membership_index <- rownames(sub_input_raw_data) %in% rownames(input_data_membership)
    
    input_matrix_membership <- input_matrix[input_data_membership_index,input_data_membership_index]
    # input_matrix_membership[input_matrix_membership > 1] <- 0
    
    matrix_2_net <- as.matrix(input_matrix_membership)
    net <- graph_from_adjacency_matrix(matrix_2_net)
    net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T))
    
    deg <- degree(net, mode="all")
    # table(deg) # there should be no deg== 0.
    # membership_table[membership_table$membership == membership_id,]$deg_avg <- round(mean(deg),2)
    
    # Diameter (longest geodesic distance)
    diam <- get_diameter(net, directed=T)
    # diam
    membership_table[membership_table$membership == membership_id,]$diam_length <- length(diam)
    
    # Assortativity
    membership_table[membership_table$membership == membership_id,]$assortativity <- assortativity_degree(net, directed=F)
    
    # Transitivity
    membership_table[membership_table$membership == membership_id,]$transitivity <- transitivity(net, type="global")  # net is treated as an undirected network
    
    # Density
    # The proportion of present edges from all possible ties.
    membership_table[membership_table$membership == membership_id,]$edge_density <- edge_density(net, loops=F)
    
    # centralization on degree
    membership_table[membership_table$membership == membership_id,]$centr_degree <- centr_degree(net, mode="in", normalized=T)$centralization
    
    # centralization on Closeness (centrality based on distance to others in the graph)
    membership_table[membership_table$membership == membership_id,]$centr_clo <- centr_clo(net, mode="all", normalized=T)$centralization
    
    # centralization on Eigenvector (centrality proportional to the sum of connection centralities)
    # Values of the first eigenvector of the graph adjacency matrix
    membership_table[membership_table$membership == membership_id,]$eigen_centrality <- eigen_centrality(net, directed=T, weights=NA)$value
    membership_table[membership_table$membership == membership_id,]$centr_eigen<- centr_eigen(net, directed=T, normalized=T)$centralization
  }
  # head(membership_table)
  write.csv(membership_table, file = paste0("./", sample_id,"_hamming", "/Europe_covid_data_cluster_level_info_", sample_id, ".csv"), row.names=FALSE)
  # }
# }
