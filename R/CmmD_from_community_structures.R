CmmD_from_community_structures <- function(nodelist= NULL, community_structures, resolution_start, resolution_end, interval, distmethod = "hamming", threads){
  start_time <- Sys.time()
  #libraries needed:
  require("AnnotationDbi")
  require("igraph")
  require("stringr")
  require("parallelDist")
  
  if(class(resolution_end)!= 'numeric'){
    stop("ERROR: Resolution parameter must be a number")
  }
  if(class(resolution_start)!= 'numeric'){
    stop("ERROR: Resolution parameter must be a number")
  }
  if(class(interval)!= 'numeric'){
    stop("ERROR: Interval value must be a number")
  }
  
  message(paste0("Resolution parameter starts at: ",resolution_start))
  message(paste0("Resolution parameter ends at: ",resolution_end))
  # message(paste0("Looking for communties found in more than: ",needed_analysis," analysis"))
  
  resolution_interval <- seq(from = resolution_start, to = resolution_end, by = interval)
  
  message(paste0("Reading MolTi output files. Calculating Gene/Community matrix"))
  
  output_files <- community_structures
  alllists <- list()
  for(i in 1:length(output_files)){
    red <- readLines(output_files[i])
    cluster_ids <- grep("Cluster",red)
    lista <- list()
    for(j in 1:length(cluster_ids)){
      st <- cluster_ids[j]
      if(j == length(cluster_ids)){
        en <- length(red)
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-length(current_cluster)]
      }
      else{
        en <- cluster_ids[j+1]
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-c(length(current_cluster),length(current_cluster)-1)]
      }
      lista[[j]] <- current_cluster2[2:length(current_cluster2)]
      names(lista)[j] <- paste0("Cluster_",j)
    }
    kaz <- output_files[i]
    assign(paste0("com_",kaz),value=lista)
    alllists[[i]] <- lista
  }
  names(alllists) <- output_files
  tamano_alllists <- length(alllists)
  allgenes <- unique(unlist(alllists))
  
  if(length(nodelist)>0){
    inter_nodes <- intersect(allgenes,nodelist)
    allgenes <- inter_nodes
  }
  
  message(paste0("Files red. Calculating Gene/Community matrix"))
  res_matrix <- matrix(ncol= tamano_alllists+1, nrow= length(allgenes))
  rownames(res_matrix) <- allgenes
  colnames(res_matrix) <- c(output_files,"Pattern")
  ######Generating the Gene/Community matrix.
  for(i in 1:length(allgenes)){
    gen <- rownames(res_matrix)[i]
    for(j in 1:tamano_alllists){
      searched <- unlist(lapply(alllists[[j]], function(x) gen %in% x))
      comunidad <- unname(which(searched == TRUE))
      res_matrix[i,j] <- comunidad
    }
    res_matrix[i,"Pattern"] <- paste0(res_matrix[i,1:(ncol(res_matrix)-1)],collapse="_")
    percentage <- round((i/length(allgenes)),digits = 4)*100
    porcentajes <- seq(from= 0, to= 100, by = 5)
    progresos <- paste0("Progress: ",porcentajes,"%")
    porcentajes[1] <- 1
    if((percentage %in% porcentajes)==TRUE){
      cual_percentage <- which(porcentajes==percentage)
      to_post <- paste0("Progress: ",percentage,"%")
      message(to_post)
    }
  }
  
  message(paste0("Gene/Community matrix calculated, calculating Hamming distances for all gene pairs. This process may take a while: It takes about 14 min with an Intel Xeon E-2124 processor"))
  
  genes_same_communities <- split(rownames(res_matrix),res_matrix[,"Pattern"])
  final_res_matrix_length <- ncol(res_matrix) - 1
  distance_matrix <- parDist(res_matrix[,1:final_res_matrix_length],method= distmethod ,threads= threads,diag= T)
  final_output <- list(res_matrix[,1:final_res_matrix_length], genes_same_communities,distance_matrix)
  names(final_output) <- c("gene_community_matrix","l_constant","distance_matrix")
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  message(paste0("Run Time: ",diff_time))
  return(final_output)
}
