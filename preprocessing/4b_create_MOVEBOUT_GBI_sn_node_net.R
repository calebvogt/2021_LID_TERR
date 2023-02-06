## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

# load  & clean data ---------------------------------------------------------------
dd <- paste(getwd(), "data/", sep = "/")
filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
social_data = lapply(filenames, fread) ## READ IN ALL FILES
meta <- read.csv("data/metadata.csv")

# clean social data for triaged mice from social interaction bouts. 
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME

  df2 <- df %>% 
    select(!(V1)) %>% 
    filter(day %in% (1:20))
  
  print(aa)
  social_data[[aa]] <- df2
}

# Create nets  ------------------------------------------------------------
trial_net_stats_list <- list()
trial_vertex_stats_list <- list()
aa=1
for(aa in 1:length(social_data)) { ## Loop through trials 
  df <- social_data[[aa]] 
  trial <- print(unique(df$trial))
  
  net_stats_list <- list()
  vertex_stats_list <- list()
  bb=1
  for(bb in 1:20) { # loop through days
    focals <- meta %>%  ## select your focals for this run
      # filter(sex == "M") %>% ## choose males
      # filter(sex == "F") %>% ## choose females
      filter(if(bb %in% 1:5) {drop==1} # choose early animals for days 1-5, or adjust to have them for the whole thing. 
             else {drop %in% 1:2}) ## choose all animals for days 6-20, i.e. drop 1 and drop 2 mice
      
    focal_names <- unique(focals$strain_sex_name)
    
    gbi <- df %>%  # select mouse columns
      filter(day == bb) %>% #comment out if you want the full 10 day network and run from here below
      dplyr::select(matches(focal_names)) # choose your columns of mice
    
    colnames(gbi) <- gsub("C57-M-","",colnames(gbi)) # remove strain-sex info from colnames of strain you are working with. 
    colnames(gbi) <- gsub("C57-F-","",colnames(gbi))
    ids <- colnames(gbi) ## get ids from remaining cols
    
    ## some future version of this script might want to save the daily undir_matrix to a list for additional analyses., i.e. mantel tests. 
    undir_matrix <- get_network(association_data = gbi, ## asnipe: turn gbi into undirected weighted adjacency matrix
                                data_format = "GBI",
                                association_index = "SRI") ## SRI = prop. of events shared b/w 2 mice over all events they participated in. 0 indicates no time spent together, 1 indicates all time spent together
    
    g <- graph.adjacency(undir_matrix,  ## create igraph object from undirected matrix
                         mode="undirected", ## undirected. https://www.baeldung.com/cs/graphs-directed-vs-undirected-graph
                         weighted = T, ## set to T. pulls SRI values from cells in the undir_matrix 
                         diag = F) 
    g <- simplify(g) ## removes multiple edges and loop edges
    
    # create vertex attributes
    V(g)$day <- bb
    V(g)$sex <- as.character(meta$sex[match(V(g)$name,meta$name)])
    V(g)$label <- as.character(meta$code[match(V(g)$name,meta$name)])
    V(g)$num_terr <- as.character(meta$num_terr[match(V(g)$name,meta$name)])
    V(g)$color = V(g)$sex #assign the "Sex" attribute as the vertex color
    V(g)$color = gsub("F","gold",V(g)$color) 
    V(g)$color = gsub("M","lightblue",V(g)$color) 
    V(g)$drop = as.character(meta$drop[match(V(g)$name,meta$name)])
    
    ## adjust node shapes based on drop status
    # V(g)$shape = V(g)$drop ## associate node shape with drop 1 or drop 2 status
    # V(g)$shape = gsub("1","circle",V(g)$shape) 
    # V(g)$shape = gsub("2","square",V(g)$shape)
    
    ## adjust node shapes based on drop status
    add_shape("triangle") ## add shapes
    mytriangle <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }
      
      symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
              stars=cbind(vertex.size, vertex.size, vertex.size),
              add=TRUE, inches=FALSE)
    } # triangle vertex shape
    add_shape("triangle", clip=shapes("circle")$clip, plot=mytriangle) # clips as a circle
    V(g)$shape = V(g)$num_terr
    V(g)$shape = gsub("0","circle",V(g)$shape) 
    V(g)$shape = gsub("1","triangle",V(g)$shape)
    V(g)$shape = gsub("2","square",V(g)$shape)
    
    ## get vertex edge strength measures
    V(g)$node_edge_strength <- graph.strength(g) ## node edge strength = sum of all edge weights for a single node (between 0-1)
    V(g)$node_degree_centrality <- degree(g, mode="all") ## node degree centrality  = number of connections to other nodes
    V(g)$node_eigen_centrality <- eigen_centrality(g)$vector ## node eigen centrality = measure of being connected to the well connected. not clear if this should be used or centr_eigen$vector
    V(g)$node_betweeness_centrality <- betweenness(g, directed = FALSE, weights = NA) ## node betweeness centrality = # of geodesics (shortest paths) going through a node
    V(g)$node_closeness_centrality <- closeness(g, mode = "all", normalized = TRUE, weights = NA) ## node closeness centrality = measure of # steps required to access every other node. 
    V(g)$node_page_rank <- page_rank(g)$vector ## 
    V(g)$node_authority_score <- authority_score(g)$vector ## 
    
    ## OppSex edge strength
    m_vertices <- V(g)[sex=="M"] ## get male vertex objects
    f_vertices <- V(g)[sex=="F"] ## get female vertex objects
    m_vertices_index <- which(V(g)$sex=="M") ## get male vertex index
    f_vertices_index <- which(V(g)$sex=="F") ## get female vertex index
    
    ## get males' female edge strength (i.e. only female weights)
    male_mf_edges_list <- list()
    cc=m_vertices[1] ## vertex object
    for(cc in m_vertices) {
      # print(cc)
      vertex_index <- which(V(g)$sex=="F") ## pull out the female index from the vertex object
      mf_edge <- E(g)[cc %--% vertex_index] ## get the edge objects using the vertex indices
      stats <- data.frame(matrix(ncol = 1,nrow = 1)) # create empty dataframe
      stats$focal <- V(g)[cc]$name ## get vertex objet name
      stats$mf_vertex_strength <- sum(E(g)[mf_edge]$weight)
      male_mf_edges_list[[cc]] <- stats ## write stats to list
    }
    male_mf_edges <- do.call("rbind", male_mf_edges_list)
    male_mf_edges <- select(male_mf_edges, focal, mf_vertex_strength)
    
    ## female mf edge weights
    female_mf_edges_list <- list()
    dd=f_vertices[1] ## vertex object
    for(dd in f_vertices) {
      # print(dd)
      vertex_index <- which(V(g)$sex=="M") ## pull out the male index from the vertex object
      mf_edge <- E(g)[dd %--% vertex_index] ## get the edge objects using the vertex indices
      stats <- data.frame(matrix(ncol = 1,nrow = 1)) # create empty dataframe
      stats$focal <- V(g)[dd]$name ## get vertex objet name
      stats$mf_vertex_strength <- sum(E(g)[mf_edge]$weight)
      female_mf_edges_list[[dd]] <- stats ## write stats to list
    }
    female_mf_edges <- do.call("rbind", female_mf_edges_list)
    female_mf_edges <- select(female_mf_edges, focal, mf_vertex_strength)
    
    oppsex_edge_strength <- rbind(male_mf_edges,female_mf_edges)
    V(g)$node_oppsex_edge_strength <- as.character(oppsex_edge_strength$mf_vertex_strength[match(V(g)$name,oppsex_edge_strength$focal)]) ## add oppsex edge strength to each mouse
    
    vertex_stats_list[[bb]] <- igraph::as_data_frame(g, "vertices") # write daily vertex stats to list
    
    ## Edge measures
    # edge_betweenness(g, directed = FALSE, weights = NA) ## edge betweeness = # of geodesics (shortest paths) passing through a focal edge
    
    ## create net stats
    net_stats <- data.frame(matrix(ncol = 1,nrow = 1)) # create empty dataframe
    graph_centrality = centr_degree(g, mode = "all", loops = T, normalized = T) # https://igraph.org/r/doc/centr_degree.html
    net_stats$trial <- trial
    net_stats$day <- bb
    net_stats$net_centrality <- graph_centrality$centralization ## net level centrality score based on node level centrality. https://igraph.org/r/doc/centralize.html
    net_stats$net_components_num <- components(g)$no ## components = # of network isolates disconnected from other isolates
    net_stats$net_eigen_centrality <- centr_eigen(g, directed =F, scale = T, normalized = T)$value ## https://igraph.org/r/doc/centr_eigen.html
    net_stats$net_mean_dist <- mean_distance(g, directed = FALSE) ## mean distance = avg. # of edges between any two nodes in network. 
    net_stats$net_edge_density <- edge_density(g, loops = FALSE) ## edge density = # of edges : # of all possible edges (ratio)
    net_stats$net_transitivity <- transitivity(g) ## transitivity = probability that adjacent nodes of a network are connected (aka clustering coefficient)
    net_stats$net_modularity_infomap <- modularity(cluster_infomap(g)) ## infomap method attempts to map the flow of information in a network, and the different clusters in which information may get remain for longer periods. Similar to walktrap, but not necessarily maximizing modularity, but rather the so-called "map equation".
    net_stats$net_modularity_infomap_group_n <- length(cluster_infomap(g))
    net_stats$net_modularity_fast_greedy <- modularity(cluster_fast_greedy(g))
    net_stats$net_modularity_fast_greedy_group_n <- length(cluster_fast_greedy(g))
    net_stats[1] <- NULL #clean net_stats
    net_stats_list[[bb]] <- net_stats ## write daily net stats to list
    
    ## graph
    svg(file=paste0("output/",trial,"_Night_", bb, "_sn.svg"), bg = "transparent")
    # png(file=paste0("output/",trial,"_Day_", i, "_sn.png"), bg = "white")
    par(mar = c(0.4, 0.1, 2, 0.1))
    plot(g,
         layout = layout.fruchterman.reingold,
         vertex.size = 20, 
         # vertex.size = scales::rescale(V(g)$node_edge_strength, to = c(8,30)), #rescale vertex strength to reasonable min/max size
         vertex.shape=V(g)$shape,
         # vertex.label= NA, ## Remove vertex labels
         edge.width = scales::rescale(E(g)$weight, to = c(0.5,15)), # rescale edge weight to reasonable min/max size
         edge.color = "darkgray",
         edge.curved = 0.2,
    )
    title(paste0("Night ", bb), cex.main=3) 
    dev.off()
  }
  
  ##  get vertex stats and triage dead mice
  vertex_stats <- do.call("rbind", vertex_stats_list)
  trial_vertex_stats_list[[aa]] <- vertex_stats %>% 
    arrange(name, day)
  
  ## get net stats
  trial_net_stats_list[[aa]] <- do.call("rbind", net_stats_list)
}
alltrial_net_stats <- do.call("rbind", trial_net_stats_list)
alltrial_vertex_stats <- do.call("rbind", trial_vertex_stats_list)

write.csv(alltrial_net_stats, "data/ALLTRIAL_SNA_net_stats.csv")
write.csv(alltrial_vertex_stats, "data/ALLTRIAL_SNA_node_stats.csv")
