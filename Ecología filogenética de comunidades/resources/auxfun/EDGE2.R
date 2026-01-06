####### GE2 calculation #########

# generate 1,000,000 GE2 values distributed across all Red List categories 
# and sample to a distribution of GE2 values to each RL category to capture uncertainty
# this distribution of GE2 is then used to derive the GE2 for each species for each iteration of EDGE2 calculation
# a random value of GE2, corresponding to the appropriate RL category, can be assigned from the output to each species

pext.vals <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")) , pext = rev(c(0.97, 0.97/2, 0.97/4,0.97/8,0.97/16)))
GE.2.calc <- function(pext){
  iucn <- sample(1:5, size=1000000, replace=TRUE)
  data <- data.frame(species=c(1:1000000), pext=pext$pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext.tap <- c(0, pext$pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4; rank.quu <- rank^5
  model <- lm(pext.tap ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2; data$rank.cub <- data$rank^3; data$rank.qu <- data$rank^4; data$rank.quu <- data$rank^5
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  pext.LC <- data.frame(RL.cat = "LC", pext =data$rank.pext[data$pext == pext.tap[2]])
  pext.NT <- data.frame(RL.cat = "NT", pext =data$rank.pext[data$pext == pext.tap[3]])
  pext.VU <- data.frame(RL.cat = "VU", pext =data$rank.pext[data$pext == pext.tap[4]])
  pext.EN <- data.frame(RL.cat = "EN", pext =data$rank.pext[data$pext == pext.tap[5]])
  pext.CR <- data.frame(RL.cat = "CR", pext =data$rank.pext[data$pext == pext.tap[6]])
  pext.obj <- rbind(pext.CR,pext.EN, pext.VU, pext.NT, pext.LC)
  pext.obj.sample <- NULL
  for(i in unique(pext.obj$RL.cat)){
    a <- pext.obj$pext[pext.obj$RL.cat == i]
    if(median(a) < pext$pext[pext$rl.cat == i]){
      while(median(a) < pext$pext[pext$rl.cat == i]){
        a <- a[-sample(c(1:length(a)),50,prob = rev(a))]
      }
    }else{
      while(median(a) > pext$pext[pext$rl.cat == i]){
        a <- a[-sample(c(1:length(a)),50,prob = a)]
      }
    }
    pext.obj.sample <- rbind(pext.obj.sample,data.frame(RL.cat = i, pext = a))
  }
  return(pext.obj.sample)
}

####### EDGE2 calculation #########

# provide phylogenetic tree and dataframe with two columns: 
# the first comprising species names, the second comprising their associated GE2 scores (between 0 and 1)
# function returns three objects: 
# 1. dataframe with terminal branch length, GE2, ED2 and EDGE2 scores for each species
# 2. expected PD loss tree
# 3. PD and expected PD loss in MY for the clade

# remove excess pext, and reorders to same order as tree$tip.label
into_order <- function(tree, pext){
  new_pext <- pext[match(tree$tip.label, pext$species),]
  return (new_pext)
}

# order tree components
reorder_tree <- function(tree, ordering){
  tree@edge.length <- tree@edge.length[ordering]
  tree@edge <- tree@edge[ordering,]
  return(tree)
}

# EDGE2 Function
EDGE2_mod <- function(tree, pext){
  
  require(phylobase)
  require(data.table)
  
  names(pext) <- c("species","pext")
  
  N_species <- length(tree$tip.label)
  N_nodes <- tree$Nnode
  N_tot <- N_species + N_nodes
  
  # ensure extinction probabilities are given in same order as tree$tip.label.
  if (!identical(tree$tip.label, pext$species)){
    pext <- into_order(tree, pext)
  }
  
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  
  tree_dat <- data.frame(Species = as.character(tree$tip.label),
                         TBL = NA, 
                         pext = pext$pext, ED = NA, EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)
  
  tree <- as(tree, "phylo4")
  root <- rootNode(tree)
  nodes <- c(root, descendants(tree, root, "all"))
  
  # reorder tree components more instinctively, such that nodes are easier to find
  ord <- order(nodes)
  tree <- reorder_tree(tree, ord)
  nodes <- nodes[ord]
  
  tree_dat$TBL <- tree@edge.length[1:N_species]
  
  node_data <- data.frame(Node = 1:N_tot, Pext = rep(1, N_tot), Edge_Sum = NA)
  node_data[1:N_species, 2] <- pext[,2]
  
  # assign the product of its descendant tips to each node
  for (i in c(1:length(tree@label), N_tot:(root+1))){         # for each node, beginning with tips
    anc <- tree@edge[i,1]                                   # find ancestor of node
    node_data[anc, 2] <- node_data[anc, 2]*node_data[i,2]   # muliply ancestor value by node "pext" 
  }
  
  # multiply each edge.length by each pext caluclated above
  for(i in 1:length(nodes)){
    tree@edge.length[i] <- tree@edge.length[i]*node_data[i,2]
  }
  save(tree, file ="tree.rda")
  
  if (is.na(tree@edge.length[root])){
    tree@edge.length[root] <- 0
  }
  node_data$Edge_Sum[root] <- tree@edge.length[root]
  
  # for each internal node, summate ancesteral edgelengths
  for (i in (root+1):N_tot){
    ans <- tree@edge[i,1]
    node_data$Edge_Sum[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }
  
  # for each tip, summate ancesteral edgelengths to find EDGE2 score
  for (i in 1:N_species){
    ans <- tree@edge[i,1]
    tree_dat$EDGE[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }  
  
  tree_dat$ED <- tree_dat$EDGE / tree_dat$pext
  # reorder tree
  tree <- reorder_tree(tree, order(ord))
  
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}

# Esta función no estaba incluída en el archivo origina porque EDGE está diseñado para usarse con un
# único árbol, pero ahora es posible generar múltiples probabilidades de extinción para múltiples árboles
# -Dan

random_GE2 <- function(data, GE2, reps) {
  tmp.data <- data
  new.cols <- matrix(NA, nrow = nrow(data), ncol = reps)
  
  # PRE-PROCESAMIENTO: Organizar GE2 para no buscar en cada fila
  lookup <- split(GE2[, 2], GE2[, 1]) 
  all_vals <- GE2[, 2] # Pool completo para NE y DD
  
  cats_validas <- c("LC", "NT", "VU", "EN", "CR")
  cats_unk <- c("NE", "DD")
  
  # Extraemos las categorías de los datos una sola vez
  cats <- data[, 2] 
  
  for (n in 1:reps) {
    pext <- rep(NA, nrow(data)) # Vector vacío para llenar
    
    # 1. Caso: Categorías conocidas (LC, NT, etc.)
    cats_presentes <- intersect(unique(cats), cats_validas)
    
    for(cat in cats_presentes) {
      # Encuentra todas las filas que son de esta categoría
      idx <- which(cats == cat)
      
      if(!is.null(lookup[[cat]])) {
        pext[idx] <- sample(lookup[[cat]], length(idx), replace = TRUE)
      }
    }
    
    # 2. Caso: Desconocidas (NE, DD)
    idx_unk <- which(cats %in% cats_unk)
    if(length(idx_unk) > 0) {
      pext[idx_unk] <- sample(all_vals, length(idx_unk), replace = TRUE)
    }
    
    new.cols[, n] <- pext
  }
  
  colnames(new.cols) <- paste("pext_", 1:reps, sep = "")
  
  tmp.data <- cbind(tmp.data, as.data.frame(new.cols))
  return(tmp.data)
}