# Load the required packages.
library(foreach)
library(doParallel)
library(Matrix)

# Generate a random network. 
# Input: adjacency matrix adj_matrix, number of rewiring attempts per edge Q, and number of random networks to generate n_random.
generate_random_networks <- function(adj_matrix, n_random = 20, Q = 10, nThreads = 1) {
  setup_parallel <- function(nThreads) {
    if (foreach::getDoParWorkers() > 1) {
      if (nThreads > 1 & foreach::getDoParWorkers() != nThreads) {
        message("Parallel backend already registered and is inconsistent with nThreads.")
      }
      return(foreach::getDoParWorkers())
    } else {
      if (nThreads > 1) {
        cl <- parallel::makeCluster(nThreads, outfile = "")
        doParallel::registerDoParallel(cl)
        return(cl)
      }
    }
    return(NULL)
  }
  
  cl <- setup_parallel(nThreads)
  
  edges <- which(adj_matrix == 1, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], ]
  
  random_networks <- foreach::foreach(i = 1:n_random, .packages = c("Matrix")) %dopar% {
    set.seed(i) 
    
    random_adj <- adj_matrix
    
    num_swaps <- Q * nrow(edges)
    
    current_edges <- which(random_adj == 1, arr.ind = TRUE)
    current_edges <- current_edges[current_edges[, 1] < current_edges[, 2], ]
    
    for (j in 1:num_swaps) {
      selected_indices <- sample(1:nrow(current_edges), 2)
      edge1 <- current_edges[selected_indices[1], ]
      edge2 <- current_edges[selected_indices[2], ]
      
      A <- edge1[1]
      B <- edge1[2]
      C <- edge2[1]
      D <- edge2[2]
      
      if (A != D && C != B && A != C && B != D && 
          random_adj[A, D] == 0 && random_adj[C, B] == 0) {
        
        random_adj[A, B] <- random_adj[B, A] <- 0
        random_adj[C, D] <- random_adj[D, C] <- 0
        
        random_adj[A, D] <- random_adj[D, A] <- 1
        random_adj[C, B] <- random_adj[B, C] <- 1

        current_edges[selected_indices[1], ] <- c(min(A, D), max(A, D))
        current_edges[selected_indices[2], ] <- c(min(C, B), max(C, B))
      }
    }
    
    return(random_adj)
  }
  
  if (!is.null(cl) && nThreads > 1 && foreach::getDoParWorkers() > 1) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  
  return(random_networks)
}

# Construct the adjacency matrix.
build_adjacency_matrix <- function(net) {
  genes <- unique(c(as.character(net$gene1), as.character(net$gene2)))
  adj.mat <- matrix(0, nrow=length(genes), ncol=length(genes), dimnames=list(genes, genes))
  
  for (i in seq_len(nrow(net))) {
    gene1 <- net$gene1[i]
    gene2 <- net$gene2[i]
    adj.mat[gene1, gene2] <- 1
    adj.mat[gene2, gene1] <- 1
  }
  
  return(adj.mat)
}

adj.mat <- build_adjacency_matrix(net = TAIR_net)

random_TAIR <- generate_random_networks(adj.mat, Q = 100, n_random = 100, nThreads = 10)

# Define the smoothRWR function.
smoothRWR <- function(pcc, network, gamma = 0.7, nThreads = 1) {
  if (!(is(network, "matrix") || is(network, "sparseMatrix"))) {
    stop("'network' must be an adjacency matrix")
  }
  
  rownames(pcc) <- toupper(rownames(pcc))
  colnames(pcc) <- toupper(colnames(pcc))
  rownames(network) <- toupper(rownames(network))
  colnames(network) <- toupper(colnames(network))
  
  gene.overlap <- intersect(rownames(network), rownames(pcc))
  
  ep.data <- matrix(rep(0, nrow(network) * dim(pcc)[2]), nrow = nrow(network))
  rownames(ep.data) <- rownames(network)
  colnames(ep.data) <- colnames(pcc)
  
  ep.data[gene.overlap, ] <- pcc[gene.overlap, ]
  
  Anorm <- network / Matrix::rowSums(network)
  
  eye <- diag(dim(network)[1])
  
  AA <- Matrix::t(eye - gamma * Anorm)
  BB <- (1 - gamma) * ep.data
  
  smooth.data <- solve(AA, BB)
  
  smooth.pcc <- pcc
  smooth.pcc[gene.overlap, ] <- as.matrix(smooth.data[gene.overlap, ])
  
  return(smooth.pcc)
}

# Main function.
calculate_p_values_t_test <- function(seed_genes, pcc_values, original_network, random_networks, gamma = 0.7, nThreads = 1) {
  pcc_matrix <- matrix(pcc_values, ncol = 1, dimnames = list(seed_genes, "PCC"))
  
  original_adj <- build_adjacency_matrix(original_network)
  actual_result <- smoothRWR(pcc = pcc_matrix, network = original_adj, gamma = gamma, nThreads = nThreads)
  actual_result <- as.numeric(actual_result[rownames(actual_result), 1])  # 提取平滑后的 PCC 值
  
  random_results <- matrix(NA, nrow = length(seed_genes), ncol = length(random_networks))
  rownames(random_results) <- seed_genes
  
  if (foreach::getDoParWorkers() > 1) {
    message("Parallel backend already registered. Using existing backend.")
  } else {
    if (nThreads > 1) {
      cl <- parallel::makeCluster(nThreads, outfile = "")
      doParallel::registerDoParallel(cl)
      on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      })
    }
  }
  
  random_results <- foreach(
    i = seq_along(random_networks),
    .combine = cbind,
    .packages = c("Matrix", "methods"),
    .export = c("smoothRWR", "build_adjacency_matrix", "pcc_matrix", "gamma")
  ) %dopar% {
    current_network <- random_networks[[i]]
    result <- smoothRWR(pcc = pcc_matrix, network = current_network, gamma = gamma, nThreads = 1)
    result <- as.numeric(result[rownames(result), 1])  # 提取平滑后的 PCC 值
    return(result)
  }
  
  pvals <- sapply(1:length(seed_genes), function(i) {
    if (length(unique(random_results[i, ])) == 1) {
      return(NA)
    }
    t_test_result <- t.test(x = random_results[i, ], mu = actual_result[i], alternative = "two.sided")
    return(t_test_result$p.value)
  })
  
  names(pvals) <- seed_genes
  
  adjusted_pvals <- p.adjust(pvals, method = "BH")
  
  return(data.frame(
    gene = seed_genes,
    actual_probability = actual_result,
    pval = pvals,
    adjusted_pval = adjusted_pvals
  ))
}


# Example data.
gene_list <- TAIR_pcc

original_network <- TAIR_net

random_networks <- random_TAIR

TAIR_rwr_ <- calculate_p_values_t_test(
  seed_genes = gene_list$gene,
  pcc_values = gene_list$pcc,
  original_network = original_network,
  random_networks = random_networks,
  gamma = 0.3,
  nThreads = 5
)
