#' sc-rna seq flux calculation
#'
#' @param num_cell the number of cell types or clusters
#' @param fraction fraction of each cell types. Fractions need to sum up to 1
#' @param fluxscore calculated metabolic activity score(MRAS)
#' @param medium Medium profile
#' @import Seurat
#' @import utils
#'
#' @return Flux score for single cell data
#' @export
compute_sc_flux <- function(num_cell, fraction, fluxscore, medium) {
  # ========== 输入验证 ==========
  if (sum(fraction) != 1) {
    stop("Sum of fractions must be equal to 1")
  }
  if (length(fraction) != num_cell) {
    stop("Number of cell clusters does not match with the length of fraction")
  }

  # ========== 加载内部数据 ==========
  Hgem <- METAFlux:::Hgem
  mat <- Hgem$S
  reaction_name <- Hgem$Reaction
  names(reaction_name) <- NULL
  A_combined <- METAFlux:::A_combined

  # ========== 构建基础矩阵 ==========
  message("Preparing for TME S matrix.....")
  
  # 使用稀疏矩阵
  D <- Matrix::Matrix(0, nrow = 8378, ncol = 13082, sparse = TRUE)
  candi_list <- list(mat, D, D)
  
  # 构建矩阵选择器
  matrix_construct <- diag(1, nrow = num_cell, ncol = num_cell)
  matrix_construct[matrix_construct == 0] <- 2
  
  # 预分配列表
  celltype_matrix <- vector("list", num_cell)
  A_matrix <- vector("list", num_cell)
  
  # 构建矩阵
  for (i in seq_len(num_cell)) {
    A_matrix[[i]] <- A_combined
    celltype_matrix[[i]] <- do.call(cbind, candi_list[matrix_construct[, i]])
  }
  
  # 构建最终的S矩阵
  final_s <- rbind(
    do.call(cbind, A_matrix),
    do.call(rbind, celltype_matrix)
  )
  
  # 构建外部介质矩阵
  whole3 <- rbind(
    Matrix::Diagonal(1648, -1),
    Matrix::Matrix(0, nrow = num_cell * nrow(mat), ncol = 1648, sparse = TRUE)
  )
  
  final_s <- cbind(final_s, whole3)
  
  # 构建反应名称
  exchange_idx <- which(Hgem$pathway == "Exchange/demand reactions")
  external <- paste("external_medium", reaction_name[exchange_idx])
  reaction_name[exchange_idx] <- paste0("internal_medium ", reaction_name[exchange_idx])
  
  # 向量化构建细胞类型反应名称
  cell_reaction <- lapply(seq_len(num_cell), function(i) {
    paste(paste("celltype", i), reaction_name)
  })
  construct_reaction_names <- c(unlist(cell_reaction), external)
  
  message("S matrix completed......")

  # ========== 关键优化：预计算所有不变量 ==========
  message("Pre-computing invariants...")
  
  n_bootstraps <- ncol(fluxscore) / num_cell
  Seq <- seq(1, ncol(fluxscore), by = num_cell)
  
  # 1. P 矩阵（fraction 不变，所以 P 不变）
  P_diag <- c(rep(fraction, each = 13082), rep(1, 1648))
  P <- Matrix::Diagonal(ncol(final_s), P_diag)
  
  # 2. fraction_finals（最耗时的矩阵乘法，只算一次）
  fraction_finals <- final_s %*% P
  
  # 3. q 向量（不变）
  q_positions <- seq(13015, ncol(final_s), by = 13082)
  q <- rep(0, ncol(final_s))
  q[q_positions] <- -10000 * fraction
  
  # 4. A 矩阵（不变）
  P1 <- Matrix::Diagonal(ncol(final_s), 1)
  A <- rbind(fraction_finals, P1)
  
  # 5. 索引向量（不变）
  rev_indices <- rep(Hgem$rev == 1, num_cell)
  non_rev_indices <- rep(Hgem$rev == 0, num_cell)
  tail_idx <- tail(seq_len(ncol(final_s)), 1648)
  
  # 6. 常量向量（不变）
  LB_rep <- rep(Hgem$LB, num_cell)
  zero_vec_final_s <- rep(0, nrow(final_s))
  neg_one_vec <- rep(-1, 1648)
  one_vec <- rep(1, 1648)
  
  # 7. medium 匹配索引（最关键的优化！）
  message("Pre-computing medium matches...")
  medium_matches <- unlist(lapply(medium$reaction_name, function(x) {
    intersect(
      which(stringi::stri_detect_fixed(construct_reaction_names, x)),
      tail_idx
    )
  }), use.names = FALSE)
  
  # 8. OSQP 设置（保持默认设置）
  settings <- osqp::osqpSettings(
    max_iter = 1000000L,
    eps_abs = 1e-04,
    eps_rel = 1e-04,
    verbose = FALSE,
    adaptive_rho_interval = 50
  )

  # ========== 主循环：使用预计算的常量，每次创建新模型 ==========
  message("Computing metabolic flux...")
  flux_vector <- vector("list", n_bootstraps)
  pb <- utils::txtProgressBar(0, n_bootstraps, style = 3)
  
  for (i in seq_len(n_bootstraps)) {
    utils::setTxtProgressBar(pb, i)
    
    t <- Seq[i]
    
    # 提取当前 bootstrap 的分数
    score <- fluxscore[, c(t:(t + num_cell - 1)), drop = FALSE]
    
    # 计算 ras（唯一变化的部分）
    ras <- c(as.vector(unlist(score[, seq_len(num_cell)])), one_vec)
    
    # 构建边界（使用预计算的常量和索引）
    origlb <- c(LB_rep, neg_one_vec)
    origlb[rev_indices] <- -ras[rev_indices]
    origlb[non_rev_indices] <- 0
    origlb[tail_idx] <- 0
    origlb[medium_matches] <- -1  # 使用预计算的匹配索引（关键优化）
    
    origub <- ras
    
    l <- c(zero_vec_final_s, origlb)
    u <- c(zero_vec_final_s, origub)
    
    # 创建并求解模型（使用预计算的 P, q, A 矩阵）
    model <- osqp::osqp(P, q, A, l, u, settings)
    res <- model$Solve()
    
    # 存储结果
    flux_vector[[i]] <- res$x
  }
  
  close(pb)

  # ========== 合并结果 ==========
  flux_matrix <- as.data.frame(do.call(cbind, flux_vector))
  rownames(flux_matrix) <- construct_reaction_names

  return(flux_matrix)
}
