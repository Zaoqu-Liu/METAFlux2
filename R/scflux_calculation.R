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
  # 输入验证
  if (sum(fraction) != 1) {
    stop("Sum of fractions must be equal to 1")
  }
  if (length(fraction) != num_cell) {
    stop("Number of cell clusters does not match with the length of fraction")
  }

  # 加载内部数据
  Hgem <- METAFlux:::Hgem
  mat <- Hgem$S
  reaction_name <- Hgem$Reaction
  names(reaction_name) <- NULL
  A_combined <- METAFlux:::A_combined

  # 优化：使用Matrix包的稀疏矩阵
  D <- Matrix::Matrix(0, nrow = 8378, ncol = 13082, sparse = TRUE)

  # 构建候选矩阵列表
  candi_list <- list(mat, D, D)

  # 构建矩阵选择器
  matrix_construct <- diag(1, nrow = num_cell, ncol = num_cell)
  matrix_construct[matrix_construct == 0] <- 2

  # 预分配列表
  celltype_matrix <- vector("list", num_cell)
  A_matrix <- vector("list", num_cell)

  # 优化：预先构建矩阵
  for (i in seq_len(num_cell)) {
    A_matrix[[i]] <- A_combined
    celltype_matrix[[i]] <- do.call(cbind, candi_list[matrix_construct[, i]])
  }

  message("Preparing for TME S matrix.....")

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
  reaction_name[exchange_idx] <- paste0(
    "internal_medium ",
    reaction_name[exchange_idx]
  )

  # 优化：向量化构建细胞类型反应名称
  cell_reaction <- lapply(seq_len(num_cell), function(i) {
    paste(paste("celltype", i), reaction_name)
  })

  construct_reaction_names <- c(unlist(cell_reaction), external)

  # 预计算P1矩阵
  P1 <- Matrix::Diagonal(ncol(final_s), 1)

  message("S matrix completed......")

  # 预分配flux_vector
  n_bootstraps <- ncol(fluxscore) / num_cell
  flux_vector <- vector("list", n_bootstraps)

  message("Compute metabolic flux......")

  # 优化：预计算一些不变的部分
  Seq <- seq(1, ncol(fluxscore), by = num_cell)
  pb <- utils::txtProgressBar(0, length(Seq), style = 3)

  # 预计算medium匹配索引
  tail_idx <- tail(seq_len(ncol(final_s)), 1648)

  # 预计算q的位置
  q_positions <- seq(13015, ncol(final_s), by = 13082)

  # 预设osqp设置
  settings <- osqp::osqpSettings(
    max_iter = 1000000L,
    eps_abs = 1e-04,
    eps_rel = 1e-04,
    verbose = FALSE,
    adaptive_rho_interval = 50
  )

  for (t in Seq) {
    utils::setTxtProgressBar(pb, match(t, Seq))

    # 提取当前bootstrap的分数
    score <- fluxscore[, c(t:(t + num_cell - 1)), drop = FALSE]

    # 构建P矩阵
    P_diag <- c(
      rep(fraction, each = 13082),
      rep(1, 1648)
    )
    P <- Matrix::Diagonal(ncol(final_s), P_diag)

    # 计算加权的S矩阵
    fraction_finals <- final_s %*% P

    # 构建q向量
    q <- rep(0, ncol(final_s))
    q[q_positions] <- -10000 * fraction

    # 构建A矩阵
    A <- rbind(fraction_finals, P1)

    # 构建ras向量
    ras <- c(as.vector(unlist(score[, seq_len(num_cell)])), rep(1, 1648))

    # 构建边界
    origlb <- c(rep(Hgem$LB, num_cell), rep(-1, 1648))

    # 设置可逆反应的下界
    rev_indices <- rep(Hgem$rev == 1, num_cell)
    origlb[rev_indices] <- -ras[rev_indices]

    # 设置不可逆反应的下界
    non_rev_indices <- rep(Hgem$rev == 0, num_cell)
    origlb[non_rev_indices] <- 0

    # 设置上界
    origub <- ras

    # 设置tail部分的下界
    origlb[tail_idx] <- 0

    # 优化：预先计算medium匹配
    matches <- unlist(lapply(medium$reaction_name, function(x) {
      intersect(
        which(stringi::stri_detect_fixed(construct_reaction_names, x)),
        tail_idx
      )
    }))

    origlb[matches] <- -1

    # 构建最终约束
    l <- c(rep(0, nrow(final_s)), origlb)
    u <- c(rep(0, nrow(final_s)), origub)

    # 创建并求解模型
    model <- osqp::osqp(P, q, A, l, u, settings)
    res <- model$Solve()

    # 存储结果
    flux_vector[[match(t, Seq)]] <- res$x
  }

  close(pb)

  # 合并结果
  flux_matrix <- as.data.frame(do.call(cbind, flux_vector))
  rownames(flux_matrix) <- construct_reaction_names

  return(flux_matrix)
}
