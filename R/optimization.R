#' Final optimization step for flux calculation
#'
#' @param mras metabolic reaction activity scores
#' @param medium input medium file which indicates the nutrients available in the medium.
#' We provide 2 general mediums if you have no prior knowledge about your medium: cell line medium and human blood medium if prior knowledge is not available.
#' Please see tutorial for more details.
#'
#' @return Calculated fluxes
#' @export
compute_flux <- function(mras, medium) {
  message("Setting up for optimization.....")

  # 加载内部数据
  Hgem <- METAFlux:::Hgem

  # 优化：预计算不变的部分
  n_reactions <- ncol(Hgem$S)
  n_metabolites <- nrow(Hgem$S)

  P <- Matrix::Diagonal(n_reactions, 1)
  q <- rep(0, n_reactions)
  q[which(Hgem$Obj == 1)] <- -10000

  A <- rbind(Hgem$S, P)

  # 预计算exchange reactions索引
  exchange_idx <- which(Hgem$pathway == "Exchange/demand reactions")
  medium_idx <- which(Hgem$Reaction %in% medium$reaction_name)
  rev_idx <- which(Hgem$rev == 1)

  # 预设置osqp参数（不变部分）
  settings <- osqp::osqpSettings(
    max_iter = 1000000L,
    eps_abs = 1e-04,
    eps_rel = 1e-04,
    adaptive_rho_interval = 50,
    verbose = FALSE
  )

  # 预分配结果列表
  n_samples <- ncol(mras)
  flux_vector <- vector("list", n_samples)

  message("Computing bulk RNA-seq flux.....")

  # 优化：使用更简洁的进度条
  pb <- utils::txtProgressBar(0, n_samples, style = 3)

  for (i in seq_len(n_samples)) {
    utils::setTxtProgressBar(pb, i)

    # 计算边界
    origlb <- Hgem$LB
    origlb[rev_idx] <- -mras[rev_idx, i]
    origlb[!rev_idx] <- 0
    origlb <- origlb[, 1]

    origub <- mras[, i]

    # 设置exchange reactions边界
    origlb[exchange_idx] <- 0
    origlb[medium_idx] <- -1

    # 构建约束
    l <- c(rep(0, n_metabolites), origlb)
    u <- c(rep(0, n_metabolites), origub)

    # 创建并求解模型
    model <- osqp::osqp(P, q, A, l, u, settings)
    res <- model$Solve()

    # 存储结果
    flux_vector[[i]] <- res$x
  }

  close(pb)

  # 优化：使用do.call一次性合并
  flux_matrix <- do.call(cbind, flux_vector)
  colnames(flux_matrix) <- colnames(mras)
  rownames(flux_matrix) <- Hgem$Reaction

  return(flux_matrix)
}
