#' Generate bootstrap index data
#'
#' @param celltype colnames of single cell data.Colnames should be labeled as cell type or cluster.
#' @param n  number of bootstrap
#' @import dplyr
#' @noRd
generate_boots <- function(celltype, n) {
  dt <- data.frame(cluster = celltype, id = seq_along(celltype))

  # 优化：预分配矩阵，避免多次cbind
  n_cells <- length(celltype)
  index <- matrix(0, nrow = n_cells, ncol = n)

  # 优化：使用split预先分组，避免每次循环都group_by
  split_indices <- split(dt$id, dt$cluster)

  for (i in seq_len(n)) {
    # 对每个cluster进行有放回抽样
    sampled_ids <- unlist(lapply(split_indices, function(ids) {
      sample(ids, length(ids), replace = TRUE)
    }), use.names = FALSE)
    index[, i] <- sampled_ids
  }

  return(index)
}




#
#' Calculate mean expression for one bootstrap
#'
#' @param i index
#' @param myseurat single cell Seurat object.METAFlux will calculate on "data" slot
#' @param samples generated bootstrap index data
#' @param myident Seurat idents.This will be a character string indicating the grouping of the seurat object
#' @noRd
get_ave_exp <- function(i, myseurat, samples, myident) {
  # 提取当前bootstrap的索引
  sample_idx <- samples[, i]

  meta.data <- myseurat@meta.data[sample_idx, , drop = FALSE]
  sample <- myseurat@assays$RNA@counts[, sample_idx, drop = FALSE]

  # === BUG FIX: 确保列名和行名匹配 ===
  # 当 sample_idx 包含重复索引时，meta.data 的行名会被 R 自动添加后缀 (.1, .2)
  # 但 sample (matrix) 的列名不会改变，导致 CreateSeuratObject 失败
  # 解决方案：重置行名和列名为一致的序列
  new_names <- paste0("cell_", seq_len(ncol(sample)))
  rownames(meta.data) <- new_names
  colnames(sample) <- new_names
  # === END BUG FIX ===

  # 创建Seurat对象
  SeuratObject <- suppressWarnings(
    Seurat::CreateSeuratObject(count = sample, meta.data = meta.data)
  )
  SeuratObject <- Seurat::NormalizeData(SeuratObject, verbose = FALSE)
  ave <- Seurat::AverageExpression(
    SeuratObject,
    group.by = myident,
    return.seurat = TRUE
  )[["RNA"]]@data

  return(ave)
}


#' Calculate bootstrap mean expression for single cell data
#'
#' @param n_bootstrap number of bootstrap
#' @param seed random seed
#' @param myseurat Seurat object. METAFlux will calculate on "data" slot
#' @param myident Seurat idents for grouping.This will be a character string indicating the grouping of the seurat object
#'
#' @return mean expression data
#' @export
calculate_avg_exp <- function(myseurat, myident, n_bootstrap, seed) {
  set.seed(seed)

  # 生成bootstrap索引
  samples <- generate_boots(myseurat@meta.data[, myident], n_bootstrap)

  # 优化：使用lapply并预分配，避免多次cbind
  exp <- lapply(seq_len(n_bootstrap), get_ave_exp,
    myseurat = myseurat,
    samples = samples,
    myident = myident
  )

  # 合并结果
  exp <- do.call(cbind, exp)

  return(exp)
}
