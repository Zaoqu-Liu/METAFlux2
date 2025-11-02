#' Calculate reaction score for isoenzymes (or relationship only)
#'
#' @param x index
#' @param data input gene expression data
#' @param list isoenzymes reaction list
#' @param gene_num file of number of times one gene has involved in all pathways
#'
#' @return Calculated isoenzyme list scores
#' @noRd
calculate_iso_score <- function(x, data, list, gene_num) {
  data_feature <- list[[x]][list[[x]] %in% rownames(data)]

  if (length(data_feature) > 0) {
    vec <- gene_num[data_feature, ]$V1
    expr <- as.matrix(data[data_feature, , drop = FALSE])

    # 优化：向量化操作，避免sweep
    norm <- colSums(expr / vec, na.rm = TRUE)
    return(norm)
  } else {
    # 优化：直接返回NA向量
    return(rep(NA_real_, ncol(data)))
  }
}




#' Calculate reaction score for enzyme complex(and relationship only)
#'
#' @param x index
#' @param data input gene expression data
#' @param list enzyme complex reaction list
#' @param gene_num file of number of times one gene has involved in all pathways
#'
#' @return Calculated enzyme complex list scores
#' @noRd
calculate_simple_comeplex_score <- function(x, data, list, gene_num) {
  data_feature <- list[[x]][list[[x]] %in% rownames(data)]

  if (length(data_feature) > 0) {
    vec <- gene_num[data_feature, ]$V1
    expr <- as.matrix(data[data_feature, , drop = FALSE])

    # 优化：向量化的最小值计算
    norm <- apply(expr / vec, 2, min, na.rm = TRUE)
    return(norm)
  } else {
    return(rep(NA_real_, ncol(data)))
  }
}


#' calculate reaction score for complicated complex(with both 'and' and 'or' relationship)
#'
#' @param x index
#' @param data input gene expression data
#' @param gene_num file of number of times one gene has involved in all pathways
#'
#' @return Calculated complicated enzyme complex list scores
#' @noRd
calculate_multi_comp <- function(x, data, gene_num) {
  reaction <- multi_comp[[x]]

  # 提取括号内容
  bracket_contents <- stringr::str_extract_all(reaction, "\\([^()]+\\)")[[1]]

  # 检测是否包含'or'
  com <- stringi::stri_detect_fixed(bracket_contents, "or")

  # 去除括号
  c <- gsub("\\)", "", gsub("\\(", "", bracket_contents))

  if (unique(com) == TRUE) {
    # 处理包含'or'的情况
    newiso <- lapply(c, function(x) {
      trimws(unlist(strsplit(x, "or", fixed = TRUE)))
    })

    sum_score <- do.call(rbind, lapply(
      seq_along(newiso),
      calculate_iso_score,
      data = data,
      list = newiso,
      gene_num = gene_num
    ))

    # 提取不包含'or'的特征
    split_parts <- unlist(strsplit(reaction, "and", fixed = TRUE))
    feature <- trimws(split_parts[!stringi::stri_detect_fixed(split_parts, "or")])

    data_feature <- feature[feature %in% rownames(data)]

    if (length(data_feature) > 0) {
      vec <- gene_num[data_feature, ]$V1
      expr <- as.matrix(data[data_feature, , drop = FALSE])
      whole_score <- rbind(expr / vec, sum_score)
    } else {
      whole_score <- sum_score
    }

    norm <- apply(whole_score, 2, min, na.rm = TRUE)
  } else if (unique(com) == FALSE) {
    # 处理包含'and'的情况
    newiso <- lapply(c, function(x) {
      trimws(unlist(strsplit(x, "and", fixed = TRUE)))
    })

    sum_score <- do.call(rbind, lapply(
      seq_along(newiso),
      calculate_simple_comeplex_score,
      data = data,
      list = newiso,
      gene_num = gene_num
    ))

    # 提取不包含'and'的特征
    split_parts <- unlist(strsplit(reaction, "or", fixed = TRUE))
    feature <- trimws(split_parts[!stringi::stri_detect_fixed(split_parts, "and")])

    data_feature <- feature[feature %in% rownames(data)]

    if (length(data_feature) > 0) {
      vec <- gene_num[data_feature, ]$V1
      expr <- as.matrix(data[data_feature, , drop = FALSE])
      upper_score <- expr / vec
      whole_score <- rbind(upper_score, sum_score)
    } else {
      whole_score <- sum_score
    }

    norm <- colSums(whole_score, na.rm = TRUE)
  }

  return(norm)
}



#' Normalize scores
#'
#' @param x index
#' @param ... additional parameters to `max()`
#' @noRd
stdize <- function(x, ...) {
  x / max(x, ...)
}


#' Calculate metabolic reaction scores (MRAS) for 13082 reactions
#'
#' @param data gene expression data.1.The gene expression matrix should be gene by sample matrix where row names are human gene names (gene symbols),
#' and column names should be sample names. Please note that METAFlux does not support other gene IDs.
#' 2.The input gene expression matrix should be normalized (e.g., log-transformed, etc.) before using METAFlux.
#' METAflux will not perform any normalization on expression data.
#' 3.Gene expression data cannot have negative values.
#' @export
calculate_reaction_score <- function(data) {
  # 输入验证
  if (sum(data < 0) > 0) {
    stop("Expression data needs to be all positive")
  }

  # 确保特征存在
  features <- rownames(data)
  if (sum(features %in% rownames(gene_num)) == 0) {
    stop("Requested gene names cannot be found. Rownames of input data should be human gene names. Please check the rownames of input data.")
  }

  message(paste0(
    round(sum(features %in% rownames(gene_num)) / 3625 * 100, 3),
    "% metabolic related genes were found......"
  ))

  # 加载内部数据
  gene_num <- METAFlux:::gene_num
  Hgem <- METAFlux:::Hgem
  iso <- METAFlux:::iso
  multi_comp <- METAFlux:::multi_comp
  simple_comp <- METAFlux:::simple_comp

  message("Computing metabolic reaction activity scores......")

  # 优化：使用vapply或并行计算（如果列表很大）
  # 计算三种类型的分数
  n_iso <- length(iso)
  n_simple <- length(simple_comp)
  n_multi <- length(multi_comp)

  # 预分配矩阵
  core <- matrix(NA_real_, nrow = n_iso, ncol = ncol(data))
  core2 <- matrix(NA_real_, nrow = n_simple, ncol = ncol(data))
  core3 <- matrix(NA_real_, nrow = n_multi, ncol = ncol(data))

  # 计算iso分数
  for (i in seq_len(n_iso)) {
    core[i, ] <- calculate_iso_score(i, data = data, list = iso, gene_num = gene_num)
  }

  # 计算simple_comp分数
  for (i in seq_len(n_simple)) {
    core2[i, ] <- calculate_simple_comeplex_score(i, data, list = simple_comp, gene_num = gene_num)
  }

  # 计算multi_comp分数
  for (i in seq_len(n_multi)) {
    core3[i, ] <- calculate_multi_comp(i, data = data, gene_num = gene_num)
  }

  message("Preparing for score matrix......")

  # 设置行名
  rownames(core) <- names(iso)
  rownames(core2) <- names(simple_comp)
  rownames(core3) <- names(multi_comp)

  # 合并分数矩阵
  big_score_matrix <- rbind(core, core2, core3)

  # 优化：向量化标准化，避免apply
  max_vals <- apply(big_score_matrix, 2, max, na.rm = TRUE)
  big_score_matrix <- sweep(big_score_matrix, 2, max_vals, "/")

  # 处理NaN和NA
  big_score_matrix[is.na(big_score_matrix)] <- 0

  # 准备最终数据框
  empty_helper <- data.frame(reaction = Hgem$Reaction, stringsAsFactors = FALSE)

  # 转换为数据框进行合并
  big_score_df <- as.data.frame(big_score_matrix)
  big_score_df$reaction <- rownames(big_score_df)

  Final_df <- merge(
    empty_helper,
    big_score_df,
    all.x = TRUE,
    by = "reaction"
  )

  # 处理NA
  Final_df[is.na(Final_df)] <- 1

  # 设置行名并删除reaction列
  rownames(Final_df) <- Final_df$reaction
  Final_df$reaction <- NULL

  # 重新排序
  Final_df <- Final_df[Hgem$Reaction, , drop = FALSE]

  # 验证
  if (all.equal(rownames(Final_df), Hgem$Reaction)) {
    message("Metabolic reaction activity scores successfully calculated \n")
  } else {
    message("Calculation not reliable. Check input data format \n")
  }

  # 将biomass reaction设为0
  Final_df[which(Hgem$LB == 0 & Hgem$UB == 0), ] <- 0

  return(Final_df)
}
