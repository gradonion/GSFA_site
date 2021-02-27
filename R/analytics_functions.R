ZZT_mean_diff <- function(A_mat, B_mat){
  if (nrow(A_mat) != nrow(B_mat)){
    stop("Please make sure the two matrices have the same number of rows (samples).")
  }
  zzt_diff <- A_mat %*% t(A_mat) - B_mat %*% t(B_mat)
  return(sqrt(sum(zzt_diff^2)) / nrow(zzt_diff))
}

factor_matrix_regression <- function(factors, G_mat){
  # factors: nxk, G_mat: nxm
  pval_tb <- matrix(nrow = ncol(factors), ncol = ncol(G_mat))
  beta_reg_tb <- pval_tb
  if (is.null(colnames(G_mat))){
    colnames(G_mat) <- 1:ncol(G_mat)
  }
  colnames(pval_tb) <- paste0("pval-", colnames(G_mat))
  colnames(beta_reg_tb) <- paste0("beta_reg-", colnames(G_mat))
  for (i in 1:ncol(factors)){
    for (j in 1:ncol(G_mat)){
      regress_mat <- data.frame(G = G_mat[, j], Z = factors[, i])
      lm.summary <- summary(lm(Z ~ G, data = regress_mat))
      pval_tb[i, j] <- lm.summary$coefficients[2, 4]
      beta_reg_tb[i, j] <- lm.summary$coefficients[2, 1]
    }
  }
  return(list(pval = pval_tb, beta = beta_reg_tb))
}

make_gibbs_res_tb <- function(gibbs_obj, G, compute_pve = TRUE){
  names(gibbs_obj) <- sub(pattern = "[.]", replacement = "_", x = names(gibbs_obj))
  K <- ncol(gibbs_obj$Z_pm)
  res_tb <- data.frame(index = 1:K, pi = gibbs_obj$pi_pm[, 1])
  if (compute_pve){
    sum_var <- rep(0, K)
    for (i in 1:K) {
      mat_tmp <- outer(gibbs_obj$Z_pm[, i], gibbs_obj$W_pm[, i])
      sum_var[i] <- sum(apply(mat_tmp, 2, var))
    }
    res_tb$sum_var <- sum_var
  }
  if (is.null(ncol(G))){
    if (length(G) != nrow(gibbs_obj$Z_pm)){
      stop("Number of samples in genotype and in the gibbs object do not match!")
    }
    res_tb$beta_pm <- gibbs_obj$beta_pm[, 1]
    res_tb$beta_reg <- rep(NA, K)
    res_tb$pval <- rep(NA, K)
    for (i in 1:K){
      lm.summary <- summary(lm(gibbs_obj$Z_pm[, i] ~ G))
      res_tb$beta_reg[i] <- lm.summary$coefficients[2, 1]
      res_tb$pval[i] <- lm.summary$coefficients[2, 4]
    }
  } else {
    if (nrow(G) != nrow(gibbs_obj$Z_pm)){
      stop("Number of samples in genotype and in the gibbs object do not match!")
    }
    res_regress <- factor_matrix_regression(gibbs_obj$Z_pm, G)
    tmp_beta_pm <- t(gibbs_obj$beta_pm)
    colnames(tmp_beta_pm) <- paste0("beta_pm-", 1:ncol(tmp_beta_pm))
    res_tb <- cbind(res_tb, tmp_beta_pm,
                    res_regress$beta, res_regress$pval)
  }
  return(res_tb)
}

make_flash_res_tb <- function(flashier_obj, G){
  K <- flashier_obj$n.factors
  flash_res_tb <- data.frame(matrix(nrow = K, ncol = 3))
  names(flash_res_tb) <- c('pve', 'pi1', 'lfsr_0.05')
  flash_res_tb$pve <- flashier_obj$pve
  for (i in 1:K){
    flash_res_tb$pi1[i] <- flashier_obj$fitted.g[[2]][[i]]$pi[2]
    factor_lfsr <- flashier_obj$loadings.lfsr[[2]][, i]
    flash_res_tb$lfsr_0.05[i] <- sum(factor_lfsr < 0.05) / length(factor_lfsr)
  }
  if (is.null(ncol(G))){
    flash_res_tb$beta_reg <- rep(NA, K)
    flash_res_tb$pval <- rep(NA, K)
    for (i in 1:K){
      flash.summary <- summary(lm(flashier_obj$loadings.pm[[1]][, i] ~ G))
      flash_res_tb$beta_reg[i] <- flash.summary$coefficients[2, 1]
      flash_res_tb$pval[i] <- flash.summary$coefficients[2, 4]
    }
  } else {
    res_regress <- factor_matrix_regression(flashier_obj$loadings.pm[[1]], G)
    flash_res_tb <- cbind(flash_res_tb, res_regress$beta, res_regress$pval)
  }
  return(flash_res_tb)
}

plot_pval_heatmap <- function(heatmap_matrix, factor_annot = NULL, snp_annot = NULL,
                              row_title = "Factor", column_title = "KO Perturbations"){
  # 'heatmap_matrix' has factors in the rows and SNPs in the columns
  # col_fun <- circlize::colorRamp2(c(0, 3, 15), c("blue3", "white", "firebrick"))
  col_fun <- circlize::colorRamp2(c(0, 3, 15), c("blue3", "white", "firebrick"))
  main_lgd <- Legend(col_fun = col_fun, title = "-log10(Factor~Guide\np value)", at = seq(0, 15, 3))
  lgd_list <- list(main = main_lgd)
  
  if (is.null(snp_annot)){
    column_annot <- NULL
    column_lgd <- NULL
  } else {
    snp_col_fun <- circlize::colorRamp2(c(-0.46, 0, 0.46), c("turquoise2", "white", "orange2"))
    column_annot <- columnAnnotation(beta = snp_annot,
                                     col = list(beta = snp_col_fun),
                                     show_annotation_name = T,
                                     annotation_name_side = "left",
                                     show_legend = F)
    column_lgd <- Legend(col_fun = snp_col_fun, title = "PLIER~SNP beta", at = seq(-0.4, 0.4, 0.2))
    lgd_list[["column"]] <- column_lgd
  }
  
  if (is.null(factor_annot)){
    row_annot <- NULL
    row_lgd <- NULL
  } else {
    # factor_col_fun <- circlize::colorRamp2(c(0, 0.7, 1), c("gold", "lightgreen", "seagreen"))
    factor_col_fun <- circlize::colorRamp2(c(0, 1), c("palegoldenrod", "turquoise4"))
    row_annot <- rowAnnotation(pi1 = factor_annot,
                               col = list(pi1 = factor_col_fun),
                               show_annotation_name = T,
                               show_legend = F)
    row_lgd <- Legend(col_fun = factor_col_fun, title = "Density (pi1)", at = seq(0, 1, 0.2))
    lgd_list[["row"]] <- row_lgd
  }
  
  # pd <- packLegend(main_lgd, column_lgd, row_lgd, direction = "vertical")
  ht <- Heatmap(-log10(heatmap_matrix), col = col_fun, name = "-log10(p value)",
                row_title = row_title, column_title = column_title,
                cluster_rows = F, cluster_columns = F,
                show_heatmap_legend = F,
                top_annotation = column_annot) + row_annot
  draw(ht, annotation_legend_list = lgd_list)
}

source("/project2/xinhe/yifan/GTEx/scripts/qqplot_uniform.R")
summ_pvalues <- function(pvalues, title_text = NULL){
  requireNamespace("gridExtra", quietly = TRUE)
  # distribution histogram
  plot1 <- histogram(pvalues, col = 'grey', type = "count",
                     xlim = c(0, 1), breaks = 50,
                     main = "p value distribution", xlab = "P value", ylab = "Count")
  # uniform qq-plot
  plot2 <- qqunif.plot(pvalues, main = "p value qq-plot")
  gridExtra::grid.arrange(plot1, plot2, ncol = 2, top = title_text)
}

plot_pairwise.corr_heatmap <- function(input_mat_1, input_mat_2 = NULL,
                                       name_1 = NULL, name_2 = NULL,
                                       corr_type = "pearson",
                                       return_corr = FALSE){
  # Please store samples in the columns of 'input_mat'.
  if (is.null(input_mat_2)){
    input_mat_2 <- input_mat_1
    name_2 <- name_1
  }
  stopifnot(nrow(input_mat_1) == nrow(input_mat_2))
  corr_mat <- matrix(nrow = ncol(input_mat_1), ncol = ncol(input_mat_2))
  
  for (i in 1:ncol(input_mat_1)){
    for (j in 1:ncol(input_mat_2)){
      vec_1 <- input_mat_1[, i]
      vec_2 <- input_mat_2[, j]
      if (corr_type == "pearson"){
        corr_mat[i, j] <- cor(vec_1, vec_2, method = "pearson")
      } else if (corr_type == "jaccard"){
        corr_mat[i, j] <- sum(vec_1 * vec_2) / sum(vec_1 + vec_2 > 0)
      } else {
        stop("Please provide a valid method to compute pairwise correlation.")
      }
    }
  }
  
  if (is.null(colnames(input_mat_1))){
    rownames(corr_mat) <- 1:ncol(input_mat_1)
  } else {
    rownames(corr_mat) <- colnames(input_mat_1)
  }
  if (is.null(colnames(input_mat_2))){
    colnames(corr_mat) <- 1:ncol(input_mat_2)
  } else {
    colnames(corr_mat) <- colnames(input_mat_2)
  }

  if (corr_type == "pearson"){
    ht <- Heatmap(corr_mat,
                  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  name = "Pearson Correlation",
                  row_title = name_1, column_title = name_2,
                  cluster_rows = F, cluster_columns = F,
                  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
  }
  if (corr_type == "jaccard"){
    ht <- Heatmap(corr_mat,
                  # col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "steelblue1", "steelblue4")),
                  col = circlize::colorRamp2(breaks = c(0, 0.5, 1),
                                             colors = c("black", "purple", "gold")),
                  name = "Jaccard Index",
                  row_title = name_1, column_title = name_2,
                  cluster_rows = F, cluster_columns = F,
                  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
  }
  draw(ht)
  if (return_corr){
    return(corr_mat)
  }
}

plot_PIP_hist_grid <- function(PIP_mat, cutoff = 0.5){
  plot.lst <- list()
  for (i in 1:ncol(PIP_mat)){
    pi1 <- mean(PIP_mat[, i] > cutoff)
    plot_df <- data.frame(prob = PIP_mat[, i])
    p <- ggplot(plot_df, aes_string("prob")) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept = cutoff, color = "red") +
      labs(title = paste0("Factor ", i, ": ", signif(pi1, digits = 2))) +
      theme(axis.title = element_blank())
    plot.lst[[i]] <- p
  }
  args <- c(plot.lst, list(nrow = 4,
                           left = "Count", bottom = "PIP"))
  do.call(grid.arrange, args)
}

print_enrich_tb <- function(enrich_list, qvalue_cutoff = 0.05, FC_cutoff = 2,
                            type = "per_factor"){
  signif_num <- rep(0, length(enrich_list))
  for (i in 1:length(enrich_list)){
    enrich_df <- enrich_list[[i]]@result
    tg_ratio <- enrich_df %>% pull(GeneRatio)
    tg_mat <- do.call(rbind, strsplit(tg_ratio, split = "/"))
    enrich_df$tg_1 <- as.numeric(tg_mat[, 1])
    enrich_df$tg_2 <- as.numeric(tg_mat[, 2])
    
    bg_ratio <- enrich_df %>% pull(BgRatio)
    bg_mat <- do.call(rbind, strsplit(bg_ratio, split = "/"))
    enrich_df$bg_1 <- as.numeric(bg_mat[, 1])
    enrich_df$bg_2 <- as.numeric(bg_mat[, 2])
    
    enrich_df <- enrich_df %>% mutate(FoldChange = (tg_1/tg_2) / (bg_1/bg_2))
    signif_tb <- enrich_df %>% filter(qvalue < qvalue_cutoff) %>%
      filter(FoldChange >= FC_cutoff) %>%
      select(ID, Description, GeneRatio, BgRatio, FoldChange, pvalue, qvalue, GS_size) %>%
      arrange(-FoldChange)
    signif_tb$FoldChange <- signif(signif_tb$FoldChange, digits = 3)
    signif_tb$pvalue <- format(signif_tb$pvalue, digits = 3)
    signif_tb$qvalue <- format(signif_tb$qvalue, digits = 3)
    
    signif_num[i] <- nrow(signif_tb)
    if (nrow(signif_tb) > 0){
      if (type == "per_factor"){
        caption_text <- paste("Factor", i, ":", nrow(signif_tb), "significant GO terms")
      }
      if (type == "per_marker"){
        caption_text <- paste(names(enrich_list)[i], ":", nrow(signif_tb), "significant GO terms")
      }
      print(kable(signif_tb, caption = caption_text) %>%
              kable_styling() %>%
              scroll_box(width = '100%', height = '500px'))
      cat("\n")
      cat("------------")
      cat("\n")
    }
  }
  return(signif_num)
}

permute_pval <- function(Z_pm, G_mat, num_perm = 1000){
  # G_mat: sample by condition
  # Z_pm: sample by factor
  pval_list <- list()
  for (i in 1:num_perm){
    if (i %% 100 == 0){
      print(i)
    }
    rand_order <- sample(1:nrow(Z_pm), size = nrow(Z_pm))
    rand_conditions <- G_mat[rand_order, ]
    pval_list[[i]] <- factor_matrix_regression(Z_pm, rand_conditions)$pval
  }
  return(pval_list)
}

plot_pval_list_grid <- function(pval_list, plot_by = "factor"){
  plot.lst <- list()
  if (plot_by == "factor"){
    K <- nrow(pval_list[[1]])
    for (k in 1:K){
      pval_vec <- unlist(lapply(pval_list, function(x){ x[k, ] }))
      pval_df <- data.frame(pval = pval_vec)
      plot.lst[[k]] <- ggplot(pval_df, aes_string("pval")) +
        geom_histogram(bins = 100) +
        labs(title = paste("Factor", k)) +
        theme(axis.title = element_blank())
    }
  } else if (plot_by == "marker"){
    K <- ncol(pval_list[[1]])
    if (is.null(colnames(pval_list[[1]]))){
      colnames(pval_list[[1]]) <- 1:K
    }
    for (k in 1:K){
      pval_vec <- unlist(lapply(pval_list, function(x){ x[k, ] }))
      pval_df <- data.frame(pval = pval_vec)
      plot.lst[[k]] <- ggplot(pval_df, aes_string("pval")) +
        geom_histogram(bins = 100) +
        labs(title = paste("Condition", colnames(pval_list[[1]])[k])) +
        theme(axis.title = element_blank())
    }
  } else {
    stop("Please provide a value for \'plot_by\', (coulde be either \'factor\' or \'marker\'.)")
  }
  args <- c(plot.lst, list(nrow = 5, left = "Count", bottom = "P value"))
  do.call(grid.arrange, args)
}

paired_pval_ranked_scatterplot <- function(pval_vec_1, pval_vec_2,
                                           name_1, name_2, zoom_in_y = NULL){
  paired_df <- data.frame(V1 = pval_vec_1[order(pval_vec_1)],
                          V2 = pval_vec_2[order(pval_vec_2)])
  paired_df <- reshape2::melt(paired_df, value.name = "pval")
  paired_df <- paired_df %>%
    mutate(neg_logp = -log10(pval)) %>%
    mutate(rank = rep(1:length(pval_vec_1), 2)) %>%
    rename(type = variable) %>%
    mutate(type = factor(type)) %>%
    mutate(type = recode(type, V1 = name_1, V2 = name_2))
  p1 <- ggplot(paired_df, aes(x = rank, y = neg_logp, color = type)) +
    geom_point(size = 0.8, alpha = 0.6) +
    labs(y = "-log10(P value)",
         title = paste(name_1, "vs", name_2, "Factor~KO Associations")) +
    theme(legend.title = element_blank())

  if (is.numeric(zoom_in_y)){
    p1 <- p1 +
      geom_hline(yintercept = zoom_in_y, color = "grey", linetype = "dashed") +
      theme(legend.position = "none")
    p2 <- ggplot(paired_df, aes(x = rank, y = neg_logp, color = type)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_y_continuous(limits = c(0, zoom_in_y)) +
      labs(y = "-log10(P value)",
           title = paste0("Zoomed in (Y-axis truncated at ", zoom_in_y, ")")) +
      theme(legend.title = element_blank())
    grid.arrange(p1, p2, nrow = 1)
  } else {
    print(p1)
  }
}
