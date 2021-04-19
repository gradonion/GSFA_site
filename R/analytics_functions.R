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
    # res_tb$beta_pm <- gibbs_obj$beta_pm[, 1]
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
    # tmp_beta_pm <- t(gibbs_obj$beta_pm)
    # colnames(tmp_beta_pm) <- paste0("beta_pm-", 1:ncol(tmp_beta_pm))
    res_tb <- cbind(res_tb, res_regress$beta, res_regress$pval)
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

dotplot_beta_PIP <- function(beta_pip_matrix, beta_pm_matrix,
                             marker_names, reorder_markers = marker_names,
                             exclude_offset = TRUE, return_dataframe = FALSE){
  # Both 'beta_pip_matrix' and 'beta_pm_matrix' should be factor by guide/marker matrices
  if (exclude_offset){
    beta_pip_matrix <- beta_pip_matrix[, -ncol(beta_pip_matrix)]
    beta_pm_matrix <- beta_pm_matrix[, -ncol(beta_pm_matrix)]
  }
  rownames(beta_pip_matrix) <- 1:nrow(beta_pip_matrix)
  colnames(beta_pip_matrix) <- marker_names
  beta_pip_df <- as.data.frame(beta_pip_matrix)
  beta_pip_df$Factor <- paste0("Factor ", 1:nrow(beta_pip_df))
  beta_pip_plot_df <- reshape2::melt(beta_pip_df, value.name = "PIP")
  
  rownames(beta_pm_matrix) <- 1:nrow(beta_pm_matrix)
  colnames(beta_pm_matrix) <- marker_names
  beta_pm_df <- as.data.frame(beta_pm_matrix)
  beta_pm_df$Factor <- paste0("Factor ", 1:nrow(beta_pm_df))
  beta_pm_plot_df <- reshape2::melt(beta_pm_df, id.var = "Factor",
                                    variable.name = "Perturbation",
                                    value.name = "Estimated effect size")
  # beta_pm_plot_df$PIP <- beta_pip_plot_df$PIP
  beta_pm_plot_df <- beta_pm_plot_df %>%
    mutate(PIP = beta_pip_plot_df$PIP,
           Factor = factor(Factor, levels = paste0("Factor ", nrow(beta_pip_df):1)),
           Perturbation = factor(Perturbation, levels = reorder_markers))
  plot_out <- ggplot(beta_pm_plot_df) +
    geom_point(aes(x = Perturbation, y = Factor,
                   size = PIP, color = `Estimated effect size`)) +
    scale_color_gradient2(low = "purple3", mid = "grey90", high = "darkorange1") +
    # scale_color_gradientn(colors = c("purple", "grey90", "darkorange1"),
    #                       values = scales::rescale(c(-0.6, 0, 0.6))) +
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12))
  print(plot_out)
  if (return_dataframe){
    return(beta_pm_plot_df)
  }
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

plot_pairwise.corr_heatmap <- function(input_mat_1, input_mat_2 = NULL,
                                       name_1 = NULL, name_2 = NULL,
                                       corr_type = c("pearson", "jaccard", "prop_overlap"),
                                       return_corr = FALSE,
                                       label_size = 8,
                                       color_vec = NULL){
  # Please store samples in the columns of 'input_mat'.
  if (is.null(input_mat_2)){
    input_mat_2 <- input_mat_1
  }
  if (is.null(name_2)){
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
        corr_mat[i, j] <- sum(vec_1 * vec_2) / sum(vec_1)
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
    legend_name <- "Pearson Correlation"
    if (is.null(color_vec)){
      colormap <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                       colors = c("blue", "white", "red"))
    } else {
      colormap <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                       colors = color_vec)
    }
    
  }
  if (corr_type == "jaccard" | corr_type == "prop_overlap"){
    legend_name <- ifelse(corr_type == "jaccard",
                          "Jaccard Index", "% of Shared Non-Zero Genes")
    if (is.null(color_vec)){
      colormap <- circlize::colorRamp2(breaks = c(0, 0.5, 1),
                                       colors = c("black", "purple", "gold"))
    } else {
      colormap <- circlize::colorRamp2(breaks = c(0, 0.5, 1),
                                       colors = color_vec)
    }
  }
  
  ht <- Heatmap(corr_mat,
                col = colormap,
                name = legend_name,
                row_title = name_1, column_title = name_2,
                cluster_rows = F, cluster_columns = F,
                row_names_gp = gpar(fontsize = label_size),
                column_names_gp = gpar(fontsize = label_size))
  draw(ht)
  if (return_corr){
    return(corr_mat)
  }
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
                                           name_1, name_2,
                                           fdr_cutoff = 0.05,
                                           zoom_in_start = NULL, zoom_in_end = NULL,
                                           title_text = "Factor~KO Associations"){
  df1 <- data.frame(pval = pval_vec_1,
                    fdr = p.adjust(pval_vec_1, method = "fdr"),
                    type = name_1)
  df1 <- df1 %>% arrange(pval) %>%
    mutate(rank = 1:nrow(df1),
           pass_fdr = fdr < fdr_cutoff)
  df2 <- data.frame(pval = pval_vec_2,
                    fdr = p.adjust(pval_vec_2, method = "fdr"),
                    type = name_2)
  df2 <- df2 %>% arrange(pval) %>%
    mutate(rank = 1:nrow(df2),
           pass_fdr = fdr < fdr_cutoff)
  paired_df <- rbind(df1, df2) %>%
    mutate(neg_logp = -log10(pval),
           type = factor(type, levels = c(name_1, name_2)))

  if (is.numeric(zoom_in_start) | is.numeric(zoom_in_end)){
    if (is.null(zoom_in_end)){
      paired_df <- paired_df %>% filter(rank >= zoom_in_start)
    }
    if (is.null(zoom_in_start)){
      paired_df <- paired_df %>% filter(rank <= zoom_in_end)
    }
    paired_df <- paired_df %>% filter(rank >= zoom_in_start & rank <= zoom_in_end)
  }
  p1 <- ggplot(paired_df, aes(x = rank, y = neg_logp, color = type, shape = pass_fdr)) +
    geom_point(size = 1.2, alpha = 1) +
    scale_shape_manual(values = c(3, 17)) +
    labs(title = paste(name_1, "vs", name_2, title_text),
         y = "-log10(P value)",
         color = "Method",
         shape = paste0("FDR < ", fdr_cutoff)) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.5)),
           shape = guide_legend(override.aes = list(alpha = 1, size = 1.5)))
  print(p1)
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

print_enrich_ORA_tb <- function(enrich_list,
                                fdr_cutoff = 0.05, FC_cutoff = 2,
                                print_index = NULL,
                                enrich_type = "GO terms",
                                list_type = c("per_factor", "per_marker"),
                                convert_genes = FALSE,
                                gene_map = NULL){
  if (is.null(print_index)){
    print_index <- 1:length(enrich_list)
  }
  signif_num <- rep(0, length(enrich_list))
  for (i in 1:length(enrich_list)){
    if (is.null(enrich_list[[i]])){ next }
    signif_tb <- enrich_list[[i]] %>%
      filter(enrichmentRatio >= FC_cutoff, FDR <= fdr_cutoff)
    if (nrow(signif_tb) == 0){ next }
    signif_tb <- signif_tb %>%
      mutate(GeneRatio = paste(overlap, size, sep = "/")) %>%
      mutate(GeneSet = paste0("[", geneSet, "](", link, ")")) %>%
      select(GeneSet, description, enrichmentRatio, pValue, FDR, GeneRatio, userId) %>%
      arrange(-enrichmentRatio) %>%
      dplyr::rename(enrichRatio = enrichmentRatio, geneIDs = userId)
    signif_tb$enrichRatio <- signif(signif_tb$enrichRatio, digits = 3)
    signif_tb$pValue <- format(signif_tb$pValue, digits = 3)
    signif_tb$FDR <- format(signif_tb$FDR, digits = 3)
    signif_num[i] <- nrow(signif_tb)
    
    if (i %in% print_index){
      if (convert_genes){
        stopifnot(is.data.frame(gene_map) & c("ID", "Name") %in% names(gene_map))
        signif_tb <- signif_tb %>% rowwise() %>%
          mutate(geneSymbols = convert_IDs_to_Symbols(geneIDs, gene_map)) %>%
          select(-geneIDs)
      }
      if (list_type == "per_factor"){
        caption_text <- paste("Factor", i, ":", nrow(signif_tb), "significant", enrich_type)
      }
      if (list_type == "per_marker"){
        caption_text <- paste(names(enrich_list)[i], ":", nrow(signif_tb), "significant", enrich_type)
      }
      cat(caption_text)
      cat("\n")
      print(signif_tb %>%
              kable(escape = FALSE, format = "html") %>%
              kable_styling() %>%
              column_spec(column = 2, width = "12em; display: inline-block;") %>%
              column_spec(column = 7, width = "50em; display: inline-block;") %>%
              scroll_box(width = "100%", height = '400px'))
      cat("\n")
      cat("------------")
      cat("\n")
    }
  }
  return(signif_num)
}

print_enrich_GSEA_tb <- function(enrich_list,
                                 fdr_cutoff = 0.1,
                                 print_index = NULL,  
                                 enrich_type = "GO terms",
                                 list_type = c("per_factor", "per_marker"),
                                 convert_genes = FALSE,
                                 gene_map = NULL){
  if (is.null(print_index)){
    print_index <- 1:length(enrich_list)
  }
  signif_num <- rep(0, length(enrich_list))
  for (i in 1:length(enrich_list)){
    if (is.null(enrich_list[[i]])){ next }
    signif_tb <- enrich_list[[i]] %>%
      dplyr::filter(FDR <= fdr_cutoff)
    if (nrow(signif_tb) == 0){ next }
    signif_tb <- signif_tb %>%
      mutate(GeneSet = paste0("[", geneSet, "](", link, ")")) %>%
      select(GeneSet, description, normalizedEnrichmentScore, pValue, FDR, size, leadingEdgeNum, userId) %>%
      arrange(-normalizedEnrichmentScore) %>%
      dplyr::rename(NES = normalizedEnrichmentScore, numLeadGenes = leadingEdgeNum, geneIDs = userId)
    signif_tb$NES <- signif(signif_tb$NES, digits = 3)
    signif_tb$pValue <- format(signif_tb$pValue, digits = 3)
    signif_tb$FDR <- format(signif_tb$FDR, digits = 3)
    signif_num[i] <- nrow(signif_tb)
    
    if (i %in% print_index){
      if (convert_genes){
        stopifnot(is.data.frame(gene_map) & c("ID", "Name") %in% names(gene_map))
        signif_tb <- signif_tb %>% rowwise() %>%
          mutate(geneSymbols = convert_IDs_to_Symbols(geneIDs, gene_map)) %>%
          select(-geneIDs)
      }
      if (list_type == "per_factor"){
        caption_text <- paste("Factor", i, ":", nrow(signif_tb), "significant", enrich_type)
      }
      if (list_type == "per_marker"){
        caption_text <- paste(names(enrich_list)[i], ":", nrow(signif_tb), "significant", enrich_type)
      }
      cat(caption_text)
      cat("\n")
      print(signif_tb %>%
              mutate(NES = cell_spec(NES, color = ifelse(NES > 0, "firebrick", "forestgreen"))) %>%
              kable(escape = FALSE, format = "html") %>%
              kable_styling(full_width = F) %>%
              column_spec(column = 2, width = "12em; display: inline-block;") %>%
              column_spec(column = 8, width = "50em; display: inline-block;") %>%
              scroll_box(width = "100%", height = '400px'))
      cat("\n")
      cat("------------")
      cat("\n")
    }
  }
  return(signif_num)
}

convert_IDs_to_Symbols <- function(gene_id_string, gene_map){
  gene_ids <- strsplit(gene_id_string, split = ";")[[1]]
  gene_symbols <- gene_map %>% filter(ID %in% gene_ids) %>% pull(Name)
  gene_symbol_string <- paste(gene_symbols, collapse = "; ")
  return(gene_symbol_string)
}

