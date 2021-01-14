DE_test <- function(expression_mat, condition, gene_names = NULL, test.use = "t.test"){
  # expression_mat: feature by sample matrix
  # condition: 0-1 vector indicating the treatment of each sample (0: control)
  stopifnot(ncol(expression_mat) == length(condition))
  DE_res <- data.frame(matrix(nrow = nrow(expression_mat),
                              ncol = 4))
  names(DE_res) <- c("mean_1", "mean_0", "p_val", "fdr")
  for (i in 1:nrow(DE_res)){
    group_0 <- expression_mat[i, condition == 0]
    group_1 <- expression_mat[i, condition == 1]
    if (test.use == "wilcox"){
      test_res <- wilcox.test(group_1, group_0, alternative = "two.sided")
    }
    if (test.use == "t.test"){
      test_res <- t.test(group_1, group_0, alternative = "two.sided")
    }
    DE_res$mean_0[i] <- mean(group_0)
    DE_res$mean_1[i] <- mean(group_1)
    DE_res$p_val[i] <- test_res$p.value
  }
  DE_res$fdr <- p.adjust(DE_res$p_val, method = "fdr")
  if (is.null(gene_names)){
    gene_names <- rownames(expression_mat)
  }
  DE_res$gene_ID <- gene_names
  DE_res <- DE_res %>% arrange(fdr, p_val)
  return(DE_res)
}

gridplot_dge_volcano <- function(DE_list, cutoff, title_text){
  # Traditional DGE volcano plots (logFC vs -log p-value) in a grid
  plot.lst <- list()
  for (m in names(DE_list)){
    DE_df <- DE_list[[m]]
    stopifnot(c("avg_logFC", "p_val", "fdr") %in% colnames(DE_df))
    DE_df <- DE_df %>% mutate(pass_fdr = fdr < cutoff,
                              neglog_pval = -log10(p_val))
    p <- ggplot(DE_df,
                aes_string(x = "avg_logFC", y = "neglog_pval", color = "pass_fdr")) +
      geom_point(size = 1) +
      scale_color_manual(values = c("black", "red")) +
      labs(x = "Average logFC", y = "-log10(P value)",
           title = paste0(m, " (", sum(DE_df$pass_fdr), " genes)")) +
      theme(legend.position = "none")
    plot.lst[[m]] <- p
  }
  args <- c(plot.lst, list(ncol = floor(sqrt(length(DE_list))),
                           top = paste0(title_text, "\n(Red: FDR < ", cutoff, ")")))
  do.call(grid.arrange, args)
}

merge_lists <- function(DE_list_1, DE_list_2){
  merged_list <- list()
  common_names <- intersect(names(DE_list_1), names(DE_list_2))
  for (m in common_names){
    merged_list[[m]] <- inner_join(DE_list_1[[m]], DE_list_2[[m]], by = "gene_ID")
  }
  return(merged_list)
}

DGE_overlap_num <- function(lfsr_mat, DE_list,
                            lfsr_cutoff = 0.05, dge_cutoff = 0.05,
                            lfsr_name = "GSFA", dge_name = "DGE"){
  markers <- colnames(lfsr_mat)
  overlap_df <- as.data.frame(matrix(nrow = length(markers), ncol = 4))
  names(overlap_df) <- c("KO", "olap_num", "dge_num", "lfsr_num")
  overlap_df$KO <- markers
  for (i in 1:length(markers)){
    m <- markers[i]
    DE_df <- DE_list[[m]]
    stopifnot(c("fdr", "gene_ID") %in% names(DE_df))
    lfsr_df <- data.frame(gene_ID = rownames(lfsr_mat),
                          lfsr = lfsr_mat[, m])
    DE_df <- inner_join(DE_df, lfsr_df, by = "gene_ID")
    overlap_df$olap_num[i] <- DE_df %>% filter(fdr < dge_cutoff, lfsr < lfsr_cutoff) %>% nrow()
    overlap_df$dge_num[i] <- DE_df %>% filter(fdr < dge_cutoff) %>% nrow()
    overlap_df$lfsr_num[i] <- DE_df %>% filter(lfsr < lfsr_cutoff) %>% nrow()
  }
  print(knitr::kable(overlap_df,
                     caption = paste0(dge_name, " FDR cutoff: ", dge_cutoff, ";\n",
                                      lfsr_name, " LFSR cutoff: ", lfsr_cutoff)) %>%
          kable_styling() %>% scroll_box(width = '100%'))
  return(overlap_df)
}

compute_beta_dot_W <- function(lfsr_mat, gibbs_PM){
  DE_list <- list()
  for (m in 1:length(colnames(lfsr_mat))){
    DE_df <- data.frame(gene_ID = rownames(lfsr_mat),
                        lfsr = lfsr_mat[, m],
                        beta_W = NA)
    for (i in 1:nrow(DE_df)){
      gene_id <- DE_df$gene_ID[i]
      DE_df$beta_W[i] <- sum(gibbs_PM$W_pm[i, ] * gibbs_PM$beta_pm[m, ])
      DE_df <- DE_df %>% mutate(pass_lfsr = lfsr < 0.05)
    }
    marker <- colnames(lfsr_mat)[m]
    DE_list[[marker]] <- DE_df
  }
  return(DE_list)
}

gridplot_betaW_lfsr <- function(DE_list, lfsr_cutoff = 0.05, title_text = "GSFA"){
  # Plot GSFA estimated effect size (beta dotprod W) vs LFSR
  # GSFA version of "volcano plot"
  plot.lst <- list()
  for (m in names(DE_list)){
    DE_df <- DE_list[[m]]
    DE_df <- DE_df %>% mutate(neglog_lfsr = -log10(lfsr + 1e-3)) %>%
      mutate(neglog_pval = -log10(p_val)) %>%
      mutate(neglog_pval = ifelse(neglog_pval > 10, 10, neglog_pval))
    p <- ggplot(DE_df, aes_string(x = "beta_W", y = "neglog_lfsr", color = "neglog_pval")) +
      geom_point(size = 0.6) +
      # scale_y_continuous(limits = c(-0.3, 0.3)) +
      scale_color_gradientn(colors = c("lightblue", "yellow", "red"),
                            values = scales::rescale(c(0, 3, 10)),
                            limits = c(0, 10)) +
      labs(title = paste0(m, " (", sum(DE_df$lfsr < lfsr_cutoff), " genes)"),
           x = "Estimated Beta * W",
           y = "-log10(LFSR)",
           color = "-log10(p-value)") +
      theme(axis.text = element_text(size = 11),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10))
    if (min(DE_df$lfsr) < 0.05){
      p <- p + geom_hline(yintercept = -log10(lfsr_cutoff),
                          color = "red", linetype = "dashed")
    }
    plot.lst[[m]] <- p
  }
  args <- c(plot.lst, list(ncol = floor(sqrt(length(DE_list))),
                           top = paste0(title_text, "\n(colored by p-value from t-test)")))
  do.call(grid.arrange, args)
}
