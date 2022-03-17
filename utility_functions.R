compute_tsne <- function(M) {
  tsne <- Rtsne::Rtsne(M)
  df <- data.frame(tsne$Y, row.names = rownames(M), check.names = F)
  colnames(df) <- c('tSNE1', 'tSNE2')
  return(df)
}


# M: matrix of samples (on rows) vs features (on columns)
# factors: annotation vector for the samples 
plot_pca <- function(M, factors) {
  pca <- stats::prcomp(M, rank. = 10) # go checkout https://www.youtube.com/watch?v=HMOI_lkzW08
  df <- data.frame(pca$x, check.names = FALSE, row.names = rownames(pca$x))
  df$group <- factors
  
  var_exp <- round(diag(cov(pca$x))/sum(diag(cov(pca$x))) * 100, 1)
  ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = group), size = 5, alpha = 0.5) +
    theme_bw(base_size = 8) +
    labs(x = paste0('PC1 (',var_exp[['PC1']],'%)'),
         y = paste0('PC2 (',var_exp[['PC2']],'%)')) 
}

get_factor_specific_variables <- function(M, factors)  {
  res <- do.call(rbind, pbapply::pblapply(colnames(M), function(x) {
    l <- split(M[,x], factors)
    do.call(rbind, lapply(names(l), function(i) {
      g1 <- l[[i]]
      g2 <- as.vector(unlist(l[setdiff(names(l), i)]))
      t <- wilcox.test(g1, g2, 
                       alternative = 'greater')
      data.frame("variable" = x, "ref_cl" = i, "pval" = t$p.value)
    }))
  }))
  res$padj <- p.adjust(res$pval, method = 'BH')
  return(data.table::as.data.table(res))
}

get_SNF_clusters <- function(matrix.list) {
  require(SNFtool)
  common <- Reduce(intersect, lapply(matrix.list, rownames))
  message(date(), ' => getting distance matrix')
  W <- lapply(matrix.list, function(x) {
    d <- SNFtool::dist2(as.matrix(x[common,]),as.matrix(x[common,]))^(1/2)
    am <- affinityMatrix(d, K = 20, sigma = 0.5)
  })
  FN <- SNFtool::SNF(W)
  # get best clustering k
  best <- SNFtool::estimateNumberOfClustersGivenGraph(FN, NUMC = 3:7)
  message(date(), ' => clustering')
  res <- do.call(cbind, lapply(names(best), function(x) {
    group = spectralClustering(FN,best[[x]])    # the final subtypes information
    df <- data.frame('snf' = group, row.names = common)
    colnames(df)[1] <- x
    return(df)
  }))
  return(res)
}

plot_cluster_correspondence <- function(df1, df2) {
  df1 <- data.frame(df1, check.names = F)
  df2 <- data.frame(df2, check.names = F)
  if(sum(is.na(df1[,1])) > 0) {
    warning("Converting NA values to 'Undefined'")
    df1[is.na(df1[,1]),1] <- 'Undefined'
  }
  if(sum(is.na(df2[,1])) > 0) {
    warning("Converting NA values to 'Undefined'")
    df2[is.na(df2[,1]),1] <- 'Undefined'
  }
  df1[,1] <- as.factor(df1[,1])
  df2[,1] <- as.factor(df2[,1])
  
  dt <- data.table(merge(df1[,1,drop=F], df2[,1,drop = F], by = 'row.names'))
  labels <- colnames(dt[,2:3])
  colnames(dt) <- c('rn', 'g1', 'g2')
  dt1 <- dt[order(g1), c('rn', 'g1')]
  dt2 <- dt[order(g2), c('rn', 'g2')]
  dt1$r1 <- 1:nrow(dt1)
  dt2$r2 <- 1:nrow(dt2)
  dt <- merge(dt1, dt2, by = 'rn')[order(rn)]
  ami <- aricode::AMI(dt$g1, dt$g2)
  label_pos_left <- -0.05
  label_pos_right <- 1.05
  # plot segments 
  ggplot(dt[order(g1)]) + 
    geom_point(aes(x = 0, y = r1, color = g1))  +
    geom_point(aes(x = 1, y = r2, color = g2)) + 
    geom_segment(aes(x = 0, xend = 1, y = r1, yend = r2, color = g1), alpha = 0.25) +
    geom_label(data = dt[,median(r1), by = g1], aes(x = label_pos_left, y = V1, label = g1, color = g1), hjust = 1) +
    geom_segment(data = dt[,list('max' = max(r1), 'min' = min(r1), 'median' = median(r1)), by = g1], 
                 aes(x = label_pos_left, y = median, xend = 0, yend = max, color = g1)) + 
    geom_segment(data = dt[,list('max' = max(r1), 'min' = min(r1), 'median' = median(r1)), by = g1], 
                 aes(x = label_pos_left, y = median, xend = 0, yend = min, color = g1)) + 
    geom_label(data = dt[,median(r2), by = g2], aes(x = label_pos_right, y = V1, label = g2, color = g2), hjust = 0) +
    geom_segment(data = dt[,list('max' = max(r2), 'min' = min(r2), 'median' = median(r2)), by = g2], 
                 aes(x = label_pos_right, y = median, xend = 1, yend = max, color = g2)) + 
    geom_segment(data = dt[,list('max' = max(r2), 'min' = min(r2), 'median' = median(r2)), by = g2], 
                 aes(x = label_pos_right, y = median, xend = 1, yend = min, color = g2)) + 
    annotate('label', x = 0, y = max(dt$r1)+1, label = labels[1], vjust = 0) + 
    annotate('label', x = 1, y = max(dt$r2)+1, label = labels[2], vjust = 0) + 
    ggtitle(label = paste0(labels, collapse = " <-> "), 
            subtitle = paste0("Adjusted Mutual Information: ",round(ami,2))) + 
    theme_minimal() + 
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), 
          legend.position = 'none',
          plot.margin = unit(c(0.1, 1.5, 0, 1.5), units = 'in'),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5)) + 
    coord_cartesian(clip = 'off') 
}

# M: matrix of samples (in the rows) vs features (in the columns)
# y: sample labels 
train_random_forest <- function(M, y) {
  require(ranger)
  df <- data.frame(M)
  df$y <- y
  train_samples <- sample(rownames(M), round(0.7 * nrow(M)))
  test_samples <- setdiff(rownames(M), train_samples)
  training <- df[train_samples,]
  testing <- df[test_samples,]
  fit <- ranger::ranger(seed = 42, formula = y ~ ., data = training, 
                        importance = 'permutation', 
                        num.trees = 500, num.threads = 2)
  return(list('fit' = fit, 'test_data' = testing))
}

repeat_modelling <- function(M, y, reps = 10) {
  do.call(rbind, pbapply::pblapply(1:reps, function(i) {
    res <- train_random_forest(M, y)
    pred <- stats::predict(res$fit, res$test_data)
    data.table('cor' = cor(res$test_data$y, pred$predictions), 
               'R2' = caret::R2(res$test_data$y, pred$predictions),
               'i' = i)
  }))
} 

