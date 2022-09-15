#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install(c('cowplot', 'survminer', 'clValid', 'aricode', 'pbapply', 'ggpubr', 
 #                      'omicade4', 'SNFtool', 'Rtsne', 'ggplot2', 'ranger', 'survival',
  #                     'data.table', 'caret'))

library(data.table)
library(ggplot2)
library(omicade4)
library(SNFtool)
library(ranger)
library(Rtsne)
library(cowplot)
library(survival)
library(survminer)
library(clValid)
library(caret)
library(aricode) 
library(pbapply)
library(ggpubr)
library(cowplot)

ggplot2::theme_set(ggpubr::theme_pubclean())

#Import data
dat <- readRDS('./TCGA.pan_gastrointestinal.RDS')
colData <- dat$colData

#Description of how I got this data 
#mut/gex/meth: https://github.com/BIMSBbioinfo/uyar_et_al_multiomics_deeplearning/blob/main/src/download.tcga.R
#cnv: https://github.com/BIMSBbioinfo/uyar_et_al_multiomics_deeplearning/blob/main/src/download.tcga.firehose.R

sapply(dat, dim)
dat$colData[1:5, 1:10]
table(colData$project)
table(colData$ajcc_pathologic_tumor_stage)
table(colData$site_of_resection_or_biopsy)
hist(colData$age_at_initial_pathologic_diagnosis)

head(dat$gex[1:5, 1:10])

dim(dat$gex)
cor(dat$gex)
pheatmap::pheatmap(cor(dat$gex))
pheatmap::pheatmap(cor(dat$meth))

# t-distributed stochastic neighbor embedding
# M: matrix of samples (on rows) vs features (on columns)
compute_tsne <- function(M) {
  tsne <- Rtsne::Rtsne(M)
  df <- data.frame(tsne$Y, row.names = rownames(M), check.names = F)
  colnames(df) <- c('tSNE1', 'tSNE2')
  return(df)
}

df <- compute_tsne(t(dat$gex))

ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point()

df$project <- colData[match(rownames(df), samples)]$project
p1 <- ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = project), size = 3) + 
  labs(title = 'GEX')

df$subtype <- colData[match(rownames(df), samples)]$subtype
p2 <- ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = subtype), size = 3)
cowplot::plot_grid(p1, p2)

df.cnv <- compute_tsne(t(dat$cnv))
df.cnv$project <- colData[match(rownames(df.cnv), samples)]$project
p3 <- ggplot(df.cnv, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = project), size = 3) +
  labs(title = 'CNV')

cowplot::plot_grid(p1, p3)

df.mut <- compute_tsne(t(dat$mut))
df.mut$project <- colData[match(rownames(df.mut), samples)]$project
p4 <- ggplot(df.mut, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = project), size = 3) +
  labs(title = 'MUT')

df.meth <- compute_tsne(t(dat$meth))
df.meth$project <- colData[match(rownames(df.meth), samples)]$project
p5 <- ggplot(df.meth, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = project), size = 3) +
  labs(title = 'METH')

#Integrate layers using MCIA
integration <- omicade4::mcia(df.list = dat[c('mut', 'cnv', 'gex', 'meth')], 
                              cia.nf = 10)
latent_factors <- integration$mcoa$SynVar
rownames(latent_factors) <- gsub("\\.", "-", rownames(latent_factors))

df.int <- compute_tsne(latent_factors)
df.int$project <- colData[match(rownames(df.int), samples)]$project
ggplot(df.int, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = project), size = 3) +
  labs(title = 'Integrated')

#Exploring clinical variables with latent factors 
age <- colData[match(rownames(latent_factors), samples)]$age_at_initial_pathologic_diagnosis
sort(abs(apply(latent_factors, 2, function(x) {
  stats::cor(x, age, use = 'complete', method = 'pearson')
})), decreasing = T)[1:10]

ggpubr::ggscatter(data.frame('Var1' = latent_factors$SynVar1, 
                             'age' = age), 
                  x = 'age', y = 'Var1', add = 'reg.line', 
                  cor.coef = T)

tumor_stage <- colData[match(rownames(latent_factors), samples)]$ajcc_pathologic_tumor_stage

dt <- data.table(cbind(data.frame(latent_factors, check.names = F), 
            data.frame('stage' = tumor_stage)))

mdt <- melt.data.table(dt[!is.na(stage)], id.vars = 'stage')

ggboxplot(mdt, x = 'stage', y = 'value', color = 'stage', 
          add = 'jitter', facet.by = 'variable', nrow = 3) + 
  stat_compare_means(label.y = 4)

ggboxplot(mdt[variable == 'SynVar4'], x = 'stage', y = 'value', color = 'stage', 
          add = 'jitter', facet.by = 'variable', nrow = 3) + 
  stat_compare_means(label.y = 4)


#Clustering
#1. kmeans
ggplot(df.int, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = project), size = 3) +
  labs(title = 'Integrated')

clustering_experiment <- clValid::clValid(latent_factors, nClust = 2:5, 
                                          clMethods = 'kmeans', 
                                          validation = c('internal', 'stability'))

clValid::optimalScores(clustering_experiment)
cl <- clustering_experiment@clusterObjs$kmeans$`5`$cluster

df.int$cluster <- paste0("subtype", cl[rownames(df.int)])
ggplot(df.int, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = cluster), size = 3) +
  labs(title = 'Integrated')


# find cluster specific features
diff <- get_factor_specific_variables(latent_factors, df.int$cluster)
diff[,.SD[which.min(padj)],by = c('ref_cl')]

df <- cbind(df.int, latent_factors[,c('SynVar1'),drop = F])
p1 <- ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = SynVar1), size = 3) +
  labs(title = 'Integrated') + 
  scale_color_gradient(low = 'gray', high = 'red')
p2 <- ggplot(df, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = cluster)) +
  labs(title = 'Integrated')

cowplot::plot_grid(p1, p2)

#2. snf
# simple snf clustering function
snf <- get_SNF_clusters(lapply(dat[1:4], t))
df.int$snf_cluster <- snf[rownames(df.int), 4]
ggplot(df.int, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(color = as.factor(snf_cluster))) +
  labs(title = 'Integrated') 

aricode::AMI(df.int$project, df.int$cluster)
aricode::AMI(df.int$project, df.int$snf_cluster)
aricode::AMI(df.int$cluster, df.int$snf_cluster)

p1 <- plot_cluster_correspondence(df.int[,'project',drop=F], df.int[,'cluster',drop=F])
p2 <- plot_cluster_correspondence(df.int[,'project',drop=F], df.int[,'snf_cluster',drop=F])
cowplot::plot_grid(p1, p2)

#survival analysis
df.int$status <- colData[match(rownames(df.int), samples)]$PFI
df.int$time <- colData[match(rownames(df.int), samples)]$PFI.time
surv_object <-  survival::Surv(time = df.int$time, event = df.int$status)
fit1 <- survminer::surv_fit(surv_object ~ cluster, df.int)
survminer::ggsurvplot(fit1, df.int, pval = TRUE, risk.table = T, surv.median.line = 'hv')

fit2 <- survminer::surv_fit(surv_object ~ snf_cluster, df.int)
survminer::ggsurvplot(fit2, df.int, pval = TRUE, risk.table = T, surv.median.line = 'hv')

# classification: predicting subtype labels 
y <- colData$subtype
M <- latent_factors[colData$samples,]
res <- train_random_forest(M, as.factor(y))
pred <- stats::predict(res$fit, res$test_data)
table(res$test_data$y, pred$predictions)
caret::confusionMatrix(table(res$test_data$y, pred$predictions))

# regression analysis: predict age using latent factors
M <- latent_factors[colData[!is.na(age_at_initial_pathologic_diagnosis)]$samples,]
y <- colData[!is.na(age_at_initial_pathologic_diagnosis)]$age_at_initial_pathologic_diagnosis
res <- train_random_forest(M, y)
pred <- stats::predict(res$fit, res$test_data)
df <- data.frame('known' = res$test_data$y, 'predicted' = pred$predictions)
ggscatter(df, 
          x = 'known', y = 'predicted', add = 'reg.line', 
          cor.coef = T)
stats::cor(df$known, df$predicted)

#regression analysis: trametinib drug response prediction in cell lines 
dat <- readRDS('./CCLE.drug_response.RDS')
sapply(dat, dim)
head(dat$drugResponse)
hist(dat$drugResponse$value)

y <- dat$drugResponse[match(colnames(dat$mut), sample_id)]$value
mut_cnv <- mcia(dat[c('mut', 'cnv')], cia.nf = 10)
res <- train_random_forest(mut_cnv$mcoa$SynVar, y)
pred <- stats::predict(res$fit, res$test_data)
df <- data.frame('known' = res$test_data$y, 'predicted' = pred$predictions)
ggscatter(df, 
          x = 'known', y = 'predicted', add = 'reg.line', 
          cor.coef = T)

mut_cnv_gex <- mcia(dat[c('mut', 'cnv', 'gex')], cia.nf = 10)
res <- train_random_forest(mut_cnv_gex$mcoa$SynVar, y)
pred <- stats::predict(res$fit, res$test_data)
df <- data.frame('known' = res$test_data$y, 'predicted' = pred$predictions)
ggscatter(df, 
          x = 'known', y = 'predicted', add = 'reg.line', 
          cor.coef = T)

mut_cnv.stats <- repeat_modelling(mut_cnv$mcoa$SynVar, y, reps = 50)
mut_cnv_gex.stats <- repeat_modelling(mut_cnv_gex$mcoa$SynVar, y, reps = 50)

mut_cnv.stats$type <- 'Genome'
mut_cnv_gex.stats$type <- 'Genome+Transcriptome'

dt <- rbind(mut_cnv.stats, mut_cnv_gex.stats)
ggboxplot(dt, x = 'type', y = 'cor', add = 'jitter', color = 'type') + 
  stat_compare_means()









