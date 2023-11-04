if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("devtools")    # only if devtools not yet installed
BiocManager::install("pachterlab/sleuth")
devtools::install_github("DillonHammill/HeatmapR")

library(sleuth)
library(data.table)
library(corrplot)
library(psych)
library(ggplot2)


# path to kallisto data
base_dir <- "/Users/yasnezhok/Desktop/Kallisto/Samples"
sample_id <- dir(file.path(base_dir))
s2c <- read.table(file.path("/Users/yasnezhok/Desktop/Kallisto/info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)                 
kal_dirs <- file.path(base_dir, sample_id)
s2c <- dplyr::mutate(s2c, path = kal_dirs)


# create sleuth model 
so <- sleuth_prep(s2c, extra_bootstrap_summary = T)

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')

# model - likehood ratio test
so <- sleuth_lrt(so, 'reduced', 'full')

# show our model
models(so)

# matrix of differentially expressed TEs
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval < 0.05)
head(sleuth_significant, 22)

# significant vs non significant TEs
table(sleuth_table[,"qval"] < 0.05)

# how different is the expression of the most significant TEs
mean(sleuth_significant$test_stat[1:22] / sd(sleuth_significant$test_stat))

# some statistical tests 
shapiro.test(sleuth_significant$test_stat)
wilcox.test(sleuth_significant$test_stat, (sleuth_significant$test_stat[1:22]))
# write results
write.csv(sleuth_significant, file = "/Users/yasnezhok/Desktop/sleuth_significant.csv")
# table of results
kallisto_table(so, use_filtered=T, normalized=T, include_covariates=T)

# matrix of results
matrix = sleuth_to_matrix(so, "obs_norm", "tpm")
write.csv(matrix, file = "/Users/yasnezhok/Desktop/res_matrix.csv")

# correlation matrix
filtered.cor <- function(x){    
  num_var <- sapply(x, function(x) is.numeric(x))    
  cor_mat <- cor(x[, num_var], method = 'spearman')   
  diag(cor_mat) <- 1    
  return(cor_mat)}
data2 = data1[c(1,5,2,8,7,4,3,6)]
str(data1)
data1 = as.data.frame(matrix)
str(data1)
data2 = data1[c(3,6,7,4,2,8,1,5)]
setnames(data2, old=c( "data.ERR5881995", "data.ERR5881996", "data.ERR5881997", "data.ERR5882004", "data.ERR5882006", "data.ERR5882007", "data.SRR11802228", "data.SRR11802260", "data.SRR16145899", "data.SRR19548554", "data.SRR19548555", "data.SRR19548556", "data.SRR19548559", "data.SRR19763997", "data.SRR19763998", "data.SRR19763999", "data.SRR6781219", "data.SRR6781220", "data.SRR6781222", "data.SRR6781223", "data.SRR6781224"), new=c ("data.ERR5881995", "data.ERR5881996", "data.ERR5881997", "data.ERR5882004", "data.ERR5882006", "data.ERR5882007", "data.SRR11802228", "data.SRR11802260", "data.SRR16145899", "data.SRR19548554", "data.SRR19548555", "data.SRR19548556", "data.SRR19548559", "data.SRR19763997", "data.SRR19763998", "data.SRR19763999", "data.SRR6781219", "data.SRR6781220", "data.SRR6781222", "data.SRR6781223", "data.SRR6781224"))
data2 = data1[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
setnames(data2, old=c( "data.ERR5881995", "data.ERR5881996", "data.ERR5881997", "data.ERR5882004", "data.ERR5882006", "data.ERR5882007", "data.SRR11802228", "data.SRR11802260", "data.SRR16145899", "data.SRR19548554", "data.SRR19548555", "data.SRR19548556", "data.SRR19548559", "data.SRR19763997", "data.SRR19763998", "data.SRR19763999", "data.SRR6781219", "data.SRR6781220", "data.SRR6781222", "data.SRR6781223", "data.SRR6781224"), new=c ("data.ERR5881995", "data.ERR5881996", "data.ERR5881997", "data.ERR5882004", "data.ERR5882006", "data.ERR5882007", "data.SRR11802228", "data.SRR11802260", "data.SRR16145899", "data.SRR19548554", "data.SRR19548555", "data.SRR19548556", "data.SRR19548559", "data.SRR19763997", "data.SRR19763998", "data.SRR19763999", "data.SRR6781219", "data.SRR6781220", "data.SRR6781222", "data.SRR6781223", "data.SRR6781224"))


cor_matrix = filtered.cor(data2)
cor.plot(data2)

# p-value of correlation
r = corr.test(data1, method = "spearman")
corr.p(r$r, 21)

# vizualization of correlation matrix
corrplot.mixed(cor_matrix, tl.cex = 0.5)

write.csv(cor_matrix, file = "/Users/yasnezhok/Desktop/cor_matrix.csv")

# boxplot of expression of particular 
plot_bootstrap(so, "DF0000623_4_LTR9B", units = "est_counts", color_by = "condition", divide_groups=F ) + theme_bw() 
plot_bootstrap(so, "DF0000623_4_LTR9B", units = "est_counts", color_by = "condition")
plot1 <- plot_bootstrap(so, "DF0000623_4_LTR9B", units = "est_counts", color_by = "condition")
ggsave("DF0000623_4_LTR9B.png")

# bootstrap values used in boxplots 
penelope_bootstrap <- get_bootstrap_summary(so, "DF0000623_4_LTR9B", units = "est_counts")
write.csv(penelope_bootstrap, file = "/Users/yasnezhok/Desktop/penelope_bootstrap.csv")

# heatmap of 22 most significant TEs
transcripts = sleuth_significant$target_id[1:22]
plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log", cluster_transcripts = FALSE) # +theme(axis.text.x = element_text(size=20))
hp <- plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log")

transcripts

heatmap <- plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log", cluster_transcripts = FALSE)
heat_map_save("hp.png")
png(file="hpp.png")
heatmap("hp")
dev.off()

# pca
p <- plot_pca(so, pc_x=1L, pc_y=2L, units="est_counts",text_labels = T, color_by = "condition", use_filtered = T) + theme_bw()
  
# % of variance in PC
plot_pc_variance(so,units='est_counts') + theme_bw()

# loadings of PC
plot_loadings(so, use_filtered = TRUE, sample = NULL, pc_input = "PC1",
              units = "est_counts", pc_count=5, scale = F,
              pca_loading_abs = TRUE) + theme_bw() + theme(axis.text.x=element_text(angle=70, hjust=1))
sleuth_live(so)

#  PCA values
mat = sleuth:::spread_abundance_by(
  abund = so$obs_norm_filt,
  var = "est_counts",
  which_order = so$sample_to_covariates$sample)
pca_res <- prcomp(mat,scale. = F,center = F)


#  tpm and est_counts for a particular TE
tmp <- so$obs_raw %>% dplyr::filter(target_id == "DF0000405.4")
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp

