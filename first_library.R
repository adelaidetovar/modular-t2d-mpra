##############################
### First library analysis ###
##############################

librarian::shelf(tidyverse, factoextra, data.table, MPRAnalyze,
                 reshape2, ggtext, UpSetR, patchwork, bedtoolsr,
                 glmnet, GGally, mvtnorm)

# mk work dir
work_dir <- "/path/to/workdir/"

if (!dir.exists(work_dir)){
  dir.create(work_dir)
} else {
  print("Dir already exists!")
}

# set fig dir
fig_dir <- paste0(work_dir, "figs/")

# data dir
data_dir <- paste0(work_dir, "data/")

# out dir
out_dir <- paste0(work_dir, "output/")


# functions
rnaSum <- function(rna_df){
  rna_df %>% mutate(sum_rep1 = rowSums(select(., ends_with("_rep1"))),
                    sum_rep2 = rowSums(select(., ends_with("_rep2"))),
                    sum_rep3 = rowSums(select(., ends_with("_rep3"))))
}

rnaLimit <- function(rnaSumA){
  rbind(rnaSumA[,c(1,(length(rnaSumA[1,])-2):length(rnaSumA[1,]))])
}

dnaSum <- function(dna_df){
  dna_df %>% mutate(sum_d = rowSums(across(2:length(dna_df[1,]))))
}

dnaLimit <- function(dnaSumA){
  rbind(dnaSumA[,c(1,length(dnaSumA[1,]))])
}

mpraNorm <- function(df){
  df %>%
    mutate(dna_norm = sum_d/(sum(sum_d))*1000000) %>%
    mutate_at(., vars(starts_with("sum")),
              funs(norm = ./(sum(.))*1000000)) %>%
    mutate_at(., vars(7:9),
              funs(dna_norm = . / dna_norm))
}

getSigQuant <- function(df) {
  df %>% mutate(fdr = p.adjust(df$pval.zscore, method = "BH"),
                sig = case_when(fdr <= 0.05 ~ "fdr_05",
                                fdr >= 0.05 & df$pval.zscore <= 0.05 ~ "nominal_sig",
                                df$pval.zscore > 0.05 ~ "not_sig")) %>%
    filter(fdr != 0)
}

getSigInserts <- function(df) {
  rownames(df[df$sig=="fdr_05",])
}

#############################
### Make annotation files ###
#############################

# read in Yoshi's annotations
oligo_annot <- read.delim(paste0(data_dir, 'kyono_library_annotations.txt'))
oligo_annot$basename <- gsub("_ref|_alt", "", oligo_annot$name)
oligo_annot$allele_pair <- paste0(oligo_annot$proxyalt, oligo_annot$proxyref)

# read in regulatory data
promoter = read.table(paste0(data_dir, "Islets.all_promoter.bed"))
enhancer = read.table(paste0(data_dir, "Islets.all_enhancer.bed"))
atac = read.table(paste0(data_dir, "Islets.atac_peaks.bed"))
tss = read.table(paste0(data_dir, "gencode_pc_5kb.bed"))

# get regulatory annotations
reg_annot <- oligo_annot
prom_annot <- bt.intersect(reg_annot[,c(6:8)],promoter,wao=TRUE)
enh_annot <- bt.intersect(reg_annot[,c(6:8)],enhancer,wao=TRUE)
atac_annot <- bt.intersect(reg_annot[,c(6:8)],atac,wao=TRUE)
tss_annot <- bt.intersect(reg_annot[,c(6:8)],tss,wao=TRUE)

reg_annot$prom <- prom_annot$V7
reg_annot$enh <- enh_annot$V7
reg_annot$atac <- atac_annot$V7
reg_annot$tss <- tss_annot$V7

reg_annot <- reg_annot %>%
  mutate(state = case_when(prom == enh ~ "Other",
                          prom > enh ~ "Promoter",
                          prom < enh ~ "Enhancer"),
        tss_rel_pos = ifelse(tss > 0, "Proximal", "Distal"))

write.table(reg_annot, paste0(out_dir, "kyono_annot_plus_regulatoryinfo.txt"), sep = "\t",
            row.names = F, quote = F)

###########################
### Variant-trait table ###
###########################

var_df <- oligo_annot %>%
  select(proxy_rsid, index_rsid, trait) %>%
  distinct(proxy_rsid, .keep_all = TRUE)

var_tab <- as.data.frame(table(var_df$trait))
colnames(var_tab) <- c("Trait", "Count")
var_tab$Trait <- as.character(var_tab$Trait)
var_tab$Trait[3:11] <- c("T2D", "T2D,2hrGluadjBMI,FGluadjBMI","T2D,2hrGluadjBMI,FGluadjBMI,FInsadjBMI",
                         "T2D,2hrGluadjBMI,FGluadjBMI,FInsadjBMI,HbA1c","T2D,2hrGluadjBMI,FInsadjBMI",
                         "T2D,FGluadjBMI", "T2D,FGluadjBMI,HbA1c","T2D,FInsadjBMI","T2D,HbA1c")

total_var <- sum(var_tab$Count)
var_tab$Prop <- var_tab$Count/total_var
var_tab$Type <- "Variant"

var_tab %>% ggplot(aes(x = Type, y = Count, fill = Trait, label = paste(Trait, Count, sep = "\n"))) +
  geom_col() + scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_linedraw(base_size = 12) +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none") +
  geom_text_repel(size = 4,
                  direction = "x",
                  hjust = 0)

###################################################
### Initial data load - from STARRseq pipeline  ###
###################################################

# RNA counts from MPRA
rna_counts <- list()
rna_counts[[1]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_UP.minDNA10.minBarcode2.rna_counts.tsv'))
rna_counts[[2]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_dwn.minDNA10.minBarcode2.rna_counts.tsv'))
rna_counts[[3]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_UP.minDNA10.minBarcode2.rna_counts.tsv'))
rna_counts[[4]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_dwn.minDNA10.minBarcode2.rna_counts.tsv'))
names(rna_counts) <- c("UpINS","DwnINS","UpSCP1","DwnSCP1")

# DNA counts from MPRA
dna_counts <- list()
dna_counts[[1]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_UP.minDNA10.minBarcode2.dna_counts.tsv'))
dna_counts[[2]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_dwn.minDNA10.minBarcode2.dna_counts.tsv'))
dna_counts[[3]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_UP.minDNA10.minBarcode2.dna_counts.tsv'))
dna_counts[[4]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_dwn.minDNA10.minBarcode2.dna_counts.tsv'))
names(dna_counts) <- c("UpINS","DwnINS","UpSCP1","DwnSCP1")

# Collect common inserts across all configurations
# sum RNA counts for all barcodes used for a specific insert
sum_rna <- list()
sum_dna <- list()
for (i in 1:4){
  sum_rna[[i]] <- rnaSum(rna_counts[[i]])
  sum_rna[[i]] <- rnaLimit(sum_rna[[i]])
  sum_dna[[i]] <- dnaSum(dna_counts[[i]])
  sum_dna[[i]] <- dnaLimit(sum_dna[[i]])
}

merge_counts <- list()
for (i in 1:4){
  merge_counts[[i]] <- merge(sum_rna[[i]],sum_dna[[i]],by="refname")
}

merge_counts[[1]] <- mpraNorm(merge_counts[[1]])
merge_counts[[1]] <- merge_counts[[1]][,c(1,11:13)]
colnames(merge_counts[[1]]) <- c("refname","UpINS_rep1",
                                 "UpINS_rep2","UpINS_rep3")

merge_counts[[2]] <- mpraNorm(merge_counts[[2]])
merge_counts[[2]] <- merge_counts[[2]][,c(1,11:13)]
colnames(merge_counts[[2]]) <- c("refname","DwnINS_rep1",
                                 "DwnINS_rep2","DwnINS_rep3")

merge_counts[[3]] <- mpraNorm(merge_counts[[3]])
merge_counts[[3]] <- merge_counts[[3]][,c(1,11:13)]
colnames(merge_counts[[3]]) <- c("refname","UpSCP1_rep1",
                                 "UpSCP1_rep2","UpSCP1_rep3")

merge_counts[[4]] <- mpraNorm(merge_counts[[4]])
merge_counts[[4]] <- merge_counts[[4]][,c(1,11:13)]
colnames(merge_counts[[4]]) <- c("refname","DwnSCP1_rep1",
                                 "DwnSCP1_rep2","DwnSCP1_rep3")

all_counts <- merge(merge(merge(merge_counts[[1]], merge_counts[[2]], by="refname"),
                         merge_counts[[3]], by="refname"),
                   merge_counts[[4]], by="refname")

rownames(all_counts) <- all_counts[,c(1)]
all_counts <- all_counts[,c(2:13)]

# save list of all common inserts across four configs
common_inserts <- rownames(all_counts)

#####################################
### Principal components analysis ###
#####################################

pca_counts <- as.data.frame(t(all_counts))
pca_counts$config <- c("INS-Up","INS-Up","INS-Up",
                       "INS-Down","INS-Down","INS-Down",
                       "SCP1-Up","SCP1-Up","SCP1-Up",
                       "SCP1-Down","SCP1-Down","SCP1-Down")

# generate PCA for all counts
res_pca <- prcomp(pca_counts[,c(1:11656)])
pca_loadings <- as.data.frame(res_pca$rotation)
top_pc1 <- pca_loadings %>%
  arrange(-abs(PC1)) %>%
  slice(1:50)

top_pc2 <- pca_loadings %>%
  arrange(-abs(PC2)) %>%
  slice(1:50)

pca_data = as.data.frame(res_pca$x)
pca_data$indiv <- rownames(pca_data)

# add metadata info
promoter = c("INS","INS","INS","INS","INS","INS",
             "SCP1","SCP1","SCP1","SCP1","SCP1","SCP1")
config = c("Upstream","Upstream","Upstream",
          "Downstream","Downstream","Downstream",
          "Upstream","Upstream","Upstream",
          "Downstream","Downstream","Downstream")
indiv = rownames(pca_data)
col_data <- data.frame(config, promoter, indiv)

# add plot labels to PCA df
pca_data$full_name <- paste(pca_data$config, pca_data$promoter, sep = " + ")
ix_label <- c(1, 4, 7, 10)
pca_data$label <- ""
pca_data$label[ix_label] <- pca_data$full_name[ix_label]

# plot PCA
pca_data=merge(pca_data,col_data)
percent_var = round(100*(res_pca$sdev^2/sum(res_pca$sdev^2)),digits=2)
pca_plot <- ggplot(data=pca_data,aes(x = PC1,y=PC2,fill=interaction(config,promoter),
                                     shape=config,label=label)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,24),
                     name="Insert Position") +
  scale_fill_manual(values=c("#0072B2", "#56B4E9", "#E69F00", "#F0E442"),
                    name="Promoter") +
  geom_text_repel(family = "Helvetica", min.segment.length = 0,
                  box.padding = 1) +
  theme_linedraw(base_size=12,
                 base_family = "Helvetica") +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position = "none") +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  guides(fill = guide_legend(override.aes=list(shape=22,color=NA)))

ggsave(plot=pca_plot,filename=paste0(fig_dir,'pca_inserts_minDNA10.minBarcode2.png'),
        unit="in",width=4.5,height=4.5,dpi=600)

# plot PCA correlogram with metadata
####
# check correlation btwn pcs and counts
pca_df <- as.data.frame(res_pca$x)
all_counts <- as.data.frame(t(pca_counts))

pca_df$rna_prop <- 200000/pca_df$rna_count

pca_df$rna_count <- c(sum(sum_rna[[1]]$sum_rep1), sum(sum_rna[[1]]$sum_rep2), sum(sum_rna[[1]]$sum_rep3),
                      sum(sum_rna[[2]]$sum_rep1), sum(sum_rna[[2]]$sum_rep2), sum(sum_rna[[2]]$sum_rep3),
                      sum(sum_rna[[3]]$sum_rep1), sum(sum_rna[[3]]$sum_rep2), sum(sum_rna[[3]]$sum_rep3),
                      sum(sum_rna[[4]]$sum_rep1), sum(sum_rna[[4]]$sum_rep2), sum(sum_rna[[4]]$sum_rep3))
pca_df$dna_count <- c(rep(sum(sum_dna[[1]]$sum_d),3),rep(sum(sum_dna[[2]]$sum_d),3),
                      rep(sum(sum_dna[[3]]$sum_d),3),rep(sum(sum_dna[[4]]$sum_d),3),)
pca_df$count_ratio <- log2(pca_df$rna_count/pca_df$dna_count)

# Matrix of plots
p1 <- ggpairs(pca_df[,c(1:5,13:15)], lower = list(continuous = wrap("points", alpha = 0.5)),
              upper = list(continuous = wrap(ggally_cor, method = "spearman"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
# Correlation matrix plot
p2 <- ggcorr(pca_df[,c(1:5,13:15)], label = TRUE, label_round = 2, method = c("everything","spearman"))

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill
colors <- c(colors[1:24],colors[27],colors[25:27]) # one color missing

p <- 8
# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

png(filename = paste0(fig_dir, "pca_correlogram_spearman.png"),
    units = "in", height = 6, width = 6, res = 600)
p1
dev.off()

# regress out PC1/PC2 and replot PCA

pheno <- data.frame(
  config = pca_counts$config,
  PC1 = pca_df$PC1,
  PC2 = pca_df$PC2
)

clean_pca_counts <- as.data.frame(t(pca_counts[,c(1:11656)]))

mod <- with(pheno, model.matrix(~ PC1))
Hat <- solve(t(mod) %*% mod) %*% t(mod)
ty <- t(clean_pca_counts)
ty[is.na(ty)] <- 0
beta <- (Hat %*% ty)
cleany_pc1 <- clean_pca_counts - t(as.matrix(mod) %*% beta)

mod <- with(pheno, model.matrix(~ PC2))
Hat <- solve(t(mod) %*% mod) %*% t(mod)
ty <- t(clean_pca_counts)
ty[is.na(ty)] <- 0
beta <- (Hat %*% ty)
cleany_pc2 <- clean_pca_counts - t(as.matrix(mod) %*% beta)

tcleany_pc1 <- as.data.frame(t(cleany_pc1))
tcleany_pc2 <- as.data.frame(t(cleany_pc2))
tcleany_pc1$config <- tcleany_pc2$config <- pheno$config


clean_pc1_res <- prcomp(tcleany_pc1[,c(1:11656)])
clean_pc2_res <- prcomp(tcleany_pc2[,c(1:11656)])

autoplot(clean_pc1_res, data = tcleany_pc1, color = "config")
autoplot(clean_pc2_res, data = tcleany_pc2, color = "config")

pc1_data <- as.data.frame(clean_pc1_res$x)
pc1_data$indiv <- rownames(pc1_data)

# add metadata info
promoter = c("INS","INS","INS","INS","INS","INS",
             "SCP1","SCP1","SCP1","SCP1","SCP1","SCP1")
config = c("Upstream","Upstream","Upstream",
           "Downstream","Downstream","Downstream",
           "Upstream","Upstream","Upstream",
           "Downstream","Downstream","Downstream")
indiv = rownames(pc1_data)
col_data <- data.frame(config, promoter, indiv)

# add plot labels to PCA df
pc1_data$full_name <- paste(col_data$config, col_data$promoter, sep = " + ")
ix_label <- c(1, 4, 7, 10)
pc1_data$label <- ""
pc1_data$label[ix_label] <- pc1_data$full_name[ix_label]

# plot PCA
pc1_data=merge(pc1_data,col_data)
percent_var = round(100*(clean_pc1_res$sdev^2/sum(clean_pc1_res$sdev^2)),digits=2)
pca_plot <- ggplot(data=pc1_data,aes(x = PC1,y=PC2,fill=interaction(config,promoter),
                                     shape=config,label=label)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,24),
                     name="Insert Position") +
  scale_fill_manual(values=c("#0072B2", "#56B4E9", "#E69F00", "#F0E442"),
                    name="Promoter") +
  geom_text_repel(family = "Helvetica", min.segment.length = 0,
                  box.padding = 1) +
  theme_linedraw(base_size=12,
                 base_family = "Helvetica") +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position = "none") +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  guides(fill = guide_legend(override.aes=list(shape=22,color=NA)))

ggsave(plot=pca_plot,filename=paste0(fig_dir, "pca_plot_PC1.png"),
       unit="in",width=4.5,height=4.5,dpi=600)

pc2_data <- as.data.frame(clean_pc2_res$x)
pc2_data$indiv <- rownames(pc2_data)

# add metadata info
promoter = c("INS","INS","INS","INS","INS","INS",
             "SCP1","SCP1","SCP1","SCP1","SCP1","SCP1")
config = c("Upstream","Upstream","Upstream",
           "Downstream","Downstream","Downstream",
           "Upstream","Upstream","Upstream",
           "Downstream","Downstream","Downstream")
indiv = rownames(pc2_data)
col_data <- data.frame(config, promoter, indiv)

# add plot labels to PCA df
pc2_data$full_name <- paste(col_data$config, col_data$promoter, sep = " + ")
ix_label <- c(1, 4, 7, 10)
pc2_data$label <- ""
pc2_data$label[ix_label] <- pc2_data$full_name[ix_label]

# plot PCA
pc2_data=merge(pc2_data,col_data)
percent_var = round(100*(clean_pc2_res$sdev^2/sum(clean_pc2_res$sdev^2)),digits=2)
pca_plot <- ggplot(data=pc2_data,aes(x = PC1,y=PC2,fill=interaction(config,promoter),
                                     shape=config,label=label)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,24),
                     name="Insert Position") +
  scale_fill_manual(values=c("#0072B2", "#56B4E9", "#E69F00", "#F0E442"),
                    name="Promoter") +
  geom_text_repel(family = "Helvetica", min.segment.length = 0,
                  box.padding = 1) +
  theme_linedraw(base_size=12,
                 base_family = "Helvetica") +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position = "none") +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  guides(fill = guide_legend(override.aes=list(shape=22,color=NA)))

ggsave(plot=pca_plot,filename=paste0(fig_dir, "pca_plot_PC2.png"),
       unit="in",width=4.5,height=4.5,dpi=600)

############
### ECDF ###
############

indiv_res <- list()
indiv_res[[1]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_UP.minDNA10.minBarcode2.type_norm.mpra.tsv'))
indiv_res[[2]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_dwn.minDNA10.minBarcode2.type_norm.mpra.tsv'))
indiv_res[[3]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_UP.minDNA10.minBarcode2.type_norm.mpra.tsv'))
indiv_res[[4]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_dwn.minDNA10.minBarcode2.type_norm.mpra.tsv'))
names(indiv_res) <- c("UpINS","DwnINS","UpSCP1","DwnSCP1")
indiv_res <- lapply(indiv_res, function(x) x %>%
         rownames_to_column(var = "refname") %>%
         filter(refname %in% common_inserts))
for (i in 1:4){
colnames(indiv_res[[i]])[2:6] <- paste(colnames(indiv_res[[i]])[2:6], names(indiv_res)[i], sep = "_")
}
indiv_df <- merge(merge(merge(indiv_res[[1]],indiv_res[[2]],by="refname"),
                        indiv_res[[3]], by="refname"),
                  indiv_res[[4]], by = "refname")

ec_df <- indiv_df[,c(2,7,12,17)]
ec_df <- melt(ec_df)

ecdf_plot <- ggplot(ec_df, aes(x = log2(value), color = variable)) +
  stat_ecdf(size = 0.5) + theme_bw(base_size = 12) +
  scale_color_manual(values = c("#F0E442", "#E69F00","#56B4E9", "#0072B2"),
                     name = "Configuration",
                     labels = c("*SCP1* + Up","*SCP1* + Down","*INS* + Up","*INS* + Down")) +
  labs(y = "Proportion", x = expression(log[2]~"(RNA/DNA)")) +
  theme(legend.text = element_markdown(),
        text = element_text(family = "Helvetica"))
  

ggsave(plot = ecdf_plot, paste0(fig_dir, "ecdf.png"), units = "in", height = 5, width = 6)

###############################
### Calculating proportions ###
###############################

length(which(indiv_res[[1]]$pval.zscore<0.05))/length(indiv_res[[1]]$pval.zscore)
length(which(indiv_res[[2]]$pval.zscore<0.05))/length(indiv_res[[2]]$pval.zscore)
length(which(indiv_res[[3]]$pval.zscore<0.05))/length(indiv_res[[3]]$pval.zscore)
length(which(indiv_res[[4]]$pval.zscore<0.05))/length(indiv_res[[4]]$pval.zscore)

###################
### UpSet plots ###
###################

indiv_sig <- list()
for (i in 1:4){
  indiv_sig[[i]] <- getSigQuant(indiv_res[[i]])
}

list_sig <- list(UpINS = getSigInserts(indiv_sig[[1]]),
                 DownINS = getSigInserts(indiv_sig[[2]]),
                 UpSCP1 = getSigInserts(indiv_sig[[3]]),
                 DownSCP1 = getSigInserts(indiv_sig[[4]]))

all_sig <- Reduce(union, list(list_sig[[1]],
                              list_sig[[2]],
                              list_sig[[3]],
                              list_sig[[4]]))

sig_one <- all_sig %in% list_sig[[1]]
sig_two <- all_sig %in% list_sig[[2]]
sig_three <- all_sig %in% list_sig[[3]]
sig_four <- all_sig %in% list_sig[[4]]

sig_df = as.data.frame(cbind(sig_four, sig_three, sig_two, sig_one))
rownames(sig_df) = all_sig
colnames(sig_df) = c("A","B","C","D") # corresponds top DownSCP1, UpSCP1, DownINS, UpINS, in that order

sig_df_upset <- upset(sig_df,
                      intersect = c("A", "B",
                                    "C", "D"),
                      base_annotations = list(
                        'Intersection Size' = intersection_size(
                          counts = FALSE,
                          text = list(size = 5)
                        )
                      ),
                      queries = list(
                        upset_query(set = "A", fill = "#0072B2"),
                        upset_query(set = "B", fill = "#56B4E9"),
                        upset_query(set = "C", fill = "#E69F00"),
                        upset_query(set = "D", fill = "#F0E442")
                      ),
                      themes = upset_default_themes(
                        text = element_text(family = "Helvetica", size = 15)
                      ),
                      sort_sets = FALSE)

intersect_df <- as.data.frame(table(factor(sig_df_upset[[2]]$data$intersection)))

intersect_totals <- intersect_df %>%
  group_by(Var1) %>%
  summarize(total = sum(Freq)) %>%
  mutate(percent = round(((total*100)/11656), digits = 2),
         label = paste0(total, "\n", percent, "%"))

intersect_df <- intersect_df %>%
  left_join(intersect_totals, by = "Var1")
intersect_bar <- intersect_df %>%
  mutate(condition = case_when(Var1 == 2 ~ "UpSCP",
                               Var1 == 1 ~ "DownSCP",
                               Var1 == 4 ~ "UpINS",
                               Var1 == 3 ~ "DownINS",
                               TRUE ~ "other"),
         condition = factor(condition, levels = c("UpSCP", "DownSCP", "UpINS", "DownINS", "other")),
         Var1 = reorder(Var1, -Freq)) %>%
  ggplot(aes(x = Var1, y = Freq, fill = condition)) +
  geom_bar(stat = "identity", linewidth = 0.25, color = "black") +
  geom_text(aes(label = paste0(total, "\n(", percent, "%", ")"), 
                y = total + 200),
            size = 4, family = "Helvetica",
            angle = 45, lineheight = 0.75) +
  scale_fill_manual(values = c("#F0E442",
                               "#E69F00",
                               "#56B4E9",
                               "#0072B2",
                               "#8E8E8E")) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0,0,0,0,"pt")) +
  labs(x = NULL, y = "Intersection Size") + coord_cartesian(clip = "off")

intersection_order <- levels(reorder(intersect_df$Var1, -intersect_df$Freq))
replace_intersect <- c("3" = "DownINS", "4" = "UpINS",
                       "1" = "DownSCP", "2" = "UpSCP")
intersection_order <- as.factor(gsub("3", "DownINS",
                                     gsub("4", "UpINS",
                                          gsub("1", "DownSCP",
                                               gsub("2", "UpSCP", intersection_order)))))
intersection_order <- factor(intersection_order, levels = intersection_order)

sig_mat <- sig_df_upset[[4]]$data %>%
  mutate(start = str_sub(as.character(intersection), 1, 1),
         end = str_sub(as.character(intersection), nchar(as.character(intersection)), nchar(as.character(intersection))),
         intersection_recode = str_replace_all(intersection, replace_intersect),
         condition = case_when(!(intersection_recode %in% c("UpINS", "DownINS","UpSCP", "DownSCP")) ~ "other",
                               TRUE ~ intersection_recode),
         condition = factor(condition, levels = c("DownSCP", "UpSCP", "DownINS", "UpINS")),
         group = str_replace_all(group, replace_intersect),
         group = factor(group, levels = c("DownSCP", "UpSCP", "DownINS", "UpINS")),
         start = str_replace_all(start, replace_intersect),
         end = str_replace_all(end, replace_intersect),
         intersection_recode = factor(intersection_recode, levels = intersection_order)) %>%
  ggplot(aes(x = intersection_recode, y = group, fill = condition,
             size = value)) +
  geom_point(fill = "#E4E4E2", color = "#C0C0C0", size = 4, pch = 21, stroke = 0.5) +
  geom_segment(aes(x = intersection_recode, xend = intersection_recode, y = start, yend = end),
               linewidth = 0.75) +
  geom_point(aes(stroke = value * 0.75), pch = 21) + scale_fill_manual(values = c("#E69F00",
                                                                                  "#F0E442",
                                                                                  "#0072B2",
                                                                                  "#56B4E9",
                                                                                  "#000000")) +
  scale_size_manual(values = c(0, 5)) + theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black", size = 14),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  labs(x = "Intersection", y = NULL) + scale_y_discrete(labels = c("Downstream + SCP1",
                                                                   "Upstream + SCP1",
                                                                   "Downstream + INS",
                                                                   "Upstream + INS"))


design <- "
11
22
"

active_upset_plot <- free(intersect_bar, type = "space", side = "l") + sig_mat + plot_layout(design = design)
ggsave(active_upset_plot, filename = paste0(fig_dir, "upset_plot.png", height = 5.5, width = 9, units = "in", dpi = 600)

#####################
### Joint testing ###
#####################

## Subset RNA, DNA counts; create new annot file
rna_subset <- list()
for (i in 1:4){
  rna_subset[[i]] <- subset(rna_counts[[i]],rna_counts[[i]]$refname %in% common_inserts)
}

dna_subset <- list()
for (i in 1:4){
  dna_subset[[i]] <- subset(dna_counts[[i]], dna_counts[[i]]$refname %in% common_inserts)
}

# separate out barcodes for each refname in annot file
rna_annot <- list()
for (i in 1:4){
  rna_annot[[i]] <- data.frame("bcnum_rep" = colnames(rna_subset[[i]][2:length(colnames(rna_subset[[i]]))]))
  rna_annot[[i]]$barcode <- sub("_.*","",rna_annot[[i]]$bcnum_rep)
  rna_annot[[i]]$rep <- sub(".*_","",rna_annot[[i]]$bcnum_rep)
}

# add metadata for annot file
rna_annot[[1]]$pos <- rep("up",length(rna_annot[[1]][1]))
rna_annot[[2]]$pos <- rep("down",length(rna_annot[[2]][1]))
rna_annot[[3]]$pos <- rep("up",length(rna_annot[[3]][1]))
rna_annot[[4]]$pos <- rep("down",length(rna_annot[[4]][1]))

rna_annot[[1]]$prom <- rep("INS",length(rna_annot[[1]][1]))
rna_annot[[2]]$prom <- rep("INS",length(rna_annot[[2]][1]))
rna_annot[[3]]$prom <- rep("SCP1",length(rna_annot[[3]][1]))
rna_annot[[4]]$prom <- rep("SCP1",length(rna_annot[[4]][1]))

# create joint identifier
for (i in 1:4){
  rna_annot[[i]]$bcnum_rep_pos_prom <- paste(rna_annot[[i]]$bcnum_rep,rna_annot[[i]]$pos,
                                             rna_annot[[i]]$prom,sep="_")
}

# make DNA annot file
dna_annot <- list()
for (i in 1:4){
  dna_annot[[i]] <- data.frame("bcnum" = colnames(dna_subset[[i]][2:length(colnames(dna_subset[[i]]))]))
  dna_annot[[i]]$barcode <- sub("_.*","",dna_annot[[i]]$bcnum)
  dna_annot[[i]]$lib <- rep("dna",length(dna_annot[[i]][1]))
}

dna_annot[[1]]$pos <- rep("up",length(dna_annot[[1]][1]))
dna_annot[[2]]$pos <- rep("down",length(dna_annot[[2]][1]))
dna_annot[[3]]$pos <- rep("up",length(dna_annot[[3]][1]))
dna_annot[[4]]$pos <- rep("down",length(dna_annot[[4]][1]))

dna_annot[[1]]$prom <- rep("INS",length(dna_annot[[1]][1]))
dna_annot[[2]]$prom <- rep("INS",length(dna_annot[[2]][1]))
dna_annot[[3]]$prom <- rep("SCP1",length(dna_annot[[3]][1]))
dna_annot[[4]]$prom <- rep("SCP1",length(dna_annot[[4]][1]))

for (i in 1:4){
  dna_annot[[i]]$bcnum_pos_prom <- paste(dna_annot[[i]]$bcnum,dna_annot[[i]]$pos,
                                         dna_annot[[i]]$prom,sep="_")
}

# create new rna counts file with updated annotations (joint)
new_rnacounts <- rna_counts
for (i in 1:4){
  colnames(new_rnacounts[[i]])[2:length(colnames(new_rnacounts[[i]]))] <- rna_annot[[i]]$bcnum_rep_pos_prom
}

new_dnacounts <- dna_counts
for (i in 1:4){
  colnames(new_dnacounts[[i]])[2:length(colnames(new_dnacounts[[i]]))] <- dna_annot[[i]]$bcnum_pos_prom
}

# create new joint dataframe across all constructs
new_rnadf <- merge(merge(merge(new_rnacounts[[1]],new_rnacounts[[2]],by="refname"),
                         new_rnacounts[[3]], by="refname"),new_rnacounts[[4]],by="refname")
new_dnadf <- merge(merge(merge(new_dnacounts[[1]],new_dnacounts[[2]],by="refname"),
                         new_dnacounts[[3]],by="refname"),new_dnacounts[[4]],by="refname")

# create new annotation df across all constructs
new_rnaannot <- rbind(rna_annot[[1]],rna_annot[[2]],rna_annot[[3]],rna_annot[[4]])
new_dnaannot <- rbind(dna_annot[[1]],dna_annot[[2]],dna_annot[[3]],dna_annot[[4]])

# joint analysis in MPRAnalyze
dna = new_dnadf
rownames(dna) <- dna$refname
dna <- dna[,-1]
dna <- as.matrix(dna)
rna = new_rnadf
rownames(rna) <- rna$refname
rna <- rna[,-1]
rna <- as.matrix(rna)

dna_annot = new_dnaannot
rownames(dna_annot) <- dna_annot$bcnum_pos_prom
dna_annot <- dna_annot[,c(1:5)]
dna_annot[,c(1,2,4,5)] <- lapply(dna_annot[,c(1,2,4,5)], factor)
rna_annot = new_rnaannot
rownames(rna_annot) <- rna_annot$bcnum_rep_pos_prom
rna_annot <- rna_annot[,c(1:5)]
rna_annot[,c(1:5)] <- lapply(rna_annot[,c(1:5)], factor)
rna_annot$config <- paste(rna_annot$pos,rna_annot$prom, sep="_")
dna_annot$config <- paste(dna_annot$pos,dna_annot$prom, sep="_")

write.table(dna, paste0(out_dir, "combined.dna_counts.tsv"), sep="\t", row.names = F)
write.table(rna, paste0(out_dir, "combined.rna_counts.tsv"), sep="\t", row.names = F)
write.table(dna_annot, paste0(out_dir, "combined.dna_annots.tsv"), sep="\t", row.names = F)
write.table(rna_annot, paste0(out_dir, "combined.rna_annots.tsv"), sep="\t", row.names = F)

obj = MpraObject(dnaCounts = dna, rnaCounts = rna, dnaAnnot = dna_annot, rnaAnnot = rna_annot)
obj <- estimateDepthFactors(obj, which.lib = "dna",
                            depth.estimator = "uq")
obj <- estimateDepthFactors(obj, lib.factor = "rep",
                            which.lib = "rna",
                            depth.estimator = "uq")

# fit comparative analysis model
comp = analyzeComparative(obj, dnaDesign = ~ barcode, rnaDesign = ~ rep + pos + prom,
                          reducedDesign = ~ rep, fit.se=TRUE)

# fxn to get significant results
getSig <- function(df) {
  df %>% mutate(sig = case_when(df$fdr <= 0.05 ~ "fdr_05",
                                df$fdr >= 0.05 ~ "not_sig")) %>%
    filter(fdr != 0)
}

# Wald test results for promoter and position
wald_res <- list()
wald_res[[1]] <- testCoefficient(comp, "pos", "up")
wald_res[[2]] <- testCoefficient(comp, "prom", "SCP1")
wald_res <- lapply(wald_res, getSig)
wald_res <- lapply(wald_res, function(x) {x %>% rownames_to_column(var = "refname")})

# fxn to get marginal results
getPermis <- function(df) {
  df %>% mutate(sig = case_when(df$fdr <= 0.1 ~ "fdr_1",
                                df$fdr > 0.1 ~ "not_sig"))
}

# make plotting df from Wald testing results
wald_plot_df <- merge(wald_res[[1]], wald_res[[2]], by = "refname")
wald_plot_df <- wald_plot_df %>%
  mutate(sig = as.factor(if_else(sig.x == "fdr_05" & sig.y == "fdr_05", "both",
                                 if_else(sig.x == "fdr_05" & sig.y != "fdr_05", "pos",
                                         if_else(sig.x != "fdr_05" & sig.y == "fdr_05", "prom",
                                                 "neither")))))
# adjust levels of which significance they show
wald_plot_df$sig <- factor(wald_plot_df$sig, levels = c("pos", "prom", "both","neither"))
table(wald_plot_df$sig)
# change shading of points based on significance level
wald_plot_df <- wald_plot_df %>%
  mutate(alpha = if_else(sig == "neither", FALSE, TRUE))

# alter logFC just for plotting -- calculated at base e, should be base 2
wald_plot_df <- wald_plot_df %>%
  mutate(corrected_logFC.x = log2(exp(logFC.x)),
         corrected_logFC.y = log2(exp(logFC.y)))

# main scatter for wald plot
wald_plot <- wald_plot_df %>% ggplot(aes(x = corrected_logFC.y, y = corrected_logFC.x, color = sig, alpha = alpha)) +
  geom_point(size = 1) + theme_linedraw(base_size = 14) +
  scale_alpha_discrete(range = c(0.1, 0.5), breaks = c(0.1, 0.5)) +
  scale_color_manual(values = c("#ffb14e", "#ea5f94", "#0000ff", "#e3e3e3"), name = "Significant\nEffect",
                     labels = c("Position", "Promoter", "Both", "Neither")) +
  theme(legend.position = c(0.1, 0.875),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(3, units = "pt"),
        plot.margin = margin(0,0,0,0,unit="cm"),
        text = element_text(family = "Helvetica")) +
  labs(x = expression("log"[2]*(FC)*","~"Promoter"), y = expression("log"[2]*(FC)*","~"Position"))

wald_dens1 <- wald_plot_df %>% filter(sig.y == "fdr_05") %>% ggplot(aes(x = corrected_logFC.y)) +
  geom_density(fill = "#ea5f94", alpha = 0.5) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0,unit="cm"))

wald_dens2 <- wald_plot_df %>% filter(sig.x == "fdr_05") %>% ggplot(aes(x = corrected_logFC.x)) +
  geom_density(fill = "#ffb14e", alpha = 0.5) +
  theme_void() +
  coord_flip() +
  theme(plot.margin = margin(0,0,0,0,unit="cm"))

full_wald_plot <- wald_dens1 + plot_spacer() +
  wald_plot + wald_dens2 +
  patchwork::plot_layout(
    ncol = 2, 
    nrow = 2, 
    heights = c(0.75, 4),
    widths = c(4, 0.75)
  )

ggsave(full_wald_plot,
       filename = paste0(fig_dir, "wald_plot.png"), units = "in",
       height = 7, width = 7)

#########################################
### Prep for annotation/motif overlap ###
#########################################

# add a bias statistic - aka signed z-score
wald_res <- lapply(wald_res, function(x) {x %>% mutate(zscore = ifelse(logFC < 0, statistic*-1, statistic))})

# now convert insert names to other convention for ease of using existing pipeline
wald_v2 <- lapply(wald_res, function(x) {x %>% mutate(basename = gsub("_ref|_alt", "", refname))})

# remove duplicate inserts used for allelic bias, keeping random Wald statistic
wald_v2 <- lapply(wald_v2, function(x) {x %>% arrange(basename) %>% distinct(., basename, .keep_all = TRUE)})
wald_v2 <- lapply(wald_v2, function(x) {merge(x, reg_annot[,c(1,6:8)], by.x = "refname", by.y = "name")})
wald_v2 <- lapply(wald_v2, function(x) {x %>% mutate(name = paste("TCs", chr, start, end, sep = "_"))})

# remove 12 duplicate inserts that have overlapping positions
wald_v2 <- lapply(wald_v2, function(x) {x %>% distinct(., name, .keep_all = TRUE)})
wald_v2 <- lapply(wald_v2, function(x) {x %>% column_to_rownames(., var = "name")})
write.table(wald_v2[[1]], paste0(out_dir, "wald_pos_zscore.tsv"), sep = "\t", quote = F)
write.table(wald_v2[[2]], paste0(out_dir, "wald_prom_zscore"), sep = "\t", quote = F)

#####################
### LASSO testing ###
#####################

# Perform overlap and calculate z-scores outside of R
# Load in scores with annotations

# load data

pos_lasso <- read.table(file = paste0(data_dir, "gwas_pos.zscores.tsv"), header = TRUE)
pos_lasso_negp <- read.table(file = paste0(data_dir, "gwas_pos.neg-log-p.zscores.tsv"), header = TRUE)

prom_lasso <- read.table(file = paste0(data_dir, "gwas_prom.zscores.tsv"), header=TRUE)
prom_lasso_negp <- read.table(file = paste0(data_dir, "gwas_prom.neg-log-p.zscores.tsv"), header = TRUE)

# establish chrom states in data

chrom_states <- c('X8_Genic_enhancer', 'X9_Active_enhancer_1', 'stretchEnhancer',
                  'typicalEnhancer', 'X10_Active_enhancer_2', 'X11_Weak_enhancer',
                  'X14_Bivalent_poised_TSS', 'X16_Repressed_polycomb', 'X17_Weak_repressed_polycomb',
                  'X18_Quiescent_low_signal', 'X1_Active_TSS', 'X2_Weak_TSS', 
                  'X3_Flanking_TSS', 'X5_Strong_transcription', 'X6_Weak_transcription',
                  'X8_Genic_enhancer', 'X9_Active_enhancer', 'stretchEnhancer', 'typicalEnhancer', 'atac_seq')

trimmed_motifs <- scan(paste0(data_dir, 'motif_list.encode_1KG.trimmed.txt'), character())
trimmed_motifs <- c('zscore', trimmed_motifs)

####
# pos full set + chrom state

df <- dplyr::select(pos_lasso, -which(grepl('disc', names(pos_lasso))))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit1 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit1)
l.lasso.min1 <- cvfit1$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef1 <- as.data.frame(as.matrix(cvfit1$glmnet.fit$beta))
all_coef1$mark <- row.names(all_coef1) # select s5 for plots

lasso.model1 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min1*0.8)
coef1 <- as.data.frame(as.matrix(lasso.model1$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef1 <- as.data.frame(as.matrix(coef(cvfit1, s="lambda.min")))

# change rownames into another column
coef_parsed1 <- setDT(coef1, keep.rownames = TRUE)[]
colnames(coef_parsed1) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef1 <- coef_parsed1[coef_parsed1$coef != 0, ]
nonzero_coef1 <- nonzero_coef1[order(nonzero_coef1$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed1, file=paste0(out_dir, "pos.lasso-zscore.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef1, file=paste0(out_dir, "pos.lasso-zscore.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef1, file=paste0(out_dir, "pos.lasso-zscore.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# pos full set - chrom state

df <- dplyr::select(pos_lasso, -which(grepl('disc', names(pos_lasso))))
df <- dplyr::select(df, -which(names(df) %in% chrom_states))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit2 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit2)
l.lasso.min2 <- cvfit2$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef2 <- as.data.frame(as.matrix(cvfit2$glmnet.fit$beta))
all_coef2$mark <- row.names(all_coef2) # select s5 for plots

lasso.model2 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min2*0.8)
coef2 <- as.data.frame(as.matrix(lasso.model2$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef2 <- as.data.frame(as.matrix(coef(cvfit2, s="lambda.min")))

# change rownames into another column
coef_parsed2 <- setDT(coef2, keep.rownames = TRUE)[]
colnames(coef_parsed2) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef2 <- coef_parsed2[coef_parsed2$coef != 0, ]
nonzero_coef2 <- nonzero_coef2[order(nonzero_coef2$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed2, file=paste0(out_dir, "pos.lasso-zscore.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef2, file=paste0(out_dir, "pos.lasso-zscore.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef2, file=paste0(out_dir, "pos.lasso-zscore.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# pos red set

df <- dplyr::select(pos_lasso, which(names(pos_lasso) %in% trimmed_motifs))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit3 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit3)
l.lasso.min3 <- cvfit3$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef3 <- as.data.frame(as.matrix(cvfit3$glmnet.fit$beta))
all_coef3$mark <- row.names(all_coef3) # select s5 for plots

lasso.model3 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min3*0.8)
coef3 <- as.data.frame(as.matrix(lasso.model3$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef3 <- as.data.frame(as.matrix(coef(cvfit3, s="lambda.min")))

# change rownames into another column
coef_parsed3 <- setDT(coef3, keep.rownames = TRUE)[]
colnames(coef_parsed3) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef3 <- coef_parsed3[coef_parsed3$coef != 0, ]
nonzero_coef3 <- nonzero_coef3[order(nonzero_coef3$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed3, file=paste0(out_dir, "pos.lasso-zscore.red.coeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef3, file=paste0(out_dir, "pos.lasso-zscore.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef3, file=paste0(out_dir, "pos.lasso-zscore.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom full set + chrom state

df <- dplyr::select(prom_lasso, -which(grepl('disc', names(prom_lasso))))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit4 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit4)
l.lasso.min4 <- cvfit4$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef4 <- as.data.frame(as.matrix(cvfit4$glmnet.fit$beta))
all_coef4$mark <- row.names(all_coef4) # select s5 for plots

lasso.model4 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min4*0.8)
coef4 <- as.data.frame(as.matrix(lasso.model4$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef4 <- as.data.frame(as.matrix(coef(cvfit4, s="lambda.min")))

# change rownames into another column
coef_parsed4 <- setDT(coef4, keep.rownames = TRUE)[]
colnames(coef_parsed4) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef4 <- coef_parsed4[coef_parsed4$coef != 0, ]
nonzero_coef4 <- nonzero_coef4[order(nonzero_coef4$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed4, file=paste0(out_dir, "prom.lasso-zscore.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef4, file=paste0(out_dir, "prom.lasso-zscore.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef4, file=paste0(out_dir, "prom.lasso-zscore.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom full set - chrom state

df <- dplyr::select(prom_lasso, -which(grepl('disc', names(prom_lasso))))
df <- dplyr::select(df, -which(names(df) %in% chrom_states))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit5 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit5)
l.lasso.min5 <- cvfit5$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef5 <- as.data.frame(as.matrix(cvfit5$glmnet.fit$beta))
all_coef5$mark <- row.names(all_coef5) # select s5 for plots

lasso.model5 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min5*0.8)
coef5 <- as.data.frame(as.matrix(lasso.model5$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef5 <- as.data.frame(as.matrix(coef(cvfit5, s="lambda.min")))

# change rownames into another column
coef_parsed5 <- setDT(coef5, keep.rownames = TRUE)[]
colnames(coef_parsed5) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef5 <- coef_parsed5[coef_parsed5$coef != 0, ]
nonzero_coef5 <- nonzero_coef5[order(nonzero_coef5$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed5, file=paste0(out_dir, "prom.lasso-zscore.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef5, file=paste0(out_dir, "prom.lasso-zscore.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef5, file=paste0(out_dir, "prom.lasso-zscore.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom red set

df <- dplyr::select(prom_lasso, which(names(prom_lasso) %in% trimmed_motifs))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit6 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit6)
l.lasso.min6 <- cvfit6$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef6 <- as.data.frame(as.matrix(cvfit6$glmnet.fit$beta))
all_coef6$mark <- row.names(all_coef6) # select s5 for plots

lasso.model6 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min6*0.8)
coef6 <- as.data.frame(as.matrix(lasso.model6$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef6 <- as.data.frame(as.matrix(coef(cvfit6, s="lambda.min")))

# change rownames into another column
coef_parsed6 <- setDT(coef6, keep.rownames = TRUE)[]
colnames(coef_parsed6) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef6 <- coef_parsed6[coef_parsed6$coef != 0, ]
nonzero_coef6 <- nonzero_coef6[order(nonzero_coef6$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed6, file=paste0(out_dir, "prom.lasso-zscore.red.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef6, file=paste0(out_dir, "prom.lasso-zscore.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef6, file=paste0(out_dir, "prom.lasso-zscore.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

####
# pos full set + chrom state - neg pval

df <- dplyr::select(pos_lasso_negp, -which(grepl('disc', names(pos_lasso_negp))))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit7 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit7)
l.lasso.min7 <- cvfit7$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef7 <- as.data.frame(as.matrix(cvfit7$glmnet.fit$beta))
all_coef7$mark <- row.names(all_coef7) # select s5 for plots

lasso.model7 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min7*0.8)
coef7 <- as.data.frame(as.matrix(lasso.model7$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef7 <- as.data.frame(as.matrix(coef(cvfit7, s="lambda.min")))

# change rownames into another column
coef_parsed7 <- setDT(coef7, keep.rownames = TRUE)[]
colnames(coef_parsed7) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef7 <- coef_parsed7[coef_parsed7$coef != 0, ]
nonzero_coef7 <- nonzero_coef7[order(nonzero_coef7$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed7, file=paste0(out_dir, "pos.lasso-zscore-negp.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef7, file=paste0(out_dir, "pos.lasso-zscore-negp.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef7, file=paste0(out_dir, "pos.lasso-zscore-negp.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# pos full set - chrom state

df <- dplyr::select(pos_lasso_negp, -which(grepl('disc', names(pos_lasso_negp))))
df <- dplyr::select(df, -which(names(df) %in% chrom_states))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit8 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit8)
l.lasso.min8 <- cvfit8$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef8 <- as.data.frame(as.matrix(cvfit8$glmnet.fit$beta))
all_coef8$mark <- row.names(all_coef8) # select s5 for plots

lasso.model8 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min8*0.8)
coef8 <- as.data.frame(as.matrix(lasso.model8$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef8 <- as.data.frame(as.matrix(coef(cvfit8, s="lambda.min")))

# change rownames into another column
coef_parsed8 <- setDT(coef8, keep.rownames = TRUE)[]
colnames(coef_parsed8) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef8 <- coef_parsed8[coef_parsed8$coef != 0, ]
nonzero_coef8 <- nonzero_coef8[order(nonzero_coef8$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed8, file=paste0(out_dir, "pos.lasso-zscore-negp.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef8, file=paste0(out_dir, "pos.lasso-zscore-negp.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef8, file=paste0(out_dir, "pos.lasso-zscore-negp.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# pos red set

df <- dplyr::select(pos_lasso_negp, which(names(pos_lasso_negp) %in% trimmed_motifs))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit9 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit9)
l.lasso.min9 <- cvfit9$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef9 <- as.data.frame(as.matrix(cvfit9$glmnet.fit$beta))
all_coef9$mark <- row.names(all_coef9) # select s5 for plots

lasso.model9 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min9*0.8)
coef9 <- as.data.frame(as.matrix(lasso.model9$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef9 <- as.data.frame(as.matrix(coef(cvfit9, s="lambda.min")))

# change rownames into another column
coef_parsed9 <- setDT(coef9, keep.rownames = TRUE)[]
colnames(coef_parsed9) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef9 <- coef_parsed9[coef_parsed9$coef != 0, ]
nonzero_coef9 <- nonzero_coef9[order(nonzero_coef9$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed9, file=paste0(out_dir, "pos.lasso-zscore-negp.red.coeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef9, file=paste0(out_dir, "pos.lasso-zscore-negp.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef9, file=paste0(out_dir, "pos.lasso-zscore-negp.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom full set + chrom state

df <- dplyr::select(prom_lasso_negp, -which(grepl('disc', names(prom_lasso_negp))))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit10 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit10)
l.lasso.min10 <- cvfit10$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef10 <- as.data.frame(as.matrix(cvfit10$glmnet.fit$beta))
all_coef10$mark <- row.names(all_coef10) # select s5 for plots

lasso.model10 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min10*0.8)
coef10 <- as.data.frame(as.matrix(lasso.model10$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef10 <- as.data.frame(as.matrix(coef(cvfit10, s="lambda.min")))

# change rownames into another column
coef_parsed10 <- setDT(coef10, keep.rownames = TRUE)[]
colnames(coef_parsed10) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef10 <- coef_parsed10[coef_parsed10$coef != 0, ]
nonzero_coef10 <- nonzero_coef10[order(nonzero_coef10$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed10, file=paste0(out_dir, "prom.lasso-zscore-negp.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef10, file=paste0(out_dir, "prom.lasso-zscore-negp.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef10, file=paste0(out_dir, "prom.lasso-zscore-negp.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom full set - chrom state

df <- dplyr::select(prom_lasso_negp, -which(grepl('disc', names(prom_lasso_negp))))
df <- dplyr::select(df, -which(names(df) %in% chrom_states))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit11 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit11)
l.lasso.min11 <- cvfit11$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef11 <- as.data.frame(as.matrix(cvfit11$glmnet.fit$beta))
all_coef11$mark <- row.names(all_coef11) # select s5 for plots

lasso.model11 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min11*0.8)
coef11 <- as.data.frame(as.matrix(lasso.model11$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef11 <- as.data.frame(as.matrix(coef(cvfit11, s="lambda.min")))

# change rownames into another column
coef_parsed11 <- setDT(coef11, keep.rownames = TRUE)[]
colnames(coef_parsed11) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef11 <- coef_parsed11[coef_parsed11$coef != 0, ]
nonzero_coef11 <- nonzero_coef11[order(nonzero_coef11$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed11, file=paste0(out_dir, "prom.lasso-zscore-negp.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef11, file=paste0(out_dir, "prom.lasso-zscore-negp.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef11, file=paste0(out_dir, "prom.lasso-zscore-negp.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom red set

df <- dplyr::select(prom_lasso_negp, which(names(prom_lasso_negp) %in% trimmed_motifs))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit12 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit12)
l.lasso.min12 <- cvfit12$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef12 <- as.data.frame(as.matrix(cvfit12$glmnet.fit$beta))
all_coef12$mark <- row.names(all_coef12) # select s5 for plots

lasso.model12 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min12*0.8)
coef12 <- as.data.frame(as.matrix(lasso.model12$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef12 <- as.data.frame(as.matrix(coef(cvfit12, s="lambda.min")))

# change rownames into another column
coef_parsed12 <- setDT(coef12, keep.rownames = TRUE)[]
colnames(coef_parsed12) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef12 <- coef_parsed12[coef_parsed12$coef != 0, ]
nonzero_coef12 <- nonzero_coef12[order(nonzero_coef12$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed12, file=paste0(out_dir, "prom.lasso-zscore-negp.red.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef12, file=paste0(out_dir, "prom.lasso-zscore-negp.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef12, file=paste0(out_dir, "prom.lasso-zscore-negp.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# pos red set + chrom states

df <- dplyr::select(pos_lasso, which(names(pos_lasso) %in% c(trimmed_motifs, chrom_states)))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit13 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit13)
l.lasso.min13 <- cvfit13$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef13 <- as.data.frame(as.matrix(cvfit13$glmnet.fit$beta))
all_coef13$mark <- row.names(all_coef13) # select s5 for plots

lasso.model13 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min13*0.8)
coef13 <- as.data.frame(as.matrix(lasso.model13$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef13 <- as.data.frame(as.matrix(coef(cvfit13, s="lambda.min")))

# change rownames into another column
coef_parsed13 <- setDT(coef13, keep.rownames = TRUE)[]
colnames(coef_parsed13) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef13 <- coef_parsed13[coef_parsed13$coef != 0, ]
nonzero_coef13 <- nonzero_coef13[order(nonzero_coef13$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed13, file=paste0(out_dir, "pos.lasso-zscore.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef13, file=paste0(out_dir, "pos.lasso-zscore.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef13, file=paste0(out_dir, "pos.lasso-zscore.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom red set + chrom states

df <- dplyr::select(prom_lasso, which(names(prom_lasso) %in% c(trimmed_motifs, chrom_states)))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit14 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit14)
l.lasso.min14 <- cvfit14$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef14 <- as.data.frame(as.matrix(cvfit14$glmnet.fit$beta))
all_coef14$mark <- row.names(all_coef14) # select s5 for plots

lasso.model14 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min14*0.8)
coef14 <- as.data.frame(as.matrix(lasso.model14$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef14 <- as.data.frame(as.matrix(coef(cvfit14, s="lambda.min")))

# change rownames into another column
coef_parsed14 <- setDT(coef14, keep.rownames = TRUE)[]
colnames(coef_parsed14) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef14 <- coef_parsed14[coef_parsed14$coef != 0, ]
nonzero_coef14 <- nonzero_coef14[order(nonzero_coef14$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed14, file=paste0(out_dir, "prom.lasso-zscore.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef14, file=paste0(out_dir, "prom.lasso-zscore.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef14, file=paste0(out_dir, "prom.lasso-zscore.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# pos red set + chrom states

df <- dplyr::select(pos_lasso_negp, which(names(pos_lasso_negp) %in% c(trimmed_motifs, chrom_states)))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit15 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit15)
l.lasso.min15 <- cvfit15$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef15 <- as.data.frame(as.matrix(cvfit15$glmnet.fit$beta))
all_coef15$mark <- row.names(all_coef15) # select s5 for plots

lasso.model15 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min15*0.8)
coef15 <- as.data.frame(as.matrix(lasso.model15$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef15 <- as.data.frame(as.matrix(coef(cvfit15, s="lambda.min")))

# change rownames into another column
coef_parsed15 <- setDT(coef15, keep.rownames = TRUE)[]
colnames(coef_parsed15) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef15 <- coef_parsed15[coef_parsed15$coef != 0, ]
nonzero_coef15 <- nonzero_coef15[order(nonzero_coef15$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed15, file=paste0(out_dir, "pos.lasso-zscore-negp.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef15, file=paste0(out_dir, "pos.lasso-zscore-negp.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef15, file=paste0(out_dir, "pos.lasso-zscore-negp.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

###
# prom red set + chrom states

df <- dplyr::select(prom_lasso_negp, which(names(prom_lasso_negp) %in% c(trimmed_motifs, chrom_states)))
annots <- colnames(df)

set.seed(330)

# create model matrix and zscore
x <- model.matrix(zscore~., df)[,-1]
y <- df$zscore

# fit a lasso regression model to x and y
cvfit16 <- cv.glmnet(x, y, alpha = 1)
plot(cvfit16)
l.lasso.min16 <- cvfit16$lambda.min
#abline(log(l.lasso.min*0.8))
all_coef16 <- as.data.frame(as.matrix(cvfit16$glmnet.fit$beta))
all_coef16$mark <- row.names(all_coef16) # select s5 for plots

lasso.model16 <- glmnet(x, y, alpha = 1, lambda = l.lasso.min16*0.8)
coef16 <- as.data.frame(as.matrix(lasso.model16$beta))

# generate dataframe with matrix format for lambda that gives minimum MSE for fitted model 
coef16 <- as.data.frame(as.matrix(coef(cvfit16, s="lambda.min")))

# change rownames into another column
coef_parsed16 <- setDT(coef16, keep.rownames = TRUE)[]
colnames(coef_parsed16) <- c('marks', 'coef')

# create dataframe with sorted nonzero coefficients
nonzero_coef16 <- coef_parsed16[coef_parsed16$coef != 0, ]
nonzero_coef16 <- nonzero_coef16[order(nonzero_coef16$coef, decreasing=TRUE), ]

# create final output
# write.table(coef_parsed16, file=paste0(out_dir, "prom.lasso-zscore-negp.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef16, file=paste0(out_dir, "prom.lasso-zscore-negp.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
# write.table(all_coef16, file=paste0(out_dir, "prom.lasso-zscore-negp.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

# plot LASSO results
tf_info = read.table("/home/albanus/data/motifs/motif_library_information/pwmInfo.dat", header = T,row.names = NULL, fill = T)

## position full

pos_full_df <- nonzero_coef7[nonzero_coef7$marks != "(Intercept)",]

# several names to replace
pos_full_df$marks[pos_full_df$marks=="MA0060"] <- "NFYA"
pos_full_df$marks[pos_full_df$marks=="MA0502"] <- "NFYB"
pos_full_df$marks[pos_full_df$marks=="J1793"] <- "BATF3"
pos_full_df$marks[pos_full_df$marks=="J2733"] <- "HOXC13"
pos_full_df$marks[pos_full_df$marks=="J2843"] <- "ISX"
pos_full_df$marks[pos_full_df$marks=="J1078"] <- "POU3F1"
pos_full_df$marks[pos_full_df$marks=="J1398"] <- "TBX21"
pos_full_df$marks[pos_full_df$marks=="J0418"] <- "ELF1"
pos_full_df$marks[pos_full_df$marks=="J1093"] <- "POU3F2"
pos_full_df$marks[pos_full_df$marks=="J3168"] <- "PRRX1"
pos_full_df$marks[pos_full_df$marks=="J2763"] <- "HOXD13"
pos_full_df$marks[pos_full_df$marks=="MA0078"] <- "SOX17"
pos_full_df$marks[pos_full_df$marks=="MA0495"] <- "MAFF"
pos_full_df$marks[pos_full_df$marks=="J2023"] <- "FOXB1"
pos_full_df$marks[pos_full_df$marks=="J1113"] <- "POU3F4"
pos_full_df$marks[pos_full_df$marks=="MA0498"] <- "MEIS1"
pos_full_df$marks[pos_full_df$marks=="MA0073"] <- "RREB1"


pos_full_plot <- pos_full_df %>% arrange(marks, coef) %>% mutate(abs_coef = abs(coef)) %>% 
  ggplot(aes(x = reorder(marks, coef), y = coef)) +
  geom_bar(stat = "identity", fill = "#ffb14e", color = "#FFFFFF") +
  geom_hline(yintercept = 0.0) +
  theme_linedraw(base_size = 13) +
  scale_y_continuous(limits = c(-0.03, 0.22), breaks = seq(-0.03, 0.22, by = 0.03),
                     labels = c("-0.03", "0.0", "0.03", "0.06", "0.09", "0.12", "0.15", "0.18", "0.21")) +
  scale_x_discrete(position = "bottom") +
  theme(text = element_text(family = "Helvetica")) +
  coord_flip() +
  labs(y = "Coefficient",
       x = "TF Motif or Genomic Annotation\n")

ggsave(pos_full_plot,
       filename = paste0(fig_dir, "pos_full_plot.png"), units = "in",
       height = 8, width = 5)


## position red + chrom states

pos_rc_df <- nonzero_coef15[nonzero_coef15$marks != "(Intercept)",]

pos_rc_df_filt <- pos_rc_df
pos_rc_df_filt$coef_abs <- abs(pos_rc_df_filt$coef)
pos_rc_df_filt <- pos_rc_df_filt %>% arrange(-coef_abs)
pos_rc_df_filt <- pos_rc_df_filt[c(1:25),]

pos_rc_plot <- pos_rc_df_filt %>% arrange(marks, coef) %>% mutate(abs_coef = abs(coef)) %>% 
  ggplot(aes(x = reorder(marks, coef), y = coef)) +
  geom_bar(stat = "identity", fill = "#ffb14e", color = "#FFFFFF") +
  geom_hline(yintercept = 0.0) +
  theme_linedraw(base_size = 13) +
  scale_y_continuous(limits = c(-0.05, 0.33), breaks = seq(-0.05, 0.30, by = 0.05),
                     labels = c("-0.05", "0.0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30")) +
  scale_x_discrete(position = "bottom") +
  theme(text = element_text(family = "Helvetica")) +
  coord_flip() +
  labs(y = "Coefficient",
       x = "TF Motif or Genomic Annotation\n")

ggsave(pos_rc_plot,
       filename = paste0(fig_dir, "pos_rc_plot.png"), units = "in",
       height = 6.5, width = 5)

## promoter full

prom_full_df <- nonzero_coef10[nonzero_coef10$marks != "(Intercept)",]

# update some annotations that are not in readable format
prom_full_df$marks[prom_full_df$marks=="MA0101"] <- "RELA"
prom_full_df$marks[prom_full_df$marks=="J3918"] <- "SRY"

prom_full_plot <- prom_full_df %>% arrange(marks, coef) %>% mutate(abs_coef = abs(coef)) %>% 
  ggplot(aes(x = reorder(marks, coef), y = coef)) +
  geom_bar(stat = "identity", fill = "#ea5f94", color = "#FFFFFF") +
  geom_hline(yintercept = 0.0) +
  theme_linedraw() +
  scale_y_continuous(limits = c(-0.12, 0.03), breaks = seq(-0.12, 0.03, by = 0.03),
                     labels = c("-0.12", "-0.09", "-0.06", "-0.03", "0.0", "0.03")) +
  scale_x_discrete(position = "bottom") +
  theme(text = element_text(family = "Helvetica"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  coord_flip() +
  labs(y = "Coefficient",
       x = "TF Motif")

ggsave(prom_full_plot,
       filename = paste0(fig_dir, "prom_full_plot.png"), units = "in",
       height = 4.8, width = 3.75)

## promoter red + chrom states

prom_rc_df <- nonzero_coef16[nonzero_coef16$marks != "(Intercept)",]

prom_rc_plot <- prom_rc_df %>% arrange(marks, coef) %>% mutate(abs_coef = abs(coef)) %>% 
  ggplot(aes(x = reorder(marks, coef), y = coef)) +
  geom_bar(stat = "identity", fill = "#ea5f94", color = "#FFFFFF") +
  geom_hline(yintercept = 0.0) +
  theme_linedraw(base_size = 11) +
  scale_y_continuous(limits = c(-0.12, 0), breaks = seq(-0.12, 0, by = 0.03), labels = c("-0.12", "-0.09", "-0.06", "-0.03", "0.0")) +
  scale_x_discrete(position = "bottom") +
  theme(text = element_text(family = "Helvetica")) +
  coord_flip() +
  labs(y = "Coefficient",
       x = "TF Motif or Genomic Annotation\n")

ggsave(prom_rc_plot,
       filename = paste0(fig_dir, "prom_rc_plot.png"), units = "in",
       height = 4, width = 3)
