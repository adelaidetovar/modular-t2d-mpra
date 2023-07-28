#########################
### Main analysis ###
#########################

mpra_packages <- c("tidyverse", "factoextra", "data.table", "MPRAnalyze", "reshape2",
                   "ggtext", "UpSetR")
lapply(mpra_packages, require, character.only=TRUE)

# mk work dir
work_dir <- "/lab/work/tovar/2023_07_mpra-manuscript/"

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

#########################
### Initial data load ###
#########################

# Helper functions
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

# RNA counts from MPRA
rna_counts <- list()
rna_counts[[3]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_UP.minDNA10.minBarcode2.rna_counts.tsv'))
rna_counts[[4]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_dwn.minDNA10.minBarcode2.rna_counts.tsv'))
rna_counts[[1]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_UP.minDNA10.minBarcode2.rna_counts.tsv'))
rna_counts[[2]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_dwn.minDNA10.minBarcode2.rna_counts.tsv'))
names(rna_counts) <- c("UpINS","DwnINS","UpSCP1","DwnSCP1")

# DNA counts from MPRA

dna_counts <- list()
dna_counts[[3]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_UP.minDNA10.minBarcode2.dna_counts.tsv'))
dna_counts[[4]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_dwn.minDNA10.minBarcode2.dna_counts.tsv'))
dna_counts[[1]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_UP.minDNA10.minBarcode2.dna_counts.tsv'))
dna_counts[[2]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_dwn.minDNA10.minBarcode2.dna_counts.tsv'))
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

# plot PCA
pca_data=merge(pca_data,col_data)
percent_var = round(100*(res_pca$sdev^2/sum(res_pca$sdev^2)),digits=2)
pca_plot <- ggplot(data=pca_data,aes(x = PC1,y=PC2,fill=promoter,shape=config)) +
    geom_point(size=3) +
    scale_shape_manual(values=c(21,24),
                      name="Insert Position") +
    scale_fill_manual(values=c("#FA8775","#774BE5"),
                      name="Promoter") +
    theme_linedraw(base_size=8) +
    theme(plot.title=element_text(hjust=0.5,face="bold")) +
    ggtitle("Principal Components Analysis") +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    guides(fill = guide_legend(override.aes=list(shape=22,color=NA)))

ggsave(plot=pca_plot,filename=paste0(fig_dir,'pca_inserts_minDNA10.minBarcode2.png'),
        unit="in",width=4,height=3,dpi=300)

############
### ECDF ###
############

indiv_res <- list()
indiv_res[[3]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_UP.minDNA10.minBarcode2.type_norm.mpra.tsv'))
indiv_res[[4]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_SCP1_promoter_dwn.minDNA10.minBarcode2.type_norm.mpra.tsv'))
indiv_res[[1]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_UP.minDNA10.minBarcode2.type_norm.mpra.tsv'))
indiv_res[[2]] <- read.delim(paste0(data_dir, 'enh_and_oth_and_pro.config_INS1_promoter_dwn.minDNA10.minBarcode2.type_norm.mpra.tsv'))
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
  theme(legend.text = element_markdown())

ggsave(plot = ecdf_plot, paste0(fig_dir, "ecdf.png"), units = "in", height = 3, width = 5)

###################
### UpSet plots ###
###################

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


indiv_sig <- list()
for (i in 1:4){
  indiv_sig[[i]] <- getSigQuant(indiv_res[[i]])
}

list_sig <- list(UpINS = getSigInserts(indiv_sig[[1]]),
                 DownINS = getSigInserts(indiv_sig[[2]]),
                 UpSCP1 = getSigInserts(indiv_sig[[3]]),
                 DownSCP1 = getSigInserts(indiv_sig[[4]]))

upset_plot <- upset(fromList(list_sig), sets = c("UpSCP1","DownSCP1","UpINS", "DownINS"),
                    keep.order = TRUE, set_size.scale_max = 8500, set_size.show = TRUE,
                    sets.bar.color = c("#F0E442", "#E69F00", "#56B4E9", "#0072B2"),
                    text.scale = 0.75)

# change ratio of plot grid to prevent labels from cutting off
upset_plot$mb.ratio <- c(0.6, 0.4)

png(filename = paste0(fig_dir, "upset_plot.png"),
    width = 8, height = 4, units = "in", res = 300)
upset_plot
dev.off()

#####################
### Joint testing ###
#####################

## Subset RNA, DNA counts; create new annot file
subRna <- list()
for (i in 1:4){
  subRna[[i]] <- subset(rnaCounts[[i]],rnaCounts[[i]]$refname %in% commonInserts)
}

subDna <- list()
for (i in 1:4){
  subDna[[i]] <- subset(dnaCounts[[i]], dnaCounts[[i]]$refname %in% commonInserts)
}

rnaAnnot <- list()
for (i in 1:4){
  rnaAnnot[[i]] <- data.frame("bcnum_rep" = colnames(subRna[[i]][2:length(colnames(subRna[[i]]))]))
  rnaAnnot[[i]]$barcode <- sub("_.*","",rnaAnnot[[i]]$bcnum_rep)
  rnaAnnot[[i]]$rep <- sub(".*_","",rnaAnnot[[i]]$bcnum_rep)
}

rnaAnnot[[1]]$pos <- rep("up",length(rnaAnnot[[1]][1]))
rnaAnnot[[2]]$pos <- rep("down",length(rnaAnnot[[2]][1]))
rnaAnnot[[3]]$pos <- rep("up",length(rnaAnnot[[3]][1]))
rnaAnnot[[4]]$pos <- rep("down",length(rnaAnnot[[4]][1]))

rnaAnnot[[1]]$prom <- rep("SCP1",length(rnaAnnot[[1]][1]))
rnaAnnot[[2]]$prom <- rep("SCP1",length(rnaAnnot[[2]][1]))
rnaAnnot[[3]]$prom <- rep("INS",length(rnaAnnot[[3]][1]))
rnaAnnot[[4]]$prom <- rep("INS",length(rnaAnnot[[4]][1]))

for (i in 1:4){
  rnaAnnot[[i]]$bcnum_rep_pos_prom <- paste(rnaAnnot[[i]]$bcnum_rep,rnaAnnot[[i]]$pos,rnaAnnot[[i]]$prom,sep="_")
}

dnaAnnot <- list()
for (i in 1:4){
  dnaAnnot[[i]] <- data.frame("bcnum" = colnames(subDna[[i]][2:length(colnames(subDna[[i]]))]))
  dnaAnnot[[i]]$barcode <- sub("_.*","",dnaAnnot[[i]]$bcnum)
  dnaAnnot[[i]]$lib <- rep("dna",length(dnaAnnot[[i]][1]))
}

dnaAnnot[[1]]$pos <- rep("up",length(dnaAnnot[[1]][1]))
dnaAnnot[[2]]$pos <- rep("down",length(dnaAnnot[[2]][1]))
dnaAnnot[[3]]$pos <- rep("up",length(dnaAnnot[[3]][1]))
dnaAnnot[[4]]$pos <- rep("down",length(dnaAnnot[[4]][1]))

dnaAnnot[[1]]$prom <- rep("SCP1",length(dnaAnnot[[1]][1]))
dnaAnnot[[2]]$prom <- rep("SCP1",length(dnaAnnot[[2]][1]))
dnaAnnot[[3]]$prom <- rep("INS",length(dnaAnnot[[3]][1]))
dnaAnnot[[4]]$prom <- rep("INS",length(dnaAnnot[[4]][1]))

for (i in 1:4){
  dnaAnnot[[i]]$bcnum_pos_prom <- paste(dnaAnnot[[i]]$bcnum,dnaAnnot[[i]]$pos,dnaAnnot[[i]]$prom,sep="_")
}

head(rnaCounts[[1]])
rnaCountsNew <- rnaCounts
for (i in 1:4){
  colnames(rnaCountsNew[[i]])[2:length(colnames(rnaCountsNew[[i]]))] <- rnaAnnot[[i]]$bcnum_rep_pos_prom
}

dnaCountsNew <- dnaCounts
for (i in 1:4){
  colnames(dnaCountsNew[[i]])[2:length(colnames(dnaCountsNew[[i]]))] <- dnaAnnot[[i]]$bcnum_pos_prom
}

newRna <- merge(merge(merge(rnaCountsNew[[1]],rnaCountsNew[[2]],by="refname"),rnaCountsNew[[3]], by="refname"),rnaCountsNew[[4]],by="refname")
newDna <- merge(merge(merge(dnaCountsNew[[1]],dnaCountsNew[[2]],by="refname"),dnaCountsNew[[3]],by="refname"),dnaCountsNew[[4]],by="refname")

newRnaAnnot <- rbind(rnaAnnot[[1]],rnaAnnot[[2]],rnaAnnot[[3]],rnaAnnot[[4]])
newDnaAnnot <- rbind(dnaAnnot[[1]],dnaAnnot[[2]],dnaAnnot[[3]],dnaAnnot[[4]])


#MPRAnalyze
dna = newDna
rownames(dna) <- dna$refname
dna <- dna[,-1]
dna <- as.matrix(dna)
rna = newRna
rownames(rna) <- rna$refname
rna <- rna[,-1]
rna <- as.matrix(rna)

dna_annot = newDnaAnnot
rownames(dna_annot) <- dna_annot$bcnum_pos_prom
dna_annot <- dna_annot[,c(1:5)]
dna_annot[,c(1,2,4,5)] <- lapply(dna_annot[,c(1,2,4,5)], factor)
rna_annot = newRnaAnnot
rownames(rna_annot) <- rna_annot$bcnum_rep_pos_prom
rna_annot <- rna_annot[,c(1:5)]
rna_annot[,c(1:5)] <- lapply(rna_annot[,c(1:5)], factor)
rna_annot$config <- paste(rna_annot$pos,rna_annot$prom, sep="_")
dna_annot$config <- paste(dna_annot$pos,dna_annot$prom, sep="_")

write.table(dna, "combined.dna_counts.tsv", sep="\t")
write.table(rna, "combined.rna_counts.tsv", sep="\t")
write.table(dna_annot, "combined.dna_annots.tsv", sep="\t")
write.table(rna_annot, "combined.rna_annots.tsv", sep="\t")

obj = MpraObject(dnaCounts = dna, rnaCounts = rna, dnaAnnot = dna_annot, rnaAnnot = rna_annot)
obj <- estimateDepthFactors(obj, which.lib = "dna",
                            depth.estimator = "uq")
obj <- estimateDepthFactors(obj, lib.factor = "rep",
                            which.lib = "rna",
                            depth.estimator = "uq")

obj = analyzeQuantification(obj = obj, dnaDesign = ~ barcode, rnaDesign = ~ rep + config)
output <- testEmpirical(obj)

comp = analyzeComparative(obj, dnaDesign = ~ barcode, rnaDesign = ~ rep + interaction(pos,prom),
                          reducedDesign = ~ rep, fit.se=TRUE)
alpha <- getAlpha(comp,by.factor="all")
results <- testLrt(comp)
results <- testCoefficient(comp_two, interaction(pos,prom))

comp_two = analyzeComparative(obj, dnaDesign = ~ barcode, rnaDesign = ~ rep + prom,
                              reducedDesign = ~ rep, fit.se=TRUE)
results <- testCoefficient(comp_two, "config", "SCP1")

alpha <- getAlpha(comp_two,by.factor="all")
results <- testLrt(comp_two)
results <- testCoefficient(comp_two, "prom",contrast = "SCP1")


library(tidyverse)
library(MPRAnalyze)

getSigQuant <- function(df) {
  df %>% mutate(fdr = p.adjust(df$pval.zscore, method = "BH"),
                sig = case_when(fdr <= 0.05 ~ "fdr_05",
                                fdr >= 0.05 & df$pval.zscore <= 0.05 ~ "nominal_sig",
                                df$pval.zscore > 0.05 ~ "not_sig")) %>%
    filter(fdr != 0)
}

getSig <- function(df) {
  df %>% mutate(sig = case_when(df$fdr <= 0.05 ~ "fdr_05",
                                df$fdr >= 0.05 ~ "not_sig")) %>%
    filter(fdr != 0)
}

getPermis <- function(df) {
  df %>% mutate(sig = case_when(df$fdr <= 0.1 ~ "fdr_1",
                                df$fdr > 0.1 ~ "not_sig"))
}

load("../output/combined_v2/quant.rds")
quant_v2 <- testEmpirical(obj)
quant_v2 <- getSigQuant(quant_v2)
length(which(quant_v2$fdr < 0.05))
write.table(quant_v2, "../output/combined_v2/quant.tsv", sep = "\t", quote = F)

quantInserts05 <- rownames(quant_v2[quant_v2$fdr < 0.05,])
quantInserts10 <- rownames(quant_v2[quant_v2$fdr < 0.1,])

p <- quant_v2 %>% count(sig) %>% ggplot(aes(fill = sig, y = n, x = '')) +
  geom_bar(position = "fill", stat = "identity") + labs(y = "Proportion", x = NULL) +
  theme_bw(base_size = 12) + scale_fill_discrete(name = "Transcriptional\nActivity",
                                                 labels = c("FDR < 0.05",
                                                            "Nominal",
                                                            "Not Sig.")) +
  geom_text(x = 1, y = (7484/11656)/2, label = "7484") +
  geom_text(x = 1, y = (7484/11656) + (828/11656)/2, label = "828") +
  geom_text(x = 1, y = (7484+828)/11656 + (3344/11656)/2 , label = "3344") +
  theme(axis.ticks.x = element_blank())
ggsave(p, filename = "../figs/joint_quant.png", units = "in",
       height = 5, width = 5)

LRT_v2 <- list()
load("../output/combined_v2/comp.rds")
LRT_v2[[1]] <- testLrt(comp)
write.table(LRT_v2[[1]], "../output/combined_v2/comp_LRT.tsv", sep = "\t", quote = F)

# pos test
load("../output/combined_v2/comp_two.rds")
LRT_v2[[2]] <- testLrt(comp_two)
write.table(LRT_v2[[2]], "../output/combined_v2/comp_two_LRT.tsv", sep = "\t", quote = F)

# prom test
load("../output/combined_v2/comp_three.rds")
LRT_v2[[3]] <- testLrt(comp_three)
write.table(LRT_v2[[3]], "../output/combined_v2/comp_three_LRT.tsv", sep = "\t", quote = F)
names(LRT_v2) <- c("comp", "comp_two", "comp_three")

LRT_v2 <- lapply(LRT_v2, getSig)

# Wald test results for promoter and position
wald_v2 <- list()
wald_v2[[1]] <- testCoefficient(comp, "pos", "up")
wald_v2[[2]] <- testCoefficient(comp, "prom", "SCP1")
wald_v2 <- lapply(wald_v2, getSig)


p <- ggplot(wald_v2[[1]],aes(x = logFC,y = -log10(pval), col = sig)) +
  geom_point(alpha = 0.5,size = 1) + scale_color_manual(values = c("red", "gray", "black")) +
  xlab("logFC") + ylab("-log10 p-value") + theme_linedraw(base_size = 12) + 
  theme(legend.position = "none") + scale_y_continuous(limits = c(-2,35))

ggsave(filename = "../figs/jointwald_pos_ashr_volcano_lfdr_nol.png", plot = p, units = "in", width = 4, height = 4)

which(wald_v2[[]])

# Position volcano plot

plot1 <- ggplot(wald_v2[[2]],aes(x = -logFC,y = -log10(pval), col = sig)) +
  geom_point(alpha = 0.5,size = 1) + theme_linedraw(base_size = 14) +
  scale_color_manual(values = c("#ea5f94", "#e3e3e3"))  +
  labs(x = expression("log"[2]*(FC)*","~"Promoter"), y = expression("log"[10]*(p-value))) + 
  theme(legend.position = "none", plot.margin = margin(0,0,0,0,unit="cm")) +
  geom_hline(yintercept = -log10(max(wald_v2[[2]]$pval[wald_v2[[2]]$sig=="fdr_05"])), linetype = "dashed") +
  scale_y_continuous(limits = c(-0.1,8.5), expand = c(0,0))

ggsave(filename = "../../figs/wald_combined_promoter_vol.png", plot = plot1, units = "in", width = 7, height = 7)

dens <- wald_v2[[2]] %>% filter(sig == "fdr_05") %>% ggplot(aes(x = -logFC)) +
  geom_density(fill = "#ea5f94", alpha = 0.5) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0,unit="cm"))

p <- dens + plot1 + 
  plot_layout(
    ncol = 1, 
    nrow = 2, 
    heights = c(0.75, 4)
  )

ggsave(p, filename = "../../figs/wald_results_promoter_vol.png", units = "in",
       height = 7, width = 7)

wald_v2 <- wald_v3
wald_v2 <- map(wald_v2, ~ data.frame(.) %>%
                 mutate(ref_name = rownames(.)))
wald_plot_df <- merge(wald_v2[[1]], wald_v2[[2]], by = "ref_name")
wald_plot_df <- wald_plot_df %>%
  mutate(sig = as.factor(if_else(sig.x == "fdr_05" & sig.y == "fdr_05", "both",
                                 if_else(sig.x == "fdr_05" & sig.y != "fdr_05", "pos",
                                         if_else(sig.x != "fdr_05" & sig.y == "fdr_05", "prom",
                                                 "neither")))))
wald_plot_df$sig <- factor(wald_plot_df$sig, levels = c("pos", "prom", "both","neither"))
table(wald_plot_df$sig)
wald_plot_df <- wald_plot_df %>%
  mutate(alpha = if_else(sig == "neither", FALSE, TRUE))

plot1 <- wald_plot_df %>% ggplot(aes(x = logFC.x, y = logFC.y, color = sig, alpha = alpha)) +
  geom_point(size = 1) + theme_linedraw(base_size = 14) +
  scale_alpha_discrete(range = c(0.1, 0.5), breaks = c(0.1, 0.5)) +
  scale_color_manual(values = c("#ffb14e", "#ea5f94", "#0000ff", "#e3e3e3"), name = "Significant\nEffect",
                     labels = c("Position", "Promoter", "Both", "Neither")) +
  theme(legend.position = c(0.1, 0.875),
        legend.background = element_rect(color = "black", size = 0.25),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(3, units = "pt"),
        plot.margin = margin(0,0,0,0,unit="cm")) +
  labs(x = expression("log"[2]*(FC)*","~"Position"), y = expression("log"[2]*(FC)*","~"Promoter"))

dens1 <- wald_plot_df %>% filter(sig.x == "fdr_05") %>% ggplot(aes(x = logFC.x)) +
  geom_density(fill = "#ffb14e", alpha = 0.5) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0,unit="cm"))

dens2 <- wald_plot_df %>% filter(sig.y == "fdr_05") %>% ggplot(aes(x = logFC.y)) +
  geom_density(fill = "#ea5f94", alpha = 0.5) +
  theme_void() +
  coord_flip() +
  theme(plot.margin = margin(0,0,0,0,unit="cm"))

p <- dens1 + plot1 +
  plot_layout(
    ncol = 1, 
    nrow = 2, 
    heights = c(0.75, 4)
  )

ggsave(p, filename = "../../figs/wald_results_comparison_cage.png", units = "in",
       height = 7, width = 7)


wald_plot_df <- wald_plot_df %>%
  mutate(up = if_else(sig.x == "fdr_05" & logFC.x > 0, TRUE, FALSE),
         down = if_else(sig.x == "fdr_05" & logFC.x < 0, TRUE, FALSE),
         ins = if_else(sig.y == "fdr_05" & logFC.y < 0, TRUE, FALSE),
         scp = if_else(sig.y == "fdr_05" & logFC.y > 0, TRUE, FALSE))

wald_plot_df <- wald_plot_df %>%
  mutate(sig_any = if_else(sig.x == "fdr_05" | sig.y == "fdr_05", TRUE, FALSE))

wald_plot_df <- wald_plot_df %>%
  mutate(pos = case_when(sig.x == "fdr_05" & logFC.x > 0 ~ "up",
                         sig.x == "fdr_05" & logFC.x < 0 ~ "down",
                         TRUE ~ "NA"),
         prom = case_when(sig.y == "fdr_05" & logFC.y > 0 ~ "scp",
                          sig.y == "fdr_05" & logFC.y < 0 ~ "ins",
                          TRUE ~ "NA"))

count_df <- data.frame("variable"=c("down", "up","ins","scp"),"class"=c("Position","Position","Promoter","Promoter"),"count"=c(308,395,512,186))

count_df$variable <- as.factor(count_df$variable)
count_df$variable <- factor(count_df$variable, levels = c("down", "up", "scp", "ins"))

p <- count_df %>% group_by(class) %>% ggplot(aes(fill = variable, y = count, x = class)) +
  geom_bar(position = "fill", stat = "identity") + labs(x = NULL, y = "Proportion") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  scale_fill_manual(values = c("#ffb14e", "#D17D00", "#ea5f94", "#CA0068")) +
  geom_text(x = 2, y = (512/698)/2, label = "INS", fontface = "italic") +
  geom_text(x = 2, y = (186/698)/2 + (512/698), label = "SCP1", fontface = "italic") +
  geom_text(x = 1, y = (395/702)/2, label = "Down") +
  geom_text(x = 1, y = (395/702) + (308/702)/2, label = "Up") +
  scale_y_continuous(expand = c(0,0))
ggsave(p, filename = "../figs/wald_results_comparison_proportion.png", units = "in",
       height = 3, width = 3)
