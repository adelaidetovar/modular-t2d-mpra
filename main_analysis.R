#########################
### Main analysis ###
#########################

mpra_packages <- c("tidyverse", "factoextra", "data.table", "MPRAnalyze", "reshape2",
                   "ggtext", "UpSetR", "patchwork", "bedtoolsr", "glmnet")
lapply(mpra_packages, require, character.only=TRUE)

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

#############################
### Make annotation files ###
#############################

# read in Yoshi's annotations

oligo_annot <- read.delim(paste0(data_dir, 'oligo_annot.txt'))
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

write.table(reg_annot, paste0(out_dir, "oligo_annot_plus_regulatoryinfo.txt"), sep = "\t",
            row.names = F, quote = F)

############################
### Variant-trait figure ###
############################

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
  theme(legend.text = element_markdown(),
        text = element_text(family = "Helvetica"))
  

ggsave(plot = ecdf_plot, paste0(fig_dir, "ecdf.png"), units = "in", height = 5, width = 6)

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
                    text.scale = 1.75,
                    point.size = 2.8,
                    line.size = 1)

# change ratio of plot grid to prevent labels from cutting off
upset_plot$mb.ratio <- c(0.6, 0.4)

png(filename = paste0(fig_dir, "upset_plot.png"),
    width = 12.5, height = 5, units = "in", res = 300)
upset_plot
dev.off()

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

rna_annot[[1]]$prom <- rep("SCP1",length(rna_annot[[1]][1]))
rna_annot[[2]]$prom <- rep("SCP1",length(rna_annot[[2]][1]))
rna_annot[[3]]$prom <- rep("INS",length(rna_annot[[3]][1]))
rna_annot[[4]]$prom <- rep("INS",length(rna_annot[[4]][1]))

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

dna_annot[[1]]$prom <- rep("SCP1",length(dna_annot[[1]][1]))
dna_annot[[2]]$prom <- rep("SCP1",length(dna_annot[[2]][1]))
dna_annot[[3]]$prom <- rep("INS",length(dna_annot[[3]][1]))
dna_annot[[4]]$prom <- rep("INS",length(dna_annot[[4]][1]))

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

write.table(dna, paste0(out_dir, "combined.dna_counts.tsv"), sep="\t")
write.table(rna, paste0(out_dir, "combined.rna_counts.tsv"), sep="\t")
write.table(dna_annot, paste0(out_dir, "combined.dna_annots.tsv"), sep="\t")
write.table(rna_annot, paste0(out_dir, "combined.rna_annots.tsv"), sep="\t")

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
# change shading of points based on significant level
wald_plot_df <- wald_plot_df %>%
  mutate(alpha = if_else(sig == "neither", FALSE, TRUE))

# main scatter for wald plot
wald_plot <- wald_plot_df %>% ggplot(aes(x = logFC.y, y = logFC.x, color = sig, alpha = alpha)) +
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

wald_dens1 <- wald_plot_df %>% filter(sig.y == "fdr_05") %>% ggplot(aes(x = logFC.y)) +
  geom_density(fill = "#ea5f94", alpha = 0.5) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0,unit="cm"))

wald_dens2 <- wald_plot_df %>% filter(sig.x == "fdr_05") %>% ggplot(aes(x = logFC.x)) +
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
write.table(coef_parsed1, file=paste0(out_dir, "pos.lasso-zscore.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef1, file=paste0(out_dir, "pos.lasso-zscore.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef1, file=paste0(out_dir, "pos.lasso-zscore.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed2, file=paste0(out_dir, "pos.lasso-zscore.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef2, file=paste0(out_dir, "pos.lasso-zscore.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef2, file=paste0(out_dir, "pos.lasso-zscore.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed3, file=paste0(out_dir, "pos.lasso-zscore.red.coeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef3, file=paste0(out_dir, "pos.lasso-zscore.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef3, file=paste0(out_dir, "pos.lasso-zscore.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed4, file=paste0(out_dir, "prom.lasso-zscore.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef4, file=paste0(out_dir, "prom.lasso-zscore.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef4, file=paste0(out_dir, "prom.lasso-zscore.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed5, file=paste0(out_dir, "prom.lasso-zscore.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef5, file=paste0(out_dir, "prom.lasso-zscore.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef5, file=paste0(out_dir, "prom.lasso-zscore.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed6, file=paste0(out_dir, "prom.lasso-zscore.red.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef6, file=paste0(out_dir, "prom.lasso-zscore.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef6, file=paste0(out_dir, "prom.lasso-zscore.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed7, file=paste0(out_dir, "pos.lasso-zscore-negp.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef7, file=paste0(out_dir, "pos.lasso-zscore-negp.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef7, file=paste0(out_dir, "pos.lasso-zscore-negp.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed8, file=paste0(out_dir, "pos.lasso-zscore-negp.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef8, file=paste0(out_dir, "pos.lasso-zscore-negp.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef8, file=paste0(out_dir, "pos.lasso-zscore-negp.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed9, file=paste0(out_dir, "pos.lasso-zscore-negp.red.coeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef9, file=paste0(out_dir, "pos.lasso-zscore-negp.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef9, file=paste0(out_dir, "pos.lasso-zscore-negp.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed10, file=paste0(out_dir, "prom.lasso-zscore-negp.full.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef10, file=paste0(out_dir, "prom.lasso-zscore-negp.full.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef10, file=paste0(out_dir, "prom.lasso-zscore-negp.full.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed11, file=paste0(out_dir, "prom.lasso-zscore-negp.mid.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef11, file=paste0(out_dir, "prom.lasso-zscore-negp.mid.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef11, file=paste0(out_dir, "prom.lasso-zscore-negp.mid.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed12, file=paste0(out_dir, "prom.lasso-zscore-negp.red.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef12, file=paste0(out_dir, "prom.lasso-zscore-negp.red.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef12, file=paste0(out_dir, "prom.lasso-zscore-negp.red.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed13, file=paste0(out_dir, "pos.lasso-zscore.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef13, file=paste0(out_dir, "pos.lasso-zscore.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef13, file=paste0(out_dir, "pos.lasso-zscore.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed14, file=paste0(out_dir, "prom.lasso-zscore.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef14, file=paste0(out_dir, "prom.lasso-zscore.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef14, file=paste0(out_dir, "prom.lasso-zscore.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed15, file=paste0(out_dir, "pos.lasso-zscore-negp.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef15, file=paste0(out_dir, "pos.lasso-zscore-negp.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef15, file=paste0(out_dir, "pos.lasso-zscore-negp.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

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
write.table(coef_parsed16, file=paste0(out_dir, "prom.lasso-zscore-negp.redchrom.finalcoeff.tsv"), sep='\t', row.names=FALSE)
write.table(nonzero_coef16, file=paste0(out_dir, "prom.lasso-zscore-negp.redchrom.nonzerocoeff.tsv"), sep='\t', row.names=FALSE)
write.table(all_coef16, file=paste0(out_dir, "prom.lasso-zscore-negp.redchrom.allcoeff.tsv"), sep='\t', row.names=FALSE)

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


############################
### Second MPRA - 832/13 ###
############################

# read in dictionaries
adPath = paste0(data_dir, 'ADIPOQ.txt')
apPath = paste0(data_dir, 'APOA2.txt')
inPath = paste0(data_dir, 'INS.txt')
myPath = paste0(data_dir, 'MYBPC2.txt')
miPath = paste0(data_dir, 'minP.txt')
scPath = paste0(data_dir, 'SCP1.txt')

adSA = fread(adPath)
apSA = fread(apPath)
inSA = fread(inPath)
mySA = fread(myPath)
miSA = fread(miPath)
scSA = fread(scPath)

# reduce and remove duplicates
adSA <- adSA %>% distinct(readgroupid, refname)
apSA <- apSA %>% distinct(readgroupid, refname)
inSA <- inSA %>% distinct(readgroupid, refname)
mySA <- mySA %>% distinct(readgroupid, refname)
miSA <- miSA %>% distinct(readgroupid, refname)
scSA <- scSA %>% distinct(readgroupid, refname)

adSA$prom = "ADIPOQ"
apSA$prom = "APOA2"
inSA$prom = "INS"
mySA$prom = "MYBPC2"
miSA$prom = "minP"
scSA$prom = "SCP1"

subassembly = bind_rows(list(adSA, apSA, inSA, mySA, miSA, scSA))
rm(mySA)
rm(scSA)
rm(adSA)
rm(apSA)
rm(inSA)
rm(miSA)

# read in counts
countsPath = paste0(data_dir, '832_count_table.txt')

oligoCounts = fread(countsPath)

fullCounts = merge(oligoCounts, subassembly[,c(1,4,5)], by = "V1")
rm(oligoCounts)
rm(subassembly)
fullCounts$oligo_id <- substr(fullCounts$refname, 1, nchar(fullCounts$refname)-6)

oligonames <- read.delim(paste0(data_dir, "oligonames.tsv"))
oligonames$allele <- substr(oligonames$full_name, nchar(oligonames$full_name), nchar(oligonames$full_name))

formatCounts = merge(fullCounts, oligonames, by = "oligo_id")

colnames(formatCounts)[3:8] <- c("rna1", "rna2", "rna3", "rna4", "rna5", "dna")

sumCounts <- formatCounts %>%
  group_by(oligo_id, prom) %>%
  summarize(oligo_id = oligo_id,
            rna1 = sum(rna1),
            rna2 = sum(rna2),
            rna3 = sum(rna3),
            rna4 = sum(rna4),
            rna5 = sum(rna5),
            dna = sum(dna),
            base_name = base_name,
            prom = prom) %>%
  distinct(base_name, .keep_all = TRUE)

# filter for rows where all have at least 5 counts
filtSumCounts <- sumCounts %>% filter(rna1 > 5 & rna2 > 5 & rna3 > 5 & rna4 > 5 & rna5 > 5 & dna > 5)

# make design matrix
design = data.frame("material" = c(rep("RNA",5),rep("DNA",1)))
design$replicate = c("cDNA1","cDNA2","cDNA3","cDNA4","cDNA5","plasmid1")
rownames(design) = names(filtSumCounts)[3:8]

# build and execute DESeq2 call and retrieve summary statistics
deseq_object = DESeqDataSetFromMatrix(countData=filtSumCounts[3:8],
                                      colData = design,
                                      design = ~material)
deseq_object =  DESeq(deseq_object, fitType = 'local', minReplicatesForReplace=Inf)

# format contrasts, then bind into final dataframe
results_material = results(deseq_object, contrast = c("material","RNA","DNA"),
                           cooksCutoff = FALSE, independentFiltering = FALSE)
names(results_material) = paste0(names(results_material),"_", "expr")

results_DESeq2 = cbind(filtSumCounts, as.data.frame(results_material))

write.table(results_DESeq2, paste0(out_dir, "832-lib2-expr-filt.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

# now perform allelic analysis with linear contrast
design_contrast = data.frame("material" = factor(rep(c(rep("RNA",5),rep("DNA",1)),2)),
                             "allele" = factor(c(rep("ref",6),rep("alt",6))),
                             "sample" = factor(rep(c("cDNA1","cDNA2","cDNA3","cDNA4","cDNA5","plasmid1"),2)),
                             "dnaAllele" = c(0,0,0,0,0,0,0,0,0,0,0,1),
                             "rnaAllele" = c(0,0,0,0,0,0,1,1,1,1,1,0))
rownames(design_contrast) = names(allelicCounts)[3:14]

deseq_object_contrast = DESeqDataSetFromMatrix(countData=allelicCounts[,3:14],
                                               colData = design_contrast,
                                               design = ~0 + sample + dnaAllele + rnaAllele)
deseq_object_contrast = DESeq(deseq_object_contrast)

results_expr_contrast = results(deseq_object_contrast, contrast = list(c("samplecDNA1", "samplecDNA2","samplecDNA3",
                                                                         "samplecDNA4","samplecDNA5"),
                                                                       c("sampleplasmid1")))
results_allele_contrast = results(deseq_object_contrast, contrast = list("dnaAllele","rnaAllele"))

names(results_expr_contrast) = paste0(names(results_expr_contrast),"_", "expr")
names(results_allele_contrast) = paste0(names(results_allele_contrast),"_", "allele")

results_DESeq2_contrast = cbind(allelicCounts,
                                as.data.frame(results_expr_contrast),
                                as.data.frame(results_allele_contrast))

write.table(out_dir, paste0(out_dir, "832-lib2-allelic-filt.tsv"), quote = F, sep = '\t', row.names = F, col.names = T)

#############################
### Second MPRA - LHCN-M2 ###
#############################

# read in dictionaries
inPath = paste0(data_dir, 'INS.txt')
myPath = paste0(data_dir, 'MYBPC2.txt')
scPath = paste0(data_dir, 'SCP1.txt')

inSA = fread(inPath)
mySA = fread(myPath)
scSA = fread(scPath)

# reduce and remove duplicates
inSA <- inSA %>% distinct(readgroupid, refname)
mySA <- mySA %>% distinct(readgroupid, refname)
scSA <- scSA %>% distinct(readgroupid, refname)

inSA$prom = "INS"
mySA$prom = "MYBPC2"
scSA$prom = "SCP1"

subassembly = bind_rows(list(inSA, mySA, scSA))
rm(mySA)
rm(scSA)
rm(inSA)

# read in counts
countsPath = paste0(data_dir, 'lhcn_count_table.txt')

oligoCounts = fread(countsPath)

fullCounts = merge(oligoCounts[,c(1:5,18:21)], subassembly[,c(1,4,5)], by = "V1")
rm(oligoCounts)
rm(subassembly)
fullCounts$oligo_id <- substr(fullCounts$refname, 1, nchar(fullCounts$refname)-6)

oligonames <- read.delim(paste0(data_dir, "oligonames.tsv"))
oligonames$allele <- substr(oligonames$full_name, nchar(oligonames$full_name), nchar(oligonames$full_name))

formatCounts = merge(fullCounts, oligonames, by = "oligo_id")

colnames(formatCounts)[3:10] <- c("rna1", "rna2", "rna3", "rna4", "dna1", "dna2", "dna3", "dna4")

sumCounts <- formatCounts %>%
  group_by(oligo_id, prom) %>%
  summarize(oligo_id = oligo_id,
            rna1 = sum(rna1),
            rna2 = sum(rna2),
            rna3 = sum(rna3),
            rna4 = sum(rna4),
            dna1 = sum(dna1),
            dna2 = sum(dna2),
            dna3 = sum(dna3),
            dna4 = sum(dna4),
            base_name = base_name,
            prom = prom) %>%
  distinct(base_name, .keep_all = TRUE)

# filter for rows where all have at least 5 counts
filtSumCounts <- sumCounts %>% filter(rna1 > 5 & rna2 > 5 & rna3 > 5 & rna4 > 5 &
                                        dna1 > 5 & dna2 > 5 & dna3 > 5 & dna4 > 5)

# make design matrix
design = data.frame("material" = c(rep("RNA",4),rep("DNA",4)))
design$replicate = c("cDNA1","cDNA2","cDNA3","cDNA4","gDNA1","gDNA2","gDNA3","gDNA4")
rownames(design) = names(filtSumCounts)[3:10]

# build and execute DESeq2 call and retrieve summary statistics
deseq_object = DESeqDataSetFromMatrix(countData=filtSumCounts[3:10],
                                      colData = design,
                                      design = ~material)
deseq_object =  DESeq(deseq_object, fitType = 'local', minReplicatesForReplace=Inf)

# format contrasts, then bind into final dataframe
results_material = results(deseq_object, contrast = c("material","RNA","DNA"),
                           cooksCutoff = FALSE, independentFiltering = FALSE)
names(results_material) = paste0(names(results_material),"_", "expr")

results_DESeq2 = cbind(filtSumCounts, as.data.frame(results_material))

write.table(results_DESeq2, paste0(out_dir, "lhcn-lib2-expr-filt.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

# now perform allelic analysis with linear contrast
design_contrast = data.frame("material" = factor(rep(c(rep("RNA",5),rep("DNA",1)),2)),
                             "allele" = factor(c(rep("ref",6),rep("alt",6))),
                             "sample" = factor(rep(c("cDNA1","cDNA2","cDNA3","cDNA4","cDNA5","plasmid1"),2)),
                             "dnaAllele" = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1),
                             "rnaAllele" = c(0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0))
rownames(design_contrast) = names(allelicCounts)[3:18]

deseq_object_contrast = DESeqDataSetFromMatrix(countData=allelicCounts[,3:18],
                                               colData = design_contrast,
                                               design = ~0 + sample + dnaAllele + rnaAllele)
deseq_object_contrast = DESeq(deseq_object_contrast)

results_expr_contrast = results(deseq_object_contrast, contrast = list(c("samplecDNA1", "samplecDNA2","samplecDNA3",
                                                                         "samplecDNA4"),
                                                                       c("samplegDNA1", "samplegDNA2","samplegDNA3",
                                                                         "samplegDNA4")))
results_allele_contrast = results(deseq_object_contrast, contrast = list("dnaAllele","rnaAllele"))

names(results_expr_contrast) = paste0(names(results_expr_contrast),"_", "expr")
names(results_allele_contrast) = paste0(names(results_allele_contrast),"_", "allele")

results_DESeq2_contrast = cbind(allelicCounts,
                                as.data.frame(results_expr_contrast),
                                as.data.frame(results_allele_contrast))

write.table(out_dir, paste0(out_dir, "lhcn-lib2-allelic-filt.tsv"), quote = F, sep = '\t', row.names = F, col.names = T)
