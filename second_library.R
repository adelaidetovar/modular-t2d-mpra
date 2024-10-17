###############################
### Second library analysis ###
###############################

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

############################
### Second MPRA - 832/13 ###
############################

# read in dictionaries
dict_path = file.path('/path/to/tovar_library_barcode_pairing.txt')

dict = fread(dict_path)
# fix neg control barcodes b/c they're reverse complements
neg_oligos <- dict %>%
  filter(str_detect(refname, "^NA"))
barcodes <- Biostrings::DNAStringSet(neg_oligos$barcode)
barcodes <- Biostrings::reverseComplement(barcodes)
barcodes <- as.character(barcodes)
neg_oligos$barcode <- barcodes
dict <- dict %>%
  filter(!str_detect(refname, "^NA"))
dict <- bind_rows(dict, neg_oligos)

cts_path = file.path('/path/to/tovar_library_counts.txt')
raw_cts = fread(cts_path)

merged_raw_cts = merge(dict, raw_cts, by = "barcode")

filtered_cts = merged_raw_cts %>%
  filter(rep1 + rep2 + rep3 + rep4 + rep5 > 0) %>%
  filter(dna > 2)

summ_filt_cts = filtered_cts %>%
  group_by(refname, config) %>%
  reframe(refname = refname,
          config = config,
          across(dna:rep5, sum)) %>%
  distinct(refname, config, .keep_all = TRUE)

# format counts for MPRAnalyze
reform_filt_cts = filtered_cts %>%
  group_by(refname, config) %>%
  mutate(bcnum = 1:n()) %>%
  ungroup() %>%
  mutate(bcnum = paste0("bc",bcnum))

gather_filt_cts = reform_filt_cts %>% gather(rep, value, -c(refname, bcnum, barcode, config)) %>%
  mutate(refname_full = paste(config, refname, sep = "_"),
         bcnum_rep = paste(bcnum, rep, sep = "_"))

dna_filt_cts_pivot = gather_filt_cts %>%
  filter(rep == "dna") %>%
  select(refname_full, bcnum, value) %>%
  pivot_wider(names_from = bcnum, values_from = value, values_fill = 0)

rna_filt_cts_pivot = gather_filt_cts %>%
  filter(rep != "dna") %>%
  select(refname_full, bcnum_rep, value) %>%
  pivot_wider(names_from = bcnum_rep, values_from = value, values_fill = 0)

dna_filt_annot = data.frame("bcnum" = colnames(dna_filt_cts_pivot)[-1])
dna_filt_annot$barcode = dna_filt_annot$bcnum
dna_filt_annot$lib = "dna"

rna_filt_annot = data.frame("bcnum_rep" = colnames(rna_filt_cts_pivot)[-1])
rna_filt_annot <- rna_filt_annot %>%
  separate(bcnum_rep, into = c("barcode", "rep"), remove = FALSE, sep = "_")

### run MPRAnalyze
dna = as.data.frame(dna_filt_cts_pivot)
rownames(dna) = dna$refname_full
dna = as.matrix(dna[,-1])

rna = as.data.frame(rna_filt_cts_pivot)
rownames(rna) = rna$refname_full
rna = as.matrix(rna[,-1])

dna_annot = dna_filt_annot
rna_annot = rna_filt_annot

obj_test = MpraObject(dnaCounts = dna, rnaCounts = rna, dnaAnnot = dna_annot, rnaAnnot = rna_annot)

obj_test <- estimateDepthFactors(obj_test, which.lib = "dna",
                                 depth.estimator = "uq")
obj_test <- estimateDepthFactors(obj_test, lib.factor = "rep",
                                 which.lib = "rna",
                                 depth.estimator = "uq")

obj_test = analyzeQuantification(obj = obj_test, dnaDesign = ~ barcode, rnaDesign = ~rep)
output <- testEmpirical(obj_test)
alpha_out <- getAlpha(obj_test)
output$alpha_orig <- alpha_out$alpha
output$alpha_log <- log2(exp(output$alpha_orig))
output <- output %>%
  rownames_to_column(var = "refname_full") %>%
  separate(refname_full, into = c("prom", "rsid", "chr", "pos_hg38", "ref", "alt", "allele", "type", "site"),
           sep = "_", remove = FALSE)

output <- output %>%
  mutate(site = case_when(type == "left" ~ "left",
                          type == "right" ~ "right",
                          TRUE ~ site),
         type = case_when(is.na(type) | type == "left" | type == "right" ~ "orig",
                          TRUE ~ type),
         site = case_when(is.na(site) ~ "ctr",
                          TRUE ~ site))

write.table(output, paste0(out_dir, "tovar_library.res.tsv"), sep='\t', row.names=FALSE)

out_filt <- output %>%
  dplyr::filter(pval.zscore < 0.05)

neg_filt <- output %>%
  dplyr::filter(pval.zscore > 0.2 & rsid == "NA")

rsid_filt <- out_filt %>%
  dplyr::filter(rsid != "NA") %>%
  distinct(rsid, site, prom, .keep_all = TRUE) %>%
  mutate(min_name = paste(rsid, site, prom, sep = "_"))
  
summ_filt_cts <- summ_filt_cts %>%
  separate(refname, into = c("rsid", "chr", "pos_hg38", "ref", "alt", "allele", "type", "site"),
           sep = "_", remove = FALSE) %>%
  mutate(site = case_when(type == "left" ~ "left",
                          type == "right" ~ "right",
                          TRUE ~ site),
         type = case_when(is.na(type) | type == "left" | type == "right" ~ "orig",
                          TRUE ~ type),
         site = case_when(is.na(site) ~ "ctr",
                          TRUE ~ site))

summ_filt_pairs <- summ_filt_cts %>%
  mutate(min_name = paste(rsid, site, config, sep = "_")) %>%
  filter(min_name %in% unique(rsid_filt$min_name))

summ_filt_pairs$category <- "test_vars"
neg_filt_pairs <- summ_filt_cts %>%
  mutate(min_name = paste(rsid, site, config, sep = "_")) %>%
  filter(paste(config, refname, sep = "_") %in% neg_filt$refname_full) %>%
  mutate(category = "neg_ctrls")

summ_filt_pairs <- bind_rows(summ_filt_pairs, neg_filt_pairs)

tpm_filt_pairs <- summ_filt_pairs %>%
  pivot_longer(starts_with(c("rep","dna")), 
               names_to = "sample", 
               values_to = "counts") %>%
  ungroup() %>%
  select(refname, rsid, chr, pos_hg38, ref, alt, allele, type, site, config, min_name, sample, counts) %>%
  group_by(sample) %>%
  mutate(tpm = (counts / (sum(counts)) * 1e6)) %>%
  select(refname, rsid, chr, pos_hg38, ref, alt, allele, type, site, config, min_name, sample, tpm) %>%
  pivot_wider(names_from = sample, values_from = tpm)

tpm_filt_pairs <- tpm_filt_pairs %>%
  mutate(ratio1 = log2((rep1/dna)+1),
         ratio2 = log2((rep2/dna)+1),
         ratio3 = log2((rep3/dna)+1),
         ratio4 = log2((rep4/dna)+1),
         ratio5 = log2((rep5/dna)+1))

ins_filt_pairs <- tpm_filt_pairs %>%
  filter(config == "ins")

scp_filt_cts <- summ_filt_cts %>%
  filter(config == "scp")

scp_filt_tpm <- scp_filt_cts %>%
  pivot_longer(starts_with(c("rep","dna")), 
               names_to = "sample", 
               values_to = "counts") %>%
  ungroup() %>%
  select(refname, rsid, chr, pos_hg38, ref, alt, allele, type, site, config, sample, counts) %>%
  group_by(sample) %>%
  mutate(tpm = (counts / (sum(counts)) * 1e6)) %>%
  select(refname, rsid, chr, pos_hg38, ref, alt, allele, type, site, config, sample, tpm) %>%
  pivot_wider(names_from = sample, values_from = tpm)

scp_filt_tpm <- scp_filt_tpm %>%
  mutate(ratio1 = log2((rep1/dna)+1),
         ratio2 = log2((rep2/dna)+1),
         ratio3 = log2((rep3/dna)+1),
         ratio4 = log2((rep4/dna)+1),
         ratio5 = log2((rep5/dna)+1))

ins_filt_melt <- melt(ins_filt_pairs,
                      id.vars = c("refname", "config", "allele", "type",
                                  "site"),
                      measure.vars = c("ratio1", "ratio2", "ratio3",
                                       "ratio4", "ratio5"))

scp_filt_melt <- melt(scp_filt_tpm,
                      id.vars = c("refname", "config", "allele", "type",
                                  "site"),
                      measure.vars = c("ratio1", "ratio2", "ratio3",
                                       "ratio4", "ratio5"))


neg_controls <- c("NA_chr19_3586640_NA_NA_NA",
                  "NA_chr5_117093713_NA_NA_NA",
                  "NA_chr11_116285064_NA_NA_NA")

neg_filt_melt <- ins_filt_melt %>%
  filter(refname %in% neg_controls) %>%
  group_by(variable) %>%
  reframe(refname = refname,
          config = config,
          allele = allele,
          type = type,
          site = site,
          variable = variable,
          value = mean(value)) %>%
  distinct(variable, .keep_all = TRUE)

neg_scp_filt_melt <- scp_filt_melt %>%
  filter(refname %in% neg_controls) %>%
  group_by(variable) %>%
  reframe(refname = refname,
          config = config,
          allele = allele,
          type = type,
          site = site,
          variable = variable,
          value = mean(value)) %>%
  distinct(variable, .keep_all = TRUE)
  

rs163_scp_melt <- scp_filt_melt %>%
  filter(refname %in% c("rs1635852_chr7_28149792_T_C_R_del_left", "rs1635852_chr7_28149792_T_C_R_left",
                        "rs1635852_chr7_28149792_T_C_R_shf_left","rs1635852_chr7_28149792_T_C_A_left"))

rs163_scp_melt <- left_join(rs163_scp_melt, neg_scp_filt_melt[,c(1,7)], by = "variable") %>%
  mutate(plot_value = value.x - value.y)

rs163_scp_melt <- rs163_scp_melt %>%
  mutate(allele = case_when(type %in% c("del","shf") ~ "",
                            TRUE ~ allele),
         allele = factor(allele, levels = c("R","","A")),
         type = factor(type, levels = c("orig","del","shf")))

# rs1635852
rs163_melt <- ins_filt_melt %>%
  filter(refname %in% c("rs1635852_chr7_28149792_T_C_R_del_left", "rs1635852_chr7_28149792_T_C_R_left",
                        "rs1635852_chr7_28149792_T_C_R_shf_left","rs1635852_chr7_28149792_T_C_A_left"))

rs163_melt <- left_join(rs163_melt, neg_filt_melt[,c(1,7)], by = "variable") %>%
  mutate(plot_value = value.x - value.y)
rs163_melt <- rs163_melt %>%
  mutate(allele = case_when(type %in% c("del","shf") ~ "",
                            TRUE ~ allele),
         allele = factor(allele, levels = c("R","","A")),
         type = factor(type, levels = c("orig","del","shf")))

rs163_melt %>% rstatix::wilcox_test(value.x ~ refname, exact = TRUE)

rs163_scp_melt %>%
  ggplot(aes(x = allele, y = plot_value, color = interaction(type,allele))) +
  geom_boxplot(aes(linetype = type, fill = interaction(type,allele)),
               outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_point(aes(group = interaction(allele,type)),
             pch = 21, size = 2.5, color = "black",
             fill = "white", position = position_dodge(width = 0.75),
             stroke = 0.75) +
  theme_bw(base_size = 15) + labs(x = "Allele", y = expression(paste("log"[2], "(RNA/DNA), rs1635852"))) +
  scale_color_manual(values = c("#ffffff","#000000","#000000",
                                "#ffffff","#0000ff", "#0000ff")) +
  scale_fill_manual(values = c("#ff0000", "#ffffff", "#ffffff",
                               "#0000ff", "#ffffff", "#ffffff")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"),
                        labels = c("Original", "Deleted", "Shuffled")) +
  scale_x_discrete(labels = c("ref" = "Ref\n(Risk, T)", "none" = "", "alt" = "Alt\n(Non-risk, C)")) +
  guides(linetype = guide_legend(
    override.aes = (list(shape = 21, fill = c("black", "white", "white"),
                         linetype = c("solid", "solid", "dashed"),
                         color = c("white", "black", "black"))),
    title = "Oligo Type"
  ),
  color = "none",
  fill = "none") +
  theme(legend.position = c(0.8, 0.15),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="black"),
        text=element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "#818181")

ggsave(filename = paste0(fig_dir, "rs1635852_plot.png"), plot = rs118_plot,
       dpi = 600, device = ragg::agg_png(), units = "in", width = 5, height = 5)

# rs11819995
rs118_scp_melt <- scp_filt_melt %>%
  filter(refname %in% c("rs11819995_chr11_128519496_C_T_A_right", "rs11819995_chr11_128519496_C_T_A_shf_right",
                        "rs11819995_chr11_128519496_C_T_A_del_right","rs11819995_chr11_128519496_C_T_R_right",
                        "rs11819995_chr11_128519496_C_T_R_shf_right","rs11819995_chr11_128519496_C_T_R_del_right"))
  

rs118_melt <- ins_filt_melt %>%
  filter(refname %in% c("rs11819995_chr11_128519496_C_T_A_right", "rs11819995_chr11_128519496_C_T_A_shf_right",
                        "rs11819995_chr11_128519496_C_T_A_del_right","rs11819995_chr11_128519496_C_T_R_right",
                        "rs11819995_chr11_128519496_C_T_R_shf_right","rs11819995_chr11_128519496_C_T_R_del_right"))

rs118_scp_melt <- left_join(rs118_scp_melt, neg_scp_filt_melt[,c(1,7)], by = "variable") %>%
  mutate(plot_value = value.x - value.y)
rs118_scp_melt <- rs118_scp_melt %>%
  mutate(allele = factor(allele, levels = c("R","A")),
         type = factor(type, levels = c("orig","del","shf")))

rs118_scp_melt %>% rstatix::wilcox_test(value.x ~ refname, exact = TRUE)


rs118_melt <- left_join(rs118_melt, neg_filt_melt[,c(1,7)], by = "variable") %>%
  mutate(plot_value = value.x - value.y)
rs118_melt <- rs118_melt %>%
  mutate(allele = factor(allele, levels = c("R","A")),
         type = factor(type, levels = c("orig","del","shf")))

rs118_melt %>% rstatix::wilcox_test(value.x ~ refname, exact = TRUE)

rs118_scp_melt %>%
  ggplot(aes(x = allele, y = plot_value, color = interaction(type,allele))) +
  geom_boxplot(aes(linetype = type, fill = interaction(type,allele)),
               outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_point(aes(group = interaction(allele,type)),
             pch = 21, size = 2.5, color = "black",
             fill = "white", position = position_dodge(width = 0.75),
             stroke = 0.75) +
  theme_bw(base_size = 15) + labs(x = "Allele", y = expression(paste("log"[2],"(RNA/DNA), rs11819995"))) +
  scale_color_manual(values = c("#ffffff","#ff0000","#ff0000",
                                "#ffffff","#0000ff", "#0000ff")) +
  scale_fill_manual(values = c("#ff0000", "#ffffff", "#ffffff",
                               "#0000ff", "#ffffff", "#ffffff")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"),
                        labels = c("Original", "Deleted", "Shuffled")) +
  scale_x_discrete(labels = c("ref" = "Ref\n(Risk, T)", "none" = "", "alt" = "Alt\n(Non-risk, C)")) +
  guides(linetype = guide_legend(
    override.aes = (list(shape = 21, fill = c("black", "white", "white"),
                         linetype = c("solid", "solid", "dashed"),
                         color = c("white", "black", "black"))),
    title = "Oligo Type"
  ),
  color = "none",
  fill = "none") +
  theme(legend.position = c(0.825, 0.17),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="black"),
        text=element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "#818181") +
  scale_y_continuous(limits = c(-2.15,3.1))

ggsave(filename = paste0(fig_dir, "rs11819995_plot.png"), plot = rs118_plot,
       dpi = 600, device = ragg::agg_png(), units = "in", width = 5, height = 5)