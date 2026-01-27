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

cts_path = file.path('/path/to/tovar_library_counts.txt')
raw_cts = fread(cts_path)

raw_tpm = raw_cts %>%
  pivot_longer(starts_with(c("rep","dna")),
               names_to = "sample",
               values_to = "counts") %>%
  group_by(sample) %>%
  mutate(tpm = (counts / (sum(counts)) * 1e6)) %>%
  ungroup() %>%
  select(-counts) %>%
  pivot_wider(names_from = sample, values_from = tpm)

merged_tpm = merge(raw_tpm, dict, by = "barcode") %>%
  group_by(refname, config) %>%
  reframe(refname = refname,
          config = config,
          across(c(rep1:dna), sum)) %>%
  distinct(refname, config, .keep_all = TRUE)

merged_ratio = merged_tpm %>%
  separate(refname, into = c("rsid", "chr", "pos_hg38", "ref", "alt", "allele", "type", "site"),
           sep = "_", remove = FALSE) %>%
  mutate(site = case_when(type == "left" ~ "left",
                          type == "right" ~ "right",
                          TRUE ~ site),
         type = case_when(is.na(type) | type == "left" | type == "right" ~ "orig",
                          TRUE ~ type),
         site = case_when(is.na(site) ~ "ctr",
                          TRUE ~ site)) %>%
  mutate(ratio1 = log2((rep1/dna)),
         ratio2 = log2((rep2/dna)),
         ratio3 = log2((rep3/dna)),
         ratio4 = log2((rep4/dna)),
         ratio5 = log2((rep5/dna)),
         ratio1 = case_when(is.infinite(ratio1) ~ 0,
                            TRUE ~ ratio1),
         ratio2 = case_when(is.infinite(ratio2) ~ 0,
                            TRUE ~ ratio2),
         ratio3 = case_when(is.infinite(ratio3) ~ 0,
                            TRUE ~ ratio3),
         ratio4 = case_when(is.infinite(ratio4) ~ 0,
                            TRUE ~ ratio4),
         ratio5 = case_when(is.infinite(ratio5) ~ 0,
                            TRUE ~ ratio5))

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

summ_filt_tpm = summ_filt_cts %>%
  separate(refname, into = c("rsid", "chr", "pos_hg38", "ref", "alt", "allele", "type", "site"),
           sep = "_", remove = FALSE) %>%
  mutate(site = case_when(type == "left" ~ "left",
                          type == "right" ~ "right",
                          TRUE ~ site),
         type = case_when(is.na(type) | type == "left" | type == "right" ~ "orig",
                          TRUE ~ type),
         site = case_when(is.na(site) ~ "ctr",
                          TRUE ~ site)) %>%
  pivot_longer(starts_with(c("rep","dna")), 
               names_to = "sample", 
               values_to = "counts") %>%
  ungroup() %>%
  select(refname, rsid, chr, pos_hg38, ref, alt, allele, type, site, config, sample, counts) %>%
  group_by(sample) %>%
  mutate(tpm = (counts / (sum(counts)) * 1e6)) %>%
  select(refname, rsid, chr, pos_hg38, ref, alt, allele, type, site, config, sample, tpm) %>%
  pivot_wider(names_from = sample, values_from = tpm)

nonzero_oligos <- summ_filt_cts %>%
  filter(rowSums(across(rep1:rep5, ~ .x > 0)) >= 3) %>%
  mutate(refname_full = paste(config, refname, sep = "_"))

neg_control_oligos <- nonzero_oligos %>%
  filter(substr(refname_full, 5, 5) == "N") %>%
  mutate(ratio1 = log2(rep1/dna),
         ratio2 = log2(rep2/dna),
         ratio3 = log2(rep3/dna),
         ratio4 = log2(rep4/dna),
         ratio5 = log2(rep5/dna),
         avg_ratio = (ratio1 + ratio2 + ratio3 + ratio4 + ratio5)/5) %>%
  filter(avg_ratio >= -1 & avg_ratio <= 0.5)

# format counts for MPRAnalyze
reform_filt_cts = filtered_cts %>%
  mutate(refname_full = paste(config, refname, sep = "_")) %>%
  filter(refname_full %in% nonzero_oligos$refname_full) %>%
  select(-c(refname_full)) %>%
  group_by(refname, config) %>%
  mutate(bcnum = 1:n()) %>%
  ungroup() %>%
  mutate(bcnum = paste0("bc",bcnum))

gather_filt_cts = reform_filt_cts %>%
  gather(rep, value, -c(refname, bcnum, barcode, config)) %>%
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

obj_test = MpraObject(dnaCounts = dna, rnaCounts = rna, dnaAnnot = dna_annot, rnaAnnot = rna_annot, controls = neg_control_df$refname_full)

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
  dplyr::filter(pval.mad < 0.05) %>%
  mutate(min_name = paste(rsid, site, prom, sep = "_"))

rsid_filt <- out_filt %>%
  dplyr::filter(rsid != "NA") %>%
  distinct(rsid, site, prom, .keep_all = TRUE) %>%
  mutate(min_name = paste(rsid, site, prom, sep = "_"))

tpm_filt_pairs <- ratio_merge %>%
  mutate(refname_full = paste(config, refname, sep = "_")) %>%
  filter(refname_full %in% output$refname_full) %>%
  mutate(min_name = paste(rsid, site, config, sep = "_")) %>%
  filter(min_name %in% out_filt$min_name)

ins_filt_pairs <- tpm_filt_pairs %>%
  filter(config == "ins")

ins_filt_long <- ins_filt_pairs %>%
  filter(rsid != "NA") %>%
  pivot_longer(cols = matches("ratio"),
               names_to = "replicate",
               values_to = "activity") %>%
  mutate(allele_type = paste(allele, type, sep = "_"))

## test contrasts for active oligo sets

test_contrasts <- function(data) {
  contrasts <- list(
    ref_vs_alt = list(ref = "R_orig", comparison = "A_orig"),
    ref_vs_ref_del = list(ref = "R_orig", comparison = "R_del"),
    ref_vs_ref_shf = list(ref = "R_orig", comparison = "R_shf"),
    ref_vs_alt_del = list(ref = "R_orig", comparison = "A_del"),
    ref_vs_alt_shf = list(ref = "R_orig", comparison = "A_shf"),
    alt_vs_alt_del = list(ref = "A_orig", comparison = "A_del"),
    alt_vs_alt_shf = list(ref = "A_orig", comparison = "A_shf"),
    alt_vs_ref_del = list(ref = "A_orig", comparison = "R_del"),
    alt_vs_ref_shf = list(ref = "A_orig", comparison = "R_shf")
  )
  
  results <- map_dfr(names(contrasts), function(contrast_name) {
    ref_group <- contrasts[[contrast_name]]$ref
    comp_group <- contrasts[[contrast_name]]$comparison
    
    ref_data <- data %>%
      filter(allele_type == ref_group) %>%
      pull(activity)
    
    comp_data <- data %>%
      filter(allele_type == comp_group) %>%
      pull(activity)
    
    if(length(ref_data) == 0 | length(comp_data) == 0){
      return(tibble(contrast = contrast_name, p_value = NA,
                    estimate = NA, mean_intact = NA, mean_perturbed = NA))
    }
    
    tryCatch({
      test_result <- wilcox.test(comp_data, ref_data, exact = FALSE)
      
      tibble(
        contrast = contrast_name,
        p_value = test_result$p.value,
        estimate = mean(comp_data, na.rm = TRUE) - mean(ref_data, na.rm = TRUE),
        mean_intact = mean(ref_data, na.rm = TRUE),
        mean_perturbed = mean(comp_data, na.rm = TRUE)
      )
    }, error = function(e) {
      tibble(contrast = contrast_name, p_value = NA, estimate = NA,
             mean_intact = NA, mean_perturbed = NA)
    })
  })
  
  return(results)
}

contrast_results <- ins_filt_long %>%
  group_by(min_name) %>%
  nest() %>%
  mutate(tests = map(data, test_contrasts)) %>%
  select(-data) %>%
  unnest(tests)

significant_sets <- contrast_results %>%
  filter(!is.na(p_value)) %>%
  group_by(min_name) %>%
  summarize(
    any_sig = any(p_value < 0.05, na.rm = TRUE),
    n_sig = sum(p_value < 0.05, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    sig_contrast = paste(contrast[p_value < 0.05], collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(any_sig)

## plotting

neg_control_df <- summ_filt_tpm %>%
  filter(rowSums(across(rep1:rep5, ~ .x > 0)) >= 5) %>%
  mutate(refname_full = paste(config, refname, sep = "_")) %>%
  filter(substr(refname_full, 5, 5) == "N") %>%
  mutate(ratio1 = log2((rep1/dna)),
         ratio2 = log2((rep2/dna)),
         ratio3 = log2((rep3/dna)),
         ratio4 = log2((rep4/dna)),
         ratio5 = log2((rep5/dna)),
         ratio1 = case_when(is.infinite(ratio1) ~ 0,
                            TRUE ~ ratio1),
         ratio2 = case_when(is.infinite(ratio2) ~ 0,
                            TRUE ~ ratio2),
         ratio3 = case_when(is.infinite(ratio3) ~ 0,
                            TRUE ~ ratio3),
         ratio4 = case_when(is.infinite(ratio4) ~ 0,
                            TRUE ~ ratio4),
         ratio5 = case_when(is.infinite(ratio5) ~ 0,
                            TRUE ~ ratio5),
         avg_ratio = (ratio1+ ratio2 + ratio3 + ratio4 + ratio5)/5) %>%
  filter(avg_ratio >= 0 & avg_ratio <= 1) %>%
  group_by(config) %>%
  reframe(ratio1 = median(ratio1),
          ratio2 = median(ratio2),
          ratio3 = median(ratio3),
          ratio4 = median(ratio4),
          ratio5 = median(ratio5)) %>%
  mutate(refname = "neg_ctrl") %>%
  pivot_longer(!c(config,refname))

ratio_plots <- ratio_merge %>%
  select(c(refname:config, ratio1:ratio5)) %>%
  pivot_longer(!c(refname:config)) %>%
  left_join(., neg_control_df[,c(1,3,4)], by = c("name", "config")) %>%
  mutate(norm_value = value.x - value.y) %>%
  select(c(refname:config, name, norm_value)) %>%
  pivot_wider(values_from = norm_value) %>%
  mutate(refname_full = paste(config, refname, sep = "_"))

rs163_ins_df <- ratio_plots %>%
  select(c(refname:config, ratio1:ratio5)) %>%
  filter(rsid == "rs1635852" & site == "left" & config == "ins") %>%
  pivot_longer(!c(refname:config)) %>%
  mutate(allele = case_when(type %in% c("del", "shf") ~ "",
                            TRUE ~ allele),
         allele = factor(allele, levels = c("R", "", "A")),
         type = factor(type, levels = c("orig", "del", "shf")))

rs163_ins_res <- rs163_ins_df %>%
  rstatix::wilcox_test(value ~ refname)
  
rs163_ins_plot <- rs163_ins_df %>%
  ggplot(aes(x = allele, y = value, color = interaction(type,allele))) +
  geom_boxplot(aes(linetype = type, fill = interaction(type,allele)),
               outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_point(aes(group = interaction(allele,type)),
             pch = 21, size = 3, color = "black",
             fill = "white", position = position_dodge(width = 0.75),
             stroke = 1) +
  theme_bw(base_size = 15) + labs(x = "Allele", y = expression(paste("log"[2], "(RNA/DNA), rs1635852"))) +
  scale_color_manual(values = c("#ffffff","#000000","#000000",
                                "#ffffff","#0000ff", "#0000ff")) +
  scale_fill_manual(values = c("#ff0000", "#ffffff", "#ffffff",
                               "#0000ff", "#ffffff", "#ffffff")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"),
                        labels = c("Original", "Deleted", "Shuffled")) +
  scale_x_discrete(labels = c("R" = "Ref\n(Risk, T)", "none" = "", "A" = "Alt\n(Non-risk, C)")) +
  guides(linetype = guide_legend(
    override.aes = (list(shape = 21, fill = c("black", "white", "white"),
                         linetype = c("solid", "solid", "dashed"),
                         color = c("white", "black", "black"))),
    title = "Oligo Type"
  ),
  color = "none",
  fill = "none") +
  theme(legend.position = c(0.875, 0.175),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="black"),
        text=element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "#818181") +
  scale_y_continuous(limits = c(-1, 1.75))

ggsave(plot = rs163_ins_plot, filename = paste0(fig_dir, "rs1635852-ins.png"),
       dpi = 600, units = "in", width = 6.25, height = 5, device = ragg::agg_png())

rs118_ins_df <- ratio_plots %>%
  select(c(refname:config, ratio1:ratio5)) %>%
  filter(rsid == "rs11819995" & site == "right" & config == "ins") %>%
  pivot_longer(!c(refname:config)) %>%
  mutate(allele = factor(allele, levels = c("R", "", "A")),
         type = factor(type, levels = c("orig", "del", "shf")))

rs118_ins_res <- rs118_ins_df %>%
  rstatix::wilcox_test(value ~ refname)

rs118_ins_plot <- rs118_ins_df %>%
  ggplot(aes(x = allele, y = value, color = interaction(type,allele))) +
  geom_boxplot(aes(linetype = type, fill = interaction(type,allele)),
               outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_point(aes(group = interaction(allele,type)),
             pch = 21, size = 3, color = "black",
             fill = "white", position = position_dodge(width = 0.75),
             stroke = 0.75) +
  theme_bw(base_size = 15) + labs(x = "Allele", y = expression(paste("log"[2],"(RNA/DNA), rs11819995"))) +
  scale_color_manual(values = c("#ffffff","#ff0000","#ff0000",
                                "#ffffff","#0000ff", "#0000ff")) +
  scale_fill_manual(values = c("#ff0000", "#ffffff", "#ffffff",
                               "#0000ff", "#ffffff", "#ffffff")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"),
                        labels = c("Original", "Deleted", "Shuffled")) +
  scale_x_discrete(labels = c("R" = "Ref\n(Non-risk, C)", "none" = "", "A" = "Alt\n(Risk, T)")) +
  guides(linetype = guide_legend(
    override.aes = (list(shape = 21, fill = c("black", "white", "white"),
                         linetype = c("solid", "solid", "dashed"),
                         color = c("white", "black", "black"))),
    title = "Oligo Type"
  ),
  color = "none",
  fill = "none") +
  theme(legend.position = c(0.875, 0.175),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="black"),
        text=element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "#818181") +
  scale_y_continuous(limits = c(-1.25, 3))

ggsave(plot = rs118_ins_plot, filename = paste0(fig_dir, "rs11819995-ins_plot.png"),
       dpi = 600, units = "in", width = 6.25, height = 5, device = ragg::agg_png())

rs163_df <- ratio_plots %>%
  select(c(refname:config, ratio1:ratio5)) %>%
  filter(rsid == "rs1635852" & site == "left") %>%
  pivot_longer(!c(refname:config)) %>%
  mutate(allele = case_when(type %in% c("del", "shf") ~ "",
                            TRUE ~ allele),
         allele = factor(allele, levels = c("R", "", "A")),
         type = factor(type, levels = c("orig", "del", "shf")))

rs163_scp_res <- rs163_df %>%
  filter(config == "scp") %>%
  rstatix::wilcox_test(value ~ refname)

prom_label <- c("ins" = "INS", "scp" = "SCP1")

rs163_joint_plot <- rs163_df %>%
  ggplot(aes(x = allele, y = value, color = interaction(type,allele))) +
  facet_wrap(. ~ config, labeller = labeller(config = prom_label)) +
  geom_boxplot(aes(linetype = type, fill = interaction(type,allele)),
               outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_point(aes(group = interaction(allele,type)),
             pch = 21, size = 2.5, color = "black",
             fill = "white", position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.35),
             stroke = 0.75) +
  theme_bw(base_size = 15) + labs(x = "Allele", y = expression(paste("log"[2], "(RNA/DNA), rs1635852"))) +
  scale_color_manual(values = c("#ffffff","#000000","#000000",
                                "#ffffff","#0000ff", "#0000ff")) +
  scale_fill_manual(values = c("#ff0000", "#ffffff", "#ffffff",
                               "#0000ff", "#ffffff", "#ffffff")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"),
                        labels = c("Original", "Deleted", "Shuffled")) +
  scale_x_discrete(labels = c("R" = "Ref\n(Risk, T)", "none" = "", "A" = "Alt\n(Non-risk, C)")) +
  guides(linetype = guide_legend(
    override.aes = (list(shape = 21, fill = c("black", "white", "white"),
                         linetype = c("solid", "solid", "dashed"),
                         color = c("white", "black", "black"))),
    title = "Oligo Type"
  ),
  color = "none",
  fill = "none") +
  theme(legend.position = c(0.1, 0.22),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="black"),
        text=element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "#818181") +
  scale_y_continuous(breaks = seq(-5, 1, by = 1), labels = seq(-5, 1, by = 1),
                     limits = c(-5.6, 1.75))

ggsave(plot = rs163_joint_plot, filename = paste0(fig_dir, "rs1635852-joint_plot.png"),
       dpi = 600, units = "in", width = 8, height = 5, device = ragg::agg_png())

rs118_df <- ratio_merge %>%
  select(c(refname:config, ratio1:ratio5)) %>%
  filter(rsid == "rs11819995" & site == "right") %>%
  pivot_longer(!c(refname:config)) %>%
  mutate(allele = factor(allele, levels = c("R", "", "A")),
         type = factor(type, levels = c("orig", "del", "shf")))

rs118_scp_res <- rs118_df %>%
  filter(config == "scp") %>%
  rstatix::wilcox_test(value ~ refname)

rs118_joint_plot <- rs118_df %>%
  ggplot(aes(x = allele, y = value, color = interaction(type,allele))) +
  facet_wrap(. ~ config, labeller = labeller(config = prom_label)) +
  geom_boxplot(aes(linetype = type, fill = interaction(type,allele)),
               outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_point(aes(group = interaction(allele,type)),
             pch = 21, size = 2.5, color = "black",
             fill = "white", position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.35),
             stroke = 0.75) +
  theme_bw(base_size = 15) + labs(x = "Allele", y = expression(paste("log"[2],"(RNA/DNA), rs11819995"))) +
  scale_color_manual(values = c("#ffffff","#ff0000","#ff0000",
                                "#ffffff","#0000ff", "#0000ff")) +
  scale_fill_manual(values = c("#ff0000", "#ffffff", "#ffffff",
                               "#0000ff", "#ffffff", "#ffffff")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"),
                        labels = c("Original", "Deleted", "Shuffled")) +
  scale_x_discrete(labels = c("R" = "Ref\n(Non-risk, C)", "none" = "", "A" = "Alt\n(Risk, T)")) +
  guides(linetype = guide_legend(
    override.aes = (list(shape = 21, fill = c("black", "white", "white"),
                         linetype = c("solid", "solid", "dashed"),
                         color = c("white", "black", "black"))),
    title = "Oligo Type"
  ),
  color = "none",
  fill = "none") +
  theme(legend.position = c(0.1, 0.22),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="black"),
        text=element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "#818181") +
  scale_y_continuous(limits = c(-6.5, 3))

ggsave(plot = rs118_joint_plot, filename = paste0(fig_dir, "rs11819995-joint_plot.png"),
       dpi = 600, units = "in", width = 10, height = 5, device = ragg::agg_png())
