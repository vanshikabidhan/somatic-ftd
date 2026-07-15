options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(ggeffects)
library(broom.mixed)
library(DHARMa)
library(car)
library(ggrepel)
library(performance)
library(marginaleffects)
library(patchwork)

#######################################################################
# Analyse TARDBP
f <- "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files/gene/TARDBP_Coding_variants.tsv"
df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "ID"

df$ID <- factor(df$ID)
df$Sex <- factor(df$Sex, levels = c("M", "F"))
df$Condition <- factor(df$Condition, levels = c("FTLD-TDP C", "Control"))
df$Recruitment_Site <- factor(df$Recruitment_Site)
df$Pool <- factor(df$Pool)
df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
stopifnot(all(df$n_count >= 0))
stopifnot(all(df$n_count <= df$total_cells))
# center age 
df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
df$Hemisphere <- trimws(df$Hemisphere)
df$Hemisphere[df$Hemisphere == ""] <- NA
df$Hemisphere[df$ID == "A3_P"] <- "Right"
df$Hemisphere[df$ID == "O1_P"] <- "Right"
df$Hemisphere <- factor(df$Hemisphere)


# Fit full model 
fit_full <- glmmTMB(
  n_count ~ Condition + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + (1 | Pool) + (1| Recruitment_Site) ,
  family = nbinom2(link = "log"),
  data = df
)

# Check overdispersion
print(check_overdispersion(fit_full)) # No overdispersion detected
summary(fit_full)
# Strong effect of age at death on the mutation burden, Beta = -0.0073532, P_value = 0.000194

# Run interaction 
fit_interaction <- glmmTMB(
  n_count ~ Condition * AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + 
    (1 | Pool) + (1| Recruitment_Site),
  family = nbinom2(link = "log"),
  data = df
)
summary(fit_interaction)
print(check_overdispersion(fit_interaction))
anova(fit_full, fit_interaction) 
# LRT: χ²(1) = 0.99, p = 0.32

#######################################################
# Run TARDBP model based on disease stratification:
df_P <- df[df$Group=="P",]
df_C <- df[df$Group=="C",]

fit_Patient_TARDBP <- glmmTMB(
  n_count ~ AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + 
    (1 | Pool) + (1 | Recruitment_Site)  ,
  family = nbinom2(link = "log"),
  data = df_P,
  control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS"))
)

fit_Control_TARDBP <- glmmTMB(
  n_count ~ AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + 
    (1 | Pool) + (1 | Recruitment_Site)  ,
  family = nbinom2(link = "log"),
  data = df_C,
  control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS"))
)

check_overdispersion(fit_Patient_TARDBP)
check_overdispersion(fit_Control_TARDBP)

coefficients <- c("SexF", "AAD_centered")
fits <- list(Control = fit_Control_TARDBP, Patient = fit_Patient_TARDBP)

for (coeff in coefficients) {
  effects <- do.call(rbind, lapply(names(fits), function(grp) {
    summ <- summary(fits[[grp]])$coefficients$cond
    ci   <- confint(fits[[grp]], parm = coeff, method = "Wald")
    if (coeff %in% rownames(summ)) {
      data.frame(
        group   = grp,
        beta    = summ[coeff, "Estimate"],
        RR      = exp(summ[coeff, "Estimate"]),
        CI_low  = exp(ci[coeff, "2.5 %"]),
        CI_high = exp(ci[coeff, "97.5 %"]),
        p       = summ[coeff, "Pr(>|z|)"]
      )
    } else {
      data.frame(group = grp, beta = NA, RR = NA, CI_low = NA, CI_high = NA, p = NA)
    }
  }))
  
  effects$sig <- ifelse(effects$p < 0.05, "p < 0.05", "ns")
  cat("\n", coeff, "\n"); print(effects)
}
#Controls: β = −0.014, RR = 0.986 per year, 95% CI: 0.975–0.996, p = 0.008
#Patients: β = −0.006, RR = 0.994 per year, 95% CI: 0.989–0.999, p = 0.018
###########################################################################

### Analyse per gene
base_dir <- "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files/gene"
files <- list.files(base_dir, pattern = "_Coding_variants\\.tsv$", full.names = TRUE)

# Check overdispersion per gene using Poisson fit
for (f in files) {
  print(f)
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "ID"
  df$ID <- factor(df$ID)
  df$Sex <- factor(df$Sex, levels = c("M", "F"))
  df$Condition <- factor(df$Condition, levels = c("FTLD-TDP C", "Control"))
  df$Recruitment_Site <- factor(df$Recruitment_Site)
  df$Pool <- factor(df$Pool)
  
  df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
  df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
  df$Hemisphere <- trimws(df$Hemisphere)
  df$Hemisphere[df$Hemisphere == ""] <- NA
  df$Hemisphere[df$ID == "A3_P"] <- "Right"
  df$Hemisphere[df$ID == "O1_P"] <- "Right"
  df$Hemisphere <- factor(df$Hemisphere)
  
  
  fit <- glmmTMB(
    n_count ~ Condition + AAD_centered + Sex + median_DP + Hemisphere + 
      offset(log(total_cells)) +
      (1 | Pool) + (1 | Recruitment_Site) ,
    family = poisson(link = "log"),
    data = df
  )
  
  print(check_overdispersion(fit))
}

# Overdispersion detected from GRN (ratio = 4.272, p-value = < 0.001)  

#Per-gene model
fits <- list()
for (f in files) {
  gene <- sub("_Coding_variants\\.tsv$", "", basename(f)) 
  gene <- sub("\\.tsv$", "", gene)
  print(gene)
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "ID"
  
  df$ID <- factor(df$ID)
  df$Sex <- factor(df$Sex, levels = c("M", "F"))
  df$Condition <- factor(df$Condition, levels = c("FTLD-TDP C", "Control"))
  df$Recruitment_Site <- factor(df$Recruitment_Site)
  df$Pool <- factor(df$Pool)
  df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
  stopifnot(all(df$n_count >= 0))
  stopifnot(all(df$n_count <= df$total_cells))
  df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
  df$Hemisphere <- trimws(df$Hemisphere)
  df$Hemisphere[df$Hemisphere == ""] <- NA
  df$Hemisphere[df$ID == "A3_P"] <- "Right"
  df$Hemisphere[df$ID == "O1_P"] <- "Right"
  df$Hemisphere <- factor(df$Hemisphere)
  
  fit <- glmmTMB(
    n_count ~ Condition + AAD_centered + Sex + median_DP + Hemisphere + (1 | Pool) + (1 | Recruitment_Site) +
      offset(log(total_cells)) ,
    family = nbinom2(link = "log"),
    data = df, 
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
 
  obj_name <- paste0("fit_nb_", gene)
  assign(obj_name, fit, envir = .GlobalEnv)   
  fits[[gene]] <- fit   
  print(check_overdispersion(fit))
  print(summary(fit))
}


# Of note: GRN has the lowest dispersion, and the smallest n_count values (no. of somatic variants per individual)
# Extract the coefficients and p-value
coefficients <- c("ConditionControl", "SexF", "AAD_centered")
for (coeff in coefficients) {
  print(coeff)
  effects <- sapply(fits, function(model) {
    summ <- summary(model)$coefficients$cond
    ci <- confint(model, parm = coeff, method = "Wald")  # or method = "profile"
    
    if (coeff %in% rownames(summ)) {
      c(
        beta    = summ[coeff, "Estimate"],
        RR      = exp(summ[coeff, "Estimate"]),
        CI_low  = exp(ci[coeff, "2.5 %"]),
        CI_high = exp(ci[coeff, "97.5 %"]),
        p       = summ[coeff, "Pr(>|z|)"]
      )
    } else {
      c(beta = NA, RR = NA, CI_low = NA, CI_high = NA, p = NA)
    }
  })
  
  effects <- as.data.frame(t(effects))
  effects$gene <- rownames(effects)
  effects$p_adj <- p.adjust(effects$p, method = "bonferroni")
  effects$sig <- ifelse(effects$p_adj < 0.05, "p < 0.05", "ns")
  print(effects)
}

# For condition - check with marginal means per condition
condition_emm <- lapply(names(fits), function(g) {
  emm <- emmeans(fits[[g]], ~ Condition, type = "response")
  cont <- as.data.frame(contrast(emm, method = "pairwise"))
  cont$gene <- g
  cont
})
condition_df <- do.call(rbind, condition_emm)
condition_df$p_adj <- p.adjust(condition_df$p.value, method = "BH")
condition_df$sig <- ifelse(condition_df$p_adj < 0.05, "p < 0.05", "ns")
print(condition_df)
# No differences/effect observed between conditions or sexes.

# Marginal Mean for Age (Continuous predictor)
for (g in names(fits)) {
  cat("\n====", g, "====\n")
  # Age trend
  emm_age <- emtrends(fits[[g]], ~ 1, var = "AAD_centered")
  cat("\nAge trend:\n")
  print(emm_age)
}
# gives slope (trend) per group if interaction exists

test(emtrends(fits[['TARDBP']], ~ 1, var = "AAD_centered"))
# TARDBP shows 0.73% lower burden per year corresponding to a 7.1% lower mutation burden for additional decade of age

#######################################################################
# Forest plot to visualize the effect of age at death on mutation burden per gene
p <- ggplot(effects, aes(x = RR, y = reorder(gene, RR), color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.8) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                height = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("p < 0.05" = "Red", "ns" = "#457B9D")) +
  labs(
    x = "Incidence rate ratio per year (95% CI)",
    y = NULL,
    color = NULL
  ) +
  geom_text(
    data = effects[effects$sig == "p < 0.05", ],
    aes(x = effects[effects$sig == "p < 0.05", "RR"],
        label = sprintf("IRR=%.3f (%.1f%% %s/year)\np_adj=%.3f",
                        RR,
                        abs((RR - 1) * 100),
                        ifelse(RR < 1, "decrease", "increase"),
                        p_adj)),
    vjust = -0.5, size = 4, color = "Red"
  ) +
  coord_cartesian(clip = "off") +
  theme_set(theme_classic(base_size = 14, base_family = "DejaVu Sans")) + 
  theme(
    # axis labels
    axis.title.x     = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y     = element_text(size = 18, margin = margin(r = 10)),
    # tick labels
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 16, face = "italic"),
    # tick marks
    axis.ticks       = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    # axis line width
    axis.line        = element_line(linewidth = 1.1),
    # transparent background
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "none"
  )
p

ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/panels/panel4/figure3_nolegend.svg", plot = p, width = 9.5, height = 7, units = "in", bg = "transparent")
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/panels/panel4/figure3_nolegend.png", plot = p, width = 9.5, height = 7, units = "in", dpi = 300, bg = "transparent")

###########################################################################
## Perform a NB GLM on a joined dataframe to compare burden across genes
genes <- c('TARDBP', 'TET2', 'OPTN', 'TBK1', 'GRN', 'UNC13A',  'TMEM106B')
# Number of base pairs for which variants are called across the gene 
gene_len_df <- data.frame(
  gene = genes,
  gene_length = c(TARDBP = 1099, OPTN = 1733, TBK1 = 1729, GRN = 1532, UNC13A = 4426, TET2= 4291, TMEM106B = 806)[genes])

dd_list <- lapply(files, function(f) {
  gene <- sub("_Coding_variants\\.tsv$", "", basename(f))
  gene <- sub("\\.tsv$", "", gene)
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "ID"
  # factors
  df$ID <- factor(df$ID)
  df$Sex <- factor(df$Sex, levels = c("M", "F"))
  df$Condition <- factor(df$Condition, levels = c("FTLD-TDP C", "Control"))
  df$Recruitment_Site <- factor(df$Recruitment_Site)
  df$Pool <- factor(df$Pool)
  df$gene <- factor(gene)
  df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
  stopifnot(all(df$n_count >= 0))
  stopifnot(all(df$n_count <= df$total_cells))
  df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
  df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
  df$Hemisphere <- trimws(df$Hemisphere)
  df$Hemisphere[df$Hemisphere == ""] <- NA
  df$Hemisphere[df$ID == "A3_P"] <- "Right"
  df$Hemisphere[df$ID == "O1_P"] <- "Right"
  df$Hemisphere <- factor(df$Hemisphere)
  df
})

dd <- do.call(rbind, dd_list)
dd$gene <- factor(dd$gene)
dd$gene <- relevel(dd$gene, ref = "TARDBP") #
dd <- merge(dd, gene_len_df, by = "gene", all.x = TRUE)
dd$gene <- relevel(factor(dd$gene), ref = "TARDBP")

fit_joint_nb <- glmmTMB(n_count ~ Condition + gene + AAD_centered + Sex + median_DP + Hemisphere + 
                          offset(log(total_cells)) + offset(log(gene_length)) +
                          (1 | Pool) +  (1 | ID) + (1 | Recruitment_Site) , family = nbinom2(link = "log"), data = dd) #
summary(fit_joint_nb)
print(check_overdispersion(fit_joint_nb))

sim_res <- simulateResiduals(fit_joint_nb, n = 1000)
plot(sim_res)
testDispersion(sim_res) # Underdispersion detected = 0.35111, p-value = 0.002; Outlier test; p= 2e-05

testZeroInflation(sim_res)
plotResiduals(sim_res, dd$gene) # Seems to be driven by GRN and UNC13A (Leven Test for homogeneity of variance)
plotResiduals(sim_res, dd$AAD_centered)
plotResiduals(sim_res, dd$Condition)
plotResiduals(sim_res, dd$AAD_centered, quantreg = TRUE)

## Factors causing Dharma misfit
# 1. Zero Inflation
check_zeroinflation(fit_joint_nb) # Underfitting zeroes
# 2. Check if it's condition-driven  -> NO
plotResiduals(sim_res, form = dd$Condition)
# 3. Check if it's gene-driven -> Yes
plotResiduals(sim_res, form = dd$gene)
#Check for collinearity between the median depth and gene
boxplot(median_DP ~ gene, data = dd)
# The vif for median depth is not extremely large -> collinearity is not an issue
vif(lm(n_count ~ gene + median_DP + Condition + AAD_centered + Sex + Hemisphere, data = dd))

# Use a gene specific dispersio parameter
fit_joint_nb_test2 <- glmmTMB(
  n_count ~ Condition + gene + AAD_centered + Sex + median_DP + Hemisphere + 
    offset(log(total_cells)) + offset(log(gene_length)) +
    (1 | Pool) + (1| ID) + (1 | Recruitment_Site),
  family = nbinom2(link="log"),
  dispformula = ~gene,  # gene-specific dispersion
  data = dd
)
summary(fit_joint_nb_test2)
# Model convergence issues for TBK1 

fit_joint_nb_test3 <- glmmTMB(
  n_count ~ Condition + gene + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + offset(log(gene_length)) +
    (1 | Pool) + (1 | Recruitment_Site) + (1|ID),
  family = nbinom2(link="log"),
  dispformula = ~gene,# gene-specific dispersion
  control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS")),
  data = dd
)
summary(fit_joint_nb_test3) #Model converges well

sim_res <- simulateResiduals(fit_joint_nb_test3, n = 1000)
plot(sim_res)
testDispersion(sim_res) # No dispersion: 0.71307, p-value = 0.606
testZeroInflation(sim_res) #ratioObsSim = 1.0834, p-value = 1
plotResiduals(sim_res, dd$gene) 
anova(fit_joint_nb,fit_joint_nb_test3) # LRT: χ²(6) = 313.78, p = < 2.2e-16

### Forest plot: Differences in mutation burden of each gene, compared to TARDBP
fit_joint <- fit_joint_nb_test3
summary(fit_joint)
# Plot the average burden per gene
emm_joint <- emmeans(fit_joint,~ gene, type   = "response",at  = list(total_cells  = 1, gene_length  = 1))
pairs(emm_joint, adjust = "Bonferroni")

# Summary Table for mutation rate per gene
# get pairwise contrasts vs TARDBP (or all pairs)
contrasts_joint <- contrast(emm_joint, method = "trt.vs.ctrl", 
                            ref = which(emm_joint@grid$gene == "TARDBP"),
                            reverse = TRUE) %>%
  summary(infer = TRUE, type = "response", adjust = "bonferroni")

emm_table <- as.data.frame(emm_joint) %>%
  rename(
    mutation_rate = response,
    LCL = asymp.LCL,
    UCL = asymp.UCL
  ) %>%
  mutate(
    rate_scaled = mutation_rate * 1e4,
    LCL_scaled  = LCL * 1e4,
    UCL_scaled  = UCL * 1e4,
    rate_ci = sprintf("%.2f (%.2f–%.2f)", rate_scaled, LCL_scaled, UCL_scaled)
  ) %>%
  left_join(
    contrasts_joint %>%
      as.data.frame() %>%
      mutate(gene = sub(" / TARDBP|TARDBP / ", "", contrast)) %>%
      select(gene, ratio, asymp.LCL, asymp.UCL, p.value) %>%
      rename(fold_vs_TARDBP = ratio,
             fold_LCL = asymp.LCL,
             fold_UCL = asymp.UCL),
    by = "gene"
  ) %>%
  mutate(
    fold_ci = sprintf("%.2f (%.2f–%.2f)", fold_vs_TARDBP, fold_LCL, fold_UCL),
    p_txt   = format.pval(p.value, digits = 2, eps = 1e-16)
  ) %>%
  arrange(desc(mutation_rate)) %>%
  select(gene, rate_ci, fold_ci, p_txt) %>%
  rename(
    "Gene"                                           = gene,
    "Mutation rate per bp per cell (×10⁻⁴, 95% CI)" = rate_ci,
    "Fold change vs. TARDBP (95% CI)"                = fold_ci,
    "P value"                                        = p_txt
  )
emm_table


#########################################
# =============================================================================
# Forest plot — TARDBP-referenced ratios (per gene)
# =============================================================================
sig_colors <- c("p < 0.001" = "#1a6faf", "p < 0.05" = "#56a0d3", "ns" = "grey60")
res_emm_gene <- contrast(emm_joint,
                         method  = "trt.vs.ctrl",
                         ref     = which(emm_joint@grid$gene == "TARDBP"),
                         reverse = FALSE) %>%
  summary(infer = TRUE, type = "response", adjust = "none") %>%
  as.data.frame() %>%
  mutate(
    gene  = sub("/.*", "", contrast) %>% trimws(),
    rr    = ratio,
    lo    = asymp.LCL,
    hi    = asymp.UCL,
    rr_ci = sprintf("%.2f (%.2f–%.2f)", rr, lo, hi),
    p_txt = format.pval(p.value, digits = 2),
    sig   = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"
    ),
    sig = factor(sig, levels = c("p < 0.001", "p < 0.05", "ns"))
  ) %>%
  arrange(rr) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

p_forest_gene <- ggplot(res_emm_gene, aes(x = rr, y = gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick", alpha = 0.8, linewidth = 1.1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = sig),
                 height = 0.18, linewidth = 0.5) +
  geom_point(aes(color = sig), size = 2.2) +
  scale_color_manual(
    values = sig_colors, name = NULL,
    guide  = guide_legend(override.aes = list(
      linetype = 1, linewidth = 0.5, shape = 16, size = 2.2
    ))
  ) +
  scale_x_continuous(breaks = c(0.5, 0.75, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Mutation Burden Ratio (relative to TARDBP)", y = NULL) +
  theme_classic(base_size = 11, base_family = "DejaVu Sans") +
  theme(
    axis.line.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    axis.title.x      = element_text(size = 16, margin = margin(t = 10)),
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 14, face = "italic"),
    axis.line         = element_line(linewidth = 1.1),
    axis.ticks        = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    legend.position   = c(1.25, -0.05),
    legend.text       = element_text(size = 16),
    legend.background = element_rect(fill = "transparent", colour = NA),
    plot.margin       = margin(5.5, 2, 5.5, 5.5),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA)
  )

p_table_gene <- ggplot(res_emm_gene, aes(y = gene)) +
  geom_text(aes(x = 0, label = rr_ci), hjust = 0, size = 5.5) +
  #annotate("text", x = 0, y = Inf, label = "Mutation Burden Ratio (95% CI)",
           # hjust = 0, vjust = 1, fontface = "bold", size = 5) +
  scale_x_continuous(limits = c(0, 1.4), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

p_combined_gene <- (p_forest_gene | p_table_gene) +
  plot_layout(widths = c(1.4, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

p_combined_gene
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/burden/forest_plot_gene.svg", plot = p_combined_gene, width = 9, height = 6, units = "in", bg = "transparent")
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/burden/forest_plot_gene.png", plot = p_combined_gene, width = 9, height = 6, units = "in", dpi = 300, bg = "transparent")


###########################################
# Plot the Mean mutational burden with CI's
## Plot with data points
dd$burden_norm <- dd$n_count/(dd$total_cells * dd$gene_length)
emm_df <- as.data.frame(emm_joint)
gene_order <- c('GRN', 'TET2', 'UNC13A', 'OPTN', 'TBK1', 'TMEM106B', 'TARDBP')
emm_df <- emm_df %>% mutate(gene_num = match(gene, gene_order))

pairs_df <- as.data.frame(pairs(emm_joint, adjust = "bonferroni")) %>%
  filter(grepl("TARDBP", contrast)) %>%
  filter(p.value < 0.05) %>%
  mutate(
    gene1 = gsub(" /.*", "", contrast),
    gene2 = gsub(".* / ", "", contrast),
    x1    = match(gene1, gene_order),
    x2    = match(gene2, gene_order),
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  )

y_min  <- max(emm_df$asymp.UCL)
y_step <- max(emm_df$asymp.UCL) * 0.25

pairs_df <- pairs_df %>%
  arrange(desc(x2)) %>%
  mutate(y_pos   = y_min + row_number() * y_step,
         tip_len = y_step * 0.15)  # vertical tip length

space = 0.00015
p <- ggplot() +
  geom_col(data = emm_df,
           aes(x = gene_num, y = response),
           fill = "grey80", colour = "black", linewidth = 0.5, width = 0.6) +
  geom_jitter(data = dd %>% mutate(gene_num = match(gene, gene_order)),
              aes(x = gene_num, y = burden_norm),
              width = 0.15, size = 1.5, alpha = 0.5, colour = "black") +
  geom_errorbar(data = emm_df,
                aes(x = gene_num, ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.15, linewidth = 1, colour = "Red") +
  geom_point(data = emm_df,
             aes(x = gene_num, y = response),
             size = 2, colour = "Red") +
  # horizontal bracket
  geom_segment(data = pairs_df,
               aes(x = x1, xend = x2, y = y_pos + space, yend = y_pos + space),
               inherit.aes = FALSE) +
  # left vertical tip
  geom_segment(data = pairs_df,
               aes(x = x1, xend = x1, y = y_pos + space, yend = y_pos  + space - tip_len),
               inherit.aes = FALSE) +
  # right vertical tip
  geom_segment(data = pairs_df,
               aes(x = x2, xend = x2, y = y_pos  + space, yend = y_pos  + space- tip_len),
               inherit.aes = FALSE) +
  # significance label
  geom_text(data = pairs_df,
            aes(x = (x1 + x2) / 2, y = y_pos + space, label = sig_label),
            inherit.aes = FALSE, size = 4, vjust = -0.3) +
  scale_x_continuous(
    breaks = seq_along(gene_order),
    labels = gene_order
  ) +
  # scale_y_continuous(
  #   labels = scales::scientific,
  #   expand = expansion(mult = c(0, 0.2))
  # ) +
  scale_y_continuous(
    labels = function(x) x * 1e4,
    breaks = c(0, 2e-4, 4e-4, 6e-4),
    limits = c(0, 65e-5),
    expand = expansion(mult = c(0, 0)),
    name = expression("Length-normalised mutational burden (×10"^-4*")")
  ) +
  theme_classic(base_size = 14, base_family = "DejaVu Sans") +
  labs(x = "Gene") +
  theme(
    axis.title.x      = element_text(size = 18, margin = margin(t = 1)),
    axis.title.y      = element_text(size = 18, margin = margin(r = 10)),
    axis.text.x       = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 16),
    axis.ticks        = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    axis.line         = element_line(linewidth = 1.1),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA)
  )
p
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/panels/panel4/figure2.svg", plot = p, width = 10.5, height = 7.5, units = "in", bg = "transparent")
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/panels/panel4/figure2.png", plot = p, width = 10.5, height = 7.5, units = "in", dpi = 300, bg = "transparent")

######################################
#Plot the pair differences in mutational burden
pairs_result <- as.data.frame(pairs(emm_joint, adjust = "Bonferroni")) #Gives ratio
pairs_result$pval_adj <- ifelse(pairs_result$p.value < 0.05, "p < 0.05", "ns")

p <- ggplot(pairs_result, aes(x = ratio, 
                         y = reorder(contrast, ratio),
                         color = pval_adj,
                         alpha = pval_adj)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ratio - 1.96*SE, 
                     xmax = ratio + 1.96*SE), 
                 height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", 
             color = "firebrick", alpha = 1) +
  scale_color_manual(values = c("p < 0.05" = "#1a6faf", "ns" = "grey60")) +
  scale_alpha_manual(values = c("p < 0.05" = 1, 
                                "ns" = 0.6)) +
  # scale_y_discrete(labels = function(x) parse(text = paste0("italic(", x, ")")))+
  labs(
    x = "Mutation Burden Ratio (Gene1 / Gene2)",
    y = NULL,
    #title = "Pairwise somatic mutation rates between genes"
  ) +
  theme_classic(base_size = 14, base_family = "DejaVu Sans")+
    theme(
      axis.title.x      = element_text(size = 16, margin = margin(t = 10)),
      axis.title.y      = element_text(size = 16, margin = margin(r = 10)),
      axis.text.x       = element_text(size = 14),
      axis.text.y       = element_text(size = 14),
      axis.ticks        = element_line(linewidth = 1.1),
      axis.ticks.length = unit(6, "pt"),
      axis.line         = element_line(linewidth = 1.1),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      panel.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )
  
p
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/burden/burden_pairwise_gene.svg", plot = p, width = 9, height = 6, units = "in", bg = "transparent")
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/burden/burden_pairwise_gene.png", plot = p, width = 9, height = 6, units = "in", dpi = 300, bg = "transparent")


############ Interactions #######################
# 1. Interaction to confirm whether no patient-vs-control differences are observed across the genes (without condition as a main effect) 
fit_interaction1 <- glmmTMB(
  n_count ~ gene*Condition + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + offset(log(gene_length)) +
    (1 | Pool) + (1 | Recruitment_Site) , # Removed ID as OPTN did not converge even with BFGS optimizer
  family = nbinom2,
  dispformula = ~gene,
  data = dd
)
summary(fit_interaction1)

# Plot relative burden between patients and controls for each gene -> Highlights TARDBP with highest mean burden, 
# However, no significant differences between condition groups.
## 1) Get EMMs once (so rates + p-values are consistent)
emm <- emmeans(
  fit_interaction1, ~ Condition | gene,
  type = "response",
  at = list(total_cells = 1, gene_length = 1),
  adjust = "bonferroni" 
)

rate_df <- as.data.frame(summary(emm))     
p_df <- as.data.frame(pairs(emm)) %>%      
  dplyr::select(gene, p.value) %>%
  distinct() %>%
  mutate(p.adj = p.adjust(p.value, method = "bonferroni"))
p_df

## 2) Convert to wider format
wide_rates <- rate_df %>%
  dplyr::select(gene, Condition, response) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = response)

plot_df <- wide_rates %>%
  left_join(p_df, by = "gene") %>%
  mutate(
    mlog10p = -log10(p.value),
    mean_burden = rowMeans(across(where(is.numeric) & !any_of(c("p.value","mlog10p"))), na.rm = TRUE)
  )

# Compute shared limits and breaks for the axes
shared_min <- min(plot_df$Control, plot_df$`FTLD-TDP C`, na.rm = TRUE)
shared_max <- max(plot_df$Control, plot_df$`FTLD-TDP C`, na.rm = TRUE)
# shared_breaks <- scales::trans_breaks("log10", function(x) 10^x, n = 6)(
#   c(shared_min, shared_max), 
# )
shared_breaks <- seq(5e-5, shared_max, by = 2.5e-5)

# Label position at upper end of equality line
x_label   <- max(shared_max, na.rm = TRUE) 
eq_label_y <- x_label
p <- ggplot(plot_df, aes(x = Control, y = `FTLD-TDP C`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey70", linewidth = 0.5) +
  geom_point(size = 3, color = "grey70", alpha = 0.9) +
  geom_text_repel(
    aes(label = gene),
    size               = 4.5,
    box.padding        = 0.5,
    point.padding      = 0.3,
    min.segment.length = 0.2,
    segment.color      = "grey60",
    max.overlaps       = Inf
  ) +
  annotate("text", x = x_label, y = eq_label_y,
           label    = "Equal\nburden",
           color    = "grey50",
           hjust    = 0, vjust=-0.2, size = 3.5,
           fontface = "italic") +
  labs(
    x     = expression("Mutation burden in Controls (×10"^{-5}*")"),
    y     = expression("Mutation burden in FTLD-TDP Cases (×10"^{-5}*")"),
  ) +
  scale_x_continuous(
    breaks = shared_breaks,
    labels = function(x) format(x * 1e5, digits = 2),
    limits = c(shared_min, shared_max)
  ) +
  scale_y_continuous(
    breaks = shared_breaks,
    labels = function(x) format(x * 1e5, digits = 2),
    limits = c(shared_min, shared_max)
  )+
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14, base_family = "DejaVu Sans") +
  theme(
    axis.title.x      = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 18, margin = margin(r = 10)),
    axis.text.x       = element_text(size = 16, hjust = 1),
    axis.text.y       = element_text(size = 16),
    axis.ticks        = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    axis.line         = element_line(linewidth = 1.1),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA),
  )
p <- p +
  annotate("text",
           x        = shared_min,
           y        = shared_max,
           label    = "All pairwise p > 0.05\n(Bonferroni-adjusted)",
           hjust    = 0,
           vjust    = 1,
           size     = 5,
           color    = "grey40",
           fontface = "italic")
p
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/burden/burden_patientvscontrol.svg", plot = p, width = 10, height = 8, units = "in", bg = "transparent")
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/burden/burden_patientvscontrol.png", plot = p, width = 10, height = 8, units = "in", dpi = 300, bg = "transparent")

########### Age at death vs Condition trend for TARDBP ##################
df <- read.table("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files/gene/TARDBP_Coding_variants.tsv", header=T, sep="\t")
colnames(df)[colnames(df) == "X"] <- "ID"
df$ID <- as.factor(df$ID)
df$Sex <- factor(df$Sex, levels=c("M", "F"))
df$Condition <- factor(df$Condition, levels=c("FTLD-TDP C", "Control"))
df$Recruitment_Site <- as.factor(df$Recruitment_Site)
df$Pool <- as.factor(df$Pool)
df$AAD_centered <- scale(df$Age_at_Death, center = TRUE, scale = FALSE)
df$Hemisphere <- trimws(df$Hemisphere)
df$Hemisphere[df$Hemisphere == ""] <- NA
df$Hemisphere[df$ID == "A3_P"] <- "Right"
df$Hemisphere[df$ID == "O1_P"] <- "Right"
df$Hemisphere <- factor(df$Hemisphere)

fit_nb <- glmmTMB(
  n_count ~ Condition + Age_at_Death + Sex + median_DP+ Hemisphere + 
    (1 | Pool) + (1| Recruitment_Site) +
    offset(log(total_cells)) ,
  family = nbinom2(link = "log"),
  control = glmmTMBControl(optimizer = optim,optArgs = list(method = "BFGS")),
  data = df
)
summary(fit_nb)

# Fit a regression line
# Predict
df$burden_rate <- df$n_count / df$total_cells
newdat <- expand.grid(
  Age_at_Death     = seq(min(df$Age_at_Death),
                         max(df$Age_at_Death),
                         length.out = 100),
  Condition        = c("FTLD-TDP C", "Control"),
  Sex              = "M" ,
  Hemisphere       = "Left", 
  median_DP        = median(df$median_DP),
  total_cells      = median(df$total_cells),
  Pool             = df$Pool[1],              
  Recruitment_Site = df$Recruitment_Site[1]   
)

pred <- predict(fit_nb, newdata = newdat,
                type = "link", se.fit = TRUE, re.form = NA)

newdat$fit_link  <- pred$fit
newdat$se_link   <- pred$se.fit

# These are counts — need to divide by total_cells for rate
newdat$pred_rate  <- exp(newdat$fit_link) / newdat$total_cells
newdat$conf.low   <- exp(newdat$fit_link - 1.96 * newdat$se_link) / newdat$total_cells
newdat$conf.high  <- exp(newdat$fit_link + 1.96 * newdat$se_link) / newdat$total_cells

## Plot
p <- ggplot(df, aes(Age_at_Death, burden_rate, fill = Condition)) +
  geom_point(alpha = 0.9, size = 3.5, shape = 21, colour = "black", stroke = 0.5) +
  geom_ribbon(
    data = newdat,
    aes(x = Age_at_Death, ymin = conf.low, ymax = conf.high,
        fill = Condition),
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = newdat,
    aes(y = pred_rate, color = Condition),
    linewidth = 1.2
  ) +
  scale_color_manual(values = c(
    "Control"    = "#92C7FF",
    "FTLD-TDP C" = "#F44542"
  )) +
  scale_fill_manual(values = c(
    "Control"    = "#92C7FF",
    "FTLD-TDP C" = "#F44542"
  )) +
  labs(
    x     = "Age at death (years)",
    y     = expression(italic("TARDBP") ~ "mutational burden"),
    color = "Condition",
    fill  = "Condition"
  ) +
  theme_classic(base_size = 14, base_family = "DejaVu Sans") +
  theme(
    axis.title.x      = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 18, margin = margin(r = 10)),
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 16),
    axis.ticks        = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    axis.line         = element_line(linewidth = 1.1),
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 14),
    #legend.position   = "bottom_left",
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA)
  )
p
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/panels/panel4/figure4_legend.svg", plot = p, width = 9.5, height = 7, units = "in", bg = "transparent")
ggsave("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/panels/panel4/figure4_legend.png", plot = p, width = 9.5, height = 7, units = "in", dpi = 300, bg = "transparent")

############################################################################
# Comparison for all genes together (19,102 variants)
f <- '/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files/All_Coding_variants.tsv'
df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df %>% rename(ID = X)
df$ID <- factor(df$ID)
df$Sex <- factor(df$Sex, levels = c("M", "F"))
df$Condition <- factor(df$Condition, levels = c("FTLD-TDP C", "Control"))
df$Recruitment_Site <- factor(df$Recruitment_Site)
df$Pool <- factor(df$Pool)
df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
stopifnot(all(df$n_count >= 0))
stopifnot(all(df$n_count <= df$total_cells))
# center age 
df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
df$Hemisphere <- trimws(df$Hemisphere)
df$Hemisphere[df$Hemisphere == ""] <- NA
df$Hemisphere[df$ID == "A3_P"] <- "Right"
df$Hemisphere[df$ID == "O1_P"] <- "Right"
df$Hemisphere <- factor(df$Hemisphere)

fit_all <- glmmTMB(
  n_count ~ Condition + AAD_centered + Sex + median_DP + Hemisphere + 
    offset(log(total_cells)) + 
    (1 | Pool) + (1 | Recruitment_Site)  ,
  family = poisson(link = "log"),
  data = df
)

print(check_overdispersion(fit_all)) # Poisson works (no overdispersion) 
sim_res_all <- simulateResiduals(fit_all, n = 1000)
plot(sim_res_all)

summary(fit_all)
# No effect of condition: No difference in the overall mutation burden between patients and controls
