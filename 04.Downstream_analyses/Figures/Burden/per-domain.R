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
library(tidyr)
library(patchwork)
library(performance)
library(DHARMa)
library(tibble)
library(tidyverse)

# =============================================================================
# Load Data
# =============================================================================
df <- read.table("/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files/tardp_domains/TARDBP_Coding_variants_perdomain.tsv", header=T, sep="\t")
colnames(df)[colnames(df) == "X"] <- "ID"
df$ID <- as.factor(df$ID)
df$Sex <- factor(df$Sex, levels=c("M", "F"))
df$Condition <- factor(df$Condition, levels=c("FTLD-TDP C", "Control"))
df$Recruitment_Site <- as.factor(df$Recruitment_Site)
df$Pool <- as.factor(df$Pool)
df$domain <- as.factor(df$domain)
df$non_mutated <- df$total_cells - df$n_count
df$AAD_centered <- as.numeric(scale(df$Age_at_Death, center = TRUE, scale = FALSE))
df$Hemisphere <- trimws(df$Hemisphere)
df$Hemisphere[df$Hemisphere == ""] <- NA
df$Hemisphere[df$ID == "A3_P"] <- "Right"
df$Hemisphere[df$ID == "O1_P"] <- "Right"
df$Hemisphere <- factor(df$Hemisphere)
# =============================================================================
# Model fit
# =============================================================================
fit_poisson_domain <- glmmTMB(
  n_count ~ Condition + domain + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) +
    offset(log(domain_length)) +
    (1 | ID)+
    (1 | Pool) +
    (1 | Recruitment_Site),
  family = poisson(link = "log"),
  data = df
)
check_overdispersion(fit_poisson_domain)
#Overdispersion detected. Dispersion ratio = 1.557, p-value = < 0.001

fit_nb_domain <- glmmTMB(
  n_count ~ Condition + domain + AAD_centered + Sex + median_DP +  Hemisphere +
    offset(log(total_cells)) +
    offset(log(domain_length)) +
    (1 | ID) +
    (1 | Pool)+ 
    (1 | Recruitment_Site),
  family = nbinom2(link = "log"),
  data = df
)
check_overdispersion(fit_nb_domain) # No overdispersion
summary(fit_nb_domain)

fit_nb_domain_noRS <- glmmTMB(
  n_count ~ Condition + domain + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) +
    offset(log(domain_length)) +
    (1 | Pool)+ 
    (1 | ID),
  family = nbinom2(link = "log"),
  data = df
)
check_overdispersion(fit_nb_domain_noRS) # Underdispersion detected
summary(fit_nb_domain_noRS)

fit_nb_domain_noPoolID <- glmmTMB(
  n_count ~ Condition + domain + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) +
    offset(log(domain_length)) +
    (1 | ID) +
    (1 | Recruitment_Site),
  family = nbinom2(link = "log"),
  data = df
)
check_overdispersion(fit_nb_domain_noPoolID) # Underdispersion detected
summary(fit_nb_domain_noPoolID)

# Chosen model -> fit_nb_domain
sim <- simulateResiduals(fit_nb_domain)
plot(sim)

# =============================================================================
# Plotting functions
# =============================================================================
# Output path
OUT_DIR <- "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/"

save_plot <- function(p, name, width, height) {
  for (ext in c("svg", "png")) {
    path <- file.path(OUT_DIR, paste0(name, ".", ext))
    ggsave(path, plot = p, width = width, height = height,
           units = "in", dpi = if (ext == "png") 300 else 100,
           bg = "transparent")
  }
}

# Shared theme elements
theme_custom <- function(base_size = 14, base_family = "DejaVu Sans") {
  theme_classic(base_size = base_size) +
    theme(
      axis.title.x      = element_text(size = 18, margin = margin(t = 10)),
      axis.title.y      = element_text(size = 18, margin = margin(r = 10)),
      axis.text.x       = element_text(size = 16),
      axis.text.y       = element_text(size = 16),
      axis.ticks        = element_line(linewidth = 1.1),
      axis.ticks.length = unit(6, "pt"),
      axis.line         = element_line(linewidth = 1.1),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      panel.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA)
    )
}

sig_colors <- c("p < 0.001" = "#1a6faf", "p < 0.05" = "#56a0d3", "ns" = "grey60")
# =============================================================================
# Marginal means — bar + jitter plot
# =============================================================================
fit_sel      <- fit_nb_domain
domain_levels <- c("NTD", "RRM2", "RRM1", "CTD")

emm    <- emmeans(fit_sel, ~ domain, type = "response", offset = log(1) + log(1))
emm_df <- as.data.frame(summary(emm, infer = TRUE, type = "response", adjust="Bonferroni")) %>%
  mutate(domain = factor(domain, levels = domain_levels))

df_domain <- df %>%
  filter(domain %in% domain_levels) %>%
  mutate(burden_norm = n_count / (total_cells * domain_length))

domain_order <- c(NTD = 1, RRM2 = 2, RRM1 = 3, CTD = 4)
pairs_res <- pairs(emm, adjust = "Bonferroni") %>%
  as_tibble() %>%
  mutate(contrast = as.character(contrast)) %>%
  separate(contrast, into = c("g1", "g2"), sep = " / ") %>%
  mutate(
    g1 = trimws(g1),
    g2 = trimws(g2),
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    x1 = domain_order[g1],
    x2 = domain_order[g2]
  )

tip_len  <- 2e-05
pairs_df <- pairs_res %>%
  select(g1, g2, sig_label, x1, x2) %>%
  mutate(span = abs(x2 - x1)) %>%
  arrange(span) %>%
  mutate(y_pos = seq(9.8e-04, by = 1e-04, length.out = n()),
         y_pos = case_when(
           g1 == "NTD"  & g2 == "RRM2" ~ y_pos[which(g1 == "RRM1" & g2 == "RRM2")],
           g1 == "RRM1" & g2 == "RRM2" ~ y_pos[which(g1 == "NTD"  & g2 == "RRM2")],
           TRUE ~ y_pos
         )) %>%
  select(-span)

p_bar <- ggplot() +
  geom_col(data = emm_df,
           aes(x = factor(domain, levels = domain_levels), y = response),
           fill = "grey80", colour = "black", linewidth = 0.5, width = 0.6) +
  geom_jitter(data = df_domain,
              aes(x = factor(domain, levels = domain_levels), y = burden_norm),
              width = 0.15, size = 1.5, alpha = 0.5, colour = "black") +
  geom_errorbar(data = emm_df,
                aes(x = factor(domain, levels = domain_levels),
                    ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.15, linewidth = 1, colour = "#A32D2D") +
  geom_point(data = emm_df,
             aes(x = factor(domain, levels = domain_levels), y = response),
             size = 2, colour = "#A32D2D") +
  # Significance brackets
  geom_segment(data = pairs_df,
               aes(x = x1, xend = x2, y = y_pos, yend = y_pos),
               inherit.aes = FALSE) +
  geom_segment(data = pairs_df,
               aes(x = x1, xend = x1, y = y_pos, yend = y_pos - tip_len),
               inherit.aes = FALSE) +
  geom_segment(data = pairs_df,
               aes(x = x2, xend = x2, y = y_pos, yend = y_pos - tip_len),
               inherit.aes = FALSE) +
  geom_text(data = pairs_df,
            aes(x = (x1 + x2) / 2, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 5, vjust = -0.3) +
  scale_x_discrete(limits = domain_levels) +
  scale_y_continuous(
    labels = function(x) x * 1e4,
    breaks = seq(0, 6e-3, by = 3e-4),
    limits = c(0, 1.6e-3),
    expand = expansion(mult = c(0, 0)),
    name = expression("Length-normalised mutational burden (×10"^-4*")")
  ) +
  labs(x = "TDP-43 domain") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_bar
save_plot(p_bar, "panels/panel4/figure1", width = 9.5, height = 7)

# =============================================================================
# Forest plot — CTD-referenced ratios
# =============================================================================
emm_joint <- emmeans(fit_sel, ~ domain, type = "response", at = list(total_cells = 1, domain_length = 1))
res_emm <- contrast(emm_joint,
                    method  = "trt.vs.ctrl",
                    ref     = which(emm_joint@grid$domain == "CTD"),
                    reverse = FALSE) %>%
  summary(infer = TRUE, type = "response", adjust = "Bonferroni") %>%
  as.data.frame() %>%
  mutate(
    domain = sub("/.*", "", contrast) %>% trimws(),  # take numerator (now domain, not CTD)
    rr     = ratio,
    lo     = asymp.LCL,
    hi     = asymp.UCL,
    rr_ci  = sprintf("%.2f (%.2f–%.2f)", rr, lo, hi),
    p_txt  = format.pval(p.value, digits = 2),
    sig    = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"
    ),
    sig = factor(sig, levels = c("p < 0.001", "p < 0.05", "ns"))
  ) %>%
  arrange(rr) %>%
  mutate(domain = factor(domain, levels = unique(domain)))

p_forest <- ggplot(res_emm, aes(x = rr, y = domain)) +
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
  scale_x_continuous(breaks = c(1, 1.5, 2, 2.5, 3),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Mutation Burden Ratio (relative to CTD)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.line.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    axis.title.x      = element_text(size = 16, margin = margin(t = 20)),
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 14),
    axis.line         = element_line(linewidth = 1.1),
    axis.ticks        = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    legend.position   = c(1.15, 0.05),
    legend.text       = element_text(size = 16),
    legend.background = element_rect(fill = "transparent", colour = NA),
    plot.margin       = margin(5.5, 2, 5.5, 5.5),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA)
  )

p_table <- ggplot(res_emm, aes(y = domain)) +
  geom_text(aes(x = 0, label = rr_ci), hjust = 0, size = 5.5) +
  annotate("text", x = 0, y = Inf, label = "Mutation Burden\n Ratio (95% CI)",
           hjust = 0, vjust = 1.2, fontface = "bold", size = 5) +
  scale_x_continuous(limits = c(0, 1.4), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

p_combined <- (p_forest | p_table) +
  plot_layout(widths = c(1.4, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))
p_combined

save_plot(p_combined, "supplementary/burden/forest_plot_domain", width = 9, height = 6)

# =============================================================================
# Pairwise comparisons — all domain pairs
# =============================================================================
pairs_result <- as.data.frame(pairs(emm_joint, adjust = "Bonferroni")) %>%
  mutate(pval_adj = ifelse(p.value < 0.05, "p < 0.05", "ns"))

pairs_ci <- confint(pairs(emm_joint), type = "response")

pairs_result <- pairs_result %>%
  left_join(pairs_ci[, c("contrast", "asymp.LCL", "asymp.UCL")], by = "contrast")

p_pairs <- ggplot(pairs_result, aes(x = ratio,
                                    y = reorder(contrast, ratio),
                                    color = pval_adj,
                                    alpha = pval_adj)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "firebrick", alpha = 0.8, linewidth = 1.1) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL), height = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(limits = c(0.2, 1.6), breaks = seq(0, 1.6, by = 0.4)) +
  scale_color_manual(values = c("p < 0.05" = "#1a6faf", "ns" = "grey60")) +
  scale_alpha_manual(values = c("p < 0.05" = 1,         "ns" = 0.6)) +
  labs(x = "Mutation Burden Ratio (Domain1 / Domain2)", y = NULL) +
  theme_custom() +
  theme(
    axis.text.x  = element_text(hjust = 1),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
p_pairs
save_plot(p_pairs, "supplementary/burden/domain_pairwise_comparison", width = 9, height = 6)

# =============================================================================
# Interaction tests (LRT)
# =============================================================================
fit_nb_domain_x_condition <- glmmTMB(
  n_count ~ Condition * domain + AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + offset(log(domain_length)) +
    (1 | Pool) + (1 | ID) + (1| Recruitment_Site),
  family = nbinom2(link = "log"), data = df
)

fit_nb_domain_x_aad <- glmmTMB(
  n_count ~ Condition + domain * AAD_centered + Sex + median_DP + Hemisphere +
    offset(log(total_cells)) + offset(log(domain_length)) +
    (1 | Pool) + (1 | ID) + (1| Recruitment_Site),
  family = nbinom2(link = "log"), data = df
)

anova(fit_nb_domain, fit_nb_domain_x_condition) # χ²(3) = 2.83, p = 0.418
anova(fit_nb_domain, fit_nb_domain_x_aad)       # χ²(3) = 3.13, p = 0.372

# =============================================================================
# Results summary table
# =============================================================================
# as.data.frame(summary(emm)) %>%
#   mutate(estimate_ci = sprintf("%.2f (%.2f–%.2f)",
#                                response * 1e4,
#                                asymp.LCL * 1e4,
#                                asymp.UCL * 1e4))
result_df <- as.data.frame(emm) %>%
  rename(mutation_rate = response, LCL = asymp.LCL, UCL = asymp.UCL) %>%
  mutate(
    rate_ci = sprintf("%.2f (%.2f–%.2f)",
                      mutation_rate * 1e4, LCL * 1e4, UCL * 1e4)
  ) %>%
  left_join(
    contrast(emm, method = "trt.vs.ctrl", ref = "CTD") %>%
      summary(infer = c(TRUE, TRUE), type = "response", adjust = "bonferroni") %>%
      as.data.frame() %>%
      mutate(domain = sub(" / CTD|CTD / ", "", contrast)) %>%
      select(domain, ratio, asymp.LCL, asymp.UCL, p.value) %>%
      rename(fold_vs_ref = ratio,
             fold_LCL = asymp.LCL,
             fold_UCL = asymp.UCL),
    by = "domain"
  ) %>%
  mutate(
    fold_ci = sprintf("%.2f (%.2f–%.2f)", fold_vs_ref, fold_LCL, fold_UCL),
    p_txt   = format.pval(p.value, digits = 2, eps = 1e-16)
  ) %>%
  arrange(desc(mutation_rate)) %>%
  select(domain, rate_ci, fold_ci, p_txt) %>%
  rename(
    "Domain" = domain,
    "Mutation burden (×10⁻⁴, 95% CI)" = rate_ci,
    "Fold change vs. CTD (95% CI)" = fold_ci,
    "P value" = p_txt
  )
result_df