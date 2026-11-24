# Patient vs control burden per gene: full data + downsampled sets (scatter plot): Figure 7E (Supplementary)
#
# gene x Condition model (= fit_interaction1, gene-model.R), EMMs per Condition | gene at
# total_cells = 1, gene_length = 1; contrasts Bonferroni-adjusted across genes.
#
# In the downsampled sets some genes' NB dispersion goes to the Poisson limit and the site variance to ~0,
# so the Hessian is not positive-definite although all fixed effects are estimable. Fits are kept if they
# converge with finite fixed-effect SEs; otherwise the model is simplified. 

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
})

args     <- commandArgs(trailingOnly = TRUE)
base     <- "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files"
full_dir <- if (length(args) >= 1) args[1] else file.path(base, "gene")
sets_dir <- if (length(args) >= 2) args[2] else file.path(base, "gene_halves")
out_dir  <- if (length(args) >= 3) args[3] else file.path(sets_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene_len  <- c(TARDBP = 1099, OPTN = 1733, TBK1 = 1841, GRN = 1471, UNC13A = 4426, TET2 = 4291, TMEM106B = 806)
gene_cols <- c(TARDBP = "#117733", TMEM106B = "#AA4499", TBK1 = "#CC6677", OPTN = "#88CCEE",
               UNC13A = "#DDCC77", TET2 = "#332288", GRN = "#44AA99")
set_names <- list.files(sets_dir, pattern = "^set_[0-9]+$")
set_names <- set_names[order(as.integer(sub("set_", "", set_names)))]
datasets  <- c(Full = full_dir, setNames(file.path(sets_dir, set_names), sub("set_", "Set ", set_names)))

# same preprocessing as dd in gene-model.R; age centred within each gene file
build_dd <- function(dir) {
  files <- list.files(dir, pattern = "_Coding_variants\\.tsv$", full.names = TRUE)
  dd <- dplyr::bind_rows(lapply(files, function(f) {
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "ID"
    df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
    stopifnot(all(df$n_count >= 0), all(df$n_count <= df$total_cells))
    df$gene <- sub("_Coding_variants\\.tsv$", "", basename(f))
    df$AAD_centered <- df$Age_at_Death - mean(df$Age_at_Death, na.rm = TRUE)
    df$Hemisphere <- trimws(df$Hemisphere)
    df$Hemisphere[df$Hemisphere == ""] <- NA
    df$Hemisphere[df$ID %in% c("A3_P", "O1_P")] <- "Right"
    df
  }))
  dplyr::mutate(dd,
                ID = factor(ID), Sex = factor(Sex, levels = c("M", "F")),
                Condition = factor(Condition, levels = c("FTLD-TDP C", "Control")),
                Recruitment_Site = factor(Recruitment_Site), Pool = factor(Pool), Hemisphere = factor(Hemisphere),
                gene = relevel(factor(gene), ref = "TARDBP"), gene_length = unname(gene_len[as.character(gene)]))
}

ctrls <- list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

# NaN SEs of dispersion parameters (Poisson-limit genes) only trigger a warning in summary()
fe_ok <- function(f) !is.null(f) && f$fit$convergence == 0 &&
  all(is.finite(suppressWarnings(summary(f))$coefficients$cond[, "Std. Error"]))

# first spec x optimiser giving an acceptable fit
fit_first <- function(fixed, specs, data) {
  for (sp in specs) for (ctl in ctrls) {
    form <- as.formula(paste(paste(deparse(fixed, width.cutoff = 500L), collapse = ""), "+", sp$re))
    f <- tryCatch(suppressWarnings(glmmTMB(form, data = data, family = sp$fam, dispformula = sp$disp,
                                           control = ctl)), error = function(e) NULL)
    if (fe_ok(f)) { attr(f, "model") <- sp$lab; attr(f, "pdHess") <- isTRUE(f$sdr$pdHess); return(f) }
  }
  NULL
}

re_full <- "(1 | Pool) + (1 | Recruitment_Site)"
spec_int <- list(
  list(lab = "NB, dispersion ~gene (original)",        re = re_full,      fam = nbinom2, disp = ~gene),
  list(lab = "NB, dispersion ~gene, no site effect",  re = "(1 | Pool)", fam = nbinom2, disp = ~gene),
  list(lab = "NB, common dispersion, no site effect", re = "(1 | Pool)", fam = nbinom2, disp = ~1),
  list(lab = "Poisson, no site effect",               re = "(1 | Pool)", fam = poisson, disp = ~1))

f_int  <- n_count ~ gene * Condition + AAD_centered + Sex + median_DP + Hemisphere +
  offset(log(total_cells)) + offset(log(gene_length))

# ---- condition: interaction model ----
condition_by_gene <- function(dd, label) {
  fit <- fit_first(f_int, spec_int, dd)
  if (is.null(fit)) { warning(label, ": no acceptable interaction fit"); return(NULL) }
  emm  <- emmeans(fit, ~ Condition | gene, type = "response", at = list(total_cells = 1, gene_length = 1))
  s    <- as.data.frame(summary(emm))
  pick <- function(d, opts) d[[intersect(opts, names(d))[1]]]
  rates <- data.frame(gene = s$gene, Condition = ifelse(s$Condition == "Control", "C", "P"),
                      est = pick(s, c("response", "rate"))) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = est, names_prefix = "rate_")
  pr <- as.data.frame(confint(pairs(emm)))                 # ratio = FTLD-TDP C / Control
  pv <- as.data.frame(pairs(emm))
  ctr <- data.frame(gene = pr$gene, RR = pick(pr, c("ratio", "estimate")),
                    CI_low = pick(pr, c("asymp.LCL", "lower.CL")), CI_high = pick(pr, c("asymp.UCL", "upper.CL")),
                    p = pv$p.value[match(pr$gene, pv$gene)])
  ctr$p_adj <- p.adjust(ctr$p, "bonferroni")
  dplyr::left_join(rates, ctr, by = "gene") %>%
    dplyr::mutate(dataset = label, model = attr(fit, "model"), pdHess = attr(fit, "pdHess"))
}

cond <- list()
for (d in names(datasets)) {
  cat("--", d, "\n")
  cond[[d]] <- condition_by_gene(build_dd(datasets[[d]]), d)
}
cond <- dplyr::bind_rows(cond) %>%
  dplyr::mutate(dataset = factor(dataset, levels = names(datasets)), is_full = dataset == "Full")
stopifnot(nrow(cond) > 0)

write.table(cond, file.path(out_dir, "condition_full_and_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(dplyr::distinct(cond, dataset, model, pdHess))
print(dplyr::filter(cond, p_adj < 0.05) %>% dplyr::select(dataset, gene, RR, p_adj))

# ---- scatter: patient vs control burden ----

full <- dplyr::filter(cond, is_full); sets <- dplyr::filter(cond, !is_full)
lims <- range(c(cond$rate_C, cond$rate_P), na.rm = TRUE)
brks <- seq(5e-5, lims[2], by = 2.5e-5)

# annotation replaces the legend; only claims "all p > 0.05" if that holds in every data set
n_sig <- sum(cond$p_adj < 0.05, na.rm = TRUE)
sig_label <- if (n_sig == 0) "All pairwise p > 0.05\n(Bonferroni-adjusted;\nfull and downsampled data)" else
  sprintf("%d of %d gene-level comparisons\np < 0.05 (Bonferroni-adjusted)", n_sig, sum(!is.na(cond$p_adj)))

p <- ggplot(mapping = aes(x = rate_C, y = rate_P, colour = gene, fill = gene)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey70", linewidth = 0.5) +
  geom_point(data = sets, size = 3, alpha = 0.5, stroke = 0) +                       # downsampled sets
  geom_point(data = full, shape = 21, size = 4, colour = "black", stroke = 0.6) +  # full data
  geom_text_repel(
    data               = full,
    aes(label = gene),
    size               = 5,
    fontface           = "italic",
    box.padding        = 0.5,
    point.padding      = 0.3,
    min.segment.length = 0.2,
    segment.colour     = "grey60",
    max.overlaps       = Inf
  ) +
  annotate("text", x = lims[2], y = lims[2],
           label    = "Equal\nburden",
           colour   = "grey50",
           hjust    = 0, vjust = -0.2, size = 3.5,
           fontface = "italic") +
  annotate("text", x = lims[1], y = lims[2],
           label    = sig_label,
           hjust    = 0, vjust = 1, size = 5,
           colour   = "grey40",
           fontface = "italic") +
  scale_colour_manual(values = gene_cols, guide = "none") +
  scale_fill_manual(values = gene_cols, guide = "none") +
  scale_x_continuous(
    breaks = brks,
    labels = function(x) format(x * 1e5, digits = 2),
    limits = lims
  ) +
  scale_y_continuous(
    breaks = brks,
    labels = function(x) format(x * 1e5, digits = 2),
    limits = lims
  ) +
  labs(
    x = expression("Mutation burden in controls (×10"^{-5}*")"),
    y = expression("Mutation burden in FTLD-TDP type C patients (×10"^{-5}*")")
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14, base_family = "DejaVu Sans") +
  theme(
    axis.title.x      = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y      = element_text(size = 18, margin = margin(r = 10)),
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 16),
    axis.ticks        = element_line(linewidth = 1.1),
    axis.ticks.length = unit(6, "pt"),
    axis.line         = element_line(linewidth = 1.1),
    legend.position   = "none",
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent", colour = NA)
  )
p
for (ext in c("svg", "png"))
  ggsave(file.path(out_dir, paste0("burden_patientvscontrol_full_and_sets.", ext)), p,
         width = 10, height = 8, units = "in", dpi = 300, bg = "transparent")

cat("Results written to", out_dir, "\n")
