# TARDBP burden vs age at death: full data and the six downsampled sets (forest plot).
# Figure 7F
# Per-gene model as for Fig. 4C (gene-model.R):
#   n_count ~ Condition + AAD_centered + Sex + median_DP + Hemisphere + offset(log(total_cells))
#             + (1|Pool) + (1|Recruitment_Site)
# NB first, Poisson if the NB dispersion goes to the Poisson limit; a fit is kept if it converges with
# finite fixed-effect SEs. 
# Sets 1-3 and 4-6 share their patients (half A and B) and all sets share the controls.


options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(glmmTMB); library(dplyr); library(ggplot2) })

args     <- commandArgs(trailingOnly = TRUE)
base     <- "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/Burden_files"
full_dir <- if (length(args) >= 1) args[1] else file.path(base, "gene")
sets_dir <- if (length(args) >= 2) args[2] else file.path(base, "gene_halves")
out_dir  <- if (length(args) >= 3) args[3] else file.path(sets_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set_names <- list.files(sets_dir, pattern = "^set_[0-9]+$")
set_names <- set_names[order(as.integer(sub("set_", "", set_names)))]
datasets  <- c(Full = full_dir, setNames(file.path(sets_dir, set_names), sub("set_", "Set ", set_names)))
half      <- c(Full = "Full data", setNames(ifelse(seq_along(set_names) <= 3, "Half A", "Half B"),
                                            sub("set_", "Set ", set_names)))

prep <- function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "ID"
  df <- df[!is.na(df$total_cells) & df$total_cells > 0, ]
  stopifnot(all(df$n_count >= 0), all(df$n_count <= df$total_cells))
  df$Hemisphere <- trimws(df$Hemisphere)
  df$Hemisphere[df$Hemisphere == ""] <- NA
  df$Hemisphere[df$ID %in% c("A3_P", "O1_P")] <- "Right"
  transform(df, ID = factor(ID), Sex = factor(Sex, levels = c("M", "F")),
            Condition = factor(Condition, levels = c("FTLD-TDP C", "Control")),
            Recruitment_Site = factor(Recruitment_Site), Pool = factor(Pool), Hemisphere = factor(Hemisphere),
            AAD_centered = Age_at_Death - mean(Age_at_Death, na.rm = TRUE))
}

form  <- n_count ~ Condition + AAD_centered + Sex + median_DP + Hemisphere + offset(log(total_cells)) +
  (1 | Pool) + (1 | Recruitment_Site)
ctrls <- list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

fe_ok <- function(f) !is.null(f) && f$fit$convergence == 0 &&
  all(is.finite(suppressWarnings(summary(f))$coefficients$cond[, "Std. Error"]))

fit_tardbp <- function(dir) {
  d <- prep(file.path(dir, "TARDBP_Coding_variants.tsv"))
  for (fam in c("nbinom2", "poisson")) for (ctl in ctrls) {
    f <- tryCatch(suppressWarnings(glmmTMB(form, data = d, family = get(fam)(), control = ctl)),
                  error = function(e) NULL)
    if (fe_ok(f)) { attr(f, "model") <- fam; return(f) }
  }
  NULL
}

res <- dplyr::bind_rows(lapply(names(datasets), function(ds) {
  fit <- fit_tardbp(datasets[[ds]])
  if (is.null(fit)) return(data.frame(dataset = ds, RR = NA, CI_low = NA, CI_high = NA, p = NA, model = NA))
  cf <- suppressWarnings(summary(fit))$coefficients$cond["AAD_centered", ]
  b <- cf[["Estimate"]]; se <- cf[["Std. Error"]]
  data.frame(dataset = ds, RR = exp(b), CI_low = exp(b - 1.96 * se), CI_high = exp(b + 1.96 * se),
             p = cf[["Pr(>|z|)"]], model = attr(fit, "model"))
})) %>%
  dplyr::mutate(p_bonf7 = pmin(p * 7, 1), half = unname(half[dataset]),
                pct_per_year = (RR - 1) * 100,
                dataset = factor(dataset, levels = rev(names(datasets))))          # Full at the top

write.table(res, file.path(out_dir, "TARDBP_age_full_and_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 3)

p <- ggplot(res, aes(x = RR, y = dataset, colour = dataset == "Full")) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey45"), guide = "none") +   # full data in black
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "n.s.", `TRUE` = "p < 0.05"), name = NULL, drop = FALSE) +
  labs(x = expression("exp(β) per year of age at death (95% CI) for " * italic("TARDBP")), y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14, base_family = "DejaVu Sans") +        
  theme(axis.title.x      = element_text(size = 18, margin = margin(t = 10)),
        axis.title.y      = element_text(size = 18, margin = margin(r = 10)),
        axis.text.x       = element_text(size = 16),
        axis.text.y       = element_text(size = 16),
        axis.ticks        = element_line(linewidth = 1.1),
        axis.ticks.length = unit(6, "pt"),
        axis.line         = element_line(linewidth = 1.1),
        legend.text       = element_text(size = 14),
        plot.background   = element_rect(fill = "transparent", colour = NA),
        panel.background  = element_rect(fill = "transparent", colour = NA))
p
for (ext in c("svg", "png"))
  ggsave(file.path(out_dir, paste0("TARDBP_age_full_and_sets.", ext)), p, width = 9.5, height = 7,
         units = "in", dpi = 300, bg = "transparent")