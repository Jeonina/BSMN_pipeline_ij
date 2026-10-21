# Reproduce Bae et al. 2022 (Science, doi:10.1126/science.abm6222) style figures
# for the Hypermutability cohort: Fig. 1 (A, B, C, E, F), Fig. 3A and Table 1.
#
# Run:  Rscript analysis/paper_figures.R  [results_dir]  [phenotype_xlsx]  [qc_mode]
#   qc_mode = detailed (default): technical outliers labelled by cause
#           = pending: technical outliers shown neutrally as "QC pending"
#             (for presentations before the QC findings are discussed);
#             outputs go to paper_figures_qcpending/
#
# Inputs in results_dir (copied from the server):
#   cohort_audit.tsv     scripts/cohort_audit.py output (per-sample QC + counts)
#   somatic_variants.tsv bcftools +split-vep table, no header:
#                        sample chrom pos ref alt AF AD SYMBOL Consequence IMPACT
#                        gnomADg_AF gnomADe_AF Existing_variation
#   cohort.SBS96.all     SigProfilerMatrixGenerator SBS96 matrix
#   sample_flags.tsv     hand-curated technical flags (sample, flag, note);
#                        flag = contamination | read_end_artifact
# Output: results_dir/paper_figures/*.png, *.tsv
#
# Fig. 1D (haplotype phasing) and 1G (aneuploidy) are not reproduced: no phasing
# or aneuploidy calls exist for this cohort.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(readxl); library(stringr)
  library(ggplot2); library(patchwork); library(scales); library(ggrepel); library(forcats)
  library(gridExtra); library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else
  "~/Desktop/BaeLab/Project_Data/BSMN/20261020_200WGS_results"
pheno_xlsx <- if (length(args) >= 2) args[2] else
  "~/Desktop/BaeLab/Project_Data/BSMN/Hypermutability_Corrected_Barcode_List_20250722_JHS_08_25_2025.xlsx"
qc_mode <- if (length(args) >= 3) args[3] else "detailed"
res_dir <- path.expand(res_dir); pheno_xlsx <- path.expand(pheno_xlsx)
out_dir <- file.path(res_dir, if (qc_mode == "pending") "paper_figures_qcpending" else "paper_figures")
dir.create(out_dir, showWarnings = FALSE)

# ---- palette (validated categorical slots, fixed order) --------------------------
dx_levels <- c("Control", "Bipolar", "Alzheimer's", "Autism")
dx_cols   <- c(Control = "#2a78d6", Bipolar = "#eb6834", "Alzheimer's" = "#1baf7a", Autism = "#eda100")
flag_cols <- c("Non-hypermutable" = "grey65", Hypermutable = "#111111",
               Contamination = "#e34948", "Read-end artifact" = "#6250d6")
if (qc_mode == "pending") flag_cols <- c("Non-hypermutable" = "grey65", Hypermutable = "#111111",
                                         "Above fence, QC pending" = "white")
clean_label <- if (qc_mode == "pending") "samples passing QC" else "non-hypermutable samples"
sbs6_cols <- c("C>A" = "#03BCEE", "C>G" = "#010101", "C>T" = "#E32926",
               "T>A" = "#CAC9C9", "T>C" = "#A1CE63", "T>G" = "#EBC6C4")   # SigProfiler convention
theme_set(theme_classic(base_size = 11) +
            theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
                  plot.title = element_text(face = "bold", size = 12),
                  plot.tag = element_text(face = "bold", size = 14)))
save_png <- function(p, name, w, h) ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300, bg = "white")

# ---- load ----------------------------------------------------------------------
audit <- read_tsv(file.path(res_dir, "cohort_audit.tsv"), show_col_types = FALSE) |>
  filter(filtered_ok %in% c(TRUE, "True"))
vars <- read_tsv(file.path(res_dir, "somatic_variants.tsv"), col_names = c(
  "sample", "chrom", "pos", "ref", "alt", "vaf", "ad", "symbol", "consequence", "impact",
  "gnomadg", "gnomade", "existing"), col_types = cols(.default = "c")) |>
  mutate(pos = as.integer(pos), vaf = as.numeric(vaf),
         key = paste(chrom, pos, ref, alt, sep = ":"))
flags <- read_tsv(file.path(res_dir, "sample_flags.tsv"), show_col_types = FALSE)

pheno <- read_excel(pheno_xlsx, sheet = "Phenotype") |>
  transmute(donor = BrNum, age = AgeDeath, sex = Sex,
            dx = recode(PrimaryDx, Alzheimer = "Alzheimer's"))
sinfo <- read_excel(pheno_xlsx, sheet = "SampleInfo") |>
  transmute(sample = `Sample ID`, donor = str_extract(`Br&Region`, "^[^_]+"),
            region = str_extract(`Br&Region`, "[^_]+$"))

samples <- audit |>
  select(sample, final_snv, est_depth, contamination) |>
  inner_join(sinfo, by = "sample") |> inner_join(pheno, by = "donor") |>
  left_join(flags |> select(sample, flag), by = "sample") |>
  mutate(dx = factor(dx, dx_levels))

# ---- per-donor burden (union of both regions) and hypermutability call ----------
donor_counts <- vars |> inner_join(sinfo, by = "sample") |>
  group_by(donor) |> summarise(n_union = n_distinct(key), .groups = "drop")
donors <- samples |>
  group_by(donor, age, sex, dx) |>
  summarise(n_regions = n(), tech_flag = paste(na.omit(unique(flag)), collapse = ","),
            depth = mean(est_depth), .groups = "drop") |>
  left_join(donor_counts, by = "donor") |>
  mutate(n_union = coalesce(n_union, 0L)) |>
  filter(n_regions == 2)   # donors with both regions only (a region still being reprocessed is excluded)
q <- quantile(donors$n_union, c(.25, .75)); fence <- unname(q[2] + 3 * (q[2] - q[1]))   # upper outer fence
donors <- donors |> mutate(
  status = case_when(
    n_union > fence & str_detect(tech_flag, "contamination") ~ "Contamination",
    n_union > fence & str_detect(tech_flag, "read_end")      ~ "Read-end artifact",
    n_union > fence                                           ~ "Hypermutable",
    TRUE                                                      ~ "Non-hypermutable"),
  status = if (qc_mode == "pending")
    if_else(status %in% c("Contamination", "Read-end artifact"), "Above fence, QC pending", status) else status,
  status = factor(status, names(flag_cols)))
hyper_donors <- donors$donor[donors$status == "Hypermutable"]
samples <- samples |> mutate(
  group = case_when(!is.na(flag) ~ "Technical",
                    donor %in% hyper_donors ~ "Hypermutable",
                    TRUE ~ "Non-hypermutable"))
clean_samples <- samples$sample[samples$group == "Non-hypermutable"]
write_tsv(donors, file.path(out_dir, "donor_burden.tsv"))
message(sprintf("donors=%d  fence=%.1f  hypermutable=%s", nrow(donors), fence,
                paste(hyper_donors, collapse = ",")))

# ---- Fig. 1A: per-donor burden by diagnosis -------------------------------------
med_dx <- donors |> filter(status == "Non-hypermutable") |> group_by(dx) |>
  summarise(med = median(n_union), .groups = "drop")
med_ctrl <- med_dx$med[med_dx$dx == "Control"]
f1a <- ggplot(donors, aes(dx, n_union)) +
  geom_hline(yintercept = fence, linetype = "dotted") +
  geom_crossbar(data = med_dx, aes(dx, med, ymin = med, ymax = med, colour = dx),
                inherit.aes = FALSE, width = 0.6, linewidth = 0.8, show.legend = FALSE) +
  annotate("text", x = 2.5, y = fence * 1.1, label = sprintf("fence = %.0f", fence), size = 3) +
  geom_hline(yintercept = med_ctrl, linewidth = 0.4) +
  geom_quasirandom(aes(fill = status), shape = 21, size = 2.3, width = 0.25, stroke = 0.3) +
  geom_text_repel(data = filter(donors, n_union > fence, qc_mode != "pending" | status == "Hypermutable"),
                  aes(label = donor), size = 3, min.segment.length = 0) +
  scale_y_log10() + scale_colour_manual(values = dx_cols) + scale_fill_manual(values = flag_cols, drop = FALSE) +
  labs(tag = "A", x = NULL, y = "Somatic SNVs per donor\n(union of DLPFC + caudate, log)", fill = NULL,
       title = sprintf("Burden per donor (n = %d)", nrow(donors))) +
  theme(legend.position = "bottom")

# ---- Fig. 1B: combined SBS96 spectrum of non-hypermutable samples ---------------
sbs <- read_tsv(file.path(res_dir, "cohort.SBS96.all"), show_col_types = FALSE)
sbs_long <- sbs |> pivot_longer(-MutationType, names_to = "sample", values_to = "n") |>
  mutate(sub = str_sub(MutationType, 3, 5),
         ctx = paste0(str_sub(MutationType, 1, 1), str_sub(MutationType, 3, 3), str_sub(MutationType, 7, 7)))
pooled <- sbs_long |> filter(sample %in% clean_samples) |> group_by(MutationType, sub, ctx) |>
  summarise(n = sum(n), .groups = "drop") |> mutate(frac = n / sum(n)) |> arrange(sub, ctx) |>
  mutate(MutationType = factor(MutationType, MutationType))
f1b <- ggplot(pooled, aes(MutationType, frac, fill = sub)) + geom_col(width = 0.75) +
  facet_grid(~sub, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = sbs6_cols, guide = "none") + scale_y_continuous(labels = percent) +
  scale_x_discrete(labels = function(x) paste0(str_sub(x, 1, 1), str_sub(x, 3, 3), str_sub(x, 7, 7))) +
  labs(tag = "B", x = NULL, y = "Fraction of SNVs",
       title = sprintf("Combined spectrum, %d %s (%s SNVs)", length(clean_samples), clean_label, comma(sum(pooled$n)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, family = "mono"),
        panel.spacing.x = unit(1, "pt"))

# ---- Fig. 1C: substitution spectra by group --------------------------------------
sbs6 <- sbs_long |> inner_join(samples |> select(sample, dx, group, donor), by = "sample") |>
  filter(group != "Technical") |>
  mutate(set = if_else(group == "Hypermutable", paste0(donor, " (hyper)"), as.character(dx))) |>
  group_by(set, sub) |> summarise(n = sum(n), .groups = "drop") |>
  group_by(set) |> mutate(frac = n / sum(n), tot = sum(n)) |> ungroup() |>
  mutate(set = factor(set, c(dx_levels, sort(setdiff(unique(set), dx_levels)))))
f1c <- ggplot(sbs6, aes(set, frac, fill = sub)) + geom_col(width = 0.75) +
  geom_text(data = distinct(sbs6, set, tot), aes(set, 1.03, label = comma(tot)), inherit.aes = FALSE, size = 2.8) +
  scale_fill_manual(values = sbs6_cols) + scale_y_continuous(labels = percent, breaks = seq(0, 1, .25)) +
  labs(tag = "C", x = NULL, y = "Fraction of SNVs", fill = NULL, title = "Substitution spectra") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- Fig. 1E: burden vs age, hypermutable proportion by age bin ------------------
age_bins <- donors |> mutate(bin = cut(age, c(0, 40, 60, Inf), labels = c("<40", "40-60", ">60"))) |>
  group_by(bin) |> summarise(k = sum(status == "Hypermutable"), n = n(), .groups = "drop") |>
  rowwise() |> mutate(p = k / n, lo = binom.test(k, n)$conf.int[1], hi = binom.test(k, n)$conf.int[2]) |> ungroup()
ct <- with(filter(donors, status == "Non-hypermutable"), cor.test(age, n_union, method = "spearman", exact = FALSE))
f1e_top <- ggplot(age_bins, aes(bin, p)) + geom_col(fill = "#9ec5f4", width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  geom_text(aes(label = paste0(k, "/", n)), vjust = -0.6, size = 3, nudge_x = 0.22) +
  scale_y_continuous(labels = percent) + labs(tag = "E", x = "Age at death", y = "Hypermutable\ndonors")
f1e_bot <- ggplot(donors, aes(age, n_union)) +
  geom_hline(yintercept = fence, linetype = "dotted") +
  geom_point(aes(colour = dx, shape = sex, alpha = status %in% c("Non-hypermutable", "Hypermutable")), size = 2, stroke = 0.9) +
  geom_text_repel(data = filter(donors, status == "Hypermutable"), aes(label = donor), size = 3) +
  scale_y_log10() + scale_colour_manual(values = dx_cols) + scale_shape_manual(values = c(F = 4, M = 1)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(x = "Age at death (years)", y = "SNVs per donor (log)", colour = NULL, shape = NULL,
       subtitle = sprintf("Non-hypermutable donors: Spearman rho = %.2f, p = %.2g (faded = %s)",
                          ct$estimate, ct$p.value, if (qc_mode == "pending") "QC pending" else "technical outliers"))
f1e <- f1e_top / f1e_bot + plot_layout(heights = c(1, 2.2))

# ---- Fig. 1F: VAF distributions --------------------------------------------------
vaf_df <- vars |> inner_join(samples |> select(sample, dx, group, donor), by = "sample") |>
  filter(group != "Technical") |>
  mutate(set = if_else(group == "Hypermutable", paste0(donor, " (hyper)"), as.character(dx)))
vaf_n <- vaf_df |> count(set)
vaf_df <- vaf_df |> left_join(vaf_n, by = "set") |>
  mutate(lab = factor(paste0(set, " (", comma(n), ")")),
         lab = fct_relevel(lab, grep(paste(dx_levels, collapse = "|"), levels(lab), value = TRUE)))
f1f <- ggplot(vaf_df, aes(vaf, lab)) +
  geom_violin(fill = "grey90", colour = "grey50", scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_x_log10(labels = percent, breaks = c(.02, .05, .1, .2, .4)) +
  labs(tag = "F", x = "Variant allele fraction (log)", y = NULL, title = "VAF distributions (mutation counts)")

fig1 <- (f1a | f1e) / (f1b) / (f1c | f1f) + plot_layout(heights = c(1.25, 0.8, 1))
save_png(fig1, "Fig1_panels_ABCEF.png", 13, 14)
save_png(f1a, "Fig1A_burden_per_donor.png", 6.5, 5)
save_png(f1b, "Fig1B_combined_spectrum.png", 11, 3.2)
save_png(f1c, "Fig1C_substitution_spectra.png", 6, 4)
save_png(f1e, "Fig1E_age.png", 6.5, 6)
save_png(f1f, "Fig1F_vaf.png", 6.5, 4)

# ---- Fig. 3A: genomic distribution and VEP impact (non-hypermutable samples) -----
loc_of <- function(csq) {
  first <- str_split_fixed(csq, "&", 2)[, 1]
  case_when(
    str_detect(csq, "stop_gained|frameshift|start_lost|stop_lost|splice_acceptor|splice_donor") ~ "Splice/stop (HIGH)",
    str_detect(csq, "missense") ~ "Missense",
    str_detect(csq, "synonymous") ~ "Synonymous",
    str_detect(csq, "splice_region") ~ "Splice region",
    str_detect(csq, "UTR") ~ "UTR",
    str_detect(csq, "non_coding_transcript_exon") ~ "ncRNA exon",
    first == "intron_variant" ~ "Intron",
    str_detect(first, "upstream|downstream") ~ "Up/downstream",
    first == "intergenic_variant" ~ "Intergenic",
    TRUE ~ "Other")
}
loc_levels <- c("Intergenic", "Up/downstream", "Intron", "ncRNA exon", "UTR", "Splice region",
                "Synonymous", "Missense", "Splice/stop (HIGH)", "Other")
loc_cols <- setNames(c("#cde2fb", "#9ec5f4", "#5598e7", "#1baf7a", "#eda100", "#e87ba4",
                       "#8c8c8c", "#eb6834", "#e34948", "#f2f2f2"), loc_levels)
v3 <- vars |> filter(sample %in% clean_samples) |>
  inner_join(samples |> select(sample, dx), by = "sample") |>
  mutate(loc = factor(loc_of(consequence), loc_levels),
         impact = factor(impact, c("HIGH", "MODERATE", "LOW", "MODIFIER")))
tot3 <- count(v3, dx)
f3a_loc <- ggplot(v3, aes(dx, fill = fct_rev(loc))) + geom_bar(position = "fill", width = 0.7) +
  geom_text(data = tot3, aes(dx, 1.04, label = comma(n)), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = loc_cols, breaks = loc_levels) + scale_y_continuous(labels = percent) +
  labs(tag = "A", x = NULL, y = "Fraction of SNVs", fill = "Location (VEP)", title = "Genomic distribution")
imp_tab <- v3 |> count(dx, impact) |> group_by(dx) |> mutate(frac = n / sum(n)) |> ungroup() |>
  filter(impact != "MODIFIER")
f3a_imp <- ggplot(imp_tab, aes(dx, frac, fill = impact)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.75), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = c(HIGH = "#e34948", MODERATE = "#eb6834", LOW = "#eda100")) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, .15))) +
  labs(x = NULL, y = "Fraction of SNVs", fill = "VEP impact",
       title = "Predicted impact (MODIFIER omitted)")
fig3a <- f3a_loc | f3a_imp
save_png(fig3a, "Fig3A_genomic_distribution_vep.png", 12, 4.8)
write_tsv(v3 |> count(dx, loc) |> group_by(dx) |> mutate(frac = n / sum(n)), file.path(out_dir, "Fig3A_table.tsv"))

# ---- Table 1: protein-altering variants in CH / cancer genes, and all HIGH --------
# Edit this list as needed; the paper used Bailey et al. 2018, Cancer Gene Census,
# SFARI (ASD) and published SCZ gene lists.
gene_set <- c(
  # clonal hematopoiesis
  "DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "JAK2", "SF3B1", "SRSF2", "U2AF1", "CBL", "GNB1", "GNAS",
  "IDH1", "IDH2", "BCOR", "BCORL1", "STAG2", "ZRSR2", "KRAS", "NRAS",
  # recurrent cancer / brain somatic drivers
  "PIK3CA", "PTEN", "AKT3", "MTOR", "TSC1", "TSC2", "BRAF", "EGFR", "NF1", "ATRX", "CIC", "FUBP1", "PDGFRA",
  "H3-3A", "H3F3A", "SETD2", "ARID1A", "CTNNB1", "APC", "NOTCH1")
tab1 <- vars |> inner_join(samples |> select(sample, donor, region, dx, age, group), by = "sample") |>
  filter(impact == "HIGH" | (impact == "MODERATE" & symbol %in% gene_set)) |>
  mutate(Mutation = sprintf("%s:%s %s>%s", chrom, comma(pos), ref, alt),
         Hypermutable = if_else(group == "Hypermutable", "Yes", "No"),
         Gene = if_else(symbol == ".", "(lncRNA, unnamed)", symbol),
         `Gene set` = if_else(symbol %in% gene_set, "CH / cancer", "-"),
         Consequence = str_replace_all(consequence, "_variant|&non_coding_transcript", "") |> str_replace_all("_", " ")) |>
  arrange(Mutation, region) |>
  transmute(Mutation, Donor = donor, Region = region, Phenotype = as.character(dx), Age = round(age),
            Hypermutable, VAF = percent(vaf, 0.1), Consequence, Gene, `Gene set`)
write_tsv(tab1, file.path(out_dir, "Table1_damaging_variants.tsv"))
tg <- tableGrob(tab1, rows = NULL, theme = ttheme_minimal(
  base_size = 10, core = list(fg_params = list(hjust = 0, x = 0.04)),
  colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold"))))
ggsave(file.path(out_dir, "Table1_damaging_variants.png"), arrangeGrob(tg,
       top = grid::textGrob("Table 1. Protein-altering somatic SNVs in CH/cancer genes and all HIGH-impact SNVs",
                            x = 0.02, hjust = 0, gp = grid::gpar(fontface = "bold"))),
       width = 12, height = 0.6 + 0.32 * (nrow(tab1) + 2), dpi = 300, bg = "white")

message("written: ", out_dir)
