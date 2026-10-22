# Sample identity QC figures (scripts/sample_identity.py analyze output).
#   A  identity heatmap: chr22 allele-fraction r, samples ordered by donor then region;
#      the two regions of a donor form 2x2 blocks on the diagonal, a bright cell off them = swap/duplicate
#   B  r distribution: same-donor pairs vs different-donor pairs, with the call threshold
#   C  sex: chrX heterozygosity vs chrY/chr22 depth, coloured by metadata sex
# Usage: Rscript analysis/sample_identity_figures.R <identity_dir> <sample_info.tsv> [min_r]
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
id_dir <- path.expand(if (length(args) >= 1) args[1] else "~/Desktop/BaeLab/Project_Data/BSMN/20261020_200WGS_results/qc_identity")
info_path <- path.expand(if (length(args) >= 2) args[2] else file.path(dirname(id_dir), "sample_info.tsv"))
min_r <- if (length(args) >= 3) as.numeric(args[3]) else 0.8

sex_cols <- c(M = "#2a78d6", F = "#eb6834")
pair_cols <- c("Same donor" = "#2a78d6", "Different donors" = "#eb6834")
theme_set(theme_classic(base_size = 11) +
            theme(plot.title = element_text(face = "bold", size = 12),
                  plot.tag = element_text(face = "bold", size = 14)))
save_png <- function(p, name, w, h) ggsave(file.path(id_dir, name), p, width = w, height = h, dpi = 300, bg = "white")

info <- read_tsv(info_path, show_col_types = FALSE) |>
  separate(br_region, c("donor", "region"), sep = "_", remove = FALSE)
mat <- read_tsv(file.path(id_dir, "identity_matrix.tsv"), show_col_types = FALSE, name_repair = "minimal")
names(mat)[1] <- "a"
r_long <- mat |> pivot_longer(-a, names_to = "b", values_to = "r") |>
  mutate(r = as.numeric(r)) |>
  left_join(info |> select(a = br_region, donor_a = donor), by = "a") |>
  left_join(info |> select(b = br_region, donor_b = donor), by = "b")
pairs <- read_tsv(file.path(id_dir, "identity_pairs.tsv"), show_col_types = FALSE)
sex <- read_tsv(file.path(id_dir, "identity_sex.tsv"), show_col_types = FALSE)

# ---- A: heatmap -----------------------------------------------------------------
lvl <- info |> filter(br_region %in% mat$a) |> arrange(donor, region) |> pull(br_region)
fa <- ggplot(r_long |> mutate(a = factor(a, lvl), b = factor(b, rev(lvl))), aes(a, b, fill = r)) +
  geom_raster() +
  scale_fill_gradient(low = "#f4f4f2", high = "#0d366b", limits = c(min(r_long$r, na.rm = TRUE), 1),
                      name = "r") +
  labs(tag = "A", x = NULL, y = NULL,
       title = sprintf("Genotype correlation, chr22 (n = %d samples)", length(lvl)),
       subtitle = "Ordered by donor, then region: same-donor pairs sit in 2x2 blocks on the diagonal") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey80"))

# ---- B: r distribution ------------------------------------------------------------
dist <- r_long |> filter(as.character(a) < as.character(b)) |>
  mutate(kind = factor(if_else(donor_a == donor_b, "Same donor", "Different donors"), names(pair_cols)))
n_kind <- count(dist, kind)
fb <- ggplot(dist, aes(r, fill = kind)) +
  geom_histogram(bins = 80, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = min_r, linetype = "dashed", colour = "grey35") +
  annotate("text", x = min_r, y = Inf, label = sprintf(" call threshold r = %.1f", min_r),
           hjust = 0, vjust = 1.5, size = 3, colour = "grey25") +
  scale_fill_manual(values = pair_cols, labels = \(k) sprintf("%s (n = %s)", k, format(n_kind$n[match(k, n_kind$kind)], big.mark = ",")),
                    name = NULL) +
  scale_y_sqrt() +
  labs(tag = "B", x = "Genotype correlation r (chr22 allele fractions)", y = "Sample pairs (sqrt scale)",
       title = "Same-donor vs different-donor pairs") +
  theme(legend.position = "top", legend.justification = "left")

# ---- C: sex -------------------------------------------------------------------------
bad_sex <- sex |> filter(status != "OK")
fc <- ggplot(sex, aes(x_het_rate, chrY_chr22_depth, colour = sex_meta)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_point(data = bad_sex, shape = 21, size = 4.5, stroke = 1, colour = "black") +
  geom_text_repel(data = bad_sex, aes(label = br_region), colour = "black", size = 3, min.segment.length = 0) +
  scale_colour_manual(values = sex_cols, name = "Metadata sex") +
  labs(tag = "C", x = "chrX heterozygous-site fraction (non-PAR)", y = "chrY / chr22 depth (MQ >= 30)",
       title = sprintf("Sex check: %d / %d match metadata", sum(sex$status == "OK"), nrow(sex))) +
  theme(legend.position = "top", legend.justification = "left")

# ---- summary + output ---------------------------------------------------------------------
same <- dist |> filter(kind == "Same donor"); diff <- dist |> filter(kind == "Different donors")
message(sprintf("same-donor r: min %.3f median %.3f | different-donor r: median %.3f max %.3f",
                min(same$r), median(same$r), median(diff$r), max(diff$r)))
message("pair status: ", paste(names(table(pairs$status)), table(pairs$status), sep = "=", collapse = "  "))
message("sex status : ", paste(names(table(sex$status)), table(sex$status), sep = "=", collapse = "  "))

save_png(fa, "identity_A_heatmap.png", 7.5, 7)
save_png(fb, "identity_B_r_distribution.png", 7, 4.5)
save_png(fc, "identity_C_sex.png", 6, 5)
save_png((fa | (fb / fc)) + plot_layout(widths = c(1.15, 1)), "identity_panels.png", 15, 8.5)
message("figures -> ", id_dir)
