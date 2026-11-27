library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# --- 1. INPUT / OUTPUT Snakemake ---
file_blast   <- snakemake@input[["blast"]]
file_output  <- snakemake@output[["plot"]]
name_sample <- paste(snakemake@wildcards[["species"]], snakemake@wildcards[["sex"]], sep = "_")

# Security if the TSV file is empty 
if (!file.exists(file_blast) || file.info(file_blast)$size == 0) {
  p_vide <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Blast file empty:", name_sample)) +
    theme_void()
  ggsave(file_output, p_vide, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# --- 2. read et extraction of the best species of each contig ---
data_detected <- read.delim(file_blast, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    sample = name_sample,
    contig      = coalesce(str_extract(qseqid, "^NODE_\\d+"), qseqid),
    coverage    = as.numeric(str_extract(qseqid, "(?<=cov_)[0-9.]+")),
    species      = sseqid %>% str_remove("_[^_]+$") %>% str_replace_all("_", " ")
  ) %>%
  arrange(desc(bitscore), desc(pident)) %>%
  group_by(sample, qseqid, contig, coverage) %>%
  slice(1) %>%
  ungroup() %>%
  filter(bitscore >= 1000) %>%
  arrange(desc(coverage)) %>%
  mutate(
    pct_reads = (coverage / sum(coverage)) * 100,
    species    = factor(species, levels = unique(species)),
    contig    = factor(contig, levels = unique(contig))
  ) %>%
  arrange(species, desc(coverage)) %>%
  mutate(
    y_mid = rev(cumsum(rev(pct_reads))) - (pct_reads / 2)
  )

# Security if any hit has superior to bitscore >= 1000
if (nrow(data_detected) == 0) {
  p_vide <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Any hit BLAST (bitscore >= 1000) :", name_sample)) +
    theme_void()
  ggsave(file_output, p_vide, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# if more than 8 species
palettes_species <- if (nlevels(data_detected$species) <= 8) {
  scale_fill_brewer(palette = "Set2")
} else {
  scale_fill_viridis_d(option = "D")
}

# --- 3. Graphic A : Raw coverage + pident per contig ---
p_cov_raw <- ggplot(data_detected, aes(x = sample, y = coverage, fill = species, group = contig)) +
  geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.55) +
  geom_text(
    aes(label = sprintf("%s\n%.1f× (%.1f%%)", contig, coverage, pident)),
    position = position_dodge(width = 0.7),
    vjust = -0.25, lineheight = 0.95, fontface = "bold", size = 3.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  palettes_species +
  labs(
    title    = "Raw coverage and identity per contig",
    subtitle = "Coverage (×) and BLAST identity of the best-matching species",
    x        = NULL,
    y        = "Mean coverage (×)",
    fill     = "Species detected"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text     = element_text(face = "italic")
  )

# --- 4. Graphic B : Part of relative observed ---
p_part_obs <- ggplot(data_detected, aes(x = sample, y = pct_reads, fill = species, group = contig)) +
  geom_col(width = 0.4, color = "white", linewidth = 0.4) +
  geom_text(
    data = filter(data_detected, pct_reads >= 5),
    aes(y = y_mid, label = sprintf("%.1f %%", pct_reads)),
    color = "white", fontface = "bold", size = 3.8
  ) +
  geom_label_repel(
    data = filter(data_detected, pct_reads < 5),
    aes(y = y_mid, label = sprintf("%.2f %%", pct_reads)),
    color = "white", fontface = "bold", size = 3.5,
    nudge_x = 0.35, show.legend = FALSE,
    min.segment.length = 0, segment.color = "#FC8D62", segment.size = 0.8
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, " %"),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  palettes_species +
  labs(
    title    = "Observed relative proportion",
    subtitle = "Coverage proportion between detected contigs",
    x        = NULL,
    y        = "Share of coverage (%)",
    fill     = "Species detected"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 5. Assembly of plot et save ---
figure_detection <- p_cov_raw + p_part_obs + plot_layout(widths = c(1.2, 1))
ggsave(file_output, figure_detection, width = 10, height = 7, dpi = 300)