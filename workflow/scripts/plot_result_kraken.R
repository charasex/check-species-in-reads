library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# --- 1. INPUT / OUTPUT Snakemake ---
file_kraken  <- snakemake@input[["report"]]
file_output  <- snakemake@output[["plot"]]
name_sample <- paste(snakemake@wildcards[["species"]], snakemake@wildcards[["sex"]], sep = "_")

# Security if report kraken2 is empty
if (!file.exists(file_kraken) || file.info(file_kraken)$size == 0) {
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("File empty:", name_sample)) +
    theme_void()
  ggsave(file_output, p_empty, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# --- 2. read of kraken report and extraction of 2 principles species ---
data_detected <- read.delim(
  file_kraken, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
  col.names = c("pct", "kmers_clade", "kmers_taxon", "rang", "taxid", "name")
) %>%
  mutate(
    specie = str_trim(name),
    rang   = str_trim(rang)
  ) %>%
  filter(rang == "S", kmers_clade > 0) %>%
  slice_max(order_by = kmers_clade, n = 2, with_ties = FALSE) %>%
  arrange(desc(kmers_clade)) %>%
  mutate(
    sample      = name_sample,
    kmers       = kmers_clade,
    pct_kmers   = (kmers / sum(kmers)) * 100,
    specie      = factor(specie, levels = specie),
    y_mid       = rev(cumsum(rev(pct_kmers))) - (pct_kmers / 2)
  )

# Security if any specie is detected
if (nrow(data_detected) == 0) {
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Any specie detected:", name_sample)) +
    theme_void()
  ggsave(file_output, p_empty, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# --- 3. Graphic A : Number absolute of k-mer detected ---
p_kmers_raw <- ggplot(data_detected, aes(x = sample, y = kmers, fill = specie)) +
  geom_col(position = position_dodge(width = 0.6, preserve = "single"), width = 0.45) +
  geom_text(
    aes(label = format(kmers, big.mark = " ")),
    position = position_dodge(width = 0.6),
    vjust = -0.35, fontface = "bold", size = 3.8
  ) +
  scale_y_continuous(
    labels = function(x) format(x, big.mark = " "),
    expand = expansion(mult = c(0, 0.22))
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Absolute number of k-mers for the best species",
    subtitle = "K-mers assigned to the Species rank (Kraken2)",
    x        = NULL,
    y        = "Number of assigned k-mers",
    fill     = "Species detected"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text     = element_text(face = "italic")
  )

# --- 4. Graphic B : Proportion relative observed ---
p_part_obs <- ggplot(data_detected, aes(x = sample, y = pct_kmers, fill = specie)) +
  geom_col(width = 0.4) +
  geom_text(
    data = filter(data_detected, pct_kmers >= 5),
    aes(y = y_mid, label = sprintf("%.1f %%", pct_kmers)),
    color = "white", fontface = "bold", size = 3.8
  ) +
  geom_label_repel(
    data = filter(data_detected, pct_kmers < 5),
    aes(y = y_mid, label = sprintf("%.2f %%", pct_kmers)),
    color = "white", fontface = "bold", size = 3.5,
    nudge_x = 0.35, show.legend = FALSE,
    min.segment.length = 0, segment.color = "#FC8D62", segment.size = 0.8
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, " %"),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Observed relative proportion",
    subtitle = "Proportion of k-mers among the detected species",
    x        = NULL,
    y        = "Share of k-mers (%)",
    fill     = "Species detected"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 5. Assembly of plot et save ---
figure_kraken <- p_kmers_raw + p_part_obs + plot_layout(widths = c(1.2, 1))
ggsave(file_output, figure_kraken, width = 10, height = 7, dpi = 300)