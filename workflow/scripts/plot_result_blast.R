library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# --- 1. Entrées / sorties Snakemake ---
fichier_blast   <- snakemake@input[["blast"]]
fichier_sortie  <- snakemake@output[["plot"]]
nom_echantillon <- paste(snakemake@wildcards[["species"]], snakemake@wildcards[["sex"]], sep = "_")

# Sécurité si le fichier TSV est vide (0 octet)
if (!file.exists(fichier_blast) || file.info(fichier_blast)$size == 0) {
  p_vide <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Fichier BLAST vide :", nom_echantillon)) +
    theme_void()
  ggsave(fichier_sortie, p_vide, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# --- 2. Lecture et extraction de la meilleure espèce de chaque contig ---
donnees_detectees <- read.delim(fichier_blast, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    echantillon = nom_echantillon,
    contig      = coalesce(str_extract(qseqid, "^NODE_\\d+"), qseqid),
    coverage    = as.numeric(str_extract(qseqid, "(?<=cov_)[0-9.]+")),
    espece      = sseqid %>% str_remove("_[^_]+$") %>% str_replace_all("_", " ")
  ) %>%
  arrange(desc(bitscore), desc(pident)) %>%
  group_by(echantillon, qseqid, contig, coverage) %>%
  slice(1) %>%
  ungroup() %>%
  filter(bitscore >= 1000) %>%
  arrange(desc(coverage)) %>%
  mutate(
    pct_reads = (coverage / sum(coverage)) * 100,
    espece    = factor(espece, levels = unique(espece)),
    contig    = factor(contig, levels = unique(contig))
  ) %>%
  arrange(espece, desc(coverage)) %>%
  mutate(
    y_mid = rev(cumsum(rev(pct_reads))) - (pct_reads / 2)
  )

# Sécurité si aucun hit ne passe bitscore >= 800
if (nrow(donnees_detectees) == 0) {
  p_vide <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Aucun hit BLAST (bitscore >= 800) :", nom_echantillon)) +
    theme_void()
  ggsave(fichier_sortie, p_vide, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# Palette extensible si plus de 8 espèces sont détectées (limite de Set2)
palette_especes <- if (nlevels(donnees_detectees$espece) <= 8) {
  scale_fill_brewer(palette = "Set2")
} else {
  scale_fill_viridis_d(option = "D")
}

# --- 3. Graphique A : Couverture brute + pident par contig ---
p_cov_brute <- ggplot(donnees_detectees, aes(x = echantillon, y = coverage, fill = espece, group = contig)) +
  geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.55) +
  geom_text(
    aes(label = sprintf("%s\n%.1f× (%.1f%%)", contig, coverage, pident)),
    position = position_dodge(width = 0.7),
    vjust = -0.25, lineheight = 0.95, fontface = "bold", size = 3.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  palette_especes +
  labs(
    title    = "Couverture brute et identité par contig",
    subtitle = "Profondeur de couverture (×) et identité BLAST de la meilleure espèce",
    x        = NULL,
    y        = "Couverture moyenne (×)",
    fill     = "Espèce détectée"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text     = element_text(face = "italic")
  )

# --- 4. Graphique B : Part relative observée ---
p_part_obs <- ggplot(donnees_detectees, aes(x = echantillon, y = pct_reads, fill = espece, group = contig)) +
  geom_col(width = 0.4, color = "white", linewidth = 0.4) +
  geom_text(
    data = filter(donnees_detectees, pct_reads >= 5),
    aes(y = y_mid, label = sprintf("%.1f %%", pct_reads)),
    color = "white", fontface = "bold", size = 3.8
  ) +
  geom_label_repel(
    data = filter(donnees_detectees, pct_reads < 5),
    aes(y = y_mid, label = sprintf("%.2f %%", pct_reads)),
    color = "white", fontface = "bold", size = 3.5,
    nudge_x = 0.35, show.legend = FALSE,
    min.segment.length = 0, segment.color = "#FC8D62", segment.size = 0.8
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, " %"),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  palette_especes +
  labs(
    title    = "Proportion relative observée",
    subtitle = "Part de couverture entre les contigs détectés",
    x        = NULL,
    y        = "Part de la couverture (%)",
    fill     = "Espèce détectée"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 5. Assemblage et sauvegarde ---
figure_detection <- p_cov_brute + p_part_obs + plot_layout(widths = c(1.2, 1))
ggsave(fichier_sortie, figure_detection, width = 10, height = 7, dpi = 300)