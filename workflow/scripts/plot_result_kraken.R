library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# --- 1. Gestion des entrées/sorties (Snakemake ou exécution manuelle) ---
if (exists("snakemake")) {
  fichier_kraken  <- snakemake@input[[1]]
  fichier_sortie  <- snakemake@output[[1]]
  nom_echantillon <- if (!is.null(snakemake@wildcards[["sample"]])) {
    snakemake@wildcards[["sample"]]
  } else {
    tools::file_path_sans_ext(basename(fichier_kraken))
  }
} else {
  # Mode test local dans RStudio
  fichier_kraken  <- "mix_5_30_conf015.report"
  fichier_sortie  <- "detection_test_contamination_kraken.png"
  nom_echantillon <- tools::file_path_sans_ext(basename(fichier_kraken))
}

# --- 2. Lecture du rapport Kraken2 et extraction des 2 principales espèces ---
donnees_detectees <- read.delim(
  fichier_kraken, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
  col.names = c("pct", "kmers_clade", "kmers_taxon", "rang", "taxid", "nom")
) %>%
  mutate(
    espece = str_trim(nom),
    rang   = str_trim(rang)
  ) %>%
  filter(rang == "S", kmers_clade > 0) %>%
  slice_max(order_by = kmers_clade, n = 2, with_ties = FALSE) %>%
  arrange(desc(kmers_clade)) %>%
  mutate(
    echantillon = nom_echantillon,
    kmers       = kmers_clade,
    pct_kmers   = (kmers / sum(kmers)) * 100,
    espece      = factor(espece, levels = espece),
    # Position verticale exacte du centre de chaque segment empilé
    y_mid       = rev(cumsum(rev(pct_kmers))) - (pct_kmers / 2)
  )

# Sécurité si aucune espèce n'est détectée dans un échantillon vide
if (nrow(donnees_detectees) == 0) {
  p_vide <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Aucune espèce détectée :", nom_echantillon)) +
    theme_void()
  ggsave(fichier_sortie, p_vide, width = 10, height = 7, dpi = 300)
  quit(save = "no", status = 0)
}

# --- 3. Graphique A : Nombre absolu de k-mers détectés ---
p_kmers_brute <- ggplot(donnees_detectees, aes(x = echantillon, y = kmers, fill = espece)) +
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
    title    = "Nombre absolu de k-mers des meilleures espèces",
    subtitle = "K-mers assignés au rang Espèce (Kraken2)",
    x        = NULL,
    y        = "Nombre de k-mers assignés",
    fill     = "Espèce détectée"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text     = element_text(face = "italic")
  )

# --- 4. Graphique B : Proportion relative observée ---
p_part_obs <- ggplot(donnees_detectees, aes(x = echantillon, y = pct_kmers, fill = espece)) +
  geom_col(width = 0.4) +
  # Étiquette centrée pour les proportions >= 5 %
  geom_text(
    data = filter(donnees_detectees, pct_kmers >= 5),
    aes(y = y_mid, label = sprintf("%.1f %%", pct_kmers)),
    color = "white", fontface = "bold", size = 3.8
  ) +
  # Étiquette orange déportée avec flèche pour les faibles proportions (< 5 %)
  geom_label_repel(
    data = filter(donnees_detectees, pct_kmers < 5),
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
    title    = "Proportion relative observée",
    subtitle = "Part des k-mers entre les espèces détectées",
    x        = NULL,
    y        = "Part des k-mers (%)",
    fill     = "Espèce détectée"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 5. Assemblage et sauvegarde ---
figure_kraken <- p_kmers_brute + p_part_obs + plot_layout(widths = c(1.2, 1))
ggsave(fichier_sortie, figure_kraken, width = 10, height = 7, dpi = 300)