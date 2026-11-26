library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

# --- 1. Fichier BLAST unique (avec en-têtes) ---
fichier_blast <- "blast_output_0-30.tsv"  # Remplace par le nom exact de ton fichier
nom_echantillon <- tools::file_path_sans_ext(basename(fichier_blast))

# --- 2. Lecture et extraction des meilleurs contigs ---
donnees_detectees <- read.delim(fichier_blast, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    echantillon = nom_echantillon,
    contig      = str_extract(qseqid, "^NODE_\\d+"),
    coverage    = as.numeric(str_extract(qseqid, "(?<=cov_)[0-9.]+")),
    espece      = sseqid %>% str_remove("_[^_]+$") %>% str_replace_all("_", " ")
  ) %>%
  arrange(desc(bitscore), desc(pident)) %>%
  group_by(echantillon, qseqid, contig, coverage) %>%
  slice(1) %>%
  ungroup() %>%
  filter(bitscore >= 800) %>%                 # Filtre des faibles bitscores
  slice_max(order_by = coverage, n = 2) %>%    # Prend au maximum les 2 meilleurs contigs
  arrange(desc(coverage)) %>%                  # Le plus couvert en 1er (cible), le 2e ensuite (contaminant)
  mutate(
    role      = if_else(row_number() == 1, "Contig principal (Cible)", "2e contig détecté (Contaminant)"),
    pct_reads = (coverage / sum(coverage)) * 100
  )

# --- 3. Graphique A : Couverture brute + pident en dessous ---
p_cov_brute <- ggplot(donnees_detectees, aes(x = echantillon, y = coverage, fill = espece)) +
  geom_col(position = position_dodge(width = 0.6, preserve = "single"), width = 0.45) +
  geom_text(
    aes(label = sprintf("%.1f×\n(%.1f %% id.)", coverage, pident)),
    position = position_dodge(width = 0.6),
    vjust = -0.25, lineheight = 0.95, fontface = "bold", size = 3.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Couverture brute et identité des meilleurs contigs",
    subtitle = "Profondeur de couverture (×) et pourcentage d'identité BLAST (pident)",
    x        = NULL, y = "Couverture moyenne (×)", fill = "Espèce détectée"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text = element_text(face = "italic")
  )

# --- 4. Graphique B : Part relative observée ---
p_part_obs <- ggplot(donnees_detectees, aes(x = echantillon, y = pct_reads, fill = espece)) +
  geom_col(width = 0.4, color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = ifelse(pct_reads >= 3, sprintf("%.1f %%", pct_reads), "")),
    position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 3.8
  ) +
  scale_y_continuous(labels = function(x) paste0(x, " %"), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Proportion relative observée",
    subtitle = "Part de couverture entre les espèces détectées",
    x        = NULL, y = "Part de la couverture (%)", fill = "Espèce détectée"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# library(ggrepel)
# 
# # Calcul de la position verticale exacte de chaque segment empilé
# donnees_detectees <- donnees_detectees %>%
#   mutate(
#     espece = factor(espece, levels = espece),
#     y_mid  = rev(cumsum(rev(pct_reads))) - (pct_reads / 2)
#   )
# 
# # --- 4. Graphique B : Part relative observée ---
# p_part_obs <- ggplot(donnees_detectees, aes(x = echantillon, y = pct_reads, fill = espece)) +
#   # Suppression de color = "white" pour laisser apparaître la fine bande orange
#   geom_col(width = 0.4) +
#   # Pourcentage au centre pour les espèces majoritaires (>= 5 %)
#   geom_text(
#     data = filter(donnees_detectees, pct_reads >= 5),
#     aes(y = y_mid, label = sprintf("%.1f %%", pct_reads)),
#     color = "white", fontface = "bold", size = 3.8
#   ) +
#   # Étiquette orange déportée avec flèche pour les faibles proportions (< 5 %)
#   geom_label_repel(
#     data = filter(donnees_detectees, pct_reads < 5),
#     aes(y = y_mid, label = sprintf("%.2f %%", pct_reads)),
#     color = "white", fontface = "bold", size = 3.5,
#     nudge_x = 0.35, show.legend = FALSE,
#     min.segment.length = 0, segment.color = "#FC8D62", segment.size = 0.8
#   ) +
#   scale_y_continuous(
#     labels = function(x) paste0(x, " %"),
#     expand = expansion(mult = c(0.01, 0.04))
#   ) +
#   scale_fill_brewer(palette = "Set2") +
#   labs(
#     title    = "Proportion relative observée",
#     subtitle = "Part de couverture entre les espèces détectées",
#     x        = NULL, y = "Part de la couverture (%)", fill = "Espèce détectée"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(legend.position = "none")

# --- 5. Assemblage final ---
figure_detection <- p_cov_brute + p_part_obs + plot_layout(widths = c(1.2, 1))
print(figure_detection)

# Sauvegarde optionnelle :
ggsave("detection_contamination_blast.png", figure_detection, width = 10, height = 7, dpi = 300)