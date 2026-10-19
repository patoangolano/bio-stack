# Figuras das arvores ML (IQ-TREE) com fosseis destacados e UFBoot nos nos.
#   cd r; Rscript ../scripts/rapinantes/fosseis/p08_figuras_arvores.R
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
})
BASE <- normalizePath(file.path("..", "runs", "rapinantes", "fosseis"), winslash = "/")
FIG <- file.path(BASE, "figuras")
dir.create(FIG, showWarnings = FALSE)
EXT <- c(K19 = "Eudromia_elegans", M11 = "Tinamidae", M05 = "Tinamidae", M16 = "Anhimidae",
         VZ26 = "Anhimidae", MK22 = "Anhimidae", MC03 = "Apsaravis", C25 = "Gansus_yumenensis",
         LZ07 = "Lithornis", D21 = "Meleagris_gallopavo")
CORES <- c(Accipitriformes = "#b2182b", Strigiformes = "#762a83", Falconiformes = "#e08214",
           Cariamiformes = "#8c510a", Psittaciformes = "#1b7837", Passeriformes = "#2166ac",
           Coraciimorphae = "#5aae61", outros = "#7f7f7f", fossil = "#000000")

desenhar <- function(id, modo) {
  f <- file.path(BASE, "arvores", id, paste0(modo, ".treefile"))
  if (!file.exists(f)) return(invisible(NULL))
  tr <- read.tree(f)
  tr <- root(tr, EXT[[id]], resolve.root = TRUE)
  tr <- ladderize(tr)
  # LZ07 nao tem suporte nos nos (UFBoot interrompido)
  ufb <- if (is.null(tr$node.label)) rep("", Nnode(tr)) else sub("^.*/", "", tr$node.label)
  ufb[!grepl("^[0-9.]+$", ufb)] <- ""
  tr$node.label <- ufb
  tx <- read.delim(file.path(BASE, "matrizes", id, "taxa.tsv"), stringsAsFactors = FALSE)
  meta <- data.frame(label = tx$terminal, clado = tx$clado, tipo = tx$tipo)
  meta$rot <- ifelse(meta$tipo == "fossil", paste0("† ", gsub("_", " ", meta$label)), gsub("_", " ", meta$label))
  n <- Ntip(tr)
  p <- ggtree(tr, linewidth = 0.3) %<+% meta +
    geom_tiplab(aes(label = rot, color = clado, fontface = ifelse(tipo == "fossil", "bold", "plain")),
                size = 2.2, show.legend = FALSE) +
    geom_nodelab(aes(label = label), size = 1.7, hjust = 1.2, vjust = -0.4, color = "grey30") +
    scale_color_manual(values = CORES, na.value = "grey40") +
    theme_tree2() + xlab("substituicoes por carater (Mk)") +
    ggtitle(sprintf("%s, %s: ML (IQ-TREE 3, Mk+ASC+G4 particionado); nos = UFBoot (%%), n = %d terminais", id,
                    ifelse(modo == "livre", "morfologia sozinha", "atuais presos a Stiller et al. 2024"), n)) +
    theme(plot.title = element_text(size = 7))
  p <- p + xlim(0, max(node.depth.edgelength(tr)) * 1.6)
  ggsave(file.path(FIG, sprintf("arvore_%s_%s.png", id, modo)), p, width = 7, height = max(3, 0.13 * n + 1),
         dpi = 220, limitsize = FALSE)
}
for (id in c("K19", "M11", "M05", "VZ26", "MK22", "M16", "MC03", "C25", "D21", "LZ07")) {
  for (modo in c("livre", "esqueleto")) {
    tryCatch(desenhar(id, modo), error = function(e) message(id, " ", modo, ": ", conditionMessage(e)))
  }
}
cat("ggtree", as.character(packageVersion("ggtree")), "\n")
