# 05_parcimonia.R
# Linha de base por parcimônia (Sankoff, phangorn): número mínimo de passos e estados
# MPR nos nós focais, sob custos iguais (ganho = perda) e sob ganho 2x mais caro que a
# perda (aproxima o raciocínio de Jarvis et al. 2014, em que a rapinância é cara de ganhar).
# Pontas polimórficas entram como ambíguas. Só árvores sem comprimento relevante: a
# topologia de cada T* (versão datada sem fósseis e com F1).
# Execução, a partir de r/:  Rscript ../runs/rapinantes/ancestral/scripts/05_parcimonia.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
arv <- readRDS(file.path(base, "arvores", "arvores.rds"))
mat <- read.delim(file.path(base, "caracteres", "matriz_mk.tsv"), colClasses = "character")
mat$chave <- ifelse(mat$tipo == "fossil", mat$ponta, mat$family)
fams <- read.delim(file.path(base, "caracteres", "familias_pontas.tsv"))

grupos <- list(
  Neognathae_raiz = NULL,
  Telluraves = fams$family[fams$clade_mag7 == "Telluraves"],
  Australaves = fams$family[fams$order %in% c("Cariamiformes", "Falconiformes", "Psittaciformes", "Passeriformes")],
  Afroaves = fams$family[fams$order %in% c("Accipitriformes", "Strigiformes", "Coliiformes", "Leptosomatiformes",
                                           "Trogoniformes", "Bucerotiformes", "Coraciiformes", "Piciformes")],
  Eufalconimorphae = fams$family[fams$order %in% c("Falconiformes", "Psittaciformes", "Passeriformes")],
  Accipitrimorphae = fams$family[fams$order == "Accipitriformes"]
)
custos <- list(iguais = matrix(c(0, 1, 1, 0), 2, 2), ganho_2x = matrix(c(0, 2, 1, 0), 2, 2))
# phangorn::sankoff lê o custo como [destino, origem]: conferido num exemplo de brinquedo
# ((A=1,B=1),C=0), em que só esta orientação faz a raiz 1 (uma perda) vencer com ganho = 2
dimnames(custos$iguais) <- dimnames(custos$ganho_2x) <- list(c("0", "1"), c("0", "1"))

linhas <- list()
for (nome in grep("_(E_datada|F1)$", names(arv), value = TRUE)) {
  tr <- arv[[nome]]
  tr$edge.length <- NULL
  for (char in c("C1_T50", "C4_bico_estrito", "C5_pes")) {
    cod <- setNames(mat[[char]], mat$chave)[tr$tip.label]
    cod[cod == "0&1"] <- "?"
    dat <- phyDat(as.list(cod), type = "USER", levels = c("0", "1"), ambiguity = "?")
    for (cn in names(custos)) {
      passos <- sankoff(tr, dat, cost = custos[[cn]])
      anc <- ancestral.pars(tr, dat, type = "MPR", cost = custos[[cn]])
      for (g in names(grupos)) {
        no <- if (is.null(grupos[[g]])) Ntip(tr) + 1 else getMRCA(tr, grupos[[g]])
        p <- anc[[no]][1, ]
        linhas[[length(linhas) + 1]] <- data.frame(arvore = nome, caractere = char, custo = cn,
                                                   passos = passos, no = g,
                                                   MPR = paste(c("0", "1")[p > 0], collapse = "/"),
                                                   frac1 = p[2] / sum(p))
      }
    }
  }
}
res <- do.call(rbind, linhas)
write.table(res, file.path(base, "resultados", "parcimonia_mpr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w <- reshape(res[, c("arvore", "caractere", "custo", "passos", "no", "MPR")],
             idvar = c("arvore", "caractere", "custo", "passos"), timevar = "no", direction = "wide")
names(w) <- sub("^MPR\\.", "", names(w))
print(w[order(w$caractere, w$custo, w$arvore), ], row.names = FALSE)
