# 02_arvores.R
# Árvores em nível de família para a reconstrução de estado ancestral do hábito
# rapinante em Telluraves, sob seis topologias:
#   T1   Stiller et al. 2024 (Nature), árvore principal 63K ASTRAL, datada (MCMCtree)
#   T1b  Stiller et al. 2024, árvore ASTRAL 94K (mesma topologia de Telluraves em 128K/153K/159K)
#   T1c  Stiller et al. 2024, árvore ASTRAL de éxons (posições 1+2 do códon), "clado das rapinantes"
#   T2   Prum et al. 2015 (Nature), espinha dorsal da árvore datada de 198 espécies
#   T3   Jarvis et al. 2014 (Science), TENT ExaML datada (48 espécies)
#   T4   Kuhl et al. 2021 (MBE), Fig. 2 (árvore datada em nível de ordem), codificada à mão
# T1, T1b, T1c: poda direta das árvores de espécies de Stiller para uma espécie por família.
# T2, T3, T4: híbridas -- espinha dorsal entre ordens da fonte; relações entre famílias
# dentro de cada unidade (ordem ou clado-externo) tomadas de T1.
# Comprimentos: (a) idades de Stiller 2024 -- nós cujo conjunto de pontas existe em T1
# recebem a idade de T1, os demais são interpolados por espaçamento uniforme (BLADJ);
# (b) Grafen (rho = 1), escalado para a idade da raiz de T1.
# Fósseis (F1, F2) enxertados só nas árvores datadas.
#
# Execução, a partir de r/:  Rscript ../runs/rapinantes/ancestral/scripts/02_arvores.R

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(treeio)
})

base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
dl <- file.path(base, "downloads", "erda")
out <- file.path(base, "arvores")
dir.create(out, showWarnings = FALSE)

fams <- read.delim(file.path(base, "caracteres", "familias_pontas.tsv"))
rep2fam <- setNames(fams$family, fams$ponta)
ord_of <- setNames(fams$order, fams$family)
clade_of <- setNames(fams$clade_mag7, fams$family)
stopifnot(!anyDuplicated(fams$ponta))

# ---------------------------------------------------------------- T1 datada
st <- read.mcmctree(file.path(dl, "04_timetree__main_alternative_mod.tre"))@phylo
stopifnot(all(fams$ponta %in% st$tip.label))
idade_coroa <- function(tr, tips) {
  h <- max(node.depth.edgelength(tr))
  nd <- getMRCA(tr, tips)
  h - node.depth.edgelength(tr)[nd]
}
ann <- read.delim(file.path(dl, "07_scripts_figures__Fig1__annotations.txt"))
coroa_falconidae <- idade_coroa(st, ann$taxa[ann$family == "Falconidae"])
coroa_cariamidae <- idade_coroa(st, ann$taxa[ann$family == "Cariamidae"])

T1 <- keep.tip(st, fams$ponta)
T1$tip.label <- unname(rep2fam[T1$tip.label])
T1$node.label <- NULL
stopifnot(is.ultrametric(T1, tol = 1e-6), is.binary(T1))
raiz_T1 <- max(node.depth.edgelength(T1))

# ---------------------------------------------------------------- utilidades
galloanserae <- fams$family[fams$clade_mag7 == "Galloanseres"]

poda_astral <- function(arquivo) {
  t <- read.tree(file.path(dl, arquivo))
  stopifnot(all(fams$ponta %in% t$tip.label))
  t <- keep.tip(t, fams$ponta)
  t$tip.label <- unname(rep2fam[t$tip.label])
  t <- root(unroot(t), outgroup = galloanserae, resolve.root = TRUE)
  t$edge.length <- NULL
  t$node.label <- NULL
  if (!is.binary(t)) t <- multi2di(t, random = FALSE)
  t
}

idades_no <- function(tr) {
  h <- max(node.depth.edgelength(tr))
  h - node.depth.edgelength(tr)
}

# BLADJ: idade de T1 para nós cujo conjunto de pontas é clado em T1; demais nós
# espaçados uniformemente entre o ancestral fixo mais próximo e o descendente fixo mais antigo.
bladj <- function(topo, ref) {
  ref_age <- idades_no(ref)
  nt <- Ntip(topo); nn <- topo$Nnode
  age <- rep(NA_real_, nt + nn)
  age[1:nt] <- 0
  fixo <- rep(FALSE, nt + nn); fixo[1:nt] <- TRUE
  desc <- prop.part(topo)
  for (i in seq_len(nn)) {
    tips <- topo$tip.label[desc[[i]]]
    if (length(tips) == Ntip(ref) || is.monophyletic(ref, tips)) {
      age[nt + i] <- ref_age[getMRCA(ref, tips)]
      fixo[nt + i] <- TRUE
    }
  }
  raiz <- nt + 1
  stopifnot(fixo[raiz])
  filhos <- function(n) topo$edge[topo$edge[, 1] == n, 2]
  # L: idade do descendente fixo mais antigo alcançado só por nós não fixos; d: profundidade da cadeia
  L <- rep(NA_real_, nt + nn); d <- rep(0L, nt + nn)
  # ordem dos nós internos por profundidade em número de ramos: raiz primeiro
  tu <- topo; tu$edge.length <- rep(1, nrow(tu$edge))
  prof <- node.depth.edgelength(tu)
  internos <- (nt + 1):(nt + nn)
  ordem <- internos[order(prof[internos])]
  for (n in rev(ordem)) {
    ch <- filhos(n)
    Lc <- sapply(ch, function(c) if (fixo[c]) age[c] else L[c])
    dc <- sapply(ch, function(c) if (fixo[c]) 0L else d[c])
    L[n] <- max(Lc)
    d[n] <- 1L + max(dc)
  }
  for (n in ordem) {
    for (c in filhos(n)) {
      if (!fixo[c]) {
        P <- age[n]
        age[c] <- P - (P - L[c]) / (d[c] + 1)
      }
    }
  }
  bl <- age[topo$edge[, 1]] - age[topo$edge[, 2]]
  if (any(bl <= 0)) stop("ramo de comprimento <= 0 após BLADJ: ", sum(bl <= 0))
  topo$edge.length <- bl
  attr(topo, "nos_fixos") <- sum(fixo[(nt + 1):(nt + nn)])
  attr(topo, "nos_interpolados") <- sum(!fixo[(nt + 1):(nt + nn)])
  topo
}

grafen <- function(topo) {
  t <- compute.brlen(topo, method = "Grafen", power = 1)
  t$edge.length <- t$edge.length * raiz_T1 / max(node.depth.edgelength(t))
  t
}

# ---------------------------------------------------------------- unidades para as híbridas
unidade <- function(f) {
  o <- ord_of[[f]]; cl <- clade_of[[f]]
  if (cl == "Galloanseres") return("GALLOANSERAE")
  if (cl == "Aequornithes") return("AEQUORNITHES")
  if (o == "Strisores") return("STRISORES")
  if (o == "Charadriiformes") return("CHARADRIIFORMES")
  toupper(o)
}
un <- sapply(fams$family, unidade)
sub_unidade <- function(u) {
  fs <- names(un)[un == u]
  if (length(fs) == 1) return(fs)
  s <- keep.tip(T1, fs); s$edge.length <- NULL
  sub(";$", "", write.tree(s))
}
monta <- function(espinha) {
  tok <- unique(regmatches(espinha, gregexpr("[A-Z]{4,}", espinha))[[1]])
  stopifnot(setequal(tok, unique(un)))
  for (u in tok) espinha <- gsub(paste0("\\b", u, "\\b"), sub_unidade(u), espinha, perl = TRUE)
  t <- read.tree(text = espinha)
  stopifnot(setequal(t$tip.label, fams$family), is.binary(t))
  t
}

tell_prum <- "(ACCIPITRIFORMES,((CARIAMIFORMES,(FALCONIFORMES,(PSITTACIFORMES,PASSERIFORMES))),(STRIGIFORMES,(COLIIFORMES,(LEPTOSOMATIFORMES,(TROGONIFORMES,(BUCEROTIFORMES,(CORACIIFORMES,PICIFORMES))))))))"
tell_jk   <- "((ACCIPITRIFORMES,(STRIGIFORMES,(COLIIFORMES,(LEPTOSOMATIFORMES,(TROGONIFORMES,(BUCEROTIFORMES,(CORACIIFORMES,PICIFORMES))))))),(CARIAMIFORMES,(FALCONIFORMES,(PSITTACIFORMES,PASSERIFORMES))))"
espinhas <- list(
  # Prum 2015: (Strisores,(Columbaves,(Gruiformes,((Aequornithes+Phaeth),(Charadrii+Mirandornithes)),(Opisthocomus,Telluraves)))) podada
  T2 = sprintf("(GALLOANSERAE,(STRISORES,((((PTEROCLIFORMES,MESITORNITHIFORMES),COLUMBIFORMES),((CUCULIFORMES,OTIDIFORMES),MUSOPHAGIFORMES)),(%s,(AEQUORNITHES,CHARADRIIFORMES)))));", tell_prum),
  # Jarvis 2014 TENT: Columbea vs Passerea; Telluraves + (Aequornithes+Phaeth); + (Charadrii+Grui, Opistho); + (Otidimorphae, Strisores)
  T3 = sprintf("(GALLOANSERAE,(((PTEROCLIFORMES,MESITORNITHIFORMES),COLUMBIFORMES),(((%s,AEQUORNITHES),CHARADRIIFORMES),(((MUSOPHAGIFORMES,OTIDIFORMES),CUCULIFORMES),STRISORES))));", tell_jk),
  # Kuhl 2021 Fig. 2: (Mirandornithes,(Basal Landbirds,(Aquatic&Semiaquatic, Telluraves)))
  T4 = sprintf("(GALLOANSERAE,((((((PTEROCLIFORMES,MESITORNITHIFORMES),(COLUMBIFORMES,CUCULIFORMES)),(MUSOPHAGIFORMES,OTIDIFORMES)),STRISORES),CHARADRIIFORMES),(AEQUORNITHES,%s)));", tell_jk)
)

topos <- list(
  T1 = { t <- T1; t$edge.length <- NULL; t },
  T1b = poda_astral("03_species_trees__94K.tre"),
  T1c = poda_astral("03_species_trees__exon_c12.tre"),
  T2 = monta(espinhas$T2),
  T3 = monta(espinhas$T3),
  T4 = monta(espinhas$T4)
)

# ---------------------------------------------------------------- fósseis
fos <- read.delim(file.path(base, "caracteres", "fosseis_codificacao_manual.tsv"))
psitpas <- fams$family[fams$order %in% c("Psittaciformes", "Passeriformes")]
alvo_no <- function(tr, pos) {
  switch(pos,
    "tronco:Strigiformes" = list(no = getMRCA(tr, c("Tytonidae", "Strigidae")), min_idade = 0),
    "tronco:Psittacopasserae" = list(no = getMRCA(tr, psitpas), min_idade = 0),
    "tronco:Falconidae" = list(no = which(tr$tip.label == "Falconidae"), min_idade = coroa_falconidae),
    "tronco:Cariamidae" = list(no = which(tr$tip.label == "Cariamidae"), min_idade = coroa_cariamidae),
    stop("posição desconhecida: ", pos)
  )
}
enxerta <- function(tr, cenario) {
  log <- data.frame()
  for (i in seq_len(nrow(fos))) {
    pos <- fos[[paste0("posicao_", cenario)]][i]
    a <- alvo_no(tr, pos)
    ida <- idades_no(tr)
    pai <- tr$edge[tr$edge[, 2] == a$no, 1]
    c_idade <- max(ida[a$no], a$min_idade)
    f <- fos$idade_Ma[i]
    piso <- max(c_idade, f, fos$idade_minima_linhagem_Ma[i])
    if (ida[pai] <= piso) stop("fóssil ", fos$ponta[i], " mais antigo que o nó-pai do ramo-alvo")
    anc <- (ida[pai] + piso) / 2
    # bind.tip: 'position' é medido do nó 'where' para trás (rumo à raiz)
    tr <- bind.tip(tr, fos$ponta[i], edge.length = anc - f, where = a$no, position = anc - ida[a$no])
    log <- rbind(log, data.frame(cenario = cenario, fossil = fos$ponta[i], posicao = pos,
                                 idade_ponta = f, idade_enxerto = round(anc, 3),
                                 idade_no_pai = round(ida[pai], 3), idade_no_filho = round(ida[a$no], 3)))
  }
  # conferência: a idade da ponta fóssil reproduz a declarada
  ida <- idades_no(tr)
  for (i in seq_len(nrow(fos))) {
    k <- which(tr$tip.label == fos$ponta[i])
    stopifnot(abs(ida[k] - fos$idade_Ma[i]) < 1e-6)
  }
  attr(tr, "log_enxerto") <- log
  tr
}

# ---------------------------------------------------------------- monta tudo
arv <- list(); resumo <- data.frame(); enx <- data.frame()
for (nm in names(topos)) {
  tp <- topos[[nm]]
  dat <- if (nm == "T1") T1 else bladj(tp, T1)
  arv[[paste0(nm, "_E_datada")]] <- dat
  arv[[paste0(nm, "_E_grafen")]] <- grafen(tp)
  for (cen in c("F1", "F2")) {
    f <- enxerta(dat, cen)
    enx <- rbind(enx, cbind(topologia = nm, attr(f, "log_enxerto")))
    attr(f, "log_enxerto") <- NULL
    arv[[paste0(nm, "_", cen)]] <- f
  }
  resumo <- rbind(resumo, data.frame(
    topologia = nm, n_pontas = Ntip(tp),
    RF_vs_T1 = dist.topo(unroot(tp), unroot(topos$T1)),
    nos_fixos_T1 = if (nm == "T1") tp$Nnode else attr(dat, "nos_fixos"),
    nos_interpolados = if (nm == "T1") 0 else attr(dat, "nos_interpolados"),
    raiz_Ma = round(max(node.depth.edgelength(dat)), 3)
  ))
}
for (nm in names(arv)) write.tree(arv[[nm]], file.path(out, paste0(nm, ".tre")))
saveRDS(arv, file.path(out, "arvores.rds"))
write.table(resumo, file.path(out, "resumo_topologias.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(enx, file.path(out, "enxertos_fosseis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(sprintf("coroa Falconidae (Stiller 2024, Herpetotheres-Falco) = %.3f Ma", coroa_falconidae),
             sprintf("coroa Cariamidae (Stiller 2024, Cariama-Chunga) = %.3f Ma", coroa_cariamidae),
             sprintf("raiz T1 (Galloanserae vs Neoaves) = %.3f Ma", raiz_T1)),
           file.path(out, "idades_referencia.txt"))

# espinha de Telluraves por topologia (ordens), para conferência
colapsa <- function(tr) {
  t2 <- tr; t2$edge.length <- NULL
  tel <- fams$family[fams$clade_mag7 == "Telluraves"]
  t2 <- keep.tip(t2, tel)
  keep <- c()
  for (o in unique(ord_of[t2$tip.label])) {
    sp <- t2$tip.label[ord_of[t2$tip.label] == o]
    keep <- c(keep, sp[1])
  }
  t3 <- keep.tip(t2, keep); t3$tip.label <- ord_of[t3$tip.label]
  write.tree(ladderize(t3))
}
esp <- sapply(topos, colapsa)
writeLines(paste(names(esp), esp, sep = "\t"), file.path(out, "espinha_telluraves_por_topologia.tsv"))
print(resumo)
print(enx)
cat(paste(names(esp), esp, sep = "  "), sep = "\n")
