# 04_resumo_figuras.R
# Tabelas-resumo e figuras a partir das saídas de 03_mk.R.
# Execução, a partir de r/:  Rscript ../runs/rapinantes/ancestral/scripts/04_resumo_figuras.R

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
})

base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
res_dir <- file.path(base, "resultados")
fig_dir <- file.path(base, "figuras")
dir.create(fig_dir, showWarnings = FALSE)

mg <- read.delim(file.path(res_dir, "marginais_nos.tsv"))
aj <- read.delim(file.path(res_dir, "ajustes_modelos.tsv"))
sm <- read.delim(file.path(res_dir, "simmap_transicoes.tsv"))
sn <- read.delim(file.path(res_dir, "simmap_estado_nos.tsv"))
arv <- readRDS(file.path(base, "arvores", "arvores.rds"))
mat <- read.delim(file.path(base, "caracteres", "matriz_mk.tsv"), colClasses = "character")
mat$chave <- ifelse(mat$tipo == "fossil", mat$ponta, mat$family)
fams <- read.delim(file.path(base, "caracteres", "familias_pontas.tsv"))

TOPO <- c("T1", "T1b", "T1c", "T2", "T3", "T4")
NOS <- c("Neognathae_raiz", "Telluraves", "Australaves", "Afroaves", "Eufalconimorphae", "Accipitrimorphae",
         "pai_Accipitrimorphae", "Strigiformes", "pai_Strigiformes")
fmt <- function(p, mono = TRUE) ifelse(is.na(p), "NA", ifelse(mono, sprintf("%.2f", p), sprintf("(%.2f)", p)))

# P(rapinante): binários -> P1; C6 -> P(estado 2 = síndrome completa)
mg$Prap <- ifelse(mg$caractere == "C6_sindrome", mg$P2, mg$P1)
mg$Pqualquer <- ifelse(mg$caractere == "C6_sindrome", mg$P1 + mg$P2, mg$P1)

# ---------------------------------------------------------------- tabela principal
sel <- subset(mg, caractere == "C1_T50" & ramos == "datada" & fosseis == "nenhum")
linhas <- list()
for (r in c("fitzjohn", "zero")) for (t in TOPO) {
  s <- subset(sel, topologia == t & raiz == r)
  melhor <- unique(s$modelo[s$melhor])
  a <- subset(aj, caractere == "C1_T50" & ramos == "datada" & fosseis == "nenhum" & topologia == t & raiz == r)
  pesos <- paste(sprintf("%s=%.2f", a$modelo, a$peso), collapse = "; ")
  d <- data.frame(raiz = r, topologia = t, melhor_modelo = melhor, pesos_AICc = pesos)
  for (n in NOS) {
    x <- subset(s, modelo == "media_AICc" & no == n)
    d[[n]] <- fmt(x$Prap, x$monofiletico)
  }
  x <- subset(s, modelo == "media_AICc" & no == "pai_Accipitrimorphae")
  d$identidade_pai_Accipitrimorphae <- x$identidade
  x <- subset(s, modelo == "media_AICc" & no == "pai_Strigiformes")
  d$identidade_pai_Strigiformes <- x$identidade
  linhas[[length(linhas) + 1]] <- d
}
tp <- do.call(rbind, linhas)
write.table(tp, file.path(res_dir, "tabela_principal_C1_T50.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# P1 por modelo (não só a média), C1_T50, árvores datadas sem fósseis
pm <- subset(mg, caractere == "C1_T50" & ramos == "datada" & fosseis == "nenhum" & no %in% NOS)
pm_w <- reshape(pm[, c("topologia", "raiz", "modelo", "no", "Prap")], idvar = c("topologia", "raiz", "modelo"),
                timevar = "no", direction = "wide")
names(pm_w) <- sub("^Prap\\.", "", names(pm_w))
pm_w <- pm_w[order(pm_w$raiz, match(pm_w$topologia, TOPO), pm_w$modelo), c("raiz", "topologia", "modelo", NOS)]
write.table(pm_w, file.path(res_dir, "tabela_C1_T50_por_modelo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA")

# ---------------------------------------------------------------- sensibilidade: caracteres, ramos, fósseis
sens <- subset(mg, modelo == "media_AICc" & no %in% c("Telluraves", "Afroaves", "Australaves", "pai_Accipitrimorphae",
                                                      "Accipitrimorphae", "Strigiformes"))
sens$cenario <- paste(sens$ramos, sens$fosseis, sep = "/")
sw <- reshape(sens[, c("caractere", "raiz", "cenario", "topologia", "no", "Prap")],
              idvar = c("caractere", "raiz", "cenario", "topologia"), timevar = "no", direction = "wide")
names(sw) <- sub("^Prap\\.", "", names(sw))
sw <- sw[order(sw$caractere, sw$raiz, sw$cenario, match(sw$topologia, TOPO)), ]
write.table(sw, file.path(res_dir, "tabela_sensibilidade_media_AICc.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")

# melhor modelo por cenário
best <- subset(aj, dAICc == 0)
bt <- aggregate(modelo ~ caractere + raiz + ramos + fosseis, data = best, FUN = function(v) {
  tb <- table(v); paste(sprintf("%s:%d", names(tb), tb), collapse = " ")
})
write.table(bt, file.path(res_dir, "tabela_melhor_modelo_contagem.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- transições
tr_tab <- subset(sm, (escopo == "arvore_toda" | escopo == "coroa_Telluraves"))
tr_tab$resumo <- sprintf("%.1f [%.0f-%.0f]", tr_tab$media, tr_tab$q025, tr_tab$q975)
write.table(tr_tab[, c("topologia", "ramos", "fosseis", "caractere", "raiz", "modelo_mapa", "semente", "n_mapas",
                       "escopo", "transicao", "media", "mediana", "q025", "q975", "resumo")],
            file.path(res_dir, "tabela_transicoes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- figura 1: espinhas de Telluraves
rapt <- c("Accipitriformes", "Strigiformes", "Falconiformes")
png(file.path(fig_dir, "fig1_espinhas_telluraves.png"), width = 2400, height = 1600, res = 200)
par(mfrow = c(2, 3), mar = c(3, 1, 3, 1))
leg <- c(T1 = "T1 Stiller 2024 (63K ASTRAL)", T1b = "T1b Stiller 2024 (94K-159K ASTRAL)",
         T1c = "T1c Stiller 2024 (éxons c12)", T2 = "T2 Prum 2015", T3 = "T3 Jarvis 2014 (TENT)",
         T4 = "T4 Kuhl 2021 (3'-UTR)")
for (t in TOPO) {
  tr <- arv[[paste0(t, "_E_datada")]]
  tel <- fams$family[fams$clade_mag7 == "Telluraves"]
  tt <- keep.tip(tr, tel)
  ordem <- setNames(fams$order, fams$family)
  keep <- sapply(unique(ordem[tt$tip.label]), function(o) tt$tip.label[ordem[tt$tip.label] == o][1])
  t3 <- keep.tip(tt, keep); t3$tip.label <- unname(ordem[t3$tip.label])
  t3 <- ladderize(t3)
  cores <- ifelse(t3$tip.label %in% rapt, "firebrick", "grey20")
  plot(t3, tip.color = cores, cex = 0.9, main = leg[t], edge.width = 1.5, label.offset = 0.3)
  axisPhylo(cex.axis = 0.7)
}
dev.off()

# ---------------------------------------------------------------- figura 2: mapa de calor P(rapinante)
hm <- subset(mg, caractere %in% c("C1_T50", "C4_bico_estrito", "C5_pes", "C6_sindrome") & ramos == "datada" &
               fosseis == "nenhum" & modelo == "media_AICc" &
               no %in% c("Telluraves", "Australaves", "Afroaves", "Eufalconimorphae", "Accipitrimorphae",
                         "pai_Accipitrimorphae"))
hm$no <- factor(hm$no, levels = c("Telluraves", "Australaves", "Afroaves", "Eufalconimorphae",
                                  "pai_Accipitrimorphae", "Accipitrimorphae"))
hm$topologia <- factor(hm$topologia, levels = rev(TOPO))
hm$rot <- ifelse(hm$monofiletico, sprintf("%.2f", hm$Prap), "n.m.")
hm$caractere <- factor(hm$caractere, levels = c("C1_T50", "C4_bico_estrito", "C5_pes", "C6_sindrome"),
                       labels = c("C1 trófico (AVONET, T50)", "C4 bico (estrito)", "C5 pés", "C6 síndrome = 2"))
g <- ggplot(hm, aes(no, topologia, fill = Prap)) + geom_tile(colour = "white") +
  geom_text(aes(label = rot), size = 2.6) +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = "P(rapinante)") +
  facet_grid(raiz ~ caractere) + theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL,
       title = "Probabilidade marginal do estado rapinante (média ponderada por AICc; árvores datadas, sem fósseis)",
       subtitle = "linhas do painel: prior de raiz (fitzjohn / zero); n.m. = nó não monofilético na topologia; n = 202 famílias")
ggsave(file.path(fig_dir, "fig2_mapa_calor_P_rapinante.png"), g, width = 12, height = 5.5, dpi = 200)

# figura 2b: efeito dos fósseis (raiz zero), P no nó de Telluraves
fz <- subset(mg, ramos == "datada" & raiz == "zero" & modelo == "media_AICc" & no == "Telluraves" &
               caractere %in% c("C1_T50", "C3_elton_T50", "C4_bico_estrito", "C5_pes", "C6_sindrome"))
fq <- read.delim(file.path(res_dir, "marginais_nos_fosQ.tsv"))
fq <- subset(fq, raiz == "zero" & modelo == "media_AICc" & no == "Telluraves" & caractere %in% c("C1_T50_fosQ", "C3_elton_T50_fosQ"))
fq$Prap <- fq$P1
fq$fosseis <- paste0(fq$fosseis, "-dieta?")
fq$caractere <- sub("_fosQ$", "", fq$caractere)
fz <- rbind(fz[, c("topologia", "fosseis", "caractere", "Prap")], fq[, c("topologia", "fosseis", "caractere", "Prap")])
fz$fosseis <- factor(fz$fosseis, levels = c("nenhum", "F1", "F1-dieta?", "F2", "F2-dieta?"))
fz$topologia <- factor(fz$topologia, levels = rev(TOPO))
g <- ggplot(fz, aes(fosseis, topologia, fill = Prap)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", Prap)), size = 2.6) +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = "P(rapinante)") +
  facet_wrap(~ caractere, nrow = 1) + theme_minimal(base_size = 9) +
  labs(x = "cenário de fósseis", y = NULL,
       title = "Nó de Telluraves, raiz fixada em 0: efeito dos fósseis (F1, F2) e da dieta fóssil tratada como desconhecida",
       subtitle = "C6: P(síndrome completa, estado 2); média ponderada por AICc; árvores datadas (idades de Stiller 2024)")
ggsave(file.path(fig_dir, "fig2b_efeito_fosseis_telluraves.png"), g, width = 12, height = 4, dpi = 200)

# ---------------------------------------------------------------- figura 3: T1 com marginais (C1_T50, média AICc)
desenha_T <- function(nome, char, raiz, arquivo) {
  job <- Filter(function(r) r$arvore == nome && r$caractere == char && r$raiz == raiz,
                lapply(list.files(file.path(res_dir, "jobs"), full.names = TRUE), readRDS))[[1]]
  tr <- arv[[nome]]
  cod <- setNames(mat[[char]], mat$chave)[tr$tip.label]
  # recalcula as marginais nos nós internos com o melhor modelo (para a figura)
  st <- if (char == "C6_sindrome") c("0", "1", "2") else c("0", "1")
  x <- matrix(0, length(cod), length(st), dimnames = list(names(cod), st))
  for (i in seq_along(cod)) { s <- strsplit(cod[i], "&", fixed = TRUE)[[1]]; x[i, s] <- 1 / length(s) }
  # Q do melhor modelo reconstruída das taxas guardadas em 03_mk.R (mesma MLE, sem reotimizar)
  txt <- job$ajustes$taxas[job$ajustes$modelo == job$melhor]
  rates <- as.numeric(regmatches(txt, gregexpr("[0-9.]+(e-?[0-9]+)?", sub(" [(].*$", "", txt)))[[1]])
  mods <- list(ER = c(0, 1, 1, 0), ARD = c(0, 1, 2, 0), IRR_perda = c(0, 0, 1, 0), IRR_ganho = c(0, 1, 0, 0))
  dm <- matrix(mods[[job$melhor]], 2, 2, byrow = TRUE)
  if (grepl("Q de ", txt)) dm <- matrix(mods[[sub("^.*Q de ([A-Za-z_]+)[)].*$", "\\1", txt)]], 2, 2, byrow = TRUE)
  Q <- matrix(0, 2, 2, dimnames = list(st, st))
  for (k in seq_along(rates)) Q[dm == k] <- rates[k]
  diag(Q) <- -rowSums(Q)
  pi <- if (raiz == "fitzjohn") "fitzjohn" else setNames(c(1, 0), st)
  fit <- fitMk(tr, x, fixedQ = Q, pi = pi)
  a <- ancr(fit)$ace
  png(arquivo, width = 2200, height = 4200, res = 220)
  cols <- setNames(c("grey85", "firebrick"), st)
  plotTree(tr, fsize = 0.42, ftype = "i", lwd = 1, offset = 0.6)
  tc <- t(sapply(cod, function(cc) { v <- setNames(c(0, 0), st); s <- strsplit(cc, "&")[[1]]; v[s] <- 1 / length(s); v }))
  tiplabels(pie = tc, piecol = cols, cex = 0.28)
  foco <- c(getMRCA(tr, fams$family[fams$clade_mag7 == "Telluraves"]),
            getMRCA(tr, fams$family[fams$order %in% c("Accipitriformes")]),
            getMRCA(tr, fams$family[fams$order %in% c("Strigiformes")]),
            getMRCA(tr, fams$family[fams$order %in% c("Falconiformes", "Psittaciformes", "Passeriformes")]),
            getMRCA(tr, fams$family[fams$order %in% c("Cariamiformes", "Falconiformes", "Psittaciformes", "Passeriformes")]),
            Ntip(tr) + 1)
  internos <- as.numeric(rownames(a))
  nodelabels(node = internos, pie = a, piecol = cols, cex = 0.22)
  nodelabels(node = foco, pie = a[as.character(foco), ], piecol = cols, cex = 0.7)
  legend("topleft", legend = c("não rapinante (0)", "rapinante (1)"), fill = cols, bty = "n", cex = 0.8)
  title(sprintf("%s | %s | raiz %s | melhor modelo %s (logL %.2f)", nome, char, raiz, job$melhor,
                job$ajustes$logL[job$ajustes$modelo == job$melhor]), cex.main = 0.8)
  dev.off()
}
desenha_T("T1_E_datada", "C1_T50", "fitzjohn", file.path(fig_dir, "fig3_T1_C1_T50_raiz_fitzjohn.png"))
desenha_T("T1_E_datada", "C1_T50", "zero", file.path(fig_dir, "fig3b_T1_C1_T50_raiz_zero.png"))

# ---------------------------------------------------------------- figura 4: ganhos e perdas
gp <- subset(sm, caractere %in% c("C1_T50", "C4_bico_estrito", "C5_pes") & escopo == "coroa_Telluraves" &
               fosseis == "nenhum" & transicao %in% c("0->1", "1->0"))
gp$topologia <- factor(gp$topologia, levels = TOPO)
gp$transicao <- factor(gp$transicao, labels = c("ganhos (0->1)", "perdas (1->0)"))
g <- ggplot(gp, aes(topologia, media, colour = modelo_mapa, shape = raiz)) +
  geom_pointrange(aes(ymin = q025, ymax = q975), position = position_dodge(width = 0.6), size = 0.3) +
  facet_grid(transicao ~ caractere, scales = "free_y") + theme_bw(base_size = 9) +
  labs(y = "número de transições na coroa de Telluraves (média e IC 95%, 1000 mapas)", x = NULL,
       colour = "modelo do mapa", shape = "raiz")
ggsave(file.path(fig_dir, "fig4_ganhos_perdas_telluraves.png"), g, width = 11, height = 6, dpi = 200)

print(tp)
