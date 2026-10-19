# 03_mk.R
# Modelos Mk (ER, SYM, ARD, irreversível só-perda, irreversível só-ganho) com AICc,
# probabilidades marginais (empirical Bayes; phytools::ancr) nos nós focais,
# média ponderada por AICc, e mapeamento estocástico (make.simmap, Q fixa na MLE).
#
# Duas condições de raiz (Neognathae = Galloanserae vs Neoaves):
#   fitzjohn  prior de raiz de FitzJohn et al. (2009)
#   zero      raiz fixada no estado 0 (não rapinante), justificada pelo grupo externo
#             (Palaeognathae e Galloanserae sem nenhuma família rapinante); com ela o
#             modelo só-perda fica inviável sempre que houver ponta no estado 1.
# Para caractere binário SYM = ER (a matriz simétrica 2x2 tem uma única taxa); o SYM
# só é ajustado para o caractere de 3 estados (C6).
#
# Entradas: arvores/arvores.rds (02_arvores.R), caracteres/matriz_mk.tsv (01_caracteres.py)
# Saídas:   resultados/ajustes_modelos.tsv, resultados/marginais_nos.tsv,
#           resultados/simmap_transicoes.tsv, resultados/simmap_estado_nos.tsv,
#           resultados/jobs/*.rds, logs/03_mk_sessionInfo.txt
#
# Ajuste: castor::fit_mk (16 partidas) para a MLE; reavaliação no phytools com a Q fixa;
# se um modelo ficou abaixo de um modelo aninhado nele (ótimo local), recebe a Q do aninhado.
# Aleatoriedade: partidas do fit_mk (set.seed(SEMENTE_BASE + 10*id) antes dos ajustes) e
# make.simmap (set.seed(SEMENTE_BASE + 10*id + j) imediatamente antes; j = 1 melhor modelo,
# 2 = ER), sempre dentro do processo-filho;
# RNGkind padrão (Mersenne-Twister / Inversion / Rejection). As sementes vão para a tabela.
# Núcleos: no máximo 4 (cluster PSOCK).
#
# Execução, a partir de r/:
#   Rscript ../runs/rapinantes/ancestral/scripts/03_mk.R [n_mapas] [regex_arvore] [regex_caractere]

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(parallel)
  library(castor)
})

SEMENTE_BASE <- 20261017L
NUCLEOS <- 4L
args <- commandArgs(trailingOnly = TRUE)
N_MAPAS <- if (length(args) >= 1) as.integer(args[1]) else 1000L
CAR_SIMMAP <- c("C1_T50", "C4_bico_estrito", "C5_pes", "C6_sindrome")

base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
res_dir <- file.path(base, "resultados")
job_dir <- file.path(res_dir, "jobs")
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base, "logs"), showWarnings = FALSE)

arv <- readRDS(file.path(base, "arvores", "arvores.rds"))
mat <- read.delim(file.path(base, "caracteres", "matriz_mk.tsv"), colClasses = "character")
# pontas da árvore: nome da família (viventes) ou FOSSIL_* (fósseis)
mat$chave <- ifelse(mat$tipo == "fossil", mat$ponta, mat$family)
stopifnot(!anyDuplicated(mat$chave))
fams <- read.delim(file.path(base, "caracteres", "familias_pontas.tsv"))

caracteres <- c("C1_T50", "C1_T25", "C1_T75", "C1_bin50", "C2_aquat_T50", "C3_elton_T50",
                "C4_bico_estrito", "C4b_bico_amplo", "C5_pes", "C6_sindrome")

modelos_bin <- list(
  ER = matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE),
  ARD = matrix(c(0, 1, 2, 0), 2, 2, byrow = TRUE),
  IRR_perda = matrix(c(0, 0, 1, 0), 2, 2, byrow = TRUE),   # só 1 -> 0
  IRR_ganho = matrix(c(0, 1, 0, 0), 2, 2, byrow = TRUE)    # só 0 -> 1
)
modelos_tri <- list(
  ER = matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3, byrow = TRUE),
  SYM = matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3, byrow = TRUE),
  ARD = matrix(c(0, 1, 2, 3, 0, 4, 5, 6, 0), 3, 3, byrow = TRUE),
  IRR_perda = matrix(c(0, 0, 0, 1, 0, 0, 2, 3, 0), 3, 3, byrow = TRUE),  # só para estado menor
  IRR_ganho = matrix(c(0, 1, 2, 0, 0, 3, 0, 0, 0), 3, 3, byrow = TRUE)   # só para estado maior
)

grupos <- list(
  Neognathae_raiz = fams$family,
  Neoaves = fams$family[fams$clade_mag7 != "Galloanseres"],
  Telluraves = fams$family[fams$clade_mag7 == "Telluraves"],
  Australaves = fams$family[fams$order %in% c("Cariamiformes", "Falconiformes", "Psittaciformes", "Passeriformes")],
  Afroaves = fams$family[fams$order %in% c("Accipitriformes", "Strigiformes", "Coliiformes", "Leptosomatiformes",
                                           "Trogoniformes", "Bucerotiformes", "Coraciiformes", "Piciformes")],
  Eufalconimorphae = fams$family[fams$order %in% c("Falconiformes", "Psittaciformes", "Passeriformes")],
  Accipitrimorphae = fams$family[fams$order == "Accipitriformes"],
  Strigiformes = fams$family[fams$order == "Strigiformes"],
  Coraciimorphae = fams$family[fams$order %in% c("Coliiformes", "Leptosomatiformes", "Trogoniformes",
                                                 "Bucerotiformes", "Coraciiformes", "Piciformes")],
  Accipitrimorphae_Strigiformes = fams$family[fams$order %in% c("Accipitriformes", "Strigiformes")],
  Strigiformes_Coraciimorphae = fams$family[fams$order %in% c("Strigiformes", "Coliiformes", "Leptosomatiformes",
                                                              "Trogoniformes", "Bucerotiformes", "Coraciiformes",
                                                              "Piciformes")],
  Telluraves_sem_Accipitrimorphae = fams$family[fams$clade_mag7 == "Telluraves" & fams$order != "Accipitriformes"]
)
NOS_REPORTADOS <- c("Neognathae_raiz", "Neoaves", "Telluraves", "Australaves", "Afroaves", "Eufalconimorphae",
                    "Accipitrimorphae", "Strigiformes", "Coraciimorphae")

nos_focais <- function(tr) {
  fos <- grep("^FOSSIL_", tr$tip.label, value = TRUE)
  ext <- if (length(fos)) drop.tip(tr, fos) else tr
  viv <- function(n) setdiff(extract.clade(tr, n)$tip.label, fos)
  ident <- function(no) {
    tips <- viv(no)
    hit <- names(grupos)[sapply(grupos, function(g) setequal(g, tips))]
    if (length(hit)) hit[1] else paste0("outro(", length(tips), " familias)")
  }
  out <- list()
  for (g in NOS_REPORTADOS) {
    no <- if (g == "Neognathae_raiz") Ntip(tr) + 1 else getMRCA(tr, grupos[[g]])
    out[[g]] <- list(no = no, mono = (g == "Neognathae_raiz") || is.monophyletic(ext, grupos[[g]]), ident = g)
  }
  # nó-pai (o primeiro ancestral com conjunto de viventes diferente; pula nós de enxerto fóssil)
  for (g in c("Accipitrimorphae", "Strigiformes")) {
    filho <- out[[g]]$no
    base_set <- viv(filho)
    no <- filho
    repeat {
      no <- tr$edge[tr$edge[, 2] == no, 1]
      if (!setequal(viv(no), base_set)) break
    }
    out[[paste0("pai_", g)]] <- list(no = no, mono = TRUE, ident = ident(no))
  }
  out
}

prior_tips <- function(tr, char) {
  cod <- setNames(mat[[char]], mat$chave)[tr$tip.label]
  stopifnot(!anyNA(cod))
  estados <- if (char == "C6_sindrome") c("0", "1", "2") else c("0", "1")
  m <- matrix(0, length(cod), length(estados), dimnames = list(names(cod), estados))
  for (i in seq_along(cod)) {
    s <- strsplit(cod[i], "&", fixed = TRUE)[[1]]
    m[i, s] <- 1 / length(s)
  }
  m
}

ajusta <- function(tr, x, dm, pi) {
  # MLE pelo castor::fit_mk (16 partidas; muito mais rápido que o otimizador do phytools e,
  # conferido em T1, chega ao mesmo ótimo), depois o modelo é reavaliado no phytools com
  # a Q fixa (fitMk(fixedQ=)), de onde saem logL, ancr e make.simmap.
  k <- max(dm)
  st <- colnames(x)
  rp <- if (identical(pi, "fitzjohn")) "likelihoods" else as.numeric(pi)
  idx <- unname(dm); storage.mode(idx) <- "integer"
  cf <- tryCatch(fit_mk(tr, Nstates = ncol(x), tip_priors = unname(x), rate_model = idx, root_prior = rp,
                        Ntrials = 16, Nthreads = 1), error = function(e) NULL)
  if (is.null(cf) || !isTRUE(cf$success) && !is.finite(cf$loglikelihood)) return(NULL)
  if (!is.finite(cf$loglikelihood) || cf$loglikelihood < -1e10) return(NULL)
  Q <- cf$transition_matrix; dimnames(Q) <- list(st, st)
  fp <- tryCatch(fitMk(tr, x, fixedQ = Q, pi = pi), error = function(e) NULL)
  if (is.null(fp)) return(NULL)
  ll <- as.numeric(logLik(fp))
  if (!is.finite(ll) || ll < -1e10) return(NULL)
  taxas <- sapply(seq_len(k), function(i) Q[which(dm == i)[1]])
  list(fit = fp, Q = Q, k = k, logL = ll, logL_castor = cf$loglikelihood, taxas = taxas)
}

bind_fill <- function(lst) {
  lst <- Filter(function(d) !is.null(d) && nrow(d) > 0, lst)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols]
  }))
}

roda_tarefa <- function(id, nome_arv, char, raiz, faz_simmap) {
  tr <- arv[[nome_arv]]
  x <- prior_tips(tr, char)
  st <- colnames(x)
  mods <- if (char == "C6_sindrome") modelos_tri else modelos_bin
  pi <- if (raiz == "fitzjohn") "fitzjohn" else setNames(c(1, rep(0, length(st) - 1)), st)
  n <- Ntip(tr)
  ajustes <- list(); tab <- data.frame()
  set.seed(SEMENTE_BASE + 10L * id)  # partidas aleatórias do castor::fit_mk
  for (m in names(mods)) {
    dm <- mods[[m]]; dimnames(dm) <- list(st, st)
    k <- max(dm)
    aj <- ajusta(tr, x, dm, pi)
    if (is.null(aj)) {
      tab <- rbind(tab, data.frame(modelo = m, k = k, logL = -Inf, logL_castor = NA, AICc = Inf, taxas = "inviavel"))
      next
    }
    aicc <- -2 * aj$logL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
    ajustes[[m]] <- aj$fit
    tab <- rbind(tab, data.frame(modelo = m, k = k, logL = aj$logL, logL_castor = aj$logL_castor, AICc = aicc,
                                 taxas = paste(signif(aj$taxas, 4), collapse = ";")))
  }
  # reparo de aninhamento: a MLE de um modelo não pode ficar abaixo da de um modelo aninhado nele
  aninhados <- if (length(st) == 2) list(ARD = c("ER", "IRR_perda", "IRR_ganho"))
               else list(SYM = "ER", ARD = c("ER", "SYM", "IRR_perda", "IRR_ganho"))
  tab$reparo <- ""
  for (M in names(aninhados)) {
    for (N in aninhados[[M]]) {
      iM <- which(tab$modelo == M); iN <- which(tab$modelo == N)
      if (is.null(ajustes[[N]])) next
      if (is.null(ajustes[[M]]) || tab$logL[iM] < tab$logL[iN] - 1e-6) {
        ajustes[[M]] <- ajustes[[N]]
        k <- tab$k[iM]
        tab$logL[iM] <- tab$logL[iN]
        tab$AICc[iM] <- -2 * tab$logL[iN] + 2 * k + 2 * k * (k + 1) / (n - k - 1)
        tab$taxas[iM] <- paste0(tab$taxas[iN], " (Q de ", N, ")")
        tab$reparo[iM] <- paste0("recebeu a MLE de ", N)
      }
    }
  }
  tab$dAICc <- tab$AICc - min(tab$AICc)
  tab$peso <- ifelse(is.finite(tab$AICc), exp(-0.5 * tab$dAICc), 0)
  tab$peso <- tab$peso / sum(tab$peso)
  melhor <- tab$modelo[which.min(tab$AICc)]

  nos <- nos_focais(tr)
  anc_list <- list()
  for (m in names(ajustes)) {
    a <- tryCatch(ancr(ajustes[[m]])$ace, error = function(e) NULL)
    if (!is.null(a)) anc_list[[m]] <- a
  }
  ok_m <- intersect(tab$modelo[tab$peso > 0], names(anc_list))
  w <- tab$peso[match(ok_m, tab$modelo)]; w <- w / sum(w)
  anc_list[["media_AICc"]] <- Reduce(`+`, Map(function(a, ww) a * ww, anc_list[ok_m], w))
  marg <- list()
  for (m in names(anc_list)) {
    a <- anc_list[[m]]
    for (g in names(nos)) {
      lin <- a[as.character(nos[[g]]$no), ]
      d <- data.frame(modelo = m, melhor = (m == melhor), no = g, identidade = nos[[g]]$ident,
                      monofiletico = nos[[g]]$mono)
      for (s in st) d[[paste0("P", s)]] <- unname(lin[s])
      marg[[length(marg) + 1]] <- d
    }
  }
  marg <- bind_fill(marg)

  sm <- list(); sn <- list(); sementes <- c()
  if (faz_simmap) {
    # melhor modelo nas duas raízes; ER como referência só com a raiz de FitzJohn
    alvo <- if (raiz == "fitzjohn") unique(c(melhor, "ER")) else melhor
    tel_no <- getMRCA(tr, grupos$Telluraves)
    rotulos <- as.vector(outer(st, st, paste, sep = "->"))
    for (j in seq_along(alvo)) {
      m <- alvo[j]
      if (is.null(ajustes[[m]])) next
      Qm <- unclass(as.Qmatrix(ajustes[[m]]))[seq_along(st), seq_along(st)]
      dimnames(Qm) <- list(st, st)
      semente <- SEMENTE_BASE + 10L * id + j
      set.seed(semente)
      mapas <- tryCatch(make.simmap(tr, x, Q = Qm, pi = pi, nsim = N_MAPAS, message = FALSE),
                        error = function(e) NULL)
      if (is.null(mapas)) next
      sementes <- c(sementes, setNames(semente, m))
      tot <- countSimmap(mapas, states = st, message = FALSE)
      tot <- if (is.list(tot) && !is.null(tot$Tr)) tot$Tr else tot
      conta_tel <- t(sapply(mapas, function(mp) {
        as.numeric(countSimmap(extract.clade.simmap(mp, tel_no), states = st, message = FALSE)$Tr)
      }))
      colnames(conta_tel) <- rotulos
      rs <- function(v) data.frame(media = mean(v), mediana = median(v),
                                   q025 = unname(quantile(v, 0.025)), q975 = unname(quantile(v, 0.975)))
      for (cn in colnames(tot)) {
        if (cn == "N") next
        p <- strsplit(cn, ",", fixed = TRUE)[[1]]
        if (length(p) == 2 && p[1] == p[2]) next
        sm[[length(sm) + 1]] <- cbind(data.frame(modelo_mapa = m, semente = semente, escopo = "arvore_toda",
                                                 transicao = sub(",", "->", cn)), rs(tot[, cn]))
      }
      for (cn in rotulos) {
        p <- strsplit(cn, "->", fixed = TRUE)[[1]]
        if (p[1] == p[2]) next
        sm[[length(sm) + 1]] <- cbind(data.frame(modelo_mapa = m, semente = semente, escopo = "coroa_Telluraves",
                                                 transicao = cn), rs(conta_tel[, cn]))
      }
      # estado amostrado nos nós focais (frequência entre os mapas)
      for (g in names(nos)) {
        no <- nos[[g]]$no
        e <- sapply(mapas, function(mp) {
          k <- which(mp$edge[, 1] == no)[1]
          names(mp$maps[[k]])[1]
        })
        fr <- table(factor(e, levels = st)) / length(e)
        d <- data.frame(modelo_mapa = m, semente = semente, no = g, identidade = nos[[g]]$ident)
        for (s in st) d[[paste0("F", s)]] <- as.numeric(fr[s])
        sn[[length(sn) + 1]] <- d
      }
    }
  }
  res <- list(id = id, arvore = nome_arv, caractere = char, raiz = raiz, ajustes = tab, melhor = melhor,
              marginais = marg, simmap = bind_fill(sm), simmap_nos = bind_fill(sn),
              sementes = sementes, n_mapas = N_MAPAS)
  saveRDS(res, file.path(job_dir, sprintf("job_%04d.rds", id)))
  res
}

tarefas <- expand.grid(arvore = names(arv), caractere = caracteres, raiz = c("fitzjohn", "zero"),
                       stringsAsFactors = FALSE)
tarefas$id <- seq_len(nrow(tarefas))
tarefas$simmap <- tarefas$caractere %in% CAR_SIMMAP & !grepl("grafen", tarefas$arvore)
# filtros opcionais para teste: argumento 2 = regex de árvore, 3 = regex de caractere
if (length(args) >= 2) tarefas <- tarefas[grepl(args[2], tarefas$arvore), ]
if (length(args) >= 3) tarefas <- tarefas[grepl(args[3], tarefas$caractere), ]
# tarefas com mapeamento primeiro (mais longas) para equilibrar a carga
tarefas <- tarefas[order(!tarefas$simmap, tarefas$caractere != "C6_sindrome"), ]
# retomada: tarefa cujo resultado já está em jobs/ não é refeita (a semente é função do id,
# então refazer daria o mesmo resultado)
feitos <- as.integer(sub("^job_0*([0-9]+)[.]rds$", "\\1", list.files(job_dir, pattern = "^job_[0-9]+[.]rds$")))
todas <- tarefas
tarefas <- tarefas[!tarefas$id %in% feitos, ]
cat(sprintf("%d tarefas no total, %d já feitas, %d a rodar (%d com mapeamento), %d mapas, %d núcleos\n",
            nrow(todas), nrow(todas) - nrow(tarefas), nrow(tarefas), sum(tarefas$simmap), N_MAPAS, NUCLEOS))

cl <- makeCluster(NUCLEOS)
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(ape); library(phytools); library(castor) })))
clusterExport(cl, c("arv", "mat", "fams", "grupos", "NOS_REPORTADOS", "modelos_bin", "modelos_tri",
                    "nos_focais", "prior_tips", "ajusta", "bind_fill", "roda_tarefa",
                    "SEMENTE_BASE", "N_MAPAS", "job_dir"))
t0 <- Sys.time()
# chunk.size = 1: cada tarefa é uma unidade de escalonamento; sem isso o parLapplyLB
# entrega blocos contíguos e um núcleo fica com todos os mapeamentos pesados
res <- parLapplyLB(cl, seq_len(nrow(tarefas)), function(i, tf) {
  tryCatch(roda_tarefa(tf$id[i], tf$arvore[i], tf$caractere[i], tf$raiz[i], tf$simmap[i]),
           error = function(e) list(id = tf$id[i], erro = conditionMessage(e)))
}, tf = tarefas, chunk.size = 1)
stopCluster(cl)
cat("tempo:", format(Sys.time() - t0), "\n")

erros <- Filter(function(r) !is.null(r$erro), res)
for (e in erros) cat("ERRO tarefa", e$id, ":", e$erro, "\n")
# agrega a partir dos arquivos de jobs/ (inclui as tarefas de rodadas anteriores)
arqs <- file.path(job_dir, sprintf("job_%04d.rds", todas$id))
ok <- lapply(arqs[file.exists(arqs)], readRDS)
tarefas <- todas

meta <- function(r) {
  p <- strsplit(r$arvore, "_")[[1]]
  data.frame(arvore = r$arvore, topologia = p[1], ramos = if (p[2] == "E") p[3] else "datada",
             fosseis = if (p[2] == "E") "nenhum" else p[2], caractere = r$caractere, raiz = r$raiz)
}
junta <- function(campo) bind_fill(lapply(ok, function(r) {
  d <- r[[campo]]
  if (is.null(d) || !nrow(d)) return(NULL)
  cbind(meta(r)[rep(1, nrow(d)), ], d, row.names = NULL)
}))
write.table(junta("ajustes"), file.path(res_dir, "ajustes_modelos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(junta("marginais"), file.path(res_dir, "marginais_nos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sm <- junta("simmap"); if (!is.null(sm)) sm$n_mapas <- N_MAPAS
write.table(sm, file.path(res_dir, "simmap_transicoes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(junta("simmap_nos"), file.path(res_dir, "simmap_estado_nos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "",
             sprintf("SEMENTE_BASE=%d; semente = SEMENTE_BASE + 10*id + j (j=1 melhor modelo, j=2 ER); N_MAPAS=%d; NUCLEOS=%d",
                     SEMENTE_BASE, N_MAPAS, NUCLEOS),
             sprintf("tarefas=%d ok=%d erros=%d", nrow(tarefas), length(ok), length(erros))),
           file.path(base, "logs", "03_mk_sessionInfo.txt"))
cat("ok:", length(ok), " erros:", length(erros), "\n")
