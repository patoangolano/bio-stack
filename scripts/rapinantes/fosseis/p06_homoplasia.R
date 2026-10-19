# Passo 7 (2d): homoplasia por carater na arvore morfologica x na arvore
# molecular-restrita, e teste de permutacao para os caracteres "rapinantes".
#
#   cd r; Rscript ../scripts/rapinantes/fosseis/p06_homoplasia.R
#
# Arvores (IQ-TREE, passo 3):  livre.treefile (morfologia sozinha) e
# esqueleto.treefile (atuais presos a Stiller et al. 2024, fosseis livres).
# Duas versoes de cada comparacao: "completa" (com fosseis) e "atuais" (fosseis
# podados, dados restritos aos atuais). IC e IR de Fitch (phangorn), carater
# nao ordenado, '?' = qualquer estado. So entram nos testes os caracteres
# informativos para parcimonia no conjunto de taxa usado.
#
# Classe "rapinante" (pela descricao publicada do carater): bico em gancho,
# cera, sulco tomial, bico alto e comprimido; unguais, halux, dedos do pe e
# seus tendoes; tarsometatarso e hipotarso. Descricao so existe para VZ26
# (Mayr 2026, File S1), MK22 (Mayr & Kitchener 2022, SD) e C25 (Chen et al.
# 2025, S1; so cintura escapular e membro anterior).
#
# Teste: diferenca de medias (rapinante - resto) de IC na arvore molecular e de
# deltaIC = IC_mol - IC_morf; p unilateral (rapinante menor) por 10 000
# permutacoes dos rotulos, semente fixa; no conjunto agregado a permutacao e
# estratificada por matriz.

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})
set.seed(20261017)
BASE <- normalizePath(file.path("..", "runs", "rapinantes", "fosseis"), winslash = "/")
OUT <- file.path(BASE, "homoplasia")
dir.create(OUT, showWarnings = FALSE)
IDS <- c("VZ26", "MK22", "C25", "D21", "K19", "M11", "M05", "M16", "MC03", "LZ07")
NPERM <- 10000

# Ingles (Mayr; Chen et al.), latim anatomico (Livezey & Zusi) e portugues
# (Degrange et al.). "bico" = gancho do bico, cera, tomio/crista tomial,
# bico alto e comprimido; "pe_garra" = unguais, halux, dedos do pe e seus
# tendoes; "tarsometatarso" = tarsometatarso, hipotarso, metatarsais.
classe_rapinante <- function(d) {
  d <- tolower(d)
  asa <- grepl(paste0("alul|digiti major|digiti minor|majoris|minoris|major digit|minor digit|manus|manual|",
                      "carpometacarp|wing|\\bradi|ulna|humer|úmer|coracoid|\\basa\\b|álula|\\balae\\b|",
                      "antebrach|remig|stern|furcul|scapul"), d)
  bico <- (grepl("hook|uncinat|adunc|gancho", d) &
             grepl("beak|bill|rhamph|maxill|bico|ranfoteca|premaxil|pré-maxil|culmen", d)) |
    grepl("\\bcere\\b|\\bcera\\b|tomi(al|um|a)\\b|tomial", d) |
    grepl("upper beak dorsoventrally very deep", d)
  pe <- grepl(paste0("ungual|ungu(es|is)|claw|garra|hallu|hálux|halux|\\btoes?\\b|\\bdedos?\\b|",
                     "digit[a-z]* pedis|digitorum pedis|pedal|falang|phalan"), d) & !asa
  tmt <- grepl("tarsometatars|hypotars|hipotars|metatars", d) & !asa
  ifelse(bico, "bico", ifelse(pe, "pe_garra", ifelse(tmt, "tarsometatarso", "outro")))
}

ler_arv <- function(f) {
  t <- read.tree(f)
  t$node.label <- NULL
  t
}

informativo <- function(dat) {
  # informativo para parcimonia: >= 2 estados com >= 2 taxa cada
  m <- as.character(dat)
  apply(m, 2, function(x) {
    x <- x[x != "?"]
    tb <- table(x)
    sum(tb >= 2) >= 2
  })
}

linhas <- list()
ens <- list()
for (id in IDS) {
  pa <- file.path(BASE, "arvores", id)
  if (!file.exists(file.path(pa, "livre.treefile"))) next
  col <- read.delim(file.path(BASE, "matrizes", id, "colunas.tsv"), quote = "", stringsAsFactors = FALSE)
  tx <- read.delim(file.path(BASE, "matrizes", id, "taxa.tsv"), stringsAsFactors = FALSE)
  dat <- read.phyDat(file.path(BASE, "matrizes", id, paste0(id, ".phy")), format = "phylip",
                     type = "USER", levels = as.character(0:9), ambiguity = "?")
  atuais <- tx$terminal[tx$tipo == "atual"]
  t_morf <- ler_arv(file.path(pa, "livre.treefile"))
  t_mol <- if (file.exists(file.path(pa, "esqueleto.treefile"))) ler_arv(file.path(pa, "esqueleto.treefile")) else NULL
  if (is.null(t_mol)) next
  # arvore alternativa: esqueleto com Falconidae junto de Accipitridae (busca ML restrita, passo 5)
  f_alt <- file.path(BASE, "au", id, "atuais_Falconidae", "H_Falconidae_com_Accipitridae.treefile")
  t_alt <- if (file.exists(f_alt)) ler_arv(f_alt) else NULL
  versoes <- list(completa = list(dat = dat, morf = t_morf, mol = t_mol, alt = t_alt))
  if (length(atuais) >= 6) {
    da <- subset(dat, atuais)
    versoes$atuais <- list(dat = da, morf = keep.tip(t_morf, atuais), mol = keep.tip(t_mol, atuais),
                           alt = if (is.null(t_alt)) NULL else keep.tip(t_alt, atuais))
  }
  for (v in names(versoes)) {
    x <- versoes[[v]]
    # sitewise=TRUE devolve um valor por padrao de sitio; expandir para colunas
    idx <- attr(x$dat, "index")
    ci_m <- CI(x$morf, x$dat, sitewise = TRUE)[idx]
    ci_l <- CI(x$mol, x$dat, sitewise = TRUE)[idx]
    ri_m <- RI(x$morf, x$dat, sitewise = TRUE)[idx]
    ri_l <- RI(x$mol, x$dat, sitewise = TRUE)[idx]
    st_m <- parsimony(x$morf, x$dat, site = "site")[idx]
    st_l <- parsimony(x$mol, x$dat, site = "site")[idx]
    st_a <- if (is.null(x$alt)) rep(NA_real_, length(idx)) else parsimony(x$alt, x$dat, site = "site")[idx]
    inf <- informativo(x$dat)
    ens[[length(ens) + 1]] <- data.frame(
      matriz = id, versao = v, ntax = length(x$morf$tip.label), n_informativos = sum(inf),
      passos_morf = parsimony(x$morf, x$dat), passos_mol = parsimony(x$mol, x$dat),
      CI_morf = CI(x$morf, x$dat), CI_mol = CI(x$mol, x$dat),
      RI_morf = RI(x$morf, x$dat), RI_mol = RI(x$mol, x$dat))
    desc <- col$descricao
    desc[is.na(desc)] <- ""
    cls <- if (any(nchar(desc) > 0)) classe_rapinante(desc) else rep(NA_character_, nrow(col))
    linhas[[length(linhas) + 1]] <- data.frame(
      matriz = id, versao = v, coluna = col$coluna, carater_original = col$carater_original,
      classe = cls, informativo = inf, CI_morf = ci_m, CI_mol = ci_l, RI_morf = ri_m, RI_mol = ri_l,
      passos_morf = st_m, passos_mol = st_l, passos_FalcAcc = st_a, descricao = substr(desc, 1, 160))
  }
}
car <- do.call(rbind, linhas)
ensemble <- do.call(rbind, ens)
write.table(car, file.path(OUT, "caracteres_CI_RI.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ensemble, file.path(OUT, "ensemble_CI_RI.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------- teste de permutacao ----------------
teste <- function(df, var, estrato = NULL) {
  r <- df$classe != "outro"
  obs <- mean(df[[var]][r]) - mean(df[[var]][!r])
  perm <- replicate(NPERM, {
    rr <- if (is.null(estrato)) sample(r) else ave(r, df[[estrato]], FUN = sample)
    mean(df[[var]][rr]) - mean(df[[var]][!rr])
  })
  c(n_rap = sum(r), n_resto = sum(!r), media_rap = mean(df[[var]][r]), media_resto = mean(df[[var]][!r]),
    dif = obs, p_unilateral_menor = (1 + sum(perm <= obs)) / (NPERM + 1))
}
res <- list()
base <- subset(car, !is.na(classe) & informativo & !is.na(CI_mol) & !is.na(CI_morf))
base$dCI <- base$CI_mol - base$CI_morf
base$dRI <- base$RI_mol - base$RI_morf
for (v in unique(base$versao)) {
  for (id in unique(base$matriz)) {
    d <- subset(base, matriz == id & versao == v)
    if (sum(d$classe != "outro") < 3 || sum(d$classe == "outro") < 3) next
    for (var in c("CI_mol", "CI_morf", "dCI", "RI_mol", "dRI")) {
      dd <- d[!is.na(d[[var]]), ]
      res[[length(res) + 1]] <- data.frame(matriz = id, versao = v, variavel = var, t(teste(dd, var)))
    }
  }
  ok <- names(which(tapply(base$classe[base$versao == v] != "outro", base$matriz[base$versao == v], sum) >= 3))
  d <- subset(base, versao == v & matriz %in% ok)
  if (length(unique(d$matriz)) >= 2) {
    for (var in c("CI_mol", "CI_morf", "dCI", "RI_mol", "dRI")) {
      dd <- d[!is.na(d[[var]]), ]
      res[[length(res) + 1]] <- data.frame(matriz = paste(sort(unique(dd$matriz)), collapse = "+"), versao = v,
                                           variavel = var, t(teste(dd, var, estrato = "matriz")))
    }
  }
}
perm <- do.call(rbind, res)

# ------- uniao falcao-gaviao: que caracteres pagam menos passos com Falconidae junto de Accipitridae -------
fa <- subset(car, informativo & !is.na(passos_FalcAcc))
fa$ganho <- fa$passos_mol - fa$passos_FalcAcc   # > 0: carater favorece Falconidae + Accipitridae
uniao <- do.call(rbind, lapply(split(fa, list(fa$matriz, fa$versao), drop = TRUE), function(d) {
  r <- !is.na(d$classe) & d$classe != "outro"
  fav <- d$ganho > 0
  ft <- if (any(!is.na(d$classe)) && sum(r) > 0 && sum(fav) > 0) fisher.test(table(factor(r, c(FALSE, TRUE)),
                                                                              factor(fav, c(FALSE, TRUE)))) else NULL
  data.frame(matriz = d$matriz[1], versao = d$versao[1], n_informativos = nrow(d),
             n_favorecem_FalcAcc = sum(fav), n_contra = sum(d$ganho < 0), saldo_passos = sum(d$ganho),
             n_rapinantes = if (any(!is.na(d$classe))) sum(r) else NA,
             rapinantes_que_favorecem = if (any(!is.na(d$classe))) sum(r & fav) else NA,
             OR_Fisher = if (is.null(ft)) NA else unname(ft$estimate), p_Fisher = if (is.null(ft)) NA else ft$p.value)
}))
write.table(uniao, file.path(OUT, "uniao_falcao_gaviao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fa[fa$ganho > 0, c("matriz", "versao", "carater_original", "classe", "ganho", "CI_mol", "descricao")],
            file.path(OUT, "caracteres_que_unem_falcao_gaviao.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(uniao, digits = 3)
write.table(perm, file.path(OUT, "permutacao_rapinantes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(ensemble, digits = 3)
print(perm, digits = 3)
cat("phangorn", as.character(packageVersion("phangorn")), "ape", as.character(packageVersion("ape")), "\n",
    file = file.path(OUT, "versoes_R.txt"))
