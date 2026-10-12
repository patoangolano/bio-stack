# bio-stack — ambiente R do projeto (renv), Windows.
#
#   Rscript r/instalar.R
#
# Idempotente. Sem r/renv.lock: cria o projeto renv, instala os pacotes abaixo
# e grava o lock. Com r/renv.lock: RESTAURA as versoes travadas (renv::restore)
# e so instala, e regrava o lock, o que tiver sido acrescentado as listas
# depois. O lock e o que a secao de Metodos cita, do mesmo jeito que o uv.lock
# do lado Python — instalar "a mais nova" e regravar o lock em silencio
# desfaria exatamente isso.
#
# Por que renv e nao a biblioteca do usuario: o R 4.6.1 desta maquina estava
# vazio (30 pacotes base). Biblioteca global acumula versao de projeto em
# projeto e ninguem sabe depois qual phytools fez a figura. O renv guarda
# biblioteca por projeto, com cache compartilhado, e o lockfile fixa versao e
# repositorio.
#
# Conjuntos escolhidos pela pesquisa de 2026-10-11 (docs/PESQUISA-STACK-2026-10.md):
#   filogenia comparativa e ASR de carater — ape, phangorn, phytools, castor,
#     geiger, corHMM, TreeDist; figura de arvore — ggtree, treeio, tidytree
#   expressao diferencial (referencia para conferir o pydeseq2) — DESeq2,
#     edgeR, limma, sva, apeglm, ashr, IHW, qvalue

local({
  aqui <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "."), winslash = "/")
  setwd(aqui)
})

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c("ape", "phangorn", "phytools", "castor", "geiger", "corHMM", "TreeDist",
          "tidytree", "ashr", "BiocManager",
          "phylter",  # 2026-10-11: outliers por gene e por especie em filogenomica
          # 2026-10-12: morfometria geometrica e fosseis. geomorph faz Procrustes
          # de pontos, curvas e superficies, ANOVA de forma e sinal filogenetico;
          # paleobioDB consulta ocorrencia fossil pela API da Paleobiology Database.
          "geomorph", "paleobioDB")
bioc <- c("ggtree", "treeio", "DESeq2", "edgeR", "limma", "sva", "apeglm", "IHW", "qvalue",
          # 2026-10-11: figura anotada, efeito aleatorio no DE (dream), importacao salmon/kallisto
          "ggtreeExtra", "variancePartition", "tximport")

# Pacote que so existe no GitHub. Vetor NOMEADO: o nome e como o pacote entra no
# renv.lock e e por ele que a checagem de `faltando` compara; o valor e a
# especificacao que o renv resolve. Sem o nome, `setdiff` compararia
# "nclark-lab/RERconverge" contra "RERconverge" e reinstalaria a cada execucao.
#
# RERconverge (2026-10-12): correlaciona taxa evolutiva relativa de gene com
# evolucao de traco ao longo de uma filogenia. E a ponte entre morfologia medida
# e sequencia — o que liga este bloco ao projeto de atavismo em Archosauria.
github <- c(RERconverge = "nclark-lab/RERconverge")

# renv no diretorio de usuario: a biblioteca do sistema (Program Files) exige admin.
lib_usuario <- Sys.getenv("R_LIBS_USER")
dir.create(lib_usuario, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_usuario, .libPaths()))
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", lib = lib_usuario)
}

instalar_e_travar <- function(pacotes) {
  alvo <- ifelse(pacotes %in% bioc, paste0("bioc::", pacotes), pacotes)
  # Pacote de GitHub precisa da especificacao "dono/repo", nao do nome simples:
  # `renv::install("RERconverge")` iria ao CRAN e nao acharia nada.
  de_github <- pacotes %in% names(github)
  alvo[de_github] <- unname(github[pacotes[de_github]])
  renv::install(alvo, prompt = FALSE)
  renv::snapshot(type = "all", prompt = FALSE)
}

if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE, bioconductor = "3.23", restart = FALSE)
  renv::settings$bioconductor.version("3.23")
  # "all": o ambiente e do projeto inteiro, nao de scripts soltos em r/. No modo
  # implicito o renv so conta pacote chamado por library() em arquivo daqui, marca
  # os outros como "nao usados" e reclama de projeto fora de sincronia a cada inicio.
  renv::settings$snapshot.type("all")
  instalar_e_travar(c(cran, bioc, names(github)))
} else {
  renv::load(".")
  renv::restore(prompt = FALSE)
  faltando <- setdiff(c(cran, bioc, names(github)), names(renv::lockfile_read("renv.lock")$Packages))
  if (length(faltando) > 0) {
    cat("Acrescentados as listas depois do lock:", faltando, "\n")
    instalar_e_travar(faltando)
  }
}

cat("\nVersoes:\n")
for (p in c(cran, bioc, names(github))) cat(sprintf("  %-12s %s\n", p, as.character(packageVersion(p))))
