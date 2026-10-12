# Ambiente R do projeto (renv)

R 4.6.1 (Windows) + Bioconductor 3.23, biblioteca do projeto em `r/renv/`
(fora do git), versões travadas em `r/renv.lock` — o equivalente do
`uv.lock` para o lado R. Montado em 2026-10-11: 186 pacotes instalados em
88 s, todos binários do CRAN/Bioconductor exceto três compilados (IHW,
KEGGREST, ggfun); o lock travou 202 — esses 186, os 15 recomendados que
vêm com o próprio R (Matrix, MASS, mgcv, nlme…) e o renv, porque o snapshot
`type = "all"` registra tudo o que a biblioteca usa. No mesmo dia entraram
phylter, ggtreeExtra, variancePartition e tximport (35 pacotes com as
dependências, nenhuma versão travada mudou): o `renv.lock` trava 237.

```powershell
Rscript r/instalar.R          # sem lock: monta e trava; com lock: restaura as versões travadas
cd r; Rscript seu_script.R    # o .Rprofile ativa o renv ao iniciar o R aqui
```

| Para quê | Pacotes (versão travada) |
|---|---|
| filogenia comparativa, ASR de caráter | ape 5.8.1, phangorn 2.12.1, phytools 2.5.2, castor 1.8.7, geiger 2.0.12, corHMM 2.8, TreeDist 2.15.0; phylter 0.9.12 (gene ou espécie fora do padrão numa filogenômica) |
| figura de árvore | ggtree 4.2.0, treeio 1.36.1, tidytree 0.4.8, ggtreeExtra 1.22.1 (camadas anotadas ao lado da árvore) |
| expressão diferencial (referência para conferir o pydeseq2) | DESeq2 1.52.0, edgeR 4.10.5, limma 3.68.5, sva 3.60.0, apeglm 1.34.0, ashr 2.2.63, IHW 1.40.0, qvalue 2.44.0; variancePartition 1.42.0 (`dream`: efeito aleatório de animal, ninhada, doador); tximport 1.40.0 (salmon/kallisto) |

Teste rápido, a partir de `r/`: ler a árvore de espécies do teste do
pipeline de filogenia e compará-la com a verdadeira dá
`TreeDist::RobinsonFoulds(...) = 0` (conferido na montagem).

Duas armadilhas:

- A árvore do wASTRAL traz `:nan` em boa parte dos ramos **internos** (os
  terminais vêm sem comprimento nenhum), e o comprimento que existe não é
  unidade coalescente: o modo híbrido mistura coalescente com substituição
  (CU + 2×SU). O `ape::read.tree` **aceita em silêncio** e põe `NaN` em
  `edge.length` — e tirar o `:nan` do texto antes não resolve: o ape grava
  `NaN` em todo ramo sem comprimento (20 dos 22 na árvore de teste, com ou
  sem o `gsub`; conferido em 2026-10-11). Para topologia (RF,
  figura), descarte o comprimento: `tr <- read.tree(f); tr$edge.length <- NULL`.
  Para análise que use comprimento, reestime na topologia fixa com o IQ-TREE
  (receita em `docs/FILOGENIA.md`, "Como ler o resultado").
- `TreeTools` (dependência do TreeDist) e `tidytree` (importado por treeio e
  ggtree, logo também pelo ggtreeExtra) declaram, cada um, a classe S4
  `phylo`. Com os dois na mesma sessão o R emite, por `message()`, `Found
  more than one class "phylo" in cache; using the first, from namespace
  'TreeTools'` — é inofensivo (RF, `as.treedata`, `drop.tip` e `ggtree` dão
  o mesmo resultado nas duas ordens de carga). O phytools não tem papel
  nisso. Para silenciar, `suppressMessages()`.

Para acrescentar pacote: inclua o nome nas listas `cran`/`bioc` de
`instalar.R` e rode o script — ele restaura o que está travado e instala e
trava só o que é novo. (`renv::install` + `renv::snapshot(type = "all")`
dentro de `r/` faz o mesmo à mão.)
