# Pipeline de filogenia — IQ-TREE 3

`pipelines/filogenia` leva locos não alinhados (um FASTA por loco) até a
árvore de máxima verossimilhança por loco, a supermatriz particionada, os
fatores de concordância e a árvore de espécies — e devolve o parágrafo de
Métodos com as versões que de fato rodaram. Montado e validado em
2026-10-11.

```
FASTA por loco
  └─ PREPARA_LOCO ........ saneia nomes, valida alfabeto, descarta loco com < 4 táxons
      └─ MAFFT_ALIGN ..... L-INS-i até 200 seqs (--auto acima)     [ou FAMSA_ALIGN]
          └─ CLIPKIT ..... smart-gap                                [ou TRIMAL -automated1]
              ├─ IQTREE_GENES ........ ModelFinder (MFP) + SH-aLRT 1000 + aBayes + UFBoot2 1000 --bnni, por loco
              │    └─ ASTER_WASTRAL .. árvore de espécies ponderada pelo aBayes (≥ 3 locos)
              └─ CONCATENA ........... supermatriz + partições NEXUS + ocupação táxon × loco
                   └─ IQTREE_CONCAT .. MFP+MERGE (rclusterf), mesmos suportes (com --amostragem_concat: MFP por loco, sem fusão)
                        └─ IQTREE_CONCORDANCIA ... gCF + sCFL (100 quartetos)
RESUMO ... resumo.tsv, resumo.json, metodos.md (versões do versoes.yml, DOIs conferidos)
```

## Por que um pipeline local

A regra da stack é nf-core antes de pipeline caseiro. Em 2026-10-11 o
catálogo do nf-core não tinha pipeline publicado para esse caminho:
`genephylomodeler` está em desenvolvimento sem release, `multiplesequencealign`
(1.1.1) para no alinhamento, `reportho` (1.1.0) parte de uma proteína de
consulta (ID UniProt ou FASTA) e devolve listas de ortólogos, sem árvore (a
1.1.0 removeu MSA e filogenia), e `phyloplace` (2.3.0) faz posicionamento,
não inferência. Então o fluxo é
local, mas **toda ferramenta que tem módulo nf-core vem do módulo nf-core**
(`mafft/align`, `famsa/align`, `clipkit`, `trimal`, `iqtree`), instalado com
`nf-core modules install` e rastreado em `modules.json`. Só o que não existe
no nf-core é módulo local (`aster_wastral`, `iqtree_concordancia`,
`prepara_loco`, `concatena`, `resumo`).

O módulo `iqtree` tem um **patch local registrado** pelo `nf-core modules
patch` em `modules/nf-core/iqtree/iqtree.diff`: IQ-TREE 3.1.3 → 3.1.4 e
`-nt AUTO` → `-T task.cpus`. Os dois motivos estão abaixo.

## Como rodar

Pelo invólucro, do Windows (é o caminho normal):

```powershell
uv run python scripts/filogenia/filogenia.py --input data/proc/meus_locos --outgroup Homo_sapiens
```

Ele escreve o comando em `runs/filogenia/<id>/executar.sh`, roda no WSL com
`bash -l` (as três armadilhas do CLAUDE.md), desenha as figuras e grava uma
linha de proveniência em `runs/provenance.jsonl`. Opções: `--alinhador
famsa|nenhum`, `--limpeza trimal|nenhuma`, `--tipo-seq DNA|AA`,
`--nome-taxon primeiro|completo|colchete`, `--cpus-concatenado N`,
`--retomar <id>`, e qualquer argumento extra do Nextflow depois de `--`.
`--input` aceita caminho do Windows (`C:\...`), do Git Bash (`/c/...`), do
Explorer (`\\wsl.localhost\Ubuntu-24.04\...`) e do WSL (`/home/...`,
`~/...`). Um CSV com caminhos do Windows é reescrito em
`runs/filogenia/<id>/entrada.csv` com os caminhos já do WSL.

Direto no WSL:

```bash
nextflow run /mnt/c/Users/Manec/bio-stack/pipelines/filogenia -profile docker \
    --input /caminho/locos --outgroup Homo_sapiens --outdir resultados -w ~/nf-work/filogenia/work
```

Teste (5 locos simulados, ~1 minuto depois que as imagens estão no cache;
na primeira vez o Docker baixa ~0,3 GB, que ocupam ~1 GB no disco):

```powershell
uv run --script tests/check_filogenia.py
```

## Entrada

- **Diretório** com um FASTA por loco (`.fa .fasta .fas .fna .faa .aln .fsa`,
  opcionalmente `.gz`, sem distinguir maiúscula — `loco.FASTA` entra; o que
  não é FASTA fica de fora com aviso no log; arquivo oculto (`.*`) e
  subpasta ficam de fora sem aviso), **um FASTA** (aceito com qualquer
  extensão, que entra no nome do loco: `ENAM.txt` vira `ENAM_txt`), um
  **glob**, ou
  um **CSV** com colunas
  `loco,fasta` — caminho relativo no CSV é relativo ao próprio CSV, e `~/`
  vale (o pipeline expande; o `file()` do Nextflow não, medido). O nome do
  loco vem do arquivo (ou da coluna `loco` do CSV) e é saneado para
  `[A-Za-z0-9_]` — vira prefixo de
  arquivo, charset do NEXUS (que não aceita hífen) e chave do resumo; dois
  locos que colidem depois do saneamento são erro, antes de qualquer etapa.
- **Nome do táxon** é o que casa as linhas da supermatriz entre locos.
  `--nome_taxon primeiro` (padrão) pega o primeiro token do cabeçalho;
  `completo` usa o cabeçalho inteiro; `colchete` pega o organismo nos dois
  formatos do NCBI — o `[organism=…]` do NCBI Datasets (`>NP_114095.2 ENAM
  [organism=Homo sapiens] [GeneID=10117]` → `Homo_sapiens`) ou o último
  colchete sem chave do RefSeq clássico (`>XP_015134201.2 enamelin [Gallus
  gallus]` → `Gallus_gallus`). Colchete só com outra chave (`[GeneID=…]`,
  `[isoform=…]`) é erro: pegá-lo daria um nome por gene, e nenhum táxon
  casaria entre locos. Caracteres que o IQ-TREE recusaria viram `_`.
- **Um registro por táxon por loco.** Nome repetido depois do saneamento é
  erro, com a mensagem dizendo quais cabeçalhos colidiram.
- **Tipo** detectado (DNA se ≥ 90 % dos caracteres informativos forem
  ACGTUN) ou forçado com `--tipo_seq`, que **não traduz**: forçar `AA` sobre
  DNA é recusado — A, C, G, T e N também são aminoácidos, e o IQ-TREE
  ajustaria um modelo de proteína a nucleotídeos sem reclamar. Para a árvore
  de proteína, traduza antes (`seqkit translate -T <tabela> --trim`). `*` de
  proteína sai da sequência desalinhada e vira gap na já alinhada (as
  colunas ficam); `*` interno, ou registro que fica vazio e sai, gera aviso
  em `resumo.json → avisos_preparo` e pendência no `metodos.md` (pseudogene
  ou erro de tradução — o stderr da tarefa o Nextflow não mostra).
- Loco com menos de `--min_taxa` (padrão 4: abaixo disso não há ramo
  interno, e o UFBoot não faz sentido — o IQ-TREE só avisa) sai com aviso
  e aparece em `resumo.json → locos_fora`.

## Parâmetros

| Parâmetro | Padrão | O que faz |
|---|---|---|
| `--alinhador` | `mafft` | `mafft`, `famsa` ou `nenhum` (entrada já alinhada) |
| `--limpeza` | `clipkit` | `clipkit`, `trimal` ou `nenhuma` |
| `--clipkit_modo` | `smart-gap` | modo do ClipKIT |
| `--modelo` | `MFP` | `-m` do IQ-TREE; na supermatriz vira `MFP+MERGE` (só `MFP` com `--amostragem_concat`); modelo fixo (ex.: `GTR+F+G4`) não roda o ModelFinder |
| `--ufboot` / `--bnni` | `1000` / `true` | réplicas UFBoot2 (0 desliga; entre 1 e 999 é recusado) e otimização NNI das réplicas |
| `--alrt` | `1000` | réplicas SH-aLRT (0 desliga) |
| `--abayes` | `true` | suporte aBayes; é o peso que o wASTRAL usa |
| `--amostragem_concat` | vazio | reamostragem do UFBoot na supermatriz: vazio = sítios (sem `--sampling`); `GENE` ou `GENESITE` reamostram locos e desligam o MERGE |
| `--seed` | `12345` | semente de IQ-TREE, sCFL e wASTRAL |
| `--outgroup` | — | táxon(s) para enraizar a apresentação |
| `--min_taxa` | `4` | mínimo de táxons por loco |
| `--cpus_concatenado` | `1` | threads da supermatriz (> 1 perde reprodução bit a bit) |
| `--concatenar` / `--concordancia` | `true` | supermatriz e gCF/sCFL (≥ 2 locos) |
| `--arvore_especies` / `--min_locos_especies` | `true` / `3` | wASTRAL |

## Saída (`runs/filogenia/<id>/`)

| Caminho | Conteúdo |
|---|---|
| `resumo.tsv` | por loco: táxons, comprimento antes/depois da limpeza, sítios informativos, modelo, lnL, UFBoot médio, fração de ramos fortes |
| `resumo.json` | tudo acima + supermatriz, concordância, árvore de espécies, versões, parâmetros |
| `metodos.md` | parágrafo de Métodos em português (cada frase só se a etapa rodou), referências com DOI e pendências: versão não registrada, supermatriz com mais de uma thread, limpeza agressiva, dado alterado no preparo, outgroup ausente de todos os locos, ocupação táxon × loco abaixo de 50 % e nenhum táxon em mais de um loco — os três últimos denunciam `--nome_taxon` errado, que roda com exit 0 |
| `figuras/` | PNG + SVG de cada árvore, enraizada no outgroup (ponto médio se ele faltar), ramos fracos em cinza com o valor; o suporte é lido por bipartição antes de enraizar e reatribuído depois. A supermatriz tem legenda própria (julgar por gCF/sCFL) e a do wASTRAL sai sem escala |
| `entrada.csv` | só com CSV de entrada: o CSV que rodou, com os caminhos convertidos para o WSL |
| `01_entrada/` … `07_arvore_especies/` | intermediários de cada etapa |
| `pipeline_info/` | `versoes.yml`, `parametros.json`, `bco.json` (BioCompute Object), trace, timeline, relatório, DAG |
| `ro-crate-metadata.json`, `nextflow.config`, `main.nf`, `README.md` | Workflow Run RO-Crate do nf-prov: os metadados do crate; a configuração resolvida da execução (gerada pelo nf-prov, não é cópia do `nextflow.config` do pipeline); cópia do `main.nf` que rodou (a `mainEntity` do crate, que pode diferir do `pipelines/filogenia/main.nf` atual); e cópia do `README.md` do pipeline (não descreve a execução). **Fica inválido enquanto o `manifest` não tiver `license`** — o repositório ainda não tem licença; escolha uma e ponha em `manifest.license` (SPDX) no `nextflow.config` |
| `executar.sh` | o comando exato que rodou |

## Reprodutibilidade — o que foi medido

Duas execuções do zero (work dirs separados, sem `-resume`) do perfil de
teste produzem árvores, alinhamentos e tabelas idênticos byte a byte; o que
difere é relógio, tempo de CPU e caminho (números abaixo). Para chegar nisso
foram necessárias três medidas, cada uma medida antes:

1. **MAFFT `--threadit 0`.** O refinamento iterativo multithread do MAFFT
   não é determinístico: `--thread 8` no L-INS-i deu 2 alinhamentos
   distintos em 3 rodadas do mesmo loco; com `--threadit 0`, 1 em 3. A fase
   progressiva continua paralela.
2. **IQ-TREE `-T 1` em vez de `-T AUTO`.** Com `AUTO` o número de threads
   depende do tempo medido na hora. Com semente fixa, `-T 1` reproduziu o
   `.treefile` byte a byte em 3 rodadas. Com `-T 2`, no loco_A do teste
   (build Linux), topologia, lnL e suportes se mantiveram entre rodadas e o
   comprimento de ramo variou na 8ª casa decimal (até 5e-8); em relação ao
   `-T 1`, o SH-aLRT mudou em 6 de 9 ramos (até 1,3 ponto), com UFBoot e
   aBayes iguais. A pesquisa de 2026-10-11 mediu no build
   Windows, em `example.phy`, UFBoot e SH-aLRT variando 1 a 2 pontos entre
   rodadas com `-T 2` e `-T 4`. Como as árvores de genes já rodam em paralelo
   entre locos, 1 thread por loco não custa nada; `--cpus_concatenado > 1`
   faz o `metodos.md` declarar essa variação.
3. **IQ-TREE 3.1.4**, não o 3.1.3 do módulo. O 3.1.4 (2026-09-10) corrige
   um `setRootNode()` que abortava quando o outgroup falta numa partição
   (iqtree3#207) — exatamente o caso de uma supermatriz em que algum loco não
   tem o grupo externo. O teste tem um loco assim (`loco_D`).

FAMSA (`-t 8`) deu alinhamento idêntico em 3 rodadas. O wASTRAL (semente
fixa) repete o arquivo byte a byte em 3 rodadas com o mesmo `-t` (medido com
1, 2, 4 e 8) e não precisou de ajuste: o módulo passa `-t ${task.cpus}`, fixo
na configuração (2), e não escolhido na hora como o `-T AUTO` do IQ-TREE.
Mas o `especies.tre` muda com o `-t`, sem mudar topologia nem suporte: com
`-t 2` o ramo tartarugas + arcossauros sai `0.944625:0.94446`; com `-t 1`,
`4` e `8`, `0.944625:nan`. Quem mudar o `cpus` do ASTER_WASTRAL obtém outro
arquivo. Remedido em 2026-10-11 com a
configuração padrão atual (aBayes ligado, supermatriz sem `--sampling`): duas
execuções do zero deram 52 de 70 arquivos idênticos — todas as árvores,
alinhamentos, esquemas de partição, fatores de concordância, `resumo.json` e
`metodos.md`. Os 18 que diferem são `.iqtree` e `.log` (data, hora, tempo de
CPU, host) e os relatórios de execução que embutem relógio ou o caminho do
outdir (`trace`, `timeline`, relatório HTML, `bco.json`, o RO-Crate e o
`nextflow.config` que ele copia). Uma curiosidade: o `.iqtree` da supermatriz
traz na coluna w-AICc um número subnormal (2,2e-314 numa execução, 4,1e-314
na outra) — lixo de memória no relatório do IQ-TREE, sem efeito em nada.

**`--sampling` muda a verossimilhança das partições — e, com MERGE, o que o
UFBoot sorteia.** Medido na supermatriz de teste, IQ-TREE 3.1.4, mesma
semente:

- Passar `--sampling`, com qualquer valor, desliga a remoção de sequências
  vazias por partição (`remove_empty_seq = false`, em `utils/tools.cpp`), e
  a verossimilhança das partições muda: com `--sampling SITE` ou `GENESITE`
  o MFP+MERGE (rclusterf) parou em 3 partições (BIC 48 147,6); sem a opção,
  fundiu as 4 em 1 (BIC 48 167,7). Os BIC não se comparam entre si — a
  verossimilhança é calculada de outro jeito. `SITE` explícito não equivale
  a omitir a opção, e com `--bnni` aborta (`ERROR: Wrong -bsam`): por isso
  `SITE` não é valor aceito, e vazio já é reamostragem de sítios.
- Depois do MERGE o UFBoot sorteia os blocos **já fundidos**, não os locos
  (`SuperAlignment::createBootstrapAlignment`). Com tudo fundido em 1,
  `GENESITE` equivale a `SITE` e `GENE` dá 100 em todo ramo, seja qual for o
  dado — com 8 locos simulados em conflito, `GENE` deu 100/100 com MERGE e
  63/63 sem. Por isso, com `--amostragem_concat`, a supermatriz roda MFP por
  loco, sem fusão, e a unidade sorteada é o loco.

O `metodos.md` registra o esquema que de fato saiu, não o esperado.

## Validação

Os locos de teste (`assets/test/`) foram simulados com o AliSim do IQ-TREE
3.1.4 sobre `arvore_verdadeira.nwk` (12 amniotas: aves, crocodilianos,
tartarugas, *Anolis*, *Homo* como externo), com indels, e depois
desalinhados. Cada um exercita um caminho: `loco_B` tem cabeçalho com
descrição (saneamento de nome), `loco_C` não tem dois táxons (preenchimento
da supermatriz), `loco_D` não tem o outgroup, `loco_E` tem 3 táxons (tem de
ser descartado). `tests/check_filogenia.py` roda o perfil de teste e
confere contra a árvore verdadeira. Resultado em 2026-10-11, configuração
padrão (aBayes ligado, supermatriz sem `--sampling`), 66 s:

| Verificação | Resultado |
|---|---|
| árvore concatenada vs. verdadeira | RF 0/18 |
| árvore de espécies (wASTRAL) vs. verdadeira | RF 0/18 |
| árvores de genes A, B, C, D | RF 0/18, 0/18, 2/14, 0/16 |
| `loco_E` descartado | sim, `poucos_taxa` |
| gCF médio / sCFL médio | 97,2 / 64,1 |
| esquema de partição | 4 partições de entrada → 1 após MERGE (HKY+G4) |
| `metodos.md` | IQ-TREE v3.1.4; pendências: nenhuma |

O `loco_C` errar a posição de *Chrysemys* (SH-aLRT 76,7 / aBayes 0,70 /
UFBoot 60 no ramo errado) é o esperado para 751 sítios com dois táxons
faltando — e a figura mostra o ramo em cinza, com os três valores.

## Escolhas que vieram da pesquisa (medidas, não opinião)

- **aBayes para o wASTRAL.** O wASTRAL lê só o **primeiro** número de um
  rótulo `a/b` — com o rótulo padrão do IQ-TREE ponderaria pelo SH-aLRT sem
  avisar, e aborta se o rótulo começar com `/` (medido no v1.25). O módulo
  extrai o campo certo; e o campo é o aBayes porque o artigo do wASTRAL
  (Zhang & Mirarab 2022) mostrou melhor acurácia com o aBayes do IQ-TREE —
  UFBoot não foi testado lá. Não colapse ramos antes: o peso faz isso melhor.
- **Reamostragem de locos é opcional** (`--amostragem_concat GENE|GENESITE`):
  `GENE` sorteia locos inteiros em vez de sítios (Nei et al. 2001, citado no
  tutorial avançado do IQ-TREE); `GENESITE` sorteia locos e depois sítios
  dentro deles, o que, segundo o mesmo tutorial, "pode ajudar a reduzir
  falsos positivos" (ramo errado com 100 %). Medido aqui, com 8 locos
  simulados em conflito 4 contra 4 e uma partição por loco: reamostrando
  sítios o ramo em conflito teve 97; com `GENE` ou `GENESITE`, 63. Só tem
  sentido com uma partição por loco — por isso liga MFP sem MERGE. O padrão
  é vazio: sítios, com MERGE.
- **Limpeza leve.** Filtro automático agressivo piora árvore de gene único
  (Tan et al. 2015); no benchmark do ClipKIT, "não aparar" ficou entre os
  melhores. Loco com mais de 20 % dos sítios removidos vira pendência no
  `metodos.md` (`resumo.json → locos_limpeza_agressiva`). Evite todos os
  modos `kpi*` (`kpi`, `kpic` e as variantes `-gappy`/`-smart-gap`) quando
  o comprimento de ramo importa (datação, saturação, dN/dS): todos removem
  singletons, e `kpi`/`kpi-*` também os sítios constantes. Medido com a
  topologia verdadeira fixa (GTR+F+G4, ClipKIT 2.11.4): o `kpic` encurta os
  terminais (0,65 a 0,74 vezes o verdadeiro) e alonga os internos (1,6 a 2,0
  vezes); o `kpi` alonga todos. O `+ASC` não corrige, porque também faltam
  os singletons.
- **Nomes de táxon.** O IQ-TREE estende o nome palavra a palavra, até 4,
  só enquanto houver duplicata **no arquivo**: num loco com *Gallus gallus*
  e *Gallus varius* saem `Gallus_gallus` e `Gallus_varius`; num loco sem
  duplicata, `Gallus gallus` vira `Gallus` (medido no 3.1.4; `alignment.cpp`).
  `: , ( ) [ ] ; '` viram `_` com só um WARNING. Dentro de um loco não há
  colisão — o risco é o mesmo táxon ganhar nomes diferentes em locos
  diferentes, e a supermatriz, o gCF e o ASTER deixarem de casar as linhas.
  Por isso o `PREPARA_LOCO` saneia antes, com a mesma regra em todos os
  locos, e registra o mapa em `01_entrada/*.entrada.json`.
- **Sequências idênticas.** O IQ-TREE mantém duas e recoloca as demais com
  ramo zero; isso entra no gCF e no ASTER como resolução arbitrária.
  Deduplique por espécie antes se acontecer muito.

## Como ler o resultado

- **Árvores de genes:** ramo confiável é **SH-aLRT ≥ 80 e UFBoot ≥ 95**. A
  FAQ do IQ-TREE restringe esse critério a árvore de gene único; UFBoot não
  é bootstrap padrão.
- **Supermatriz:** com muitos sítios o suporte satura em 100 até em ramo que
  metade dos genes contradiz. Julgue cada ramo pelo gCF (Minh, Hahn &
  Lanfear 2020) e pelo sCFL (Mo et al. 2023); a figura da supermatriz tem legenda própria dizendo
  isso. Suporte não é concordância: UFBoot 100 com gCF 30 é a supermatriz
  inteira preferindo um ramo que só 30 % dos genes contêm.
- **Discordância** entre concatenada e wASTRAL é informação, não erro do
  pipeline — reporte as duas. A causa pode ser ILS, paralogia, erro de
  estimativa das árvores de genes ou erro sistemático. gCF ≈ sCF com
  gDF1 ≈ gDF2 é compatível com ILS; gCF ≪ sCF aponta sinal fraco por loco
  (árvores de genes mal estimadas), não conflito biológico (Lanfear & Hahn
  2024); gDF1 ≫ gDF2 pede introgressão, paralogia ou atração de ramos
  longos. Sob ILS a concatenação pode convergir com suporte alto para a
  árvore errada, dentro e fora da zona de anomalia (Kubatko & Degnan 2007;
  Mendes & Hahn 2018) — zona de anomalia é a região em que a árvore de gene
  mais frequente difere da de espécies, não o nome dessa inconsistência.
- **Comprimento do wASTRAL não é unidade coalescente.** O pipeline roda o
  modo híbrido (`-B`, `--mode 1`), cujo comprimento mistura unidade
  coalescente e substituição (CU + 2×SU — tutorial do ASTER; Zhang &
  Mirarab 2022). Terminais vêm sem comprimento, e em qualquer `--mode` todo
  ramo interno com q1 = 1 — nenhum quarteto discordante — sai `:nan`: o
  comprimento de máxima verossimilhança não tem limite (7 dos 9 internos no
  teste, que tem só 4 árvores de genes). Por isso a figura vira cladograma.
  Reporte esses ramos como "comprimento não estimável (q1 = 1)", nunca como
  0 — em unidade coalescente, 0 é discordância máxima. O `ape::read.tree`
  aceita o `nan` em silêncio e põe `NaN` em `edge.length`, mesmo depois de
  tirar o `:nan` do texto. Só topologia (RF, figura): `tr$edge.length <- NULL`.
  Unidade coalescente: `wastral --mode 4` (sem peso; `-u 2` põe `q1` e
  `CULength` no rótulo) — o `:nan` dos ramos com q1 = 1 continua. Substituições por
  sítio: reestime na topologia fixa, da raiz de `runs/filogenia/<id>/` —
  `sed 's/:nan//g' 07_arvore_especies/especies.tre > 07_arvore_especies/especies.sem_nan.tre`
  e `iqtree3 -s 05_concatenado/supermatriz.fa -p 05_concatenado/concatenado.best_model.nex -te 07_arvore_especies/especies.sem_nan.tre -T 1 --seed 12345 --prefix 07_arvore_especies/especies.sl`;
  a árvore sai em `07_arvore_especies/especies.sl.treefile` (com o `:nan` cru
  o IQ-TREE aborta em "Numerical underflow" — medido).

O agente `phylo-analyst` tem essas regras; o comando `/filogenia` é a ordem
das etapas.

## Versões e containers

| Ferramenta | Versão | Container |
|---|---|---|
| IQ-TREE | 3.1.4 | `quay.io/biocontainers/iqtree:3.1.4--h8471819_0` |
| MAFFT | 7.520 | `quay.io/biocontainers/mulled-v2-12eba4a0…:42570789…-0` (módulo nf-core) |
| FAMSA | 2.4.1 | `quay.io/biocontainers/famsa:2.4.1--h9ee0642_0` |
| ClipKIT | 2.11.4 | `quay.io/biocontainers/clipkit:2.11.4--pyhdfd78af_0` |
| trimAl | 1.5.0 | `quay.io/biocontainers/trimal:1.5.0--h9948957_2` |
| ASTER (wASTRAL) | 1.25 (v1.25.3.8) | `quay.io/biocontainers/aster:1.25--hc52dbad_1` |
| Python (scripts) | 3.12.12 | `quay.io/biocontainers/python:3.12.12@sha256:4cc84261…` |
| Nextflow | 26.04.6 | WSL, `~/.local/bin` |

Os módulos nf-core fixam MAFFT 7.520 e ClipKIT 2.11.4; o bioconda tem 7.525 e
2.14.0. Ficam os do módulo até o nf-core atualizar — trocar por patch só se
houver motivo medido, como no IQ-TREE.

## Limites e próximos passos

- **Alinhamento consciente de códon** não está no fluxo: para testes de
  seleção (HyPhy, PAML) o alinhamento tem de respeitar o quadro de leitura.
  MACSE 2.07, HyPhy 2.5.103 e PAML 4.10.10 já estão na camada pixi do WSL
  (`bootstrap/pixi-global.toml`); o nf-core tem `macse/refinealignment` e
  `hyphy/*` para quando isso virar etapa.
- **Ortologia** é responsabilidade da entrada. Para montar locos a partir de
  genomas, o OrthoFinder vem antes deste pipeline — a 3.1.5 está na camada
  pixi (`orthofinder -f <dir de proteomas>`), e o nf-core tem módulo.
- **Datação** (LSD2 no próprio IQ-TREE, MCMCtree) e **reconstrução
  ancestral** (`--ancestral`) não estão ligadas; são extensões naturais do
  `IQTREE_CONCAT`.
- **TreeShrink** para ramos terminais anômalos (parálogo, contaminação)
  também não.

## Referências

Todas conferidas no Crossref em 2026-10-11.

- Wong TKF, Ly-Trong N, Ren H, et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Mol Biol Evol* 43(5):msag117. doi:10.1093/molbev/msag117
- Kalyaanamoorthy S et al. (2017) ModelFinder. *Nat Methods* 14:587–589. doi:10.1038/nmeth.4285
- Hoang DT et al. (2018) UFBoot2. *Mol Biol Evol* 35:518–522. doi:10.1093/molbev/msx281
- Guindon S et al. (2010) PhyML 3.0 / SH-aLRT. *Syst Biol* 59:307–321. doi:10.1093/sysbio/syq010
- Chernomor O, von Haeseler A, Minh BQ (2016) Partition models. *Syst Biol* 65:997–1008. doi:10.1093/sysbio/syw037
- Minh BQ, Hahn MW, Lanfear R (2020) Concordance factors. *Mol Biol Evol* 37:2727–2733. doi:10.1093/molbev/msaa106
- Mo YK et al. (2023) Updated site concordance factors. *Bioinformatics* 39(1):btac741. doi:10.1093/bioinformatics/btac741
- Katoh K, Standley DM (2013) MAFFT v7. *Mol Biol Evol* 30:772–780. doi:10.1093/molbev/mst010
- Deorowicz S et al. (2016) FAMSA. *Sci Rep* 6:33964. doi:10.1038/srep33964
- Steenwyk JL et al. (2020) ClipKIT. *PLoS Biol* 18:e3001007. doi:10.1371/journal.pbio.3001007
- Capella-Gutiérrez S et al. (2009) trimAl. *Bioinformatics* 25:1972–1973. doi:10.1093/bioinformatics/btp348
- Zhang C, Mirarab S (2022) wASTRAL. *Mol Biol Evol* 39:msac215. doi:10.1093/molbev/msac215
- Zhang C, Nielsen R, Mirarab S (2025) ASTER. *Mol Biol Evol* 42(8):msaf172. doi:10.1093/molbev/msaf172
- Anisimova M et al. (2011) aBayes / survey de suporte. *Syst Biol* 60:685–699. doi:10.1093/sysbio/syr041
- Lanfear R et al. (2014) rcluster (fusão de partições). *BMC Evol Biol* 14:82. doi:10.1186/1471-2148-14-82
- Lanfear R et al. (2017) PartitionFinder 2 — o rclusterf rápido que o `MFP+MERGE` usa. *Mol Biol Evol* 34:772–773. doi:10.1093/molbev/msw260
- Kubatko LS, Degnan JH (2007) Inconsistência da concatenação sob coalescência. *Syst Biol* 56:17–24. doi:10.1080/10635150601146041
- Mendes FK, Hahn MW (2018) Why concatenation fails near the anomaly zone. *Syst Biol* 67:158–169. doi:10.1093/sysbio/syx063
- Lanfear R, Hahn MW (2024) The meaning and measure of concordance factors in phylogenomics. *Mol Biol Evol* 41:msae214. doi:10.1093/molbev/msae214
- Tan G et al. (2015) Limpeza automática e árvores de gene único. *Syst Biol* 64:778–791. doi:10.1093/sysbio/syv033
- Ly-Trong N et al. (2022) AliSim (dados de teste). *Mol Biol Evol* 39:msac092. doi:10.1093/molbev/msac092
- Di Tommaso P et al. (2017) Nextflow. *Nat Biotechnol* 35:316–319. doi:10.1038/nbt.3820
- Ewels PA et al. (2020) nf-core. *Nat Biotechnol* 38:276–278. doi:10.1038/s41587-020-0439-x
