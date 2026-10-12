# bio-stack

Stack de bioinformática e biologia sintética computacional para Claude Code,
com aceleração via NVIDIA BioNeMo NIM.

Montada para trabalho acadêmico: análise de dados ômicos, biologia estrutural
computacional, desenho de construtos in silico, engenharia de pipelines.

---

## O diagnóstico que define esta stack

Antes de recomendar qualquer coisa, verifiquei a máquina:

```
GPU        NVIDIA GeForce RTX 5070 Laptop · 8 GB VRAM · CUDA 13.4 · sm_120
Presentes  nvidia-smi, wsl (sem distro)
Ausentes   git, node, npm, npx, uv, python (só o stub da Store), docker, conda
```

Esse é o retrato de antes do bootstrap. Depois dele, `docker` **continua**
ausente no PATH do Windows — e isso é projeto, não lacuna pendente: o daemon é
o `docker-ce` nativo dentro da `Ubuntu-24.04`, sob systemd, sem Docker Desktop
no Windows. Container é chamado atravessando a fronteira, com
`wsl -d Ubuntu-24.04 -- bash -lc "docker ..."`, e `-profile docker` é o perfil
padrão do Nextflow aqui — não precisa de conda nem singularity.

Duas consequências que mudam tudo:

**1. Nada está instalado.** O sistema de plugins do Claude Code precisa de
`git`; o BioNeMo Agent Toolkit precisa de `node`; os servidores MCP em Python
precisam de `uv`. Por isso a stack começa com uma camada de bootstrap — não é
burocracia, é pré-requisito.

**2. 8 GB de VRAM não rodam modelo de estrutura.** ESMFold quer ~14 GB,
Evo 2 quer ~40 GB, Parabricks exige ≥16 GB. A solução não é otimizar: é usar
os **endpoints hospedados** da NVIDIA, onde a sua licença AI Enterprise já
cobre o compute. A GPU local fica com o que ela faz bem — ESM-2 650M,
DNABERT-2, RDKit, nvMolKit, pré-processamento.

---

## Instalação

```powershell
# 1. Ferramentas de base (uma vez; reinicie o terminal depois)
.\bootstrap\00_bootstrap_windows.ps1

# 2. Chaves de API (valida a chave NVIDIA ao vivo)
.\bootstrap\01_setup_keys.ps1

# 3. Ambiente Python
uv sync --python 3.11 --locked
```

O passo 3 é o caminho canônico: `pyproject.toml` declara `gpu` como grupo em
`[dependency-groups]` — ligado por `default-groups` junto com `dev`,
`evolucao`, `estatistica`, `pesquisa` e `notebooks`, então o `uv sync`
instala todos — e aponta `torch` para o índice `pytorch-cu128` em
`[[tool.uv.index]]`. Não é extra: `uv sync --extra gpu` falha com
"Extra `gpu` is not defined".
Esse apontamento não é preferência — não existe wheel oficial de PyTorch para
sm_120 abaixo de cu128, e um wheel cu121 instala sem reclamar e só falha no
primeiro kernel. O `--python 3.11` não é detalhe: sem `.python-version` no
repositório, o `uv sync` pelado num clone sem `.venv` usa o primeiro Python
que satisfaz `>=3.11` (nesta máquina, o CPython 3.13 gerenciado pelo uv), e
como o `uv.lock` bifurca em 3.12 o ambiente sai outro (pandas 3, anndata
0.13). `bootstrap/02_python_env.ps1` roda exatamente esse comando, para quem
prefere o script.

Reinicie o Claude Code, abra esta pasta como projeto, e siga
[`bootstrap/03_install_plugins.md`](bootstrap/03_install_plugins.md) para as
skills. Depois rode `/bio-doctor` para conferir tudo.

---

## O que a stack tem

### Servidores MCP

| Servidor | Cobre |
|---|---|
| **`nvidia-bio`** — código próprio | NIMs de biologia: dobramento, cofolding, docking, MSA, Evo 2 |
| **`nvidia-nemo`** — código próprio | NeMo Retriever: busca semântica na sua biblioteca de PDFs, embeddings, reranking |
| `biocontext` | STRING, Open Targets, Reactome, UniProt, HPA, KEGG, AlphaFold DB, Ensembl, InterPro (14 bases) |
| `gget` | busca de gene, BLAST, BLAT, MUSCLE, PDB, Enrichr, ARCHS4, COSMIC, CELLxGENE |
| `biomcp` | ensaios clínicos, ClinVar, MyVariant, OncoKB, openFDA (~30 fontes) |
| `tooluniverse` | metacatálogo: busca e executa ferramentas de outros servidores por descrição |
| **`claude-science`** — código próprio | opera o daemon local do Claude Science pela API REST dele (varreduras longas, 24 conectores) |

Mais os conectores já ligados na sua conta: PubMed, bioRxiv, BioRender,
Synthesize Bio, Hugging Face.

### Por que um servidor MCP próprio para a NVIDIA

O BioNeMo Agent Toolkit da NVIDIA entrega **skills** — excelentes para ensinar
o procedimento científico. Mas cada chamada vira um `curl` improvisado pelo
modelo, sem retry, sem o polling assíncrono que o NVCF exige, e com a resposta
caindo inteira no contexto. Um mmCIF passa de 100 mil tokens.

`servers/nvidia_bio_mcp/server.py` resolve isso: trata o 202 + `nvcf-reqid`,
valida sequência antes de gastar uma chamada, traduz os códigos de erro
(403 = licença não cobre, 404 = slug mudou), **grava estruturas em disco e
devolve só o caminho**, e registra tudo para reprodutibilidade.

Os dois convivem: as skills dizem *o que* fazer, o MCP *executa* sem estourar
o contexto.

### Agentes

`structure-bio` · `omics-analyst` · `seq-designer` · `variant-curator` ·
`pipeline-engineer` · `lit-scout` · `evo-devo` · `phylo-analyst`

Cada um escrito para este hardware — o `structure-bio` sabe que não deve
propor AlphaFold local, o `pipeline-engineer` sabe que Nextflow precisa de WSL2.

### Comandos

`/bio-doctor` · `/nim` · `/nvidia` · `/fold` · `/metodos` · `/atavismo` · `/filogenia`

### Pipeline de filogenia (IQ-TREE 3)

`pipelines/filogenia` — Nextflow DSL2 com módulos nf-core, roda no WSL com
Docker: locos não alinhados → MAFFT/FAMSA → ClipKIT/trimAl → IQ-TREE 3 por
loco (ModelFinder + SH-aLRT + aBayes + UFBoot2) → supermatriz MFP+MERGE →
gCF/sCFL → wASTRAL → `resumo.tsv`, figuras e o parágrafo de Métodos com as
versões que rodaram. Árvores, alinhamentos, partições, fatores de
concordância, `resumo.json` e `metodos.md` saem idênticos byte a byte entre
duas execuções do zero (medido no perfil de teste, configuração padrão; com
`--cpus_concatenado` > 1 isso se perde) — 52 de 70 arquivos; os outros
carimbam data, tempo de CPU, host ou caminho. Validado contra árvore
simulada conhecida (RF = 0 na concatenada e na de espécies). Guia: [`docs/FILOGENIA.md`](docs/FILOGENIA.md).

```powershell
uv run python scripts/filogenia/filogenia.py --input <dir com 1 FASTA por loco> --outgroup <taxon>
```

### CLIs de bioinformática no WSL (pixi)

37 ambientes pixi (um por ferramenta, ou grupo coeso), versão exata em
`bootstrap/pixi-global.toml` e build + sha256 em
`bootstrap/pixi-global.snapshot.tsv`: filogenia (IQ-TREE, MAFFT, FAMSA,
MUSCLE 5, ClipKIT, trimAl, RAxML-NG, FastTree, ASTER), evolução molecular
(HyPhy, PAML, MACSE, PRANK), homologia e ortologia (BLAST+, DIAMOND, HMMER,
MMseqs2, OrthoFinder 3.1.5, miniprot, minimap2), transcrição (subread 2.1.1
com featureCounts, StringTie 3.0.3), genomas e dados
(samtools/bcftools 1.24, seqkit, csvtk, taxonkit — com a base do NCBI que o
bootstrap baixa, data em `~/.taxonkit/BAIXADO_EM` —, NCBI Datasets, EDirect, gffread, AGAT, utilitários UCSC,
jq, DuckDB, parallel) e nf-test. `bash bootstrap/05_wsl_ferramentas.sh` (no WSL) reaplica tudo.
O levantamento que escolheu cada uma está em
[`docs/PESQUISA-STACK-2026-10.md`](docs/PESQUISA-STACK-2026-10.md).

### R (renv)

`r/` tem R 4.6.1 + Bioconductor 3.23, com biblioteca do projeto e versões
travadas em `r/renv.lock` (237 pacotes): filogenia comparativa e ASR de
caráter (ape, phangorn, phytools, castor, geiger, corHMM, TreeDist,
phylter), figura de árvore (ggtree, treeio, ggtreeExtra) e expressão
diferencial de referência para conferir o pydeseq2 (DESeq2, edgeR, limma,
sva, apeglm, ashr, IHW, qvalue, variancePartition/`dream`, tximport).
`Rscript r/instalar.R` restaura as versões travadas; armadilhas (o `:nan`
do wASTRAL, entre outras) em [`r/README.md`](r/README.md).

### Hooks

| Hook | O que faz |
|---|---|
| `session_context` | injeta GPU, ferramentas e chaves no contexto ao abrir a sessão |
| `guard_data` | bloqueia `rm -rf` em dados brutos, push forçado, curl-pipe-bash |
| `lint_gate` | roda o linter do arquivo recém-editado e devolve os erros ao modelo — cobertura real abaixo |
| `provenance` | caderno de laboratório automático em `runs/provenance.jsonl` |
| `session_mark` | marca início, fim e compactação de sessão no mesmo log, com `kind="marco"` |

### LSP

O Claude Code **tem** cliente LSP nativo (v2.0.74+), instalável como plugin:

```
/plugin install pyright-lsp@claude-plugins-official
```

Isso cobre Python. O resto — Nextflow, Snakemake, schemas de nf-core, R —
é o que o hook `lint_gate` tenta cobrir, rodando o analisador em modo batch e
devolvendo os diagnósticos pelo canal de feedback (saída 2, texto no stderr).
As duas metades em [`docs/LSP.md`](docs/LSP.md).

### Cobertura real do `lint_gate`

Medida nesta máquina em 2026-10-10. Metade dos ramos está inerte, e o hook não
avisa:

| Arquivo | Analisador | Estado |
|---|---|---|
| `.py` | `ruff check` + `basedpyright` | **funciona** — nenhum dos dois está no PATH do Windows; `uv tool run` resolve e cacheia sob demanda |
| `.smk`, `Snakefile` | `snakemake --lint` | **funciona** pela mesma rota `uv tool run` |
| `.json` | `json.loads` embutido | **funciona** |
| `.yaml`, `.yml` | PyYAML | **inerte** — o cabeçalho PEP 723 do hook declara `dependencies = []`, então `import yaml` levanta `ImportError` e o ramo devolve lista vazia |
| `.nf`, `nextflow.config` | `nextflow lint` | **inerte** — invocado como binário cru, não por `uv tool run`; o nextflow só existe dentro da WSL |
| `.R`, `.r` | `Rscript -e lintr::lint` | **inerte** — o `Rscript` existe, o pacote `lintr` não; o erro sai no stderr e o ramo só lê o stdout |

O silêncio é total por duas razões que se somam. O guarda `have()` devolve
verdadeiro para qualquer ferramenta desde que o `uv` exista, então nunca
reprova nada; e `run()` converte `OSError` em saída −1 com stdout e stderr
vazios, que cada analisador lê como "nada a relatar". Um `.nf` truncado no meio
de um heredoc passa com saída 0:

```console
$ uv run --script .claude/hooks/lint_gate.py < payload_nf.json; echo "exit=$?"
exit=0

$ uv run --script .claude/hooks/lint_gate.py < payload_py.json; echo "exit=$?"
[lint_gate] 6 diagnostico(s) em quebrado.py:
  ...
  quebrado.py:4:12: F821 Undefined name `nao_existe`
  quebrado.py:4: "nao_existe" is not defined
Corrija antes de rodar a analise.
exit=2
```

Um gate que não roda é pior que gate nenhum: ele dá falsa confiança. Para o
ramo Nextflow valer algo aqui, a chamada tem de atravessar a fronteira da WSL
— `wsl -d Ubuntu-24.04 -- bash -lc "nextflow lint ..."` — e não procurar um
binário que, por projeto, não existe no Windows.

---

## As duas peças que rendem mais

**`/metodos`** lê `runs/provenance.jsonl` e escreve a seção de Métodos do seu
relatório a partir do que você *realmente* rodou — versões, parâmetros, datas
de acesso, inclusive as tentativas que falharam. É o registro que todo mundo
promete manter à mão e abandona na terceira semana.

**`lint_gate`** fecha o ciclo de correção **em Python**. Sem ele, um erro de
tipo num script de análise só aparece depois de 40 minutos de pipeline; com
ele, volta na mesma iteração. Nos ramos `.nf`, `.R` e `.yaml` o gate está
inerte hoje — leia a tabela de cobertura antes de confiar nele para um
`main.nf`.

---

## Testes

```bash
powershell -ExecutionPolicy Bypass -File tests/run_all.ps1
```

Nenhuma chamada de rede (com o cache do uv já populado), nenhum crédito da
NVIDIA consumido. Cada suíte imprime uma linha `OK`/`FALHA` por caso e fecha
com o próprio veredito; o `run_all.ps1` termina em `Toda a suite passou.` ou
em `FALHARAM: <suítes>`. Nenhum total de casos é impresso — para tê-lo,
conte as linhas `OK` da saída em vez de confiar num número escrito aqui.

`tests/check_nemo.py` substitui as chamadas de embedding e reranking por um
dublê determinístico, então exercita indexação incremental, fatiamento e o
pipeline de busca inteiro sem gastar um token de crédito.

`tests/check_hooks.py` exercita quatro dos cinco hooks (`guard_data`,
`lint_gate`, `provenance`, `session_context`; o `session_mark` não tem
teste) pelo contrato real do Claude Code
— JSON no stdin, código de saída como veredito. Cobre o `guard_data` nos dois
sentidos (inclusive os comandos que ele **deve** deixar passar, como
`nextflow run` e `blastp`), a redação de chaves de API no log de proveniência,
e a truncagem de campos gigantes.

`tests/smoke_design.py` fica **fora** de `run_all.ps1`: e o unico teste que
chama os NIM de verdade e gasta credito. Rode de proposito, com
`uv run --script tests/smoke_design.py`, para confirmar que os endpoints de
desenho continuam de pe e que os formatos de resposta nao mudaram.

`tests/check_server.py` carrega o servidor MCP e verifica registro de ferramentas,
validação de sequência, mensagens de erro de credencial, e — o mais importante —
que coordenadas atômicas **nunca** voltam inline para o contexto.

Rode depois de qualquer alteração nos hooks ou no servidor.

`tests/check_filogenia_figura.py` cobre, sem WSL e em segundos, a figura de
árvore, o invólucro e as bordas do utilitário do pipeline de filogenia
(cabeçalho do NCBI Datasets, Métodos com um loco só): o suporte segue a bipartição
quando a árvore é reenraizada (a Biopython sozinha desloca 3 de 4 rótulos — o
teste confere isso como controle), rótulo com campo vazio ou com outro número
de campos, o `:nan` do wASTRAL, e a tradução de caminhos para o WSL (Git Bash,
UNC, `~`, CSV).

Duas suítes dependem do WSL e ficam fora do `run_all.ps1` padrão (a
`check_filogenia.py` baixa imagens Docker na primeira vez); entram com
`run_all.ps1 -Completo` e pelo `uv run pytest`:

- `tests/check_wsl_ferramentas.py` — o manifesto pixi e o instalado
  concordam, e cada um dos 37 ambientes tem ao menos um binário executado
  de verdade (41 sondas, não todos os binários expostos), sob `timeout 60` e
  com stdin fechado; a versão fixada é conferida quando o binário a
  imprime. Exceções: do ambiente `kent` rodam só `twoBitToFa`, `liftOver` e
  `mafsInRegion`, que não imprimem versão (vale a linha de uso); o gotree
  também não imprime, e a sonda roda `gotree stats` numa árvore de 4 folhas;
  o PRANK do pacote 251117 se apresenta como `v.250331` (confere-se só
  `PRANK v.`); o mmseqs imprime `18.8cc5c` e o trimAl 1.5.1 `v1.5.rev1`, e
  as sondas conferem `18` e `1.5`.
- `tests/check_filogenia.py` — roda o perfil de teste do pipeline e exige
  RF = 0 entre a árvore inferida e a árvore verdadeira da simulação.

---

## Escopo

Análise computacional acadêmica: genômica, transcriptômica, estrutura de
proteínas, desenho de construtos para organismos de laboratório de ensino,
engenharia de pipelines.

Projetos envolvendo agentes de risco elevado, toxinas, ou aumento de virulência
têm via institucional própria — comitê de biossegurança e CTNBio. Não é o
objeto desta stack, e não é uma questão de configuração.

Na prática, o trabalho de graduação e pós em bioinformática não esbarra nisso:
dobrar uma proteína, montar um plasmídeo repórter, rodar um DESeq2 ou desenhar
primers são tarefas comuns e diretas. O que esta stack adiciona é rastreabilidade
— cada análise fica registrada com parâmetro e data, o que serve tanto para a
sua banca quanto para deixar o contexto do trabalho legível.

---

## Documentação

| Arquivo | Conteúdo |
|---|---|
| [`CLAUDE.md`](CLAUDE.md) | roteamento de agentes, restrições de hardware, padrões |
| [`docs/NVIDIA.md`](docs/NVIDIA.md) | endpoints NIM, protocolo assíncrono, orçamento de VRAM, custo |
| [`docs/STACK.md`](docs/STACK.md) | catálogo completo do ecossistema com links |
| [`docs/CHAVES.md`](docs/CHAVES.md) | inventário de chaves de API, medido ao vivo |
| [`docs/AUTONOMIA.md`](docs/AUTONOMIA.md) | auditoria e roteiro para ambiente autônomo |
| [`docs/LSP.md`](docs/LSP.md) | camada de análise estática e o que o Claude Code não faz |
| [`docs/FILOGENIA.md`](docs/FILOGENIA.md) | pipeline IQ-TREE 3: uso, parâmetros, reprodutibilidade medida, validação |
| [`r/README.md`](r/README.md) | ambiente R (renv): pacotes travados, como restaurar, armadilhas |
| [`docs/PESQUISA-STACK-2026-10.md`](docs/PESQUISA-STACK-2026-10.md) | levantamento verificado de ferramentas, linguagens e recursos |
| [`bootstrap/03_install_plugins.md`](bootstrap/03_install_plugins.md) | skills e plugins, com comandos |

---

## Estrutura

```
bio-stack/
├── CLAUDE.md                    roteamento e contexto do projeto
├── .mcp.json                    7 servidores MCP
├── .claude/
│   ├── settings.json            hooks, permissões, variáveis
│   ├── agents/                  8 agentes especializados
│   ├── commands/                7 slash commands
│   └── hooks/                   5 hooks Python
├── servers/
│   ├── nvidia_bio_mcp/          NIMs de biologia (código próprio)
│   └── nvidia_nemo_mcp/         NeMo Retriever (código próprio)
├── pipelines/filogenia/         Nextflow: locos → IQ-TREE 3 → Métodos
├── scripts/filogenia/           invólucro Windows→WSL e figuras de árvore
├── r/                           R 4.6 + Bioconductor (renv): filogenia comparativa, DE de referência
├── bootstrap/                   instalação em camadas (+ pixi-global.toml no WSL)
├── tests/                       suítes check_*.py; o run_all.ps1 padrão roda sem rede e sem crédito (exceções em Testes)
├── editor/                      configuração do VS Code
├── docs/                        referência
├── data/{raw,proc}/             dados (raw é protegido por hook)
├── runs/                        saídas + provenance.jsonl
└── notebooks/                   análise exploratória
```
