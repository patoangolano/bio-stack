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
.\bootstrap\02_python_env.ps1 -Gpu -Dev
```

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
`pipeline-engineer` · `lit-scout`

Cada um escrito para este hardware — o `structure-bio` sabe que não deve
propor AlphaFold local, o `pipeline-engineer` sabe que Nextflow precisa de WSL2.

### Comandos

`/bio-doctor` · `/nim` · `/nvidia` · `/fold` · `/metodos`

### Hooks

| Hook | O que faz |
|---|---|
| `session_context` | injeta GPU, ferramentas e chaves no contexto ao abrir a sessão |
| `guard_data` | bloqueia `rm -rf` em dados brutos, push forçado, curl-pipe-bash |
| `lint_gate` | roda ruff/basedpyright/nextflow lint e devolve os erros ao modelo |
| `provenance` | caderno de laboratório automático em `runs/provenance.jsonl` |

### LSP

O Claude Code **tem** cliente LSP nativo (v2.0.74+), instalável como plugin:

```
/plugin install pyright-lsp@claude-plugins-official
```

Isso cobre Python. O que o LSP não cobre — Nextflow, Snakemake, schemas de
nf-core, R — fica com o hook `lint_gate`, que roda os linters e devolve os
diagnósticos pelo canal de feedback. As duas metades em
[`docs/LSP.md`](docs/LSP.md).

---

## As duas peças que rendem mais

**`/metodos`** lê `runs/provenance.jsonl` e escreve a seção de Métodos do seu
relatório a partir do que você *realmente* rodou — versões, parâmetros, datas
de acesso, inclusive as tentativas que falharam. É o registro que todo mundo
promete manter à mão e abandona na terceira semana.

**`lint_gate`** fecha o ciclo de correção. Sem ele, um erro de tipo num script
de análise só aparece depois de 40 minutos de pipeline. Com ele, volta na mesma
iteração.

---

## Testes

```bash
powershell -ExecutionPolicy Bypass -File tests/run_all.ps1
```

103 testes, nenhuma chamada de rede, nenhum crédito da NVIDIA consumido.

`tests/test_nemo.py` substitui as chamadas de embedding e reranking por um
dublê determinístico, então exercita indexação incremental, fatiamento e o
pipeline de busca inteiro sem gastar um token de crédito.

`tests/test_hooks.py` exercita os quatro hooks pelo contrato real do Claude Code
— JSON no stdin, código de saída como veredito. Cobre os 15 casos do
`guard_data` (inclusive os que ele **deve** deixar passar, como `nextflow run` e
`blastp`), a redação de chaves de API no log de proveniência, e a truncagem de
campos gigantes.

`tests/smoke_design.py` fica **fora** de `run_all.ps1`: e o unico teste que
chama os NIM de verdade e gasta credito. Rode de proposito, com
`uv run --script tests/smoke_design.py`, para confirmar que os endpoints de
desenho continuam de pe e que os formatos de resposta nao mudaram.

`tests/test_server.py` carrega o servidor MCP e verifica registro de ferramentas,
validação de sequência, mensagens de erro de credencial, e — o mais importante —
que coordenadas atômicas **nunca** voltam inline para o contexto.

Rode depois de qualquer alteração nos hooks ou no servidor.

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
| [`bootstrap/03_install_plugins.md`](bootstrap/03_install_plugins.md) | skills e plugins, com comandos |

---

## Estrutura

```
bio-stack/
├── CLAUDE.md                    roteamento e contexto do projeto
├── .mcp.json                    5 servidores MCP
├── .claude/
│   ├── settings.json            hooks, permissões, variáveis
│   ├── agents/                  6 agentes especializados
│   ├── commands/                5 slash commands
│   └── hooks/                   4 hooks Python
├── servers/
│   ├── nvidia_bio_mcp/          NIMs de biologia (código próprio)
│   └── nvidia_nemo_mcp/         NeMo Retriever (código próprio)
├── bootstrap/                   instalação em 3 camadas
├── tests/                       103 testes, sem rede e sem custo
├── editor/                      configuração do VS Code
├── docs/                        referência
├── data/{raw,proc}/             dados (raw é protegido por hook)
├── runs/                        saídas + provenance.jsonl
└── notebooks/                   análise exploratória
```
