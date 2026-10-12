# bio-stack

Ambiente de bioinformática e biologia sintética computacional para trabalho
acadêmico, com aceleração via NVIDIA BioNeMo NIM.

## Contexto do projeto

Trabalho de graduação/pós: análise de dados ômicos, biologia estrutural
computacional, desenho de construtos in silico, engenharia de pipelines.
Organismos e sistemas de laboratório de ensino e pesquisa básica.

O que **não** é o objeto deste projeto e deve ser encaminhado ao comitê de
biossegurança da instituição: trabalho com agentes de risco elevado, toxinas,
ou qualquer coisa que aumente virulência ou transmissibilidade. Esse tipo de
projeto tem via institucional própria e não passa por aqui.

## Hardware — restrição que molda todas as decisões

```
GPU:    NVIDIA GeForce RTX 5070 Laptop
VRAM:   8 GB (8151 MiB)
CUDA:   13.4 (driver 617.42)
Arch:   Blackwell, sm_120 (compute 12.0)  → PyTorch precisa de wheels cu128+
Torch:  2.11.0+cu128 · torch.version.cuda 12.8 · cuda.is_available() True
SO:     Windows 11 · WSL2 com Ubuntu-24.04 (20 núcleos, 15 GB RAM; disco: ver abaixo)
```

Disco: dentro do WSL o `df -h /` mostra ~1 TB (1007G) porque o `ext4.vhdx`
é um VHD dinâmico de 1 TiB virtual guardado em C:. O espaço que conta é o
livre em C: — 98 GB de 476 GB em 2026-10-12, dividido com o Windows. Meça
com `df -h /c` no Git Bash antes de baixar genomas.

A exigência de wheel cu128+ não é preferência: o driver 617.42 expõe CUDA 13.4,
mas não existe wheel oficial de PyTorch para sm_120 abaixo de cu128 — um wheel
cu121 instala sem reclamar e só falha no primeiro kernel. O `.venv` é montado
com `uv sync --python 3.11 --locked` (sem o `--python`, num clone novo o uv
escolhe o Python mais novo e o lock bifurca em 3.12): `gpu`, `dev`,
`evolucao`, `estatistica`, `pesquisa` e
`notebooks` são grupos ligados por `default-groups`, e o grupo `gpu` puxa
torch do índice `pytorch-cu128` declarado em `pyproject.toml`. `gpu` não é
extra — `uv sync --extra gpu` falha.

**Não proponha rodar localmente**: AlphaFold2, ESMFold, OpenFold2/3, Boltz-2,
Evo 2, RFdiffusion, Parabricks. Nenhum cabe em 8 GB. Use as ferramentas
`mcp__nvidia-bio__*`, que batem nos endpoints hospedados da NVIDIA.

**Cabe local**: ESM-2 650M, DNABERT-2, scanpy, pydeseq2, RDKit, biopython,
nvMolKit, inferência de modelos pequenos.

**Vive no WSL, já instalado e pronto** (Ubuntu 24.04.5 LTS, nome de distro
`Ubuntu-24.04`, reinstalada em 2026-10-10, systemd como PID 1): `nextflow`
26.04.6, `nf-core` 4.1.0, `bedtools` 2.31.1, `docker` 29.9.0, Java openjdk
21.0.12.1, e a **camada pixi** (2026-10-11): 37 ambientes de CLIs,
declaradas em `bootstrap/pixi-global.toml` — IQ-TREE 3.1.4, MAFFT, FAMSA,
MUSCLE 5, ClipKIT, trimAl, RAxML-NG, FastTree, ASTER/wASTRAL, HyPhy, PAML
(codeml), MACSE, PRANK, BLAST+ 2.17.0, DIAMOND, HMMER, MMseqs2, OrthoFinder 3.1.5, subread (featureCounts) 2.1.1, StringTie 3.0.3, miniprot,
minimap2, samtools/bcftools 1.24, seqkit 2.14.0, csvtk, taxonkit, NCBI
Datasets, EDirect, gffread, AGAT, utilitários UCSC (twoBitToFa, liftOver,
mafsInRegion, mafSpeciesSubset, bigWigAverageOverBed), jq, DuckDB, GNU
parallel, nf-test. Nada disso
precisa ser instalado e nada disso está no PATH do Windows — não procure lá
e não peça `wsl --install`.

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "samtools flagstat amostra.bam"
```

Onde cada binário mora decide em qual armadilha ele cai: as CLIs do pixi
estão em `~/.pixi/bin`, que vem antes de `/usr/bin` no PATH (por isso
`samtools` responde 1.24 e não o 1.19 do apt — em shell de login; fora
dele responde o 1.19, ver armadilha 3); bedtools veio do apt;
nextflow é binário solto em `~/.local/bin`; nf-core veio por pipx; docker-ce
veio do repositório oficial do Docker.

Pipeline de filogenia (IQ-TREE 3, ver `docs/FILOGENIA.md`), utilizável daqui
também — o invólucro é Python puro e chama o WSL sozinho:

```powershell
uv run python scripts/filogenia/filogenia.py --input <dir com 1 FASTA por loco> --outgroup <taxon>
```

O projeto aparece lá como `/mnt/c/Users/Manec/bio-stack`. Prefira trabalhar em
`~/` dentro do WSL quando o I/O importar — `/mnt/c` atravessa a fronteira de
sistema de arquivos e é ordens de grandeza mais lento para muitos arquivos
pequenos, que é exatamente o padrão de um pipeline Nextflow.

### Três armadilhas medidas ao chamar o WSL de fora, e elas se somam

Na prática as três defesas são necessárias ao mesmo tempo: script em arquivo,
`MSYS_NO_PATHCONV=1`, e `bash -l`.

1. **`wsl.exe` engole `$var`.** Chamado de fora (por `subprocess`, por exemplo),
   o `wsl.exe` perde as sequências `$var` na tradução da linha de comando — um
   laço `for t in ...; do echo "$t"; done` executa o número certo de vezes e
   imprime vazio em todas, com exit 0 e sem erro nenhum. Escreva comandos sem
   cifrão, ou ponha o script num arquivo e execute o arquivo. Ver
   `_wsl_presentes()` em `.claude/hooks/session_context.py` (o hook do
   Codex não sonda o WSL).
2. **Um shell MSYS reescreve o caminho antes do `wsl.exe` ver.** Passando
   `/mnt/c/Users/...` para `wsl.exe` de dentro de um Git Bash, o retorno é
   `bash: C:/Program Files/Git/mnt/c/Users/...: No such file or directory` — o
   MSYS traduziu o caminho POSIX para caminho Windows e colou a raiz da
   instalação do Git na frente. Contorno: `export MSYS_NO_PATHCONV=1` antes da
   chamada. Pelo PowerShell o problema não existe.
3. **Execute com `bash -l`, não `bash`.** O `-l` não é estilo: o PATH de
   `~/.local/bin` e dos shims do pipx está em `~/.profile`, lido só em shell de
   **login**. O `~/.bashrc` do Ubuntu não substitui — ele abre com
   `case $- in *i*) ;; *) return;; esac` e desiste em shell não-interativa, que
   é exatamente o caso de `bash arquivo.sh`. O efeito é uma falha **parcial**:
   nextflow, nf-core e seqkit dão `command not found` enquanto samtools,
   bcftools, bedtools e tabix seguem funcionando porque estão em `/usr/bin`
   — e samtools, bcftools e tabix respondem, sem aviso, o 1.19 do apt em vez
   do 1.24 do pixi, porque o `~/.pixi/bin` também só entra no login.
   Tem cara de instalação incompleta e manda investigar a coisa errada — é a
   mais traiçoeira das três.

### Docker é nativo na distro, não é Docker Desktop

`docker-ce` 29.9.0 roda **dentro** da Ubuntu-24.04, com daemon sob systemd
(`systemctl enable --now docker`), socket em `/var/run/docker.sock` e o usuário
`manec` no grupo `docker`. O Docker Desktop não está instalado no Windows e não
será. Consequência que evita um diagnóstico errado: **não existe `docker` no
PATH do Windows, e isso é o esperado, não uma lacuna** — toda chamada passa pelo
WSL.

Verificado de ponta a ponta, não só com `docker info`:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "cd ~/nf-smoke && nextflow run nf-core/demo -profile test,docker --outdir resultados"
```

Resultado: 8 processos em container (FASTQC ×3, SEQTK_TRIM ×3, COWPY, MULTIQC),
`Pipeline completed successfully`, relatório MultiQC de 2 452 732 bytes.
Então `-profile docker` é o padrão a usar aqui — não precisa de conda nem
singularity.

## Roteamento de agentes

| A pergunta é sobre… | Delegue para |
|---|---|
| estrutura 3D, docking, afinidade, pLDDT | `structure-bio` |
| desenho de backbone/binder, geração de molécula | `structure-bio` |
| matriz de contagem, DEG, scRNA-seq, enriquecimento | `omics-analyst` |
| plasmídeo, primer, códon, Golden Gate, CRISPR guide | `seq-designer` |
| variante, VCF, ACMG, gnomAD, farmacogenômica | `variant-curator` |
| árvore filogenética, alinhamento múltiplo, IQ-TREE, suporte de ramo, concordância, árvore de espécies | `phylo-analyst` |
| Nextflow, nf-core, Snakemake, cluster, container | `pipeline-engineer` |
| revisão de literatura, estado da arte, precedente | `lit-scout` |
| "onde nos meus PDFs fala sobre X" | `mcp__nvidia-nemo__library_search` |

Perguntas rápidas de lookup (o que é o gene X, qual o UniProt de Y) resolva
direto com `mcp__gget__*` ou `mcp__biocontext` — não delegue.

## Ordem de preferência das ferramentas

1. **Base de dados antes de predição.** Estrutura experimental no PDB >
   AlphaFold DB > predição nova. Anotação no UniProt > inferência.
   Predição custa tempo e crédito e é sempre menos confiável que medida.
2. **Endpoint hospedado antes de container local.** Nesta máquina não há
   escolha, mas mesmo com GPU maior o hospedado poupa dias de setup.
3. **nf-core antes de pipeline caseiro.**

## Economia de contexto

Arquivos biológicos são enormes. Regras:

- **Nunca** leia um `.cif`, `.pdb`, `.bam`, `.vcf` de genoma inteiro, ou FASTQ
  com Read. Um mmCIF típico passa de 100 mil tokens.
- Escreva um script curto que extrai o que interessa e leia a **saída do script**.
- As ferramentas `mcp__nvidia-bio__*` já gravam em `runs/nim/` e devolvem só
  o caminho — é de propósito.
- Para VCF e BAM, filtre com bcftools/samtools antes de olhar.

## Reprodutibilidade

Toda operação é registrada em `runs/provenance.jsonl` pelo hook
`.codex/hooks/provenance.py`, ligado em `PostToolUse` no `.codex/hooks.json`.
Isso existe para virar a seção de Métodos depois — rode `/metodos`.

Consequência prática: quando o usuário perguntar "o que eu rodei mesmo?",
consulte o log, não a memória da conversa.

## Padrões de código

- Python 3.11, `ruff` para lint e formato, type hints em funções públicas.
- Análise em script ou notebook versionado, nunca comando solto no terminal.
- Seeds fixadas em tudo que é estocástico — UMAP e Leiden mudam de resultado
  entre execuções sem isso. **Nunca use `np.random.seed()`**: é estado global de
  processo, e qualquer biblioteca importada depois pode consumir números do
  mesmo gerador e deslocar a sequência sem que nada falhe. O lint recusa
  (`ruff`, regra NPY002, ligada em `[tool.ruff.lint]` do `pyproject.toml`).
  Use `random_state=` no parâmetro da biblioteca
  (`sc.tl.leiden(..., random_state=0)`) e `rng = np.random.default_rng(0)` para
  aleatoriedade sua.
- Figuras com eixos rotulados, unidade, e n declarado.
- Caminhos relativos ao projeto; dados brutos em `data/raw/` ficam somente-leitura.

## Comandos

Espelhados em `.agents/skills/source-command-*`, usáveis daqui:

- `/bio-doctor` — diagnóstico do ambiente
- `/filogenia` — locos → IQ-TREE 3 (ML por loco, supermatriz, concordância, wASTRAL) → Métodos

Só do lado Claude (`.claude/commands/`, sem espelho aqui):

- `/nim` — estado dos endpoints NVIDIA
- `/fold` — predição de estrutura
- `/atavismo` — testa se um caráter foi perdido por gene ou por regulação
  (delega ao agente `evo-devo`, que também só existe em `.claude/agents/`)
- `/metodos` — gera a seção de Métodos a partir do log
- `/nvidia` — audita tudo que a chave NVIDIA destrava

## Estrutura

```
bootstrap/   instalação (rodar uma vez); 05_wsl_ferramentas.sh + pixi-global.toml = CLIs do WSL
servers/     servidor MCP da NVIDIA (código próprio)
pipelines/   pipelines Nextflow próprios (filogenia/ — ver docs/FILOGENIA.md)
scripts/     análises e invólucros (filogenia/filogenia.py roda o pipeline do Windows)
r/           R 4.6 + Bioconductor (renv): filogenia comparativa, DE de referência
.codex/      agents, hooks, config.toml — o lado não-Claude
.agents/     skills compartilhadas (fonte dos comandos espelhados)
.claude/     equivalentes do lado Claude: agents, commands, hooks, settings
data/raw/    dados brutos — protegidos por hook, nunca sobrescrever
data/proc/   dados processados
runs/        saídas de análise + provenance.jsonl
docs/        catálogo da stack, referência NVIDIA, revisões de literatura
notebooks/   análise exploratória
```
