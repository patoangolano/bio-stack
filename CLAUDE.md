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
CUDA:   13.5 (driver 621.01, medido em 2026-10-23; era 13.4 no 617.42)
Arch:   Blackwell, sm_120 (compute 12.0)  → PyTorch precisa de wheels cu128+
Torch:  2.11.0+cu128 · torch.version.cuda 12.8 · cuda.is_available() True
SO:     Windows 11 · WSL2 com Ubuntu-24.04 (20 núcleos, 15 GB RAM; disco: ver abaixo)
```

Disco: dentro do WSL o `df -h /` mostra ~1 TB (1007G) porque o `ext4.vhdx`
é um VHD dinâmico de 1 TiB virtual guardado em C:. O espaço que conta é o
livre em C: — 98 GB de 476 GB em 2026-10-12, dividido com o Windows. Meça
com `df -h /c` no Git Bash antes de baixar genomas.

O par `cu128` / `torch.version.cuda 12.8` é o que torna a exigência de wheels
cu128+ verificável: o driver 621.01 expõe CUDA 13.5, mas não existe wheel
oficial de PyTorch para sm_120 abaixo de cu128 — um wheel cu121 instala sem
reclamar e só falha no primeiro kernel.

**Não proponha rodar localmente**: AlphaFold2, ESMFold, OpenFold2/3, Boltz-2,
Evo 2, RFdiffusion, Parabricks. Nenhum cabe em 8 GB. Use as ferramentas
`mcp__nvidia-bio__*`, que batem nos endpoints hospedados da NVIDIA.

**Cabe local**: ESM-2 650M, DNABERT-2, scanpy, pydeseq2, RDKit, biopython,
nvMolKit, inferência de modelos pequenos. Desde 2026-10-12 também TotalSegmentator,
nnInteractive, DeepLabCut, SLEAP e BioEmu — todos em ambiente isolado com
`torch cu128` próprio, cada um verificado com kernel real em `sm_120`, não só
com `cuda.is_available()`.

**Vive no WSL** (Ubuntu 24.04.5 LTS, distro `Ubuntu-24.04`, reinstalada em
2026-10-10, 2026-10-16 e 2026-10-23): `nextflow` 26.04.6, `nf-core` 4.1.0,
`bedtools` 2.31.1, `docker` 29.9.2, Java 21 — a **camada base**, reconstruída
por `bootstrap/04b_wsl_base.sh` (como root) — e, desde 2026-10-11, a **camada
pixi**: 43 ambientes, 39 de CLIs (um por ferramenta, ou grupo coeso; `rclone` e
`pigz`, do backup, desde 2026-10-23) e 4 de biofísica (GROMACS, PyMOL, OpenMM, OpenSim), declarados em
`bootstrap/pixi-global.toml` (versão exata) com build e sha256 em
`bootstrap/pixi-global.snapshot.tsv`. Entre elas: `iqtree3` 3.1.4, `mafft`,
`famsa`, `muscle`, `clipkit`, `trimal`, `raxml-ng`, `FastTree`, `wastral`,
`hyphy` 2.5.103, `codeml` (PAML 4.10.10), `macse`, `prank`, BLAST+ 2.17.0,
`diamond`, `hmmer`, `mmseqs`, `orthofinder` 3.1.5, `featureCounts` (subread 2.1.1), `stringtie` 3.0.3, `miniprot`, `minimap2`, `samtools`/`bcftools`
1.24, `seqkit` 2.14.0, `csvtk`, `taxonkit`, `datasets`/`dataformat`
(NCBI), `efetch`/`esearch`, `gffread`, AGAT, utilitários UCSC (`twoBitToFa`,
`liftOver`, `mafsInRegion`, `mafSpeciesSubset`, `bigWigAverageOverBed`),
`jq`, `duckdb`, `parallel`, `nf-test`, `rclone`, `pigz`.
Não procure nenhuma delas no PATH do Windows. Reconstruir a distro inteira,
de uma vez e fora da árvore de processos da sessão (são horas, e job do Bash
tool morre de SIGHUP quando a sessão fecha):

```powershell
.\bootstrap\tarefa_longa.ps1 -Nome reconstruir-wsl -Script bootstrap\09_wsl_reconstruir.sh
```

O `09` faz, na ordem e sem nenhum passo que espere teclado, o mesmo que os
comandos abaixo, mais `wsl --set-default Ubuntu-24.04` e a conferência. Log e
código de saída em `runs/bootstrap/`. A mão, do Windows:

```bash
wsl --install -d Ubuntu-24.04 --no-launch --web-download
wsl -d Ubuntu-24.04 -u root -- bash /mnt/c/Users/Manec/bio-stack/bootstrap/04b_wsl_base.sh
wsl --terminate Ubuntu-24.04
wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/05_wsl_ferramentas.sh
wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/06_wsl_envs_biomec.sh
wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/07_wsl_claude_science.sh
```

Conferir: `uv run --script tests/check_wsl_ferramentas.py`. O `04b` cria o
usuário `manec` **sem senha** — o `-u root` do `wsl` dispensa senha, e o `sudo`
só passa a valer depois de `wsl -d Ubuntu-24.04 -u root passwd manec`. O `07`
deixa o daemon do Claude Science de pé na 8765, mas o login da conta é manual:
`wsl -d Ubuntu-24.04 -- bash -lc "claude-science url"` dá um link de uso único.
Nenhum script pede `sudo` interativo: o `06` travou 19 min num deles em
2026-10-16, esperando a senha de um usuário que não tinha senha.

O lado Windows tem o par equivalente. O `bootstrap/08_windows_uv.sh` refaz o
`.venv`, os `envs/` isolados e as seis uv tools, com kernel real em `sm_120`
em cada uma. O `bootstrap/restaura_drive.sh` traz de volta do backup no Google
Drive o que o git não guarda (`runs/`, `data/`, `~/rapinantes`, `~/nf-work`),
conferindo o md5 antes de extrair e sem sobrescrever nada que já exista. Os
dois rodam por `tarefa_longa.ps1`. Foi assim que a stack voltou em 2026-10-23,
depois de um `wsl --unregister` e de um re-clone do repositório feitos a mão.

Onde cada uma mora importa para a armadilha do PATH, logo abaixo: as do
pixi estão em `~/.pixi/bin`; bedtools veio do apt e está em `/usr/bin`;
nextflow é binário em `~/.local/bin`; nf-core veio por pipx. O apt **não**
instala samtools, bcftools nem tabix desde a reconstrução de 2026-10-16 — o
htslib 1.19 do Ubuntu tem CVEs abertas no decodificador CRAM, e quem responde
é o 1.24 do pixi. Isso vale **só em shell de login**: fora dele (`bash -c`,
cron, unidade systemd) o PATH não tem o `~/.pixi/bin` e `samtools` dá
`command not found`. Mais um motivo para a armadilha 3.

```bash
wsl -d Ubuntu-24.04 -- bash -lc "samtools flagstat amostra.bam"
```

O projeto aparece lá como `/mnt/c/Users/Manec/bio-stack`. Prefira trabalhar em
`~/` dentro do WSL quando o I/O importar — `/mnt/c` atravessa a fronteira de
sistema de arquivos e é ordens de grandeza mais lento para muitos arquivos
pequenos, que é exatamente o padrão de um pipeline Nextflow.

> **Três armadilhas medidas, e elas se somam.** Na prática você precisa das três
> defesas ao mesmo tempo: script em arquivo, `MSYS_NO_PATHCONV=1`, e `bash -l`.
>
> 1. **`wsl.exe` engole `$var`.** Chamado via `subprocess` do Python, o `wsl.exe`
>    perde as sequências `$var` na tradução da linha de comando — um laço
>    `for t in ...; do echo "$t"; done` executa o número certo de vezes e imprime
>    vazio em todas, com exit 0 e sem erro. Escreva comandos sem cifrão, ou ponha
>    o script num arquivo e execute o arquivo. Ver `_wsl_presentes()` em
>    `.claude/hooks/session_context.py`.
> 2. **Git Bash reescreve o caminho antes do `wsl.exe` ver.** Pelo Bash tool,
>    passar `/mnt/c/Users/...` para `wsl.exe` volta como
>    `bash: C:/Program Files/Git/mnt/c/Users/...: No such file or directory` — o
>    MSYS traduziu o caminho POSIX para caminho Windows e colou a raiz da
>    instalação do Git na frente. Contorno: `export MSYS_NO_PATHCONV=1` antes da
>    chamada.
> 3. **Execute com `bash -l`, não `bash`.** O `-l` não é detalhe de estilo: o
>    PATH de `~/.local/bin` e dos shims do pipx está em `~/.profile`, que só é
>    lido em shell de **login**. O `~/.bashrc` do Ubuntu não serve de substituto
>    — ele abre com `case $- in *i*) ;; *) return;; esac` e desiste em shell
>    não-interativa, que é exatamente o caso de `bash arquivo.sh`. O efeito é uma
>    falha **parcial**: nextflow, nf-core e todas as CLIs do pixi (samtools,
>    seqkit, tblastn…) dão `command not found` enquanto bedtools, Java e docker
>    seguem funcionando porque estão em `/usr/bin`. Isso tem cara de instalação
>    incompleta e manda o leitor investigar a coisa errada. É a mais traiçoeira
>    das três.

**Docker é nativo dentro da distro** (instalado em 2026-10-10, refeito em
2026-10-16 pelo `04b`): `docker-ce` 29.9.1 do repositório oficial do Docker,
daemon sob systemd
(`systemctl enable --now docker`), socket em `/var/run/docker.sock`, usuário
`manec` no grupo `docker`. Não há Docker Desktop no Windows e não haverá — logo
**não existe `docker` no PATH do Windows, e isso é o esperado, não lacuna**. Quem
procura o binário do lado Windows conclui errado que falta instalar. Ele entrou
uma vez, em 2026-10-21, e saiu em 2026-10-23 (`winget uninstall --id
Docker.DockerDesktop`): virou a distro padrão do WSL, e a WSL integration dele
disputa o `/var/run/docker.sock` com o docker-ce da distro.

Verificado com o pipeline inteiro, não só com `docker info`: o `info` prova
apenas que o daemon responde, enquanto o pipeline exercita pull de imagem, bind
mount do diretório de trabalho e propagação de código de saída — é nesses três
que a configuração costuma falhar depois de um `info` verde.

```bash
wsl -d Ubuntu-24.04 -- bash -lc "cd ~/nf-smoke && nextflow run nf-core/demo -profile test,docker --outdir resultados"
```

Resultado: 8 processos em container (FASTQC ×3, SEQTK_TRIM ×3, COWPY, MULTIQC),
`Pipeline completed successfully`, relatório MultiQC de 2 452 732 bytes
(2,45 MB) gerado. Então `-profile docker` é o padrão a usar aqui — não precisa
de conda nem singularity.

## Roteamento de agentes

| A pergunta é sobre… | Delegue para |
|---|---|
| estrutura 3D, docking, afinidade, pLDDT | `structure-bio` |
| desenho de backbone/binder, geração de molécula | `structure-bio` |
| matriz de contagem, DEG, scRNA-seq, enriquecimento | `omics-analyst` |
| plasmídeo, primer, códon, Golden Gate, CRISPR guide | `seq-designer` |
| variante, VCF, ACMG, gnomAD, farmacogenômica | `variant-curator` |
| perda gênica, pseudogene, estado ancestral, atavismo, evo-devo | `evo-devo` |
| árvore filogenética, alinhamento múltiplo, IQ-TREE, suporte de ramo, concordância, árvore de espécies | `phylo-analyst` |
| Nextflow, nf-core, Snakemake, cluster, container | `pipeline-engineer` |
| revisão de literatura, estado da arte, precedente | `lit-scout` |
| "onde nos meus PDFs fala sobre X" | `mcp__nvidia-nemo__library_search` |
| varredura longa que não deve ocupar este contexto | `mcp__claude-science__science_ask` |

O `claude-science` opera o daemon local pela API REST dele — ver
`docs/ECOSSISTEMA.md`. Mande para lá o que é longo ou precisa dos 24 conectores
nativos de biologia do Claude Science; pergunta trivial se responde aqui. O uso
sai da assinatura Max, e o servidor **recusa** se a conta cair para cobrança por
token. Peça saída enxuta no texto do pedido ("devolva só a tabela", "não imprima
sequência") — a resposta volta para este contexto.

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

Toda operação é registrada em `runs/provenance.jsonl` pelo hook `provenance.py`.
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
  (`ruff` regra NPY002). Use `random_state=` no parâmetro da biblioteca
  (`sc.tl.leiden(..., random_state=0)`) e `rng = np.random.default_rng(0)` para
  aleatoriedade sua.
- Figuras com eixos rotulados, unidade, e n declarado.
- Caminhos relativos ao projeto; dados brutos em `data/raw/` ficam somente-leitura.

## Comandos

- `/bio-doctor` — diagnóstico do ambiente
- `/nim` — estado dos endpoints NVIDIA
- `/fold` — predição de estrutura
- `/atavismo` — testa se um caráter foi perdido por gene ou por regulação
- `/filogenia` — locos → IQ-TREE 3 (ML por loco, supermatriz, concordância, wASTRAL) → Métodos
- `/metodos` — gera a seção de Métodos a partir do log
- `/nvidia` — audita tudo que a chave NVIDIA destrava

## Biofísica, biomecânica, morfologia e osteologia

Camada de 2026-10-12. O catálogo completo, com as armadilhas medidas, está em
[`docs/STACK.md`](docs/STACK.md) — leia antes de mexer. O resumo operacional:

**Onde cada coisa vive, e por quê.** A regra única é que nada negocie `numpy`
nem `scipy` com o venv do projeto, onde as duas estão travadas em versão exata
por peso estatístico.

| Camada | O que está lá | Como chamar |
|---|---|---|
| venv do projeto (Windows) | MDAnalysis, MuJoCo, febio-python, trimesh, PyVista, Open3D, gmsh, morphosource | `uv run python` |
| `uv tool` (isolado) | TotalSegmentator, nnInteractive, DeepLabCut, SLEAP | `TotalSegmentator`, `dlc`, `sleap-train`, … |
| `envs/<nome>/.venv` | BioEmu, BoneJ_Headless | `envs/bioemu/.venv/Scripts/python.exe` |
| pixi global (WSL) | GROMACS, PyMOL, OpenMM, OpenSim 4.6 | `gmx`, `pymol`, `python-openmm`, `python-opensim` |
| pixi em `~/envs` (WSL) | MyoConverter (OpenSim 4.4.1), opencap-processing (OpenSim 4.5) | `pixi run --manifest-path ~/envs/<nome>/pixi.toml python` |
| aplicativo gráfico | 3D Slicer 5.12.4 + SlicerMorph, Fiji + BoneJ 7.2.2, Blender 5.2.1 | MCP `slicer`, `blender`, `chimerax` |

**Há três OpenSim na máquina e isso é de propósito.** 4.6 no pixi global
(com Moco dentro), 4.4.1 no `envs/myoconverter`, 4.5 build `py311np123` no
`envs/opencap`. Cada upstream exige a sua, e as versões de numpy que elas
arrastam são incompatíveis entre si. Não unifique.

**Os MCP de ferramenta gráfica sobem com o programa fechado.** `slicer`,
`blender` e `chimerax` respondem ao handshake e listam ferramentas mesmo sem o
programa aberto; só a chamada falha. Um `mcp list` verde não quer dizer que a
ponte está aberta — no Slicer é preciso iniciar o módulo Web Server, no Blender
clicar Start MCP Server.

**ChimeraX 1.12 e FEBio 4.13 entraram em 2026-10-12**, em
`Program Files/ChimeraX 1.12` e `Program Files/FEBioStudio`. Duas consequências
que mudam como se usa:

- **Não rode `remotecontrol rest start` na mão.** O `chimerax-mcp` tem
  `auto_launch=True`: ele acha o executável e sobe o ChimeraX com o REST na
  8765 na primeira chamada. Verificado abrindo a lisozima 1AKI do PDB.
- **`febio-python` 0.2.1 não lê o `.xplt` do FEBio 4.13** (`Incorrect XPLIT
  file version: 53`). O `.feb` e o `run()` funcionam; para número em Python,
  declare `<Output><logfile>` no `.feb` e leia o texto delimitado com pandas. O
  plot binário serve só para a GUI do FEBio Studio. Detalhe e exemplo em
  `docs/STACK.md`.
- `C:/Program Files/FEBioStudio/bin` está no PATH do usuário — sem isso o
  `febio_python.feb.run.run()` não acha o `febio4.exe`.

**Reaplicar:** `uv sync` (Windows) · `Rscript r/instalar.R` (precisa de
Rtools45) · `bash bootstrap/05_wsl_ferramentas.sh` e
`bash bootstrap/06_wsl_envs_biomec.sh` (no WSL; o `05` traz GROMACS, PyMOL,
OpenMM e OpenSim desde 2026-10-16) · `bash vendor/clonar.sh`.

## Estrutura

```
bootstrap/   instalação (rodar uma vez)
servers/     servidor MCP da NVIDIA (código próprio)
pipelines/   pipelines Nextflow próprios (filogenia/ — ver docs/FILOGENIA.md)
scripts/     análises e invólucros (filogenia/filogenia.py roda o pipeline do Windows)
r/           R 4.6 + Bioconductor (renv): filogenia comparativa, DE de referência
envs/        ambientes isolados: requirements.in/pixi.toml versionados, .venv fora do git
vendor/      clones de terceiros que não são pacote (reconstrói: vendor/clonar.sh)
.claude/     agents, commands, hooks, settings
data/raw/    dados brutos — protegidos por hook, nunca sobrescrever
data/proc/   dados processados
runs/        saídas de análise + provenance.jsonl
docs/        catálogo da stack, referência NVIDIA, revisões de literatura
notebooks/   análise exploratória
```
