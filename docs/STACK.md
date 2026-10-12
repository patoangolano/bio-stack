# Catálogo do ecossistema

Levantamento feito em setembro de 2026. Marcado com **★** o que está instalado
ou configurado nesta stack; o resto é o mapa para expandir.

---

## Servidores MCP

### Instalados

| Servidor | Fonte | Cobre |
|---|---|---|
| ★ `nvidia-bio` | **código deste projeto** (`servers/`) | NIM de biologia da NVIDIA com polling assíncrono e gravação em disco |
| ★ `biocontext` | [biocontext-ai/knowledgebase-mcp](https://github.com/biocontext-ai/knowledgebase-mcp) | 14 bases: STRING, Open Targets, Reactome, UniProt, HPA, KEGG, AlphaFold DB, Ensembl, InterPro, EuropePMC, PRIDE, PanglaoDB, bioRxiv, Antibody Registry — **exige `--with "bibtexparser<2"`**, ver nota abaixo |
| ★ `gget` | [longevity-genie/gget-mcp](https://github.com/longevity-genie/gget-mcp) | gget do Pachter Lab: busca de gene, BLAST, BLAT, MUSCLE, PDB, AlphaFold, Enrichr, ARCHS4, COSMIC, CELLxGENE |
| ★ `biomcp` | [genomoncology/biomcp](https://github.com/genomoncology/biomcp) | ~30 fontes biomédicas: ensaios clínicos, ClinVar, MyVariant, OncoKB, openFDA, PubMed |
| ★ `tooluniverse` | [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse) | catálogo de ferramentas científicas do MIMS (Harvard) |
| ★ `claude-science` | **código deste projeto** (`servers/`) | opera o daemon local do Claude Science pela API REST dele |
| ★ `slicer` | [zhaoyouj/mcp-slicer](https://pypi.org/project/mcp-slicer/) | 3D Slicer: lista nós da cena, executa Python dentro do Slicer, captura tela. **Exige CPython 3.13** — por isso `uvx --python 3.13` |
| ★ `blender` | [ahujasid/mcp-for-blender](https://github.com/ahujasid/mcp-for-blender) | Blender: executa Python na cena, limpeza e retopologia de malha. 9 ferramentas |
| ★ `chimerax` | [chimerax-mcp](https://pypi.org/project/chimerax-mcp/) | ChimeraX pela REST API: abre/salva estrutura, sequência, script. 39 ferramentas |
| ★ `totalsegmentator` | vem dentro do próprio TotalSegmentator (`totalseg_mcp`) | segmentação de TC: classes, tarefas, modalidade, fase de contraste. 8 ferramentas |

> **Por que o `bibtexparser<2` no `biocontext`.** O pacote importa `scholarly`,
> que faz `from bibtexparser.bibdatabase import BibDatabase`. O bibtexparser 2.x
> removeu esse módulo, então o `uvx` resolvia a versão nova e o servidor morria
> no import com `ModuleNotFoundError` — aparecia no Claude Code apenas como
> `CONNECTION_CLOSED`, sem pista da causa. Diagnosticado em 2026-09-12 rodando
> o handshake stdio na mão. Com o pin, o servidor sobe (`BioContextAI 2.13.1`).

> **Os três MCP de ferramenta gráfica sobem sem o programa aberto.** `slicer`,
> `blender` e `chimerax` são pontes: o processo MCP inicia e responde ao
> handshake mesmo com o programa fechado, e só a *chamada de ferramenta* falha.
> Medido em 2026-10-12 com handshake stdio na mão — os quatro responderam
> `tools/list` (3, 9, 39 e 8 ferramentas). Consequência prática: um `mcp list`
> verde não significa que a ponte do outro lado está aberta. Para cada um:
> no Slicer, abrir o módulo **Web Server** e iniciar o servidor; no Blender, o
> addon já está habilitado e persistido (`blender_mcp.py` em
> `AppData/Roaming/Blender Foundation/Blender/5.2/scripts/addons`), falta
> clicar **Start MCP Server**; no ChimeraX, o programa precisa estar instalado
> — ver a seção de biofísica abaixo.

> **`blender-mcp` virou `mcp-for-blender`.** O nome antigo ainda funciona, como
> invólucro de compatibilidade que imprime aviso de depreciação a cada início.
> O `.mcp.json` aqui usa o nome novo.

### Já conectados nesta conta Claude

PubMed · bioRxiv/medRxiv · BioRender · Synthesize Bio (GEM-1) · Hugging Face Hub

### Vale adicionar conforme a necessidade

| Servidor | Fonte | Quando |
|---|---|---|
| Seqera MCP | [docs.seqera.io](https://docs.seqera.io/platform-cloud/seqera-mcp/overview) | trabalha com Nextflow/nf-core a sério |
| `ChatSpatial` | [cafferychen777/ChatSpatial](https://github.com/cafferychen777/ChatSpatial) | transcriptômica espacial |
| `plant-genomics-mcp` | [musharna/plant-genomics-mcp](https://github.com/musharna/plant-genomics-mcp) | genomas de plantas (TAIR, 12 espécies) |
| `roda-mcp` | [awslabs/mcp](https://github.com/awslabs/mcp/tree/main/src/roda-mcp-server) | datasets públicos na AWS (1.100+) |
| `Helix GenoSphere` | registro de conectores | genômica humana + dados clínicos longitudinais |
| `Boltz API` | registro de conectores | alternativa hospedada ao Boltz-2 |
| `LatchBio` | registro de conectores | lançar workflows na nuvem |
| `Synapse.org` | registro de conectores | dados de consórcios (Sage Bionetworks) |

### Ferramentas de linha de comando — ★ camada pixi no WSL (2026-10-11)

37 ambientes pixi (um por ferramenta, ou grupo coeso), versão exata em
[`bootstrap/pixi-global.toml`](../bootstrap/pixi-global.toml), build e
sha256 em `bootstrap/pixi-global.snapshot.tsv`, conferidas por
`tests/check_wsl_ferramentas.py`. Chame direto pelo WSL
(`wsl -d Ubuntu-24.04 -- bash -lc "blastn -version"`).

| Família | Ferramentas |
|---|---|
| filogenia | ★ IQ-TREE 3.1.4, MAFFT 7.525, FAMSA 2.4.1, MUSCLE 5.3, ClipKIT 2.14.0, trimAl 1.5.1, RAxML-NG 2.0.3, FastTree 2.2.0, ASTER 1.25 (wASTRAL), gotree, goalign |
| evolução molecular | ★ HyPhy 2.5.103, PAML 4.10.10 (codeml, MCMCtree), MACSE 2.07, PRANK |
| homologia e ortologia | ★ BLAST+ 2.17.0, DIAMOND 2.2.8, HMMER 3.4, MMseqs2 18, OrthoFinder 3.1.5, miniprot 0.18, minimap2 2.31 |
| transcrição | ★ subread 2.1.1 (featureCounts), StringTie 3.0.3 — leitura de RNA-seq sobre janela não anotada (relíquia transcrita ou não) |
| genomas e dados | ★ samtools/bcftools 1.24, seqkit 2.14.0, csvtk 0.38.0, taxonkit 0.20.0 (base do NCBI baixada pelo bootstrap, data em `~/.taxonkit/BAIXADO_EM`), NCBI Datasets 18.38.0, EDirect 26.0, gffread 0.12.9, AGAT 1.7.0, UCSC kent v482 (twoBitToFa, faToTwoBit, liftOver, mafsInRegion, mafSpeciesSubset, bigWigAverageOverBed, bigWigToWig), jq 1.8.2, DuckDB 1.5.6, GNU parallel |
| workflow | ★ nf-test 0.9.5 |

> **Não use a organização `bio-mcp`** (`bio-mcp-blast`, `bio-mcp-samtools`…)
> que este catálogo sugeria antes: sem release, e uma camada MCP por cima de
> CLI que o modelo já chama direto pelo WSL só acrescenta ponto de falha.
> Verificado na pesquisa de 2026-10-11.

O levantamento completo — o que entrou, o que ficou como receita e o que
foi descartado, com fonte — está em
[`PESQUISA-STACK-2026-10.md`](PESQUISA-STACK-2026-10.md).

### R — ★ biblioteca do projeto em `r/` (renv, 2026-10-11)

R 4.6.1 + Bioconductor 3.23, versões travadas em `r/renv.lock` (237
pacotes). Restaurar: `Rscript r/instalar.R`. Detalhes em
[`r/README.md`](../r/README.md).

| Para quê | Pacotes |
|---|---|
| filogenia comparativa, ASR de caráter | ape, phangorn, phytools, castor, geiger, corHMM, TreeDist, phylter |
| figura de árvore | ggtree, treeio, tidytree, ggtreeExtra |
| DE de referência (confere o pydeseq2) | DESeq2, edgeR, limma, sva, apeglm, ashr, IHW, qvalue, variancePartition (`dream`), tximport |
| ★ morfometria geométrica e fósseis (2026-10-12) | geomorph 4.1.1, paleobioDB 1.0.1, RERconverge 0.3.0 |

> **Rtools45 passou a ser dependência do ambiente R.** `RERconverge` compila
> C++ e a instalação morria em `'make' not found`; o renv faz rollback
> atômico, então as 42 bibliotecas da leva caíam junto — parecia erro do
> geomorph. Resolvido com `winget install RProject.Rtools` (4.5.6768). A
> página do Rtools45 diz que ele serve "R 4.5 and later, incluindo o R-devel
> que vai ser 4.7.0", então é o toolchain correto para o R 4.6.1 desta
> máquina; o CRAN não publica um `rtools46`. O instalador põe
> `C:\rtools45\usr\bin` no PATH e o R acha o `make` sem configuração.

### Registros para consultar

- [BioContextAI Registry](https://github.com/biocontext-ai/registry/) — catálogo curado de MCPs biomédicos
- [MCPmed](https://github.com/MCPmed) — implementações de referência (GEO, STRING, UCSC Cell Browser, PLSDB)
- [awesome-mcp-servers #bio](https://github.com/punkpeye/awesome-mcp-servers#bio)
- [BioinfoMCP](https://github.com/florensiawidjaja/BioinfoMCP) — gera MCP automaticamente a partir da documentação de uma ferramenta CLI

---

---

## ★ Biofísica, biomecânica, morfologia e osteologia (2026-10-12)

Camada instalada a partir de `Downloads/ferramentas-bio-claude-code.md`. A
divisão por ambiente não é arbitrária: ela existe para que **nada** negocie
`numpy` e `scipy` com o venv do projeto, onde as duas estão travadas em versão
exata porque sustentam afirmação de Métodos (ver `pyproject.toml`).

### Grupos do `pyproject.toml` — venv do projeto, Windows

Entram num `uv sync` pelado (estão em `default-groups`). Todos com wheel
`win_amd64` ou puro-Python para o CPython 3.11.17.

| Grupo | Pacotes | Para quê |
|---|---|---|
| `biofisica` | MDAnalysis 2.10.0 | ler trajetória de dinâmica molecular |
| `biomecanica` | MuJoCo 3.15.0, febio-python 0.2.1 | motor de física; escrever e ler o `.feb` do FEBio |
| `morfologia` | trimesh 5.1.1, PyVista 0.49.1, Open3D 0.20.0, gmsh 4.15.2, morphosource 1.2.0 | limpar/medir/visualizar malha, gerar malha de EF, buscar mídia no MorphoSource |

A resolução foi conferida antes de sincronizar: `numpy==2.4.6`,
`scipy==1.17.1`, `statsmodels==0.15.0`, `scanpy==1.11.5`, `pydeseq2==0.5.4` e
`torch==2.11.0+cu128` ficaram onde estavam. O wheel do `gmsh` traz o SDK, então
o binário `gmsh` aparece no PATH da venv sem instalação separada.

### Camada pixi no WSL — o que é conda-only ou CLI nativa

| Ambiente | Versão | Exposto como |
|---|---|---|
| `gromacs` | 2026.3 | `gmx` |
| `pymol` | pymol-open-source 3.1.0 | `pymol` |
| `openmm` | OpenMM 8.6.1 + python 3.11.14 | `python-openmm` |
| `opensim` | OpenSim 4.6 + python 3.11 | `python-opensim` |

OpenMM enxerga a GPU: `Platform` lista `['Reference', 'CPU', 'CUDA']`. O
**Moco já vem dentro do OpenSim 4.6** — `MocoStudy`, `MocoTrack` e
`MocoInverse` respondem no ambiente `opensim`, então não há nada separado a
instalar para controle ótimo.

OpenMM e OpenSim são bibliotecas, não CLIs, e por isso o ambiente expõe um
`python` próprio em vez de um binário de ferramenta. É deliberado: script que
faz `import openmm` precisa do módulo no **seu** interpretador.

### Ferramentas isoladas (`uv tool`) — cada uma com sua própria árvore de torch

| Ferramenta | Versão | torch | Executáveis |
|---|---|---|---|
| TotalSegmentator | 2.18.0 | 2.11.0+cu128 | `TotalSegmentator`, `totalseg_*` (16), `totalseg_mcp` |
| nnInteractive (servidor) | 0.3.2 | 2.8.0+cu128 | `nninteractive-slicer-server` |
| DeepLabCut | 3.0.2 | 2.11.0+cu128 | `dlc` |
| SLEAP | 1.6.5 + sleap-nn 0.3.3 | 2.11.0+cu128 | `sleap`, `sleap-train`, `sleap-track`, … (8) |
| chimerax-mcp | 0.1.1 | — | `chimerax-mcp` |

> **`cuda.is_available() == True` não prova nada; kernel prova.** O CLAUDE.md
> já registra que um wheel errado instala calado e só falha no primeiro kernel.
> Cada um dos quatro ambientes acima rodou um `matmul` real em `device='cuda'`
> e sincronizou — `sm_120`, RTX 5070 Laptop. Só depois disso foram declarados
> instalados.

> **Três armadilhas medidas na instalação dessas ferramentas.**
>
> 1. **`--index` do cu128 não basta.** `uv tool install totalsegmentator
>    --index https://download.pytorch.org/whl/cu128 --index-strategy
>    unsafe-best-match` instalou `torch==2.14.1+cpu`: aquele índice também
>    serve build de CPU, e ele ganhou. O que resolve é `--torch-backend cu128`.
> 2. **Um extra chamado `cuda128` entregou CPU.** O SLEAP declara os extras
>    `nn-cuda128`, `nn-cuda118`, `nn-cuda130`. Instalar `sleap[nn-cuda128]`
>    trouxe `torch==2.14.1+cpu` — o extra só marca a intenção e não carrega
>    índice. `sleap[nn-cuda128] --torch-backend cu128` é que traz o
>    `2.11.0+cu128`. É a armadilha do CLAUDE.md numa forma nova: o nome promete
>    CUDA e o `AssertionError: Torch not compiled with CUDA enabled` só aparece
>    no primeiro kernel.
> 3. **`totalseg_mcp` vem quebrado de fábrica.** O pacote não declara
>    `uvicorn` nem `mcp`, e o executável morre em `ModuleNotFoundError`.
>    Corrigido com `--with uvicorn --with mcp` na instalação da ferramenta.

### Venvs dedicadas (`envs/<nome>/`) — pins que brigam com o projeto

O `requirements.in` de cada uma é versionado e explica o porquê; o
`requirements.txt` ao lado é o `uv pip freeze` do que caiu. A `.venv` fica fora
do git.

| Ambiente | Conteúdo | Por que não é um grupo do `pyproject.toml` |
|---|---|---|
| `envs/bioemu` | BioEmu 1.4.1, torch 2.11.0+cu128 (103 pacotes) | traz a própria árvore de torch; e `uv tool install` recusa, porque o pacote não declara console script (`No executables are provided`) — o uso é `python -m bioemu.sample` |
| `envs/bonej-headless` | pynrrd, pydicom, nibabel, SimpleITK 2.5.6, tifffile, matplotlib | pilota o Fiji; o `requirements.txt` do upstream pede `numpy==1.26.0` e o projeto trava `2.4.6` |

BioEmu prevê ensemble de equilíbrio **em representação de backbone apenas** —
não devolve cadeia lateral. Licença MIT.

> **Os pins do BoneJ_Headless são de 2023 e não instalam.** `matplotlib==3.4.3`
> não tem wheel para CPython 3.11, cai para compilação e para em "Unable to
> find a compatible Visual Studio installation". A venv usa piso de versão, e o
> desvio está escrito em `envs/bonej-headless/requirements.in` — nenhuma
> estatística depende dessas versões, diferente do que vale no `pyproject.toml`.

### Ambientes pixi de biomecânica no WSL (`envs/*/pixi.toml`)

Manifesto versionado no repo; o ambiente é construído em `~/envs/<nome>` dentro
do WSL, não em `/mnt/c` — um ambiente pixi são dezenas de milhares de arquivos
pequenos, que é o padrão de acesso em que `/mnt/c` mais sofre.

Aplicar: `bash bootstrap/06_wsl_envs_biomec.sh` (no WSL).

| Ambiente | Conteúdo |
|---|---|
| `envs/myoconverter` | python 3.10, OpenSim 4.4.1 (`py310np121`), myoconverter do git |
| `envs/opencap` | python 3.11, OpenSim 4.5 (`py311np123`), scipy 1.10.0, casadi |

> **Os dois upstreams exigem versão exata do OpenSim, e cada um exige uma
> diferente.** O MyoConverter precisa de 4.4.1 com numpy 1.21; o
> opencap-processing manda `opensim=4.5=py311np123` no README, e é esse build
> que faz o `scipy==1.10.0` do `requirements.txt` dele resolver — com o 4.6,
> que traz numpy 2.4.6, a resolução morre em "Because scipy==1.10.0 depends on
> numpy>=1.19.5,<1.27.0 and numpy==2.4.6". Nenhum dos dois cabe no ambiente
> `opensim` do pixi global, que tem 4.6. Daí três OpenSim na máquina, cada um
> com seu dono.

> **O MyoConverter se contradiz consigo mesmo.** O `conda_env.yml` dele pede
> `pyvista=0.41.1`; o `pyproject.toml` do próprio pacote exige
> `pyvista>=0.40.0,<0.41.0`. Quem manda na instalação é o `pyproject`, e o
> conda-forge não empacota 0.40.x — então não há pin conda que sirva. O
> manifesto aqui é mínimo de propósito: o conda entrega só o interpretador e o
> OpenSim, e o resto vem do PyPI com as faixas do próprio pacote. **Aviso do
> upstream, que muda o que se pode afirmar: modelo convertido pode exigir
> ajuste antes de servir para pesquisa.**

> **O `evdev` do opencap custou três conclusões erradas antes da certa.** O
> `requirements.txt` do upstream traz `maskpass==0.3.6`, que depende de
> `evdev`, que compila extensão C. A sequência, toda medida em 2026-10-12:
>
> 1. O erro manda `apt-get install linux-headers-$(uname -r)`. Não funciona: o
>    kernel é `6.18.33.2-microsoft-standard-WSL2` e a Ubuntu não empacota
>    headers para ele. `linux-libc-dev` (6.8.0) passa dessa etapa.
> 2. Em seguida: `error: [Errno 2] No such file or directory: 'gcc'`. O build
>    roda com o PATH do ambiente pixi, onde o gcc do sistema não aparece.
>    Resolvido com o metapacote `c-compiler` **dentro** do ambiente.
> 3. Em seguida: `error: 'BUS_AMD_SFH' undeclared`, mais
>    `KEY_ALL_APPLICATIONS`, `KEY_LINK_PHONE` e `KEY_REFRESH_RATE_TOGGLE`. O
>    `ecodes.c` pré-gerado do evdev referencia constantes de `linux/input.h`
>    posteriores ao 6.8. Fixar `evdev<2` **não resolve** — o 1.9.3 falha nas
>    mesmas constantes, então não existe pin de evdev que sirva.
> 4. O que resolve é `kernel-headers_linux-64 >=6.12` do conda-forge, no
>    próprio ambiente: o gcc do `c-compiler` lê o sysroot do ambiente, não o
>    `/usr/include` da distro.
>
> Moral para o próximo ambiente deste tipo: toolchain e header pertencem ao
> **manifesto**, não ao `apt`. Só assim o ambiente reproduz em outra máquina. O
> `linux-libc-dev` continua instalado e o passo segue no
> `bootstrap/06_wsl_envs_biomec.sh`, mas ele não é o que conserta o evdev.

### Aplicativos gráficos

| Programa | Versão | Onde | Como entrou |
|---|---|---|---|
| 3D Slicer | 5.12.4 (rev 34645) | `AppData/Local/slicer.org/3D Slicer 5.12.4` | download direto, sha512 conferido, instalador NSIS em `/S` |
| SlicerMorph | git `0934a26` | extensão do Slicer | pacote do `slicer-packages.kitware.com`, sha512 conferido |
| Fiji + BoneJ | Fiji stable; BoneJ 7.2.2 + bonej-plus 0.0.24 | `C:/Users/Manec/Apps/Fiji.app` | zip portátil; update site adicionada por CLI |
| Blender | 5.2.1 LTS | `Program Files/Blender Foundation/Blender 5.2` | já estava instalado |
| ChimeraX | 1.12 (2026-06-12) | `Program Files/ChimeraX 1.12` | instalado por Мухамед em 2026-10-12, depois do aceite de licença |
| FEBio Studio + solver | Studio 3.2, **FEBio 4.13.0** | `Program Files/FEBioStudio` | instalado por Мухамед em 2026-10-12, depois do cadastro |

O SlicerMorph veio com as dependências que ele declara —
`SegmentEditorExtraEffects`, `SurfaceMarkup` e `Sandbox`. Dos 197 módulos que o
Slicer carrega, respondem **ALPACA** (landmarking automático), **GPA** (análise
de Procrustes generalizada), `MarkupEditor`, `PlaceLandmarkGrid`,
`SlicerMorphSampleData` e **`MorphoSourceImport`** — importação direta do
MorphoSource dentro do Slicer. *MALPACA* não é módulo separado na versão atual:
é o modo multi-template dentro do ALPACA.

> **Instalar extensão do Slicer pelo gestor, em modo headless, pendura para
> sempre.** `downloadAndInstallExtensionByName` é assíncrono e precisa do laço
> de eventos do Qt, que `--no-main-window --python-script` não roda. Duas
> tentativas ficaram presas em 12–14 MB de RSS, sem nenhuma saída. O que
> funciona é baixar o pacote e chamar `installExtension(arquivo)`, que é
> síncrono. O endpoint do pacote é
> `https://slicer-packages.kitware.com/api/v1/item/<extension_id>/download` — o
> `extensions.slicer.org/api/v1/...` devolve o HTML da SPA (1061 bytes
> idênticos para qualquer id, que é como o erro se apresenta). O
> `extension_id` e o `meta.sha512` de cada extensão estão em
> `.../Extensions-34645/ExtensionsMetadataFromServer.json`, que o próprio
> Slicer grava ao atualizar os metadados.

> **O Slicer no Windows não despeja `stdout` antes de sair, e `slicer.util.exit`
> corta o buffer.** Script headless que precise reportar resultado tem de
> escrever em arquivo; `print` se perde. Foi isso que fez as duas primeiras
> tentativas parecerem silenciosas em vez de penduradas.

> **O Fiji do winget está quebrado.** `Fiji.ImageJ` aponta para
> `fiji-latest-portable-nojava.zip` no arquivo morto e devolve **403**. O zip
> atual é `fiji-stable-win64-jdk.zip` (521 754 144 bytes), e o nome tem
> `win64`, não `windows`. Sendo portátil, instalar é descompactar. A update
> site do BoneJ entra sem GUI:
> `ImageJ-win64.exe --update add-update-site BoneJ https://sites.imagej.net/BoneJ/`
> seguido de `--update update`.

> **O Blender 5.2.1 já estava na máquina e o `winget upgrade` falha com 1602.**
> O desinstalador silencioso da versão antiga recusa. É inofensivo: 5.2.1 está
> muito acima do 3.0+ que a ponte MCP pede. O addon `blender_mcp.py` foi
> instalado com `uvx mcp-for-blender install-addon` e habilitado e persistido
> com `addon_enable` + `save_userpref`. A pasta de addons do usuário não
> existia porque o Blender nunca havia rodado — foi criada pelo próprio
> Blender, com `bpy.utils.user_resource('SCRIPTS', path='addons', create=True)`.

> **O `chimerax-mcp` sobe o ChimeraX sozinho — não precisa de
> `remotecontrol rest start` na mão.** O cliente tem `auto_launch=True` por
> padrão: na primeira chamada que não encontra REST em `127.0.0.1:8765`, ele
> acha o executável e lança `ChimeraX --cmd "remotecontrol rest start port
> 8765"`. A descoberta varre `%ProgramFiles%` por pasta começando com
> "chimerax" e testa `<pasta>/bin/ChimeraX.exe` — encontra
> `C:/Program Files/ChimeraX 1.12/bin/ChimeraX.exe` sem configuração.
> `CHIMERAX_BIN` sobrepõe, se um dia a instalação sair do lugar. Verificado em
> 2026-10-12 com ida e volta real: abriu a lisozima **1AKI** do PDB,
> `info models` devolveu `model id #1 type AtomicStructure name 1aki` e a
> enumeração de átomos voltou por REST.

> **`febio-python` 0.2.1 NÃO lê o `.xplt` do FEBio 4.13, e isso muda o fluxo de
> pós-processamento.** O parser recusa com `ValueError: Incorrect XPLIT file
> version: 53, expected version: 5 or [8 or 49]` — o formato binário de plot
> andou e a biblioteca ficou atrás. O que **funciona**, medido em 2026-10-12 no
> mesmo modelo: `Feb(arquivo)` lê o `.feb`, `febio_python.feb.run.run()` chama
> o solver, e a saída numérica sai por `<Output><logfile>` em texto delimitado,
> que entra em `pandas` direto. Então: **plot binário só para a GUI do FEBio
> Studio; para número em Python, declare `<logfile>` no `.feb`.**
>
> ```xml
> <logfile>
>   <node_data    file="deslocamento.txt" data="ux;uy;uz" delim=","/>
>   <element_data file="tensao.txt"       data="sx;sy;sz" delim=","/>
> </logfile>
> ```
>
> O logfile intercala um cabeçalho `*Step = n` antes de cada bloco, então o
> parse é um `re.split` por passo e um `read_csv` por bloco. Validado contra
> física conhecida: cubo unitário neo-Hookeano (E=1, ν=0,3) em tração
> prescrita `uz=0,2` convergiu em 10 passos, 3 iterações por passo,
> `N O R M A L   T E R M I N A T I O N`; os nós do topo chegaram a `uz=0,2`
> exatos com contração lateral de 0,0386779723 — Poisson, como manda.

> **`C:/Program Files/FEBioStudio/bin` entrou no PATH do usuário.** Sem isso
> `febio_python.feb.run.run()` falha: quando o argumento `executable` é omitido
> ele resolve o nome `febio4.exe` pelo PATH, e o instalador do FEBio não mexe
> no PATH. Acrescentado em 2026-10-12 via
> `[Environment]::SetEnvironmentVariable('Path', …, 'User')`.

### Clones de terceiros (`vendor/`)

Reconstruir: `bash vendor/clonar.sh`. Não são versionados aqui — são clones.

| Clone | Licença | Nota |
|---|---|---|
| `vendor/BoneJ_Headless` | BSD-3-Clause | scripts que pilotam o Fiji+BoneJ pela API do ImageJ; não existe no PyPI nem no conda |
| `vendor/opencap-processing` | Apache-2.0 | scripts que importam `opensim`; rodam no ambiente `~/envs/opencap` do WSL |

O repositório inteiro tem 415 MB e quase tudo é imagem de exemplo em `ROIs/` e
`Examples/`; o código são ~20 KB. O checkout é esparso (332 KB em disco). Para
trazer os exemplos:
`git -C vendor/BoneJ_Headless sparse-checkout add Examples ROIs`.

### Já estava na stack — não reinstalar

| Item da lista | Onde já está |
|---|---|
| Marketplace life-sciences | 15 plugins instalados (PubMed, BioRender, ChEMBL, ClinicalTrials, bioRxiv, Open Targets, Owkin, Synapse, 10x Genomics, nextflow-development, scvi-tools, single-cell-rna-qc, …) |
| **K-Dense Scientific Skills** | **`.agents/skills/` — as 177 skills do upstream, todas as 177 presentes** |
| **NVIDIA BioNeMo Agent Toolkit** | **`.agents/skills/` — as 37 restantes: `boltz2-nim`, `openfold2-nim`, `openfold3-nim`, `evo2-nim`, `rfdiffusion-nim`, `proteinmpnn-nim`, `diffdock-nim`, `genmol-nim`, `molmim-nim`, `msa-search-nim`, `codonfm-*`, `complexa-*`, `cuequivariance`** |
| ToolUniverse | MCP `tooluniverse` no `.mcp.json` |
| Claude Science | MCP `claude-science` (`servers/claude_science_mcp/`) |
| Boltz-2, OpenFold3, Evo 2 | `mcp__nvidia-bio__*` nos endpoints hospedados |
| OpenSim Moco | dentro do OpenSim 4.6 (`MocoStudy`/`MocoTrack`/`MocoInverse`) |
| oVert (openVertebrate) | projeto de dados, não software — acessado pelo MorphoSource |
| Paleobiology Database | pacote R `paleobioDB` |

> **O K-Dense não entrou como marketplace de propósito, e isso é acerto, não
> lacuna.** `.agents/skills/` tem 214 diretórios; os 177 do
> `K-Dense-AI/scientific-agent-skills` estão **todos** lá (interseção 177,
> nenhum só no upstream), e os 37 extras são o BioNeMo Agent Toolkit.
> Acrescentar o marketplace faria as duas fontes somarem — a armadilha que
> `.claude/marketplace/plugins/clonar.sh` já documenta ("o bioinfo-213 reportou
> 426 skills para 213 em disco"). O repositório também mudou de nome:
> `K-Dense-AI/claude-scientific-skills` responde **301** para
> `K-Dense-AI/scientific-agent-skills`, e não traz
> `.claude-plugin/marketplace.json`, então `marketplace add` recusaria de todo
> jeito. Para quem tentar: o clone precisa de
> `CLAUDE_CODE_PLUGIN_GIT_TIMEOUT_MS` acima de 120 s (são 474 MB) e a URL tem
> de ser HTTPS, porque o `marketplace add` tenta `git@github.com` e a chave de
> host do GitHub não está conhecida nesta máquina.

### Falta a mão de Мухамед — bloqueado por licença ou conta

Nada aqui é limitação técnica; é um aceite que não é meu para dar.

| Item | O que falta | Depois |
|---|---|---|
| 3D Slicer → Web Server | passo de sessão, não de instalação: abrir o módulo **Web Server** e iniciar | o MCP `slicer` responde com 3 ferramentas |
| Blender → Start MCP Server | addon já habilitado e persistido; falta clicar | o MCP `blender` responde com 9 ferramentas |
| App Claude Science | só macOS e Linux; **não roda no Windows** | o MCP `claude-science` já cobre o daemon por REST |
| `molecule-mcp` (ChatMol) | o README exige Claude Desktop, que não existe aqui | `pymol` 3.1.0 está no WSL e o MCP `chimerax` cobre o outro lado |
| DeepLabCut SuperAnimal | os pesos pré-treinados são **somente pesquisa, sem uso comercial** | o pacote está instalado; a restrição é dos pesos, não do código |

## Skills

| Coleção | Tamanho | Destaque |
|---|---|---|
| ★ [NVIDIA-BioNeMo/bionemo-agent-toolkit](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit) | ~15 skills + workflows | a camada NVIDIA; NIMs, Parabricks, nvMolKit, cuEquivariance |
| ★ [anthropics/life-sciences](https://github.com/anthropics/life-sciences) | 5 skills + 8 MCPs | oficial Anthropic; scRNA-QC, Nextflow, scvi-tools |
| ★ [ClawBio/ClawBio](https://github.com/ClawBio/ClawBio) | 97 skills | local-first, sem API externa; Galaxy, PRS, farmacogenômica |
| [GoekeLab/awesome-genomic-skills](https://github.com/GoekeLab/awesome-genomic-skills) | lista curada | o mapa do ecossistema inteiro |
| [google-deepmind/science-skills](https://github.com/google-deepmind/science-skills) | ~36 | genômica, estrutural, quimioinformática, literatura |
| [openai/plugins — life-science-research](https://github.com/openai/plugins/tree/main/plugins/life-science-research) | 50 | genética humana, genômica funcional, vias |
| [GPTomics/bioSkills](https://github.com/GPTomics/bioSkills) | ponta a ponta | RNA-seq, variantes, ChIP, scRNA, espacial, Hi-C, proteômica |
| [swaruplab/operon](https://github.com/swaruplab/operon) | 180+ protocolos | protocolos detalhados em formato SKILL.md |
| [mims-harvard/ToolUniverse](https://github.com/mims-harvard/ToolUniverse) | 1.000+ APIs | 68 skills prontas sobre modelos e bases |
| [aws-samples/sample-kiro-power-life-sciences](https://github.com/aws-samples/sample-kiro-power-life-sciences) | 24 MCPs | 100+ bases via AWS |
| [HKU-BAL/Clair-skills](https://github.com/HKU-BAL/Clair-skills) | focada | suíte Clair de chamada de variantes |
| [jaechang-hits/SciAgent-Skills](https://github.com/jaechang-hits/SciAgent-Skills) | 197 | multi-domínio científico |
| [K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills) | 135 | ciência geral |
| [InternScience/Awesome-Scientific-Skills](https://github.com/InternScience/Awesome-Scientific-Skills) | curada | pesquisa científica ampla |

> Instalar tudo é contraproducente: cada skill ocupa espaço de descoberta e
> aumenta a chance de o modelo escolher a ferramenta errada. Comece com as três
> marcadas e adicione conforme a dor aparecer.

---

## Agentes desta stack

Definidos em `.claude/agents/`, escritos para este hardware e este contexto.

| Agente | Domínio |
|---|---|
| `structure-bio` | dobramento, docking, afinidade, interpretação de pLDDT/PAE |
| `omics-analyst` | RNA-seq bulk, single-cell, enriquecimento funcional |
| `seq-designer` | plasmídeos, primers, códons, Golden Gate, guias CRISPR |
| `variant-curator` | anotação de variantes, ACMG, gnomAD, farmacogenômica |
| `pipeline-engineer` | Nextflow, nf-core, Snakemake, cluster, containers |
| `lit-scout` | revisão de literatura, mineração de bases |
| `evo-devo` | perda gênica, estado ancestral, polarização por grupo externo, atavismo |
| `phylo-analyst` | alinhamento, IQ-TREE 3, suporte de ramo, concordância, árvore de espécies — ver [`FILOGENIA.md`](FILOGENIA.md) |

---

## Hooks

| Hook | Evento | O que faz |
|---|---|---|
| `session_context.py` | SessionStart | injeta GPU, ferramentas presentes, chaves, histórico |
| `guard_data.py` | PreToolUse (Bash) | bloqueia remoção recursiva em dados brutos, push forçado, curl-pipe-bash |
| `lint_gate.py` | PostToolUse (Write/Edit) | roda ruff/basedpyright/nextflow lint e devolve diagnósticos |
| `provenance.py` | PostToolUse (tudo) | grava caderno de laboratório em `runs/provenance.jsonl` |

---

## Benchmarks

Para avaliar se a stack está ajudando de fato:

| Benchmark | Fonte | Mede |
|---|---|---|
| BioAgent Bench | [bioagent-bench](https://github.com/bioagent-bench/bioagent-bench) | 10 pipelines ponta a ponta com artefato verificável |
| BixBench | [Future-House/BixBench](https://github.com/Future-House/BixBench) | tarefas reais de biologia computacional |
| CompBioBench | [Genentech/compbiobench-runner](https://github.com/Genentech/compbiobench-runner) | 100 questões: single-cell, epigenômica, genômica |
| LAB-Bench | [futurehouse/lab-bench](https://huggingface.co/datasets/futurehouse/lab-bench) | 2.400+ questões de biologia prática |
| SkillsBench | [benchflow-ai/skillsbench](https://github.com/benchflow-ai/skillsbench) | se as skills ajudam mesmo, 86 tarefas |
| BiomniBench | [phylobio/BiomniBench-DA](https://huggingface.co/datasets/phylobio/BiomniBench-DA) | 100 tarefas de análise tiradas de artigos |

---

## Bibliotecas GPU da NVIDIA

| Biblioteca | Para quê | Roda em 8 GB? |
|---|---|---|
| [nvMolKit](https://github.com/NVIDIA-BioNeMo/nvMolKit) | quimioinformática acelerada (fingerprints, conformers) | sim |
| [bionemo-recipes](https://github.com/NVIDIA-BioNeMo/bionemo-recipes) | treinar/adaptar modelos de descoberta de fármacos | depende |
| [bionemo-framework](https://github.com/NVIDIA-BioNeMo/bionemo-framework) | treino de ESM-2, Evo2, CodonFM, Llama3 | não (treino) |
| [Proteina-Complexa](https://github.com/NVIDIA-BioNeMo/Proteina-Complexa) | desenho de binder por flow matching | não |
| [genmol](https://github.com/NVIDIA-BioNeMo/genmol) | geração de moléculas por difusão discreta | provavelmente |
| [CodonFM](https://github.com/NVIDIA-BioNeMo) | modelo de RNA em nível de códon (130M sequências) | inferência sim |
| [RNAPro](https://github.com/NVIDIA-BioNeMo) | estrutura 3D de RNA | não |
| Parabricks | genômica GPU (fq2bam, DeepVariant) | não, ≥16 GB |
| RAPIDS / rapids-singlecell | scRNA-seq 10–50× mais rápido | sim, **só em WSL2** |
