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
VRAM:   8 GB
CUDA:   13.4 (driver 616.92)
Arch:   Blackwell, sm_120  → PyTorch precisa de wheels cu128+
SO:     Windows 11 · WSL2 presente (distro docker-desktop), sem Ubuntu
```

**Não proponha rodar localmente**: AlphaFold2, ESMFold, OpenFold2/3, Boltz-2,
Evo 2, RFdiffusion, Parabricks. Nenhum cabe em 8 GB. Use as ferramentas
`mcp__nvidia-bio__*`, que batem nos endpoints hospedados da NVIDIA.

**Cabe local**: ESM-2 650M, DNABERT-2, scanpy, pydeseq2, RDKit, biopython,
nvMolKit, inferência de modelos pequenos.

**Precisa de WSL2**: nextflow, samtools/bcftools, rapids-singlecell, pysam.
O kernel do WSL2 **já está instalado** — Docker Desktop roda sobre ele
(`wsl -l -v` mostra `docker-desktop`, versão 2). Falta só a distro:
`wsl --install -d Ubuntu-24.04`, que não exige reinício nem habilitar recurso
do Windows. Peça isso antes, em vez de tentar no PowerShell.

## Roteamento de agentes

| A pergunta é sobre… | Delegue para |
|---|---|
| estrutura 3D, docking, afinidade, pLDDT | `structure-bio` |
| desenho de backbone/binder, geração de molécula | `structure-bio` |
| matriz de contagem, DEG, scRNA-seq, enriquecimento | `omics-analyst` |
| plasmídeo, primer, códon, Golden Gate, CRISPR guide | `seq-designer` |
| variante, VCF, ACMG, gnomAD, farmacogenômica | `variant-curator` |
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
- `/metodos` — gera a seção de Métodos a partir do log
- `/nvidia` — audita tudo que a chave NVIDIA destrava

## Estrutura

```
bootstrap/   instalação (rodar uma vez)
servers/     servidor MCP da NVIDIA (código próprio)
.claude/     agents, commands, hooks, settings
data/raw/    dados brutos — protegidos por hook, nunca sobrescrever
data/proc/   dados processados
runs/        saídas de análise + provenance.jsonl
docs/        catálogo da stack, referência NVIDIA, revisões de literatura
notebooks/   análise exploratória
```
