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

> **Por que o `bibtexparser<2` no `biocontext`.** O pacote importa `scholarly`,
> que faz `from bibtexparser.bibdatabase import BibDatabase`. O bibtexparser 2.x
> removeu esse módulo, então o `uvx` resolvia a versão nova e o servidor morria
> no import com `ModuleNotFoundError` — aparecia no Claude Code apenas como
> `CONNECTION_CLOSED`, sem pista da causa. Diagnosticado em 2026-09-12 rodando
> o handshake stdio na mão. Com o pin, o servidor sobe (`BioContextAI 2.13.1`).

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

### Ferramentas de linha de comando (exigem WSL2 + bioconda)

Organização [bio-mcp](https://github.com/bio-mcp): `bio-mcp-blast`,
`bio-mcp-samtools`, `bio-mcp-bcftools`, `bio-mcp-seqkit`, `bio-mcp-bwa`,
`bio-mcp-bedtools`, `bio-mcp-fastqc`, `bio-mcp-interpro`, `bio-mcp-amber`,
`bio-mcp-evo2`, `bio-mcp-queue` (fila para jobs longos).

Config de exemplo, depois de instalar as ferramentas no WSL:

```jsonc
{
  "mcpServers": {
    "blast": { "command": "wsl", "args": ["-e", "uvx", "bio-mcp-blast"] },
    "samtools": { "command": "wsl", "args": ["-e", "uvx", "bio-mcp-samtools"] }
  }
}
```

### Registros para consultar

- [BioContextAI Registry](https://github.com/biocontext-ai/registry/) — catálogo curado de MCPs biomédicos
- [MCPmed](https://github.com/MCPmed) — implementações de referência (GEO, STRING, UCSC Cell Browser, PLSDB)
- [awesome-mcp-servers #bio](https://github.com/punkpeye/awesome-mcp-servers#bio)
- [BioinfoMCP](https://github.com/florensiawidjaja/BioinfoMCP) — gera MCP automaticamente a partir da documentação de uma ferramenta CLI

---

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
