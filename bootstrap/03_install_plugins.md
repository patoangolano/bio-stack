# Camada 3 — Skills e plugins

Rode estes comandos **dentro do Claude Code**, com o projeto `bio-stack` aberto.
Exigem `git` e `node` no PATH (camada 0).

---

## 3.1 BioNeMo Agent Toolkit — NVIDIA

O pacote central da camada NVIDIA. Skills para dobramento, docking, química
generativa, genômica e desenho de proteína, escritas pela própria NVIDIA.

```
/plugin marketplace add NVIDIA-BioNeMo/bionemo-agent-toolkit
/plugin install bionemo-agent-toolkit@bionemo-agent-toolkit
```

Alternativa pela CLI `skills` (fora do Claude Code):

```bash
npx skills add NVIDIA-BioNeMo/bionemo-agent-toolkit --agent claude-code
npx skills add NVIDIA-BioNeMo/bionemo-agent-toolkit --list
npx skills add NVIDIA-BioNeMo/bionemo-agent-toolkit --skill boltz2-nim --yes
```

**Skills incluídas:**

| Categoria | Skills |
|---|---|
| NIM | boltz2, diffdock, evo2, genmol, molmim, msa-search, openfold2, openfold3, proteinmpnn, rfdiffusion |
| Modelos abertos | proteina-complexa (desenho de binder), kermt (propriedade molecular) |
| Bibliotecas | parabricks (genômica GPU), nvmolkit (quimioinformática), cuequivariance |
| Workflows | pipeline de descoberta de fármacos, desenho de binder, estrutura com MSA |

**Variável exigida:** `NGC_API_KEY`. A mesma chave do `NVIDIA_API_KEY` serve —
defina as duas apontando para o mesmo valor:

```powershell
$k = [Environment]::GetEnvironmentVariable('NVIDIA_API_KEY','User')
[Environment]::SetEnvironmentVariable('NGC_API_KEY', $k, 'User')
```

> **Skills vs. o servidor MCP deste projeto.** As skills da NVIDIA ensinam o
> procedimento científico; o servidor `nvidia-bio` em `servers/` trata o
> protocolo (polling assíncrono, retry, gravação em disco). Os dois convivem —
> use as skills para saber *o que* fazer e as ferramentas MCP para *executar*
> sem estourar o contexto.

---

## 3.2 Life Sciences — Anthropic

Bundle oficial: skills + conectores MCP de parceiros.

```
/plugin marketplace add https://github.com/anthropics/life-sciences.git
/plugin install pubmed@life-sciences
/plugin install single-cell-rna-qc@life-sciences
/plugin install nextflow-development@life-sciences
/plugin install scvi-tools@life-sciences
```

Skills: QC de scRNA-seq, dados de instrumento → Allotrope, desenvolvimento
Nextflow (rnaseq/sarek/atacseq), scvi-tools, seleção de problema científico.

MCP de parceiros: PubMed (sem auth), BioRender, Synapse.org, Scholar Gateway
(Wiley), Consensus, Cortellis, AdisInsight, 10x Genomics Cloud. Os autenticados
pedem conta — configure em `/plugin` → Manage plugins → Configure.

---

## 3.3 BioMCP — dados biomédicos

```powershell
uv tool install biomcp-python
biomcp health --apis-only
```

O `.mcp.json` deste projeto chama `biomcp run --mode stdio` — `biomcp serve`
não existe como subcomando.

Depois, no Claude Code:

```
/plugin marketplace add genomoncology/biomcp
/plugin install biomcp@biomcp
```

O `.mcp.json` deste projeto já aponta para `biomcp serve` — se instalar pelo
`uv tool`, o servidor sobe sozinho na próxima sessão.

Cobre ~30 fontes: PubMed, ClinVar, ClinicalTrials.gov, MyVariant, OncoKB,
Reactome, openFDA. Entidades: gene, variant, article, trial, drug, disease,
pathway, protein, adverse-event, pgx.

---

## 3.4 ClawBio — skills locais de genômica

97 skills reproduzíveis, sem depender de API externa.

```
/plugin marketplace add ClawBio/ClawBio
/plugin install clawbio
```

Cobre: chamada de variantes, RNA-seq, genética de populações, PRS, ancestralidade,
farmacogenômica, single-cell, integração com Galaxy (8.000+ ferramentas).

> O modo MCP do ClawBio foi **descontinuado na v0.7.0**. Use como skills.

---

## 3.5 Coleções amplas (opcional)

Instalar tudo é contraproducente: cada skill ocupa espaço de descoberta e
aumenta a chance de o modelo escolher a ferramenta errada. Escolha no máximo uma.

**Atenção ao método de instalação** — os repositórios abaixo não são todos
marketplaces de plugin, e `/plugin marketplace add` falha nos que não são.

### Via marketplace de plugin

```
/plugin marketplace add mims-harvard/ToolUniverse
```

[`ToolUniverse`](https://github.com/mims-harvard/ToolUniverse) — 1.000+ APIs e
modelos científicos, com 68 skills prontas.

### Via CLI `skills` (padrão Agent Skills)

```bash
npx skills add GPTomics/bioSkills --agent claude-code
npx skills add google-deepmind/science-skills --agent claude-code
```

| Repositório | Cobre |
|---|---|
| [`GPTomics/bioSkills`](https://github.com/GPTomics/bioSkills) | pipelines ponta a ponta: RNA-seq, variantes, ChIP, scRNA, espacial, Hi-C |
| [`google-deepmind/science-skills`](https://github.com/google-deepmind/science-skills) | ~36 skills: genômica, estrutural, quimioinformática, literatura |

### Só para leitura, não instalável

[`GoekeLab/awesome-genomic-skills`](https://github.com/GoekeLab/awesome-genomic-skills)
é uma lista curada do ecossistema — o mapa, não o território. Vale ler antes de
decidir o que instalar.

[`swaruplab/operon`](https://github.com/swaruplab/operon) é um aplicativo
desktop separado com 180+ protocolos, fora do fluxo de plugins do Claude Code.

---

## 3.6 Verificação

```
/plugin
/bio-doctor
```

`/bio-doctor` roda o diagnóstico completo e diz o que ainda falta.
