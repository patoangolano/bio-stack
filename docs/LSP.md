# Camada LSP

## Duas metades que se completam

O Claude Code **tem cliente LSP nativo** (a partir da v2.0.74), ativado por
plugins de code intelligence do marketplace oficial. Isso dá ao modelo o que
um editor tem: tipos, definições, referências, renomeação, diagnósticos ao vivo.

Mas o LSP cobre *linguagens*. Bioinformática vive de coisas que não têm
language server: Nextflow, Snakemake, samplesheets CSV, schemas de nf-core.

Então a stack usa as duas metades:

| Metade | Quem fornece | Cobre |
|---|---|---|
| **LSP nativo** | plugin `pyright-lsp` | Python: tipos, go-to-definition, referências |
| **Ponte de linters** | hook `lint_gate.py` | ruff, Nextflow, Snakemake, R, JSON/YAML |

---

## Metade 1 — LSP nativo

```
/plugin marketplace add anthropics/claude-plugins-official
/plugin install pyright-lsp@claude-plugins-official
```

Python é a linguagem desta stack, então `pyright-lsp` é o que importa. O
marketplace oficial cobre ~12 implementações (TypeScript, Rust, Go, Java,
C/C++, entre outras) caso o projeto cresça para fora de Python.

Depois de instalar, o modelo passa a resolver símbolos em vez de procurar
texto: ele acha a definição de `DeseqDataSet` na sua venv em vez de grepar por
`class DeseqDataSet`.

> **Nota:** houve um bug conhecido (issue #16219) em que a configuração
> `lspServers` do marketplace não era propagada para o plugin instalado.
> Se o LSP não responder depois de instalar, é o primeiro lugar para olhar.

---

## Metade 2 — ponte de linters

`.claude/hooks/lint_gate.py` é um `PostToolUse` que dispara em `Write|Edit|NotebookEdit`:

```
Claude edita pipeline.nf
   ↓
hook detecta .nf
   ↓
nextflow lint
   ↓
diagnósticos → stderr → exit 2
   ↓
Claude recebe os erros e corrige
```

Cobertura por extensão:

| Extensão | Analisador | Já coberto pelo LSP? |
|---|---|---|
| `.py` | `ruff check` + `basedpyright --outputjson` | tipos sim; regras do ruff não |
| `.nf`, `nextflow.config` | `nextflow lint` | **não existe LSP instalável aqui** |
| `.smk`, `Snakefile` | `snakemake --lint` | **não** |
| `.R` | `Rscript -e lintr::lint` | parcial (`languageserver`) |
| `.json`, `.yaml` | validação de sintaxe | não |

O ganho real está nas linhas sem LSP. Um erro de sintaxe num `.nf` que só
apareceria depois de 40 minutos de pipeline volta na mesma iteração.

### Como os linters são encontrados

O hook resolve cada ferramenta **por PATH ou por `uvx`**:

```python
def tool_cmd(tool, *args):
    if shutil.which(tool):
        return [tool, *args]
    return [UVX, "tool", "run", "--quiet", tool, *args]
```

Isso não é detalhe de implementação — é o que faz a ponte funcionar. Os
linters instalados por `02_python_env.ps1 -Dev` vão para `.venv/Scripts/`, que
**não está no PATH do processo do hook**. Uma versão anterior deste hook
resolvia só por `shutil.which` e ficava completamente inerte numa instalação
que seguiu a própria documentação — silenciosamente, sem erro nenhum.

### Calibragem

`pyproject.toml` na raiz fixa `typeCheckingMode = "standard"` para o
basedpyright. O padrão dele é `recommended`, mais estrito que o Pylance, e
devolveria uma parede de diagnósticos em qualquer script científico normal
(`np.ndarray` sem parâmetros de tipo, `pd.DataFrame` não tipado, etc.).

---

## Instalação dos analisadores

```powershell
# Recomendado: uv tool install poe no ~/.local/bin, que ESTA no PATH
uv tool install ruff
uv tool install basedpyright

# Nextflow e Snakemake (dentro do WSL)
wsl -e bash -c "curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/"
uv tool install snakemake

# R
Rscript -e 'install.packages(c("lintr","languageserver"), repos="https://cloud.r-project.org")'
```

`uv tool install` em vez de `uv pip install`: o primeiro põe o binário em
`~/.local/bin` (no PATH), o segundo esconde dentro da venv. Com `uvx` disponível
o hook funciona nos dois casos, mas o PATH direto evita o custo de resolução a
cada chamada.

---

## Editor — VS Code

O LSP do Claude Code e o do seu editor são independentes: instalar um não
configura o outro.

| Extensão | ID | Para quê |
|---|---|---|
| Python + Pylance | `ms-python.python`, `ms-python.vscode-pylance` | base |
| Ruff | `charliermarsh.ruff` | lint/format em tempo real |
| Nextflow | `nextflow.nextflow` | LS oficial: completar, diagnóstico, DAG |
| R | `REditorSupport.r` | Bioconductor |
| Jupyter | `ms-toolsai.jupyter` | notebooks |
| Rainbow CSV | `mechatroner.rainbow-csv` | ler samplesheet sem enlouquecer |
| Even Better TOML | `tamasfe.even-better-toml` | `pyproject.toml` |

Copie `editor/vscode-settings.json` para `.vscode/settings.json`.

O language server do Nextflow é um JAR Java autônomo e funciona com qualquer
cliente LSP (Neovim, Helix, Emacs) — só a integração com VS Code é oficialmente
suportada, mas o protocolo é padrão.

---

## O que não existe

- Não há language server de FASTA, GenBank, VCF ou PDB. São formatos de dados,
  não linguagens. Validação desses arquivos se faz com `seqkit`,
  `bcftools view -h`, `gt` — não com LSP.
- Não há LSP de "biologia". Sequência de DNA não tem gramática para analisar
  desse jeito.
