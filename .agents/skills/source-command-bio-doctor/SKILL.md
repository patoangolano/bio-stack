---
name: "source-command-bio-doctor"
description: "Diagnostica a bio-stack — ferramentas, chaves, GPU, servidores MCP e endpoints NVIDIA — e diz exatamente o que consertar."
---

# source-command-bio-doctor

Use this skill when the user asks to run the migrated source command `bio-doctor`.

## Command Template

Faça um diagnóstico completo da bio-stack e entregue um relatório acionável.

## 1. Ferramentas de base

```bash
nvidia-smi --query-gpu=name,memory.total,memory.used,compute_cap,driver_version --format=csv
```

Depois verifique a presença de: `git`, `node`, `npm`, `uv`, `uvx`, `python`,
`R`. No PowerShell (`nextflow`, `docker` e `samtools` ficam de fora de
propósito: só existem no WSL, e MISS deles aqui é o esperado — a presença
real é conferida na seção 1b):

```powershell
'git','node','npm','uv','uvx','python','R' | ForEach-Object { $c = Get-Command $_ -ErrorAction SilentlyContinue; if ($c) { "OK   $_" } else { "MISS $_" } }
```

## 1b. Camada de CLIs no WSL (pixi)

As ferramentas de bioinformática vivem no WSL, num ambiente pixi por
ferramenta, declaradas em `bootstrap/pixi-global.toml`. Sem cifrão na linha
de comando e com `bash -lc`:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "pixi --version && pixi global list"
wsl -d Ubuntu-24.04 -- bash -lc "command -v nextflow docker iqtree3 mafft hyphy codeml blastn diamond miniprot datasets nf-test samtools"
```

Ferramenta faltando: `wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/05_wsl_ferramentas.sh`.
`samtools`/`bcftools` respondem pelo pixi (1.24), não pelo apt (1.19).

## 1c. Pipeline de filogenia

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "cd /mnt/c/Users/Manec/bio-stack/pipelines/filogenia && nextflow lint main.nf nextflow.config conf modules/local"
uv run --script tests/check_filogenia.py
```

## 2. Chaves de API

Verifique se estão definidas (sem imprimir o valor):

```powershell
'NVIDIA_API_KEY','NCBI_API_KEY' | ForEach-Object { $v = [Environment]::GetEnvironmentVariable($_, 'User'); if ($v) { "OK   $_ ($($v.Length) chars)" } else { "MISS $_" } }
```

> O `S2_API_KEY` saiu desta sonda em 2026-10-10: era repassado ao `biomcp`, que
> não o lê. Ver `docs/CHAVES.md` seção 2. **Não reponha** sem achar um
> consumidor real.

## 3. Camada NVIDIA

Chame `nim_config` e depois `nim_probe_endpoints`. Reporte quais modelos
respondem e quais deram 404 (slug mudou) ou 403 (licença não cobre).

## 4. Ambiente Python

Se `.venv` existir, liste os pacotes-chave:

```powershell
.\.venv\Scripts\python.exe -c "import importlib.metadata as m
for p in ['biopython','scanpy','pydeseq2','gget','rdkit','torch']:
    try: print(f'{p}: {m.version(p)}')
    except m.PackageNotFoundError: print(f'{p}: AUSENTE')"
```

> Consulte os **metadados da distribuição**, não o módulo importável: o pacote
> `biopython` importa como `Bio`, e checar por `find_spec("biopython")` reporta
> "AUSENTE" mesmo com ele instalado.

E se torch estiver instalado, confirme CUDA:

```powershell
.\.venv\Scripts\python.exe -c "import torch; print('cuda', torch.cuda.is_available(), torch.version.cuda if torch.cuda.is_available() else '')"
```

## 5. Relatório

Entregue uma tabela com três colunas: **componente | estado | o que fazer**.

Ordene por impacto: o que bloqueia mais trabalho primeiro. Para cada item
quebrado, dê o comando exato de correção, não uma descrição.

Se tudo estiver verde, diga isso em uma linha e liste os agentes disponíveis
(`structure-bio`, `omics-analyst`, `seq-designer`, `variant-curator`,
`evo-devo`, `phylo-analyst`, `pipeline-engineer`, `lit-scout`) com um
exemplo de uso de cada.

Não conserte nada sem perguntar — este comando é só diagnóstico.
