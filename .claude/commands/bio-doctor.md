---
description: Diagnostica a bio-stack — ferramentas, chaves, GPU, servidores MCP e endpoints NVIDIA — e diz exatamente o que consertar.
allowed-tools: Bash, Read, Glob, mcp__nvidia-bio__nim_config, mcp__nvidia-bio__nim_probe_endpoints
---

Faça um diagnóstico completo da bio-stack e entregue um relatório acionável.

## 1. Ferramentas de base

```bash
nvidia-smi --query-gpu=name,memory.total,memory.used,compute_cap,driver_version --format=csv
```

Depois verifique a presença de: `git`, `node`, `npm`, `uv`, `uvx`, `python`,
`nextflow`, `docker`, `samtools`, `R`. No PowerShell:

```powershell
'git','node','npm','uv','uvx','python','nextflow','docker','samtools','R' | ForEach-Object { $c = Get-Command $_ -ErrorAction SilentlyContinue; if ($c) { "OK   $_" } else { "MISS $_" } }
```

## 2. Chaves de API

Verifique se estão definidas (sem imprimir o valor):

```powershell
'NVIDIA_API_KEY','NCBI_API_KEY','S2_API_KEY' | ForEach-Object { $v = [Environment]::GetEnvironmentVariable($_, 'User'); if ($v) { "OK   $_ ($($v.Length) chars)" } else { "MISS $_" } }
```

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
`pipeline-engineer`, `lit-scout`) com um exemplo de uso de cada.

Não conserte nada sem perguntar — este comando é só diagnóstico.
