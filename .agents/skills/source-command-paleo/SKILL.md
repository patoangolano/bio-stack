---
name: source-command-paleo
description: Executa análises de paleontologia computacional, datação de pontas (FBD/PBDB), matrizes morfológicas (Mk/Mkv +ASC), paleoproteômica (ZooMS/desamidação), DNA antigo (aDNA) e reconstrução de proteínas ancestrais (ASR + AltAll + Boltz-2/OpenFold3).
---

Alvo paleontológico (táxons, matriz, sequência ou fóssil): $ARGUMENTS

Delegue ao agente `paleo-analyst` ou execute diretamente o fluxo abaixo usando
`scripts/paleo/paleo_toolkit.py` e as ferramentas do WSL (`Ubuntu-26.04`).

## Passo 1 — Classificar a janela tafonômica e o Extant Phylogenetic Bracket (EPB)

1. Identifique a idade aproximada do material ou clado:
   - **0 – ~1,5 Ma**: DNA antigo (`aDNA` / `sedaDNA`) viável (`adapterremoval`, `bwa aln`, `damageprofiler`, `mapDamage`, `angsd`).
   - **0 – ~48 Ma**: Paleoproteômica intracristalina (esmalte dentário: `AMELX`, `ENAM`, `AMBN`, `MMP20`, `ODAM`; osso/dentina: `COL1A1`, `COL1A2`).
   - **> 50 Ma**: Sem DNA endógeno; use **Extant Phylogenetic Bracket (Witmer 1995)**, **ASR + AltAll** (`iqtree3 --ancestral` + `paleo_toolkit.py altall-asr`), **matrizes morfológicas discretas** (`mk-audit`) e **Fossilized Birth-Death (FBD)** (`rb` / `beast` / `mb`).
2. Declare os dois clados viventes que flanqueiam o táxon extinto (ex.: Crocodylia + Aves para Dinosauria/Archosauria) e o nível de inferência de Witmer (I, II ou III).

## Passo 2 — Calibrar idades estratigráficas (FAD/LAD) na Paleobiology Database (PBDB)

```powershell
uv run python scripts/paleo/paleo_toolkit.py pbdb-range --taxa "Tyrannosaurus rex,Archaeopteryx lithographica,Spinosaurus aegyptiacus,Alligator mississippiensis"
```

## Passo 3 — Auditar matriz morfológica (Mk/Mkv +ASC)

```powershell
uv run python scripts/paleo/paleo_toolkit.py mk-audit --input <matriz.nex>
```

## Passo 4 — Reconstrução de Sequência Ancestral (ASR) com controle AltAll (Eick et al. 2017)

```powershell
uv run python scripts/paleo/paleo_toolkit.py altall-asr --state-file runs/paleo/asr.state --node Node1 --pp-threshold 0.80 --out-fasta runs/paleo/ancestral_pair.fa
```

## Passo 5 — Controle de autenticidade em Paleoproteômica / ZooMS ou aDNA

```powershell
uv run python scripts/paleo/paleo_toolkit.py zooms-qc --sequence "<sequencia_AMELX_ou_COL1A1>"
```
