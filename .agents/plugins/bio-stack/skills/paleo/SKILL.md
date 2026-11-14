---
name: paleo
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

Nunca invente datas de calibração. Extraia `firstapp_max_ma` (FAD) e
`lastapp_min_ma` (LAD) diretamente da API v1.2 da PBDB:

```powershell
uv run python scripts/paleo/paleo_toolkit.py pbdb-range --taxa "Tyrannosaurus rex,Archaeopteryx lithographica,Spinosaurus aegyptiacus,Alligator mississippiensis"
```

Isso grava a consulta em `runs/provenance.jsonl` e devolve os blocos prontos
para:
- `RevBayes` (`dnUniform(lad_min, fad_max)`)
- `PAML MCMCTree` (limite mínimo suave em escala de 100 Ma)
- `IQ-TREE 3 LSD2` (`--date`)

## Passo 3 — Auditar matriz morfológica (Mk/Mkv +ASC) ou alinhamento molecular

Se houver matriz osteológica discreta (NEXUS ou FASTA):

```powershell
uv run python scripts/paleo/paleo_toolkit.py mk-audit --input <matriz.nex>
```

- Se `invariant_characters == 0`, use obrigatoriamente `+ASC` (`MK+FQ+ASC+G4`
  no `iqtree3` ou `coding="variable"` no `rb` / `mb`) para corrigir o viés de
  aquisição de Lewis (2001).

## Passo 4 — Reconstrução de Sequência Ancestral (ASR) com controle AltAll (Eick et al. 2017)

Quando reconstruir uma proteína ancestral para testar função ou estabilidade:

1. Rode `iqtree3` com `--ancestral` no WSL:
   ```powershell
   wsl -d Ubuntu-26.04 -- bash -lc "iqtree3 -s alinhamento.fa -m MFP --ancestral -pre runs/paleo/asr"
   ```
2. Gere o par `ML` + `AltAll` (troca todos os sítios com $PP < 0.80$ pelo 2º
   estado mais provável):
   ```powershell
   uv run python scripts/paleo/paleo_toolkit.py altall-asr --state-file runs/paleo/asr.state --node Node1 --pp-threshold 0.80 --out-fasta runs/paleo/ancestral_pair.fa
   ```
3. Dobre **ambas** as sequências (`ML` e `AltAll`) no `Boltz-2` ou `OpenFold3`
   (`mcp__nvidia-bio__fold_protein`) e compare RMSD e pLDDT. Se houver ligante,
   faça cross-check físico da pose com `smina` ou `vina` no WSL.

## Passo 5 — Controle de autenticidade em Paleoproteômica / ZooMS ou aDNA

- Para **Paleoproteômica / ZooMS**:
  ```powershell
  uv run python scripts/paleo/paleo_toolkit.py zooms-qc --sequence "<sequencia_AMELX_ou_COL1A1>"
  ```
  Reporte os peptídeos trípticos, massas monoisotópicas $[M+H]^+$ / $[M+2H]^{2+}$
  e o deslocamento exato de desamidação diagenética ($+0.984016\text{ Da}$ por N/Q).
- Para **aDNA** (no WSL `Ubuntu-26.04`):
  - Trimming e colapso: `adapterremoval --file1 R1.fq.gz --file2 R2.fq.gz --collapse`
  - Alinhamento sem semente: `bwa aln -l 1024 -n 0.01 -o 2 ref.facollapsed.fq > aln.sai`
  - Perfil de desaminação C→T / G→A: `damageprofiler -i aln.bam -r ref.fa -o dp_out/` ou `mapDamage -i aln.bam -r ref.fa`
  - Verossimilhança de genótipo: `angsd -bam bam.list -GL 2 -doMajorMinor 1 -doMaf 1`

## Passo 6 — Entregar

Entregue a tabela estratigráfica (PBDB `taxon_no` e FAD/LAD), os controles de
precisão aplicados (`+ASC`, `AltAll`, desamidação/dano C→T, `smina`/`vina`), e
um parágrafo explícito sobre o que o registro tafonômico não resolve. Tudo fica
registrado em `runs/provenance.jsonl` para `/metodos`.
