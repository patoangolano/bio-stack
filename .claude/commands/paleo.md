---
description: Executa análises de paleontologia computacional, datação de pontas (FBD/PBDB), matrizes morfológicas (Mk/Mkv +ASC), paleoproteômica (ZooMS/desamidação), DNA antigo (aDNA) e reconstrução de proteínas ancestrais (ASR + AltAll + Boltz-2/OpenFold3).
argument-hint: <táxons, matriz morfológica, sequência ancestral ou fóssil>
---

Alvo paleontológico (táxons, matriz, sequência ou fóssil): $ARGUMENTS

Delegue ao agente `paleo-analyst` ou execute diretamente o fluxo abaixo usando
`scripts/paleo/paleo_toolkit.py` e as ferramentas do WSL (`Ubuntu-26.04`).

1. **Calibração estratigráfica na Paleobiology Database (PBDB API v1.2)**:
   `uv run python scripts/paleo/paleo_toolkit.py pbdb-range --taxa "..."`
2. **Auditoria de matriz morfológica discreta (Lewis 2001 Mk/Mkv +ASC)**:
   `uv run python scripts/paleo/paleo_toolkit.py mk-audit --input <matriz.nex>`
3. **Reconstrução ancestral com controle de viés de consenso AltAll (Eick et al. 2017)**:
   `uv run python scripts/paleo/paleo_toolkit.py altall-asr --state-file <asr.state> --node <Node> --pp-threshold 0.80`
4. **Digestão in silico e desamidação diagenética (ZooMS / Paleoproteômica)**:
   `uv run python scripts/paleo/paleo_toolkit.py zooms-qc --sequence "<seq>"`
