# Expansão de Paleontologia Computacional e Engenharia de Precisão da `bio-stack` (2026-11)

Levantamento bibliográfico primário, auditoria de precisão e instalação
verificada de ferramentas, modelos de IA, bases de dados e protocolos
estatísticos para projetos de **paleontologia computacional**, **paleoproteômica**,
**DNA antigo (aDNA)**, **filogenia em tempo profundo (*tip-dating* / FBD)** e
**paleobiomecânica 3D**.

---

## 1. Mapa de Precisão por Janela Tafonômica Temporal

O primeiro erro metodológico em paleobiologia computacional é aplicar a
ferramenta errada para a escala de tempo geológico do espécime. A `bio-stack`
passa a operar sob três regimes tafonômicos estritos:

| Janela Temporal | Biomolécula / Sinal Físico Sobrevivente | Ferramentas Instaladas & Verificadas na `bio-stack` | Controles de Precisão Obrigatórios |
| :--- | :--- | :--- | :--- |
| **0 – ~1,5 Ma** *(Pleistoceno; até ~2 Ma em permafrost polar)* | **DNA antigo (aDNA)** e DNA ambiental sedimentar (**sedaDNA**): fragmentos ultracurtos ($30\text{–}70\text{ pb}$) com desaminação hidrolítica de citosina ($\text{C}\to\text{U}\to\text{T}$) nas pontas 5' de fita simples. | `adapterremoval` (`3.0.2`), `bwa` (`0.7.19`), `damageprofiler` (`1.1`), `mapDamage2` (`2.2.3`), `angsd` (`0.940`), `vsearch` (`2.32.0`), `samtools`/`bcftools` (`1.24`). | **1)** Semeadura desativada no alinhamento (`bwa aln -l 1024 -n 0.01 -o 2`); **2)** Curva exponencial de dano $\text{C}\to\text{T}$ 5' e $\text{G}\to\text{A}$ 3' ($\delta_s, \delta_d, \lambda$) via `damageprofiler`/`mapDamage`; **3)** Verossimilhança de genótipo (`angsd -GL 2`) em vez de chamada dura diploide em cobertura $<5\times$. |
| **0 – ~48 Ma** *(Cenozoico / Eoceno Médio — esmalte intracristalino e casca de ovo calcítica)* | **Paleoproteômica e ZooMS**: proteínas intracristalinas do esmalte (`AMELX`, `AMBN`, `ENAM`, `MMP20`, `ODAM`) e colágeno ósseo/dentinário (`COL1A1`, `COL1A2`, `COL2A1`, `ALB`, `AHSG`). | `scripts/paleo/paleo_toolkit.py zooms-qc`, `ESM-2 650M` (local), `ProstT5` / `Casanovo` / `InstaNovo+` (arquiteturas *de novo*), `mafft` (`7.525`), `prank` (`251117`), `iqtree3` (`3.1.4`), `rb` (`1.4.2`). | **1)** Deslocamento exato de desamidação não-enzimática ($\text{Asn}\to\text{Asp}$ e $\text{Gln}\to\text{Glu}$, $\Delta m = +0.984016\text{ Da}$) e hidroxiprolina ($+15.994915\text{ Da}$); **2)** Clivagem diagenética preferencial em $\text{Asp-Pro}$ ($\text{D-P}$); **3)** Descarte de peptídeos com $0\%$ de desamidação como contaminação moderna. |
| **> 50 Ma** *(Mesozoico / Paleozoico — dinossauros não-avianos, arcossauros basais, sinápsidos)* | **Sem DNA endógeno e sem peptídeo longo livre de biomineral fechado.** O sinal vem de **matrizes osteológicas discretas**, **Reconstrução de Sequência Ancestral (ASR)** nos clados viventes (*Extant Phylogenetic Bracket* — Witmer 1995) e **micro-CT / FEA 3D**. | `scripts/paleo/paleo_toolkit.py` (`pbdb-range`, `mk-audit`, `altall-asr`), `rb` (`RevBayes 1.4.2`), `beast` (`BEAST 2.7.7` + `beagle-lib 4.2.0`), `mb` (`MrBayes 3.2.9`), `iqtree3` (`3.1.4`), `Boltz-2` / `OpenFold3` (NIM), `smina` (`2020.12.10`), `vina` (`1.1.2`), `SlicerMorph`, `nnInteractive`, `gmsh`, `febio-python`, `geomorph` (`4.1.1`). | **1)** Correção de viés de aquisição de Lewis (2001) (`Mkv` / `+ASC`) em matrizes morfológicas; **2)** Calibração estratigráfica `[LAD_min, FAD_max]` extraída da API v1.2 da **PBDB** sob processo **Fossilized Birth-Death (FBD)**; **3)** Protocolo **AltAll** (Eick et al. 2017) em toda proteína ancestral antes do dobramento 3D; **4)** Controle físico `smina`/`vina` para qualquer pose do `diffdock`. |

---

## 2. Literatura Primária Recente (2024–2026) que Fundamenta as Atualizações

1. **Sobrevivência de Paleoproteínas do Esmalte até o Eoceno (48 Ma) e Filogenia Molecular Profunda**:
   - **Gatti et al. (2026)** *Ancient amino acids from fossil teeth reveal deep-time evolutionary relationships.* *Communications Biology* — DOI: `10.1038/s42003-026-09716-6` | PMID: `41720999`. Demonstra recuperação de sinal de aminoácidos intracristalinos em dentes fósseis de até **48 milhões de anos** (Eoceno de Messel), estendendo o alcance tafonômico do esmalte muito além do limite do colágeno ósseo.
   - **Patramanis et al. (2026)** *Phylogenetic inference using ancient proteins.* *Genome Biology and Evolution* — DOI: `10.1093/gbe/evag035` | PMID: `41744466`. Avalia sistematicamente o poder de resolução filogenética do conjunto canônico de **12 paleoproteínas** (`AMELX`, `AMELY`, `ENAM`, `AMBN`, `ODAM`, `MMP20`, `COL1A1`, `COL1A2`, `COL2A1`, `ALB`, `AHSG`, `SERPINC1`), que coincide diretamente com o módulo estrutural de esmalte de `docs/PROJETO-ARCHOSAURIA.md`.
2. **Modelos de IA para Sequenciamento *De Novo* de Peptídeos Antigos e Homologia Remota**:
   - **Casanovo (`v4.2+` / `v5.x`)**: Yilmaz et al. (2024) *Nature Communications* 15:6427 (DOI: `10.1038/s41467-024-49731-x`, PMID: `39080256`) e Straub et al. (2025/2026) *J. Proteome Res.* (DOI: `10.1021/acs.jproteome.5c00706`, PMID: `41468557`). *Transformer* encoder-decoder para MS/MS livre de banco de dados com suporte a PTMs diagenéticas ($\text{N/Q}+0.9840\text{ Da}$, $\text{P/M}+15.9949\text{ Da}$).
   - **InstaNovo / InstaNovo+**: Eloff, Kalogeropoulos et al. (2025) *Nature Machine Intelligence* 7:565–579 (DOI: `10.1038/s42256-025-01019-5`, PMID: `41847501`). Combina decodificação autorregressiva com **difusão discreta multinomial** (`InstaNovo+`), elevando o *recall* em peptídeos não-trípticos gerados por hidrólise diagenética espontânea.
   - **$\pi$-PrimeNovo**: Zhang et al. (2025) *Nature Communications* 16:267 (DOI: `10.1038/s41467-024-55021-3`, PMID: `39747823`). Arquitetura não-autorregressiva (CTC) com controle exato de massa (`PMC` $\le 0.01\text{ Da}$), resolvendo ambiguidades isobáricas ($\text{GG} = \text{N}$, $\text{GA}/\text{AG} = \text{Q}/\text{K}$).
   - **Modanovo**: Klaproth-Andrade et al. (2025) *Molecular & Cellular Proteomics* (DOI: `10.1016/j.mcpro.2025.101501`, PMID: `41453616`). Decodificação robusta sob múltiplas modificações pós-traducionais/diagenéticas.
   - **ProstT5 & Foldseek 3Di**: Heinzinger et al. (2024) *NAR Genomics and Bioinformatics* 6(4):lqae150 (DOI: `10.1093/nargab/lqae150`, PMID: `39633723`). Traduz sequências 1D diretamente no alfabeto estrutural 3Di para busca de homologia remota em linhagens extintas sem genoma próximo.
3. **Autenticidade Diagenética vs. Idade Térmica**:
   - **Ramsøe et al. (2020)** `deamiDATE` (*J. Archaeol. Sci.* 115:105080, DOI: `10.1016/j.jas.2020.105080`), **Schroeter & Cleland (2016)** (*Rapid Commun. Mass Spectrom.* 30:597–602, DOI: `10.1002/rcm.7475`, PMID: `26842581`) e **Demarchi et al. (2016)** (*eLife* 5:e17092, DOI: `10.7554/eLife.17092`, PMID: `27673330`). A desamidação de $\text{Asn}/\text{Gln}$ depende do resíduo $i+1$ ($\text{NG}/\text{QG}$ rápidos vs. $\text{NV}/\text{QP}$ lentos) e da ligação mineral; serve como **filtro de autenticidade endógena**, nunca como cronômetro absoluto isolado.
   - **Briggs et al. (2007)** (*PNAS* 104:14616–14621, DOI: `10.1073/pnas.0704665104`), **Jónsson et al. (2013)** `mapDamage2` (*Bioinformatics* 29:1682–1684, DOI: `10.1093/bioinformatics/btt193`) e **Neukamm et al. (2021)** `DamageProfiler` (*Bioinformatics* 37:3652–3653, DOI: `10.1093/bioinformatics/btab190`). Parametrização formal de $\delta_s$, $\delta_d$, $\lambda$ e $\nu$ em aDNA.
4. **Inferência Filogenética FBD, Viés de Aquisição Morfológico (`+ASC`) e Viés de Consenso em ASR (`AltAll`)**:
   - **Heath, Huelsenbeck & Stadler (2014)** (*PNAS* 111:E2957–E2966, DOI: `10.1073/pnas.1319091111`), **Gavryushkina et al. (2014)** (*PLOS Comput. Biol.* 10:e1003919) e **Höhna et al. (2016)** `RevBayes` (*Syst. Biol.* 65:726–736, DOI: `10.1093/sysbio/syw021`).
   - **Lewis (2001)** (*Syst. Biol.* 50:913–925, DOI: `10.1080/106351501753462876`, PMID: `12116640`). Modelo condicional `Mkv` (`L_i / (1 - L_const)`) obrigatório para matrizes morfológicas sem sítios constantes.
   - **Eick, Bridgham, Anderson, Harms & Thornton (2017)** *Robustness of Reconstructed Ancestral Protein Functions to Statistical Uncertainty.* *Molecular Biology and Evolution* 34(2):247–261 — DOI: `10.1093/molbev/msw223` | PMID: `27795231`. Define o protocolo **AltAll** para eliminar o viés de estabilidade térmica artificial de proteínas ancestrais ML.
5. **Segmentação Volumétrica 3D de Micro-CT Fóssil e Paleobiomecânica (FEA)**:
   - **`nnInteractive`**: Isensee et al. (2025) *arXiv:2503.08373* (DOI: `10.48550/arXiv.2503.08373`). Segmentação volumétrica 3D interativa por pontos/laços 3D, superior ao `TotalSegmentator` (treinado em humano) quando há matriz calcítica/pirítica aderida ao osso fóssil.
   - **`SlicerMorph` / `ALPACA` / `MALPACA`**: Rolfe et al. (2021) (*Methods Ecol. Evol.* 12:1816–1824, DOI: `10.1111/2041-210X.13669`), Porto et al. (2021) (DOI: `10.1111/2041-210X.13689`) e Zhang et al. (2022) `MALPACA` (*PLOS ONE* 17:e0278035, DOI: `10.1371/journal.pone.0278035`). Registro deformável *multi-template* para transferência de semilandmarks em clados com alta disparidade morfológica.
   - **Propriedades Elásticas Reais de Osso de Archosauria para FEA (`FEBio`)**:
     - *Alligator mississippiensis* (ortotropia cortical completa por ultrassom): **Zapata et al. (2010)** *Bone* 46(3):860–867 (DOI: `10.1016/j.bone.2009.11.010`, PMID: `19922820`) — $E_1 = 8.8\text{ GPa}$, $E_2 = 10.7\text{ GPa}$, $E_3 = 20.5\text{ GPa}$; $G_{12} = 3.35\text{ GPa}$, $G_{13} = 4.96\text{ GPa}$, $G_{23} = 5.92\text{ GPa}$; $\nu_{12} = 0.40$, $\nu_{13} = 0.14$, $\nu_{23} = 0.14$.
     - *Struthio camelus* (validação craniana *ex vivo*): **Cuff, Bright & Rayfield (2015)** *PeerJ* 3:e1294 (DOI: `10.7717/peerj.1294`, PMID: `26500813`) — $E_{\text{cortical}} = 13.5\text{–}17.3\text{ GPa}$, $\nu = 0.35$.

---

## 3. Ferramentas Instaladas e Verificadas no WSL (`Ubuntu-26.04` — Camada `pixi global`)

Onze novos ambientes foram adicionados a [`bootstrap/pixi-global.toml`](../bootstrap/pixi-global.toml)
(totalizando **54 ambientes** e **64 pacotes** travados com `sha256` em
[`bootstrap/pixi-global.snapshot.tsv`](../bootstrap/pixi-global.snapshot.tsv) e
auditados por [`tests/check_wsl_ferramentas.py`](../tests/check_wsl_ferramentas.py)):

| Ambiente `pixi` | Versão Fixada | Binários Expostos | Papel na `bio-stack` |
| :--- | :--- | :--- | :--- |
| ★ `revbayes` | `1.4.2` | `rb`, `rb-mpi` | Inferência Bayesiana sob processo **Fossilized Birth-Death (FBD)**, *total-evidence tip-dating* e modelos morfológicos particionados `Mkv+G`. |
| ★ `beast2` | `2.7.7` (+ `beagle-lib 4.2.0`) | `beast`, `beauti`, `treeannotator`, `loganalyser`, `logcombiner`, `densitree`, `packagemanager` | Relógio molecular relaxado, *Sampled Ancestors* e filodinâmica temporal. |
| ★ `mrbayes` | `3.2.9` | `mb`, `mb-mpi` | Inferência Bayesiana combinada (molecular + morfológica discreta `Standard` com `lset coding=variable`). |
| ★ `damageprofiler` | `1.1` | `damageprofiler` | Quantificação rápida de desaminação pós-morte 5' $\text{C}\to\text{T}$ e 3' $\text{G}\to\text{A}$ e distribuição de comprimentos em BAMs de aDNA. |
| ★ `mapdamage2` | `2.2.3` | `mapDamage` | Estimação MCMC dos parâmetros de Briggs-Johnson ($\delta_s, \delta_d, \lambda, \nu$) e reescalonamento Phred (`--rescale`) de bases deaminadas. |
| ★ `adapterremoval` | `3.0.2` | `adapterremoval`, `adapterremoval3` | Trimming e colapso (`--collapse`) de leituras *paired-end* sobrepostas em fragmentos curtos de aDNA ($30\text{–}70\text{ pb}$). |
| ★ `bwa` | `0.7.19` | `bwa` | Alinhamento de aDNA com semente desativada (`bwa aln -l 1024 -n 0.01 -o 2`) tolerante a dano terminal $\text{C}\to\text{T}$. |
| ★ `angsd` | `0.940` | `angsd`, `realSFS`, `thetaStat`, `NGSadmix`, `contamination`, `ibs`, `haploToPlink` | Genômica populacional em aDNA de baixa cobertura via verossimilhanças de genótipo (`-GL 1`/`-GL 2`) sem chamada dura de variantes. |
| ★ `smina` | `2020.12.10` | `smina` | **Upgrade de precisão estrutural**: docking físico clássico e minimização/rescoring local (`--score_only`, `--minimize`) para auditar poses geradas pelo modelo de difusão `diffdock` (NIM). |
| ★ `autodock-vina` | `1.1.2` | `vina`, `vina_split` | **Upgrade de precisão estrutural**: função de energia empírica/física de referência para cross-check de ligantes em bolsões ancestrais ou modernos. |
| ★ `vsearch` | `2.32.0` | `vsearch` | Detecção de quimeras *de novo*, desreplicação e agrupamento para **sedaDNA** (DNA antigo sedimentar) e metabarcoding. |

---

## 4. Novo Toolkit Python Auditável: `scripts/paleo/paleo_toolkit.py`

Script autocontido (integrado ao log de reprodutibilidade `runs/provenance.jsonl`
para geração automática de Métodos via `/metodos`) com quatro subcomandos:

1. **`pbdb-range` — Calibração Estratigráfica Direta na Paleobiology Database (PBDB API v1.2)**:
   ```powershell
   uv run python scripts/paleo/paleo_toolkit.py pbdb-range --taxa "Tyrannosaurus rex,Archaeopteryx lithographica,Spinosaurus aegyptiacus,Alligator mississippiensis"
   ```
   Retorna `taxon_no`, status (`extinct`/`extant`), contagem de ocorrências,
   intervalos $FAD = [FAD_{\min}, FAD_{\max}]$ e $LAD = [LAD_{\min}, LAD_{\max}]$
   em Ma, e escreve os comandos prontos para `RevBayes` (`dnUniform`),
   `MCMCTree` e `IQ-TREE 3 LSD2`.
2. **`mk-audit` — Auditor de Viés de Aquisição de Lewis (2001) para Matrizes Morfológicas**:
   ```powershell
   uv run python scripts/paleo/paleo_toolkit.py mk-audit --input runs/paleo/archosauria_morpho_demo.fa
   ```
   Conta caracteres invariáveis, autapomórficos e parcimônia-informativos,
   calcula a porcentagem de dados ausentes (`?`/`-`) por táxon fóssil e emite a
   parametrização exata (`MK+FQ+ASC+G4`, `coding="variable"`).
3. **`altall-asr` — Gerador do Par Ancestral `ML` + `AltAll` (Eick et al. 2017)**:
   ```powershell
   uv run python scripts/paleo/paleo_toolkit.py altall-asr --state-file runs/paleo/amelx_asr.state --node Node2 --pp-threshold 0.80 --out-fasta runs/paleo/amelx_ancestral_pair.fa
   ```
   Lê a tabela `.state` produzida por `iqtree3 --ancestral` e exporta um FASTA
   contendo a sequência ancestral de máxima verossimilhança (`ML`) e a sequência
   `AltAll` (que substitui todos os sítios com $PP_1 < 0.80$ pelo segundo
   aminoácido mais provável), pronta para dobramento comparativo no `Boltz-2` /
   `OpenFold3`.
4. **`zooms-qc` — Digestão Tríptica *In Silico* e Relógio Diagenético de Desamidação**:
   ```powershell
   uv run python scripts/paleo/paleo_toolkit.py zooms-qc --sequence "MPLPPHPGHPGYINFSYEVLTPLK..."
   ```
   Aplica a regra de Keil (clivagem após K/R exceto antes de P) e calcula as
   massas monoisotópicas exatas $[M+H]^+$ e $[M+2H]^{2+}$ intactas e após
   desamidação total de $\text{Asn}/\text{Gln}$ ($+0.984016\text{ Da}$ por sítio).

---

## 5. Novos Agentes, Skills e Referências Adicionados

- **Agente `paleo-analyst`**:
  - `.agents/agents/paleo-analyst.md` (Antigravity)
  - `.agents/plugins/bio-stack/agents/paleo-analyst.md` (Plugin)
  - `.claude/agents/paleo-analyst.md` (Claude Code)
  - `.codex/agents/paleo-analyst.toml` (Codex)
- **Comando / Skill `/paleo`**:
  - `.agents/skills/paleo/SKILL.md`
  - `.agents/plugins/bio-stack/skills/paleo/SKILL.md`
  - `.agents/skills/source-command-paleo/SKILL.md`
  - `.claude/commands/paleo.md`
- **Contrato de API da Paleobiology Database (PBDB v1.2) + Open Tree of Life (OTOL v3)**:
  - `.agents/skills/database-lookup/references/pbdb.md` (registrado em `.agents/skills/database-lookup/SKILL.md`).
- **Correção estrutural nos subagentes do Antigravity**:
  - Removido o campo inválido `tools: [..., call_mcp_tool, ...]` do frontmatter de todos os agentes em `.agents/agents/*.md` e `.agents/plugins/bio-stack/agents/*.md`, que causava falha na instanciação de subagentes (`unknown component: tool "call_mcp_tool" not found in registry`).
