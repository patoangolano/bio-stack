# Braço paleontológico computacional: fósseis-chave e convergência morfológica dos rapinantes

2026-10-17/18. Tudo reproduzível a partir de `scripts/rapinantes/fosseis/`
(ordem: `p01` → `f02` → `p01 mapear` → `p02` → `f03` → `p03` → `p04 gerar/gerar_atuais` → `f04`/`f05`
→ `f06`/`f07`/`p04 lz07_spr`/`f08` → `p06` (R) → `p07` → `p08` (R) → `p09` → `p10`). Saídas em
`runs/rapinantes/fosseis/` (figuras em `figuras/`: `arvore_<ID>_{livre,esqueleto}.png`,
`au_fosseis.png`, `homoplasia_classes.png`, `morfoespaco_<ID>.png`); tabelas compiladas em `runs/rapinantes/fosseis/tabelas_resumo.md`;
proveniência dos dados em `runs/rapinantes/fosseis/LEIAME_dados_brutos.md` (o hook `guard_data`
bloqueou escrever o LEIAME dentro de `data/raw/`; os arquivos baixados estão em
`data/raw/rapinantes/morfologia/`, sha256 em `runs/rapinantes/fosseis/log/sha256_brutos.txt`).

## 1. Matrizes usadas (DOI conferido na Crossref/DataCite)

| ID | Fonte | Táxons × caracteres usados | Rapinantes atuais | Fósseis de interesse | Lista de caracteres |
|---|---|---|---|---|---|
| K19 | Ksepka, Grande & Mayr 2019, Curr Biol 29:657, 10.1016/j.cub.2018.12.040 (NEXUS via G. T. Lloyd) | 62 × 146 (34 atuais) | *Falco, Buteo, Sagittarius, Cathartes, Tyto, Bubo* (+ *Cariama*) | Halcyornithidae (*Pseudasturides, Cyrilavis*) + 26 Psittacopasserae/Coraciimorphae paleógenos | não obtida |
| M11 | Mayr 2011, J Syst Palaeontol 9:159, 10.1080/14772019.2010.505252 (Lloyd) | 33 × 100 (famílias) | Accipitridae, Sagittariidae, Cathartidae, Falconidae, Strigiformes | Messelasturidae, Halcyornithidae, Sandcoleidae | não obtida |
| M05 | Mayr 2005, J Vertebr Paleontol 25:635, 10.1671/0272-4634(2005)025[0635:tpoapp]2.0.co;2 (Lloyd) | 23 × 110 | idem + Cariamidae | *Messelastur, Tynskya*, Pseudasturidae | não obtida |
| VZ26 | Mayr 2026, Vertebr Zool 76:121, 10.3897/vz.76.e182505 (versão emendada de Mayr & Kitchener 2023) | 42 × 105 | só Falconidae e Strigiformes (sem Accipitridae, sem Cariamidae) | *Messelastur, Tynskya*, 6 halcionitídeos, zigodactilídeos | sim (File S1) |
| MK22 | Mayr & Kitchener 2022, J Vertebr Paleontol 41:e2083515, 10.1080/02724634.2021.2083515 (figshare 10.6084/m9.figshare.20267628.v2) | 10 × 42 (3 atuais) | Falconidae (+ *Cariama*) | ***Masillaraptor*, *Danielsraptor***, *Bathornis*, *Dynamopterus*, 3 fororracídeos | sim |
| M16 | Mayr 2016, J Paleontol 90:357, 10.1017/jpa.2016.45 (Lloyd) | 13 × 40 (5 atuais) | Falconidae, Cathartidae (+ *Cariama*) | *Bathornis, Paracrax*, Idiornithidae, fororracídeos | não obtida |
| D21 | Degrange et al. 2021, Pap Palaeontol, 10.1002/spp2.1361 (MorphoBank P3896, CC0) | 72 × 204 | 27 Accipitridae, 8 Falconidae, 7 Cathartidae, *Pandion, Sagittarius* | *Dryornis pampeanus*, *Pampagyps imperator* (Cathartidae) | sim (em português) |
| MC03 | Mayr & Clarke 2003, Cladistics 19:527, 10.1016/j.cladistics.2003.10.003 (Lloyd) | 46 × 148 (famílias de Neornithes) | Accipitridae, Falconidae, Sagittariidae, Cathartidae, Strigiformes | só grupos externos mesozoicos | não obtida |
| C25 | Chen et al. 2025, Integr Org Biol 7, 10.1093/iob/obaf029 | 77 × 202, cintura escapular e asa, codificada em boa parte em TC | *Pandion, Elanus, Coragyps, Micrastur, Tyto, Ninox* | — | sim |
| LZ07 | Livezey & Zusi 2007, Zool J Linn Soc 149:1, 10.1111/j.1096-3642.2006.00293.x (Lloyd = MorphoBank P510, 100% das células iguais) | 154 Neornithes × 2325 caracteres variáveis (de 2954) | *Accipiter, Gyps, Pandion, Falco, Polyborus, Sagittarius, Cathartes*, 4 corujas | nenhum rapinante fóssil | sim (MorphoBank) |

Não existe matriz baixável com **corujas-tronco** (*Ogygoptynx, Berruornis*, Protostrigidae,
*Primoptynx, Sophiornis, Ypresiglaux*), ***Horusornis*** ou **Pandionidae fósseis** (busca na coleção
de Lloyd, na listagem completa da API pública do MorphoBank, figshare, Europe PMC e web).
***Danielsavis* não é rapinante** (Danielsavidae, Galloanseres); o masilarraptorídeo de
Walton-on-the-Naze é ***Danielsraptor***. *Masillaraptor buchheimi* (Li et al. 2025,
10.7302/27584) e as listas de caracteres de K19/M11/M05/MC03 ficaram atrás de paywall ou anti-robô.

## 2. Métodos (resumo; comandos exatos em `runs/rapinantes/fosseis/log/comandos_*.txt`)

- **Codificação**: polimorfismo e “−” viram “?”; estados de cada caráter renumerados 0..k−1;
  caracteres invariantes após isso saem; ordenados pelos autores entram não ordenados (decisões
  por matriz em `matrizes/<ID>/decisoes.json`).
- **ML**: IQ-TREE 3.1.4 (pixi, WSL), partições por número de estados (k = 2, 3, ≥4), cada uma
  **MK+ASC+G4** (`-p`, comprimentos proporcionais); 1000 UFBoot com `--bnni` e 1000 SH-aLRT;
  `--seed 20261017`; ≤ 4 threads. LZ07: UFBoot interrompido (tempo) e análise com esqueleto
  sem UFBoot; MK22: com 3 atuais não há esqueleto a impor.
- **Esqueleto molecular**: árvore principal de **Stiller et al. 2024** (Nature 629,
  10.1038/s41586-024-07323-1; `63K.tre`, ASTRAL em 63 430 regiões intergênicas; repositório
  10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a), podada a um representante por terminal
  (família pelo NCBI Taxonomy local via `taxonkit`; escolhas manuais justificadas em
  `p01_taxonomia.py`). Terminais que dividem representante viram politomia. Restrição `-g` só
  nos atuais; fósseis livres. Em Telluraves essa árvore dá ((Accipitriformes + Strigiformes) +
  Coraciimorphae) e (Cariamiformes + (Falconiformes + Psittacopasserae)).
- **Posição de fóssil**: em cada uma das 1000 árvores UFBoot, os outros fósseis são podados e
  registra-se o conjunto de atuais irmão do fóssil; “grupo total de X” = fração em que esse
  conjunto está contido em X (fóssil no tronco ou dentro de X).
- **AU**: para cada hipótese “fóssil irmão de X”, busca ML restrita (esqueleto + enxerto);
  as melhores árvores, mais a ML com fósseis livres, entram em `-z … -zb 10000 -zw -au`.
  Para os atuais, o mesmo com o Falconidae deslocado para junto de Accipitridae (a “Falconiformes”
  clássica) e, em D21/LZ07, Cathartidae junto de Ciconiidae; o teste final é só entre variantes
  do esqueleto (a árvore morfológica livre domina qualquer conjunto e rejeita todas).
- **Homoplasia**: IC e IR de Fitch por caráter (phangorn 2.12.1, R 4.6.1) na árvore ML
  morfológica × na ML com esqueleto, com e sem fósseis; classe “rapinante” pela descrição do
  caráter (gancho/cera/tômio; unguais, hálux, dedos do pé; tarsometatarso/hipotarso);
  permutação de rótulos, 10 000, `set.seed(20261017)`, estratificada por matriz no agregado.
  Passos por caráter na árvore molecular × na alternativa falcão-gavião.
- **Morfoespaço**: MORD (Lloyd 2016) e Gower só em caracteres comparáveis (≥ 15 por par;
  táxon com mais pares faltantes sai até a matriz ficar completa), PCoA; nos atuais de
  Telluraves, PERMANOVA (clado molecular; guilda rapinante = Accipitridae, Pandionidae,
  Sagittariidae, Falconidae, Strigiformes; termo sequencial guilda | clado), Mantel parcial e
  distância média entre rapinantes de clados diferentes contra sorteio do mesmo número de táxons
  (9999 permutações, `np.random.default_rng(20261017)`).

## 3. Resultados

### 3a. Morfologia sozinha: que grupos “tradicionais” aparecem

UFBoot (%) do grupo nas 1000 árvores da análise livre (podadas aos atuais); † = está na árvore ML.
LZ07: só presença na árvore ML — o refinamento das UFBoot (154 táxons × 2325 caracteres) foi
interrompido em 100/1000 réplicas por tempo; a busca ML estava concluída (lnL ≈ −103 015).

| grupo | K19 | M11 | M05 | MC03 | D21 | C25 | VZ26 | LZ07 |
|---|---|---|---|---|---|---|---|---|
| Accipitridae + Falconidae | **81,7** † | 26,5 | **92,1** † | **81,5** † | 1,3 | 0,8 | — | não (*Pandion* entra junto de Falconidae) |
| Falconiformes s.l. (+ Sagittariidae, Pandionidae, Cathartidae) | 19,4 | 0,1 | 19,4 | 3,4 | 65,3 † (= todos os rapinantes da matriz) | 0,0 | — | sim † |
| “Raptores” (falconiformes s.l. + corujas) | **78,7** † | 39,8 † | **87,1** † | **89,8** † | — | 0,0 | 86,6 † (só Falconidae + corujas) | sim † |
| Strigiformes + Caprimulgiformes tradicionais | — | 0,0 | 0,0 | 0,0 | — | 0,7 | — | não |
| Cathartidae + Ciconiidae | — | — | 7,8 | 0,0 | 0,4 | 0,6 | — | não |
| Accipitriformes (molecular, com Cathartidae) | 0,0 | 0,0 | 0,0 | 0,1 | 1,5 | 12,8 | — | não |
| Hieraves (Accipitriformes + Strigiformes) | 0,0 | 0,0 | 0,0 | 0,3 | — | 0,3 | — | não |
| Eufalconimorphae / Australaves | 0,0 / 0,0 | 0,0 / — | — | 0,0 / 0,0 | — | 0,0 / 0,0 | 0,0 / — | não / não |

A ML morfológica recupera só 18–42% dos *splits* do esqueleto molecular (K19 13/31, M11 6/27,
M05 3/17, MC03 8/40, C25 14/72, LZ07 55/134; D21 14/21 porque o esqueleto dela é quase todo
politomia). Na ML de LZ07 a ordem é (((Accipitridae, (*Pandion*, Falconidae)), *Sagittarius*),
*Cathartes*) + corujas, com *Cariama* fora e os caprimulgiformes junto de passeriformes e
*Colius* — nenhuma matriz junta corujas e Caprimulgiformes.
Nenhuma matriz recupera Accipitriformes molecular, Hieraves, Australaves ou Eufalconimorphae
(≤ 1,5%, exceto C25: Accipitriformes 12,8%).

### 3b. Fósseis sobre o esqueleto molecular (UFBoot da posição)

| fóssil (matriz) | irmão na ML | posições mais frequentes nas UFBoot | grupo total |
|---|---|---|---|
| *Masillaraptor* e *Danielsraptor* (MK22; com 3 atuais o esqueleto não restringe nada) | Falconidae | irmãos de Falconidae 77,8%; fora de Falconidae + *Cariama* 17,9%; com *Cariama* 4,3% | pan-Falconiformes 77,8% |
| Messelasturidae (M11) | Falconidae | Falconidae 36,2; Accipitridae 25,3; Coliidae 22,0 | pan-Eufalconimorphae 39,6; pan-Hieraves 37,1; pan-Psittaciformes 1,8 |
| *Messelastur* (M05) | Psittacidae | Psittacidae 66,5; Falconidae 31,3 | pan-Psittaciformes 66,5; pan-Falconiformes 31,3 |
| *Messelastur* (VZ26) | tronco de Psittacopasserae | tronco de Psittacopasserae 41,8; Passeriformes 30,0; Falconidae 25,7 | pan-Psittacopasserae 73,3; pan-Eufalconimorphae 99,0 |
| Halcyornithidae (K19: *Pseudasturides, Cyrilavis*) | Passeriformes | Psittaciformes 41,9; Passeriformes 33,6; tronco 23,7 | pan-Psittacopasserae 99,2 |
| Halcyornithidae (M11) | Coliidae | Coliidae 36,9; Psittacidae 28,9; Falconidae 21,3 | pan-Eufalconimorphae 58,8 |
| Pseudasturidae (M05) | Psittacidae | Psittacidae 92,6 | pan-Psittaciformes 92,6 |
| fororracídeos (M16, MK22) | *Cariama* | *Cariama* 100 | pan-Cariamiformes 100 |
| *Bathornis* (M16) | Opisthocomidae | Opisthocomidae 44,5; *Cariama* 42,5 | pan-Cariamiformes 42,5 (livre: *Cariama* 97%) |
| *Paracrax*, *Elaphrocnemus* (M16) | Opisthocomidae | Opisthocomidae ~75; *Cariama* 13 | instável (matriz de 13 táxons) |
| *Dryornis*, *Pampagyps* (D21) | dentro de Cathartidae (perto dos condores) | — | pan-Cathartidae 100 |

Sem esqueleto, a morfologia puxa os messelasturídeos para junto dos rapinantes: irmão das
corujas em M05 (68%), dentro de (Accipitridae, Falconidae, Strigiformes) em M11 (63%). Com o
esqueleto, essa atração se desfaz e eles caem em Eufalconimorphae (Falconidae, Psittacopasserae
ou o tronco deles), exceto em M11, que divide o sinal entre Falconidae e Accipitridae.

### 3c. Testes AU (p-AU; rejeição a 5% em negrito)

- **Masillaraptoridae (MK22)**: melhor = tronco de Falconidae (p = 1); com *Cariama* /
  pan-Cariamiformes p = 0,10; tronco de (Falconidae + Cariamiformes) p = 0,25. Nada rejeitado.
  **Tronco de Accipitriformes e de Psittacopasserae não são testáveis**: a matriz não tem esses
  táxons.
- **Messelasturidae**: M11 — melhor irmão de Falconidae; Accipitridae 0,58, Psittaciformes 0,40,
  Accipitriformes 0,36, Psittacopasserae 0,35, Strigiformes 0,14, Passeriformes 0,086: nada
  rejeitado. M05 (*Messelastur*) — melhor Psittaciformes; corujas 0,54, Falconidae 0,41,
  Hieraves 0,28, Accipitridae 0,19; **tronco de Australaves 0,037**; Accipitriformes 0,050.
  VZ26 — melhor tronco de Psittacopasserae; nada rejeitado (mín. 0,07).
- **Halcyornithidae**: K19 — Passeriformes, Psittaciformes e Psittacopasserae empatados
  (p = 1; 0,54; 0,53); **Accipitriformes 0,007, Strigiformes 0,008, Hieraves 0,048,
  Australaves 0,044**; Falconiformes 0,16 não rejeitado. M05 (Pseudasturidae): **Accipitridae
  0,012, Accipitriformes 0,022, Cariamidae 0,027**.
- **Cathartidae fósseis (D21)**: *Dryornis* dentro de Cathartidae p = 0,58; **Falconidae
  8×10⁻⁶, Accipitridae 2×10⁻⁵, Ciconiidae 1×10⁻⁴, Pandionidae 0,032**. *Pampagyps*: **Accipitridae
  0,021, Pandionidae 0,015**; Ciconiidae 0,058.
- ***Bathornis*/*Elaphrocnemus* (M16)**: nada rejeitado (13 táxons × 40 caracteres).
- **Atuais — Falconidae junto de Accipitridae contra a posição molecular**: a morfologia
  prefere a “Falconiformes” clássica em K19 (ΔlnL 45,3; posição molecular **p = 0,004**),
  M11 (63,6; **0,0017**), M05 (48,8; **0,0001**), MC03 (60,1; **0,0003**); não em D21 (5,3;
  p = 0,43) nem em C25, a matriz de cintura e asa codificada em TC (2,9; p = 0,46). Em LZ07
  (alternativa feita por SPR sobre a árvore restrita, sem reotimizar — o que joga contra ela) o
  ganho é de 595 lnL e a posição molecular cai com **p = 1×10⁻³⁹**.
  Cathartidae com Ciconiidae: D21 +1,0 lnL (p = 0,41); LZ07 −908 lnL (**p = 8×10⁻⁵⁶** contra a
  hipótese antiga urubu-cegonha).

### 3d. Homoplasia

- IC do conjunto cai ao trocar a árvore morfológica pela molecular-restrita: K19 0,286 → 0,271,
  M11 0,335 → 0,285, M05 0,367 → 0,313, MC03 0,228 → 0,197, D21 0,285 → 0,278, C25 0,099 → 0,089,
  VZ26 0,328 → 0,315; LZ07 0,217 → 0,189 (só atuais; 15 152 → 17 354 passos). São de 3% (D21) a
  18% (M11) de passos a mais.
- Caracteres informativos que pagam menos passos com Falconidae junto de Accipitridae (contra os
  que pagam mais), só atuais: K19 17 × 3, M11 18 × 2, M05 18 × 3, MC03 18 × 3; D21 25 × 19;
  C25 21 × 21 (saldo 0); **LZ07 161 × 30** (saldo de 131 passos). Em C25 os 21 são todos de esterno, coracoide, úmero e
  carpometacarpo — a asa não carrega o sinal da união.
- Classe “rapinante” (bico, garra/pé, tarsometatarso) contra o resto: em VZ26 e D21 os
  caracteres rapinantes têm IC **maior**, não menor, nas duas árvores (p unilateral “menor”
  0,97–0,99); a perda de IC ao passar para a árvore molecular é um pouco maior neles (VZ26
  ΔIC −0,044 × −0,012; p = 0,11). Essas duas matrizes quase não têm Accipitridae (VZ26) ou têm
  só 4 caracteres rapinantes informativos (D21).
- **LZ07 é a matriz que responde** (lista de caracteres + falcões, gaviões, *Pandion*,
  *Sagittarius*, urubu e corujas; só atuais): os 221 caracteres rapinantes informativos têm IC
  menor na árvore molecular (0,426 × 0,483; **p = 0,013**) e perdem mais IC ao sair da árvore
  morfológica (ΔIC −0,090 × −0,053; **p = 0,0008**; ΔIR **p = 0,003**); na árvore morfológica a
  diferença some (p = 0,22) — é a topologia molecular que os torna homoplásticos. Por classe:
  pé/garra (90 caracteres) ΔIC −0,104, tarsometatarso (113) −0,089, bico (24) −0,044, resto (1739)
  −0,053. Entre os 161 caracteres que favorecem Falconidae + Accipitridae, os rapinantes estão
  super-representados (29 de 227 × 132 de 1739; OR = 1,78, **Fisher p = 0,010**); pé/garra dá 19
  dos 29: proporção das falanges dos dedos II e IV, unguais, hálux, tendões flexores
  (*m. flexor hallucis longus*), canais do hipotarso, almofadas plantares; no bico, cera e
  serrilha tomial. Agregado VZ26 + D21 + LZ07 (permutação estratificada): ΔIC p = 0,0007.

### 3e. Morfoespaço (MORD, atuais de Telluraves)

| matriz | n (rapinantes) | PERMANOVA clado R² | guilda R² | guilda \| clado | Mantel guilda \| clado | MORD rapinantes de clados diferentes × outros pares entre clados (p) |
|---|---|---|---|---|---|---|
| K19 | 32 (5) | 0,61 (p = 1e-4) | 0,17 (1e-4) | 0,009 (0,80) | r = 0,11 (0,068) | 0,195 × 0,344 (1e-4) |
| M11 | 18 (4) | 0,59 (3e-4) | 0,36 (4e-4) | 0,018 (0,80) | r = 0,42 (0,0025) | 0,148 × 0,310 (0,004) |
| M05 | 10 (4) | 0,76 (0,005) | 0,26 (0,04) | 0,045 (0,67) | r = 0,22 (0,091) | 0,206 × 0,346 (0,014) |
| MC03 | 11 (4) | 0,75 (0,02) | 0,26 (0,03) | 0,047 (0,67) | r = 0,24 (0,074) | 0,126 × 0,219 (0,006) |
| C25 | 25 (5) | 0,49 (1e-4) | 0,13 (0,001) | 0,023 (0,69) | r = 0,07 (0,16) | 0,252 × 0,325 (0,002) |
| D21 | 44 (37) | 0,21 (1e-4) | 0,52 (1e-4) | 0,51 (1e-4)† | r = 0,91 (1e-4)† | 0,270 × 0,445 (1e-4) |
| VZ26 | 15 (2) | 0,43 (0,008) | 0,19 (0,009) | — | r = 0,02 (0,44) | 0,152 × 0,316 (0,022) |
| LZ07 | 47 (10) | 0,70 (1e-4) | 0,28 (1e-4) | 0,017 (0,019) | r = 0,53 (1e-4) | 0,122 × 0,142 (0,038) |

† Em D21 “guilda” separa Cathartidae do resto de Accipitriformes; o número que importa é a
última coluna: falcões ficam mais perto dos gaviões (0,27) do que dos urubus (0,45). Em todas
as matrizes os rapinantes de clados diferentes ficam mais próximos entre si do que pares
sorteados de clados diferentes (p de 1e-4 a 0,038), mas a guilda quase não explica variância
além do clado na PERMANOVA (R² ≤ 0,05; só LZ07 chega a p = 0,019 com R² = 0,017), porque guilda
e clado são quase colineares (Falconidae e Strigiformes são inteiramente rapinantes). O Mantel
parcial guilda | clado é significativo em M11 (r = 0,42) e LZ07 (r = 0,53).

## 4. Micro-TC (MorphoSource, API pública, 2026-10-17)

Inventário completo em `runs/rapinantes/fosseis/microct/inventario_morphosource.tsv` (56 consultas:
48 nomes de fósseis, 8 famílias atuais; 604 mídias; tamanho em MB; URL
`https://www.morphosource.org/concern/media/<id>`; resumo por consulta em `resumo_morphosource.tsv`).

- **Nenhuma mídia** (aberta ou restrita) de *Masillaraptor*, *Danielsraptor*, *Messelastur*,
  *Tynskya*, halcionitídeos, corujas paleogênicas, *Bathornis*, Idiornithidae, *Horusornis* ou
  *Pelargopappus*.
- Fororracídeos e batornitídeo **só restritos**: *Titanis walleri* (UF:VP, 9 mídias incl. 3 séries
  de TC; 123–3921 MB por mídia), *Andalgalornis* (FMNH:PA:14357, malha de luz estruturada, CC BY-NC), *Paracrax*
  sp. (AMNH:FARB:42998, 6 séries de TC, InC).
- Abertos de interesse: Teratornithidae (*Teratornis merriami* LACM, 4 malhas 7–35 MB, InC-NC;
  *Cathartornis gracilis* UCMP:V:12598, série de TC do tarsometatarso, 18 GB, CNE), *Breagyps clarki*
  (LACM/UCMP, 3 malhas + 1 série de TC de 7,7 GB), *Gymnogyps amplus* (LACM, 2 crânios em malha);
  *Argentavis* (endocasto) é restrito.
- Atuais abertos: Accipitridae 113 de 179 mídias (44 volumes, 69 malhas; esqueletos inteiros
  NHMUK de *Buteo buteo*, *Nisaetus*, *Torgos*, *Parabuteo*, *Accipiter nisus*; tarsometatarsos
  UCMP de *Aquila*/*Buteo*), Falconidae 24/83, Pandionidae 13/32 (esqueleto NHMUK), Sagittariidae
  15/15 (esqueleto NHMUK:ZOO:1898.5.7.1), Cathartidae 26/61 (esqueleto de *Vultur*), Strigidae 56/128,
  Tytonidae 25/43, Cariamidae 8/40.
- A demonstração (curvatura de ungual) não foi feita: não há fóssil relevante aberto, e os atuais
  abertos são volumes de esqueleto inteiro (GB) sem ungual segmentado.

## 5. Limitações e o que estes dados não dizem

- **Não posicionam corujas-tronco nem *Horusornis*/Pandionidae fósseis**: não há matriz.
- **Não testam *Masillaraptor* como tronco de Accipitriformes ou de Psittacopasserae**: a única
  matriz com masilarraptorídeos (MK22) tem 3 atuais. O apoio a “tronco de Falconidae” (78% UFBoot)
  não exclui Cariamiformes (p-AU 0,10).
- AU que não rejeita não é apoio: M11, VZ26, MK22 e M16 têm de 10 a 42 táxons e de 40 a 105
  caracteres; quase nenhuma posição de messelasturídeo é rejeitada.
- As posições dependem da árvore de Stiller et al. 2024 (ASTRAL 63K). O nó Hieraves × resto de
  Afroaves muda entre análises do próprio artigo; com outro esqueleto, “irmão de Accipitridae” e
  “irmão de Hieraves” podem trocar de lugar.
- Matrizes do mesmo grupo de autores (Mayr: M05, M11, M16, MK22, VZ26; Ksepka: K19) podem
  compartilhar caracteres (não conferido caráter a caráter); os resultados não somam evidência
  independente.
- Mk com polimorfismo como “?”, caracteres ordenados tratados como não ordenados, mesma árvore
  para caracteres de esqueleto, plumagem e músculo.
- A classe “rapinante” sai de palavras-chave nas descrições; só VZ26, D21, C25 e LZ07 têm
  descrição. K19, M11, M05 e MC03 — as que mais unem falcão e gavião — não têm lista de
  caracteres acessível, então não dá para dizer **quais** caracteres carregam a convergência nelas.
- PERMANOVA e Mantel tratam táxons como independentes; não há correção filogenética.
- O inventário do MorphoSource é por nome; mídia sem taxonomia preenchida não aparece.

## 6. Conclusão

1. **A morfologia une os rapinantes por convergência, e o sinal está no pé.** Sozinha, ela monta
   “Raptores” (falcões + gaviões + corujas) em K19, M05, MC03, M11 e LZ07 e Accipitridae +
   Falconidae com 82–92% de UFBoot; nenhuma matriz recupera Accipitriformes molecular, Hieraves,
   Australaves ou corujas + Caprimulgiformes. Sobre o esqueleto molecular, cinco das sete matrizes
   testáveis rejeitam a posição molecular de Falconidae a favor da “Falconiformes” clássica
   (p-AU de 0,004 a 1×10⁻³⁹); em LZ07 os caracteres de pé/garra, tarsometatarso e bico perdem mais
   IC na árvore molecular (p = 0,0008) e estão super-representados entre os que unem falcão e
   gavião (p = 0,010). Cintura escapular e asa (C25, codificada em TC) não carregam o sinal.
2. ***Masillaraptor* + *Danielsraptor***: tronco de Falconidae (UFBoot 78%), sem rejeitar
   Cariamiformes (p = 0,10); não há matriz pública que permita testar tronco de Accipitriformes ou
   de Psittacopasserae.
3. **Messelasturidae**: o que os punha junto das corujas e gaviões é a mesma convergência — com o
   esqueleto molecular caem em Eufalconimorphae em M05 (pan-Psittaciformes 66%) e VZ26 (tronco de
   Psittacopasserae 42%, pan-Psittacopasserae 73%); M11 divide entre Falconidae (36%) e
   Accipitridae (25%). O AU não rejeita quase nada. Halcionitídeos: pan-Psittacopasserae 99% (K19), com Accipitriformes e
   Strigiformes rejeitados (p ≤ 0,008).
4. **Cariamiformes e Cathartidae fósseis**: fororracídeos com *Cariama* (100%); *Bathornis* e
   idiornitídeos instáveis (matriz de 13 táxons); *Dryornis* e *Pampagyps* são Cathartidae de
   coroa (100%), com falcões, gaviões e cegonhas rejeitados para *Dryornis* (p ≤ 1×10⁻⁴).
5. **Corujas-tronco, *Horusornis* e Pandionidae fósseis ficam sem resposta**: não há matriz, e o
   MorphoSource não tem TC aberta de nenhum dos fósseis-chave (só de teratornitídeos e
   catartídeos pleistocênicos, e restrita para *Titanis* e *Paracrax*).
