# Braço computacional — filogenômica com conjunto de dados novo (nota interina)

Escrita em 2026-10-18, 02:20 (horário de Brasília), com a execução do pipeline
`20261018-020210-rapinantes_AA_mafft_clipkit` em andamento. **Nenhum resultado
filogenético está nesta versão**: as árvores de genes, o wASTRAL, o gCF/sCFL e a
supermatriz ainda não terminaram (seção 4). Os resultados entram depois, em
seções datadas próprias.

Tabelas desta etapa: `runs/rapinantes/filogenomica/` (`taxons.tsv`, `genes.tsv`,
`ocupacao_taxon.tsv`, `auditoria_refseq.tsv`, `auditoria_miniprot.tsv`,
`auditoria_miniprot_resumo.tsv`, `validacao_miniprot_vs_anotacao_GB.tsv`).
Scripts: `scripts/rapinantes/filogenomica/f01*`–`f08*`. Locos de entrada (AA):
`data/proc/rapinantes/filogenomica/locos_aa/` (300 FASTA, sha256 do conjunto
`b9dd5609247645ed…`); CDS correspondentes no WSL, `~/rapinantes/locos/cds_final/`.

## 1. Conjunto de dados

**57 táxons, 300 locos nucleares codificantes, 95,0 % das células táxon × loco
preenchidas** (média 54,2 táxons por loco; mínimo 49 = 85 % de 57).

### 1.1 Montagens

- **46 táxons RefSeq** (anotação NCBI Eukaryotic Genome Annotation; acessões
  `GCF_` em `taxons.tsv`). Um único caso trocado: *Phaethon aethereus*
  (GCF_964289735.1, anotação de 2025) não tem dados no NCBI Orthologs; entrou
  *P. lepturus* (GCF_000687285.1). Todos os outros 45 foram conferidos com um
  gene-sonda (BSG) antes do download (`ortologos_disponiveis.tsv`).
- **11 táxons de montagens GenBank sem anotação RefSeq**, com CDS extraído por
  `miniprot` 0.18 (proteína do parente RefSeq mais próximo como consulta):
  *Sarcoramphus papa* (GCA_037962945.1, consulta *Gymnogyps*), *Sagittarius
  serpentarius* (GCA_034782735.1), *Pandion haliaetus* (GCA_026930845.1),
  *Elanus axillaris* (GCA_026212715.1), *Gypaetus barbatus* (GCA_028022735.1),
  *Leptodon cayanensis* (GCA_055772235.1), *Aegypius monachus*
  (GCA_982185325.1), *Pithecophaga jefferyi* (GCA_054553165.1) — os oito com
  consulta *Aquila chrysaetos* —, *Herpetotheres cachinnans* (GCA_013399355.1)
  e *Ibycter americanus* (GCA_035590275.1), consulta *Falco peregrinus*, e
  *Ninox strenua* (GCA_982290615.1), consulta *Strix aluco*. Download: 11
  montagens, 4,1 GB comprimidos, em `~/rapinantes/raw` (teto pedido: ~10
  genomas, ≤ 15 GB).

### 1.2 Seleção dos locos

1. Relatório NCBI Gene do *Gallus gallus* (28 263 genes) → **12 677 candidatos**:
   codificante, com grupo "NCBI Ortholog", símbolo oficial (fora `LOC`),
   cromossomo nomeado (fora W, MT, Un) e fora de famílias multigênicas ou com
   muitas retrocópias (olfativos, MHC, imunoglobulinas/TCR, histonas,
   queratinas, ZNF/TRIM, defensinas, quimiocinas, CYP/UGT, proteínas
   ribossomais etc.; regex em `f02d_candidatos.py`).
2. Permutação fixa dos candidatos com `rng = np.random.default_rng(0)`; os
   primeiros 1 200 foram consultados no NCBI Orthologs (`datasets download gene
   --ortholog <46 tax IDs> --include cds,protein,product-report`, datasets
   18.38). **1 180 baixados**; 20 (um lote, `lote_008_p0`) falharam quatro
   vezes com erro HTTP/2 do servidor e saíram por motivo técnico, não
   biológico.
3. Por gene e táxon: um registro, o CDS mais longo (os relatórios de produto
   das aves não trazem RefSeq Select); táxon com mais de um gene no grupo de
   ortólogos sai do loco; CDS tem de ser múltiplo de 3, sem stop interno, e
   ter 0,5×–1,5× a mediana do loco. **686 genes aprovados** (≥ 40 dos 46
   táxons RefSeq, CDS do galo ≥ 600 pb, proteína de consulta disponível). Fora:
   417 com < 40 táxons RefSeq, 74 com CDS do galo ausente ou < 600 pb, 3 sem
   consulta.
4. Depois do miniprot, na mesma ordem da permutação, ficaram os **primeiros 300
   com ≥ 49/57 táxons** (o pool foi usado até a posição 554 de 1 200).
   Nos 300 locos, 197 registros RefSeq saíram por CDS inválido (stop interno ou
   fora de quadro) e 82 por comprimento fora de 0,5×–1,5× a mediana; nenhum
   por duplicata no grupo de ortólogos.

Tamanho: CDS do galo com mediana de 1 557 pb (600–14 604); soma das medianas de
proteína por loco, 201 989 aa (média 673 aa por loco), antes da limpeza.

**Ortologia.** RefSeq: grupos 1:1 do NCBI Orthologs, mais o filtro de um gene
por táxon. miniprot: alinhamento Rank = 1 da proteína de consulta no genoma
inteiro, com os filtros da §1.3; não há melhor hit recíproco contra o proteoma
inteiro. Retrotradução **não** entra na ortologia: é só o passo planejado para
a matriz DNA 1+2 (§4), que reconstrói os códons a partir do alinhamento de AA.

### 1.3 Auditoria do miniprot

miniprot 0.18 (`-M 2`, `--gff`), CDS extraído com gffread 0.12.9. Registro sai
do loco com frameshift, stop interno, cobertura da consulta < 50 % ou
identidade < 70 %; isso vira dado faltante, **não** é tratado como lesão
genética (nenhuma inferência de pseudogene é feita aqui).

| Táxon | Contig N50 | ok / 686 | frameshift | stop interno | outros | ocupação final (300) |
|---|---|---|---|---|---|---|
| *Sarcoramphus papa* | 40,9 Mb | 668 | 14 | 3 | 1 | 294 (98,0 %) |
| *Pithecophaga jefferyi* | 45,2 Mb | 642 | 30 | 13 | 1 | 294 (98,0 %) |
| *Leptodon cayanensis* | 143 kb | 649 | 24 | 12 | 1 | 292 (97,3 %) |
| *Gypaetus barbatus* | 18,1 Mb | 640 | 36 | 9 | 1 | 287 (95,7 %) |
| *Aegypius monachus* | 34,5 Mb | 633 | 41 | 11 | 1 | 287 (95,7 %) |
| *Ninox strenua* | 141 kb | 638 | 27 | 16 | 5 | 283 (94,3 %) |
| *Ibycter americanus* | 45 kb | 632 | 42 | 11 | 1 | 281 (93,7 %) |
| *Herpetotheres cachinnans* | 62 kb | 625 | 48 | 9 | 4 | 279 (93,0 %) |
| *Sagittarius serpentarius* | 34 kb | 608 | 61 | 13 | 4 | 278 (92,7 %) |
| *Elanus axillaris* | 30 kb | 617 | 48 | 18 | 3 | 278 (92,7 %) |
| *Pandion haliaetus* | 4,9 kb | 526 | 122 | 36 | 2 | 237 (79,0 %) |

"Outros" = cobertura < 50 % + identidade < 70 %. Frameshift + stop por táxon:
2,5 % (*Sarcoramphus*) a 23,0 % (*Pandion*). Mesmo nas montagens de
cromossomo (*Gypaetus*, contig N50 18 Mb) ficam 6,6 %: parte dessas lesões é
artefato de alinhamento ou diferença real de modelo de gene em relação à
consulta, não só erro de montagem — por isso o registro sai em vez de ser
corrigido. Identidade média com a consulta nos registros aceitos: 0,961
(*Pandion*) a 0,994 (*Sarcoramphus*). Um segundo alinhamento quase empatado
(≥ 0,99 do escore) apareceu em 0–11 locos por táxon (11 no *Leptodon*); o loco
fica com o Rank = 1.

`-M 2` (metade da amostragem de modimizadores no índice) foi adotado depois de
o WSL reiniciar durante o miniprot do *Pandion* (≈ 00:52 no relógio do WSL) com o padrão `-M 1`
(pico medido de 10,2–11,5 GB num WSL de 15 GB dividido com outra análise; a
causa do reinício não foi apurada). No *Sagittarius*, `-M 2` deu pico de
6,5 GB e 684 dos 686 alinhamentos Rank = 1 no mesmo lugar que `-M 1`. As 11
montagens foram rodadas de novo, todas com `-M 2`. Sobrou uma cópia
descompactada temporária (1,25 GB) em
`~/rapinantes/raw/GCA_026930845.1/genoma.fna`, que o hook de proteção de
`raw/` não deixa apagar por aqui.

**Controle na mesma montagem.** *Gypaetus* e *Sarcoramphus* têm anotação NCBI
de GenBank (GB_2026_01_30, GB_2025_12_03). A proteína do miniprot foi
comparada ao proteoma dessa anotação (DIAMOND blastp `--ultra-sensitive`, `-k 1`):

| Táxon | registros | com hit | ident. ≥ 98 % e cob. ≥ 90 % | identidade média |
|---|---|---|---|---|
| *Gypaetus barbatus* | 640 | 640 | 591 (92,3 %) | 99,42 % |
| *Sarcoramphus papa* | 668 | 668 | 621 (93,0 %) | 99,75 % |

Nos 6 casos do *Gypaetus* abaixo de 90 % contra a anotação (BCL2L1, DECR2,
ERO1B, KCNQ2, NFU1, PEX11A; 42,6–89,4 %), a proteína do miniprot tem 94,5–99,1 % de
identidade com a consulta de *Aquila*: o melhor hit da anotação é outra
isoforma ou outro membro da família, e o registro do miniprot é compatível com
o ortólogo. O controle mede a fidelidade da extração, não a ortologia.

### 1.4 Ocupação

RefSeq: média 95,4 %, com 7 táxons abaixo de 85 % — *Merops nubicus* 74,3 %,
*Antrostomus carolinensis* 78,3 %, *Eurypyga helias* 79,0 %, *Buceros
rhinoceros* 80,3 %, *Fulmarus glacialis* 83,0 %, *Acanthisitta chloris* 84,0 %
e *Apaloderma vittatum* 84,7 %, todos anotações de 2014 sobre montagens em
scaffold (N50 de contig 15–29 kb). miniprot: média 93,6 % (*Pandion* 79,0 % a
*Sarcoramphus*/*Pithecophaga* 98,0 %). Nenhum táxon fica abaixo de 50 %.

## 2. Amostragem diante dos nós contestados

| Grupo | Táxons (n) | Observação |
|---|---|---|
| Cathartidae | *Gymnogyps* (RefSeq), *Sarcoramphus* (miniprot) — 2 | quebra o ramo longo dos urubus; sem *Cathartes*/*Coragyps* |
| Sagittariidae | *Sagittarius* — 1 | miniprot |
| Pandionidae | *Pandion* — 1 | miniprot; a montagem mais fragmentada (79 % de ocupação) |
| Accipitridae | 10: Elaninae (*Elanus*), Gypaetinae (*Gypaetus*), Perninae (*Leptodon*), Aegypiinae (*Aegypius*), Circaetinae (*Pithecophaga*), Harpiinae (*Harpia*), Aquilinae (*Aquila*), Accipitrinae (*Astur*), Haliaeetinae (*Haliaeetus albicilla*), Buteoninae (*Buteo*) | Circinae e Milvinae de fora (teto de genomas) |
| Falconidae | *Herpetotheres* (Herpetotherinae), *Ibycter* (Polyborinae), *Falco peregrinus*, *F. naumanni* — 4 | sem *Micrastur* nem *Microhierax* |
| Strigiformes | *Tyto* (Tytonidae), *Strix* (Striginae), *Athene* (Surniinae), *Ninox* (Ninoxinae) — 4 | sem Asioninae, *Bubo*, *Phodilus* |
| Cariamiformes | *Cariama* — 1 | sem *Chunga*; ramo longo e sem quebra |
| Coraciimorphae | *Colius*, *Leptosomus*, *Apaloderma*, *Buceros*, *Merops*, *Dryobates*, *Indicator*, *Pogoniulus* — 8 | uma espécie por ordem, exceto Piciformes (3) |
| Psittacopasserae | *Strigops*, *Nestor*, *Melopsittacus*; *Acanthisitta*, *Chiroxiphia* (suboscine), *Taeniopygia*, *Corvus* (oscines) — 7 | |
| Hipóteses clássicas | *Antrostomus* (Caprimulgidae, H5), *Calypte*, *Apus* (Strisores), *Ciconia boyciana* (H6), *Rissa* (Charadriiformes) | um representante por hipótese |
| Resto de Neoaves | *Columba*, *Tauraco*, *Chlamydotis*, *Cuculus*, *Opisthocomus*, *Grus*, *Eurypyga*, *Phaethon lepturus*, *Gavia*, *Aptenodytes*, *Fulmarus*, *Pelecanus* | contexto, não objeto |
| Raiz | *Gallus gallus* (grupo externo da apresentação), *Anas platyrhynchos* | Galloanserae |

Accipitridae tem 10 táxons, acima do "~6 por família": a pergunta pedia uma
espécie por subfamília, e as subfamílias venceram o teto.

## 3. Parâmetros do pipeline (execução AA principal)

Comando exato em `runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit/executar.sh`:

```
uv run python scripts/filogenia/filogenia.py --input data/proc/rapinantes/filogenomica/locos_aa \
  --outgroup Gallus_gallus --tipo-seq AA --cpus-concatenado 8 --nome rapinantes_AA_mafft_clipkit \
  -- -c /mnt/c/Users/Manec/bio-stack/scripts/rapinantes/filogenomica/limite_10cpus.config
```

- **Alinhamento:** MAFFT L-INS-i (`--localpair --maxiterate 1000 --threadit 0`;
  todos os locos têm ≤ 200 sequências), módulo nf-core.
- **Limpeza:** ClipKIT, modo `smart-gap`.
- **Árvore por loco:** IQ-TREE 3, `-m MFP` (ModelFinder, BIC), `-B 1000 --bnni`
  (UFBoot2), `--alrt 1000` (SH-aLRT), `--abayes`, `--seed 12345`,
  `-o Gallus_gallus`, `-st AA`, `-T 1`, `-mem 1G`.
- **Supermatriz:** partição por loco, `-m MFP+MERGE` (rclusterf), mesmos
  suportes, `--seed 12345`, **8 threads** (`--cpus-concatenado 8`: a supermatriz
  AA de ~200 mil sítios não roda em tempo útil com 1; UFBoot e SH-aLRT podem
  variar 1–2 pontos entre execuções, e o `metodos.md` vai declarar isso), `-mem 9G`.
- **Concordância:** gCF (árvores de genes) e sCFL (100 quartetos) sobre a
  árvore da supermatriz, com o esquema de partição que ela escolheu.
- **Árvore de espécies:** wASTRAL, modo híbrido, ponderado pelo aBayes,
  semente 12345.
- **Recursos:** `limite_10cpus.config` limita o executor local a 10 CPUs e
  13 GB (outra análise pesada divide o WSL) e alarga os tempos máximos; não
  muda nenhum argumento das ferramentas.

As versões que de fato rodaram vão ser copiadas do `metodos.md` da execução
quando ela terminar — não estão escritas aqui à mão.

## 4. Status em 2026-10-18, 02:20

| Etapa | Estado |
|---|---|
| Amostragem, ortólogos, miniprot, auditoria, validação | **feito** |
| PREPARA_LOCO, MAFFT | feito (300/300) |
| ClipKIT | em andamento (192/300) |
| Árvores de genes (IQ-TREE, AA) | **pendente** — rodando; 0/300 terminadas |
| Supermatriz AA (MFP+MERGE + suportes) | **pendente** |
| gCF / sCFL | **pendente** (depende da supermatriz) |
| wASTRAL | **pendente** (depende das 300 árvores de genes) |
| Matriz DNA 1+2 | **pendente** — `f05_dna12.py` pronto: códons tirados do CDS na ordem do alinhamento de AA do MAFFT, a mesma máscara de colunas do ClipKIT, 3ª posição removida; roda com `--alinhador nenhum --limpeza nenhuma` |
| AU / KH / SH (H1–H7, S1–S3; `f07_hipoteses.py`, `f07b_topologia.sh`) | **pendente** — depende das supermatrizes AA e DNA 1+2 |
| χ² gDF1 × gDF2 e sDF1 × sDF2 por ramo (`f08_nos_chave.py`) | **pendente** |
| Sensibilidade AA (FAMSA, trimAl, FAMSA + trimAl) | **pendente** — por custo, sem UFBoot/SH-aLRT (só aBayes), o que vai estar declarado |

Previsão, não medida: uma árvore de gene AA levou 21 min (PIGT, 576 aa) e
36 min (NUBP1, 323 aa) com 46 táxons em 1 thread; com até 9 em paralelo, as
300 levam ~10–15 h. A supermatriz AA (~200 mil sítios antes da limpeza, 57
táxons) passa da memória concedida e roda no modo de economia do IQ-TREE;
estimativa de 12 h ou mais.

## 5. O que este conjunto pode e não pode testar

**Pode testar**, com os limites da §2:

- monofilia de Accipitrimorphae e de Accipitriformes s.s. (dois Cathartidae
  quebram o ramo dos urubus);
- a resolução das três linhagens de Afroaves — Accipitrimorphae,
  Strigiformes, Coraciimorphae (H1/H1b/H2/H3) —, por concatenação, por wASTRAL
  e pelos fatores de concordância;
- as posições de *Sagittarius* e *Pandion* (S1–S3) e a raiz de Accipitridae
  (*Elanus* contra o resto);
- Falconiformes s.l. (H4), corujas + bacuraus (H5), Cathartidae + cegonhas
  (H6), *Cariama* irmão de Eufalconimorphae ou só de Falconidae (H7);
- a posição de *Herpetotheres* em Falconidae;
- assimetria gDF1 × gDF2 nos ramos curtos (ILS sozinha produz as duas
  alternativas na mesma frequência).

**Não pode testar, ou testa mal:**

- **Independência estatística das análises de 2014–2024.** A independência é
  de anotação e de ortologia, não de dado bruto: várias montagens são as
  mesmas de Jarvis et al. 2014 (*Cariama*, *Merops*, *Buceros*, *Apaloderma*,
  *Nestor*, *Acanthisitta*, *Eurypyga*, *Fulmarus*, *Phaethon lepturus*,
  *Antrostomus*), e *Herpetotheres* é montagem do B10K. Erro dessas montagens
  é compartilhado.
- **Sinal não codificante.** Só há éxons. A dependência do tipo de dado em
  Neoaves (Jarvis 2014; Reddy et al. 2017) fica fora; íntron, UCE e região
  intergênica não entram.
- **Ramos de Cariamiformes, Sagittariidae e Pandionidae** sem quebra: um táxon
  cada (sem *Chunga*, só uma espécie em cada família monotípica). Atração de
  ramos longos nesses três não pode ser descartada pela amostragem.
- **Detalhe interno de Coraciimorphae, Strigidae e Accipitridae.** Uma espécie
  por ordem em Coraciimorphae (exceto Piciformes); sem Asioninae nem *Bubo*;
  sem Circinae nem Milvinae.
- **H5 e H6 com um representante cada** (*Antrostomus*, com 78 % de ocupação;
  *Ciconia*): rejeição ou não rejeição vale para esses dois táxons.
- **Neoaves fora de Telluraves.** A amostragem fora do clado é contexto, rala
  demais para afirmar a raiz de Neoaves ou o irmão de Telluraves.
- **Tempo de divergência e direção de caracteres.** Precisam de datação
  calibrada e de reconstrução ancestral (outros braços do projeto).
- **Causa da discordância.** gDF1 × gDF2 separa ILS simétrica de uma
  alternativa dominante, mas não distingue introgressão de paralogia ou de
  atração de ramos longos. 300 locos codificantes têm pouco poder por loco
  para entrenós de fração de milhão de anos; gCF baixo com sCF mais alto vai
  indicar árvore de gene mal estimada, não conflito biológico.

## Resultados — árvores de gene e wASTRAL (2026-10-18, 08:50)

Execução `20261018-020210-rapinantes_AA_mafft_clipkit`: as 300 árvores de genes
(AA) e o wASTRAL terminaram às 08:41; a supermatriz começou às 08:40 e **não
entra nesta seção**. gCF e quartetos abaixo foram calculados sobre a
**topologia do wASTRAL** (`f08a_genes_astral.sh`: IQ-TREE 3.1.4 `--gcf`;
wASTRAL v1.25.3.8 `--mode 4 -C -u 2`, não ponderado), não sobre a supermatriz.
Tabelas: `runs/rapinantes/filogenomica/nos_chave_genes_astral_AA.tsv` (nós
pré-definidos), `ramos_astral_AA.tsv` (todos os 54 ramos),
`suporte_arvores_genes_AA.tsv`, `ramos_terminais_anomalos_AA.tsv`. Figura:
`runs/rapinantes/filogenomica/figuras/wastral_AA_300locos.png`.

### Qualidade das árvores de genes

- **Resolução por loco baixa.** Fração de ramos internos com SH-aLRT ≥ 80 **e**
  UFBoot ≥ 95, por árvore: mediana 0,119 (IQR 0,051–0,207; máximo 0,451). Em
  57 táxons, uma árvore de gene típica tem ~6 de 54 ramos confiáveis.
- **Limpeza agressiva.** 51 dos 300 locos perderam > 20 % das colunas no
  ClipKIT (mediana 7,6 %; máximo 41,9 %, CCDC69). Pela Regra 4, isso pede
  também uma execução com `--limpeza nenhuma` — pendente.
- **Ramos terminais anômalos** (> 5× a mediana do táxon e > 0,3 subst./sítio):
  150 em 73 locos, 1,12 % dos registros RefSeq contra 0,10 % dos miniprot.
  Dez passam de 1 subst./sítio: ATF6B/*Rissa* (7,0), STK10/*Fulmarus* (4,4),
  SNRK/*Pogoniulus* (4,0), SUMF1/*Opisthocomus* (3,9), NKIRAS1/*Athene* (2,2),
  CLCN4/*Pandion* (1,6; o único miniprot), ASAP1/*Phaethon* (1,5),
  P2RY12/*Acanthisitta* (1,5), POC1A/*Dryobates* (1,4), WBP1L/*Rissa* (1,2).
  São suspeitos de parálogo, modelo de gene quimérico ou erro de extração.
  Estão dentro da supermatriz que está rodando; uma reanálise sem eles
  (poda à maneira do TreeShrink) fica pendente. Como são ramos terminais,
  pesam pouco nos nós profundos, mas entram no gCF dos clados desses táxons.

### Topologia do wASTRAL

Sustentado (localPP = 1,00, q1 ≥ 0,52):

- **Accipitrimorphae** monofilético (Cathartidae + Accipitriformes s.s.).
- **Sagittarius** irmão de **Pandion + Accipitridae**; *Elanus* irmão do resto
  de Accipitridae; dentro de Accipitridae: (*Gypaetus*, *Leptodon*),
  (*Aegypius*, *Pithecophaga*), *Aquila*, *Harpia* (PP 0,675) e (*Astur*,
  (*Haliaeetus*, *Buteo*)).
- **Falconiformes**, com *Herpetotheres* irmão de (*Ibycter*, *Falco*);
  **Strigiformes**, com *Tyto* irmão de Strigidae; Psittacopasserae;
  Passeriformes; Piciformes + *Merops* (PP 0,99).

Não sustentado — o esqueleto de Telluraves:

- *Leptosomus* sai de Coraciimorphae e fica irmão de (Accipitrimorphae +
  Strigiformes) (PP 0,560). Afroaves, Australaves, Coraciimorphae e
  Eufalconimorphae **não aparecem** na árvore.
- (Accipitrimorphae + Strigiformes + *Leptosomus*) é irmão de (Falconidae +
  *Cariama*) (PP 0,919); o outro clado de Telluraves é Psittacopasserae +
  Coraciimorphae sem *Leptosomus* (PP 0,950).
- Fora de Telluraves, o esqueleto de Neoaves também é instável (PP 0,43–0,93;
  Strisores parte: *Antrostomus* longe de *Calypte* + *Apus*), o que é esperado
  com essa amostragem (§5).

### Nós contestados

q1/q2/q3: frequência dos quartetos na topologia do wASTRAL e nas duas
alternativas NNI (não ponderado); p(q2 = q3): χ² sobre as contagens efetivas
f2 e f3. gCF/gDF1/gDF2/gDFP em %, gene a gene sobre a topologia do wASTRAL.
"Genes fortes": árvores de genes com o clado (restrito aos táxons do gene) a
SH-aLRT ≥ 80 e UFBoot ≥ 95, a favor / contra (ramo forte incompatível).

| Nó | no wASTRAL | localPP | q1 / q2 / q3 | p(q2 = q3) | gCF / gDF1 / gDF2 / gDFP | genes fortes a favor / contra |
|---|---|---|---|---|---|---|
| Telluraves | sim | 0,939 | 0,380 / 0,360 / 0,259 | 0,028 | 0,7 / 0 / 0 / 99,3 | 0 / 23 |
| Afroaves | **não** | — | — | — | — | 0 / 23 |
| Australaves | **não** | — | — | — | — | 0 / 18 |
| Accipitrimorphae | sim | 1,000 | 0,425 / 0,280 / 0,295 | 0,74 | 7,7 / 0,3 / 0,3 / 91,7 | 1 / 11 |
| Accipitriformes s.s. | sim | 1,000 | 0,643 / 0,174 / 0,184 | 0,78 | 39,6 / 0,7 / 2,2 / 57,6 | 28 / 9 |
| *Pandion* + Accipitridae (S1) | sim | 1,000 | 0,658 / 0,153 / 0,189 | 0,35 | 45,8 / 3,1 / 4,4 / 46,7 | 58 / 10 |
| *Sagittarius* + Accipitridae (S2) | não | — | — | — | — | 11 / 42 |
| *Sagittarius* + *Pandion* (S3) | não | — | — | — | — | 0 / 33 |
| Strigiformes | sim | 1,000 | 0,602 / 0,189 / 0,209 | 0,59 | 33,8 / 0 / 0,3 / 65,9 | 26 / 5 |
| Coraciimorphae | **não** | — | — | — | — | 0 / 16 |
| Strigiformes + Coraciimorphae (H1/H3) | não | — | — | — | — | 0 / 18 |
| Strigiformes + Accipitrimorphae (H2) | sim | 0,736 | 0,356 / 0,310 / 0,334 | 0,59 | 1,0 / 1,0 / 1,7 / 96,3 | 0 / 11 |
| Accipitrimorphae + Coraciimorphae (H1b) | não | — | — | — | — | 0 / 24 |
| (Accip. + Strig.) + *Leptosomus* | sim | 0,560 | 0,365 / 0,348 / 0,288 | 0,19 | 0,7 / 0 / 0,7 / 98,7 | — |
| Falconidae + Accipitrimorphae (H4) | não | — | — | — | — | 0 / 12 |
| Falconiformes | sim | 1,000 | 0,810 / 0,085 / 0,105 | 0,49 | 64,6 / 0,4 / 0,4 / 34,6 | 126 / 4 |
| *Cariama* + Falconiformes (H7b) | sim | 0,729 | 0,361 / 0,311 / 0,327 | 0,75 | 5,4 / 0 / 0 / 94,6 | 23 / 8 |
| Eufalconimorphae (H7a) | **não** | — | — | — | — | 0 / 15 |
| Psittacopasserae + Coraciimorphae s/ *Leptosomus* | sim | 0,950 | 0,385 / 0,304 / 0,311 | 0,89 | 0 / 0 / 0 / 100 | — |
| Strigiformes + Caprimulgidae (H5) | não | — | — | — | — | 7 / 6 |
| Cathartidae + Ciconiiformes (H6) | não | — | — | — | — | 0 / 11 |

### Leitura (Regras 1 e 2)

1. **O que as árvores de genes sustentam:** os nós dentro de ordem e de
   família. Accipitriformes s.s., *Sagittarius* irmão de *Pandion* +
   Accipitridae, Falconiformes e Strigiformes têm localPP 1,00 e q1 de 0,60 a
   0,81, com mais genes fortes a favor que contra (Falconiformes 126/4;
   *Pandion* + Accipitridae 58/10, contra 11/42 para *Sagittarius* +
   Accipitridae). Accipitrimorphae é o caso de fronteira: localPP 1,00, mas q1
   0,425, gCF 7,7 % e 1 gene forte a favor contra 11 — a monofilia aparece no
   somatório dos quartetos, quase nunca inteira numa árvore de gene.
2. **O que elas não resolvem:** o esqueleto de Telluraves. Nos ramos que
   decidem a posição das linhagens de rapinantes (H2, *Leptosomus*, H7b, os
   dois clados-irmãos de Telluraves), q1 fica entre 0,356 e 0,385, perto do
   1/3 de uma politomia; localPP entre 0,56 e 0,95 (só um ramo chega a 0,950);
   gCF de 0 a 5 %, com gDFP de 95 a 100 %. Quase nenhuma árvore de gene contém
   esses clados inteiros, nem as alternativas NNI: o χ² gDF1 × gDF2 fica sem
   contagem para testar, e o teste útil é o de quartetos. **Isso é sinal fraco
   por loco (árvores de genes mal resolvidas, mediana de 12 % de ramos
   fortes), não evidência de conflito biológico.** A ausência de Afroaves,
   Australaves e Coraciimorphae no wASTRAL não contradiz Jarvis 2014, Prum
   2015 ou Stiller 2024: com localPP < 0,95 e q1 ≈ 1/3, este conjunto não
   escolhe entre as alternativas.
3. **ILS × sinal sistemático nos ramos curtos:** em nenhum dos ramos
   contestados as duas alternativas diferem (p(q2 = q3) de 0,19 a 0,89) —
   compatível com ILS simétrica, mas também com falta de poder. A exceção é
   Telluraves (p = 0,028: q2 0,360 contra q3 0,259), onde uma alternativa
   domina; com q1 0,380 quase empatado com q2, a monofilia de Telluraves nas
   árvores de genes é frágil.
4. **Hipóteses clássicas nas árvores de genes:** H4 (Falconiformes s.l.), H6
   (Cathartidae + cegonhas) e H1b não aparecem no wASTRAL, e nenhum gene as
   sustenta com ramo forte (0 a favor). H5 (corujas + bacurau) também não
   aparece, mas tem 7 genes fortes a favor e 6 contra — sinal de poucos genes,
   que o wASTRAL não acompanha. O teste formal dessas hipóteses é o AU sobre a
   supermatriz (pendente).
5. **Status depois desta seção:** supermatriz AA (MFP+MERGE, 8 threads) em
   andamento desde 08:40, com gCF/sCFL em seguida; árvores de genes da matriz
   DNA 1+2 rodando em paralelo em 2 CPUs (`--concatenar false`; a supermatriz
   DNA 1+2 entra depois, por `-resume`). Os nós profundos vão ser julgados de
   novo pelo gCF e pelo sCFL da supermatriz e pelo teste AU.

### DNA 1+2 — árvores de gene e wASTRAL (2026-10-18, 11:30)

Execução `20261018-084205-rapinantes_DNA12`, fase de árvores de gene
(`--concatenar false`, 2 CPUs; comando em `executar_fase1_genes.sh`). Matriz:
códons do CDS colocados no alinhamento de AA do MAFFT da execução principal,
mesma máscara do ClipKIT, 3ª posição fora (`f05_dna12.py`): 300 locos,
430 732 sítios; a tradução de cada CDS conferiu com a sequência de AA
alinhada em todos os registros. Pipeline com `--alinhador nenhum --limpeza
nenhuma --tipo-seq DNA`, mesmos suportes e semente. Tabelas:
`nos_chave_genes_astral_DNA12.tsv`, `ramos_astral_DNA12.tsv`; figura
`figuras/wastral_DNA12_300locos.png`.

- Resolução por loco igual à do AA: mediana de 0,125 dos ramos internos com
  SH-aLRT ≥ 80 e UFBoot ≥ 95 (IQR 0,057–0,208). Ramos terminais anômalos: 35
  em 27 locos (0,26 % dos registros RefSeq, 0,03 % dos miniprot), contra 150
  no AA; o limiar absoluto de 0,3 subst./sítio não tem o mesmo sentido em AA
  e em DNA 1+2, então a diferença não diz qual matriz tem mais parálogos.
- Mesmos nós sustentados que no AA (localPP 1,00): Accipitrimorphae (q1 0,435;
  gCF 8,3), Accipitriformes s.s. (0,656; 38,1), *Pandion* + Accipitridae
  (0,693; 48,5; genes fortes 65 a favor / 16 contra, contra 10/50 para
  *Sagittarius* + Accipitridae), Accipitridae com *Elanus* primeiro,
  Strigiformes, Falconiformes com *Herpetotheres* primeiro.
- **Telluraves não sai monofilético:** *Calypte* + *Apus* entram como irmãos
  de Psittacopasserae (PP 0,669). Afroaves, Australaves, Coraciimorphae,
  Eufalconimorphae, H1/H1b/H2/H3 não aparecem. As linhagens de rapinantes
  formam, de novo, um agrupamento só: ((Strigiformes + *Leptosomus*) [PP
  0,606] + Accipitrimorphae) [0,897] + (Falconidae + *Cariama*) [0,608], esse
  último nó com PP 0,664. Todos com q1 entre 0,338 e 0,360.
- **Leitura.** AA e DNA 1+2 concordam em reunir Accipitrimorphae,
  Strigiformes, *Leptosomus*, Falconidae e *Cariama* num mesmo clado, mas
  com localPP < 0,95 e q1 ≈ 1/3 em todos os ramos que fazem isso. As duas
  matrizes vêm dos mesmos genes e do mesmo alinhamento: a concordância entre
  elas **não** é replicação independente. Esse agrupamento lembra o antigo
  "grupo das aves de rapina" e é exatamente o tipo de sinal que convergência
  em genes codificantes ou atração entre ramos curtos produziria; com este
  poder, não há como separá-lo de ruído. Não deve ser lido como apoio a
  Falconiformes s.l. nem a um clado de rapinantes.
- Assimetria de alternativas (χ² q2 × q3): só no ramo Telluraves + *Calypte*/
  *Apus* (p = 0,034; q1 0,380, q2 0,358, q3 0,261 — o mesmo padrão do ramo de
  Telluraves no AA), em Psittacopasserae (p = 0,009) e em (*Ibycter* +
  *Falco*) (p = 0,002); em nenhum ramo que envolve linhagens de rapinantes.

### Status e aviso operacional (2026-10-18, 19:55)

- Supermatriz AA (`20261018-020210-rapinantes_AA_mafft_clipkit`): ModelFinder
  com fusão terminou em 5 h 48 min (300 → 48 partições; 3,9 GB de RAM, sem
  modo de economia); busca ML convergiu na iteração ~103 (lnL −2 230 208); o
  refinamento das 1000 réplicas do UFBoot com `--bnni` está a ~30 s por
  réplica (100 em 50 min), o que projeta o fim da supermatriz para a
  madrugada de 2026-10-19; gCF/sCFL vêm depois. Execução não tocada.
- **A retomada da execução FAMSA (tarefa bp7kj6ggo) não reaproveitou o
  cache.** O log mostra 300 PREPARA_LOCO, 300 FAMSA_ALIGN e 300 CLIPKIT
  "Submitted" e nenhum "Cached", e as árvores de genes estão sendo refeitas:
  o `bin/filogenia_util.py` editado entra como arquivo de entrada do
  PREPARA_LOCO, o hash dele muda e tudo a jusante roda de novo (~7 h com 8
  CPUs, pela primeira passada). Como FAMSA, ClipKIT e IQ-TREE `-T 1` com
  semente fixa são determinísticos, as 300 árvores devem sair iguais às da
  primeira passada; o custo é tempo. O mesmo vale para a fase 2 da DNA 1+2
  (`--retomar`): as árvores de genes vão rodar de novo.
- **Disco C: cheio (2026-10-18, ~20:10): 1,6 GB livres de 476 GB** (eram 94 GB
  de manhã). Este braço ocupa ~11 GB dentro do WSL (`~/rapinantes`, dos quais
  4,1 GB de montagens e 1,25 GB de cópia temporária em `raw/`) e ~0,4 GB em
  `runs/filogenia`; o `ext4.vhdx` está com 55 GB e `AppData/Local/Temp` com
  17 GB, nenhum dos grandes itens de lá é deste braço. Com o disco nesse
  ponto, a supermatriz AA corre risco de morrer ao gravar a saída final.
  Nada foi apagado por aqui (o hook de `raw/` bloqueia, e os outros itens não
  são deste braço).

## Resultados — supermatriz AA, concordância e testes de topologia (2026-10-19)

A supermatriz da execução `20261018-020210-rapinantes_AA_mafft_clipkit` morreu por
SIGHUP às 21:09 de 2026-10-18, no refinamento dos UFBoot. Foi refeita com os 48 modelos
que o ModelFinder tinha escolhido, fixados, e terminou às 16:27 de 2026-10-19. A
história completa está no `RETOMADA_MANUAL.md` da execução e na seção 4.7 do
`docs/metodos_rapinantes.md`.

- **Matriz:** 57 táxons, 215 366 sítios (99 768 informativos por parcimônia), 48
  grupos de partições.
- **Resultado:** lnL −2 230 233,545; UFBoot médio de 89,4, com 72 % dos 54 ramos
  internos em ≥ 95.
- **Comparação com a execução interrompida:** a topologia é idêntica à da melhor
  árvore da busca interrompida (Robinson-Foulds 0).

Tabelas:

- `runs/rapinantes/filogenomica/nos_chave_AA.tsv` (`f08_nos_chave.py`);
- `scfl_alternativas_AA.tsv` (`f08d`/`f08e`): o sCF de cada alternativa NNI, medido
  na árvore restrita que a força;
- `splits_genes.tsv` (`f08f`): as bipartições concorrentes contadas nas 300 árvores
  de genes das três execuções;
- `topologia_AA/` (`f07_hipoteses.py`, `f07b_topologia.sh`): as buscas restritas e o
  `au.iqtree`.

Figura: `runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit/figuras/concatenado.png`.

### Topologia da supermatriz

Telluraves = Australaves + Afroaves, com a ordem interna que Jarvis et al. 2014
(TENT) e a concatenação do B10K dão.

- **Australaves:** (*Cariama*, (Falconiformes, Psittacopasserae)). Australaves tem
  SH-aLRT/aBayes/UFBoot 100/1/96; Eufalconimorphae, 99/1/95.
- **Afroaves:** (Accipitrimorphae, (Strigiformes, Coraciimorphae)). Afroaves tem
  98,6/1/93; Strigiformes + Coraciimorphae, 93,8/1/90. *Leptosomus* é o primeiro ramo
  de Coraciimorphae (99,8/1/53). A ordem interna de Cavitaves é fraca (UFBoot 47–55).
- **Accipitrimorphae:** (Cathartidae, (*Sagittarius*, (*Pandion*, Accipitridae))). Em
  Accipitridae, *Elanus* sai primeiro, seguido de (*Gypaetus*, *Leptodon*), depois
  (*Aegypius*, *Pithecophaga*), depois *Aquila*, depois *Harpia*, e por fim (*Astur*,
  (*Buteo*, *Haliaeetus*)). Todos os ramos têm 100/1/100. É a espinha de Catanach et
  al. 2024.
- **Falconidae:** *Herpetotheres* irmão de (*Ibycter*, *Falco*), com 100/1/100.
- **Strigiformes:** *Tyto* irmão de Strigidae (100/1/100). Dentro de Strigidae,
  (*Ninox*, *Strix*) tem só 76,3/1/79; Salter et al. 2020 põem *Ninox* junto de
  *Athene*.
- **Grupo-irmão de Telluraves:** o clado aquático (Aequornithes, Phaethontimorphae,
  *Grus*, *Rissa*), com UFBoot 63. Está aberto, como esperado da amostragem.

### Concordância nos nós profundos

| Nó | SH-aLRT/aBayes/UFBoot | gCF (N/gN) | sCF / sDF1 / sDF2 | p(sDF1 = sDF2) | wASTRAL (localPP) |
|---|---|---|---|---|---|
| Telluraves | 100/1/100 | 0,7 % (2/300) | 45,7 / 19,6 / 34,7 | 1,7 × 10⁻⁶ | sim (0,939) |
| Australaves | 100/1/96 | 0,4 % (1/258) | 53,1 / 20,9 / 26,0 | 0,14 | não |
| Eufalconimorphae | 99/1/95 | 2,3 % (6/258) | 41,8 / 28,8 / 29,5 | 0,81 | não |
| Afroaves | 98,6/1/93 | 0 (0/300) | 33,7 / 34,2 / 32,1 | 0,56 | não |
| Strigiformes + Coraciimorphae | 93,8/1/90 | 0,3 % (1/300) | 43,1 / 24,5 / 32,5 | 0,0095 | não |
| Accipitrimorphae | 100/1/100 | 7,7 % (23/300) | 51,7 / 24,7 / 23,7 | 0,76 | sim (1,00) |
| Psittacopasserae | 100/1/98 | 7,7 % (23/299) | 51,9 / 26,1 / 22,1 | 0,044 | sim (0,999) |
| Accipitriformes s.s. | 100/1/100 | 39,6 % (110/278) | 74,9 / 10,4 / 14,8 | 0,0098 | sim (1,00) |
| *Pandion* + Accipitridae | 100/1/100 | 45,8 % (104/227) | 72,5 / 11,6 / 15,9 | 0,040 | sim (1,00) |
| Strigiformes | 100/1/100 | 33,8 % (100/296) | 67,5 / 16,7 / 15,8 | 0,62 | sim (1,00) |
| Falconiformes | 100/1/100 | 66,0 % (184/279) | 84,1 / 5,9 / 10,1 | 3,9 × 10⁻⁶ | sim (1,00) |

**Quais são as alternativas.** O `.cf.stat` não dá o nome das alternativas NNI.
Elas foram identificadas por dois caminhos:

- **sDF:** o sCFL de cada árvore restrita (`f08d`/`f08e`) reproduz o sDF do ramo
  correspondente da ML.
- **gDF:** pelo número de árvores de genes que contêm cada bipartição (`f08f`).

**Eufalconimorphae.** O ramo separa quatro linhagens: Falconiformes,
Psittacopasserae, *Cariama* e o resto.

- O sCF de *Cariama* + Falconiformes (árvore H7b) é 29,6 e o de *Cariama* +
  Psittacopasserae é 28,6: as duas alternativas empatam nos sítios.
- Nas árvores de genes, a alternativa *Cariama* + Falconiformes é o gDF2: 14 árvores,
  contra 6 concordantes e 0 para *Cariama* + Psittacopasserae (p = 0,0002).
- Contada solta, sem condicionar nas quatro linhagens, a bipartição *Cariama* +
  Falconiformes aparece em 41–45 das 300 árvores de genes nas três execuções. Em 23
  delas o ramo é forte (MAFFT e DNA 1+2). Eufalconimorphae aparece em 6–10 árvores, e
  *Cariama* + Psittacopasserae em 3–6.

**Base de Afroaves.** O ramo separa Strigiformes, Coraciimorphae, Accipitrimorphae e o
resto. Os sítios se dividem em 43,1 % para corujas + Coraciimorphae, 32,5 % para
Hieraves (sCF da árvore H2: 32,6) e 24,5 % para Accipitrimorphae + Coraciimorphae
(sCF da H1b: 24,2). As duas minoritárias são desiguais (p = 0,0095), com excesso para
Hieraves. As árvores de genes quase nunca contêm qualquer um desses três clados
(0–3 de 300).

### Testes de topologia (KH, SH, AU; RELL 10 000)

| Árvore | Hipótese | ΔlnL | p-KH | p-SH | p-AU |
|---|---|---|---|---|---|
| ML | — | 0 | 0,887 | 1 | 1 |
| H2 | Strigiformes + Accipitrimorphae (Hieraves) | 34,8 | 0,113 | 0,857 | **0,130** |
| H7b | *Cariama* + Falconiformes | 55,1 | 0,020 | 0,777 | **0,019** |
| H1b | Accipitrimorphae + Coraciimorphae | 88,1 | 0,0004 | 0,579 | < 0,001 |
| wASTRAL | árvore de espécies inteira | 455,0 | 0,0005 | 0,008 | 0,0007 |
| H4 | Falconidae + Accipitrimorphae | 391,7 | 0 | 0,005 | < 0,001 |
| H4b | Falconidae + Accipitriformes s.s. | 604,1 | 0 | 0 | < 0,001 |
| H5 | Strigiformes + Caprimulgidae | 418,0 | 0 | 0,003 | < 0,001 |
| H6 | Cathartidae + Ciconiiformes | 1544,3 | 0 | 0 | < 0,001 |
| S2 | *Sagittarius* + Accipitridae | 1095,7 | 0 | 0 | < 0,001 |
| S3 | *Sagittarius* + *Pandion* | 1084,4 | 0 | 0 | < 0,001 |

H1, H3, H7a e S1 já estão na árvore ML e não têm busca. Os p-AU muito pequenos saem do
ajuste do bootstrap multiescala (H1b dá 4 × 10⁻⁸⁹ com ΔlnL 88). Por isso estão escritos
só como < 0,001.

### Leitura

1. **Concatenação e coalescente divergem justamente onde a literatura diverge.** A
   supermatriz dá Australaves, Eufalconimorphae e Afroaves com Accipitrimorphae na
   base. O wASTRAL dos mesmos 300 locos dá *Cariama* + Falconiformes, junta todas as
   linhagens de rapinantes num clado e não tem Afroaves nem Australaves. O teste AU
   rejeita a árvore do wASTRAL inteira (p = 0,0007) e rejeita *Cariama* +
   Falconiformes (p = 0,019).
2. **O sinal *Cariama* + Falconiformes está nas árvores de genes, não nos sítios.**
   Nas árvores de genes ele é forte e assimétrico: 41–45 contra 3–6 para a outra
   alternativa, nas três execuções, com MAFFT, com FAMSA e em DNA 1+2. No sCFL, que
   pontua sítios sobre quartetos da supermatriz com o modelo particionado, as duas
   alternativas empatam (29,6 × 28,6) e Eufalconimorphae vence (41,8). O padrão é o
   que §5 previa para éxons: o mesmo artefato de dado codificante descrito por Jarvis
   2014, Kuhl 2021 e nos éxons do B10K (relatório v1). Nas árvores de genes, ele não
   se separa de atração do ramo longo de *Cariama*, um táxon só.
3. **A base de Afroaves continua aberta.** A concatenação prefere corujas +
   Coraciimorphae (UFBoot 90; 43 % dos sítios), mas o AU não rejeita Hieraves
   (ΔlnL 34,8; p = 0,13). O sCF de Afroaves é 33,7, contra 34,2 e 32,1 nas
   alternativas, o que significa que nenhum sítio decide se Accipitrimorphae fica
   dentro ou fora do clado.

   As duas alternativas minoritárias são desiguais, com excesso para Hieraves
   (p = 0,0095), na mesma direção do B10K e do excesso de retroelementos de Gatesy &
   Springer 2022. Isso é compatível com reticulação entre a linhagem das corujas e a
   dos gaviões, mas não a demonstra.
4. **O que fica decidido.** Este conjunto sustenta dez resultados:
   - Accipitrimorphae;
   - Accipitriformes s.s.;
   - *Sagittarius* irmão de *Pandion* + Accipitridae, com S2 e S3 rejeitadas por mais
     de 1000 unidades de lnL;
   - *Elanus* na base de Accipitridae;
   - Falconiformes, com *Herpetotheres* primeiro;
   - Strigiformes;
   - a rejeição de Falconiformes s.l. (H4, H4b);
   - a rejeição de corujas + bacuraus (H5);
   - a rejeição de urubus + cegonhas (H6);
   - Australaves com *Cariama* fora de Eufalconimorphae, pelo AU, e não pelas árvores
     de genes.
5. **gCF baixo com sCF mais alto.** Nos nós profundos, o gCF fica entre 0 e 8 %, com
   gDFP de 84 a 100 %. O sCF fica entre 34 e 53 %. Como §5 previa, isso aponta árvores
   de genes mal estimadas, não conflito biológico. A exceção é o excesso *Cariama* +
   Falconiformes do item 2, que é sistemático.

## Sensibilidade ao alinhador (FAMSA) — wASTRAL (2026-10-19)

Execução `20261018-112050-rapinantes_AA_famsa_clipkit`: as mesmas proteínas, com FAMSA
2.4.1 e ClipKIT, árvores de genes só com aBayes e sem supermatriz. Terminou às 09:37.
Os clados foram lidos com `f07a_clados_wastral.py` e os grupos-irmãos com
`f07b_irmaos_wastral.py`, na árvore enraizada em *Gallus*. A tabela dá a localPP.

| Clado | FAMSA | MAFFT | DNA 1+2 |
|---|---|---|---|
| Accipitriformes, Accipitrimorphae, Strigiformes, Falconiformes | 1,00 | 1,00 | 1,00 |
| Strigiformes + Accipitrimorphae | 0,69 | 0,74 | não (*Leptosomus* + Strigiformes 0,61) |
| + *Leptosomus* | 0,59 | 0,56 | 0,90 |
| Falconiformes + *Cariama* | 0,77 | 0,73 | 0,61 |
| todos os rapinantes + *Cariama* + *Leptosomus* | 0,83 | 0,92 | 0,66 |
| Australaves / Eufalconimorphae | não | não | não |
| Psittacopasserae | 0,91 | 1,00 | 1,00 |

O alinhador não muda nada nos nós de rapinantes. Muda *Apus* + *Calypte*: irmãos de
Psittacopasserae no FAMSA (0,59) e no DNA 1+2 (0,67), na base de Neoaves no MAFFT
(1,00). Nas árvores de genes do FAMSA, *Cariama* + Falconiformes aparece em 43 das 300,
contra 6 para Eufalconimorphae e 6 para *Cariama* + Psittacopasserae
(`splits_genes.tsv`). É o mesmo padrão do MAFFT.

### Status (2026-10-19, 17:10)

| Etapa | Estado |
|---|---|
| Árvores de genes e wASTRAL: AA MAFFT, AA FAMSA, DNA 1+2 | feito |
| Supermatriz AA MAFFT (48 modelos fixos), gCF, sCFL | feito |
| χ² gDF1 × gDF2 e sDF1 × sDF2; nomes das alternativas | feito |
| AU/KH/SH (13 hipóteses + wASTRAL) na supermatriz AA | feito |
| Supermatriz DNA 1+2, supermatriz FAMSA, trimAl | **pendente** |
| Reanálise sem os ramos terminais anômalos | **pendente** |
