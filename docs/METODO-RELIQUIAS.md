# Método — relíquias de genes do esmalte em amniotas sem dente

Especificação do método dos scripts `scripts/atavismo/s24`–`s28`, escrita em
2026-10-17 seguindo a skill `methodology-architect` (problema → candidatos →
teste de estresse → crítica externa) e a escada de evidência da skill
`atavismo`. Resultados ficam em `docs/PESQUISA-2026-10-12-FENOTIPO-E-ENSEMBLE.md`
(§A9–A11); aqui fica só o **como** e o **quanto se pode afirmar**.

## 1. Problema

**Pergunta.** Depois da perda funcional da dentição (aves: origem única,
~116 Ma, Meredith 2014; tartarugas: origem independente), o que resta dos genes
efetores do esmalte no genoma, e o rastro é graduado por linhagem?

**Quem usa.** Genômica comparativa de perda de caráter (PhyloG2P); o índice de
retenção é candidato a fenótipo contínuo para RERconverge (§4 do projeto).

**Restrições.** Genomas públicos e anotação RefSeq; 8 GB de VRAM; sem dado de
leitura bruta processado localmente; n de espécies pequeno e não independente
(filogenia).

**Critério de sucesso.** Cada afirmação de perda carrega o nível de evidência
(§2) e o controle que a sustenta; nenhum zero vale como ausência sem controle
positivo na mesma espécie e no mesmo limiar.

## 2. Escada de evidência (skill `atavismo`)

| Nível | Critério | Onde o pipeline está |
|---|---|---|
| 1 | ausente da anotação RefSeq | `s26` PARTE 1, `anotacao_painel.json` |
| 2 | tBLASTn não acha, **com controle positivo no mesmo limiar** | `s26` (LPL), `s27` (âncoras por intervalo) |
| 3 | relíquia com **lesão identificável** e coordenada | coordenadas sim (`s25`, `s27_reciproco`); **lesões não** — §4, C1 |
| 4 | lesão **compartilhada** por espécies que divergiram depois — data a perda | não feito — §4, C2 |

## 3. Pipeline atual e seus controles

1. **Intervalo sintênico por gene** (`s26` para o cluster ENAM–AMBN–AMTN;
   `s27_painel_intervalos.py` para AMELX, MMP20, DSPP, ACP4). Âncoras: vizinho
   nomeado mais próximo de cada lado, exigido nas espécies sem dente e no
   *Anolis*; fallback de uma âncora quando o par se rompe.
2. **tBLASTn** das proteínas do jacaré (evalue 10, SEG ligado). Métrica:
   bitscore, não e-value (o e-value escala com o tamanho do banco; as janelas
   vão de 170 kb a 1,8 Mb).
3. **Competição local** (`s27_indice.py`): um trecho do alvo é descartado se
   outra proteína do intervalo alinhar ≥ 50% dele com bitscore maior.
4. **Melhor hit recíproco** contra o proteoma inteiro do jacaré
   (`s27_reciproco.sh`; `s25` para o cluster no frango).
5. **Índices**: melhor HSP ÷ âncoras (o do `s26`) e soma de trechos recíprocos ÷
   o próprio jacaré.

**Controles e o que cada um descarta**

| Controle | Descarta |
|---|---|
| âncora alinha nas sem dente | intervalo errado; método cego na espécie |
| LPL/DSPP negativo (22–30 bits) | define o piso de ruído |
| *Anolis*, o com dente mais distante | "o índice mede só distância filogenética" |
| tartaruga, perda independente | artefato específico de aves |
| recíproco no proteoma inteiro | parálogo de fora do intervalo |
| Ψ-AMEL conhecido (`s27_calibra_amelx.sh`) | **mede o ponto cego**: proteína de crocodiliano não acha pseudogene de ave divergido |

**O que o pipeline não permite afirmar** (skill `rigor-estatistico`): n = 6
espécies, não independentes; "com dente e sem dente não se sobrepõem" é
descrição, não teste. Nenhum p-valor foi calculado e nenhum deve ser escrito.

## 4. Candidatos de melhoria (gerados sem filtro, ranqueados por parcimônia)

| # | Candidato | Mecanismo | Ganho | Custo | Risco |
|---|---|---|---|---|---|
| C0 | *baseline simples*: só anotação | nível 1 | — | zero | perda falsa por pipeline de anotação |
| **C1** | **miniprot** (instalado, nunca usado) | alinhamento proteína→genoma com splicing; marca **frameshift** e **stop** | nível 3; mais sensível que tBLASTn em gene com íntron | baixo | frameshift falso por erro de montagem ou por alinhamento em região de baixa identidade |
| **C2** | lesões compartilhadas entre aves | posição da lesão na coordenada da proteína do jacaré, comparada entre *Gallus* e *Strigops* | nível 4; data a perda antes de Galloanserae × Neoaves | baixo, depois de C1 | coincidência em hotspot de indel; precisa de ≥ 2 lesões coincidentes |
| C3 | maquinário do `s27` no cluster do `s26` | competição + recíproco + soma para ENAM/AMBN/AMTN | método único para os 9 genes | baixo | — |
| C4 | mais espécies (tartarugas, aves, *Crocodylus*, escamados) | o mesmo pipeline | n de 6 para ~14; permite estatística descritiva com filogenia | médio (API) | anotação heterogênea entre montagens |
| C5 | consulta de ave, em nucleotídeo | relíquia de uma ave → outras aves por blastn | fecha o ponto cego medido no Ψ-AMEL | baixo | só funciona depois de a primeira relíquia de ave existir |
| C6 | *ambicioso*: RELAX (HyPhy) e datação por ω | seleção relaxada nos ramos sem dente; tempo de inativação por ω no ramo de transição | evidência independente de perda de restrição e data | alto (alinhamento de códons de fragmentos) | fragmentos curtos dão ω instável |
| C7 | *ambicioso*: potência por simulação na própria árvore | modo diagnóstico do TRACCER; simulação para RERconverge | decide se o §4 do projeto é viável | alto | — |

**Ordem de execução**: C1 → C2 → C5 → C3 → C4; C6 e C7 depois de haver
alinhamentos de relíquias em ≥ 4 espécies sem dente.

## 5. Teste de estresse de C1 (o candidato principal)

- **Argumento contra**: em relíquia de 116 Ma a identidade é baixa; o miniprot
  pode inventar frameshift para costurar alinhamento ruim.
- **Experimento que derruba**: rodar o miniprot nas espécies **com dente**, onde
  os genes são funcionais. Se aparecer frameshift ou stop ali na mesma taxa que
  nas aves, a lesão é artefato do alinhador, não do genoma.
- **Caso degenerado**: se o miniprot só reencontrar os HSPs do tBLASTn sem
  lesão nenhuma, ele não acrescenta nível de evidência — só sensibilidade.
- **Objeção do revisor**: erro de montagem. Resposta: lesão compartilhada por
  duas espécies em montagens independentes (C2) não é erro de montagem.

## 6. Crítica externa

Pedida ao agente `evo-devo` em 2026-10-17 (fase de *red team* da skill
`methodology-architect`); a resposta entra aqui como §7.

## 7. Crítica externa (agente `evo-devo`, 2026-10-17) e o que já foi feito com ela

| # | Crítica | Estado |
|---|---|---|
| a | **Piso de ruído no índice**: 3 genes a ~25 bits ÷ ~170 dão ~0,4 sem relíquia; *Gallus* marca 0,61. O "com dente" inclui a própria consulta; o contraste informativo é *Anolis* × *Chrysemys*, n = 1 × 1 | **Confirmado nos dados**: o IBSP (ósseo, funcional nas aves, tão desordenado quanto ENAM/AMBN) recupera só 23% (*Gallus*) e 37% (*Strigops*) do jacaré no `s27`. O índice em bits mede parte de cegueira do método para proteína desordenada. **As conclusões passam a se apoiar em lesões (nível 3), não no índice** |
| b | Faltam controle de composição (SPP1/IBSP/DMP1/MEPE), nulo de decaimento neutro por linhagem, segunda consulta | pendente — §8 |
| c | Calibrar o miniprot em **identidade equivalente**, não só em crocodiliano | **Feito** (`s29`/`s30`): *Python* 0,36, *Podarcis* 0,38, *Anolis* 0,39 de identidade no ENAM — mesma faixa do kakapo (0,35) — e zero lesão. Em ~12 mil aa alinhados em 6 espécies com dente, 1 lesão (AMBN do *Anolis*): ~0,1 por mil aa, contra 10–21 por mil nas relíquias |
| c | Coincidência na coordenada da proteína do jacaré **não prova** homologia da lesão; exige alinhamento ave × ave em nucleotídeo, mesmo indel, fora de homopolímero; GRCg7b × GRCg7w; paleognata | paleognatas incluídas (`s30`); o teste em nucleotídeo é o próximo passo |
| d | Afirmações excessivas em §A10/§A11 | corrigidas em `PESQUISA-2026-10-12` §A12 |
| e | Simples: multiz100way do hg38 + `mafsInRegion` (já instalado). Ambicioso: TOGA com o jacaré como referência | pendente — §8 |

## 8. Próximos passos, em ordem de ganho de precisão

1. **Nível 4 em nucleotídeo** (s31): relíquia de ENAM/AMBN de uma ave → blastn nas outras aves → MAFFT ave × ave + CDS do jacaré → indels e stops compartilhados por coluna, fora de homopolímero. Gabarito: as lesões compartilhadas de Meredith 2013 (PMID 23342979).
2. **Nulo por janela**: consultas embaralhadas que preservam composição (`rng = np.random.default_rng(0)`), mesmo tBLASTn; índice reportado só acima do percentil 95 do nulo.
3. **Controle de composição**: SPP1, IBSP, MEPE pelo mesmo pipeline; o "decaimento" do esmalte só conta na parte que excede o dos SCPP ósseos.
4. **Segunda consulta** (*Anolis*) para tirar a circularidade do jacaré.
5. multiz100way / TOGA como validação independente.

## 9. Resultado dos pendentes (2026-10-17, segunda rodada)

**Nível 4 em nucleotídeo** (`s31_aves_nt.sh`, `s31_nivel4.py`). A relíquia de
ENAM do kakapo, como consulta em nucleotídeo (dc-megablast), acha o trecho
homólogo em pato, mandarim e nas três paleognatas (66–73% de identidade); frango
e tinamu: nada. Alinhamento múltiplo com três crocodilianos (quadro e
polaridade) em MAFFT L-INS-i **e** MUSCLE 5; conta só lesão presente nos dois,
fora de homopolímero, com o mesmo tipo, posição e tamanho. Datação = menor clado
que contém as aves com a lesão, **exigindo** que toda ave que cobre a posição
dentro desse clado também a tenha.

| Clado | Lesões congruentes |
|---|---|
| Palaeognathae | 4 stops (códons 158, 389, 472, 539) |
| Novaeratitae (*Apteryx* + *Dromaius*) | 3 stops + deleção de 10 nt |
| Psittacopasserae (*Strigops* + *Taeniopygia*) | 1 stop |
| Neornithes | **nenhuma** — 3 candidatas incongruentes com a árvore (homoplasia; 2 são TGG→TGA, 1 usa stops diferentes) |

Controle: zero lesão nos crocodilianos nos dois alinhadores. Homoplasia medida:
3 de 12 stops compartilhados contradizem a árvore — stops compartilhados valem
menos que indel compartilhado.

**Nulo por janela e segunda consulta** (`s33_indice_nulo.py`, `s33_alinhar.sh`):
20 embaralhadas por proteína (semente 20261017); conta só o excesso sobre o
percentil 95 do nulo (26–32 bits em todas as janelas).

| Consulta | Com dente | Sem dente |
|---|---|---|
| jacaré | crocodilianos 12,3–15,0; escamados **1,45**–2,45 | tartarugas 0–**2,20**; aves 0–0,71 |
| *Anolis* | crocodilianos **2,29**–2,84; escamados 4,2–6,6 | tartarugas 0–0,57; aves **0–0,10** |

Com a consulta do jacaré, **há sobreposição** (tartaruga-verde 2,20 > píton
1,57): a separação limpa do s26 era artefato de amostragem de espécies. Com a
consulta do *Anolis*, crocodilianos, tartarugas e aves estão à **mesma
distância** filogenética da consulta (Archelosauria × Lepidosauria), e a
separação reaparece: crocodilianos 2,29–2,84 contra aves 0–0,10. Controle de
taxa por linhagem: com a mesma consulta, o LPL funcional rende 162–172 bits
nas aves e 149–151 nos crocodilianos.

**Ainda pendente**: controle de composição com SPP1/IBSP/MEPE pela consulta do
*Anolis*; multiz100way/TOGA; ODAM/ODAPH (agente `evo-devo` rodando).

**Controle de composição** (`s34_composicao.py`, `s34_alinhar.sh`): IBSP do
*Anolis* — sialoproteína óssea, desordenada como ENAM/AMBN, funcional nas aves —
contra os intervalos do DSPP do s27 (onde o IBSP é âncora), mesmo nulo do s33.
Excesso sobre o p95: crocodilianos 52–57, aves 39–41 (**~73%** do nível dos
crocodilianos), *Chrysemys* 49. O ENAM, pela mesma consulta: crocodilianos
229–254, aves 0–1,7 (**~0%**). O tBLASTn perde ~¼ da sensibilidade nas aves
para essa classe de proteína; a queda do ENAM não cabe nisso.

**Adiado, com motivo**: multiz100way do hg38 exige o MAF do chr4 (vários GB) ou
uma sessão do Table Browser; TOGA exige cadeias lastz genoma × genoma. Os dois
validariam por alinhamento de genoma inteiro o que o s29/s31 já mostram por
outro caminho; ficam para quando houver espaço em disco e tempo de CPU para isso.

## 10. Busca no genoma inteiro (2026-10-17, terceira rodada)

`s35_genomas.sh` (frango, kakapo, jacaré em `~/genomas` no WSL),
`s35_busca_genoma.sh/.py`, `s35_sintenia_odaph.py`, `s35_regua_odaph.sh`.

**O e-value não serve em genoma inteiro.** A 1e-2, só a proteína do próprio
jacaré se achou; o mesmo bitscore que dava e = 3e-6 numa janela passa de 1e-2
num banco de 2,3 Gb. Decisão por bitscore, como no s26.

**ODAPH**: no jacaré, as consultas de outras espécies (lagartixa, humano,
gambá) **não acham o ODAPH do próprio jacaré** no genoma inteiro. Nas duas aves
o melhor hit da ODAPH do jacaré é recíproco a ela — frango NC_052546.1:3,02 Mb
(chr15), kakapo NC_044279.2:31,6 Mb — mas a 32–33 bits, e os dois loci **não
têm um vizinho em comum** (±400 kb; pontos aleatórios do controle também zero):
não são o mesmo locus em duas linhagens; é recíproco por acaso no piso.

**Régua**: a ODAPH do jacaré acha a ODAPH funcional do *Anolis* com 36,6 bits e
a do *Podarcis* com 26,9 — a ~280 Ma, um gene **funcional** cai na faixa do
ruído. Busca de sequência no genoma inteiro não tem poder para o ODAPH a
distância das aves. A evidência que vale é **posicional** (s32): em 7 aves o
CDKL2 encosta no RCHY1, sem o espaço onde ODAPH e THAP6 ficam no jacaré.

**ODAM**: a consulta de *Xenopus* acha dezenas de trechos do chr1 do jacaré a
34–38 bits; o recíproco de cada um é uma proteína qualquer (FERM, TANC1,
IFT140, ...). Composição de baixa complexidade, não ODAM. Sem referência
positiva em saurópsido com dente, o ODAM sai do painel.

## 11. Validação por alinhamento de genoma inteiro — multiz100way (`s36_multiz`)

O Table Browser do UCSC bloqueia acesso automatizado (Cloudflare Turnstile) e o
PECAN de amniotas do Ensembl só traz primatas na região; foi o MAF do chr4
inteiro (4,6 GB), recortado nos 29 éxons codificantes de ENAM/AMBN/AMTN do hg38
(73 blocos).

- **Frameshift pelo multiz: descartado.** O controle reprova — o ENAM funcional
  do jacaré sai com 54 contra o humano.
- **Cobertura (presença)**: jacaré 53 blocos (ENAM 86%), *Anolis* 29 (74%);
  *Chelonia* 38 (56%), *Chrysemys* 19 (63%), *Pelodiscus* 0; das **14 aves** do
  multiz100way, **9 sem nenhuma base alinhada** (galGal4, melGal1, anaPla1,
  taeGut2, melUnd1, colLiv1, geoFor1, pseHum1, araMac1) e 5 com fragmentos
  (amaVit1 17 blocos, ENAM 41%; falcões, *Ficedula* e *Zonotrichia* 3–4).

Método independente (cadeias lastz + multiz, humano como referência), mesmo
padrão: répteis com dente > tartarugas > aves, e o frango sem nada — coerente com
o miniprot (s29), que também não alinhou o ENAM do frango.
