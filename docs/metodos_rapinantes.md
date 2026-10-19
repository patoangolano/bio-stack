# Métodos — posição filogenética das aves de rapina em Telluraves

Gerado em 2026-10-18 a partir de `runs/provenance.jsonl` (1.529 registros entre
2026-10-18T01:04Z e 07:01Z: 569 comandos de shell, 482 execuções de análise, 189
edições de código, 11 consultas a bases, 6 chamadas ao servidor do Claude Science,
267 outras e 5 marcos de sessão), dos scripts em `scripts/rapinantes/` e
`runs/rapinantes/ancestral/scripts/`, dos logs de cada braço e das versões lidas no
ambiente. Duas linhas antigas do log não são JSON válido e foram ignoradas. Os
horários estão em UTC. Os relógios do Windows e do WSL estão em America/Sao_Paulo
(UTC−3), e os horários locais dos logs foram convertidos. Todos os acessos a bases
de dados ocorreram em 2026-10-18 (UTC), ou seja, na noite de 2026-10-17 no horário
local.

O trabalho teve quatro braços computacionais: (1) reanálise das árvores de genes
publicadas pelo B10K e síntese das árvores depositadas na Open Tree of Life; (2)
posicionamento de fósseis em matrizes morfológicas publicadas; (3) reconstrução do
estado ancestral do hábito rapinante sob topologias alternativas; e (4) um conjunto
filogenômico novo, com 57 táxons e 300 locos. As seções 4.7–4.11 do braço 4 foram
acrescentadas em 2026-10-19, depois do fim das execuções, a partir dos logs e dos
scripts que elas citam. Uma síntese de literatura (NotebookLM, agentes do Claude Science e buscas
na Europe PMC e na web) orientou a escolha das hipóteses. Ela está em
`research_notes/Posição filogenética das aves de rapina/` e não é descrita aqui.

Os rótulos de topologia dos braços 1 e 3 são independentes. No braço 1, T1–T8 são
oito topologias fixas de Telluraves pontuadas pelo ASTRAL. No braço 3, T1, T1b, T1c
e T2–T4 são seis árvores-fonte (Stiller, Prum, Jarvis, Kuhl).

## Registro e reprodutibilidade

Todas as chamadas de ferramenta foram registradas automaticamente, com parâmetros e
horário UTC, em `runs/provenance.jsonl`. Cada etapa corresponde a um script
versionado, e a ordem e o comando de execução estão no cabeçalho de cada um. As
sementes foram fixas: 12345 no ASTER, no IQ-TREE do braço 1 e do pipeline, e na
amostragem de quartetos da Open Tree (`random.Random(12345)`); 20261017 no IQ-TREE
e nos testes de permutação do braço 2 e em todo o braço 3; 0 na permutação dos
genes candidatos do braço 4 (`numpy.random.default_rng(0)`). Não usamos estado
global de gerador em Python. Em R, `set.seed` foi chamado dentro de cada processo,
antes de cada etapa estocástica.

## Ambiente computacional

As análises rodaram num notebook com processador Intel Core Ultra 7 255HX, sob
Windows 11 (build 26300) e WSL2 (Ubuntu 24.04, 20 núcleos lógicos e 15 GB de RAM).
A GPU não foi usada. As ferramentas de linha de comando do WSL foram instaladas por
pixi a partir de `bootstrap/pixi-global.toml`, com build e sha256 em
`bootstrap/pixi-global.snapshot.tsv`: IQ-TREE 3.1.4 (bioconda, h8471819_0); ASTER
1.25 (hc52dbad_1), que traz o ASTRAL-IV v1.25.4.8 e o wASTRAL v1.25.3.8; gotree
0.5.2; miniprot 0.18-r281; gffread 0.12.9; DIAMOND 2.2.8; SeqKit 2.14.0; TaxonKit
0.20.0, com a taxonomia do NCBI baixada em 2026-10-16T03:29Z; e NCBI Datasets
(`datasets`/`dataformat`) 18.38.0. Os scripts que rodaram no WSL sem ambiente
virtual usaram só a biblioteca padrão do Python 3.12.3.

No Windows, o Python 3.11.17 era gerenciado por uv, com numpy 2.4.6, pandas 2.3.3,
scipy 1.17.1, DendroPy 5.1.0, Biopython 1.88, matplotlib 3.11.2, requests 2.34.2,
openpyxl 3.1.5 e o cliente `morphosource` 1.2.0. O R 4.6.1 (Rscript.exe, Windows)
usou as bibliotecas do renv do projeto (`r/renv.lock`, renv 1.3.1): ape 5.8-1,
phytools 2.5-2, castor 1.8.7, phangorn 2.12.1, treeio 1.36.1, ggtree 4.2.0 e
ggplot2 4.0.3. O pipeline Nextflow do braço 4 rodou em contêineres Docker, com
versões próprias (seção 4.6).

Os braços rodaram em paralelo e dividiram a máquina. Por isso cada um teve um teto
de CPU: 8 threads no ASTER do braço 1; 4 (às vezes 2) threads no IQ-TREE do braço
2; 4 processos R no braço 3; 8 threads no miniprot; e 10 CPUs e 13 GB no pipeline
do braço 4. Essa concorrência coincidiu com três reinicializações da VM do WSL,
cuja causa não foi apurada, e com a morte por falta de memória descrita na seção
1.9.

## 1. Reanálise das árvores de genes do B10K e síntese da literatura depositada

### 1.1 Dados

As árvores de genes e as árvores de espécies de Stiller et al. (2024) vieram do
repositório oficial do estudo (doi:10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a;
arquivo ERDA da Universidade de Copenhague), por links `sid.erda.dk/share_redirect`.
O download ocorreu em 2026-10-18 entre 01:58 e 02:00 UTC. O índice
`published-files.json` lista 3.340 arquivos (397 GB). Baixamos 0,81 GB e nenhum
alinhamento:

- as árvores de máxima verossimilhança dos 63.430 locos intergênicos de 1 kb, com
  suporte aLRT, e a versão oficial colapsada em aLRT < 0,95;
- as árvores dos 14.972 éxons (1ª e 2ª posições do códon), dos 44.846 íntrons e das
  4.985 UCEs;
- as 11 árvores de espécies publicadas (63K ASTRAL, 63K concatenada RAxML-NG, 63K
  GOBS, éxons, íntrons, UCE, 94K, 128K, 153K e 159K);
- o teste de politomia do artigo, a árvore datada por MCMCtree e os dados-fonte da
  Fig. 2.

Conferimos os tamanhos contra o índice e guardamos o sha256 de cada arquivo em
`runs/rapinantes/reanalise/b10k/brutos_manifesto/`. O ERDA não publica md5.

### 1.2 Preparação e poda

O script `preparar_b10k.py` traduziu os nomes de folha para os 363 nomes da árvore
63K. Os éxons usam um código de 6 letras (GENESP), e 6 exceções foram resolvidas à
mão. As UCEs foram traduzidas pelo `uce_name_map.txt`. O rótulo de suporte
`)/0.999` foi reescrito como `)0.999`, porque o wASTRAL lê só o primeiro número e
aborta quando o rótulo começa com barra.

Como no artigo, ramos com aLRT < 0,95 foram colapsados. Nos éxons, íntrons e UCEs
usamos `gotree collapse support -s 0.95`. No intergênico usamos o arquivo colapsado
oficial.

As árvores foram então podadas, uma a uma, para 42 táxons focais
(`b10k/taxa_focais.tsv`):

- os 15 rapinantes do B10K;
- 2 Cariamiformes, 3 Psittaciformes e 3 Passeriformes;
- 13 Coraciimorphae;
- 6 grupos externos, com *Gallus gallus* como raiz de apresentação.

A poda usou `podar_arvores.py`, porque o `gotree prune` 0.5.2 aborta o arquivo
inteiro quando uma árvore fica com duas folhas. Árvores com menos de 4 táxons
focais foram descartadas. Ficaram 62.706 de 63.430 intergênicas, 14.271 de 14.972
éxons, 40.933 de 44.846 íntrons e 4.919 de 4.985 UCEs, com média de 39,2; 33,7;
28,3 e 39,1 táxons focais por árvore.

### 1.3 Frequências de quarteto

Rodamos o ASTER 1.25 (Zhang, Nielsen & Mirarab 2025) para cada tipo de dado
(`rodar_aster_b10k.sh`). Todas as execuções usaram `-u 3`, que grava q1–q3, f1–f3,
pp1–pp3 e o número efetivo de genes no rótulo e escreve `freqQuad.csv`. Todas
usaram também `--seed 12345`, `-t 8` e `--root Gallus_gallus`. As quatro análises
foram:

- **(a) pontuação da topologia publicada.** O ASTRAL-IV pontuou (`-C`) a topologia
  63K publicada, podada aos 42 táxons e sem comprimentos nem suportes, nas árvores
  colapsadas. Assim todos os tipos de dado foram medidos nos mesmos ramos.
- **(b) busca livre.** O ASTRAL-IV, sem peso, buscou livremente nas árvores
  colapsadas.
- **(c) wASTRAL ponderado.** O wASTRAL híbrido (`-L --mode 1`; Zhang & Mirarab
  2022), ponderado pelo aLRT, buscou livremente nas árvores não colapsadas.
- **(d) pontuação com 363 táxons.** O ASTRAL-IV pontuou a topologia 63K completa
  nas árvores colapsadas sem poda. Isso só foi concluído para as UCEs (seção 1.9).

### 1.4 Topologias alternativas e escore bruto

A frequência q que o ASTRAL imprime é uma média por quarteto de espécies. O
objetivo da busca, porém, é a soma bruta de quartetos. Por isso pontuamos oito
topologias fixas de 42 táxons nas árvores colapsadas de cada tipo de dado
(`pontuar_alternativas.sh`; ASTRAL-IV `-C -t 2 -u 1 --seed 12345`) e registramos o
escore total:

- T1, a 63K;
- T2, Cariamiformes + Falconiformes;
- T3, Cariamiformes + Psittacopasserae;
- T4, Strigiformes + Coraciimorphae;
- T5, Accipitrimorphae + Coraciimorphae;
- T6, T4 e T2 juntas;
- T7, Accipitrimorphae irmã do resto de Telluraves;
- T8, Accipitrimorphae + Strigiformes irmãs do resto de Telluraves.

Para separar o efeito do colapso dos ramos (`teste_ponderacao.sh`), repetimos a
pontuação de T1, T2 e T4 (`-u 2`) nas 1.209 árvores UCE que contêm os 42 táxons
focais, uma vez colapsadas e outra sem colapsar.

### 1.5 Testes por ramo

O script `sintese_b10k.py` aplicou três testes de qui-quadrado sobre as contagens
efetivas f1–f3 e o número efetivo de genes (EN) que o ASTER grava:

- o teste de politomia (Sayyari & Mirarab 2018), com 2 g.l., de f1–f3 contra EN/3;
- a igualdade das duas resoluções mais frequentes, com 1 g.l.;
- a simetria das duas resoluções menos frequentes, com 1 g.l.

A localPP é a do ASTER (Sayyari & Mirarab 2016). Os valores de p são nominais, sem
correção para testes múltiplos. O ASTER não tem teste de politomia embutido, então
o teste do artigo (ASTRAL 5.15.1, `-t 10`) foi lido do arquivo publicado.

### 1.6 Concordância

Calculamos o gCF (Minh, Hahn & Lanfear 2020) no IQ-TREE 3.1.4 (Wong et al. 2026) com
`iqtree3 -t <63K podada> --gcf <árvores não colapsadas podadas> -T 1 --seed 12345`
(`gcf_b10k.sh`). Não calculamos sCF nem sCFL, porque os alinhamentos não foram
baixados.

Também lemos a resolução enraizada dos quatro grupos (Accipitrimorphae,
Strigiformes, Coraciimorphae e Australaves) em cada uma das 11 árvores publicadas e
em cada árvore desta reanálise.

### 1.7 Open Tree of Life e síntese Aves 1.6

Consultamos a API v3 da Open Tree of Life (`api.opentreeoflife.org/v3`) em
2026-10-18, entre 02:18 e 02:21 UTC. A síntese era a opentree16.1 (criada em
2025-12-20), com taxonomia OTT 3.7draft3, 1.931 estudos e 2.064 árvores.

O script `opentree_nos.py` fez quatro coisas:

1. resolveu 51 exemplares no `tnrs/match_names`, no contexto Birds;
2. obteve a árvore induzida;
3. para cada um de 24 nós-chave, chamou `tree_of_life/mrca` e `node_info`;
4. conferiu se os táxons que deviam ficar fora do nó estavam fora.

Quando não estavam, o nó foi registrado como ausente da síntese. As citações vieram
de `studies/find_studies`.

A síntese Aves 1.6 (McTavish et al. 2025) foi baixada do ramo `main` do repositório
GitHub McTavishLab/AvesData em 2026-10-18 às 02:30 UTC (`aves_nos.py`). Os arquivos
foram `annotations.json`, `labelled_supertree_ottnames.tre`, `rank_collection.json`
e `summary_dated_clements.nex`. O commit não foi registrado.

Para a meta-análise (`opentree_votos.py`), baixamos da phylesystem as 28 árvores
filogenéticas que sustentam, conflitam ou resolvem algum nó-chave. As árvores
taxonômicas OTT, eBird e Clements ficaram de fora.

Cada folha recebeu um grupo pela taxonomia OTT (`taxonomy/subtree`), com duas
correções medidas:

- *Cerchneis* foi incluído em Falconidae;
- "Bubonidae" foi tratado como parte de Strigidae.

Folhas fósseis sem subgrupo foram excluídas.

Por árvore e por pergunta, registramos três coisas:

1. o critério de clado estrito;
2. a fração de 2.000 quartetos amostrados que sustenta cada resolução;
3. a topologia enraizada dos quatro grupos, em 2.000 amostras de 5 táxons.

### 1.8 Idades

As idades vieram da API REST do TimeTree 5 (Kumar et al. 2022), versão 20260820,
consultada em 2026-10-18 às 02:30 UTC (`timetree_pares.py`). Foram nove pares de
táxons. Os taxids do NCBI vieram do E-utilities com o filtro `Aves[Subtree]`. Para
cada par, n, mediana e quartis foram recalculados de `studies.time_estimates`.

Para comparação, lemos a idade do MRCA dos mesmos exemplares em duas árvores
datadas (`idades_arvores_datadas.py`):

- a árvore MCMCtree do B10K (HPD 95 %);
- a MCC da Aves 1.6.

### 1.9 Execução e falhas

As análises do ASTER rodaram no WSL entre 02:06 e 06:40 UTC. O journal do WSL
registrou reinicializações da VM às 03:52, 04:07 e 04:12 UTC, com o miniprot do
braço 4 e outras análises pesadas rodando em paralelo. A pontuação da topologia 63K
com os 363 táxons e as 63.430 árvores intergênicas (`-u 3`) coincidiu com as duas
últimas reinicializações.

Às 04:14 UTC relançamos essa pontuação como unidade transitória do systemd
(`systemd-run`), para que ela não dependesse da sessão do `wsl.exe`. Às 04:17 UTC o
OOM killer do kernel a matou, com pico de 14,2 GB de RAM mais 3,9 GB de swap. Ela
não foi repetida. Por isso:

- os valores intergênicos com 363 táxons citados vêm dos dados-fonte da Fig. 2a de
  Stiller et al. (2024);
- a pontuação com 363 táxons de éxons e íntrons, prevista na mesma unidade, foi
  abandonada.

O script pula etapas que já têm `freqQuad.csv`. A etapa restante (busca livre do
ASTRAL-IV e do wASTRAL no intergênico) foi relançada às 05:05 UTC numa unidade com
teto de memória (`-p MemoryMax=9G`) e terminou às 06:40 UTC.

## 2. Fósseis e matrizes morfológicas

### 2.1 Matrizes

Baixamos as matrizes de quatro fontes: a coleção NEXUS de G. T. Lloyd, o MorphoBank
(API pública), o figshare e os suplementos dos artigos. Isso ocorreu em 2026-10-18,
até 02:54 UTC. As URLs, as licenças e os DOIs estão em
`runs/rapinantes/fosseis/LEIAME_dados_brutos.md`, e os sha256 dos 28 arquivos estão
em `runs/rapinantes/fosseis/log/sha256_brutos.txt`. Os arquivos ficam em
`data/raw/rapinantes/morfologia/`.

Analisamos dez matrizes:

| ID | Fonte | Táxons × caracteres usados | Grupo externo |
|---|---|---|---|
| K19 | Ksepka, Grande & Mayr 2019 (via Lloyd) | 62 × 146 | *Eudromia elegans* |
| M11 | Mayr 2011 (via Lloyd) | 33 × 100 | Tinamidae |
| M05 | Mayr 2005 (via Lloyd) | 23 × 110 | Tinamidae |
| VZ26 | Mayr 2026 | 42 × 105 | Anhimidae |
| MK22 | Mayr & Kitchener 2022 (figshare) | 10 × 42 | Anhimidae |
| M16 | Mayr 2016 (via Lloyd) | 13 × 40 | Anhimidae |
| D21 | Degrange et al. 2021 (MorphoBank P3896) | 72 × 204 | *Meleagris gallopavo* |
| MC03 | Mayr & Clarke 2003 (via Lloyd) | 46 × 148 | *Apsaravis* |
| C25 | Chen et al. 2025 | 77 × 202 | *Gansus yumenensis* |
| LZ07 | Livezey & Zusi 2007 (Lloyd = MorphoBank P510) | 154 Neornithes × 2.325 caracteres variáveis | *Lithornis* |

A versão de Lloyd e a do MorphoBank P510 da matriz LZ07 têm 100 % das células
iguais. Três matrizes foram baixadas e não analisadas: Ksepka et al. 2017, contida
na K19; Mayr et al. 2003; e Degrange et al. 2015, que é interna a Phorusrhacidae.

Não encontramos matriz publicada com corujas-tronco, *Horusornis* ou Pandionidae
fósseis. A descrição de *Masillaraptor buchheimi* (Li et al. 2025) e as listas de
caracteres de K19, M11, M05 e MC03 ficaram atrás de paywall ou de verificação
anti-robô. O LEIAME não pôde ser gravado em `data/raw/`, porque o hook de proteção
bloqueia escrita lá (ver Pendências).

### 2.2 Codificação

O script `p02_preparar.py` aplicou as seguintes regras, registradas por matriz em
`matrizes/<ID>/decisoes.json`:

- polimorfismos e caracteres inaplicáveis ("−") viraram "?";
- os estados foram renumerados 0..k−1;
- caracteres invariantes após a recodificação foram removidos;
- caracteres que os autores ordenaram entraram não ordenados.

A família de cada terminal veio da taxonomia do NCBI via TaxonKit (`f02_taxonkit.sh`),
com escolhas manuais justificadas em `p01_taxonomia.py`.

### 2.3 Máxima verossimilhança

Usamos o IQ-TREE 3.1.4 com `-st MORPH` e partições por número de estados (k = 2, 3
e ≥ 4; partição com menos de 3 caracteres fundida à vizinha). Cada partição recebeu
o modelo MK+ASC+G4, com comprimentos de ramo proporcionais entre partições (`-p`).
Os suportes foram 1.000 réplicas de UFBoot com `--bnni` (Hoang et al. 2018) e 1.000
réplicas de SH-aLRT (Guindon et al. 2010). Também usamos `-wbt`, `--seed 20261017` e
`-T 4`; a D21 e a análise livre da LZ07 rodaram com `-T 2`. Os comandos exatos estão em
`log/comandos_iqtree.txt`.

Cada matriz foi analisada duas vezes:

- **análise livre**: a morfologia sozinha;
- **análise com esqueleto**: os táxons atuais presos por `-g` à árvore principal de
  Stiller et al. (2024), com os fósseis livres.

O esqueleto foi a árvore `63K.tre` (ASTRAL, 63.430 regiões intergênicas), podada a
um representante por terminal morfológico. Terminais que dividem o mesmo
representante viraram politomia. A MK22 tem só três táxons atuais, então não
recebeu esqueleto.

Na matriz LZ07, a análise livre foi interrompida em 2026-10-18 às 03:39 UTC. A
busca de ML já estava concluída (lnL −103.014,6), mas o refinamento do UFBoot
estava em 100 de 1.000 réplicas (`f09_lz07_corte.sh`). Por isso usamos só a árvore
de ML dessa análise. A análise com esqueleto da LZ07 rodou sem UFBoot.

### 2.4 Posição dos fósseis

Em cada uma das 1.000 árvores UFBoot, podamos os demais fósseis e registramos o
conjunto de táxons atuais irmão do fóssil (`p03_resumo_arvores.py`, DendroPy).
"Grupo total de X" é a fração de árvores em que esse conjunto está contido em X.

### 2.5 Testes AU

Para cada hipótese "fóssil irmão de X", o script `p04_au.py` gerou uma árvore de
restrição: o esqueleto, com o fóssil enxertado no ramo-tronco de X e os outros
fósseis livres. Para cada hipótese, rodamos uma busca de ML restrita (`-g`,
`--seed 20261017`).

As melhores árvores de todas as hipóteses, mais a árvore de ML livre, entraram num
único teste AU no IQ-TREE: `-te <árvore livre> -z <hipóteses> -n 0 -zb 10000 -zw
-au`, com 10.000 reamostragens RELL (`log/comandos_au.txt`).

Para os táxons atuais, comparamos o esqueleto molecular com duas variantes:

- Falconidae deslocado para junto de Accipitridae, a "Falconiformes" clássica;
- em D21 e LZ07, também Cathartidae junto de Ciconiidae.

O teste final comparou só variantes do esqueleto, porque a árvore morfológica livre
domina qualquer conjunto. Na LZ07, as alternativas foram feitas por SPR sobre a
árvore restrita, sem reotimização (`p04_au.py lz07_spr`, `f08_lz07_au.sh`), o que
joga contra elas.

### 2.6 Homoplasia

Calculamos o IC e o IR de Fitch por caráter, não ordenado, com "?" = qualquer
estado (phangorn 2.12.1; `p06_homoplasia.R`). As comparações foram entre a árvore
de ML morfológica e a de ML com esqueleto, com e sem fósseis.

A classe "rapinante" veio das descrições publicadas dos caracteres:

- bico em gancho, cera e tômio;
- unguais, hálux, dedos e tendões do pé;
- tarsometatarso e hipotarso.

Essas descrições só existem para VZ26, MK22, C25, D21 e LZ07.

Testamos a diferença de médias do IC, e do ΔIC entre as árvores, entre caracteres
rapinantes e os demais. Usamos 10.000 permutações de rótulos (`set.seed(20261017)`),
estratificadas por matriz no conjunto agregado, com p unilateral. Também contamos
os passos por caráter na árvore molecular e na alternativa falcão + gavião. As
proporções foram comparadas por teste exato de Fisher.

### 2.7 Morfoespaço

O script `p07_morfoespaco.py` calculou as distâncias MORD e Gower por par de
táxons, só nos caracteres codificados nos dois. Pares com menos de 15 caracteres
comparáveis contaram como faltantes. Os táxons foram removidos um a um, começando
pelo que tinha mais pares faltantes, até a matriz de distâncias ficar completa.
Depois fizemos uma PCoA.

Nos táxons atuais de Telluraves, rodamos quatro testes:

- PERMANOVA por clado molecular, por guilda (rapinante ou não) e com o termo
  sequencial guilda | clado;
- Mantel parcial;
- a distância média entre rapinantes de clados diferentes, contra sorteios do mesmo
  número de táxons.

Os testes usaram 9.999 permutações, `numpy.random.default_rng(20261017)` e
implementação própria em numpy.

### 2.8 Inventário de micro-TC

O MorphoSource foi consultado pela API pública, com o cliente `morphosource` 1.2.0,
em 2026-10-18, entre 02:31 e 02:33 UTC (`p05_morphosource.py`). Foram 56 consultas
(48 nomes de fósseis e 8 famílias atuais), que devolveram 604 mídias, registradas em
`microct/inventario_morphosource.tsv`. Nenhuma mídia foi baixada.

## 3. Reconstrução do estado ancestral

### 3.1 Árvores

Construímos árvores em nível de família, com 202 famílias e uma espécie
representante por família, tiradas das 363 espécies de Stiller et al. (2024) e
montadas em `02_arvores.R`. Galloanserae foi a raiz. Os grupos externos foram
Strisores, Aequornithes, Charadriiformes e Columbaves.

Usamos seis topologias:

- T1: árvore datada 63K por MCMCtree (`04_timetree/main_alternative_mod.tre`);
- T1b: ASTRAL 94K;
- T1c: ASTRAL de éxons (1ª + 2ª posições);
- T2: Prum et al. (2015);
- T3: Jarvis et al. (2014);
- T4: Kuhl et al. (2021).

T1, T1b e T1c são podas diretas das árvores de Stiller, baixadas do ERDA (`erda.ku.dk`)
em 2026-10-18 entre 01:57 e 02:04 UTC, com md5 em `downloads/MANIFESTO_downloads.tsv`.

T2 e T3 vieram das cópias redistribuídas no arquivo de Stiller, na pasta da Extended
Data Fig. 6a: `Avian-TimeTree.tre`, com 198 espécies, e
`Chronogram01.TENT.ExAML.tre`, com 48. Não foram conferidas contra os repositórios
originais.

O `TREES.tar.gz` de Kuhl et al. (Dryad doi:10.5061/dryad.ngf1vhhpx) não pôde ser
baixado, porque a API exige token e a página serve um desafio de prova de trabalho
em JavaScript. Por isso T4 foi codificada à mão a partir da Fig. 2 do artigo, em
nível de ordem (PMC7783168, `msaa191f2.jpg`). A Fig. 3A, em nível de família, não
tem resolução legível no PMC.

T2, T3 e T4 são híbridas. A espinha entre ordens veio da fonte, e as relações dentro
de cada ordem e dentro dos grupos externos vieram de T1.

Cada topologia recebeu dois conjuntos de comprimentos de ramo:

- **(a) idades de T1.** Nós cujo conjunto de famílias forma clado em T1 receberam a
  idade de T1. Os demais foram espaçados uniformemente, à maneira do BLADJ (Webb,
  Ackerly & Kembel 2008). A raiz ficou em 88,43 Ma.
- **(b) Grafen.** Usamos ρ = 1 (Grafen 1989), escalado para a mesma idade de raiz.

Nas árvores datadas enxertamos quatro fósseis como pontas, em dois cenários (F1 e
F2), com idades, posições e codificações em `caracteres/fosseis_codificacao_manual.tsv`
e `arvores/enxertos_fosseis.tsv`:

- *Primoptynx poliotauros* (Mayr, Gingerich & Smith 2020), no tronco de Strigiformes;
- Masillaraptoridae (Mayr & Kitchener 2022), no tronco de Falconidae;
- Messelasturidae (Mayr 2021, 2025), no tronco de Psittacopasserae (F1) ou de
  Strigiformes (F2);
- Phorusrhacidae, representados por *Andalgalornis* (Degrange et al. 2010; Tambussi
  et al. 2012; Angst et al. 2013), no tronco de Cariamidae.

O ponto de enxerto ficou a meio caminho entre o nó-pai do ramo-alvo e o maior valor
entre a idade do nó-filho, a idade da ponta e a idade mínima da linhagem.

Ao todo foram 24 árvores: 6 topologias × (datada, Grafen, F1, F2).

### 3.2 Caracteres

A dieta veio do AVONET (Tobias et al. 2022), arquivo 34480856 do figshare
(doi:10.6084/m9.figshare.16586228.v7), com md5 conferido contra o declarado. Também
usamos o EltonTraits 1.0 (Wilman et al. 2014; `BirdFuncDat.txt`). Os dois foram
baixados em 2026-10-18 às 02:11–02:14 UTC.

O script `01_caracteres.py` (uv, pandas) pareou as espécies de Stiller com a
taxonomia eBird/Clements do AVONET: 312 pareamentos por nome exato, 15 por tabela de
correspondência, 2 por gênero e 1 por sinônimo manual.

Os dez caracteres foram:

- **C1 (dieta de vertebrados).** É a proporção p de espécies da família com
  `Trophic.Niche` em {Vertivore, Scavenger}, codificada pelos limiares T50
  (0,10/0,50; principal), T25 (0,05/0,25), T75 (0,25/0,75) e bin50 (sem
  polimorfismo).
- **C2.** C1 + "Aquatic predator", pelo limiar T50.
- **C3.** EltonTraits, com Vend + Vect + Vunk + Scav ≥ 50 %, pelo limiar T50.
- **C4 e C4b.** Bico em gancho, nas versões estrita e ampla.
- **C5.** Pés rapinantes.
- **C6.** Síndrome de 3 estados, derivada de C1_T50, C4 e C5.

Um polimorfismo entrou como prior de ponta (½, ½). As células de bico e pés foram
codificadas pelo analista. Só parte delas tem fonte primária específica (Fowler et
al. 2009; Portugal et al. 2016). A fonte de cada célula está em
`caracteres/caracteres_familia.tsv`.

### 3.3 Modelos, raiz e seleção

Nos caracteres binários ajustamos quatro modelos Mk: ER, ARD, irreversível
só-perda (1→0) e irreversível só-ganho (0→1). No C6 ajustamos ER, SYM, ARD,
só-perda e só-ganho. A MLE veio de `castor::fit_mk` (Louca & Doebeli 2018), com 16
partidas aleatórias (`Ntrials = 16`; o valor era 8 nos testes e foi elevado antes
da rodada final). Depois cada modelo foi reavaliado no `phytools::fitMk`, com Q fixa
(Revell 2024). Se um modelo ficou abaixo de outro aninhado nele, recebeu a Q do
aninhado.

Usamos duas condições de raiz:

- (i) o prior de FitzJohn, Maddison & Otto (2009);
- (ii) a raiz de Neognathae fixada no estado 0, justificada pelo grupo externo:
  nenhuma família de Palaeognathae ou de Galloanserae é rapinante. Com essa raiz,
  o modelo só-perda fica inviável.

A seleção usou AICc, com n = número de pontas. As probabilidades marginais
(`phytools::ancr`) foram calculadas nos nós focais e promediadas pelos pesos de
AICc. As sementes foram `set.seed(20261017 + 10·id)` antes dos ajustes, em cada
processo-filho, com o RNG padrão do R.

### 3.4 Mapeamento estocástico

Rodamos o `phytools::make.simmap` (Huelsenbeck, Nielsen & Bollback 2003) com 1.000
mapas por condição e Q fixa na MLE. Ele rodou para o melhor modelo, nas duas
raízes, e para o ER como referência, na raiz de FitzJohn. Os caracteres foram
C1_T50, C4, C5 e C6, nas 18 árvores datadas: 144 tarefas com mapeamento. A semente
foi `set.seed(20261017 + 10·id + j)`, com j = 1 para o melhor modelo e j = 2 para
o ER, e cada valor está em `resultados/tabela_transicoes.tsv`.

### 3.5 Sensibilidade e parcimônia

Os testes de sensibilidade estão embutidos no desenho:

- os limiares T25, T75 e bin50;
- o caráter C2, com predadores aquáticos;
- o caráter C3, do EltonTraits;
- o bico amplo (C4b);
- as árvores de Grafen;
- os cenários de fósseis F1 e F2.

O script `03b_mk_fosQ.R` repetiu C1_T50, C3 e C6 nas 12 árvores com fósseis, com a
dieta dos fósseis codificada como desconhecida. Foram 72 tarefas, sem mapeamento.

Como referência, rodamos uma parcimônia de Sankoff (`phangorn::sankoff`,
`ancestral.pars` com `type = "MPR"`) com custos iguais e com ganho = 2× perda. A
orientação da matriz de custo foi conferida num exemplo de brinquedo.

### 3.6 Execução

O `03_mk.R` montou 480 tarefas: 24 árvores × 10 caracteres × 2 raízes. Elas rodaram
num cluster PSOCK de 4 processos (`parLapplyLB`), com as tarefas de mapeamento na
frente da fila.

A primeira rodada começou em 2026-10-18 às 03:10 UTC e foi interrompida às 04:01
UTC, com 325 de 480 tarefas prontas (`logs/03_mk_run.log`). O motivo foi
desequilíbrio de carga: sem `chunk.size`, o `parLapplyLB` entregava blocos
contíguos e um núcleo ficava com todos os mapeamentos pesados.

O script foi alterado em dois pontos: `chunk.size = 1`, para que cada tarefa fosse
uma unidade de escalonamento, e a retomada, que pula tarefas cujo `.rds` já existe.
Antes de retomar, conferimos se havia `.rds` corrompido. A rodada foi retomada às
04:02 UTC e terminou as 155 tarefas restantes em 2,32 h, com 480 de 480 tarefas
concluídas e 0 erros (`logs/03_mk_run2.log`).

Como a semente é função do id da tarefa, refazer qualquer tarefa daria o mesmo
resultado. As 72 tarefas de dieta fóssil desconhecida rodaram em seguida, das 06:21
às 06:23 UTC, com 0 erros (`logs/03b_mk_fosQ_run.log`). O `sessionInfo()` de cada
rodada está em `logs/`.

## 4. Conjunto filogenômico novo (57 táxons, 300 locos)

### 4.1 Táxons e montagens

As montagens disponíveis por clado foram levantadas com
`datasets summary genome taxon` (`f01*.sh`). O conjunto final tem 57 táxons, com
nome, acessão, clado, nível de montagem e N50 em
`runs/rapinantes/filogenomica/taxons.tsv`.

Desses, 46 têm anotação RefSeq (NCBI Eukaryotic Genome Annotation). Houve uma troca:
*Phaethon aethereus* (GCF_964289735.1) não tinha dados no NCBI Orthologs e foi
substituído por *P. lepturus* (GCF_000687285.1). Os outros 45 foram conferidos antes
do download com o gene-sonda BSG (`f02c`).

Onze táxons vieram de montagens GenBank sem anotação RefSeq, das quais extraímos os
CDS com o miniprot. As 11 montagens (4,1 GB comprimidos) foram baixadas com
`datasets download genome accession --include genome` em 2026-10-18, entre 02:17 e
03:04 UTC. Os táxons e as proteínas de consulta foram:

| Táxon | Montagem | Consulta (parente RefSeq) |
|---|---|---|
| *Sarcoramphus papa* | GCA_037962945.1 | *Gymnogyps californianus* |
| *Sagittarius serpentarius* | GCA_034782735.1 | *Aquila chrysaetos* |
| *Pandion haliaetus* | GCA_026930845.1 | *Aquila chrysaetos* |
| *Elanus axillaris* | GCA_026212715.1 | *Aquila chrysaetos* |
| *Gypaetus barbatus* | GCA_028022735.1 | *Aquila chrysaetos* |
| *Leptodon cayanensis* | GCA_055772235.1 | *Aquila chrysaetos* |
| *Aegypius monachus* | GCA_982185325.1 | *Aquila chrysaetos* |
| *Pithecophaga jefferyi* | GCA_054553165.1 | *Aquila chrysaetos* |
| *Herpetotheres cachinnans* | GCA_013399355.1 | *Falco peregrinus* |
| *Ibycter americanus* | GCA_035590275.1 | *Falco peregrinus* |
| *Ninox strenua* | GCA_982290615.1 | *Strix aluco* |

### 4.2 Seleção de locos e ortologia

O relatório NCBI Gene do *Gallus gallus* foi baixado com
`datasets summary gene taxon 9031` em 2026-10-18 às 02:03 UTC e tinha 28.263 genes.
O script `f02d_candidatos.py` reduziu esse conjunto a 12.677 candidatos com os
seguintes critérios:

- codificante, com grupo "NCBI Ortholog";
- símbolo oficial, fora dos `LOC`;
- cromossomo nomeado, fora de W, MT e Un;
- fora de famílias multigênicas e de genes com muitas retrocópias (expressões
  regulares no script).

Os candidatos foram embaralhados numa permutação fixa (`default_rng(0)`). Os
primeiros 1.200 foram pedidos ao NCBI Orthologs com `datasets download gene gene-id
--ortholog <46 tax IDs> --include cds,protein,product-report`, em lotes de 100
genes, entre 02:16 e 03:48 UTC (`f02e`). Os lotes que o servidor cortou com erro
HTTP/2 foram refeitos em sublotes menores (`f02e2`, `f02e3` e `f02g`; sublotes de
20 genes, até 4 tentativas). Um sublote de 20 genes falhou quatro vezes e saiu do
pool por motivo técnico. Ficaram 1.180 genes.

O script `f02f_locos_refseq.py` aplicou as regras por gene e táxon:

- um registro por táxon, o CDS mais longo;
- um táxon com mais de um gene no grupo de ortólogos saiu do loco;
- o CDS tinha de ser múltiplo de 3, sem stop interno, e ter entre 0,5× e 1,5× a
  mediana do loco.

Um loco foi aprovado com três condições: pelo menos 40 dos 46 táxons RefSeq, CDS do
galo com pelo menos 600 pb e proteína de consulta disponível. Passaram 686 genes.
Saíram 417 por terem menos de 40 táxons, 74 por CDS do galo ausente ou curto e 3 por
falta de consulta. A auditoria está em `auditoria_refseq.tsv`.

A ortologia dos táxons RefSeq vem dos grupos 1:1 do NCBI Orthologs, mais o filtro
de um gene por táxon. A dos táxons miniprot é o alinhamento de Rank 1 da proteína
de consulta no genoma, sem teste de melhor hit recíproco.

### 4.3 Extração por miniprot e auditoria

Rodamos o miniprot 0.18-r281 (Li 2023) com `-M 2 -t 8 --gff`, uma montagem por vez
(`f03b_miniprot.sh`). Removemos as linhas `##PAF` do GFF3 e extraímos o CDS com o
gffread 0.12.9 (`-x`).

O valor `-M 2` substituiu o padrão (`-M 1`) depois que a VM do WSL reiniciou
durante o miniprot do *Pandion*, às 03:52 UTC, com pico medido de 10,2–11,5 GB
numa VM dividida com outras análises. No *Sagittarius*, `-M 2` baixou o pico para
6,5 GB e manteve 684 de 686 alinhamentos de Rank 1 no mesmo lugar que `-M 1`. As
11 montagens foram rodadas de novo, todas com `-M 2`, entre 04:55 e 05:00 UTC, com
pico de 4,7–9,1 GB por montagem.

O script `f03c_miniprot_cds.py` aplicou os filtros. O registro saiu do loco, e
virou dado faltante, nos seguintes casos:

- frameshift ou stop interno;
- cobertura da consulta < 50 %;
- identidade < 70 %;
- CDS que não traduz limpo.

Nenhum desses casos foi tratado como lesão genética. Um segundo alinhamento com
escore ≥ 0,99 do primeiro foi registrado como ambíguo, e o loco ficou com o Rank 1.
O resultado por táxon está em `auditoria_miniprot_resumo.tsv`: de 526 registros
aceitos em 686 (*Pandion*) a 668 em 686 (*Sarcoramphus*).

### 4.4 Controle na mesma montagem

*Gypaetus* e *Sarcoramphus* têm anotação NCBI de GenBank (GB_2026_01_30 e
GB_2025_12_03), cujos proteomas baixamos (`f03d`). As proteínas extraídas pelo
miniprot foram traduzidas com `seqkit translate -T 1` e comparadas a esses
proteomas com o DIAMOND 2.2.8 (`blastp --ultra-sensitive -k 1`; `f03e_validacao.sh`).
O critério foi identidade ≥ 98 % com cobertura ≥ 90 %. O resultado está em
`validacao_miniprot_vs_anotacao_GB.tsv`.

### 4.5 Matrizes de entrada

O script `f04_locos_finais.py` percorreu os genes na ordem da permutação e ficou
com os primeiros 300 que tinham pelo menos 49 dos 57 táxons (85 %). Registros do
miniprot com comprimento fora de 0,5×–1,5× a mediana do loco saíram. O pool foi
usado até a posição 554 de 1.200.

A matriz final tem 57 táxons × 300 locos, com 95,0 % das células preenchidas. Os
detalhes estão em `genes.tsv` e `ocupacao_taxon.tsv`. As proteínas por loco estão em
`data/proc/rapinantes/filogenomica/locos_aa/`.

A matriz de nucleotídeos das posições 1 e 2 do códon foi gerada por
`f05_dna12.py` às 05:26 UTC, a partir dos alinhamentos de aminoácidos e das máscaras
do ClipKIT da execução descrita abaixo, e está em `locos_dna12/`. O script:

1. substituiu cada resíduo pelo códon do CDS, depois de conferir que a tradução
   coincide com a sequência alinhada;
2. aplicou a mesma máscara de colunas do ClipKIT;
3. removeu a 3ª posição do códon.

A matriz DNA 1+2 foi analisada na execução `20261018-084205-rapinantes_DNA12`
(seção 4.8).

### 4.6 Pipeline de inferência (parâmetros como lançados)

A inferência usou o pipeline local `pipelines/filogenia` (Nextflow DSL2 com módulos
nf-core; Di Tommaso et al. 2017; Ewels et al. 2020). O pipeline foi lançado em
2026-10-18 às 05:02:10 UTC por `scripts/filogenia/filogenia.py`, na execução
`runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit/` (comando em
`executar.sh` e parâmetros em `pipeline_info/parametros.json`). A execução usou o
Nextflow 26.04.6, o Docker 29.9.1 (`-profile docker`) e o diretório de trabalho
`~/nf-work/filogenia/work`.

A configuração `scripts/rapinantes/filogenomica/limite_10cpus.config`, passada com
`-c`, limitou o executor local a 10 CPUs e 13 GB e baixou os pedidos de memória ao
uso medido. Ela não altera nenhum argumento das ferramentas. As versões abaixo são
as dos contêineres, e não as do pixi. Foram lidas do `.command.env` das tarefas e
da imagem registrada no `trace.tsv`, e conferidas em 2026-10-19 contra o
`pipeline_info/versoes.yml` de cada execução.

1. **Preparação.** PREPARA_LOCO, em Python 3.12.12, com `--tipo AA --nome-taxon
   primeiro`. O mínimo foi de 4 táxons por loco.
2. **Alinhamento.** MAFFT 7.520 (contêiner mulled com pigz 2.8; Katoh & Standley
   2013), em L-INS-i (`--localpair --maxiterate 1000 --threadit 0`, `--thread 4`).
   Todos os locos têm até 200 sequências. O `--threadit 0` torna o refinamento
   iterativo determinístico.
3. **Limpeza.** ClipKIT 2.11.4 (Steenwyk et al. 2020), com `--mode smart-gap
   --threads 2`.
4. **Árvore por loco.** IQ-TREE 3.1.4 com `-m MFP -B 1000 --bnni --alrt 1000
   --abayes --seed 12345 -o Gallus_gallus -st AA -T 1 -mem 1GB`. As fontes do
   método são o ModelFinder (Kalyaanamoorthy et al. 2017), com o critério padrão
   (BIC); o UFBoot2 (Hoang et al. 2018); o SH-aLRT (Guindon et al. 2010); e o
   aBayes (Anisimova et al. 2011). A memória foi de 1 GB por tarefa, contra 2 GB no
   padrão do pipeline, e o tempo máximo foi de 48 h.
5. **Supermatriz.** Uma partição por loco (`-spp`), com `-m MFP+MERGE` (fusão de
   partições; Chernomor, von Haeseler & Minh 2016; Lanfear et al. 2014, 2017), `-B
   1000 --bnni --alrt 1000 --abayes --seed 12345 -o Gallus_gallus -st AA` e `-T 8`
   (`--cpus-concatenado 8`).
   - Com 8 threads, a reprodução é de topologia e lnL. O UFBoot e o SH-aLRT podem
     variar 1–2 pontos entre execuções.
   - O módulo omite `-mem` com modelo particionado, porque o IQ-TREE o recusa nesse
     caso. O teto de memória é o do contêiner, 9 GB, e o tempo máximo é de 144 h.
   - Esta execução foi interrompida depois do ModelFinder. A árvore que entra nos
     resultados vem de uma segunda corrida com o mesmo esquema e os mesmos modelos
     fixados (seção 4.7).
6. **Concordância.** O gCF foi calculado das árvores de genes sobre a árvore da
   supermatriz (`-T 1`). O sCFL usou 100 quartetos (`--scfl 100`; Mo et al. 2023),
   com o modelo por partição escolhido na supermatriz (`-p`) e `--seed 12345`. Na
   execução AA com MAFFT, o sCFL rodou fora do Nextflow com `-T 8` (o módulo usa
   `-T 1`; seção 4.7).
7. **Árvore de espécies.** ASTER 1.25.3.8 (contêiner hc52dbad_1). O wASTRAL usou o
   preset de ponderação por aBayes (`-B`; Zhang & Mirarab 2022), `-t 2`, `--root
   Gallus_gallus` e `--seed 12345`, com mínimo de 3 locos. Só o campo aBayes do
   rótulo "SH-aLRT/aBayes/UFBoot" do IQ-TREE é passado ao wASTRAL.

### 4.7 Interrupção da supermatriz e nova corrida com modelos fixos

As árvores de genes e o wASTRAL da execução AA com MAFFT terminaram em 2026-10-18,
às 08:41 (horário de Brasília). O ModelFinder da supermatriz agrupou as 300 partições
em 48 grupos (BIC 4467486,907; lnL -2230458,528; 535 parâmetros), e a busca de árvore
terminou às 17:16. Às 21:09, durante o refinamento dos UFBoot, o processo morreu por
SIGHUP, quando a sessão que tinha lançado o Nextflow foi fechada.

A retomada pelo checkpoint foi abandonada. Com `-m MFP+MERGE`, o IQ-TREE 3.1.4 voltou
a testar as partições isoladas com o conjunto completo de modelos (~32 h previstas).
O cache `concatenado.model.gz` só guardava, para elas, os dois modelos da fase de
agrupamento.

A supermatriz que entra nos resultados foi refeita do zero, sem a busca de modelo:

- `modelos_48.py` reconstruiu do `concatenado.model.gz` o modelo de menor BIC de cada
  um dos 48 grupos do `best_scheme.nex`. A soma reproduz o lnL e os graus de
  liberdade do ModelFinder.
- A corrida usou `iqtree -p modelos_fixos.nex -B 1000 --bnni --alrt 1000 --abayes
  --seed 12345 -o Gallus_gallus -st AA` na mesma imagem do pipeline.
- Rodou com `-T 8` até a iteração 76 da busca e com `-T 18` daí em diante, retomando
  do próprio checkpoint. Antes da troca, um teste numa cópia confirmou que a retomada
  restaura a árvore, os 1000 UFBoot e os parâmetros sem recalcular nada.
- A corrida terminou em 2026-10-19, às 16:27, com lnL -2230233,545. A topologia é
  idêntica à da melhor árvore da busca interrompida (distância de Robinson-Foulds 0).

A publicação, a concordância e o resumo foram replicados fora do Nextflow por
`pos_concat_020210.sh`, com os contêineres e argumentos do pipeline. A única
diferença é o sCFL com `-T 8`. Os scripts estão em
`scripts/rapinantes/filogenomica/retomada_020210/`, e o registro detalhado em
`runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit/RETOMADA_MANUAL.md`. O
`metodos.md` gerado pelo RESUMO descreve o comando original (MFP+MERGE) e deve ser
lido com essa ressalva.

### 4.8 Matriz DNA 1+2

Execução `20261018-084205-rapinantes_DNA12`, com a matriz da seção 4.5 e `--alinhador
nenhum --limpeza nenhuma --tipo_seq DNA --concatenar false` (2 CPUs,
`limite_2cpus.config`). Foram feitas árvores de genes no IQ-TREE 3.1.4, com ModelFinder
de nucleotídeos e os mesmos suportes e semente da seção 4.6, e wASTRAL. Esta matriz
não tem supermatriz nem concordância. Como vem dos mesmos genes e do mesmo
alinhamento da execução AA, ela testa o efeito do tipo de dado, não replica a análise
de forma independente.

### 4.9 Sensibilidade ao alinhador

Execução `20261018-112050-rapinantes_AA_famsa_clipkit`: as mesmas proteínas, alinhadas
com FAMSA 2.4.1 (argumentos padrão do módulo nf-core) e limpas com ClipKIT
`smart-gap`. Por custo, as árvores de genes rodaram só com aBayes (`--ufboot 0 --alrt
0`), sem supermatriz (`--concatenar false`, 8 CPUs, `limite_8cpus.config`). O wASTRAL
usou os mesmos argumentos. A limpeza com trimAl não foi executada.

Nas três árvores de espécies, a monofilia e a probabilidade local dos clados
relevantes foram lidas, com a raiz em *Gallus*, por `f07a_clados_wastral.py`. O
grupo-irmão dos clados focais foi lido por `f07b_irmaos_wastral.py`.

### 4.10 Nós-chave, concordância e assimetria

`f08_nos_chave.py` lê cada nó-chave como bipartição, de modo que a raiz não importa,
e reúne quatro coisas:

- os suportes da supermatriz;
- gCF, gDF1, gDF2 e gDFP;
- sCF, sDF1 e sDF2 (sCFL);
- a probabilidade local do wASTRAL.

Os ramos do `.cf.stat` são casados com a árvore `.cf.tree` pelo rótulo e pelo
comprimento. A assimetria entre as duas alternativas NNI foi testada por
qui-quadrado (`scipy.stats.chisquare`) sobre gDF1_N × gDF2_N e sDF1_N × sDF2_N. Sob
ILS pura, as duas têm a mesma frequência esperada (Minh, Hahn & Lanfear 2020).

O script também conta as árvores de genes que contêm o clado, restrito aos táxons do
gene, com SH-aLRT ≥ 80 e UFBoot ≥ 95, e as que têm um ramo nesse critério
incompatível com ele.

Para saber qual alternativa NNI corresponde a gDF1 e qual a gDF2, as bipartições
concorrentes foram contadas diretamente nas árvores de genes. Para sDF1 e sDF2, o
sCFL foi calculado nas árvores restritas da seção 4.11 que forçam cada alternativa.

### 4.11 Testes de topologia

`f07_hipoteses.py` escreve, para cada hipótese, duas árvores:

- a árvore de restrição do IQ-TREE (`-g`), multifurcada, com os clados exigidos;
- uma árvore inicial compatível, obtida por SPR sobre a árvore ML da supermatriz: o
  clado X é podado e reenxertado como irmão de Y.

As 13 hipóteses são H1, H1b, H2, H3, H4, H4b, H5, H6, H7a, H7b, S1, S2 e S3, descritas
em `runs/rapinantes/filogenomica/topologia_AA/hipoteses.tsv`. Quando a árvore ML já
satisfaz uma hipótese, a melhor árvore restrita é a própria ML e a busca não roda.

`f07b_topologia.sh` faz as buscas restritas sobre a supermatriz AA:

- comando `iqtree3 -s supermatriz.fa -p concatenado.best_model.nex -g <restrição> -t
  <inicial> --fast --seed 12345 -T 18`;
- IQ-TREE 3.1.4 do pixi, a mesma versão do contêiner;
- o modelo particionado vem com os parâmetros estimados na supermatriz.

As árvores candidatas são a ML, a árvore do wASTRAL e as melhores árvores restritas.
Elas foram comparadas com `-z candidatos.nwk -n 0 -zb 10000 -au`, isto é, KH, SH e
AU [CITAÇÃO A CONFIRMAR] com 10 000 réplicas RELL. A execução foi lançada como tarefa
agendada por `lanca_topologia.ps1` e `roda_topologia.sh`.

## Uso de inteligência artificial

O desenho das análises, a escrita dos scripts, a execução e a redação deste texto
foram feitos com assistência do Claude Code (Anthropic), incluindo agentes
especializados (`phylo-analyst`, `evo-devo`), sob direção do autor. Nenhum dado
biológico foi gerado por IA. As árvores de genes e de espécies, as matrizes
morfológicas, os traços ecológicos e as montagens são dados públicos dos
repositórios citados. As exceções são a árvore T4, transcrita à mão de uma figura
publicada, e as células de bico e pés do braço 3, codificadas pelo analista. As
duas estão sinalizadas acima.

## Referências citadas nesta seção

Os DOIs abaixo foram conferidos na Crossref, na DataCite ou na Europe PMC durante o
trabalho, e estão registrados em `docs/FILOGENIA.md` e em `research_notes/`. Software
sem DOI conferido no repositório está marcado nas Pendências.

- Angst D et al. 2013. *PLoS ONE* 8:e80357. doi:10.1371/journal.pone.0080357
- Anisimova M et al. 2011. *Syst Biol* 60:685–699. doi:10.1093/sysbio/syr041
- Chen et al. 2025. *Integr Org Biol* 7. doi:10.1093/iob/obaf029
- Chernomor O, von Haeseler A, Minh BQ. 2016. *Syst Biol* 65:997–1008. doi:10.1093/sysbio/syw037
- Degrange FJ et al. 2010. *PLoS ONE* 5:e11856. doi:10.1371/journal.pone.0011856
- Degrange FJ et al. 2021. *Pap Palaeontol*. doi:10.1002/spp2.1361 (MorphoBank P3896, doi:10.7934/P3896)
- Di Tommaso P et al. 2017. *Nat Biotechnol* 35:316–319. doi:10.1038/nbt.3820
- Ewels PA et al. 2020. *Nat Biotechnol* 38:276–278. doi:10.1038/s41587-020-0439-x
- FitzJohn RG, Maddison WP, Otto SP. 2009. *Syst Biol* 58:595–611. doi:10.1093/sysbio/syp067
- Fowler DW, Freedman EA, Scannella JB. 2009. *PLoS ONE* 4:e7999. doi:10.1371/journal.pone.0007999
- Grafen A. 1989. *Phil Trans R Soc B* 326:119–157. doi:10.1098/rstb.1989.0106
- Guindon S et al. 2010. *Syst Biol* 59:307–321. doi:10.1093/sysbio/syq010
- Hoang DT et al. 2018. *Mol Biol Evol* 35:518–522. doi:10.1093/molbev/msx281
- Huelsenbeck JP, Nielsen R, Bollback JP. 2003. *Syst Biol* 52:131–158. doi:10.1080/10635150390192780
- Jarvis ED et al. 2014. *Science* 346:1320–1331. doi:10.1126/science.1253451
- Kalyaanamoorthy S et al. 2017. *Nat Methods* 14:587–589. doi:10.1038/nmeth.4285
- Katoh K, Standley DM. 2013. *Mol Biol Evol* 30:772–780. doi:10.1093/molbev/mst010
- Ksepka DT, Grande L, Mayr G. 2019. *Curr Biol* 29:657–663. doi:10.1016/j.cub.2018.12.040
- Kuhl H et al. 2021. *Mol Biol Evol* 38:108–127. doi:10.1093/molbev/msaa191 (Dryad doi:10.5061/dryad.ngf1vhhpx)
- Kumar S et al. 2022. *Mol Biol Evol* 39(8):msac174. doi:10.1093/molbev/msac174
- Lanfear R et al. 2014. *BMC Evol Biol* 14:82. doi:10.1186/1471-2148-14-82
- Lanfear R et al. 2017. *Mol Biol Evol* 34:772–773. doi:10.1093/molbev/msw260
- Li H. 2023. *Bioinformatics* 39(1):btad014. doi:10.1093/bioinformatics/btad014
- Livezey BC, Zusi RL. 2007. *Zool J Linn Soc* 149:1. doi:10.1111/j.1096-3642.2006.00293.x (MorphoBank P510, doi:10.7934/P510)
- Louca S, Doebeli M. 2018. *Bioinformatics* 34:1053–1055. doi:10.1093/bioinformatics/btx701
- Mayr G. 2005. *J Vertebr Paleontol* 25:635–645. doi:10.1671/0272-4634(2005)025[0635:TPOAPP]2.0.CO;2
- Mayr G. 2011. *J Syst Palaeontol* 9:159–171. doi:10.1080/14772019.2010.505252
- Mayr G. 2016. *J Paleontol* 90:357. doi:10.1017/jpa.2016.45
- Mayr G. 2021. *PalZ* 95:337–357. doi:10.1007/s12542-020-00541-8
- Mayr G. 2025. *J Morphol* 286:e70070. doi:10.1002/jmor.70070
- Mayr G. 2026. *Vertebr Zool* 76:121. doi:10.3897/vz.76.e182505
- Mayr G, Clarke J. 2003. *Cladistics* 19:527. doi:10.1016/j.cladistics.2003.10.003
- Mayr G, Gingerich PD, Smith T. 2020. *J Vertebr Paleontol* 40:e1769116. doi:10.1080/02724634.2020.1769116
- Mayr G, Kitchener AC. 2022. *J Vertebr Paleontol* 41(6):e2083515. doi:10.1080/02724634.2021.2083515 (figshare doi:10.6084/m9.figshare.20267628.v2)
- McTavish EJ et al. 2025. *PNAS* 122(18):e2409658122. doi:10.1073/pnas.2409658122
- Minh BQ, Hahn MW, Lanfear R. 2020. *Mol Biol Evol* 37:2727–2733. doi:10.1093/molbev/msaa106
- Mo YK et al. 2023. *Bioinformatics* 39(1):btac741. doi:10.1093/bioinformatics/btac741
- Portugal SJ et al. 2016. *Curr Biol* 26:R58–R59. doi:10.1016/j.cub.2015.12.004
- Prum RO et al. 2015. *Nature* 526:569–573. doi:10.1038/nature15697
- Redelings BD, Holder MT. 2017. *PeerJ* 5:e3058. doi:10.7717/peerj.3058
- Revell LJ. 2024. *PeerJ* 12:e16505. doi:10.7717/peerj.16505
- Sayyari E, Mirarab S. 2016. *Mol Biol Evol* 33:1654. doi:10.1093/molbev/msw079
- Sayyari E, Mirarab S. 2018. *Genes* 9:132. doi:10.3390/genes9030132
- Steenwyk JL et al. 2020. *PLoS Biol* 18:e3001007. doi:10.1371/journal.pbio.3001007
- Stiller J et al. 2024. *Nature* 629:851–860. doi:10.1038/s41586-024-07323-1 (dados: doi:10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a)
- Tambussi CP et al. 2012. *PLoS ONE* 7:e37701. doi:10.1371/journal.pone.0037701
- Tobias JA et al. 2022. *Ecol Lett* 25:581–597. doi:10.1111/ele.13898 (figshare doi:10.6084/m9.figshare.16586228.v7)
- Webb CO, Ackerly DD, Kembel SW. 2008. *Bioinformatics* 24:2098–2100. doi:10.1093/bioinformatics/btn358
- Wilman H et al. 2014. *Ecology* 95:2027. doi:10.1890/13-1917.1
- Wong TKF et al. 2026. *Mol Biol Evol* 43(5):msag117. doi:10.1093/molbev/msag117
- Zhang C, Mirarab S. 2022. *Mol Biol Evol* 39:msac215. doi:10.1093/molbev/msac215
- Zhang C, Nielsen R, Mirarab S. 2025. *Mol Biol Evol* 42(8):msaf172. doi:10.1093/molbev/msaf172

## Pendências

1. **Braço 4: conferência feita em 2026-10-19.** As versões da seção 4.6 batem com
   o `pipeline_info/versoes.yml` das três execuções. Os argumentos do wASTRAL foram
   conferidos contra o `.command.sh` real das três. Os da supermatriz são os da
   corrida com modelos fixos (seção 4.7).
2. **Análises do braço 4 ainda não rodadas.** Ainda faltam quatro análises:
   - a supermatriz e a concordância da matriz DNA 1+2;
   - a limpeza com trimAl (trimAl 1.5.0, configurado e não usado);
   - a supermatriz da execução FAMSA;
   - uma reanálise sem os ramos terminais anômalos (poda à maneira do TreeShrink).
3. **Citações a confirmar.** Faltam DOIs conferidos no repositório para ape,
   phangorn, DendroPy, DIAMOND, gffread, SeqKit, TaxonKit, NCBI Datasets, o modelo
   Mk para morfologia, o teste AU, a distância MORD, a distância de Gower, a
   PERMANOVA, o AICc, a parcimônia de Sankoff, R e Python. Todas estão hoje sem
   referência no texto: [CITAÇÃO A CONFIRMAR].
4. **Árvores T2–T4 do braço 3.** A T4 (Kuhl et al. 2021) foi codificada à mão de uma
   figura em nível de ordem. Ela deve ser conferida por uma segunda pessoa ou contra
   o `TREES.tar.gz` do Dryad, quando o download for possível. As árvores de Prum
   (T2) e de Jarvis (T3) vieram da redistribuição no arquivo de Stiller e não foram
   conferidas contra os repositórios originais (suplemento da *Nature*; GigaDB
   10.5524/101041).
5. **Commit da Aves 1.6.** Os arquivos foram baixados do ramo `main` do
   McTavishLab/AvesData em 2026-10-18, sem registro do commit. É preciso registrar o
   commit correspondente.
6. **LEIAME das matrizes morfológicas.** O hook `guard_data` bloqueou a gravação em
   `data/raw/`. É preciso copiar à mão `runs/rapinantes/fosseis/LEIAME_dados_brutos.md`
   para `data/raw/rapinantes/morfologia/LEIAME.md`.
7. **Fontes não obtidas no braço 2.** A descrição de *Masillaraptor buchheimi* (Li
   et al. 2025, doi:10.7302/27584) e as listas de caracteres de K19, M11, M05 e MC03
   ficaram atrás de paywall ou de verificação anti-robô. Valem uma tentativa por
   acesso institucional, porque sem as listas não se identificam os caracteres
   convergentes nessas matrizes.
8. **Divergência na nota de pesquisa do braço 4.** A
   `research_notes/.../computacional_filogenomica.md` (§3) diz que a supermatriz
   recebe `-mem 9G`. O módulo do pipeline omite `-mem` com modelo particionado, e o
   teto real é o do contêiner, 9 GB. Este texto segue o código. É preciso corrigir a
   nota ou confirmar no `.command.sh` quando a tarefa for submetida.
9. **Arquivo temporário.** Sobrou uma cópia descompactada de 1,25 GB em
   `~/rapinantes/raw/GCA_026930845.1/genoma.fna`, no WSL, que o hook de proteção de
   `raw/` não deixa apagar por ferramenta. É preciso removê-la à mão.
10. **Declaração de uso de IA.** O formato depende da revista ou da instituição de
    destino, e o texto atual é genérico.
