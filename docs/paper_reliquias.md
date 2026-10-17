# Relíquias dos genes do esmalte em aves e tartarugas: localização no GRCg7b, lesões compartilhadas e os limites da detecção por similaridade

% Matheus Lima
% [AFILIAÇÃO A PREENCHER]
% Manuscrito de trabalho, 17 de outubro de 2026. Não revisado por pares.

:::resumo Resumo
Aves e tartarugas perderam os dentes de forma independente, e os genes da matriz do esmalte deixaram restos nos dois genomas. Revisitamos esses restos em 18 montagens RefSeq de amniotas (3 crocodilianos, 3 escamados, 4 tartarugas, 8 aves) com uma escada de evidência explícita: ausência na anotação; ausência com controle positivo; lesão em quadro com coordenada; e lesão compartilhada que data a perda. No frango, relíquias codificantes de ENAM e AMBN ficam no cromossomo Z (GRCg7b, cerca de 55,16 Mb), no vão sintênico entre NRG1 e LPL, na mesma ordem do jacaré, sem anotação no NCBI nem no Ensembl. O índice de retenção por similaridade (bitscore de tBLASTn), comum em genômica comparativa de perda, mostrou-se confundido. Com o jacaré como consulta e um nulo por janela, espécies com e sem dente se sobrepõem (tartaruga-verde 2,20 contra píton 1,57). Com uma consulta equidistante (*Anolis*), a separação reaparece (crocodilianos 2,29–2,84; tartarugas 0–0,57; aves 0–0,10), mas o controle de composição mostra que o tBLASTn perde cerca de um quarto da sensibilidade nas aves para proteínas desordenadas. Pelo miniprot, as espécies com dente tiveram 1 lesão em 12.160 aminoácidos alinhados (0,08 por mil), contra 53 lesões em 3.518 aminoácidos (15,1 por mil) nas relíquias com lesão. Num trecho de 2,4 kb do ENAM, recuperamos lesões compartilhadas por Palaeognathae (4 stops), Novaeratitae (3 stops e uma deleção de 10 nt) e Psittacopasserae (1 stop), e nenhuma congruente em escala de Neornithes; 3 dos 11 stops compartilhados contradizem a árvore. Esses resultados replicam e localizam achados anteriores (Meredith et al. 2013) e mostram que índices de similaridade só servem como fenótipo contínuo depois de calibrados contra nulo, distância da consulta e composição.

**Palavras-chave:** pseudogene; perda gênica; esmalte; ENAM; AMBN; edentulismo; Archelosauria; genômica comparativa.
:::

:::resumo Abstract
Birds and turtles lost their teeth independently, and enamel matrix genes left remnants in both genomes. We revisited these remnants in 18 RefSeq amniote assemblies (3 crocodylians, 3 squamates, 4 turtles, 8 birds) using an explicit evidence ladder: absence from annotation; absence with a positive control; an in-frame lesion with coordinates; and a shared lesion that dates the loss. In chicken, coding remnants of ENAM and AMBN lie on the Z chromosome (GRCg7b, about 55.16 Mb), inside the syntenic gap between NRG1 and LPL, in the same order as in alligator, and annotated in neither NCBI nor Ensembl. A similarity-based retention index (tBLASTn bitscore), common in comparative gene-loss genomics, proved confounded. With alligator as query and a per-window null, toothed and toothless species overlap (green turtle 2.20 vs python 1.57). With an equidistant query (*Anolis*), separation reappears (crocodylians 2.29–2.84; turtles 0–0.57; birds 0–0.10), but a composition control shows that tBLASTn loses about a quarter of its sensitivity in birds for disordered proteins. With miniprot, toothed species showed 1 lesion in 12,160 aligned amino acids (0.08 per thousand), against 53 lesions in 3,518 amino acids (15.1 per thousand) in the remnants that carry lesions. In a 2.4-kb ENAM segment we recovered lesions shared by Palaeognathae (4 stops), Novaeratitae (3 stops and a 10-nt deletion) and Psittacopasserae (1 stop), and none congruent at the scale of Neornithes; 3 of 11 shared stops conflict with the tree. These results replicate and localize earlier findings (Meredith et al. 2013) and show that similarity indices are usable as continuous phenotypes only after calibration against a null, query distance and composition.

**Keywords:** pseudogene; gene loss; enamel; ENAM; AMBN; edentulism; Archelosauria; comparative genomics.
:::

## 1. Introdução

A ausência de dentes evoluiu várias vezes entre os vertebrados, inclusive em aves, tartarugas e algumas linhagens de mamíferos (Meredith et al. 2014). Nas aves modernas, mutações inativadoras compartilhadas em genes expressos no esmalte e na dentina indicam que o ancestral comum de Neornithes já não tinha dentes mineralizados, e a perda do esmalte foi estimada em cerca de 116 milhões de anos (Ma) (Meredith et al. 2014). Nas tartarugas a perda é independente: os crocodilianos, grupo-irmão das aves dentro de Archosauria, mantêm dentes, e a perda nas tartarugas ocorreu no Jurássico (Meredith et al. 2013).

Os genes da matriz do esmalte, sobretudo amelogenina (AMEL), ameloblastina (AMBN) e enamelina (ENAM), perderam a função nessas linhagens e passaram a decair. Sire et al. (2008) encontraram o pseudogene de AMEL da galinha e concluíram que AMBN e ENAM tinham desaparecido após rearranjos cromossômicos. Meredith et al. (2013) recuperaram cópias inativadas de AMBN, AMEL e ENAM em Galloanserae, Passeriformes e Psittaciformes, com frameshifts compartilhados anteriores à divisão basal de Neognathae, e restos de AMBN e ENAM na tartaruga-pintada (*Chrysemys picta*) com a sintenia original. Em mamíferos sem dente ou sem esmalte, o mesmo raciocínio estendeu a lista de genes específicos do dente a C4orf26, hoje ODAPH (Springer et al. 2016), e a ODAM (Springer et al. 2019).

Que a perda ocorreu já se sabe. O problema que tratamos aqui é de detecção. Uma busca por similaridade que não encontra um gene não distingue entre o gene ter sumido e o método ter ficado cego, e a cegueira cresce com a distância entre a consulta e o alvo, com a baixa complexidade da proteína e com o próprio decaimento. Esse problema afeta qualquer abordagem que associe perda gênica a perda de caráter em escala filogenética (Prudent et al. 2016; Kowalczyk et al. 2019; Saputra et al. 2021; Redlich et al. 2024), porque a variável de entrada (presente, ausente, ou um índice contínuo de retenção) herda os vieses da busca.

Tivemos três objetivos. (i) Localizar as relíquias de ENAM, AMBN e AMTN numa montagem atual do frango (GRCg7b) e em mais 16 espécies, com coordenadas. (ii) Medir quanto um índice de retenção por similaridade é confundido pela distância à consulta e pela composição da proteína. (iii) Graduar cada afirmação de perda numa escada de evidência explícita (Tabela 1), cujo nível mais alto é a lesão compartilhada que data o decaimento numa árvore.

**Tabela 1.** Escada de evidência usada no trabalho.

| Nível | Critério | Exemplo neste trabalho |
|---|---|---|
| 1 | gene ausente da anotação RefSeq | aves: 0 de 9 genes do painel anotados |
| 2 | busca não encontra, com controle positivo no mesmo limiar | ODAPH nas aves; âncoras encontradas |
| 3 | relíquia com lesão (frameshift ou stop) e coordenada | ENAM do kakapo, 6 frameshifts e 8 stops |
| 4 | lesão compartilhada por linhagens que divergiram depois | 4 stops do ENAM em Palaeognathae |

## 2. Material e métodos

O texto completo, com versões, parâmetros, datas de acesso e as tentativas que mudaram o rumo, está em `docs/metodos.md` no repositório do projeto. Resumimos aqui o que é necessário para ler os resultados.

**Espécies.** Usamos 18 montagens RefSeq do NCBI (Tabela 2), acessadas entre 12 e 17 de outubro de 2026: três crocodilianos e três escamados (com dente), quatro tartarugas e oito aves (sem dente). As aves incluem quatro paleognatas (*Struthio*, *Apteryx*, *Dromaius*, *Tinamus*) e quatro neognatas (*Gallus*, *Anas*, *Strigops*, *Taeniopygia*).

**Tabela 2.** Montagens analisadas.

| Grupo | Espécie | Montagem |
|---|---|---|
| Crocodylia | *Alligator mississippiensis* | GCF_030867095.1 |
| | *Gavialis gangeticus* | GCF_001723915.1 |
| | *Crocodylus porosus* | GCF_001723895.1 |
| Squamata | *Anolis carolinensis* | GCF_035594765.1 |
| | *Python bivittatus* | GCF_059951505.1 |
| | *Podarcis muralis* | GCF_964188315.1 |
| Testudines | *Chrysemys picta* | GCF_011386835.1 |
| | *Chelonia mydas* | GCF_015237465.2 |
| | *Gopherus evgoodei* | GCF_007399415.2 |
| | *Pelodiscus sinensis* | GCF_049634645.1 |
| Aves | *Gallus gallus* | GCF_016699485.2 (GRCg7b) |
| | *Anas platyrhynchos* | GCF_047663525.1 |
| | *Strigops habroptila* | GCF_004027225.2 |
| | *Taeniopygia guttata* | GCF_048771995.1 |
| | *Dromaius novaehollandiae* | GCF_036370855.1 |
| | *Struthio camelus* | GCF_040807025.1 |
| | *Apteryx rowi* | GCF_003343035.1 |
| | *Tinamus guttatus* | GCF_059628265.1 |

**Intervalos sintênicos.** No jacaré, o cluster fica no chr3 na ordem NRG1 → ENAM → AMBN → AMTN → LPL → PSD3. Em cada espécie extraímos o intervalo entre as âncoras NRG1 e LPL, com 60 kb de margem. Para AMELX, MMP20, DSPP e ACP4, que estão em outros cromossomos, cada gene recebeu o próprio intervalo, delimitado pelo vizinho nomeado mais próximo de cada lado.

**Busca e índices.** As proteínas do jacaré foram buscadas nos intervalos com tBLASTn (BLAST+ 2.17.0; Camacho et al. 2009; e-value 10, filtro SEG). Decidimos pelo bitscore, que não depende do tamanho do banco. Um trecho só contou se (i) nenhuma outra proteína do intervalo o alinhou melhor (competição local) e (ii) o melhor hit dele no proteoma inteiro do jacaré (48.040 proteínas) foi a própria consulta (recíproco). O índice de retenção foi a soma dos bitscores de ENAM, AMBN e AMTN dividida pelo bitscore do LPL na mesma espécie. LPL é um gene funcional no mesmo intervalo e serve de controle positivo.

**Nulo por janela, segunda consulta e composição.** Para cada proteína geramos 20 versões embaralhadas com a mesma composição (semente 20261017) e as buscamos no mesmo intervalo. No índice só entrou o excesso do bitscore sobre o percentil 95 desse nulo. Repetimos tudo com as proteínas do *Anolis*, que está à mesma distância de crocodilianos, tartarugas e aves (cerca de 281 Ma; Kumar et al. 2022). Como controle de composição usamos a IBSP do *Anolis*, uma proteína óssea desordenada como ENAM e AMBN que é funcional nas aves.

**Lesões (nível 3).** As proteínas do jacaré foram alinhadas aos intervalos com miniprot 0.18-r281 (Li 2023), que modela íntrons e marca frameshifts e stops em quadro. A taxa de lesão do alinhador foi calibrada nas seis espécies com dente, incluindo três escamados cuja identidade ao jacaré no ENAM (0,36–0,39) está na mesma faixa das relíquias.

**Lesões compartilhadas (nível 4).** A relíquia de ENAM do kakapo foi usada como consulta em nucleotídeo (dc-megablast) nas outras aves. Os trechos homólogos e o CDS de três crocodilianos foram alinhados com MAFFT L-INS-i (Katoh e Standley 2013) e com MUSCLE 5 (Edgar 2022). Uma lesão contou quando teve o mesmo tipo, posição e tamanho em duas ou mais aves, apareceu nos dois alinhadores e não estava em homopolímero (4 ou mais bases iguais). A lesão foi datada no menor clado que contém as aves que a têm, desde que toda ave desse clado que cobre a posição também a tivesse.

**ODAPH, ODAM, recombinação e alinhamento de genoma inteiro.** Identificamos o ODAPH do jacaré por posição e estrutura de éxons (Liang et al. 2021) e o buscamos nos intervalos e nos genomas inteiros de frango, kakapo e jacaré. Medimos a sensibilidade dessa busca contra ODAPH funcionais de lagartos. A taxa de recombinação na janela foi lida do mapa de ligação do frango de Elferink et al. (2010), em bins de cerca de 500 kb (taxa de macho). A validação independente usou o alinhamento multiz100way do hg38 (Blanchette et al. 2004), recortado nos 29 éxons codificantes humanos de ENAM, AMBN e AMTN.

**Estatística.** Não calculamos teste de hipótese. As espécies compartilham filogenia e o número de origens independentes de perda é dois, então todas as comparações são descritivas.

## 3. Resultados

### 3.1 Relíquias de ENAM e AMBN no cromossomo Z do frango

No chrZ do frango (GRCg7b), o intervalo NRG1–LPL tem 105 kb, contra 893 kb no jacaré (Figura 1). Antes de olhar os alvos, conferimos os controles. As proteínas do jacaré se encontraram no próprio genoma com 99,9–100% de identidade. O LPL do jacaré, controle positivo, foi encontrado no intervalo do frango com e = 1,16 × 10^−57^ (64,9% de identidade em 154 aa). O DSPP, controle negativo de outro cromossomo, deu 22,7 bits (e = 8,6), o que define o piso de ruído.

Contra esse piso, AMBN rendeu e = 1,40 × 10^−9^ (54,8% em 31 aa) e ENAM e = 1,08 × 10^−7^ (41,2% em 68 aa). Quatro dos oito HSPs de ENAM se agrupam em cerca de 1,8 kb. As duas regiões passaram nos testes de cópia única e de melhor hit recíproco (ENAM 42,7 bits, com o concorrente seguinte a e = 1,7; AMBN 44,7 bits, com o seguinte a e = 4,2). AMTN não saiu do ruído (e = 1,6). A única região com similaridade nucleotídica forte ao span do AMTN (72,8% em 158 pb) não codifica proteína do cluster quando traduzida e foi classificada como elemento não codificante conservado ou repetição.

As relíquias ficam em chrZ:55.155.480–55.157.241 (ENAM) e 55.168.855–55.168.992 (AMBN), nos vãos entre lncRNAs anotados pelo Ensembl, e não dentro deles. A ordem ENAM → AMBN é a mesma do jacaré. Nenhuma das duas regiões tem gene, pseudogene ou lncRNA anotado no NCBI ou no Ensembl.

![**Figura 1.** Sintenia do cluster ENAM–AMBN–AMTN entre o chr3 do jacaré (rAllMis1) e o chrZ do frango (GRCg7b). As âncoras NRG1, LPL e PSD3 estão ligadas por linhas tracejadas. No frango, as relíquias de ENAM e AMBN (ΨENAM, ΨAMBN) ficam no vão NRG1–LPL, na ordem do jacaré; AMTN não foi detectado. Cada traço está em escala própria, proporcional às coordenadas da montagem.](runs/relatorios/paper_fig/fig1_sintenia.png)

### 3.2 Anotação e sintenia nas demais espécies

Nas sete espécies do primeiro painel, o número de genes do esmalte e da dentina anotados (de nove) foi 7 no jacaré, 6 no *Anolis*, 5 no gavial, 1 na tartaruga-pintada (só ENAM, no chr6) e 0 no frango, no periquito e no kakapo. O intervalo NRG1–LPL encolhe de 893 kb no jacaré para 337 kb na tartaruga-pintada, 293 kb no kakapo e 105 kb no frango. No kakapo, como no frango, o intervalo fica no cromossomo Z. Como o Z é ancestral das aves, isso é uma observação única e não duas.

### 3.3 O índice de retenção por similaridade é confundido

Com seis espécies, o índice bruto separou os grupos (com dente 2,51–13,02; sem dente 0,61–1,84). A separação não resistiu à ampliação e aos controles (Figura 2, Tabela 3).

Com o jacaré como consulta e o nulo por janela (percentil 95 entre 26 e 32 bits em todas as janelas), os grupos se sobrepõem: a tartaruga-verde (2,20) supera a píton (1,57) e a lagartixa-dos-muros (1,45). A consulta do jacaré favorece os parentes próximos dele, e o *Anolis*, escamado com dente, fica a cerca de 281 Ma do jacaré, mais longe que tartarugas (257 Ma) e aves (243 Ma) (Kumar et al. 2022).

Com o *Anolis* como consulta, crocodilianos, tartarugas e aves ficam à mesma distância da consulta, e a separação reaparece: crocodilianos 2,29–2,84, tartarugas 0–0,57, aves 0–0,10. O LPL funcional, pela mesma consulta, rende mais nas aves (162–172 bits) que nos crocodilianos (149–151), o que exclui uma taxa de evolução maior nas aves como explicação geral.

A composição pesa. A IBSP, desordenada como ENAM e AMBN mas funcional nas aves, teve excesso sobre o nulo de 52–57 bits nos crocodilianos e 39–41 nas aves, ou seja, as aves ficaram com cerca de 73% do nível dos crocodilianos. O ENAM, pela mesma consulta, teve 229–254 bits nos crocodilianos e 0–1,7 nas aves (cerca de 0%). O tBLASTn perde cerca de um quarto da sensibilidade nas aves para essa classe de proteína, mas essa perda não explica a queda do ENAM.

**Tabela 3.** Índice de retenção com nulo por janela (faixa por grupo; a própria espécie da consulta fica fora da faixa).

| Consulta | Crocodylia | Squamata | Testudines | Aves |
|---|---|---|---|---|
| *Alligator* | 12,34–14,11 | 1,45–2,45 | 0–2,20 | 0–0,71 |
| *Anolis* | 2,29–2,84 | 4,16–6,59 | 0–0,57 | 0–0,10 |

![**Figura 2.** Índice de retenção (excesso de bitscore sobre o percentil 95 de 20 consultas embaralhadas, ENAM + AMBN + AMTN, dividido pelo LPL) nas 18 espécies, com a consulta do jacaré (esquerda) e do *Anolis* (direita). Círculos cheios: com dente; vazados: sem dente. A linha da própria espécie da consulta está marcada. Com a consulta do jacaré os grupos se sobrepõem; com a do *Anolis*, equidistante de crocodilianos, tartarugas e aves, eles se separam.](runs/relatorios/paper_fig/fig2_indice.png)

### 3.4 Lesões em quadro (nível 3)

Nas seis espécies com dente, o miniprot alinhou 12.160 aminoácidos dos genes do painel com uma única lesão, um frameshift no AMBN do *Anolis*: 0,08 lesão por mil aminoácidos (Figura 3). Os três escamados alinharam o ENAM com identidade de 0,36 a 0,39 e zero lesão. Esse valor está na mesma faixa do kakapo (0,35), o que afasta a explicação de que lesões surgem apenas porque o alinhamento ficou ruim.

Nas espécies sem dente, os alinhamentos com lesão somaram 53 lesões em 3.518 aminoácidos (15,1 por mil). Somando também os alinhamentos sem lesão, foram 53 lesões em 4.350 aminoácidos (12,2 por mil) (Tabela 4). O ENAM tem lesões no kakapo e em duas tartarugas, e o AMBN em duas tartarugas e na ema. No frango a relíquia está degradada demais para o miniprot: o ENAM não alinha e o AMBN alinha só 10 aminoácidos. O AMTN não alinhou em nenhuma espécie sem dente, exceto um trecho de 20 aminoácidos no *Gopherus*, sem lesão.

**Tabela 4.** Alinhamentos do miniprot com lesão nas espécies sem dente.

| Gene | Espécie | Cobertura | Identidade | Frameshifts | Stops |
|---|---|---|---|---|---|
| ENAM | *Strigops habroptila* | 95% | 0,35 | 6 | 8 |
| ENAM | *Chrysemys picta* | 83% | 0,50 | 11 | 4 |
| ENAM | *Chelonia mydas* | 66% | 0,55 | 8 | 3 |
| AMBN | *Chrysemys picta* | 100% | 0,47 | 3 | 1 |
| AMBN | *Chelonia mydas* | 100% | 0,59 | 2 | 1 |
| AMBN | *Dromaius novaehollandiae* | 74% | 0,51 | 1 | 3 |
| DSPP | *Chrysemys picta* | 48% | 0,42 | 2 | 0 |

O DSPP entra na tabela, mas não é interpretado: o controle positivo dele falhou no gavial e no *Anolis* (seção 3.6).

![**Figura 3.** Lesões (frameshift + stop) por mil aminoácidos alinhados pelo miniprot. Abaixo da linha, todos os alinhamentos das espécies com dente (1 lesão em 12.160 aa); acima, os alinhamentos das espécies sem dente que têm lesão. O rótulo de cada barra dá lesões / aminoácidos alinhados.](runs/relatorios/paper_fig/fig3_lesoes.png)

### 3.5 Lesões compartilhadas no ENAM (nível 4)

A relíquia de ENAM do kakapo, usada como consulta em nucleotídeo, encontrou o trecho homólogo no pato, no mandarim e nas três paleognatas (66–73% de identidade). No frango e no tinamu não encontrou nada. No alinhamento de 2,4 kb com três crocodilianos, os crocodilianos tiveram zero lesão nos dois alinhadores. Das 12 lesões compartilhadas por duas ou mais aves e robustas aos dois alinhadores, 9 foram congruentes com a árvore (Tabela 5, Figura 4):

- **Palaeognathae** (*Struthio*, *Apteryx*, *Dromaius*): 4 stops, nos códons 158, 389, 472 e 539 do jacaré.
- **Novaeratitae** (*Apteryx* e *Dromaius*): 3 stops (códons 189, 347 e 455) e uma deleção de 10 nt (códon 200).
- **Psittacopasserae** (*Strigops* e *Taeniopygia*): 1 stop (códon 791).
- **Neornithes**: nenhuma lesão congruente.

As outras três, todas stops, juntam paleognatas e neognatas sem incluir todas as aves que cobrem a posição, e por isso são incongruentes com a árvore. Duas estão em códons TGG, em que uma única transição G→A gera stop. Em duas delas as aves usam códons de parada diferentes (TAG e TGA; TAA e TGA), o que indica origens independentes. A homoplasia medida (3 de 11 stops compartilhados) mostra que um stop compartilhado vale menos como evidência que um indel compartilhado.

**Tabela 5.** Lesões do ENAM compartilhadas por duas ou mais aves (posição no CDS do jacaré).

| Códon | Tipo | Aves com a lesão | Leitura |
|---|---|---|---|
| 158 | stop (CAA→TGA) | *Apteryx*, *Dromaius*, *Struthio* | Palaeognathae |
| 189 | stop (TCC→TGA) | *Apteryx*, *Dromaius* | Novaeratitae |
| 200 | deleção de 10 nt | *Apteryx*, *Dromaius* | Novaeratitae |
| 227 | stop (GCA→TAG/TGA) | *Anas*, *Apteryx*, *Dromaius* | incongruente |
| 347 | stop (AGC→TGA) | *Apteryx*, *Dromaius* | Novaeratitae |
| 350 | stop (TGG→TGA) | *Struthio*, *Strigops* | incongruente |
| 389 | stop (TTT→TAA) | *Apteryx*, *Dromaius*, *Struthio* | Palaeognathae |
| 403 | stop (TGG→TAA/TGA) | *Dromaius*, *Struthio*, *Strigops* | incongruente |
| 455 | stop (CAG→TAG) | *Apteryx*, *Dromaius* | Novaeratitae |
| 472 | stop (CAA→TAG) | *Apteryx*, *Dromaius*, *Struthio* | Palaeognathae |
| 539 | stop (TGG→TGA) | *Apteryx*, *Dromaius*, *Struthio* | Palaeognathae |
| 791 | stop (TGG→TGA) | *Strigops*, *Taeniopygia* | Psittacopasserae |

![**Figura 4.** Lesões do ENAM compartilhadas e congruentes, mapeadas no cladograma das seis aves com trecho homólogo. Círculos marcam o menor clado que contém todas as aves com a lesão. Nenhuma lesão congruente foi encontrada no nó de Neornithes; as três candidatas desse nível são homoplasias. Os crocodilianos (gene funcional) não têm lesão.](runs/relatorios/paper_fig/fig4_arvore.png)

### 3.6 O resto do painel e a calibração com pseudogene conhecido

Para MMP20, a fração recuperada do sinal do jacaré (soma de bits dos trechos recíprocos) foi 89% no gavial, 57% no *Anolis*, 30% na tartaruga-pintada, 9% no kakapo e 6% no frango. As âncoras do intervalo (DCUN1D5 e TMEM123) renderam 75–124 bits nas aves e na tartaruga, então o intervalo estava correto. O MMP20 da tartaruga-pintada, porém, alinhou 91% pelo miniprot, com identidade de 0,55 e sem lesão. Os 30% medem perda de identidade e não demonstram pseudogene.

AMELX, DSPP e ACP4 ficaram sem conclusão. No AMELX, nem o *Anolis* saiu do ruído. No DSPP, o controle falhou no gavial e no *Anolis*. No ACP4, a sintenia não se estabeleceu nas aves. O pseudogene de AMEL da galinha (Ψ-AMEL; Sire et al. 2008) serviu de calibração. Ele está dentro do intervalo MSL3–MID1, com 100% de identidade, em NC_052532.1:124.410.142–124.413.840 (fita −, GRCg7b), a primeira coordenada dele numa montagem atual. O tBLASTn do AMELX do jacaré não deixou HSP nenhum ali e, sem o filtro SEG, rendeu 25,8 bits, indistinguíveis do ruído. Ou seja, a proteína de crocodiliano não detecta um pseudogene de ave conhecido. Um zero do tBLASTn não vale como ausência sem controle positivo na mesma espécie.

### 3.7 ODAPH e ODAM

No jacaré, o ODAPH é o LOC106740521 (XP_019349100.2, anotado como "uncharacterized protein"), entre CDKL2 e THAP6 (NC_081825.1:164.672.404–164.676.187), com três éxons codificantes e o motivo SSEE descrito por Liang et al. (2021). Ele está intacto em *Crocodylus*, *Anolis* e *Podarcis*. Nas duas últimas espécies o locus está anotado como lncRNA, um artefato de anotação. Nas aves e nas tartarugas, o ODAPH fica no nível 2: as mesmas consultas o encontram nas espécies com dente (e ≤ 3 × 10^−6^) e as âncoras são encontradas em todas (e ≤ 10^−23^). Em sete aves que não o frango, CDKL2 encosta em RCHY1 (vão de 0,15 a 2,4 kb): o espaço onde ODAPH e THAP6 ficam no jacaré não existe.

A busca no genoma inteiro não tem poder para decidir se o gene foi deslocado ou perdido. Os melhores hits no frango e no kakapo (32–33 bits, recíprocos ao ODAPH do jacaré) não são o mesmo locus nas duas aves, pois não têm nenhum vizinho em comum em ±400 kb. E o ODAPH do jacaré encontra o ODAPH funcional do *Anolis* com apenas 36,6 bits e o do *Podarcis* com 26,9, dentro da faixa de ruído. A evidência sobre o ODAPH nas aves é posicional, e não de sequência.

O ODAM não foi encontrado em nenhuma das 18 espécies, nem nas com dente. Sem referência positiva em saurópsido, a escada não começa, e o ODAM saiu do painel.

### 3.8 Recombinação na janela da relíquia

No mapa de Elferink et al. (2010), a janela NRG1–LPL do frango (galGal3, chrZ:52.836.615–53.408.327) tem taxa ponderada de 2,27 cM/Mb em machos. A média do Z é 3,00 e a autossômica 2,97, e a janela fica no percentil 54 entre os 122 bins do Z. Ela não é um deserto de recombinação. A taxa varia bastante ao longo de cada cromossomo do frango e se correlaciona com o conteúdo de GC (Groenen et al. 2009), e por isso a resolução importa. A relíquia, porém, cai perto da fronteira entre um bin de 0,68 cM/Mb (em que as duas populações do mapa divergem: 2,72 e 0,85) e um de 3,94 cM/Mb. Nessa resolução, uma taxa localmente baixa não fica excluída.

### 3.9 Validação por alinhamento de genoma inteiro

No multiz100way do hg38 (Tabela 6), o jacaré tem 86% do CDS humano de ENAM alinhado e o *Anolis* 74%. Nas tartarugas são 56% (*Chelonia*) e 63% (*Chrysemys*); o *Pelodiscus* não tem nenhuma base alinhada. Das 14 aves do alinhamento, 9 (incluindo frango, peru, pato, mandarim e periquito) não têm nenhuma base alinhada em nenhum dos três genes. As outras cinco têm fragmentos: papagaio-de-porto-rico (*Amazona vittata*), ENAM 41%; *Ficedula* e *Zonotrichia*, ENAM 16–17%; os dois falcões, AMBN 64%. É o mesmo padrão por um método independente (cadeias lastz e multiz, com o humano como referência). A contagem de frameshifts pelo multiz foi descartada porque reprovou no controle: o ENAM funcional do jacaré teve 54 contra o humano.

**Tabela 6.** Fração do CDS humano com base alinhada no multiz100way (hg38).

| Montagem | Grupo | ENAM | AMBN | AMTN |
|---|---|---|---|---|
| allMis1 (*Alligator*) | com dente | 86% | 85% | 80% |
| anoCar2 (*Anolis*) | com dente | 74% | 38% | 0% |
| cheMyd1 (*Chelonia*) | sem dente | 56% | 84% | 7% |
| chrPic2 (*Chrysemys*) | sem dente | 63% | 0% | 0% |
| pelSin1 (*Pelodiscus*) | sem dente | 0% | 0% | 0% |
| amaVit1 (*Amazona*) | ave | 41% | 0% | 0% |
| falChe1, falPer1 (falcões) | ave | 0% | 64% | 0% |
| ficAlb2, zonAlb1 | ave | 16–17% | 0% | 0% |
| 9 outras aves | ave | 0% | 0% | 0% |

## 4. Discussão

### 4.1 O que é replicação e o que é novo

Os resultados centrais replicam a literatura. Meredith et al. (2013) já tinham recuperado cópias inativadas de ENAM e AMBN em aves, com frameshifts compartilhados, e restos dos dois genes na tartaruga-pintada com a sintenia original. A conclusão de Sire et al. (2008), de que ENAM e AMBN desapareceram do frango, não se sustenta no GRCg7b: os restos codificantes estão no vão sintênico do chrZ, na ordem ancestral. O que este trabalho acrescenta é mais modesto e mais verificável: coordenadas numa montagem atual (incluindo a do Ψ-AMEL), quatro tartarugas e quatro paleognatas no mesmo procedimento, e um conjunto de controles que dizem quanto cada afirmação vale.

### 4.2 Índices de similaridade como fenótipo

O resultado metodológico mais útil é negativo. Um índice de retenção calculado por bitscore, do tipo que se usaria como fenótipo contínuo num teste de associação entre taxas de evolução e caráter (Kowalczyk et al. 2019), mede ao mesmo tempo a perda, a distância entre consulta e alvo e a cegueira do alinhador para proteínas desordenadas. Com seis espécies e a consulta do jacaré, o índice separou os grupos; com 18 espécies e um nulo por janela, os grupos se sobrepuseram. O contraste realmente informativo exigiu uma consulta equidistante das linhagens comparadas, um controle de taxa (LPL) e um controle de composição (IBSP). Recomendamos esses três controles, mais o nulo por janela, antes de usar qualquer índice de similaridade como variável em métodos de associação filogenética (Prudent et al. 2016; Saputra et al. 2021; Redlich et al. 2024). Mesmo com eles, duas origens independentes de perda são pouco para um teste.

### 4.3 Por que quase não há lesão compartilhada por todas as aves

No trecho de 2,4 kb do ENAM, as lesões robustas datam o decaimento até o ancestral das paleognatas e o de Psittacopasserae, e não até o de Neornithes. Isso não contradiz a perda única dos dentes no ancestral das aves. A função pode ter caído primeiro, por outro gene ou por outra região do ENAM, e as lesões detectáveis acumulam-se depois, ao longo de cada ramo. A janela para lesões compartilhadas por todas as aves é curta. A perda do esmalte foi estimada em cerca de 116 Ma (Meredith et al. 2014), e a divergência entre paleognatas e neognatas em cerca de 108 Ma (IC 102–114; Kumar et al. 2022). As duas estimativas vêm de fontes diferentes e os intervalos se sobrepõem, mas a ordem de grandeza indica um ramo de poucos milhões de anos.

Meredith et al. (2013) relatam frameshifts compartilhados anteriores à divisão basal de Neognathae, que não recuperamos. Três fatores explicam a diferença sem que um dos resultados invalide o outro. Primeiro, Galloanserae está representado aqui só pelo pato, porque o frango não teve trecho homólogo em nucleotídeo. Segundo, nossos critérios são estritos: exigimos os dois alinhadores, excluímos homopolímeros e exigimos congruência com a árvore. Terceiro, o trecho analisado é outro. "Não recuperado" não é "ausente".

### 4.4 ODAPH: posição sem sequência

Em mamíferos sem dente, o ODAPH (C4orf26) está inativado (Springer et al. 2016), o que torna esperado o decaimento dele nas aves. Aqui a sequência não decide: a régua mostra que a proteína do jacaré mal encontra o ODAPH funcional de um lagarto. O que sustenta a ausência nas aves é que o espaço do gene entre CDKL2 e THAP6 não existe em sete aves. Isso é compatível tanto com perda quanto com deslocamento antes de Neornithes, e a distinção exige um método que não dependa de similaridade de proteína a 250–280 Ma, como alinhamento de genoma inteiro com projeção de genes.

### 4.5 Limitações

- O número de espécies é pequeno e elas não são independentes. Nenhuma comparação entre grupos foi testada estatisticamente.
- A relíquia do frango está degradada demais para o miniprot, e o nível 4 não inclui Galloanserae além do pato.
- As lesões podem refletir erros de montagem. Lesões compartilhadas entre montagens independentes reduzem esse risco, mas o nível 3 numa única espécie não o elimina.
- O alinhamento multiz100way usa montagens antigas (galGal4, por exemplo) e o humano como referência, a cerca de 319 Ma do jacaré (Kumar et al. 2022).
- A taxa de recombinação vem de bins de 500 kb em galGal3, sem conversão de coordenadas da relíquia.
- O AMTN não foi detectado em nenhuma ave. Isso é nível 2 na melhor das hipóteses, com o ponto cego medido no Ψ-AMEL.

### 4.6 Próximos passos

Três passos aumentam a precisão sem exigir dados novos: projeção de genes por alinhamento de genoma inteiro com o jacaré como referência (TOGA), para o ODAPH e o AMTN; o nível 4 para o AMBN e para as tartarugas; e um teste de seleção relaxada (RELAX, HyPhy) nos ramos sem dente, onde houver alinhamento de códons suficiente. Montagens de leitura longa de Galloanserae permitiriam incluir o frango no nível 4.

## 5. Conclusão

As relíquias de ENAM e AMBN persistem no cromossomo Z do frango, no vão sintênico entre NRG1 e LPL, e com lesões em quadro no kakapo, na ema e em duas tartarugas. As lesões robustas datam o decaimento até o ancestral das paleognatas e o de Psittacopasserae. O índice de retenção por similaridade, sem nulo, sem consulta equidistante e sem controle de composição, não distingue perda de cegueira do método. A afirmação de perda só vale com o nível de evidência e o controle que a sustentam.

## Disponibilidade de dados e código

Todos os genomas são públicos (NCBI RefSeq e UCSC). Está no repositório do projeto (github.com/patoangolano/bio-stack):

- scripts de análise: `scripts/atavismo/s24`–`s36`;
- saídas tabulares: `runs/atavismo/`;
- registro de proveniência de cada operação: `runs/provenance.jsonl`;
- figuras e PDF: `scripts/paper_figuras.py` e `scripts/paper_pdf.py`;
- métodos detalhados: `docs/metodos.md`.

## Declaração de uso de inteligência artificial

O desenho das análises, a escrita dos scripts, a execução e a redação deste texto foram feitos com assistência do Claude Code (Anthropic, modelo Claude Opus 5.5), sob direção e revisão do autor. Uma crítica externa do método foi produzida por um agente de IA e conferida contra os dados. Nenhum dado biológico foi gerado por IA, e todas as referências foram conferidas no PubMed.

## Referências

- Blanchette M, Kent WJ, Riemer C, et al. 2004. Aligning multiple genomic sequences with the threaded blockset aligner. *Genome Research* 14(4):708–715. doi:10.1101/gr.1933104
- Camacho C, Coulouris G, Avagyan V, et al. 2009. BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. doi:10.1186/1471-2105-10-421
- Edgar RC. 2022. Muscle5: high-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. *Nature Communications* 13(1):6968. doi:10.1038/s41467-022-34630-w
- Elferink MG, van As P, Veenendaal T, Crooijmans RPMA, Groenen MAM. 2010. Regional differences in recombination hotspots between two chicken populations. *BMC Genetics* 11:11. doi:10.1186/1471-2156-11-11
- Groenen MAM, Wahlberg P, Foglio M, et al. 2009. A high-density SNP-based linkage map of the chicken genome reveals sequence features correlated with recombination rate. *Genome Research* 19(3):510–519. doi:10.1101/gr.086538.108
- Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30(4):772–780. doi:10.1093/molbev/mst010
- Kowalczyk A, Meyer WK, Partha R, Mao W, Clark NL, Chikina M. 2019. RERconverge: an R package for associating evolutionary rates with convergent traits. *Bioinformatics* 35(22):4815–4817. doi:10.1093/bioinformatics/btz468
- Kumar S, Suleski M, Craig JM, et al. 2022. TimeTree 5: an expanded resource for species divergence times. *Molecular Biology and Evolution* 39(8):msac174. doi:10.1093/molbev/msac174
- Li H. 2023. Protein-to-genome alignment with miniprot. *Bioinformatics* 39(1):btad014. doi:10.1093/bioinformatics/btad014
- Liang T, Hu Y, Kawasaki K, et al. 2021. Odontogenesis-associated phosphoprotein truncation blocks ameloblast transition into maturation in *Odaph*^C41\*/C41\*^ mice. *Scientific Reports* 11(1):1132. doi:10.1038/s41598-020-80912-y
- Meredith RW, Gatesy J, Springer MS. 2013. Molecular decay of enamel matrix protein genes in turtles and other edentulous amniotes. *BMC Evolutionary Biology* 13:20. doi:10.1186/1471-2148-13-20
- Meredith RW, Zhang G, Gilbert MTP, Jarvis ED, Springer MS. 2014. Evidence for a single loss of mineralized teeth in the common avian ancestor. *Science* 346(6215):1254390. doi:10.1126/science.1254390
- Prudent X, Parra G, Schwede P, Roscito JG, Hiller M. 2016. Controlling for phylogenetic relatedness and evolutionary rates improves the discovery of associations between species' phenotypic and genomic differences. *Molecular Biology and Evolution* 33(8):2135–2150. doi:10.1093/molbev/msw098
- Redlich R, Kowalczyk A, Tene M, et al. 2024. RERconverge expansion: using relative evolutionary rates to study complex categorical trait evolution. *Molecular Biology and Evolution* 41(11):msae210. doi:10.1093/molbev/msae210
- Saputra E, Kowalczyk A, Cusick L, Clark N, Chikina M. 2021. Phylogenetic permulations: a statistically rigorous approach to measure confidence in associations in a phylogenetic context. *Molecular Biology and Evolution* 38(7):3004–3021. doi:10.1093/molbev/msab068
- Sire JY, Delgado SC, Girondot M. 2008. Hen's teeth with enamel cap: from dream to impossibility. *BMC Evolutionary Biology* 8:246. doi:10.1186/1471-2148-8-246
- Springer MS, Starrett J, Morin PA, et al. 2016. Inactivation of C4orf26 in toothless placental mammals. *Molecular Phylogenetics and Evolution* 95:34–45. doi:10.1016/j.ympev.2015.11.002
- Springer MS, Emerling CA, Gatesy J, et al. 2019. Odontogenic ameloblast-associated (ODAM) is inactivated in toothless/enamelless placental mammals and toothed whales. *BMC Evolutionary Biology* 19(1):31. doi:10.1186/s12862-019-1359-6
