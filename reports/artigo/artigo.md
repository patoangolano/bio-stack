---
titulo: POSIÇÃO FILOGENÉTICA DAS AVES DE RAPINA EM TELLURAVES
subtitulo: síntese da evidência genômica, morfológica e fóssil e um conjunto filogenômico novo
titulo_en: Phylogenetic position of birds of prey within Telluraves
subtitulo_en: a synthesis of genomic, morphological and fossil evidence and a new phylogenomic dataset
autor: Matheus Lima
nota_autor: [Formação acadêmica e vínculo institucional a preencher]. E-mail: limatheus458@gmail.com.
---

#* RESUMO

As linhagens tradicionais de aves de rapina (Cathartidae, Accipitriformes, Falconiformes e Strigiformes) pertencem a Telluraves, mas a ordem dos ramos profundos desse clado segue em disputa. Este trabalho combinou quatro abordagens computacionais: a reanálise de 128.233 árvores de genes do projeto B10K e de 28 árvores-fonte da Open Tree of Life; o posicionamento de fósseis em dez matrizes morfológicas sob arcabouço molecular; a reconstrução do estado ancestral do hábito rapinante em 202 famílias sob seis topologias; e um conjunto filogenômico novo, com 57 táxons e 300 locos codificantes. Australaves, com Falconiformes como grupo-irmão de Psittacopasserae, e Accipitrimorphae, na ordem (Cathartidae, (Sagittariidae, (Pandionidae, Accipitridae))), resistiram à troca de tipo de dado, de estimador e de amostragem. A base de Afroaves não resistiu: Hieraves e Strigiformes + Coraciimorphae receberam 35,0 % e 34,6 % dos quartetos do B10K, e a resolução vencedora mudou com o estimador. No conjunto novo, a supermatriz recuperou Australaves (UFBoot 96), e o teste AU rejeitou Cariamiformes + Falconiformes (p = 0,019), mas não Hieraves (p = 0,13); as árvores de genes de éxons agruparam falcões e seriema. A morfologia uniu as rapinantes por caracteres convergentes de pé e garra, e os pés rapinantes surgiram três vezes em todas as topologias. Um ancestral predador de Telluraves só apareceu quando a dieta dos fósseis, inferida da própria anatomia rapinante, entrou na análise. A base de Afroaves deve ser tratada como politomia não resolvida.

Palavras-chave: filogenômica; Telluraves; aves de rapina; fatores de concordância; convergência morfológica.

#* ABSTRACT

The traditional lineages of birds of prey (Cathartidae, Accipitriformes, Falconiformes and Strigiformes) belong to Telluraves, but the order of the deep branches of this clade remains contested. This study combined four computational approaches: a reanalysis of 128,233 gene trees from the B10K project and of 28 source trees from the Open Tree of Life; the placement of fossils in ten morphological matrices under a molecular backbone; ancestral state reconstruction of the raptorial habit in 202 families under six topologies; and a new phylogenomic dataset of 57 taxa and 300 coding loci. Australaves, with Falconiformes sister to Psittacopasserae, and Accipitrimorphae, in the order (Cathartidae, (Sagittariidae, (Pandionidae, Accipitridae))), withstood changes of data type, estimator and sampling. The base of Afroaves did not: Hieraves and Strigiformes + Coraciimorphae received 35.0% and 34.6% of the B10K quartets, and the winning resolution changed with the estimator. In the new dataset, the supermatrix recovered Australaves (UFBoot 96), and the AU test rejected Cariamiformes + Falconiformes (p = 0.019) but not Hieraves (p = 0.13); exon gene trees grouped falcons with the seriema. Morphology united raptors through convergent foot and claw characters, and raptorial feet arose three times under every topology. A predatory ancestor of Telluraves was recovered only when fossil diets, inferred from the raptorial anatomy itself, entered the analysis. The base of Afroaves should be treated as an unresolved polytomy.

Keywords: phylogenomics; Telluraves; birds of prey; concordance factors; morphological convergence.

# 1 INTRODUÇÃO

Durante mais de um século, gaviões, falcões, urubus e o secretário formaram a ordem Falconiformes, e as corujas, a ordem Strigiformes, vizinha nas classificações. As análises cladísticas reforçaram essa vizinhança: com 148 caracteres, Strigiformes saiu como grupo-irmão de Falconidae + Accipitridae [@mayrclarke2003], e com 2.954 caracteres, como grupo-irmão de Falconiformes, com Pandionidae junto de Falconidae [@livezey2007]. Os dados moleculares desfizeram esse arranjo. Cinco genes nucleares já separavam o secretário e os acipitrídeos dos falcões [@ericson2006], e 19 locos consolidaram os falcões junto de papagaios e passeriformes [@hackett2008]. Inserções de retroposons, marcadores quase livres de homoplasia, confirmaram a posição dos falcões [@suh2011], e os genomas completos de 48 espécies [@jarvis2014], os 259 locos de captura de 198 espécies [@prum2015], as regiões 3′ não traduzidas de 429 espécies [@kuhl2021] e os 363 genomas do B10K [@stiller2024] dispersaram as rapinantes por Telluraves, as aves terrestres centrais de Neoaves.

Os nomes de clado usados aqui seguem a literatura de origem. Psittacopasserae (papagaios + passeriformes) e Eufalconimorphae (falcões + Psittacopasserae) foram propostos por @suh2011; Telluraves e Australaves (Cariamiformes + Eufalconimorphae), por @yuri2013; Afroaves vem do cenário biogeográfico de @ericson2012; Hieraves (Accipitriformes + Strigiformes), de @wu2024. Accipitrimorphae reúne Cathartidae e Accipitriformes no sentido estrito (Sagittariidae, Pandionidae e Accipitridae), e Coraciimorphae, Coliiformes e Cavitaves (Leptosomiformes, Trogoniformes, Bucerotiformes, Coraciiformes e Piciformes).

Parte do problema está resolvida. Telluraves é recuperado por quase toda análise molecular de grande escala desde 2008 e tem 64 inserções de retroposons sem conflito [@suh2015]. O grupo-irmão de Telluraves, ao contrário, mudou a cada grande estudo: Aequornithes + Phaethontimorphae [@jarvis2014; @kuhl2021], o hoatzin [@prum2015], Elementaves [@stiller2024] ou uma dicotomia com Aquaterraves [@wu2024], esta última sem sustentação numa reanálise [@springer2024]. Outras análises colocam Telluraves como grupo-irmão de todo o restante de Neoaves [@braun2021] ou reduzem a base de Neoaves a uma politomia [@houde2019; @suh2016]. A causa é cronológica. As divergências profundas de Telluraves se comprimem em pouco mais de 1 milhão de anos logo após o limite Cretáceo–Paleógeno, e os retroposons estimam de 59 % a 81 % de separação incompleta de linhagens (ILS) nos três ramos mais profundos do clado [@suh2015]. Um coliiforme-tronco do Paleoceno inicial indica que as linhagens de Telluraves já tinham divergido nesse intervalo [@ksepka2017], e a presença de Telluraves de tipo moderno no Cretáceo é considerada improvável [@mayr2014].

A pergunta mais disputada é a posição das corujas. Entre 2008 e 2026, quatro arranjos receberam suporte alto em pelo menos uma análise: Accipitrimorphae como grupo-irmão do restante de Afroaves, com as corujas junto de Coraciimorphae, sobretudo nas concatenações [@hackett2008; @jarvis2014; @reddy2017]; Hieraves como grupo-irmão do restante de Afroaves [@stiller2024]; Accipitrimorphae, só ou com as corujas, como grupo-irmão dos demais Telluraves [@prum2015; @houde2019; @gupta2025; @sapci2026]; e Hieraves como grupo-irmão de Australaves [@wu2024]. Reanálises de retroelementos e de árvores de genes sugerem reticulação antiga nessa região da árvore [@gatesy2022; @sapci2026].

A topologia também decide uma pergunta ecológica. @jarvis2014 observaram rapinantes nos ramos mais profundos de Australaves e de Afroaves e propuseram que o ancestral de Telluraves teria sido um predador. Essa leitura sustenta uma definição filogenética de rapinantes como as ordens Accipitriformes, Strigiformes, Cathartiformes, Falconiformes e Cariamiformes [@mcclure2019], e um mapeamento estocástico recente concluiu que o ancestral das aves terrestres era vertívoro [@mcclurejohnson2025].

Este trabalho pesou a evidência pela independência, e não pela contagem de árvores publicadas, que compartilham dados. Um nó foi tratado como resolvido apenas quando sobreviveu à troca de tipo de dado, de estimador, de número de locos e de amostragem de táxons. O objetivo foi estabelecer a posição filogenética das aves de rapina em Telluraves com quatro abordagens computacionais: (1) a reanálise das árvores de genes do B10K e das árvores depositadas na Open Tree of Life; (2) o posicionamento de fósseis e a medida de convergência em matrizes morfológicas publicadas; (3) a reconstrução do estado ancestral do hábito rapinante sob topologias alternativas; e (4) um conjunto filogenômico novo, independente da anotação e da ortologia do B10K.

# 2 MATERIAL E MÉTODOS

## 2.1 Reanálise das árvores de genes do B10K e da literatura depositada

As árvores de genes e de espécies de @stiller2024 vieram do repositório oficial do estudo [@stiller2024dados], em 18 de outubro de 2026. Foram usadas as árvores de máxima verossimilhança de 63.430 locos intergênicos de 1 kb, 14.972 éxons (1ª e 2ª posições do códon), 44.846 íntrons e 4.985 elementos ultraconservados (UCEs), além das 11 árvores de espécies publicadas. Como no artigo original, ramos com aLRT < 0,95 foram colapsados com o gotree 0.5.2 [@lemoine2021]. As árvores foram então podadas para 42 táxons focais: os 15 rapinantes do B10K, 2 Cariamiformes, 3 Psittaciformes, 3 Passeriformes, 13 Coraciimorphae e 6 grupos externos, com *Gallus gallus* como raiz de apresentação. Árvores com menos de quatro táxons focais foram descartadas, e ficaram 62.706 árvores intergênicas, 14.271 de éxons, 40.933 de íntrons e 4.919 de UCEs.

As frequências de quarteto foram calculadas no ASTER 1.25 [@zhang2025] com semente 12345. Para cada tipo de dado, o ASTRAL-IV pontuou a topologia publicada do B10K podada aos 42 táxons e fez uma busca livre nas árvores colapsadas; o wASTRAL híbrido, ponderado pelo aLRT [@zhangmirarab2022], buscou livremente nas árvores não colapsadas. A pontuação com os 363 táxons só foi concluída para as UCEs. Como a frequência normalizada q é uma média por quarteto de espécies e a busca otimiza a soma bruta de quartetos, oito topologias fixas de 42 táxons foram pontuadas para registrar o escore bruto. O efeito do colapso foi isolado em 1.209 árvores de UCE que contêm os 42 táxons, pontuadas com e sem colapso.

Em cada ramo foram aplicados três testes de qui-quadrado sobre as contagens efetivas de quartetos: o teste de politomia, com 2 graus de liberdade [@sayyari2018]; a igualdade das duas resoluções mais frequentes; e a simetria das duas resoluções menos frequentes, ambos com 1 grau de liberdade. A probabilidade local posterior (localPP) é a do ASTER [@sayyari2016]. O fator de concordância de genes (gCF) [@minh2020] foi calculado no IQ-TREE 3.1.4 [@wong2026] sobre as árvores não colapsadas.

A síntese 16.1 da Open Tree of Life [@hinchliff2015; @redelings2017] foi consultada pela API v3 para 24 nós-chave, e a síntese Aves 1.6 [@mctavish2025], pelo repositório público dos autores. As 28 árvores filogenéticas que sustentam, conflitam ou resolvem algum nó-chave foram baixadas, e em cada uma foram registrados o critério de clado estrito e a fração de 2.000 quartetos amostrados que sustenta cada resolução. As idades vieram da árvore datada do B10K, da árvore de consenso da Aves 1.6 e do TimeTree 5 [@kumar2022].

## 2.2 Posicionamento de fósseis e homoplasia em matrizes morfológicas

Foram analisadas dez matrizes morfológicas publicadas (Tabela 1). Polimorfismos e caracteres inaplicáveis foram recodificados como dado faltante, caracteres invariantes foram removidos e caracteres ordenados pelos autores entraram como não ordenados.

::tabela Tabela 1 – Matrizes morfológicas analisadas
| Matriz | Táxons × caracteres | Grupo externo |
|---|---|---|
| @ksepka2019 | 62 × 146 | *Eudromia elegans* |
| @mayr2011 | 33 × 100 | Tinamidae |
| @mayr2005 | 23 × 110 | Tinamidae |
| @mayr2026 | 42 × 105 | Anhimidae |
| @mayrkitchener2022 | 10 × 42 | Anhimidae |
| @mayr2016 | 13 × 40 | Anhimidae |
| @degrange2021 | 72 × 204 | *Meleagris gallopavo* |
| @mayrclarke2003 | 46 × 148 | *Apsaravis* |
| @chen2025 | 77 × 202 | *Gansus yumenensis* |
| @livezey2007 | 154 × 2.325 variáveis | *Lithornis* |
::fonte Fonte: elaborada pelo autor (2026).

As árvores foram inferidas no IQ-TREE 3.1.4 com o modelo Mk [@lewis2001], na forma MK+ASC+G4, com partições por número de estados e comprimentos de ramo proporcionais entre partições. O suporte veio de 1.000 réplicas de UFBoot com refinamento NNI [@hoang2018] e 1.000 réplicas de SH-aLRT [@guindon2010], com semente 20261017. Cada matriz foi analisada duas vezes: sozinha e com os táxons atuais presos à árvore principal do B10K, com os fósseis livres. A matriz de Mayr e Kitchener, com três táxons atuais, não recebeu arcabouço. Na matriz de Livezey e Zusi, o UFBoot da análise livre foi interrompido em 100 de 1.000 réplicas, e só a árvore de máxima verossimilhança foi usada.

A posição de cada fóssil foi lida nas 1.000 árvores UFBoot, depois da poda dos demais fósseis. Hipóteses de posição foram comparadas pelo teste AU [@shimodaira2002] com 10.000 reamostragens RELL [@kishino1990], sobre buscas restritas. Para os táxons atuais, o arcabouço molecular foi comparado com a Falconiformes clássica (Falconidae junto de Accipitridae) e, em duas matrizes, com Cathartidae junto de Ciconiidae.

A homoplasia foi medida pelo índice de consistência [@kluge1969] e pelo índice de retenção [@farris1989] de Fitch, por caráter, no phangorn 2.12.1 [@schliep2011]. Os caracteres rapinantes (bico em gancho, cera, tômio, unguais, hálux, dedos e tendões do pé, tarsometatarso e hipotarso) foram identificados pelas descrições publicadas, disponíveis para cinco matrizes. As diferenças entre classes foram testadas com 10.000 permutações de rótulos e pelo teste exato de Fisher. O morfoespaço usou as distâncias MORD [@lloyd2016] e de Gower [@gower1971], a PERMANOVA [@anderson2001] e o teste de Mantel parcial [@smouse1986], com 9.999 permutações. O repositório MorphoSource foi consultado (56 consultas, 604 mídias) para inventariar tomografias de fósseis rapinantes.

## 2.3 Reconstrução do estado ancestral

As árvores tiveram 202 famílias, com uma espécie por família tirada das 363 do B10K, e raiz em Galloanserae. Foram usadas seis topologias: a árvore principal datada do B10K (T1), a árvore de 94 mil locos (T1b), a árvore de éxons (T1c) e as espinhas de @prum2015 (T2), @jarvis2014 (T3) e @kuhl2021 (T4). A T4 foi transcrita à mão da figura em nível de ordem do artigo, porque o arquivo de árvores depositado não pôde ser baixado. Nas topologias T2 a T4, as relações dentro de cada ordem vieram de T1. Os comprimentos de ramo vieram das idades de T1, com espaçamento uniforme à maneira do BLADJ [@webb2008] e raiz em 88,43 milhões de anos (Ma), ou do método de Grafen com ρ = 1 [@grafen1989].

Nas árvores datadas foram enxertados quatro fósseis em dois cenários: *Primoptynx poliotauros* [@mayretal2020] no tronco de Strigiformes; Masillaraptoridae [@mayrkitchener2022] no tronco de Falconidae; Messelasturidae [@mayr2021; @mayr2025] no tronco de Psittacopasserae (cenário F1) ou de Strigiformes (F2); e Phorusrhacidae, representados por *Andalgalornis* [@degrange2010; @tambussi2012; @angst2013], no tronco de Cariamidae. Ao todo foram 24 árvores.

O caráter principal (C1) foi a proporção de espécies vertívoras ou necrófagas por família no AVONET [@tobias2022], com limiares de 10 % e 50 %. A sensibilidade foi medida com outros limiares, com a inclusão de predadores aquáticos (C2), com o EltonTraits 1.0 (C3) [@wilman2014], com o bico em gancho (C4), com os pés rapinantes (C5) e com uma síndrome de três estados (C6). As células de bico e pés foram codificadas pelo autor, com fonte primária específica só para parte delas [@fowler2009; @portugal2016].

Foram ajustados modelos Mk de taxas iguais (ER), de todas as taxas diferentes (ARD) e irreversíveis, por máxima verossimilhança no castor [@louca2018] com 16 partidas aleatórias, e reavaliados no phytools [@revell2024]. Duas condições de raiz foram comparadas: o prior de @fitzjohn2009 e a raiz de Neognathae fixada no estado não rapinante, justificada pelos grupos externos. Os modelos foram promediados pelo critério de informação de Akaike corrigido (AICc) [@hurvich1989], e cada condição recebeu 1.000 mapas estocásticos [@huelsenbeck2003]. Os caracteres C1, C3 e C6 foram repetidos nas árvores com fósseis com a dieta dos fósseis codificada como desconhecida. Uma parcimônia de Sankoff [@sankoff1975] serviu de referência. Ao todo rodaram 480 tarefas de ajuste, mais 72 com dieta fóssil desconhecida, com sementes derivadas de 20261017.

## 2.4 Conjunto filogenômico novo

### 2.4.1 Táxons e locos

O conjunto teve 57 táxons, com dez subfamílias de Accipitridae, as três subfamílias de Falconidae, quatro linhagens de corujas e oito ordens de Coraciimorphae. Quarenta e seis táxons vieram de anotações RefSeq do NCBI, obtidas com o NCBI Datasets [@oleary2024]. Onze vieram de montagens GenBank sem anotação, das quais as sequências codificantes (CDS) foram extraídas com o miniprot 0.18 [@li2023] a partir da proteína do parente anotado mais próximo: *Sarcoramphus papa*, *Sagittarius serpentarius*, *Pandion haliaetus*, *Elanus axillaris*, *Gypaetus barbatus*, *Leptodon cayanensis*, *Aegypius monachus*, *Pithecophaga jefferyi*, *Herpetotheres cachinnans*, *Ibycter americanus* e *Ninox strenua*. A taxonomia foi conferida com o TaxonKit [@shen2021].

Os 28.263 genes do galo foram reduzidos a 12.677 candidatos codificantes com ortólogo anotado, fora de famílias multigênicas, e embaralhados numa permutação de semente fixa. Dos primeiros 1.200, consultados no NCBI Orthologs, 686 passaram nos critérios de quadro de leitura, comprimento e ocupação. Os primeiros 300 com pelo menos 49 dos 57 táxons formaram o conjunto, com 95,0 % das células táxon × loco preenchidas e média de 54,2 táxons por loco.

Na extração pelo miniprot, registros com mudança de quadro, códon de parada interno, cobertura menor que 50 % ou identidade menor que 70 % viraram dado faltante, de 2,5 % (*Sarcoramphus*) a 23,0 % (*Pandion*). Os CDS foram extraídos com o gffread [@pertea2020]. Nas duas espécies com anotação GenBank independente, as proteínas extraídas foram traduzidas com o SeqKit [@shen2024] e comparadas às anotadas com o DIAMOND [@buchfink2021]. Coincidiram (identidade ≥ 98 % e cobertura ≥ 90 %) 92,3 % dos registros de *Gypaetus* (identidade média de 99,42 %) e 93,0 % dos de *Sarcoramphus* (99,75 %). Esse controle mede a fidelidade da extração, não a ortologia.

### 2.4.2 Inferência

A inferência usou um pipeline Nextflow [@ditommaso2017] com módulos nf-core [@ewels2020] em contêineres Docker. Cada loco foi alinhado com o MAFFT 7.520 em L-INS-i [@katoh2013] e limpo com o ClipKIT 2.11.4 no modo *smart-gap* [@steenwyk2020]. As árvores de genes de aminoácidos foram inferidas no IQ-TREE 3.1.4, com seleção de modelo pelo ModelFinder [@kalyaanamoorthy2017], 1.000 réplicas de UFBoot2, 1.000 de SH-aLRT, suporte aBayes [@anisimova2011] e semente 12345.

A supermatriz teve uma partição por loco, agrupadas pelo ModelFinder com fusão de partições [@chernomor2016; @lanfear2014; @lanfear2017] em 48 grupos. A primeira corrida foi interrompida durante o refinamento do UFBoot; a supermatriz foi refeita com os 48 modelos fixados, terminou com lnL de −2.230.233,545 e reproduziu a topologia da melhor árvore da corrida interrompida (distância de Robinson-Foulds 0). Sobre ela foram calculados o gCF das 300 árvores de genes e o fator de concordância de sítios por verossimilhança (sCFL) com 100 quartetos [@mo2023]. A árvore de espécies veio do wASTRAL ponderado pelo suporte aBayes [@zhangmirarab2022]. A assimetria entre as duas alternativas de cada nó foi testada por qui-quadrado sobre as contagens de árvores de genes e de sítios, já que sob ILS pura as duas têm a mesma frequência esperada [@minh2020].

Duas matrizes de sensibilidade tiveram árvores de genes e wASTRAL, sem supermatriz: as mesmas proteínas alinhadas com o FAMSA 2.4.1 [@deorowicz2016], com suporte aBayes, e as 1ª e 2ª posições do códon tiradas do alinhamento de aminoácidos. Treze hipóteses e a árvore do wASTRAL foram comparadas com a árvore da supermatriz pelo teste AU, com buscas restritas sobre a supermatriz e 10.000 reamostragens RELL.

## 2.5 Reprodutibilidade

Todas as chamadas de ferramenta foram registradas com parâmetros e horário num log de proveniência. As sementes foram fixas: 12345 no ASTER e no IQ-TREE dos braços 1 e 4; 20261017 no braço morfológico e no estado ancestral; 0 na permutação dos genes candidatos. As análises em Python 3.11.17 [@python2026] usaram numpy 2.4.6, pandas 2.3.3, scipy 1.17.1 e DendroPy 5.1.0 [@moreno2024]; as análises em R 4.6.1 [@rcore2026] usaram ape 5.8-1 [@paradis2019], phytools 2.5-2, castor 1.8.7 e phangorn 2.12.1. Os valores de p são nominais, sem correção para testes múltiplos.

# 3 RESULTADOS

## 3.1 Panorama por nó

A Tabela 2 resume o veredito para cada nó, combinando a literatura e as quatro análises deste trabalho. Dois nós rapinantes, Australaves e Accipitrimorphae, sobreviveram a todas as variações testadas. A base de Afroaves e o grupo-irmão de Telluraves não sobreviveram.

::tabela Tabela 2 – Veredito por nó para a posição das aves de rapina em Telluraves
| Nó | Veredito | Evidência central | Confiança |
|---|---|---|---|
| Telluraves | monofilético | 64 retroposons sem conflito; Open Tree 18 a favor e 0 contra; supermatriz UFBoot 100 | alta |
| Grupo-irmão de Telluraves | não resolvido | quatro propostas diferentes nos estudos de 2014 a 2024 | baixa |
| Australaves: Cariamiformes + (Falconiformes + Psittacopasserae) | sustentado | 176 retroposons; q de 0,381 a 0,441 em dado não codificante; Open Tree 15 a 0; supermatriz UFBoot 96; AU rejeita Cariamiformes + Falconiformes (p = 0,019) | alta |
| Psittacopasserae | sustentado | q de 0,375 a 0,405; Open Tree 17 a 0 | alta |
| Accipitrimorphae | sustentado | 192 retroposons; q de 0,892 a 0,969; gCF de 82 % a 95 %; Open Tree 16 a 0; supermatriz UFBoot 100 | muito alta |
| Sagittariidae + (Pandionidae + Accipitridae) | sustentado | q de 0,980 a 0,995; gCF de 94,8 % a 98,3 %; AU rejeita as alternativas por mais de 1.000 unidades de lnL | muito alta |
| Strigiformes: Tytonidae + Strigidae | sustentado | q de 0,997 a 0,999; gCF de 99,5 % a 99,8 % | muito alta |
| Base de Afroaves | não resolvido | quartetos do B10K 35,0 %, 34,6 % e 30,5 %; sCF 43,1 %, 32,5 % e 24,5 %; AU não rejeita Hieraves (p = 0,13) | baixa |
| Afroaves monofilético | provável | Open Tree 14 a 4; supermatriz UFBoot 93, mas sCF 33,7 % contra 34,2 % e 32,1 % | moderada |
| Falconidae: Herpetotherinae + (Polyborinae + Falconinae) | sustentado | 13 árvores-fonte, fração de quartetos 1,00; supermatriz UFBoot 100 | alta |
::fonte Fonte: elaborada pelo autor (2026), com dados de @suh2015, @stiller2024, da Open Tree of Life 16.1 e deste trabalho.
::nota Nota: q, frequência de quartetos da resolução indicada nas árvores de genes do B10K; gCF, fator de concordância de genes; sCF, fator de concordância de sítios no conjunto novo; contagens da Open Tree, árvores-fonte filogenéticas a favor e contra.

## 3.2 Telluraves e seu grupo-irmão

Nenhuma das 18 árvores-fonte filogenéticas da Open Tree ligadas a Telluraves o contradiz, e na Aves 1.6 a contagem é de 13 a 0. Na árvore datada do B10K, a coroa de Telluraves tem 63,8 Ma (intervalo de 95 % de 63,0 a 65,0), Afroaves 63,6 Ma, o nó Accipitrimorphae + Strigiformes 63,4 Ma, Australaves 63,1 Ma e Falconiformes + Psittacopasserae 62,6 Ma: internós sucessivos de 0,2 a 0,7 Ma. As medianas do TimeTree para os mesmos nós (de 65,8 a 71,2 Ma) não respeitam a topologia, pois vêm de conjuntos diferentes de estudos para cada par, e não foram usadas como estimativa.

## 3.3 Australaves e a posição dos falcões

Nos dados não codificantes do B10K, Falconiformes + Psittacopasserae venceu as duas alternativas: q = 0,381 contra 0,351 (Cariamiformes + Falconiformes) e 0,269 (Cariamiformes + Psittacopasserae) no intergênico, 0,428 contra 0,325 e 0,248 nos íntrons, e 0,406 contra 0,306 no wASTRAL ponderado. O escore bruto favoreceu o mesmo lado por 4,4 milhões de quartetos no intergênico. Nos éxons sem ponderação, a ordem se inverteu, com Cariamiformes + Falconiformes em 0,420 contra 0,316 (p = 2 × 10^−43^), e nas UCEs também (0,395 contra 0,362 com 42 táxons; 0,409 contra 0,355 com 363). Em todos os tipos de dado e estimadores, as duas resoluções minoritárias foram desiguais, sempre com Cariamiformes + Falconiformes acima de Cariamiformes + Psittacopasserae (p de simetria de 0,002 a 10^−133^). Das 11 árvores de espécies publicadas pelo B10K, todas têm Falconiformes + Psittacopasserae, exceto a de éxons, em que Australaves não é monofilético. As 15 árvores-fonte informativas da Open Tree têm Falconiformes + Psittacopasserae, sem nenhuma contra.

Dentro de Falconidae, Herpetotherinae (*Herpetotheres* + *Micrastur*) é grupo-irmão de Polyborinae + Falconinae em todas as 13 árvores-fonte informativas, com fração de quartetos 1,00, o mesmo arranjo de @fuchs2015.

## 3.4 Accipitrimorphae e Strigiformes

Na reanálise do B10K, Cathartidae como grupo-irmão de Accipitriformes no sentido estrito teve q de 0,892 a 0,969 em íntrons, intergênico e UCEs (gCF de 82,1 % a 95,1 %), e Pandionidae + Accipitridae, com Sagittariidae fora, teve q de 0,980 a 0,995 (gCF de 94,8 % a 98,3 %). As resoluções alternativas foram simétricas (p de simetria de 0,03 a 0,63). Os éxons deram os menores valores (0,556 e 0,673), sem mudar a resolução. Nenhuma das 12 árvores-fonte informativas contradiz esses nós. A ordem (Cathartidae, (Sagittariidae, (Pandionidae, Accipitridae))) coincide com a das UCEs de 134 genomas [@catanach2024], e a arquitetura cromossômica a confirma por outro caminho: Accipitridae e Pandionidae compartilham um plano genômico rearranjado, enquanto Cathartidae conserva a arquitetura ancestral [@canesin2024]. Em Strigiformes, Strigidae com *Tyto* fora teve q de 0,997 a 0,999 e gCF de 99,5 % a 99,8 %.

## 3.5 A base de Afroaves

No nó correspondente da árvore principal do B10K, Hieraves e corujas + Coraciimorphae receberam 35,0 % e 34,6 % dos quartetos, e a terceira resolução, 30,5 % [@stiller2024]. Ao pontuar a topologia publicada, Hieraves venceu em todos os tipos de dado: q = 0,374 contra 0,333 no intergênico, 0,377 contra 0,336 nos íntrons, 0,419 contra 0,304 nos éxons e 0,403 contra 0,316 nas UCEs (0,424 contra 0,303 com os 363 táxons). A busca livre do ASTRAL-IV, porém, devolveu corujas + Coraciimorphae no intergênico e nos íntrons. O escore bruto preferiu esse lado por 723.012 quartetos no intergênico e 174.401 nos íntrons, e o wASTRAL ponderado empatou (intergênico, 0,348 contra 0,343, p = 0,14; íntrons, 0,349 contra 0,345, p = 0,39). Nas 1.209 árvores de UCE com os 42 táxons, q e escore bruto concordaram nas árvores binárias e se separaram depois do colapso: a inversão veio do colapso dos ramos, não da falta de táxons.

Na Open Tree, Hieraves aparece por construção, porque o nó entra pela árvore do B10K, de maior prioridade, apesar de 4 árvores-fonte a favor e 14 contra (4 contra 9 na Aves 1.6). Pela fração de quartetos, corujas + Coraciimorphae foi majoritário em 14 árvores-fonte (10 de inferência primária) e Hieraves em 4 (3 primárias). A monofilia de Afroaves também vacilou: quatro árvores-fonte não a têm, e nas UCEs com 363 táxons Coraciimorphae + Australaves venceu Afroaves (q = 0,394 contra 0,342). Nos dados intergênicos e de íntrons, contudo, Afroaves venceu (q = 0,360 e 0,370), e a topologia com Accipitrimorphae como grupo-irmão de todo Telluraves perdeu 10,2 milhões de quartetos no intergênico.

As duas resoluções minoritárias da base de Afroaves foram desiguais em todos os tipos de dado (p de simetria de 0,007 a 10^−35^). O sinal por locus foi mínimo: o gCF do ramo Hieraves foi de 2,5 % no intergênico, 6,6 % nos íntrons, 3,9 % nas UCEs e 1,6 % nos éxons, com 87 % a 98 % das árvores decisivas discordando por polifilia. A Tabela 3 mostra que a literatura reproduz essa oscilação.

::tabela Tabela 3 – A base de Afroaves nos principais estudos filogenômicos
| Estudo | Dados | Método | Resultado |
|---|---|---|---|
| @hackett2008 | 19 locos, 169 espécies | ML | Accipitrimorphae na base; corujas + Coliiformes |
| @mccormack2013 | 1.541 UCEs, 33 táxons | Bayes e ML | (Accipitrimorphae, (corujas, Cavitaves)) |
| @jarvis2014 | 48 genomas | ExaML | Accipitrimorphae na base; corujas + Coraciimorphae |
| @prum2015 | 259 locos, 198 espécies | RAxML e ExaBayes | Accipitrimorphae na base de Telluraves |
| @suh2015 | 2.118 retroposons | parcimônia | corujas + Cavitaves; corujas + Accipitrimorphae quase tão forte |
| @reddy2017 | 235 espécies | RAxML | Accipitrimorphae na base; corujas + Coraciimorphae |
| @houde2019 | indels, 48 táxons | ASTRAL | corujas + Accipitrimorphae na base de Telluraves |
| @kuhl2021 | 3′-UTR, 429 espécies | IQ-TREE | Accipitrimorphae na base; corujas + (Coliiformes, Cavitaves) |
| @gatesy2022 | 2.113 retroelementos | ASTRAL | corujas + Cavitaves; excesso de corujas + Accipitrimorphae |
| @stiller2024 | 63.430 locos intergênicos | ASTRAL | Hieraves na base de Afroaves |
| @stiller2024 | idem | RAxML-NG | Accipitrimorphae na base; corujas + Coraciimorphae |
| @wu2024 | 124 genomas | NJst | Hieraves + Australaves |
| @gupta2025 | 62.413 árvores de genes | coalescente | Accipitrimorphae na base de Telluraves |
| @sapci2026 | árvores de genes do B10K | wASTRAL | Accipitrimorphae na base de Telluraves; sem Afroaves |
::fonte Fonte: elaborada pelo autor (2026), com base nos estudos citados.
::nota Nota: Os resultados de @hackett2008, @mccormack2013, @prum2015 e @reddy2017 foram lidos da compilação de @gatesy2022.

## 3.6 Fósseis, convergência e morfoespaço

Sem arcabouço, a morfologia montou um clado de rapinantes (falconiformes no sentido amplo mais corujas) com UFBoot de 78,7 % na matriz de @ksepka2019, 87,1 % na de @mayr2005 e 89,8 % na de @mayrclarke2003, e Accipitridae + Falconidae com 81,5 % a 92,1 %. Nenhuma matriz recuperou Accipitriformes no sentido molecular (até 1,5 %, exceto 12,8 % na de @chen2025), Hieraves (até 0,3 %), Australaves ou Eufalconimorphae (0,0 %). As árvores morfológicas recuperaram de 18 % a 42 % das bipartições do arcabouço molecular.

Sobre o arcabouço, cinco das sete matrizes testáveis rejeitaram a posição molecular de Falconidae em favor da Falconiformes clássica: @ksepka2019 (ΔlnL 45,3; p = 0,004), @mayr2011 (p = 0,0017), @mayr2005 (p = 0,0001), @mayrclarke2003 (p = 0,0003) e @livezey2007 (ganho de 595 unidades de lnL; p = 1 × 10^−39^). As exceções foram a matriz de catartídeos de @degrange2021 (p = 0,43) e a de cintura escapular e asa de @chen2025 (p = 0,46). Na matriz de Livezey e Zusi, os 221 caracteres rapinantes informativos tiveram índice de consistência menor na árvore molecular (0,426 contra 0,483; p = 0,013) e perderam mais consistência ao sair da árvore morfológica (ΔIC de −0,090 contra −0,053; p = 0,0008), sem diferença na própria árvore morfológica (p = 0,22). Entre os 161 caracteres que custam menos passos com Falconidae junto de Accipitridae, os rapinantes estiveram super-representados (razão de chances 1,78; teste de Fisher, p = 0,010), e 19 dos 29 eram de pé e garra: falanges dos dedos II e IV, unguais, hálux, tendão do flexor longo do hálux, canais do hipotarso e almofadas plantares.

No morfoespaço, rapinantes de clados moleculares diferentes ficaram mais próximos entre si do que pares sorteados em todas as oito matrizes analisadas (p de 10^−4^ a 0,038). A guilda explicou pouca variância além do clado na PERMANOVA (R^2^ ≤ 0,05), porque guilda e clado são quase colineares, mas o Mantel parcial guilda | clado foi significativo em @mayr2011 (r = 0,42) e em @livezey2007 (r = 0,53).

Na matriz de @mayrkitchener2022, *Masillaraptor* e *Danielsraptor* saíram como grupo-irmão de Falconidae em 77,8 % das réplicas UFBoot, e o teste AU não rejeitou nem a posição junto de *Cariama* (p = 0,10) nem o tronco de Falconidae + Cariamiformes (p = 0,25). Sem arcabouço, os messelasturídeos caíram junto das corujas (68 % na matriz de 2005) ou dentro de (Accipitridae, Falconidae, Strigiformes) (63 % na de 2011). Com o arcabouço, a atração se desfez: pan-Psittaciformes em 66,5 % na matriz de 2005, e tronco de Psittacopasserae em 41,8 % e pan-Eufalconimorphae em 99,0 % na matriz emendada de @mayr2026. O teste AU quase nada rejeitou para os messelasturídeos; a exceção foi o tronco de Australaves (p = 0,037) na matriz de 2005. Os halcionitídeos saíram pan-Psittacopasserae em 99,2 % na matriz de Ksepka, Grande e Mayr, com Accipitriformes (p = 0,007), Strigiformes (0,008), Hieraves (0,048) e Australaves (0,044) rejeitados, e os fororracídeos ficaram com *Cariama* em 100 %. *Dryornis pampeanus* e *Pampagyps imperator* caíram dentro da coroa de Cathartidae em 100 % das réplicas, e a hipótese que aproxima urubus e cegonhas custou 908 unidades de lnL na matriz de Livezey e Zusi (AU, p = 8 × 10^−56^). O inventário do MorphoSource não encontrou tomografia, aberta ou restrita, de *Masillaraptor*, *Danielsraptor*, *Messelastur*, *Tynskya*, halcionitídeos, corujas paleogênicas, *Bathornis* ou *Horusornis*.

## 3.7 Estado ancestral do hábito rapinante

Garras grandes, recurvas e preênseis surgiram três vezes, no tronco de Pandionidae + Accipitridae, no de Falconidae e no de Strigiformes, sem nenhuma perda, nas seis topologias, nas duas raízes, com e sem fósseis. Para a dieta, com a raiz ancorada no estado dos grupos externos, a probabilidade de Telluraves ser vertívoro foi de 0,00 a 0,03 em todas as topologias. Com o prior de FitzJohn, ela subiu para 0,82 a 0,90, mas o mesmo valor apareceu na raiz de Neognathae, e a diferença de log-verossimilhança entre as duas raízes foi de no máximo 1,5 unidade.

Com os fósseis enxertados, a probabilidade de Telluraves vertívoro subiu para 0,93 a 0,99 mesmo com a raiz não rapinante, sob um modelo ARD com perda de 73 a 80 vezes mais rápida que o ganho e de 8 a 10 perdas dentro da coroa. Com a dieta dos fósseis tratada como desconhecida, o efeito desapareceu (0,00). Com Hieraves (T1), houve uma única origem da dieta de vertebrados no ancestral de Accipitrimorphae + Strigiformes (P = 0,73 nesse nó); com as corujas junto de Coraciimorphae (T3, T4) ou Accipitrimorphae na base de Telluraves (T1b, T2), as origens foram separadas (P = 0,00 no nó-pai). Em todas as topologias, o ganho de Falconidae ficou no ramo da família, e não em Eufalconimorphae. O bico em gancho surgiu de 4 a 6 vezes só com os viventes. A topologia de éxons do B10K reuniu todas as origens numa só.

## 3.8 O conjunto filogenômico novo

A supermatriz de aminoácidos teve 215.366 sítios em 48 grupos de partições (Figura 1A). Telluraves, Accipitrimorphae, *Pandion* + Accipitridae e Falconiformes tiveram SH-aLRT, aBayes e UFBoot de 100, 1 e 100; Australaves, 100, 1 e 96; Eufalconimorphae, 99, 1 e 95; Afroaves, 98,6, 1 e 93; e corujas + Coraciimorphae, com Accipitrimorphae na base de Afroaves, 93,8, 1 e 90 (Tabela 4). Dentro de Accipitridae, *Elanus* saiu primeiro, depois (*Gypaetus*, *Leptodon*), (*Aegypius*, *Pithecophaga*), *Aquila*, *Harpia* e (*Astur*, (*Buteo*, *Haliaeetus*)), a mesma espinha de @catanach2024. *Herpetotheres* saiu como grupo-irmão de *Ibycter* + *Falco* com 100, 1 e 100.

::figura reports/artigo/figuras/figura1_arvores.png|14.5|Figura 1 – Árvores do conjunto filogenômico novo, com 57 táxons e 300 locos: (A) supermatriz de máxima verossimilhança; (B) árvore de espécies do wASTRAL
::fonte Fonte: elaborada pelo autor (2026).
::nota Nota: (A) supermatriz de aminoácidos, com UFBoot nos nós citados no texto; (B) wASTRAL ponderado por aBayes nas 300 árvores de genes de aminoácidos alinhadas com MAFFT, com probabilidade local. Cladogramas enraizados em *Gallus gallus*. Azul, Accipitrimorphae; laranja, Falconiformes; verde-água, Strigiformes; preto, Cariamiformes; cinza-escuro, demais Telluraves; cinza-claro, táxons fora de Telluraves.

Os fatores de concordância de sítios mostraram onde está o conflito (Figura 2). Accipitrimorphae teve sCF de 51,7 %, contra 24,7 % e 23,7 % para as alternativas, e Falconiformes, 84,1 %. Eufalconimorphae teve 41,8 %, contra 29,5 % para Cariamiformes + Falconiformes e 28,8 % para Cariamiformes + Psittacopasserae: as alternativas empataram (p = 0,81). Na base de Afroaves, as três resoluções receberam 43,1 % (corujas + Coraciimorphae), 32,5 % (Hieraves) e 24,5 % (Accipitrimorphae + Coraciimorphae), com as minoritárias desiguais (p = 0,0095) e o excesso do lado de Hieraves. Os sítios não decidiram a monofilia de Afroaves: 33,7 % contra 34,2 % e 32,1 %. Telluraves também mostrou alternativas desiguais (34,7 % e 19,6 %; p = 1,7 × 10^−6^), coerente com a incerteza sobre seu grupo-irmão.

::figura reports/artigo/figuras/figura2_scf.png|16|Figura 2 – Fatores de concordância de sítios (sCFL) de dez nós da árvore da supermatriz
::fonte Fonte: elaborada pelo autor (2026).
::nota Nota: Cada barra divide 100 % dos quartetos de sítios informativos entre o nó (azul) e as duas alternativas por troca de vizinho mais próximo (cinza). Sob separação incompleta de linhagens pura, as duas alternativas têm a mesma frequência esperada; o valor de p sob cada nome vem do qui-quadrado de simetria entre elas, nominal e sem correção. A linha tracejada marca 1/3, o valor esperado numa politomia.

As árvores de genes inverteram o lado de Australaves. A bipartição Cariamiformes + Falconiformes apareceu em 41 a 45 das 300 árvores, com ramo forte em 23, nas três matrizes, contra 6 a 10 para Eufalconimorphae e 3 a 6 para Cariamiformes + Psittacopasserae. Por isso o wASTRAL juntou Falconiformes e *Cariama* (localPP de 0,61 a 0,77), não recuperou Australaves e montou um clado com todas as linhagens rapinantes (localPP de 0,66 a 0,92; Figura 1B). Nos aminoácidos, o wASTRAL deu Hieraves (0,74 com MAFFT; 0,69 com FAMSA); nos códons 1 + 2, corujas + *Leptosomus* (0,61); e não recuperou Afroaves em nenhuma das três matrizes. As árvores de genes foram pouco resolvidas, com mediana de 12 % de ramos internos com SH-aLRT ≥ 80 e UFBoot ≥ 95. O gCF dos nós profundos ficou entre 0 % e 8 %; o de Accipitrimorphae foi de 7,7 %, porque 92 % das árvores decisivas não continham nem esse ramo nem suas alternativas.

::tabela Tabela 4 – Nós profundos no conjunto filogenômico novo
| Nó na supermatriz | Supermatriz | sCF: nó / alternativas | wASTRAL | Alternativa (p do teste AU) |
|---|---|---|---|---|
| Telluraves | 100/1/100 | 45,7 / 34,7 / 19,6 | 0,94 / – / – | – |
| Australaves | 100/1/96 | 53,1 / 26,0 / 20,9 | – / – / – | – |
| Eufalconimorphae | 99/1/95 | 41,8 / 29,5 / 28,8 | – / – / – | Cariamiformes + Falconiformes (0,019) |
| Afroaves | 98,6/1/93 | 33,7 / 34,2 / 32,1 | – / – / – | – |
| Strigiformes + Coraciimorphae | 93,8/1/90 | 43,1 / 32,5 / 24,5 | – / – / – | Hieraves (0,13); Accipitrimorphae + Coraciimorphae (< 0,001) |
| Accipitrimorphae | 100/1/100 | 51,7 / 24,7 / 23,7 | 1,00 / 1,00 / 1,00 | – |
| *Pandion* + Accipitridae | 100/1/100 | 72,5 / 15,9 / 11,6 | 1,00 / 1,00 / 1,00 | *Sagittarius* + Accipitridae; *Sagittarius* + *Pandion* (< 0,001) |
| Falconiformes | 100/1/100 | 84,1 / 10,1 / 5,9 | 1,00 / 1,00 / 1,00 | Falconiformes no sentido amplo (< 0,001) |
::fonte Fonte: elaborada pelo autor (2026).
::nota Nota: Supermatriz, SH-aLRT / aBayes / UFBoot; sCF, porcentagem de sítios a favor do nó e das duas alternativas, em ordem decrescente; wASTRAL, probabilidade local com MAFFT / FAMSA / DNA 1 + 2, com traço quando o nó não está na árvore.

O teste AU sobre a supermatriz rejeitou Cariamiformes + Falconiformes (ΔlnL 55; p = 0,019) e a árvore inteira do wASTRAL (p = 0,0007), mas não Hieraves (ΔlnL 34,8; p = 0,13). Rejeitou ainda Accipitrimorphae + Coraciimorphae (p < 0,001) e, por diferenças de 392 a 1.544 unidades de lnL (p < 0,001), a Falconiformes no sentido amplo em duas versões, corujas + bacurau, urubus + cegonha e as duas alternativas para *Sagittarius*.

# 4 DISCUSSÃO

## 4.1 O sinal de “aves de rapina” nos éxons é efeito de estimador e de tipo de dado

O agrupamento de falcões e seriemas tem histórico conhecido em dado codificante. Nas árvores de códons e de aminoácidos de @kuhl2021, Falconiformes e Cariamiformes formaram um agrupamento de aves de rapina, ao qual a árvore de aminoácidos juntou Strigiformes, enquanto as 3′-UTR recuperaram Australaves. A árvore de éxons de @jarvis2014 reuniu rapinantes diurnos e seriemas com 40 % de bootstrap, efeito que os autores atribuíram à conversão gênica enviesada para GC. No B10K, a árvore de éxons difere da principal em 38 de 360 ramos, contra 6 ou 7 nos outros tipos de dado, e 30,83 % dos locos de éxon correm risco de saturação, contra 0,07 % dos intergênicos [@stiller2024]; a recodificação RY desloca falcões e seriemas para perto de papagaios e passeriformes [@braun2021].

O conjunto novo reproduziu esse padrão dentro de um único conjunto de dados. Os sítios da supermatriz sustentaram Australaves, e as duas alternativas empataram. O excesso para Cariamiformes + Falconiformes apareceu na estimativa das árvores de genes, uma a uma, e não nos sítios. Somado no wASTRAL, esse excesso montou um clado de rapinantes que o teste AU rejeitou. A leitura mais consistente é que Cariamiformes + Falconiformes reflete viés sistemático de composição e saturação em sequência codificante, e não história. Com um único cariamiforme na amostra, porém, este conjunto não separa o viés de composição da atração do ramo longo de *Cariama*. A morfologia mantém a alternativa viva, porque o mosaico de caracteres dos masillaraptorídeos lembra tanto falconídeos quanto cariamiformes [@mayrkitchener2022].

## 4.2 A base de Afroaves está no limite de resolução dos dados

Nenhuma das três resoluções da base de Afroaves sobreviveu às quatro variações testadas. Duas armadilhas de leitura explicam parte da confusão na literatura. A localPP supõe resoluções minoritárias iguais; com a assimetria observada e dezenas de milhares de genes efetivos, ela vai a 0 ou a 1 com diferenças de q de 0,004. Uma localPP de 1 nesse ramo não indica resolução, e o bootstrap 100 da concatenação, no ramo em que concatenação e coalescente discordam com quartetos quase empatados, tampouco indica concordância. A segunda armadilha é o colapso dos ramos das árvores de genes, que separou a frequência normalizada do escore bruto e inverteu o lado vencedor.

A assimetria persistente entre as resoluções minoritárias aponta para mais do que ILS, porque sob ILS pura elas teriam a mesma frequência esperada [@minh2020]. Duas linhas independentes apontam na mesma direção: nos retroelementos, houve excesso de suporte para corujas + Accipitriformes entre as resoluções minoritárias, e uma rede com a linhagem das corujas como híbrida foi favorecida [@gatesy2022]; nas árvores de genes do B10K, o método Phlag encontrou blocos genômicos com distribuição anômala nesse ramo [@sapci2026]. No conjunto novo, o excesso também ficou do lado de Hieraves. A preponderância de análises independentes inclina para corujas + Coraciimorphae, com Accipitrimorphae na base de Afroaves, enquanto Hieraves conta com a análise coalescente mais bem amostrada [@stiller2024] e com @wu2024. Hieraves não é um clado demonstrado. A representação mais fiel aos dados é uma politomia mole entre Accipitrimorphae, Strigiformes e Coraciimorphae, provavelmente marcada por reticulação antiga.

## 4.3 Falconiformes era uma guilda, e o ancestral predador não está demonstrado

Como os falcões ficam em Australaves e os gaviões em Accipitrimorphae em todas as topologias que sobreviveram ao exame, o aparato de caça dos dois evoluiu de forma independente. Os pés rapinantes surgiram pelo menos três vezes contando as corujas, e a velha ordem Falconiformes era uma guilda, não uma linhagem. A homoplasia medida na matriz de Livezey e Zusi concentra-se justamente nos caracteres de pé e garra, e é a topologia molecular que os torna homoplásticos. O mesmo mecanismo explica o comportamento dos fósseis rapinantes do Eoceno: sem arcabouço, os messelasturídeos são atraídos pelas corujas e pelos falconiformes; com arcabouço, a atração se desfaz, coerente com a cautela de @mayr2021 sobre sua afinidade com Psittacopasserae.

Se a predação de vertebrados é herança comum de Telluraves é outra pergunta. Os viventes respondem que não: com a raiz ancorada nos grupos externos, a probabilidade de um ancestral vertívoro ficou perto de zero, e o resultado positivo sob o prior de FitzJohn reapareceu na raiz de Neognathae, onde nenhum grupo externo o sustenta. Os fósseis invertem a resposta, mas por um caminho parcialmente circular, porque sua dieta é inferida do bico e das garras. O ancestral predador de Telluraves é plausível, mas não demonstrado. A lista de ordens de @mcclure2019 é robusta a todas as topologias; a narrativa ancestral, não. Testes de irreversibilidade por reconstrução ancestral são frágeis a esse tipo de premissa [@goldberg2008].

O registro fóssil restringe parte das datações. Uma falange ungueal de águia-pescadora do Oligoceno inferior da Alemanha [@mayr2006] exige que Pandionidae e Accipitridae já estivessem separados no início do Oligoceno, e @kuhl2021 situam os fósseis de Pandioninae entre 32 e 34 Ma. As corujas têm o registro mais antigo entre as rapinantes de coroa, do Paleoceno médio, e *Primoptynx poliotauros*, do Eoceno inferior, tinha garras aumentadas como as de um acipitrídeo [@mayretal2020]. *Ypresiglaux michaeldanielsi* indica que as corujas desenvolveram adaptações rapinantes antes das especializações da vida noturna [@mayrkitchener2023]. Nenhum desses fósseis de coruja, porém, está numa matriz pública, e por isso eles não puderam ser posicionados aqui.

## 4.4 Limitações

O conjunto novo tem quatro limites. Só tem éxons, o tipo de dado mais sujeito ao artefato discutido acima. *Cariama*, *Sagittarius* e *Pandion* têm um táxon cada, e a atração de ramos longos neles não pode ser descartada. A independência de dado bruto é parcial, porque várias montagens são as mesmas de @jarvis2014 e *Herpetotheres* é montagem do B10K. As árvores de genes são pouco resolvidas, e o gCF baixo dos nós profundos reflete falta de resolução por locus, não conflito. Ficaram pendentes as supermatrizes das matrizes DNA 1 + 2 e FAMSA, a limpeza com trimAl e uma reanálise sem os 150 ramos terminais anômalos achados nas árvores de aminoácidos.

A reconstrução do estado ancestral foi feita por família. *Pandion* saiu como não vertívoro porque o AVONET o classifica como predador aquático, e incluir essa categoria faz o resultado depender da espinha de Neoaves. A topologia T4 foi transcrita à mão de uma figura em nível de ordem, e as células de bico e pés foram codificadas pelo autor, com fonte primária só para parte delas. No braço morfológico, as matrizes de Mayr compartilham autoria e provavelmente caracteres, e não somam evidência independente.

# 5 CONSIDERAÇÕES FINAIS

A posição das aves de rapina em Telluraves se resolve em dois níveis. No primeiro, os fatos estão estabelecidos: Falconiformes é grupo-irmão de Psittacopasserae, Cariamiformes fecha Australaves, e Accipitrimorphae é monofilético na ordem (Cathartidae, (Sagittariidae, (Pandionidae, Accipitridae))). Esses nós sobreviveram à troca de tipo de dado, de estimador, de número de locos e de amostragem, e bastam para concluir que os pés rapinantes evoluíram pelo menos três vezes e que a antiga Falconiformes era uma guilda.

No segundo nível, a base de Afroaves não é um problema de falta de dados. De 63 mil a 159 mil locos, retroposons e indels já foram aplicados, e o resultado oscila com o estimador porque três linhagens se separaram em cerca de 0,2 Ma, com discordância assimétrica compatível com reticulação. Mais locos do mesmo tipo dificilmente decidirão esse nó. Modelos que admitam rede, o mapeamento do sinal por região genômica e amostragem mais densa de corujas e Coraciimorphae são os caminhos com maior chance de decidi-lo. Na nomenclatura, Hieraves e mesmo Afroaves devem ser usados como hipóteses, não como clados estabelecidos.

O conjunto de 300 locos deste trabalho não desempatou o que a literatura deixa aberto, mas mostrou, dentro de um só conjunto de dados, o mecanismo do conflito: a concatenação e o teste AU dão Australaves, enquanto as árvores de éxons, estimadas uma a uma, puxam falcões e seriema para junto e, no somatório coalescente, montam um clado de rapinantes que os sítios não sustentam. Um clado de aves de rapina que apareça em árvores de genes codificantes deve ser lido primeiro como efeito do estimador e do tipo de dado. Na paleontologia, o gargalo é o acesso: não há matriz com corujas-tronco nem tomografia aberta de um único rapinante do Paleógeno, embora sejam justamente as corujas do Paleoceno e do Eoceno as que poderiam indicar de que lado da tricotomia a linhagem nasceu.

#* REFERÊNCIAS

::referencias

#* DISPONIBILIDADE DE DADOS E CÓDIGO

Os scripts, o log de proveniência e as tabelas de resultados estão no repositório público https://github.com/patoangolano/bio-stack. Uma cópia dedicada ao projeto, em https://github.com/patoangolano/filogenia-rapinantes, está privada na data desta versão. Os dados brutos são públicos e foram obtidos dos repositórios citados nos Métodos.

#* DECLARAÇÃO DE USO DE INTELIGÊNCIA ARTIFICIAL

O desenho das análises, a escrita dos scripts, a execução e a redação deste texto foram feitos com assistência do Claude Code (Anthropic), sob direção do autor. Nenhum dado biológico foi gerado por inteligência artificial. As árvores de genes e de espécies, as matrizes morfológicas, os traços ecológicos e as montagens genômicas são dados públicos dos repositórios citados. As exceções são a topologia T4, transcrita à mão de uma figura publicada, e as células de bico e pés da reconstrução do estado ancestral, codificadas pelo autor. Os DOIs das referências foram resolvidos no doi.org durante a preparação do texto.
