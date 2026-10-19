# Estado ancestral do hábito rapinante em Telluraves sob topologias alternativas

Agente evo-devo · 2026-10-18 · tudo in silico, dados públicos · diretório de trabalho
`C:\Users\Manec\bio-stack\runs\rapinantes\ancestral\` (árvores, caracteres, scripts, resultados, figuras).

## 0. Pergunta polarizada e resposta curta

Jarvis et al. (2014) escrevem no resumo: *"we infer the common ancestor of core landbirds to have been an
apex predator"*. A forma polarizada da pergunta: o hábito rapinante é **ancestral** em Telluraves e foi
**perdido** em Coraciimorphae, Psittacopasserae, Cariamidae etc., ou é **derivado** e foi **ganho** de modo
independente em Accipitrimorphae, Strigiformes e Falconidae? E a resposta muda com os nós instáveis
(posição de Strigiformes, monofilia de Afroaves, posição de Telluraves em Neoaves)?

Resposta, nível por nível:

- **Pés rapinantes: derivados, 3 ganhos independentes** (Pandionidae+Accipitridae, Falconidae, Strigiformes),
  zero perdas, em todas as 6 topologias, 2 raízes, com e sem fósseis. **[robusto]**
- **Dieta de vertebrados, só viventes:** ancestral em Telluraves **só** sob modelos de perda dominante
  (só-perda/ARD) com raiz de FitzJohn — e então o estado "rapinante" se estende à raiz de Neognathae
  (Galloanserae + Neoaves, 88 Ma), o que nenhum grupo externo sustenta. Com a raiz ancorada no estado do grupo
  externo, P(Telluraves rapinante) = 0,00–0,03 em todas as topologias. **[não sustentado por viventes]**
- **Com fósseis rapinantes do Paleógeno** (Masillaraptoridae no tronco de Falconidae, Phorusrhacidae no
  tronco de Cariamidae, coruja-tronco *Primoptynx*): P(Telluraves rapinante) = 0,93–0,99 mesmo com a raiz em 0,
  em T1, T1b, T2, T3, T4. Esse efeito **desaparece por inteiro** (0,00) quando a dieta dos fósseis é tratada
  como desconhecida. **[condicional: depende da dieta inferida da morfologia fóssil]**
- **Topologia** pesa em dois pontos. (i) Hieraves: com Strigiformes irmã de Accipitrimorphae (T1, Stiller 2024)
  há uma única origem no ancestral comum delas (P = 0,73); com corujas junto de Coraciimorphae (T2, T3, T4) ou
  Accipitriformes na base de Telluraves (T1b, T2) são origens separadas (P no nó-pai = 0,00). (ii) A topologia
  de éxons de Stiller (T1c, "clado rapinante", provável artefato) colapsa tudo numa origem e zera o sinal
  ancestral até com fósseis (0,03–0,06). **[condicional à topologia]**

## 1. Árvores

Todas em nível de família: 202 famílias, uma espécie representante por família, das 363 espécies de Stiller
et al. (2024). Grupos externos pedidos: Strisores (7), Aequornithes (16), Charadriiformes (16), Columbaves (6) e
Galloanserae (8, raiz). Telluraves: 149 famílias (120 de Passeriformes). Fora: Palaeognathae, Mirandornithes,
Gruiformes, Opisthocomiformes, Phaethontimorphae.

| Código | Fonte | Arquivo | Construção | Espinha de Telluraves (ordens) | RF vs T1 | nós interpolados |
|---|---|---|---|---|---|---|
| T1 | Stiller et al. 2024, ASTRAL 63K intergênico, datada (MCMCtree) | ERDA `04_timetree/main_alternative_mod.tre` | poda direta | (Acc+Strig) + Coraciimorphae; Cariam + (Falc + (Psitt + Pass)) | 0 | 0 |
| T1b | Stiller 2024, ASTRAL 94K (Telluraves igual em 128K/153K/159K) | ERDA `03_species_trees/94K.tre` | poda direta, raiz em Galloanserae | Acc + resto; Colii + (Strig + Cavitaves) | 8 | 4 |
| T1c | Stiller 2024, ASTRAL de éxons (1ª+2ª posições) | ERDA `03_species_trees/exon_c12.tre` | poda direta | Falc + (Lepto + (Cariam + (Acc + Strig))) — "clado rapinante" | 56 | 28 |
| T2 | Prum et al. 2015 | `Avian-TimeTree.tre` (198 spp.; cópia no ERDA de Stiller, `EDF6a/`) | híbrida* | Acc + (Australaves + (Strig + Coraciimorphae)) | 10 | 5 |
| T3 | Jarvis et al. 2014, TENT ExaML | `Chronogram01.TENT.ExAML.tre` (48 spp.; ERDA, `EDF6a/`) | híbrida* | Afroaves = Acc + (Strig + Coraciimorphae) | 12 | 6 |
| T4 | Kuhl et al. 2021, 3'-UTR | **codificada à mão** da Fig. 2 (PMC7783168, `msaa191f2.jpg`) | híbrida* | igual a T3; Telluraves irmã de Aequornithes+Phaethontimorphae | 14 | 7 |

\* Híbrida: espinha entre ordens (e entre grupos externos) da fonte; relações entre famílias dentro de cada
ordem, e dentro de Strisores, Aequornithes, Charadriiformes e Galloanserae, tomadas de T1. Isso isola o efeito
dos nós entre ordens, que é o que varia entre os estudos.

- O `TREES.tar.gz` de Kuhl et al. existe no Dryad (doi:10.5061/dryad.ngf1vhhpx) mas não foi baixável: a API
  exige token e a página serve um desafio de prova de trabalho em JavaScript. Por isso T4 saiu da figura (nível
  de ordem, legível; a Fig. 3A em nível de família tem resolução insuficiente no PMC).
- As árvores de Prum e Jarvis vêm da redistribuição no arquivo de Stiller (pasta da Extended Data Fig. 6a). A
  topologia de Jarvis bate com o resumo publicado. Não foram conferidas contra os repositórios originais
  (suplemento da Nature, GigaDB 10.5524/101041).
- Comprimentos: **(a) idades de Stiller 2024** — um nó cujo conjunto de famílias forma clado em T1 recebe a
  idade de T1; os demais são espaçados uniformemente entre o ancestral fixo e o descendente fixo mais antigo
  (lógica do BLADJ; Webb et al. 2008). Raiz = 88,43 Ma. Os nós profundos de Telluraves caem todos entre 63,1 e
  63,8 Ma. **(b) Grafen (ρ = 1)**, escalado para a mesma idade de raiz. Com 120 passeriformes, o Grafen deixa
  Afroaves com 9,7 Ma e Australaves com 55 Ma: é uma perturbação forte, não uma alternativa realista.
- Arquivos: `arvores/T*_E_datada.tre`, `T*_E_grafen.tre`, `T*_F1.tre`, `T*_F2.tre`, `arvores.rds`,
  `resumo_topologias.tsv`, `espinha_telluraves_por_topologia.tsv`; figura `figuras/fig1_espinhas_telluraves.png`.

### Fósseis enxertados como pontas (só nas árvores datadas)

O ponto de enxerto fica a meio caminho entre o nó-pai do ramo-alvo e o maior valor entre a idade do nó-filho,
a idade da ponta e a idade mínima da linhagem.

| Ponta | Posição F1 | Posição F2 | Idade da ponta | Enxerto em T1 | Trófico / bico / pés | Referências |
|---|---|---|---|---|---|---|
| *Primoptynx poliotauros* (coruja-tronco) | tronco de Strigiformes | idem | 55,0 Ma (Willwood, 54,5–55) | 59,5 Ma | 1 inferido / ? / 1 | Mayr, Gingerich & Smith 2020 (pé tipo Accipitridae, garras D-I e D-II maiores) |
| Masillaraptoridae (*Danielsraptor*, *Masillaraptor*) | tronco de Falconidae (coroa 31,97 Ma) | idem | 55,0 Ma (London Clay, ~55, aproximada) | 58,8 Ma | 1 inferido / 1 / ? | Mayr & Kitchener 2022 ("long-legged, raptorial birds", tronco de Falconiformes) |
| Messelasturidae (*Messelastur*, *Tynskya*) | tronco de Psittacopasserae | tronco de Strigiformes | 55,0 Ma | 62,1 Ma (F1) / 57,6 Ma (F2) | ? / 1 / ? | Mayr 2021 (crânio "raptor-like", pés semi-zigodáctilos; com arcabouço molecular, irmão de Psittaciformes+Passeriformes); F2 por Mayr 2005 e Mayr 2025 (caracteres partilhados com Strigiformes/Falconiformes) |
| Phorusrhacidae (*Andalgalornis* como representante) | tronco de Cariamidae (coroa 7,71 Ma) | idem | 6,0 Ma; idade mínima da linhagem 42 Ma (*Eleutherornis*, Luteciano superior) | 52,6 Ma | 1 / 1 / 0 (inferência) | Degrange et al. 2010; Tambussi et al. 2012; Angst et al. 2013 |

Ogygoptynx (Paleoceno; Rich & Bohaska 1976) ficou de fora como ponta: só há o tarsometatarso. Tabela de
entrada: `caracteres/fosseis_codificacao_manual.tsv`; enxertos: `arvores/enxertos_fosseis.tsv`.

## 2. Caracteres

São 202 famílias viventes mais 4 fósseis. A fonte de cada célula está em `caracteres/caracteres_familia.tsv`.
Um polimorfismo entra como prior de ponta (½, ½). A matriz que vai para o Mk é `caracteres/matriz_mk.tsv`.

- **Correspondência de taxonomia:** espécie de Stiller → espécie eBird/Clements do AVONET (312 por nome exato,
  15 por tabela de correspondência, 2 por gênero, 1 por sinônimo manual). Família de Stiller → conjunto de
  famílias Clements. Exemplo: Ramphastidae de Stiller = Ramphastidae + Capitonidae + Semnornithidae + Lybiidae +
  Megalaimidae; Bucerotidae inclui Bucorvidae. Oito famílias de Stiller que desmembram famílias Clements
  (Dendrocolaptidae, Scleruridae, Tachurididae…) herdam a proporção da família-mãe; todas são insetívoras ou
  granívoras.
- **C1 trófico (primário).** AVONET (Tobias et al. 2022; figshare 10.6084/m9.figshare.16586228 v7, arquivo
  34480856, md5 confere). Proporção p das espécies com `Trophic.Niche` ∈ {Vertivore, Scavenger}; cada nicho exige
  ≥ 60 % da dieta, pela definição do AVONET.
  - Limiares (inferior/superior): **T50** 0,10/0,50 (1 se p ≥ 0,5; 0&1 se 0,10 ≤ p < 0,5); T25 0,05/0,25;
    T75 0,25/0,75; bin50 sem polimorfismo.
  - Em T50 ficam com 1: Accipitridae (0,82), Cathartidae (1,0), Falconidae (0,60), Sagittariidae (1,0),
    Tytonidae (1,0). Strigidae (0,31) fica 0&1.
  - **Pandionidae fica 0**: no AVONET a águia-pescadora é "Aquatic predator". É uma consequência da definição,
    não da biologia, e é o que o C2 corrige.
- **C2:** C1 + "Aquatic predator" (T50). Corrige *Pandion*, mas põe 1 em boa parte de Aequornithes e
  Charadriiformes (piscivoria).
- **C3:** EltonTraits 1.0 (Wilman et al. 2014; `esapubs.org/archive/ecol/E095/178/BirdFuncDat.txt`).
  Proporção de espécies com Vend+Vect+Vunk+Scav ≥ 50 %, limiar T50.
- **C4 bico em gancho rapinante (estrito):** gancho usado para despedaçar presa ou carcaça de vertebrado. 1 em
  Accipitridae, Pandionidae, Sagittariidae, Cathartidae, Falconidae, Strigidae, Tytonidae, Cariamidae, Laniidae
  e Stercorariidae; 0&1 em Artamidae, Vangidae e Malaconotidae. **C4b amplo:** gancho morfológico de qualquer
  função; acrescenta Psittaciformes, Procellariiformes, Pelecanidae, Phalacrocoracidae, Fregatidae e
  Balaenicipitidae, e Laridae como 0&1.
- **C5 pés rapinantes** (garras grandes, recurvas, preênseis): 1 em Accipitridae, Pandionidae, Falconidae,
  Strigidae e Tytonidae (Fowler et al. 2009). 0 em Sagittariidae (mata por chute; Portugal et al. 2016),
  Cathartidae, Cariamidae e Psittaciformes.
- **C6 síndrome** (3 estados): 0 = nenhum componente; 1 = um ou dois de (C1_T50, C4, C5); 2 = os três.
- **Limite de fonte:** as células de bico e pés marcadas `ANALISTA` são codificação do analista a partir de
  conhecimento ornitológico geral, sem verificação célula a célula em fonte primária. Só Fowler 2009 e
  Portugal 2016 sustentam células específicas. As dietas fósseis "1" são inferências de morfologia feitas pelos
  autores citados.

## 3. Modelos

- **Mk binário:** ER, ARD, irreversível só-perda (1→0) e irreversível só-ganho (0→1). Em caráter binário, SYM = ER
  (uma única taxa), por isso SYM só entra no C6. **Mk de 3 estados (C6):** ER, SYM, ARD, só-perda (só para estado
  menor) e só-ganho.
- **Ajuste:** MLE por `castor::fit_mk` com 16 partidas, depois reavaliada no `phytools::fitMk(fixedQ=)`.
  - Em T1, castor e phytools chegam ao mesmo ótimo: logL −18,8046 no ARD e −35,684 no ARD do C6.
  - Reparo de aninhamento: se o ARD (ou o SYM) ficou abaixo de um modelo aninhado nele, recebe a Q do aninhado.
    Isso aconteceu, por exemplo, com o ARD de T1, cujo ótimo está na fronteira q01 = 0.
- **Seleção:** AICc com n = número de pontas; marginais em média ponderada pelos pesos de AICc.
- **Duas raízes:**
  - (i) prior de FitzJohn et al. (2009);
  - (ii) raiz de Neognathae fixada em 0, justificada pelo grupo externo: nenhuma família de Palaeognathae ou
    Galloanserae é rapinante. Com essa raiz o modelo só-perda é inviável.
- **Marginais:** `phytools::ancr` nos nós Neognathae, Neoaves, Telluraves, Australaves, Afroaves,
  Eufalconimorphae, Accipitrimorphae, Strigiformes, Coraciimorphae e no **nó-pai** de Accipitrimorphae e de
  Strigiformes. A identidade desse nó-pai é reportada por topologia.
- **Mapeamento estocástico:** `make.simmap`, **1000 mapas** por condição, Q fixa na MLE. Rodou para o melhor
  modelo (as duas raízes) e para o ER como referência (raiz de FitzJohn), nos caracteres C1_T50, C4, C5 e C6, em
  18 árvores datadas (com e sem fósseis).
  - Sementes: `set.seed(20261017 + 10·id + j)` dentro de cada processo-filho, RNG padrão; o valor de cada mapa
    está em `tabela_transicoes.tsv`. Nenhum `np.random.seed`.
  - 4 núcleos (cluster PSOCK, `chunk.size = 1`).
- **Total:** 480 tarefas (24 árvores × 10 caracteres × 2 raízes) + 72 da sensibilidade "dieta fóssil
  desconhecida"; 0 erros.
- **Parcimônia de referência:** Sankoff (phangorn), custos iguais e ganho = 2× perda. A orientação da matriz de
  custo foi conferida num exemplo de brinquedo.
- **Software:** R 4.6.1, ape 5.8.1, phytools 2.5.2, castor 1.8.7, phangorn 2.12.1; o renv do projeto não precisou
  de instalação nova.

## 4. Resultados

### 4.1 Tabela principal — C1_T50 (AVONET), árvores datadas, média ponderada por AICc, P(estado rapinante)

Célula: **raiz FitzJohn / raiz 0 / raiz 0 + fósseis F1**. Um * marca nó não monofilético na topologia (o valor
é do MRCA).

| Topologia | Telluraves | Australaves | Afroaves | Eufalconimorphae | Accipitrimorphae | nó-pai de Accipitrimorphae (identidade) |
|---|---|---|---|---|---|---|
| T1 Stiller 63K | 0,82 / 0,00 / 0,93 | 0,82 / 0,00 / 0,93 | 0,82 / 0,00 / 0,93 | 0,82 / 0,00 / 0,93 | 0,94 / 0,82 / 1,00 | 0,92 / 0,73 / 1,00 (Acc+Strig) |
| T1b Stiller 94K | 0,90 / 0,00 / 0,99 | 0,90 / 0,00 / 0,99 | *0,90 / 0,00 / 0,99 | 0,90 / 0,00 / 0,99 | 0,92 / 0,38 / 0,99 | 0,90 / 0,00 / 0,99 (Telluraves) |
| T1c Stiller éxons | 0,59 / 0,03 / 0,04 | *0,59 / 0,03 / 0,04 | *0,59 / 0,03 / 0,04 | *0,59 / 0,03 / 0,04 | 0,95 / 0,88 / 1,00 | 0,93 / 0,83 / 1,00 (Acc+Strig) |
| T2 Prum 2015 | 0,88 / 0,00 / 0,99 | 0,88 / 0,00 / 0,99 | *0,88 / 0,00 / 0,99 | 0,88 / 0,00 / 0,99 | 0,91 / 0,38 / 0,99 | 0,88 / 0,00 / 0,99 (Telluraves) |
| T3 Jarvis 2014 | 0,86 / 0,00 / 0,98 | 0,86 / 0,00 / 0,98 | 0,86 / 0,00 / 0,98 | 0,86 / 0,00 / 0,98 | 0,89 / 0,36 / 0,99 | 0,86 / 0,00 / 0,98 (Afroaves) |
| T4 Kuhl 2021 | 0,88 / 0,00 / 0,99 | 0,88 / 0,00 / 0,99 | 0,88 / 0,00 / 0,99 | 0,88 / 0,00 / 0,99 | 0,91 / 0,35 / 0,99 | 0,88 / 0,00 / 0,99 (Afroaves) |
| *raiz Neognathae* | 0,82–0,90 (T1c 0,58) / 0 fixo / 0 fixo | | | | | |

**Por modelo** (`resultados/tabela_C1_T50_por_modelo.tsv`), mesmo padrão nas 6 topologias, P(Telluraves = 1):

| Modelo | Raiz FitzJohn | Raiz 0 |
|---|---|---|
| ER | 0,00 | 0,00 |
| só-ganho | 0,00 | 0,00 |
| ARD | 1,00 (T1c 0,04) | 0,00 (T1c 0,04) |
| só-perda | 1,00 | inviável |

- **Pesos de AICc com raiz de FitzJohn:** só-perda 0,58–0,66, ARD 0,22–0,32, ER 0,07–0,15.
- **Com raiz 0:** ER/ARD dividem o peso (T1: ARD 0,51; T1c: ARD 0,77).
- **Diferença de logL entre as raízes** (melhor modelo de cada uma), só viventes: 0,11 (T1), 1,25–1,54 (T1b, T2,
  T3, T4) e −0,03 (T1c). Os dados mal distinguem "raiz rapinante" de "raiz não rapinante".
- **Leitura:** o "ancestral rapinante" sob FitzJohn é o mesmo número na raiz de Neognathae e em Telluraves.
  O modelo põe o estado na raiz da árvore inteira, não numa origem localizada em Telluraves.

### 4.2 Efeito dos fósseis e da codificação da dieta fóssil (nó Telluraves, raiz 0, média AICc)

| Caráter | Sem fósseis | F1 | F1, dieta fóssil "?" | F2 | F2, dieta fóssil "?" |
|---|---|---|---|---|---|
| C1_T50 | 0,00 (T1c 0,03) | 0,93–0,99 (T1c 0,04) | **0,00** (T1c 0,03) | 0,93–0,99 (T1c 0,04) | 0,00 |
| C3 EltonTraits | 0,00 | 0,82–0,85 (T1c 0,06) | **0,00** | 0,82–0,85 | 0,00 |
| C4 bico estrito | 0,00 | 0,94–0,99 (T1c 0,00) | — | **0,18 (T1), 0,63–0,78** (T1c 0,00) | — |
| C5 pés | 0,00 | 0,00 | — | 0,00 | — |
| C6 P(≥1 componente) | 0,00 | 0,82–0,99 (T1c 0,00) | 0,82–0,99* | 0,06 (T1), 0,37–0,55 | 0,06–0,55* |
| C6 P(síndrome completa) | 0,00 | 0,00 (T3: 0,89†) | | 0,00–0,48 | |

\* No C6 com dieta fóssil "?" o bico e os pés observados dos fósseis continuam informando.
† Ótimo instável: o ARD de 6 parâmetros tem duas soluções com ΔlogL 0,08, uma pondo o estado 1 e a outra o
estado 2 em Telluraves; P(≥1) = 0,98 nas duas.

- **Modelos com fósseis (C1, raiz 0):** o ARD leva peso 0,94–0,99, com perda de 73 a 80× a taxa de ganho
  (q01 ≈ 0,00022; q10 ≈ 0,016–0,018 por Ma). No ER, Telluraves continua 0.
- **O ancestral rapinante com fósseis exige duas coisas ao mesmo tempo:** (a) dietas fósseis codificadas como
  rapinantes e (b) um modelo assimétrico de perda dominante.
- **O bico ancestral depende de onde fica Messelasturidae:** tronco de Psittacopasserae (F1) contra tronco de
  Strigiformes (F2).

### 4.3 Ganhos e perdas na coroa de Telluraves (1000 mapas; média [IC 95 %])

| Caráter / condição | T1 | T1b | T1c | T2 | T3 | T4 |
|---|---|---|---|---|---|---|
| C1, sem fósseis, FitzJohn, só-perda (melhor) | 0 G / 13,4 [7–21] P | 0 / 14,0 [7–21] | 0 / 11,2 [4–21] | 0 / 13,8 [7–21] | 0 / 14,0 [7–21] | 0 / 13,8 [7–21] |
| C1, sem fósseis, ER (referência) | 2,6 [2–5] G / 0,9 [0–2] P | 4,2 [3–5] / 0,5 | 3,1 [2–5] / 0,8 | 4,2 [3–6] / 0,5 | 4,2 [3–6] / 0,5 | 4,3 [3–6] / 0,5 |
| C1, sem fósseis, raiz 0, melhor | ARD 2,3 [2–4] G / 1,3 [1–2] P | ER 4,2 / 0,5 | ARD 1,3 [0–3] / 3,2 [1–5] | ER 4,2 / 0,5 | ER 4,3 / 0,5 | ER 4,3 / 0,5 |
| C1, F1, raiz 0, ARD | 0,2 [0–1] G / 8,7 [5–17] P | 0,2 / 10,0 [6–18] | 1,0 [0–2] / 3,4 [3–5] | 0,2 / 8,9 [5–17] | 0,2 / 8,9 [6–16] | 0,2 / 8,8 [5–17] |
| C4 bico, sem fósseis (só-ganho/ER) | 4,1 [4–6] G / 0 P | 5,7 [5–7] / 0 | 2,6 [2–4] / 0,5 | 5,7 / 0 | 5,7 / 0 | 5,7 / 0 |
| C4 bico, F1, ARD | 1,3 [1–3] G / 4,5 [2–8] P | 1,3 / 5,6 [3–9] | ER 3,1 / 1,0 | 1,3 / 4,6 | 1,3 / 4,7 | 1,3 / 4,9 |
| C5 pés, todas as condições com raiz 0 | 3,0–3,2 [3–4] G / 0 P | idem | idem | idem | idem | idem |

Tabelas completas: `tabela_transicoes.tsv`, também para a árvore toda e para o C6. Figura: `fig4`.

### 4.4 Sensibilidade (nó Telluraves, média AICc; sem fósseis, a não ser que indicado)

- **Limiar trófico** (T25, T75, bin50), raiz 0: 0,00 em todas as topologias (T1c ≤ 0,04).
  - Raiz FitzJohn: T25 vai de 0,28 (T1) a 0,90; T75 de 0,26 (T1) a 0,79; bin50 de 0,86 a 0,92. O limiar
    mexe com a raiz de FitzJohn, não com a raiz ancorada.
- **EltonTraits (C3):** raiz FitzJohn 0,78–0,83; raiz 0, 0,00. Concorda com o AVONET.
- **Predador aquático incluído (C2), raiz 0:** T3 0,93, T4 0,54, T2 0,09, T1/T1b/T1c 0,00. Depende da espinha
  de Neoaves: em Jarvis e Kuhl, Telluraves é irmã de Aequornithes (piscívoras), e o predador "de vertebrados
  incluindo peixe" vira ancestral do clado conjunto.
- **Bico em gancho amplo (C4b, inclui papagaios e aves marinhas), raiz 0:** T4 0,96, T3 0,93, T2 0,58,
  T1/T1b 0,02, T1c 0,00. É o caráter mais dependente de topologia.
- **Grafen:** raiz 0, 0,01–0,10 em todos os caracteres. Raiz FitzJohn: de 0,36 a 0,99 (T3/T4 0,99), um efeito
  das idades absurdas de Afroaves sob Grafen.
- **Parcimônia (Sankoff):**
  - C1 com custos iguais: Telluraves = 0 em todas as topologias só com viventes (3–4 passos); com F1, 0/1 em
    T2, T3 e T4.
  - Ganho = 2× perda: Telluraves ambíguo (0/1) só com viventes e 1 com fósseis (T1c continua 0/1).
  - C4 com custos iguais: Telluraves = 1 em T2, T3 e T4 (a espinha de corujas com Coraciimorphae favorece o
    gancho ancestral), 0/1 em T1 e T1b, 0 em T1c.
  - C5: 0 em toda parte com custos iguais.
  - Parcimônia e Mk divergem no bico: o Mk põe as mudanças nos ramos terminais longos.

## 5. Leitura polarizada: em que ramo está a mudança

1. **Pés rapinantes:** três ganhos, nos ramos-tronco de Pandionidae+Accipitridae (depois da separação de
   Cathartidae e Sagittariidae), de Falconidae e de Strigiformes. O caráter não existe no ancestral de nenhum
   nó profundo. A convergência é o estado mais bem sustentado do trabalho.
2. **Dieta de vertebrados, só viventes e raiz ancorada:**
   - Em T1 o ganho fica no ramo-tronco de Hieraves (Accipitrimorphae + Strigiformes, ~63,4 Ma), com P = 0,73
     no ancestral comum.
   - Em T2, T3, T4 e T1b ficam ganhos separados no tronco de Accipitrimorphae e no de Strigiformes (P no
     nó-pai = 0,00), mais o de Falconidae.
   - Em todas as topologias o ganho de Falconidae fica no ramo próprio da família, não em Eufalconimorphae
     (0,00).
3. **"Ancestral rapinante" de Telluraves (Jarvis 2014):** reproduz-se só se se aceitar ao mesmo tempo uma raiz
   que torna Galloanserae + Neoaves rapinantes, ou fósseis do Paleógeno codificados como predadores com um
   modelo de perda ~75× mais rápida que o ganho. No segundo caso a mudança cai no **ramo-tronco de Telluraves**
   (que termina em 63,8 Ma; o nó-pai varia com a topologia), seguida de 8–10 perdas dentro da coroa.
4. **Bico em gancho:** só com viventes, 4–6 ganhos independentes. Com Messelasturidae no tronco de
   Psittacopasserae (F1), vira ancestral de Telluraves e se perde 4–6 vezes (Coraciimorphae, Psittacopasserae
   etc.). Com Messelasturidae na linhagem das corujas (F2), a probabilidade cai para 0,18–0,78.

## 6. O que estes dados não podem dizer

- **A ordem das mudanças dentro de um ramo.** Os ramos internos de Telluraves duram < 1 Ma (63,1–63,8 Ma);
  dizer se a dieta precedeu o bico, ou o bico os pés, num mesmo ramo-tronco está fora do alcance do Mk e de
  qualquer outro método. O C6 não testa sequência.
- **A trajetória histórica real.** As marginais são probabilidades condicionais a um modelo homogêneo no tempo,
  numa árvore fixa e sem incerteza de datação. Não houve integração sobre árvores nem sobre idades.
- **Discriminar a raiz.** A diferença de log-verossimilhança entre "raiz rapinante" e "raiz não rapinante" é
  ≤ 1,5 unidade só com viventes. O resultado de FitzJohn é uma propriedade do modelo, não um sinal dos dados.
  Goldberg & Igić (2008) mostram por que testes de irreversibilidade por ASR são frágeis.
- **A dieta de fósseis.** É inferida da forma do bico e das garras. Usar isso como evidência de que o ancestral
  era predador é parcialmente circular, e a sensibilidade "dieta ?" mostra que **todo** o efeito dos fósseis
  sobre o C1/C3 vem dessa inferência. As posições de Messelasturidae (incerta pelos próprios autores) e de
  *Primoptynx* (dentro do tronco de corujas, sem resolução) e as idades aproximadas da London Clay não foram
  integradas como incerteza: entram como dois cenários (F1, F2).
- **A amostragem e a codificação.**
  - Nível de família: polimorfismos internos (caracarás, *Pernis*, *Gypohierax*, pequenos Strigidae
    insetívoros) viram proporção, e Pandionidae sai 0 na definição AVONET.
  - Bico e pés são codificação do analista, com fontes primárias só para parte das células.
  - Passeriformes (120 pontas quase todas 0) dominam a estimação das taxas.
  - Subfamílias de Accipitridae e Falconidae não foram separadas: não mudariam os nós profundos, mas
    mudariam a contagem de perdas dentro das famílias.
- **As árvores T2–T4.**
  - São híbridas: só a espinha entre ordens é da fonte.
  - T4 foi codificada à mão de uma figura em nível de ordem.
  - T2 e T3 vêm da redistribuição de Stiller.
  - Nenhuma considera discordância gene-espécie, que é justamente o que torna instáveis os nós de corujas e
    Accipitriformes (ver a variação entre T1, T1b e T1c dentro do próprio B10K).
- **O C6.** O ARD de 3 estados tem superfície plana (soluções com ΔlogL < 0,1 e estados ancestrais
  diferentes). Use P(≥1), não P(estado 2).

## 7. Arquivos (`C:\Users\Manec\bio-stack\runs\rapinantes\ancestral\`)

- `downloads/MANIFESTO_downloads.tsv`: URL e md5 de cada insumo. `downloads/`: ERDA de Stiller, AVONET,
  EltonTraits, figuras de Kuhl.
- `caracteres/`:
  - `caracteres_familia.tsv`: proporções e codificações, com fonte por célula.
  - `matriz_mk.tsv`.
  - `morfologia_codificacao_manual.tsv` e `fosseis_codificacao_manual.tsv`.
  - `familias_pontas.tsv` e `pareamento_especies_stiller_avonet.tsv`.
- `arvores/`: 24 Newick (6 topologias × datada / Grafen / F1 / F2), `arvores.rds`, `resumo_topologias.tsv`,
  `enxertos_fosseis.tsv`.
- `scripts/`:
  - `01_caracteres.py` (uv) e `02_arvores.R`.
  - `03_mk.R` (modelos, ASR, simmap) e `03b_mk_fosQ.R` (dieta fóssil desconhecida).
  - `04_resumo_figuras.R` e `05_parcimonia.R`.
- `resultados/`:
  - `tabela_principal_C1_T50.tsv`, `tabela_C1_T50_por_modelo.tsv`.
  - `tabela_sensibilidade_media_AICc.tsv`, `tabela_melhor_modelo_contagem.tsv`, `tabela_transicoes.tsv`.
  - `ajustes_modelos*.tsv`, `marginais_nos*.tsv`, `simmap_*.tsv`, `parcimonia_mpr.tsv`.
  - `jobs/` (480 RDS) e `jobs_fosQ/` (72).
- `figuras/`:
  - `fig1_espinhas_telluraves.png`, `fig2_mapa_calor_P_rapinante.png`, `fig2b_efeito_fosseis_telluraves.png`.
  - `fig3_T1_C1_T50_raiz_fitzjohn.png`, `fig3b_T1_C1_T50_raiz_zero.png`.
  - `fig4_ganhos_perdas_telluraves.png`.
- `logs/`: `03_mk_sessionInfo.txt` e logs das rodadas. A primeira rodada foi interrompida por desequilíbrio de
  carga e retomada sem refazer as tarefas prontas; a semente é função do id, então o resultado é o mesmo.

## 8. Referências (DOI conferido no Crossref ou na Europe PMC; sem DOI quando não verificado)

- Angst D, Buffetaut E, Lécuyer C, Amiot R (2013) PLoS ONE 8:e80357. doi:10.1371/journal.pone.0080357
- Degrange FJ, Tambussi CP, Moreno K, Witmer LM, Wroe S (2010) PLoS ONE 5:e11856. doi:10.1371/journal.pone.0011856
- Degrange FJ et al. (2015) J Vertebr Paleontol 35:e912656. doi:10.1080/02724634.2014.912656
- FitzJohn RG, Maddison WP, Otto SP (2009) Syst Biol 58:595–611. doi:10.1093/sysbio/syp067
- Fowler DW, Freedman EA, Scannella JB (2009) PLoS ONE 4:e7999. doi:10.1371/journal.pone.0007999
- Goldberg EE, Igić B (2008) Evolution 62:2727–2741. doi:10.1111/j.1558-5646.2008.00505.x
- Grafen A (1989) Phil Trans R Soc B 326:119–157. doi:10.1098/rstb.1989.0106
- Huelsenbeck JP, Nielsen R, Bollback JP (2003) Syst Biol 52:131–158. doi:10.1080/10635150390192780
- Jarvis ED et al. (2014) Science 346:1320–1331. doi:10.1126/science.1253451 · dados: Jarvis et al. (2015) GigaScience 4. doi:10.1186/s13742-014-0038-1; GigaDB 10.5524/101041
- Kuhl H et al. (2021) Mol Biol Evol 38:108–127. doi:10.1093/molbev/msaa191 · Dryad doi:10.5061/dryad.ngf1vhhpx
- Louca S, Doebeli M (2018) Bioinformatics 34:1053–1055. doi:10.1093/bioinformatics/btx701
- Mayr G (2005) J Vertebr Paleontol 25:635–645. doi:10.1671/0272-4634(2005)025[0635:TPOAPP]2.0.CO;2
- Mayr G (2011) J Syst Palaeontol 9:159–171 (Messelasturidae + Halcyornithidae, ?Pan-Psittaciformes)
- Mayr G (2021) PalZ 95:337–357. doi:10.1007/s12542-020-00541-8
- Mayr G (2025) J Morphol 286:e70070. doi:10.1002/jmor.70070
- Mayr G, Gingerich PD, Smith T (2020) J Vertebr Paleontol 40:e1769116. doi:10.1080/02724634.2020.1769116
- Mayr G, Kitchener AC (2022) J Vertebr Paleontol 41(6):e2083515. doi:10.1080/02724634.2021.2083515
- Portugal SJ, Murn CP, Sparkes EL, Daley MA (2016) Curr Biol 26:R58–R59. doi:10.1016/j.cub.2015.12.004
- Prum RO et al. (2015) Nature 526:569–573. doi:10.1038/nature15697
- Revell LJ (2024) PeerJ 12:e16505. doi:10.7717/peerj.16505
- Rich PV, Bohaska DJ (1976) Smithson Contrib Paleobiol 27:87–93. doi:10.5479/si.00810266.27.87
- Stiller J et al. (2024) Nature 629:851–860. doi:10.1038/s41586-024-07323-1 · ERDA: https://erda.ku.dk/archives/341f72708302f1d0c461ad616e783b86/published-archive.html
- Tambussi CP, de Mendoza R, Degrange FJ, Picasso MB (2012) PLoS ONE 7:e37701. doi:10.1371/journal.pone.0037701
- Tobias JA et al. (2022) Ecol Lett 25:581–597. doi:10.1111/ele.13898 · figshare doi:10.6084/m9.figshare.16586228.v7
- Webb CO, Ackerly DD, Kembel SW (2008) Bioinformatics 24:2098–2100. doi:10.1093/bioinformatics/btn358
- Wilman H et al. (2014) Ecology 95:2027. doi:10.1890/13-1917.1
