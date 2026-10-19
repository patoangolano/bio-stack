# Reanálise de dados publicados e síntese quantitativa: posição das aves de rapina em Telluraves

Agente `phylo-analyst`, 2026-10-17/18. Escopo: reanálise das árvores de genes do B10K (Stiller et al. 2024), meta-análise de topologias sobre as árvores-fonte da Open Tree of Life (síntese 16.1) e da síntese Aves 1.6 (McTavish et al. 2025), e idades no TimeTree 5. O conjunto novo a partir de genomas do NCBI é de outro agente e não foi tocado aqui.

Todos os números abaixo saem de arquivos em `C:\Users\Manec\bio-stack\runs\rapinantes\reanalise\`, gerados pelos scripts em `C:\Users\Manec\bio-stack\scripts\rapinantes\reanalise\`. Todo DOI/URL foi conferido por ferramenta nesta sessão (Crossref, Europe PMC, resolução do DOI ou a própria API).

Abreviações: A = Accipitrimorphae (Cathartidae + Sagittariidae + Pandionidae + Accipitridae); S = Strigiformes; C = Coraciimorphae (Coliiformes, Leptosomiformes, Trogoniformes, Bucerotiformes, Coraciiformes, Piciformes); U = Australaves (Cariamiformes + Falconiformes + Psittacopasserae); K = Cariamiformes; F = Falconiformes (Falconidae); PsP = Psittacopasserae. "Hieraves" = A + S. q = frequência de quarteto em volta do ramo; EN = número efetivo de árvores de genes; localPP = probabilidade local do ASTRAL.

---

## 1. Resumo

1. **Interior de Accipitrimorphae e de Strigiformes: resolvido.** Cathartidae é irmã de Accipitriformes s.s. (Sagittariidae + (Pandionidae + Accipitridae)), e Strigidae é monofilético com Tytonidae fora: q = 0,89–0,995 nas árvores de genes intergênicas, de íntrons e de UCEs do B10K (éxons 0,56–0,85), gCF 82–99,8 %, todas as 11 árvores publicadas do B10K, 12 de 12 árvores-fonte informativas da Open Tree (Strigidae: 14 de 14) e nenhum conflito na Open Tree 16.1 nem na Aves 1.6. Falconidae: Herpetotheriinae irmã de Polyborinae + Falconinae em 13 de 13 árvores-fonte, fração de quartetos 1,00 em todas (o B10K não tem como testar).
2. **Australaves com Eufalconimorphae (Falconiformes irmã de Psittacopasserae) se sustenta, com conflito estruturado.** Dados não codificantes do B10K dão q = 0,381 (intergênico) e 0,428 (íntrons) contra 0,351 e 0,325 para Cariamiformes + Falconiformes; 15 de 15 árvores-fonte; todas as árvores publicadas do B10K menos a de éxons. Mas éxons e UCEs, sem peso, preferem Cariamiformes + Falconiformes (q = 0,420 e 0,395–0,409), e em todos os tipos de dado essa alternativa é muito mais frequente que a outra (p de simetria até 10⁻¹³³): um sinal sistemático que junta falcões e seriemas em dado codificante/conservado, o mesmo agrupamento de "aves de rapina" que Kuhl et al. (2021) descrevem nas suas árvores de sequência codificante. Psittacopasserae: q = 0,375–0,405; 18 de 18 árvores-fonte.
3. **A base de Afroaves não está resolvida.** Hieraves (Accipitrimorphae + Strigiformes) × Strigiformes + Coraciimorphae: no B10K completo, 0,3495 × 0,3456 (artigo). Nesta reanálise com 42 táxons, o q normalizado prefere Hieraves em todos os tipos de dado (0,374–0,419 × 0,304–0,336), mas o escore bruto de quartetos — o que a busca do ASTRAL otimiza — prefere Strigiformes + Coraciimorphae no intergênico (+723 012 quartetos) e nos íntrons (+174 401), e a árvore livre do ASTRAL-IV sai com S+C nos dois; o wASTRAL ponderado empata (intergênico S+C 0,348 × A+S 0,343, p = 0,14). A inversão q × escore vem do colapso dos ramos das árvores de genes (testado). A concatenada do artigo dá S+C (BS 100), e os conjuntos maiores do B10K (94K a 159K locos) também. Nas árvores-fonte: S+C 14 (10 de inferência primária) × Hieraves 4 (3 primárias, todas de 2024–2025). A Open Tree e a Aves 1.6 mostram Hieraves só porque a árvore do B10K tem a maior prioridade na síntese (4 árvores a favor, 14 contra).
4. **Nem a monofilia de Afroaves é segura.** A árvore UCE publicada do B10K e a pontuação UCE com 363 táxons põem Hieraves fora de Afroaves, ((A,S),(C,U)); as árvores 94K a 159K põem Accipitrimorphae irmã do resto de Telluraves, como Prum et al. 2015 e Feng et al. 2020; as árvores de Wu et al. 2024 e Claramunt et al. 2025 depositadas na Open Tree têm ((A,S),(C,U)). Quatro de 18 árvores-fonte informativas não têm Afroaves.
5. **Discordância assimétrica e gCF muito baixo.** Nos ramos contestados, as duas resoluções minoritárias têm frequências desiguais (sob ILS pura seriam iguais); na base de Afroaves o gCF é 1,6–6,6 %, com 87–98 % das árvores de genes decisivas discordando por polifilia (árvores de 1 kb raramente recuperam os grupos inteiros). Reticulação antiga (Şapcı et al. 2026), erro de estimativa das árvores de genes e viés de composição continuam todos possíveis; estes dados não separam um do outro.
6. **Idades.** Na árvore datada do B10K, as divergências da base de Telluraves caem todas entre 62,2 e 63,4 Ma: Accipitridae–Cathartidae 62,2 [61,3–63,0]; Falconidae–Psittaciformes 62,6; Falconidae–Cariamidae 63,1; Strigiformes–Accipitriformes 63,4 — logo depois do K–Pg, e um intervalo de cerca de 1 Ma para três ou quatro divergências é o cenário em que se espera muita ILS. As medianas do TimeTree são mais antigas (65,8–71,2 Ma) e nem respeitam a topologia (Falconidae–Psittaciformes 71,2 > Falconidae–Cariamidae 63,2).

---

## 2. O que foi feito

### 2.1 B10K (Stiller et al. 2024): só árvores de genes e árvores publicadas

- Repositório oficial achado pelo "Data availability" do artigo (texto completo PMC11111414 via Europe PMC): DOI 10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a, que resolve para o arquivo ERDA da Univ. de Copenhague (https://erda.ku.dk/archives/341f72708302f1d0c461ad616e783b86/published-archive.html). O índice `published-files.json` lista 3 340 arquivos, 397 GB no total (289 GB são janelas de 10 kb do alinhamento de genoma inteiro).
- Baixados 0,81 GB, nenhum alinhamento: árvores de genes dos 63 430 locos intergênicos de 1 kb (com suporte aLRT e a versão colapsada em aLRT < 0,95 que o artigo usou no ASTRAL), 14 972 éxons (posições 1+2 do códon), 44 846 íntrons e 4 985 UCEs; as 11 árvores de espécies publicadas (63K ASTRAL, 63K concatenada RAxML-NG, 63K GOBS, éxons, íntrons, UCE, 94K, 128K, 153K, 159K); o teste de politomia do artigo; a árvore datada MCMCtree; e a planilha de dados-fonte da Fig. 2. Tamanhos conferidos contra o índice do ERDA; sha256 em `b10k/brutos_manifesto/sha256.txt` (o ERDA não publica md5). Brutos no WSL em `~/rapinantes/reanalise/b10k/brutos`.
- Nomes de folha traduzidos para os 363 nomes da árvore 63K (`preparar_b10k.py`): éxons usam código GENESP de 6 letras (6 exceções resolvidas à mão: cormorões com gênero atual Urile/Nannopterum e os dois *Pomato-*; a troca possível entre os dois *Pomato-* não altera nenhum ramo fora de Passeriformes); UCEs pelo `uce_name_map.txt`. Rótulo de suporte ")/0,999" convertido para ")0,999" (o wASTRAL lê só o primeiro número e aborta em rótulo começado por "/").
- Conjunto focal de 42 táxons (`b10k/taxa_focais.tsv`): todos os 15 rapinantes do B10K (*Cathartes aura*; *Sagittarius serpentarius*; *Pandion haliaetus*; 5 Accipitridae: *Aquila chrysaetos*, *Circaetus pectoralis*, *Haliaeetus albicilla*, *H. leucocephalus*, *Spizaetus tyrannus*; *Falco cherrug*, *F. peregrinus*, *Herpetotheres cachinnans*; *Tyto alba*, *Strix occidentalis*, *Glaucidium brasilianum*, *Ciccaba nigrolineata*), 2 Cariamiformes, 3 Psittaciformes, 3 Passeriformes (Acanthisittidae, suboscine, oscine), 13 Coraciimorphae (2 Coliiformes, *Leptosomus*, 2 Trogoniformes, 2 Bucerotiformes, 3 Coraciiformes, 3 Piciformes) e 6 externos (*Opisthocomus*, *Charadrius*, *Egretta*, *Chaetura*, *Columba*, *Gallus*, este como raiz de apresentação). O B10K não tem *Caracara*, *Micrastur*, nenhum Polyborinae nem *Phodilus*: a base de Falconidae não é testável nele.
- Poda árvore a árvore (`podar_arvores.py`; o `gotree prune` 0.5.2 aborta o arquivo inteiro quando uma árvore fica com 2 folhas). Árvores com < 4 focais saem: ficaram 62 706/63 430 intergênicas, 14 271/14 972 éxons, 40 933/44 846 íntrons, 4 919/4 985 UCEs (média de 39,2 / 33,7 / 28,3 / 39,1 focais por árvore).
- ASTER v1.25 (pixi): por tipo de dado, (a) ASTRAL-IV v1.25.4.8 sem peso nas árvores colapsadas em aLRT < 0,95, busca livre; (b) ASTRAL-IV pontuando (`-C`) a topologia 63K publicada podada aos focais, para que todos os tipos de dado sejam medidos nos mesmos ramos; (c) wASTRAL v1.25.3.8 híbrido ponderado pelo aLRT (`-L --mode 1`) nas árvores não colapsadas; (d) ASTRAL-IV pontuando a topologia 63K com os 363 táxons e as árvores colapsadas completas — só para as UCEs: no intergênico (63 430 árvores × 363 táxons, `-u 3`) o ASTRAL-IV passou de 14 GB e foi morto pelo OOM killer do WSL; os valores com 363 táxons do intergênico vêm dos dados-fonte da Fig. 2a do artigo. Sempre `-u 3` (q1–q3, f1–f3, pp1–pp3, EN no rótulo e o `freqQuad.csv` com as quadripartições explícitas), `--seed 12345`, 8 threads. Ajuda da versão instalada conferida antes (`-u 2` = detalhado; `-u 3` = detalhado + freqQuad).
- Testes sobre f1–f3 e EN (`sintese_b10k.py`): politomia (Sayyari & Mirarab 2018; qui-quadrado 2 g.l. contra EN/3); "1ª × 2ª" (qui-quadrado 1 g.l. entre as duas resoluções mais frequentes, o teste que Stiller et al. usaram, com BH, para dizer que o ramo é estatisticamente indistinguível); simetria das duas menos frequentes (qui-quadrado 1 g.l.; sob ILS pura as duas alternativas têm a mesma frequência esperada). O ASTER não tem teste de politomia embutido (o `-t 10` era do ASTRAL-III em Java); o do artigo foi lido do arquivo deles.
- Escore total de quartetos de 8 topologias fixas (ASTRAL-IV `-C`, `pontuar_alternativas.sh`), porque a frequência q que o ASTRAL imprime é média por quarteto de espécies e o objetivo do ASTRAL é a soma bruta de quartetos; com árvores de genes colapsadas as duas podem apontar para lados opostos (aconteceu; o mecanismo foi testado, ver 3.1).
- gCF (IQ-TREE 3.1.4, `--gcf`, sem alinhamento: só gCF; sCF/sCFL exigiriam os alinhamentos, que não foram baixados por instrução) das árvores de genes não colapsadas sobre a topologia 63K podada.
- Topologia enraizada dos quatro grupos A, S, C, U em cada árvore publicada e em cada árvore desta reanálise (15 possíveis; `raiz_telluraves`).

### 2.2 Open Tree of Life (API v3) e síntese Aves 1.6

- `tree_of_life/about`: síntese **opentree16.1** (2025-12-20), OTT 3.7draft3, 1 931 estudos / 2 064 árvores no total.
- `tnrs/match_names` (contexto Birds) para 51 exemplares; `tree_of_life/induced_subtree`; para cada um de 24 nós-chave, `tree_of_life/mrca` dos táxons que o definem + verificação de que os táxons que deviam ficar fora estão fora (se não estão, o nó **não existe** na síntese e o MRCA devolvido é outro nó); `tree_of_life/node_info` (supported_by, conflicts_with, resolves, partial_path_of, terminal); `studies/find_studies` para citação e ano. Scripts `opentree_nos.py`.
- Síntese de aves: confirmado por Europe PMC/Crossref que McTavish et al. 2025 existe (PNAS 122(18):e2409658122, doi:10.1073/pnas.2409658122; "Aves 1.3"). O repositório que o artigo indica (github.com/McTavishLab/AvesData) já tem **Aves 1.6** (commits de jan–ago/2026); usei a 1.6: `annotations.json` (suporte/conflito por nó), `labelled_supertree_ottnames.tre`, `rank_collection.json` (prioridade das 338 árvores-fonte) e a árvore datada `summary_dated_clements.nex` (taxonomia Clements 2025). Script `aves_nos.py`.
- Meta-análise das árvores-fonte (`opentree_votos.py`): todas as 28 árvores filogenéticas que sustentam, conflitam ou resolvem algum nó-chave (exclui as taxonômicas OTT, eBird e Clements) foram baixadas da phylesystem (`study/<id>/tree/<id>.tre?tip_label=ot:ottId`). Cada folha recebe um grupo pela taxonomia OTT (`taxonomy/subtree`); folha dentro de Falconidae/Strigiformes/Accipitriformes sem subgrupo atribuído (gêneros fósseis) fica fora em vez de virar "externo". Duas correções medidas na OTT 3.7: os quiriquiris (*Falco sparverius* etc.) estão no gênero *Cerchneis* e há um nó "Bubonidae" separado de Strigidae; sem as correções, *F. sparverius* virava grupo externo e Strigidae saía não monofilético em 4 árvores grandes. Por árvore e por pergunta: (i) critério de clado estrito; (ii) fração de 2 000 quartetos amostrados (um táxon de cada grupo + um externo) que sustenta cada resolução — robusto a um táxon fora do lugar; (iii) topologia enraizada dos quatro grupos A, S, C, U por 2 000 amostras de 5 táxons.

### 2.3 TimeTree 5

- API REST `https://timetree.temple.edu/api/pairwise/<taxid>/<taxid>/json` respondeu (versão 20260820). Taxids do NCBI por E-utilities com filtro `Aves[Subtree]` (Sagittariidae tem homônimo fora das aves: sem o filtro saem 2 taxids). Para cada par: idade "estimada" do TimeTree com IC, mediana, e — recalculados a partir de `studies.time_estimates` — n, quartis e amplitude. `timetree_pares.py`.
- Para comparação, idade do mesmo nó (MRCA de dois exemplares) na árvore datada do B10K (MCMCtree, HPD 95 %) e na Aves 1.6 (MCC de 100 árvores datadas; HPD 95 %). `idades_arvores_datadas.py`.

---

## 3. Resultados

### 3.1 B10K — frequências de quarteto nos ramos contestados

Linhas "ASTRAL-IV (63K)": árvores colapsadas em aLRT < 0,95, pontuando a topologia 63K (os mesmos ramos em todos os tipos de dado). Linhas "wASTRAL": árvores não colapsadas ponderadas pelo aLRT, busca livre; a resolução e a localPP são as da árvore que o wASTRAL achou. "p 1ª×2ª": qui-quadrado entre as duas resoluções mais frequentes; "p simetria": entre as duas menos frequentes. Não há linha do intergênico com 363 táxons (memória; ver Limitações): os valores do artigo para a base de Afroaves estão no texto abaixo.

**Base de Afroaves** — q de A+S / S+C / A+C

| Dado | Táxons, estimador | EN | q A+S | q S+C | q A+C | Resolução da árvore (localPP) | p 1ª×2ª | p simetria |
|---|---|---|---|---|---|---|---|---|
| intergênico 63K | 42 táxons, ASTRAL-IV (63K) | 61 087 | 0,374 | 0,333 | 0,293 | A+S (1,00) | 6.5e-34 | 3.7e-35 |
| intergênico 63K | 42 táxons, wASTRAL (aLRT), livre | 59 495 | 0,343 | 0,348 | 0,309 | S+C (1,00) | 0,140 | 9.8e-25 |
| éxons C12 | 42 táxons, ASTRAL-IV (63K) | 13 634 | 0,419 | 0,304 | 0,277 | A+S (1,00) | 1.9e-55 | 2.1e-5 |
| íntrons | 42 táxons, ASTRAL-IV (63K) | 37 125 | 0,377 | 0,336 | 0,287 | A+S (1,00) | 2.5e-20 | 1.2e-33 |
| íntrons | 42 táxons, wASTRAL (aLRT), livre | 31 892 | 0,345 | 0,349 | 0,306 | S+C (1,00) | 0,394 | 2.2e-18 |
| UCE | 42 táxons, ASTRAL-IV (63K) | 4 822 | 0,403 | 0,316 | 0,281 | A+S (1,00) | 1.1e-12 | 0,002 |
| UCE | 42 táxons, wASTRAL (aLRT), livre | 4 688 | 0,360 | 0,336 | 0,304 | A+S (1,00) | 0,042 | 0,007 |
| UCE | 363 táxons, ASTRAL-IV (63K) | 4 835 | 0,424 | 0,303 | 0,273 | A+S (1,00) | 2.8e-23 | 0,006 |


**Monofilia de Afroaves (na topologia 63K)** — q de (A+S)+C / C+U / (A+S)+U

| Dado | Táxons, estimador | EN | q (A+S)+C | q C+U | q (A+S)+U | Resolução da árvore (localPP) | p 1ª×2ª | p simetria |
|---|---|---|---|---|---|---|---|---|
| intergênico 63K | 42 táxons, ASTRAL-IV (63K) | 61 738 | 0,360 | 0,333 | 0,307 | (A+S)+C (1,00) | 2.4e-16 | 6.1e-16 |
| éxons C12 | 42 táxons, ASTRAL-IV (63K) | 13 868 | 0,319 | 0,344 | 0,337 | (A+S)+C (0,00) | 0,312 | 0,009 |
| íntrons | 42 táxons, ASTRAL-IV (63K) | 38 581 | 0,370 | 0,327 | 0,303 | (A+S)+C (1,00) | 2.5e-24 | 3.4e-9 |
| UCE | 42 táxons, ASTRAL-IV (63K) | 4 864 | 0,354 | 0,347 | 0,300 | (A+S)+C (0,93) | 0,555 | 4.4e-5 |
| UCE | 42 táxons, wASTRAL (aLRT), livre | 4 695 | 0,358 | 0,338 | 0,303 | (A+S)+C (1,00) | 0,097 | 0,003 |
| UCE | 363 táxons, ASTRAL-IV (63K) | 4 881 | 0,342 | 0,394 | 0,265 | (A+S)+C (0,00) | 2.4e-5 | 4.2e-12 |


**Base de Australaves** — q de F+PsP / K+F / K+PsP

| Dado | Táxons, estimador | EN | q F+PsP | q K+F | q K+PsP | Resolução da árvore (localPP) | p 1ª×2ª | p simetria |
|---|---|---|---|---|---|---|---|---|
| intergênico 63K | 42 táxons, ASTRAL-IV (63K) | 55 219 | 0,381 | 0,351 | 0,269 | F+PsP (1,00) | 4.3e-16 | 4.6e-133 |
| intergênico 63K | 42 táxons, wASTRAL (aLRT), livre | 58 243 | 0,406 | 0,306 | 0,288 | F+PsP (1,00) | 5.0e-177 | 4.9e-9 |
| éxons C12 | 42 táxons, ASTRAL-IV (63K) | 12 843 | 0,316 | 0,420 | 0,264 | F+PsP (0,00) | 1.8e-43 | 2.2e-14 |
| íntrons | 42 táxons, ASTRAL-IV (63K) | 31 435 | 0,428 | 0,325 | 0,248 | F+PsP (1,00) | 1.0e-98 | 6.5e-73 |
| íntrons | 42 táxons, wASTRAL (aLRT), livre | 30 396 | 0,441 | 0,286 | 0,273 | F+PsP (1,00) | 1.7e-219 | 0,002 |
| UCE | 42 táxons, ASTRAL-IV (63K) | 4 535 | 0,362 | 0,395 | 0,244 | F+PsP (0,00) | 0,011 | 1.5e-24 |
| UCE | 42 táxons, wASTRAL (aLRT), livre | 4 642 | 0,404 | 0,323 | 0,274 | F+PsP (1,00) | 9.7e-11 | 1.5e-5 |
| UCE | 363 táxons, ASTRAL-IV (63K) | 4 573 | 0,355 | 0,409 | 0,236 | F+PsP (0,00) | 3.0e-5 | 1.2e-25 |


**Ramos resolvidos** — q da resolução de referência (a mais frequente em todas as análises) e gCF da topologia 63K podada:

| Ramo | Intergênico (ASTRAL-IV / wASTRAL) | Íntrons | UCE (42 / 363 táxons) | Éxons C12 | gCF % (intergênico / íntrons / UCE / éxons) |
|---|---|---|---|---|---|
| Psittacopasserae (× F+Psi × F+Pas) | 0,375 / 0,387 | 0,405 | 0,404 / 0,400 | 0,384 | 16,2 / 21,8 / 19,6 / 7,9 |
| Cathartidae irmã de Accipitriformes s.s. | 0,949 / — | 0,892 | 0,969 / 0,968 | 0,556 | 91,9 / 82,1 / 95,1 / 29,9 |
| Pandionidae + Accipitridae | 0,995 / — | 0,980 | 0,991 / 0,990 | 0,673 | 98,3 / 94,8 / 98,3 / 46,7 |
| Strigidae, *Tyto* fora | 0,999 / — | 0,999 | 0,997 / 0,997 | 0,853 | 99,7 / 99,8 / 99,5 / 69,6 |


**Leitura dos ramos resolvidos** (faixa entre intergênico, íntrons e UCE com 42 táxons; os éxons C12, sem a 3ª posição do códon e com poucos sítios variáveis por loco, dão sempre o valor mais baixo):
- Cathartidae irmã de Accipitriformes s.s. (Sagittariidae + (Pandionidae + Accipitridae)): q = 0,892–0,969; éxons 0,556. Duas alternativas simétricas (p de simetria 0,03–0,63). Concordância alta, sem sinal de conflito estruturado.
- Pandionidae + Accipitridae, com Sagittariidae fora: q = 0,980–0,995; éxons 0,673.
- Strigidae monofilético, *Tyto* fora: q = 0,997–0,999; éxons 0,853.
- Psittacopasserae: q = 0,375–0,405 contra 0,307–0,348 das alternativas; p 1ª×2ª ≤ 10⁻⁵ em todos os tipos. Resolvido, mas com muita discordância (ramo curto em unidades coalescentes) e com assimetria (F+Psi > F+Pas) em todos os tipos.

**Base de Australaves** (Eufalconimorphae, F+PsP, contra Cariamiformes + Falconiformes, K+F):
- Dados não codificantes favorecem F+PsP: intergênico q = 0,381 × 0,351 (K+F) × 0,269 (K+PsP), p 1ª×2ª = 4×10⁻¹⁶; íntrons 0,428 × 0,325 × 0,248 (p = 10⁻⁹⁸); wASTRAL no intergênico 0,406 × 0,306 e nos íntrons 0,441 × 0,286.
- Éxons e UCEs, sem peso, favorecem **K+F**: éxons 0,420 × 0,316 (F+PsP), p = 2×10⁻⁴³; UCE 0,395 × 0,362 (p = 0,011) com 42 táxons e 0,409 × 0,355 (p = 3×10⁻⁵) com os 363 táxons. Nos éxons o escore bruto de quartetos também prefere K+F (T2 − T1 = +235 066 quartetos); nas UCEs o escore bruto prefere F+PsP (+249 472) e o wASTRAL ponderado também (0,404 × 0,323).
- Em todos os tipos de dado e estimadores as duas resoluções minoritárias são desiguais, sempre com K+F > K+PsP (p de simetria de 0,002 a 10⁻¹³³). ILS sozinha daria K+F ≈ K+PsP. Há um sinal sistemático que junta falcões e seriemas, mais forte em dado codificante e conservado. É o mesmo padrão que Kuhl et al. (2021) descrevem no texto: nas árvores de sequência codificante (CODON, CODON12 e aminoácidos) "Falconiformes e Cariamiformes foram afastados de Passeriformes e Psittaciformes para formar um agrupamento de aves de rapina", enquanto as 3′-UTR recuperam Australaves com Eufalconimorphae. As 11 árvores publicadas do B10K, inclusive a concatenada, têm F+PsP (exceto a de éxons, onde Australaves nem é monofilético: Cariamiformes cai junto de Hieraves).

**Base de Afroaves** (o ramo contestado):
- 63K, 363 táxons, artigo original (dados-fonte da Fig. 2a, nó 26): A+S 0,3495 × S+C 0,3456 × A+C 0,3049; os autores dão as duas primeiras como indistinguíveis (qui-quadrado, p ajustado por BH = 0,130) e a terceira como significativamente menor. Com subamostragem de 1 a 50 espécies por ordem (Fig. 2c, dados-fonte), a ordem se inverte: S+C 0,356–0,362 × A+S 0,336–0,341.
- Esta reanálise, 42 táxons, pontuando a topologia 63K: A+S vence em todos os tipos de dado — intergênico 0,374 × 0,333 (S+C) × 0,293 (A+C); íntrons 0,377 × 0,336 × 0,287; éxons 0,419 × 0,304 × 0,277; UCE 0,403 × 0,316 × 0,281 (363 táxons: 0,424 × 0,303 × 0,273). A ordem A+S > S+C > A+C é a mesma do artigo, mas com distância maior entre as duas primeiras.
- Mas o resultado muda com o estimador: nos íntrons, a busca livre do ASTRAL-IV devolve **S+C** (o escore bruto de quartetos prefere S+C mesmo com q normalizado preferindo A+S), e o wASTRAL dá empate (S+C 0,349 × A+S 0,345; p = 0,39). Nas UCEs, o wASTRAL ainda dá A+S, por pouco (0,360 × 0,336; p = 0,04). No intergênico, a árvore da busca livre do ASTRAL-IV com os 42 táxons é (((C,S),A),U): **S+C**, com localPP 0 para esse ramo (q de A+S 0,374 × S+C 0,333) — a mesma topologia da concatenada do artigo, e o escore final (3 178 321 971) praticamente igual ao da topologia 63K com S+C trocado (3 178 318 155). O wASTRAL ponderado no intergênico também sai com S+C, num empate: S+C 0,348 × A+S 0,343 × A+C 0,309 (p 1ª×2ª = 0,14), e localPP 1 para S+C mesmo assim — exemplo de localPP inflada pela assimetria das minoritárias (ver Limitações).
- Em todos os casos, A+C é a menos frequente e as duas minoritárias são desiguais (p de simetria de 0,007 a 10⁻³⁵): o padrão "duas resoluções quase empatadas e uma claramente menor" que Stiller et al. já apontavam como incompatível com ILS pura.

**Monofilia de Afroaves** (Accipitrimorphae dentro de Afroaves ou fora):
- Intergênico (42 táxons): Afroaves q = 0,360 × C+U 0,333 × (A+S)+U 0,307 (p = 2×10⁻¹⁶). Íntrons: 0,370 × 0,327 × 0,303.
- UCE: empate com 42 táxons (0,354 × 0,347; p = 0,56) e **C+U vence com 363 táxons** (0,394 × 0,342; p = 2×10⁻⁵) — o que reproduz a árvore UCE publicada pelo B10K, ((A,S),(C,U)), sem Afroaves. Éxons: C+U 0,344 × (A+S)+U 0,337 × Afroaves 0,319.
- Escore bruto de quartetos (42 táxons): nas UCEs, a topologia 63K (((A,S),C),U) vence a de Hieraves irmã de C+U por 825 665 quartetos e a de Accipitrimorphae irmã de Telluraves por 1 674 036; nos éxons, Hieraves irmã de C+U vence a 63K por 1 639 862. No intergênico e nos íntrons, o escore bruto inverte o que o q diz: com 42 táxons, a topologia com **S+C** vence a 63K (A+S) por 723 012 quartetos no intergênico e por 174 401 nos íntrons, enquanto o q normalizado dá A+S 0,374 × 0,333 e 0,377 × 0,336. Nos mesmos dois tipos de dado, a 63K vence com folga as topologias sem Afroaves (Accipitrimorphae irmã de Telluraves: −10,2 milhões no intergênico, −4,3 milhões nos íntrons; Hieraves irmã de C+U: −15,4 e −6,0 milhões) e a de Cariamiformes + Falconiformes (−4,4 e −2,4 milhões). Em termos relativos as diferenças são pequenas (o escore total do intergênico com 42 táxons é 3,18 × 10⁹), mas só mudam os quartetos em volta do ramo alterado.

**Por que q e escore bruto discordam (testado).** O q e a localPP que o ASTRAL imprime são médias por quarteto de espécies em volta do ramo (Sayyari & Mirarab 2016): cada combinação de quatro espécies pesa igual, dividida pelo número de árvores de genes que a resolvem. O objetivo da busca soma os quartetos resolvidos de todas as árvores. Quando as árvores de genes têm ramos colapsados (e dado faltante), quartetos resolvidos por poucas árvores pesam no q tanto quanto os resolvidos por muitas, e as duas medidas podem divergir. Teste com as 1 209 árvores UCE que têm os 42 táxons focais (sem dado faltante), base de Australaves (`teste_ponderacao.sh`): com as árvores binárias (sem colapsar), q e escore bruto concordam — F+PsP 0,391 × K+F 0,337, e a topologia 63K vence a de K+F por 72 348 quartetos; com as mesmas árvores colapsadas em aLRT < 0,95, o escore bruto continua preferindo F+PsP (+60 606), mas o q vira para K+F (0,406 × 0,358; localPP de F+PsP 9 × 10⁻⁶). A inversão vem do colapso, não de falta de táxons. Isto vale também para os números publicados (o ASTRAL do artigo usou as mesmas árvores colapsadas). Na base de Afroaves, a subamostragem da Fig. 2c do artigo (S+C > A+S com poucas espécies por ordem) e o escore bruto com 42 táxons aqui (S+C) apontam para o mesmo lado; o q normalizado com 42 táxons e a árvore principal com 363 táxons apontam para o outro.

### 3.2 B10K — o que as 11 árvores publicadas dizem dos mesmos ramos

| Árvore publicada (363 táxons) | Base de Afroaves | Topologia enraizada A, S, C, U | Base de Australaves | Rótulo do ramo A+S ou S+C |
|---|---|---|---|---|
| 63K ASTRAL (principal) | A+S | (((A,S),C),U) | F+PsP | localPP 1 |
| 63K GOBS (bootstrap de locos) | A+S | (((A,S),C),U) | F+PsP | 76 |
| 63K concatenada RAxML-NG | **S+C** | (((C,S),A),U) | F+PsP | BS 100 |
| Íntrons ASTRAL | A+S | (((A,S),C),U) | F+PsP | 0 (rótulo publicado) |
| Éxons C12 ASTRAL | A+S | C e U não monofiléticos; clados A+S e A+S+K | U quebrado | 0,55 |
| UCE ASTRAL | A+S | **((A,S),(C,U))** — Afroaves não monofilético | F+PsP | 1 |
| 94K, 128K, 153K, 159K ASTRAL | **S+C** | S aninhado em C; A irmã de (S+C)+U | F+PsP | 1 |

O teste de politomia do próprio artigo (ASTRAL 5.15.1, `-t 10`, 63K) dá p = 0 para todos os ramos aqui discutidos; os únicos ramos com p > 0 em toda a árvore são a posição de Opisthocomiformes (p = 0,057) e um ramo de Elementaves (p = 0,00018). Ou seja: o ASTRAL rejeita a politomia "dura" na base de Afroaves com 63K genes, mas isso não decide entre as duas resoluções mais frequentes (ver 3.1).

### 3.3 Concatenada × coalescente e o que a discordância indica (Regra 2)

- **Onde concatenada e espécies discordam.** Nos mesmos 63 430 locos, a árvore de espécies (ASTRAL-MP) tem A+S com localPP 1 (bootstrap de locos GOBS 76) e a concatenada (RAxML-NG) tem S+C com bootstrap 100, isto é, (((C,S),A),U), a topologia de Hackett 2008 e Jarvis 2014. Em todos os outros ramos desta nota as duas concordam. Discordância concentrada num ramo só, com q quase empatadas, é o caso em que a concatenação pode convergir com suporte alto para qualquer um dos lados (Kubatko & Degnan 2007; Mendes & Hahn 2018) e em que o ASTRAL depende de cada árvore de gene estar bem estimada. Bootstrap 100 ali não é concordância.
- **Assimetria das alternativas.** Na base de Afroaves, S+C (a segunda) é sempre muito mais frequente que A+C (a terceira): 0,333 × 0,293 no intergênico, 0,336 × 0,287 nos íntrons. Sob ILS pura, as duas minoritárias têm a mesma frequência esperada. O excesso de S+C cabe em introgressão antiga entre linhagens de Strigiformes e de Coraciimorphae (ou entre Accipitrimorphae e o resto de Afroaves), em erro sistemático das árvores de genes (atração de ramos longos, composição), ou em paralogia oculta. Stiller et al. testaram hibridização com CoalHMM sem detectá-la, mas observam que o método supõe ILS e só pega sinal forte; Şapcı et al. (2026), com o Phlag sobre as mesmas árvores de genes, acharam blocos genômicos cujas árvores de genes se desviam da distribuição esperada justamente neste ramo, e o wASTRAL só desses blocos move Accipitriformes para junto de Strigiformes — "consistente com uma ou mais reticulações na base das aves terrestres". Na base de Australaves a assimetria é ainda maior (K+F ≫ K+PsP) e segue o tipo de dado (forte em éxons e UCEs), o que pesa mais para viés sistemático de dado codificante/conservado do que para reticulação.
- **gCF.** Na topologia 63K podada, o ramo A+S tem gCF = 2,5 % no intergênico (gDF1 1,0 %, gDF2 0,6 %, gDFP 95,9 %; 62 263 árvores decisivas); nos outros tipos, 6,6 % (íntrons), 3,9 % (UCE), 1,6 % (éxons), com gDFP (árvores decisivas mas discordantes por polifilia) de 87 a 98 %. Ou seja, quase nenhuma árvore de gene de 1 kb recupera os quatro grupos inteiros em volta do ramo: o sinal por loco nessa profundidade é fraco e o erro de estimativa das árvores de genes é grande. É o cenário em que gCF ≪ sCF (sCF não calculado aqui: sem alinhamentos) aponta árvores de genes mal estimadas, não conflito biológico (Minh, Hahn & Lanfear 2020; Lanfear & Hahn 2024). Para o ramo F+PsP: gCF 15,9 % (íntrons), 9,9 % (UCE), 2,4 % (éxons), 8,8 % no intergênico, com as duas discordantes desiguais (gDF1 2,8 % × gDF2 7,4 %). Nos ramos resolvidos, o gCF sobe junto com o q: Pandionidae + Accipitridae 94,8 % (íntrons) e 98,3 % (UCE); Strigidae 99,8 % e 99,5 %.
- **Leitura.** A base de Afroaves está no limite de resolução destes dados: o resultado troca de lado com o estimador (q normalizado × escore bruto × ponderação por suporte), com o tipo de dado, com a amostragem de táxons (Fig. 2c do artigo; 363 × 42 táxons aqui) e com o número de locos (63K dá A+S; 94K a 159K, S+C). Nenhuma das três resoluções tem sustentação que sobreviva a essas quatro variações.

### 3.4 Open Tree 16.1 e Aves 1.6: quem sustenta e quem conflita

Contagens só de árvores filogenéticas (as taxonômicas OTT, eBird e Clements saem da conta). "Ausente" = o clado não existe na síntese; o MRCA dos táxons que o definiriam é outro nó, e as contagens daquele nó não se aplicam.

| Nó | Open Tree 16.1: sustentam / conflitam | Aves 1.6: sustentam / conflitam | Árvore-fonte que põe o nó na Aves 1.6 (prioridade) |
|---|---|---|---|
| Telluraves | 18 / 0 | 13 / 0 | Stiller 2024 (19) |
| Afroaves (A+S+C) | 13 / 5 | 8 / 5 | Stiller 2024 (19) |
| A+S (Hieraves) | 4 / 14 | 4 / 9 | Stiller 2024 (19) |
| S+C | ausente | ausente | — |
| A+C | ausente | ausente | — |
| Coraciimorphae | 11 / 4 | 10 / 2 | Stiller 2024 (19) |
| Accipitrimorphae | 16 / 0 | 12 / 0 | Stiller 2024 (19) |
| Accipitriformes s.s. (Sag+Pan+Acd) | 12 / 0 | 10 / 0 | Stiller 2024 (19) |
| Pandionidae + Accipitridae | 12 / 0 | 10 / 0 | Stiller 2024 (19) |
| Sagittariidae + Accipitridae | ausente | ausente | — |
| Australaves | 15 / 0 | 12 / 0 | Stiller 2024 (19) |
| Eufalconimorphae (F+PsP) | 15 / 0 | 12 / 0 | Stiller 2024 (19) |
| Falconiformes + Cariamiformes | ausente | ausente | — |
| Psittacopasserae | 17 / 0 | 13 / 0 | Stiller 2024 (19) |
| Falconidae | 14 / 0 | só caminho parcial em 13 árvores | — |
| Herpetotheriinae (*Herpetotheres*+*Micrastur*) | 12 / 0 | 9 / 0 | Zhao 2025 (78) |
| Falconinae + Polyborinae | 14 / 0 | 11 / 0 | Zhao 2025 (78) |
| Strigiformes | 15 / 0 | 12 / 0 | Stiller 2024 (19) |
| Strigidae | 15 / 0 | 13 / 0 | Stiller 2024 (19) |
| Tytonidae | 13 / 0 | 11 / 0 | Salter 2020 (62) |

Sustentam A+S na Open Tree 16.1: Wu et al. 2024, Stiller et al. 2024, Claramunt et al. 2025, Zhao et al. 2025. Conflitam: Hackett 2008, Jetz 2012, D'Emic & Benson 2013, Jarvis 2014, Price 2014, Suh 2015, Prum 2015, Burleigh 2015, Claramunt & Cracraft 2015, Ksepka & Phillips 2015, Reddy 2017, Kimball 2019 (2 árvores), Feng 2020.

Como ler: a síntese da Open Tree é um "veto reverso" por prioridade — o nó entra pela árvore de maior prioridade que o contém, e árvore de prioridade menor que conflita não o remove (McTavish et al. 2025; Redelings & Holder 2017). Na Aves 1.6 a árvore de maior prioridade entre as de família para cima é a de Stiller et al. 2024 (posição 19 de 338; as 18 acima são árvores de grupos menores e da taxonomia eBird). Por isso a síntese tem A+S com 4 árvores a favor e 14 contra: o nó vem do B10K, não da maioria. Contagem de árvores na Open Tree não é peso de evidência — as árvores não são independentes (Jetz 2012 usa o esqueleto de Hackett 2008; Kimball 2019 é superárvore; Burleigh 2015 é supermatriz de dados do GenBank que incluem os de Hackett; Claramunt 2025 é árvore datada de 9 071 espécies; Feng 2020 e Stiller 2024 usam os mesmos 363 genomas).

### 3.5 Meta-análise das árvores-fonte (fração de quartetos e topologia enraizada)

Topologia enraizada dos quatro grupos A, S, C, U em cada árvore-fonte informativa (fração de 2 000 amostras de 5 táxons; em todas a topologia majoritária tem ≥ 0,82):

| Topologia (A, S, C, U) | Nome | Árvores-fonte (fração) | Primárias* |
|---|---|---|---|
| (((C,S),A),U) | Afroaves; A na base de Afroaves; S irmã de C | Hackett 2008 (0,96), Jetz 2012 (0,96), Jarvis 2014 (1,00), Price 2014 (0,96), Suh 2015 (0,82), Burleigh 2015 (0,97), Claramunt & Cracraft 2015 (0,96), Ksepka & Phillips 2015 (1,00), Reddy 2017 (0,94), Kimball 2019 (0,96; 0,97) | 7 |
| (((A,S),C),U) | Afroaves com Hieraves | Stiller 2024 (0,98), Zhao 2025 (0,97) | 2 |
| (((C,S),U),A) | A irmã do resto de Telluraves | Prum 2015 (0,95), Feng 2020 (0,97) | 2 |
| ((A,S),(C,U)) | Hieraves irmã de C+U | Wu 2024 (0,97), Claramunt 2025 (0,98) | 1 |
| (((A,C),S),U) | — | D'Emic & Benson 2013 (1,00) | 0 |

\* Primárias = inferência a partir de dados próprios; ficam fora Jetz 2012 (esqueleto de Hackett), Kimball 2019 (superárvore), Hugall & Stuart-Fox 2012 e D'Emic & Benson 2013 (árvores compiladas para estudo comparativo), Ksepka & Phillips 2015 (estudo de datação) e Claramunt et al. 2025 (árvore datada de 9 071 espécies). Os mapeamentos de tipo estão em `tabela_sintese.py → DERIVADAS`.

Fração de quartetos por pergunta (resolução majoritária por árvore; todas / só primárias):
- Base de Afroaves: S+C 14 / 10; A+S 4 / 3; A+C 3 / 1. Os 3 votos A+C vêm de Hugall 2012 e D'Emic 2013 (compiladas) e de Fuchs et al. 2011 (23 folhas, estudo de *Micrastur*, 0,67) — não são evidência sobre a base de Afroaves.
- Monofilia de Afroaves: monofilético 14 / 9; A (ou Hieraves) fora de Afroaves 4 / 3 (Prum 2015, Feng 2020, Wu 2024, Claramunt 2025).
- Base de Australaves: F+PsP em 15 / 11 árvores; nenhuma árvore com K+F ou K+PsP majoritário (frações de K+F ≤ 0,09, só nas superárvores/supermatrizes grandes).
- Psittacopasserae: 18 / 12, nenhuma contra.
- Cathartidae irmã de Accipitriformes s.s.: 12 / 8; Pandionidae + Accipitridae: 12 / 8; nenhuma contra.
- Falconidae: Herpetotheriinae irmã de Polyborinae+Falconinae em 13 / 8 árvores, fração 1,00 em todas; nenhuma contra. Fuchs et al. 2015 não entra porque, na árvore depositada, não há grupo externo a Falconidae.
- Strigiformes: Strigidae monofilético em relação a Tytonidae em 14 / 9; Hugall & Stuart-Fox 2012 (compilada) dá 0,61.

### 3.6 Idades

| Par (nó) | TimeTree: mediana (Q1–Q3; amplitude; n) | TimeTree: idade estimada [IC] | B10K MCMCtree [HPD 95 %] | Aves 1.6 MCC [HPD 95 %] |
|---|---|---|---|---|
| Accipitridae–Cathartidae | 66,0 (60,1–89,5; 38,8–109,3; 13) | 71,2 [60,8–81,5] | 62,2 [61,3–63,0] | 55,5 [44,4–62,8] |
| Accipitridae–Sagittariidae | 60,3 (40,2–76,0; 32,5–97,3; 11) | 59,0 [47,2–70,7] | 54,1 [33,8–60,4] | 50,4 [40,2–62,2] |
| Accipitridae–Pandionidae | 48,9 (30,6–58,3; 19,2–83,3; 17) | 45,8 [37,0–54,6] | 32,3 [18,3–49,1] | 41,6 [26,9–50,9] |
| Falconidae–Psittaciformes | 71,2 (63,0–88,9; 58,8–116,0; 23) | 74,2 [68,8–79,7] | 62,6 [61,7–63,8] | 59,8 [56,5–63,6] |
| Falconidae–Cariamidae | 63,2 (61,2–86,9; 49,0–133,9; 9) | 67,3 [58,6–76,1] | 63,1 [62,3–64,3] | 62,3 [58,5–65,6] |
| Tytonidae–Strigidae | 61,4 (51,5–77,9; 39,0–96,2; 15) | 65,7 [55,7–75,7] | 55,6 [52,2–58,9] | 51,8 [43,9–61,7] |
| Strigiformes–Accipitriformes | 65,8 (62,9–86,0; 57,2–110,9; 8) | 66,2 [62,0–70,4] | 63,4 [62,7–63,9] | 58,5 [48,7–63,9] |
| *Herpetotheres*–*Falco* | 38,3 (31,0–50,8; 21,7–65,4; 13) | 39,5 [33,0–46,0] | 32,0 [25,5–39,6] | 38,3 [30,4–46,0] |
| *Caracara*–*Falco* | 22,7 (TimeTree; recalculada 22,0 com um registro de 0 Ma de Ksepka & Phillips 2015; 8) | 25,6 [20,5–30,7] | ausente do B10K | 27,6 [18,6–34,7] |

Idades em Ma. Exemplares nas árvores datadas: *Aquila chrysaetos*, *Cathartes aura*, *Sagittarius*, *Pandion*, *Falco peregrinus*, *Nestor notabilis*, *Cariama cristata*, *Tyto alba*, *Strix occidentalis*, *Herpetotheres cachinnans*, *Caracara plancus*.

Três observações que valem para o texto: (1) no TimeTree, Falconidae–Psittaciformes (mediana 71,2) sai **mais antigo** que Falconidae–Cariamidae (63,2), o que contradiz Eufalconimorphae (F mais próximo de PsP que de K): as medianas do TimeTree vêm de conjuntos diferentes de estudos para cada par e não precisam obedecer a uma topologia; no B10K, com uma topologia só, a ordem é a esperada (62,6 < 63,1). (2) Para Strigiformes–Accipitriformes o "nó NCBI" que o TimeTree usa é Telluraves (a taxonomia NCBI não tem Afroaves nem Hieraves), então a mediana mistura estimativas do nó Hieraves, de Afroaves e de Telluraves conforme a topologia de cada estudo. (3) A amplitude no TimeTree é larga (p.ex. 38,8–109,3 Ma para Accipitridae–Cathartidae) e inclui estimativas pré-2010 com calibrações e dados muito diferentes.

---

## 4. Tabela-síntese: nós-chave × fontes

Siglas como no topo. q = frequência de quarteto da resolução indicada (ASTRAL-IV, árvores colapsadas em aLRT < 0,95, pontuando a topologia 63K; 42 táxons salvo indicação). "OT" = Open Tree 16.1; "Aves" = Aves 1.6; contagens só de árvores filogenéticas, sustentam/conflitam. "Fontes" = árvores-fonte da Open Tree com maioria de quartetos para cada resolução (todas / só inferência primária). Idade = TimeTree mediana (n) | B10K MCMCtree [HPD 95 %].

| Nó-chave | B10K quartetos por tipo de dado (q) | B10K: árvores publicadas (11) | OT / Aves | Árvores-fonte (todas / primárias) | Idade (Ma) |
|---|---|---|---|---|---|
| **Base de Afroaves**: A+S × S+C × A+C | intergênico 0,374 × 0,333 × 0,293; íntrons 0,377 × 0,336 × 0,287; éxons 0,419 × 0,304 × 0,277; UCE 0,403 × 0,316 × 0,281; artigo (363 táx.) 0,3495 × 0,3456 × 0,3049. Estimador muda o lado: a busca livre do ASTRAL-IV dá S+C no intergênico e nos íntrons; wASTRAL empata — intergênico A+S 0,343 × S+C 0,348 (p = 0,14), íntrons 0,345 × 0,349 (p = 0,39); escore bruto (42 táx.) prefere S+C no intergênico (+723 012) e nos íntrons (+174 401) | A+S: 63K ASTRAL (localPP 1; GOBS 76), íntrons, éxons (0,55), UCE. S+C: 63K concatenada (BS 100), 94K, 128K, 153K, 159K | A+S: OT 4/14; Aves 4/9 (entra por Stiller 2024). S+C e A+C: ausentes | S+C 14/10; A+S 4/3; A+C 3/1 (compiladas e um estudo de *Micrastur*) | Strigiformes–Accipitriformes 65,8 (8) \| 63,4 [62,7–63,9] |
| **Monofilia de Afroaves** (A dentro × A ou A+S fora) | intergênico Afroaves 0,360 × C+U 0,333 × (A+S)+U 0,307; íntrons 0,370 × 0,327; UCE 0,354 × 0,347 (p = 0,56) e, com 363 táx., C+U 0,394 × 0,342; éxons C+U 0,344 × Afroaves 0,319 | Afroaves: 63K, concatenada, íntrons. Fora: UCE ((A,S),(C,U)); 94K–159K (A irmã do resto de Telluraves) | Afroaves: OT 13/5; Aves 8/5 | monofilético 14/9; fora 4/3 (Prum 2015, Feng 2020, Wu 2024, Claramunt 2025) | — |
| **Topologia enraizada A, S, C, U** | busca livre, 42 táxons: UCE (((A,S),C),U); íntrons (((C,S),A),U); éxons: C e U não monofiléticos; intergênico (((C,S),A),U) (ASTRAL-IV e wASTRAL) | (((A,S),C),U): 63K, GOBS, íntrons; (((C,S),A),U): concatenada; ((A,S),(C,U)): UCE; S dentro de C, A fora: 94K–159K | — | (((C,S),A),U) 11/7; (((A,S),C),U) 2/2; (((C,S),U),A) 2/2; ((A,S),(C,U)) 2/1 | — |
| **Base de Australaves**: F+PsP × K+F × K+PsP | intergênico 0,381 × 0,351 × 0,269; íntrons 0,428 × 0,325 × 0,248; **éxons 0,316 × 0,420** × 0,264; **UCE 0,362 × 0,395** × 0,244 (363 táx.: 0,355 × 0,409); wASTRAL intergênico 0,406 × 0,306 e UCE 0,404 × 0,323; K+F > K+PsP em todos (p de simetria ≤ 0,002); escore bruto: F+PsP no intergênico (+4,4 mi), íntrons (+2,4 mi) e UCE (+0,25 mi); K+F nos éxons (+0,24 mi) | F+PsP em todas, exceto éxons (U não monofilético; K junto de A+S) | F+PsP: OT 15/0; Aves 12/0. K+F: ausente | F+PsP 15/11; nenhuma contra | Falconidae–Cariamidae 63,2 (9) \| 63,1 [62,3–64,3] |
| **Psittacopasserae** (Psi+Pas × F+Psi × F+Pas) | intergênico 0,375 × 0,319 × 0,307; íntrons 0,405; éxons 0,384; UCE 0,404 | todas | OT 17/0; Aves 13/0 | 18/12; nenhuma contra | Falconidae–Psittaciformes 71,2 (23) \| 62,6 [61,7–63,8] |
| **Cathartidae irmã de Accipitriformes s.s.** | intergênico 0,949; íntrons 0,892; UCE 0,969; éxons 0,556 | todas | OT 12/0; Aves 10/0 | 12/8; nenhuma contra | Accipitridae–Cathartidae 66,0 (13) \| 62,2 [61,3–63,0] |
| **Pandionidae + Accipitridae** (Sagittariidae fora) | intergênico 0,995; íntrons 0,980; UCE 0,991; éxons 0,673 | todas | OT 12/0; Aves 10/0 | 12/8; nenhuma contra | Accipitridae–Pandionidae 48,9 (17) \| 32,3 [18,3–49,1]; Accipitridae–Sagittariidae 60,3 (11) \| 54,1 [33,8–60,4] |
| **Strigidae monofilético, Tytonidae fora** | intergênico 0,999; íntrons 0,999; UCE 0,997; éxons 0,853 | todas | Strigidae OT 15/0, Aves 13/0; Tytonidae OT 13/0, Aves 11/0 | 14/9; nenhuma contra | Tytonidae–Strigidae 61,4 (15) \| 55,6 [52,2–58,9] |
| **Falconidae: Herpetotheriinae irmã de Polyborinae + Falconinae** | não testável (B10K só tem *Herpetotheres* e *Falco*) | — | Herpetotheriinae OT 12/0, Aves 9/0; Falconinae+Polyborinae OT 14/0, Aves 11/0 | 13/8, fração 1,00 em todas | *Herpetotheres*–*Falco* 38,3 (13) \| 32,0 [25,5–39,6]; *Caracara*–*Falco* 22,7 (8) \| — |

Versão completa (107 colunas, todas as análises e todos os testes): `runs/rapinantes/reanalise/tabela_sintese.tsv`.

---

## 5. Métodos (parágrafo para o trabalho)

Reanalisamos as árvores de genes publicadas por Stiller et al. (2024) para 363 espécies — 63 430 locos intergênicos de 1 kb, 14 972 éxons (1ª e 2ª posições do códon), 44 846 íntrons e 4 985 UCEs —, inferidas pelos autores por máxima verossimilhança (ParGenes/RAxML-NG) com suporte por aLRT paramétrico, e obtidas do repositório oficial do estudo (doi:10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a). Como no artigo, ramos com aLRT < 0,95 foram colapsados; as árvores foram então podadas para 42 táxons (os 15 rapinantes amostrados, representantes de todas as ordens de Telluraves e seis grupos externos), descartando as que ficaram com menos de quatro táxons. As frequências de quarteto, a probabilidade local e o número efetivo de genes de cada ramo foram obtidos com ASTRAL-IV v1.25.4.8 (pacote ASTER v1.25; Zhang, Nielsen & Mirarab 2025), pontuando a topologia publicada (`-C -u 3`) e em busca livre, e com o wASTRAL v1.25.3.8 híbrido ponderado pelo aLRT (Zhang & Mirarab 2022); para as UCEs, a topologia publicada foi também pontuada com os 363 táxons. Em cada ramo testamos a hipótese de politomia (Sayyari & Mirarab 2018), a igualdade das duas resoluções mais frequentes e a simetria das duas menos frequentes (qui-quadrado sobre as contagens efetivas), e comparamos o escore total de quartetos de oito topologias alternativas fixas; para separar o efeito do colapso dos ramos, repetimos a pontuação com as 1 209 árvores UCE que contêm os 42 táxons, colapsadas e não colapsadas. Fatores de concordância de genes (gCF; Minh, Hahn & Lanfear 2020) foram calculados no IQ-TREE 3.1.4. Para a síntese, consultamos a API v3 da Open Tree of Life (síntese opentree16.1, taxonomia OTT 3.7draft3) e a síntese Aves 1.6 (McTavish et al. 2025), listando as árvores-fonte que sustentam e que conflitam com cada nó, e classificamos cada uma das 28 árvores-fonte ligadas a esses nós pela fração de 2 000 quartetos amostrados que sustenta cada resolução e pela topologia enraizada de Accipitrimorphae, Strigiformes, Coraciimorphae e Australaves. Idades de divergência vieram da API do TimeTree 5 (Kumar et al. 2022; versão de 2026-08-20) e das árvores datadas do B10K e da Aves 1.6. Todas as etapas estocásticas usaram semente fixa (12345); scripts e saídas estão em `scripts/rapinantes/reanalise` e `runs/rapinantes/reanalise`.

---

## 6. Limitações

- **Não houve inferência de árvore de gene nesta etapa.** Os quartetos vêm das árvores de genes do B10K tal como publicadas (RAxML-NG + aLRT paramétrico). Erro de estimativa delas entra inteiro aqui. Os alinhamentos não foram baixados (instrução: só árvores), então não há sCF/sCFL, não há reestimativa de árvore de gene, nem teste de saturação, composição ou modelo. O gCF calculado é baixíssimo nos ramos profundos (ver 3.3) porque poucas árvores de genes recuperam os grupos inteiros; isso limita o que o gCF diz ali.
- **Reanálise por pipeline (`filogenia.py`) não rodou.** O caminho alternativo (Prum 2015, Jarvis 2014, Kuhl 2021) só valia se o B10K fosse inacessível, e ele foi acessível. Prum, Jarvis e outras entram pela Open Tree, como árvores publicadas, não como dados reanalisados; Kuhl et al. 2021 não está na Open Tree e entra só pelo que o texto do artigo diz.
- **Conjunto focal pequeno por grupo.** 42 táxons, com 1–5 por família de rapinante e 3 Passeriformes. Stiller et al. mostraram (Fig. 2c, dados-fonte lidos aqui) que a base de Afroaves muda com a amostragem de táxons. Por isso a pontuação com 363 táxons foi feita para as UCEs, e para o intergênico usei os números com 363 táxons do próprio artigo.
- **q normalizado × escore bruto.** A frequência q e a localPP do ASTRAL são médias por quarteto de espécies; o objetivo da busca é a soma bruta de quartetos. Com árvores de genes colapsadas, as duas podem indicar resoluções diferentes no mesmo ramo — aconteceu na base de Afroaves (intergênico e íntrons) e na de Australaves (UCEs), e o teste com árvores completas mostrou que a inversão vem do colapso. Colapsar em aLRT < 0,95 é o procedimento do artigo; não refiz a análise com outro limiar. Reportei as duas medidas; nenhuma é "a certa".
- **localPP sob violação do modelo.** A localPP supõe que as duas resoluções minoritárias têm a mesma frequência. Com assimetria forte (p de simetria até 10⁻³⁵ na base de Afroaves) e EN de dezenas de milhares, a localPP vai a 0 ou 1 com diferenças de q de 0,004. LocalPP = 1 aqui não é resolução.
- **Os p-valores são nominais**, sem correção para múltiplos testes, e tratam as contagens efetivas como independentes. Só os extremos (p < 10⁻⁵ ou p > 0,05) foram interpretados.
- **Open Tree e Aves 1.6 não são evidência independente.** O nó entra pela árvore de maior prioridade (veto reverso), e as árvores-fonte compartilham dados (ver 3.4). A contagem de árvores a favor e contra é descrição da literatura depositada, não peso de evidência. Entre as 28 árvores ligadas a estes nós não estão, por exemplo, Kuhl et al. 2021 (3′-UTR) nem as reanálises de 2025–2026 das árvores do B10K.
- **Classificação das árvores-fonte por taxonomia OTT.** Folhas sem ott_id ficam fora; táxons fósseis sem subgrupo ficam fora. Ruído residual de 2–6 % de quartetos "A+U" aparece em várias árvores e não muda nenhuma maioria.
- **B10K não amostra Polyborinae, *Micrastur* nem *Phodilus*.** Base de Falconidae e Tytonidae/Strigidae com mais de um Tytonidae vêm só da Open Tree e do TimeTree.
- **Idades.** A mediana do TimeTree mistura estudos com topologias, calibrações e dados diferentes, e o nó NCBI usado para Strigiformes–Accipitriformes é Telluraves. Um registro de 0 Ma (Ksepka & Phillips 2015, *Caracara*–*Falco*) é erro do banco; a mediana publicada pelo TimeTree (22,65) parece já o excluir, a recalculada aqui (22,04) o inclui.
- **Memória e quedas do WSL.** A pontuação da topologia 63K com 363 táxons e as 63 430 árvores intergênicas (`-u 3`) passou de 14 GB e foi morta pelo OOM killer (journal do WSL, 2026-10-18 01:17); as duas tentativas anteriores coincidiram com quedas da VM do WSL (01:07 e 01:12), e houve outra às 00:52 com outras análises pesadas rodando em paralelo. Por isso não há números intergênicos com 363 táxons recalculados aqui (usei os dados-fonte do artigo), nem pontuação com 363 táxons para éxons e íntrons. O script ficou retomável (pula etapas com `freqQuad.csv` pronto) e as execuções seguintes rodaram como unidades transitórias do systemd com teto de memória de 9 GB.

---

## 7. O que estes dados não dizem

Nenhuma destas análises diz **quando** as linhagens divergiram: as idades da seção 3.6 são compiladas de estudos datados (TimeTree) ou lidas de duas árvores datadas publicadas, não estimadas aqui, e herdam as calibrações e as topologias daqueles estudos — inclusive a topologia de Afroaves que este documento mostra estar em aberto. Também não dizem **em que direção** mudou caráter algum: dizer que a predação diurna, a garra ou o bico em gancho é ancestral de Telluraves e foi perdido em Coraciimorphae e Psittacopasserae, ou que surgiu duas ou três vezes, exige reconstrução de estado ancestral sobre árvores alternativas — e o resultado depende justamente da base de Afroaves e de Australaves, que aqui não se resolve. Não dizem que Hieraves (Accipitrimorphae + Strigiformes) é um clado: na melhor análise disponível (63K, 363 táxons) as duas resoluções mais frequentes são estatisticamente indistinguíveis, a concatenação dos mesmos dados dá a outra, e conjuntos maiores do próprio B10K (94K a 159K locos) também dão a outra. Não dizem que a convergência ecológica dos rapinantes foi ou não independente em falcões e gaviões: isso depende de Australaves × Afroaves, que é sólido (Falconiformes com Psittacopasserae), mas a parte "quantas origens" depende também da base de Afroaves e de Telluraves. E, como os quartetos vêm de árvores de genes já estimadas, não separam ILS de erro de estimativa das árvores de genes, de introgressão antiga ou de viés de composição: a assimetria das resoluções minoritárias aponta que ILS pura não explica tudo, mas não aponta qual das outras causas.

---

## 8. Fontes (DOI/URL conferidos nesta sessão)

- Stiller J et al. 2024. Complexity of avian evolution revealed by family-level genomes. *Nature* 629:851–860. doi:10.1038/s41586-024-07323-1 (PMC11111414). Dados: doi:10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a → https://erda.ku.dk/archives/341f72708302f1d0c461ad616e783b86/published-archive.html; downloads por https://sid.erda.dk/share_redirect/ENhZODU9YE/…
- McTavish EJ et al. 2025. A complete and dynamic tree of birds. *PNAS* 122(18):e2409658122. doi:10.1073/pnas.2409658122. Dados: https://github.com/McTavishLab/AvesData (Tree_versions/Aves_1.6).
- Open Tree of Life, API v3: https://api.opentreeoflife.org/v3 (síntese opentree16.1, OTT 3.7draft3). Redelings BD, Holder MT 2017. *PeerJ* 5:e3058. doi:10.7717/peerj.3058.
- Kumar S et al. 2022. TimeTree 5. *Mol Biol Evol* 39:msac174. doi:10.1093/molbev/msac174. API: https://timetree.temple.edu/api (versão 20260820).
- Zhang C, Nielsen R, Mirarab S 2025. ASTER. *Mol Biol Evol* 42:msaf172. doi:10.1093/molbev/msaf172. Zhang C, Mirarab S 2022. wASTRAL. *Mol Biol Evol* 39:msac215. doi:10.1093/molbev/msac215.
- Sayyari E, Mirarab S 2016. localPP. *Mol Biol Evol* 33:1654. doi:10.1093/molbev/msw079. Sayyari E, Mirarab S 2018. Teste de politomia. *Genes* 9:132. doi:10.3390/genes9030132.
- Minh BQ, Hahn MW, Lanfear R 2020. gCF. *Mol Biol Evol* 37:2727. doi:10.1093/molbev/msaa106. Lanfear R, Hahn MW 2024. *Mol Biol Evol* 41:msae214. doi:10.1093/molbev/msae214.
- Kubatko LS, Degnan JH 2007. *Syst Biol* 56:17. doi:10.1080/10635150601146041. Mendes FK, Hahn MW 2018. *Syst Biol* 67:158. doi:10.1093/sysbio/syx063.
- Şapcı AOB, Arasti S, Braun EL, Mirarab S 2026. Phlag. *Bioinformatics* 42:btag273. doi:10.1093/bioinformatics/btag273 (PMC13340173).
- Kuhl H et al. 2021 (online 2020). 3′-UTRs. *Mol Biol Evol* 38:108. doi:10.1093/molbev/msaa191 (PMC7783168).
- Árvores-fonte (DOIs pelo `studies/find_studies` e Crossref): Hackett 2008 doi:10.1126/science.1157704; Jarvis 2014 doi:10.1126/science.1253451; Prum 2015 doi:10.1038/nature15697; Suh 2015 doi:10.1371/journal.pbio.1002224; Reddy 2017 doi:10.1093/sysbio/syx041; Kimball 2019 doi:10.3390/d11070109; Feng 2020 doi:10.1038/s41586-020-2873-9; Wu 2024 doi:10.1073/pnas.2319696121; Claramunt 2025 doi:10.1016/j.cub.2025.07.004; Zhao 2025 doi:10.1093/sysbio/syaf080 (Syst Biol, fascículo de 2026); Fuchs 2015 doi:10.1016/j.ympev.2014.08.010; demais em `runs/rapinantes/reanalise/opentree/estudos.tsv`.

## 9. Arquivos

Resultados (Windows): `C:\Users\Manec\bio-stack\runs\rapinantes\reanalise\`
- `tabela_sintese.tsv` — tabela-síntese completa (todas as análises, todos os testes); `tabelas_md.md` — tabelas de quartetos em Markdown.
- `b10k\quartetos_ramos.tsv` — q1–q3, pp, EN, localPP, p de politomia, p 1ª×2ª e p de simetria por ramo e análise; `b10k\arvores_publicadas_ramos.tsv` — resolução e topologia enraizada A, S, C, U nas 11 árvores publicadas e nas árvores desta reanálise; `b10k\gcf_ramos.tsv` e `b10k\gcf\` — gCF; `b10k\topologias\` — 8 topologias alternativas e `escores_*.tsv` (escore total de quartetos); `b10k\aster\` — árvores anotadas do ASTER, `freqQuad.csv` e logs; `b10k\publicados\` — árvores de espécies, teste de politomia, árvore datada e dados-fonte da Fig. 2 do artigo; `b10k\brutos_manifesto\` — URLs, tamanhos, sha256 e mapas de nomes; `b10k\taxa_focais.tsv`.
- `opentree\nos.tsv`, `opentree\node_info_*.json`, `opentree\estudos.tsv` — Open Tree 16.1; `opentree\votos_por_arvore.tsv`, `opentree\votos_resumo.json`, `opentree\arvores_fonte\` (28 árvores-fonte baixadas) — meta-análise; `opentree\aves_1.6\nos_aves16.tsv` e arquivos da Aves 1.6.
- `timetree\timetree_pares.tsv`, `timetree\pairwise_*.json`, `timetree\idades_arvores_datadas.tsv`.
- `figuras\quartetos_ramos.png|svg` (q por ramo, tipo de dado e estimador) e `figuras\opentree_afroaves.png|svg` (base de Afroaves nas 21 árvores-fonte informativas).
- `versoes.txt` — versões lidas das ferramentas.

Scripts: `C:\Users\Manec\bio-stack\scripts\rapinantes\reanalise\` — `preparar_b10k.py`, `podar_arvores.py`, `rodar_aster_b10k.sh`, `pontuar_alternativas.sh`, `teste_ponderacao.sh`, `gcf_b10k.sh`, `copiar_resultados_b10k.sh`, `sintese_b10k.py`, `opentree_nos.py`, `opentree_votos.py`, `aves_nos.py`, `timetree_pares.py`, `idades_arvores_datadas.py`, `tabela_sintese.py`, `tabelas_md.py`, `figura_quartetos.py`.

Dados brutos e árvores de genes traduzidas/podadas (WSL, não copiados para o projeto): `~/rapinantes/reanalise/b10k/{brutos,genes,aster,gcf}`.
