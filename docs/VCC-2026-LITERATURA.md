# VCC 2026 — Levantamento de literatura: o que as equipes de topo fizeram

**Data da busca:** 2026-10-22
**Pergunta de pesquisa:** para a fase final do Virtual Cell Challenge 2026 (predição
zero-shot de resposta a CRISPRi em três contextos celulares não vistos, 300 alvos,
18.533 genes, prazo 2026-11-05T23:59Z), existe um salto de patamar documentado que um
modelo de transferência estatística em pseudobulk está deixando na mesa?

**Convenção deste documento.** Cada afirmação numérica carrega a fonte primária. Blocos
marcados **[INFERÊNCIA]** são raciocínio meu sobre fontes, não achado de fonte. Onde a
busca não encontrou nada, está escrito **"não encontrei fonte"** — isso é resultado, não
lacuna a preencher.

---

## 0. Método da busca

### Bases consultadas e queries, literalmente

| Base | Query | Resultados úteis |
|---|---|---|
| `library_search` (acervo local) | "Virtual Cell Challenge predição de resposta transcricional a perturbação CRISPRi em contexto celular novo" | **0 — a ferramenta devolveu erro de execução.** O acervo local não foi consultado nesta revisão. |
| WebSearch | `Virtual Cell Challenge 2025 winners method Arc Institute results` | 9 links, 3 primários |
| WebSearch | `Ahlmann-Eltze deep learning predictions gene perturbation effects do not yet outperform simple linear baselines` | 10 links, 1 primário (Nature Methods) |
| WebSearch | `Virtual Cell Challenge 2026 zero-shot unseen cell context baseline score leaderboard arcinstitute` | 9 links, 3 primários |
| WebSearch | `Kernfeld benchmark perturbation prediction transcriptome "no better than" baseline causal network expression forecasting` | 9 links |
| WebSearch | `Csendes perturbation prediction single cell benchmark "mean" baseline critique evaluation metrics misleading` | 9 links, 3 primários |
| WebSearch | `cross-cell-line transferability perturbation response Perturb-seq predicts which cell line transfer basal expression similarity not lineage` | 9 links, **3 primários diretamente na pergunta 5** |
| WebSearch | `perturbation prediction single cell "under-dispersed" OR "variance compression" calibrate dispersion negative binomial sampling to recover differential expression Jaccard significant genes` | 9 links, 2 primários |
| WebSearch | `arXiv 2026 "Virtual Cell Challenge" zero-shot unseen cell context method preprint CRISPRi 300 perturbations six metrics` | 9 links, 0 preprint de método 2026 |
| WebSearch | `turbine.ai "regression model" "Virtual Cell Challenge" top leaderboard blog` | 10 links, 2 relatos de equipe |
| WebSearch | `"Virtual Cell Challenge" 2025 final leaderboard score "overall score" baseline State model 0.4 xTrimoSCPerturb score` | 9 links, placar final |
| Europe PMC (`bc_get_europepmc_articles`) | `title:Virtual Cell Challenge Toward a Turing test for the virtual cell` | 8 hits, 1 alvo (PMID 40578317) |
| Europe PMC | `title:Predicting cellular responses to perturbation across diverse contexts with State` | 36 hits, 1 alvo (PMID 42673963) |
| Europe PMC | `title:Mechanisms Matter: Transportability of Cellular Perturbation Effects` | 4 hits, 1 alvo (PPR1210489) |
| bioRxiv API (`bc_get_biorxiv_preprint_details`) | DOIs `10.64898/2026.05.08.723625`, `10.64898/2026.07.24.740459`, `10.64898/2026.07.23.740433`, `10.1101/2024.12.23.630036` | 4/4 recuperados |
| PubMed (`article_getter`) | PMID `42648290` | 1 |

**Nota de método.** bioRxiv e cell.com devolvem HTTP 403 para busca programática; os
preprints foram recuperados pela API do bioRxiv e pelo Europe PMC. A página do placar
(`virtualcellchallenge.org/2025`, `/leaderboard`, `/evaluation`) é renderizada por
JavaScript e devolve só "Loading…" — os números do placar vieram do índice de busca da
própria página, não da página renderizada. Isso está marcado onde aparece.

---

## 1. O que ganhou o Virtual Cell Challenge de 2025, e como

Fonte primária única e completa: o relatório pós-competição do Arc Institute,
[Virtual Cell Challenge 2025 Wrap-Up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up).
Tudo nesta seção vem de lá salvo onde indicado.

### Escala e tarefa de 2025

- Mais de 5.000 inscritos de 114 países; mais de 1.200 equipes submeteram; mais de 300
  fizeram submissão final.
- Dado: ~300.000 células de H1 (células-tronco embrionárias humanas), 300 perturbações
  CRISPRi, química 10x Flex, >50.000 UMI/célula, ~1.000 células por perturbação.
- Métricas de 2025: **três** — PDS (perturbation discrimination), DES (differential
  expression) e MAE. O Prêmio Generalista usou **sete**, somando Pearson delta, Spearman
  de log-fold change, AUPRC e Spearman de effect size, as quatro vindas da avaliação do
  modelo State.
- Split (relato de participante, não do Arc): 150 perturbações de treino, 50 de
  validação, 100 de teste final —
  [yraj.in, 2025-12-17](https://yraj.in/ai/virtual-cell/perturbation-modeling/single-cell/2025/12/17/virtual-cell-challenge.html).

**Diferença estrutural com 2026, que condiciona toda a leitura abaixo:** em 2025 havia
150 perturbações de treino **no próprio contexto-alvo** (H1). Em 2026 não há conjunto de
treino nenhum e o contexto é inédito. Nenhum dos métodos abaixo foi avaliado no regime de
2026.

### Placar final de 2025

| Posição | Equipe | Modelo | Overall |
|---|---|---|---|
| 1 | BM_xTVC (BioMap Research) | xTrimoSCPerturb | 34,0 |
| 2 | XLearning Lab (Sichuan Univ.) | — | 32,8 |
| 3 | Outlier (UChicago/Dartmouth/HKU) | TransPert | 32,6 |

Números do placar em `virtualcellchallenge.org/2025`, **lidos pelo índice de busca da
página, não pela página renderizada** — a página devolve "Loading…" a cliente sem JS. As
colunas são rank, team, model, overall, des, pds, mae. Os valores por métrica não foram
recuperados. Prêmios: 1º $100.000, 2º $50.000, 3º $25.000, Generalista $100.000 (Altos
Labs, modelo `go-with-the-flow`) — estes do relatório do Arc.

### 1º lugar — xTrimoSCPerturb (BM_xTVC)

Arquitetura e truques, do relatório do Arc:

- Modelo pré-treinado de célula única baseado em arquitetura **scFoundation melhorada**,
  combinado com **embeddings de modelo de proteína** e conjuntos públicos de perturbação
  curados.
- Encoder otimizado para **recuperar sinal biológico de zeros técnicos**.
- Decoder com **cross-attention desemaranhada**.
- **Features explícitas estatísticas:** frequência de DEG no conjunto de treino e nível
  médio de expressão.
- **Perda fundida:** PDS + DES + MAE com peso pequeno.
- Treino em **pseudobulk**, com as saídas **repetidas várias vezes para satisfazer as
  exigências da avaliação DES**.
- A própria equipe registrou que abordagens puramente de IA não superaram
  consistentemente baselines estatísticos em DES e MAE.

Artigo completo prometido, **não publicado até 2026-10-22** — não encontrei fonte com o
método em detalhe.

### 2º lugar — XLearning Lab

- Geração condicional sobre **pseudobulk**, não célula única.
- Rede totalmente conectada sobre: expressão de controle agregada + **embeddings ESM-2 do
  gene perturbado** + indicador de contagem de UMI.
- **Aprendizado residual** — prediz o delta da perturbação, não o valor absoluto.
- Perda MSE com objetivos auxiliares.
- Treinou **exclusivamente em Perturb-seq público**, incluindo um pequeno conjunto de H1
  do PerturbAtlas.
- Priorizou PDS, que a equipe julgou valer **~2× o peso de DES**.

### 3º lugar — TransPert (equipe Outlier)

O mais próximo do nosso método atual, e o mais instrutivo:

- Framework **estatístico**, sem rede profunda.
- Entradas **somente de nível de resumo**: perfis de perturbação em pseudobulk e
  resultados de teste de Wilcoxon.
- Agrega resumos de perturbação por gene de **múltiplas linhagens de referência** com
  **agregação ciente de similaridade**.
- Prediz **quais genes têm chance de ser significativamente afetados** no contexto-alvo.
- Aplica **escala linear global**, com fatores escolhidos por validação cruzada no dado
  de treino do desafio, **para otimizar PDS**.
- Treinou exclusivamente em Perturb-seq público de linhagens celulares.

Os mesmos autores publicaram a análise de escala do PDS discutida na seção 4
([arXiv:2511.16954](https://arxiv.org/abs/2511.16954)).

### Prêmio Generalista — Altos Labs, `go-with-the-flow`

- **Flow matching** direto no espaço de expressão gênica; campo de velocidade
  parametrizado por **U-Net** customizada.
- Pré-treino em **~7 milhões** de células de alta qualidade de fontes públicas (incluindo
  o dado de treino do desafio) e screens internos de perturbação da Altos — **nenhum em
  H1**. Depois fine-tuning num conjunto menor específico de H1 e das perturbações do teste.
- Selecionado por **maior rank médio nas sete métricas** entre as 50 entradas finais com
  melhores pontuações nas métricas do desafio.

O preprint correspondente é **PRiMeFlow**,
[arXiv:2604.13986](https://arxiv.org/abs/2604.13986) (Yan, Wu, Ji, Agrahar, Wershof,
Nassar, Sadria, Eksi, Trifonov, Ibarra, Felgueira, Osiński, Stark) — **preprint, não
passou por revisão por pares**. Diz que a arquitetura PRiMeFlow foi a base do modelo
vencedor do Prêmio Generalista. Ablações relevantes, do próprio preprint:

- Espaço de expressão gênica **bate** espaço latente de PCA (ajuste significativamente
  mais próximo da distribuição real).
- **U-Net bate MLP** — o MLP "struggles to achieve decent fits" nas métricas
  distribucionais e de pseudobulk.
- **Acoplamento de transporte ótimo não acrescenta** desempenho de forma relevante.
- O preprint **não reivindica generalização zero-shot**; avalia transferência de
  covariável e lista contexto inédito como trabalho futuro.
- O preprint **não reporta nenhum número de PDS/DES/MAE do desafio** — só "outstanding
  performance" em H1. Verificado lendo os primeiros 100.000 caracteres do HTML.

### A equipe que chegou ao topo com regressão — relato útil por ser próximo de nós

Turbine AI, [blog, 2025-12-09](https://blog.turbine.ai/p/how-did-we-get-a-regression-model)
(post de empresa, não revisado por pares):

- Ficou **15º no geral**, e liderou o placar por um período. Modelo descrito como "um
  modelo simples descrito pela primeira vez em 1970", **sem pré-treino**, operando em
  **clusters de pseudobulk**, não células individuais.
- **O poder preditivo veio quase todo de Perturb-seq externo** (Replogle et al., Nadig
  et al.) usado **como feature de entrada**, não como dado de treino. Semanas de limpeza
  e alinhamento.
- Argumento central: o gargalo é dado, não tamanho de modelo; o detalhe de célula única
  acrescenta pouco além do número efetivo de pontos por condição.
- Write-up técnico ("Mean Predictors") existe em
  `documents.turbine.ai/view/57dc116a…/vcc-writeup-MP`, mas é um visualizador de PDF em
  JS — **não consegui extrair o conteúdo**. Os números de ablação e de métrica desse
  write-up ficam pendentes.

### Lições que o Arc tirou — a parte mais acionável do relatório

Citando curto e atribuindo: o Arc afirma que modelos **híbridos**, combinando deep
learning com features estatísticas, superaram redes neurais puras, e que o desenho da
função de perda importou tanto quanto a arquitetura. Também: features multimodais,
sobretudo embeddings de proteína, agregaram valor em todas as equipes de topo; nenhuma
métrica isolada captura qualidade de modelo, e otimizar uma frequentemente prejudica
outras.

**Resposta curta à pergunta 1:** o salto de patamar de 2025 não veio de arquitetura. Veio
de (a) treinar em pseudobulk, (b) aprender o resíduo/delta, (c) injetar Perturb-seq
público como feature, (d) embeddings de proteína para o gene perturbado, (e) perda
desenhada para as métricas pontuadas e (f) escala calibrada por validação cruzada. Três
dos quatro premiados são essencialmente isso. O quarto (Altos) é generativo e ganhou um
prêmio de consistência entre métricas, não o primeiro lugar.

---

## 2. O baseline oficial e o que ele dimensiona

**Pontuação do `baseline.py` de 2026 (sorteia 400 controles, escala por gene): não
encontrei fonte.** Nem o Arc nem terceiros publicaram o valor. O anúncio de 2026 não traz
código nem pontuação de baseline
([Arc, 2026](https://arcinstitute.org/news/virtual-cell-challenge-2026)), e busca
direcionada por `baseline.py` não devolveu nada.

O que **tem** fonte, e serve para dimensionar:

1. **A escala de 2025 era normalizada pela média de perturbação.** O placar de 2025
   declara que cada uma das três métricas é normalizada em relação ao baseline de média
   de perturbação, para medir ganho sobre o baseline ingênuo (índice de busca de
   `virtualcellchallenge.org/2025`). Então o **1º lugar, 34,0, significa ~34% do caminho
   entre o baseline de média e o teto** — em um ano em que havia 150 perturbações de
   treino no próprio contexto.
2. **O teto é réplica real.** Em 2026, as seis métricas são escaladas entre a média do
   contexto celular e dado de um **experimento de réplica real**
   ([Arc, 2026](https://arcinstitute.org/news/virtual-cell-challenge-2026)). O
   `cell-eval` implementa esse teto como um auto-split de dado real com correção de
   Spearman-Brown `r' = 2r/(1+r)`, aplicada **só** às métricas da lista `SB_METRICS`; as
   demais vêm `NaN` ([README do cell-eval](https://github.com/ArcInstitute/cell-eval)).
3. **Em 2025 o MAE saiu da competição por si.** O Arc registra que quase todos os modelos
   pontuaram **pior que o baseline** em MAE, o que efetivamente removeu o MAE da
   otimização competitiva. Em 2026 há duas métricas de erro (`mse`, `nmae`) de seis.
4. **O baseline de 2025 fornecido pelo Arc era o State**, não um sorteio de controles. Um
   blog de terceiros reporta DES 0,22 / PDS 0,62 / MAE 0,04 para esse baseline
   ([pebblous.ai](https://blog.pebblous.ai/report/nvidia-virtual-cell-challenge-2026-05/en/))
   — **fonte secundária, não confirmada no Arc; não usar como número de referência.**

**[INFERÊNCIA]** Dois cálculos de dimensionamento, meus, não de fonte:

- Se o campeão de 2025 ficou em ~0,34 **com** dado de treino no contexto, o teto
  realista para 2026 **sem** nenhum é mais baixo, e os 0,472 do 1º lugar na validação de
  2026 que medimos são surpreendentemente altos — provavelmente porque a validação tem
  placar vivo e permite ajuste iterativo contra ela, o que o teste final não permite.
- Nossa melhora de ~4,5% de MSE sobre a resposta média não se converte linearmente em
  `s`. `s = (u−b)/(r−b)` para `mse` depende de quanto a réplica real `r` bate a média `b`
  em MSE, número que **não temos**. Medir `r` localmente (duas metades dos controles, ou
  o auto-split do `cell-eval`) é o único jeito de saber se 4,5% vale 0,05 ou 0,15 nessa
  métrica. **Isso é a medição mais barata e mais informativa que falta.**

---

## 3. Predição de resposta em contexto celular NÃO visto — os dois lados, com números

Esta é a pergunta onde a literatura **discorda**, e o resumo honesto é que o lado crítico
está mais bem documentado que o lado otimista.

### Afirmação A — "modelos profundos não batem baselines simples em predição de perturbação"

**Sustentam:**

- **Ahlmann-Eltze, Huber & Anders**, *Nature Methods* 2025;22:1657–1661,
  [doi:10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6)
  (**revisado por pares**; preprint
  [bioRxiv 2024.09.16.613342](https://www.biorxiv.org/content/10.1101/2024.09.16.613342v1)).
  Cinco modelos fundacionais (scGPT, scFoundation, scBERT, Geneformer, UCE) mais GEARS e
  CPA contra baselines deliberadamente simples, em perturbação simples e dupla.
  **Nenhum superou os baselines.** Métrica: distância L2 no perfil dos 1.000 genes mais
  expressos. Baselines para perturbação dupla: "sem efeito" e modelo aditivo. Código
  público com arquivo no Zenodo.
- **Csendes et al.**, *BMC Genomics* 2025,
  [doi:10.1186/s12864-025-11600-2](https://doi.org/10.1186/s12864-025-11600-2)
  (**revisado por pares**). O baseline mais simples possível — média dos exemplos de
  treino — **superou scGPT e scFoundation**. Avaliação em pseudobulk, com Pearson em
  espaço de expressão e em espaço de expressão diferencial. Os autores registram que usar
  dado de célula única não deu vantagem óbvia sobre pseudobulk nessas linhagens de baixa
  heterogeneidade.
- **Kernfeld**, tese de doutorado, Johns Hopkins,
  [*Assessing Reliability of Causal Models of Transcription*](https://jscholarship.library.jhu.edu/server/api/core/bitstreams/d18fbb3c-9703-4296-b3e8-e930324fd5cf/content).
  Baselines triviais **tipicamente com erro de predição menor** que modelos dedicados;
  proporções de falsa descoberta muito acima do esperado. **Tese, não artigo revisado por
  pares;** li o resumo, não as tabelas. Citada por TxPert e PerturBench para a afirmação
  de que a média/mediana do treino prediz melhor que GEARS.
- **Li et al.**, *Science Advances* 2026,
  [doi:10.1126/sciadv.aed3414](https://doi.org/10.1126/sciadv.aed3414)
  (**revisado por pares**; preprint
  [bioRxiv 2024.12.23.630036](https://www.biorxiv.org/content/10.1101/2024.12.23.630036)).
  13 métodos e baselines, **25 conjuntos de dados, 24 métricas, 3 tarefas**. Na **tarefa
  3, transferência entre tipos celulares, nenhum método generalizou consistentemente**, e
  o resultado foi dominado por consistência biológica. Os modelos compartilham **viés
  conservador**; modelos fundacionais com fine-tuning **comprimem variância** e
  subestimam sinergia em combinações. PerturbNet foi o melhor em recuperar assinaturas de
  DE nas tarefas 1 e 2.
- **Qi & Chapfuwa** (Microsoft), *Mechanisms Matter: Transportability of Cellular
  Perturbation Effects*, bioRxiv 2026-05-12,
  [doi:10.64898/2026.05.08.723625](https://doi.org/10.64898/2026.05.08.723625)
  (**preprint, não revisado**). Seis modelos profundos e cinco baselines simples, em
  Perturb-seq real e semi-sintético com estrutura causal conhecida. **Sob split
  cross-context o desempenho cai substancialmente, muitas vezes ao nível dos baselines
  simples.** E, decisivo: **mesmo em dado sintético com estrutura causal totalmente
  especificada, nenhum modelo generaliza entre contextos com mecanismos causais
  diferentes.** Código em `github.com/microsoft/mechanisms_matter`.

### Afirmação B — "há ganho real de modelo em contexto novo"

**Sustentam:**

- **Adduri et al. (State)**, *Cell* 2026;189(19):5914–5931.e20,
  [doi:10.1016/j.cell.2026.07.052](https://doi.org/10.1016/j.cell.2026.07.052),
  PMID 42673963 (**revisado por pares**). State melhorou **discriminação de efeitos em
  mais de 30%** em conjuntos grandes e identificou genes diferencialmente expressos com
  acurácia significativamente melhor que baselines. Embeddings de célula treinados em
  dado observacional de **167 milhões de células** permitem **identificar perturbações
  fortes em contextos onde nenhuma perturbação foi observada no treino**.
  **Leitura crítica:** a afirmação zero-shot do State é *identificar perturbações fortes*
  — ou seja, magnitude/triagem — não predizer o perfil completo de resposta. É mais fraca
  do que o que o VCC 2026 pede. Primeira autoria e financiamento são do Arc Institute,
  que organiza o desafio e desenha as métricas: conflito de interesse estrutural a
  declarar.
- **Molina & Zhang** (AITHYRA/Academia Austríaca de Ciências), bioRxiv 2026-07-27,
  [doi:10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459)
  (**preprint, não revisado**). Reportam que alinhar priors biológicos à geometria da
  resposta **permite que modelos lineares ou MLP simples superem arquiteturas
  estado-da-arte** em vários cenários de generalização, **incluindo linhagens celulares
  não vistas**. Note que isso sustenta a afirmação B num sentido torto: o ganho existe,
  mas vem do prior, não da profundidade — o que é uma terceira posição, não o lado
  otimista clássico.
- **TxPert**, [arXiv:2505.14919](https://arxiv.org/abs/2505.14919) (**preprint, sem
  venue revisado listado**). Testa explicitamente três cenários OOD: perturbação simples
  nova, dupla nova e **linhagem celular nova**, usando grafos de conhecimento gene-gene.
  A página de resumo **não lista baselines nem números**; não extraí resultados
  quantitativos e não os afirmo.

### Estado da afirmação

**Contestado, com o peso da evidência revisada por pares no lado crítico.** Duas
publicações revisadas (Ahlmann-Eltze *Nat Methods*; Csendes *BMC Genomics*) e uma terceira
com cobertura de 25 conjuntos (Li *Sci Adv*) dizem que o baseline empata ou ganha, e a
*Sci Adv* diz isso **especificamente na tarefa de transferência entre tipos celulares**.
Do lado otimista há uma publicação revisada (State, *Cell*) com conflito de interesse
estrutural e uma reivindicação zero-shot mais fraca que a tarefa, e dois preprints.

**O que isso significa para 14 dias, direto:** não encontrei nenhuma fonte que mostre um
modelo profundo batendo baselines lineares/de média no regime **simultâneo** de contexto
novo e perturbação sem evidência medida. A hipótese de que existe um salto de patamar
arquitetural a ser capturado em duas semanas **não é sustentada pela literatura**. O salto
documentado está em prior e em calibração, não em profundidade.

---

## 4. Como as métricas de DE são ganhas, e a tensão com o encolhimento ótimo para MSE

### A tensão existe e está documentada — mas de forma fragmentada

**Sustentam que métricas de erro não ponderadas premiam colapso/sub-dispersão:**

- **Mejia et al.**, *Diversity by Design*, [arXiv:2506.22641](https://arxiv.org/abs/2506.22641)
  (2025-06-27, **preprint**). Afirmam que deltas referenciados ao controle e métricas de
  erro não ponderadas **premiam mode collapse** sempre que o controle é viesado ou o sinal
  biológico é esparso. Propõem métricas cientes de DEG — **WMSE** e **R²_w(Δ)** — e
  baselines negativo e positivo para calibrá-las; com essa calibração **o baseline de
  média pontua em desempenho nulo**. E propõem **WMSE como função de perda**, reportando
  que reduz mode collapse. A página de resumo **não traz números**; os valores ficam
  pendentes de leitura do texto completo.
- **Li et al.**, *Sci Adv* 2026 (acima): modelos fundacionais com fine-tuning
  **comprimem variância** e subestimam efeitos, e os autores atribuem isso ao fato de
  **preservação de distribuição não fazer parte do objetivo de fine-tuning**. Esta é a
  citação limpa para "encolhimento ótimo para erro quadrático produz predição
  sub-dispersa".
- **Qi & Chapfuwa** 2026 (acima): adaptam o **Vendi score** ao cenário de perturbação
  para quantificar diversidade efetiva de resposta e **diagnosticar mode collapse, um
  modo de falha invisível às métricas por perturbação**.

**Sustentam que o lado da emissão é mecanicamente decisivo:**

- **README do `cell-eval`** ([ArcInstitute/cell-eval](https://github.com/ArcInstitute/cell-eval)):
  células duplicadas **afundam** as métricas de recuperação (recall, overlap, AUC) e
  fazem as **métricas de DE com gate de FDR super-chamarem**. Ou seja: o número de
  células distintas emitidas e a dispersão entre elas entram direto em `jac` e `reach`,
  independentemente de o delta médio estar certo.
- **BM_xTVC** (relatório do Arc): treinou em pseudobulk e **repetiu as saídas várias
  vezes para satisfazer as exigências da avaliação DES**. O campeão tratou a emissão como
  problema separado do modelo.
- **TransPert** (relatório do Arc): **escala linear global** com fator escolhido por
  validação cruzada **para otimizar PDS**, mais um passo dedicado a predizer **quais
  genes** são significativamente afetados. Dois estágios: direção/conjunto, e depois
  magnitude.

**Sustentam que a magnitude era a alavanca do PDS em 2025, e que isso mudou:**

- **Liu, Zhang, Du, Zhao & Wang**, [arXiv:2511.16954](https://arxiv.org/abs/2511.16954)
  (2025-11-21, **preprint**; mesma equipe do TransPert, 3º lugar). O PDS é **altamente
  sensível à escolha da medida de distância e à escala dos efeitos preditos**; PDS
  baseado em ℓ1 e ℓ2 se comporta de forma muito diferente de medidas baseadas em cosseno
  **mesmo após casamento de norma**. O relatório do Arc resume o achado da equipe: à
  medida que a magnitude da predição cresce, o PDS ℓ1 se aproxima de similaridade de
  cosseno baseada em sinal.
- **Turbine AI** (blog, acima): as equipes descobriram que, sob PDS, **a magnitude das
  predições importa mais que a direção**, o que convidou a gaming. E registram que
  escalar transformações subia PDS sem ser biologicamente mais correto.
- **Dall'Olio** ([substack, parte 4](https://gmdbioinformatics.substack.com/p/arc-virtual-cell-challenge-part-4),
  comentário de terceiro): descreve um write-up de participante que otimizou DES notando
  que **genes altamente expressos no controle têm chance alta de serem diferencialmente
  expressos**, e predizendo o conjunto de DE a partir da expressão basal. O autor chama
  isso de "leaderboard hacking".

**Alguém documentou como resolver a tensão numa submissão única?**

**Não encontrei fonte** que trate o problema como otimização explícita de compromisso
entre erro e recuperação de conjunto de DE numa única submissão. O que existe são duas
soluções de engenharia usadas por vencedores — **perda fundida multi-métrica**
(BM_xTVC: PDS + DES + MAE com peso pequeno) e **calibração de escala em estágio separado
por validação cruzada** (TransPert) — e uma proposta metodológica de perda ponderada por
DEG (Mejia et al., WMSE) que **não foi avaliada no VCC**. Há também um trabalho de 2026
que organiza o problema sem resolvê-lo: **Schäfer et al.**, *Towards Principled
Evaluation of Single-Cell Perturbation Prediction Models*, bioRxiv 2026-07-27,
[doi:10.64898/2026.07.23.740433](https://doi.org/10.64898/2026.07.23.740433)
(**preprint**), que decompõe protocolos de avaliação em representação, métrica,
transformação de score e estratégia de reporte, e entrega o pacote `scPertEval` com
implementações de referência sobre sete conjuntos públicos.

**[INFERÊNCIA]** A consequência operacional para 2026 é específica e, acho, subexplorada:
a métrica de 2026 é `pds_cosine`, **não** o PDS ℓ1 de 2025. Cosseno é invariante a escala
global. Se isso vale literalmente, **a alavanca de magnitude que dominou 2025 deixou de
exister no PDS** — e o encolhimento passa a ser disputado só entre `mse`/`nmae` (que
querem encolher) e `jac`/`reach`/`fid` (que querem dispersão realista). Isso é mais limpo
que 2025: são dois grupos com direções opostas, não três. Não encontrei fonte que afirme
isso sobre o preset `vcc2026`; é dedução a partir do nome da métrica e do achado de Liu
et al. Mede-se em uma tarde no avaliador local, escalando a predição por um fator `k` e
vendo quais métricas se movem.

---

## 5. O que prediz transferência entre linhagens, se não é proximidade

Boa notícia: **a literatura concorda com a medição de vocês**, e há três fontes que
propõem substitutos concretos para "proximidade de linhagem".

### A afirmação negativa: similaridade distribucional não é o que governa

- **Qi & Chapfuwa** 2026, [doi:10.64898/2026.05.08.723625](https://doi.org/10.64898/2026.05.08.723625)
  (**preprint**). Usando teoria de **transportabilidade causal**, mostram que
  generalização entre contextos é governada por **mecanismos causais compartilhados, não
  por mera similaridade distribucional**. O substituto empírico que propõem é medível:
  **menor concordância cross-context do impacto transcriptômico amplo acompanha gaps de
  transferência maiores**. Isto é, não se mede a distância entre os estados basais — se
  mede a concordância entre as **respostas** nos dois contextos, nas perturbações que os
  dois têm em comum.
- **Li et al.**, *Sci Adv* 2026: na tarefa de transferência entre tipos celulares,
  nenhum método generalizou consistentemente, e o resultado foi dominado por
  **consistência biológica**, não por capacidade de modelo.

### A afirmação positiva mais útil: decomposição da resposta

**Molina & Zhang** 2026, [doi:10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459)
(**preprint, não revisado**) — aplicado a **quatro screens Perturb-seq de CRISPRi em
múltiplas linhagens**. Decompõem a resposta transcricional em quatro componentes:
**global**, **específico da perturbação**, **específico da linhagem** e **interação
perturbação × linhagem**. Os achados, nas palavras deles reescritas:

- O **componente global é de baixa dimensão**, reflete programas recorrentes de
  **proliferação e resposta a estresse**, e **pode ser inferido da expressão de
  controle**.
- Os componentes **específicos de perturbação e de linhagem são de alta dimensão e NÃO
  podem ser recuperados da expressão de controle**.
- Por isso desenvolvem modelos alinhados por componente, mapeando priors biológicos —
  notadamente **coessencialidade gênica** — sobre a geometria das respostas observadas.
- Esse alinhamento **permite que modelos lineares ou MLP superem arquiteturas
  estado-da-arte** em múltiplos cenários de generalização, incluindo **linhagens não
  vistas** e combinações de perturbações não vistas.
- Conclusão deles: generalização depende **de casar a informação biológica ao componente
  certo da resposta**, não de complexidade de modelo.

Resumo de resultados secundário sobre a mesma preprint, do índice de busca (não
verificado no texto completo, marcar como tal): coessencialidade funciona como sistema de
coordenadas funcional externo que alinha as bases de resposta de doador e alvo apesar de
sobreposição limitada entre os programas de expressão; **o grau de alinhamento das bases
de perturbação não se correlaciona com desempenho do modelo**; e PCs de expressão de
controle capturam bem a direção do template mas pouco do resíduo após remoção do template.

### Transferência com confiabilidade estimada em âncoras

**PerturbMap**, [arXiv:2607.28090](https://arxiv.org/abs/2607.28090) (Cui, Liu & Sun,
2026-07-30, **preprint**). Exatamente o problema de vocês: a perturbação foi medida em
outros contextos e falta no contexto de destino.

- Combina uma **base local de baixo rank** no contexto de destino com **propostas
  aceitas** que transportam as respostas medidas nos contextos doadores.
- Os transportes usam **experts de ridge doador→destino** treinados em perturbações
  pareadas.
- **O peso de cada proposta vem de confiabilidade de rota estimada em âncoras de
  validação** — não de similaridade de linhagem, não de profundidade do doador.
- Números: **+4,1% de MSE de efeito completo** sobre a base local de baixo rank, em
  Perturb-CITE-seq de melanoma. Recuperação top-10 do contraparte no mesmo destino por
  cosseno sobe de **74,5% para 80,5%**. Bate baselines FedAvg, resposta zero, cópia crua,
  cópia calibrada e controles afins com embaralhamento de identidade.

**[INFERÊNCIA]** Dois casamentos com os números de vocês, meus:

1. O ganho de +4,1% de MSE do PerturbMap sobre a base local está na **mesma ordem** dos
   ~4,5% de vocês sobre a resposta média. Isso sugere que ~4-5% é o patamar do que
   transferência estatística bem-feita entrega neste problema, e que esperar um salto
   de patamar por reorganização do mesmo mecanismo é otimismo. O salto, se vier, vem de
   tratar métricas diferentes de modo diferente (seção 4) e de melhorar o **componente
   global**, não a parcela específica.
2. A decomposição de Molina & Zhang é a explicação exata do achado de vocês de que a
   resposta média vale 0,252 e o desvio específico 0,045: a média é o **componente
   global** (baixa dimensão, inferível do controle, transfere — correlação mediana +0,470
   entre contextos) e o desvio específico são os **componentes de alta dimensão** que a
   preprint diz não serem recuperáveis da expressão de controle. A razão ~5× que vocês
   mediram é o preço empírico dessa distinção.

---

## 6. O que mudaríamos em 14 dias

Cinco itens, ordenados por ganho esperado sobre esforço. Cada um com a evidência que o
sustenta e o que é inferência.

### 1. Medir `r` localmente e otimizar o encolhimento **por grupo de métrica**, não um β único

Construir o teto de réplica em casa (duas metades dos controles / auto-split do
`cell-eval`), varrer um fator de escala `k` sobre a predição e registrar as seis métricas
escaladas separadamente. Hoje o β é escolhido por MSE; `jac`, `reach` e `fid` querem o
oposto, e `pds_cosine` provavelmente não quer nada.

**Evidência:** Liu et al. [arXiv:2511.16954](https://arxiv.org/abs/2511.16954) — PDS é
altamente sensível a métrica de distância e escala, e versões ℓ1/ℓ2 diferem de cosseno
mesmo após casamento de norma. Li et al. *Sci Adv* 2026 — encolhimento por objetivo de
erro comprime variância porque preservação de distribuição não está no objetivo. TransPert
(3º lugar, Arc wrap-up) fez exatamente isso: escala linear global com fator por validação
cruzada, em estágio separado. `cell-eval` README — o teto é auto-split com correção de
Spearman-Brown, logo é reproduzível localmente.
**[INFERÊNCIA]:** que `pds_cosine` seja invariante a escala e remova a alavanca de
magnitude de 2025. Verificável na primeira varredura de `k`.
**Esforço:** baixo. Uma varredura, o avaliador local já está no caminho crítico.

### 2. Calibrar a emissão estocástica pelo número real de genes significativos

Vocês já sabem que emitir células idênticas dá ~9.701 genes "significativos" contra
~2.739 reais. O item é fechar o laço: tratar a dispersão da emissão como **um escalar a
ajustar** contra `jac`/`reach`/`fid` no avaliador local, com o número real de DE do
contexto como alvo, em vez de herdar o nível de ruído do `baseline.py`.

**Evidência:** README do `cell-eval` — células duplicadas super-chamam as métricas de DE
com gate de FDR e afundam as métricas de recuperação. BM_xTVC, o campeão de 2025, repetiu
as saídas várias vezes especificamente para satisfazer a avaliação DES (Arc wrap-up).
Mejia et al. [arXiv:2506.22641](https://arxiv.org/abs/2506.22641) — métricas de erro não
ponderadas premiam colapso quando o sinal é esparso. Qi & Chapfuwa 2026 — colapso de modo
é invisível às métricas por perturbação, precisa de medida de diversidade dedicada.
**Esforço:** baixo. Um parâmetro, o gerador CSR já existe.

### 3. Mover o esforço para o **componente global** da resposta, não para o desvio específico

Parar de tentar melhorar os 0,045 do desvio específico e investir em predizer melhor a
resposta média do contexto a partir da expressão de controle — que é onde estão os 0,252,
e que entra nos 300 alvos, não nos 157.

**Evidência:** Molina & Zhang 2026
[doi:10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459) — o
componente global é de **baixa dimensão**, reflete proliferação e estresse, e **é
inferível da expressão de controle**; os componentes específicos de perturbação e de
linhagem são de **alta dimensão e não são** recuperáveis do controle. A validação
biológica que vocês já têm (DDIT3/CHOP, BTG2, HSPA6, DUSP10 — assinatura canônica de
estresse) é literalmente o componente global descrito nessa preprint.
**Esforço:** médio. É modelagem nova, mas de baixa dimensão por construção, e a
preprint diz que linear/MLP basta quando o prior está alinhado.

### 4. Escolher e **pesar** doadores por concordância medida, não por profundidade nem linhagem

Trocar "um melhor doador" por mistura de doadores com peso de confiabilidade estimado nas
perturbações que doador e alvo compartilham — e usar os quatro alvos reservados
(FUBP1, HMGXB4, TIGD5, ZBTB7B) como âncoras de calibração, não só de validação.

**Evidência:** Qi & Chapfuwa 2026 — transportabilidade é governada por mecanismo causal
compartilhado, e o substituto empírico é a **concordância cross-context do impacto
transcriptômico**, que acompanha o gap de transferência. PerturbMap
[arXiv:2607.28090](https://arxiv.org/abs/2607.28090) — peso por **confiabilidade de rota
estimada em âncoras de validação**, com +4,1% de MSE sobre a base local e recuperação
top-10 de 74,5% → 80,5%. TransPert (3º lugar) — **agregação ciente de similaridade** de
múltiplas linhagens de referência, não doador único. Li et al. *Sci Adv* 2026 — na tarefa
de transferência, consistência biológica domina o resultado.
**Ressalva honesta:** quatro âncoras é pouco para estimar confiabilidade de rota. A
versão viável em 14 dias é usar as perturbações compartilhadas entre doadores como
âncoras internas, não os quatro reservados.
**Esforço:** médio.

### 5. Cobrir os 143 alvos sem evidência com coessencialidade, não com embedding cru

Para os alvos sem evidência medida, usar coessencialidade (DepMap) **projetada sobre a
geometria da resposta observada** — base funcional, não feature solta — com ESM-2/STRING
como complemento.

**Evidência:** Molina & Zhang 2026 — coessencialidade mapeada sobre a geometria da
resposta é o que faz linear/MLP superar estado-da-arte em linhagens não vistas; o ganho
vem do **alinhamento**, não da feature. XLearning Lab (2º lugar) usou embeddings ESM-2 do
gene perturbado; BM_xTVC (1º) usou embeddings de proteína; o Arc registra que features
multimodais, sobretudo embeddings de proteína, agregaram valor nas equipes de topo.
**Contra-evidência de vocês, que eu não ignoro:** a medição de que as cinco knowledge
sources do PIE rendem 0,010 de AUPRC total para os 143 sem evidência. A diferença entre
isso e a reivindicação da preprint é exatamente **alinhamento à geometria da resposta
versus feature concatenada** — se o alinhamento não for implementado, o ganho medido por
vocês é o que se deve esperar.
**Esforço:** o mais alto dos cinco, e o de retorno mais incerto. Último por isso.

---

## 7. Lacunas declaradas

- **Pontuação do `baseline.py` de 2026:** não encontrei fonte.
- **Valores por métrica (des/pds/mae) do placar final de 2025:** não recuperados; só o
  `overall`, e esse pelo índice de busca.
- **Método completo do campeão de 2025 (xTrimoSCPerturb):** artigo prometido, não
  publicado em 2026-10-22.
- **Números de PDS/DES/MAE do PRiMeFlow no desafio:** não constam no preprint.
- **Write-up técnico "Mean Predictors" da Turbine:** visualizador JS, conteúdo não
  extraído.
- **Números de Mejia et al. (WMSE):** só o resumo foi lido; os valores ficam pendentes.
- **Tabelas da tese de Kernfeld:** só o resumo foi lido.
- **TxPert:** resultados quantitativos não extraídos; nada afirmado sobre eles.
- **Acervo local de PDFs:** `library_search` falhou com erro de execução; não consultado.
- **Estudo que resolva explicitamente a tensão MSE × conjunto de DE numa submissão
  única:** não existe, pelo que encontrei.

---

## 8. Referências

Ordenadas por tipo e, dentro de cada tipo, por relevância para este documento.

### Revisado por pares

1. Roohani YH, Hua TJ, Tung P-Y, *et al.* **Virtual Cell Challenge: Toward a Turing test
   for the virtual cell.** *Cell* 2025;188(13):3370–3374.
   [doi:10.1016/j.cell.2025.06.008](https://doi.org/10.1016/j.cell.2025.06.008) ·
   PMID 40578317
2. Arc Virtual Cell Initiative Team, Goodarzi H. **Virtual Cell Challenge 2026:
   Benchmarking zero-shot generalization across cellular contexts.** *Cell* 2026
   (2026-08-26). PMID 42648290
3. Adduri AK, Gautam D, Bevilacqua B, *et al.* **Predicting cellular responses to
   perturbation across diverse contexts with State.** *Cell* 2026;189(19):5914–5931.e20.
   [doi:10.1016/j.cell.2026.07.052](https://doi.org/10.1016/j.cell.2026.07.052) ·
   PMID 42673963
4. Ahlmann-Eltze C, Huber W, Anders S. **Deep-learning-based gene perturbation effect
   prediction does not yet outperform simple linear baselines.** *Nature Methods*
   2025;22:1657–1661.
   [doi:10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6)
5. Li L, You Y, Liao W, *et al.* **A systematic comparison of single-cell perturbation
   response prediction models.** *Science Advances* 2026.
   [doi:10.1126/sciadv.aed3414](https://doi.org/10.1126/sciadv.aed3414)
6. Csendes G, *et al.* **Benchmarking foundation cell models for post-perturbation RNA-seq
   prediction.** *BMC Genomics* 2025.
   [doi:10.1186/s12864-025-11600-2](https://doi.org/10.1186/s12864-025-11600-2)

### Preprints (não revisados por pares)

7. Molina A, Zhang X. **Perturbation response decomposition enables biologically aligned
   generalization to unseen perturbations and cellular contexts.** bioRxiv 2026-07-27.
   [doi:10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459)
8. Qi S-a, Chapfuwa P. **Mechanisms Matter: Transportability of Cellular Perturbation
   Effects.** bioRxiv 2026-05-12.
   [doi:10.64898/2026.05.08.723625](https://doi.org/10.64898/2026.05.08.723625) ·
   código: `github.com/microsoft/mechanisms_matter`
9. Liu Q, Zhang Q, Du J-H, Zhao S, Wang J. **Effects of Distance Metrics and Scaling on
   the Perturbation Discrimination Score.** arXiv:2511.16954, 2025-11-21.
   [arxiv.org/abs/2511.16954](https://arxiv.org/abs/2511.16954)
10. Yan Z, Wu Y, Ji MX, *et al.* **PRiMeFlow: Capturing Complex Expression Heterogeneity
    in Perturbation Response Modelling.** arXiv:2604.13986.
    [arxiv.org/abs/2604.13986](https://arxiv.org/abs/2604.13986)
11. Mejia GM, Miller HE, Leblanc FJA, Wang B, Swain B, de Lima Camillo LP. **Diversity by
    Design: Addressing Mode Collapse Improves scRNA-seq Perturbation Modeling on
    Well-Calibrated Metrics.** arXiv:2506.22641, 2025-06-27.
    [arxiv.org/abs/2506.22641](https://arxiv.org/abs/2506.22641)
12. Cui P, Liu Y, Sun W. **PerturbMap: Cross-Context Transfer of Single-Cell Perturbation
    Responses.** arXiv:2607.28090, 2026-07-30.
    [arxiv.org/abs/2607.28090](https://arxiv.org/abs/2607.28090)
13. Schäfer PSL, Reid K, Boldyga Z, Aksu ED, Hakem H, Saez-Rodriguez J. **Towards
    Principled Evaluation of Single-Cell Perturbation Prediction Models.** bioRxiv
    2026-07-27.
    [doi:10.64898/2026.07.23.740433](https://doi.org/10.64898/2026.07.23.740433)
14. Wenkel F, Tu W, Masschelein C, *et al.* **TxPert: Leveraging Biochemical Relationships
    for Out-of-Distribution Transcriptomic Perturbation Prediction.** arXiv:2505.14919,
    2025-05-20. [arxiv.org/abs/2505.14919](https://arxiv.org/abs/2505.14919)
15. Jiang W, Liu Y, Cai Y, *et al.* **What Makes a Representation Good for Single-Cell
    Perturbation Prediction?** arXiv:2605.19343, 2026-05-19.
    [arxiv.org/abs/2605.19343](https://arxiv.org/abs/2605.19343)

### Tese

16. Kernfeld E. **Assessing Reliability of Causal Models of Transcription.** Tese de
    doutorado, Johns Hopkins University.
    [jscholarship.library.jhu.edu](https://jscholarship.library.jhu.edu/server/api/core/bitstreams/d18fbb3c-9703-4296-b3e8-e930324fd5cf/content)

### Documentação oficial do desafio

17. Arc Institute. **Virtual Cell Challenge 2025 Wrap-Up: Winners and Reflections.**
    [arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up)
    — fonte primária única mais importante deste documento.
18. Arc Institute. **The 2026 Virtual Cell Challenge: predicting perturbation responses in
    cell contexts a model has never seen.**
    [arcinstitute.org/news/virtual-cell-challenge-2026](https://arcinstitute.org/news/virtual-cell-challenge-2026)
19. Arc Institute. **Generating high-quality Perturb-seq data for the 2026 Virtual Cell
    Challenge.**
    [arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026)
    — seis linhagens, >10.000 perturbações perfiladas por linhagem, curadas para ≥80% de
    knockdown mediano on-target, mediana de 20.000 UMI/célula, >800 bilhões de reads por
    linhagem na Ultima UG100, ~50.000 células de controle NTC por experimento.
20. ArcInstitute/cell-eval (GitHub).
    [github.com/ArcInstitute/cell-eval](https://github.com/ArcInstitute/cell-eval)
21. Placar final de 2025. [virtualcellchallenge.org/2025](https://virtualcellchallenge.org/2025)
    — página renderizada por JS; valores lidos pelo índice de busca.

### Relatos de equipe e comentário (não revisados, úteis como relato primário de método)

22. Turbine AI. **How did we get a regression model to the top of the Virtual Cell
    Challenge.** 2025-12-09.
    [blog.turbine.ai/p/how-did-we-get-a-regression-model](https://blog.turbine.ai/p/how-did-we-get-a-regression-model)
23. Dall'Olio GM. **Arc Virtual Cell Challenge — part 4: Winners announced!**
    [gmdbioinformatics.substack.com](https://gmdbioinformatics.substack.com/p/arc-virtual-cell-challenge-part-4)
24. Kilroy N. **Arc Institute's Virtual Cell Challenge: Are We Modelling Cells or Just
    Their Fingerprints?** Kiin Bio Weekly, 2025-12-17.
    [newsletter.kiin.bio](https://newsletter.kiin.bio/p/arc-institutes-virtual-cell-challenge)
25. **Compute, Confidence, and a Bad Pivot: My Lessons from Arc Institute's Virtual Cell
    Challenge 2025.** 2025-12-17.
    [yraj.in](https://yraj.in/ai/virtual-cell/perturbation-modeling/single-cell/2025/12/17/virtual-cell-challenge.html)
26. Pebblous. **Virtual Cell Challenge — When Hybrid AI and Curated Data Beat Pure
    Scale.**
    [blog.pebblous.ai](https://blog.pebblous.ai/report/nvidia-virtual-cell-challenge-2026-05/en/)
    — fonte secundária; único lugar com números do baseline State de 2025, **não
    confirmados**.
