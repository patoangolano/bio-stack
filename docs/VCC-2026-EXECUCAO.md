<!-- Execução das etapas 1 a 4 de docs/VCC-2026-ROTA.md, medida em 2026-10-22. -->

# VCC 2026 — execução das etapas 1 a 4

**Data:** 2026-10-22 · **Prazo do desafio:** 2026-11-05T23:59Z

O memorando `VCC-2026-ROTA.md` é o registro de decisão. Este é o registro do
que foi construído e, principalmente, do que o dado contradisse. Onde os dois
divergem, vale este.

---

## 1. O achado que reorganiza a estratégia: profundidade por perturbação

Medido sobre as primeiras 600 linhas de cada repo (`fdr < 0,05` entre os genes
testados):

| dataset | genes significativos / pert (mediana) | média | p90 | `\|delta_p\|` mediano | nossos alvos no split de treino |
|---|---|---|---|---|---|
| Orion | **2** | 7,9 | 11 | 0,0054 | **148** |
| Replogle | 46 | 428 | 1.464 | 0,0389 | 0 |
| **arc_vcc_25** | **2.016** | **2.754** | 6.958 | 0,0068 | 12 |
| Jiang | 13 | 243 | 962 | 0,0046 | 6 |

A referência real do desafio tem ~2.739 genes significativos por perturbação.
**Só o `arc_vcc_25` está nesse regime** — é o dado da própria competição, mesmo
ensaio e mesma profundidade.

Orion, que carrega 148 dos 157 alvos com evidência, tem mediana de **2**. A
consequência não é que Orion seja inútil; é que o plano original, que o tratava
como cavalo de batalha por causa da sobreposição de alvos, escolheu o doador
pelo critério errado.

## 2. A qualidade do alvo é o teto; a do doador modula dentro dele

Vinte e quatro direções de transferência medidas, em 11 pares de contexto
(dentro e **entre** repos — pares entre repos não estavam no plano e o dado
exigiu).

Contraste com **doador fixo** em `orion.hct116`:

| alvo | mse/mse₀ | cosseno |
|---|---|---|
| replogle.rpe1 | **0,9540** | +0,178 |
| replogle.hepg2 | 0,9552 | +0,164 |
| replogle.k562 | 0,9702 | +0,137 |
| replogle.jurkat | 0,9808 | +0,112 |
| arc_vcc_25.ARC_H1 (n=67) | 0,9948 | +0,070 |
| orion.hek293t | **0,9984** | +0,034 |

Contraste com **alvo fixo** em `orion.hct116`: o melhor doador não passa de
0,9624. Nenhuma escolha de doador salva um alvo raso.

A leitura: `mse/mse₀` mede fração da variância do alvo explicada. Alvo ruidoso
tem `mse₀` inflado por ruído, então a razão tende a 1 por construção. O nosso
alvo é o ensaio da competição — profundo. Logo a direção análoga à nossa é
**doador raso → alvo profundo**, cujo melhor valor medido é **0,954**, e não o
0,9984 de Orion→Orion.

## 3. O encolhimento vale mais que o modelo

Em **todos** os 11 pares, nas duas direções:

```
beta_por_gene  >  beta_global  >  zero (resposta média)  >>  copia_doador (beta=1)
```

`copia_doador` pontua de **1,2 a 2,8 vezes pior que não predizer nada**. Como o
zero da escala do `vcc2026` é a resposta média real do contexto, copiar o doador
cru colocaria a submissão abaixo de zero. A magnitude do encolhimento é a
decisão de engenharia mais importante do pipeline, acima de qualquer escolha de
arquitetura.

Resumo das 24 direções com `beta_por_gene`:

- **24 de 24 batem o zero.** Consistente, e por margem pequena.
- `mse/mse₀`: mediana 0,9661 · melhor 0,9040 (`replogle.hepg2→rpe1`) · pior 0,9986
- cosseno: mediana 0,153 · melhor 0,265 · pior 0,032

## 4. Proximidade de linhagem não prediz transferência — reproduzido

O memorando já marcava isso como refutado, com base no `test_unseen` do PIE.
Medição independente, nossa, concorda: o melhor par de todos é
**`hepg2 → rpe1`** (0,9040), e RPE1 é a linhagem mais divergente do conjunto.
Dentro do Replogle, os quatro melhores alvos são RPE1 e HepG2.

E o pior par de todos é `orion.hct116 → orion.hek293t`, duas linhagens do
**mesmo estudo e mesmo ensaio**, com a maior amostra disponível (3.277
perturbações). Não é distância de contexto; é profundidade.

## 5. Correções de semântica dos arrays PIE

Tudo medido, não lido em documentação.

| afirmação | o que o dado diz |
|---|---|
| `fold_changes` é log2FC | **Razão linear.** Mediana entre os testados: 0,9974 a 1,0006 nos quatro repos |
| zero em `fold_changes` é razão zero | **É ausência de medida.** `fold_changes == 0` ⟺ `tested == False`, por igualdade de máscara |
| `delta_p` é discretizado em bins | **Contínuo.** 197.588 valores únicos em 200 mil amostras, passo na resolução do float32. O passo 0,00979 é da cabeça de predição, não do alvo |
| `ctrl_means` é log1p do pseudobulk | **Média de log1p por célula.** `sum(expm1(ctrl_means))` = 8.615, não 10.000 — desigualdade de Jensen |
| `arc_vcc_25` tem um contexto | **Dois:** `ARC_H1` (146 perts) e `ARC_H1_VAL` (44), com **zero** perturbações em comum |

A última fecha uma porta: `ARC_H1 → ARC_H1_VAL` seria a mesma linhagem e daria
o **teto de reprodutibilidade** — quanto se pode acertar predizendo uma medida
nova da mesma perturbação no mesmo contexto. Com conjuntos disjuntos, esse teto
não é mensurável com o dado disponível.

Consequência de desenho: a variável transferida é `log2(fold_changes)`
mascarado por `tested`, não um delta em log. A razão é adimensional e se aplica
ao **nosso** `b_ctrl`, o que cancela a divergência de estimador em vez de
propagá-la. Também dispensa `ctrl_means`, que o `arc_vcc_25` não publica.

## 6. Universo de DE: a definição medida

Duas definições de "CPM médio no controle acima de 5" circulam e dão universos
diferentes. Os números oficiais são D 9.633 · E 9.743 · F 9.537.

| definição | D | E | F | reproduz |
|---|---|---|---|---|
| CPM do pseudobulk | 9.565 | 9.744 | 9.499 | não |
| **média do CPM por célula** | **9.633** | **9.743** | **9.537** | **sim, nos três** |

Diferença de 68 genes em D e 38 em F. Cada um deles seria predito-sem-teste ou
testado-sem-predição, o que bate direto em `jaccard` e `reach`.

E o universo é muito mais específico de contexto do que parece: só **7.772
genes** estão no universo dos três contextos, contra **11.593** em pelo menos
um. Um terço da união depende de contexto.

## 7. Evidência doadora: 157, não 161

| contagem | valor |
|---|---|
| alvos medidos em algum repo (`pert_to_id`) | 161 |
| alvos no **split de treino** (`PIE_splits/replogle_xdataset/train.json`) | **157** |
| por repo no treino | orion 148 · arc_vcc_25 12 · jiang 6 · replogle 0 |
| medidos mas reservados para validação | 4: **FUBP1, HMGXB4, TIGD5, ZBTB7B** |
| sem nenhuma evidência | 143 |

Os per-repo batem com o memorando. Treinar nos 4 reservados seria vazamento;
eles entram só como checagem.

## 8. Gate de expressão

CRISPRi reprime transcrição: alvo que o contexto não transcreve não tem o que
derrubar. Com corte em 2% das células:

| contexto | alvos ativos | sem alvo transcrito |
|---|---|---|
| D | 298 | 2 |
| E | 299 | 1 |
| F | 299 | 1 |

Menos do que se poderia temer — os 300 alvos são quase todos transcritos nos
três contextos.

## 9. Imputação dos 271 genes fora do `gene_text`

Os 18.533 genes do painel contra as 31.684 chaves de `gene_text`: **18.262
cobertos, 271 ausentes** — o número do memorando, reproduzido.

**Decisão, tomada antes e não depois:** os 271 **permanecem** no eixo de saída,
porque a submissão exige os 18.533. Perdem apenas a *feature* de texto,
imputada pela média da coluna. Não são gateados a zero: expressão de controle
existe para todos os 18.533, e a evidência cruzada de outro dataset não depende
de `gene_text`. Alias não é aplicado aqui porque o `datamodule` da Arc também
não aplica — o número de cobertura tem de ser o mesmo que o PIE veria.

As outras quatro fontes, indexadas por perturbação, cobrem: `esm2` 300/300,
`depmap_gene_effect` 300/300, `string_space` 300/300, `perturbation_text`
161/300.

---

## 10. O que foi construído

| artefato | etapa | o que faz |
|---|---|---|
| `scripts/vcc/baixar_pie_arrays.py` | — | baixa 4,47 GB de arrays PIE, retomável. Exclui `lfc_true` (float64, 2,36 GB só no Orion) e Tahoe (droga-dose) |
| `scripts/vcc/eixo_18262.py` | 1 | eixo canônico, cobertura das 5 fontes, mapa de índice para cada repo, gate de expressão, evidência ciente do split |
| `scripts/vcc/ctrl_bulk.py` | 2 | `b_ctrl`, erro padrão por jackknife delete-1 em forma fechada, universo de DE nas duas definições, distribuição de biblioteca |
| `scripts/vcc/transfer_fit.py` | 3 e 4 | `validar` mede 24 direções de transferência; `ajustar` aplica a política a D/E/F |

Saídas em `runs/vcc/eixo/`, `runs/vcc/ctrl_bulk/`, `runs/vcc/transfer/`.

O jackknife delete-1 tem forma fechada sobre uma soma — `P` sem a célula *i* é
`P − x_i` e o total é `T − t_i` —, então não há 18.400 recálculos: é um produto
externo por bloco de genes, corrigido nas entradas não nulas. 37 s para os três
contextos.

**Formato de `runs/vcc/transfer/delta_previsto.npz`**, que o `emissor_400`
consome:

```
genes    (18533,)  object  símbolos na ordem de gene_names.csv
alvos    (300,)    object  símbolos na ordem de pert_counts.csv
delta_D  (300, 18533) float32   delta em espaço log1p de CP50k
delta_E  (300, 18533) float32
delta_F  (300, 18533) float32
```

Convenção: `b_previsto = b_ctrl + delta`. Delta zero significa predizer o
controle. **Atenção:** predizer o controle pontua *abaixo* de zero, não zero —
ver a seção 15, que corrige um erro desta mesma frase na primeira versão deste
documento. Todos os 300 alvos carregam a resposta média estimada do contexto.

## 11. Conversão de espaço, resolvida

O memorando listava como questão aberta "como mapear `delta_p_pred` para o que
cada métrica espera". A conversão é fechada:

- `b = log1p(5e4 · p)`, logo `expm1(b) = 5e4 · p`
- a proporção perturbada é `p' ∝ p · 2^log2r`, **renormalizada para somar 1**
- `delta = log1p(expm1(b_ctrl) · 2^log2r · renorm) − b_ctrl`

A renormalização não é detalhe: sem ela, muitos genes se movendo deslocam todo
o eixo. É composicional, e esquecê-la enviesa as métricas de magnitude.

## 12. Baixa dimensão: o modelo diagonal não era suficiente

O modelo da seção 2 é **diagonal** — um coeficiente por gene, nenhuma interação.
Respostas de Perturb-seq não ocupam 18.533 dimensões independentes. Adicionado
um segundo modelo: SVD randomizada da resposta do doador, `k = 32` programas,
um coeficiente γ por programa em vez de um β por gene (32 parâmetros contra
18.533). A base vem **só do doador** — usar a do alvo seria vazamento e, para
D/E/F, impossível.

Escada completa, 24 direções:

| modelo | mse/mse₀ mediano | melhor | cosseno mediano |
|---|---|---|---|
| zero (resposta média) | 1,0000 | 1,0000 | — |
| copia doador (β=1) | 1,6203 | 1,1960 | 0,112 |
| β global | 0,9830 | 0,9473 | 0,112 |
| β por gene encolhido | 0,9661 | 0,9040 | 0,153 |
| γ por programa (k=32) | 0,9624 | 0,8748 | 0,153 |
| **média dos dois** | **0,9552** | **0,8721** | **0,178** |

Vencedor por direção, em 24: média dos dois **20**, γ por programa 3, β global 1.

Os dois modelos são **complementares, não redundantes** — é o que a média
vencendo em 20 de 24 demonstra. A diagonal acerta gene a gene, a de baixa
dimensão acerta a estrutura. A melhora sobre o zero sai de 3,4% para 4,5% e o
cosseno de 0,153 para 0,178.

Um detalhe que merece registro: o β por gene **destrói** o ordenamento de
magnitude (spearman de |efeito| cai de 0,112 no β global para 0,038), enquanto
melhora cosseno e MSE. A baixa dimensão recupera parte (0,078). Não comparar
esse número com os 0,741 do `Rhaister-O` citados no memorando: aquele é
calculado entre perturbações, este é entre genes dentro de uma perturbação. São
quantidades diferentes.

Na aplicação a D/E/F a parte de baixa dimensão vale só nos **6.177 genes** onde
a base foi medida; fora dela, a diagonal responde inteira.

## 13. Os 143 alvos sem evidência: resultado negativo, medido

`esm2`, `depmap_gene_effect` e `string_space` cobrem **300 de 300** alvos. Se
as cargas de programa fossem previsíveis a partir delas, os 143 alvos sem
evidência sairiam do piso — dobrando os alvos que contribuem.

Regra de decisão, **fixada antes de rodar**: os 143 entram apenas se o cosseno
fora da amostra superar 0,05 em validação cruzada de 5 dobras.

| contexto | melhor α | cosseno fora da amostra | R² | veredito |
|---|---|---|---|---|
| D | 1.000 | +0,0092 | −0,231 | reprovado |
| E | 1.000 | +0,0249 | −0,199 | reprovado |
| F | 1.000 | +0,0027 | −0,236 | reprovado |

O R² é **negativo em todos os seis valores de α e nos três contextos**: a
predição por ridge é pior que prever a carga média. Não é um pico estreito mal
calibrado, é ausência de sinal — e concorda com o 0,010 de AUPRC que o próprio
artigo do PIE mediu para as cinco fontes em regime de contexto e perturbação
ambos novos.

**Os 143 não recebem desvio específico de perturbação.** Recebem a resposta
média do contexto (seção 15), que é o zero da escala — e não o controle, que
pontuaria abaixo de zero.

Ressalva honesta: o rótulo de treino vem, em 148 dos 157 casos, do Orion, com
mediana de dois genes significativos por perturbação. Um mapa perfeito de
features para cargas ainda mostraria correlação baixa contra rótulo tão
ruidoso. O resultado negativo é sobre *esta* combinação de features e rótulo,
não uma prova de que as fontes sejam inúteis em princípio. Mas é a combinação
que temos.

## 14. Expectativa honesta, e o que ela não autoriza

Para os 157 alvos com evidência, o análogo medido mais próximo é
`orion.hct116 → replogle.rpe1`, agora com o modelo combinado: **mse/mse₀ ≈
0,955 mediano**, cosseno ≈ 0,18. Para os 143 sem evidência, a predição é zero,
que pontua zero.

Ponderado pelos 157 de 300, isso é da ordem de **2 a 3% de melhora de MSE**
sobre a resposta média — pequeno, positivo e medido. Não é um resultado que
autorize otimismo sobre colocação; é um resultado que autoriza submeter algo
estritamente melhor que o trivial, e saber por quanto.

Das três frentes que a primeira versão deste documento listou, duas foram
resolvidas no mesmo dia:

1. ~~Estrutura de baixa dimensão~~ — **feita** (seção 12). Ganhou 1,1 ponto
   percentual de MSE e 0,025 de cosseno.
2. ~~Os 143 sem evidência~~ — **testada e reprovada** (seção 13). Fecha a
   linha: não há sinal nas fontes de conhecimento contra este rótulo.
3. **Um avaliador local funcionando** (etapa 5) — segue sendo o gargalo. Sem
   ele, nenhuma mudança é verificável antes de gastar uma das 2 submissões
   diárias, com placar oculto até 05/11.

O que resta com retorno plausível, em ordem:

- **Rótulo melhor para os 157.** O limite não é o modelo, é que 148 dos 157
  vêm de um ensaio com mediana de 2 genes significativos. Os parquet de `de/`
  do Orion (3,4 GB para hct116) trazem a estatística por perturbação em tabela
  longa e podem permitir uma reestimativa menos moderada que o `fold_changes`.
- **γ por gene e por programa ao mesmo tempo**, em vez da média 50/50 fixa. O
  peso ótimo é estimável nos pares medidos, e hoje é um palpite redondo.
- **Pesar os programas pela sua própria transferibilidade.** Hoje γ é estimado
  por componente mas a base é a do doador; uma base compartilhada entre todos
  os doadores poderia transferir melhor.
- **Os 4 alvos reservados são conservadorismo meu, não necessidade.** FUBP1,
  HMGXB4, TIGD5 e ZBTB7B estão no split de *validação* do PIE e eu os excluí da
  evidência doadora. Mas aquele split existe para que as métricas reportadas do
  PIE sejam honestas — não tem relação com a competição. Usar a resposta medida
  deles em `ARC_H1_VAL` para predizer D/E/F é dado público, não vazamento: os
  contextos nem se cruzam. Incluí-los daria 4 alvos a mais com desvio
  específico, 161 em vez de 157. Ganho de 1,3%, custo de uma reemissão de 400 s
  — fica registrado como escolha, não como regra.

---

Vinte e quatro direções medidas, vinte e quatro acima do zero, e a maior
descoberta foi que o doador com 148 dos nossos alvos tem mediana de dois genes
significativos por perturbação.

---

## 15. A correção que vale mais que todo o resto: a resposta média do contexto

A primeira versão deste documento continha um erro meu, e ele atingia 143 dos
300 alvos por inteiro.

**O erro.** A escala do desafio é `s = (u − b) / (r − b)`, e `b` é o valor da
métrica alcançado por uma predição que devolve, para toda perturbação, a
**média real das 300 respostas daquele contexto**. Predizer essa média pontua
exatamente zero. Predizer o **controle** pontua *abaixo* de zero, porque o
controle está mais longe de cada resposta real do que a média delas está — a
média carrega o componente comum a todas as 300.

O `transfer_fit` prediz o desvio **centrado** de cada perturbação, centrado na
média do contexto doador, que era subtraída e descartada. Na emissão,
`delta = 0` virou "prediga o controle". Para os 143 alvos sem evidência isso
era a predição inteira, e eu havia escrito que "pontua zero". Pontua negativo.

**O que o dado diz.** A resposta média é mensurável nas 24 direções do
`validar`, porque `media_por_gene` de doador e alvo ficam gravadas em cada
`beta_*.npz`. `scripts/vcc/resposta_media.py` mede:

| | valor |
|---|---|
| correlação da resposta média entre contextos | mediana **+0,470** (min +0,107, max +0,839) |
| coeficiente de transferência | mediana **+0,408** |
| ganho sobre predizer o controle | mediana **+0,252** |
| direções com ganho positivo | **24 de 24** |
| magnitude típica | 0,084 em log2 |

Ganho de 0,252 do erro quadrático. Compare com o desvio específico de
perturbação, que rende 0,045 (mse/mse₀ 0,955). **A parcela que eu estava
jogando fora é cerca de cinco vezes maior que a que eu estava calculando** — e
vale para os 300 alvos, não só para os 157.

**Uma armadilha dentro da correção.** A primeira medição da resposta média deu
magnitude RMS 0,88, contra mediana de 0,033. A distribuição é de cauda pesada:
`tested` cobre 51% dos genes no Orion e 87% no arc_vcc_25, então um gene
testado numa **única** perturbação herdava aquele log2 como se fosse a média do
contexto, e a cauda batia no limite de ±5. O Replogle, que testa 100% dos
genes, não tinha o problema (RMS 0,18, máximo 0,94) — o que identifica a causa.

A correção é encolher a média pela cobertura, `n / (n + 20)`, a mesma lógica do
β. E ela **melhorou** a medição em vez de piorar, que é a assinatura de ter
removido ruído e não sinal:

| | antes de encolher | depois |
|---|---|---|
| correlação mediana | +0,436 | **+0,470** |
| ganho sobre o controle | +0,164 | **+0,252** |
| magnitude típica | 0,135 | 0,084 |

**Como entra na predição.** A resposta média vale para **todos os 300**,
inclusive os gateados por expressão: é o efeito comum do procedimento de
perturbação — entrega de guia, seleção, atividade de Cas —, pelo qual todos os
300 passaram, inclusive quem não transcreve o próprio alvo. O desvio específico
é que depende de evidência e de o alvo estar transcrito.

O coeficiente 0,408 foi medido doador-único contra alvo-único. Agrupar doadores
reduz o ruído do preditor, então o ótimo para a média agrupada seria **maior**.
Usar o valor medido erra por baixo de propósito: a predição fica entre o
controle e a média real, que ainda bate o controle, em vez de passar da média.

**A verificação que convenceu.** Os genes que o deslocamento mais move em D:

```
DDIT3  +0,123   CHOP, resposta integrada a estresse
BTG2   +0,125   antiproliferativo responsivo a p53
HSPA6  +0,171   choque térmico
DUSP10 +0,120   fosfatase de MAPK
IRF8   +0,242   fator mielóide mestre  (D é mielóide)
INPP5D +0,128   SHIP1, mielóide
SIRPB1 +0,145   mielóide
CXCR4  +0,136   quimiocina
```

`DDIT3` e `BTG2` juntos são a assinatura canônica de resposta genérica a
perturbação, e `IRF8`/`INPP5D`/`SIRPB1` são certos para um contexto mielóide.
Todos os 18 genes mais movidos estão no universo de DE de D, e todos positivos.
Isso é biologia, não artefato numérico.

**Decomposição final.** Nos 157 alvos com evidência, o desvio específico
responde por 84% da norma do delta (p10 0,72, p90 0,93) e a resposta média
pelo resto. Nos 143 sem evidência, o delta é só a resposta média, e é idêntico
entre eles — verificado.

**Por que isto está documentado como erro e não corrigido em silêncio.** Porque
a frase errada já estava escrita aqui, e alguém poderia lê-la. E porque o
mecanismo do erro — centrar os dados, modelar o resíduo, e esquecer que a
constante subtraída era a referência da métrica — é o tipo de coisa que volta.

---

Vinte e quatro direções medidas, vinte e quatro acima do zero em ambas as
parcelas, e o maior ganho do dia veio de uma constante que eu havia subtraído e
jogado fora.

---

## 16. O gate duplo: o mesmo erro, reintroduzido na camada de baixo

O `emissor_400` tem um gate próprio, ligado por padrão, que **zera o delta**
nos pares cujo alvo não é transcrito no contexto. Lido isoladamente, faz
sentido: CRISPRi reprime transcrição, e sem transcrito não há o que derrubar.

Mas zerar o delta é predizer o **controle**, que pontua abaixo de zero. É o
erro da seção 15 exatamente, uma camada mais embaixo.

O gate pertence ao `transfer_fit`, onde ele suprime o **desvio específico** e
preserva a resposta média — que vale para todos os 300 porque é o efeito comum
do procedimento, pelo qual até quem não transcreve o alvo passou. No emissor,
ele apaga as duas parcelas.

Quatro pares foram atingidos em `submissao_v2.h5ad`: `D/SEMA4F`, `D/ZNF264`,
`E/ZNF383`, `F/RAB26` — 4 de 900, 0,44%. Pequeno, e na direção errada.

**Correção:** emitir com `--sem-gate`. O gate do `transfer_fit` já fez o
trabalho certo, e desligar o do emissor evita que ele desfaça.

A lição não é sobre os 4 pares. É que a mesma decisão — "este par não tem
resposta, zere" — aparece em duas camadas com significados diferentes, e a
camada de baixo não sabe que a de cima já tratou o caso. Quando duas camadas
podem gatear a mesma coisa, só uma deve ter a permissão.

---

## 17. Calibração de magnitude: o encolhimento não esmagou a amplitude

Encolhimento ótimo para MSE produz predição sistematicamente **sub-dispersa** —
é a regressão à média. Isso minimiza `mse` e `nmae` mas pode afundar as
métricas de DE (`jac`, `reach`, `fid`), que premiam acertar o *conjunto* de
genes significativos e portanto precisam de amplitude realista. E o memorando
diz que é DE que separa o top 10.

A tensão é real, então foi medida. Genes por alvo com `|delta|` acima de um
limiar, **dentro do universo de DE**:

| limiar | D | E | F |
|---|---|---|---|
| 0,001 | 8.756 | 8.786 | 8.488 |
| 0,005 | 5.224 | 5.310 | 4.911 |
| **0,010** | **2.640** | **2.672** | **2.459** |
| 0,020 | 832 | 890 | 746 |
| 0,050 | 141 | 154 | 142 |
| 0,100 | 24 | 25 | 23 |

A referência real tem **~2.739** genes significativos por perturbação. No
limiar de 0,010 a predição produz 2.459 a 2.672 — a mesma ordem. A
distribuição de magnitude está na faixa certa.

`|delta|` no universo: p50 0,0054–0,0062 · p90 0,0258 · p99 0,0827. (O valor de
0,0027 citado antes era sobre os 18.533 genes, incluindo os não expressos, e
subestima o que importa.)

**Se o limiar de detecção do Wilcoxon com 400 contra 18.400 células é da ordem
de 0,01, a contagem sai certa. Esse "se" é o que o `avaliador_proxy` existe
para resolver**, e é a razão de a etapa 5 ser o gargalo e não um item de
conveniência.

### Os extremos são os próprios alvos, e são negativos

Os 14 maiores `|delta|` em D: 7 são o par alvo→alvo, e **todos negativos** —
`HMGB2` −2,91, `ANP32E` −1,40, `PA2G4` −0,91, `ANP32B` −0,89, `NUDT5` −0,89,
`EMP3` −0,84, `ASPM` −0,80. Knockdown, que é o que CRISPRi faz. O gene-alvo é
excluído de todas as métricas, então isso não pontua — serve como verificação
de que o sinal não está invertido.

Os não-autólogos entre os extremos são coerentes: `SIN3B` sobe `MT-ND3`,
`MT-ATP6`, `MT-CO2` e `DDIT3` enquanto desce `PCLAF`; `MTA1` desce `NCAPD3`;
`STAT6` desce `CENPN`. Mitocondrial, estresse e ciclo celular.

**O que é fraco:** só **18 dos 300** alvos têm seu maior `|delta|` no próprio
gene. Em dado real de CRISPRi, o gene derrubado costuma ser o mais alterado.
Isso não custa pontuação — o alvo é excluído —, mas é mais uma medida da rasura
do Orion, de onde vem a evidência de 148 dos 157.

---

## 18. Etapa 5 — o `avaliador_proxy` roda, e o número que o limita

Ambiente em `~/envs/avaliador-proxy/.venv` no WSL: `cell-eval2` 0.19.0, `pdex`
0.3.0, `numba-mwu` 0.2.0, numpy 2.4.6, scipy 1.17.1, polars 1.44.2, numba
0.68.0, anndata 0.13.4, Python 3.12.3. Detalhe e armadilhas em
[`AVALIADOR-PROXY.md`](AVALIADOR-PROXY.md); números em
`runs/vcc/avaliador/verificacao.json`.

**Os oito campos de identidade batem.** `todos_os_oito_peers_batem: true`,
`rule_mismatches: []`. Essa era a coisa com maior chance de quebrar na hora da
submissão — comparação fatal, sem perdão — e está provada localmente. Das 22
variantes testadas, 7 dos 8 campos se movem com alguma; 1 não se move com
nenhuma.

**Memória, medida, e é ela que manda:**

| | valor |
|---|---|
| pico de RSS | **10,855 GiB** (no passo `baseline`) |
| células medidas | 31.200 (32 perturbações × 400 + 18.400 controles) |
| custo por mil células | 0,3303 GiB, mais 0,55 GiB constante |
| **perturbações que cabem** | **48**, com os 18.400 controles, na VM de 15 GiB |
| painel oficial inteiro | 138.400 células → **46,3 GiB projetados** |

Os 46,3 GiB confirmam o memorando por medição e não por estimativa: o painel
completo não cabe nem no host de 31,4 GiB. **Avaliar em subconjunto não é
otimização, é condição de existência.**

**Tempo, e onde ele está:**

| passo | wall | RSS |
|---|---|---|
| `proxy` | 16,0 s | 3,69 GiB |
| `baseline` | 187,5 s | 10,86 GiB |
| `prep-real-bundle` | **840,9 s** | 8,27 GiB |
| `run` | 191,5 s | 8,21 GiB |
| `score` | 2,0 s | 0,31 GiB |

Total 1.237,9 s. Mas `prep-real-bundle` é **reaproveitável enquanto o lado real
não muda**, então o ciclo de iteração de modelo é `run` + `score` ≈ **3,2
minutos**. Isso é um ciclo utilizável, e é o que faltava.

**O que a pontuação do proxy diz, e o que não diz.** A corrida usou efeitos
reais do `arc_vcc_25` como "predição" sobre as células de D — é quase-oráculo,
serve para provar a máquina, não para estimar a nossa nota:

```
de_wilcoxon_lfc_nmae                      0,285
expr_mse_unbiased_capped_norm             0,541
de_wilcoxon_direction_fidelity_yield_raw  0,230
de_wilcoxon_direction_reach_raw           0,681
de_wilcoxon_sig_jaccard                   0,089
pds_cosine                                1,000
avg_score (from_replicate)                0,529
```

Duas leituras com cuidado. `pds_cosine = 1,000` é **degenerado**: `pds` é um
ranque entre as perturbações do painel, e com 32 em vez de 300 ele infla — não
leia isso como sinal. E `avg_score` tem `from_baseline` nulo de propósito: é a
marca de uma submissão matriculada contra bundle de competição, e o
`cell-eval2` anula a média do baseline para que ninguém leia número não-oficial
sob o rótulo oficial.

O `jaccard` de 0,089 num quase-oráculo é a medida mais interessante da tabela:
bate a faixa que o memorando reporta para o top 10 da validação (`jac` escalado
≤ 0,09). Se até o oráculo fica ali, DE é difícil para todo mundo — o que
reforça que é DE que separa.


## 19. A curva de `k_desvio`, reproduzida numa cadeia inteiramente reconstruída

O código perdido no re-clone de 2026-10-23 voltou por replay dos transcripts, e
a primeira coisa que ele tinha de provar era que voltou **certo**. A varredura de
`k_desvio` é o teste que prova três coisas de uma vez e nenhuma isoladamente: o
`transfer_fit` recuperado produz o mesmo delta, o gerador de transporte novo
emite o que o `emissor_400` perdido emitia, e o avaliador remontado no WSL
pontua igual. Se qualquer um dos três divergisse, a curva não fecharia.

Seis braços contra um bundle real novo (`proxy-vcc-D-orion-k40-s0`, 40 alvos do
painel, orion/hct116 gerando o lado real, 34.400 células por braço):

```
membro        k=1     k=2     k=3    [k=4]    k=5     k=6
fidelity     0,062   0,452   0,781   0,834   0,823   0,802
jaccard      0,084   0,401   0,532   0,554   0,551   0,541
reach        0,343   0,620   0,715   0,756   0,779   0,790
lfc_nmae     0,101   0,238   0,351   0,425   0,439   0,401
expr_mse     0,232   0,385   0,447   0,395   0,179   0,000
pds_cosine   0,995   0,995   0,995   0,994   0,994   0,993
AVG_SCORE    0,303   0,515   0,637   0,660   0,628   0,588
registrado   0,275   0,460   0,590   0,647   0,627   0,587
```

Ótimo em k=4, confirmado. k=5 e k=6 batem o registro em **0,0005 e 0,0008**;
k=4 em 0,0127.

**O que a curva diz de novo, e não estava no registro: os seis membros têm
ótimos DIFERENTES.** `expr_mse` pica em k=3 e cai para 0,395 em k=4, 0,179 em
k=5 e 0 em k=6 — este último é corte de escala, então o `avg_score` de k=6 está
otimista, como o registro já avisava. `lfc_nmae` pica em k=5. `reach` sobe
monotonicamente até k=6. `fidelity` e `jaccard` picam em k=4.

Então **k=4 é um compromisso e não um ótimo compartilhado**, e a assimetria
aponta onde atacar: o único membro que PERDE indo de k=3 para k=4 é o
`expr_mse`, justamente um dos dois de erro que concentram 60% do que falta para
o oráculo. `pds_cosine` não se move em toda a varredura (0,995 → 0,993),
confirmando que é invariante a escala.

Duas coisas que o caminho ensinou. `arc_vcc_25` não serve como gerador do lado
real com rótulos do painel: tem 12 alvos do VCC e só 11 passam o gate de 200
genes moventes, de 32 pedidos; o orion tem 148. A profundidade rasa dele
(mediana de 2 genes significativos por perturbação) não importa aqui, porque ele
gera o lado REAL e não doa efeito à predição — o que o avaliador precisa é de
estrutura de DE com os rótulos certos, e é o número de alvos que manda. E script
com CRLF faz o bash do WSL abortar com `set: pipefail: invalid option name`,
porque ele lê `set -euo pipefail\r`; os `.sh` escritos no Windows saem com CRLF, e
isso está fixado em `.gitattributes`.

Custo medido, que é o que torna a varredura viável: `prep-real-bundle` custa
527 s e depende SÓ do lado real, então o `avaliador_varredura.sh` o constrói uma
vez e reusa; cada braço depois é recorte 7 s + run 107 s + score 1,3 s. Pico de
RSS: baseline 11,7 GiB, bundle 8,8 GiB, run 8,8 GiB — não cabe nos ~10 GiB
livres do host Windows e cabe nos 14 GiB dedicados da VM do WSL. É por isso que
este ambiente mora lá, e é a única peça do projeto de VCC em que o WSL ganha: a
emissão da submissão roda a 34 MB/s num processo só, limitada por I/O de disco,
igual dos dois lados.

---

## 20. O encolhimento POR GENE: o primeiro ganho em que nenhum membro perde

A seção 19 deixou o diagnóstico pronto. Um escalar tem **um** grau de liberdade
e a discordância entre os seis membros tem **dois**: os membros de direção querem
que o gene sub-movido ande, e `expr_mse` quer que o gene bem determinado não seja
esticado. Mas isso não é desacordo sobre o NÍVEL — é sobre QUAIS GENES. São
conjuntos diferentes, e um fator por gene os separa.

### A medida: atenuação por gene contra verdade retida

`scripts/vcc/encolhimento.py` mede, em cada direção doador → alvo:

```
a_g = sd_p(previsto[.,g]) / sd_p(real[.,g])        k_g = 1 / a_g
```

com `previsto` reproduzindo o preditor de produção (0,5 diagonal + 0,5
programas, ver `transfer_fit.ajustar`). Só direção com **alvo profundo** entra —
alvo em (replogle, arc_vcc_25) com ≥ 400 perturbações, o mesmo gate do
`politica_de_beta` e pela mesma razão: a variância real do alvo é o denominador,
e num alvo raso ela é dominada por ruído de amostragem, o que a infla e
deflaciona `k_g`. Mediria a profundidade do ensaio e não a atenuação do preditor.
Sobram 16 direções.

**Quatro medidas deram direito de aplicar:**

1. **A atenuação GLOBAL cai em 3,40 a 8,69, mediana 5,615.** Isso corrobora o
   `k_desvio = 4` por um instrumento independente — o escalar veio de varredura
   no avaliador, este número vem de verdade retida — e explica a diferença sem
   apelar a erro: casar variância pede 5 a 6, que é exatamente onde `reach` (6) e
   `lfc_nmae` (5) picam na tabela da seção 19, e `expr_mse` pede 3 porque erro
   quadrático é minimizado pela média a posteriori, que é ENCOLHIDA e não
   calibrada. Os dois instrumentos concordam; o que difere é o funcional.

2. **O fator por gene varia 6 a 7 vezes** entre p05 e p95 numa direção (0,5× a
   3,5× o global). Um escalar joga tudo isso fora.

3. **E transfere.** Spearman de `k_g` em **36 pares TOTALMENTE DISJUNTOS** — sem
   doador em comum, sem alvo em comum — dá mediana **+0,543**, faixa +0,425 a
   +0,683. Com um lado em comum sobe a +0,629. Par disjunto é a única evidência
   de que a forma é propriedade do GENE: num par que compartilha um lado, a
   correlação poderia vir do lado compartilhado.

4. **Não é artefato de cobertura.** `rho(k_g, n_válido) = -0,092` contra
   `rho(k_g, var_prevista) = -0,806`. O que dirige o fator é quanto o preditor
   MOVE o gene: ele concentra movimento em poucos genes e sub-move o resto. Se o
   driver fosse a contagem de entradas válidas, o fator seria sobre a qualidade
   da medida e não sobre a predição, e não teria por que transferir.

### Duas verificações que o resultado exigia antes de ser legível

**O peso não é `1/sd_prev` disfarçado.** Se fosse, aplicá-lo igualaria a
magnitude de todo gene e destruiria a estrutura de que `pds_cosine` e `expr_mse`
vivem. Medido contra o `sd_prev` da predição em D/E/F: `r(log w, log sd_prev) =
-0,44`, R² **0,19**, e inclinação **-0,26** — a tautologia "iguale todo gene"
seria -1,0. `sd(log w) = 0,368`, dos quais **0,332 sobram** depois de remover a
componente de `sd_prev`. E a ordem de magnitude por gene sobrevive:
`spearman(sd_prev, sd_aplicado)` = +0,98 em alpha=0,5 e +0,91 em alpha=1.

Por que a correlação dentro de uma direção era -0,806 e aqui é -0,44: dentro de
uma direção, `k_g = sd_real/sd_prev` tem acoplamento mecânico com o próprio
`sd_prev` daquela direção. Agrupando 16 direções e medindo contra o `sd_prev` de
OUTRO preditor, a parte mecânica cancela e sobra a propriedade do gene — o que o
rho de +0,543 nos pares disjuntos já implicava.

**Não há vazamento do orion.** O proxy gera a verdade a partir do orion, e 4 das
16 direções têm orion como doador. Se o peso codificasse a magnitude por gene do
orion, ele encolheria justamente os genes que a verdade do proxy move, e o
veredito viria pessimista por construção. Tirar as quatro direções move o peso em
Spearman **+0,9936** (r em log +0,9939), com a razão `peso_tudo/peso_sem` entre
0,938 e 1,074 nos percentis 5–95, e o k global de 5,615 para 5,538. As doze
direções replogle→replogle dominam porque carregam ~1.000 perturbações cada
contra 480–608 do orion.

**Limite declarado.** 6.290 dos 18.533 genes do painel têm medida — 63,6% do
universo de DE de D e **53,8% da variância do delta**. O resto fica em peso 1,0 e
herda o escalar, não é imputado. Então `alpha = 1` remodela **metade da massa**,
e não "toda a predição reescalada por gene".

### A aplicação: um eixo só, com o escalar como caso especial exato

`varredura_k.py --alphas` emite `k_g = k_desvio * c(alpha) * peso_g**alpha`.

`c(alpha)` resolve `sum (c * peso**alpha * d)² == sum d²`, em forma fechada. Sem
ele, subir `alpha` aumentaria a magnitude TOTAL do desvio ao mesmo tempo que a
redistribui — a cauda direita do peso pesa mais na soma de quadrados que a
esquerda — e a varredura mediria a soma dos dois efeitos sem poder separá-los.
Com `c`, `alpha` é redistribuição pura e o nível continua sendo `k_desvio`.

Em `alpha = 0` isso dá `c = 1` e `fator = 1` exatamente, então o braço reproduz o
escalar **bit a bit** — o mesmo papel que `k=1` reproduzir o `delta_previsto` tem
na varredura escalar. Conferido nos dois lugares: por `np.array_equal` no teste, e
na emissão real, onde `delta_kd4_a0.npz` produz 682.368.320 não nulos, o mesmo
número do `kd4` da seção 19.

### O resultado

```
membro        a=0    a=0,25   a=0,5   a=0,75   a=1     a=1,5
fidelity     0,8340  0,8350  0,8361  0,8362  0,8353  0,8299
jaccard      0,5542  0,5585  0,5622  0,5634  0,5640  0,5581
reach        0,7557  0,7607  0,7625  0,7634  0,7632  0,7532
lfc_nmae     0,4249  0,4417  0,4520  0,4576  0,4602  0,4563
expr_mse     0,3953  0,4516  0,4788  0,4931  0,4994  0,4944
pds_cosine   0,9939  0,9951  0,9951  0,9951  0,9951  0,9951
AVG_SCORE    0,6597  0,6738  0,6811  0,6848  0,6862  0,6812
delta vs a=0 +0,0000 +0,0141 +0,0215 +0,0251 +0,0266 +0,0215
```

**Nenhum membro perde.** É isso que o escalar não conseguia, porque mover `k`
trocava membro por membro — a dimensão que o fator por gene usa estava
simplesmente inexplorada.

Os dois membros de ERRO ganham mais, e os dois **passam o melhor que qualquer `k`
escalar alcançou**: `expr_mse` vai a 0,4994 contra os 0,4466 de k=3, e
`lfc_nmae` a 0,4602 contra os 0,4390 de k=5. Era onde a seção 19 dizia que o
ganho estava.

`pds_cosine` SOBE (0,9939 → 0,9951) em vez de cair, que é o que a verificação de
`spearman(sd_prev, sd_aplicado) = +0,91` previa. Se o peso tivesse achatado a
estrutura por gene, ele desabaria.

**O ótimo está entre parênteses e cai exatamente em alpha=1.** A curva é
unimodal: 0,6597 → 0,6738 → 0,6811 → 0,6848 → **0,6862** → 0,6812. O braço de
alpha=1,5 existe por isso — "melhor dos testados" é mais fraco que "o ótimo está
aqui", a mesma lição que a seção 19 registrou sobre o grid de `k` parar em 4
subindo. E `alpha = 1` significa "aplique a forma medida uma vez, inteira": o
ótimo cair ali, e não em 0,5 nem em 2, é evidência de que a medida está certa —
forma forte demais poria o ótimo abaixo de 1, fraca demais acima.

---

## 21. Os dois auditores, que CHAMAM o `cell-eval2` em vez de reimplementá-lo

Os dois riscos abertos do registro, fechados com o motor que vai pontuar. A regra
em comum: `moments.jackknife_correction`, `metrics.delta._across_pert_budget`,
`_row_weights`, `_exclusion_cols`, `de_compute.compute_de`, `de.prep_de_side` e
`de.rank_de_side` são CHAMADOS de `envs/cell-eval2`, nunca portados. Uma
reimplementação concordaria com o original até o dia em que o original mudasse, e
o erro apareceria como um número plausível e errado. Por isso os dois módulos
rodam no python daquele ambiente, e por isso `pytest` entrou no `requirements.in`
dele: sem pytest lá os testes só conseguiriam PULAR, e um teste que sempre pula
não é um teste.

### `orcamento.py` — a razão `orçamento/pedido` do `expr_mse`

`expr_mse_unbiased_capped_norm` subtrai da distância uma correção de amostragem
`C_p`, mas a correção é **limitada** (issue #348): o avaliador compara o que a
submissão pede com um orçamento derivado da própria submissão, e escala por
`r = min(1, orçamento/pedido)`. Abaixo de
`PRED_CORRECTION_BUDGET_FLAG_RATIO = 0,7` registra aviso.

Medido na v6:

```
contexto   orçamento   pedido   razão    r
D            0,9274    0,2337   3,969   1,0
E            0,9127    0,2501   3,649   1,0
F            0,8961    0,2388   3,753   1,0
```

Todos os 300 alvos resolvem o próprio gene nos três contextos — taxa de câmbio
única no `_row_weights`, nada farmável. Contra o registro: v4 (k=1) dava **0,179**,
abaixo do flag; v5 (k=4) dava 2,854. O risco está fechado, e a mesma escolha de
`k_desvio` que resolveu `fidelity` resolveu isto de lado: subir a dispersão ENTRE
perturbações sobe o orçamento e não sobe `C_p`, que depende do espalhamento
DENTRO do grupo.

**O certificado é de um lado só, e é o lado que importa.** O pedido do avaliador
usa `C_p` já limitado por `min(C_p, PRED_TRACE_CAP_K * C_p_real)`, e `C_p_real` é
do lado REAL, retido em D/E/F. Mas limitar só pode REDUZIR o pedido, e reduzir o
pedido AUMENTA a razão. Logo `razão_sem_limite ≤ razão_do_avaliador`, e
`razão ≥ 1` certifica `r = 1` sem o lado real. Abaixo de 0,7 o veredito é "em
risco" e não "sinalizado", porque o limite ainda pode resgatar — e o relatório
diz qual dos dois.

**Memória: um grupo por vez, e é exato.** O jackknife é `O(n·G)` DENSO em blocos
e a submissão tem 682 milhões de não nulos por contexto, 5,5 GB em CSR. Varrer
grupo a grupo dá o MESMO vetor, verificado idêntico **bit a bit** contra a chamada
com todos os grupos juntos. Cada grupo é ~18 MB. O `chunk` fica no padrão 512
embora o próprio docstring do pacote meça 32 como mais rápido E menor: o objetivo
é reproduzir o número que o avaliador vai calcular, e ele usa o padrão.
Fidelidade antes de velocidade quando o produto é uma auditoria.

### `auditoria.py` — e o proxy de limiar estava errado por 4,8× para o lado contrário

A `varredura_k.py` conta `|delta| > 0,01` e chama o número de aproximação
declarada: 5.532 em k=4, "2,02× os 2.739 reais". Mas `|delta|` não é o teste de
DE. O motor de verdade dá **mediana 1.154, 0,42×**.

Confirmação independente de que o instrumento está certo: o gate de CPM dá
**9.633** genes em D, exatamente o universo registrado em `ctrl_bulk`.

E a medida é EXATA, não proxy, por três propriedades do preset:

1. **`control_source: real`** — o lado de controle do teste vem do dado REAL, e
   os controles de D/E/F estão publicados em `controls-final`. Não há lado retido
   nesta conta: são nossas 400 células preditas contra os 18.400 controles reais.
2. **`fdr_scope: per_pert`** — `de._bh_per_target` refaz o BH dentro de cada alvo
   sobre os p-valores dele, então o `nsig` de um alvo não depende de quais outros
   entraram. Amostrar alvos é EXATO para os amostrados, e é o que torna a
   auditoria viável: 900 testes de Wilcoxon contra 18.400 células são horas.
3. **`filter_gene_min_cpm_cell: 5,0` com UM conjunto de genes para todos os
   alvos** (#351) — o gate é a média aritmética do CPM por célula da REFERÊNCIA,
   que é o controle real. Universo fixo por contexto, independente da nossa
   predição. Um gate que dependesse da predição poderia ser movido por ela.

**E a divisão é BINÁRIA** (contexto D, amostra de 40 alvos):

```
classe            n   mediana nsig   zeros   faixa
com desvio       23       1.845         0    540 - 4.018
sem evidência    17           0         9    0 - 3
```

Os 23 alvos com desvio específico dão 0,67× os reais — ordem certa. Os 143 de 300
que carregam **só a resposta média** produzem praticamente ZERO.

Mecanismo medido: a parcela média move cada gene **1,6% RMS em log2** (razão de
proporção p05 0,985 / p95 1,015) e a parcela de desvio é **12× maior** (rms 0,093
contra 0,0078 em espaço de delta). Um deslocamento de 1,5% por gene fica abaixo
do piso de detecção de 400 células contra 18.400 controles.

**E a aposta que isso representa é RUIM, medido.** `sig_jaccard` pontua
vazio-contra-vazio como 1 e vazio-contra-cheio como 0, então emitir zero só ganha
se "sem evidência doadora" for proxy de "sem DEG real". Não é: nos dois contextos
profundos do Replogle, estar no split de treino prediz resposta MENOR (mediana
220 contra 567 em rpe1, p = 1,3e-14; 52 contra 127 em k562, p = 1,5e-4) mas a
fração com ZERO genes significativos é praticamente IGUAL (0,054 contra 0,039;
0,105 contra 0,104). O zero vem de ausência de DADO, não de ausência de RESPOSTA.
Com a taxa base de ~36,5% de perturbações realmente sem DEG (73 de 200), cerca de
**91 dos 143** são realmente não-vazios e hoje pontuam ~0 em quatro dos seis
membros.

**Isso reabre `k_media`, que nunca foi varrido porque "o proxy é cego a ele".**
Verdade para o proxy, cuja verdade é Orion transplantado e não contém a resposta
média de D — FALSO para o `auditoria.py`, que usa os controles REAIS e portanto vê
a parcela média. O que o auditoria mede é quantos genes um dado `k_media` produz;
o que ele NÃO mede é se são os genes certos, e para isso o que existe é a
transferência medida da média (seção 15: correlação mediana +0,4736, ganho
+0,2520, 24 de 24 direções positivas). É decisão de quem assina a submissão,
porque muda a aposta em 143 dos 300 alvos.

### Dois consertos que o caminho cobrou, os dois calados

O gate de reuso de bundle do `avaliador_varredura.sh` procurava
`$OUT/bundle/bundle.json` e `$OUT/bundle/splits`, e o `prep-real-bundle` do
cell-eval2 0.19.0 não escreve nenhum dos dois — ele escreve `manifest.json`,
`config.yaml`, `baseline_agg.csv` e três `anchor_*`. Na varredura de k isso nunca
apareceu porque os seis braços vinham numa invocação só e o bundle era construído
nela. Com um braço por chamada o gate passou a decidir, e decidia sempre
reconstruir: 567 s de `prep-real-bundle` mais 127 s de `baseline` por braço.
Medido antes e depois: **16 min no primeiro braço, 6 min no seguinte**. O gate
agora é `manifest.json`, que é `real_bundle.MANIFEST` e o único arquivo que
precisa existir para o bundle ser legível.

E `gerar_submissao.py` não tinha o idioma `sys.path.insert(0, parent)` dos irmãos
`varredura_k.py` e `pool_media.py`, e morria em `ModuleNotFoundError: No module
named 'conversao'` sob `-I` — que é o modo que o projeto exige para script que lê
dado baixado.
