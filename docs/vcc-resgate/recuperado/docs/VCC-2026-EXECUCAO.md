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

