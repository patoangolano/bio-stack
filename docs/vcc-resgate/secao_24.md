
## 24. Depois de k=4: duas ideias mortas por medida, um mecanismo achado, e F4 resolvido pela metade

Com a magnitude calibrada, a decomposição do que falta diz onde procurar. Contra o
oráculo de 0,970, a v5 em 0,647 deixa 0,323, e ele não está distribuído:

| membro | k=4 | oráculo | disponível em `avg_score` |
|---|---:|---:|---:|
| `de_lfc_nmae` | 0,427 | 1,024 | **+0,0995** |
| `expr_mse` | 0,392 | 0,947 | **+0,0925** |
| `sig_jaccard` | 0,567 | 1,032 | +0,0775 |
| `direction_reach` | 0,722 | 0,989 | +0,0445 |
| `direction_fidelity` | 0,783 | 0,828 | +0,0075 |
| `pds_cosine` | 0,989 | 1,000 | +0,0018 |

`direction_fidelity` e `pds_cosine` estão **fechados**: 0,009 somados. **60% do que
resta está nas duas métricas de erro**, que querem acurácia da média predita, não
dispersão. Magnitude esgotada.

### Duas ideias que pareciam boas e morreram em cinco minutos cada

**Teto na saída.** A §21 propôs k=4 com teto, porque `expr_mse` virou em k=3 e o
`LIMITE_LOG2R = 5,0` limita o doador antes do ajuste, não a saída depois de
multiplicar por `k`. Medido, o desenho estava errado duas vezes:

1. **A cauda negativa já satura de graça.** `delta = log1p(p·2^log2r·renorm) − b_ctrl`
   tem piso em `−b_ctrl`: predizer um gene em zero é a repressão mais forte
   possível e não há como ir além. O mínimo medido é **−3,2894 = exatamente
   −`b_ctrl`**, e só **2 entradas** chegam a 99% do piso. Uma delas é
   `HMGB2`/`HMGB2` — o alvo derrubando a si mesmo, e o gene-alvo é excluído de
   todas as métricas. O `log2r` chega a **−23,3** e isso é inofensivo.
2. **A cauda positiva não move um MSE de 2,89 M entradas.** 22 entradas passam de
   |delta| = 2. O que cresceu foi o corpo: |delta| > 0,5 vai de 37 (k=1) a 4.966
   (k=4). Cortar 22 outliers de 2,89 M não recupera os 0,050 que k=3→k=4 custou
   ao `expr_mse`.

**`k` por região de gene.** `expr_mse` soma sobre `G_p = G − {g_p}` — **todos** os
18.533 genes menos o alvo — enquanto `lfc_nmae`, `jaccard`, `reach` e `fidelity`
vivem no universo de DE (9.633 em D). Os conjuntos **diferem**, então aplicar `k`
alto dentro do universo e baixo fora pareceria almoço grátis: ~8.900 genes
contribuem para o `expr_mse` e para nenhuma métrica de DE.

Medido, a massa fora do universo é **0,61% / 0,89% / 0,91%** do delta² total.
Gene fora do universo tem `b_ctrl` minúsculo, então a conversão composicional já
lhe dá delta minúsculo — rms 0,0055 contra 0,0675 dentro. Zerar os 8.900
inteiramente mexeria menos de 1% no numerador. Ideia morta.

### O mecanismo de por que `expr_mse` SUBIU com `k`

A §21 registrou que a inferência da literatura errou nisso: ela previu que `mse`
quereria encolher, e `mse` subiu de 0,228 a 0,442 até k=3. A causa está na
fórmula, lida no `cell-eval2`:

```
expr_mse_unbiased_capped_p = ( ||mu_pred − mu_real||²
                               − r · min(tr Σ_pred/n_pred, k·tr Σ_real/n_real)
                               − tr Σ_real/n_real ) / |G_p|

r = min(1, B_pred / Σ_q w_q min(tr Σ_pred/n_pred, k·tr Σ_real,q/n_real))
B_pred = espalhamento ENTRE perturbações das médias preditas
```

A métrica devolve crédito pelo nosso próprio ruído de amostragem — mas **escalado
por `r`**, e `r` é a razão entre o espalhamento que a submissão de fato tem e a
correção que ela reclama. É a defesa do issue #348 contra farmar a correção.

Medido, `B_pred` exato dos nossos deltas contra a variância jackknife das 400
células emitidas (6 grupos de D, `tr Σ_pred/n_pred` ≈ 14,97 por grupo):

| versão | `B_pred` | claim | razão | `r` | flag de 0,7 |
|---|---:|---:|---:|---:|---|
| v4 (k=1) | 0,0438 | 0,2444 | **0,179** | 0,179 | **dispara** |
| v5 (k=4) | 0,6975 | 0,2444 | **2,854** | 1,000 | limpa |

Em k=1 recebíamos **17,9%** do crédito a que tínhamos direito, e os outros 82,1%
do nosso próprio ruído de emissão iam para o numerador como penalidade pura. `r`
chega a 1 em `k = √(0,2444/0,0438) = 2,36`, e `expr_mse` tem pico em k=3 — o
ponto medido imediatamente acima. A curva tem explicação de primeiros princípios,
não só forma.

E o `cell-eval2` **loga um aviso** abaixo de razão 0,7: *"its per-cell scatter
largely CANCELS in the pseudobulk, which is not what an i.i.d. sample of a
predicted cell population looks like — honest arms measure ~1.0 ... the
submission is worth a look"*. A v4 teria disparado isso. A v5 não.

`B_pred` não é alavanca: o docstring do `_across_pert_budget` diz que injetar
variação centrada **sobe o numerador e o orçamento na mesma medida**, então
comprar orçamento é empate na melhor hipótese. Fechado por construção.

### F4, medido: dois defeitos acoplados

**O gate serve duas parcelas com requisitos diferentes.** O laço do `ajustar`
abre com `usaveis = alvos ∩ disponiveis ∩ split; if not usaveis: continue`, e a
contribuição da **média** estava dentro desse `continue`. O gate é correto para o
desvio, que precisa do nosso alvo medido naquele contexto. A média não precisa de
alvo nosso nenhum: `media_do_contexto()` percorre **todas** as perturbações do
doador. Consequência: os **4 contextos do Replogle entravam com peso zero**,
porque nenhum dos nossos 300 alvos está no split de treino do Replogle — e são os
de melhor transferência medida de média do conjunto inteiro.

**O peso contradizia a medida.** Peso efetivo contra transferência medida:

| contexto | peso | share | ganho como doador |
|---|---:|---:|---|
| `arc_vcc_25.ARC_H1` | 2016 | **93,2%** | **+0,0000** (r = −0,0001, 8.800 genes) |
| `jiang` (11 ctx) | 143 | 6,6% | não medido |
| `orion.hct116` | 2 | **0,1%** | **+0,3917** (6 direções) |
| `orion.hek293t` | 2 | 0,1% | +0,6589 (1 direção) |
| `replogle` (4 ctx) | 0 | 0,0% | +0,195 a +0,434 |

O zero do `arc_vcc_25` não é artefato: 8.800 genes comparáveis — mais que os 6.287
do par que transfere melhor — e magnitudes equivalentes (0,034 contra 0,043).

### O teste, e o erro que ele me fez cometer primeiro

`scripts/vcc/pool_media.py` deixa um contexto de fora, agrupa a média dos outros
com cada esquema, e mede o ganho ao predizer a média do retido. Métrica idêntica
à do `resposta_media.py`. Isso fecha **sem** D/E/F, o que importa porque o
`avaliador_proxy` é cego à parcela da média.

Primeiro resultado, corte de 40 perturbações, 8 contextos:

```
uniforme       +0,2176    <- vencedor
amostra        +0,2025
concordancia   +0,1685
profundidade   +0,1660
atual          +0,1002
```

Eu ia implementar peso igual. Então olhei as duas linhas que mais se parecem com
o nosso alvo — os dois contextos do ensaio da própria competição — e
`profundidade` ganhava **5×** delas (+0,5248 e +0,5503 contra +0,1038 e +0,1713).
Com `ARC_H1` retido, `profundidade` põe **91% do peso em `ARC_H1_VAL`**: o mesmo
estudo. **Atalho de mesmo-estudo, indisponível para D/E/F.**

Refeito com `--sem-mesmo-repo`, que retira também os irmãos do repo do contexto
retido — o análogo honesto da nossa situação:

```
profundidade   mediana +0,0983   media +0,0859
uniforme       mediana +0,0579   media +0,0873
atual          mediana +0,0524   media +0,0618
concordancia   mediana +0,0403   media +0,0463
amostra        mediana +0,0231   media +0,0448
```

**O vencedor inverte, e o sinal fica misto:** mediana diz profundidade, média diz
peso igual, e as médias diferem em **0,0014**. Todos os ganhos também caem por um
fator de 2 a 3 — a maior parte da transferência aparente era mesmo-estudo.

Minha hipótese original — pesar por concordância medida, que é o lever 4 da
revisão de literatura — **perdeu em todas as quatro configurações** (cortes de 10,
40, 200, e sem mesmo estudo). Era a ideia mais atraente da página e é a que a
medida descarta com mais clareza.

### O que foi mudado, e o que não

**Mudado: o gate.** A média passa a ter laço próprio, sobre todo contexto com ≥ 40
perturbações, sem o gate do desvio. É a única parte robusta: `profundidade` sem
gate bate `atual` com gate em **4 de 4** configurações (+0,0983/+0,0524;
+0,1660/+0,1002; +0,0100/+0,0044; +0,2816/+0,1058).

**Não mudado: o peso.** Sinal misto da ordem de 0,0014 não escolhe vencedor. A
regra que eu dei ao agente na varredura de `k` — *"se der sinais mistos de
magnitude comparável ao ruído, diga que é ruído em vez de escolher um vencedor"* —
vale para mim.

O corte de 40 perturbações foi fixado **antes** de medir, com razão declarada (a
média de um contexto com poucas perturbações é ruído; no limite de uma, ela *é*
aquela perturbação). Ele cai uma perturbação acima do maior contexto do Jiang
(39), o que é perto demais para passar sem teste — daí os cortes de 10 e 200. Em
10, os 30 contextos do Jiang entram e **todos os ganhos colapsam** (mediana 0,043
contra 0,218); em 200 o `arc_vcc_25` sai e `profundidade` vence, mas esse corte
remove justamente a patologia que torna profundidade ruim, então não é evidência a
favor dela.

### Em aberto, e registrado porque incomoda

Com o gate fora, `profundidade` dá **95,5%** do agrupamento aos dois contextos do
`arc_vcc_25` (2 × 2016 contra 4 × 46 do Replogle e 2 × 2 do Orion) — concentração
**maior** que os 93,2% de antes, e justamente no doador cuja média é ortogonal à
do Orion. Nas duas linhas de Orion do teste sem mesmo estudo, `profundidade` mede
+0,0031 e +0,0009 contra +0,0329 e +0,0038 de peso igual. Qual contexto se parece
com D/E/F é desconhecido, e é por isso que nenhum dos dois pode ser escolhido com
a evidência que existe.

### Teste nulo: a mudança não vazou para o desvio

A correção toca só a média, e o desvio é o que o proxy mede. Então o desvio tem de
ficar intacto — e se não ficar, vazou. Medido:

```
D/E/F   desvio: max|dif| 1,05e-04 contra rms 0,0325
```

0,3% relativo, o que é pequeno mas não é zero. A diferença é **constante por
alvo**: desvio-padrão dentro da linha **3,68e-09** contra média de linha de até
1,05e-04, e exatamente **156 de 300** linhas têm constante não nula — as 156 com
desvio específico em D; as 144 sem evidência têm zero exato. Removida a constante,
o resíduo é **2,27e-07 contra rms 0,0354, seis partes por milhão**.

A constante é a renormalização composicional, que muda porque o total mudou. Não é
vazamento. E revela uma imprecisão pequena do `varredura_k.py`: a decomposição
atribui a constante da renormalização à parcela do **desvio**, então `k × desvio`
escala a constante também. A 1e-4 contra rms 0,035 isso é 0,3% e não muda nenhuma
conclusão, mas é real e fica escrito em vez de implícito.
