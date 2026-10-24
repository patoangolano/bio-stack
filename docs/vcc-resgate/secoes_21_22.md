
## 21. A varredura de `k_desvio`: o buraco era calibração, e fechou

A seção 18 nomeou o buraco com número: `direction_fidelity` em **0,022 de um
teto de 0,828** — 3% do alcançável — contra `pds_cosine` em 99%. O modelo
acertava o padrão e errava a direção por perturbação. Direção por perturbação é,
por construção, a parcela do **desvio específico**, que o β de mínimos quadrados
encolhe a ~0,22.

A hipótese era que isso fosse uma troca: mínimos quadrados é ótimo para erro
quadrático e sub-disperso para tudo que depende de acertar o sinal, então subir a
magnitude deveria comprar `fidelity` e `jaccard` ao custo de `expr_mse`.

`scripts/vcc/varredura_k.py` decompõe o delta em suas duas parcelas — resposta
média e desvio específico — em espaço de **log2 da razão**, porque em espaço de
delta elas não se somam (`log1p` de proporção renormalizada não é linear). A
decomposição é exata e dispensa reajustar o modelo: os 143 alvos sem evidência
carregam **apenas** a média e são idênticos entre si, então o delta de qualquer um
deles *é* a parcela média. O script afirma essa identidade em tempo de execução e
falha se ela não valer.

### O que foi medido

Onze braços no mesmo bundle `proxy-D-orion-k40`, mesmos 40 alvos, oito campos de
identidade de par batendo em todos.

| membro | k=1 | k=1,5 | k=2 | k=3 | **k=4** | k=5 | k=6 | k=8 | oráculo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `direction_fidelity` | 0,022 | 0,129 | 0,283 | 0,589 | 0,783 | **0,841** | 0,826 | 0,790 | 0,828 |
| `sig_jaccard` | 0,040 | 0,173 | 0,317 | 0,499 | 0,567 | **0,582** | 0,576 | 0,546 | 1,032 |
| `direction_reach` | 0,267 | 0,441 | 0,551 | 0,667 | 0,722 | 0,750 | 0,765 | **0,770** | 0,989 |
| `de_lfc_nmae` | 0,102 | 0,172 | 0,239 | 0,353 | 0,427 | **0,438** | 0,378 | 0,090 | 1,024 |
| `expr_mse` | 0,228 | 0,315 | 0,380 | **0,442** | 0,392 | 0,165 | 0,000 | 0,000 | 0,947 |
| `pds_cosine` | **0,990** | 0,990 | 0,990 | 0,990 | 0,989 | 0,988 | 0,978 | 0,957 | 1,000 |
| **`avg_score`** | 0,275 | 0,370 | 0,460 | 0,590 | **0,647** | 0,627 | 0,587 | 0,526 | 0,970 |

### A hipótese estava certa no diagnóstico e errada no mecanismo

**Não há troca até k=3 — tudo sobe junto.** De k=1 a k=3: `fidelity` ×27,
`jaccard` ×12,6, `reach` ×2,5, `lfc_nmae` ×3,5, **e `expr_mse` também**
(0,228 → 0,442). O β não estava trocando direção por magnitude. Estava
sub-disperso para os seis ao mesmo tempo, **inclusive para o membro que ele
teoricamente otimiza**.

A troca existe, mas começa longe de 1: `expr_mse` vira em k=3,
`fidelity`/`jaccard`/`lfc_nmae` em k=5, `reach` ainda sobe em k=8. O `avg_score`
concilia em k=4.

`fidelity` em k=5 dá **0,841 contra 0,828 do oráculo** — passa o quase-oráculo.
Não é paradoxo: o oráculo carrega ruído de amostragem de 400 células, e uma
predição com desvio amplificado tem direção mais limpa que uma medição ruidosa.

### Por que a objeção do proxy é mais fraca aqui que na seção 18

A seção 18 registrou que este instrumento é **cego à parcela da média**: a verdade
dele é Orion/HCT116 transplantado, que não contém a resposta média de D. Foi por
isso que a correção do viés de seleção (seção 19) moveu `avg_score` em +0,0028,
dentro do ruído.

`k_desvio` escala só o **desvio específico**, que é exatamente a parcela que a
verdade do proxy **contém**. É a pergunta para a qual este instrumento é o certo.
Ganho de k=1 para k=4: **+0,372**, que é 133× o +0,0028 anterior.

### Três coisas que a tabela esconde

**`expr_mse` em 0,0000 para k=6 e k=8 é corte de escala, não empate com o
controle.** O bruto vai a 1,337 e 3,092 contra ~1,0 do baseline — fica **pior que
o baseline** e a escala corta em zero. O `avg_score` desses dois está otimista; a
queda real é mais funda. Isso importa para a assimetria do risco: errar `k` para
baixo custa pouco e devagar, errar para cima cai num penhasco.

**k=1 é a v4, verificado numericamente.** `delta_kd1.npz` reproduz
`delta_previsto.npz` a 7e-8 relativo (máx. da diferença 9,3e-10 contra rms de
0,013, nos três contextos). Se a decomposição não fechasse, a varredura mediria o
erro dela e não `k` — era o único check que podia invalidar a tabela inteira.

**A cauda de magnitude cresceu 53×, e é uma entrada.** `escala_de_gene_maxima` vai
de 3,81 na v4 para **199,17** na v5. Medido entrada por entrada: é sempre o mesmo
par, `SIN3B` → `DDIT3`, e apenas 46 a 54 das 5.559.900 entradas passam de 10×; o
p99,99 da escala é 3,5 a 3,6. O mecanismo é coerente — SIN3B é andaime do
corepressor SIN3/HDAC, e perder repressão desreprime CHOP — e a magnitude é
exatamente `k` vezes o maior log2r do ajuste (1,9 × 4 = 7,6). Nada de contrato é
violado: 41,74% do teto de não nulos, contagem máxima por célula 66.829 contra
teto de 1.000.000. O `LIMITE_LOG2R = 5,0` do `transfer_fit` limita o log2r do
**doador** antes do ajuste, não a saída depois de multiplicar por `k` — então a
cauda não está sob limite nenhum. O proxy pontuou k=4 como vencedor **com** essa
cauda presente, logo ela está precificada; o experimento natural seguinte é k=4
com teto de saída, que pode recuperar parte do `expr_mse` que k=4 perdeu.

### O que isto autoriza, e o que não

**Autoriza `k_desvio = 4`.** Curva unimodal e suave, monótona dos dois lados do
pico, sem sinal misto, num instrumento adequado à pergunta.

**Não autoriza k=4 como ótimo do painel.** Sem réplica e sem intervalo de
confiança; e o ótimo de uma soma de seis membros depende das pontas `b` e `r`, que
no painel são as de D/E/F e não as do proxy. O que transfere é a **direção**
(k ≫ 1) e a **ordem de magnitude** (4, não 1 nem 40) — não a terceira casa. O
platô de 3 a 5 é largo (0,590 / 0,647 / 0,627), então a escolha não está no fio da
navalha; e como a queda acima de 5 é um penhasco e abaixo de 3 é uma rampa, ficar
em 4 erra para o lado barato.

`runs/vcc/transfer/delta_calibrado_kd4.npz` é o delta promovido a artefato de
primeira classe — cópia imutável de `varredura_k/delta_kd4.npz`, que uma próxima
corrida da varredura sobrescreveria sem avisar. `submissao_v5.h5ad` carrega esse
arquivo.

## 22. Pendências de revisão: duas fechadas, duas medidas e abertas

### F7 — literais obsoletos, fechado na raiz

O código afirmava `correlacao mediana 0,436, ganho mediano 0,164` num comentário e
`faixa_medida: [0.041, 0.788]` no resumo. Os valores medidos em
`resposta_media.json` são **0,474** e **0,252**; a faixa real em
`validacao_loco.json` é **[0,2178, 0,8735]** com mediana **0,4235**.

A revisão que encontrou a divergência propôs `[0,0415, 0,9056]` — **também
obsoleto**, de uma corrida intermediária. Dois valores errados para o mesmo campo
em duas gerações é diagnóstico de classe, não de instância: corrigir o literal
reproduz o defeito na corrida seguinte. `transfer_fit.py` passa a **ler** esses
números dos arquivos que já tem em mão (`medida` e `validacao` estão no escopo
onde o resumo é montado).

### F8 — a mediana de dois é a média dos dois, registrado

`beta_mediano = np.median(np.stack(...), axis=0)` sobre os betas medidos. Com dois
doadores com direção de referência (`orion` e `replogle`), a mediana **é** a média
deles. A aritmética está certa — a média de dois betas é a resposta defensável, e
o comentário original já dizia que o objetivo era não escolher arbitrariamente
entre eles. O problema é a etiqueta: "mediana" promete robustez a outlier que
`n=2` não tem. O resumo agora registra quantos betas entraram e diz explicitamente
qual das duas coisas aconteceu.

### F4 — a multiplicidade de contexto, medida: pior do que a revisão disse

A ponderação por profundidade acumula `profundidade[repo]` **uma vez por
contexto** do doador. A revisão apontou que isso dá ao Jiang (11 contextos × 13)
peso 143 contra os 4 do Orion (2 × 2), e chamou de acidente de contagem de
contexto. O peso efetivo medido:

| repo | contextos usáveis | nossos alvos no split | profundidade | peso | share |
|---|---:|---:|---:|---:|---:|
| `arc_vcc_25` | 1 (`ARC_H1`) | 12 | 2016 | 2016 | **93,2%** |
| `jiang` | 11 | 25 somados | 13 | 143 | 6,6% |
| `orion` | 2 (`hct116`, `hek293t`) | **158** | 2 | 4 | **0,2%** |
| `replogle` | 0 | 0 | 46 | 0 | 0,0% |

O agrupamento da resposta média é **93% um único contexto doador**, e o Orion —
que carrega 146 dos 157 alvos com evidência — entra com **0,2%**. A multiplicidade
de contexto do Jiang é a menor parte do problema.

Duas leituras defensáveis, e nenhuma medida:

- **A favor do peso por profundidade.** A média de um ensaio raso é atenuada para
  zero por ruído de medida, e pior: o ruído do Orion **não é de média zero** — é
  exatamente o viés de seleção do `tested` da seção 19, que tem corr = −0,878 em
  `orion.hct116`. Pesar o Orion para baixo é defesa contra o viés residual que o
  `FRACAO_MINIMA_COBERTURA = 0,9` não pegou.
- **Contra.** `media_do_contexto()` percorre **todas** as perturbações do contexto
  doador, não só nossos alvos. A média de milhares de perturbações do Orion é um
  estimador de variância muito menor que a de 146 do `ARC_H1`. Se o viés residual
  for pequeno, 93% para uma fonte é desperdício de informação.

O instrumento que decide existe: `resposta_media.py` mede transferência da média
entre as 24 direções, e já mostrou que **proximidade de linhagem não prevê
transferência** (melhor par `hepg2→rpe1`; pior `orion.hct116→orion.hek293t`).
Pesar doador por **concordância medida** em vez de profundidade declarada é o
lever 4 da revisão de literatura. Fica aberto por sequenciamento, não por
dificuldade: o ganho de k=1 para k=4 é +0,372 medido, F4 é zero medido. Uma
mudança por vez, cada uma com sua própria medição.

### F5 — centragem inconsistente entre `validar` e `ajustar`, aberto

`validar` e `ajustar` não centram o `log2r` do mesmo jeito. O efeito é atenuante,
então o resultado atual é conservador — **certo por acidente, não por desenho**.
Não foi tocado: mexer nisso move o β, que move a calibração de `k` que acabou de
ser medida. A ordem é medir F5 contra o mesmo bundle, não corrigir e esperar que
k=4 continue valendo.

### `se_jackknife` — calculado, consumido por ninguém

`ctrl_bulk.py` calcula erro padrão jackknife delete-1 em forma fechada (mediana
3,95e-03). É o código mais arriscado do arquivo e **não alimenta a submissão**.
Fica como está: apagar perde a medida, usar sem desenho adiciona um parâmetro sem
evidência.
