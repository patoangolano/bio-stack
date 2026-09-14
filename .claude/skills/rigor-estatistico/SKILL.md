---
name: rigor-estatistico
description: Regras de correção estatística para expressão diferencial e single-cell nesta stack — o que padj significa em cada ferramenta, quando um p-valor é NaN, por que nunca se filtra log2FC encolhido, e como fixar semente. Use ao rodar pydeseq2 ou scanpy, ao contar DEGs, ao escrever a seção de Métodos, ou ao revisar um resultado de expressão.
---

# Rigor estatístico — o que as ferramentas desta stack realmente fazem

Tudo aqui foi lido do código **instalado nesta máquina**, não da documentação.
Documentação descreve intenção; o código descreve o que roda. Onde as duas
divergem, o código ganha.

Versões verificadas em 2026-09-14: `pydeseq2 0.5.4`, `scanpy 1.11.5`,
`statsmodels 0.15.0`, `scipy 1.17.1`, `numpy 2.4.6`. Travadas em `uv.lock`.
Se essas versões mudarem, **reverifique antes de reusar qualquer frase abaixo** —
cada afirmação traz o comando que a reproduz.

---

## 1. `padj = NaN` significa "filtrado fora", nunca "não significativo"

No `pydeseq2` instalado, duas coisas produzem `NaN` em `padj`:

**Filtro independente.** `_independent_filtering` cria o quadro de resultado
inteiro como `NaN` e preenche só as linhas que passaram do corte de `baseMean`:

```
.venv/Lib/site-packages/pydeseq2/ds.py:510   result = pd.DataFrame(np.nan, ...)
.venv/Lib/site-packages/pydeseq2/ds.py:517   result.loc[use, i] = false_discovery_control(U2, method="bh")
```

**Outlier de Cook.** Gene com amostra influente demais tem o p-valor zerado:

```
.venv/Lib/site-packages/pydeseq2/ds.py:554   self.p_values[self.dds.cooks_outlier()] = np.nan
```

### Por que isso morde

`res[res["padj"] < 0.05]` avalia `False` em `NaN`. A **contagem de DEG sai
certa**. O que sai errado é o denominador: se você escrever "de 18.000 genes,
420 foram diferencialmente expressos", os 18.000 estão errados — milhares nunca
foram testados. Esse é o número que um revisor confere primeiro, e o erro é
invisível a menos que você conte os `NaN`.

### O que fazer

```python
res = stat_res.results_df
n_testados = int(res["padj"].notna().sum())     # entraram no teste
n_filtrados = int(res["padj"].isna().sum())     # filtro independente + Cook
deg = res[res["padj"] < 0.05]
# Métodos: "de N genes detectados, n_testados entraram no teste após filtro
# independente; n_filtrados foram removidos."
```

Reproduzir a evidência:

```bash
grep -n "false_discovery_control\|np.nan\|cooks_outlier" .venv/Lib/site-packages/pydeseq2/ds.py
```

---

## 2. `padj` do pydeseq2 vem do `scipy`, não do `statsmodels`

```
.venv/Lib/site-packages/pydeseq2/ds.py:10    from scipy.stats import false_discovery_control
.venv/Lib/site-packages/pydeseq2/ds.py:544   self.padj.loc[~self.p_values.isna()] = false_discovery_control(...)
```

Duas consequências:

- **Nunca rode `statsmodels.multipletests` sobre a saída do pydeseq2.** `padj`
  já está corrigido; corrigir de novo é correção dupla, e o resultado não
  significa nada.
- A família corrigida **depende do corte de filtro escolhido**. O pydeseq2
  testa 50 quantis de `baseMean`, aplica BH dentro de cada um, e escolhe o
  corte por ajuste lowess sobre a contagem de rejeições. Ou seja: o conjunto de
  genes sobre o qual o FDR foi controlado não é "todos os genes" — é o conjunto
  que sobreviveu ao corte selecionado. Isso precisa aparecer nos Métodos.

**Frase de Métodos, verificável:**

> Valores de p ajustados foram calculados pelo PyDESeq2 v0.5.4 por FDR de
> Benjamini–Hochberg via `scipy.stats.false_discovery_control`, aplicado após
> filtragem independente sobre 50 quantis de `baseMean` com ajuste lowess para
> seleção do limiar (α = 0,05, filtro de outliers de Cook ativo).

---

## 3. log2FC encolhido não é critério de significância

```
.venv/Lib/site-packages/pydeseq2/ds.py   def lfc_shrink(self, coeff, adapt=True):
    """LFC shrinkage with an apeGLM prior.
    Shrinks LFCs using a heavy-tailed Cauchy prior, leaving p-values unchanged."""
```

`lfc_shrink()` reescreve `log2FoldChange` e `lfcSE`, mas **não toca em
`pvalue` nem `padj`**, e deixa `stat` como a estatística de Wald **não
encolhida**. Depois do encolhimento o quadro é deliberadamente misto:
`log2FoldChange / lfcSE != stat`.

Portanto:

- Filtrar por `abs(log2FoldChange) > 1` **depois** de encolher não é um teste.
  É um corte arbitrário sobre uma estimativa bayesiana, apresentado como se
  carregasse significância.
- Recalcular z a partir dessas duas colunas dá um número que não é o teste que
  produziu o p-valor.

**Se você quer mesmo testar um tamanho de efeito**, o pydeseq2 faz isso direito
— refaz o teste de Wald contra o novo nulo:

```python
stat_res = DeseqStats(dds, contrast=["condicao", "tratado", "controle"],
                      lfc_null=1.0, alt_hypothesis="greaterAbs")
```

Isso devolve `pvalue` e `padj` válidos para |LFC| > 1. Note que muda a hipótese
testada: resultados não são comparáveis a uma rodada com `lfc_null=0`. Diga qual
das duas produziu cada figura.

---

## 4. O FDR do scanpy é por grupo, não pela família toda

Em `scanpy 1.11.5`, a correção está **dentro** do laço de grupos:

```
.venv/Lib/site-packages/scanpy/tools/_rank_genes_groups.py:442   for group_index, scores, pvals in generate_test_results:
.venv/Lib/site-packages/scanpy/tools/_rank_genes_groups.py:464       if corr_method == "benjamini-hochberg":
.venv/Lib/site-packages/scanpy/tools/_rank_genes_groups.py:468           _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
```

A família corrigida é **um grupo × seus genes**. Não há correção sobre a família
conjunta grupos × genes. Com 12 clusters, você fez 12 correções independentes e
o FDR global é maior que o α declarado.

Isso não é bug — é a escolha do scanpy, e é defensável para marcadores
exploratórios. O que não é defensável é omitir.

**Frase de Métodos:**

> A correção de FDR foi aplicada independentemente dentro de cada grupo pelo
> scanpy v1.11.5 (`statsmodels`, `fdr_bh`, α = 0,05); não houve correção sobre
> a família conjunta grupos × genes.

### `method='logreg'` não tem p-valor nenhum

Com regressão logística, `pvals` é `None` e **nada** é corrigido. Não existe
`pvals_adj` para reportar. Se você usou `logreg`, os "marcadores" são um
ranqueamento por peso de modelo, não um teste — descreva assim.

### Análise de sensibilidade (opcional, honesta)

Se quiser reportar o FDR conjunto, faça como *sensibilidade*, não como
substituição:

```python
from statsmodels.stats.multitest import multipletests
todos_p = np.concatenate([r["pvals"] for r in resultados_por_grupo])
_, padj_familia, _, _ = multipletests(todos_p, alpha=0.05, method="fdr_bh")
# "Sob correção conjunta, X dos Y marcadores permanecem abaixo de 0,05."
```

---

## 5. Semente: nunca `np.random.seed()`

`np.random.seed()` mexe no gerador **global do processo**. Qualquer biblioteca
importada depois pode consumir números do mesmo gerador e deslocar a sequência —
sem erro, sem aviso. O script parece determinístico e não é.

O lint recusa: regra `NPY002` está ativa em `pyproject.toml`.

```python
# errado
np.random.seed(42)

# certo — gerador próprio, isolado
rng = np.random.default_rng(42)

# certo — parâmetro da biblioteca
sc.tl.leiden(adata, random_state=0, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata, random_state=0)
sc.pp.pca(adata, random_state=0)
```

Para Leiden, fixe também `flavor`, `n_iterations` e `directed`: o padrão do
scanpy já mudou entre versões, e `random_state` sozinho não protege de uma
mudança de algoritmo. Registre os quatro nos Métodos.

---

## 6. Checklist antes de declarar um resultado

Rode mentalmente, na ordem:

| # | Pergunta | Onde falha silenciosamente |
|---|---|---|
| 1 | Qual é o n por grupo? | `< 3` torna a estimativa de dispersão do DESeq2 não confiável |
| 2 | O desenho inclui os confundidores? | `design="~condicao"` quando havia lote é viés, não ruído |
| 3 | Quantos genes entraram no teste? | `padj.notna().sum()`, não o total detectado |
| 4 | Estou filtrando LFC encolhido? | Ver §3 |
| 5 | A correção cobre a família que eu digo cobrir? | Ver §4 |
| 6 | Células são réplicas? | `rank_genes_groups` trata célula como réplica independente. Para comparar **condições** (não clusters), agregue em pseudobulk por amostra e use pydeseq2 |
| 7 | A semente está fixa e registrada? | Ver §5 |
| 8 | A árvore git estava limpa? | `runs/provenance.jsonl` grava `commit` e `arvore_limpa`; árvore suja = SHA é só dica |

O item 6 é o erro mais comum em single-cell e o mais caro: com 5.000 células de
3 camundongos, o n não é 5.000 — é 3. Tratar célula como réplica infla o poder
estatístico por ordens de grandeza e produz listas de marcadores cheias de
falso-positivo com `padj` minúsculo.

---

## 7. O que esta stack ainda não tem

Honestidade sobre limites, para não afirmar mais do que dá:

- **Não há dataset.** `data/raw/` e `data/proc/` estão vazios. Tudo acima é
  maquinário validado contra nenhuma entrada. Antes de escrever qualquer frase
  de Métodos, nomeie a acessão (GEO/ArrayExpress/SRA) e o arquivo exato.
- **QC a montante decide o resultado de DE** mais do que a escolha do teste —
  dublês, RNA ambiente, cortes de mitocondrial. Nada disso está automatizado
  aqui.
- **Não há validação de identificador.** Símbolo de gene inventado por um modelo
  passa despercebido. Confira contra o índice da própria matriz antes de citar
  qualquer gene em texto.
