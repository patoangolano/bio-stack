<!-- runs/vcc/submissao/LEIA-ME.md — qual arquivo é o atual, e por que os outros não são. -->

# Qual submissão é a válida

Duas candidatas, as duas aprovadas nos dois portões. **A escolha é do Мухамед e
não é automática** — a seção "v5 ou v6" abaixo diz exatamente o que pesa de cada
lado.

| arquivo | delta que carrega | contrato | numérico |
|---|---|---|---|
| **`submissao_v6.h5ad`** | `delta_calibrado_kd4_gate.npz` | **28/28** | **0,998 / 0,999 / 1,001** |
| **`submissao_v5.h5ad`** | `delta_calibrado_kd4.npz` | **28/28** | **0,998 / 0,999 / 1,001** |

Nada aqui foi submetido ao desafio. Submeter passa por `vcc_submit` do MCP em
`servers/vcc_mcp/`, que recusa amostra na fase final.

## Antes de enviar qualquer arquivo, rode os dois

```bash
./.venv/Scripts/python.exe -I scripts/vcc/validar_submissao.py <arquivo>
./.venv/Scripts/python.exe -I scripts/vcc/validar_numerico.py <arquivo> --delta <npz>
```

O primeiro testa o contrato: estrutura do CSR, rótulos, eixo, tetos. O segundo
testa se o arquivo **carrega o delta que você pensa que ele carrega**.

**O `--delta` muda com a versão, e passar o par errado agora REPROVA.** Desde
2026-10-23 o `validar_numerico` confere procedência por identidade antes de
qualquer aritmética, lendo `delta_sha256` (ou o nome, nas emissões anteriores ao
campo) do `medicao_<nome>.json` que o emissor grava. Isso existe porque a
regressão **não distinguia v5 de v6**: as duas diferem só na parcela da média,
têm correlação +0,996, e a inclinação cruzada dava 0,994 a 0,998 — dentro da
janela. Detalhe em `docs/VCC-2026-EXECUCAO.md` §26.

Os pares corretos:

```
submissao_v6.h5ad  ->  runs/vcc/transfer/delta_calibrado_kd4_gate.npz
submissao_v5.h5ad  ->  runs/vcc/transfer/delta_calibrado_kd4.npz
submissao_v4.h5ad  ->  runs/vcc/transfer/delta_previsto_v4_backup.npz
```

A v4 **não** confere contra `delta_previsto.npz`: esse arquivo foi reescrito pela
correção do gate do agrupamento (§24) e já não é o que a v4 carrega.

## v5 ou v6

As duas carregam o mesmo `k_desvio = 4`. Diferem só no agrupamento da resposta
média — correlação **+0,996**, rms da diferença 0,0061 contra 0,0675.

**A favor da v6.** O agrupamento da média da v5 tinha um defeito medido: um gate
desenhado para o desvio excluía da média os 4 contextos do Replogle, que têm a
melhor transferência medida de média do conjunto. Corrigido, o esquema ganha em
**4 de 4** configurações de teste (`pool_media_sem_mesmo_repo.json`, ganho mediano
+0,0983 contra +0,0524).

**A favor da v5.** O ganho de `k=4` foi medido no instrumento **certo para a
pergunta** — `avg_score` +0,372 contra 133× o ruído, e o `avaliador_proxy` contém
a parcela do desvio que `k` escala. O ganho do gate foi medido num instrumento
**analógico** (transferência entre contextos doadores que existem), em unidade que
não converte para `avg_score`, sobre uma diferença de 9% em rms. O proxy é **cego**
à parcela que a v6 mudou e não pode arbitrar.

Em aberto, e é o que impede uma recomendação forte: com o gate corrigido, o peso
por profundidade dá **95,5%** do agrupamento aos dois contextos do `arc_vcc_25` —
concentração maior que os 93,2% de antes, e justamente no doador cuja média é
medida como ortogonal à do Orion. O teste sem mesmo-estudo dá sinal misto entre
profundidade e peso igual (diferença de 0,0014, ruído), então o peso não foi
tocado.

## Estado de cada arquivo

| arquivo | contrato | numérico | por quê |
|---|---|---|---|
| `submissao_v6.h5ad` | 28/28 | **aprovado** | k=4 + agrupamento da média com o gate corrigido |
| `submissao_v5.h5ad` | 28/28 | **aprovado** | k=4; +0,372 de `avg_score` sobre v4 |
| `submissao_v4.h5ad` | 28/28 | aprovado | desvio sub-disperso (k=1); disparava o flag de 0,7 do `expr_mse` |
| `submissao_v3.h5ad` | 28/28 | não conferido | resposta média contaminada por viés de seleção |
| `submissao_v2.h5ad` | 28/28 | não conferido | gate do emissor zerou 4 pares |
| `submissao_v1.h5ad` | 28/28 | **+0,809 reprovado** | 143 alvos sem evidência ficaram no controle |
| `stub.h5ad` | **28/28** | **−0,001 reprovado** | delta zero; existe para testar o escritor |

A linha do `stub` é a razão de o validador numérico existir: passa nos 28 itens
sendo o **pior arquivo que este pipeline consegue produzir**. Estrutura válida não
é predição válida — e, desde §26, delta parecido não é delta certo.

## Por que v1 a v4 ficaram atrás

- **v1** — o zero da escala é a resposta **média real do contexto**, não o
  controle. Delta zero pontua abaixo de zero, e 143 dos 300 alvos estavam assim.
- **v2** — o gate do emissor zerava o delta dos pares com alvo não transcrito,
  apagando junto a resposta média. 4 de 900 pares.
- **v3** — a resposta média era 87% viés de seleção do `tested`: gene pouco
  coberto só é testado quando sua contagem sai alta.
- **v4** — correta em tudo acima, mas com o desvio no β de mínimos quadrados, que
  o encolhe a ~0,22. Pontuava `direction_fidelity` em 0,022 de um teto de 0,828.
  Não é erro de código: é o β sendo ótimo para a métrica errada. E a razão de
  orçamento do `expr_mse` dava **0,179**, abaixo do limiar de 0,7 em que o
  `cell-eval2` registra que a submissão "merece um olhar".

Detalhe em `docs/VCC-2026-EXECUCAO.md`, seções 15, 16, 19, 21 e 24 a 26.

## Um detalhe da v5 e da v6 que não é defeito, mas é extremo

`escala_de_gene_maxima` é **199,17** (v5) e **230,16** (v6), contra 3,81 na v4. É
sempre **uma** entrada de 5.559.900 — `SIN3B` → `DDIT3` — e só 46 a 55 passam de
10×; o p99,99 é 3,6. Coerente (perder o corepressor SIN3/HDAC desreprime CHOP) e
dentro de todo teto de contrato. O aumento de v5 para v6 é a correção operando: o
que entrou no agrupamento foi o Replogle, rastreio de genes essenciais com
assinatura de estresse forte.

## Por que nenhuma foi apagada

Não há pressão de disco, e a progressão v1 → v6 é o registro de cinco correções
que a seção de Métodos vai precisar citar. Apagar o erro apaga a evidência de que
ele foi encontrado — e os dois validadores, agora com checagem de procedência,
impedem que o arquivo errado seja enviado por engano.
