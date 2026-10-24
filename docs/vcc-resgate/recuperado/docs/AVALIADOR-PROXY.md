# `avaliador_proxy` — o avaliador local do Virtual Cell Challenge 2026

Etapa 5 de [`VCC-2026-ROTA.md`](VCC-2026-ROTA.md). Monta e roda a sequência de
pontuação do `cell-eval2` 0.19.0 de ponta a ponta sobre dado nosso, para que
cada mudança no modelo tenha um número antes de virar submissão.

Por que ele existe: o placar do desafio fica oculto até 2026-11-05 e aceita **2
submissões por dia**. Sem avaliador local, toda decisão de modelagem entre hoje
e o prazo é aposta cega.

---

## 1. Montar

Uma linha, do Windows:

```bash
wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/envs/avaliador-proxy/montar.sh
```

O venv vai para `~/envs/avaliador-proxy/.venv` **dentro da distro**, não em
`/mnt/c`. As versões exatas estão em
[`envs/avaliador-proxy/requirements.in`](../envs/avaliador-proxy/requirements.in),
que também explica cada escolha; o `requirements.txt` ao lado registra o que
entrou de fato, transitivas incluídas.

O `montar.sh` termina com uma conferência que **importa de verdade** em vez de
chamar `--help`: resolve o preset, resolve o backend de DE e calcula o gate de
CPM num exemplo de duas células. O `--help` não importa `polars`, não resolve
backend e não toca o gate — passa num ambiente que não pontuaria nada.

### O que ficou de fora, e por quê

**`gpudge` não entra.** É o motor de DE de GPU que a Arc usou nos bundles
oficiais `-r2`, e arrasta `torch>=2.5`. Três razões:

1. `de.backend` tem de ser fixado em `pdex` de todo modo (seção 3). Com o
   `backend: auto` do preset, o `cell-eval2` **levanta exceção** num host onde
   há dispositivo CUDA visível e nenhum motor de DE de GPU instalado — calcular
   em CPU calado ali moveria todo número de DE sem sinal.
2. `torch` só tem wheel para `sm_120` (Blackwell) em `cu128+`, e o `gpudge` não
   declara índice: um wheel de índice default instala calado e falha no primeiro
   kernel.
3. Os 8 GB de VRAM da RTX 5070 Laptop não sustentam DE de 18 533 genes × 138 400
   células.

**Consequência declarada:** `pdex` e `gpudge` não são idênticos bit a bit. O
README do `gpudge` promete Pearson > 0,999999999 e o código do `cell-eval2` diz
que os números de DE diferem ~1e-5. Nossa pontuação local é aproximação da do
servidor, e essa é a maior fonte de erro conhecida — maior que fp32 vs fp64 em
`pds_cosine`.

---

## 2. Rodar

```bash
wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/vcc/avaliador_corrida.sh
```

Variáveis de ambiente, todas com default:
`CONTEXTO=D`, `N_PERTS=32`, `N_CTRL=18400`, `SEMENTE=0`,
`TRABALHO=$HOME/vcc-avaliador`.

O script faz cinco coisas: monta o par (real, predito) e roda os quatro passos
do `cell-eval2`, cronometrando cada um com `/usr/bin/time -f '%M %e'`.

### A sequência exata, verbatim

Os quatro passos, e os quatro com o **mesmo preset**, o **mesmo `--pert-col`** e
os **mesmos `--set`**. Divergência em qualquer um dos três aborta no passo 4.

```bash
# passo 0 — o par (real, predito). Não faz parte do cell-eval2.
python -I scripts/vcc/avaliador_proxy.py \
    --contexto D --n-perts 32 --n-ctrl 18400 --semente 0 --saida "$PROXY"

# passo 1 — a ponta 0 da escala: o braço de resposta genérica.
cell-eval2 baseline \
    -ar "$PROXY/real.h5ad" \
    --preset vcc2026 --pert-col target_gene --set de.backend=pdex \
    --save-pred "$OUT/baseline_pred.h5ad" \
    -o "$OUT/baseline"

# passo 2 — as duas pontas e a identidade que as liga.
cell-eval2 prep-real-bundle \
    --real "$PROXY/real.h5ad" --baseline "$OUT/baseline_pred.h5ad" \
    --preset vcc2026 --pert-col target_gene --set de.backend=pdex \
    --id proxy-D-k32 -o "$OUT/bundle" --force

# passo 3 — a submissão.
cell-eval2 run \
    -ap "$PROXY/pred.h5ad" -ar "$PROXY/real.h5ad" \
    --preset vcc2026 --pert-col target_gene --set de.backend=pdex \
    -o "$OUT/run"

# passo 4 — a pontuação.
cell-eval2 score \
    --user-agg "$OUT/run/agg_results.csv" --real-bundle "$OUT/bundle" \
    -o "$OUT/scored.csv"
```

O passo 2 é o caro — cinco corridas completas de métrica, uma por split da
âncora, mais a do braço de baseline — e o único reaproveitável: **enquanto o
lado real não muda, o bundle vale**. O ciclo de iteração de modelo é
passo 3 → passo 4, não os quatro.

### O que uma corrida devolve

```
metric,from_baseline,from_replicate,...
de_wilcoxon_lfc_nmae,...
expr_mse_unbiased_capped_norm,...
de_wilcoxon_direction_fidelity_yield_raw,...
de_wilcoxon_direction_reach_raw,...
de_wilcoxon_sig_jaccard,...
pds_cosine,...
avg_score,,<a pontuação>
```

`from_baseline` **nulo** em `avg_score` é a marca de uma submissão
**matriculada** contra um bundle de competição: o `score` anula a média do
baseline de propósito, para que ninguém leia um número não-oficial sob o rótulo
que todo mundo lê. Se `from_baseline` vier preenchido em `avg_score`, a corrida
foi rebaixada a diagnóstico — procure a razão antes de acreditar no número.

Os números medidos estão em
[`runs/vcc/avaliador/verificacao.json`](../runs/vcc/avaliador/verificacao.json).

---

## 3. Os oito campos de identidade — e os quatro que não são oito

O `score --real-bundle` é uma afirmação: *isto é uma submissão de competição*.
Por isso compara oito campos entre o `manifest.json` do bundle e o
`run_meta.json` da corrida, e **qualquer divergência é fatal**. Não há waiver:
`--allow-config-mismatch` ao lado de `--real-bundle` é **erro de uso**, não
escape.

Medido com [`scripts/vcc/avaliador_peers.py`](../scripts/vcc/avaliador_peers.py),
que chama `baseline.build_run_meta`, `real_bundle.check_submission` e
`score.expect_from_run_meta` — as mesmas funções que o `run` e o `score`
chamam — e portanto dá o mesmo veredito sem a conta de CPU.

| campo | o que ele é | o que o faz divergir |
|---|---|---|
| `config_digest` | sha256 sobre a config de pontuação **requisitada**, mais a estatística de agregação e a normalização resolvida de cada métrica | **quase tudo.** Medido: `de.p_adj_threshold`, `de.min_abs_log2fc`, `de.fdr_scope`, `de.backend`, `filter.filter_gene_min_cpm_cell`, `bulk_target_sum`, `max_counts_per_cell`, `discrimination.exclude_target_gene`, `discrimination.tie_policy`, `--profile`, `--input-type`. **Não** divergem: `num_threads`, `cache_real`, `cache_pred`, `cache_strict`, `outdir`, `gather_threads`, `allow_fractional_counts` (ver a seção 3.1) |
| `comparator` | o espaço em que as métricas de expressão são computadas: `bulk_lognorm` ou `lognorm` | asimetria de `input_type` entre os dois lados, ou `version: v1`. Sob `vcc2026` (counts nos dois lados, `autodetect_input_type: false`) é a constante `bulk_lognorm` — só sai dela quem mudar `--input-type` ou `--version` |
| `source_fingerprint` | hash do lado **real**, de **conteúdo** sob `cache_strict` (que o preset liga) | qualquer byte diferente na matriz do real. Medido com dois `real.h5ad` de mesma forma, mesmos rótulos, mesmos genes e sementes diferentes: diverge. **É o campo que impede pontuar uma submissão contra o bundle de outro dataset** |
| `source_fingerprint_strict` | se o campo acima é hash de conteúdo ou de metadado | `--no-cache-strict`. ⚠️ Não basta que os dois lados concordem: `expect_from_run_meta` exige **True**, e recusa por nome quando é falso — dois artefatos ambos não-estritos concordariam entre si e seriam recusados de todo modo |
| `resolved_device` | `cpu` ou `cuda`, com `device: auto` resolvido no host | **o host.** `device="auto"` e `device="cpu"` dão o mesmo valor numa máquina sem cupy — medido: ambos `'cpu'` — então escrever um ou outro não muda nada aqui. O que muda é a máquina: um host com cupy resolveria `cuda` e o bundle feito nele não casaria com uma corrida feita nesta. Único dos oito que nenhuma variante de config move |
| `resolved_de_backend` | o motor de DE, com `auto` resolvido no host | `--set de.backend=<outro>`; e, mais sutil, **fornecer `--de-real`**: `_de_backend_used` devolve falso quando a tabela do lado real é fornecida, e o campo vira `null` em vez de `'pdex'` |
| `input_type_real_effective` | `counts` ou `lognorm` para o lado real, depois da resolução | `--input-type lognorm`. Sob `vcc2026` é constante, porque o preset desliga `autodetect_input_type` |
| `de_real_fingerprint` | hash da tabela de DE do lado real **quando fornecida**; `null` quando computada | `run --de-real <parquet>`. No lado do bundle é **sempre** `null`: o `prep-real-bundle` **recusa** `--de-real` de saída, porque o braço de baseline usaria a tabela fornecida enquanto a âncora recomputa a sua, e as duas pontas de `de_wilcoxon_lfc_nmae` sairiam de tabelas diferentes |

### 3.1. Os outros quatro portões, que não estão nos oito

`score.expect_from_run_meta` e `check_submission` acrescentam quatro verificações
que abortam do mesmo jeito. A primeira é a armadilha de verdade.

1. **`allow_fractional_counts` move a identidade semântica da âncora e NÃO move
   o `config_digest`.** Ele está em `DIGEST_EXEMPT_FIELDS` e está em
   `anchor._SEMANTIC_FIELDS`. Medido: com `allow_fractional_counts=true`, os
   **oito peers batem todos** e `check_submission` **aceita** — e
   `anchor_semantic_identity` diverge. É o único campo encontrado que é
   invisível aos oito e fatal no caminho do bundle. O vizinho dele na mesma
   tupla, `max_counts_per_cell`, move os dois; a isenção é exceção, não regra.
2. **`anchor_semantic_identity`** — hash sobre `_SEMANTIC_FIELDS` (`version`,
   `input_type`, `target_sum`, `bulk_target_sum`, `allow_discrete`,
   `autodetect_input_type`, `allow_fractional_counts`, `max_counts_per_cell`,
   `validate_input`, `pert_col`, `control`) mais `target_gene_map` e os campos de
   `discrimination`. Praticamente toda variante que move o `config_digest` move
   este também, o que torna a armadilha 1 a única assimetria que importa.
3. **`anchor_metric_names`** — o conjunto de métricas que a âncora cobre. Move
   com `--profile`.
4. **`de_pred_fingerprint`**, checado de **um lado só**, contra a constante
   `null`. Medido: com `run --de-pred <parquet>` os oito peers batem todos e a
   submissão é **recusada** de todo modo. O opt-in é
   `score --diagnostic-supplied-de-pred`, e ele compra a corrida, não o rótulo:
   o resultado não é matriculado.

### 3.2. E duas coisas que abortam antes de chegar aqui

- **`--pert-col` errado morre mais cedo e melhor.** Medido: com
  `--pert-col target` (o default do preset) contra um painel cuja coluna é
  `target_gene`, a leitura do `obs` aborta com
  `perturbation column 'target' missing from the real side's obs` **antes** de
  qualquer `run_meta.json` existir. Falhar antes é mais barato que divergir
  depois.
- **`rule_digest`.** O manifesto grava um de três estados, decidido uma vez na
  construção: igual ao `competition_digest()` deste build (bundle de
  **competição**, matricula), `null` (bundle **diagnóstico**, e
  `rule_mismatches` lista as razões), ou presente e diferente (bundle feito sob
  uma regra que mudou — o `score` **recusa**). Confira
  `rule_mismatches: []` no manifesto antes de confiar num número.

### 3.3. A conclusão operacional

**Um bundle local nunca certifica contra um bundle oficial da Arc, e não é
defeito.** Os bundles oficiais `-r2` foram feitos com `gpudge`, logo o
`resolved_de_backend` deles é `'gpudge'` e o `config_digest` carrega esse valor;
o nosso é `'pdex'`. Os dois lados do nosso bundle são nossos, então eles casam
entre si — e é isso que dá o ciclo de iteração. O que o avaliador local mede é
**diferença entre duas versões do nosso modelo**, não o número que o servidor
vai imprimir.

Receita prática para não perder cota:

1. Monte o bundle **uma vez** por lado real e guarde-o.
2. Rode `run` + `score` a cada iteração do modelo.
3. Antes de submeter, confira que `rule_mismatches` está vazio e que
   `from_baseline` de `avg_score` veio **nulo**.

---

## 4. O que é real no proxy e o que é proxy

[`scripts/vcc/avaliador_proxy.py`](../scripts/vcc/avaliador_proxy.py) monta o par
que a sequência consome, porque os bundles de D/E/F estão ocultos.

**Real:**

- **As células.** Todas saem das 18 400 células de controle entregues em
  `data/raw/vcc/controls-final/context_{D,E,F}.h5ad` — contagens medidas.
  Distribuição de biblioteca, esparsidade (5 678 nnz/célula medido em D) e
  correlação entre genes vêm do dado.
- **O efeito de perturbação.** Razões de dobra de `fold_changes.npy` de um dos
  datasets preprocessados do PIE, medidas pela Arc. **A escolha da fonte decide
  quantas perturbações podem levar rótulo real do painel**, e a diferença é
  grande — medido em 2026-10-22:

  | fonte / contexto | linhas | alvos entre os 300 do VCC | genes que cruzam nosso eixo |
  |---|---:|---:|---:|
  | `arc_vcc_25` / `ARC_H1` | 146 | **12** | 12 286 de 18 533 |
  | `arc_vcc_25` / `ARC_H1_VAL` | 44 | 5 | idem |
  | `orion` / `hct116` | 8 987 | **148** | 14 690 de 18 533 |
  | `orion` / `hek293t` | 4 297 | 14 | idem |

  Os números 12 e 148 são os mesmos que o memorando da rota tirou do
  `train.json` do PIE, por um caminho independente. Para recortar uma submissão
  real contra o lado proxy os rótulos têm de existir na submissão, logo
  `--fonte-efeitos orion --so-alvos-vcc` é a única combinação que passa de uma
  dúzia de perturbações.

**Proxy, declarado:**

- **O emparelhamento contexto-efeito.** O efeito medido em `ARC_H1` (H1
  embrionária) é aplicado sobre células de D/E/F, que são outro tipo celular.
  Não é predição biológica — é um gerador de dado com a estatística certa para o
  avaliador morder.
- **O eixo de genes.** `arc_vcc_25` tem 12 286 genes, todos presentes nos nossos
  18 533; os outros 6 247 entram com razão de dobra 1, assim como o que a Arc
  marcou `tested = False`.
- **O lado predito.** O mesmo mecanismo com o efeito encolhido por um expoente
  (`--encolhimento 0.5`) e com jitter lognormal por gene
  (`--jitter 0.08`): direção majoritariamente certa, magnitude errada. É o regime
  em que um modelo real cai, e o único que produz as seis métricas fora dos
  extremos. O encolhimento sozinho seria o real reescalado, e `delta_pearson`
  e `pds_cosine` não distinguiriam isso de um modelo perfeito com escala errada.

### Armadilhas medidas na construção do proxy

1. **O lado predito PRECISA ter linhas `non-targeting`, e isso contraria o
   contrato de submissão.** `cell_eval2.io.validate_pair` exige conjunto de
   perturbações **idêntico** nos dois lados e aborta com
   `perturbation sets differ: 4 pred vs 5 real; real-only=['non-targeting']`. O
   arquivo `.vcc` que se envia ao desafio **não pode** ter linha
   `non-targeting`; o h5ad que chega ao `run` **tem** — o passo de preparação da
   Arc reinjeta os controles reais na predição, e é esse artefato que o avaliador
   lê. O conteúdo dessas linhas é inerte no caminho de CPU: com
   `control_source: real` (o que o preset fixa), `run._pred_de_input` **descarta**
   as linhas de controle da predição e dobra o pool do lado **real** como
   referência.
2. **A fonte de efeitos limita quantos rótulos do painel o proxy pode ter.** Ver
   a tabela acima. Sem `--so-alvos-vcc` a seleção completa com nocautes que não
   são alvos do painel, o que não prejudica o avaliador — ele mede a **forma** do
   dado, não a identidade dos genes — mas **impede recortar uma submissão real
   contra esse lado**, porque os rótulos não existem na submissão. Medido: com
   `arc_vcc_25/ARC_H1 --so-alvos-vcc` e o filtro de ≥200 genes moventes, só **11**
   dos 12 passam; com `orion/hct116`, 40 passam folgado.
3. **Com `--n-ctrl 18400` não sobra célula-fonte.** O painel entrega 18 400
   células de controle e **não** entrega as 138 400 de verdade-fundamental, então
   as células perturbadas têm de sair das mesmas 18 400. Daí `--fonte`:
   `disjunta` tira as fontes do que o pool de controle não usou (exige
   `n_ctrl < 18400`), `compartilhada` (o default) particiona as 18 400 inteiras
   de modo que o pool de controle possa ser o painel completo. Em
   `compartilhada`, cada célula perturbada é uma cópia escalada de uma célula
   que também está na referência — 400 de 18 400, ou seja 2,2% da média da
   referência. Em qualquer dos dois, os pools por perturbação são disjuntos
   **entre si**: a âncora divide o lado real em metades e o `pds_cosine` ordena
   cada perfil contra todos os outros, e duas perturbações que reamostram as
   mesmas fontes inflam as duas coisas.
4. **`pds_cosine` sai saturado por duas razões independentes, e a segunda é
   corrigível.** A primeira é o tamanho do recorte: com 4 perturbações deu
   `1,000` e com 32 também, porque discriminar entre poucas é fácil. A segunda é
   o vazamento de identidade de célula descrito na seção 5b — e essa é a que
   dominava. Com 40 perturbações e `--pred-pool-unico`, um controle puro tira
   `0,419` bruto (perto do acaso) em vez de `0,999`. Então: **sempre
   `--pred-pool-unico` para braços de referência**, e leia `pds_cosine` sabendo
   que o ranque é entre o número de perturbações do recorte, não entre 300.

---

## 5. Memória e tempo — o número que decide o desenho

**Medido** no WSL (VM de 15 GiB, 20 núcleos), pico de RSS do processo e dos
filhos por passo, com `/usr/bin/time -f '%M %e'`. Corrida `k32`: contexto D,
32 perturbações, 18 400 controles, 31 200 células no lado real, 5 698
nnz/célula.

| passo | pico RSS | tempo |
|---|---:|---:|
| `proxy` (monta o par) | 3,69 GiB | 16 s |
| `baseline` | **10,86 GiB** | 187 s |
| `prep-real-bundle` | 8,27 GiB | **841 s** |
| `run` | 8,21 GiB | 192 s |
| `score` | 0,31 GiB | 2 s |
| **total** | **10,86 GiB** | **1 238 s** (20,6 min) |

Dois fatos que o perfil por passo entrega e que uma medida agregada esconderia:

1. **O pico está no `baseline`, não no passo mais lento.** O `prep-real-bundle`
   custa 841 s — 68% do tempo — e ainda assim pica 2,6 GiB *abaixo* do
   `baseline`, porque a âncora trabalha em metades do lado real, cinco vezes.
2. **O ciclo de iteração é `run` + `score` = 194 s.** Os 1 028 s de `baseline` +
   `prep-real-bundle` são pagos uma vez por lado real. Enquanto o lado real não
   muda, o bundle vale, e testar uma versão nova do modelo custa 3,2 min.

**Dois pontos medidos, e a reta entre eles.** O pico do passo `baseline` foi
medido duas vezes, em corridas independentes:

| células no lado real | perturbações | pico RSS |
|---:|---:|---:|
| 31 200 | 32 | 10,86 GiB |
| 34 400 | 40 | 11,70 GiB |

A reta pelos dois dá **0,263 GiB por mil células** com intercepto **+2,64 GiB**.
O intercepto positivo diz que o crescimento é **compatível com linear** nessa
faixa — não superlinear — e o termo fixo de 2,6 GiB é o que o pool de 18 400
controles e a maquinaria de DE custam independentemente do número de
perturbações. Com 2 GiB de margem reservada, a reta projeta **~52 perturbações**
com os 18 400 controles, e **39 GiB** para o painel inteiro de um contexto
(138 400 células) — mais que o dobro da VM.

⚠️ **Projeção não é medida, e o ponto de OOM não foi medido.** O que está
**observado** é que 40 perturbações completam a sequência inteira com 3,3 GiB de
folga. O teto de ~52 é extrapolação de dois pontos. Uma sondagem que ia mapear o
pico em K = 8, 48 e 64 foi morta com SIGKILL antes de gravar a primeira linha, e
o `dmesg` não registrou nada que permitisse atribuir a causa com certeza; não foi
repetida porque pontuar a submissão real (seção 5b) valia mais que fechar a
curva.

**O que o número manda fazer:** avaliar em subconjunto de perturbações não é
otimização, é condição de existência. O painel oficial de um contexto tem
138 400 células de verdade-fundamental por lado; nenhuma engenharia local põe
isso num processo de 15 GiB. O ciclo de iteração roda em dezenas de
perturbações, e o modelo completo só é materializado na emissão da submissão
(etapa 6), que é streaming e nunca materializa a matriz.

Quais métricas sobrevivem a qual recorte:

- `fid` (`de_wilcoxon_direction_fidelity_yield_raw`), `reach`
  (`..._direction_reach_raw`) e `jac` (`de_wilcoxon_sig_jaccard`) exigem células
  emitidas — não há atalho por pseudobulk.
- `pds_cosine` e `expr_mse_unbiased_capped_norm` saem de pseudobulk, mas
  `pds_cosine` **depende do número de perturbações no painel** (é um ranque
  entre elas), então um recorte pequeno o infla.
- `de_wilcoxon_lfc_nmae` sai do gate de DE, que é por perturbação.

---

## 5b. O primeiro número real do projeto

**Medido** em 2026-10-22. Lado real: contexto D, **40 alvos reais do painel VCC**,
efeito medido em `orion/hct116`, 18 400 controles, 34 400 células. Um único
bundle (`proxy-D-orion-k40`, `rule_mismatches: []`) serve os cinco braços — os
dois lados reais saem byte a byte idênticos (sha256 conferido), então a diferença
entre braços é diferença de **predição** e nada mais. Os cinco braços foram
matriculados: oito peers batendo e `from_baseline` nulo em `avg_score`.

| braço | o que é | `avg_score` |
|---|---|---:|
| oráculo | o efeito real do Orion, só com ruído de amostragem | **+0,970** |
| **nossa (v4)** | recorte de `submissao_v4.h5ad`, 40 alvos do contexto D | **+0,275** |
| nossa (v3) | o mesmo recorte de `submissao_v3.h5ad` | +0,272 |
| controle puro | delta zero em todo gene | **−0,033** |

A v4 fica **30,7%** do caminho entre o piso do controle e o teto do oráculo, e
**+0,308 acima do controle**. Por membro, com a diferença v4−v3 na última coluna
— e ela é a primeira medida de iteração de modelo do projeto, porque os dois
braços passaram pelo **mesmo bundle** e pelo mesmo lado real:

| membro | controle | v3 | **v4** | oráculo | v4−v3 |
|---|---:|---:|---:|---:|---:|
| `pds_cosine` | −0,107 | 0,989 | **0,990** | 1,000 | +0,0012 |
| `expr_mse_unbiased_capped_norm` | 0,000 | 0,211 | **0,228** | 0,947 | **+0,0171** |
| `de_wilcoxon_lfc_nmae` | −0,050 | 0,098 | **0,102** | 1,024 | +0,0034 |
| `de_wilcoxon_direction_reach_raw` | −0,005 | 0,262 | **0,267** | 0,989 | +0,0057 |
| `de_wilcoxon_sig_jaccard` | −0,019 | 0,044 | **0,040** | 1,032 | −0,0047 |
| `de_wilcoxon_direction_fidelity_yield_raw` | −0,018 | 0,027 | **0,022** | 0,828 | −0,0057 |
| **`avg_score`** | **−0,033** | 0,272 | **0,275** | 0,970 | **+0,0028** |

Onde o modelo ganha e onde não ganha: `direction_reach` e `expr_mse` compram a
maior parte do ganho sobre o controle; **`direction_fidelity` é onde ele quase
não sai do controle** (+0,039 de um teto de +0,846, ou seja 5% do alcançável) — a
direção do efeito, perturbação por perturbação, é o que falta.

### ⚠️ O que o v4−v3 de +0,0028 significa, e o que não significa

**Não significa que a correção de viés do `tested` foi inútil.** Significa que
**este proxy é cego para ela**, e isso era previsível pela construção, não é
descoberta nova: a correção mexeu na parcela de **resposta média do contexto**, e
a "verdade" deste proxy é o efeito do Orion em HCT116 transplantado — ela contém
o desvio específico da perturbação e **não contém** a resposta média de D. O
proxy mede a parcela que não mudou e trata como ruído a parcela que mudou.

O que o número **sustenta**, e é afirmação útil: **a correção não custou nada
mensurável.** Zerar a média de gene com cobertura abaixo de 0,9 derrubou a
magnitude do deslocamento de 0,063 para 0,019 — uma mudança grande, do tipo que
podia ter quebrado a parcela específica — e a pontuação na parcela que o proxy
*sabe* medir ficou igual ou melhor em quatro dos seis membros. Para a pergunta
"a correção estava certa?", o diagnóstico do próprio revisor é melhor evidência
que este proxy: fração positiva de ~1,0 para 0,553, e o Replogle (que testa 100%
dos genes) com corr(cobertura, média) de **+0,040** contra −0,878 no Orion.

Os sinais por membro são mistos (+0,017 em `expr_mse`, −0,006 em `fidelity`) com
soma quase nula, o que num recorte de 40 perturbações e sem réplica é
indistinguível de ruído. **Não há intervalo de confiança aqui** — uma réplica
exigiria reconstruir o bundle com outra semente, 10 min por réplica, e não foi
feito.

### O vazamento que eu encontrei medindo isso, e que invalidou a primeira leitura

A primeira rodada deu **+0,188** para o controle puro, acima de zero, contra a
premissa de que ele pontua abaixo. A causa estava na minha construção do proxy,
não no modelo nem no avaliador: o lado real e o lado predito de uma mesma
perturbação sorteavam da **mesma partição de ~460 células-fonte**, e duas
amostras com reposição do mesmo pool se parecem mais entre si do que amostras de
outro pool. O `pds_cosine` lê esse rastro de identidade de célula como acerto.

Medido, com 40 perturbações e predição de controle puro:

| configuração | `pds_cosine` bruto | `avg_score` |
|---|---:|---:|
| partição por perturbação | 0,99936 | **+0,188** |
| pool único (`--pred-pool-unico`) | 0,41859 | **−0,033** |

O rastro valia **+1,106 de `pds_cosine` bruto** e **+0,221 de `avg_score`**, e era
ele — não biologia — que punha o controle acima de zero. Removido, o controle cai
abaixo de zero, confirmando a premissa de forma medida. Os números da tabela
principal acima já são os do regime sem rastro para as duas referências.

Duas consequências:

1. **Para braços de referência emitidos pelo `avaliador_proxy.py`, use sempre
   `--pred-pool-unico`.** Uma submissão real não precisa: o `emissor_400` sorteia
   do pool de controle inteiro, logo já está no regime sem rastro. Por isso o
   braço `nossa` não foi refeito — ele já estava correto.
2. **O `pds_cosine` da nossa predição, +0,990, é discriminação de verdade.** No
   mesmo recorte e com o mesmo bundle, o controle sem rastro tira −0,107. Então os
   deltas por perturbação do `transfer_fit` são distintivos o suficiente para
   casar cada perfil previsto com o seu próprio alvo. ⚠️ Com a ressalva de que é
   um ranque entre **40** perturbações, não entre 300 — com 300 o problema é mais
   difícil, e quanto mais difícil não está medido.

### O que este número não é

⚠️ **Não é uma estimativa do placar.** A "verdade" deste proxy é o efeito do
Orion em HCT116 transplantado para células de D. Ela contém o desvio específico
da perturbação e **não contém** a resposta média do contexto D. O delta do nosso
modelo tem as duas parcelas. Então o proxy pontua a parcela específica e é cego
para a parcela média — é a mesma razão pela qual o v4−v3 saiu em zero (acima).
Isso é raciocínio sobre a construção, não medição: **+0,275 é leitura justa da
parcela específica e nenhuma leitura da parcela média**, e o quanto a média vale
no painel real este instrumento não tem como dizer.

⚠️ **O instrumento está limitado pela VERDADE, não pela memória.** Vale dizer
porque muda a decisão de para onde levar esforço: ampliar o recorte de 40 para
300 perturbações — o que uma A100 de 167 GiB tornaria possível — compraria escala
e tornaria o `pds_cosine` interpretável em escala de painel, mas **não** compraria
validade, porque a verdade continuaria sendo Orion transplantado e só 148 dos 300
alvos têm efeito medido lá. O que destrava a parcela média é verdade-fundamental
de D, que está oculta até 2026-11-05.

Além disso, o número local é aproximação do servidor por `pdex` vs `gpudge`
(~1e-5 em DE) e as pontas `b` e `r` são as nossas, não as de D.

### Tempo e memória desta corrida

| passo | pico RSS | tempo |
|---|---:|---:|
| `proxy` ×2 (oráculo e controle) | 3,96 GiB | 21 s |
| `recorte` (fatia a submissão de 4 GB) | 6,24 GiB | 15 s |
| `baseline` | **11,70 GiB** | 120 s |
| `prep-real-bundle` | 8,81 GiB | 605 s |
| `run` (cada braço) | ~8,8 GiB | 129–191 s |
| `score` (cada braço) | 0,31 GiB | ~2 s |

**O bundle foi reaproveitado cinco vezes.** `baseline` + `prep-real-bundle`
custaram 725 s uma vez; cada braço de predição depois custou `recorte` + `run` +
`score` ≈ 209 s. Pontuar a v4 contra o mesmo lado real da v3 custou 3,5 min, e é
por isso que os dois números são comparáveis: mesmas duas pontas de escala,
mesmo `source_fingerprint`, oito campos batendo nas duas corridas.

**Segundo ponto medido de memória:** 34 400 células → 11,70 GiB, contra 31 200 →
10,86 GiB na corrida `k32`. É o segundo ponto que fecha a reta da seção 5
(0,263 GiB por mil células, intercepto +2,64 GiB). **Observado:** 40
perturbações completam a sequência inteira com 3,3 GiB de folga.

Vale notar que o `prep-real-bundle` custou **605 s com 40 perturbações** contra
841 s com 32 na corrida `k32`. O custo da âncora não escala com células — escala
com o número de genes que o gate de DE admite, que depende da fonte de efeitos.

---

## 5c. Varredura de `k_desvio` — o avaliador escolhendo um hiperparâmetro

**Medido** em 2026-10-23. Oito valores de `k_desvio` (a escala aplicada à parcela
de desvio específico do delta, em log2-razão, com a parcela média intacta),
emitidos pelo `emissor_400.py` de produção e pontuados contra o **mesmo bundle**
`proxy-D-orion-k40`. Os 11 braços deste bundle estão matriculados.

`k=1` não foi re-emitido: `delta_kd1.npz` reproduz `delta_previsto.npz` a **7e-8
relativo** (max|dif| 9,3e-10 contra rms|x| 0,013, conferido nos três contextos), e
`idx_ctx` do emissor é o índice dentro da tupla emitida — logo o bloco de D é
idêntico com `--contextos D` ou `D,E,F`, e a v4 já pontuada **é** o braço `k=1`.

| membro | k=1 | k=1,5 | k=2 | k=3 | **k=4** | k=5 | k=6 | k=8 | oráculo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `fidelity` | 0,022 | 0,129 | 0,283 | 0,589 | 0,783 | **0,841** | 0,826 | 0,790 | 0,828 |
| `jaccard` | 0,040 | 0,173 | 0,317 | 0,499 | 0,567 | **0,582** | 0,576 | 0,546 | 1,032 |
| `reach` | 0,267 | 0,441 | 0,551 | 0,667 | 0,722 | 0,750 | 0,765 | **0,770** | 0,989 |
| `lfc_nmae` | 0,102 | 0,172 | 0,239 | 0,353 | 0,427 | **0,438** | 0,378 | 0,090 | 1,024 |
| `expr_mse` | 0,228 | 0,315 | 0,380 | **0,442** | 0,392 | 0,165 | 0,000 | 0,000 | 0,947 |
| `pds_cosine` | **0,990** | 0,990 | 0,990 | 0,990 | 0,989 | 0,988 | 0,978 | 0,957 | 1,000 |
| **`avg_score`** | 0,275 | 0,370 | 0,460 | 0,590 | **0,647** | 0,627 | 0,587 | 0,526 | 0,970 |

**O ótimo é k = 4, e está entre parênteses:** 0,590 → **0,647** → 0,627 → 0,587 →
0,526. Curva unimodal e suave, com o pico entre 4 e 5. O ganho sobre `k=1` é
**+0,372 de `avg_score`** — 133 vezes o +0,0028 que separou v3 de v4, ou seja
muito acima de qualquer plausível piso de ruído deste recorte.

### A hipótese estava certa no diagnóstico e errada no mecanismo

A hipótese era que o encolhimento mata `direction` **para comprar** `expr_mse`,
isto é, uma troca. **Não é troca até k=3: tudo sobe junto.** De k=1 a k=3,
`fidelity` ×27, `jaccard` ×12,6, `reach` ×2,5, `lfc_nmae` ×3,5 — e `expr_mse`
também, de 0,228 para 0,442. O β de mínimos quadrados não estava trocando
direção por magnitude; estava **sub-disperso para os seis membros ao mesmo
tempo**, inclusive para o que ele teoricamente otimiza.

A troca existe, mas começa muito acima de 1: o ponto de virada de `expr_mse` é
k=3, o de `fidelity`/`jaccard`/`lfc_nmae` é k=5, e `reach` ainda sobe em k=8. O
`avg_score` concilia isso em k=4.

**`fidelity` em k=5 chega a 0,841 contra 0,828 do oráculo** — passa do
quase-oráculo. Não é paradoxo: o oráculo carrega o ruído de amostragem de 400
células, enquanto uma predição com desvio amplificado tem direção mais limpa que
uma medição ruidosa. Diz que `fidelity` deixou de ser o gargalo: o buraco que a
seção 5b nomeou (0,022 de um teto de 0,828, 3% do alcançável) fecha em 102% com
um único escalar.

⚠️ **`expr_mse` em 0,0000 para k=6 e k=8 é clamp, não empate com o controle.** O
valor bruto vai a 1,337 e 3,092 contra ~1,0 do baseline de resposta genérica, ou
seja fica **pior que o baseline**, e a escala corta em 0. O `avg_score` desses
dois k está portanto otimista — a queda real é mais funda que a tabela mostra.

### O que esta varredura autoriza, e o que não

**Autoriza:** `k_desvio = 4` em vez de 1, no mesmo instrumento, com o mesmo
bundle e os oito campos de identidade batendo em todos os braços. A curva é
monótona até o pico e monótona depois dele; não há sinal misto a interpretar.

**Não autoriza ler 0,647 como estimativa de placar**, pelas mesmas razões da
seção 5b — a verdade do proxy é Orion transplantado e não contém a resposta média
de D. **Mas aqui a objeção é mais fraca que lá**, e vale dizer por quê: `k` só
escala a parcela de **desvio específico**, que é exatamente a parcela que a
verdade do proxy **contém**. É a pergunta para a qual este instrumento é o
instrumento certo — ao contrário da correção de viés da média, para a qual ele
era cego.

⚠️ **Não autoriza k=4 como ótimo do painel real.** Dois motivos declarados:
40 perturbações sem réplica, logo sem intervalo de confiança; e o ótimo de uma
soma de seis membros depende das pontas `b` e `r`, que no painel real são as de
D e não as nossas. O que transfere com confiança é a **direção** (k ≫ 1) e a
**ordem de magnitude** (4, não 1 nem 40), não a terceira casa.

---

## 6. O universo de DE — a definição que reproduz o número oficial

O preset fixa `filter.filter_gene_min_cpm_cell: 5.0`, e a quantidade que o gate
compara é a **média, por gene, do CPM por célula** `x_ig / L_i · 10⁶` **sobre as
células de referência** — não o CPM do pseudobulk. Confirmado no código
(`de_compute._ref_cpm_from_cells`, cujo docstring chama a quantidade de
invariante a qualquer reescala por célula) e testado no `montar.sh` com um
exemplo de duas células, onde as duas definições dão
`[2,5e5; 7,5e5]` contra `[3,33e5; 6,67e5]`.

É a segunda definição que reproduz o número oficial: medido nesta campanha,
**D 9 633 · E 9 743 · F 9 537** genes de 18 533, contra 9 565/9 744/9 499 pela
definição de pseudobulk. O `ctrl_bulk` da etapa 2 usa a correta.

O gate é **somente de referência** (issue #351 do upstream): o `cell-eval2`
assume o gate até no caminho do `gpudge`, porque o gate nativo do `gpudge`
mantém um par (alvo, gene) quando a média do **alvo** OU a da referência passa
do limiar, e esse OR é vazamento.

---

## 7. Arquivos

| arquivo | o que é |
|---|---|
| `envs/avaliador-proxy/requirements.in` | versões exatas, com a razão de cada uma |
| `envs/avaliador-proxy/requirements.txt` | o que entrou de fato (`uv pip freeze`) |
| `envs/avaliador-proxy/montar.sh` | monta o venv e confere com import real |
| `scripts/vcc/avaliador_proxy.py` | monta o par (real, predito) |
| `scripts/vcc/avaliador_recorte.py` | recorta uma submissão real para o formato pontuável |
| `scripts/vcc/avaliador_corrida.sh` | os quatro passos, cronometrados |
| `scripts/vcc/avaliador_peers.py` | mede o que faz divergir cada campo de identidade |
| `scripts/vcc/avaliador_verificacao.py` | junta tudo em `runs/vcc/avaliador/verificacao.json` |
| `runs/vcc/avaliador/verificacao.json` | os números medidos |

O `.venv` fica fora do git (`envs/**/.venv/` no `.gitignore`), e
`runs/vcc/**` também — o `verificacao.json` é cópia de trabalho, e o que
reconstrói tudo é o `requirements.in` mais os scripts.

---

## 8. Em aberto

- **`gpudge` ou `pdex` no servidor do teste?** Os bundles `-r2` de validação
  foram `gpudge`. Para `-r3`/`-r4` e para D/E/F, não confirmado. A trava dos
  oito campos torna `gpudge` o mais provável, e nossa pontuação local é
  aproximação (~1e-5 em DE, fp32 vs fp64 em `pds_cosine` e nas de expressão).
- **`b` e `r` de D/E/F.** Estão nos bundles ocultos. As faixas publicadas são de
  A/B/C, e o `avaliador_proxy` produz as nossas, que não são as deles.
- **Os controles entregues são os de pontuação?** Indício forte
  (`manifest.json`: 138 400 células, `panel_id vcc2026-test-1`), texto ambíguo.
  Se forem outros, o `source_fingerprint` do bundle oficial não é o do nosso
  lado real — o que já é verdade pelo backend, então não muda nada aqui.
- **`resolved_device` num host com cupy.** Não testável nesta máquina sem
  instalar `cupy`, o que traria o problema do wheel `sm_120`. Documentado
  analiticamente na seção 3. É o único dos oito campos que nenhuma variante de
  config move.
- **O ponto exato de OOM.** Dois pontos medidos (31 200 → 10,86 GiB; 34 400 →
  11,70 GiB), reta compatível com linear, teto projetado em ~52 perturbações. O
  teto não foi **observado**: uma sondagem em K = 8/48/64 foi morta com SIGKILL
  sem gravar linha, e não foi repetida porque pontuar a submissão real valia mais.
- **Quanto o proxy penaliza a parcela de resposta média do contexto.** Sabe-se a
  direção (penaliza) e por quê (a verdade do proxy não contém essa parcela), mas
  não a magnitude. Medi-la exigiria um lado real que contenha a resposta média de
  D — isto é, a verdade-fundamental oculta.
- **`pds_cosine` com 300 perturbações em vez de 40.** O ranque fica mais difícil;
  quanto, não medido. Localmente não cabe em 15 GiB; numa A100 de 167 GiB
  caberia, mas ver a ressalva da seção 5b — o que falta é verdade, não memória.
- **Se a correção de viés do `tested` ajudou.** O proxy mede +0,0028 de v3 para
  v4, indistinguível de ruído, **e é cego para a parcela que mudou**. Medir isso
  exige verdade-fundamental de D. O que o proxy sustenta é o limite inferior: a
  correção não piorou a parcela específica.
- **Intervalo de confiança de qualquer um dos números.** Sem réplica. Uma réplica
  por semente custaria ~12 min (bundle novo) e não foi feita. É o que falta para
  saber se o `k=4` da seção 5c é 4 ou 4±1.
- **Se `k_desvio = 4` é o ótimo do PAINEL, não só do proxy.** O ótimo de uma soma
  de seis membros depende das pontas `b` e `r`, que no painel real são as de D.
  Transfere com confiança a direção (k ≫ 1) e a ordem de magnitude (4, não 1 nem
  40); não a terceira casa.
