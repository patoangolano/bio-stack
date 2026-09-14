# Ecossistema da stack — o que vale adotar, o que não vale, e o que já foi feito

Auditoria de 2026-09-14. Sucessora de [`AUTONOMIA.md`](AUTONOMIA.md), que
perguntou "o que falta para a stack rodar sozinha?". Esta pergunta é outra:
**quais MCPs, skills, plugins, comandos, hooks e LSPs realmente melhoram esta
stack** — e quais só parecem melhorar.

## Como isto foi apurado

Dois workflows de pesquisa, 197 agentes no total, com verificação adversarial de
3 votos por alegação (cada verificador instruído a *refutar*, com refutação como
padrão na dúvida). Depois, verificação local minha de tudo que era verificável
nesta máquina.

O resultado honesto da verificação:

| | alegações extraídas | sobreviveram |
|---|---|---|
| Workflow 1 — ecossistema MCP e mecânica do Claude Code | 124 → 25 testadas | 20 |
| Workflow 2 — reprodutibilidade e rigor estatístico | 24 | **4** |

Seis dos oito sub-ângulos do segundo workflow **não produziram nada aproveitável**.
Isso não é resultado negativo — é ausência de cobertura, e está registrado como
tal na §4. O que sobreviveu, sobreviveu porque foi conferido contra o código
instalado nesta venv, não contra documentação.

> Regra que emergiu e vale manter: **documentação descreve intenção, código
> descreve o que roda.** Três das quatro descobertas mais úteis vieram de ler
> `.venv/Lib/site-packages/`, e contradiziam a leitura natural da doc.

---

## 1. O que já foi executado

Três commits. Cada item abaixo foi verificado na máquina, não presumido.

### 1.1 `git init` — o caderno de laboratório passou a ter história

`runs/provenance.jsonl` registrava 1.500+ operações sem que houvesse como saber
o que o código *era* quando cada uma rodou. `.gitattributes` traz
`merge=union` nesse arquivo: append-only, duas sessões registrando em paralelo
não geram conflito.

### 1.2 O guard de dados tinha um buraco do tamanho do problema

`guard_data.py` casava só `Bash|PowerShell`. **`Write` e `Edit` sobre
`data/raw/` passavam sem checagem nenhuma** — o caminho mais curto para destruir
dado bruto era o único desprotegido. E o `CLAUDE.md` afirmava que `data/raw`
ficava somente-leitura; `attrib` mostrou que não havia atributo nenhum.

Corrigido, com 11 casos em `tests/test_guard_write.py` que rodam o hook como
subprocesso de verdade. Durante o próprio teste o guard bloqueou meu comando —
prova empírica de que `exit 2` bloqueia nesta versão.

### 1.3 O log só registrava sucesso

`PostToolUse` não dispara quando a ferramenta falha. Consequência: **uma seção
de Métodos gerada desse log descreve um experimento que não aconteceu**, e o
viés é invisível porque nada no arquivo denuncia a falta.

Pior: a exceção que pulava os servidores `nvidia-*` (porque "eles já registram
sozinhos") não valia para falha. Todo `_log_provenance()` em
`servers/nvidia_bio_mcp/server.py` roda **depois** de obter o resultado — se o
POST levanta exceção, a linha nunca é alcançada. Foi assim que o 422 crônico do
`evo2_score` passou meses sem aparecer em lugar nenhum.

Agora `PostToolUseFailure` também registra, com campo `ok` e a mensagem de erro.

### 1.4 Resultado agora está amarrado a código

Cada linha do log carrega o SHA do `HEAD`, lido direto de `.git/HEAD` (um
subprocesso por chamada custaria minutos em 900 operações). Os limites de sessão
gravam `git status`: **árvore suja significa que o SHA é só uma dica**, porque o
código que rodou não está em nenhum commit.

`session_mark.py` marca também compactação de contexto — sem isso `/metodos` não
tem como saber que houve descontinuidade no meio de uma análise.

### 1.5 O ambiente não estava travado

`pyproject.toml` não tinha `[project]`, e `uv.lock` tinha **52 bytes e zero
pacotes**. `uv sync` noutra máquina instalava nada. Agora: 173 pacotes travados,
com versão exata nos que carregam peso estatístico — trocar `pydeseq2` ou
`scanpy` muda o que a seção de Métodos *afirma*, então não podem flutuar.

### 1.6 Skill `rigor-estatistico` — quatro fatos que mudam conclusão

Lidos do código instalado. Cada um traz o comando que reproduz a evidência.

| Fato | Onde | Por que morde |
|---|---|---|
| `padj = NaN` significa "filtrado fora", nunca "não significativo" | `pydeseq2/ds.py:510,517,554` | A contagem de DEG sai **certa** e o denominador sai **errado**. "De 18.000 genes, 420 são DE" — os 18.000 nunca foram testados |
| `padj` vem de `scipy.stats.false_discovery_control`, não do statsmodels | `ds.py:10,544` | Rodar `multipletests` em cima é correção dupla |
| `lfc_shrink()` deixa o p-valor intacto | docstring do próprio método | Filtrar `abs(log2FC) > 1` depois de encolher **não é teste**. O certo é `lfc_null=1.0, alt_hypothesis="greaterAbs"` |
| O FDR do scanpy é **por grupo** | `_rank_genes_groups.py:442,468` | Com 12 clusters são 12 correções independentes; o FDR global é maior que o α declarado. Com `method='logreg'` não existe p-valor nenhum |

### 1.7 `np.random.seed()` virou erro de lint

Regra `NPY002` ativa. É estado global de processo: qualquer biblioteca importada
depois consome do mesmo gerador e desloca a sequência **sem falhar**. O
`CLAUDE.md` recomendava exatamente isso — corrigido no mesmo commit.

A regra correta, confirmada por teste: `random_state=<int>` em chamada de
biblioteca, `np.random.default_rng()` só em código próprio. O scanpy tipa
`random_state` como `_LegacyRandom` e **recusa** um `Generator` moderno — o
basedpyright rejeita a linha estaticamente, e `tests/test_determinism.py`
confirma em runtime (5/5, incluindo o controle negativo: semente diferente
*precisa* mudar o embedding, senão o teste não estaria testando nada).

### 1.8 Cache endereçado por conteúdo no servidor NIM

Inserido no `_post`, por onde toda ferramenta passa. Varredura de mutantes com
prefixo compartilhado para de recomprar MSA. Desligável com
`BIO_STACK_NIM_CACHE=0`; `nim_config` reporta entradas e bytes. Cache-hit é
registrado na proveniência — reaproveitar não pode virar "recalculei".
15 casos em `tests/test_cache.py`.

### 1.9 ToolUniverse ligado, e o motivo do plugin falhar

Medido nesta máquina: `uvx tooluniverse` sobe a v3.4.7, carrega **2.621
ferramentas** e expõe **4** — `Compact mode: Exposed 4 core discovery tools`.
O handshake MCP responde JSON limpo no stdout.

Isso derruba a objeção óbvia ("2.700 ferramentas estouram o contexto"): a
separação entre *descobrir* e *chamar* é o desenho do servidor. Ligado direto no
`.mcp.json`. O `plugin:tooluniverse` continua dando `CONNECTION_CLOSED` — é a
fiação do plugin, não o pacote.

### 1.10 LSP: o plugin oficial não serviria

`settings.json` **não aceita** `lspServers` — conferido contra o schema oficial,
que só tem o campo em `plugin.json`. E o `pyright-lsp` do marketplace chama
`pyright-langserver`, que não existe aqui; o que existe é
`basedpyright-langserver.exe`, dentro da `.venv`.

Criado `.claude/marketplace/` no próprio repositório declarando `basedpyright` e
`ruff server`. Ambos respondem LSP — verificado com handshake `initialize`
(basedpyright 1.40.1, ruff 0.16.7).

---

## 2. O que depende de você

| # | Ação | Por quê |
|---|---|---|
| 1 | `wsl --install -d Ubuntu-24.04` | Destrava `nextflow`, `samtools`, `seqkit`, `pysam`. O kernel já está lá — falta só a distro, sem reinício |
| 2 | Aprovar o marketplace `bio-stack` quando o Claude Code perguntar | Sem isso o `bio-lsp` não carrega |
| 3 | Nomear uma acessão real (GEO/SRA/ArrayExpress) | `data/raw/` e `data/proc/` estão **vazios**. Todo o maquinário estatístico acima está validado contra nenhuma entrada |
| 4 | Indexar sua bibliografia (`library_index`) | `runs/library/` continua vazio |
| 5 | Repor `S2_API_KEY` | Semantic Scholar ausente |

### Plugins que poluem o contexto

Estão carregados `small-business`, `legal`, `marketing`, `sales`,
`human-resources`, `finance`, `figma`, `postiz`, `customer-support` — dezenas de
skills e servidores MCP irrelevantes, em toda sessão. Desligar o que não é de
bioinformática é o ganho de contexto mais barato disponível.

### Benchmark

**LABBench2** (EdisonScientific — a linha do FutureHouse mudou de organização,
então "o benchmark morreu" é falso). CI declara `windows-latest`, Python
3.11–3.13, sem Docker e sem GPU. 1.912 tarefas em 15 tags. Último commit
2026-05-07, CC BY-SA 4.0.

```bash
git clone --filter=blob:none --depth 1 https://github.com/EdisonScientific/labbench2
```

Ressalva real: o `uv sync` dentro do clone monta uma venv **separada**, não
reaproveita a da bio-stack — e o dataset é gated no HuggingFace.

---

## 3. O que NÃO adotar — e por quê

Resultado negativo verificado vale tanto quanto positivo, e economiza dias.

| Candidato | Veredito |
|---|---|
| Família **Augmented-Nature** (24 servidores, um por base: UniProt, PDB, KEGG, STRING…) | Quase dormente. Todos com o mesmo carimbo de 2025-12-21 de uma varredura de LICENSE; o de UniProt teve última mudança real em ~05/2025. Nenhum publicado no npm, todos presos ao SDK MCP pré-1.0, todos auto-rotulados "Unofficial". Contagem de ferramentas no README contradiz o próprio `package.json` |
| **scmcp** (single-cell) | Último push 2025-07-10; a organização inteira sem atividade em 2026. A lacuna de single-cell é real e **continua aberta** |
| **miniwdl**, **cwltool** | Quebram no import em Windows nativo (`fcntl`, `pwd`) |
| **nf-core pipelines lint** | Aborta antes de selecionar teste — chama o Nextflow por baixo |
| **zenodo_client** | `ensure_zenodo` **não é idempotente**: cada disparo publica uma nova versão com novo DOI. O oposto do que um hook precisa |
| **multipy** (q-valor de Storey) | Morto: último commit 2020-10-27. Use `multipletests(method="fdr_tsbh")`, já instalado — e chame de BH adaptativo de dois estágios, não de "q-valor de Storey" |
| Registro do **BioContextAI** como filtro de qualidade | A alegação de que a inclusão garante manutenção foi refutada 0-3. Serve para descobrir, nunca para confiar |

O que **roda** em Windows nativo: `snakemake --lint` e `snakefmt`.

---

## 4. O que continua em aberto

Honestidade sobre os limites desta auditoria.

1. **Seis dos oito sub-ângulos do segundo workflow não produziram nada.**
   RO-Crate/PROV, deposição Zenodo/ORCID, versionamento de dados, modos de falha
   publicados, conformidade de workflow e padrões de reporte ficaram **sem
   resposta** — não com resposta negativa. Se algum importa para o TCC, precisa
   de pesquisa própria.
2. **`runs/provenance.jsonl` mudou de forma no meio.** Linhas antigas não têm
   `ok`, `commit` nem `v`. Agora há `"v": 3`, mas ninguém migrou o histórico, e
   o arquivo mistura análise real com ruído de shell sob nenhuma marcação que
   os separe. `/metodos` exige curadoria manual que não fica registrada.
3. **Pseudorreplicação não está coberta por nada.** `rank_genes_groups` trata
   célula como réplica independente. Com 5.000 células de 3 camundongos, o n é
   3, não 5.000. O caminho correto — pseudobulk por amostra e depois pydeseq2 —
   não foi avaliado.
4. **QC a montante decide o resultado de DE** mais que a escolha do teste
   (dublês, RNA ambiente, cortes de mitocondrial). Nada automatizado aqui.
5. **Não há validação de identificador.** Símbolo de gene inventado por um
   modelo passa despercebido.

---

## 5. Ordem sugerida

| # | Ação | Custo |
|---|---|---|
| 1 | Desligar os plugins fora de escopo | minutos |
| 2 | `wsl --install -d Ubuntu-24.04` | ~15 min |
| 3 | Escolher a acessão e baixar o dado | 1 dia |
| 4 | Indexar a bibliografia | ~1 h |
| 5 | Rodar `--tag seqqa2 --limit 5` do LABBench2 como linha de base | ~2 h |
| 6 | Migrar o histórico do `provenance.jsonl` para `v: 3` | ~80 linhas |
| 7 | Caminho de pseudobulk para DE entre condições | ~150 linhas |

Os itens 1–2 valem hoje. Do 3 em diante, cada um é independente — dá para parar
em qualquer ponto sem deixar nada pela metade.
