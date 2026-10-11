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
| 1 | ~~`wsl --install -d Ubuntu-24.04`~~ | **Feito.** Distro reinstalada em 2026-10-10 — ver §6 |
| 2 | Aprovar o marketplace `bio-stack` quando o Claude Code perguntar | Sem isso o `bio-lsp` não carrega |
| 2b | ~~Docker dentro do WSL~~ | **Feito.** `docker-ce` nativo na distro — ver §6 |
| 3 | Nomear uma acessão real (GEO/SRA/ArrayExpress) | `data/raw/` e `data/proc/` estão **vazios**. Todo o maquinário estatístico acima está validado contra nenhuma entrada |
| 4 | Indexar sua bibliografia (`library_index`) | `runs/library/` continua vazio |
| 5 | ~~Repor `S2_API_KEY`~~ | **Não é pendência.** Rastreada em 2026-10-10: nenhum código lê a variável, ela saiu das sondas, e `docs/CHAVES.md` §2 documenta o rastreamento |

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
| ~~**nf-core pipelines lint**~~ | **Refutação revertida.** Abortava por ausência do Nextflow, não por defeito próprio. Com o WSL instalado: 237 testes passados, 0 falhas |
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
| 2 | ~~Instalar o Docker dentro da distro~~ | **Feito.** `docker-ce` nativo, daemon sob systemd; cadeia Nextflow→container verificada |
| 3 | Escolher a acessão e baixar o dado | 1 dia |
| 4 | Indexar a bibliografia | ~1 h |
| 5 | Rodar `--tag seqqa2 --limit 5` do LABBench2 como linha de base | ~2 h |
| 6 | Migrar o histórico do `provenance.jsonl` para `v: 3` | ~80 linhas |
| 7 | Caminho de pseudobulk para DE entre condições | ~150 linhas |

Os itens 1–2 valem hoje. Do 3 em diante, cada um é independente — dá para parar
em qualquer ponto sem deixar nada pela metade.


---

## 6. Adendo — WSL2 instalado (2026-09-14, distro reinstalada em 2026-10-10)

Ubuntu 24.04.5 LTS, 20 núcleos, 15 GB RAM, 949 GB livres de 1007 GB, `systemd`
como PID 1. A distro chama-se `Ubuntu-24.04` e é a **única** entrada de
`wsl -l -v`.

O disco foi refeito: o `ext4.vhdx` de setembro desapareceu e o atual nasceu em
2026-10-10 — `/lost+found` carimba 02:23 daquele dia. É o que explica a
correção de arquitetura do Docker no fim desta seção. As **versões** da tabela
abaixo foram remedidas em 2026-10-10 e batem; as **verificações** da terceira
coluna são da rodada de 2026-09-14 e não foram repetidas no disco novo — o
`pipelines lint` de 237 testes, em particular, está por reconfirmar. Onde cada
binário mora importa:
`samtools`/`bcftools`/`bedtools`/`tabix` vieram do `apt` e estão em `/usr/bin`;
`seqkit` e `nextflow` são binário em `~/.local/bin`; `nf-core` veio por pipx.

Verificado ponta a ponta:

| Ferramenta | Versão | Verificação |
|---|---|---|
| `nextflow` | 26.04.6 | `nextflow run hello` completa os 4 processos |
| `nf-core` | 4.1.0 | `pipelines lint`: **237 passaram, 0 falharam** |
| `samtools` | 1.19.2 | — |
| `bcftools` | 1.19 | — |
| `bedtools` | 2.31.1 | — |
| `seqkit` | 2.13.0 | — |
| Java | OpenJDK 21.0.12.1 | pré-requisito do Nextflow |

Isso **reverte** uma linha da §3: `nf-core pipelines lint` não tinha defeito —
abortava por ausência do Nextflow. Verificação negativa também envelhece.

### Duas armadilhas que custaram tempo

**1. O `wsl.exe` engole `$var` quando chamado por `subprocess`.** Um laço
`for t in ...; do echo "$t"; done` executa o número certo de vezes e imprime
vazio em todas — sem erro, sem `stderr`, código de saída 0. Pelo Git Bash o
mesmo comando funciona. A sonda em `session_context.py` foi reescrita sem cifrão
nenhum (`command -v` aceita vários nomes de uma vez).

**2. O guard de dados barrou o instalador do Nextflow** por ser script encanado
direto para o shell — e estava certo. O script foi baixado, conferido (Apache
2.0, Seqera Labs) e só então executado. Depois barrou também uma *mensagem de
commit* que **mencionava** o padrão: a regra varre a string inteira do comando.
Contorno: `git commit -F arquivo`, que é melhor prática para texto longo de
qualquer forma.

### Docker: nativo na distro, verificado ponta a ponta

> **Correção de 2026-10-10 — mudança de arquitetura.** A redação original deste
> trecho descrevia `docker` 29.7.2 na distro falando com o daemon do **Docker
> Desktop** por integração WSL. Essa montagem não existe mais: o Docker Desktop
> **não está instalado** nesta máquina, nem há distro `docker-desktop` em
> `wsl -l -v`. Portanto não existe `docker` no PATH do Windows, e isso é o
> esperado — não lacuna a corrigir.

O que roda é **`docker-ce` 29.9.0 nativo dentro da Ubuntu-24.04**, instalado pelo
repositório oficial do Docker: daemon sob systemd (`systemctl enable --now
docker`), socket em `/var/run/docker.sock`, usuário `manec` no grupo `docker`,
20 CPUs e 15 GB de RAM visíveis. Nenhuma peça no Windows — um reparo no Desktop
não conserta nada aqui, porque não há Desktop no caminho.

Verificado com o pipeline inteiro, não com `docker info`:

```
wsl -d Ubuntu-24.04 -- bash -lc "cd ~/nf-smoke && \
  nextflow run nf-core/demo -profile test,docker --outdir resultados"
→ 8 processos em container: FASTQC x3, SEQTK_TRIM x3, COWPY, MULTIQC
→ Pipeline completed successfully (exit 0)
→ resultados/multiqc/multiqc_report.html, 2 452 732 bytes (2,45 MB)
```

Com isso a cadeia inteira está fechada: **Nextflow → container → ferramenta de
bioinformática → relatório**. `-profile docker` é o padrão daqui em diante;
conda e singularity ficam desnecessários.

`/mnt/c` atravessa a fronteira de sistema de arquivos e é ordens de grandeza
mais lento para muitos arquivos pequenos — o padrão exato de um pipeline
Nextflow. Trabalhe em `~/` dentro da distro quando o I/O importar.

## Seis repositórios de skills de terceiros (2026-10-11)

Nove plugins, 310 skills e 30 agentes entraram pelo marketplace local
`bio-stack`, no escopo do projeto. Nenhum deles instalou pelo caminho curto, e
as razões valem registro porque cada uma é uma armadilha diferente.

| Plugin | Origem | Skills | Agentes | Licença |
|---|---|---|---|---|
| `bio-scientist-crispr` | kexindon/claude-bio-scientist | 11 | 0 | MIT |
| `bioinfo-pipelines` | YuliaNuzhnenko/bioinformatics-agent-skills | 12 | 0 | MIT |
| `biomedical-ai` | zamushwani/biomedical-ai-skills | 16 | 0 | MIT |
| `arms-biostat` | koaeraser/ARMS-Biostat | 10 | 0 | **CC BY-NC 4.0** |
| `bioinfo-213` | Pavel-Kravchenko/Bioinformatics | 213 | 0 | **não declarada** |
| `bayesian-modeling` | choxos/BiostatAgent | 9 | 6 | MIT |
| `itc-modeling` | choxos/BiostatAgent | 6 | 7 | MIT |
| `r-tidy-modeling` | choxos/BiostatAgent | 23 | 10 | MIT |
| `clinical-trial-simulation` | choxos/BiostatAgent | 7 | 7 | MIT |

Duas licenças pedem atenção. O **ARMS-Biostat é CC BY-NC 4.0** — uso acadêmico
está dentro, uso comercial não. O **repositório do Pavel não declara licença
nenhuma**, então o padrão legal é "todos os direitos reservados": serve como
referência de leitura, não para redistribuir.

### Por que não `claude plugin marketplace add <repo>`

Nenhum dos seis entrou por esse caminho:

- quatro não têm `.claude-plugin/marketplace.json`, e o `add` responde
  `Marketplace file not found`;
- o do Pavel tem um arquivo com esse nome, mas declara só
  `{name, description, version, skills}` — sem `owner` nem `plugins`, não valida
  como marketplace;
- o do choxos valida e os plugins validam, mas carregam com erro:
  `conflicting manifests: both plugin.json and marketplace entry specify
  components`. O `marketplace.json` dele enumera agentes, comandos e skills com
  `strict: false`, e cada plugin tem o próprio `plugin.json`. Declarado no
  marketplace local como ponteiro, sem componentes, ele carrega limpo.

### Três armadilhas medidas

1. **`enable` não instala.** `claude plugin enable` liga o interruptor de um
   plugin já instalado; sem `claude plugin install` antes, o plugin não aparece
   em `claude plugin list` e nada carrega — sem erro, só ausência.

2. **O caminho da entrada não pode sair do marketplace, nem por junção.** A
   recusa é explícita: `an absolute, climbing, network-shaped,
   backslash-containing or link-traversing entry`. Um `vendor/` na raiz do
   projeto com `source: ../../vendor/x` falha, e um `mklink /J` apontando para
   lá falha igual. Por isso os clones moram dentro de
   `.claude/marketplace/plugins/`, ignorados pelo git.

3. **Enumerar skills soma com a descoberta por convenção.** Declarar as 213 do
   Pavel na entrada do marketplace *e* manter o diretório `Skills/` fez o
   `details` reportar **426**. As entradas deste marketplace não enumeram
   componentes de propósito.

E uma quarta, específica do Windows: o repositório do Pavel tem
`Legacy/.../Thumbs.db:encryptable` no índice. O dois-pontos é sintaxe de
alternate data stream no NTFS, e o git aborta a operação inteira com
`invalid path` — `clone` deixa a árvore vazia com exit 0 e sem aviso, `checkout`
falha, `sparse-checkout` falha porque a validação roda antes do filtro, e
`git archive` falha também. Só passa com `core.protectNTFS=false`, e é o que o
`extrair_pavel.sh` faz. Como o archive leva apenas `Skills/`, o arquivo
problemático nunca é materializado — e de brinde os 496 MB caem para 2,5 MB,
porque o peso era `Course/` e `Legacy/`, material de curso sem relação com
skill nenhuma.

### Custo em contexto

O `details` projeta o custo sempre-ligado de cada plugin. Somados, os nove
custam **~38,5 k tokens em toda sessão**, e **27 k são do `bioinfo-213`** — as
213 descrições de skill, das quais 29 são de algoritmo (`algo-avl-trees`,
`algo-dijkstra`) e 31 de Python, sem relação com os dois projetos. Para cortar
esse pedaço sem perder as de bioinformática:

```bash
cd .claude/marketplace/plugins/bioinfo-213
rm -rf Skills/algo-* Skills/python-*
cd ../../../.. && claude plugin uninstall bioinfo-213@bio-stack --scope project \
  && claude plugin install bioinfo-213@bio-stack --scope project
```

### Reconstrução

Os clones não são versionados. Para refazer tudo:

```bash
bash .claude/marketplace/plugins/clonar.sh
claude plugin marketplace update bio-stack
```
