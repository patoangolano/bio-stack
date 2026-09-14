# O caminho para um ambiente autônomo

Auditoria de 2026-09-12. O que está aqui ou é medido nesta máquina, ou vem de
pesquisa com fonte citada. O que é opinião está marcado como tal.

---

## O diagnóstico honesto

A stack hoje é uma **bancada bem equipada**, não um laboratório autônomo. A
diferença não é quantidade de ferramenta — é que falta o que transforma uma
sessão de trabalho em um processo que continua sem você olhando.

Três evidências de que o problema não é falta de ferramenta:

1. **164 operações no `runs/provenance.jsonl`, 123 delas `Bash`.** As ferramentas
   especializadas — os NIMs, o `gget`, o `biomcp` — respondem por 10. A stack
   tem mais capacidade do que uso.
2. **`runs/library/` está vazio.** A camada de RAG sobre a sua bibliografia
   nunca rodou. Descobrimos por quê: o modelo de embedding estava morto havia
   quatro meses (HTTP 410). Ninguém notou porque ninguém dependia dela.
3. **Três rotas quebradas silenciosamente** — `biocontext` não subia, o embedding
   do NeMo estava aposentado, e `evo2_score` dava 422 em toda chamada sem
   parâmetro opcional. Nenhuma falhava de forma visível.

O padrão é **falha silenciosa**. Um ambiente autônomo não pode ter isso, porque
não há ninguém olhando quando quebra.

---

## As seis peças que faltam

Em ordem de quanto destravam, não de dificuldade.

### 1. Versionamento — a fundação que não existe

O `.gitignore` foi escrito com cuidado. Documenta que `data/raw/**` fica fora,
que `runs/nim/**` fica fora, e que `runs/provenance.jsonl` **é versionado de
propósito, porque é o caderno de laboratório e precisa de histórico**.

Só que `git init` nunca rodou. O caderno de laboratório não tem histórico.

```bash
git init && git add -A && git commit -m "estado inicial da bio-stack"
```

Sem isso, `/metodos` reconstrói o que foi rodado mas não o que o código *era*
quando rodou — e é justamente essa a pergunta da banca.

### 2. Cache de resultado — o que paga a conta

Predição de estrutura é a operação cara. Hoje, dobrar a mesma sequência duas
vezes custa duas vezes. Um cache endereçado por conteúdo (hash de sequência +
parâmetros, apontando para o arquivo já gravado) elimina isso.

Onde entra: `servers/nvidia_bio_mcp/server.py`, antes do `_post`. O `_slug(seq)`
já existe e já é usado no nome do arquivo — metade do trabalho está feita.

Impacto direto: uma varredura de mutantes com prefixo compartilhado reaproveita
MSA em vez de recomprá-lo.

### 3. Livro-caixa de crédito

Nenhum arquivo do projeto contém a palavra `budget` ou `ledger`. Não há como
responder "quanto já gastei" nem "esta varredura cabe no que sobrou".

O `provenance.py` já intercepta toda chamada. Acrescentar um campo de custo
estimado por ferramenta (a ordem de grandeza está em `docs/NVIDIA.md`) torna a
pergunta respondível sem infraestrutura nova.

Isso é pré-requisito para autonomia: um agente que roda sozinho e não sabe o que
gasta é um agente que você não deixa rodar sozinho.

### 4. Fila de trabalhos longos

O servidor tem `NVIDIA_NIM_TIMEOUT=900`. Quinze minutos é o teto de uma predição
— e amarra a sessão inteira enquanto espera. Uma varredura de 50 sequências não
cabe numa sessão interativa.

Duas rotas:

- **Barata:** `bio-mcp-queue` (organização [bio-mcp](https://github.com/bio-mcp)),
  feita exatamente para isso. Exige WSL2 — ver item 6.
- **Própria:** `runs/queue.jsonl` + um comando `/fila` que despacha e coleta.
  Mais trabalho, mas fica no mesmo padrão de proveniência do resto.

*Opinião:* comece pela própria. A fila é ~150 linhas e você controla o formato do
registro, que é o que alimenta `/metodos`.

### 5. Ligar os benchmarks

O `docs/STACK.md` lista seis benchmarks e nenhum está conectado. Sem isso não há
resposta para "a stack está ajudando?" — só impressão.

O mais adequado ao seu caso é o **BixBench** (Future-House): tarefas reais de
biologia computacional com artefato verificável. Rodar 10 tarefas dá uma linha de
base contra a qual medir cada mudança da stack.

Isso também é o que transforma o trabalho em portfólio defensável: "montei um
ambiente" é fraco; "montei um ambiente e ele sai de 40% para 65% no BixBench" é
um resultado.

### 6. WSL2 com Ubuntu — mais barato do que a documentação dizia

O `CLAUDE.md` dizia que não havia WSL. **Está instalado**: `wsl -l -v` mostra
`docker-desktop` na versão 2, e o Docker Desktop 29.7.2 está rodando sobre ele.
Falta só a distro:

```powershell
wsl --install -d Ubuntu-24.04
```

Sem reinício, sem habilitar recurso do Windows. Destrava de uma vez: `nextflow`,
`samtools`/`bcftools`, `seqkit`, `pysam`, `rapids-singlecell`, e os servidores da
organização `bio-mcp`.

---

## Camada de aquisição de dados — o buraco silencioso

Para um loop de pesquisa fechar, o agente precisa **buscar dados**, não só
analisar os que você já baixou. Testei estas fontes sem chave nenhuma em
2026-09-12, e todas responderam `200`:

| Fonte | O que traz |
|---|---|
| GEO (via E-utilities) | matrizes de expressão publicadas |
| ENA / SRA | leituras brutas |
| BioStudies / ArrayExpress | experimentos funcionais |
| PRIDE | proteômica |
| CELLxGENE | single-cell curado |
| Zenodo | dados suplementares de artigo |
| Addgene | plasmídeos caracterizados |

Nenhuma exige credencial. Nenhuma está exposta como ferramenta nesta stack — hoje
o caminho é o modelo improvisar um `curl`, que é exatamente o padrão que o
`nvidia-bio` foi criado para evitar.

*Opinião:* um servidor MCP pequeno `bio-data`, com quatro ferramentas (buscar
série no GEO, baixar matriz, resolver acesso no ENA, buscar no Zenodo), rende
mais que qualquer skill nova. É a peça que falta para o loop fechar sem você.

---

## O que a pesquisa confirmou sobre a camada NVIDIA

Verificação adversarial de 3 votos por alegação; sobreviveram 15 de 25.

- **O catálogo BioNeMo tem exatamente 12 microserviços documentados.** Esta stack
  agora registra 13 rotas (as 12 mais o `/generate` do Evo 2), com o `esmfold`
  mantido só para explicar o 410.
- **Slugs de documentação não são caminhos de API.** `evo2` na doc é
  `arc/evo2-40b` na API. Qualquer tabela de slug não testada por chamada real é
  não confiável — inclusive as deste projeto, que por isso trazem data de
  verificação.
- **O licenciamento tem três camadas**, não duas: API hospedada gratuita com
  limite de taxa; containers auto-hospedados grátis para pesquisa até 2 nós ou
  16 GPUs sob o NVIDIA Developer Program; e produção, que exige AI Enterprise
  (~US$ 4.500/GPU/ano). O `README.md` afirma que "sua licença AI Enterprise já
  cobre o compute" — **isso provavelmente está errado**: o comportamento medido
  (Evo2-40B respondendo 200, sem nenhum 403) é o da camada gratuita. Confira em
  <https://build.nvidia.com> antes de planejar volume.
- **Auto-hospedar está fora de alcance** por três barreiras independentes, cada
  uma fatal: piso de 48 GB de VRAM nas matrizes de suporte oficiais, ausência de
  model-parallelism (uma predição nunca atravessa duas GPUs), e 700 GB–1,66 TB de
  bancos de dados. A barreira não é a sua GPU ser pequena — é que não há
  configuração de GPU consumer suportada.
- **Existe um conector first-party NVIDIA BioNeMo ↔ Claude** (confiança média,
  voto 2-1): `Settings > Customize > Compute > NVIDIA BioNeMo NIM`, com importação
  de skills do `bionemo-agent-toolkit`. Vale investigar — pode substituir parte do
  servidor próprio.

---

## Ordem sugerida

| # | Ação | Custo | Destrava |
|---|---|---|---|
| 1 | `git init` + primeiro commit | minutos | reprodutibilidade real |
| 2 | Rotacionar `NVIDIA_API_KEY`, repor `NCBI_API_KEY` | minutos | ver `docs/CHAVES.md` |
| 3 | `wsl --install -d Ubuntu-24.04` | ~15 min | nextflow, samtools, bio-mcp |
| 4 | Indexar sua bibliografia (`library_index`) | ~1 h | a camada de RAG que nunca rodou |
| 5 | Cache de resultado no `nvidia-bio` | ~150 linhas | para de pagar duas vezes |
| 6 | Livro-caixa no `provenance.py` | ~60 linhas | orçamento previsível |
| 7 | Fila de trabalhos | ~150 linhas | varredura sem babá |
| 8 | Servidor `bio-data` | ~300 linhas | o loop fecha sem você |
| 9 | BixBench como linha de base | ~2 h | mede se algo disso adiantou |

Os itens 1–4 são de ambiente e valem hoje. Do 5 em diante é código, e cada um é
uma peça independente — dá para parar em qualquer ponto sem deixar nada pela
metade.
