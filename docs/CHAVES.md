# Chaves de API — inventário medido

Tudo nesta página foi **testado ao vivo em 2026-09-12**, desta máquina, com as
credenciais desta conta. Código HTTP real, não documentação.

As linhas da `NVIDIA_API_KEY` e da `NCBI_API_KEY` foram **remedidas em
2026-10-22**, e as duas mudaram de estado desde a medição original — ver
*Estado atual* e *O incidente do `NCBI_API_KEY`*.

Reverifique com `/bio-doctor` — endpoints mudam.

---

## Estado atual

| Chave | Estado | Impacto |
|---|---|---|
| `NVIDIA_API_KEY` | **OK** · **rotacionada em 2026-10-22** | NIMs de biologia, embeddings, 3 LLMs, login no `nvcr.io` |
| `NCBI_API_KEY` | **OK** — remedida em 2026-10-22 | 36 hex, **não** começa com `nvapi-`; E-utilities devolve `200`. O incidente abaixo está fechado |
| `S2_API_KEY` | ausente — **e nenhum código lê** | ver abaixo: a variável é passada ao `biomcp`, que não a consome |
| `OPENFDA_API_KEY` | ausente | openFDA responde `200` sem chave; só sobe o teto |
| `HF_TOKEN` | ausente | download anônimo de modelo não-gated funciona (`307`) |
| `NGC_API_KEY` | ausente | **não é necessária** — `docker login nvcr.io` com a `NVIDIA_API_KEY` deu *Login Succeeded* |

### O incidente do `NCBI_API_KEY`

A variável guardava a string `nvapi-...` — a chave da NVIDIA, não uma chave do
NCBI. Duas consequências:

1. **Funcional.** O E-utilities rejeitava toda requisição com `HTTP 400`
   (`{"error":"API key invalid"}`). Sem chave nenhuma a mesma chamada devolve
   `200`. Ou seja: a chave errada era *pior* que nenhuma chave.
2. **Segurança.** O `api_key` do E-utilities vai na **query string**. A chave da
   NVIDIA foi transmitida em claro para servidores do NIH, que a ecoam de volta
   na resposta de erro e a registram em log de acesso.

O valor inválido foi removido. **Rotacione a chave NVIDIA** em
<https://build.nvidia.com> — considere-a exposta.

### Fechamento, 2026-10-22

Os dois lados do incidente estão resolvidos, e os dois foram medidos:

- **`NCBI_API_KEY` reposta.** 36 caracteres hexadecimais, não começa com
  `nvapi-`, e `einfo.fcgi` com ela na query string devolve `200`. A seção
  *"reposição (prioridade alta)"* abaixo já não é pendência.
- **`NVIDIA_API_KEY` rotacionada.** Chave nova gravada no ambiente de usuário,
  70 caracteres. Autenticação confirmada por sondagem: `boltz2` devolve `422`
  (erro de validação de corpo vazio), não `401` nem `403`.

**A rotação não restaura modelo aposentado, e isso foi medido.** `esmfold`,
`alphafold2` e `alphafold2_multimer` devolvem `410` na chave nova exatamente
como na antiga. São aposentadorias upstream da NVIDIA, sem slug substituto —
não limite de licença da conta. Quem rotacionar a chave esperando recuperá-los
vai gastar tempo à toa. A superfície real e o que usar no lugar estão em
[`MODELOS-BIO.md`](MODELOS-BIO.md).

Duas consequências operacionais da rotação:

1. **Reinicie o Claude Code.** Os servidores MCP leem o ambiente na partida.
   Até reiniciar, o `nim_config` segue reportando o prefixo da chave **antiga**
   e as ferramentas `mcp__nvidia-bio__*` continuam usando ela.
2. **O hook protege o log, o transcript não.** O `provenance.py` tem um
   `redact()` com a regex `nvapi-[A-Za-z0-9_\-]{10,}` aplicada ao `input` antes
   de gravar, então `runs/provenance.jsonl` guarda `[REDACTED]`. Chave colada
   em conversa, porém, fica no transcript da sessão em texto claro. Prefira
   `bootstrap/01_setup_keys.ps1`, que lê por `SecureString`.

---

## O que é aberto (não precisa de chave)

Testado com `GET` simples, sem cabeçalho de autorização — todos `HTTP 200`:

| Domínio | Bases |
|---|---|
| Proteína / estrutura | UniProt · RCSB PDB · AlphaFold DB · InterPro · STRING |
| Química | ChEMBL · PubChem · Rhea |
| Genoma / variante | Ensembl · UCSC · gnomAD · MyVariant · Open Targets (POST GraphQL) |
| Via / anotação | Reactome · KEGG · Bgee |
| Literatura | Europe PMC · OpenAlex · Crossref |
| Dados | GEO (E-utilities) · ENA · BioStudies/ArrayExpress · PRIDE · Zenodo · CELLxGENE · Addgene |

Nenhuma dessas precisa de credencial para uso acadêmico em volume normal.
Mande um `User-Agent` com e-mail de contato — Crossref e OpenAlex dão limite
melhor no "polite pool", e é boa educação.

---

## Chaves que valem a pena pegar

Todas gratuitas, todas em minutos.

### 1. `NCBI_API_KEY` — ~~reposição (prioridade alta)~~ **já reposta**

Resolvida: medida em 2026-10-22, válida, `200` no E-utilities. Fica aqui como
referência de formato, não como tarefa.

<https://account.ncbi.nlm.nih.gov/settings/> → *API Key Management*.
Sobe o E-utilities de 3 para 10 requisições por segundo. Formato: 36 caracteres
hexadecimais. **Não** começa com `nvapi-`.

### 2. `S2_API_KEY` — Semantic Scholar (**não é necessária**)

Rastreado em 2026-10-10: a variável aparece em quatro lugares e **nenhum deles
a consome**.

| Lugar | O que faz com ela |
|---|---|
| `.mcp.json:65` | repassa como env para o servidor `biomcp` |
| `.claude/commands/bio-doctor.md:26` | só relata presente/ausente |
| `.claude/hooks/session_context.py:166` | só relata presente/ausente |
| `.codex/config.toml:28` | só declara como variável a repassar |

O único consumidor possível era o `biomcp`, e o `biomcp-python 0.7.3` instalado
tem zero ocorrências de `S2_API_KEY` ou `semanticscholar`:

```bash
grep -rc "S2_API_KEY\|semanticscholar" \
  "$(find ~/AppData/Local/uv/cache/archive-v0 -maxdepth 5 -type d -name biomcp | head -1)"
# 0
```

O `lit-scout` **nunca** dependeu dela: as ferramentas dele são PubMed,
bioRxiv, Europe PMC, Open Targets, STRING, Reactome e `nvidia-nemo`. Nenhuma
fala com o Semantic Scholar. A versão anterior desta seção chamava a chave de
"prioridade alta para `lit-scout`" — estava errado, e o `/bio-doctor` repetia
o erro reportando `MISS` como se fosse lacuna.

Se um dia um servidor que *use* o S2 entrar na stack: conta não basta, o pedido
é um formulário separado em <https://www.semanticscholar.org/product/api> →
*Request an API Key*, a chave chega por e-mail, o teto inicial é 1 req/s e o
Allen Institute não publica prazo de aprovação. Medido sem chave em 2026-09-12:
a 1ª requisição passa, a 2ª e a 3ª imediatas dão `429`.

### 3. `CROSSREF_MAILTO` / `OPENALEX_MAILTO` — **não existem na stack**

A recomendação anterior aqui era colocar um e-mail no `User-Agent` para entrar
na fila rápida do Crossref e do OpenAlex, a custo zero. Verificado em
2026-10-10: essas duas variáveis aparecem **apenas nesta página** — nenhum
servidor MCP, hook ou script as lê. Definir no ambiente não muda nada.

Para valer, o e-mail teria de ser injetado por quem faz a requisição HTTP, e
hoje ninguém na stack chama Crossref ou OpenAlex direto.

### 4. `HF_TOKEN` — Hugging Face (opcional)

<https://huggingface.co/settings/tokens>. Modelo não-gated baixa sem token.
Precisa dele para: modelos com licença aceita manualmente, repositórios
privados, e limite de download maior. O conector MCP da conta já está
autenticado (`patooriginal`) — isso cobre *busca*, não *download local*.

### 5. `OPENFDA_API_KEY` — baixa prioridade

<https://open.fda.gov/apis/authentication/>. Sem chave: 240 req/min. Com chave:
240 req/min e 120 000/dia. Só importa em varredura grande.

---

## Chaves que **não** valem a pena agora

| Serviço | Por que não |
|---|---|
| `NGC_API_KEY` | redundante — a `NVIDIA_API_KEY` já autentica no `nvcr.io` |
| Benchling / IDT / Twist | plataformas de bancada e encomenda; o trabalho aqui é in silico |
| Seqera / Tower | só depois que houver Nextflow rodando de verdade |
| Elsevier / Springer | acesso a texto completo depende de assinatura institucional; peça pela biblioteca do Senac antes de gastar tempo com chave |
| GPU serverless (Modal, Replicate, RunPod) | só se um modelo específico não tiver NIM hospedado; hoje o catálogo da NVIDIA cobre o que esta stack usa |

---

## Como gravar

```powershell
.\bootstrap\01_setup_keys.ps1
```

Lê por `SecureString`, grava no ambiente de usuário do Windows, nunca em
arquivo versionado. Reinicie o Claude Code depois — os servidores MCP herdam o
ambiente na partida.
