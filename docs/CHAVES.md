# Chaves de API — inventário medido

Tudo nesta página foi **testado ao vivo em 2026-09-12**, desta máquina, com as
credenciais desta conta. Código HTTP real, não documentação.

Reverifique com `/bio-doctor` — endpoints mudam.

---

## Estado atual

| Chave | Estado | Impacto |
|---|---|---|
| `NVIDIA_API_KEY` | **OK** | NIMs de biologia, embeddings, 3 LLMs, login no `nvcr.io` |
| `NCBI_API_KEY` | **removida em 2026-09-12** | estava com o valor da chave NVIDIA — o NCBI devolvia `400 API key invalid` |
| `S2_API_KEY` | ausente | Semantic Scholar dá `429` na 2ª requisição seguida sem chave |
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

### 1. `NCBI_API_KEY` — reposição (prioridade alta)

<https://account.ncbi.nlm.nih.gov/settings/> → *API Key Management*.
Sobe o E-utilities de 3 para 10 requisições por segundo. Formato: 36 caracteres
hexadecimais. **Não** começa com `nvapi-`.

### 2. `S2_API_KEY` — Semantic Scholar (prioridade alta para `lit-scout`)

<https://www.semanticscholar.org/product/api> → formulário, aprovação por
e-mail. Medido sem chave: a 1ª requisição passa, a 2ª e a 3ª imediatas dão
`429`. O pool anônimo é compartilhado entre o mundo inteiro — para revisão de
literatura em lote, é inutilizável.

### 3. `CROSSREF_MAILTO` / `OPENALEX_MAILTO` — não é chave

Só um e-mail no `User-Agent` ou no parâmetro `mailto=`. Move você para a fila
rápida das duas APIs. Custo zero.

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
