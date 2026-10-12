# Camada NVIDIA — referência

## As três rotas de compute

| Rota | Quando | Nesta máquina |
|---|---|---|
| **NIM hospedado** `health.api.nvidia.com` | padrão para tudo que é modelo grande | **use esta** |
| **NIM auto-hospedado** container `nvcr.io` | privacidade de dados, volume alto, offline | inviável (8 GB) |
| **GPU local** PyTorch/RAPIDS | modelos pequenos, pré-processamento | modelos ≤ 3 GB |

Com NVIDIA AI Enterprise, a mesma chave NGC serve para as três.

---

## Endpoints hospedados

Padrão de URL:

```
https://health.api.nvidia.com/v1/biology/<org>/<modelo>/<ação>
```

Autenticação:

```
Authorization: Bearer $NVIDIA_API_KEY
Content-Type: application/json
```

### Protocolo assíncrono

Predições longas retornam **202** com o header `nvcf-reqid`. O resultado sai
por polling:

```
GET https://health.api.nvidia.com/v1/status/<reqid>
```

até virar 200. O servidor `nvidia-bio` deste projeto já faz isso —
`servers/nvidia_bio_mcp/server.py`, função `_poll_nvcf`. Se você chamar por
`curl` direto, precisa implementar o loop.

### Modelos

| Chave | Caminho | Tarefa |
|---|---|---|
| `boltz2` | `biology/mit/boltz2/predict` | cofolding proteína-ligante + afinidade |
| `openfold2` | `biology/openfold/openfold2/predict-structure-from-msa-and-template` | dobramento com MSA |
| `openfold3` | `biology/openfold/openfold3/predict` | complexos biomoleculares |
| ~~`esmfold`~~ | ~~`biology/meta/esmfold`~~ | **aposentado — HTTP 410** |
| `alphafold2` | `biology/deepmind/alphafold2` | AF2 monomérico |
| `alphafold2_multimer` | `biology/deepmind/alphafold2-multimer` | complexo proteico multicadeia |
| `msa_search` | `biology/colabfold/msa-search/predict` | busca de MSA |
| `diffdock` | `biology/mit/diffdock` | docking cego |
| `rfdiffusion` | `biology/ipd/rfdiffusion/generate` | geração de backbone |
| `proteinmpnn` | `biology/ipd/proteinmpnn/predict` | desenho de sequência p/ backbone |
| `molmim` | `biology/nvidia/molmim/generate` | moléculas otimizadas por propriedade |
| `genmol` | `biology/nvidia/genmol/generate` | geração fragmento→molécula |
| `evo2` | `biology/arc/evo2-40b/forward` | embeddings por camada; `unembed` dá logits |
| `evo2_generate` | `biology/arc/evo2-40b/generate` | geração autoregressiva de DNA |

> **Os slugs mudam.** A NVIDIA renomeia quando promove versões, e esta tabela
> já esteve errada. Os caminhos acima foram verificados endpoint por endpoint
> em **2026-09-12**; `boltz2`, `rfdiffusion`, `proteinmpnn`, `molmim` e `genmol`
> têm resposta 200 de uma chamada real, os demais respondem 422 (rota existe,
> payload de teste inválido). Rode `/nim` (que chama `nim_probe_endpoints`)
> para reverificar antes de assumir que um 404 é erro seu.
>
> Três armadilhas confirmadas na medição:
>
> - `alphafold2` e `diffdock` **não têm segmento de ação**. Acrescentar
>   `/predict` ou `/generate` devolve 404.
> - `esmfold` foi desligado. `biology/meta/esmfold` devolve **410 Gone** e não
>   há slug substituto — use `boltz2`, `openfold2` ou `alphafold2`.
> - `openfold3` está sob a org `openfold`, não `nvidia`.
> - `alphafold2-multimer` existe e **não estava registrado** até 2026-09-12.
>   O campo obrigatório é `sequences` (lista), não `sequence` — é a única
>   rota hospedada para complexo proteico sem ligante.
>
> Catálogo autoritativo: <https://build.nvidia.com/explore/biology>

### Exemplo cru

```bash
curl -X POST https://health.api.nvidia.com/v1/biology/mit/boltz2/predict \
  -H "Authorization: Bearer $NVIDIA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "polymers": [{"id":"A","molecule_type":"protein","sequence":"MALWMRLLPLLALLALWGPDPAAA"}],
    "recycling_steps": 3,
    "sampling_steps": 50,
    "diffusion_samples": 1,
    "step_scale": 1.638,
    "output_format": "mmcif"
  }'
```

### Evo 2 — como pontuar variante de verdade

O endpoint tem duas rotas e é fácil usar a errada. Medido em 2026-09-12:

| Rota | Devolve | Serve para |
|---|---|---|
| `/forward` | só as camadas pedidas em `output_layers` | embeddings |
| `/forward` com `output_layers: ["unembed"]` | **logits por posição da entrada**, `(1, n, 512)` | **verossimilhança** |
| `/generate` com `enable_sampled_probs` | probabilidade dos tokens **gerados** | continuação sintética |

A armadilha: `sampled_probs` tem o comprimento de `num_tokens`, não da
sequência de entrada. Pontuar uma variante com ele mede a continuação que o
modelo inventou, não a sequência que você quer avaliar.

O caminho certo é `unembed`. O Evo 2 é byte-level: o vocabulário de 512
posições é indexado por `ord()` do caractere, então o logit de `A` fica no
índice 65. Log-softmax nos logits, gather na base realmente observada na
posição seguinte, e sai `log P(base_i+1 | prefixo)` — a sequência inteira numa
chamada só. Confirmação de que a indexação está certa: a massa de
probabilidade em `{A,C,G,T}` soma exatamente 1,0.

Validação de sanidade que vale repetir ao mexer nisso: um trecho real do `lacI`
de *E. coli* dá perplexidade **1,065**; uma sequência inventada dá **3,9**
(perto de 4,0, que é o acaso sobre quatro bases). Se a sua sequência real
pontuar perto de 4, alguma coisa está errada no pipeline, não no modelo.

Duas armadilhas do payload:

- `output_layers` é **obrigatório e não aceita lista vazia** — devolve
  `422 List should have at least 1 item`. Não existe "só os logits" por omissão.
- Nomes de camada são atributos do StripedHyena: `unembed`, `norm`,
  `embedding_layer` e `blocks.N.mlp.l3` funcionam; `logits`, `lm_head` e
  `output` devolvem 422.

A ferramenta `evo2_likelihood` do servidor MCP já implementa tudo isso e grava
o vetor por base em disco, devolvendo só o resumo.

---

## Desenho de proteína e química generativa

Quatro NIMs que estavam no catálogo mas só eram alcançáveis por `nim_invoke`
cru. Desde 2026-09-12 têm ferramenta tipada, com validação de parâmetro e
gravação em disco como o resto do servidor.

| Ferramenta MCP | NIM | O que faz |
|---|---|---|
| `design_backbone` | rfdiffusion | gera esqueleto de proteína condicionado a um alvo |
| `design_sequence` | proteinmpnn | veste um esqueleto com sequência de aminoácidos |
| `generate_molecules` | molmim | otimiza uma molécula semente por propriedade |
| `generate_from_fragment` | genmol | cresce uma molécula a partir de um fragmento |

### O fluxo de desenho é encadeado

RFdiffusion devolve **só backbone** — coordenadas sem identidade de resíduo.
Sozinho ele não responde "que proteína é essa". A sequência vem do ProteinMPNN,
e a validação vem de dobrar o resultado de volta:

```
design_backbone  ->  design_sequence  ->  fold_protein  ->  compara com o alvo
   (esqueleto)        (aminoácidos)        (estrutura)      (o desenho fecha?)
```

Pular o último passo é o erro comum: o ProteinMPNN sempre devolve *alguma*
sequência, inclusive para um backbone ruim.

### Formatos de resposta que fogem do padrão

Medidos com chamada real, não lidos da documentação:

- **RFdiffusion** devolve a estrutura na chave `output_pdb` (string), não em
  `structures`. O handler do servidor foi corrigido para reconhecê-la — antes
  caía no fallback e gravava JSON cru.
- **ProteinMPNN** devolve `mfasta` (multi-FASTA), `scores` e `probs`. O primeiro
  registro do `mfasta` é a **sequência nativa** do PDB de entrada, não um
  desenho; os desenhos vêm depois, com `seq_recovery` no cabeçalho. A matriz
  `probs` é `[n_seq][n_resíduo][21]` — 42 kB só para a crambina, com 46
  resíduos e 2 sequências. Vai para disco sempre.
- **MolMIM** devolve `molecules` como **JSON dentro de uma string**, e chama o
  SMILES de `sample`. **GenMol** devolve lista de verdade, com a chave `smiles`.
  O servidor normaliza os dois para `{"smiles", "score"}`.
- **GenMol** exige a máscara `[*{N-N}]` no fragmento indicando onde crescer.
  Sem ela devolve variações do próprio fragmento. Com `unique: true` o retorno
  vem com menos moléculas do que você pediu — é deduplicação, não erro.

### Como ler o `seq_recovery`

Fração de resíduos do desenho idênticos à sequência nativa. Numa redesenhagem
de proteína conhecida, 0,4–0,5 é típico e saudável. Perto de 1,0 significa que
o modelo copiou a entrada em vez de desenhar; perto de 0 costuma indicar
backbone ruim ou temperatura alta demais.

`sampling_temp` controla isso: 0,1 é conservador, 0,3 diversifica, acima disso
degrada.

### Escopo

Estas ferramentas são para desenho estrutural em alvos modelo — crambina,
lisozima, GFP, peças caracterizadas do Registry. Desenho de proteína com função
nova, ou dirigido a alvo de risco elevado, sai do escopo deste projeto e vai ao
comitê de biossegurança da instituição, como o resto do trabalho de bancada.

---

## NeMo Retriever — busca semântica na sua biblioteca

A mesma chave cobre o NeMo Retriever, que faz parte do NVIDIA AI Enterprise.
Host diferente dos NIMs de biologia: `integrate.api.nvidia.com`, não
`health.api.nvidia.com`.

| Operação | Endpoint | Modelo | Estado em 2026-09-12 |
|---|---|---|---|
| Embedding | `POST /v1/embeddings` | `nvidia/nemotron-3-embed-1b` | **200**, dim 2048 |
| Reranking | `POST /v1/retrieval/nvidia/reranking` | — | **404 em todas as rotas** |
| Catálogo | `GET /v1/models` | — | 200, 82 modelos listados |

> **O modelo antigo morreu.** `nvidia/llama-3.2-nv-embedqa-1b-v2` responde
> **HTTP 410**: *"reached its end of life on 2026-05-18"*. O `...-v1` aparece no
> catálogo mas devolve 404 para esta conta. O substituto que funciona é
> `nvidia/nemotron-3-embed-1b` — 2048 dimensões, aceita `input_type`
> `passage`/`query` igual ao anterior. Já trocado no `.mcp.json`.
>
> **Reranking não tem rota viva.** Foram testadas sete: `/v1/retrieval/nvidia/reranking`,
> `/v1/retrieval/nvidia/<modelo>/reranking` (quatro variantes), `/v1/ranking`,
> `/v1/rerank`, e o host alternativo `ai.api.nvidia.com`. Todas 404. Nenhum
> modelo com "rerank" no nome aparece no catálogo desta chave.
>
> Consequência: o segundo estágio do pipeline de RAG está indisponível pela
> NVIDIA. O `library_search` degrada sozinho para ordem por cosseno e avisa —
> mas `rerank_passages` levanta erro. A alternativa que cabe nesta máquina é um
> cross-encoder local (`BAAI/bge-reranker-v2-m3`, ~560 M, multilíngue, ~1,2 GB
> em fp16), que roda folgado nos 8 GB e não custa crédito.

### O parâmetro que quase todo mundo erra

`input_type` é obrigatório e assimétrico:

```
indexar documento  ->  "input_type": "passage"
buscar             ->  "input_type": "query"
```

Usar o mesmo valor nos dois lados **não dá erro** — só degrada silenciosamente a
qualidade da busca. O servidor `nvidia-nemo` trata isso internamente.

### Exemplo cru

```bash
curl -X POST https://integrate.api.nvidia.com/v1/embeddings   -H "Authorization: Bearer $NVIDIA_API_KEY"   -H "Content-Type: application/json"   -d '{
    "input": ["Batch effect em scRNA-seq de tecido criopreservado"],
    "model": "nvidia/nemotron-3-embed-1b",
    "input_type": "passage",
    "encoding_format": "float",
    "truncate": "END"
  }'
```

### Por que o reranking importa

O pipeline canônico da NVIDIA tem dois estágios porque eles fazem coisas
diferentes:

```
pergunta -> embedding -> busca por cosseno -> top-40 candidatos
         -> reranking -> top-8 na ordem certa
```

O embedding compara vetores calculados **separadamente** para pergunta e
passagem: é rápido, escala para milhões de trechos, e é grosseiro. O reranker lê
pergunta e passagem **juntas** num único forward pass: é caro, não escala, e é
bem mais preciso.

Implementação caseira de RAG quase sempre pula o segundo estágio. É a maior
diferença de qualidade por linha de código nesta stack inteira.

### Custo

Embedding e reranking são ordens de magnitude mais baratos que predição de
estrutura. Indexar 200 PDFs custa alguns milhares de tokens de embedding —
comparável a uma única predição do Boltz-2. Não economize aqui; economize em
folding.

---

## NIM auto-hospedado

Só vale com ≥ 24 GB de VRAM. Registrado aqui para quando você tiver acesso ao
cluster da universidade ou DGX Cloud.

```bash
echo $NGC_API_KEY | docker login nvcr.io -u '$oauthtoken' --password-stdin

docker run --rm --gpus all --shm-size=4g \
  -e NGC_API_KEY \
  -v ~/.cache/nim:/opt/nim/.cache \
  -p 8000:8000 \
  nvcr.io/nim/mit/boltz2:latest

curl http://localhost:8000/v1/health/ready
```

Com o container no ar, troque o modo do servidor MCP:

```jsonc
// .mcp.json
"env": { "NVIDIA_NIM_MODE": "local", "NVIDIA_NIM_LOCAL_BASE": "http://localhost:8000" }
```

Note que o caminho local **não** tem o prefixo `/v1`:
`http://localhost:8000/biology/mit/boltz2/predict`.

---

## GPU local — o que realmente roda em 8 GB

| Carga | VRAM | Roda? |
|---|---|---|
| ESM-2 t33 650M (embeddings) | 1,25 GiB fp16, medido | sim, folgado |
| ESM-2 t36 3B | ~6,5 GB fp16 | sim, batch 1 |
| DNABERT-2 117M | ~0,5 GB | sim |
| nvMolKit (fingerprints, conformers) | baixo | sim |
| rapids-singlecell | ~4 GB p/ 100k células | sim, **mas só em WSL2** |
| ESMFold | ~14 GB | não |
| Evo 2 7B | ~20 GB | não |
| Parabricks (fq2bam) | ≥ 16 GB | não |

### PyTorch em Blackwell

RTX 50xx é `sm_120`. Precisa de wheels **cu128 ou superior**:

```powershell
uv pip install --index-url https://download.pytorch.org/whl/cu128 torch torchvision
```

Instalar do PyPI padrão traz build de CPU ou um CUDA velho que não reconhece
`sm_120` — sintoma típico é `CUDA error: no kernel image is available`.

### Embeddings de proteína local

Pelo `fair-esm`, que é o caminho sem armadilha:

```python
import torch, esm

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.half().cuda().eval()
bc = alphabet.get_batch_converter()

_, _, toks = bc([("alvo", "MALWMRLLPLLALLALWGPDPAAA")])
with torch.no_grad():
    out = model(toks.cuda(), repr_layers=[33])
# O batch converter envolve a sequência em <cls> … <eos> (tokens 0 e 2); os dois
# carregam estado próprio e entram na média como se fossem resíduos.
emb = out["representations"][33][:, 1:-1, :].mean(1)   # (1, 1280)
```

Uso acadêmico: comparar variantes por distância no espaço de embedding,
agrupar homólogos, alimentar um classificador simples de propriedade.

#### Pelo `transformers`: não use `AutoModel`

> **`AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D")` ganha um pooler
> inicializado aleatoriamente.** O `AutoModel` resolve para `EsmModel`, que
> constrói `pooler.dense` por padrão — e esse peso **não existe** no checkpoint do
> ESM-2. O relatório de carga do `transformers` 5.19 marca
> `missing_keys: {'pooler.dense.weight', 'pooler.dense.bias'}`: são inicializados
> na hora, com o RNG do processo. Medido nesta máquina carregando duas vezes no
> mesmo processo: `last_hidden_state` sai bit-a-bit idêntico, `pooler_output`
> **não** (delta máximo 0,938). Quem usar `pooler_output` como embedding colhe um
> número diferente a cada execução, sem erro, sem aviso no stderr.
>
> É a mesma classe de defeito que o `CLAUDE.md` proíbe em `np.random.seed()` —
> estado global de processo mudando o resultado sem nada falhar — numa camada
> diferente. Fixar seed não resolve: o pooler continua sendo ruído, só ruído
> reprodutível.

O contorno é pedir o modelo sem a camada e fazer o pooling à mão:

```python
import torch
from transformers import AutoTokenizer, EsmModel

NOME = "facebook/esm2_t33_650M_UR50D"
SEQUENCIAS = ["MALWMRLLPLLALLALWGPDPAAA", "MKTAYIAKQRQISFVKSHFSRQ"]

tokenizador = AutoTokenizer.from_pretrained(NOME)

# add_pooling_layer=False: sem isto o transformers cria pooler.dense com pesos
# aleatórios, porque o checkpoint do ESM-2 não tem esse tensor.
modelo = EsmModel.from_pretrained(NOME, add_pooling_layer=False, dtype=torch.float16)
modelo = modelo.cuda().eval()

entrada = tokenizador(SEQUENCIAS, padding=True, return_tensors="pt")
entrada = {chave: valor.cuda() for chave, valor in entrada.items()}

with torch.no_grad():
    estados = modelo(**entrada).last_hidden_state

# O tokenizador envolve cada sequência em <cls> … <eos> e completa com <pad>.
# Os três carregam estado próprio; a média tem de cair só sobre os resíduos.
mascara = entrada["attention_mask"].clone()
mascara[:, 0] = 0                                    # <cls>
ultimo = entrada["attention_mask"].sum(1) - 1        # <eos> de cada sequência
mascara[torch.arange(mascara.shape[0], device=mascara.device), ultimo] = 0

peso = mascara.unsqueeze(-1).to(estados.dtype)
emb = (estados * peso).sum(1) / peso.sum(1)          # (2, 1280)

print(emb.shape, emb.dtype)
```

Verificado nesta máquina com `.venv/Scripts/python.exe`, `transformers` 5.19.0 e
`torch` 2.11.0+cu128:

- `EsmModel(..., add_pooling_layer=False)` carrega com `missing_keys` vazio — só
  `lm_head.*` sobra como `UNEXPECTED`, que é o esperado ao tirar a cabeça de
  linguagem. Nenhum peso inicializado na hora.
- `mascara.sum(1)` devolve `[24, 22]`, o comprimento real das duas sequências:
  `<cls>`, `<eos>` e `<pad>` ficaram fora da média.
- O embedding sai **bit-a-bit idêntico** entre duas cargas no mesmo processo.
- Pico de VRAM com `torch.cuda.max_memory_allocated()`: **1250,3 MiB**
  (1 311 052 800 bytes) para estas duas sequências em fp16 — as ativações crescem
  com o comprimento, o peso fixo do modelo é ~1,25 GiB. Folga larga nos 8151 MiB.

---

## Custo

Endpoints hospedados consomem crédito da conta. Ordem de custo por chamada:

```
molmim ≈ genmol  <  proteinmpnn  <  msa_search  <  rfdiffusion  <  diffdock
  <  boltz2  <  openfold3  ≈  alphafold2
```

Antes de varrer muitas sequências, estime o número de chamadas. Uma varredura
de 500 mutantes em Boltz-2 é cara — considere triar primeiro com embeddings do
ESM-2 local e só dobrar os 20 melhores.
