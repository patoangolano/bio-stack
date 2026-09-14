---
description: Audita tudo que a sua chave NVIDIA destrava — NIMs de biologia, NeMo Retriever, catálogo de LLMs e registry de containers — e diz o que está fora do entitlement.
argument-hint: "[bio|retriever|llm|registry|tudo]"
allowed-tools: Bash, Read, mcp__nvidia-bio__nim_config, mcp__nvidia-bio__nim_probe_endpoints, mcp__nvidia-bio__nim_list_models, mcp__nvidia-nemo__nemo_config
---

Escopo pedido: $ARGUMENTS (padrão: tudo)

Uma chave NVIDIA com AI Enterprise cobre quatro superfícies distintas. Entitlements
variam por licença — o objetivo deste comando é descobrir empiricamente o que a
chave **desta** conta alcança, não repetir o que o marketing promete.

## 1. Credencial

```bash
powershell -NoProfile -Command "if ([Environment]::GetEnvironmentVariable('NVIDIA_API_KEY','User')) { 'configurada' } else { 'AUSENTE' }"
```

Se ausente, pare aqui e instrua: `bootstrap\01_setup_keys.ps1`. Nunca peça a
chave no chat — o script lê via SecureString.

## 2. Biologia — `health.api.nvidia.com`

`nim_config` e depois `nim_probe_endpoints`. Tabela: **modelo | tarefa | estado**.

Leitura dos códigos: `4xx de validação` = rota existe (verde); `404` = slug mudou;
`403` = fora do entitlement; `401` = credencial ruim.

## 3. NeMo Retriever — `integrate.api.nvidia.com`

`nemo_config` para ver modelos e índices. Depois teste os dois endpoints com uma
chamada mínima:

```bash
curl -s -o /dev/null -w "embeddings: %{http_code}\n" -X POST \
  https://integrate.api.nvidia.com/v1/embeddings \
  -H "Authorization: Bearer $NVIDIA_API_KEY" -H "Content-Type: application/json" \
  -d '{"input":["teste"],"model":"nvidia/llama-3.2-nv-embedqa-1b-v2","input_type":"query"}'

curl -s -o /dev/null -w "reranking:  %{http_code}\n" -X POST \
  https://integrate.api.nvidia.com/v1/retrieval/nvidia/reranking \
  -H "Authorization: Bearer $NVIDIA_API_KEY" -H "Content-Type: application/json" \
  -d '{"model":"nvidia/llama-3.2-nv-rerankqa-1b-v2","query":{"text":"a"},"passages":[{"text":"b"}]}'
```

200 = disponível. Estas duas chamadas custam frações de centavo.

## 4. Catálogo de modelos

```bash
curl -s https://integrate.api.nvidia.com/v1/models \
  -H "Authorization: Bearer $NVIDIA_API_KEY" | python -c "
import json,sys
d=json.load(sys.stdin).get('data',[])
print(f'{len(d)} modelos visíveis')
bio=[m['id'] for m in d if any(k in m['id'].lower() for k in
     ('bio','prot','mol','gene','dna','rna','fold','med','embedqa','rerank'))]
print('relevantes para esta stack:')
[print(' ', m) for m in sorted(bio)]
"
```

**Este catálogo NÃO é o inventário de entitlement.** Medido em 2026-09-12:
`GET /v1/models` lista 82 modelos, mas dos 16 testados com uma chamada real só
três responderam 200 — `openai/gpt-oss-20b`, `nvidia/nemotron-3-super-120b-a12b`
e `nvidia/nemotron-3.5-lightning-30b-a3b`. Os outros devolvem
`404 Function '<uuid>': Not found for account '<conta>'`, inclusive modelos que
aparecem na listagem (`nvidia/llama-3.2-nv-embedqa-1b-v1`, `nvidia/embed-qa-4`,
`writer/palmyra-med-70b`).

A listagem é o catálogo do serviço, não o que a conta alcança. A inferência
válida é de mão única: **não aparecer** ali é sinal ruim, mas **aparecer não
garante nada**. A única prova é a chamada.

Então: depois de listar, teste com `max_tokens: 5` os modelos que interessam.

## 5. Registry de containers — `nvcr.io`

Só importa se for auto-hospedar NIM (precisa de ≥24 GB de VRAM — esta máquina
tem 8 GB, então é para o cluster da universidade).

```bash
docker --version 2>/dev/null || echo "docker ausente — só faz sentido com WSL2"
```

Se houver Docker: `echo $NVIDIA_API_KEY | docker login nvcr.io -u '$oauthtoken' --password-stdin`

## Relatório

Tabela **superfície | estado | o que dá para fazer com isso**, e ao final:

- O que está fora do entitlement e qual o caminho alternativo.
- Quais slugs do registro local divergem do catálogo real (proponha o patch em
  `servers/nvidia_bio_mcp/server.py`, não aplique sem confirmar).
- Lembrete de custo: predição de estrutura é a operação cara; embedding e
  reranking são baratos.

Não faça chamadas caras (folding, docking) só para testar — os probes acima bastam.
