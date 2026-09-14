---
description: Mostra o estado dos endpoints NIM de biologia da NVIDIA — quais respondem, quais mudaram de slug, quais sua licença cobre.
argument-hint: "[probe|list|config] [modelo]"
allowed-tools: Bash, Read, mcp__nvidia-bio__nim_config, mcp__nvidia-bio__nim_list_models, mcp__nvidia-bio__nim_probe_endpoints
---

Pedido: $ARGUMENTS

## Sem argumento — visão geral

1. `nim_config` — modo (hosted/local), endpoint base, se a chave está presente.
2. `nim_list_models` — catálogo local com a confiança de cada slug.
3. `nim_probe_endpoints` — teste ao vivo de quais rotas respondem.

Apresente como uma tabela: **modelo | tarefa | estado | observação**.

## Leitura dos códigos

- `4xx erro de validação` → **rota existe e funciona**. É o resultado esperado
  de um POST vazio. Verde.
- `404` → slug mudou. A NVIDIA renomeia quando promove versão. Procure o nome
  atual em https://build.nvidia.com/explore/biology e use `nim_invoke` com o
  caminho correto até o registro ser atualizado.
- `401` → chave inválida ou expirada. Rode `bootstrap/01_setup_keys.ps1`.
- `403` → chave válida mas sem acesso àquele modelo. Verifique o que sua
  licença NVIDIA AI Enterprise cobre; alguns NIMs exigem entitlement separado.
- `429` → limite de taxa. Espere.

## Se pedirem `local`

Modo local exige containers do `nvcr.io` rodando. Verifique:

```bash
docker ps --filter "ancestor=nvcr.io" --format "table {{.Names}}\t{{.Ports}}\t{{.Status}}"
```

E lembre o usuário do teto de VRAM: com 8 GB, nenhum NIM de estrutura sobe
localmente. O modo `hosted` não é um plano B aqui — é o plano certo.

## Custo

Endpoints hospedados consomem crédito da conta NVIDIA. Antes de sugerir uma
varredura de muitas sequências, estime quantas chamadas seriam e avise.
Predição de estrutura é a operação mais cara da stack.
