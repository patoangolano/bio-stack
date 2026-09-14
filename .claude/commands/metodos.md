---
description: Gera a seção de Métodos do seu relatório/TCC a partir do log de proveniência das análises que você realmente rodou.
argument-hint: "[período, ex: 'última semana' ou 'hoje'] [formato: abnt|vancouver|texto]"
allowed-tools: Read, Write, Glob, Grep, Bash
---

Escreva a seção de Métodos a partir do que foi **efetivamente executado**,
não do que se pretendia executar.

Escopo pedido: $ARGUMENTS

## Fontes, nesta ordem

1. `runs/provenance.jsonl` — registro append-only de toda ferramenta chamada,
   com parâmetros e horário UTC. Esta é a fonte primária.
2. `runs/nim/` — arquivos de saída gerados (estruturas, MSAs, embeddings).
3. Scripts e notebooks no projeto — a lógica de análise.
4. Versões de pacotes: `uv pip list` ou o `pipeline_info/` do nf-core.

## Como ler o log

```bash
python -c "
import json, collections
rows = [json.loads(l) for l in open('runs/provenance.jsonl', encoding='utf-8')]
print('total:', len(rows))
print(collections.Counter(r['kind'] for r in rows))
print('periodo:', rows[0]['ts'], '->', rows[-1]['ts'])
"
```

Depois filtre pelo período pedido e agrupe por tipo de operação.

## O que a seção precisa conter

Um leitor precisa conseguir **repetir** o trabalho. Então, para cada etapa:

- Software e **versão exata** (não "usamos o DESeq2", mas "pydeseq2 v0.4.12").
- Parâmetros que não são padrão — e diga explicitamente quando usou o padrão.
- Base de dados e **data de acesso** (bases mudam; UniProt de março não é
  UniProt de setembro).
- Para modelos preditivos: nome, versão do endpoint, e as métricas de confiança
  obtidas.
- Hardware, quando for relevante para tempo de execução.

## Honestidade

- Inclua as tentativas que falharam se elas mudaram o rumo da análise. "Primeiro
  tentamos X, que produziu Y, e por isso adotamos Z" é bom método, não fraqueza.
- Se o log mostra que um parâmetro foi mudado no meio do caminho, reporte o
  valor final e mencione a mudança.
- **Não invente versão de pacote.** Se a versão não está no log nem no ambiente,
  escreva `[VERSÃO A CONFIRMAR]` e liste esses pontos no final para o usuário
  preencher.
- Se uma análise usou dados gerados por IA (Synthesize Bio / GEM-1), isso
  precisa estar explícito no texto, não numa nota de rodapé.

## Saída

Grave em `docs/metodos.md`. Prosa corrida em português, tempo passado, voz
passiva ou primeira pessoa do plural (o padrão acadêmico brasileiro), com as
citações no formato pedido. Ao final do arquivo, uma lista de **pendências**:
tudo que precisa de confirmação humana.

Mostre o texto no chat também, para revisão imediata.
