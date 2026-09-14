---
description: Prediz a estrutura 3D de uma proteína, checando antes se já existe estrutura experimental ou no AlphaFold DB.
argument-hint: "<sequência, UniProt ID, ou nome do gene> [--model boltz2|esmfold|openfold3]"
allowed-tools: Read, Write, Bash, mcp__nvidia-bio__fold_protein, mcp__nvidia-bio__cofold_complex, mcp__nvidia-bio__nim_probe_endpoints, mcp__gget__gget_pdb, mcp__gget__gget_seq, mcp__gget__gget_search, mcp__gget__gget_info
---

Alvo: $ARGUMENTS

Delegue ao agente `structure-bio` se a tarefa envolver análise estrutural
profunda. Para uma predição direta, siga aqui.

## Passo 1 — Resolver a entrada

O que o usuário passou pode ser:

- **Sequência de aminoácidos** → use direto.
- **UniProt ID** (P04637) → `gget_seq` para obter a sequência.
- **Nome de gene** (TP53) → `gget_search` para achar o ID, depois `gget_seq`.
  Confirme o organismo com o usuário se for ambíguo.

## Passo 2 — Não prediga o que já existe

Antes de gastar minutos e crédito:

1. **Estrutura experimental?** `gget_pdb` exige um PDB ID; ele não busca por
   gene. Descubra os IDs com `gget info <Ensembl_ID>` ou pela Search API do RCSB,
   e só então puxe com `gget_pdb <PDB_ID>`. Mostre resolução e método, e pergunte
   se ainda quer a predição.
2. **AlphaFold DB?** Não use `gget_alphafold` — ele roda o AlphaFold local e não
   cabe em 8 GB. Consulte o banco direto:
   `https://alphafold.ebi.ac.uk/api/prediction/<UniProt_ID>` — o mmCIF traz o
   pLDDT no campo B-factor.

Se qualquer um dos dois resolver, entregue isso e **explique por que não
prediz**. Predizer o que já foi cristalizado é desperdício.

## Passo 3 — Predizer

`mcp__nvidia-bio__fold_protein` com o modelo escolhido (padrão `boltz2`).
Se der 404, rode `nim_probe_endpoints` antes de reportar falha.

Se a sequência tiver > 1500 resíduos, avise que vai demorar e pergunte se
prefere predizer por domínio.

## Passo 4 — Interpretar

Nunca entregue só o caminho do arquivo. Reporte:

- pLDDT médio e a distribuição (quais regiões estão baixas).
- Se há regiões de pLDDT < 50, diga que provavelmente são desordenadas —
  não é defeito da predição.
- ipTM/PAE se for complexo.
- Onde o arquivo foi gravado e como abrir (PyMOL/ChimeraX).

Gere também um snippet pronto de visualização:

```
# PyMOL - escala de cor do AlphaFold (azul = confiavel, laranja = desordenado)
load runs/nim/<arquivo>.cif
show cartoon
spectrum b, red_yellow_white_blue, minimum=50, maximum=90
```

Confira a escala antes: AF2/ESMFold gravam pLDDT em 0-100, mas algumas saidas do
Boltz usam 0-1. Se os valores estiverem entre 0 e 1, use `minimum=0.5, maximum=0.9`.

## Passo 5 — Registrar

A chamada já entra em `runs/provenance.jsonl` automaticamente. Confirme ao
usuário que o registro existe e que `/metodos` vai transformar isso em texto.
