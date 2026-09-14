---
description: Testa se um caráter foi perdido por perda de gene ou por mudança regulatória, e desenha a reversão in silico se houver alvo.
argument-hint: "<caráter> em <clado foco> versus <clado irmão>  (ex: 'dentes em Aves versus Crocodylia')"
allowed-tools: Read, Write, Edit, Bash, Glob, Grep, Task, WebSearch, WebFetch, mcp__gget__gget_search, mcp__gget__gget_info, mcp__gget__gget_seq, mcp__gget__gget_blast, mcp__gget__gget_blat, mcp__biocontext, mcp__nvidia-bio__evo2_score, mcp__nvidia-bio__evo2_likelihood
---

Caráter e clados: $ARGUMENTS

Delegue ao agente `evo-devo`. Ele tem as regras de evidência; este comando é a
ordem das etapas. Não pule a 2 nem a 5.

## Passo 1 — Montar o eixo, não o par

Duas espécies dão uma diferença; três ou mais com topologia conhecida dão uma
**mudança**, com direção e ramo. Antes de consultar qualquer base, fixe:

- **foco** — onde se suspeita da perda
- **irmão** — o clado que reteve o caráter
- **externo** — enraíza os dois
- **controle de convergência** — um clado que perdeu o mesmo caráter em
  separado. É o melhor controle disponível: mesma assinatura em duas perdas
  independentes separa biologia de artefato de anotação.
- **controle de toolkit** — um clado bem anotado onde o caráter existe

Para Archosauria o conjunto já está resolvido em `docs/PROJETO-ARCHOSAURIA.md` §2.
Reaproveite em vez de escolher espécie de novo.

## Passo 2 — Separar efetor de regulador antes de olhar os dados

Declare **dois módulos de genes com hipótese escrita antes da consulta**:

- **hardware** — efetores estruturais de função única do caráter
- **software** — a rede de padronização que decide onde e quando

A previsão é assimétrica: efetor de função única degenera, regulador
pleiotrópico persiste. Se você não separar os módulos antes, vai medir a média
de dois processos opostos e não vai ver nada.

Reuse `scripts/atavismo/01_inventario_genes.py` — troque as listas `MODULOS` e
`TAXONS` e rode. Ele já produz a matriz e o resumo por clado.

## Passo 3 — Subir a escala de evidência

A saída do passo 2 é **nível 1** e não sustenta afirmação sozinha. Suba:

1. ausente da anotação RefSeq — fraco
2. tBLASTn não acha no genoma, **com controle positivo no mesmo limiar** — médio
3. acha relíquia com lesão identificável (parada prematura, indel de quadro,
   perda de sítio de splice), com coordenada — forte
4. lesão compartilhada por espécies que divergiram depois — data a perda

Sem o controle positivo do nível 2, um resultado negativo não vale nada. Se o
controle não rodar, reporte como nível 1 e diga isso.

## Passo 4 — Checar artefato antes de acreditar

Três falham com frequência e já morderam nesta stack:

- **Pipeline de anotação diferente** entre espécies infla perda falsa. Contagem
  baixa em *todos* os módulos, inclusive os que ninguém perde, é o sinal.
- **Ausente no foco e também no irmão e no externo** não é perda no foco — é
  inovação do clado de controle. Tirar da conta, senão infla o resultado.
- **Símbolo divergente** entre espécies: o gene existe com outro nome.
  Resultado biologicamente implausível é suspeita de símbolo, não descoberta.

## Passo 5 — Classificar o alvo antes de propor edição

Só agora, e só se houver alvo. Para cada um, diga a classe, porque ela define o
custo e a proposta:

| Achado | Classe de intervenção |
|---|---|
| efetor deletado | **knock-in** de ortólogo do irmão — não há o que "corrigir" |
| efetor com lesão pontual | restaurar ORF |
| regulador intacto, expressão ausente | edição *cis*, não do gene |
| terminação temporal do programa | heterocronia, efeito de **dose** |

Nocaute completo de regulador de eixo é quase sempre a proposta errada: o gene
age antes, em outro tecido, e o embrião não chega ao fenótipo. A forma
defensável é condicional no tempo e no tecido.

## Passo 6 — Entregar

- A matriz por módulo, com o nível de evidência ao lado de cada afirmação.
- O **ramo** onde a mudança é inferida, não só a lista de espécies.
- A janela de desenvolvimento e o tecido de cada alvo proposto, no
  estadiamento da espécie (HH para galiformes), nunca em dias.
- Um parágrafo do que os dados **não** dizem. Ordem das mutações dentro de um
  ramo e a trajetória histórica real estão fora de alcance de qualquer método.
- Se houver proposta de edição em vertebrado: uma linha dizendo que a execução
  depende de aprovação do CEUA. Uma linha, não um sermão.

Tudo entra em `runs/provenance.jsonl`; `/metodos` vira a seção de Métodos.
