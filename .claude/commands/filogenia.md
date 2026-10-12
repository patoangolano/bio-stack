---
description: Constrói a filogenia de um conjunto de locos com IQ-TREE 3 (alinhamento, limpeza, ML por loco, supermatriz, concordância, árvore de espécies) e interpreta o resultado.
argument-hint: "<diretório ou FASTA de locos> [outgroup=Taxon] [opções: famsa, trimal, AA, colchete]"
allowed-tools: Read, Write, Edit, Bash, Glob, Grep, Task
---

Entrada pedida: $ARGUMENTS

Delegue ao agente `phylo-analyst`. Ele tem as regras de interpretação; este
comando é a ordem das etapas.

## Passo 1 — Olhar a entrada sem abrir os arquivos

Nunca leia FASTA com Read. Conte com um script:

```bash
uv run python -c "
import gzip, sys
from pathlib import Path
SUF = ('.fa', '.fasta', '.fas', '.fna', '.faa', '.aln', '.fsa')
alvo = Path(sys.argv[1])
if alvo.is_dir():
    # Mesmo filtro do pipeline para pasta: oculto (.gitkeep, ._x.fa) ele pula
    # calado; extensao fora da lista (sem distinguir maiuscula) ele avisa e ignora.
    todos = sorted(f for f in alvo.iterdir() if f.is_file() and not f.name.startswith('.'))
    arqs = [f for f in todos if f.name.lower().removesuffix('.gz').endswith(SUF)]
else:
    # Arquivo unico entra direto, com qualquer extensao (so .csv vira manifesto);
    # extensao que nao e de FASTA fica no nome do loco (ENAM.txt -> ENAM_txt).
    todos = arqs = [alvo]
for f in arqs:
    abrir = gzip.open if f.name.lower().endswith('.gz') else open
    with abrir(f, 'rt', errors='replace') as fh:
        nomes = [l[1:].split()[0] for l in fh if l.startswith('>') and l[1:].split()]
    print(f'{f.name}: {len(nomes)} seqs, {len(set(nomes))} nomes unicos')
fora = [f.name for f in todos if f not in arqs]
if fora:
    print('fora da entrada (o pipeline avisa e ignora):', ', '.join(fora))
" <caminho>
```

Confira: um FASTA por loco, um registro por táxon por loco, e os mesmos
nomes de táxon entre locos (é por nome que a supermatriz casa as linhas).
Cabeçalho do NCBI — RefSeq clássico (`>XP_... descrição [Gallus gallus]`) ou
NCBI Datasets (`>NP_... ENAM [organism=Gallus gallus] [GeneID=...]`) — pede
`--nome-taxon colchete`. Para árvore de proteína a partir de CDS, traduza
antes (`seqkit translate -T <tabela> --trim`): `--tipo-seq AA` não traduz, e
o pipeline recusa DNA forçado como AA.

## Passo 2 — Fixar o grupo externo antes de rodar

Sem grupo externo a árvore não tem direção. Se `$ARGUMENTS` não trouxe um,
proponha o mais próximo fora do grupo de interesse e diga em uma frase por
quê (Regra 3 do agente). Para Archosauria, use a tabela de papéis do
`evo-devo`.

## Passo 3 — Rodar

```bash
uv run python scripts/filogenia/filogenia.py --input <caminho> --outgroup <taxon> [--alinhador famsa] [--limpeza trimal] [--tipo-seq AA] [--nome-taxon colchete]
```

Primeira execução baixa as imagens Docker (~0,3 GB, que ocupam ~1 GB no
disco). Um conjunto de 5 locos e
12 táxons roda em ~1 minuto. Se falhar, o log está em
`runs/filogenia/<id>/nextflow_console.log`; a causa real, em
`~/nf-work/filogenia/<id>/.nextflow.log` no WSL.

## Passo 4 — Ler o resultado

Leia `runs/filogenia/<id>/resumo.tsv` e `resumo.json` (são pequenos) e olhe
as figuras em `figuras/`. Responda:

1. Quais locos ficaram de fora e por quê (`locos_fora`).
2. Concatenada e árvore de espécies concordam? Em quais ramos não?
3. Nas árvores de genes, quais ramos passam SH-aLRT ≥ 80 **e** UFBoot ≥ 95?
   Na supermatriz esse critério não vale (o suporte satura): julgue cada ramo
   pelo gCF e pelo sCFL. gCF ≈ sCF com gDF1 ≈ gDF2 é discordância real,
   compatível com ILS; gCF ≪ sCF é sinal fraco por loco (árvores de genes mal
   estimadas), não conflito biológico (Regra 2 do agente).
4. Algum loco com árvore muito diferente das outras (suspeita de parálogo)?

## Passo 5 — Sensibilidade, se um nó sustenta a conclusão

Rode de novo com `--alinhador famsa` e com `--limpeza trimal` e diga se o nó
sobrevive. Um nó que muda com o alinhador não sustenta conclusão. Se o
`metodos.md` listar loco com limpeza agressiva (> 20 % dos sítios removidos),
rode também `--limpeza nenhuma`.

## Passo 6 — Entregar

- Figura, tabela de suporte dos ramos que importam, e o parágrafo de
  `metodos.md` (não reescreva versões à mão).
- O que a árvore **não** diz, em parágrafo próprio.
- A execução já está em `runs/provenance.jsonl`; `/metodos` a incorpora.
