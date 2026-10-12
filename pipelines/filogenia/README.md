# bio-stack/filogenia

Locos não alinhados → MAFFT/FAMSA → ClipKIT/trimAl → IQ-TREE 3 por loco
(ModelFinder + UFBoot2 + SH-aLRT) → supermatriz MFP+MERGE → gCF/sCFL →
wASTRAL → resumo e parágrafo de Métodos.

```powershell
# do Windows (normal)
uv run python scripts/filogenia/filogenia.py --input <dir> --outgroup <taxon>

# teste ponta a ponta contra a árvore verdadeira
uv run --script tests/check_filogenia.py
```

```bash
# no WSL
nextflow run . -profile test,docker --outdir /tmp/filogenia_teste
nextflow lint main.nf nextflow.config conf/ modules/local/
```

Guia completo, parâmetros, o que foi medido sobre reprodutibilidade e a
validação: [`docs/FILOGENIA.md`](../../docs/FILOGENIA.md).

Layout:

```
main.nf                  workflow
nextflow.config          parametros, perfis (docker, test), relatorios
conf/                    recursos, argumentos por processo, perfil de teste
modules/nf-core/         modulos nf-core (modules.json); iqtree com patch registrado
modules/local/           o que nao existe no nf-core
bin/filogenia_util.py    preparo, concatenacao e resumo (so biblioteca padrao)
assets/test/             locos simulados com AliSim + arvore verdadeira + gerar.sh
```
