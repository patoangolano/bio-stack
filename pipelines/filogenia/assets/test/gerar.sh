#!/usr/bin/env bash
# Regera os locos de teste do pipeline de filogenia.
#
# Os locos sao simulados com o AliSim (IQ-TREE 3) sobre uma arvore conhecida,
# com indels, e depois "desalinhados" (gaps removidos). O pipeline tem de
# realinhar e recuperar a topologia: o teste compara a arvore inferida com
# arvore_verdadeira.nwk pela distancia de Robinson-Foulds.
#
# Uso (no WSL): bash pipelines/filogenia/assets/test/gerar.sh
set -euo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
IMG="quay.io/biocontainers/iqtree:3.1.4--h8471819_0"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

cp "$AQUI/arvore_verdadeira.nwk" "$TMP/"
simula() {  # nome modelo comprimento semente
    docker run --rm -u "$(id -u):$(id -g)" -v "$TMP:/w" -w /w "$IMG" \
        iqtree --alisim "$1" -t arvore_verdadeira.nwk -m "$2" --length "$3" --seed "$4" \
        --indel 0.03,0.03 --indel-size 'POW{1.7/30},POW{1.7/30}' -af fasta -quiet
}
simula loco_A 'GTR{1.2/3.5/0.8/1.1/4.2}+F{0.3/0.2/0.2/0.3}+G4{0.6}' 900 101
simula loco_B 'HKY{3.0}+F{0.25/0.25/0.25/0.25}+G4{0.8}' 1200 202
simula loco_C 'GTR{1.0/2.8/0.9/1.3/3.9}+F{0.28/0.22/0.24/0.26}+I{0.2}+G4{0.7}' 750 303
simula loco_D 'K80{2.5}+G4{0.5}' 1050 404
simula loco_E 'JC' 600 505

desalinha() {  # entrada saida [taxons a remover...]
    local ent="$1" sai="$2"; shift 2
    python3 - "$ent" "$sai" "$@" <<'PY'
import sys
ent, sai, *fora = sys.argv[1:]
regs, nome = [], None
for linha in open(ent):
    linha = linha.strip()
    if linha.startswith(">"):
        nome = linha[1:]
        regs.append([nome, ""])
    elif linha:
        regs[-1][1] += linha.replace("-", "")
with open(sai, "w") as fh:
    for nome, seq in regs:
        if nome in fora:
            continue
        cab = f"{nome} loco_B simulado" if "loco_B" in sai else nome
        fh.write(f">{cab}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")
PY
}
desalinha "$TMP/loco_A.fa" "$AQUI/locos/loco_A.fa"
# Cabecalhos com descricao: exercita o saneamento de nome (--nome_taxon primeiro).
desalinha "$TMP/loco_B.fa" "$AQUI/locos/loco_B.fa"
# Loco sem dois taxons: a supermatriz tem de preencher com gap.
desalinha "$TMP/loco_C.fa" "$AQUI/locos/loco_C.fa" Strigops_habroptila Pelodiscus_sinensis
# Loco sem o outgroup: exercita -o ausente numa particao (bug corrigido no IQ-TREE 3.1.4).
desalinha "$TMP/loco_D.fa" "$AQUI/locos/loco_D.fa" Homo_sapiens
# Loco com 3 taxons: tem de ser descartado com aviso (min_taxa = 4).
python3 - "$TMP/loco_E.fa" "$AQUI/locos/loco_E.fa" <<'PY'
import sys
ent, sai = sys.argv[1:]
manter = {"Gallus_gallus", "Alligator_mississippiensis", "Chrysemys_picta"}
regs, atual = [], None
for linha in open(ent):
    linha = linha.strip()
    if linha.startswith(">"):
        atual = linha[1:]
        regs.append([atual, ""])
    elif linha:
        regs[-1][1] += linha.replace("-", "")
with open(sai, "w") as fh:
    for nome, seq in regs:
        if nome in manter:
            fh.write(f">{nome}\n{seq}\n")
PY
echo "locos de teste gerados em $AQUI/locos"
