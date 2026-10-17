#!/usr/bin/env bash
# Nivel 3 da escada de evidencia (skill atavismo): LESAO identificavel nas
# reliquias — frameshift e stop em quadro — com miniprot, que alinha proteina
# contra genoma atravessando introns e marca as duas coisas no alinhamento.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s29_lesoes.sh
#   uv run python scripts/atavismo/s29_lesoes.py
#
# Sementes sensiveis (-k 4 -M 0 -n 1 -l 4): com o padrao (-k 6 -n 3) a reliquia
# de ENAM do kakapo nao alinha. Medido em 2026-10-17 que a sensibilidade NAO cria
# lesao falsa: ENAM do Anolis (com dente, ~280 Ma do jacare) sai 1-849 de 862 aa
# com fs=0 st=0 nas duas configuracoes. Esse e o controle negativo de lesao, e
# ele roda de novo aqui em toda especie com dente.
#
# Saida bruta (PAF com cg/cs) em runs/atavismo/lesoes/; a leitura das lesoes e o
# teste de lesao compartilhada ficam no s29_lesoes.py.
set -euo pipefail

P=/mnt/c/Users/Manec/bio-stack/runs/atavismo
OUT="$P/lesoes"
mkdir -p "$OUT"
OPC=(-j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5 -t 8)
echo "miniprot $(miniprot --version)  |  ${OPC[*]}"

alinha() {  # $1=gene  $2=faa com o alvo  $3=intervalo  $4=rotulo
    local q="$OUT/q_$1.faa"
    awk -v id=">$1_Amis" '/^>/{p=($1==id)} p' "$2" > "$q"
    [ -s "$q" ] || { echo "  $1: proteina ausente em $2" >&2; return; }
    miniprot "${OPC[@]}" "$3" "$q" 2>/dev/null > "$OUT/mp_$1_$4.paf" || true
    printf '  %-6s %-5s %3s alinhamentos\n' "$1" "$4" "$(grep -c . "$OUT/mp_$1_$4.paf" || true)"
}

# cluster do s26: um intervalo NRG1..LPL por especie
for g in ENAM AMBN AMTN; do
    for f in "$P"/lncrna/intervalo_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#intervalo_}
        alinha "$g" "$P/lncrna/jacare_proteinas.faa" "$f" "$rot"
    done
done
# painel do s27: um intervalo por gene e especie
for g in MMP20 AMELX DSPP ACP4; do
    for f in "$P"/painel/intervalo_"$g"_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#intervalo_"$g"_}
        alinha "$g" "$P/painel/painel_${g}_proteinas.faa" "$f" "$rot"
    done
done
