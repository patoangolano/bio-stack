"""
f08f — em quantas arvores de genes aparece cada biparticao concorrente dos nos profundos de Telluraves, com qualquer
suporte e com SH-aLRT >= 80 e UFBoot >= 95 (so onde a execucao tem os dois suportes). Os grupos sao restritos aos
taxons de cada gene, como no f08. Serve para nomear as alternativas NNI do gCF: gDF1/gDF2 so podem ser as biparticoes
que aparecem pelo menos tantas vezes quanto elas.

    uv run python scripts/rapinantes/filogenomica/f08f_splits_genes.py <saida.tsv> rotulo=run_id [...]
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

_spec = importlib.util.spec_from_file_location("f08", Path(__file__).with_name("f08_nos_chave.py"))
assert _spec and _spec.loader
f08 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f08)
G = f08.G
G["LEPTO"] = {"Leptosomus_discolor"}

BIPARTICOES = [
    "EUFAL",
    "CARI+FALC",
    "CARI+PSPA",
    "AUST",
    "STRIG+ACCMORPH",
    "STRIG+CORAC",
    "ACCMORPH+CORAC",
    "STRIG+LEPTO",
    "STRIG+ACCMORPH+LEPTO",
    "AFRO",
]


def main() -> None:
    saida = Path(sys.argv[1])
    linhas = []
    for par in sys.argv[2:]:
        rotulo, run_id = par.split("=", 1)
        run = f08.RAIZ / "runs" / "filogenia" / run_id
        rot = (
            json.loads((run / "resumo.json").read_text(encoding="utf-8")).get("rotulos_suporte")
            or []
        )
        i_sh = rot.index("alrt") if "alrt" in rot else None
        i_ufb = rot.index("ufboot") if "ufboot" in rot else None
        cont = {b: [0, 0, 0] for b in BIPARTICOES}  # informativos, presentes, presentes fortes
        for gf in sorted((run / "04_arvores_genes").glob("*.treefile")):
            gt, gb = f08.biparticoes(f08.ler(gf))
            lados = {lado: lab for lado, lab, _ in gb}
            for b in BIPARTICOES:
                c = f08.conj(b, gt)
                if len(c) < 2 or len(gt) - len(c) < 2:
                    continue
                cont[b][0] += 1
                k = f08.canon(c, gt)
                if k not in lados:
                    continue
                cont[b][1] += 1
                if i_sh is None or i_ufb is None:
                    continue
                v = lados[k].split("/")
                try:
                    if float(v[i_sh]) >= 80 and float(v[i_ufb]) >= 95:
                        cont[b][2] += 1
                except (ValueError, IndexError):
                    continue
        fortes_ok = i_sh is not None and i_ufb is not None
        for b, (n, p, f) in cont.items():
            linhas.append(
                {
                    "conjunto": rotulo,
                    "run_id": run_id,
                    "biparticao": b,
                    "genes_informativos": n,
                    "genes_com_biparticao": p,
                    "genes_fortes": f if fortes_ok else "",
                }
            )
    saida.parent.mkdir(parents=True, exist_ok=True)
    with saida.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    print(f"{len(linhas)} linhas -> {saida}")


if __name__ == "__main__":
    main()
