"""
Open Tree of Life (API v3): estudos-fonte que sustentam e que conflitam com os
nos-chave das aves de rapina na arvore sintetica.

    uv run python scripts/rapinantes/reanalise/opentree_nos.py

1. tnrs/match_names (contexto Birds) para os taxa exemplares;
2. tree_of_life/induced_subtree com os ott_ids (topologia sintetica dos exemplares);
3. por no-chave: tree_of_life/mrca dos taxa que o definem, conferindo se os taxa
   que deveriam ficar FORA estao fora (senao o no "nao existe" na sintese);
4. tree_of_life/node_info: supported_by, conflicts_with, resolves, partial_path_of, terminal;
5. studies/find_studies para citacao e ano de cada estudo-fonte.

Grava em runs/rapinantes/reanalise/opentree/: respostas brutas (json), nos.tsv,
estudos.tsv, induzida.tre. Somente biblioteca padrao.
"""

from __future__ import annotations

import json
import re
import time
import urllib.request
from pathlib import Path

API = "https://api.opentreeoflife.org/v3"
RAIZ = Path(__file__).resolve().parents[3]
OUT = RAIZ / "runs" / "rapinantes" / "reanalise" / "opentree"

TAXA = [
    # Accipitrimorphae
    "Cathartes aura", "Vultur gryphus", "Coragyps atratus", "Sagittarius serpentarius",
    "Pandion haliaetus", "Aquila chrysaetos", "Haliaeetus albicilla", "Elanus caeruleus",
    "Pernis apivorus", "Accipiter nisus", "Gyps fulvus", "Circaetus pectoralis",
    # Falconiformes
    "Falco peregrinus", "Falco cherrug", "Caracara plancus", "Milvago chimachima",
    "Herpetotheres cachinnans", "Micrastur semitorquatus", "Polihierax semitorquatus",
    "Spiziapteryx circumcincta",
    # Strigiformes
    "Tyto alba", "Phodilus badius", "Strix occidentalis", "Bubo bubo", "Athene noctua",
    "Glaucidium brasilianum", "Asio otus", "Ninox novaeseelandiae",
    # Cariamiformes, Psittacopasserae
    "Cariama cristata", "Chunga burmeisteri", "Nestor notabilis", "Psittacus erithacus",
    "Melopsittacus undulatus", "Acanthisitta chloris", "Passer domesticus", "Corvus corax",
    # Coraciimorphae
    "Colius striatus", "Leptosomus discolor", "Trogon melanurus", "Buceros rhinoceros",
    "Upupa epops", "Merops nubicus", "Coracias garrulus", "Picus viridis", "Galbula dea",
    # externos
    "Gallus gallus", "Columba livia", "Opisthocomus hoazin", "Charadrius vociferus",
    "Egretta garzetta", "Chaetura pelagica",
]

CORAC = ["Colius striatus", "Leptosomus discolor", "Trogon melanurus", "Buceros rhinoceros",
         "Upupa epops", "Merops nubicus", "Coracias garrulus", "Picus viridis", "Galbula dea"]
ACCIPITRI = ["Cathartes aura", "Vultur gryphus", "Coragyps atratus", "Sagittarius serpentarius",
             "Pandion haliaetus", "Aquila chrysaetos", "Haliaeetus albicilla", "Elanus caeruleus",
             "Pernis apivorus", "Accipiter nisus", "Gyps fulvus", "Circaetus pectoralis"]
ACCIPITRIDAE = ["Aquila chrysaetos", "Haliaeetus albicilla", "Elanus caeruleus", "Pernis apivorus",
                "Accipiter nisus", "Gyps fulvus", "Circaetus pectoralis"]
STRIG = ["Tyto alba", "Phodilus badius", "Strix occidentalis", "Bubo bubo", "Athene noctua",
         "Glaucidium brasilianum", "Asio otus", "Ninox novaeseelandiae"]
STRIGIDAE = ["Strix occidentalis", "Bubo bubo", "Athene noctua", "Glaucidium brasilianum",
             "Asio otus", "Ninox novaeseelandiae"]
FALCONIDAE = ["Falco peregrinus", "Falco cherrug", "Caracara plancus", "Milvago chimachima",
              "Herpetotheres cachinnans", "Micrastur semitorquatus", "Polihierax semitorquatus",
              "Spiziapteryx circumcincta"]
CARIAM = ["Cariama cristata", "Chunga burmeisteri"]
PSITTACO = ["Nestor notabilis", "Psittacus erithacus", "Melopsittacus undulatus",
            "Acanthisitta chloris", "Passer domesticus", "Corvus corax"]
EXT = ["Gallus gallus", "Columba livia", "Opisthocomus hoazin", "Charadrius vociferus",
       "Egretta garzetta", "Chaetura pelagica"]
TODOS_TELL = ACCIPITRI + FALCONIDAE + STRIG + CARIAM + PSITTACO + CORAC

# (nome, taxa que definem o no, taxa que tem de ficar fora)
NOS = [
    ("Telluraves", TODOS_TELL, EXT),
    ("Afroaves (Accipitrimorphae+Strigiformes+Coraciimorphae)", ACCIPITRI + STRIG + CORAC,
     FALCONIDAE + CARIAM + PSITTACO + EXT),
    ("Strigiformes+Coraciimorphae", STRIG + CORAC, ACCIPITRI),
    ("Accipitrimorphae+Strigiformes", ACCIPITRI + STRIG, CORAC),
    ("Accipitrimorphae+Coraciimorphae", ACCIPITRI + CORAC, STRIG),
    ("Coraciimorphae", CORAC, ACCIPITRI + STRIG),
    ("Accipitrimorphae (Cathartidae+Accipitriformes)", ACCIPITRI, STRIG + CORAC + FALCONIDAE),
    ("Accipitriformes s.s. (Sagittariidae+Pandionidae+Accipitridae)",
     ["Sagittarius serpentarius", "Pandion haliaetus"] + ACCIPITRIDAE,
     ["Cathartes aura", "Vultur gryphus", "Coragyps atratus"]),
    ("Pandionidae+Accipitridae", ["Pandion haliaetus"] + ACCIPITRIDAE,
     ["Sagittarius serpentarius", "Cathartes aura"]),
    ("Sagittariidae+Accipitridae", ["Sagittarius serpentarius"] + ACCIPITRIDAE,
     ["Pandion haliaetus", "Cathartes aura"]),
    ("Cathartidae", ["Cathartes aura", "Vultur gryphus", "Coragyps atratus"],
     ["Sagittarius serpentarius", "Pandion haliaetus"] + ACCIPITRIDAE),
    ("Accipitridae", ACCIPITRIDAE, ["Pandion haliaetus", "Sagittarius serpentarius"]),
    ("Australaves (Cariamiformes+Falconiformes+Psittacopasserae)", CARIAM + FALCONIDAE + PSITTACO,
     ACCIPITRI + STRIG + CORAC + EXT),
    ("Eufalconimorphae (Falconiformes+Psittacopasserae)", FALCONIDAE + PSITTACO, CARIAM),
    ("Falconiformes+Cariamiformes", FALCONIDAE + CARIAM, PSITTACO),
    ("Psittacopasserae", PSITTACO, FALCONIDAE + CARIAM),
    ("Falconiformes+Psittaciformes", FALCONIDAE + ["Nestor notabilis", "Psittacus erithacus",
                                                   "Melopsittacus undulatus"],
     ["Acanthisitta chloris", "Passer domesticus", "Corvus corax"]),
    ("Falconidae", FALCONIDAE, CARIAM + PSITTACO),
    ("Herpetotheriinae (Herpetotheres+Micrastur)", ["Herpetotheres cachinnans", "Micrastur semitorquatus"],
     ["Falco peregrinus", "Caracara plancus"]),
    ("Falconinae+Polyborinae (Falco+Caracara)", ["Falco peregrinus", "Falco cherrug", "Caracara plancus",
                                                 "Milvago chimachima", "Polihierax semitorquatus",
                                                 "Spiziapteryx circumcincta"],
     ["Herpetotheres cachinnans", "Micrastur semitorquatus"]),
    ("Strigiformes", STRIG, CORAC + ACCIPITRI),
    ("Strigidae", STRIGIDAE, ["Tyto alba", "Phodilus badius"]),
    ("Tytonidae", ["Tyto alba", "Phodilus badius"], STRIGIDAE),
]


def post(caminho: str, corpo: dict, tentativas: int = 4) -> dict:
    dados = json.dumps(corpo).encode()
    for k in range(tentativas):
        try:
            req = urllib.request.Request(f"{API}/{caminho}", data=dados,
                                         headers={"content-type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa: BLE001 - registra e tenta de novo; a ultima falha sobe
            if k == tentativas - 1:
                raise
            print(f"  {caminho}: {e}; nova tentativa")
            time.sleep(3 * (k + 1))
    raise RuntimeError("inalcancavel")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sobre = post("tree_of_life/about", {})
    (OUT / "about.json").write_text(json.dumps(sobre, indent=1))
    tn = post("tnrs/match_names", {"names": TAXA, "context_name": "Birds", "do_approximate_matching": False})
    (OUT / "tnrs.json").write_text(json.dumps(tn, indent=1))
    ott: dict[str, int] = {}
    for r in tn["results"]:
        if r["matches"]:
            m = r["matches"][0]["taxon"]
            ott[r["name"]] = m["ott_id"]
    faltam = [t for t in TAXA if t not in ott]
    print("sem match:", faltam)

    # quais estao na arvore sintetica (taxon pode estar "broken"/fora)
    ind = post("tree_of_life/induced_subtree", {"ott_ids": list(ott.values()), "label_format": "name"})
    (OUT / "induced_subtree.json").write_text(json.dumps(ind, indent=1))
    (OUT / "induzida.tre").write_text(ind.get("newick", "") + "\n")
    ausentes = ind.get("unknown", {}) or {}
    print("fora da sintese (unknown):", ausentes)
    nao_sint = {k for k, v in ott.items() if f"ott{v}" in ausentes}

    linhas = []
    estudos: dict[str, dict] = {}
    for nome, dentro, fora in NOS:
        ids_dentro = [ott[t] for t in dentro if t in ott and t not in nao_sint]
        mr = post("tree_of_life/mrca", {"ott_ids": ids_dentro})
        no = mr["mrca"]["node_id"]
        info = post("tree_of_life/node_info", {"node_id": no, "include_lineage": False})
        # taxa que deviam estar fora: estao sob o no?
        ids_fora = [ott[t] for t in fora if t in ott and t not in nao_sint]
        intrusos = []
        for t in fora:
            if t not in ott or t in nao_sint:
                continue
            m2 = post("tree_of_life/mrca", {"ott_ids": ids_dentro + [ott[t]]})
            if m2["mrca"]["node_id"] == no:
                intrusos.append(t)
        sup = info.get("supported_by", {}) or {}
        con = info.get("conflicts_with", {}) or {}
        res = info.get("resolves", {}) or {}
        par = info.get("partial_path_of", {}) or {}
        ter = info.get("terminal", {}) or {}
        for chave in list(sup) + list(con) + list(res) + list(par) + list(ter):
            sid = chave.split("@")[0]
            estudos.setdefault(sid, {})
        (OUT / f"node_info_{re.sub(r'[^A-Za-z0-9]+', '_', nome)[:60]}.json").write_text(
            json.dumps({"mrca": mr, "node_info": info, "intrusos": intrusos, "ids_fora": ids_fora}, indent=1))
        linhas.append({
            "no": nome, "node_id": no, "num_tips": info.get("num_tips"),
            "taxon_ott": (info.get("taxon") or {}).get("name", ""),
            "existe_na_sintese": "sim" if not intrusos else "nao (inclui: " + "; ".join(intrusos) + ")",
            "n_sustentam": len(sup), "n_conflitam": len(con), "n_resolve": len(res),
            "n_caminho_parcial": len(par), "n_terminal": len(ter),
            "sustentam": ";".join(sorted(sup)), "conflitam": ";".join(sorted(con)),
            "caminho_parcial": ";".join(sorted(par)),
        })
        print(f"{nome}: {no} sust={len(sup)} confl={len(con)} intrusos={intrusos}")

    # citacoes
    for sid in sorted(estudos):
        try:
            r = post("studies/find_studies", {"property": "ot:studyId", "value": sid, "verbose": True})
            m = (r.get("matched_studies") or [{}])[0]
        except Exception as e:  # noqa: BLE001
            m = {"erro": str(e)}
        estudos[sid] = m
    (OUT / "estudos.json").write_text(json.dumps(estudos, indent=1))

    cab = list(linhas[0].keys())
    with open(OUT / "nos.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cab) + "\n")
        for li in linhas:
            f.write("\t".join(str(li[c]) for c in cab) + "\n")
    with open(OUT / "estudos.tsv", "w", encoding="utf-8") as f:
        f.write("study_id\tano\tdoi\treferencia\n")
        for sid, m in sorted(estudos.items()):
            ref = (m.get("ot:studyPublicationReference") or m.get("erro") or "").replace("\t", " ").replace("\n", " ")
            f.write(f"{sid}\t{m.get('ot:studyYear', '')}\t{m.get('ot:studyPublication', '')}\t{ref}\n")
    print("ok:", OUT)


if __name__ == "__main__":
    main()
