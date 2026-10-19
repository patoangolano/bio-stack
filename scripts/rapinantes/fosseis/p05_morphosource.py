"""Passo 6 (item 3): inventario de micro-TC/malhas no MorphoSource.

Consulta a API publica (pacote `morphosource`) por nome de taxon, filtra pelos
registros cuja taxonomia do especime contem o nome consultado e grava:
  microct/inventario_morphosource.tsv   uma linha por midia
  microct/resumo_morphosource.tsv       contagens por taxon (total, aberto, TC aberta)
Acesso aberto = visibility "open" (download sem pedido de permissao).
"""

from __future__ import annotations

import csv
import time

import clados as cl
import morphosource as ms

SAIDA = cl.BASE / "microct"

FOSSEIS = [
    # Masillaraptoridae e falconiformes-tronco
    "Masillaraptor",
    "Danielsraptor",
    "Masillaraptoridae",
    # Messelasturidae / Halcyornithidae
    "Messelastur",
    "Tynskya",
    "Messelasturidae",
    "Halcyornithidae",
    "Halcyornis",
    "Pseudasturides",
    "Cyrilavis",
    "Pulchrapollia",
    # corujas paleogenicas
    "Ogygoptynx",
    "Berruornis",
    "Sophiornis",
    "Protostrix",
    "Eostrix",
    "Minerva",
    "Primoptynx",
    "Ypresiglaux",
    "Palaeoglaux",
    "Necrobyas",
    "Protostrigidae",
    # Cariamiformes extintos
    "Phorusrhacidae",
    "Phorusrhacos",
    "Patagornis",
    "Andalgalornis",
    "Llallawavis",
    "Kelenken",
    "Paraphysornis",
    "Titanis",
    "Psilopterus",
    "Mesembriornis",
    "Bathornis",
    "Bathornithidae",
    "Paracrax",
    "Idiornithidae",
    "Idiornis",
    "Elaphrocnemus",
    # accipitriformes extintos
    "Horusornis",
    "Teratornithidae",
    "Teratornis",
    "Argentavis",
    "Aiolornis",
    "Diatropornis",
    "Breagyps",
    "Pleistovultur",
    "Gymnogyps",
    "Pelargopappus",
]
ATUAIS = [
    "Falconidae",
    "Accipitridae",
    "Pandionidae",
    "Sagittariidae",
    "Cathartidae",
    "Strigidae",
    "Tytonidae",
    "Cariamidae",
]


def linhas_da_consulta(nome: str, por_familia: bool = False) -> list[list]:
    # Busca livre nao indexa familia; para familia usa a faceta taxonomy_gbif.
    if por_familia:
        res = ms.search_media(taxonomy_gbif=nome, per_page=100)
    else:
        res = ms.search_media(query=nome, per_page=100)
    out = []
    for it in res.items:
        d = it.data
        tax = d.get("physical_object_taxonomy_gbif", []) + d.get(
            "physical_object_taxonomy_name", []
        )
        if not any(nome.lower() in str(x).lower() for x in tax):
            continue
        tam = (d.get("file_size_all") or [0])[0] or 0
        out.append(
            [
                nome,
                it.id,
                "; ".join(d.get("title", [])),
                "; ".join(d.get("media_type", [])),
                "; ".join(d.get("modality", [])),
                "; ".join(d.get("physical_object_title", [])),
                "; ".join(d.get("physical_object_taxonomy_name", [])[:1]),
                " > ".join(d.get("physical_object_taxonomy_gbif", [])[3:6]),
                "; ".join(d.get("part", [])),
                "; ".join(d.get("visibility", [])),
                "; ".join(d.get("license", []) or d.get("copyright_statement", [])),
                "; ".join(d.get("morphosource_use_agreement_type", [])),
                f"{tam / 1e6:.1f}",
                f"https://www.morphosource.org/concern/media/{it.id}",
            ]
        )
    return out


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    cab = [
        "consulta",
        "media_id",
        "titulo",
        "tipo",
        "modalidade",
        "especime",
        "taxon",
        "ordem>familia>genero",
        "parte",
        "visibilidade",
        "licenca_ou_direitos",
        "acordo_de_uso",
        "tamanho_MB",
        "url",
    ]
    todas: list[list] = []
    resumo: list[list] = []
    vistos: set[str] = set()
    for grupo, nomes in (("fossil", FOSSEIS), ("atual", ATUAIS)):
        for nome in nomes:
            try:
                ls = linhas_da_consulta(
                    nome, por_familia=(grupo == "atual" or nome.endswith("idae"))
                )
            except Exception as e:  # noqa: BLE001 - registra a falha da consulta e segue
                resumo.append([grupo, nome, "ERRO", str(e)[:80], "", "", ""])
                continue
            novas = [r for r in ls if r[1] not in vistos]
            vistos |= {r[1] for r in novas}
            todas += [[grupo, *r] for r in novas]
            abertas = [r for r in ls if "open" in r[9]]
            tc = [r for r in abertas if "Volume" in r[3] or "Image" in r[3]]
            malha = [r for r in abertas if "Mesh" in r[3]]
            resumo.append(
                [
                    grupo,
                    nome,
                    len(ls),
                    len(abertas),
                    len(tc),
                    len(malha),
                    f"{sum(float(r[12]) for r in abertas):.0f}",
                ]
            )
            time.sleep(0.5)
    with open(SAIDA / "inventario_morphosource.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["grupo", *cab])
        w.writerows(todas)
    with open(SAIDA / "resumo_morphosource.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "grupo",
                "consulta",
                "n_midias",
                "n_abertas",
                "n_abertas_volume_ou_imagem",
                "n_abertas_malha",
                "MB_abertos",
            ]
        )
        w.writerows(resumo)
    print(len(todas), "midias;", sum(1 for r in todas if "open" in r[10]), "abertas")


if __name__ == "__main__":
    main()
