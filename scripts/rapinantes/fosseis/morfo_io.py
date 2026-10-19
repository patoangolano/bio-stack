"""Leitura das matrizes morfologicas baixadas para o braco paleontologico.

Cada matriz vira um objeto Matriz com:
  taxa      lista de nomes (sanitizados: so [A-Za-z0-9_])
  linhas    dict taxon -> lista de celulas; cada celula e um frozenset de
            estados (str de um caractere) ou None para '?'/'-'/inaplicavel
  ordenados conjunto de indices (0-based) de caracteres ordenados pelos autores
  descricoes dict indice -> texto do carater (quando a fonte publica a lista)
  fosseis   conjunto de taxa extintos
  externo   lista de taxa do grupo externo (raiz)

Polimorfismos ficam guardados como conjunto; quem exporta para o IQ-TREE
decide o que fazer com eles (aqui: viram '?', ver exportar_iqtree).
"""

from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

RAIZ = Path(__file__).resolve().parents[3]
BRUTOS = RAIZ / "data" / "raw" / "rapinantes" / "morfologia"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


@dataclass
class Matriz:
    id: str
    fonte: str
    taxa: list[str]
    linhas: dict[str, list[frozenset | None]]
    ordenados: set[int] = field(default_factory=set)
    descricoes: dict[int, str] = field(default_factory=dict)
    fosseis: set[str] = field(default_factory=set)
    externo: list[str] = field(default_factory=list)

    @property
    def nchar(self) -> int:
        return len(next(iter(self.linhas.values())))

    def estados_observados(self, j: int) -> set[str]:
        s: set[str] = set()
        for t in self.taxa:
            c = self.linhas[t][j]
            if c is not None:
                s |= c
        return s

    def estados_sem_polimorfismo(self, j: int) -> set[str]:
        """Estados vistos em celulas monomorficas (o que sobra depois de '?')."""
        s: set[str] = set()
        for t in self.taxa:
            c = self.linhas[t][j]
            if c is not None and len(c) == 1:
                s |= c
        return s


def sanear(nome: str) -> str:
    nome = nome.replace("†", "").strip().strip(".").strip()
    nome = re.sub(r"[^A-Za-z0-9_]+", "_", nome)
    return re.sub(r"_+", "_", nome).strip("_")


def _celulas(seq: str) -> list[frozenset | None]:
    out: list[frozenset | None] = []
    i = 0
    seq = seq.replace(" ", "").replace("\t", "")
    while i < len(seq):
        ch = seq[i]
        if ch in "({":
            fim = seq.index(")" if ch == "(" else "}", i)
            estados = frozenset(x for x in seq[i + 1 : fim] if x not in " ,")
            out.append(estados if estados else None)
            i = fim + 1
            continue
        out.append(None if ch in "?-" else frozenset(ch.upper()))
        i += 1
    return out


def _bloco_matrix(texto: str) -> list[tuple[str, str]]:
    # ancorado na linha: "matrix" aparece dentro de rotulos de carater (Livezey & Zusi)
    m = re.search(r"^\s*MATRIX\s*$(.*?)\n\s*;", texto, re.S | re.I | re.M)
    if m is None:
        raise ValueError("bloco MATRIX nao encontrado")
    pares = []
    for linha in m.group(1).split("\n"):
        linha = linha.strip()
        if not linha or linha.startswith("["):
            continue
        if linha.startswith("'"):
            fim = linha.index("'", 1)
            nome, resto = linha[1:fim], linha[fim + 1 :]
        else:
            nome, _, resto = linha.partition(" ")
            if not resto:
                nome, _, resto = linha.partition("\t")
        pares.append((sanear(nome), resto.strip()))
    return pares


def _ordenados_nexus(texto: str) -> set[int]:
    m = re.search(r"TYPESET[^=]*=\s*(.*?);", texto, re.S | re.I)
    if not m:
        return set()
    corpo = m.group(1)
    mo = re.search(r"ord\s*:\s*([^,;]*?)(?:,\s*unord|$)", corpo, re.I | re.S)
    if not mo:
        return set()
    idx: set[int] = set()
    for tok in re.findall(r"\d+\s*-\s*\d+|\d+", mo.group(1)):
        if "-" in tok:
            a, b = (int(x) for x in tok.split("-"))
            idx |= set(range(a - 1, b))
        else:
            idx.add(int(tok) - 1)
    return idx


def ler_nexus(caminho: Path, id_: str, fonte: str) -> Matriz:
    texto = caminho.read_text(encoding="utf-8", errors="replace")
    pares = _bloco_matrix(texto)
    linhas = {n: _celulas(s) for n, s in pares}
    ncs = {len(v) for v in linhas.values()}
    mn = re.search(r"NCHAR\s*=\s*(\d+)", texto, re.I)
    if mn is None:
        raise ValueError(f"{id_}: NCHAR ausente")
    nchar = int(mn.group(1))
    if ncs != {nchar}:
        raise ValueError(f"{id_}: comprimentos {ncs} != NCHAR {nchar}")
    m = Matriz(id_, fonte, [n for n, _ in pares], linhas, ordenados=_ordenados_nexus(texto))
    # CHARLABELS (Chen et al. 2025, formato NDE)
    cl = re.search(r"CHARLABELS(.*?);", texto, re.S | re.I)
    if cl:
        for num, rot in re.findall(r"\[(\d+)\]\s*'([^']*)'", cl.group(1)):
            m.descricoes[int(num) - 1] = rot.strip()
    return m


def _docx_tabelas(caminho: Path) -> list[list[list[str]]]:
    with zipfile.ZipFile(caminho) as z:
        raiz = ET.fromstring(z.read("word/document.xml"))
    tabelas = []
    for tbl in raiz.iter(W + "tbl"):
        linhas = []
        for tr in tbl.iter(W + "tr"):
            celulas = []
            for tc in tr.findall(W + "tc"):
                celulas.append("".join(t.text or "" for t in tc.iter(W + "t")).strip())
            linhas.append(celulas)
        tabelas.append(linhas)
    return tabelas


def _docx_paragrafos(caminho: Path) -> list[str]:
    with zipfile.ZipFile(caminho) as z:
        raiz = ET.fromstring(z.read("word/document.xml"))
    return ["".join(t.text or "" for t in p.iter(W + "t")) for p in raiz.iter(W + "p")]


def _descricoes_numeradas(paragrafos: list[str]) -> dict[int, str]:
    """'12. Texto...' -> {11: 'Texto...'}, pela ordem de aparicao. A numeracao
    impressa tem erro de digitacao na fonte (Mayr 2026, File S1: o 15 saiu
    como '16.', seguido do 16 verdadeiro); a ordem e o que casa com a matriz,
    e quem chama confere que o total bate com NCHAR."""
    out: dict[int, str] = {}
    k = 0
    for p in paragrafos:
        m = re.match(r"\s*(\d+)\.\s+(.*)", p)
        if not m:
            continue
        out[k] = m.group(2).strip()
        k += 1
    return out


def ler_vz2026() -> Matriz:
    """Mayr (2026) Vertebrate Zoology 76: 121-134, File S2 (matriz em tabela .docx)."""
    pasta = BRUTOS / "mayr_2026_vz"
    tabelas = _docx_tabelas(pasta / "File_S2_matriz.docx")
    linhas: dict[str, list] = {}
    ordem: list[str] = []
    fosseis: set[str] = set()
    for tab in tabelas:
        cab = None
        for lin in tab:
            if lin and lin[0] == "Taxa":
                cab = lin
                continue
            if cab is None or not lin or not lin[0] or len(lin) < 5:
                continue
            nome_bruto = lin[0]
            nome = sanear(nome_bruto)
            if "†" in nome_bruto:
                fosseis.add(nome)
            cel = []
            for v in lin[1 : len(cab)]:
                v = v.strip()
                if v in ("?", "", "-", "‒", "–", "—", "‐"):
                    cel.append(None)
                else:
                    cel.append(frozenset(v))
            if nome not in linhas:
                linhas[nome] = []
                ordem.append(nome)
            linhas[nome].extend(cel)
    nchars = {len(v) for v in linhas.values()}
    if len(nchars) != 1:
        raise ValueError(f"VZ26: comprimentos desiguais {nchars}")
    m = Matriz(
        "VZ26", "Mayr 2026, Vertebrate Zoology 76:121-134 (File S2)", ordem, linhas, fosseis=fosseis
    )
    desc = _descricoes_numeradas(_docx_paragrafos(pasta / "File_S1_caracteres.docx"))
    m.descricoes = desc
    if len(desc) != m.nchar:
        raise ValueError(f"VZ26: {len(desc)} descricoes para {m.nchar} caracteres")
    m.ordenados = {j for j, d in desc.items() if "scored as ordered" in d}
    m.externo = ["Anhimidae"]
    return m


def ler_mk22() -> Matriz:
    pasta = BRUTOS / "mayr_kitchener_2022"
    m = ler_nexus(
        pasta / "Mayr_and_Kitchener_Masillaraptoridae.nex",
        "MK22",
        "Mayr & Kitchener 2022, J Vertebr Paleontol 41:e2083515 (figshare 20267628)",
    )
    m.descricoes = _descricoes_numeradas(
        _docx_paragrafos(pasta / "Mayr_and_Kitchener_Supplementary_Data.docx")
    )
    if len(m.descricoes) != m.nchar:
        raise ValueError(f"MK22: {len(m.descricoes)} descricoes para {m.nchar} caracteres")
    m.fosseis = {
        "Masillaraptor",
        "Danielsraptor",
        "Bathornis",
        "Dynamopterus",
        "Psilopterus",
        "Patagornis",
        "Llallawavis",
    }
    m.externo = ["Anhimidae"]
    return m


def ler_chen25() -> Matriz:
    m = ler_nexus(
        BRUTOS / "chen_etal_2025" / "S4_matriz.nex",
        "C25",
        "Chen et al. 2025, Integr Org Biol obaf029 (Suppl. S4)",
    )
    m.fosseis = {"Gansus_yumenensis", "Ichthyornis_dispar"}
    m.externo = ["Gansus_yumenensis"]
    return m


LLOYD = {
    "K19": (
        "Ksepka_etal_2019a.nex",
        "Ksepka, Grande & Mayr 2019, Curr Biol 29:657-663 (via G. T. Lloyd)",
    ),
    "K17": (
        "Ksepka_etal_2017a.nex",
        "Ksepka, Stidham & Williamson 2017, PNAS 114:8047-8052 (via G. T. Lloyd)",
    ),
    "M11": ("Mayr_2011c.nex", "Mayr 2011, J Syst Palaeontol 9:159-171 (via G. T. Lloyd)"),
    "M05": ("Mayr_2005a.nex", "Mayr 2005, J Vertebr Paleontol 25:635-645 (via G. T. Lloyd)"),
    "M16": ("Mayr_2016a.nex", "Mayr 2016, J Paleontol 90:357-374 (via G. T. Lloyd)"),
    "MC03": (
        "Mayr_et_Clarke_2003a.nex",
        "Mayr & Clarke 2003, Cladistics 19:527-553 (via G. T. Lloyd)",
    ),
    "LZ07": (
        "Livezey_et_Zusi_2007a.nex",
        "Livezey & Zusi 2007, Zool J Linn Soc 149:1-95 (via G. T. Lloyd)",
    ),
}

FOSSEIS = {
    "K19": {
        "Cyrilavis_colburnorum",
        "Pseudasturides_macrocephalus",
        "Anneavis_anneae",
        "Tsidiiyazhi_abini",
        "Sandcoleus_copiosus",
        "Oligocolius_psittacocephalon",
        "Masillacolius_brevidactylus",
        "Chascacocolius_oscitans",
        "Selmes_absurdipes",
        "Plesiocathartes_wyomingensis",
        "Plesiocathartes_kelleri",
        "Paracoracias_occidentalis",
        "Eocoracias_brachyptera",
        "Primotrogon_wintersteini",
        "Palaeotodus_itardiensis",
        "Messelirrisor_halcyrostris",
        "Rupelramphastoides_knopfi",
        "Wieslochia_weissi",
        "Zygodactylus_luberonensis",
        "Eozygodactylus_americanus",
        "Primozygodactylus_quintus",
        "Eofringillirostrum_boudreauxi",
        "Eofringillirostrum_parvulus",
        "Morsoravis_sedilis",
        "Morsoravis_sp",
        "Pumiliornis_tessellatus",
        "Psittacopes_lepidus",
        "Eocuculus_cherpinae",
    },
    "M11": {"Halcyornithidae", "Sandcoleidae", "Messelasturidae"},
    "M05": {"Messelastur", "Tynskya", "Pseudasturidae"},
    "M16": {
        "Paracrax",
        "Bathornis",
        "Ameghinornis",
        "Dynamopterus",
        "Elaphrocnemus",
        "Psilopterus",
        "Patagornis",
        "Llallawavis",
    },
    "MC03": {"Apsaravis", "Hesperornis", "Ichthyornis"},
    "LZ07": {
        "Lithornis",
        "Aepyornithidae",
        "Dinornithiformes",
        "Sylviornis",
        "Diatryma",
        "Dromornithidae",
        "Aptornis",
        "Raphus",
        "Pezophaps",
    },
}
FOSSEIS["K17"] = FOSSEIS["K19"] & {
    "Cyrilavis_colburnorum",
    "Pseudasturides_macrocephalus",
    "Anneavis_anneae",
    "Tsidiiyazhi_abini",
    "Sandcoleus_copiosus",
    "Oligocolius_psittacocephalon",
    "Masillacolius_brevidactylus",
    "Chascacocolius_oscitans",
    "Plesiocathartes_wyomingensis",
    "Plesiocathartes_kelleri",
    "Paracoracias_occidentalis",
    "Eocoracias_brachyptera",
    "Palaeotodus_itardiensis",
    "Rupelramphastoides_knopfi",
    "Zygodactylus_luberonensis",
    "Eozygodactylus_americanus",
}

EXTERNO = {
    "K19": ["Eudromia_elegans", "Dendrocygna_bicolor"],
    "K17": ["Eudromia_elegans", "Dendrocygna_bicolor"],
    "M11": ["Tinamidae", "Anseriformes"],
    "M05": ["Tinamidae", "Galliformes"],
    "M16": ["Anhimidae"],
    "MC03": ["Apsaravis", "Hesperornis", "Ichthyornis"],
    "LZ07": [
        "Lithornis",
        "Struthio",
        "Rhea",
        "Casuarius",
        "Dromaius",
        "Apteryx",
        "Eudromia",
        "Aepyornithidae",
        "Dinornithiformes",
    ],
}

# Livezey & Zusi 2007: os 34 primeiros terminais sao nao-neornitinos (de
# Non_Archosauromorpha a Apsaravis). A pergunta e sobre Neoaves; o subconjunto
# Neornithes poupa CPU e remove a raiz mais instavel da matriz.
LZ07_NAO_NEORNITHES = 34


def ler_lloyd(id_: str) -> Matriz:
    arq, fonte = LLOYD[id_]
    m = ler_nexus(BRUTOS / "lloyd" / arq, id_, fonte)
    m.fosseis = set(FOSSEIS.get(id_, set()))
    m.externo = list(EXTERNO[id_])
    if id_ == "LZ07":
        fora = m.taxa[:LZ07_NAO_NEORNITHES]
        m.taxa = m.taxa[LZ07_NAO_NEORNITHES:]
        for t in fora:
            del m.linhas[t]
        # A mesma matriz esta no MorphoBank (P510) com a lista de caracteres;
        # conferido em 2026-10-17: 188 x 2954, 100% das celulas iguais as do Lloyd.
        m.descricoes = lista_morphobank(
            next((BRUTOS / "morphobank_P510").glob("*character_list*.txt"))
        )
        if len(m.descricoes) != m.nchar:
            raise ValueError("LZ07: lista de caracteres do MorphoBank nao bate com NCHAR")
    falta = (m.fosseis | set(m.externo)) - set(m.taxa)
    if falta:
        raise ValueError(f"{id_}: taxa declarados ausentes da matriz: {falta}")
    return m


def lista_morphobank(caminho: Path) -> dict[int, str]:
    """Lista de caracteres do MorphoBank: ' 12. nome : estado (0); estado (1).'"""
    out: dict[int, str] = {}
    for lin in caminho.read_text(encoding="utf-8", errors="replace").splitlines():
        mm = re.match(r"\s*(\d+)\.\s+(.*)", lin)
        if mm:
            out[int(mm.group(1)) - 1] = mm.group(2).strip()
    return out


def ler_degrange21() -> Matriz:
    """Degrange et al. 2021, Pap Palaeontol (Dryornis), MorphoBank P3896, CC0.
    Lista de caracteres em portugues."""
    pasta = BRUTOS / "morphobank_P3896"
    nex = next(p for p in pasta.glob("*.nex"))
    m = ler_nexus(
        nex,
        "D21",
        "Degrange et al. 2021, Pap Palaeontol 7 (doi 10.1002/spp2.1361), MorphoBank P3896",
    )
    m.descricoes = lista_morphobank(next(pasta.glob("*character_list*.txt")))
    if len(m.descricoes) != m.nchar:
        raise ValueError("D21: lista de caracteres nao bate com NCHAR")
    m.fosseis = {"Pampagyps_imperator", "Dryornis_pampeanus"}
    m.externo = [
        "Meleagris_gallopavo",
        "Acrylium_vulturinum",
        "Anseranas_semipalmata",
        "Dendrocygna_viduata",
    ]
    return m


def todas() -> dict[str, Matriz]:
    out = {k: ler_lloyd(k) for k in LLOYD}
    out["VZ26"] = ler_vz2026()
    out["MK22"] = ler_mk22()
    out["C25"] = ler_chen25()
    out["D21"] = ler_degrange21()
    return out
