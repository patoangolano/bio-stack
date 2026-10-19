# referencias.py
"""Citações e referências ABNT do artigo a partir de DOIs resolvidos no doi.org.

Lê reports/artigo/referencias.tsv (chave -> DOI ou referência manual), baixa o
CSL-JSON de cada DOI por negociação de conteúdo no doi.org (Crossref e DataCite) com
cache em reports/artigo/build/csl/, e:

- ``--verificar``: imprime autor, ano e título de cada DOI, para conferência;
- padrão: substitui as marcas de citação de reports/artigo/artigo.md, acrescenta as
  REFERÊNCIAS (NBR 6023:2018) e grava reports/artigo/build/artigo_resolvido.md e o
  recibo reports/artigo/recibos_referencias.tsv.

Marcas de citação (NBR 10520:2023, sistema autor-data):
  [@chave] ou [@a; @b]   -> (Sobrenome, ano; Sobrenome et al., ano)
  [@chave, p. 12]        -> (Sobrenome, ano, p. 12)
  @chave                 -> Sobrenome (ano)

Uso: uv run python scripts/rapinantes/artigo/referencias.py [--verificar]
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
DIR = RAIZ / "reports/artigo"
TABELA = DIR / "referencias.tsv"
ARTIGO = DIR / "artigo.md"
BUILD = DIR / "build"
CACHE = BUILD / "csl"
SAIDA_MD = BUILD / "artigo_resolvido.md"
RECIBO = DIR / "recibos_referencias.tsv"

PARTICULAS = {
    "von",
    "van",
    "de",
    "da",
    "do",
    "dos",
    "das",
    "del",
    "della",
    "der",
    "den",
    "du",
    "le",
    "la",
    "ten",
    "ter",
}
SEM_LOCAL = "[*S. l.*]"
# Campos que o Crossref registra incompletos, corrigidos pela página da editora.
CORRECOES_CSL = {
    # Crossref traz só "1"; academic.oup.com/sysbio/issue/18/1 dá "Pages 1–32" (2026-10-19).
    "10.2307/2412407": {"page": "1-32"},
}


@dataclass
class Ref:
    chave: str
    doi: str
    manual: str
    nota: str
    csl: dict = field(default_factory=dict)
    # preenchidos em preparar()
    sobrenomes: list[str] = field(default_factory=list)
    ano: str = ""
    sufixo: str = ""
    abnt: str = ""
    cit_autor_manual: str = ""
    citada: bool = False

    @property
    def ano_cit(self) -> str:
        return f"{self.ano}{self.sufixo}"


def baixar_csl(doi: str) -> dict:
    """CSL-JSON do doi.org, com cache em disco; falha alta se o DOI não resolver."""
    CACHE.mkdir(parents=True, exist_ok=True)
    arq = CACHE / (re.sub(r"[^A-Za-z0-9._-]", "_", doi) + ".json")
    if arq.exists():
        return json.loads(arq.read_text(encoding="utf-8"))
    url = "https://doi.org/" + urllib.parse.quote(doi, safe="/:;()[]")
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.citationstyles.csl+json",
            "User-Agent": "bio-stack-artigo/1.0 (referencias.py)",
        },
    )
    for tentativa in range(4):
        try:
            with urllib.request.urlopen(req, timeout=40) as resp:
                dado = json.loads(resp.read().decode("utf-8"))
            arq.write_text(json.dumps(dado, ensure_ascii=False, indent=1), encoding="utf-8")
            time.sleep(0.3)
            return dado
        except urllib.error.HTTPError as erro:
            if erro.code == 404:
                raise RuntimeError(f"DOI não encontrado no doi.org: {doi}") from erro
            espera = 2 * (tentativa + 1)
        except (urllib.error.URLError, TimeoutError):
            espera = 2 * (tentativa + 1)
        time.sleep(espera)
    raise RuntimeError(f"doi.org não respondeu após 4 tentativas: {doi}")


def limpar_html(texto: str) -> str:
    """Itálico do Crossref vira *...*; demais marcas somem; entidades são decodificadas."""
    # O Crossref às vezes cola o itálico na palavra vizinha ("Eocene<i>Messelastur</i>").
    texto = re.sub(r"<(?:i|em)>", "\x01", texto, flags=re.I)
    texto = re.sub(r"</(?:i|em)>", "\x02", texto, flags=re.I)
    texto = re.sub(r"<[^>]+>", "", texto)
    texto = html.unescape(texto)
    texto = re.sub(r"(\w)\x01", "\\1 \x01", texto)
    texto = re.sub(r"\x02(\w)", "\x02 \\1", texto)
    texto = re.sub(r"\x01\s+", "\x01", texto)
    texto = re.sub(r"\s+\x02", "\x02", texto)
    texto = texto.replace("\x02\x01", "").replace("\x01", "*").replace("\x02", "*")
    texto = re.sub(r"\s+", " ", texto).strip()
    return texto.replace("**", "")


def separar_sobrenome(familia: str) -> tuple[str, str]:
    """('von Haeseler') -> ('Haeseler', 'von'): partícula minúscula vai para depois do prenome."""
    partes = familia.split()
    part = []
    while len(partes) > 1 and partes[0] in PARTICULAS:
        part.append(partes.pop(0))
    return " ".join(partes), " ".join(part)


def iniciais(prenome: str) -> str:
    saida = []
    for bloco in re.split(r"\s+", prenome.replace(".", ". ").strip()):
        if not bloco:
            continue
        hifen = [p for p in bloco.split("-") if p]
        saida.append("-".join(p[0].upper() + "." for p in hifen))
    return " ".join(saida)


def nome_abnt(autor: dict) -> tuple[str, str]:
    """Nome na referência (SOBRENOME, I. partícula) e sobrenome para a citação."""
    if "literal" in autor:
        return autor["literal"].upper(), autor["literal"]
    familia = autor.get("family", "").strip()
    sobrenome, part = separar_sobrenome(familia)
    if autor.get("non-dropping-particle"):
        part = (autor["non-dropping-particle"] + " " + part).strip()
    if sobrenome.isupper() and len(sobrenome) > 3:
        sobrenome = sobrenome.title()
    ref = sobrenome.upper()
    prenome = autor.get("given", "")
    if prenome:
        ref += ", " + iniciais(prenome)
    if part:
        ref += " " + part
    return ref, sobrenome


def ano_csl(csl: dict) -> str:
    # O ano do fascículo impresso acompanha volume e número; "issued" do Crossref é a
    # data mais antiga (em geral a versão online), que descasaria de v. e n.
    for campo in ("published-print", "issued", "published-online", "created"):
        partes = csl.get(campo, {}).get("date-parts")
        if partes and partes[0] and partes[0][0]:
            return str(partes[0][0])
    return "[s. d.]"


def primeiro(valor) -> str:
    if isinstance(valor, list):
        return valor[0] if valor else ""
    return valor or ""


def autores_ref(ref: Ref) -> str:
    """Bloco de autoria da referência: até três nomes; quatro ou mais, o primeiro + et al."""
    if ref.manual:
        return ref.manual.split(". ")[0]
    autores = ref.csl.get("author") or ref.csl.get("editor") or []
    nomes = [nome_abnt(a)[0] for a in autores]
    if len(nomes) >= 4:
        return nomes[0] + " et al."
    return "; ".join(nomes) or "[AUTOR NÃO INFORMADO]"


def com_ponto(texto: str) -> str:
    """Fecha o elemento com ponto, sem duplicar o de uma abreviatura ou o de '?'/'!'."""
    return texto if texto.endswith((".", "?", "!")) else texto + "."


def formatar_abnt(ref: Ref) -> str:
    csl = ref.csl
    aut = com_ponto(autores_ref(ref))
    titulo = limpar_html(primeiro(csl.get("title", "")))
    sub = limpar_html(primeiro(csl.get("subtitle", "")))
    if sub and sub.lower() not in titulo.lower():
        titulo = f"{titulo}: {sub}"
    titulo = titulo.rstrip(".")
    if titulo.isupper():
        # Título todo em caixa alta no Crossref (Cladistics 1989, Evolution 2008).
        titulo = titulo.capitalize()
    titulo_p = com_ponto(titulo)
    titulo_b = f"**{titulo}**" + ("" if titulo.endswith(("?", "!")) else ".")
    ano = ref.ano_cit
    doi = f"DOI: https://doi.org/{ref.doi}."
    tipo = csl.get("type", "article-journal")
    local = csl.get("publisher-place") or SEM_LOCAL
    editora = csl.get("publisher", "")
    if tipo in ("article-journal", "article", "journal-article"):
        periodico = limpar_html(primeiro(csl.get("container-title", "")))
        partes = [f"**{periodico}**", local]
        if csl.get("volume"):
            partes.append(f"v. {csl['volume']}")
        if csl.get("issue"):
            fasciculo = re.sub(r"Supplement[_ ]?(\d+)", r"supl. \1", str(csl["issue"]))
            partes.append(f"n. {fasciculo}")
        pag = str(csl.get("page") or "").replace("–", "-").replace("--", "-")
        ini, _, fim = pag.partition("-")
        if fim == ini:
            pag = ini
        num_artigo = str(csl.get("article-number") or csl.get("number") or "")
        if pag and re.match(r"^\d", pag):
            partes.append(f"p. {pag}")
        elif pag:
            partes.append(pag)
        elif num_artigo:
            partes.append(num_artigo)
        partes.append(ano)
        return f"{aut} {titulo_p} " + ", ".join(partes) + f". {doi}"
    if tipo == "dataset":
        return (
            f"{aut} {titulo_b} {local}: {editora or '[*s. n.*]'}, {ano}. Conjunto de dados. {doi}"
        )
    if tipo in ("chapter", "book-chapter"):
        livro = limpar_html(primeiro(csl.get("container-title", "")))
        pag = csl.get("page")
        pag_txt = f" p. {pag}." if pag else ""
        return f"{aut} {titulo_p} *In*: **{livro}**. {local}: {editora}, {ano}.{pag_txt} {doi}"
    if tipo == "posted-content":
        return f"{aut} {titulo_p} {local}: {editora}, {ano}. Preprint. {doi}"
    return f"{aut} {titulo_b} {local}: {editora or '[*s. n.*]'}, {ano}. {doi}"


def carregar_tabela() -> dict[str, Ref]:
    refs: dict[str, Ref] = {}
    with TABELA.open(encoding="utf-8") as fh:
        for lin in csv.DictReader(fh, delimiter="\t"):
            r = Ref(lin["chave"], lin["doi"].strip(), lin["abnt_manual"].strip(), lin["nota"])
            if r.manual:
                r.cit_autor_manual = lin["cit_autor"].strip()
                r.ano = lin["ano"].strip()
            if r.chave in refs:
                sys.exit(f"chave duplicada: {r.chave}")
            refs[r.chave] = r
    return refs


def preparar(refs: dict[str, Ref]) -> None:
    for r in refs.values():
        if r.manual:
            r.sobrenomes = [r.cit_autor_manual]
            continue
        r.csl = baixar_csl(r.doi) | CORRECOES_CSL.get(r.doi, {})
        autores = r.csl.get("author") or r.csl.get("editor") or []
        r.sobrenomes = [nome_abnt(a)[1] for a in autores]
        r.ano = ano_csl(r.csl)


def autor_citacao(r: Ref, narrativa: bool) -> str:
    s = r.sobrenomes
    if len(s) >= 4:
        return f"{s[0]} et al."
    if len(s) == 1:
        return s[0]
    if narrativa:
        return ", ".join(s[:-1]) + " e " + s[-1]
    return "; ".join(s)


def chave_ordem(texto: str) -> str:
    sem_acento = unicodedata.normalize("NFKD", texto)
    sem_acento = "".join(c for c in sem_acento if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9 ]", "", sem_acento.lower())


def desambiguar(citadas: list[Ref]) -> None:
    """Mesmo autor de citação e mesmo ano: sufixos a, b... na ordem das referências."""
    grupos: dict[tuple[str, str], list[Ref]] = {}
    for r in citadas:
        grupos.setdefault((autor_citacao(r, False), r.ano), []).append(r)
    for membros in grupos.values():
        if len(membros) < 2:
            continue
        membros.sort(key=lambda r: chave_ordem(formatar_abnt(r) if not r.manual else r.manual))
        for i, r in enumerate(membros):
            r.sufixo = "abcdefghij"[i]


RE_PAREN = re.compile(r"\[(-?@[^\]]+)\]")
RE_NARR = re.compile(r"(?<![\w\[@/])@([A-Za-z][\w-]*\d{4}[\w-]*)")


def processar(texto: str, refs: dict[str, Ref]) -> str:
    # 1ª passada: marcar citadas
    for m in RE_PAREN.finditer(texto):
        for parte in m.group(1).split(";"):
            chave = parte.strip().lstrip("-@").split(",")[0].strip()
            if chave not in refs:
                sys.exit(f"citação sem referência: {chave}")
            refs[chave].citada = True
    for m in RE_NARR.finditer(RE_PAREN.sub("", texto)):
        if m.group(1) not in refs:
            sys.exit(f"citação sem referência: {m.group(1)}")
        refs[m.group(1)].citada = True
    citadas = [r for r in refs.values() if r.citada]
    for r in citadas:
        if not r.manual:
            r.abnt = ""
    desambiguar(citadas)

    def paren(m: re.Match) -> str:
        itens = []
        for parte in m.group(1).split(";"):
            parte = parte.strip()
            so_ano = parte.startswith("-")
            chave, _, loc = parte.lstrip("-@").partition(",")
            r = refs[chave.strip()]
            data = r.ano_cit + (f", {loc.strip()}" if loc.strip() else "")
            itens.append((None if so_ano else autor_citacao(r, False), data))
        # NBR 10520: várias obras na mesma citação vão em ordem alfabética de autoria
        if all(aut is not None for aut, _ in itens):
            itens.sort(key=lambda it: (chave_ordem(it[0] or ""), it[1]))
        # mesmo autor em sequência: (Mayr, 2005, 2011)
        saida: list[str] = []
        for aut, data in itens:
            if saida and aut is not None and saida[-1].startswith(aut + ", "):
                saida[-1] += f", {data}"
            else:
                saida.append(data if aut is None else f"{aut}, {data}")
        return "(" + "; ".join(saida) + ")"

    texto = RE_PAREN.sub(paren, texto)
    texto = RE_NARR.sub(
        lambda m: f"{autor_citacao(refs[m.group(1)], True)} ({refs[m.group(1)].ano_cit})", texto
    )
    return texto


def lista_referencias(refs: dict[str, Ref]) -> list[str]:
    """Ordem alfabética pelo bloco de autoria; o mesmo bloco segue a ordem cronológica,
    de modo que MAYR, G. sozinho vem antes de MAYR, G.; CLARKE, J."""
    citadas = [r for r in refs.values() if r.citada]
    for r in citadas:
        r.abnt = r.manual if r.manual else formatar_abnt(r)
    citadas.sort(key=lambda r: (chave_ordem(autores_ref(r)), r.ano_cit, chave_ordem(r.abnt)))
    return [r.abnt for r in citadas]


def gravar_recibo(refs: dict[str, Ref]) -> None:
    with RECIBO.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "chave",
                "citada",
                "doi",
                "consulta",
                "status",
                "ano",
                "n_autores",
                "primeiro_autor",
                "titulo_resolvido",
                "periodico_ou_editora",
                "nota",
            ]
        )
        for r in refs.values():
            if r.manual:
                w.writerow(
                    [
                        r.chave,
                        "sim" if r.citada else "não",
                        "",
                        "",
                        "manual",
                        r.ano,
                        "",
                        r.cit_autor_manual,
                        "",
                        "",
                        r.nota,
                    ]
                )
                continue
            csl = r.csl
            w.writerow(
                [
                    r.chave,
                    "sim" if r.citada else "não",
                    r.doi,
                    f"https://doi.org/{r.doi}",
                    "resolvido no doi.org (CSL-JSON)",
                    r.ano_cit,
                    len(csl.get("author") or []),
                    r.sobrenomes[0] if r.sobrenomes else "",
                    limpar_html(primeiro(csl.get("title", ""))),
                    limpar_html(primeiro(csl.get("container-title", "")))
                    or csl.get("publisher", ""),
                    r.nota,
                ]
            )


def main() -> None:
    ap = argparse.ArgumentParser(description="Citações e referências ABNT do artigo")
    ap.add_argument("--verificar", action="store_true", help="só imprime os metadados resolvidos")
    args = ap.parse_args()
    refs = carregar_tabela()
    preparar(refs)
    if args.verificar:
        for r in refs.values():
            if r.manual:
                print(f"{r.chave}\tMANUAL\t{r.cit_autor_manual}\t{r.ano}")
                continue
            t = limpar_html(primeiro(r.csl.get("title", "")))[:90]
            print(
                f"{r.chave}\t{r.ano}\t{len(r.sobrenomes)}\t{r.sobrenomes[0] if r.sobrenomes else '?'}\t{r.csl.get('type')}\t{t}"
            )
        return
    texto = ARTIGO.read_text(encoding="utf-8")
    resolvido = processar(texto, refs)
    refs_txt = "\n\n".join("::ref " + linha for linha in lista_referencias(refs))
    resolvido = resolvido.replace("::referencias", refs_txt)
    BUILD.mkdir(parents=True, exist_ok=True)
    SAIDA_MD.write_text(resolvido, encoding="utf-8")
    gravar_recibo(refs)
    nao_citadas = [r.chave for r in refs.values() if not r.citada]
    print(
        f"{sum(r.citada for r in refs.values())} referências citadas; {len(nao_citadas)} na tabela sem citação: {', '.join(nao_citadas) or '-'}"
    )
    print(SAIDA_MD.relative_to(RAIZ))
    print(RECIBO.relative_to(RAIZ))


if __name__ == "__main__":
    main()
