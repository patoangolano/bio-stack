"""Converte docs/paper_reliquias.md em PDF (reportlab, Cambria).

Le um subconjunto de markdown — o que o paper usa e nada alem:

    # titulo            % linha de frente (autor, afiliacao, data)
    ## secao            ### subsecao
    :::resumo Rotulo    ... caixa de resumo ...    :::
    paragrafo, - item, 1. item, | tabela |, ![legenda](figura.png), \\newpage
    inline: **negrito**, *italico*, `codigo`, ^sobrescrito^, \\* literal

Paragrafo que comeca com **Tabela N.** e a legenda da tabela seguinte e fica
presa a ela. Itens de lista dentro de "Referencias" saem com recuo pendente.

    uv run python scripts/paper_pdf.py
"""

from __future__ import annotations

import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Flowable,
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

RAIZ = Path(__file__).resolve().parents[1]
ENTRADA = RAIZ / "docs" / "paper_reliquias.md"
SAIDA = RAIZ / "runs" / "relatorios" / "paper_reliquias_2026-10-17.pdf"
TITULO_CURTO = "Relíquias do esmalte em amniotas sem dente"
FONTES = Path("C:/Windows/Fonts")

TINTA = colors.HexColor("#0b0b0b")
TINTA2 = colors.HexColor("#52514e")
REGUA = colors.HexColor("#9a9890")
CAIXA = colors.HexColor("#f3f2ee")

MARGEM_H, MARGEM_TOPO, MARGEM_BASE = 2.4 * cm, 2.2 * cm, 2.4 * cm
LARGURA = A4[0] - 2 * MARGEM_H


def registra_fontes() -> None:
    # Helvetica nao tem Ψ, ≤ nem →; a Cambria tem todos (conferido com fontTools).
    pdfmetrics.registerFont(TTFont("Cambria", str(FONTES / "cambria.ttc"), subfontIndex=0))
    pdfmetrics.registerFont(TTFont("Cambria-Bold", str(FONTES / "cambriab.ttf")))
    pdfmetrics.registerFont(TTFont("Cambria-Italic", str(FONTES / "cambriai.ttf")))
    pdfmetrics.registerFont(TTFont("Cambria-BoldItalic", str(FONTES / "cambriaz.ttf")))
    pdfmetrics.registerFontFamily(
        "Cambria",
        normal="Cambria",
        bold="Cambria-Bold",
        italic="Cambria-Italic",
        boldItalic="Cambria-BoldItalic",
    )
    pdfmetrics.registerFont(TTFont("Consolas", str(FONTES / "consola.ttf")))


def estilos() -> dict[str, ParagraphStyle]:
    base = ParagraphStyle(
        "corpo",
        fontName="Cambria",
        fontSize=10.5,
        leading=14.6,
        textColor=TINTA,
        alignment=TA_JUSTIFY,
        spaceAfter=6.5,
    )
    return {
        "corpo": base,
        "titulo": ParagraphStyle(
            "titulo",
            parent=base,
            fontName="Cambria-Bold",
            fontSize=17,
            leading=21.5,
            alignment=TA_CENTER,
            spaceAfter=10,
        ),
        "frente": ParagraphStyle(
            "frente",
            parent=base,
            fontSize=10,
            leading=13,
            alignment=TA_CENTER,
            textColor=TINTA2,
            spaceAfter=1.5,
        ),
        "resumo_rot": ParagraphStyle(
            "resumo_rot", parent=base, fontName="Cambria-Bold", fontSize=10.5, spaceAfter=3
        ),
        "resumo": ParagraphStyle("resumo", parent=base, fontSize=9.6, leading=13, spaceAfter=4),
        "h2": ParagraphStyle(
            "h2",
            parent=base,
            fontName="Cambria-Bold",
            fontSize=13,
            leading=16,
            alignment=TA_LEFT,
            spaceBefore=12,
            spaceAfter=5,
            keepWithNext=1,
        ),
        "h3": ParagraphStyle(
            "h3",
            parent=base,
            fontName="Cambria-BoldItalic",
            fontSize=11,
            leading=14,
            alignment=TA_LEFT,
            spaceBefore=8,
            spaceAfter=3.5,
            keepWithNext=1,
        ),
        "legenda": ParagraphStyle(
            "legenda",
            parent=base,
            fontSize=8.8,
            leading=11.6,
            textColor=TINTA2,
            spaceBefore=3,
            spaceAfter=9,
        ),
        "leg_tabela": ParagraphStyle(
            "leg_tabela",
            parent=base,
            fontSize=8.8,
            leading=11.6,
            textColor=TINTA2,
            alignment=TA_LEFT,
            spaceBefore=4,
            spaceAfter=3,
            keepWithNext=1,
        ),
        "celula": ParagraphStyle(
            "celula", parent=base, fontSize=8.3, leading=10.2, alignment=TA_LEFT, spaceAfter=0
        ),
        "item": ParagraphStyle("item", parent=base, leftIndent=14, bulletIndent=3, spaceAfter=3),
        "ref": ParagraphStyle(
            "ref",
            parent=base,
            fontSize=9,
            leading=11.8,
            alignment=TA_LEFT,
            leftIndent=14,
            firstLineIndent=-14,
            spaceAfter=3.5,
        ),
    }


def inline(texto: str) -> str:
    """Markdown inline -> marcacao de Paragraph do reportlab."""
    texto = texto.replace("\\*", "\x00")
    codigos: list[str] = []

    def guarda(m: re.Match[str]) -> str:
        codigos.append(m.group(1))
        return f"\x01{len(codigos) - 1}\x01"

    texto = re.sub(r"`([^`]+)`", guarda, texto)
    texto = texto.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    texto = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", texto)
    texto = re.sub(r"\*(.+?)\*", r"<i>\1</i>", texto)
    texto = re.sub(r"\^(.+?)\^", r"<super>\1</super>", texto)
    # "e = 1,7" e "1,16 × 10" nao quebram no meio da linha
    texto = texto.replace(" = ", " = ").replace(" × ", " × ")
    for i, c in enumerate(codigos):
        c = c.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        texto = texto.replace(f"\x01{i}\x01", f'<font name="Consolas" size="9">{c}</font>')
    return texto.replace("\x00", "*")


def texto_puro(celula: str) -> str:
    return re.sub(r"[*`^\\]", "", celula)


def tabela(linhas: list[str], est: dict[str, ParagraphStyle]) -> Table:
    celulas = [[c.strip() for c in ln.strip().strip("|").split("|")] for ln in linhas]
    celulas = [r for r in celulas if not all(re.fullmatch(r":?-{3,}:?", c) for c in r)]
    ncol = len(celulas[0])
    # Largura proporcional ao texto mais longo de cada coluna, com piso e teto,
    # para que a coluna de rotulo curto nao roube espaco da de texto.
    pesos = [min(max(max(len(texto_puro(r[j])) for r in celulas), 7), 42) for j in range(ncol)]
    larguras = [LARGURA * p / sum(pesos) for p in pesos]
    cab = ParagraphStyle("cab", parent=est["celula"], fontName="Cambria-Bold")
    dados = [
        [Paragraph(inline(c), cab if i == 0 else est["celula"]) for c in r]
        for i, r in enumerate(celulas)
    ]
    t = Table(dados, colWidths=larguras, repeatRows=1, hAlign="LEFT")
    t.setStyle(
        TableStyle(
            [
                ("LINEABOVE", (0, 0), (-1, 0), 0.9, TINTA),
                ("LINEBELOW", (0, 0), (-1, 0), 0.5, TINTA),
                ("LINEBELOW", (0, -1), (-1, -1), 0.9, TINTA),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 2.6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return t


def figura(caminho: Path, legenda: str, est: dict[str, ParagraphStyle]) -> Flowable:
    if not caminho.exists():
        raise SystemExit(f"figura ausente: {caminho} — rode scripts/paper_figuras.py")
    w, h = ImageReader(str(caminho)).getSize()
    larg = LARGURA * 0.98
    alt = larg * h / w
    teto = (A4[1] - MARGEM_TOPO - MARGEM_BASE) * 0.62
    if alt > teto:
        larg, alt = larg * teto / alt, teto
    return KeepTogether(
        [
            Spacer(1, 4),
            Image(str(caminho), width=larg, height=alt),
            Paragraph(inline(legenda), est["legenda"]),
        ]
    )


def caixa(conteudo: list[Flowable]) -> Table:
    t = Table([[conteudo]], colWidths=[LARGURA])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), CAIXA),
                ("LINEABOVE", (0, 0), (-1, 0), 0.6, REGUA),
                ("LINEBELOW", (0, -1), (-1, -1), 0.6, REGUA),
                ("LEFTPADDING", (0, 0), (-1, -1), 11),
                ("RIGHTPADDING", (0, 0), (-1, -1), 11),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    return t


def constroi(md: str, est: dict[str, ParagraphStyle]) -> tuple[list[Flowable], str]:
    historia: list[Flowable] = []
    destino: list[Flowable] = historia  # dentro de :::resumo, os blocos vao para a caixa
    buffer: list[str] = []
    legenda_tab: Paragraph | None = None
    em_refs = False
    titulo = ""
    linhas = md.splitlines()
    i = 0

    def fecha_paragrafo() -> None:
        nonlocal legenda_tab
        if not buffer:
            return
        texto = " ".join(s.strip() for s in buffer)
        buffer.clear()
        if texto.startswith("**Tabela"):
            legenda_tab = Paragraph(inline(texto), est["leg_tabela"])
        else:
            destino.append(
                Paragraph(inline(texto), est["resumo" if destino is not historia else "corpo"])
            )

    while i < len(linhas):
        ln = linhas[i]
        s = ln.strip()
        if not s:
            fecha_paragrafo()
        elif s.startswith("# "):
            fecha_paragrafo()
            titulo = s[2:].strip()
            historia.append(Paragraph(inline(titulo), est["titulo"]))
        elif s.startswith("% "):
            historia.append(Paragraph(inline(s[2:]), est["frente"]))
            if i + 1 < len(linhas) and not linhas[i + 1].startswith("%"):
                historia.append(Spacer(1, 12))
        elif s.startswith(":::resumo"):
            fecha_paragrafo()
            nova: list[Flowable] = [
                Paragraph(inline(s[len(":::resumo") :].strip()), est["resumo_rot"])
            ]
            destino = nova
        elif s == ":::":
            fecha_paragrafo()
            historia.append(caixa(destino))
            historia.append(Spacer(1, 10))
            destino = historia
        elif s.startswith("### "):
            fecha_paragrafo()
            historia.append(Paragraph(inline(s[4:]), est["h3"]))
        elif s.startswith("## "):
            fecha_paragrafo()
            em_refs = s[3:].strip().lower().startswith("referências")
            historia.append(Paragraph(inline(s[3:]), est["h2"]))
        elif s == "\\newpage":
            fecha_paragrafo()
            historia.append(PageBreak())
        elif s.startswith("|"):
            fecha_paragrafo()
            bloco = []
            while i < len(linhas) and linhas[i].strip().startswith("|"):
                bloco.append(linhas[i])
                i += 1
            t = tabela(bloco, est)
            historia.append(KeepTogether([legenda_tab, t]) if legenda_tab else t)
            historia.append(Spacer(1, 8))
            legenda_tab = None
            continue
        elif m := re.fullmatch(r"!\[(.+)\]\(([^)]+)\)", s):
            fecha_paragrafo()
            historia.append(figura(RAIZ / m.group(2), m.group(1), est))
        elif m := re.match(r"(-|\d+\.)\s+(.*)", s):
            fecha_paragrafo()
            texto = m.group(2)
            while i + 1 < len(linhas) and linhas[i + 1].startswith("  ") and linhas[i + 1].strip():
                i += 1
                texto += " " + linhas[i].strip()
            if em_refs:
                historia.append(Paragraph(inline(texto), est["ref"]))
            else:
                marca = "•" if m.group(1) == "-" else m.group(1)
                destino.append(Paragraph(inline(texto), est["item"], bulletText=marca))
        else:
            buffer.append(ln)
        i += 1
    fecha_paragrafo()
    return historia, titulo


def rodape(canvas, doc) -> None:  # noqa: ANN001 — assinatura do callback do reportlab
    canvas.saveState()
    y = MARGEM_BASE - 0.9 * cm
    canvas.setStrokeColor(REGUA)
    canvas.setLineWidth(0.4)
    canvas.line(MARGEM_H, y + 10, A4[0] - MARGEM_H, y + 10)
    canvas.setFont("Cambria-Italic", 8)
    canvas.setFillColor(TINTA2)
    canvas.drawString(MARGEM_H, y, TITULO_CURTO)
    canvas.setFont("Cambria", 8)
    canvas.drawRightString(A4[0] - MARGEM_H, y, str(doc.page))
    canvas.restoreState()


def main() -> int:
    registra_fontes()
    est = estilos()
    historia, titulo = constroi(ENTRADA.read_text(encoding="utf-8"), est)
    SAIDA.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(
        str(SAIDA),
        pagesize=A4,
        leftMargin=MARGEM_H,
        rightMargin=MARGEM_H,
        topMargin=MARGEM_TOPO,
        bottomMargin=MARGEM_BASE,
        title=titulo,
        author="Matheus Lima",
        subject="Manuscrito de trabalho — relíquias de genes do esmalte",
        creator="scripts/paper_pdf.py (reportlab)",
        lang="pt-BR",
    )
    doc.build(historia, onFirstPage=rodape, onLaterPages=rodape)
    print(f"{SAIDA}  ({doc.page} páginas)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
