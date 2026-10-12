"""Gera o relatorio em PDF do projeto Archosauria.

    uv run python scripts/relatorio_archosauria.py

Saida: runs/relatorios/archosauria_<data>.pdf

Todo numero deste relatorio foi medido nesta maquina ou lido em fonte primaria
citada. Nada e estimativa. Onde a medida nao existe, o texto diz que nao existe.

Fontes dos numeros:
  - sintenia e intervalos: runs/atavismo/lncrna/s26_resumo.json e as medidas
    diretas na API do NCBI Datasets registradas em docs/PESQUISA-2026-10-12-*.md
  - alinhamentos: runs/atavismo/lncrna/s26_*.tsv, alvo_tblastn.tsv,
    reciproco_blastx.tsv
  - literatura: DOIs conferidos nos levantamentos da mesma data
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

RAIZ = Path(__file__).resolve().parents[1]
SAIDA = RAIZ / "runs" / "relatorios"

TINTA = colors.HexColor("#1a1a1a")
CINZA = colors.HexColor("#6b6b6b")
REGUA = colors.HexColor("#c8c8c8")
FUNDO = colors.HexColor("#f2f2f2")
ACENTO = colors.HexColor("#8c2d04")

E = getSampleStyleSheet()

TITULO = ParagraphStyle("t", parent=E["Title"], fontName="Helvetica-Bold",
                        fontSize=19, leading=23, textColor=TINTA, spaceAfter=4)
SUB = ParagraphStyle("s", parent=E["Normal"], fontName="Helvetica",
                     fontSize=10.5, leading=14, textColor=CINZA, alignment=1,
                     spaceAfter=18)
H1 = ParagraphStyle("h1", parent=E["Heading1"], fontName="Helvetica-Bold",
                    fontSize=13.5, leading=17, textColor=TINTA,
                    spaceBefore=16, spaceAfter=7)
H2 = ParagraphStyle("h2", parent=E["Heading2"], fontName="Helvetica-Bold",
                    fontSize=11, leading=14, textColor=ACENTO,
                    spaceBefore=11, spaceAfter=5)
P = ParagraphStyle("p", parent=E["Normal"], fontName="Helvetica",
                   fontSize=9.6, leading=13.6, textColor=TINTA,
                   alignment=TA_JUSTIFY, spaceAfter=7)
NOTA = ParagraphStyle("n", parent=P, fontSize=8.7, leading=12,
                      textColor=CINZA, leftIndent=10, spaceAfter=7)
CAIXA = ParagraphStyle("c", parent=P, fontSize=9.4, leading=13,
                       leftIndent=9, rightIndent=9, spaceBefore=4, spaceAfter=4)
MONO = ParagraphStyle("m", parent=E["Normal"], fontName="Courier",
                      fontSize=8.2, leading=10.6, textColor=TINTA,
                      leftIndent=9, spaceAfter=7)
LEGENDA = ParagraphStyle("l", parent=E["Normal"], fontName="Helvetica-Oblique",
                         fontSize=8.3, leading=11, textColor=CINZA,
                         spaceBefore=3, spaceAfter=11)
REF = ParagraphStyle("r", parent=E["Normal"], fontName="Helvetica",
                     fontSize=8.4, leading=11.4, textColor=TINTA,
                     leftIndent=15, firstLineIndent=-15, spaceAfter=4)


def tabela(dados, larguras, destaque_linhas=()):
    t = Table(dados, colWidths=larguras, repeatRows=1, hAlign="LEFT")
    estilo = [
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.3),
        ("LEADING", (0, 0), (-1, -1), 10.6),
        ("TEXTCOLOR", (0, 0), (-1, -1), TINTA),
        ("BACKGROUND", (0, 0), (-1, 0), FUNDO),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, REGUA),
        ("LINEBELOW", (0, -1), (-1, -1), 0.7, REGUA),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
    ]
    for i in destaque_linhas:
        estilo.append(("FONTNAME", (0, i), (-1, i), "Helvetica-Bold"))
    t.setStyle(TableStyle(estilo))
    return t


def caixa(titulo, corpo):
    """Bloco destacado: um achado que muda a leitura do resto."""
    interno = [[Paragraph(f"<b>{titulo}</b><br/><br/>{corpo}", CAIXA)]]
    t = Table(interno, colWidths=[16.4 * cm], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), FUNDO),
        ("LINEBEFORE", (0, 0), (0, -1), 2.2, ACENTO),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
    ]))
    return KeepTogether([Spacer(1, 4), t, Spacer(1, 9)])


def rodape(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7.6)
    canvas.setFillColor(CINZA)
    canvas.drawString(2.2 * cm, 1.35 * cm,
                      "bio-stack  |  Projeto Archosauria  |  "
                      "Psittaciformes x Crocodylia")
    canvas.drawRightString(A4[0] - 2.2 * cm, 1.35 * cm, f"{doc.page}")
    canvas.setStrokeColor(REGUA)
    canvas.setLineWidth(0.4)
    canvas.line(2.2 * cm, 1.75 * cm, A4[0] - 2.2 * cm, 1.75 * cm)
    canvas.restoreState()


def construir() -> Path:
    SAIDA.mkdir(parents=True, exist_ok=True)
    hoje = date.today().isoformat()
    destino = SAIDA / f"archosauria_{hoje}.pdf"

    doc = SimpleDocTemplate(
        str(destino), pagesize=A4,
        leftMargin=2.2 * cm, rightMargin=2.2 * cm,
        topMargin=2.0 * cm, bottomMargin=2.3 * cm,
        title="Projeto Archosauria - perda de denticao e reliquias de esmalte",
        author="Matheus Lima", subject="Genomica comparativa de Archosauria",
    )

    h = []

    # ---------------------------------------------------------------- capa
    h.append(Paragraph("Perda de dentição em Archosauria", TITULO))
    h.append(Paragraph(
        "Relíquias de esmalte no cromossomo Z, sintenia conservada e um índice "
        "de retenção que acompanha o fenótipo<br/>"
        f"<b>Relatório de progresso — {hoje}</b>", SUB))

    h.append(Paragraph("Resumo", H1))
    h.append(Paragraph(
        "As aves perderam os dentes mineralizados uma única vez, há cerca de "
        "116 milhões de anos, no ancestral comum de Neornithes. Os seis genes de "
        "matriz de esmalte e de dentina que sustentam a odontogênese estão "
        "inativados em todas as aves examinadas e intactos no jacaré. Este "
        "projeto pergunta o que sobrou daquele cluster no genoma das aves, onde "
        "sobrou, e se o que sobrou acompanha a perda do caráter.", P))
    h.append(Paragraph(
        "O relatório estabelece quatro coisas, todas medidas e não inferidas. "
        "Primeira: a vizinhança sintênica arcossauria do cluster de esmalte "
        "(<i>NRG1 &rarr; ENAM &middot; AMBN &middot; AMTN &rarr; LPL &rarr; "
        "PSD3</i>) está conservada em ordem e em fita entre jacaré e frango, e "
        "no frango ela cai no cromossomo Z. Segunda: o intervalo que abriga os "
        "três genes no jacaré colapsou para 11,8% do tamanho no frango, e o que "
        "resta ali não contém nenhum gene codificante anotado. Terceira: ENAM e "
        "AMBN deixaram relíquias codificantes detectáveis nesse endereço, com "
        "melhor hit recíproco confirmado e cópia única, e elas não estão "
        "anotadas em nenhuma base. Quarta: um índice de retenção do cluster, "
        "medido em seis espécies, separa sem sobreposição as que têm dente das "
        "que não têm — e a separação não é explicada por distância "
        "filogenética.", P))

    h.append(caixa(
        "O que este relatório não afirma",
        "Nada aqui demonstra seleção sobre as relíquias, nem estabelece a ordem "
        "temporal entre o colapso do intervalo e a perda do caráter, nem atribui "
        "função a qualquer sequência. Essas três perguntas continuam abertas e o "
        "texto diz onde."))

    # ------------------------------------------------------- 1. a pergunta
    h.append(Paragraph("1. A pergunta, e por que ela foi reformulada", H1))
    h.append(Paragraph(
        "A formulação inicial do projeto foi: <i>quais são as mutações e "
        "mudanças que podem ter ocorrido nos últimos 250 milhões de anos que "
        "separam Psittaciformes de Crocodylia</i>. Como está, a pergunta não tem "
        "resposta — nem com computação infinita. A trajetória histórica real de "
        "substituições não é recuperável, e a ordem em que elas ocorreram dentro "
        "de um ramo não é inferível.", P))
    h.append(Paragraph(
        "O que é inferível, e ainda assim é bastante: quais diferenças "
        "moleculares são derivadas em cada linhagem, e em quais genes o padrão "
        "de substituição indica mudança de regime evolutivo desde o ancestral "
        "comum. Reconstrução de genoma ancestral é, nas palavras dos próprios "
        "autores do Cactus, <i>&ldquo;only a hypothetical reconstruction of the "
        "true history&rdquo;</i> — não um registro histórico. Delimitar isso é o "
        "que separa uma afirmação defensável de uma indefensável.", P))
    h.append(Paragraph(
        "A dentição foi escolhida como caráter focal por três razões "
        "verificadas: tem data com intervalo de confiança, tem painel de genes "
        "estabelecido na literatura, e tem controle positivo funcional — o "
        "jacaré, no qual os mesmos genes estão intactos. Nenhum conjunto de "
        "mamíferos oferece esse controle para a comparação arcossauria.", P))

    # ------------------------------------------------- 2. o que se sabia
    h.append(Paragraph("2. O que a literatura já havia estabelecido", H1))
    h.append(Paragraph(
        "Meredith e colaboradores (2014) analisaram seis genes essenciais à "
        "odontogênese em 48 genomas de aves, cobrindo quase todas as ordens "
        "viventes, mais o jacaré-americano. Encontraram mutações inativadoras "
        "compartilhadas pelas 48 aves, em genes de esmalte e de dentina, o que "
        "indica perda única no ancestral comum. A data estimada é de cerca de "
        "116 milhões de anos, com intervalo de confiança de 95% entre 105,7 e "
        "128,3 Ma.", P))

    h.append(tabela([
        ["Gene", "Papel", "48 aves", "Alligator"],
        ["DSPP", "dentina", "inativado", "intacto"],
        ["AMTN", "matriz de esmalte", "inativado", "intacto"],
        ["AMBN", "matriz de esmalte", "inativado", "intacto"],
        ["ENAM", "matriz de esmalte", "inativado", "intacto"],
        ["AMELX", "matriz de esmalte", "inativado", "intacto"],
        ["MMP20", "matriz de esmalte", "inativado", "intacto"],
    ], [2.6 * cm, 5.2 * cm, 4.3 * cm, 4.3 * cm]))
    h.append(Paragraph(
        "Tabela 1. O painel canônico de Meredith et al. (2014). Springer et al. "
        "(2019) acrescenta ODAM, inativado em placentários sem dente; ODAPH e "
        "ACP4 completam nove candidatos. O projeto rastreava quatro; os cinco "
        "restantes foram acrescentados neste ciclo.", LEGENDA))

    h.append(Paragraph(
        "A perda de dentição não é evento isolado nas aves. Nas fontes "
        "consultadas estão documentadas <b>doze origens evolutivas "
        "independentes</b> de perda de dente ou de esmalte em amniotas: "
        "Neornithes, Testudines, Pholidota, Mysticeti, o ancestral de Pilosa, "
        "Vermilingua, Cingulata, Tubulidentata, Kogiidae, Monodontidae, "
        "Tachyglossidae e Sirenia. Duas dessas origens estão representadas no "
        "conjunto de dados deste projeto — aves e tartarugas — e a segunda havia "
        "entrado originalmente apenas como grupo externo para polarização.", P))

    h.append(Paragraph("2.1 A lacuna metodológica", H2))
    h.append(Paragraph(
        "Seis consultas literais no Europe PMC não localizaram nenhum trabalho "
        "que tenha aplicado análise de <b>taxa evolutiva relativa</b> — "
        "RERconverge ou equivalente — à perda de dentição ou de esmalte. Os três "
        "trabalhos existentes sobre o tema usam detecção direta de mutação "
        "inativadora e razões dN/dS no PAML. O quase-caso mais próximo, Abramyan "
        "et al. (2025), comparou modelos de &omega; em 53 espécies com grupos "
        "fixados pelo pesquisador — que é inferência diferente de resíduo "
        "normalizado pelo fundo genômico com nulo empírico por permulação.", P))

    h.append(PageBreak())

    # --------------------------------------------- 3. o que foi medido
    h.append(Paragraph("3. O que foi medido", H1))
    h.append(Paragraph(
        "Todas as coordenadas vêm da API do NCBI Datasets e da REST do Ensembl, "
        "consultadas em 12 de outubro de 2026. Os assemblies são "
        "<font face='Courier'>rAllMis1</font> para o jacaré e "
        "<font face='Courier'>GRCg7b</font> para o frango.", P))

    h.append(Paragraph("3.1 A âncora sintênica, e por que ela precisou de quatro genes", H2))
    h.append(Paragraph(
        "Duas âncoras publicadas apontam para cromossomos diferentes no frango. "
        "Sire et al. (2008) usam SULT1E1 e IGJ e chegam ao cromossomo 4; são "
        "âncoras de sintenia <i>mamífera</i>. Este projeto usa LPL e PSD3, que "
        "são as âncoras <i>arcossaurias</i> — os vizinhos imediatos do cluster no "
        "jacaré — e chega ao cromossomo Z. As duas não se contradizem: "
        "Al-Hashimi et al. (2010) já havia reportado que o pseudogene de "
        "enamelina da galinha aparece em região não esperada pela sintenia de "
        "mamífero, sugerindo translocação em ancestral sauropsídeo.", P))
    h.append(Paragraph(
        "Para eliminar a ambiguidade, a âncora foi estendida de dois para quatro "
        "genes, acrescentando NRG1 do lado oposto. A ordem no cromossomo é a "
        "mesma nas duas espécies, e a fita também é a mesma — negativa.", P))

    h.append(tabela([
        ["Elemento", "Alligator, chr3", "Gallus, chrZ"],
        ["NRG1", "156.331.250 - 156.595.276", "54.642.186 - 55.109.784"],
        ["ENAM", "157.293.890 - 157.307.957", "ausente da anotação"],
        ["AMBN", "157.326.619 - 157.341.756", "ausente da anotação"],
        ["AMTN", "157.359.292 - 157.374.521", "ausente da anotação"],
        ["LPL", "157.488.689 - 157.512.736", "55.214.993 - 55.231.804"],
        ["PSD3", "157.601.092 - 157.780.011", "55.231.641 - 55.383.432"],
    ], [3.0 * cm, 6.7 * cm, 6.7 * cm], destaque_linhas=(1, 5, 6)))
    h.append(Paragraph(
        "Tabela 2. Sintenia conservada. A ordem NRG1 &rarr; cluster &rarr; LPL "
        "&rarr; PSD3 e a fita negativa são idênticas nas duas espécies. Isso "
        "torna o endereço do cluster no frango inequívoco.", LEGENDA))

    h.append(Paragraph("3.2 O colapso do intervalo", H2))
    h.append(tabela([
        ["Intervalo", "Alligator", "Gallus", "Razão"],
        ["NRG1 - LPL", "893.413 bp", "105.209 bp", "11,8%"],
        ["LPL - PSD3", "88.356 bp", "-797 bp (sobrepostos)", "colapso total"],
    ], [3.8 * cm, 4.2 * cm, 5.2 * cm, 3.2 * cm], destaque_linhas=(1,)))
    h.append(Paragraph(
        "Tabela 3. O intervalo que abriga os três genes de esmalte no jacaré tem, "
        "no frango, 11,8% do tamanho. São 788.204 pares de base ausentes no "
        "endereço sintênico exato.", LEGENDA))
    h.append(Paragraph(
        "O conteúdo do intervalo no frango são <b>seis lncRNAs anotados e nenhum "
        "gene codificante</b>. E a verificação negativa é dupla: nenhum dos nove "
        "genes do painel está anotado em <i>Gallus gallus</i>, nem no NCBI, em "
        "nenhum assembly, nem no Ensembl.", P))
    h.append(Paragraph(
        "Isso confirma, por outro caminho, o relato de Sire et al. (2008) de que "
        "AMBN e ENAM <i>desapareceram</i> após rearranjo cromossômico na galinha. "
        "Desapareceram, não pseudogenizaram. Deleção de DNA é o modo dominante "
        "de perda gênica, e não pseudogenização — Crombez et al. (2025) mostram "
        "em simulação que, na ausência de deleção, pseudogenes permaneceriam "
        "detectáveis por muito mais tempo do que se observa.", P))

    # ----------------------------------------------- 4. as reliquias
    h.append(Paragraph("4. As relíquias", H1))
    h.append(Paragraph(
        "Se os genes desapareceram da anotação, resta perguntar se desapareceram "
        "da sequência. O teste foi alinhamento de proteína contra genoma "
        "(<font face='Courier'>tblastn</font>), usando as proteínas do jacaré "
        "como consulta e o intervalo do frango como sujeito, com dois controles.", P))

    h.append(Paragraph("4.1 Os controles vêm primeiro", H2))
    h.append(tabela([
        ["Controle", "Resultado", "Leitura"],
        ["Sanidade: proteínas do\njacaré vs. próprio genoma",
         "ENAM 99,9%; AMBN, AMTN\ne LPL a 100%", "a janela está correta"],
        ["POSITIVO: LPL do jacaré\nvs. chrZ do frango",
         "e = 1,16e-57; 64,9% em\n154 aa; miniprot id = 0,829",
         "o método tem sensibilidade\na 116 Myr"],
        ["NEGATIVO: DSPP do jacaré\n(chr2) vs. chrZ do frango",
         "e = 8,6; 19 aa; 22,7 bits",
         "piso de ruído em ~21-24 bits"],
    ], [5.0 * cm, 5.7 * cm, 5.7 * cm]))
    h.append(Paragraph(
        "Tabela 4. Sem o controle positivo, um resultado negativo sobre o cluster "
        "seria indistinguível de falta de sensibilidade do método. LPL está "
        "anotado nas duas espécies e foi recuperado com folga, o que torna o "
        "negativo interpretável.", LEGENDA))

    h.append(Paragraph("4.2 O resultado no frango", H2))
    h.append(tabela([
        ["Gene", "Melhor e-value", "Identidade", "HSPs", "Recíproco", "Veredito"],
        ["AMBN", "1,40e-09", "54,8% em 31 aa", "3", "confirma AMBN", "relíquia"],
        ["ENAM", "1,08e-07", "41,2% em 68 aa", "8", "confirma ENAM", "relíquia"],
        ["AMTN", "1,6", "26,2% em 65 aa", "-", "-", "não detectada"],
    ], [1.9 * cm, 2.9 * cm, 3.4 * cm, 1.3 * cm, 3.4 * cm, 3.5 * cm],
        destaque_linhas=(1, 2)))
    h.append(Paragraph(
        "Tabela 5. Relíquias codificantes no intervalo do chrZ. Para a região de "
        "ENAM, oito HSPs apontam para ENAM, com o concorrente mais próximo a "
        "e = 1,7. Cada região candidata dá exatamente um hit contra a própria "
        "janela, o que exclui elemento repetitivo.", LEGENDA))

    h.append(caixa(
        "As relíquias não são os lncRNAs — e a ordem delas sobreviveu",
        "A hipótese de trabalho era que as relíquias fossem os lncRNAs anotados, "
        "com a transcrição explicando a persistência. Não são. A relíquia de ENAM "
        "está em 55.155.480-55.157.241 e a de AMBN em 55.168.855-55.168.992, "
        "ambas em sequência intergênica que nenhuma base anotou — nem como gene, "
        "nem como pseudogene, nem como lncRNA. E <b>ENAM vem antes de AMBN</b> no "
        "cromossomo, que é a mesma ordem do jacaré. Ruído não preserva ordem."))

    h.append(Paragraph(
        "Um candidato foi reprovado, e é assim que se sabe que o teste funciona. "
        "O lncRNA ENSGALG00010014557 havia dado um hit nucleotídico forte dentro "
        "do span de AMTN no jacaré — 72,8% de identidade em 158 bp, e = 4,41e-17. "
        "Traduzida, a mesma região não acerta nada: melhor e-value 1,2, em 8 "
        "aminoácidos, e o acerto é em AMBN, não em AMTN. Não é relíquia de "
        "codificante. É elemento não-codificante conservado ou repetição "
        "compartilhada.", P))

    h.append(PageBreak())

    # --------------------------------- 5. o comparativo
    h.append(Paragraph("5. O teste comparativo em seis espécies", H1))
    h.append(Paragraph(
        "Uma relíquia em uma espécie pode ser particularidade de montagem. A "
        "medida foi então repetida no conjunto da seção 2 do plano do projeto, "
        "com duas origens independentes de perda e três controles positivos.", P))

    h.append(Paragraph(
        "Primeiro, o estado da anotação. O painel de nove genes foi consultado "
        "espécie por espécie:", P))
    h.append(tabela([
        ["Espécie", "Dentição", "Genes do painel anotados"],
        ["Alligator mississippiensis", "com dente", "7 de 9"],
        ["Anolis carolinensis", "com dente", "6 de 9"],
        ["Gavialis gangeticus", "com dente", "5 de 9 (scaffolds não posicionados)"],
        ["Chrysemys picta", "sem dente", "1 de 9 — somente ENAM, no chr6"],
        ["Gallus gallus", "sem dente", "0 de 9"],
        ["Melopsittacus undulatus", "sem dente", "0 de 9"],
        ["Strigops habroptila", "sem dente", "0 de 9"],
    ], [6.3 * cm, 3.1 * cm, 7.0 * cm], destaque_linhas=(4,)))
    h.append(Paragraph(
        "Tabela 6. A tartaruga é desdentada e mantém ENAM anotado; nenhuma das "
        "três aves mantém gene algum do painel. O contraste é mensurável e "
        "aponta decaimento mais avançado nas aves.", LEGENDA))

    h.append(Paragraph("5.1 O índice de retenção", H2))
    h.append(Paragraph(
        "Comparar e-values entre espécies é inválido: o e-value escala com o "
        "tamanho do banco, e as janelas analisadas variam de 710 kb a 1,59 Mb. "
        "O <i>bitscore</i> não depende do banco. O índice abaixo é a soma dos "
        "bitscores dos três genes do cluster, dividida pelo bitscore do controle "
        "positivo da mesma espécie — o que remove também o que sobra de "
        "diferença de qualidade de montagem e de distância ao jacaré.", P))

    h.append(tabela([
        ["Espécie", "Dentição", "ENAM", "AMBN", "AMTN", "LPL (+)", "Índice"],
        ["Alligator mississippiensis", "com dente", "1590", "349", "314", "173", "13,02"],
        ["Gavialis gangeticus", "com dente", "1450", "326", "73,6", "170", "10,88"],
        ["Anolis carolinensis", "com dente", "241", "67,4", "73,6", "152", "2,51"],
        ["Chrysemys picta", "sem dente", "200", "53,9", "26,2", "152", "1,84"],
        ["Strigops habroptila", "sem dente", "106", "24,6", "37,4", "194", "0,87"],
        ["Gallus gallus", "sem dente", "42,0", "56,6", "23,5", "200", "0,61"],
    ], [5.2 * cm, 2.5 * cm, 1.6 * cm, 1.5 * cm, 1.5 * cm, 1.7 * cm, 1.7 * cm],
        destaque_linhas=(1, 2, 3, 4, 5, 6)))
    h.append(Paragraph(
        "Tabela 7. Bitscore do melhor HSP por gene, e o índice de retenção do "
        "cluster. As três espécies com dente ficam entre 2,51 e 13,02; as três "
        "sem dente, entre 0,61 e 1,84. Não há sobreposição.", LEGENDA))

    h.append(caixa(
        "Fenótipo vence filogenia — e é Anolis que mostra isso",
        "Entre as seis espécies, <i>Anolis carolinensis</i> é a mais distante do "
        "jacaré, que é a espécie de onde vêm as proteínas de consulta. Se o "
        "índice medisse apenas distância filogenética, Anolis seria o menor de "
        "todos. Ele é <b>2,51</b> — acima da tartaruga (1,84) e muito acima das "
        "duas aves (0,87 e 0,61). Anolis tem dentes. O que separa os grupos não é "
        "a distância; é o caráter."))

    h.append(Paragraph(
        "Dois detalhes finos, que só aparecem com mais de uma ave. AMBN tem "
        "bitscore 56,6 no frango e apenas 24,6 no kakapo. AMTN inverte: 23,5 no "
        "frango e 37,4 no kakapo. <b>Genes diferentes decaíram em linhagens "
        "diferentes de aves</b>, o que é o que se espera de degradação "
        "estocástica independente depois de uma perda funcional compartilhada. "
        "A perda do caráter é única; o apagamento do rastro não é.", P))

    # ---------------------------------------- 6. conclusoes
    h.append(Paragraph("6. Conclusões", H1))
    for n, (t, c) in enumerate([
        ("O endereço está estabelecido por medida.",
         "A vizinhança sintênica arcossauria do cluster de esmalte mapeia no "
         "cromossomo Z de <i>Gallus gallus</i>, entre NRG1 e LPL, com ordem e "
         "fita conservadas em relação ao jacaré. A divergência aparente com a "
         "literatura decorre de usar âncoras mamíferas em vez de arcossaurias, e "
         "não de erro de coordenada."),
        ("O colapso é grande e está no lugar exato.",
         "O intervalo caiu para 11,8% do tamanho — 788.204 pares de base — e o "
         "que resta não contém gene codificante anotado."),
        ("Duas relíquias codificantes sobreviveram e ninguém as anotou.",
         "ENAM e AMBN são detectáveis por tblastn com melhor hit recíproco "
         "confirmado, cópia única e ordem preservada. Não estão em nenhuma base, "
         "em nenhuma categoria."),
        ("A retenção do cluster acompanha a dentição em seis espécies.",
         "Separação sem sobreposição entre espécies com e sem dente, não "
         "explicada por distância filogenética, com as duas origens "
         "independentes de perda apresentando retenção reduzida."),
        ("O decaimento é graduado, não binário.",
         "Tartaruga menos degradada que as aves; entre as aves, genes diferentes "
         "decaíram em linhagens diferentes. Isso dá um fenótipo molecular "
         "contínuo onde antes havia só presença e ausência."),
    ], start=1):
        h.append(Paragraph(f"<b>{n}. {t}</b> {c}", P))

    h.append(Paragraph("6.1 O que permanece aberto", H2))
    h.append(Paragraph(
        "<b>Seleção não foi testada.</b> A persistência das relíquias tem três "
        "explicações concorrentes e nenhuma foi medida: seleção residual por "
        "função pleiotrópica; baixa taxa de deleção no ambiente genômico local; e "
        "o regime de recombinação do cromossomo Z, que não é ambiente neutro — a "
        "degeneração de cromossomo sexual começa quando o crossing-over cessa, "
        "com taxa inicialmente rápida e depois muito mais lenta. Como as "
        "relíquias se mostraram não transcritas, as duas explicações de ambiente "
        "genômico ganharam peso relativo, e o mapa de recombinação do chrZ deixou "
        "de ser item acessório.", P))
    h.append(Paragraph(
        "<b>A ordem temporal é desconhecida.</b> Não há como dizer, com estes "
        "dados, se o colapso do intervalo precedeu, acompanhou ou seguiu a perda "
        "funcional dos genes.", P))
    h.append(Paragraph(
        "<b>Taxa evolutiva relativa ainda não foi aplicada.</b> A lacuna "
        "identificada permanece aberta, e o obstáculo foi medido: com apenas as "
        "sete espécies do núcleo do projeto, o caráter é colinear com a "
        "bipartição Aves/Crocodylia e não há espaço de permulação. O pacote "
        "RERconverge 0.3.0 instalado tem <font face='Courier'>min.sp = 10</font> "
        "como padrão, e o precedente publicado opera com 30 a 190 pontas, embora "
        "exista um caso com sete. O caminho é ampliar para Amniota, onde as doze "
        "origens independentes dão réplica evolutiva, estimar o fundo com o "
        "genoma inteiro e ler os nove genes dentro do ranking.", P))

    h.append(Paragraph("6.2 Próximos passos, em ordem de custo", H2))
    h.append(tabela([
        ["#", "Passo", "Custo", "Depende de"],
        ["1", "Acrescentar AMELX, DSPP, MMP20, ODAM, ODAPH e ACP4\nao teste de relíquia nas seis espécies",
         "horas", "nada"],
        ["2", "Testar a correlação entre grau de decaimento e distância\nà borda da deleção",
         "horas", "passo 1"],
        ["3", "Conferir cobertura e qualidade da montagem do chrZ\nem 55,15-55,17 Mb",
         "horas", "nada"],
        ["4", "Obter o mapa de recombinação do chrZ e separar\nregião de alta e de baixa recombinação",
         "dias", "literatura"],
        ["5", "Ampliar para Amniota e rodar taxa evolutiva relativa\ncom permulações",
         "semanas", "download e anotação"],
    ], [0.9 * cm, 8.6 * cm, 2.3 * cm, 4.6 * cm]))
    h.append(Paragraph(
        "Tabela 8. Os três primeiros passos não exigem genoma novo nem dado "
        "novo: são execuções dos mesmos scripts com outros parâmetros.", LEGENDA))

    h.append(PageBreak())

    # ------------------------------------------- 7. metodos
    h.append(Paragraph("7. Métodos", H1))
    h.append(Paragraph(
        "<b>Dados.</b> Coordenadas e sequências obtidas da API do NCBI Datasets "
        "(v2alpha), E-utilities e REST do Ensembl, em 12 de outubro de 2026. "
        "Assemblies: <i>A. mississippiensis</i> rAllMis1 (NC_081826.1, chr3); "
        "<i>G. gallus</i> GRCg7b (NC_052572.1, chrZ); <i>C. picta</i> "
        "GCF_011386835.1; <i>S. habroptila</i> GCF_004027225.2; "
        "<i>M. undulatus</i> GCF_012275295.1; <i>G. gangeticus</i> "
        "GCA_030020295.1; <i>A. carolinensis</i> GCF_035594765.1.", P))
    h.append(Paragraph(
        "<b>Alinhamento.</b> BLAST+ 2.17.0 (<font face='Courier'>tblastn</font>, "
        "<font face='Courier'>blastx</font>, <font face='Courier'>blastn</font> "
        "nas tarefas <font face='Courier'>dc-megablast</font> e "
        "<font face='Courier'>blastn</font> com word_size 7), miniprot 0.18-r281, "
        "seqkit. Filtro de baixa complexidade ativo "
        "(<font face='Courier'>-seg yes</font> para consulta proteica, "
        "<font face='Courier'>-dust yes</font> para nucleotídica). Limiar de "
        "e-value deliberadamente frouxo em 10, para que hits marginais fossem "
        "visíveis e julgados em vez de filtrados em silêncio.", P))
    h.append(Paragraph(
        "<b>Controles.</b> Três por execução: sanidade (consulta contra o próprio "
        "genoma), positivo (LPL, anotado em todas as espécies) e negativo (DSPP, "
        "em outro cromossomo no jacaré e de baixa complexidade). Verificação de "
        "cada candidato por número de cópias na própria janela e por melhor hit "
        "recíproco contra o conjunto de proteínas de consulta.", P))
    h.append(Paragraph(
        "<b>Comparação entre espécies por bitscore</b>, não por e-value, porque "
        "o e-value escala com o tamanho do banco e as janelas diferem em até 2,2 "
        "vezes. Normalização pelo bitscore do controle positivo da mesma espécie.", P))
    h.append(Paragraph(
        "<b>Reprodutibilidade.</b> Scripts versionados em "
        "<font face='Courier'>scripts/atavismo/</font>: "
        "<font face='Courier'>s24_lncrna_relicto.py</font> e "
        "<font face='Courier'>s24_alinhar.sh</font> (teste no frango), "
        "<font face='Courier'>s25_verificar_relicto.sh</font> (verificação), "
        "<font face='Courier'>s26_relicto_comparativo.py</font> e "
        "<font face='Courier'>s26_alinhar.sh</font> (comparativo). Saídas em "
        "<font face='Courier'>runs/atavismo/lncrna/</font>. Toda operação está "
        "registrada em <font face='Courier'>runs/provenance.jsonl</font>.", P))

    h.append(Paragraph("7.1 Limitações do método", H2))
    h.append(Paragraph(
        "A evidência de relíquia é de pouca sequência: 41% de identidade em 68 "
        "aminoácidos para ENAM e 55% em 31 para AMBN. O que sustenta a afirmação "
        "não é a identidade isolada, e sim a conjunção de seis propriedades — "
        "agrupamento dos HSPs em menos de 2 kb, cópia única na própria janela, "
        "melhor hit recíproco confirmado, posição prevista <i>a priori</i> por "
        "sintenia de quatro genes, ordem ENAM&ndash;AMBN preservada, e controle "
        "negativo mudo.", P))
    h.append(Paragraph(
        "A normalização pelo bitscore do controle positivo é imperfeita: o melhor "
        "HSP de LPL fragmenta de forma diferente entre espécies por causa dos "
        "íntrons, e o frango tem o maior denominador do conjunto. O efeito é "
        "conservador — reduz o índice do frango, que já é o menor. A separação "
        "entre os grupos é grande o bastante para não depender da escolha de "
        "normalização.", P))
    h.append(Paragraph(
        "<i>Melopsittacus undulatus</i> ficou fora do comparativo porque NRG1 não "
        "está anotado nesse assembly. <i>Gavialis</i> e <i>Anolis</i> têm as "
        "âncoras em arranjo que não permite extrair um intervalo orientado; para "
        "essas duas, o alinhamento foi feito sobre a janela ampla que contém as "
        "âncoras, e o índice deve ser lido como limite inferior.", P))

    # ------------------------------------------ 8. referencias
    h.append(Paragraph("8. Referências", H1))
    for r in [
        "Abramyan J, Li M, Khansa R (2025). Reduction of tooth replacement "
        "disproportionately affects the evolution of enamel matrix proteins. "
        "<i>J Mol Evol</i>. doi:10.1007/s00239-025-10258-4",
        "Al-Hashimi N, Sire JY, Delgado S (2010). Evolutionary analysis of "
        "mammalian enamelin, the largest enamel protein. "
        "<i>Mol Biol Evol</i> 27. doi:10.1093/molbev/msq098",
        "Charlesworth D (2020). The timing of genetic degeneration of sex "
        "chromosomes. <i>Philos Trans R Soc B</i> 375. doi:10.1098/rstb.2019.0093",
        "Chikina M, Robinson JD, Clark NL (2016). Hundreds of genes experienced "
        "convergent shifts in selective pressure in marine mammals. "
        "<i>Mol Biol Evol</i> 33. doi:10.1093/molbev/msw112",
        "Crombez E et al. (2025). The subordinate role of pseudogenization to "
        "recombinative deletion following polyploidization in angiosperms. "
        "<i>Nat Commun</i> 16.",
        "Emerling CA et al. (2023). Genomic data suggest parallel dental "
        "vestigialization within the xenarthran radiation. <i>Peer Community J</i>. "
        "doi:10.24072/pcjournal.303",
        "Kowalczyk A et al. (2019). RERconverge: an R package for associating "
        "evolutionary rates with convergent traits. <i>Bioinformatics</i> 35. "
        "doi:10.1093/bioinformatics/btz468",
        "Kowalczyk A et al. (2022). Complementary evolution of coding and "
        "noncoding sequence underlies mammalian hairlessness. <i>eLife</i> 11.",
        "Meredith RW, Gatesy J, Cheng J, Springer MS (2011). Pseudogenization of "
        "the tooth gene enamelysin (MMP20) in the common ancestor of extant baleen "
        "whales. <i>Proc R Soc B</i> 278. doi:10.1098/rspb.2010.1280",
        "Meredith RW, Zhang G, Gilbert MTP, Jarvis ED, Springer MS (2014). "
        "Evidence for a single loss of mineralized teeth in the common avian "
        "ancestor. <i>Science</i> 346:1254390. doi:10.1126/science.1254390",
        "Nixon CF et al. (2024). Exploring the sequence landscape of "
        "ancestral-consensus proteins. <i>Protein Sci</i> 33:e5011. "
        "doi:10.1002/pro.5011",
        "Partha R et al. (2017). Subterranean mammals show convergent regression "
        "in ocular genes and enhancers. <i>eLife</i> 6:e25884.",
        "Partha R, Kowalczyk A, Clark NL, Chikina M (2019). Robust method for "
        "detecting convergent shifts in evolutionary rates. <i>Mol Biol Evol</i> 36. "
        "doi:10.1093/molbev/msz107",
        "Saputra E, Kowalczyk A, Cusick L, Clark NL, Chikina M (2021). "
        "Phylogenetic permulations: a statistically rigorous approach to measure "
        "confidence in associations in a phylogenetic context. "
        "<i>Mol Biol Evol</i> 38. doi:10.1093/molbev/msab068",
        "Sire JY, Delgado SC, Girondot M (2008). Hen's teeth with enamel cap: "
        "from dream to impossibility. <i>BMC Evol Biol</i> 8:246. "
        "doi:10.1186/1471-2148-8-246",
        "Springer MS et al. (2019). Odontogenic ameloblast-associated (ODAM) is "
        "inactivated in toothless/enamelless placental mammals. "
        "<i>BMC Evol Biol</i> 19. doi:10.1186/s12862-019-1359-6",
        "Williams PD, Pollock DD, Blackburne BP, Goldstein RA (2006). Assessing "
        "the accuracy of ancestral protein reconstruction methods. "
        "<i>PLoS Comput Biol</i> 2:e69. doi:10.1371/journal.pcbi.0020069",
        "Zaraisky AG et al. (2024). Loss of the ability to regenerate body "
        "appendages in vertebrates. <i>Biol Rev</i> 99.",
    ]:
        h.append(Paragraph(r, REF))

    h.append(Spacer(1, 14))
    h.append(Paragraph(
        "Este relatório foi gerado por "
        "<font face='Courier'>scripts/relatorio_archosauria.py</font>. O "
        "levantamento completo, com as armadilhas medidas e o programa de "
        "pesquisa, está em "
        "<font face='Courier'>docs/PESQUISA-2026-10-12-FENOTIPO-E-ENSEMBLE.md</font>.",
        NOTA))

    doc.build(h, onFirstPage=rodape, onLaterPages=rodape)
    return destino


if __name__ == "__main__":
    caminho = construir()
    print(f"PDF gerado: {caminho}")
    print(f"tamanho: {caminho.stat().st_size:,} bytes")
