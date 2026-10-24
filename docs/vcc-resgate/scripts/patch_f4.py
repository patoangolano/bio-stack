"""Separa o agrupamento da resposta media do gate do desvio, e usa peso igual.

Medido em `runs/vcc/transfer/pool_media.json`, deixar-um-contexto-de-fora sobre
8 contextos doadores, metrica identica a do `resposta_media.py`:

    uniforme       mediana +0,2176   <- vencedor
    amostra        mediana +0,2025
    concordancia   mediana +0,1685
    profundidade   mediana +0,1660   (so tirar o gate)
    atual          mediana +0,1002

Duas mudancas, e a ordem de tamanho importa: tirar o gate vale +0,066, trocar
profundidade por peso igual vale outros +0,052.
"""

from pathlib import Path

CONSTANTE_ANTES = "# Programas de resposta retidos na base de baixa dimensao."

CONSTANTE_DEPOIS = '''
# Minimo de perturbacoes para um contexto entrar no AGRUPAMENTO DA MEDIA. A
# media de um contexto com poucas perturbacoes e ruido: no limite de uma
# perturbacao ela E aquela perturbacao. Os 30 contextos do Jiang tem de 1 a 39
# perturbacoes e ficam todos fora; os 8 restantes tem de 44 a 8.987.
#
# O valor cai UMA perturbacao acima do maior contexto do Jiang (39), o que e
# perto demais para passar sem teste. Medido em pool_media_min10.json (30
# contextos) e pool_media_min200.json (6): `atual` perde nos tres cortes e
# `concordancia` nao vence em nenhum. Em 200 `profundidade` vence, mas esse
# corte exclui o arc_vcc_25 — a propria patologia que torna profundidade ruim —
# entao nao e evidencia a favor dela.
MIN_PERTS_MEDIA = 40
# Programas de resposta retidos na base de baixa dimensao.
'''.strip("\n")

LACO_ANTES = '''
    procedencia: dict[str, list[str]] = {}
    for repo in REPOS:
        meta = json.loads((METAS / REPOS[repo] / "meta.json").read_text(encoding="utf-8"))
        dataset = str(meta["dataset"])
        catalogo = perts_por_contexto(repo)
        genes_repo = [g for g in genes_do_repo(repo) if g in indice_gene]
        for contexto, disponiveis in sorted(catalogo.items()):
'''.strip("\n")

LACO_DEPOIS = '''
    procedencia: dict[str, list[str]] = {}

    # ========================================================================
    # AGRUPAMENTO DA RESPOSTA MEDIA — laco PROPRIO, e isso e a correcao.
    #
    # Antes a media era acumulada dentro do laco do desvio, atras do gate
    # `set(alvos) & disponiveis & split; if not usaveis: continue`. O gate e
    # correto para o DESVIO, que precisa do nosso alvo medido naquele contexto.
    # A media nao precisa de alvo nosso nenhum: `media_do_contexto()` percorre
    # TODAS as perturbacoes do doador. Um gate servindo duas parcelas com
    # requisitos diferentes excluia os 4 contextos do Replogle — os de melhor
    # transferencia medida de media do conjunto inteiro — porque zero dos
    # nossos 300 alvos esta no split de treino do Replogle.
    #
    # E o peso e IGUAL, nao profundidade. `profundidade` punha 93,2% do
    # agrupamento no arc_vcc_25.ARC_H1, cuja media e ORTOGONAL a do Orion
    # (r = -0,0001 em 8.800 genes de magnitude comparavel, logo nao e falta de
    # sobreposicao). Com estimadores ruidosos e quase intercambiaveis, peso
    # igual compra a reducao de variancia que qualquer ponderacao confiante
    # joga fora — inclusive ponderar pela concordancia medida, que perde
    # tambem (+0,1685 contra +0,2176).
    # ========================================================================
    medias_de_contexto: dict[tuple[str, str], np.ndarray] = {}
    contextos_na_media: list[str] = []
    for repo in REPOS:
        catalogo = perts_por_contexto(repo)
        genes_repo = [g for g in genes_do_repo(repo) if g in indice_gene]
        colunas_media = np.array([indice_gene[g] for g in genes_repo], dtype=np.int64)
        for contexto, disponiveis in sorted(catalogo.items()):
            if len(disponiveis) < MIN_PERTS_MEDIA:
                continue
            media = media_do_contexto(repo, contexto, genes_repo)
            medias_de_contexto[(repo, contexto)] = media
            contextos_na_media.append(f"{repo}.{contexto}")
            tem_media = media != 0.0
            media_numerador[colunas_media] += media * tem_media
            media_denominador[colunas_media] += tem_media
    print(
        f"  resposta media agrupada com peso IGUAL sobre {len(contextos_na_media)} "
        f"contextos (>= {MIN_PERTS_MEDIA} perts): {', '.join(contextos_na_media)}"
    )

    for repo in REPOS:
        meta = json.loads((METAS / REPOS[repo] / "meta.json").read_text(encoding="utf-8"))
        dataset = str(meta["dataset"])
        catalogo = perts_por_contexto(repo)
        genes_repo = [g for g in genes_do_repo(repo) if g in indice_gene]
        for contexto, disponiveis in sorted(catalogo.items()):
'''.strip("\n")

DESVIO_ANTES = '''
            print(f"  {repo}.{contexto}: {len(usaveis)} dos 300 alvos no split de treino")
            media = media_do_contexto(repo, contexto, genes_repo)
            resposta = construir(repo, contexto, usaveis, genes_repo, media_externa=media)
            # A media deste contexto entra no agrupamento ponderada pela
            # profundidade do ensaio: o arc_vcc_25 (2.016 genes significativos
            # por perturbacao) domina, e e o regime mais parecido com o nosso.
            colunas_media = np.array([indice_gene[g] for g in genes_repo], dtype=np.int64)
            tem_media = media != 0.0
            media_numerador[colunas_media] += profundidade[repo] * media * tem_media
            media_denominador[colunas_media] += profundidade[repo] * tem_media
            beta = beta_por_doador.get(repo, beta_mediano)
'''.strip("\n")

DESVIO_DEPOIS = '''
            print(f"  {repo}.{contexto}: {len(usaveis)} dos 300 alvos no split de treino")
            # Reusa a media ja calculada no laco de cima quando o contexto entrou
            # la: `orion.hct116` percorre 8.987 perturbacoes e ler duas vezes e
            # desperdicio puro. Contexto abaixo de MIN_PERTS_MEDIA nao esta no
            # cache e ainda precisa da media PROPRIA para centrar o desvio — que
            # e uso diferente de entrar no agrupamento.
            media = medias_de_contexto.get((repo, contexto))
            if media is None:
                media = media_do_contexto(repo, contexto, genes_repo)
            resposta = construir(repo, contexto, usaveis, genes_repo, media_externa=media)
            beta = beta_por_doador.get(repo, beta_mediano)
'''.strip("\n")


def main() -> None:
    alvo = Path("scripts/vcc/transfer_fit.py")
    texto = alvo.read_text(encoding="utf-8")
    trocas = (
        (CONSTANTE_ANTES, CONSTANTE_DEPOIS),
        (LACO_ANTES, LACO_DEPOIS),
        (DESVIO_ANTES, DESVIO_DEPOIS),
    )
    for antes, depois in trocas:
        if texto.count(antes) != 1:
            raise SystemExit(f"ancora nao unica ({texto.count(antes)}x): {antes[:70]!r}")
        texto = texto.replace(antes, depois)
    alvo.write_text(texto, encoding="utf-8")
    print(f"{len(trocas)} trocas aplicadas em {alvo}")


if __name__ == "__main__":
    main()
