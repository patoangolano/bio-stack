# /// script
# requires-python = ">=3.11"
# dependencies = ["biopython==1.88", "matplotlib==3.11.2"]
# ///
"""
Regressao da figura de arvore, do involucro e das bordas do utilitario do
pipeline de filogenia.

    uv run --script tests/check_filogenia_figura.py

Sem WSL, sem Docker, poucos segundos. Cobre o que ja deu resposta errada em
silencio nesta stack:

1. Reenraizar desloca o suporte para o ramo vizinho. O Newick guarda o
   suporte no no, e a Biopython so leva o comprimento junto quando reenraiza
   — medido: 3 de 4 biparticoes saiam com o rotulo de outra. O desenho tem
   de ler por biparticao antes de enraizar e reatribuir depois.
2. Rotulo de suporte com outro numero de campos era descartado e o ramo saia
   preto, como forte.
3. O ":nan" do wASTRAL ia parar no suporte do no.
4. O involucro passava caminho do Git Bash, UNC do WSL e "~" de um jeito que
   o WSL nao abria, e o outgroup cru enquanto o pipeline o saneia.
5. O --nome_taxon colchete pegava "[GeneID=...]" no cabecalho do NCBI
   Datasets, e o metodos.md descrevia supermatriz e limpeza que nao rodaram
   quando havia um loco so ou --limpeza nenhuma.
"""

from __future__ import annotations

import csv
import os
import sys
import tempfile
from io import StringIO
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    _reconf = getattr(_s, "reconfigure", None)
    if _reconf is not None:
        _reconf(encoding="utf-8", errors="replace")

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "scripts" / "filogenia"))

# Os dois modulos vem de scripts/filogenia pelo sys.path acima; o analisador
# estatico nao segue o insert, dai o ignore.
import desenhar_arvore as d  # noqa: E402  # pyright: ignore[reportMissingImports]
import filogenia as f  # noqa: E402  # pyright: ignore[reportMissingImports]
from Bio import Phylo  # noqa: E402

sys.path.insert(0, str(RAIZ / "pipelines" / "filogenia" / "bin"))
import filogenia_util as u  # noqa: E402  # pyright: ignore[reportMissingImports]

falhas: list[str] = []


def ok(cond: bool, msg: str) -> None:
    print(f"  {'OK  ' if cond else 'FALHA'} {msg}")
    if not cond:
        falhas.append(msg)


#: Suportes todos fracos e distintos: todo ramo interno ganha rotulo escrito.
NWK = (
    "(A:0.1,B:0.1,((C:0.1,D:0.1)11/0.11/11:0.2,"
    "(E:0.1,(F:0.1,G:0.1)22/0.22/22:0.1)33/0.33/33:0.1)44/0.44/44:0.3);"
)


def biparticoes(arv, rotulo) -> dict[frozenset[str], set[str]]:
    todas = frozenset(t.name for t in arv.get_terminals())
    ref = min(todas)
    out: dict[frozenset[str], set[str]] = {}
    for c in arv.find_clades():
        if c.is_terminal() or c is arv.root:
            continue
        folhas = frozenset(t.name for t in c.get_terminals())
        chave = todas - folhas if ref in folhas else folhas
        if valor := rotulo(c):
            out.setdefault(chave, set()).add(valor)
    return out


def cru(c) -> str | None:
    return c.name or (f"{c.confidence:g}" if c.confidence is not None else None)


def testa_reenraizamento(tmp: Path) -> None:
    print("\n-- suporte segue a biparticao depois de reenraizar --")
    esperado = {
        k: {d._texto(d._campos(next(iter(v)), d.ROTULOS_PADRAO))}
        for k, v in biparticoes(
            Phylo.read(StringIO(NWK), "newick"),  # pyright: ignore[reportPrivateImportUsage]
            cru,
        ).items()
    }
    ingenuo = Phylo.read(StringIO(NWK), "newick")  # pyright: ignore[reportPrivateImportUsage]
    ingenuo.root_with_outgroup("F")
    obtido_ingenuo = {
        k: {d._texto(d._campos(next(iter(v)), d.ROTULOS_PADRAO))}
        for k, v in biparticoes(ingenuo, cru).items()
    }
    errados = sum(1 for k, v in esperado.items() if obtido_ingenuo.get(k) != v)
    ok(errados > 0, f"controle: reenraizar a Biopython desloca rotulos ({errados} de 4)")

    capturado: dict[str, object] = {}
    desenhar_original = Phylo.draw  # pyright: ignore[reportPrivateImportUsage]

    def espiao(arv, **kw):
        capturado["arv"], capturado["bl"] = arv, kw["branch_labels"]
        return desenhar_original(arv, **kw)

    d.Phylo.draw = espiao
    try:
        arq = tmp / "arvore.tre"
        arq.write_text(NWK, encoding="utf-8")
        for og in ("F", "C,D", "A", ""):
            d.desenhar_arvore(arq, tmp / f"fig_{og.replace(',', '') or 'meio'}", outgroup=og)
            obtido = biparticoes(capturado["arv"], capturado["bl"])  # type: ignore[arg-type]
            ok(
                obtido == esperado,
                f"outgroup {og or '(ponto medio)'}: 4 biparticoes com o proprio suporte",
            )
    finally:
        d.Phylo.draw = desenhar_original


def testa_rotulos() -> None:
    print("\n-- leitura do rotulo de suporte --")
    pad = d.ROTULOS_PADRAO
    casos = [
        ("85/0.97/96", {"alrt": 85, "abayes": 0.97, "ufboot": 96}),
        ("85/96", {"alrt": 85, "ufboot": 96}),  # padrao do IQ-TREE sem aBayes
        ("/96", {"ufboot": 96}),  # SH-aLRT vazio: o 96 e o UFBoot
        ("/0.97/96", {"abayes": 0.97, "ufboot": 96}),  # vazio conta na posicao
        ("0.94", {"prob": 0.94}),  # wASTRAL
        ("100", {"ufboot": 100}),
    ]
    for rotulo, esperado in casos:
        obtido = d._campos(rotulo, pad)
        ok(obtido == esperado, f"{rotulo!r} -> {esperado}")
    ok(d._campos("abc/12", pad) is None, "rotulo nao numerico nao e reconhecido")
    ok(d._campos("1/2/3/4", pad) is None, "4 campos nao e reconhecido")
    ok(not d._forte({"alrt": 85, "abayes": 0.90, "ufboot": 99}), "aBayes 0,90 reprova o ramo")
    ok(d._forte({"alrt": 80, "abayes": 0.95, "ufboot": 95}), "80/0,95/95 e forte")


def testa_nan(tmp: Path) -> None:
    print("\n-- arvore do wASTRAL com :nan --")
    arq = tmp / "especies.tre"
    arq.write_text("((A,(B,C)0.94:nan)0.99:1.2,D);", encoding="utf-8")
    saidas = d.desenhar_arvore(arq, tmp / "especies", tipo="especies")
    ok(all(p.is_file() and p.stat().st_size > 0 for p in saidas), "figura sai (cladograma)")


def testa_involucro(tmp: Path) -> None:
    print("\n-- involucro: caminhos, aspas, CSV, outgroup --")
    casos = {
        "\\\\wsl.localhost\\Ubuntu-24.04\\home\\manec\\x": "/home/manec/x",
        "\\\\wsl$\\Ubuntu-24.04\\home": "/home",
        "/c/Users/Manec/x": "/mnt/c/Users/Manec/x",
        "/mnt/c/x": "/mnt/c/x",
        "/home/manec/x": "/home/manec/x",
        "~/dados": "~/dados",
    }
    if os.name == "nt":
        casos["C:\\Users\\Manec\\bio-stack"] = "/mnt/c/Users/Manec/bio-stack"
    for entrada, esperado in casos.items():
        ok(f.para_wsl(entrada) == esperado, f"para_wsl({entrada!r}) = {esperado!r}")
    try:
        f.para_wsl("\\\\wsl.localhost\\Debian\\home")
        ok(False, "caminho de outra distro e recusado")
    except SystemExit:
        ok(True, "caminho de outra distro e recusado")
    ok(f._aspas("~/a b") == "\"$HOME\"'/a b'", 'til fica fora das aspas como "$HOME"')
    ok(f._aspas("it's") == "'it'\"'\"'s'", "aspas simples escapadas")

    (tmp / "sub").mkdir()
    csv_in = tmp / "locos.csv"
    # BOM do PowerShell 5.1 e caminho relativo ao CSV, nao ao diretorio corrente.
    csv_in.write_text(
        "\ufeffloco,fasta\nA,sub/a.fa\nB,/home/manec/b.fa\nC,~/c.fa\n", encoding="utf-8"
    )
    saida = tmp / "run"
    saida.mkdir()
    destino = f.preparar_entrada(str(csv_in), saida)
    with (saida / "entrada.csv").open(encoding="utf-8", newline="") as fh:
        linhas = {r["loco"]: r["fasta"] for r in csv.DictReader(fh)}
    ok(destino.endswith("/entrada.csv"), "CSV reescrito em runs/filogenia/<id>/entrada.csv")
    ok(linhas.get("A") == f.para_wsl(tmp / "sub" / "a.fa"), "relativo resolvido contra o CSV")
    ok(
        linhas.get("B") == "/home/manec/b.fa" and linhas.get("C") == "~/c.fa",
        "caminho do WSL intacto",
    )
    ok(
        f.saneia_outgroup("Homo sapiens, Mus musculus") == "Homo_sapiens,Mus_musculus",
        "outgroup saneado com a regra do pipeline",
    )


def testa_utilitario() -> None:
    print("\n-- utilitario do pipeline: nome do taxon e Metodos nas bordas --")
    cabecalhos = {
        "XP_015134201.2 enamelin [Gallus gallus]": "Gallus_gallus",
        "NP_114095.2 ENAM [organism=Homo sapiens] [GeneID=10117] [isoform=1]": "Homo_sapiens",
        "XP_1 ENAM [organism=Gallus gallus] [GeneID=1] [isoform=X1]": "Gallus_gallus",
        "XP_2 descricao [GeneID=3]": "",  # so chave que nao e organismo: vira erro
    }
    for cab, esperado in cabecalhos.items():
        ok(u.sanear_nome(cab, "colchete") == esperado, f"colchete: {cab[:44]!r} -> {esperado!r}")

    params = {
        "alinhador": "mafft",
        "limpeza": "nenhuma",
        "modelo": "MFP",
        "alrt": 1000,
        "abayes": True,
        "ufboot": 1000,
        "bnni": True,
        "seed": 1,
    }
    versoes = {"mafft": "7.520", "iqtree": "3.1.4", "nextflow": "26.04.6"}
    um_loco = u.montar_metodos(params, versoes, {"n_locos_inferidos": 1})
    ok("Para o único loco" in um_loco, "um loco so: 'Para o único loco'")
    ok("supermatriz" not in um_loco, "um loco so: nenhuma supermatriz no texto")
    concat = {"n_particoes_entrada": 2, "n_particoes_finais": 1, "n_sitios": 900, "n_taxa": 6}
    sem_limpeza = u.montar_metodos(params, versoes, {"n_locos_inferidos": 2, "concatenado": concat})
    ok(
        "alinhamentos limpos" not in sem_limpeza
        and "Os alinhamentos foram concatenados" in sem_limpeza,
        "--limpeza nenhuma: alinhamentos sem 'limpos'",
    )
    com_aviso = u.montar_metodos(
        params, versoes, {"n_locos_inferidos": 1, "avisos_preparo": {"loco_X": ["a: '*' interno"]}}
    )
    ok(
        "avisos_preparo" in com_aviso.split("## Pendências", 1)[-1],
        "aviso do preparo vira pendencia",
    )
    # Rede de seguranca: nome de taxon errado (acesso no lugar da especie) roda
    # com exit 0; o que denuncia e a ocupacao, a falta de taxon em comum e o
    # outgroup sumido — os tres tem de virar pendencia.
    ruim = {**concat, "ocupacao_perc": 25.0, "n_taxa_multiloco": 0}
    pend = u.montar_metodos(
        {**params, "outgroup": "Homo_sapiens"},
        versoes,
        {"n_locos_inferidos": 4, "concatenado": ruim, "outgroup_ausente": True},
    ).split("## Pendências", 1)[-1]
    ok("Ocupação táxon × loco de 25.0%" in pend, "ocupacao < 50% vira pendencia")
    ok("Nenhum táxon aparece em mais de um loco" in pend, "nenhum taxon em comum vira pendencia")
    ok("não aparece em nenhum loco" in pend, "outgroup ausente vira pendencia")


def main() -> int:
    print("== figura de arvore e involucro do pipeline de filogenia ==")
    with tempfile.TemporaryDirectory() as tmp:
        testa_reenraizamento(Path(tmp))
        testa_rotulos()
        testa_nan(Path(tmp))
        testa_involucro(Path(tmp))
    testa_utilitario()
    print()
    if falhas:
        print(f"FALHOU: {len(falhas)} verificacao(oes).")
        return 1
    print("Figura, involucro e utilitario OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
