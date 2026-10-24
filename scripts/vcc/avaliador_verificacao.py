# scripts/vcc/avaliador_verificacao.py
"""Junta as corridas de verificacao do `avaliador_proxy` num unico JSON.

Etapa 5 de docs/VCC-2026-ROTA.md. Le os artefatos que `avaliador_corrida.sh`
deixou (um diretorio por corrida) mais o relatorio de
`scripts/vcc/avaliador_peers.py`, e grava `runs/vcc/avaliador/verificacao.json`.

O que entra: as versoes exatas do ambiente, a sequencia de comandos verbatim, o
pico de RSS e o tempo de parede por passo, a pontuacao das seis metricas, os
oito campos de identidade de cada lado, e a extrapolacao de quantas
perturbacoes cabem por corrida. O que NAO entra: nada que precise ser lembrado
em vez de lido — o arquivo e o que responde "o que eu rodei mesmo?".
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

RAIZ = Path(__file__).resolve().parents[2]

#: Memoria da VM do WSL, medida com `free`. O teto de uma corrida e este menos o
#: que o resto do sistema ja ocupa, e nao os 31,4 GiB do host: a VM nao cresce.
WSL_MEM_GIB = 15.0
#: Celulas do painel oficial por contexto: 18 400 de controle + 300 x 400.
CELULAS_OFICIAIS = 18_400 + 300 * 400


def _ler_tempos(caminho: Path) -> dict[str, dict[str, float]]:
    """`tempos.tsv` de `/usr/bin/time -f '%M %e'`, por passo."""
    out: dict[str, dict[str, float]] = {}
    with caminho.open(encoding="utf-8") as fh:
        for linha in csv.DictReader(fh, delimiter="\t"):
            out[linha["passo"]] = {
                "rss_gib": round(float(linha["rss_kb"]) / 1024 / 1024, 3),
                "wall_s": round(float(linha["wall_s"]), 1),
            }
    return out


def _ler_scored(caminho: Path) -> dict[str, dict[str, float | None]]:
    """`scored.csv`: uma linha por metrica mais `avg_score`.

    `from_baseline` NULO em `avg_score` e a marca de submissao MATRICULADA contra
    um bundle de competicao — o `score` anula a media do baseline de proposito,
    para que ninguem leia um numero nao-oficial sob o rotulo que todo mundo le.
    """
    out: dict[str, dict[str, float | None]] = {}
    with caminho.open(encoding="utf-8") as fh:
        for linha in csv.DictReader(fh):
            def num(v: str) -> float | None:
                return round(float(v), 6) if v not in ("", "nan") else None
            out[linha["metric"]] = {
                "from_baseline": num(linha["from_baseline"]),
                "from_replicate": num(linha["from_replicate"]),
            }
    return out


def _versoes(python_venv: Path) -> dict[str, str]:
    codigo = (
        "import importlib.metadata as m, json, sys;"
        "print(json.dumps({**{n: m.version(n) for n in ("
        "'cell-eval2','pdex','numba-mwu','numpy','scipy','polars','pyarrow','numba',"
        "'anndata','scanpy','h5py')}, 'python': sys.version.split()[0]}))"
    )
    saida = subprocess.run([str(python_venv), "-I", "-c", codigo],
                           capture_output=True, text=True, check=True, timeout=120)
    return json.loads(saida.stdout)


def _comandos_verbatim(script: Path) -> list[str]:
    """As invocacoes do `cell-eval2` como estao no script, sem reescrever nada.

    Extrai do proprio arquivo em vez de repetir a lista aqui: duas copias da
    sequencia divergem, e a que o leitor confere e a do script.
    """
    texto = script.read_text(encoding="utf-8")
    bruto = re.findall(r'cronometrar \w+ "\$VENV/bin/cell-eval2"(.*?)(?=\n\n|\n#)', texto,
                       flags=re.S)
    return [" ".join(("cell-eval2 " + t).split()).replace("\\ ", "") for t in bruto]


def _teto_perturbacoes(rss_gib: float, celulas: int, margem_gib: float,
                       base_gib: float) -> dict[str, Any]:
    """Quantas perturbacoes cabem, extrapolando LINEARMENTE no numero de celulas.

    ⚠️ ISTO E EXTRAPOLACAO DE UM PONTO, nao medicao. O unico ponto medido e
    `(celulas, rss_gib)`; tudo que esta acima dele e projecao, e o campo
    `e_extrapolacao` no retorno diz isso para quem ler o JSON sem ler este
    docstring.

    Linear e a forma certa aqui e nao uma simplificacao de conveniencia: o custo
    dominante e a matriz esparsa (nnz proporcional a celulas) mais os pseudobulks
    por perturbacao, ambos lineares. O que NAO e linear e o `pds_cosine`, que e
    O(n_pert^2) em distancias entre perfis — mas sobre 300 perturbacoes isso e uma
    matriz 300x300, 0,7 MB. Nao e o que decide.

    `base_gib` e o custo constante (interprete, numba, llvmlite, polars, pyarrow
    carregados) e vem MEDIDO do passo `score`, que roda o mesmo import e
    praticamente nenhum dado. Subtrai-lo antes de dividir evita que a
    extrapolacao multiplique o custo fixo junto com o variavel.
    """
    por_celula_gib = max(rss_gib - base_gib, 0.0) / celulas
    disponivel = WSL_MEM_GIB - margem_gib - base_gib
    celulas_cabem = int(disponivel / por_celula_gib) if por_celula_gib > 0 else 0
    return {
        "e_extrapolacao": "linear, a partir de UM ponto medido; so o ponto medido e medida",
        "rss_medido_gib": rss_gib,
        "celulas_medidas": celulas,
        "custo_constante_gib": base_gib,
        "custo_constante_de_onde": "pico do passo `score`, que importa o mesmo e nao le dado",
        "margem_reservada_gib": margem_gib,
        "gib_por_mil_celulas": round(por_celula_gib * 1000, 4),
        "celulas_que_cabem": celulas_cabem,
        "perturbacoes_que_cabem_com_18400_controles":
            max((celulas_cabem - 18_400) // 400, 0),
        "celulas_do_painel_oficial": CELULAS_OFICIAIS,
        "rss_projetado_do_painel_inteiro_gib":
            round(base_gib + por_celula_gib * CELULAS_OFICIAIS, 1),
    }


def _comparacao(diretorio: Path, bracos: list[str], margem_gib: float) -> dict[str, Any]:
    """Varios bracos de predicao pontuados contra UM bundle.

    O layout e o de `scratchpad/tres.sh`: `bundle/`, `tempos.tsv`, e um
    `run_<braco>/` + `scored_<braco>.csv` por braco. Um bundle so para todos os
    bracos nao e economia: e o que torna os numeros COMPARAVEIS. As duas pontas da
    escala (baseline e ancora) sao as mesmas, entao a diferenca entre dois bracos
    e diferenca de predicao e nada mais.
    """
    manifesto = json.loads((diretorio / "bundle" / "manifest.json").read_text(encoding="utf-8"))
    tempos = _ler_tempos(diretorio / "tempos.tsv")
    recorte_p = diretorio / "recorte.json"
    saida: dict[str, Any] = {
        "bundle": {
            "real_bundle_id": manifesto["real_bundle_id"],
            "rule_mismatches": manifesto["rule_mismatches"],
            "e_bundle_de_competicao": not manifesto["rule_mismatches"],
            "n_splits": manifesto["n_splits"],
        },
        "tempos_e_memoria_por_passo": tempos,
        "pico_rss_gib": max((v["rss_gib"] for v in tempos.values()), default=0.0),
        "passo_do_pico": max(tempos, key=lambda k: tempos[k]["rss_gib"]) if tempos else None,
        "recorte": (json.loads(recorte_p.read_text(encoding="utf-8"))
                    if recorte_p.exists() else None),
        "bracos": {},
    }
    for braco in bracos:
        meta = json.loads((diretorio / f"run_{braco}" / "run_meta.json")
                          .read_text(encoding="utf-8"))
        divergem = [c for c in ("config_digest", "comparator", "source_fingerprint",
                                "source_fingerprint_strict", "resolved_device",
                                "resolved_de_backend", "input_type_real_effective",
                                "de_real_fingerprint")
                    if manifesto.get(c) != meta.get(c)]
        saida["bracos"][braco] = {
            "pontuacao": _ler_scored(diretorio / f"scored_{braco}.csv"),
            "oito_peers_divergentes": divergem,
            "matriculada": True,
        }
        # `from_baseline` nulo em avg_score E a marca de matricula; ler o CSV e
        # mais honesto que assumir.
        pont = saida["bracos"][braco]["pontuacao"]
        saida["bracos"][braco]["matriculada"] = (
            not divergem and pont.get("avg_score", {}).get("from_baseline") is None)
    if tempos:
        celulas = None
        if saida["recorte"]:
            celulas = saida["recorte"]["celulas_total"]
        if celulas:
            saida["memoria"] = _teto_perturbacoes(
                saida["pico_rss_gib"], celulas, margem_gib,
                tempos.get("score_oraculo", {}).get("rss_gib", 0.0))
    return saida


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__ or "",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--corrida", action="append", required=True, metavar="ROTULO=DIR",
                   help="diretorio de saida de uma corrida, repetivel "
                        "(ex.: k32=/home/manec/vcc-avaliador/out-D-k32-c18400-s0)")
    p.add_argument("--proxy", action="append", required=True, metavar="ROTULO=DIR",
                   help="diretorio do proxy correspondente, repetivel")
    p.add_argument("--peers", required=True, help="JSON de scripts/vcc/avaliador_peers.py")
    p.add_argument("--python-venv", required=True, help="python do env avaliador-proxy")
    p.add_argument("--margem-gib", type=float, default=2.0,
                   help="memoria deixada para o resto do sistema na extrapolacao "
                        "[default: %(default)s]")
    p.add_argument("--comparacao", default=None, metavar="DIR",
                   help="diretorio com varios bracos de predicao pontuados contra UM "
                        "bundle (bundle/, tempos.tsv, run_<braco>/, scored_<braco>.csv)")
    p.add_argument("--bracos", default="oraculo,controle,nossa",
                   help="nomes dos bracos em --comparacao [default: %(default)s]")
    p.add_argument("--saida", default=str(RAIZ / "runs" / "vcc" / "avaliador" /
                                          "verificacao.json"))
    args = p.parse_args()

    proxies = dict(item.split("=", 1) for item in args.proxy)
    corridas: dict[str, Any] = {}
    for item in args.corrida:
        rotulo, diretorio = item.split("=", 1)
        d = Path(diretorio)
        px = Path(proxies[rotulo])
        proxy = json.loads((px / "proxy.json").read_text(encoding="utf-8"))
        manifesto = json.loads((d / "bundle" / "manifest.json").read_text(encoding="utf-8"))
        run_meta = json.loads((d / "run" / "run_meta.json").read_text(encoding="utf-8"))
        tempos = _ler_tempos(d / "tempos.tsv")
        scored = _ler_scored(d / "scored.csv")
        peers = {c: {"bundle": manifesto.get(c), "corrida": run_meta.get(c),
                     "igual": manifesto.get(c) == run_meta.get(c)}
                 for c in ("config_digest", "comparator", "source_fingerprint",
                           "source_fingerprint_strict", "resolved_device",
                           "resolved_de_backend", "input_type_real_effective",
                           "de_real_fingerprint")}
        pico = max((v["rss_gib"] for v in tempos.values()), default=0.0)
        passo_pico = max(tempos, key=lambda k: tempos[k]["rss_gib"]) if tempos else None
        corridas[rotulo] = {
            "proxy": proxy,
            "tempos_e_memoria_por_passo": tempos,
            "pico_rss_gib": pico,
            "passo_do_pico": passo_pico,
            "wall_total_s": round(sum(v["wall_s"] for v in tempos.values()), 1),
            "bundle": {
                "real_bundle_id": manifesto["real_bundle_id"],
                "rule_digest": manifesto["rule_digest"],
                "rule_mismatches": manifesto["rule_mismatches"],
                "e_bundle_de_competicao": not manifesto["rule_mismatches"],
                "members": manifesto["members"],
                "estimators": manifesto["estimators"],
                "bulk_target_sum": manifesto["bulk_target_sum"],
                "n_splits": manifesto["n_splits"],
                "base_seed": manifesto["base_seed"],
            },
            "oito_peers": peers,
            "todos_os_oito_peers_batem": all(v["igual"] for v in peers.values()),
            "de_pred_fingerprint_da_corrida": run_meta.get("de_pred_fingerprint"),
            "pontuacao": scored,
            # O custo constante sai do passo `score` quando ele existe: e o unico
            # passo que importa a pilha inteira e quase nao le dado.
            "memoria": _teto_perturbacoes(
                pico, proxy["real"]["forma"][0], args.margem_gib,
                tempos.get("score", {}).get("rss_gib", 0.0)),
        }

    relatorio = {
        "gerado_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "etapa": "5 — avaliador_proxy",
        "onde": "WSL Ubuntu-24.04, venv em ~/envs/avaliador-proxy/.venv",
        "memoria_da_vm_gib": WSL_MEM_GIB,
        "versoes": _versoes(Path(args.python_venv)),
        "sequencia_verbatim": _comandos_verbatim(
            RAIZ / "scripts" / "vcc" / "avaliador_corrida.sh"),
        "corridas": corridas,
        "campos_de_identidade": json.loads(Path(args.peers).read_text(encoding="utf-8")),
    }
    if args.comparacao:
        comp = _comparacao(
            Path(args.comparacao), [b.strip() for b in args.bracos.split(",") if b.strip()],
            args.margem_gib)
        relatorio["comparacao_de_bracos"] = comp
        # DOIS pontos medidos sao melhores que um, e a inclinacao entre eles e a
        # unica estimativa de custo marginal que nao depende de supor o intercepto.
        pontos = [(c["memoria"]["celulas_medidas"], c["memoria"]["rss_medido_gib"])
                  for c in (*corridas.values(), comp) if "memoria" in c]
        pontos = sorted(set(pontos))
        if len(pontos) >= 2:
            (c0, r0), (c1, r1) = pontos[0], pontos[-1]
            inclin = (r1 - r0) / (c1 - c0)
            intercepto = r0 - inclin * c0
            disponivel = WSL_MEM_GIB - args.margem_gib
            cabem = int((disponivel - intercepto) / inclin) if inclin > 0 else 0
            relatorio["memoria_dois_pontos"] = {
                "pontos_medidos": [{"celulas": c, "pico_rss_gib": r} for c, r in pontos],
                "gib_por_mil_celulas": round(inclin * 1000, 4),
                "intercepto_gib": round(intercepto, 3),
                # Intercepto negativo = o pico cresce MAIS rapido que linear nesta
                # faixa. A reta de dois pontos entao SUBESTIMA o custo acima dela,
                # e o teto verdadeiro e menor que o que ela projeta.
                "crescimento": ("superlinear nesta faixa (intercepto negativo)"
                                if intercepto < 0 else "compativel com linear"),
                "perturbacoes_que_cabem_com_18400_controles":
                    max((cabem - 18_400) // 400, 0),
                "rss_projetado_do_painel_inteiro_gib":
                    round(intercepto + inclin * CELULAS_OFICIAIS, 1),
                "maior_recorte_OBSERVADO_a_completar":
                    {"celulas": c1, "pico_rss_gib": r1,
                     "folga_gib": round(WSL_MEM_GIB - r1, 2)},
                "ponto_de_oom": "NAO MEDIDO — a sondagem de teto foi morta por SIGKILL",
            }
    saida = Path(args.saida)
    saida.parent.mkdir(parents=True, exist_ok=True)
    saida.write_text(json.dumps(relatorio, indent=2, ensure_ascii=False) + "\n",
                     encoding="utf-8")
    print(f"gravado: {saida}")
    for rotulo, c in corridas.items():
        print(f"  {rotulo}: pico {c['pico_rss_gib']} GiB em '{c['passo_do_pico']}', "
              f"{c['wall_total_s']} s, avg_score(replicate)="
              f"{c['pontuacao']['avg_score']['from_replicate']}, "
              f"oito peers batem={c['todos_os_oito_peers_batem']}")
    comp = relatorio.get("comparacao_de_bracos")
    if comp:
        print(f"  comparacao: pico {comp['pico_rss_gib']} GiB em "
              f"'{comp['passo_do_pico']}'")
        for braco, d in comp["bracos"].items():
            avg = d["pontuacao"]["avg_score"]["from_replicate"]
            print(f"    {braco:9s} avg_score={avg}  matriculada={d['matriculada']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
