# scripts/vcc/avaliador_peers.py
"""Mede o que faz divergir cada campo de identidade do `score --real-bundle`.

Etapa 5 de docs/VCC-2026-ROTA.md. O `score --real-bundle` compara oito campos
entre o `manifest.json` do bundle e o `run_meta.json` da corrida, e qualquer
divergencia e FATAL — `real_bundle.check_submission` levanta `ValueError` e
`--allow-config-mismatch` ao lado de `--real-bundle` e erro de uso, nao waiver.
Como o placar fica oculto ate 2026-11-05 e aceita 2 submissoes por dia,
descobrir na hora qual campo divergiu custa meio dia de cota.

⚠️ E HA MAIS QUE OITO. `score.expect_from_run_meta` acrescenta tres portoes que
nao estao na lista dos oito e que tambem abortam:

  * `source_fingerprint_strict` tem de ser **True**, nao apenas igual ao do
    bundle. Dois artefatos ambos nao-estritos concordariam entre si e seriam
    recusados de todo modo, porque o portao da ancora e o hash de CONTEUDO.
  * `anchor_semantic_identity` — hash sobre `anchor._SEMANTIC_FIELDS`, que inclui
    `allow_fractional_counts`. Esse campo esta em `DIGEST_EXEMPT_FIELDS`, logo
    **nao move o `config_digest`**: e um campo invisivel aos oito e fatal aqui.
  * `anchor_metric_names` — o conjunto de metricas que a ancora cobre.
  * e `de_pred_fingerprint`, checado de UM LADO SO, contra a constante `None`:
    uma tabela de DE do lado predito fornecida com `run --de-pred` recusa a
    submissao, mesmo com os oito batendo.

Em vez de rodar a sequencia inteira por variante (meia hora cada), este script
chama `baseline.build_run_meta` direto — a MESMA funcao que o `run` usa para
gravar o `run_meta.json`, e por construcao: o docstring dela diz que viver no
`baseline.py` em vez do `cli.py` e o que impede os dois registros de divergirem.
Depois passa o resultado por `real_bundle.check_submission` e por
`score.expect_from_run_meta`, as MESMAS funcoes que o `score` chama. Entao o
veredito aqui e o veredito de la, sem a conta de CPU.

`build_run_meta` le os dois lados BACKED e nunca materializa X — exceto pelo
`source_fingerprint` sob `cache_strict`, que e hash de conteudo e portanto le a
matriz. O preset `vcc2026` liga `cache_strict`, entao cada variante custa uma
leitura do lado real.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

from cell_eval2.baseline import build_run_meta
from cell_eval2.config import EvalConfig
from cell_eval2.real_bundle import SUBMISSION_PEERS, check_submission, read_real_bundle
from cell_eval2.score import expect_from_run_meta

#: Como a corrida de referencia e montada: o preset, a coluna do painel e o
#: backend fixado. Qualquer desvio disto e uma das variantes abaixo.
PERT_COL = "target_gene"
BACKEND = "pdex"

#: Portoes do caminho `--real-bundle` que NAO estao nos oito peers. Medidos, nao
#: inferidos: ver o docstring do modulo.
PORTOES_EXTRA = ("anchor_semantic_identity", "anchor_metric_names", "de_pred_fingerprint")


def _config_base() -> EvalConfig:
    cfg = EvalConfig.from_preset("vcc2026")
    return replace(cfg, pert_col=PERT_COL, de=replace(cfg.de, backend=BACKEND))


def _variantes(cfg: EvalConfig, real: str, pred: str,
               real_alt: str | None, de_real: str | None,
               ) -> list[tuple[str, str, EvalConfig, str, str, str | None, str | None]]:
    """(nome, mudanca, config, real, pred, de_real, de_pred). A primeira e a correta."""
    v: list[tuple[str, str, EvalConfig, str, str, str | None, str | None]] = [
        ("controle", "nada — e a corrida que o bundle espera",
         cfg, real, pred, None, None),
        ("pert_col", "--pert-col target (o default do preset) em vez de target_gene",
         replace(cfg, pert_col="target"), real, pred, None, None),
        ("de_backend_scanpy", "--set de.backend=scanpy em vez de pdex",
         replace(cfg, de=replace(cfg.de, backend="scanpy")), real, pred, None, None),
        ("de_backend_auto", "--set de.backend=auto (o default do preset, resolvido no host)",
         replace(cfg, de=replace(cfg.de, backend="auto")), real, pred, None, None),
        ("p_adj_threshold", "--set de.p_adj_threshold=0.01 em vez de 0.05",
         replace(cfg, de=replace(cfg.de, p_adj_threshold=0.01)), real, pred, None, None),
        ("min_abs_log2fc", "--set de.min_abs_log2fc=0.25 em vez de 0.0",
         replace(cfg, de=replace(cfg.de, min_abs_log2fc=0.25)), real, pred, None, None),
        ("fdr_scope", "--set de.fdr_scope=global em vez de per_pert",
         replace(cfg, de=replace(cfg.de, fdr_scope="global")), real, pred, None, None),
        ("filtro_cpm", "--set filter.filter_gene_min_cpm_cell=1.0 em vez de 5.0",
         replace(cfg, filter=replace(cfg.filter, filter_gene_min_cpm_cell=1.0)),
         real, pred, None, None),
        ("bulk_target_sum", "--set bulk_target_sum=1000000.0 em vez de 50000.0",
         replace(cfg, bulk_target_sum=1.0e6), real, pred, None, None),
        ("exclude_target_gene", "--set discrimination.exclude_target_gene=false",
         replace(cfg, discrimination=replace(cfg.discrimination, exclude_target_gene=False)),
         real, pred, None, None),
        ("tie_policy", "--set discrimination.tie_policy=position em vez de midrank",
         replace(cfg, discrimination=replace(cfg.discrimination, tie_policy="position")),
         real, pred, None, None),
        ("profile", "--profile full em vez do perfil vcc2026",
         replace(cfg, metrics="full"), real, pred, None, None),
        ("input_type_lognorm", "--input-type lognorm em vez de counts",
         replace(cfg, input_type="lognorm"), real, pred, None, None),
        ("sem_cache_strict",
         "--no-cache-strict (impressao digital de METADADO em vez de conteudo)",
         replace(cfg, cache_strict=False), real, pred, None, None),
        # ⚠️ O caso que os oito peers NAO pegam: `allow_fractional_counts` esta em
        # `DIGEST_EXEMPT_FIELDS` e nao move o config_digest, mas entra em
        # `anchor._SEMANTIC_FIELDS` e portanto move a identidade semantica da ancora.
        ("allow_fractional_counts",
         "--set allow_fractional_counts=true (ISENTO do config_digest)",
         replace(cfg, allow_fractional_counts=True), real, pred, None, None),
        # Vizinho do anterior em `_SEMANTIC_FIELDS` e NAO isento do digest: medido,
        # move os dois. E o contraste que mostra que `allow_fractional_counts` e
        # excecao, nao regra.
        ("max_counts_per_cell",
         "--set max_counts_per_cell=500000.0 em vez de 1e6",
         replace(cfg, max_counts_per_cell=5.0e5), real, pred, None, None),
        ("device_cpu_explicito", "--set device=cpu em vez de auto",
         replace(cfg, device="cpu"), real, pred, None, None),
        ("num_threads", "--set num_threads=4 (campo isento do digest)",
         replace(cfg, num_threads=4), real, pred, None, None),
        ("cache_real", "--cache-real <dir> (campo isento do digest)",
         replace(cfg, cache_real="/tmp/cache-inexistente"), real, pred, None, None),
    ]
    if real_alt is not None:
        # Mesma forma, mesmos rotulos, mesmos genes — so as CONTAGENS diferem.
        # Sob `cache_strict` isso move `source_fingerprint`; sem ele, nada moveria.
        v.append(("real_outro_conteudo",
                  "mesmo lado real em forma e rotulos, contagens diferentes (outra semente)",
                  cfg, real_alt, pred, None, None))
    if de_real is not None:
        v.append(("de_real_fornecida",
                  "run --de-real <parquet>: tabela de DE do lado real FORNECIDA",
                  cfg, real, pred, de_real, None))
        v.append(("de_pred_fornecida",
                  "run --de-pred <parquet>: tabela de DE do lado predito FORNECIDA",
                  cfg, real, pred, None, de_real))
    return v


def _diferencas(manifesto: dict[str, Any], meta: dict[str, Any],
                campos: tuple[str, ...]) -> list[str]:
    ausente = object()
    return [c for c in campos if manifesto.get(c, ausente) != meta.get(c, ausente)]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__ or "",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bundle", required=True, help="diretorio do bundle (prep-real-bundle)")
    p.add_argument("--real", required=True, help="real.h5ad da corrida")
    p.add_argument("--pred", required=True, help="pred.h5ad da corrida")
    p.add_argument("--real-alternativo", default=None,
                   help="segundo real.h5ad, mesma forma e outras contagens, para mover "
                        "source_fingerprint")
    p.add_argument("--de-real", default=None,
                   help="de_real.parquet de `run --write-degenes`, para mover "
                        "de_real_fingerprint e de_pred_fingerprint")
    p.add_argument("--saida", default=None, help="grava o relatorio como JSON aqui")
    args = p.parse_args()

    bundle = read_real_bundle(args.bundle)
    manifesto = bundle.manifest
    # A ancora do bundle: o `anchor_meta.json` e quem carrega a identidade
    # semantica, nao o manifesto (que guarda a copia sob o mesmo nome).
    espera_bundle = {
        "anchor_semantic_identity": manifesto.get("anchor_semantic_identity"),
        "anchor_metric_names": manifesto.get("anchor_metric_names"),
        "de_pred_fingerprint": None,  # nulo por construcao no lado do bundle
    }

    linhas: list[dict[str, Any]] = []
    for nome, mudanca, cfg, real, pred, de_real, de_pred in _variantes(
            _config_base(), args.real, args.pred, args.real_alternativo, args.de_real):
        try:
            meta = build_run_meta(cfg, real, pred, de_real=de_real, de_pred=de_pred)
        except Exception as exc:  # noqa: BLE001 — a mensagem exata E o resultado
            # Uma variante pode nem chegar a gravar run_meta: --pert-col errado
            # morre na leitura do obs, antes de qualquer campo existir. Falhar
            # ANTES e melhor que divergir depois, e vale registrar qual faz isso.
            linha = {"variante": nome, "mudanca": mudanca, "divergem_dos_oito": None,
                     "aborta_antes_do_run_meta": f"{type(exc).__name__}: {exc}"
                                                 .split("\n")[0][:300]}
            linhas.append(linha)
            print(f"  {nome:24s} ABORTA antes do run_meta: {type(exc).__name__}")
            continue

        oito = _diferencas(manifesto, meta, SUBMISSION_PEERS)
        extra = _diferencas(espera_bundle, meta, PORTOES_EXTRA)
        try:
            check_submission(manifesto, meta)
            veredito_oito = "aceita"
        except ValueError:
            veredito_oito = "recusada"
        try:
            expect_from_run_meta(meta)
            veredito_ancora = "ok"
        except ValueError as exc:
            veredito_ancora = f"recusa: {str(exc).split(chr(10))[0][:120]}"

        linhas.append({
            "variante": nome,
            "mudanca": mudanca,
            "divergem_dos_oito": oito,
            "divergem_dos_portoes_extra": extra,
            "veredito_check_submission": veredito_oito,
            "veredito_expect_from_run_meta": veredito_ancora,
            "valores": {c: meta.get(c) for c in (*oito, *extra)},
        })
        marca = "+extra" if extra else ""
        print(f"  {nome:24s} {veredito_oito:9s} oito={oito or '—'} {marca}")

    oito_movidos = sorted({c for ln in linhas if ln.get("divergem_dos_oito")
                           for c in ln["divergem_dos_oito"]})
    nao_movidos = [c for c in SUBMISSION_PEERS if c not in oito_movidos]
    relatorio = {
        "bundle": Path(args.bundle).name,
        "peers": list(SUBMISSION_PEERS),
        "portoes_extra": list(PORTOES_EXTRA),
        "variantes": linhas,
        "dos_oito_movidos_por_alguma_variante": oito_movidos,
        "dos_oito_que_nenhuma_variante_moveu": nao_movidos,
    }
    print(f"\n  dos oito, movidos aqui: {oito_movidos}")
    print(f"  dos oito, nao movidos aqui: {nao_movidos}")
    if args.saida:
        Path(args.saida).parent.mkdir(parents=True, exist_ok=True)
        Path(args.saida).write_text(json.dumps(relatorio, indent=2, ensure_ascii=False) + "\n",
                                    encoding="utf-8")
        print(f"  gravado: {args.saida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
