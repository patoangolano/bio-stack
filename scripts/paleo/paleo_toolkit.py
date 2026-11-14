# /// script
# requires-python = ">=3.11"
# dependencies = ["httpx>=0.27"]
# ///
"""
paleo_toolkit.py — Utilitario de precisao para projetos de Paleontologia na bio-stack.

Tres subcomandos auditaveis (todos registram proveniencia em runs/provenance.jsonl):

  1. pbdb-range    Consulta a API v1.2 da Paleobiology Database (vocab=pbdb),
                   extrai FAD/LAD (First/Last Appearance Datum em Ma) e ocorrencias,
                   e gera blocos de calibracao prontos para RevBayes (FBD),
                   PAML MCMCTree e IQ-TREE 3 (--date / LSD2).
  2. zooms-qc      Digestao in silico de paleoproteinas (esmalte e colageno:
                   AMELX, ENAM, AMBN, MMP20, COL1A1, COL1A2) com calculo exato
                   de massa monoisotopica e deslocamentos diageneticos de
                   desamidacao (Asn->Asp +0.984016 Da, Gln->Glu +0.984016 Da) e
                   hidroxilacao de prolina (+15.994915 Da).
  3. mk-audit      Audita matrizes morfologicas discretas (NEXUS ou FASTA de
                   estados 0-9/?/-) para vies de aquisicao de Lewis (2001)
                   (Mk vs Mkv/+ASC), contando sitios invariaveis, autapomorfias,
                   fracao de dados ausentes por taxon fossil e emitindo a
                   receita exata para RevBayes, MrBayes e IQ-TREE 3.
  4. altall-asr    Gera o par de sequencias ancestrais ML (estado de maior PP) e
                   AltAll (Eick et al. 2017: troca todo sitio com PP < limiar,
                   padrao 0.80, pelo 2o estado mais provavel) a partir de um
                   arquivo .state do IQ-TREE 3 (--ancestral), eliminando o vies
                   de consenso antes de dobrar no Boltz-2 / OpenFold3.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx

RAIZ = Path(__file__).resolve().parents[2]
PROV_LOG = Path(os.environ.get("BIO_STACK_PROVENANCE", RAIZ / "runs" / "provenance.jsonl"))
PBDB_BASE = "https://paleobiodb.org/data1.2"

# Massas monoisotopicas IUPAC/Unimod (Da) para paleoproteomica de alta precisao
MASSA_AA: dict[str, float] = {
    "A": 71.037114,
    "R": 156.101111,
    "N": 114.042927,
    "D": 115.026943,
    "C": 103.009185,
    "E": 129.042593,
    "Q": 128.058578,
    "G": 57.021464,
    "H": 137.058912,
    "I": 113.084064,
    "L": 113.084064,
    "K": 128.094963,
    "M": 131.040485,
    "F": 147.068414,
    "P": 97.052764,
    "S": 87.032028,
    "T": 101.047679,
    "W": 186.079313,
    "Y": 163.063320,
    "V": 99.068414,
}
MASSA_H2O = 18.010565
MASSA_PROTON = 1.00727647
DELTA_DESAMIDACAO = 0.984016  # N->D ou Q->E (NH2 -> OH)
DELTA_HIDROXIPROLINA = 15.994915  # P -> Hyp (colageno)
DELTA_OXIDACAO_M = 15.994915  # M -> Met(O)


def registrar_proveniencia(acao: str, params: dict[str, Any], resumo: dict[str, Any]) -> None:
    try:
        PROV_LOG.parent.mkdir(parents=True, exist_ok=True)
        entry = {
            "ts": datetime.now(UTC).isoformat(),
            "tool": f"paleo_toolkit:{acao}",
            "params": params,
            "summary": resumo,
        }
        with PROV_LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except OSError:
        pass


def cmd_pbdb_range(args: argparse.Namespace) -> int:
    nomes = [n.strip() for n in args.taxa.split(",") if n.strip()]
    resultados: list[dict[str, Any]] = []
    headers = {"Accept": "application/json", "User-Agent": "bio-stack-paleo-toolkit/1.0"}

    with httpx.Client(timeout=20.0, headers=headers) as client:
        for nome in nomes:
            r = client.get(
                f"{PBDB_BASE}/taxa/single.json",
                params={"name": nome, "vocab": "pbdb", "show": "app,attr,parent,ecotaph"},
            )
            r.raise_for_status()
            recs = r.json().get("records", [])
            if not recs:
                resultados.append({"query": nome, "found": False})
                continue
            rec = recs[0]
            fad_max = float(rec.get("firstapp_max_ma") or 0.0)
            fad_min = float(rec.get("firstapp_min_ma") or 0.0)
            lad_max = float(rec.get("lastapp_max_ma") or 0.0)
            lad_min = float(rec.get("lastapp_min_ma") or 0.0)
            slug = re.sub(r"[^A-Za-z0-9_]", "_", rec.get("taxon_name", nome))
            resultados.append(
                {
                    "query": nome,
                    "found": True,
                    "taxon_no": rec.get("orig_no") or rec.get("taxon_no"),
                    "taxon_name": rec.get("taxon_name"),
                    "rank": rec.get("taxon_rank"),
                    "is_extant": rec.get("is_extant"),
                    "n_occs": rec.get("n_occs", 0),
                    "fad_ma": [fad_min, fad_max],
                    "lad_ma": [lad_min, lad_max],
                    "environment": rec.get("environment"),
                    "diet": rec.get("diet"),
                    "calibrations": {
                        "revbayes_fbd_tip": f"age_{slug} ~ dnUniform({lad_min:.2f}, {fad_max:.2f})",
                        "mcmctree_soft_min": f"'>{fad_min / 100.0:.4f}'  # min {fad_min:.2f} Ma (escala 100 Ma)",
                        "iqtree_lsd2_date": f"{slug}\t-{fad_max:.2f}:-{lad_min:.2f}"
                        if rec.get("is_extant") == "extinct"
                        else f"{slug}\t0",
                    },
                }
            )

    registrar_proveniencia(
        "pbdb-range",
        {"taxa": nomes, "endpoint": f"{PBDB_BASE}/taxa/single.json"},
        {"n_queried": len(nomes), "n_found": sum(1 for x in resultados if x.get("found"))},
    )
    print(json.dumps({"access_date": datetime.now(UTC).isoformat(), "taxa": resultados}, indent=2, ensure_ascii=False))
    return 0


def digerir_tripsina(seq: str, min_len: int = 6, max_len: int = 35) -> list[str]:
    """Cliva apos K ou R, exceto antes de P (regra de Keil)."""
    s = re.sub(r"\s+", "", seq).upper()
    cortes = [0]
    for i in range(len(s) - 1):
        if s[i] in ("K", "R") and s[i + 1] != "P":
            cortes.append(i + 1)
    cortes.append(len(s))
    peps = [s[cortes[j] : cortes[j + 1]] for j in range(len(cortes) - 1)]
    return [p for p in peps if min_len <= len(p) <= max_len and set(p) <= set(MASSA_AA)]


def massa_peptideo(pep: str, n_desam: int = 0, n_hyp: int = 0) -> float:
    base = sum(MASSA_AA[aa] for aa in pep) + MASSA_H2O
    return base + n_desam * DELTA_DESAMIDACAO + n_hyp * DELTA_HIDROXIPROLINA


def cmd_zooms_qc(args: argparse.Namespace) -> int:
    seq = args.sequence.strip().upper()
    peps = digerir_tripsina(seq, min_len=args.min_len, max_len=args.max_len)
    painel: list[dict[str, Any]] = []
    for p in peps:
        n_nq = p.count("N") + p.count("Q")
        n_p = p.count("P")
        m0 = massa_peptideo(p, 0, 0)
        m_deam_all = massa_peptideo(p, n_nq, 0)
        painel.append(
            {
                "peptide": p,
                "length": len(p),
                "deamidation_sites_NQ": n_nq,
                "proline_sites_P": n_p,
                "neutral_mass_da": round(m0, 5),
                "mh_plus_1": round(m0 + MASSA_PROTON, 5),
                "mh2_plus_2": round((m0 + 2 * MASSA_PROTON) / 2.0, 5),
                "mh_plus_1_full_deamidated": round(m_deam_all + MASSA_PROTON, 5),
                "delta_full_deamidation_da": round(n_nq * DELTA_DESAMIDACAO, 5),
            }
        )
    registrar_proveniencia(
        "zooms-qc",
        {"seq_len": len(seq), "min_len": args.min_len, "max_len": args.max_len},
        {"n_peptides": len(painel), "total_nq_sites": sum(x["deamidation_sites_NQ"] for x in painel)},
    )
    print(json.dumps({"n_tryptic_peptides": len(painel), "peptides": painel}, indent=2, ensure_ascii=False))
    return 0


def cmd_mk_audit(args: argparse.Namespace) -> int:
    caminho = Path(args.input)
    texto = caminho.read_text(encoding="utf-8", errors="replace")
    sequencias: dict[str, str] = {}

    if texto.lstrip().startswith(">"):
        nome_atual = ""
        buf: list[str] = []
        for linha in texto.splitlines():
            linha = linha.strip()
            if not linha:
                continue
            if linha.startswith(">"):
                if nome_atual:
                    sequencias[nome_atual] = "".join(buf)
                nome_atual = linha[1:].split()[0]
                buf = []
            else:
                buf.append(re.sub(r"\s+", "", linha))
        if nome_atual:
            sequencias[nome_atual] = "".join(buf)
    else:
        # Bloco MATRIX de NEXUS simples
        em_matriz = False
        for linha in texto.splitlines():
            s = linha.strip()
            if not s or s.startswith("["):
                continue
            if s.upper().startswith("MATRIX"):
                em_matriz = True
                continue
            if em_matriz and s.startswith(";"):
                break
            if em_matriz:
                partes = s.split()
                if len(partes) >= 2:
                    sequencias[partes[0]] = sequencias.get(partes[0], "") + partes[1]

    if not sequencias:
        raise SystemExit("Nenhum taxon lido da matriz morfologica.")

    comprimentos = {len(v) for v in sequencias.values()}
    if len(comprimentos) != 1:
        raise SystemExit(f"Comprimentos divergentes na matriz: {comprimentos}")
    n_char = next(iter(comprimentos))
    n_taxa = len(sequencias)

    invariaveis = 0
    autapomorficos = 0
    informativos = 0
    todos_estados: set[str] = set()

    for col in range(n_char):
        obs = [s[col] for s in sequencias.values() if s[col] not in ("?", "-")]
        todos_estados.update(obs)
        contagem = {est: obs.count(est) for est in set(obs)}
        if len(contagem) <= 1:
            invariaveis += 1
        elif sum(1 for c in contagem.values() if c >= 2) >= 2:
            informativos += 1
        else:
            autapomorficos += 1

    missing_por_taxon = {
        t: round(100.0 * sum(1 for c in seq if c in ("?", "-")) / n_char, 2)
        for t, seq in sequencias.items()
    }
    k_estados = max(2, len(todos_estados))
    precisa_asc = invariaveis == 0

    relatorio = {
        "n_taxa": n_taxa,
        "n_characters": n_char,
        "max_observed_states_k": k_estados,
        "invariant_characters": invariaveis,
        "autapomorphic_characters": autapomorficos,
        "parsimony_informative_characters": informativos,
        "ascertainment_bias_ASC_required": precisa_asc,
        "missing_percent_by_taxon": missing_por_taxon,
        "recommended_models": {
            "iqtree3": "MK+FQ+ASC+G4" if precisa_asc else "MK+FQ+G4",
            "revbayes": f'Q_morph := fnJC({k_estados}); seq ~ dnPhyloCTMC(tree=psi, Q=Q_morph, type="Standard", coding="{"variable" if precisa_asc else "all"}")',
            "mrbayes": f"lset coding={'variable' if precisa_asc else 'all'} rates=gamma;",
        },
    }
    registrar_proveniencia("mk-audit", {"input": str(caminho)}, relatorio)
    print(json.dumps(relatorio, indent=2, ensure_ascii=False))
    return 0


def cmd_altall_asr(args: argparse.Namespace) -> int:
    caminho = Path(args.state_file)
    no_alvo = args.node.strip()
    limiar = float(args.pp_threshold)

    colunas_estados: list[str] = []
    ml_seq: list[str] = []
    alt_seq: list[str] = []
    sitios_ambiguos: list[dict[str, Any]] = []

    with caminho.open(encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            s = linha.strip()
            if not s or s.startswith("#"):
                continue
            partes = s.split("\t")
            if partes[0] == "Node" and partes[1] == "Site":
                # Cabecalho: Node Site State p_A p_R ...
                colunas_estados = [c.removeprefix("p_") for c in partes[3:]]
                continue
            if partes[0] != no_alvo or not colunas_estados:
                continue
            sitio = int(partes[1])
            probs = [float(x) for x in partes[3 : 3 + len(colunas_estados)]]
            pares = sorted(zip(colunas_estados, probs, strict=False), key=lambda x: x[1], reverse=True)
            est1, pp1 = pares[0]
            est2, pp2 = pares[1] if len(pares) > 1 else (est1, 0.0)
            ml_seq.append(est1)
            if pp1 < limiar and pp2 >= args.min_alt_pp:
                alt_seq.append(est2)
                sitios_ambiguos.append(
                    {
                        "site": sitio,
                        "ml_state": est1,
                        "ml_pp": round(pp1, 4),
                        "alt_state": est2,
                        "alt_pp": round(pp2, 4),
                    }
                )
            else:
                alt_seq.append(est1)

    if not ml_seq:
        raise SystemExit(f"No '{no_alvo}' nao encontrado em {caminho}.")

    seq_ml_str = "".join(ml_seq)
    seq_alt_str = "".join(alt_seq)
    resumo = {
        "node": no_alvo,
        "length": len(ml_seq),
        "pp_threshold": limiar,
        "min_alt_pp": args.min_alt_pp,
        "n_swapped_sites": len(sitios_ambiguos),
        "percent_swapped": round(100.0 * len(sitios_ambiguos) / len(ml_seq), 2),
        "ml_sequence": seq_ml_str,
        "altall_sequence": seq_alt_str,
        "swapped_sites": sitios_ambiguos,
    }
    if args.out_fasta:
        out_p = Path(args.out_fasta)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        out_p.write_text(
            f">{no_alvo}_ML\n{seq_ml_str}\n>{no_alvo}_AltAll_pp{int(limiar * 100)}\n{seq_alt_str}\n",
            encoding="utf-8",
        )
        resumo["out_fasta"] = str(out_p)

    registrar_proveniencia(
        "altall-asr",
        {"state_file": str(caminho), "node": no_alvo, "pp_threshold": limiar},
        {"length": len(ml_seq), "n_swapped_sites": len(sitios_ambiguos)},
    )
    print(json.dumps(resumo, indent=2, ensure_ascii=False))
    return 0


def main() -> int:
    sys.stdout.reconfigure(encoding="utf-8")
    p = argparse.ArgumentParser(description="Toolkit de Paleontologia da bio-stack")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_pbdb = sub.add_parser("pbdb-range", help="Consulta FAD/LAD na Paleobiology Database API v1.2")
    p_pbdb.add_argument("--taxa", required=True, help="Lista de taxons separados por virgula")

    p_zooms = sub.add_parser("zooms-qc", help="Digestao in silico e desamidacao diagenetica (ZooMS)")
    p_zooms.add_argument("--sequence", required=True, help="Sequencia proteica (ex.: AMELX ou COL1A1)")
    p_zooms.add_argument("--min-len", type=int, default=6)
    p_zooms.add_argument("--max-len", type=int, default=35)

    p_mk = sub.add_parser("mk-audit", help="Audita matriz morfologica discreta para modelo Mk/Mkv (+ASC)")
    p_mk.add_argument("--input", required=True, help="Arquivo FASTA ou NEXUS com caracteres discretos")

    p_alt = sub.add_parser(
        "altall-asr",
        help="Extrai par ancestral ML + AltAll (Eick et al. 2017) de .state do IQ-TREE 3",
    )
    p_alt.add_argument("--state-file", required=True, help="Arquivo .state gerado por iqtree3 --ancestral")
    p_alt.add_argument("--node", required=True, help="Nome do no interno (ex.: Node1)")
    p_alt.add_argument("--pp-threshold", type=float, default=0.80, help="Limiar de PP para troca AltAll (padrao: 0.80)")
    p_alt.add_argument("--min-alt-pp", type=float, default=0.10, help="PP minima do 2o estado (padrao: 0.10)")
    p_alt.add_argument("--out-fasta", default="", help="Caminho opcional para gravar FASTA com ML e AltAll")

    args = p.parse_args()
    if args.cmd == "pbdb-range":
        return cmd_pbdb_range(args)
    if args.cmd == "zooms-qc":
        return cmd_zooms_qc(args)
    if args.cmd == "mk-audit":
        return cmd_mk_audit(args)
    if args.cmd == "altall-asr":
        return cmd_altall_asr(args)
    return 1


if __name__ == "__main__":
    sys.exit(main())
