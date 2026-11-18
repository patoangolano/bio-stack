"""Consolida e verifica cruzadamente todas as frentes da revisão de energia microbiana.

Lê os artefatos estruturados das cinco frentes em
``runs/energia_microbiana_2026-11-17/``, valida invariantes de contagem e
integridade (incluindo hashes SHA-256 dos modelos locais e caches NCBI), gera a
união desduplicada de identificadores bibliográficos e grava
``resumo_integrado.json``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    """Compute the SHA-256 hex digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    """Read a UTF-8 JSON object from disk."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Esperado objeto JSON em {path}, obtido {type(data).__name__}.")
    return data


def load_tsv(path: Path) -> list[dict[str, str]]:
    """Read a UTF-8 TSV file as a list of row dicts."""
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_doi(raw: str | None) -> str | None:
    """Normalize a DOI string by stripping prefixes and lowercasing."""
    if not raw:
        return None
    cleaned = raw.strip()
    cleaned = re.sub(r"^https?://(dx\.)?doi\.org/", "", cleaned, flags=re.IGNORECASE)
    cleaned = cleaned.strip().rstrip(".")
    return cleaned.lower() or None


def normalize_pmid(raw: str | int | None) -> str | None:
    """Extract numeric PMID string if present."""
    if raw is None:
        return None
    match = re.search(r"\b(\d+)\b", str(raw))
    return match.group(1) if match else None


def normalize_pmcid(raw: str | None) -> str | None:
    """Extract normalized PMC identifier if present."""
    if not raw:
        return None
    match = re.search(r"\b(PMC\d+)\b", raw.strip(), flags=re.IGNORECASE)
    return match.group(1).upper() if match else None


def collect_bibliography(
    versatilidade: dict[str, Any],
    diversidade: dict[str, Any],
    eficiencia: dict[str, Any],
    anotacoes: dict[str, Any],
    interpretacao_md: str,
) -> dict[str, Any]:
    """Deduplicate DOIs, PMIDs, and PMCIDs across literature, annotation, and GEM fronts."""
    dois: dict[str, dict[str, Any]] = {}
    pmids: dict[str, dict[str, Any]] = {}
    pmcids: dict[str, set[str]] = {}

    def register(
        front: str,
        doi: str | None = None,
        pmid: str | None = None,
        pmcid: str | None = None,
        title: str | None = None,
        year: int | str | None = None,
    ) -> None:
        n_doi = normalize_doi(doi)
        n_pmid = normalize_pmid(pmid)
        n_pmcid = normalize_pmcid(pmcid)
        if n_doi:
            entry = dois.setdefault(
                n_doi,
                {"doi": n_doi, "fronts": set(), "pmids": set(), "title": title, "year": year},
            )
            entry["fronts"].add(front)
            if n_pmid:
                entry["pmids"].add(n_pmid)
            if not entry["title"] and title:
                entry["title"] = title
        if n_pmid:
            entry = pmids.setdefault(
                n_pmid,
                {"pmid": n_pmid, "fronts": set(), "dois": set(), "title": title, "year": year},
            )
            entry["fronts"].add(front)
            if n_doi:
                entry["dois"].add(n_doi)
            if not entry["title"] and title:
                entry["title"] = title
        if n_pmcid:
            pmcids.setdefault(n_pmcid, set()).add(front)

    for ref in versatilidade.get("references", []):
        if ref.get("type") in {"database"}:
            continue
        register(
            "versatilidade",
            doi=ref.get("doi"),
            pmid=ref.get("pmid"),
            pmcid=ref.get("pmcid"),
            title=ref.get("title"),
            year=ref.get("year"),
        )

    for src in diversidade.get("fontes", []):
        register(
            "diversidade",
            doi=src.get("doi"),
            pmid=src.get("pmid"),
            pmcid=src.get("pmcid"),
            title=src.get("titulo"),
            year=src.get("ano"),
        )

    for src in eficiencia.get("sources", []):
        register(
            "eficiencia",
            doi=src.get("doi"),
            pmid=src.get("pmid"),
            pmcid=src.get("pmcid"),
            title=src.get("title"),
            year=src.get("year"),
        )

    for rep in anotacoes.get("representatives", []):
        for ref in rep.get("references", []):
            cit = ref.get("citation", {})
            doi_val = None
            pmid_val = None
            for xref in cit.get("citationCrossReferences", []):
                if xref.get("database") == "PubMed":
                    pmid_val = xref.get("id")
                elif xref.get("database") == "DOI":
                    doi_val = xref.get("id")
            register(
                "anotacoes_uniprot",
                doi=doi_val,
                pmid=pmid_val,
                title=cit.get("title"),
                year=cit.get("publicationDate"),
            )
        for ev in rep.get("annotation_evidence", []):
            if isinstance(ev, dict) and ev.get("source") == "PubMed":
                register("anotacoes_eco", pmid=ev.get("id"))

    for rec in anotacoes.get("kegg_records", []):
        for match in re.findall(r"PMID:(\d+)", json.dumps(rec)):
            register("anotacoes_kegg", pmid=match)

    gem_refs: list[dict[str, str | None]] = []
    for line in interpretacao_md.splitlines():
        if "https://doi.org/" not in line:
            continue
        doi_match = re.search(r"https://doi\.org/([^\s—)]+)", line)
        pmid_match = re.search(r"PMID\s+(\d+)", line)
        pmc_match = re.search(r"(PMC\d+)", line)
        doi_str = doi_match.group(1) if doi_match else None
        pmid_str = pmid_match.group(1) if pmid_match else None
        pmc_str = pmc_match.group(1) if pmc_match else None
        register("modelos_gem", doi=doi_str, pmid=pmid_str, pmcid=pmc_str)
        gem_refs.append(
            {
                "line": line.lstrip("- ").strip(),
                "doi": normalize_doi(doi_str),
                "pmid": normalize_pmid(pmid_str),
                "pmcid": normalize_pmcid(pmc_str),
            }
        )

    doi_records = sorted(
        (
            {
                "doi": rec["doi"],
                "fronts": sorted(rec["fronts"]),
                "pmids": sorted(rec["pmids"]),
                "title": rec["title"],
                "year": rec["year"],
            }
            for rec in dois.values()
        ),
        key=lambda item: item["doi"],
    )
    pmid_records = sorted(
        (
            {
                "pmid": rec["pmid"],
                "fronts": sorted(rec["fronts"]),
                "dois": sorted(rec["dois"]),
                "title": rec["title"],
                "year": rec["year"],
            }
            for rec in pmids.values()
        ),
        key=lambda item: int(item["pmid"]),
    )
    return {
        "unique_doi_count": len(doi_records),
        "unique_pmid_count": len(pmid_records),
        "unique_pmcid_count": len(pmcids),
        "dois": doi_records,
        "pmids": pmid_records,
        "pmcids": {key: sorted(val) for key, val in sorted(pmcids.items())},
        "gem_literature_verified": gem_refs,
    }


def verify_genome_panel(run_dir: Path) -> dict[str, Any]:
    """Load and validate the 79-entry NCBI genome panel and historical checks."""
    gen_dir = run_dir / "genomas"
    painel = load_tsv(gen_dir / "painel_final.tsv")
    resumo = load_json(gen_dir / "painel_final_resumo.json")
    inv_main = load_json(gen_dir / "inventario_ncbi.json")
    inv_res = load_json(gen_dir / "resolucao_nomes" / "inventario_ncbi.json")
    hist = load_json(gen_dir / "linhagens_historicas.json")

    if len(painel) != 79 or resumo.get("entries") != 79:
        raise ValueError(f"Painel esperado com 79 entradas, obtido {len(painel)}.")
    distinct_acc = {row["accession"] for row in painel}
    if len(distinct_acc) != 79 or resumo.get("distinct_accessions") != 79:
        raise ValueError(f"Esperados 79 acessos distintos, obtido {len(distinct_acc)}.")

    cache_checks: list[dict[str, Any]] = []
    for item in inv_main.get("ledger", []):
        cache_path = gen_dir / "ncbi_raw" / f"{item['id']}.json"
        actual_sha = sha256_file(cache_path)
        cache_checks.append(
            {
                "id": item["id"],
                "matches": actual_sha == item["cache_sha256"],
            }
        )
    if not all(entry["matches"] for entry in cache_checks):
        raise ValueError("Divergência de SHA-256 em cache bruto de genomas/ncbi_raw.")

    res_checks: list[dict[str, Any]] = []
    for item in inv_res.get("ledger", []):
        cache_path = gen_dir / "resolucao_nomes" / "ncbi_raw" / f"{item['id']}.json"
        actual_sha = sha256_file(cache_path)
        res_checks.append(
            {
                "id": item["id"],
                "matches": actual_sha == item["cache_sha256"],
            }
        )
    if not all(entry["matches"] for entry in res_checks):
        raise ValueError("Divergência de SHA-256 em cache bruto de resolucao_nomes/ncbi_raw.")

    historical_reports = []
    for report in hist.get("data", {}).get("reports", []):
        info = report.get("assembly_info", {})
        historical_reports.append(
            {
                "accession": report.get("accession"),
                "organism": report.get("organism", {}).get("organism_name"),
                "strain": report.get("organism", {}).get("infraspecific_names", {}).get("strain"),
                "assembly_status": info.get("assembly_status"),
                "assembly_level": info.get("assembly_level"),
                "paired_accession": report.get("paired_accession"),
            }
        )

    relaxed_entries = [
        {
            "id": row["id"],
            "taxon_requested": row["taxon_pedido"],
            "strain_filter_requested": row["linhagem_filtro_pedido"],
            "selected_accession": row["accession"],
            "selected_organism": row["organism"],
            "selected_strain": row["strain"],
            "level": row["level"],
        }
        for row in painel
        if row["resolucao"].startswith("Filtro de linhagem")
    ]

    return {
        "datasets_version": inv_main.get("datasets_version"),
        "panel_entries": len(painel),
        "distinct_accessions": len(distinct_acc),
        "levels": dict(Counter(row["level"] for row in painel)),
        "groups": dict(Counter(row["grupo"] for row in painel)),
        "refseq_vs_genbank": {
            "GCF_refseq": sum(1 for row in painel if row["accession"].startswith("GCF_")),
            "GCA_genbank": sum(1 for row in painel if row["accession"].startswith("GCA_")),
        },
        "main_queries_count": inv_main.get("queries_in_run"),
        "main_returned_assemblies": inv_main.get("returned_assemblies"),
        "resolution_queries_count": inv_res.get("queries_in_run"),
        "resolution_returned_assemblies": inv_res.get("returned_assemblies"),
        "raw_cache_sha256_verified": len(cache_checks) + len(res_checks),
        "pinned_results": resumo.get("pinned_results", {}),
        "historical_accessions_checked": historical_reports,
        "strain_filter_relaxed": relaxed_entries,
    }


def verify_cobrapy_models(run_dir: Path, repo_root: Path) -> dict[str, Any]:
    """Verify local GEM hashes and summarize the COBRApy audit."""
    audit = load_json(run_dir / "modelos" / "auditoria_cobrapy.json")
    models = audit.get("models", [])
    if len(models) != 4:
        raise ValueError(f"Esperados 4 modelos auditados, obtido {len(models)}.")
    if audit.get("ranking_inferred") is not False:
        raise ValueError("auditoria_cobrapy.json deve manter ranking_inferred = false.")

    model_summaries = []
    for model in models:
        src = repo_root / model["source_path"]
        actual_sha = sha256_file(src)
        if actual_sha != model["sha256"]:
            raise ValueError(f"SHA-256 divergente para {src}.")
        model_summaries.append(
            {
                "model_id": model["model_id"],
                "source_path": model["source_path"],
                "sha256_verified": True,
                "sha256": actual_sha,
                "genes": model["genes"],
                "reactions": model["reactions"],
                "metabolites": model["metabolites"],
                "objective": model["objective"],
                "medium_import_count": len(model["medium_import_bounds"]),
                "medium_import_bounds": model["medium_import_bounds"],
                "native_bounds_optimization": model["native_bounds_optimization"],
                "closed_imports_native_maintenance": model[
                    "closed_all_boundary_imports_native_maintenance"
                ],
                "closed_imports_atpm_test": model["closed_imports_atpm_test"],
            }
        )

    return {
        "python": audit.get("python"),
        "cobra_version": audit.get("cobra_version"),
        "solver": audit.get("solver"),
        "swiglpk_version": audit.get("swiglpk_version"),
        "ranking_inferred": audit.get("ranking_inferred"),
        "models_audited": len(model_summaries),
        "models": model_summaries,
    }


def build_notebooklm_cross_check(run_dir: Path) -> dict[str, Any]:
    """Summarize NotebookLM artifacts and document convergences and corrections."""
    nlm_dir = run_dir / "notebooklm"
    fontes = load_json(nlm_dir / "fontes_notebook.json")
    auditoria = load_json(nlm_dir / "auditoria_fontes.json")
    comparacao = load_json(nlm_dir / "comparacao_modelos.json")

    if fontes.get("notebook_id") != "4cb14b1e-7be6-439c-8ec6-f24686551300":
        raise ValueError("notebook_id divergente em fontes_notebook.json.")
    if fontes.get("source_count") != 59 or len(fontes.get("sources", [])) != 59:
        raise ValueError("Esperadas 59 fontes em fontes_notebook.json.")

    return {
        "notebook_id": fontes["notebook_id"],
        "url": fontes["url"],
        "source_count": fontes["source_count"],
        "queries": [
            {
                "file": "notebooklm/auditoria_fontes.json",
                "query_id": auditoria.get("query_id"),
                "status": auditoria.get("status"),
                "answer_chars": len(auditoria.get("result", {}).get("answer", "")),
            },
            {
                "file": "notebooklm/comparacao_modelos.json",
                "query_id": comparacao.get("query_id"),
                "status": comparacao.get("status"),
                "answer_chars": len(comparacao.get("result", {}).get("answer", "")),
            },
        ],
        "convergences": [
            (
                "Separação fenotípica entre R. palustris CGA009 e TIE-1 quanto à "
                "fotoferrotrofia/eletrotrofia experimentalmente demonstradas em TIE-1."
            ),
            (
                "Inexistência de microrganismo campeão universal entre os cinco GEMs "
                "(iDT1294, iCN1361, iJN678, iHN637, iAF692) e invalidade de comparar "
                "valores ótimos de FBA sob limites padrão e BOFs distintas."
            ),
            (
                "Consumo de H2 atmosférico (~530 ppbv) por hidrogenases de alta "
                "afinidade (ex.: Huc em M. smegmatis) associado predominantemente à "
                "manutenção/persistência celular sob inanição de carbono."
            ),
            (
                "Identificação de fontes ruidosas ou fora do escopo no acervo importado "
                "automaticamente (p. ex., página institucional, sRNA de motilidade, "
                "estruturas PDB isoladas, páginas ResearchGate e sínteses automáticas "
                "sem protocolo PRISMA)."
            ),
        ],
        "corrections_by_primary_and_db_audits": [
            {
                "topic": "R. palustris CGA009 — crescimento fotoautotrófico com H2 (hupV)",
                "notebooklm_statement": (
                    "auditoria_fontes.json lista os quatro modos tróficos em CGA009 "
                    "sem registrar a lesão em hupV."
                ),
                "verified_correction": (
                    "versatilidade.md (PMID 16923881) comprova deleção de 4 pb em hupV "
                    "impedindo fotoautotrofia com H2 em CGA009 nas condições testadas "
                    "(restaurada em CGA010); CGA0092 (PMID 37078874) não reporta reparo "
                    "de hupV."
                ),
            },
            {
                "topic": "R. palustris CGA009 — presença anotada de homólogos pioABC",
                "notebooklm_statement": (
                    "auditoria_fontes.json sugere ausência de fotoferrotrofia/captação "
                    "catódica e atribui o operon pioABC especificamente a TIE-1."
                ),
                "verified_correction": (
                    "anotacoes.md (UniProt 2026_04 e KEGG) confirma ortólogos candidatos "
                    "PioA/B/C também em CGA009 (TX73_003855/003850/003845; A0AAE9XUC8, "
                    "Q6NBT3, Q6NBT4): ausência de validação fenotípica em CGA009 não "
                    "equivale a ausência genômica anotada."
                ),
            },
            {
                "topic": "Dimensões dos modelos GEM publicados vs. arquivos locais",
                "notebooklm_statement": (
                    "comparacao_modelos.json cita dimensões dos artigos originais "
                    "(iJN678: 678 genes; iCN1361: 1.361 genes / 1.263 metabólitos; "
                    "iAF692: 619 reações)."
                ),
                "verified_correction": (
                    "auditoria_cobrapy.json mede os arquivos locais efetivos "
                    "(iJN678.json: 622 genes / 863 reações / 795 metabólitos; "
                    "iCN1361.xml: 1.362 genes / 1.292 reações / 1.265 metabólitos; "
                    "iAF692.json: 692 genes / 690 reações / 628 metabólitos), reforçando "
                    "que equivalência binária entre artigo e arquivo local não pode ser "
                    "presumida sem hash/origem."
                ),
            },
        ],
    }


def reconcile_pending_items(diversidade: dict[str, Any]) -> list[dict[str, str]]:
    """Map pending items raised by subagents to their post-integration status."""
    mapping = {
        "Inventario formal de acessos de montagem sob responsabilidade do principal.": (
            "FECHADO na Frente 5: painel formal com 79 montagens distintas verificado em "
            "genomas/painel_final.tsv e genomas/INVENTARIO.md."
        ),
        "Correspondencia entre genoma e cocultura de Methanosarcina barkeri.": (
            "PARCIALMENTE FECHADO: painel_final.tsv registra M. barkeri 227 "
            "(GCF_000970065.1, filtro relaxado) e linhagens_historicas.json documenta "
            "M. barkeri Fusaro (GCF_000195895.1, status suppressed). Permanece aberto o "
            "rastreamento direto da linhagem exata do ensaio de DIET (S40) vs. montagem."
        ),
        "Linhagem e doador exato de Chloroflexus no estudo selecionado.": (
            "PARCIALMENTE FECHADO: versatilidade.md (PMID 31158281) separa o isolado "
            "termofílico de Nakabusa (H2 fotoautotrófico e quimiolitotrófico) da cepa-tipo "
            "genômica J-10-fl (GCF_000018865.1); anotacoes.md valida J-10-fl (3HP bi-cycle "
            "A9WIU3/A9WEU4) excluindo Q6QQP7 (OK-70-fl)."
        ),
        "Linhagem de Chlamydomonas no estudo S05.": (
            "ABERTO (limitação mantida): texto completo de S05 (PMID 41570075) não foi "
            "auditado quanto ao identificador exato de sublinhagem; painel genômico usa "
            "CC-503 cw92 mt+ (GCF_000002595.2) apenas como referência da espécie."
        ),
        "Correspondencia entre material catodico e genoma de Geobacter.": (
            "ABERTO (limitação mantida): G. sulfurreducens PCA (GCF_000007985.2) e "
            "G. metallireducens GS-15 (GCF_000012925.1) representam os genomas de "
            "referência, sem vínculo verificado de sublinhagem com todo ensaio catódico."
        ),
        "Avaliacao institucional de risco, pureza e disponibilidade de cada material.": (
            "FORA DO ESCOPO computacional (mantido para comitê institucional e "
            "coleções de cultura antes de qualquer trabalho experimental)."
        ),
    }
    return [
        {
            "item": item,
            "status_pos_integracao": mapping.get(
                item, "ABERTO: verificar antes de extrapolar entre linhagens."
            ),
        }
        for item in diversidade.get("pendencias", [])
    ]


def build_summary(run_dir: Path, repo_root: Path) -> dict[str, Any]:
    """Assemble and validate the integrated synthesis JSON across all five fronts."""
    agentes_dir = run_dir / "agentes"
    versatilidade = load_json(agentes_dir / "versatilidade.json")
    diversidade = load_json(agentes_dir / "diversidade.json")
    eficiencia = load_json(agentes_dir / "eficiencia.json")
    anotacoes = load_json(agentes_dir / "anotacoes.json")
    interpretacao_md = (run_dir / "modelos" / "INTERPRETACAO.md").read_text(encoding="utf-8")

    if len(versatilidade.get("candidates", [])) != 20:
        raise ValueError("Esperados 20 candidatos em versatilidade.json.")
    if len(diversidade.get("taxa", [])) != 36 or diversidade.get("n_entradas") != 36:
        raise ValueError("Esperadas 36 entradas em diversidade.json.")
    if len(eficiencia.get("quantitative_evidence", [])) != 17:
        raise ValueError("Esperados 17 registros quantitativos em eficiencia.json.")
    if len(anotacoes.get("representatives", [])) != 48:
        raise ValueError("Esperados 48 representantes em anotacoes.json.")

    genome_panel = verify_genome_panel(run_dir)
    cobrapy_audit = verify_cobrapy_models(run_dir, repo_root)
    nlm_cross_check = build_notebooklm_cross_check(run_dir)
    bibliography = collect_bibliography(
        versatilidade, diversidade, eficiencia, anotacoes, interpretacao_md
    )

    preserved_legacy = [
        "scripts/bioenergetics_efficiency_analysis.py",
        "runs/bioenergetics_efficiency_atlas.json",
        "scripts/microbios_energia/01_verificar_genomas.sh",
        "scripts/microbios_energia/02_resolver_acessos.sh",
        "scripts/microbios_energia/03_baixar_gems.py",
        "runs/microbios_energia/genomas_ncbi.tsv",
        "runs/microbios_energia/genomas_resolvidos.tsv",
    ]
    preserved_status = {rel: (repo_root / rel).exists() for rel in preserved_legacy}
    if not all(preserved_status.values()):
        raise ValueError(f"Arquivo anterior esperado ausente: {preserved_status}")

    counts = {
        "panel_entries": genome_panel["panel_entries"],
        "distinct_accessions": genome_panel["distinct_accessions"],
        "versatility_candidates": len(versatilidade["candidates"]),
        "versatility_references": len(versatilidade["references"]),
        "diversity_entries": len(diversidade["taxa"]),
        "diversity_sources": len(diversidade["fontes"]),
        "efficiency_sources": len(eficiencia["sources"]),
        "efficiency_quantitative_records": len(eficiencia["quantitative_evidence"]),
        "efficiency_conditional_rankings": len(eficiencia["conditional_rankings"]),
        "annotation_taxa": len(anotacoes["taxa"]),
        "annotation_candidates_triaged": anotacoes["counts"]["candidate_records_before_selection"],
        "annotation_representatives": len(anotacoes["representatives"]),
        "annotation_representatives_uniprot": anotacoes["counts"]["module_representatives_uniprot"],
        "annotation_representatives_kegg": anotacoes["counts"]["module_representatives_kegg"],
        "annotation_excluded_cautions": len(anotacoes["excluded_cautions"]),
        "cobrapy_models_audited": cobrapy_audit["models_audited"],
        "notebooklm_sources": nlm_cross_check["source_count"],
        "unique_dois_across_fronts": bibliography["unique_doi_count"],
        "unique_pmids_across_fronts": bibliography["unique_pmid_count"],
        "unique_pmcids_across_fronts": bibliography["unique_pmcid_count"],
    }

    return {
        "generated_at": datetime.now(UTC).isoformat(),
        "protocol": "runs/energia_microbiana_2026-11-17/PROTOCOLO.md",
        "prior_materials_audit": (
            "runs/energia_microbiana_2026-11-17/AUDITORIA_MATERIAIS_ANTERIORES.md"
        ),
        "preserved_prior_files_present": preserved_status,
        "counts": counts,
        "fronts": {
            "frente_1_versatilidade": {
                "file": "agentes/versatilidade.md",
                "json": "agentes/versatilidade.json",
                "unit_of_analysis": versatilidade.get("unit_of_analysis"),
                "critical_corrections": versatilidade.get("critical_corrections", []),
            },
            "frente_2_diversidade": {
                "file": "agentes/diversidade.md",
                "json": "agentes/diversidade.json",
                "negative_control": "Haloferax volcanii DS2 (GCF_000025685.1)",
                "exclusions": diversidade.get("exclusoes", []),
            },
            "frente_3_eficiencia": {
                "file": "agentes/eficiencia.md",
                "json": "agentes/eficiencia.json",
                "universal_winner": eficiencia["conclusion"]["universal_winner"],
                "conditional_rankings": eficiencia["conditional_rankings"],
                "pareto_summary": eficiencia["pareto"],
                "evidence_gaps": eficiencia["evidence_gaps"],
            },
            "frente_4_anotacoes": {
                "file": "agentes/anotacoes.md",
                "json": "agentes/anotacoes.json",
                "script": "agentes/auditoria_anotacoes.py",
                "uniprot_release": anotacoes["request_log"][0]["uniprot_release"],
                "validation": anotacoes["validation"],
                "excluded_cautions": [
                    {
                        "accession": item["accession"],
                        "target": item["target"],
                        "reason": item["reason_excluded_from_module_confirmation"],
                    }
                    for item in anotacoes["excluded_cautions"]
                ],
            },
            "frente_5_integracao": {
                "genome_panel": genome_panel,
                "cobrapy_audit": cobrapy_audit,
                "notebooklm_cross_check": nlm_cross_check,
            },
        },
        "pendencias_reconciliacao": reconcile_pending_items(diversidade),
        "bibliography_union": bibliography,
    }


def main() -> None:
    """Run integrated cross-verification and write resumo_integrado.json."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path("runs/energia_microbiana_2026-11-17"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("runs/energia_microbiana_2026-11-17/resumo_integrado.json"),
    )
    args = parser.parse_args()

    summary = build_summary(args.run_dir, args.repo_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary["counts"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
