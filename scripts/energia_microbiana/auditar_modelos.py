"""Audita modelos locais sem transformar crescimento padrão em ranking de espécies.

Não altera os modelos-fonte. Os limites do meio são os do arquivo, não medidas
ambientais. A análise é de consistência e rastreabilidade, não de eficiência.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from datetime import UTC, datetime
from pathlib import Path

import cobra
from cobra.io import load_json_model, read_sbml_model
from cobra.util.solver import linear_reaction_coefficients


def solve(model: cobra.Model) -> dict:
    """Return optimizer status without converting infeasibility to zero growth."""
    solution = model.optimize()
    return {
        "status": solution.status,
        "objective_value": (
            float(solution.objective_value) if solution.status == "optimal" else None
        ),
    }


def audit(path: Path) -> dict:
    """Inspect a model and perform a bounded, explicitly non-comparative audit."""
    model = load_json_model(str(path)) if path.suffix == ".json" else read_sbml_model(str(path))
    model.solver = "glpk"
    objective = {
        reaction.id: float(coefficient)
        for reaction, coefficient in linear_reaction_coefficients(model).items()
    }
    result = {
        "source_path": path.as_posix(),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_provenance": (
            "Arquivo local preexistente; origem/versionamento externo não "
            "certificados nesta auditoria."
        ),
        "model_id": model.id,
        "model_name": model.name,
        "genes": len(model.genes),
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "objective": objective,
        "objective_direction": model.objective.direction,
        "medium_import_bounds": model.medium,
        "native_bounds_optimization": solve(model),
        "warnings": [
            "Import bounds are flux capacities, not environmental concentrations.",
            "Different default media/objectives prohibit cross-model ranking.",
            "A feasible FBA solution does not prove measured growth or efficiency.",
        ],
    }
    maintenance = [
        reaction
        for reaction in model.reactions
        if reaction.id.casefold() in {"atpm", "atpmaint", "dm_atp_c_"}
        or "atp maintenance" in (reaction.name or "").casefold()
    ]
    with model:
        for reaction in model.boundary:
            if len(reaction.metabolites) != 1:
                continue
            coefficient = next(iter(reaction.metabolites.values()))
            lower, upper = reaction.bounds
            reaction.bounds = (
                (max(0.0, lower), max(0.0, upper))
                if coefficient < 0
                else (min(0.0, lower), min(0.0, upper))
            )
        result["closed_all_boundary_imports_native_maintenance"] = solve(model)
        if len(maintenance) == 1:
            atpm = maintenance[0]
            native_bounds = atpm.bounds
            atpm.lower_bound = 0.0
            model.objective = atpm
            model.objective.direction = "max"
            result["closed_imports_atpm_test"] = {
                "reaction": atpm.id,
                "equation": atpm.reaction,
                "native_bounds": native_bounds,
                "test_bounds": atpm.bounds,
                "result": solve(model),
                "interpretation": (
                    "Fluxo ATPM positivo sem importações é alerta de ciclo/"
                    "fonte interna permissiva; fluxo zero não valida todo o modelo."
                ),
            }
        else:
            result["closed_imports_atpm_test"] = {
                "status": "not_performed",
                "reason": "ATPM inexistente ou ambígua; nenhuma reação foi inventada.",
                "candidate_ids": [reaction.id for reaction in maintenance],
            }
    return result


def main() -> None:
    """Audit four public-style local models in a separate execution environment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, default=Path("data/proc/gem"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    results = []
    for name in ["iCN1361.xml", "iJN678.json", "iHN637.json", "iAF692.json"]:
        path = args.models / name
        try:
            result = audit(path)
        except Exception as exc:  # noqa: BLE001
            result = {
                "source_path": path.as_posix(),
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
            }
        results.append(result)
    output = {
        "generated_at": datetime.now(UTC).isoformat(),
        "python": platform.python_version(),
        "cobra_version": cobra.__version__,
        "solver": "glpk",
        "swiglpk_version": importlib.metadata.version("swiglpk"),
        "ranking_inferred": False,
        "models": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    for result in results:
        print(
            json.dumps(
                {key: value for key, value in result.items() if key != "medium_import_bounds"},
                ensure_ascii=False,
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
