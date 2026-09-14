# /// script
# requires-python = ">=3.10"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0"]
# ///
"""Smoke test das ferramentas de desenho e quimica, contra os NIM reais.

    uv run --script tests/smoke_design.py [caminho.pdb]

ATENCAO — ao contrario do resto da suite, este script **faz chamadas de rede e
gasta credito da NVIDIA**. Por isso ele fica de fora do `tests/run_all.ps1`.
Rode de proposito, quando quiser confirmar que os quatro endpoints de desenho
continuam de pe e que os formatos de resposta nao mudaram.

Sem argumento, baixa a crambina (1CRN, 46 residuos) do RCSB como alvo — e a
menor proteina util para este teste, o que mantem o custo baixo.

Exercita o caminho completo, do payload ao disco:

    design_backbone -> design_sequence -> generate_molecules -> generate_from_fragment

Imprime so resumo: nenhuma coordenada, nenhuma matriz de probabilidade.
"""

import asyncio
import importlib.util
import sys
import tempfile
from pathlib import Path

import httpx

ROOT = Path(__file__).resolve().parent.parent
SERVER = ROOT / "servers" / "nvidia_bio_mcp" / "server.py"

spec = importlib.util.spec_from_file_location("nvbio", SERVER)
assert spec is not None and spec.loader is not None, f"nao consegui carregar {SERVER}"
srv = importlib.util.module_from_spec(spec)
# Registrar antes de executar: @dataclass resolve anotacoes adiadas via
# sys.modules[cls.__module__].
sys.modules["nvbio"] = srv
spec.loader.exec_module(srv)


def unwrap(tool):
    """FastMCP embrulha a funcao; pega a original para chamar direto."""
    return getattr(tool, "fn", tool)


def short(path_str: str) -> str:
    p = Path(path_str)
    return f"{p.name} ({p.stat().st_size} bytes)"


def buscar_crambina(destino: Path) -> Path:
    alvo = destino / "1CRN.pdb"
    r = httpx.get("https://files.rcsb.org/download/1CRN.pdb", timeout=60, follow_redirects=True)
    r.raise_for_status()
    alvo.write_text(r.text, encoding="utf-8")
    print(f"alvo   = 1CRN baixado do RCSB ({len(r.text)} bytes)")
    return alvo


async def rodar(pdb: str) -> int:
    print(f"modo   = {srv.MODE}")
    print(f"OUT    = {srv.OUT_DIR}")
    print(f"PROV   = {srv.PROV_LOG}\n")
    falhas = 0

    print("[1/4] design_backbone (rfdiffusion)")
    try:
        bb = await unwrap(srv.design_backbone)(
            input_pdb_path=pdb, contigs="A1-46/0 20-20", diffusion_steps=15, name="smoke"
        )
        print("      formato:", bb["format"], "|", [short(f) for f in bb["files"]])
        if "ATOM" in str(bb):
            print("      REGRESSAO: coordenadas voltaram inline")
            falhas += 1
        backbone = bb["files"][0]
    except srv.NimError as e:
        print("      FALHOU:", e)
        return 1  # sem backbone nao da para seguir

    print("\n[2/4] design_sequence (proteinmpnn) sobre o backbone gerado")
    try:
        ds = await unwrap(srv.design_sequence)(
            input_pdb_path=backbone, num_sequences=3, sampling_temp=0.2, name="smoke"
        )
        print("      desenhos:", ds["num_designs"], "|", [short(f) for f in ds["files"]])
        for d in ds["designs"]:
            print(f"      score={d['score']} recovery={d['seq_recovery']} {d['sequence'][:40]}...")
        if ds.get("probs") is not None:
            print("      REGRESSAO: matriz probs voltou inline")
            falhas += 1
    except srv.NimError as e:
        print("      FALHOU:", e)
        falhas += 1

    print("\n[3/4] generate_molecules (molmim, semente = aspirina)")
    try:
        mm = await unwrap(srv.generate_molecules)(
            smiles="CC(=O)Oc1ccccc1C(=O)O", property_name="QED",
            num_molecules=5, particles=20, iterations=10, name="smoke",
        )
        print(f"      geradas {mm['num_generated']}, unicas {mm['num_unique']} ({mm['score_type']})")
        for m in mm["molecules"][:3]:
            print(f"      {m['score']:.3f}  {m['smiles']}")
    except srv.NimError as e:
        print("      FALHOU:", e)
        falhas += 1

    print("\n[4/4] generate_from_fragment (genmol, fragmento = benzeno)")
    try:
        gm = await unwrap(srv.generate_from_fragment)(
            fragment_smiles="c1ccccc1[*{20-20}]", num_molecules=5, name="smoke"
        )
        print(f"      geradas {gm['num_generated']}, unicas {gm['num_unique']}")
        for m in gm["molecules"][:3]:
            print(f"      {m['score']}  {m['smiles']}")
    except srv.NimError as e:
        print("      FALHOU:", e)
        falhas += 1

    return falhas


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        alvo = Path(sys.argv[1]) if len(sys.argv) > 1 else buscar_crambina(Path(td))
        codigo = asyncio.run(rodar(str(alvo)))
    print("\n" + ("TODOS PASSARAM" if not codigo else f"{codigo} FALHA(S)"))
    sys.exit(1 if codigo else 0)
