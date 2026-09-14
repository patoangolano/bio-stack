# /// script
# requires-python = ">=3.10"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0", "numpy>=1.26"]
# ///
"""Testa o servidor MCP nvidia-bio. Nenhuma chamada de rede e feita.

    uv run --script tests/test_server.py

Cobre: registro de ferramentas, registro de modelos, validacao de sequencia,
mensagens de erro de credencial, modo local vs hospedado, e a garantia de que
coordenadas atomicas nunca voltam inline para o contexto.
"""
import asyncio
import importlib.util
import sys
import tempfile
from pathlib import Path

SERVER = Path(__file__).resolve().parent.parent / "servers" / "nvidia_bio_mcp" / "server.py"

spec = importlib.util.spec_from_file_location("nvbio", SERVER)
srv = importlib.util.module_from_spec(spec)
# Registrar antes de executar: @dataclass resolve anotacoes adiadas via
# sys.modules[cls.__module__], que seria None num modulo nao registrado.
sys.modules["nvbio"] = srv
spec.loader.exec_module(srv)

fails = 0


def check(label, cond, detail=""):
    global fails
    fails += not cond
    print(f"{'OK   ' if cond else 'FALHA'} {label}" + (f"  -- {detail}" if not cond and detail else ""))


def unwrap(tool):
    """FastMCP embrulha a funcao; pega a original para chamar direto."""
    return getattr(tool, "fn", tool)


# --- registro de ferramentas -------------------------------------------------
tools = asyncio.run(srv.mcp.list_tools())
names = sorted(t.name for t in tools)
print(f"ferramentas registradas ({len(names)}): {', '.join(names)}\n")

expected = {"nim_config", "nim_list_models", "nim_probe_endpoints", "fold_protein",
            "cofold_complex", "msa_search", "dock_ligand", "evo2_score", "nim_invoke",
            "design_backbone", "design_sequence", "generate_molecules",
            "generate_from_fragment"}
check("todas as ferramentas registradas", expected <= set(names), f"faltam {expected - set(names)}")
check("toda ferramenta tem descricao", all(t.description for t in tools))
def schema_of(t):
    """mcp 2.x renomeou inputSchema -> input_schema."""
    return getattr(t, "input_schema", None) or getattr(t, "inputSchema", None)


check("toda ferramenta tem schema", all(schema_of(t) for t in tools))

# --- registro de modelos -----------------------------------------------------
check("boltz2 no registro", "boltz2" in srv.REGISTRY)
check("caminho do boltz2 correto", srv.REGISTRY["boltz2"].path == "biology/mit/boltz2/predict",
      srv.REGISTRY["boltz2"].path)
check("boltz2 marcado como confirmado", srv.REGISTRY["boltz2"].confidence == "confirmed")
check(">=12 modelos registrados", len(srv.REGISTRY) >= 12, str(len(srv.REGISTRY)))

listed = unwrap(srv.nim_list_models)(task="folding")
check("filtro por tarefa funciona", listed and all(m["task"] == "folding" for m in listed),
      str(listed)[:120])

# --- validacao de sequencia --------------------------------------------------
check("aceita proteina valida", srv._validate_protein("MKVLA aGG") == "MKVLAAGG",
      srv._validate_protein("MKVLA aGG"))
try:
    srv._validate_protein("MKVL0ZB")
    check("rejeita caractere invalido", False, "nao levantou")
except srv.NimError as e:
    check("rejeita caractere invalido", "invalidos" in str(e))
try:
    srv._validate_protein("")
    check("rejeita sequencia vazia", False, "nao levantou")
except srv.NimError:
    check("rejeita sequencia vazia", True)
try:
    srv._validate_nucleotide("ACGTX")
    check("rejeita nucleotideo invalido", False, "nao levantou")
except srv.NimError:
    check("rejeita nucleotideo invalido", True)
check("aceita DNA valido", srv._validate_nucleotide("acg tN") == "ACGTN")

# --- erro claro sem chave de API ---------------------------------------------
srv.API_KEY = ""
srv.MODE = "hosted"
try:
    srv._headers()
    check("erro claro sem NVIDIA_API_KEY", False, "nao levantou")
except srv.NimError as e:
    check("erro claro sem NVIDIA_API_KEY", "NVIDIA_API_KEY" in str(e) and "01_setup_keys" in str(e))

# --- modo local dispensa chave -----------------------------------------------
srv.MODE = "local"
try:
    h = srv._headers()
    check("modo local nao exige chave", "Authorization" not in h)
except srv.NimError as e:
    check("modo local nao exige chave", False, str(e))
check("url local sem prefixo /v1",
      srv._url_for("biology/mit/boltz2/predict").endswith("localhost:8000/biology/mit/boltz2/predict"),
      srv._url_for("biology/mit/boltz2/predict"))

# --- nim_invoke recusa caminho fora de biology/ ------------------------------
try:
    asyncio.run(unwrap(srv.nim_invoke)(path="admin/secrets", payload={}))
    check("nim_invoke recusa caminho arbitrario", False, "nao levantou")
except srv.NimError as e:
    check("nim_invoke recusa caminho arbitrario", "biology/" in str(e))

# --- gravacao em disco e proveniencia ----------------------------------------
with tempfile.TemporaryDirectory() as td:
    srv.OUT_DIR = Path(td) / "nim"
    srv.PROV_LOG = Path(td) / "prov.jsonl"
    res = srv._handle_structure_response(
        {"structures": ["HEADER teste\nATOM      1  N   MET A   1       0.0   0.0   0.0\n"],
         "plddt": 88.5},
        stem="teste", tool="fold_protein", params={"model": "boltz2"})
    check("grava estrutura em disco", bool(res["files"]) and Path(res["files"][0]).is_file())
    check("detecta formato pdb", res["format"] == "pdb", res["format"])
    check("extrai confianca", res["confidence"] == {"plddt": 88.5}, str(res.get("confidence")))
    check("inclui interpretacao de pLDDT", "pLDDT" in (res.get("interpretation") or ""))
    check("coordenadas NAO voltam inline", "ATOM" not in str(res))
    check("grava proveniencia", srv.PROV_LOG.is_file())

    # mmCIF deve ser detectado como cif, nao pdb
    res2 = srv._handle_structure_response({"structures": ["data_teste\n_entry.id TESTE\n"]},
                                          stem="t2", tool="fold_protein", params={})
    check("detecta formato cif", res2["format"] == "cif", res2["format"])


# --- desenho de proteina: parsing do ProteinMPNN -----------------------------
# Literal capturado de uma chamada real ao NIM (crambina 1CRN, 2 sequencias).
MFASTA = (
    ">input, score=1.3360, global_score=1.3360, fixed_chains=[], designed_chains=['A'], "
    "model_name=v_48_002, git_hash=unknown, seed=699\n"
    "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN\n"
    ">T=0.1, sample=1, score=0.7900, global_score=0.7900, seq_recovery=0.5217\n"
    "TVCCPSKKALEEYLKCLKPGTPKEECAKATGCIIIKGKTCPKNYKY\n"
    ">T=0.1, sample=2, score=0.7489, global_score=0.7489, seq_recovery=0.5217\n"
    "TVCCPSKKALEEYKKCLKDGTPKKICAKKTGCIIIKGKKCPKNYAY\n"
)

recs = srv._parse_mfasta(MFASTA)
check("mfasta rende 3 registros", len(recs) == 3, str(len(recs)))
check("primeiro registro e a nativa", recs[0]["is_native"] and not recs[1]["is_native"])
check("sequencia nativa integra",
      recs[0]["sequence"] == "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
      recs[0]["sequence"])
check("le score do desenho", recs[1]["score"] == 0.79, str(recs[1]["score"]))
check("le seq_recovery", recs[1]["seq_recovery"] == 0.5217, str(recs[1]["seq_recovery"]))
# global_score nao pode vazar para o campo score
check("nao confunde score com global_score", recs[0]["score"] == 1.336, str(recs[0]["score"]))

# --- normalizacao de moleculas: dois formatos, um resultado ------------------
# MolMIM: JSON dentro de string, SMILES na chave "sample".
mm = srv._normalize_molecules(
    {"molecules": '[{"sample": "CCO", "score": 0.9}]', "score_type": "QED"})
check("molmim: desembrulha JSON-em-string", mm == [{"smiles": "CCO", "score": 0.9}], str(mm))

# GenMol: lista de verdade, SMILES na chave "smiles".
gm = srv._normalize_molecules({"molecules": [{"smiles": "c1ccccc1", "score": 0.44}]})
check("genmol: le lista direta", gm == [{"smiles": "c1ccccc1", "score": 0.44}], str(gm))

check("aceita lista de SMILES nus",
      srv._normalize_molecules({"molecules": ["CCO"]}) == [{"smiles": "CCO", "score": None}])
try:
    srv._normalize_molecules({"molecules": "isto nao e json"})
    check("erro claro se molecules nao for JSON", False, "nao levantou")
except srv.NimError as e:
    check("erro claro se molecules nao for JSON", "JSON" in str(e))
try:
    srv._normalize_molecules({"status": "success"})
    check("erro claro se faltar molecules", False, "nao levantou")
except srv.NimError as e:
    check("erro claro se faltar molecules", "molecules" in str(e))

# --- guardas de parametro (falham antes de gastar credito) -------------------
try:
    asyncio.run(unwrap(srv.generate_from_fragment)(fragment_smiles="c1ccccc1"))
    check("genmol exige mascara de crescimento", False, "nao levantou")
except srv.NimError as e:
    check("genmol exige mascara de crescimento", "[*" in str(e))

try:
    asyncio.run(unwrap(srv.generate_molecules)(smiles="CCO", num_molecules=50, particles=30))
    check("molmim exige particles >= num_molecules", False, "nao levantou")
except srv.NimError as e:
    check("molmim exige particles >= num_molecules", "particles" in str(e))

try:
    asyncio.run(unwrap(srv.design_backbone)(input_pdb_path="/nao/existe.pdb", contigs="A1-10"))
    check("design_backbone recusa PDB inexistente", False, "nao levantou")
except srv.NimError as e:
    check("design_backbone recusa PDB inexistente", "nao encontrado" in str(e))

# --- resposta do RFdiffusion e do ProteinMPNN em disco -----------------------
with tempfile.TemporaryDirectory() as td:
    srv.OUT_DIR = Path(td) / "nim"
    srv.PROV_LOG = Path(td) / "prov.jsonl"

    # Regressao: o RFdiffusion devolve a estrutura na chave `output_pdb`, que
    # nao estava na lista do handler e caia no fallback de JSON cru.
    res = srv._handle_structure_response(
        {"output_pdb": "ATOM      1  N   GLY A   1      -0.964  -7.023 -28.786  1.00\n",
         "elapsed_ms": 2071},
        stem="rfd", tool="design_backbone", params={"contigs": "A1-46/0 20-20"})
    check("reconhece output_pdb do RFdiffusion", res["format"] == "pdb", res["format"])
    check("nao cai no fallback de JSON cru",
          not res["files"][0].endswith("_raw.json"), res["files"][0])
    check("backbone gravado em disco", Path(res["files"][0]).is_file())
    check("coordenadas do backbone NAO voltam inline", "ATOM" not in str(res))

    # ProteinMPNN: a matriz probs e enorme e nunca pode voltar inline.
    probs = [[[0.05] * 21 for _ in range(46)] for _ in range(2)]
    resm = srv._handle_mpnn_response(
        {"mfasta": MFASTA, "scores": [0.79, 0.7489], "probs": probs},
        stem="mpnn", params={"backbone": "1CRN.pdb"})
    check("mpnn grava fasta e probs", len(resm["files"]) == 2, str(resm["files"]))
    check("mpnn conta so os desenhos", resm["num_designs"] == 2, str(resm["num_designs"]))
    check("mpnn devolve as sequencias desenhadas",
          resm["designs"][0]["sequence"].startswith("TVCCPSKK"), str(resm["designs"][0]))
    check("matriz probs NAO volta inline", resm.get("probs") is None)
    check("probs no disco, nao no contexto", len(str(resm)) < 2000, f"{len(str(resm))} chars")
    check("mpnn grava proveniencia", srv.PROV_LOG.is_file())

    try:
        srv._handle_mpnn_response({"scores": []}, stem="x", params={})
        check("erro claro se faltar mfasta", False, "nao levantou")
    except srv.NimError as e:
        check("erro claro se faltar mfasta", "mfasta" in str(e))

    # sampling_temp fora de faixa: precisa de um PDB real, a leitura vem antes.
    fake = Path(td) / "fake.pdb"
    fake.write_text("ATOM      1  N   MET A   1       0.0   0.0   0.0\n", encoding="utf-8")
    try:
        asyncio.run(unwrap(srv.design_sequence)(input_pdb_path=str(fake), sampling_temp=0.0))
        check("design_sequence valida sampling_temp", False, "nao levantou")
    except srv.NimError as e:
        check("design_sequence valida sampling_temp", "sampling_temp" in str(e))

    lixo = Path(td) / "lixo.pdb"
    lixo.write_text("isto nao e um pdb\n", encoding="utf-8")
    try:
        srv._read_pdb(str(lixo))
        check("recusa arquivo sem registros ATOM", False, "nao levantou")
    except srv.NimError as e:
        check("recusa arquivo sem registros ATOM", "ATOM" in str(e))

print(f"\n{'TODOS PASSARAM' if not fails else str(fails) + ' FALHA(S)'}")
sys.exit(1 if fails else 0)
