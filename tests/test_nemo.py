# /// script
# requires-python = ">=3.10"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0", "numpy>=1.26", "pypdf>=5.0"]
# ///
"""Testa o servidor MCP nvidia-nemo. Nenhuma chamada de rede.

    uv run --script tests/test_nemo.py

As funcoes de rede (_embed, _post) sao substituidas por um dublê deterministico,
entao a suite exercita toda a logica de indexacao, fatiamento, persistencia e
busca sem consumir um unico token de credito da NVIDIA.
"""
import asyncio
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

SERVER = Path(__file__).resolve().parent.parent / "servers" / "nvidia_nemo_mcp" / "server.py"

spec = importlib.util.spec_from_file_location("nvnemo", SERVER)
srv = importlib.util.module_from_spec(spec)
sys.modules["nvnemo"] = srv
spec.loader.exec_module(srv)

fails = 0


def check(label, cond, detail=""):
    global fails
    fails += not cond
    print(f"{'OK   ' if cond else 'FALHA'} {label}" + (f"  -- {detail}" if not cond and detail else ""))


def unwrap(tool):
    return getattr(tool, "fn", tool)


# =============================================================================
# Dublê de rede: embedding deterministico por hash do texto. Textos parecidos
# nao ficam proximos, o que basta para testar o encanamento (nao a semantica).
# =============================================================================
EMBED_CALLS = []
RERANK_CALLS = []
DIM = 32


async def fake_embed(client, texts, input_type):
    EMBED_CALLS.append({"n": len(texts), "input_type": input_type})
    rng = np.zeros((len(texts), DIM), dtype=np.float32)
    for i, t in enumerate(texts):
        seed = abs(hash(t)) % (2**31)
        rng[i] = np.random.default_rng(seed).standard_normal(DIM)
    return rng / np.maximum(np.linalg.norm(rng, axis=1, keepdims=True), 1e-12)


async def fake_post(client, path, payload):
    RERANK_CALLS.append({"path": path, "model": payload.get("model"),
                         "n_passages": len(payload.get("passages", []))})
    n = len(payload.get("passages", []))
    # inverte a ordem, para ser obvio no teste se o rerank foi aplicado
    return {"rankings": [{"index": i, "logit": float(n - i)} for i in reversed(range(n))]}


srv._embed = fake_embed
srv._post = fake_post
srv.API_KEY = "nvapi-dubles-de-teste"


# =============================================================================
print("--- registro de ferramentas ---")
tools = asyncio.run(srv.mcp.list_tools())
names = sorted(t.name for t in tools)
print(f"registradas ({len(names)}): {', '.join(names)}\n")
expected = {"nemo_config", "library_index", "library_search", "embed_texts", "rerank_passages"}
check("todas registradas", expected <= set(names), f"faltam {expected - set(names)}")
check("toda ferramenta tem descricao", all(t.description for t in tools))


# =============================================================================
print("\n--- fatiamento ---")
check("texto curto vira 1 chunk", len(srv._chunk("abc " * 30)) == 1)
check("texto vazio vira 0 chunks", srv._chunk("   ") == [])

longo = ("Paragrafo sobre expressao diferencial. " * 40 + "\n\n") * 6
pedacos = srv._chunk(longo)
check("texto longo e fatiado", len(pedacos) > 1, f"{len(pedacos)} chunks")
check("chunk respeita o teto", all(len(p) <= srv.CHUNK_CHARS + 250 for p in pedacos),
      f"maior={max(len(p) for p in pedacos)}")
check("ha sobreposicao entre chunks", len(pedacos) < 2 or pedacos[0][-60:] in longo)
check("descarta fragmento minusculo", all(len(p) > 50 for p in pedacos))

check("hifen de quebra de linha e colado",
      "expressao" in srv._clean("expres-\nsao diferencial"),
      srv._clean("expres-\nsao diferencial"))
check("colapsa linhas em branco", "\n\n\n" not in srv._clean("a\n\n\n\n\nb"))


# =============================================================================
print("\n--- indexacao e busca ---")
with tempfile.TemporaryDirectory() as td:
    tmp = Path(td)
    srv.INDEX_DIR = tmp / "idx"
    srv.PROV_LOG = tmp / "prov.jsonl"

    docs = tmp / "artigos"
    docs.mkdir()
    (docs / "scrnaseq.md").write_text(
        "# Batch effect em single-cell\n\n" + ("Variacao tecnica induzida por criopreservacao "
        "em tecido congelado afeta a deteccao de genes. " * 25), encoding="utf-8")
    (docs / "deseq.txt").write_text(
        "Expressao diferencial requer contagens brutas. " * 60, encoding="utf-8")
    sub = docs / "metodos"
    sub.mkdir()
    (sub / "primers.md").write_text("Desenho de primers com Tm entre 58 e 62 graus. " * 40,
                                    encoding="utf-8")

    idx_fn = unwrap(srv.library_index)
    res = asyncio.run(idx_fn(folder=str(docs), name="teste"))
    check("indexa 3 arquivos (recursivo)", res["indexed_files"] == 3, str(res.get("indexed_files")))
    check("gera chunks", res["chunks"] > 3, str(res.get("chunks")))
    check("nenhuma falha", res["failed"] == [], str(res["failed"]))
    check("estima tokens cobrados", res["approx_tokens_billed"] > 0)
    check("embedding usa input_type=passage",
          all(c["input_type"] == "passage" for c in EMBED_CALLS), str(EMBED_CALLS))

    check("grava .npz e .jsonl",
          (srv.INDEX_DIR / "teste.npz").is_file() and (srv.INDEX_DIR / "teste.jsonl").is_file())

    # --- incremental ---
    EMBED_CALLS.clear()
    res2 = asyncio.run(idx_fn(folder=str(docs), name="teste"))
    check("segunda passada pula tudo", res2["indexed_files"] == 0, str(res2))
    check("nao chama embedding se nada mudou", EMBED_CALLS == [], str(EMBED_CALLS))
    check("reporta skipped", res2["skipped_unchanged"] == 3, str(res2.get("skipped_unchanged")))

    # --- arquivo modificado ---
    (docs / "deseq.txt").write_text("Conteudo totalmente novo sobre normalizacao. " * 60,
                                    encoding="utf-8")
    EMBED_CALLS.clear()
    res3 = asyncio.run(idx_fn(folder=str(docs), name="teste"))
    check("reindexa so o arquivo alterado", res3["indexed_files"] == 1, str(res3.get("indexed_files")))
    check("pula os dois inalterados", res3["skipped_unchanged"] == 2, str(res3.get("skipped_unchanged")))

    chunks = [json.loads(line) for line in (srv.INDEX_DIR / "teste.jsonl").read_text(encoding="utf-8").splitlines()]
    sources = {c["source"] for c in chunks}
    check("nao duplica trechos do arquivo alterado",
          len(chunks) == res3["chunks"], f"{len(chunks)} vs {res3['chunks']}")
    check("mantem os 3 arquivos no indice", len(sources) == 3, str(sources))

    # --- busca ---
    search = unwrap(srv.library_search)
    RERANK_CALLS.clear()
    EMBED_CALLS.clear()
    out = asyncio.run(search(query="o que causa variacao tecnica?", index="teste", top_k=3))
    check("devolve resultados", len(out["hits"]) == 3, str(len(out["hits"])))
    check("busca usa input_type=query",
          EMBED_CALLS and EMBED_CALLS[0]["input_type"] == "query", str(EMBED_CALLS))
    check("aplica reranking", out["reranked"] is True, str(out.get("reranked")))
    check("usa o endpoint de reranking",
          RERANK_CALLS and RERANK_CALLS[0]["path"] == "retrieval/nvidia/reranking",
          str(RERANK_CALLS))
    check("resultado traz logit do rerank", all("rerank_logit" in h for h in out["hits"]))
    check("resultado traz a fonte", all(h["source"] for h in out["hits"]))

    out_nr = asyncio.run(search(query="x", index="teste", top_k=2, rerank=False))
    check("rerank=False nao reranqueia", out_nr["reranked"] is False)

    # --- truncagem de contexto ---
    out_t = asyncio.run(search(query="x", index="teste", top_k=2, max_chars_per_hit=100))
    check("trunca texto do resultado",
          all(len(h["text"]) <= 130 for h in out_t["hits"]),
          str([len(h["text"]) for h in out_t["hits"]]))

    # --- rerank degrada com elegancia ---
    async def broken_post(client, path, payload):
        raise srv.NemoError("503 servico indisponivel")
    srv._post = broken_post
    out_deg = asyncio.run(search(query="x", index="teste", top_k=2))
    check("falha de rerank nao perde a busca", len(out_deg["hits"]) == 2, str(out_deg))
    check("avisa que o rerank falhou", "warning" in out_deg, str(out_deg.keys()))
    srv._post = fake_post

    # --- erros claros ---
    try:
        asyncio.run(search(query="x", index="inexistente"))
        check("indice inexistente da erro claro", False, "nao levantou")
    except srv.NemoError as e:
        check("indice inexistente da erro claro", "library_index" in str(e))

    try:
        asyncio.run(idx_fn(folder=str(tmp / "nao_existe")))
        check("pasta inexistente da erro claro", False, "nao levantou")
    except srv.NemoError as e:
        check("pasta inexistente da erro claro", "nao encontrada" in str(e))

    vazia = tmp / "vazia"
    vazia.mkdir()
    try:
        asyncio.run(idx_fn(folder=str(vazia)))
        check("pasta sem documentos da erro claro", False, "nao levantou")
    except srv.NemoError as e:
        check("pasta sem documentos da erro claro", "nenhum arquivo" in str(e))

    try:
        asyncio.run(idx_fn(folder=str(docs), name="t2", max_files=1))
        check("teto de arquivos protege", False, "nao levantou")
    except srv.NemoError as e:
        check("teto de arquivos protege", "teto" in str(e))

    # --- troca de modelo invalida o indice ---
    srv.EMBED_MODEL = "nvidia/outro-modelo"
    try:
        asyncio.run(search(query="x", index="teste"))
        check("modelo diferente e recusado", False, "nao levantou")
    except srv.NemoError as e:
        check("modelo diferente e recusado", "Reindexe" in str(e), str(e)[:90])
    srv.EMBED_MODEL = "nvidia/llama-3.2-nv-embedqa-1b-v2"

    # --- proveniencia ---
    check("grava proveniencia", srv.PROV_LOG.is_file())


# =============================================================================
print("\n--- credencial ausente ---")
srv.API_KEY = ""
try:
    srv._headers()
    check("erro claro sem chave", False, "nao levantou")
except srv.NemoError as e:
    check("erro claro sem chave", "NVIDIA_API_KEY" in str(e) and "01_setup_keys" in str(e))

print(f"\n{'TODOS PASSARAM' if not fails else str(fails) + ' FALHA(S)'}")
sys.exit(1 if fails else 0)
