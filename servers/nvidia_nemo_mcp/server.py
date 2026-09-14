# /// script
# requires-python = ">=3.10"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0", "numpy>=1.26", "pypdf>=5.0"]
# ///
"""
nvidia-nemo-mcp - busca semantica na sua propria biblioteca, com NeMo Retriever.

Por que isto existe
-------------------
A mesma chave NVIDIA que roda os NIM de biologia da acesso ao NeMo Retriever
(incluido no NVIDIA AI Enterprise). Isso resolve um problema que todo aluno de
bioinformatica tem e ninguem resolve: uma pasta com 200 PDFs baixados ao longo
de dois anos, e nenhuma forma de achar "aquele artigo que falava sobre o efeito
de batch em scRNA-seq de tecido congelado".

Busca por palavra-chave falha nisso porque o artigo certo pode dizer
"cryopreservation-induced technical variation" e nunca a palavra "batch". Busca
semantica acha mesmo assim.

Pipeline (o padrao de RAG da NVIDIA):
    documento -> chunks -> embedding NIM -> indice local
    pergunta  -> embedding NIM -> busca por cosseno -> top-K
              -> reranking NIM -> top-N com ordem corrigida

O reranking e o passo que a maioria das implementacoes caseiras pula, e e o que
mais melhora a precisao: o embedding e rapido mas grosseiro, o reranker le a
pergunta e a passagem juntas e reordena.

Tudo local exceto as chamadas de embedding/rerank. Seus PDFs nunca saem do disco
inteiros - so os trechos que voce pedir para buscar.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx
import numpy as np

try:
    from mcp.server.mcpserver import MCPServer as _Server  # mcp >= 2.0
except ImportError:  # pragma: no cover
    from mcp.server.fastmcp import FastMCP as _Server  # mcp 1.x

# ---------------------------------------------------------------------------
# Configuracao
# ---------------------------------------------------------------------------

API_KEY = os.environ.get("NVIDIA_API_KEY", "").strip()
BASE = os.environ.get("NVIDIA_NEMO_BASE", "https://integrate.api.nvidia.com/v1")

# llama-3.2-nv-embedqa-1b-v2 cobre 26 idiomas. Importa aqui: artigos em ingles,
# anotacoes suas em portugues, e busca crosslingual entre os dois.
EMBED_MODEL = os.environ.get("NVIDIA_EMBED_MODEL", "nvidia/llama-3.2-nv-embedqa-1b-v2")
RERANK_MODEL = os.environ.get("NVIDIA_RERANK_MODEL", "nvidia/llama-3.2-nv-rerankqa-1b-v2")

INDEX_DIR = Path(os.environ.get("BIO_STACK_INDEX", "./runs/library")).resolve()
PROV_LOG = Path(os.environ.get("BIO_STACK_PROVENANCE", "./runs/provenance.jsonl")).resolve()

# Chunk de ~1200 chars com 200 de sobreposicao. Abaixo disso o trecho perde
# contexto; acima, o embedding dilui o assunto especifico.
CHUNK_CHARS = int(os.environ.get("BIO_STACK_CHUNK", "1200"))
CHUNK_OVERLAP = 200
EMBED_BATCH = 64            # a API aceita lote; 64 e o equilibrio latencia/payload
RERANK_MAX_PASSAGES = 512   # limite documentado do endpoint
TIMEOUT = 180.0

SUPPORTED = {".pdf", ".txt", ".md", ".markdown", ".rst"}

mcp = _Server("nvidia-nemo")


class NemoError(RuntimeError):
    """Falha previsivel, com mensagem acionavel para o modelo."""


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------


def _headers() -> dict[str, str]:
    if not API_KEY:
        raise NemoError(
            "NVIDIA_API_KEY nao esta definida. Rode bootstrap/01_setup_keys.ps1 "
            "e reinicie o Claude Code."
        )
    return {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }


async def _post(client: httpx.AsyncClient, path: str, payload: dict[str, Any]) -> dict[str, Any]:
    url = f"{BASE.rstrip('/')}/{path.lstrip('/')}"
    r = await client.post(url, headers=_headers(), json=payload, timeout=TIMEOUT)

    if r.status_code == 401:
        raise NemoError("HTTP 401: NVIDIA_API_KEY invalida ou expirada.")
    if r.status_code == 403:
        raise NemoError(
            f"HTTP 403: sua chave nao tem acesso a {payload.get('model')}. "
            "NeMo Retriever faz parte do NVIDIA AI Enterprise - confira seu entitlement."
        )
    if r.status_code == 429:
        raise NemoError("HTTP 429: limite de taxa. Reduza o lote ou espere.")
    if r.status_code >= 400:
        raise NemoError(f"HTTP {r.status_code} em {url}: {r.text[:500]}")

    return r.json()


async def _embed(client: httpx.AsyncClient, texts: list[str], input_type: str) -> np.ndarray:
    """
    input_type e a parte que quase todo mundo erra: os textos do acervo sao
    'passage' e a pergunta e 'query'. Usar o mesmo para os dois degrada
    bastante a qualidade da busca, sem dar erro nenhum.
    """
    if input_type not in ("passage", "query"):
        raise NemoError(f"input_type invalido: {input_type}")

    out: list[list[float]] = []
    for i in range(0, len(texts), EMBED_BATCH):
        batch = texts[i : i + EMBED_BATCH]
        data = await _post(client, "embeddings", {
            "input": batch,
            "model": EMBED_MODEL,
            "input_type": input_type,
            "encoding_format": "float",
            "truncate": "END",
        })
        rows = sorted(data.get("data", []), key=lambda d: d.get("index", 0))
        if len(rows) != len(batch):
            raise NemoError(f"API devolveu {len(rows)} embeddings para {len(batch)} textos.")
        out.extend(r["embedding"] for r in rows)

    arr = np.asarray(out, dtype=np.float32)
    # Normalizar aqui deixa a busca por cosseno virar um produto de matriz.
    norms = np.linalg.norm(arr, axis=1, keepdims=True)
    return arr / np.maximum(norms, 1e-12)


# ---------------------------------------------------------------------------
# Extracao e fatiamento
# ---------------------------------------------------------------------------


def _read_document(path: Path) -> str:
    if path.suffix.lower() == ".pdf":
        try:
            from pypdf import PdfReader
        except ImportError as exc:  # pragma: no cover
            raise NemoError("pypdf nao instalado - nao consigo ler PDF.") from exc
        try:
            reader = PdfReader(str(path))
            return "\n\n".join((p.extract_text() or "") for p in reader.pages)
        except Exception as exc:
            raise NemoError(f"falha ao ler {path.name}: {exc}") from exc
    return path.read_text(encoding="utf-8", errors="replace")


def _clean(text: str) -> str:
    # PDFs de duas colunas quebram palavras com hifen no fim da linha.
    text = re.sub(r"(\w)-\n(\w)", r"\1\2", text)
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def _chunk(text: str) -> list[str]:
    """Fatia em pedacos com sobreposicao, cortando em fronteira de paragrafo."""
    text = _clean(text)
    if len(text) <= CHUNK_CHARS:
        return [text] if text else []

    chunks: list[str] = []
    start = 0
    while start < len(text):
        end = min(start + CHUNK_CHARS, len(text))
        if end < len(text):
            # procura uma quebra natural na ultima quarta parte do chunk
            window = text[start + (CHUNK_CHARS * 3 // 4) : end]
            for sep in ("\n\n", ". ", "\n"):
                idx = window.rfind(sep)
                if idx != -1:
                    end = start + (CHUNK_CHARS * 3 // 4) + idx + len(sep)
                    break
        piece = text[start:end].strip()
        if len(piece) > 50:  # fragmentos minusculos so poluem o indice
            chunks.append(piece)
        if end >= len(text):
            break
        start = max(end - CHUNK_OVERLAP, start + 1)
    return chunks


def _file_hash(path: Path) -> str:
    h = hashlib.sha256()
    h.update(str(path.stat().st_size).encode())
    with path.open("rb") as fh:
        h.update(fh.read(65536))
    return h.hexdigest()[:16]


# ---------------------------------------------------------------------------
# Indice em disco
# ---------------------------------------------------------------------------


@dataclass
class Index:
    vectors: np.ndarray            # (n, d) normalizados
    chunks: list[dict[str, Any]]   # {source, page_hint, text, file_hash}
    model: str

    @staticmethod
    def paths(name: str) -> tuple[Path, Path]:
        # `name` vem de parametro de ferramenta: sem sanitizar, "../../x"
        # escreveria fora de INDEX_DIR.
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", str(name)).strip("._")[:64]
        if not safe:
            raise NemoError(f"nome de indice invalido: {name!r}")
        return INDEX_DIR / f"{safe}.npz", INDEX_DIR / f"{safe}.jsonl"

    @classmethod
    def load(cls, name: str) -> Index | None:
        vpath, cpath = cls.paths(name)
        if not (vpath.is_file() and cpath.is_file()):
            return None
        with np.load(vpath) as z:
            vecs = z["vectors"]
            model = str(z["model"]) if "model" in z else EMBED_MODEL
        chunks = [json.loads(line) for line in cpath.read_text(encoding="utf-8").splitlines() if line.strip()]
        if len(chunks) != len(vecs):
            raise NemoError(
                f"indice '{name}' corrompido: {len(vecs)} vetores para {len(chunks)} trechos. "
                "Reindexe com force=true."
            )
        return cls(vectors=vecs, chunks=chunks, model=model)

    def save(self, name: str) -> None:
        INDEX_DIR.mkdir(parents=True, exist_ok=True)
        vpath, cpath = self.paths(name)
        np.savez_compressed(vpath, vectors=self.vectors, model=np.array(self.model))
        with cpath.open("w", encoding="utf-8") as fh:
            for c in self.chunks:
                fh.write(json.dumps(c, ensure_ascii=False) + "\n")


def _log(tool: str, params: dict[str, Any], result: dict[str, Any]) -> None:
    try:
        PROV_LOG.parent.mkdir(parents=True, exist_ok=True)
        entry = {
            "ts": datetime.now(UTC).isoformat(timespec="seconds"),
            "server": "nvidia-nemo-mcp",
            "tool": tool,
            "params": {k: v for k, v in params.items() if "key" not in k.lower()},
            "summary": {k: result[k] for k in ("indexed_files", "chunks", "hits") if k in result},
        }
        with PROV_LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Ferramentas
# ---------------------------------------------------------------------------


@mcp.tool()
def nemo_config() -> dict[str, Any]:
    """
    Mostra a configuracao do servidor: endpoint, modelos, se a chave esta
    presente, e quais indices ja existem. Chame primeiro quando algo falhar.
    """
    indices = []
    if INDEX_DIR.is_dir():
        for v in sorted(INDEX_DIR.glob("*.npz")):
            c = INDEX_DIR / f"{v.stem}.jsonl"
            if not c.is_file():
                continue  # .npz sem .jsonl nao e um indice pesquisavel
            n = sum(1 for _ in c.open(encoding="utf-8"))
            indices.append({"name": v.stem, "chunks": n,
                            "size_mb": round(v.stat().st_size / 1e6, 2)})

    return {
        "api_key_present": bool(API_KEY),
        "base_url": BASE,
        "embed_model": EMBED_MODEL,
        "rerank_model": RERANK_MODEL,
        "index_dir": str(INDEX_DIR),
        "chunk_chars": CHUNK_CHARS,
        "indices": indices,
        "supported_formats": sorted(SUPPORTED),
    }


@mcp.tool()
async def library_index(
    folder: str,
    name: str = "default",
    recursive: bool = True,
    force: bool = False,
    max_files: int = 500,
) -> dict[str, Any]:
    """
    Indexa uma pasta de artigos (PDF, txt, md) para busca semantica.

    Incremental por padrao: arquivos ja indexados e nao modificados sao pulados,
    entao rodar de novo depois de baixar 3 PDFs novos custa so os 3.

    Args:
        folder: pasta com os documentos.
        name: nome do indice (permite acervos separados, ex. "tcc", "metodos").
        recursive: descer em subpastas.
        force: reindexar tudo, ignorando o que ja existe.
        max_files: teto de seguranca contra apontar para a pasta errada.
    """
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise NemoError(f"pasta nao encontrada: {root}")

    pattern = "**/*" if recursive else "*"
    files = sorted(
        p for p in root.glob(pattern)
        if p.is_file() and p.suffix.lower() in SUPPORTED
    )
    if not files:
        raise NemoError(
            f"nenhum arquivo suportado em {root}. Formatos aceitos: {sorted(SUPPORTED)}"
        )
    if len(files) > max_files:
        raise NemoError(
            f"{len(files)} arquivos encontrados, acima do teto de {max_files}. "
            "Confirme que a pasta esta certa e aumente max_files se estiver."
        )

    existing = None if force else Index.load(name)
    known: dict[str, str] = {}
    kept_chunks: list[dict[str, Any]] = []
    kept_vectors: list[np.ndarray] = []

    if existing is not None:
        if existing.model != EMBED_MODEL:
            raise NemoError(
                f"o indice '{name}' foi criado com {existing.model}, mas o modelo atual e "
                f"{EMBED_MODEL}. Vetores de modelos diferentes nao sao comparaveis - "
                "reindexe com force=true."
            )
        for i, c in enumerate(existing.chunks):
            known[c["source"]] = c.get("file_hash", "")
            kept_chunks.append(c)
            kept_vectors.append(existing.vectors[i])

    new_chunks: list[dict[str, Any]] = []
    skipped, failed, indexed = 0, [], 0

    for path in files:
        rel = str(path.relative_to(root))
        try:
            fh = _file_hash(path)
        except OSError as exc:
            failed.append(f"{rel}: {exc}")
            continue

        if known.get(rel) == fh:
            skipped += 1
            continue

        # arquivo mudou: remove os trechos antigos dele
        if rel in known:
            drop = {i for i, c in enumerate(kept_chunks) if c["source"] == rel}
            kept_chunks = [c for i, c in enumerate(kept_chunks) if i not in drop]
            kept_vectors = [v for i, v in enumerate(kept_vectors) if i not in drop]

        try:
            text = _read_document(path)
        except NemoError as exc:
            failed.append(f"{rel}: {exc}")
            continue

        pieces = _chunk(text)
        if not pieces:
            failed.append(f"{rel}: sem texto extraivel (PDF escaneado? precisa de OCR)")
            continue

        for j, piece in enumerate(pieces):
            new_chunks.append({
                "source": rel, "chunk": j, "text": piece,
                "file_hash": fh, "chars": len(piece),
            })
        indexed += 1

    if not new_chunks and existing is not None:
        result = {
            "indexed_files": 0, "skipped_unchanged": skipped,
            "chunks": len(kept_chunks), "failed": failed,
            "note": "nada mudou desde a ultima indexacao.",
        }
        _log("library_index", {"folder": str(root), "name": name}, result)
        return result

    if not new_chunks:
        raise NemoError(
            "nenhum arquivo produziu texto indexavel. "
            + (f"Falhas: {failed[:5]}" if failed else "")
            + " PDFs escaneados precisam de OCR antes de indexar."
        )

    async with httpx.AsyncClient() as client:
        vectors = await _embed(client, [c["text"] for c in new_chunks], "passage")

    all_chunks = kept_chunks + new_chunks
    all_vectors = (
        np.vstack([np.asarray(kept_vectors, dtype=np.float32).reshape(-1, vectors.shape[1]), vectors])
        if kept_vectors else vectors
    )

    Index(vectors=all_vectors, chunks=all_chunks, model=EMBED_MODEL).save(name)

    total_chars = sum(c["chars"] for c in new_chunks)
    result = {
        "index": name,
        "indexed_files": indexed,
        "skipped_unchanged": skipped,
        "new_chunks": len(new_chunks),
        "chunks": len(all_chunks),
        "embed_model": EMBED_MODEL,
        "approx_tokens_billed": total_chars // 4,
        "failed": failed,
        "next_step": f"use library_search(query, index='{name}') para buscar.",
    }
    _log("library_index", {"folder": str(root), "name": name, "files": indexed}, result)
    return result


@mcp.tool()
async def library_search(
    query: str,
    index: str = "default",
    top_k: int = 8,
    candidates: int = 40,
    rerank: bool = True,
    max_chars_per_hit: int = 700,
) -> dict[str, Any]:
    """
    Busca semantica na biblioteca indexada.

    Faz o pipeline completo da NVIDIA: embedding da pergunta -> busca por
    cosseno -> reranking dos melhores candidatos. O reranking corrige a ordem
    lendo pergunta e passagem juntas, o que o embedding sozinho nao faz.

    Args:
        query: a pergunta, em linguagem natural. Frase completa funciona melhor
               que palavras-chave soltas - o modelo foi treinado para isso.
        index: qual acervo consultar.
        top_k: quantos resultados devolver.
        candidates: quantos candidatos buscar antes de reranquear.
        rerank: aplicar o reranking NIM (recomendado; custa uma chamada extra).
        max_chars_per_hit: teto de texto por resultado, para poupar contexto.
    """
    idx = Index.load(index)
    if idx is None:
        raise NemoError(
            f"indice '{index}' nao existe. Rode library_index(folder=...) primeiro. "
            f"Indices disponiveis: {[i['name'] for i in nemo_config()['indices']]}"
        )
    if idx.model != EMBED_MODEL:
        raise NemoError(
            f"indice criado com {idx.model}, modelo atual {EMBED_MODEL}. Reindexe com force=true."
        )

    candidates = max(top_k, min(candidates, len(idx.chunks), RERANK_MAX_PASSAGES))

    async with httpx.AsyncClient() as client:
        qvec = await _embed(client, [query], "query")
        scores = idx.vectors @ qvec[0]
        order = np.argsort(-scores)[:candidates]

        hits = [{
            "source": idx.chunks[i]["source"],
            "chunk": idx.chunks[i]["chunk"],
            "cosine": round(float(scores[i]), 4),
            "text": idx.chunks[i]["text"],
        } for i in order]

        reranked = False
        if rerank and len(hits) > 1:
            try:
                data = await _post(client, "retrieval/nvidia/reranking", {
                    "model": RERANK_MODEL,
                    "query": {"text": query},
                    "passages": [{"text": h["text"]} for h in hits],
                    "truncate": "END",
                })
                rankings = data.get("rankings") or []
                if rankings:
                    ordered = []
                    for r in rankings:
                        j = r.get("index")
                        if isinstance(j, int) and 0 <= j < len(hits):
                            h = dict(hits[j])
                            h["rerank_logit"] = round(float(r.get("logit", 0.0)), 4)
                            ordered.append(h)
                    if ordered:
                        hits = ordered
                        reranked = True
            except NemoError as exc:
                # Reranking e melhoria, nao requisito: se falhar, entrega a
                # ordem do cosseno e avisa, em vez de perder a busca inteira.
                hits = hits[:top_k]
                result = {
                    "hits": [_trim(h, max_chars_per_hit) for h in hits],
                    "reranked": False,
                    "warning": f"reranking falhou ({exc}); ordem por cosseno apenas.",
                }
                _log("library_search", {"query": query, "index": index}, result)
                return result

    hits = hits[:top_k]
    result = {
        "query": query,
        "index": index,
        "hits": [_trim(h, max_chars_per_hit) for h in hits],
        "reranked": reranked,
        "searched_chunks": len(idx.chunks),
        "note": "cosseno mede similaridade grosseira; rerank_logit e o criterio "
                "confiavel quando presente (maior = mais relevante).",
    }
    _log("library_search", {"query": query, "index": index}, {"hits": len(hits)})
    return result


def _trim(hit: dict[str, Any], limit: int) -> dict[str, Any]:
    h = dict(hit)
    if len(h["text"]) > limit:
        h["text"] = h["text"][:limit] + f"... [+{len(hit['text']) - limit} chars]"
    return h


@mcp.tool()
async def embed_texts(texts: list[str], input_type: str = "passage") -> dict[str, Any]:
    """
    Acesso direto ao endpoint de embedding. Util para analise propria:
    agrupar abstracts, medir similaridade entre resumos, alimentar um
    classificador.

    Os vetores sao gravados em disco (sao grandes) e o retorno traz o caminho.

    Args:
        texts: lista de textos.
        input_type: "passage" para documentos, "query" para perguntas.
    """
    if not texts:
        raise NemoError("lista de textos vazia.")
    if len(texts) > 512:
        raise NemoError(f"{len(texts)} textos excede o teto de 512 por chamada.")

    async with httpx.AsyncClient() as client:
        vecs = await _embed(client, texts, input_type)

    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    stem = hashlib.sha256("".join(texts).encode()).hexdigest()[:12]
    # Subdiretorio proprio: no mesmo nivel dos indices, nemo_config listava
    # estes arquivos como se fossem acervos pesquisaveis.
    avulsos = INDEX_DIR / "embeddings"
    avulsos.mkdir(parents=True, exist_ok=True)
    out = avulsos / f"{stem}.npz"
    np.savez_compressed(out, vectors=vecs, model=np.array(EMBED_MODEL))

    return {
        "file": str(out),
        "shape": list(vecs.shape),
        "model": EMBED_MODEL,
        "normalized": True,
        "note": "vetores normalizados (norma 1): similaridade de cosseno = produto escalar. "
                "Carregue com np.load(file)['vectors'].",
    }


@mcp.tool()
async def rerank_passages(query: str, passages: list[str], top_k: int = 10) -> dict[str, Any]:
    """
    Reordena uma lista de passagens pela relevancia a uma pergunta.

    Use quando ja tiver candidatos de outra fonte - resultados do PubMed, trechos
    que voce colou, saidas de outra busca - e quiser a ordem certa.

    Args:
        query: a pergunta.
        passages: textos candidatos (maximo 512).
        top_k: quantos devolver.
    """
    if not passages:
        raise NemoError("lista de passagens vazia.")
    if len(passages) > RERANK_MAX_PASSAGES:
        raise NemoError(f"{len(passages)} passagens excede o limite de {RERANK_MAX_PASSAGES}.")

    async with httpx.AsyncClient() as client:
        data = await _post(client, "retrieval/nvidia/reranking", {
            "model": RERANK_MODEL,
            "query": {"text": query},
            "passages": [{"text": p} for p in passages],
            "truncate": "END",
        })

    out = []
    for r in (data.get("rankings") or [])[:top_k]:
        i = r.get("index")
        if isinstance(i, int) and 0 <= i < len(passages):
            out.append({
                "original_index": i,
                "logit": round(float(r.get("logit", 0.0)), 4),
                "text": passages[i][:700],
            })

    return {"query": query, "model": RERANK_MODEL, "ranked": out, "total_input": len(passages)}


if __name__ == "__main__":
    mcp.run()
