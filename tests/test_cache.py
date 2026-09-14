# /// script
# requires-python = ">=3.10"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0", "numpy>=1.26"]
# ///
"""
Testa o cache enderecado por conteudo do servidor nvidia-bio.

O que importa verificar, em ordem de quanto doi errar:

1. **Chave estavel sob reordenacao.** Dois payloads com as mesmas chaves em
   ordem diferente precisam bater no mesmo arquivo. Sem `sort_keys`, o cache
   parece funcionar em teste e nunca acerta em producao.
2. **Payload diferente, chave diferente.** Um cache que colide devolve a
   estrutura da sequencia errada — o pior erro possivel aqui, porque o
   resultado parece plausivel.
3. **Ida e volta.** O que foi gravado e o que volta.
4. **Desligavel.** Com BIO_STACK_NIM_CACHE=0 nao le nem grava.

Nao chama a API. Testar rede aqui tornaria o teste lento e dependente de cota.

    uv run --script tests/test_cache.py
"""

from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent


def carrega(tmp: Path, cache_ligado: bool = True):
    """
    Importa o servidor com OUT_DIR apontando para um diretorio temporario.

    As constantes do modulo (OUT_DIR, CACHE_DIR, CACHE_ON) sao lidas do ambiente
    no momento do import, entao o ambiente precisa estar montado antes — e o
    modulo precisa ser recarregado a cada variacao.
    """
    os.environ["BIO_STACK_OUT"] = str(tmp)
    os.environ["BIO_STACK_PROVENANCE"] = str(tmp / "prov.jsonl")
    os.environ["BIO_STACK_NIM_CACHE"] = "1" if cache_ligado else "0"
    os.environ.setdefault("NVIDIA_API_KEY", "nvapi-teste-nao-usada")

    spec = importlib.util.spec_from_file_location(
        f"nbs_{cache_ligado}_{tmp.name}", RAIZ / "servers" / "nvidia_bio_mcp" / "server.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("nao consegui montar o spec de import do servidor nvidia-bio")
    mod = importlib.util.module_from_spec(spec)
    # Registrar antes de executar nao e opcional: o `mcp` monta dataclasses, e
    # `dataclasses._is_type` resolve anotacao via `sys.modules[cls.__module__]`.
    # Com o modulo fora do registro isso vira None e estoura AttributeError.
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    falhas = []

    def check(descricao: str, condicao: bool) -> None:
        print(f"{'ok   ' if condicao else 'FALHA'} {descricao}")
        if not condicao:
            falhas.append(descricao)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        srv = carrega(tmp)

        rota = "/biology/mit/esm2/forward"
        a = {"sequence": "MKTAYIAK", "format": "pdb", "num_models": 1}
        b = {"num_models": 1, "format": "pdb", "sequence": "MKTAYIAK"}  # mesma coisa, outra ordem
        c = {"sequence": "MKTAYIAL", "format": "pdb", "num_models": 1}  # um aminoacido diferente

        ka, kb, kc = (srv._cache_key(rota, p) for p in (a, b, c))
        check("ordem das chaves do payload nao muda a chave do cache", ka == kb)
        check("sequencia diferente gera chave diferente", ka != kc)
        check("rota diferente gera chave diferente", ka != srv._cache_key("/outra", a))
        check("chave e um sha256 completo", len(ka) == 64)

        check("cache vazio devolve None", srv._cache_read(ka) is None)

        resposta = {"pdbs": ["ATOM      1  N   MET A   1"], "plddt": 87.3}
        srv._cache_write(ka, rota, resposta)
        check("ida e volta preserva a resposta", srv._cache_read(ka) == resposta)
        check("payload reordenado acerta a mesma entrada", srv._cache_read(kb) == resposta)
        check("payload diferente continua sem entrada", srv._cache_read(kc) is None)

        stats = srv._cache_stats()
        check("nim_config reporta 1 entrada", stats.get("entradas") == 1)
        check("nim_config reporta o cache como ativo", stats.get("ativo") is True)

        envelope = (tmp / "cache" / f"{ka}.json").read_text(encoding="utf-8")
        check("envelope guarda quando foi gravado", "gravado_em" in envelope)
        check("envelope guarda a rota de origem", rota in envelope)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        srv = carrega(tmp, cache_ligado=False)
        k = srv._cache_key("/x", {"s": 1})
        srv._cache_write(k, "/x", {"r": 1})
        check("BIO_STACK_NIM_CACHE=0 nao grava", not (tmp / "cache").exists())
        check("BIO_STACK_NIM_CACHE=0 nao le", srv._cache_read(k) is None)
        check("nim_config mostra o cache desligado", srv._cache_stats() == {"ativo": False})

    total = 15
    print(f"\n{total - len(falhas)}/{total} passaram")
    return 1 if falhas else 0


if __name__ == "__main__":
    sys.exit(main())
