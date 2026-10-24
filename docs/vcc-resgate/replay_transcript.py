"""Reconstroi os arquivos perdidos da sessao 07a4432b por replay do transcript.

Executa, em ordem de timestamp e intercalando o principal com os subagentes:
  * Write  -> grava o conteudo inteiro
  * Edit   -> aplica old_string -> new_string (falha registrada, nunca silenciosa)
  * Bash com `python -I - <<'PYEOF' ... PYEOF` que toca um arquivo-alvo -> executa
    o corpo com cwd na caixa de areia (so se o original NAO terminou em erro)
  * Bash que roda um patch_*.py do scratchpad -> roda a copia resgatada em
    docs/vcc-resgate/scripts/

E confere: cada Read COMPLETO de um arquivo-alvo no transcript e um snapshot; no
instante dele, o conteudo reconstruido tem de ser igual.

Uso: python -I replay.py <pasta_de_saida>
"""

from __future__ import annotations

import json
import pathlib
import re
import subprocess
import sys
from dataclasses import dataclass, field

BASE = pathlib.Path(r"C:\Users\Manec\.claude\projects\C--Users-Manec-bio-stack")
SESSAO = "07a4432b-599e-4679-83d3-c71280471116"
RAIZ = pathlib.Path(r"C:\Users\Manec\bio-stack")
RESGATE = RAIZ / "docs" / "vcc-resgate" / "scripts"
PY = str(RAIZ / ".venv" / "Scripts" / "python.exe")

PREFIXOS_ALVO = (
    "scripts/vcc/",
    "scripts/colab/",
    "docs/VCC-2026-",
    "docs/AVALIADOR-PROXY.md",
    "servers/vcc_mcp/",
    "envs/avaliador-proxy/",
    "runs/vcc/submissao/LEIA-ME.md",
)


def normaliza(caminho: str) -> str:
    c = caminho.replace("\\", "/")
    for prefixo in (
        "C:/Users/Manec/bio-stack/",
        "c:/Users/Manec/bio-stack/",
        "/c/Users/Manec/bio-stack/",
        "/mnt/c/Users/Manec/bio-stack/",
    ):
        if c.startswith(prefixo):
            return c[len(prefixo) :]
    return c


def e_alvo(rel: str) -> bool:
    return rel.startswith(PREFIXOS_ALVO)


@dataclass(order=True)
class Evento:
    ts: str
    seq: int
    origem: str = field(compare=False)
    tipo: str = field(compare=False)
    dados: dict = field(compare=False)


def texto_do_resultado(conteudo) -> str:
    if isinstance(conteudo, str):
        return conteudo
    if isinstance(conteudo, list):
        return "".join(b.get("text", "") for b in conteudo if isinstance(b, dict))
    return ""


def coletar() -> list[Evento]:
    fontes = [("principal", BASE / f"{SESSAO}.jsonl")]
    fontes += [(p.stem, p) for p in sorted((BASE / SESSAO / "subagents").glob("*.jsonl"))]
    usos: dict[str, tuple[str, str, int, str, dict]] = {}
    resultados: dict[str, tuple[str, bool]] = {}
    ts_resultado: dict[str, str] = {}
    seq = 0
    for origem, caminho in fontes:
        for linha in caminho.open(encoding="utf-8"):
            try:
                d = json.loads(linha)
            except json.JSONDecodeError:
                continue
            c = (d.get("message") or {}).get("content")
            if not isinstance(c, list):
                continue
            for b in c:
                if not isinstance(b, dict):
                    continue
                seq += 1
                if b.get("type") == "tool_use":
                    # Cada tool_use aparece DUAS vezes no JSONL, a segunda sem
                    # timestamp. Sobrescrever com ela poe o evento em "" e o
                    # ordena antes de tudo -- inclusive antes do Write do arquivo.
                    ts = d.get("timestamp") or ""
                    if b["id"] not in usos or (ts and not usos[b["id"]][0]):
                        usos[b["id"]] = (ts, origem, seq, b.get("name", ""), b.get("input") or {})
                elif b.get("type") == "tool_result":
                    tid = b.get("tool_use_id", "")
                    if tid not in resultados:
                        resultados[tid] = (texto_do_resultado(b.get("content")), bool(b.get("is_error")))
                    if d.get("timestamp") and not ts_resultado.get(tid):
                        ts_resultado[tid] = d["timestamp"]

    eventos: list[Evento] = []
    for tid, (ts, origem, s, nome, entrada) in usos.items():
        texto, erro = resultados.get(tid, ("", False))
        base = {"entrada": entrada, "erro": erro, "resultado": texto}
        if nome in ("Write", "Edit"):
            rel = normaliza(str(entrada.get("file_path", "")))
            if e_alvo(rel):
                eventos.append(Evento(ts, s, origem, nome, {**base, "rel": rel}))
        elif nome == "Read":
            rel = normaliza(str(entrada.get("file_path", "")))
            if e_alvo(rel) and not entrada.get("offset") and not entrada.get("limit") and not erro:
                # O snapshot vale no instante do RESULTADO, nao do pedido.
                eventos.append(Evento(ts_resultado.get(tid, ts), s, origem, "Read", {**base, "rel": rel}))
        elif nome == "Bash":
            cmd = str(entrada.get("command", ""))
            eventos.append(Evento(ts, s, origem, "Bash", {**base, "cmd": cmd}))
    eventos.sort()
    return eventos


LINHA_LIDA = re.compile(r"^\s*(\d+)\t(.*)$")


def snapshot(texto: str) -> str | None:
    linhas = []
    esperado = 1
    for bruta in texto.split("\n"):
        m = LINHA_LIDA.match(bruta)
        if not m:
            if linhas:
                break
            continue
        n = int(m.group(1))
        if n != esperado:
            return None
        linhas.append(m.group(2))
        esperado += 1
    if not linhas:
        return None
    return "\n".join(linhas)


def redirecionar(codigo: str, saida: pathlib.Path) -> str:
    alvo = str(saida).replace("\\", "/")
    for forma in ("C:/Users/Manec/bio-stack", "/c/Users/Manec/bio-stack"):
        codigo = codigo.replace(forma, alvo)
    return codigo


HEREDOC = re.compile(r"<<'?(PYEOF|PY|EOF)'?\n(.*?)\n\1\s*$", re.S | re.M)
PATCH_SCRIPT = re.compile(r"scratchpad/(patch_[a-z0-9_]+\.py)")


def main() -> int:
    sys.stdout.reconfigure(encoding="utf-8")  # pyright: ignore[reportAttributeAccessIssue]
    saida = pathlib.Path(sys.argv[1])
    saida.mkdir(parents=True, exist_ok=True)
    eventos = coletar()
    estado: dict[str, str] = {}
    log: list[str] = []
    stats: dict[str, dict[str, int]] = {}

    def st(rel: str) -> dict[str, int]:
        return stats.setdefault(rel, {"write": 0, "edit_ok": 0, "edit_falhou": 0, "bash": 0, "snap_ok": 0, "snap_dif": 0})

    def gravar_caixa() -> None:
        for rel, conteudo in estado.items():
            p = saida / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(conteudo, encoding="utf-8", newline="\n")

    def ler_caixa() -> None:
        for rel in list(estado):
            p = saida / rel
            if p.exists():
                estado[rel] = p.read_text(encoding="utf-8")

    for ev in eventos:
        d = ev.dados
        if ev.tipo == "Write":
            if d["erro"]:
                continue
            estado[d["rel"]] = str(d["entrada"].get("content", ""))
            st(d["rel"])["write"] += 1
        elif ev.tipo == "Edit":
            rel = d["rel"]
            if d["erro"]:
                continue
            velho = str(d["entrada"].get("old_string", ""))
            novo = str(d["entrada"].get("new_string", ""))
            tudo = bool(d["entrada"].get("replace_all"))
            atual = estado.get(rel)
            if atual is None or velho not in atual or (not tudo and atual.count(velho) != 1):
                st(rel)["edit_falhou"] += 1
                n = 0 if atual is None else atual.count(velho)
                log.append(f"{ev.ts[:19]} EDIT FALHOU {rel} (ocorrencias {n}): {velho[:80]!r}")
                continue
            estado[rel] = atual.replace(velho, novo) if tudo else atual.replace(velho, novo, 1)
            st(rel)["edit_ok"] += 1
        elif ev.tipo == "Read":
            rel = d["rel"]
            snap = snapshot(d["resultado"])
            if snap is None or rel not in estado:
                continue
            if estado[rel].rstrip("\n") == snap.rstrip("\n"):
                st(rel)["snap_ok"] += 1
            else:
                st(rel)["snap_dif"] += 1
                a = estado[rel].rstrip("\n").split("\n")
                b = snap.rstrip("\n").split("\n")
                primeira = next((i for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y), min(len(a), len(b)))
                log.append(f"{ev.ts[:19]} SNAPSHOT DIFERE {rel}: linhas {len(a)} vs {len(b)}, 1a dif na {primeira + 1}")
                # O snapshot e a verdade daquele instante: adota.
                estado[rel] = snap + "\n"
        elif ev.tipo == "Bash":
            if d["erro"]:
                continue
            cmd = d["cmd"]
            tocados = [r for r in {normaliza(x) for x in re.findall(r"[\w./\\:-]*(?:scripts/vcc|docs/VCC-2026-|scripts/colab)[\w./-]*", cmd)} if e_alvo(r)]
            corpo = None
            m = HEREDOC.search(cmd)
            if m and "python" in cmd.split("<<")[0] and "write_text" in m.group(2):
                corpo = m.group(2)
            script = None
            mp = PATCH_SCRIPT.search(cmd)
            if mp and ("python" in cmd) and "cat > " not in cmd.split(mp.group(1))[0][-40:]:
                copia = RESGATE / mp.group(1)
                if copia.exists():
                    script = copia
            if corpo is None and script is None:
                continue
            if corpo is not None and not any(t in corpo for t in ("scripts/vcc", "docs/VCC-2026", "scripts/colab")):
                continue
            perigoso = [p for p in ("subprocess", "os.system", "shutil.rmtree", "unlink", "rmdir") if p in (corpo or "")]
            if perigoso:
                log.append(f"{ev.ts[:19]} BASH PULADO (contem {perigoso}): {cmd[:90]!r}")
                continue
            gravar_caixa()
            ambiente = {**__import__("os").environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
            if corpo is not None:
                corpo = redirecionar(corpo, saida)
                r = subprocess.run([PY, "-I", "-"], input=corpo, text=True, cwd=saida, capture_output=True, encoding="utf-8", env=ambiente)
            else:
                # Os patch_*.py resgatados apontam para o repo REAL por caminho
                # absoluto: roda uma copia redirecionada para a caixa de areia.
                assert script is not None
                copia_local = saida.parent / f"_replay_{script.name}"
                copia_local.write_text(redirecionar(script.read_text(encoding="utf-8"), saida), encoding="utf-8")
                r = subprocess.run([PY, "-I", str(copia_local)], text=True, cwd=saida, capture_output=True, encoding="utf-8", env=ambiente)
            antes = dict(estado)
            # Arquivos novos criados pelo patch tambem entram.
            for p in saida.rglob("*"):
                if p.is_file():
                    rel = p.relative_to(saida).as_posix()
                    if e_alvo(rel) and rel not in estado:
                        estado[rel] = p.read_text(encoding="utf-8")
            ler_caixa()
            mudou = [rel for rel in estado if antes.get(rel) != estado[rel]]
            for rel in mudou:
                st(rel)["bash"] += 1
            rotulo = script.name if script else "heredoc"
            if r.returncode != 0:
                log.append(f"{ev.ts[:19]} BASH {rotulo} SAIU {r.returncode}: {(r.stderr or r.stdout)[-160:]!r}")
            elif not mudou:
                log.append(f"{ev.ts[:19]} BASH {rotulo} rodou e nao mudou nada (tocava {tocados[:3]})")

    gravar_caixa()
    print(f"{'arquivo':58s} wr  ed+ ed-  sh  snap+ snap-  linhas")
    for rel in sorted(estado):
        s = st(rel)
        print(f"{rel:58s} {s['write']:2d} {s['edit_ok']:4d} {s['edit_falhou']:3d} {s['bash']:3d} {s['snap_ok']:5d} {s['snap_dif']:5d}  {estado[rel].count(chr(10)):6d}")
    print("\n--- ocorrencias ---")
    for linha in log:
        print(linha)
    return 0


if __name__ == "__main__":
    sys.exit(main())
