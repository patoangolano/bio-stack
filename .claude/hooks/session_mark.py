# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
Marcadores de limite de sessao no caderno de laboratorio.

O `provenance.py` registra operacoes, mas nao registra onde uma sessao comeca,
termina, ou tem o contexto compactado. Isso importa por dois motivos concretos:

1. **Compactacao apaga o raciocinio, nao o registro.** Quando o contexto e
   compactado no meio de uma analise, o que o modelo lembra depois e um resumo.
   Se o log nao marcar o ponto, `/metodos` nao tem como saber que houve uma
   descontinuidade — e vai narrar como continua uma analise que foi reconstruida
   de memoria resumida.

2. **Sessao e a unidade natural de "uma sentada de trabalho".** Agrupar o log
   por sessao e o que separa "rodei 40 comandos" de "rodei 40 comandos em tres
   tentativas, das quais a segunda foi descartada".

Escreve no mesmo arquivo e no mesmo formato do provenance, com kind="marco",
para que qualquer leitor do log continue funcionando sem alteracao.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path


def _log_path() -> Path:
    raw = os.environ.get("BIO_STACK_PROVENANCE", "")
    if raw and "${" not in raw:
        return Path(raw)
    return Path(__file__).resolve().parents[2] / "runs" / "provenance.jsonl"


# Traducao para portugues so na leitura humana; a chave `event` guarda o nome
# original, que e o que um script deve filtrar.
RAIZ = Path(__file__).resolve().parents[2]


def _estado_git() -> dict[str, object]:
    """
    Commit do HEAD e se havia alteracao nao commitada.

    A distincao e o que separa reprodutivel de quase-reprodutivel: um resultado
    produzido com a arvore suja nao pode ser refeito a partir do SHA, porque o
    codigo que rodou nao esta em lugar nenhum do historico.
    """
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=RAIZ, capture_output=True, text=True, timeout=8, check=False,
        )
        if sha.returncode != 0:
            return {"git": "fora de repositorio"}
        sujo = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=RAIZ, capture_output=True, text=True, timeout=8, check=False,
        )
        pendentes = [ln for ln in sujo.stdout.splitlines() if ln.strip()]
        return {
            "commit": sha.stdout.strip()[:40],
            "arvore_limpa": not pendentes,
            "arquivos_pendentes": len(pendentes),
        }
    except (OSError, subprocess.SubprocessError):
        return {"git": "indisponivel"}


DESCRICAO = {
    "SessionStart": "sessao iniciada",
    "SessionEnd": "sessao encerrada",
    "PreCompact": "contexto compactado — o que vem depois deste ponto foi "
    "reconstruido a partir de resumo, nao do transcript completo",
    "PostCompact": "contexto retomado apos compactacao",
}


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    event = payload.get("hook_event_name") or ""
    if event not in DESCRICAO:
        return 0

    entry = {
        "ts": datetime.now(UTC).isoformat(timespec="seconds"),
        "session": payload.get("session_id"),
        "kind": "marco",
        "event": event,
        "nota": DESCRICAO[event],
    }

    # Estado do repositorio no limite da sessao. Aqui `git status` cabe: roda
    # duas vezes por sessao, nao uma vez por chamada de ferramenta. E e o que
    # qualifica o SHA que o provenance.py grava em cada linha — um commit
    # limpo significa "este codigo exato produziu este resultado"; um commit
    # sujo significa "havia alteracao nao commitada, o SHA e so uma dica".
    entry.update(_estado_git())

    # SessionEnd informa por que a sessao acabou (clear, logout, saida do
    # usuario). Util para distinguir "terminou o trabalho" de "caiu".
    if payload.get("reason"):
        entry["motivo"] = str(payload["reason"])[:100]
    # PreCompact distingue compactacao automatica (estourou a janela) de manual.
    if payload.get("trigger"):
        entry["gatilho"] = str(payload["trigger"])[:50]

    log = _log_path()
    try:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except OSError:
        pass  # um marcador perdido nunca justifica derrubar a sessao

    return 0


if __name__ == "__main__":
    sys.exit(main())
