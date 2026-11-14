#!/usr/bin/env python3
"""Stop: grava marcador de fim de turno em runs/provenance.jsonl.

Suporta protocolo Claude Code / Codex e Google Antigravity (JSON decision: allow).
"""

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


def _usable_env_path(var: str) -> str | None:
    val = os.environ.get(var, "").strip()
    if not val or "$" in val or "%" in val:
        return None
    return val


def main() -> None:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        payload = {}

    is_agy = isinstance(payload, dict) and ("hookEventName" in payload or "toolCall" in payload)

    if isinstance(payload, dict) and payload.get("stop_hook_active"):
        if is_agy:
            print(json.dumps({"decision": "allow"}))
        sys.exit(0)

    root = (
        _usable_env_path("CLAUDE_PROJECT_DIR")
        or _usable_env_path("BIO_STACK_ROOT")
        or str(Path(__file__).resolve().parents[2])
    )
    log_path = Path(root) / "runs" / "provenance.jsonl"
    if not log_path.exists():
        if is_agy:
            print(json.dumps({"decision": "allow"}))
        sys.exit(0)

    entry = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "session_id": payload.get("session_id") or payload.get("conversationId"),
        "event": "turn_end",
    }
    try:
        with log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry) + "\n")
    except Exception:
        pass

    if is_agy:
        print(json.dumps({"decision": "allow"}))
    sys.exit(0)


if __name__ == "__main__":
    main()
