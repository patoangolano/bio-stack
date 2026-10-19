#!/usr/bin/env bash
set -uo pipefail
ls -la ~/.taxonkit 2>/dev/null | head
echo "Falco peregrinus" | taxonkit name2taxid 2>&1 | head -3
