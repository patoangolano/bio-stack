#!/usr/bin/env bash
# bio-stack — Claude Science no WSL (Ubuntu-24.04), como o usuario.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/07_wsl_claude_science.sh
#
# Pre-requisito: 04b_wsl_base.sh, que poe bubblewrap, socat e o uv em
# /usr/local/bin — sem os tres a sandbox do daemon nao roda os conectores.
#
# Depois deste script falta um passo que so voce faz, o login da conta no
# navegador:
#
#   wsl -d Ubuntu-24.04 -- bash -lc "claude-science url"
#
# imprime um link de uso unico, valido por ~3 min.
#
# O instalador vem de claude.ai, vai para arquivo antes de rodar, e o sha256
# dele fica no log. Ele mesmo confere o sha256 do binario contra o
# manifest.json da release (so checksum, sem assinatura) e instala em
# ~/.local/bin, sem sudo e sem mexer no ~/.profile.
#
# Quem opera esta instancia e servers/claude_science_mcp (CLAUDE_SCIENCE_TARGET=wsl
# no .mcp.json). Porta 8765 e nao 8000: o app do Windows usa a 8000, e o WSL2
# espelha as portas de localhost da distro para o Windows.
#
# Idempotente: so reinicia o daemon se ele estiver parado ou se o config mudou.
set -euo pipefail

INSTALADOR=https://claude.ai/install-claude-science.sh
PORTA=8765
CFG="$HOME/.claude-science/config.toml"

for dep in bwrap socat /usr/local/bin/uv; do
    command -v "$dep" >/dev/null || {
        echo "falta $dep; rode antes, como root:" >&2
        echo "  wsl -d Ubuntu-24.04 -u root -- bash /mnt/c/Users/Manec/bio-stack/bootstrap/04b_wsl_base.sh" >&2
        exit 1
    }
done

# --- 1. CLI ----------------------------------------------------------------------
if command -v claude-science >/dev/null; then
    echo "claude-science ja instalado: $(claude-science --version)   (atualizar: claude-science update)"
else
    inst="$(mktemp)"
    curl -fsSL --retry 3 -o "$inst" "$INSTALADOR"
    echo "instalador: sha256 $(sha256sum "$inst" | cut -d' ' -f1)"
    bash "$inst" 2>&1 | tail -3
    rm -f "$inst"
    export PATH="$HOME/.local/bin:$PATH"
    claude-science --version
fi

# --- 2. caminhos da sandbox para os conectores locais (servers/*) -------------------
mudou=0
mkdir -p "$HOME/.claude-science" "$HOME/.cache/uv"
if [ -f "$CFG" ] && grep -q '^\[sandbox\]' "$CFG"; then
    echo "config.toml ja tem [sandbox]"
else
    cat >> "$CFG" <<'TOML'
# Conectores locais rodam numa sandbox que nao enxerga quase nada do home.
# Libera so o necessario para os servidores da bio-stack, e nada alem:
#   leitura  -> o codigo dos servidores
#   gravacao -> runs/ (saidas em runs/nim e o provenance.jsonl) e o cache do uv
# data/raw/ fica de fora de proposito: e somente-leitura na politica do projeto.
[sandbox]
user_read_paths  = ["/mnt/c/Users/Manec/bio-stack/servers"]
user_write_paths = ["/mnt/c/Users/Manec/bio-stack/runs", "/home/manec/.cache/uv"]
TOML
    mudou=1
    echo "config.toml: bloco [sandbox] gravado"
fi
# Dominios somados a allowlist embutida da sandbox (a lista e aditiva; nao trava
# nada). A embutida ja cobre NCBI, Ensembl, EBI, UniProt e a API REST do UCSC;
# faltam os hosts de ARQUIVO. Sem eles cada acesso vira cartao de aprovacao na
# interface web — em 2026-10-16 uma pesquisa ficou 70 min parada esperando
# clique para o hgdownload, e morreu quando a distro desligou.
#   hgdownload.*.ucsc.edu   cadeias de liftOver, tabelas do Genome Browser
#   genome.animalgenome.org Animal QTLdb: mapas de ligacao de galinha
#   zenodo, dryad, figshare tabelas suplementares (mapas de recombinacao)
if grep -q '^\[sandbox\.network\]' "$CFG"; then
    echo "config.toml ja tem [sandbox.network]"
else
    cat >> "$CFG" <<'TOML'

[sandbox.network]
allowed_domains = [
  "hgdownload.soe.ucsc.edu",
  "hgdownload.gi.ucsc.edu",
  "genome.ucsc.edu",
  "genome.animalgenome.org",
  "zenodo.org",
  "datadryad.org",
  "figshare.com",
  "ndownloader.figshare.com",
]
TOML
    mudou=1
    echo "config.toml: bloco [sandbox.network] gravado"
fi
python3 -c "import tomllib, sys; tomllib.load(open(sys.argv[1], 'rb')); print('config.toml valido')" "$CFG"

# --- 3. daemon ---------------------------------------------------------------------
# "De pe" e de pe NA $PORTA. Quem sobe o daemon por outro caminho cai na 8000
# padrao: em 2026-10-23 o `claude-science url` do login o subiu la, e o
# wslrelay passou a segurar a 8000 do Windows, que e do app Windows.
rodando=$(claude-science status | python3 -c "
import json, sys
s = json.load(sys.stdin)
print(int(bool(s.get('running')) and s.get('port') == int(sys.argv[1])))" "$PORTA")
if [ "$rodando" = 1 ] && [ "$mudou" = 0 ]; then
    echo "daemon ja de pe na $PORTA"
else
    claude-science stop >/dev/null 2>&1 || true
    sleep 3
    claude-science serve --port "$PORTA" --detached --no-browser 2>&1 | grep -E 'daemon started|rror'
    sleep 6
fi
claude-science status | python3 -c "
import json, sys
s = json.load(sys.stdin); d = s.get('daemon') or {}
print('running', s['running'], '| porta', s.get('port'), '| versao', s.get('version'), '| sandbox', d.get('sandbox_active'))
"
echo
echo 'falta o login (uma vez): wsl -d Ubuntu-24.04 -- bash -lc "claude-science url"'
