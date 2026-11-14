#!/usr/bin/env bash
# bio-stack — aplicativos graficos de Windows (3D Slicer 5.12.4 + SlicerMorph,
# Fiji + BoneJ 7.2.2, Blender + blender_mcp.py) e ambiente R completo (Rtools45
# + r/renv com 237 pacotes).
#
# Idempotente: verifica cada componente antes de baixar ou instalar.

set -u

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
LOG="${1:-$RAIZ/runs/bootstrap/10_windows_apps_r.log}"
mkdir -p "$(dirname "$LOG")"
: > "$LOG"

registrar() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }

registrar "=== 1. Rtools45 + Blender via winget ==="
if [ -x "/c/rtools45/usr/bin/make.exe" ]; then
    registrar "OK   Rtools45 ja instalado em C:/rtools45"
else
    registrar "Instalando RProject.Rtools via winget..."
    winget.exe install --id RProject.Rtools --exact --accept-package-agreements --accept-source-agreements --disable-interactivity >>"$LOG" 2>&1 || true
    if [ -x "/c/rtools45/usr/bin/make.exe" ]; then
        registrar "OK   Rtools45 instalado"
    else
        registrar "FALHA Rtools45 nao encontrado em /c/rtools45"
    fi
fi

BLENDER_EXE=""
for b in "/c/Program Files/Blender Foundation"/Blender*/blender.exe; do
    [ -x "$b" ] && BLENDER_EXE="$b"
done
if [ -n "$BLENDER_EXE" ]; then
    registrar "OK   Blender ja instalado: $BLENDER_EXE"
else
    registrar "Instalando BlenderFoundation.Blender via winget..."
    winget.exe install --id BlenderFoundation.Blender --exact --accept-package-agreements --accept-source-agreements --disable-interactivity >>"$LOG" 2>&1 || true
    for b in "/c/Program Files/Blender Foundation"/Blender*/blender.exe; do
        [ -x "$b" ] && BLENDER_EXE="$b"
    done
    if [ -n "$BLENDER_EXE" ]; then
        registrar "OK   Blender instalado: $BLENDER_EXE"
    else
        registrar "FALHA Blender nao encontrado apos winget"
    fi
fi

if [ -n "$BLENDER_EXE" ]; then
    registrar "Configurando addon blender_mcp no Blender..."
    "$BLENDER_EXE" --background --python-expr "import bpy; p = bpy.utils.user_resource('SCRIPTS', path='addons', create=True); print('ADDONS_DIR=' + str(p))" >>"$LOG" 2>&1 || true
    /c/Users/user/.local/bin/uvx.exe mcp-for-blender install-addon >>"$LOG" 2>&1 || true
    "$BLENDER_EXE" --background --python-expr "
import bpy, addon_utils
bpy.ops.preferences.addon_enable(module='blender_mcp')
bpy.ops.wm.save_userpref()
hab = {m.__name__ for m in addon_utils.modules() if addon_utils.check(m.__name__)[1]}
print('BLENDER_MCP_HABILITADO=' + str('blender_mcp' in hab))
" >>"$LOG" 2>&1 || true
    registrar "OK   Addon blender_mcp configurado"
fi

registrar "=== 2. 3D Slicer 5.12.4 + SlicerMorph ==="
SLICER_EXE="/c/Users/user/AppData/Local/slicer.org/3D Slicer 5.12.4/Slicer.exe"
DL_SLICER="/c/Users/user/Apps/_dl-slicer"
mkdir -p "$DL_SLICER"

if [ -x "$SLICER_EXE" ]; then
    registrar "OK   3D Slicer 5.12.4 ja instalado"
else
    registrar "Baixando 3D Slicer 5.12.4 (bitstream 6aa1db04ce9de556d30112bb)..."
    curl -sL --retry 3 -o "$DL_SLICER/Slicer-5.12.4-win-amd64.exe" \
        "https://download.slicer.org/bitstream/6aa1db04ce9de556d30112bb" >>"$LOG" 2>&1
    if "$RAIZ/.venv/Scripts/python.exe" -I -c "
import hashlib, sys
esperado = '5ba320cb67f0acdaacf0a31380e9cf3b9f46d05f57a53da7e04ebcd9490251ecbe998bfefdbdab1f747e038653177868cdb0ad30986473b9d68efba0f4c6045c'
h = hashlib.sha512()
with open(r'C:\Users\user\Apps\_dl-slicer\Slicer-5.12.4-win-amd64.exe', 'rb') as f:
    for b in iter(lambda: f.read(1 << 20), b''):
        h.update(b)
sys.exit(0 if h.hexdigest() == esperado else 1)
"; then
        registrar "OK   sha512 do instalador do 3D Slicer confere; executando /S..."
        "$DL_SLICER/Slicer-5.12.4-win-amd64.exe" /S >>"$LOG" 2>&1 || true
        sleep 8
        if [ ! -x "$SLICER_EXE" ] && [ -x "/c/ProgramData/slicer.org/3D Slicer 5.12.4/Slicer.exe" ]; then
            mkdir -p "/c/Users/user/AppData/Local/slicer.org"
            cmd.exe //c "mklink /J \"C:\\Users\\user\\AppData\\Local\\slicer.org\\3D Slicer 5.12.4\" \"C:\\ProgramData\\slicer.org\\3D Slicer 5.12.4\"" >>"$LOG" 2>&1 || true
        fi
    else
        registrar "FALHA sha512 do instalador do 3D Slicer divergiu"
    fi
fi

if [ -x "$SLICER_EXE" ]; then
    registrar "Baixando e instalando SlicerMorph e dependencias de forma sincrona..."
    DL_EXT="/c/Users/user/Apps/_dl-slicerext"
    mkdir -p "$DL_EXT"
    for par in \
        "SurfaceMarkup=6aa1fd8fce9de556d3012cd4" \
        "Sandbox=6aa1e8b7ce9de556d3011d62" \
        "SegmentEditorExtraEffects=6aa1f67fce9de556d30127e2" \
        "SlicerMorph=6aa1fde7ce9de556d3012d2d"; do
        nome="${par%%=*}"
        id="${par#*=}"
        if [ ! -s "$DL_EXT/$nome.zip" ]; then
            curl -sL --retry 3 -o "$DL_EXT/$nome.zip" \
                "https://slicer-packages.kitware.com/api/v1/item/$id/download"
        fi
    done

    cat > "$DL_EXT/instalar_ext.py" <<'EOF'
import os, slicer

gestor = slicer.app.extensionsManagerModel()
gestor.interactive = False

base = r"C:\Users\user\Apps\_dl-slicerext"
ordem = ["SurfaceMarkup", "Sandbox", "SegmentEditorExtraEffects", "SlicerMorph"]

saida = []
for nome in ordem:
    caminho = os.path.join(base, nome + ".zip")
    if gestor.isExtensionInstalled(nome):
        saida.append(f"JA_INSTALADO {nome}")
        continue
    if not os.path.isfile(caminho):
        saida.append(f"SEM_PACOTE {nome}")
        continue
    ok = gestor.installExtension(caminho, False, False)
    saida.append(f"{'OK' if ok else 'FALHA'} {nome}")

saida.append("INSTALADAS=" + ",".join(sorted(gestor.installedExtensions)))
with open(r"C:\Users\user\Apps\_dl-slicerext\ext-resultado.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(saida) + "\n")

slicer.util.exit(0)
EOF
    rm -f "$DL_EXT/ext-resultado.txt"
    "$SLICER_EXE" --no-main-window --no-splash --python-script "C:/Users/user/Apps/_dl-slicerext/instalar_ext.py" >/dev/null 2>&1 || true
    cat "$DL_EXT/ext-resultado.txt" 2>/dev/null | tee -a "$LOG"
    rm -rf "$DL_SLICER"
fi

registrar "=== 3. Fiji (ImageJ) + BoneJ 7.2.2 ==="
FIJI_EXE="/c/Users/user/Apps/Fiji.app/ImageJ-win64.exe"
DL_FIJI="/c/Users/user/Apps/_dl-fiji"
if [ -x "$FIJI_EXE" ]; then
    registrar "OK   Fiji.app ja instalado"
else
    mkdir -p "$DL_FIJI"
    registrar "Baixando Fiji portable (fiji-stable-win64-jdk.zip)..."
    curl -sL --retry 3 -o "$DL_FIJI/fiji-stable-win64-jdk.zip" \
        "https://downloads.imagej.net/fiji/stable/fiji-stable-win64-jdk.zip" >>"$LOG" 2>&1
    registrar "Extraindo Fiji.app em C:/Users/user/Apps..."
    "$RAIZ/.venv/Scripts/python.exe" -I -c "
import zipfile
zipfile.ZipFile(r'C:\Users\user\Apps\_dl-fiji\fiji-stable-win64-jdk.zip').extractall(r'C:\Users\user\Apps')
" >>"$LOG" 2>&1
    rm -rf "$DL_FIJI"
fi

if [ -x "$FIJI_EXE" ]; then
    if ls /c/Users/user/Apps/Fiji.app/jars/bonej*.jar >/dev/null 2>&1; then
        registrar "OK   BoneJ ja instalado em Fiji.app/jars"
    else
        registrar "Habilitando update site BoneJ e atualizando Fiji..."
        (
            cd /c/Users/user/Apps/Fiji.app &&
            ./ImageJ-win64.exe --update add-update-site BoneJ https://sites.imagej.net/BoneJ/ >>"$LOG" 2>&1 &&
            ./ImageJ-win64.exe --update update >>"$LOG" 2>&1
        ) || true
        if ls /c/Users/user/Apps/Fiji.app/jars/bonej*.jar >/dev/null 2>&1; then
            registrar "OK   BoneJ instalado em Fiji.app/jars"
        else
            registrar "FALHA BoneJ nao apareceu em Fiji.app/jars"
        fi
    fi
fi

registrar "=== 4. Biblioteca R (r/renv - 237 pacotes) ==="
RSCRIPT="/c/Program Files/R/R-4.6.2/bin/Rscript.exe"
if [ ! -x "$RSCRIPT" ]; then
    for r in "/c/Program Files/R"/R-*/bin/Rscript.exe; do
        [ -x "$r" ] && RSCRIPT="$r"
    done
fi
if [ -x "$RSCRIPT" ]; then
    export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
    registrar "Restaurando r/renv via $RSCRIPT r/instalar.R..."
    if "$RSCRIPT" "$RAIZ/r/instalar.R" >>"$LOG" 2>&1; then
        registrar "OK   r/instalar.R concluido"
    else
        registrar "FALHA r/instalar.R saiu com codigo $?"
    fi
else
    registrar "FALHA Rscript.exe nao encontrado"
fi

registrar "=== fim 10_windows_apps_r.sh ==="
