<#
.SYNOPSIS
    Camada 0 da bio-stack: instala as dependencias de base que faltam no Windows.

.DESCRIPTION
    Diagnostico feito nesta maquina (RTX 5070 Laptop, 8 GB VRAM, CUDA 13.4):
      git    AUSENTE  -> obrigatorio para /plugin marketplace add
      node   AUSENTE  -> obrigatorio para `npx skills add` (BioNeMo Agent Toolkit)
      uv     AUSENTE  -> obrigatorio para rodar servidores MCP em Python via uvx
      python AUSENTE  -> so existe o stub da Microsoft Store
      docker AUSENTE  -> opcional (so para auto-hospedar containers NIM)

    Este script instala tudo via winget, em escopo de usuario (sem admin quando possivel).
    Rode uma vez e reinicie o terminal.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\00_bootstrap_windows.ps1
    powershell -ExecutionPolicy Bypass -File .\00_bootstrap_windows.ps1 -IncludeDocker
#>

[CmdletBinding()]
param(
    # Docker Desktop + WSL2. So necessario se voce for auto-hospedar containers NIM.
    # Com 8 GB de VRAM a maioria dos NIMs de estrutura NAO cabe local - prefira a API hospedada.
    [switch]$IncludeDocker,

    # Pula a etapa de checagem de GPU.
    [switch]$SkipGpuCheck
)

$ErrorActionPreference = 'Stop'
$script:Failures = @()

function Write-Step($msg) { Write-Host "`n=== $msg ===" -ForegroundColor Cyan }
function Write-Ok($msg)   { Write-Host "  [ok]   $msg" -ForegroundColor Green }
function Write-Warn2($msg){ Write-Host "  [aviso] $msg" -ForegroundColor Yellow }
function Write-Err($msg)  { Write-Host "  [erro] $msg" -ForegroundColor Red }

function Test-Tool($name) {
    $cmd = Get-Command $name -ErrorAction SilentlyContinue
    if (-not $cmd) { return $false }

    # O Windows planta um stub de python.exe em WindowsApps que nao e um
    # interpretador: ele so abre a Microsoft Store. Get-Command encontra o stub,
    # entao sem este descarte o script pula a instalacao do Python de verdade.
    if ($name -like 'python*' -and $cmd.Source -like '*\WindowsApps\*') {
        return $false
    }
    return $true
}

function Install-IfMissing {
    param(
        [Parameter(Mandatory)][string]$Command,
        [Parameter(Mandatory)][string]$WingetId,
        [string]$Label = $null,
        [string[]]$ExtraArgs = @()
    )
    $Label = if ($Label) { $Label } else { $Command }

    if (Test-Tool $Command) {
        Write-Ok "$Label ja instalado -> $((Get-Command $Command).Source)"
        return
    }

    Write-Host "  instalando $Label ($WingetId)..." -ForegroundColor Gray
    # nao use $args aqui: e variavel automatica do PowerShell e o splatting fica ambiguo
    $wingetArgs = @('install', '--id', $WingetId, '--exact',
                    '--accept-source-agreements', '--accept-package-agreements',
                    '--disable-interactivity') + $ExtraArgs
    try {
        & winget @wingetArgs
        # winget devolve 0 para sucesso e -1978335189 para "ja instalado"
        if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne -1978335189) {
            throw "winget saiu com codigo $LASTEXITCODE"
        }
        Write-Ok "$Label instalado"
    } catch {
        Write-Err "falha ao instalar ${Label}: $_"
        $script:Failures += $Label
    }
}

# ---------------------------------------------------------------------------
Write-Step 'Pre-requisitos'

if (-not (Test-Tool 'winget')) {
    Write-Err 'winget nao encontrado. Atualize o "Instalador de Aplicativo" pela Microsoft Store e rode de novo.'
    exit 1
}
Write-Ok "winget disponivel"

# ---------------------------------------------------------------------------
Write-Step 'Ferramentas de base'

# Git: exigido pelo sistema de plugins do Claude Code (/plugin marketplace add clona repos).
Install-IfMissing -Command 'git' -WingetId 'Git.Git' -Label 'Git'

# Node LTS: exigido pelo `npx skills add` do BioNeMo Agent Toolkit.
Install-IfMissing -Command 'node' -WingetId 'OpenJS.NodeJS.LTS' -Label 'Node.js LTS'

# Python 3.11: melhor compatibilidade com o stack de bioinfo (scanpy, pydeseq2, biopython).
# 3.13 ainda quebra algumas wheels cientificas.
Install-IfMissing -Command 'python' -WingetId 'Python.Python.3.11' -Label 'Python 3.11'

# uv: roda servidores MCP Python isolados via `uvx`, sem poluir o sistema.
Install-IfMissing -Command 'uv' -WingetId 'astral-sh.uv' -Label 'uv'

if ($IncludeDocker) {
    Write-Step 'Docker Desktop (opcional - NIM auto-hospedado)'
    Install-IfMissing -Command 'docker' -WingetId 'Docker.DockerDesktop' -Label 'Docker Desktop'
    Write-Warn2 'Docker exige WSL2 + reinicio. Depois rode: wsl --install -d Ubuntu-24.04'
}

# ---------------------------------------------------------------------------
Write-Step 'Recarregando PATH nesta sessao'

$machinePath = [Environment]::GetEnvironmentVariable('Path', 'Machine')
$userPath    = [Environment]::GetEnvironmentVariable('Path', 'User')
$env:Path    = "$machinePath;$userPath"
Write-Ok 'PATH recarregado (abra um terminal novo se algo ainda nao aparecer)'

# ---------------------------------------------------------------------------
if (-not $SkipGpuCheck) {
    Write-Step 'GPU / CUDA'
    if (Test-Tool 'nvidia-smi') {
        $q = 'name,driver_version,memory.total,compute_cap'
        $raw = (& nvidia-smi --query-gpu=$q --format=csv,noheader) -join "`n"
        Write-Host "  $raw" -ForegroundColor Gray

        $vramMb = 0
        if ($raw -match '(\d+)\s*MiB') { $vramMb = [int]$Matches[1] }

        if ($vramMb -gt 0 -and $vramMb -lt 16000) {
            Write-Warn2 "VRAM detectada: ${vramMb} MiB."
            Write-Warn2 'Abaixo de 16 GB: NIMs de estrutura (OpenFold3, Boltz-2, AlphaFold2) NAO rodam local.'
            Write-Warn2 'Use os endpoints hospedados em health.api.nvidia.com (ver docs/NVIDIA.md).'
            Write-Warn2 'Local, cabem bem: ESM-2 650M, DNABERT-2, scanpy/rapids-singlecell, nvMolKit.'
        }
    } else {
        Write-Warn2 'nvidia-smi nao encontrado - driver NVIDIA ausente ou fora do PATH.'
    }
}

# ---------------------------------------------------------------------------
Write-Step 'Resultado'

$check = @('git','node','npm','npx','python','uv','uvx')
foreach ($t in $check) {
    if (Test-Tool $t) { Write-Ok $t } else { Write-Err "$t ainda ausente" }
}

if ($script:Failures.Count -gt 0) {
    Write-Host ''
    Write-Err "Falharam: $($script:Failures -join ', ')"
    Write-Host 'Instale manualmente e rode este script de novo.' -ForegroundColor Yellow
    exit 1
}

Write-Host @'

Camada 0 pronta.

Proximos passos:
  1. Feche e reabra o Claude Code (para ele enxergar git/node/uv no PATH).
  2. Configure as chaves:      .\bootstrap\01_setup_keys.ps1
  3. Crie o ambiente Python:   .\bootstrap\02_python_env.ps1
  4. Instale os plugins:       veja bootstrap\03_install_plugins.md

'@ -ForegroundColor Green
