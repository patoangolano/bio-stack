<#
.SYNOPSIS
    Cria o ambiente Python da bio-stack com uv.

.DESCRIPTION
    Monta um .venv em bio-stack/.venv com tres camadas:

      core  - biopython, pandas, scanpy, anndata, pydeseq2, gget, httpx  (sempre)
      gpu   - PyTorch CUDA 12.8 (Blackwell sm_120) + transformers + fair-esm (--Gpu)
      dev   - ruff, basedpyright, pytest  (--Dev)

    Nota honesta sobre GPU no Windows:
      - PyTorch CUDA funciona nativo no Windows. OK.
      - RAPIDS / rapids-singlecell / cuDF NAO tem build nativo para Windows.
        Precisa de WSL2 (Ubuntu) ou Linux. O script avisa e nao tenta instalar.
      - NVIDIA Parabricks exige >= 16 GB de VRAM e Linux/Docker. Fora desta maquina;
        use DGX Cloud / cluster da universidade ou os endpoints hospedados.

.EXAMPLE
    .\02_python_env.ps1 -Gpu -Dev
#>

[CmdletBinding()]
param(
    [switch]$Gpu,
    [switch]$Dev,
    [string]$PythonVersion = '3.11',
    # Indice de wheels do PyTorch. cu128 cobre Blackwell (RTX 50xx, sm_120).
    [string]$TorchIndex = 'https://download.pytorch.org/whl/cu128'
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # .../bio-stack
$venv = Join-Path $root '.venv'

function Write-Step($m) { Write-Host "`n=== $m ===" -ForegroundColor Cyan }

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Host 'uv nao encontrado. Rode antes: .\00_bootstrap_windows.ps1' -ForegroundColor Red
    exit 1
}

Write-Step "Criando venv (Python $PythonVersion)"
uv venv --python $PythonVersion $venv
if ($LASTEXITCODE -ne 0) { throw 'uv venv falhou' }

$py = Join-Path $venv 'Scripts\python.exe'
$env:VIRTUAL_ENV = $venv

Write-Step 'Camada core'
uv pip install --python $py -r (Join-Path $PSScriptRoot 'requirements-core.txt')
if ($LASTEXITCODE -ne 0) { throw 'instalacao do core falhou' }

if ($Gpu) {
    Write-Step "Camada GPU (PyTorch via $TorchIndex)"
    # PyTorch entra separado: precisa do indice proprio da NVIDIA, nao do PyPI.
    uv pip install --python $py --index-url $TorchIndex torch torchvision
    if ($LASTEXITCODE -ne 0) {
        Write-Host '  [aviso] PyTorch CUDA falhou. Tentando indice cu126 como alternativa.' -ForegroundColor Yellow
        uv pip install --python $py --index-url 'https://download.pytorch.org/whl/cu126' torch torchvision
    }
    uv pip install --python $py -r (Join-Path $PSScriptRoot 'requirements-gpu.txt')

    Write-Step 'Verificando CUDA'
    & $py -c @'
import torch
print(f"  torch          {torch.__version__}")
print(f"  cuda disponivel {torch.cuda.is_available()}")
if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    gb = p.total_memory / 1024**3
    print(f"  gpu            {p.name}")
    print(f"  capability     sm_{p.major}{p.minor}")
    print(f"  vram           {gb:.1f} GB")
    if gb < 16:
        print("  [aviso] < 16 GB: use os NIMs hospedados para predicao de estrutura.")
else:
    print("  [aviso] build de CPU instalado - reinstale com o indice cu128.")
'@
}

if ($Dev) {
    Write-Step 'Camada dev (linters / type checker)'
    uv pip install --python $py ruff basedpyright pytest pytest-cov
}

Write-Step 'RAPIDS / GPU single-cell'
Write-Host @'
  rapids-singlecell e cuDF nao tem build nativo para Windows.
  Para usar a GPU em scRNA-seq (10-50x mais rapido que scanpy em CPU):

      wsl --install -d Ubuntu-24.04
      # dentro do WSL:
      curl -LsSf https://astral.sh/uv/install.sh | sh
      uv venv --python 3.11 && source .venv/bin/activate
      uv pip install --extra-index-url=https://pypi.nvidia.com \
          cudf-cu12 cuml-cu12 cugraph-cu12 rapids-singlecell

  Sem WSL, scanpy em CPU resolve bem datasets ate ~50k celulas.
'@ -ForegroundColor Gray

Write-Host "`nAmbiente pronto: $venv" -ForegroundColor Green
Write-Host "Ativar:  $venv\Scripts\Activate.ps1`n" -ForegroundColor Green
