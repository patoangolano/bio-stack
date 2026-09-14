<#
.SYNOPSIS
    Autentica no registry de containers da NVIDIA (nvcr.io).

.DESCRIPTION
    So e necessario para AUTO-HOSPEDAR containers NIM. Com 8 GB de VRAM nesta
    maquina, nenhum NIM de estrutura sobe local - isto e para quando voce tiver
    acesso ao cluster da universidade ou a DGX Cloud.

    A chave e lida da variavel de ambiente e passada por stdin: ela nunca
    aparece na linha de comando nem no historico do shell.

.EXAMPLE
    .\04_nvcr_login.ps1
#>
$ErrorActionPreference = 'Stop'

$key = [Environment]::GetEnvironmentVariable('NVIDIA_API_KEY', 'User')
if (-not $key) { $key = [Environment]::GetEnvironmentVariable('NGC_API_KEY', 'User') }
if (-not $key) {
    Write-Host 'Nenhuma chave configurada. Rode antes: .\01_setup_keys.ps1' -ForegroundColor Red
    exit 1
}

if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
    Write-Host 'Docker nao encontrado.' -ForegroundColor Red
    Write-Host 'Instale com: .\00_bootstrap_windows.ps1 -IncludeDocker' -ForegroundColor Yellow
    Write-Host 'Docker exige WSL2. Mas lembre: 8 GB de VRAM nao rodam NIM de estrutura.' -ForegroundColor Yellow
    exit 1
}

Write-Host 'Autenticando em nvcr.io...' -ForegroundColor Cyan
# O usuario e literalmente a string '$oauthtoken' - nao e uma variavel.
$key | docker login nvcr.io --username '$oauthtoken' --password-stdin

if ($LASTEXITCODE -eq 0) {
    Write-Host @'

Autenticado. Exemplos de pull (confira a VRAM exigida antes):

    docker pull nvcr.io/nim/mit/boltz2:latest
    docker pull nvcr.io/nim/nvidia/llama-3.2-nv-embedqa-1b-v2:latest

Depois de subir um container, aponte o servidor MCP para ele:
    NVIDIA_NIM_MODE=local  e  NVIDIA_NIM_LOCAL_BASE=http://localhost:8000

'@ -ForegroundColor Green
} else {
    Write-Host 'Login falhou. A chave e valida? Confira em build.nvidia.com' -ForegroundColor Red
    exit 1
}
