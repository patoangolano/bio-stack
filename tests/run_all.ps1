<#
.SYNOPSIS
    Roda toda a suite de testes da bio-stack.
.DESCRIPTION
    Nao faz nenhuma chamada de rede e nao consome credito da NVIDIA.
    Exige apenas `uv` no PATH - as dependencias sao resolvidas por script.
.PARAMETER Completo
    Inclui as suites que dependem do WSL: check_wsl_ferramentas.py (as CLIs
    do pixi) e check_filogenia.py (pipeline ponta a ponta com Docker).
    A segunda baixa ~0,3 GB de imagens na primeira vez (~1 GB no disco) - por
    isso fica fora do
    padrao, que promete nao usar rede.
.EXAMPLE
    .\tests\run_all.ps1
.EXAMPLE
    .\tests\run_all.ps1 -Completo
#>
param([switch]$Completo)
$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent $PSScriptRoot
$failed = @()

$suites = @('check_hooks.py', 'check_server.py', 'check_nemo.py', 'check_filogenia_figura.py')
if ($Completo) { $suites += @('check_wsl_ferramentas.py', 'check_filogenia.py') }

foreach ($suite in $suites) {
    Write-Host "`n=== $suite ===" -ForegroundColor Cyan
    & uv run --script (Join-Path $PSScriptRoot $suite)
    if ($LASTEXITCODE -ne 0) { $failed += $suite }
}

Write-Host ''
if ($failed.Count -gt 0) {
    Write-Host "FALHARAM: $($failed -join ', ')" -ForegroundColor Red
    exit 1
}
Write-Host 'Toda a suite passou.' -ForegroundColor Green
