<#
.SYNOPSIS
    Roda toda a suite de testes da bio-stack.
.DESCRIPTION
    Nao faz nenhuma chamada de rede e nao consome credito da NVIDIA.
    Exige apenas `uv` no PATH - as dependencias sao resolvidas por script.
.EXAMPLE
    .\tests\run_all.ps1
#>
$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent $PSScriptRoot
$failed = @()

foreach ($suite in @('test_hooks.py', 'test_server.py', 'test_nemo.py')) {
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
