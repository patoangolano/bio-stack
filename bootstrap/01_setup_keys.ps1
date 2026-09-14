<#
.SYNOPSIS
    Grava as chaves de API da bio-stack como variaveis de ambiente de usuario.

.DESCRIPTION
    As chaves ficam no ambiente do usuario do Windows (HKCU), NUNCA em arquivo
    versionado. O .mcp.json e os hooks leem delas por expansao ${VAR}.

    Onde conseguir cada chave:

      NVIDIA_API_KEY   (obrigatoria para a camada NVIDIA)
          https://build.nvidia.com -> canto superior direito -> "Get API Key"
          Formato: nvapi-...
          Com NVIDIA AI Enterprise a mesma chave NGC vale para:
            - endpoints hospedados  health.api.nvidia.com / integrate.api.nvidia.com
            - pull de containers    nvcr.io  (docker login -u '$oauthtoken')

      NCBI_API_KEY     (opcional, muito recomendada)
          https://account.ncbi.nlm.nih.gov/settings/ -> API Key Management
          Sobe o rate limit do E-utilities de 3 req/s para 10 req/s.

      S2_API_KEY       (opcional) https://www.semanticscholar.org/product/api
      OPENFDA_API_KEY  (opcional) https://open.fda.gov/apis/authentication/
      HF_TOKEN         (opcional) https://huggingface.co/settings/tokens

.EXAMPLE
    .\01_setup_keys.ps1
    .\01_setup_keys.ps1 -Show
#>

[CmdletBinding()]
param(
    # Apenas mostra o que ja esta configurado, sem pedir nada.
    [switch]$Show
)

$ErrorActionPreference = 'Stop'

$Keys = [ordered]@{
    'NVIDIA_API_KEY'  = @{ Required = $true;  Hint = 'build.nvidia.com -> Get API Key (nvapi-...)' }
    'NCBI_API_KEY'    = @{ Required = $false; Hint = 'account.ncbi.nlm.nih.gov/settings -> API Key' }
    'S2_API_KEY'      = @{ Required = $false; Hint = 'Semantic Scholar API' }
    'OPENFDA_API_KEY' = @{ Required = $false; Hint = 'open.fda.gov' }
    'HF_TOKEN'        = @{ Required = $false; Hint = 'huggingface.co/settings/tokens' }
}

function Get-Stored($name) {
    [Environment]::GetEnvironmentVariable($name, 'User')
}

function Format-Masked($value) {
    if ([string]::IsNullOrWhiteSpace($value)) { return '(vazio)' }
    if ($value.Length -le 12) { return ('*' * $value.Length) }
    "$($value.Substring(0,8))...$($value.Substring($value.Length-4))"
}

Write-Host "`n=== Chaves da bio-stack ===`n" -ForegroundColor Cyan

foreach ($name in $Keys.Keys) {
    $current = Get-Stored $name
    $req = if ($Keys[$name].Required) { 'obrigatoria' } else { 'opcional   ' }
    $state = if ($current) { "OK  $(Format-Masked $current)" } else { 'NAO CONFIGURADA' }
    $color = if ($current) { 'Green' } elseif ($Keys[$name].Required) { 'Red' } else { 'DarkGray' }
    Write-Host ("  {0,-16} [{1}]  {2}" -f $name, $req, $state) -ForegroundColor $color
}

if ($Show) { Write-Host ''; exit 0 }

Write-Host "`nDeixe em branco e pressione Enter para pular / manter o valor atual.`n" -ForegroundColor Gray

foreach ($name in $Keys.Keys) {
    $current = Get-Stored $name
    $label = if ($current) { "$name (atual: $(Format-Masked $current))" } else { "$name  - $($Keys[$name].Hint)" }

    # Read-Host -AsSecureString evita que a chave fique no historico do console.
    $secure = Read-Host -Prompt "  $label" -AsSecureString
    $bstr = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secure)
    try {
        $plain = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($bstr)
    } finally {
        [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($bstr)
    }

    if ([string]::IsNullOrWhiteSpace($plain)) { continue }

    $plain = $plain.Trim()

    if ($name -eq 'NVIDIA_API_KEY' -and -not $plain.StartsWith('nvapi-')) {
        Write-Host "    [aviso] chaves NVIDIA normalmente comecam com 'nvapi-'. Gravando mesmo assim." -ForegroundColor Yellow
    }

    [Environment]::SetEnvironmentVariable($name, $plain, 'User')
    Set-Item -Path "env:$name" -Value $plain
    Write-Host "    gravada." -ForegroundColor Green
}

# ---------------------------------------------------------------------------
# Validacao ao vivo do endpoint NVIDIA.
$nv = Get-Stored 'NVIDIA_API_KEY'
if ($nv) {
    Write-Host "`n=== Validando NVIDIA_API_KEY ===" -ForegroundColor Cyan
    try {
        # Endpoint de modelos: barato e nao consome credito de inferencia.
        $resp = Invoke-RestMethod -Uri 'https://integrate.api.nvidia.com/v1/models' `
                                  -Headers @{ Authorization = "Bearer $nv" } `
                                  -Method Get -TimeoutSec 25
        $n = @($resp.data).Count
        Write-Host "  [ok] chave valida - $n modelos visiveis no catalogo." -ForegroundColor Green
    } catch {
        Write-Host "  [erro] a chave foi rejeitada ou a rede falhou:" -ForegroundColor Red
        Write-Host "         $($_.Exception.Message)" -ForegroundColor Red
        Write-Host "  Confira em https://build.nvidia.com se a chave ainda esta ativa." -ForegroundColor Yellow
    }
}

Write-Host @'

Pronto. IMPORTANTE: reinicie o Claude Code para ele herdar as variaveis novas.

'@ -ForegroundColor Green
