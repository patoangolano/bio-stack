<#
.SYNOPSIS
    Materializa o ambiente Python da bio-stack a partir do uv.lock.

.DESCRIPTION
    Este script nao tem lista de pacotes, e isso e deliberado. A lista vive em
    `pyproject.toml` e o fecho transitivo com hash vive em `uv.lock`; qualquer
    segunda copia da lista viraria uma segunda fonte de verdade que divergiria
    em silencio. Era exatamente o que acontecia com os antigos
    `requirements-core.txt` / `requirements-gpu.txt`: eles resolviam scanpy,
    pydeseq2, statsmodels, scipy e numpy por piso de versao, entao um bootstrap
    rodado hoje podia instalar uma versao diferente da travada — e sao essas
    cinco que calculam FDR e padj. Metodos que citam versao ficavam
    inverificaveis.

    Camadas (todas resolvidas pelo uv, nenhuma listada aqui):

      base  - [project.dependencies]                                 (sempre)
      dev   - [dependency-groups].dev, linters em .venv\Scripts\     (sempre)
      gpu   - [dependency-groups].gpu, torch cu128 + transformers    (sempre)

    As TRES entram num sync pelado, porque `[tool.uv] default-groups` lista
    dev e gpu. Isso e deliberado e conserta uma armadilha medida: enquanto a
    camada gpu era um extra (`[project.optional-dependencies]`), um
    `uv sync --locked` sem `--extra gpu` DESINSTALAVA 26 pacotes, o
    torch==2.11.0+cu128 entre eles. Como o passo 3 do 00_bootstrap_windows.ps1
    chama este script sem flags, o caminho documentado de instalacao destruia
    a camada GPU montada no passo anterior. Verificado com --dry-run em
    2026-10-10, antes e depois.

    Para o ambiente magro, sem os ~3 GB de wheel de CUDA, use -SemGpu.
    -Gpu e -Dev sobrevivem so por compatibilidade com quem ja chamava assim;
    hoje nao fazem nada, porque as duas camadas ja sao padrao.

    O indice de wheels do PyTorch nao e mais parametro: ele esta fixado em
    `[[tool.uv.index]] pytorch-cu128` com `explicit = true`, e `torch` aponta
    para ele em `[tool.uv.sources]`. Passar o indice pela linha de comando era
    o que permitia instalar um wheel cu121 sem reclamacao nenhuma — ele instala
    e so falha no primeiro kernel, porque nao ha wheel oficial para sm_120
    abaixo de cu128.

    Nota honesta sobre GPU no Windows:
      - PyTorch CUDA funciona nativo no Windows. OK.
      - RAPIDS / rapids-singlecell / cuDF NAO tem build nativo para Windows.
        Precisa de WSL2 (Ubuntu) ou Linux. O script avisa e nao tenta instalar.
      - NVIDIA Parabricks exige >= 16 GB de VRAM e Linux/Docker. Fora desta
        maquina; use DGX Cloud / cluster da universidade ou os NIM hospedados.

    Orcamento de VRAM nesta maquina (RTX 5070 Laptop, 8 GB), medido:
      ESM-2 t33 650M      ~2,7 GB fp16   -> cabe folgado
      ESM-2 t36 3B        ~6,5 GB fp16   -> cabe apertado, batch 1
      DNABERT-2 117M      ~0,5 GB        -> cabe folgado
      ESMFold             ~14 GB         -> NAO cabe, use o NIM hospedado
      Evo 2 7B / 40B      >= 40 GB       -> NAO cabe, use o NIM hospedado

    onnxruntime-gpu estava no antigo requirements-gpu.txt e NAO entrou no extra
    `gpu` do pyproject. Decisao, nao esquecimento: nada no projeto usa ONNX, e
    sao ~250 MB de wheel que ainda disputam a mesma VRAM de 8 GB. Se algum dia
    entrar um modelo exportado para ONNX, a linha volta no pyproject — nao
    aqui.

.EXAMPLE
    .\02_python_env.ps1 -Gpu -Dev

.EXAMPLE
    .\02_python_env.ps1 -Gpu -Dev -WhatIf
    Mostra cada chamada ao uv sem executar nenhuma.
#>

[CmdletBinding(SupportsShouldProcess)]
param(
    # Nao fazem mais nada: dev e gpu sao default-groups no pyproject. Ficam
    # declarados para que `-Gpu -Dev` de quem ja chamava assim nao estoure
    # como parametro desconhecido.
    [switch]$Gpu,
    [switch]$Dev,
    # O inverso e que e util agora: pular os ~3 GB de wheel de CUDA.
    [switch]$SemGpu,
    [string]$PythonVersion = '3.11',
    # Por padrao o sync exige que o uv.lock ja cubra o pyproject (`--locked`):
    # se alguem editou o pyproject e esqueceu de travar, o script falha alto em
    # vez de resolver versao nova em silencio. -Relock troca isso por um sync
    # que atualiza o lock — e uma decisao, entao precisa ser pedida.
    [switch]$Relock
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # .../bio-stack
$venv = Join-Path $root '.venv'

function Write-Step($m) { Write-Host "`n=== $m ===" -ForegroundColor Cyan }

function Invoke-Uv {
    param(
        [string]$Alvo,
        [string[]]$Argumentos
    )
    Write-Host "  uv $($Argumentos -join ' ')" -ForegroundColor DarkGray
    if (-not $PSCmdlet.ShouldProcess($Alvo, "uv $($Argumentos -join ' ')")) { return }
    & uv @Argumentos
    if ($LASTEXITCODE -ne 0) { throw "uv $($Argumentos[0]) falhou (exit $LASTEXITCODE)" }
}

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Host 'uv nao encontrado. Rode antes: .\00_bootstrap_windows.ps1' -ForegroundColor Red
    exit 1
}

$syncArgs = @(
    'sync'
    '--project', $root
    '--python', $PythonVersion
)
if (-not $Relock) { $syncArgs += '--locked' }
if ($SemGpu) { $syncArgs += @('--no-group', 'gpu') }

# A mensagem tem que descrever o que o uv VAI fazer, nao o que as flags
# sugerem. O header antigo afirmava que dev so entrava com -Dev, e isso era
# falso mesmo antes do default-groups: o uv ja tratava dev como grupo padrao.
$camadas = @('base', 'dev')
if (-not $SemGpu) { $camadas += 'gpu' }
if ($Gpu -or $Dev) {
    Write-Host "  nota: -Gpu/-Dev nao fazem nada; dev e gpu ja sao default-groups." -ForegroundColor DarkGray
}

Write-Step "Sincronizando .venv com o uv.lock (Python $PythonVersion; camadas: $($camadas -join ', '))"
try {
    Invoke-Uv -Alvo $venv -Argumentos $syncArgs
}
catch {
    if (-not $Relock) {
        Write-Host @'
  O sync com --locked falhou. A causa normal e pyproject.toml editado sem
  travar: rode `uv lock` e confira o diff do uv.lock antes de commitar, ou
  repita este script com -Relock para deixar o uv travar sozinho.
'@ -ForegroundColor Yellow
    }
    throw
}

$py = Join-Path $venv 'Scripts\python.exe'

if ($Gpu) {
    Write-Step 'Verificando CUDA'
    # Sem acento nenhum no bloco: o python roda com o console em cp1252 e
    # travessao ou cedilha aqui sai como UnicodeEncodeError.
    $probe = @'
import torch
print(f"  torch           {torch.__version__}")
print(f"  torch.version.cuda {torch.version.cuda}")
print(f"  cuda disponivel {torch.cuda.is_available()}")
if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    gb = p.total_memory / 1024**3
    print(f"  gpu             {p.name}")
    print(f"  capability      sm_{p.major}{p.minor}")
    print(f"  vram            {gb:.1f} GB")
    if gb < 16:
        print("  [aviso] < 16 GB: use os NIM hospedados para predicao de estrutura.")
else:
    print("  [aviso] build de CPU instalado. Confira [tool.uv.sources] torch no pyproject.")
'@
    if ($PSCmdlet.ShouldProcess($py, 'python (probe de CUDA)')) {
        if (-not (Test-Path $py)) { throw "interpretador nao encontrado: $py" }
        # `python -c` com string de varias linhas nao sobrevive a chamada de
        # executavel nativo: o PowerShell quebra a string nas linhas e cada
        # pedaco chega como argumento separado, e o python morre com
        # "'(' was never closed". Arquivo resolve. `-I` porque o arquivo mora
        # no TEMP, e em modo isolado o python nao importa nada do diretorio do
        # script nem do cwd.
        $probeFile = Join-Path ([System.IO.Path]::GetTempPath()) 'bio-stack-cuda-probe.py'
        try {
            Set-Content -LiteralPath $probeFile -Value $probe -Encoding utf8
            & $py -I $probeFile
            if ($LASTEXITCODE -ne 0) { throw 'probe de CUDA falhou' }
        }
        finally {
            Remove-Item -LiteralPath $probeFile -ErrorAction SilentlyContinue
        }
    }
    else {
        Write-Host '  (pulado)' -ForegroundColor DarkGray
    }
}

Write-Step 'RAPIDS / GPU single-cell'
Write-Host @'
  rapids-singlecell e cuDF nao tem build nativo para Windows.
  Para usar a GPU em scRNA-seq (10-50x mais rapido que scanpy em CPU), dentro
  da distro WSL (aqui: Ubuntu-24.04):

      curl -LsSf https://astral.sh/uv/install.sh | sh
      uv venv --python 3.11 && source .venv/bin/activate
      uv pip install --extra-index-url=https://pypi.nvidia.com \
          cudf-cu12 cuml-cu12 cugraph-cu12 rapids-singlecell

  Sem WSL, scanpy em CPU resolve bem datasets ate ~50k celulas.
'@ -ForegroundColor Gray

Write-Host "`nAmbiente pronto: $venv" -ForegroundColor Green
Write-Host "Ativar:  $venv\Scripts\Activate.ps1`n" -ForegroundColor Green
