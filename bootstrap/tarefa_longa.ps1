<#
.SYNOPSIS
    Roda um script de horas como tarefa agendada do usuario, fora da arvore de
    processos do Claude Code.

.DESCRIPTION
    Job de horas iniciado pelo Bash tool morre de SIGHUP quando a sessao fecha
    (medido em 2026-10-18: duas execucoes de filogenia perdidas, e o WSL desligou
    sem cliente). A tarefa roda sob o Agendador, com conhost --headless: nao ha
    janela que um clique perdido feche.

    Generaliza os registradores avulsos que viviam nos scratchpads de sessao
    (registra_uv_tools.ps1, registra_wsl_pixi.ps1, lanca_backup.ps1,
    lanca_retomadas.ps1) e os .cmd que os acompanhavam. Gera um .cmd ao lado do
    log e registra `conhost --headless cmd /c <esse .cmd>`. O .cmd continua
    existindo por dois motivos medidos:
      - o redirecionamento do cmd guarda os bytes crus; o do PowerShell 5.1
        re-codifica a saida em UTF-16;
      - o %ERRORLEVEL% so reflete o comando anterior numa linha propria. Numa
        linha unica com `&`, o cmd expande a variavel antes de o comando rodar.

    Sem -Wsl o script roda no Git Bash do Windows, em shell de login; com -Wsl,
    em `wsl.exe -d Ubuntu-24.04 -- bash -l`. Caminhos absolutos em tudo: a
    tarefa nao herda o PATH do shell interativo.

    O script recebe o caminho do proprio log como $1, e o stdout/stderr dele
    vai para <Nome>.cmd.log, com o codigo de saida na ultima linha. Os dois
    ficam em runs/bootstrap/, fora do git.

    -EsperarPor segura o inicio enquanto outra tarefa de \bio-stack\ estiver
    Running. Duas reconstrucoes baixando ao mesmo tempo foi o que matou a
    primeira passada das uv tools (o torch cu128 esgotou 4 retries em 7674 s,
    2026-10-21). A espera fica dentro da tarefa, e nao em quem a registra; como
    a tarefa conta como Running enquanto espera, da para enfileirar varias de
    uma vez, cada uma esperando as anteriores.

    -Ambiente poe variaveis no .cmd, antes do script. E o jeito de mudar uma
    execucao sem editar o script, o que nao se faz com ele rodando: o bash le
    o arquivo por deslocamento conforme avanca, e uma edicao no meio troca o
    que ele executa.

.EXAMPLE
    .\bootstrap\tarefa_longa.ps1 -Nome reconstruir-windows-uv -Script bootstrap\08_windows_uv.sh

.EXAMPLE
    .\bootstrap\tarefa_longa.ps1 -Nome reconstruir-wsl -Script bootstrap\09_wsl_reconstruir.sh -EsperarPor reconstruir-windows-uv

.EXAMPLE
    .\bootstrap\tarefa_longa.ps1 -Nome backup-drive -Script bootstrap\backup_drive.sh -Wsl

.EXAMPLE
    .\bootstrap\tarefa_longa.ps1 -Nome uv-2 -Script bootstrap\08_windows_uv.sh -Ambiente @{ UV_HTTP_TIMEOUT = '300' }
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory)] [string]$Nome,
    [Parameter(Mandatory)] [string]$Script,
    [switch]$Wsl,
    [string[]]$EsperarPor = @(),
    [hashtable]$Ambiente = @{},
    [string]$Distro = 'Ubuntu-24.04'
)

$ErrorActionPreference = 'Stop'
$raiz    = Split-Path -Parent $PSScriptRoot
$pasta   = '\bio-stack\'
$logDir  = Join-Path $raiz 'runs\bootstrap'
$gitBash = 'C:\Program Files\Git\bin\bash.exe'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$alvo = if ([IO.Path]::IsPathRooted($Script)) { $Script } else { Join-Path $raiz $Script }
if (-not (Test-Path -LiteralPath $alvo)) { throw "script nao encontrado: $alvo" }
$alvo = (Resolve-Path -LiteralPath $alvo).Path

$atual = Get-ScheduledTask -TaskName $Nome -TaskPath $pasta -ErrorAction SilentlyContinue
if ($atual -and $atual.State -eq 'Running') {
    throw "a tarefa $pasta$Nome ja esta rodando; espere terminar ou pare com Stop-ScheduledTask"
}

$log    = Join-Path $logDir "$Nome.log"
$cmdLog = Join-Path $logDir "$Nome.cmd.log"
$cmd    = Join-Path $logDir "$Nome.cmd"

# C:\a\b -> /mnt/c/a/b (WSL) ou C:/a/b (Git Bash, que aceita a barra normal).
function ConverterWsl([string]$p) { '/mnt/' + $p.Substring(0, 1).ToLower() + ($p.Substring(2) -replace '\\', '/') }
function ConverterGitBash([string]$p) { $p -replace '\\', '/' }

if ($Wsl) {
    $exec = 'C:\Windows\System32\wsl.exe -d {0} -- bash -l "{1}" "{2}"' -f $Distro, (ConverterWsl $alvo), (ConverterWsl $log)
} else {
    if (-not (Test-Path -LiteralPath $gitBash)) { throw "Git Bash nao encontrado em $gitBash" }
    $exec = '"{0}" -l "{1}" "{2}"' -f $gitBash, (ConverterGitBash $alvo), (ConverterGitBash $log)
}

$linhas = @(
    '@echo off'
    ('rem gerado por bootstrap\tarefa_longa.ps1 em {0:yyyy-MM-dd HH:mm}; registre de novo em vez de editar' -f (Get-Date))
    ('echo inicio %DATE% %TIME% > "{0}"' -f $cmdLog)
)
if ($EsperarPor.Count -gt 0) {
    $lista = ($EsperarPor | ForEach-Object { "'$_'" }) -join ','
    $espera = "foreach (`$t in @($lista)) { while ((Get-ScheduledTask -TaskPath '$pasta' -TaskName `$t -ErrorAction SilentlyContinue).State -eq 'Running') { Start-Sleep -Seconds 30 } }"
    $linhas += ('echo esperando: {0} >> "{1}"' -f ($EsperarPor -join ' '), $cmdLog)
    $linhas += ('powershell.exe -NoProfile -Command "{0}"' -f $espera)
    $linhas += ('echo fim da espera %DATE% %TIME% >> "{0}"' -f $cmdLog)
}
# `set "NOME=valor"` com as aspas por fora: o valor pode ter espaco (G:\Meu
# Drive) e nao leva aspa nenhuma para dentro da variavel. No modo -Wsl a
# variavel nao atravessa para a distro sem WSLENV; aqui so o Git Bash a usa.
foreach ($chave in $Ambiente.Keys) {
    $linhas += ('set "{0}={1}"' -f $chave, $Ambiente[$chave])
    $linhas += ('echo ambiente: {0}={1} >> "{2}"' -f $chave, $Ambiente[$chave], $cmdLog)
}
$linhas += ('{0} >> "{1}" 2>&1' -f $exec, $cmdLog)
$linhas += ('echo exit=%ERRORLEVEL% >> "{0}"' -f $cmdLog)
Set-Content -LiteralPath $cmd -Value $linhas -Encoding ASCII

$acao = New-ScheduledTaskAction -Execute 'C:\Windows\System32\conhost.exe' `
    -Argument ('--headless C:\Windows\System32\cmd.exe /c "{0}"' -f $cmd) `
    -WorkingDirectory $raiz
$principal = New-ScheduledTaskPrincipal -UserId ('{0}\{1}' -f $env:USERDOMAIN, $env:USERNAME) `
    -LogonType Interactive -RunLevel Limited
# ExecutionTimeLimit zero = sem limite. O padrao (72 h) ate serviria, mas o
# padrao "parar se estiver na bateria" nao: isto e um laptop.
$config = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -StartWhenAvailable -ExecutionTimeLimit ([TimeSpan]::Zero) -MultipleInstances IgnoreNew

Register-ScheduledTask -TaskPath $pasta -TaskName $Nome -Action $acao -Principal $principal `
    -Settings $config -Description ("bio-stack: {0}" -f $alvo) -Force | Out-Null
Start-ScheduledTask -TaskPath $pasta -TaskName $Nome
Start-Sleep -Seconds 5

$tarefa = Get-ScheduledTask -TaskPath $pasta -TaskName $Nome
$info   = $tarefa | Get-ScheduledTaskInfo
'tarefa   {0}{1}: {2} (resultado {3})' -f $pasta, $Nome, $tarefa.State, $info.LastTaskResult
'script   {0}' -f $alvo
'log      {0}' -f $log
'saida    {0}' -f $cmdLog
'remover  Unregister-ScheduledTask -TaskPath ''{0}'' -TaskName {1} -Confirm:$false' -f $pasta, $Nome
