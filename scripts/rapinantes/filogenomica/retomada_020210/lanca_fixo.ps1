# lanca_fixo.ps1 — starts roda_fixo.sh as a scheduled task of the current user, outside Claude Code's process tree,
# so closing the session does not send SIGHUP down to IQ-TREE (what killed 020210 at 21:09 on 2026-10-18).
# conhost --headless keeps the console app windowless: a stray click closing a window would kill the run.
# -Threads goes to roda_fixo.sh as IQ-TREE's -T (8 = the pipeline's cpus_concatenado).
param([int]$Threads = 8)
$ErrorActionPreference = 'Stop'
$nome = 'supermatriz_fixa_020210'
$script = '/mnt/c/Users/Manec/bio-stack/scripts/rapinantes/filogenomica/retomada_020210/roda_fixo.sh'
$acao = New-ScheduledTaskAction -Execute 'C:\Windows\System32\conhost.exe' `
    -Argument "--headless C:\Windows\System32\wsl.exe -d Ubuntu-24.04 -- bash -l $script $Threads" `
    -WorkingDirectory 'C:\Users\Manec\bio-stack'
# No execution time limit, and the default "stop if on batteries" is off: this is a laptop.
$cfg = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -ExecutionTimeLimit ([TimeSpan]::Zero) -MultipleInstances IgnoreNew
Register-ScheduledTask -TaskPath '\bio-stack\' -TaskName $nome -Action $acao -Settings $cfg -Force | Out-Null
Start-ScheduledTask -TaskPath '\bio-stack\' -TaskName $nome
Start-Sleep -Seconds 30
Get-ScheduledTask -TaskPath '\bio-stack\' | ForEach-Object {
    $i = $_ | Get-ScheduledTaskInfo
    '{0}: estado={1} ultimo_resultado={2} inicio={3}' -f $_.TaskName, $_.State, $i.LastTaskResult, $i.LastRunTime
}
