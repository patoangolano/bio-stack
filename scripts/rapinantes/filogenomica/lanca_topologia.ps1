# lanca_topologia.ps1 — starts roda_topologia.sh (AU/KH/SH on the 020210 supermatrix) as a scheduled task of the
# current user, outside Claude Code's process tree, as retomada_020210/lanca_fixo.ps1 does for the supermatrix.
# -Threads goes to roda_topologia.sh as IQ-TREE's -T.
param([int]$Threads = 18)
$ErrorActionPreference = 'Stop'
$nome = 'topologia_AA_020210'
$script = '/mnt/c/Users/Manec/bio-stack/scripts/rapinantes/filogenomica/roda_topologia.sh'
$acao = New-ScheduledTaskAction -Execute 'C:\Windows\System32\conhost.exe' `
    -Argument "--headless C:\Windows\System32\wsl.exe -d Ubuntu-24.04 -- bash -l $script $Threads" `
    -WorkingDirectory 'C:\Users\Manec\bio-stack'
$cfg = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -ExecutionTimeLimit ([TimeSpan]::Zero) -MultipleInstances IgnoreNew
Register-ScheduledTask -TaskPath '\bio-stack\' -TaskName $nome -Action $acao -Settings $cfg -Force | Out-Null
Start-ScheduledTask -TaskPath '\bio-stack\' -TaskName $nome
Start-Sleep -Seconds 30
Get-ScheduledTask -TaskPath '\bio-stack\' | ForEach-Object {
    $i = $_ | Get-ScheduledTaskInfo
    '{0}: estado={1} ultimo_resultado={2} inicio={3}' -f $_.TaskName, $_.State, $i.LastTaskResult, $i.LastRunTime
}
