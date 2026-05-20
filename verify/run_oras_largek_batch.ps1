param(
    [string]$Tag = "largek",
    [string]$Parpool = "off",
    [int]$Workers = 0,
    [string]$KVals = "40 80 120",
    [string]$Degrees = "1 2 3",
    [string]$QVals = "10",
    [double]$StripOverlapExtension = 0.25,
    [int]$TimeoutSeconds = 21600
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = Resolve-Path (Join-Path $scriptDir "..")
$repoMatlab = ($repo.Path -replace "\\", "/")
$tagSafe = ($Tag -replace "[^A-Za-z0-9_]", "_")

$stdout = Join-Path $scriptDir "oras_largek_${tagSafe}_stdout.log"
$stderr = Join-Path $scriptDir "oras_largek_${tagSafe}_stderr.log"
$watchdog = Join-Path $scriptDir "oras_largek_${tagSafe}_watchdog.log"

$workerLine = ""
if ($Workers -gt 0) {
    $workerLine = "setenv('ORAS_LARGEK_WORKERS','$Workers'); "
}

$batch = "cd('$repoMatlab'); " +
    "setenv('ORAS_LARGEK_KVALS','$KVals'); " +
    "setenv('ORAS_LARGEK_DEGREES','$Degrees'); " +
    "setenv('ORAS_LARGEK_QVALS','$QVals'); " +
    "setenv('ORAS_LARGEK_PARPOOL','$Parpool'); " +
    "setenv('ORAS_LARGEK_STRIP_OVERLAP_EXTENSION','$StripOverlapExtension'); " +
    $workerLine +
    "setenv('ORAS_LARGEK_TAG','$tagSafe'); " +
    "addpath(genpath('.')); run('verify/verify_oras_largek_iterations.m');"

$launchedAt = Get-Date

"Launcher started at $($launchedAt.ToString('o')), timeout=$TimeoutSeconds, tag=$tagSafe, parpool=$Parpool, workers=$Workers" |
    Out-File -FilePath $watchdog -Encoding utf8

$argumentLine = "-nosplash -nodesktop -batch `"$batch`""

$proc = Start-Process -FilePath "matlab" `
    -ArgumentList $argumentLine `
    -WorkingDirectory $repo.Path `
    -WindowStyle Hidden `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -PassThru

"MATLAB pid=$($proc.Id)" | Out-File -FilePath $watchdog -Encoding utf8 -Append

function Get-LaunchedMatlabProcesses {
    Get-Process matlab -ErrorAction SilentlyContinue |
        Where-Object { $_.StartTime -ge $launchedAt.AddSeconds(-5) }
}

$deadline = $launchedAt.AddSeconds($TimeoutSeconds)
while ((Get-LaunchedMatlabProcesses | Measure-Object).Count -gt 0 -and (Get-Date) -lt $deadline) {
    Start-Sleep -Seconds 60
}

$remaining = @(Get-LaunchedMatlabProcesses)
if ($remaining.Count -gt 0) {
    $pidList = ($remaining | ForEach-Object { $_.Id }) -join ","
    "Timeout reached at $(Get-Date -Format o); stopping MATLAB pids=$pidList" |
        Out-File -FilePath $watchdog -Encoding utf8 -Append
    $remaining | Stop-Process -Force
} else {
    "MATLAB exited at $(Get-Date -Format o)" |
        Out-File -FilePath $watchdog -Encoding utf8 -Append
}
