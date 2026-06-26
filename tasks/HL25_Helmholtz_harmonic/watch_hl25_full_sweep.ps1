param(
    [int]$IntervalSeconds = 300,
    [int]$StaleMinutes = 30,
    [switch]$Once
)

$ErrorActionPreference = 'Continue'
$TaskDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$StatusPath = Join-Path $TaskDir 'full_sweep_watchdog_status.txt'
$CsvPath = Join-Path $TaskDir 'full_sweep_lxzz_cross_results.csv'

function Get-LatestLog {
    $logs = @()
    $logs += Get-ChildItem -LiteralPath $TaskDir -Filter 'full_sweep_corrected_*.out.log' -ErrorAction SilentlyContinue
    $logs += Get-ChildItem -LiteralPath $TaskDir -Filter 'full_sweep_*.out.log' -ErrorAction SilentlyContinue
    $logs | Sort-Object LastWriteTime -Descending | Select-Object -First 1
}

function Get-LatestErrLog {
    $logs = @()
    $logs += Get-ChildItem -LiteralPath $TaskDir -Filter 'full_sweep_corrected_*.err.log' -ErrorAction SilentlyContinue
    $logs += Get-ChildItem -LiteralPath $TaskDir -Filter 'full_sweep_*.err.log' -ErrorAction SilentlyContinue
    $logs | Sort-Object LastWriteTime -Descending | Select-Object -First 1
}

function Get-SweepProcesses {
    $script:ProcessQueryWarning = $null
    try {
        return @(Get-CimInstance Win32_Process -Filter "name='MATLAB.exe'" -ErrorAction Stop |
            Where-Object { $_.CommandLine -like '*verify_hl25_full_sweep*' })
    } catch {
        $script:ProcessQueryWarning = $_.Exception.Message
        return @()
    }
}

function Get-CsvStatusLines {
    if (!(Test-Path -LiteralPath $CsvPath)) {
        return @('CSV: missing')
    }
    $rows = Import-Csv -LiteralPath $CsvPath
    $lines = @("CSV rows: $($rows.Count)")
    $lines += ($rows | Group-Object status | Sort-Object Name |
        ForEach-Object { "  $($_.Name): $($_.Count)" })
    return $lines
}

function Get-OverallState($processes, $latestLog, $statusLines) {
    $statusText = ($statusLines -join "`n")
    $unfinished = $statusText -match 'estimated_only|running_pending'
    $failed = $statusText -match 'failed'
    if ($processes.Count -eq 0) {
        if ($unfinished) { return 'STOPPED_INCOMPLETE' }
        if ($failed) { return 'ENDED_WITH_FAILURES' }
        return 'COMPLETED'
    }
    if ($null -eq $latestLog) {
        return 'RUNNING_NO_LOG'
    }
    $age = (Get-Date) - $latestLog.LastWriteTime
    if ($age.TotalMinutes -gt $StaleMinutes) {
        return 'RUNNING_STALE_LOG'
    }
    return 'RUNNING'
}

function Write-WatchdogStatus {
    $now = Get-Date
    $processes = @(Get-SweepProcesses)
    $latestLog = Get-LatestLog
    $latestErr = Get-LatestErrLog
    $statusLines = @(Get-CsvStatusLines)
    $state = Get-OverallState $processes $latestLog $statusLines

    $report = New-Object System.Collections.Generic.List[string]
    $report.Add("HL25 full sweep watchdog")
    $report.Add("Updated: $($now.ToString('yyyy-MM-dd HH:mm:ss'))")
    $report.Add("State: $state")
    $report.Add("")
    $report.Add("Processes:")
    if ($script:ProcessQueryWarning) {
        $report.Add("  process query warning: $script:ProcessQueryWarning")
    }
    if ($processes.Count -eq 0) {
        $report.Add("  none")
    } else {
        foreach ($p in $processes) {
            $report.Add("  PID $($p.ProcessId), started $($p.CreationDate)")
        }
    }
    $report.Add("")
    if ($null -ne $latestLog) {
        $age = (Get-Date) - $latestLog.LastWriteTime
        $report.Add("Latest stdout log: $($latestLog.FullName)")
        $report.Add("Latest stdout log updated: $($latestLog.LastWriteTime) (age $([math]::Round($age.TotalMinutes, 1)) min)")
    } else {
        $report.Add("Latest stdout log: none")
    }
    if ($null -ne $latestErr) {
        $report.Add("Latest stderr log: $($latestErr.FullName) ($($latestErr.Length) bytes)")
    } else {
        $report.Add("Latest stderr log: none")
    }
    $report.Add("")
    foreach ($line in $statusLines) {
        $report.Add([string]$line)
    }
    $report.Add("")
    $report.Add("Recent stdout tail:")
    if ($null -ne $latestLog -and (Test-Path -LiteralPath $latestLog.FullName)) {
        $tail = Get-Content -LiteralPath $latestLog.FullName -Tail 35 -ErrorAction SilentlyContinue
        foreach ($line in $tail) { $report.Add($line) }
    } else {
        $report.Add("  unavailable")
    }

    $text = $report -join [Environment]::NewLine
    Set-Content -LiteralPath $StatusPath -Value $text -Encoding UTF8
    Write-Output $text
}

do {
    Write-WatchdogStatus
    if ($Once) { break }
    Start-Sleep -Seconds $IntervalSeconds
} while ($true)
