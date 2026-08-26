param(
    [Parameter(Mandatory = $true)]
    [string]$SeedExporter,
    [Parameter(Mandatory = $true)]
    [string]$MaCetspExecutable,
    [Parameter(Mandatory = $true)]
    [string]$LkhExecutable,
    [Parameter(Mandatory = $true)]
    [string]$GurobiExecutable,
    [Parameter(Mandatory = $true)]
    [string]$ResultChecker,
    [int]$Runs = 20,
    [int]$Repetitions = 1000,
    [int]$BonusRepetitions = 10000,
    [int]$PoolSize = 20,
    [int]$Iterations = 5000,
    [int]$Patience = 500,
    [int]$InitialPatience = 1500,
    [double]$MaxTime = 36000,
    [double]$FitnessWeight = 0.96,
    [double]$MinimumDistance = 5.0,
    [int]$NeighborSize = 50,
    [int]$ConstructionThreads = 1,
    [int]$SolverThreads = 1,
    [int]$ParallelJobs = 8,
    [int]$Seed = 123456789,
    [string[]]$Instances = @(),
    [switch]$Resume,
    [string]$Output = "results/mennell_hybrid/tables/mennell_protocol_runs.csv",
    [string]$Summary = "results/mennell_hybrid/tables/mennell_protocol_summary.csv",
    [string]$WorkDirectory = "build-hybrid/ma-mennell-protocol"
)

$ErrorActionPreference = "Stop"

if ($Runs -le 0 -or $Repetitions -le 0 -or $BonusRepetitions -le 0 -or
    $PoolSize -le 0 -or $Iterations -le 0 -or $Patience -le 0 -or
    $InitialPatience -le 0 -or
    $MaxTime -le 0 -or $ParallelJobs -le 0 -or
    $ConstructionThreads -le 0 -or $SolverThreads -le 0) {
    throw "Run counts, search limits, parallelism, and solver threads must be positive"
}

$repository = (Resolve-Path (Join-Path $PSScriptRoot "../../..")).Path
$workerPath = Join-Path $PSScriptRoot "run_mennell_protocol_worker.ps1"
$seedExporterPath = (Resolve-Path $SeedExporter).Path
$maCetspPath = (Resolve-Path $MaCetspExecutable).Path
$lkhPath = (Resolve-Path $LkhExecutable).Path
$gurobiPath = (Resolve-Path $GurobiExecutable).Path
$resultCheckerPath = (Resolve-Path $ResultChecker).Path
$outputPath = Join-Path $repository $Output
$summaryPath = Join-Path $repository $Summary
$workPath = Join-Path $repository $WorkDirectory
$tourDirectory = Join-Path $repository "results/mennell_hybrid/tours"
$pwshPath = (Get-Process -Id $PID).Path

New-Item -ItemType Directory -Force -Path $workPath | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $outputPath) | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $summaryPath) | Out-Null
New-Item -ItemType Directory -Force -Path $tourDirectory | Out-Null

$benchmarkDirectory = Join-Path $repository "data/mennell/reference_tours"
$inputDirectory = Join-Path $repository "data/mennell/instances"
$files = Get-ChildItem -LiteralPath $benchmarkDirectory -Filter "*.txt" |
    Sort-Object Name
if ($Instances.Count -gt 0) {
    $wanted = @{}
    foreach ($name in $Instances) {
        $wanted[$name] = $true
        if (-not $name.EndsWith(".txt")) {
            $wanted["$name.txt"] = $true
        }
    }
    $files = $files | Where-Object { $wanted.ContainsKey($_.Name) }
}
if ($files.Count -eq 0) {
    throw "No matching benchmark instances"
}

$jobs = [System.Collections.Generic.List[object]]::new()
$instanceSpecs = [System.Collections.Generic.List[object]]::new()
$precisionByInstance = @{}
foreach ($file in $files) {
    $inputPath = Join-Path $inputDirectory $file.Name
    if (-not (Test-Path -LiteralPath $inputPath)) {
        Write-Warning "Skipping $($file.Name): no processed input"
        continue
    }
    $benchmarkLine = Get-Content -LiteralPath $file.FullName -TotalCount 1
    $valueMatch = [regex]::Match($benchmarkLine, 'value\s*:\s*(.+)$')
    if (-not $valueMatch.Success) {
        throw "Benchmark value missing from $($file.FullName)"
    }
    $valueText = $valueMatch.Groups[1].Value.Trim()
    $leiHao = [double]$valueText
    $decimalIndex = $valueText.IndexOf('.')
    $publishedPrecision = if ($decimalIndex -lt 0) {
        0
    } else {
        $valueText.Length - $decimalIndex - 1
    }
    $precisionByInstance[$file.Name] = $publishedPrecision
    $repeatCount = if ($file.Name -eq "bonus1000.txt") {
        $BonusRepetitions
    } else {
        $Repetitions
    }
    $stem = [System.IO.Path]::GetFileNameWithoutExtension($file.Name)
    $instanceSpecs.Add([pscustomobject]@{
        File = $file
        InputPath = $inputPath
        LeiHao = $leiHao
        PublishedPrecision = $publishedPrecision
        Repetitions = $repeatCount
        Stem = $stem
    })
}
for ($run = 1; $run -le $Runs; ++$run) {
    foreach ($spec in $instanceSpecs) {
        $file = $spec.File
        $runDirectory = Join-Path $workPath (Join-Path $spec.Stem ("run-{0:D2}" -f $run))
        $jobs.Add([pscustomobject]@{
            Instance = $file.Name
            InputPath = $spec.InputPath
            LeiHao = $spec.LeiHao
            PublishedPrecision = $spec.PublishedPrecision
            Repetitions = $spec.Repetitions
            Run = $run
            Seed = $Seed + $run - 1
            RunDirectory = $runDirectory
            MetricsPath = Join-Path $runDirectory "metrics.json"
        })
    }
}

$results = [System.Collections.Generic.List[object]]::new()
$pending = [System.Collections.Generic.Queue[object]]::new()
foreach ($job in $jobs) {
    if ($Resume -and (Test-Path -LiteralPath $job.MetricsPath)) {
        $results.Add((Get-Content -LiteralPath $job.MetricsPath -Raw | ConvertFrom-Json))
        Remove-Item -LiteralPath (Join-Path $job.RunDirectory "failure.log") `
            -Force -ErrorAction SilentlyContinue
    } else {
        $pending.Enqueue($job)
    }
}
if (-not $Resume -and
    ($jobs | Where-Object { Test-Path -LiteralPath $_.MetricsPath } | Select-Object -First 1)) {
    throw "Existing run data found in $WorkDirectory; use -Resume or a new work directory"
}

function Export-Results {
    if ($results.Count -eq 0) {
        return
    }
    $orderedResults = $results | Sort-Object instance, @{ Expression = { [int]$_.run } }
    $orderedResults | Select-Object `
        instance, run, seed, repetitions, pool_size, iterations, patience,
        initial_patience, generated_best, ma_cetsp, final, lei_hao,
        improvement_vs_lei_hao_percent, seed_ms, ma_ms, socp_ms, ma_best_time,
        ma_total_time, ma_valid, final_valid |
        Export-Csv -LiteralPath $outputPath -NoTypeInformation

    $summaries = foreach ($group in ($orderedResults | Group-Object instance)) {
        $best = $group.Group | Sort-Object { [double]$_.final } | Select-Object -First 1
        $mean = ($group.Group | Measure-Object -Property final -Average).Average
        $publishedPrecision = $precisionByInstance[$group.Name]
        $roundingRadius = 0.5 * [math]::Pow(10.0, -$publishedPrecision)
        $comparison = if ([double]$best.final -lt
            [double]$best.lei_hao - $roundingRadius) {
            "improvement"
        } elseif ([double]$best.final -gt
            [double]$best.lei_hao + $roundingRadius) {
            "loss"
        } else {
            "tie"
        }
        $tourName = "{0}.result" -f (
            [System.IO.Path]::GetFileNameWithoutExtension($group.Name))
        $tourPath = Join-Path $tourDirectory $tourName
        Copy-Item -LiteralPath $best.final_result -Destination $tourPath -Force
        [pscustomobject]@{
            instance = $group.Name
            runs_completed = $group.Count
            best_run = $best.run
            best_seed = $best.seed
            generated_best = $best.generated_best
            ma_cetsp = $best.ma_cetsp
            best_final = $best.final
            mean_final = $mean
            lei_hao = $best.lei_hao
            published_precision = $publishedPrecision
            comparison_at_published_precision = $comparison
            improvement_vs_lei_hao_percent =
                100.0 * ([double]$best.lei_hao - [double]$best.final) /
                [double]$best.lei_hao
            final_valid = $best.final_valid
            best_tour = "results/mennell_hybrid/tours/$tourName"
        }
    }
    $summaries | Sort-Object instance |
        Export-Csv -LiteralPath $summaryPath -NoTypeInformation
}

$active = [System.Collections.Generic.List[object]]::new()
$failures = [System.Collections.Generic.List[object]]::new()
$total = $jobs.Count
$completedBeforeStart = $results.Count
$timer = [System.Diagnostics.Stopwatch]::StartNew()

Write-Host ((
    "Protocol: {0} instances x {1} runs; {2} queued, {3} resumed; " +
    "initial patience={4}, later patience={5}, parallel={6}") -f
    $instanceSpecs.Count, $Runs, $pending.Count, $completedBeforeStart,
    $InitialPatience, $Patience, $ParallelJobs)

while ($pending.Count -gt 0 -or $active.Count -gt 0) {
    while ($pending.Count -gt 0 -and $active.Count -lt $ParallelJobs) {
        $job = $pending.Dequeue()
        New-Item -ItemType Directory -Force -Path $job.RunDirectory | Out-Null
        $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
        $startInfo.FileName = $pwshPath
        $startInfo.UseShellExecute = $false
        $startInfo.CreateNoWindow = $true
        foreach ($argument in @(
            "-NoLogo", "-NoProfile", "-File", $workerPath,
            "-SeedExporter", $seedExporterPath,
            "-MaCetspExecutable", $maCetspPath,
            "-LkhExecutable", $lkhPath,
            "-GurobiExecutable", $gurobiPath,
            "-ResultChecker", $resultCheckerPath,
            "-InputPath", $job.InputPath,
            "-InstanceName", $job.Instance,
            "-RunDirectory", $job.RunDirectory,
            "-LeiHao", "$($job.LeiHao)",
            "-PublishedPrecision", "$($job.PublishedPrecision)",
            "-Run", "$($job.Run)",
            "-Seed", "$($job.Seed)",
            "-Repetitions", "$($job.Repetitions)",
            "-PoolSize", "$PoolSize",
            "-Iterations", "$Iterations",
            "-Patience", "$Patience",
            "-InitialPatience", "$InitialPatience",
            "-MaxTime", "$MaxTime",
            "-FitnessWeight", "$FitnessWeight",
            "-MinimumDistance", "$MinimumDistance",
            "-NeighborSize", "$NeighborSize",
            "-ConstructionThreads", "$ConstructionThreads",
            "-SolverThreads", "$SolverThreads"
        )) {
            $startInfo.ArgumentList.Add($argument)
        }
        $process = [System.Diagnostics.Process]::Start($startInfo)
        $active.Add([pscustomobject]@{ Job = $job; Process = $process })
    }

    $finished = @($active | Where-Object { $_.Process.HasExited })
    if ($finished.Count -eq 0) {
        Start-Sleep -Seconds 1
        continue
    }
    foreach ($entry in $finished) {
        $entry.Process.WaitForExit()
        $job = $entry.Job
        if ($entry.Process.ExitCode -eq 0 -and
            (Test-Path -LiteralPath $job.MetricsPath)) {
            $metrics = Get-Content -LiteralPath $job.MetricsPath -Raw |
                ConvertFrom-Json
            $results.Add($metrics)
            $comparison = if ([double]$metrics.final -lt [double]$metrics.lei_hao) {
                "IMPROVED"
            } else {
                "gap={0:F3}%" -f (
                    100.0 * ([double]$metrics.final - [double]$metrics.lei_hao) /
                    [double]$metrics.lei_hao)
            }
            Write-Host (
                "[{0}/{1}] {2} run {3}: {4:F6} ({5}), {6:F1} min" -f
                $results.Count, $total, $job.Instance, $job.Run,
                [double]$metrics.final, $comparison,
                ([double]$metrics.ma_ms / 60000.0))
            Export-Results
        } else {
            $failurePath = Join-Path $job.RunDirectory "failure.log"
            $message = if (Test-Path -LiteralPath $failurePath) {
                Get-Content -LiteralPath $failurePath -Raw
            } else {
                "worker exited with code $($entry.Process.ExitCode)"
            }
            $failures.Add([pscustomobject]@{
                instance = $job.Instance
                run = $job.Run
                message = $message.Trim()
            })
            Write-Warning "$($job.Instance) run $($job.Run) failed"
        }
        $active.Remove($entry) | Out-Null
        $entry.Process.Dispose()
    }
}

$timer.Stop()
Export-Results
Write-Host (
    "Completed {0}/{1} runs in {2:F2} hours; summary: {3}" -f
    $results.Count, $total, $timer.Elapsed.TotalHours, $Summary)
if ($failures.Count -gt 0) {
    $failureOutput = [System.IO.Path]::ChangeExtension($outputPath, ".failures.csv")
    $failures | Export-Csv -LiteralPath $failureOutput -NoTypeInformation
    throw "$($failures.Count) runs failed; details: $failureOutput"
}
$failureOutput = [System.IO.Path]::ChangeExtension($outputPath, ".failures.csv")
Remove-Item -LiteralPath $failureOutput -Force -ErrorAction SilentlyContinue
