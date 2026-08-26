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
    [Parameter(Mandatory = $true)]
    [string]$InputPath,
    [Parameter(Mandatory = $true)]
    [string]$InstanceName,
    [Parameter(Mandatory = $true)]
    [string]$RunDirectory,
    [Parameter(Mandatory = $true)]
    [double]$LeiHao,
    [Parameter(Mandatory = $true)]
    [int]$PublishedPrecision,
    [Parameter(Mandatory = $true)]
    [int]$Run,
    [Parameter(Mandatory = $true)]
    [int]$Seed,
    [int]$Repetitions = 1000,
    [int]$PoolSize = 20,
    [int]$Iterations = 5000,
    [int]$Patience = 500,
    [int]$InitialPatience = 1500,
    [double]$MaxTime = 36000,
    [double]$FitnessWeight = 0.96,
    [double]$MinimumDistance = 5.0,
    [int]$NeighborSize = 50,
    [int]$ConstructionThreads = 1,
    [int]$SolverThreads = 1
)

$ErrorActionPreference = "Stop"

function Read-Field([string]$Text, [string]$Name) {
    $match = [regex]::Match($Text, "(?m)^$([regex]::Escape($Name))=(.+)$")
    if (-not $match.Success) {
        throw "Output omitted '$Name'"
    }
    return $match.Groups[1].Value.Trim()
}

function Read-SummaryField([string]$Text, [string]$Name) {
    $match = [regex]::Match(
        $Text,
        "(?m)^\[SUMMARY\].*?\b$([regex]::Escape($Name)):\s*([^\s]+)")
    if (-not $match.Success) {
        throw "MA-CETSP summary omitted '$Name'"
    }
    return $match.Groups[1].Value.Trim()
}

function Invoke-Captured(
    [string]$Executable,
    [string[]]$Arguments,
    [string]$LogPath,
    [string]$Description
) {
    $lines = & $Executable @Arguments 2>&1
    $exitCode = $LASTEXITCODE
    $lines | Set-Content -LiteralPath $LogPath -Encoding utf8
    $text = $lines -join "`n"
    if ($exitCode -ne 0) {
        throw "$Description failed with exit code $exitCode"
    }
    return $text
}

try {
    New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
    $gurobiBin = Split-Path -Parent $GurobiExecutable
    $env:PATH = "$gurobiBin;$env:PATH"

    $instancePath = Join-Path $RunDirectory "instance.txt"
    $seedPath = Join-Path $RunDirectory "population.seeds"
    $maResultPath = Join-Path $RunDirectory "ma.result"
    $finalResultPath = Join-Path $RunDirectory "final.result"
    $lkhTemp = Join-Path $RunDirectory "lkh"
    $depotIndex = (Get-Content -LiteralPath $InputPath | Measure-Object -Line).Lines - 1

    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $seedOutput = Invoke-Captured $SeedExporter @(
        $InputPath,
        $instancePath,
        $seedPath,
        "$Repetitions",
        "$PoolSize",
        "$Seed",
        "$ConstructionThreads",
        "$depotIndex",
        "$MinimumDistance"
    ) (Join-Path $RunDirectory "seed.log") "seed export"
    $timer.Stop()
    $seedMilliseconds = $timer.ElapsedMilliseconds

    $timer.Restart()
    $maOutput = Invoke-Captured $MaCetspExecutable @(
        "--instance_file=$instancePath",
        "--seed_file=$seedPath",
        "--result_file=$maResultPath",
        "--lkh_executable=$LkhExecutable",
        "--lkh_temp_root=$lkhTemp",
        "--pop_size=$PoolSize",
        "--iteration=$Iterations",
        "--patience=$Patience",
        "--initial_patience=$InitialPatience",
        "--max_time=$MaxTime",
        "--fit_beta=$FitnessWeight",
        "--dist_th=$MinimumDistance",
        "--neighbor_size=$NeighborSize",
        "--solver_threads=$SolverThreads",
        "--seed=$Seed"
    ) (Join-Path $RunDirectory "ma.log") "MA-CETSP"
    $timer.Stop()
    $maMilliseconds = $timer.ElapsedMilliseconds

    $checkOutput = Invoke-Captured $ResultChecker @(
        $InputPath,
        $maResultPath,
        $GurobiExecutable,
        $finalResultPath,
        "$SolverThreads",
        $instancePath
    ) (Join-Path $RunDirectory "check.log") "result validation"

    $finalDistance = [double](Read-Field $checkOutput "final_distance")
    $metrics = [ordered]@{
        instance = $InstanceName
        run = $Run
        seed = $Seed
        repetitions = $Repetitions
        pool_size = $PoolSize
        iterations = $Iterations
        patience = $Patience
        initial_patience = $InitialPatience
        generated_best = [double](Read-Field $seedOutput "best_distance")
        ma_cetsp = [double](Read-Field $checkOutput "ma_distance")
        final = $finalDistance
        lei_hao = $LeiHao
        published_precision = $PublishedPrecision
        improvement_vs_lei_hao_percent = 100.0 * ($LeiHao - $finalDistance) / $LeiHao
        seed_ms = $seedMilliseconds
        ma_ms = $maMilliseconds
        socp_ms = [long](Read-Field $checkOutput "socp_ms")
        ma_best_time = [double](Read-SummaryField $maOutput "best_time")
        ma_total_time = [double](Read-SummaryField $maOutput "total_time")
        ma_valid = Read-Field $checkOutput "ma_valid"
        final_valid = Read-Field $checkOutput "final_valid"
        ma_result = $maResultPath
        final_result = $finalResultPath
    }
    $temporaryMetrics = Join-Path $RunDirectory "metrics.json.tmp"
    $metrics | ConvertTo-Json | Set-Content -LiteralPath $temporaryMetrics -Encoding utf8
    Move-Item -LiteralPath $temporaryMetrics `
        -Destination (Join-Path $RunDirectory "metrics.json") -Force
    Remove-Item -LiteralPath (Join-Path $RunDirectory "failure.log") `
        -Force -ErrorAction SilentlyContinue
    exit 0
} catch {
    ($_ | Out-String) | Set-Content `
        -LiteralPath (Join-Path $RunDirectory "failure.log") -Encoding utf8
    exit 1
}
