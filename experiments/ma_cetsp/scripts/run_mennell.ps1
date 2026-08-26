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
    [int]$Repetitions = 1000,
    [int]$BonusRepetitions = 10000,
    [int]$PoolSize = 40,
    [int]$Iterations = 5000,
    [double]$MaxTime = 300,
    [UInt64]$Seed = 123456789,
    [int]$MaxThreads = 0,
    [double]$MinimumDistance = 5.0,
    [string[]]$Instances = @(),
    [switch]$Resume,
    [string]$Output = "results/mennell_hybrid/tables/mennell_hybrid.csv",
    [string]$WorkDirectory = "build-hybrid/ma-mennell"
)

$ErrorActionPreference = "Stop"
$repository = (Resolve-Path (Join-Path $PSScriptRoot "../../..")).Path
$seedExporterPath = (Resolve-Path $SeedExporter).Path
$maCetspPath = (Resolve-Path $MaCetspExecutable).Path
$lkhPath = (Resolve-Path $LkhExecutable).Path
$gurobiPath = (Resolve-Path $GurobiExecutable).Path
$resultCheckerPath = (Resolve-Path $ResultChecker).Path
$outputPath = Join-Path $repository $Output
$workPath = Join-Path $repository $WorkDirectory
New-Item -ItemType Directory -Force -Path $workPath | Out-Null

$gurobiBin = Split-Path -Parent $gurobiPath
$env:PATH = "$gurobiBin;$env:PATH"

function Read-Field([string]$Text, [string]$Name) {
    $match = [regex]::Match($Text, "(?m)^$([regex]::Escape($Name))=(.+)$")
    if (-not $match.Success) {
        throw "Output omitted '$Name':`n$Text"
    }
    return $match.Groups[1].Value.Trim()
}

function Invoke-Checked(
    [string]$Executable,
    [string[]]$Arguments,
    [string]$Description
) {
    $lines = & $Executable @Arguments 2>&1
    $text = $lines -join "`n"
    if ($LASTEXITCODE -ne 0) {
        throw "$Description failed:`n$text"
    }
    return $text
}

$priorRows = Import-Csv (Join-Path $repository "results/postprocessing/mennell_postprocess.csv")
$priorByName = @{}
foreach ($row in $priorRows) {
    $priorByName[$row.instance] = $row
}

$results = [System.Collections.Generic.List[object]]::new()
$completed = @{}
if ($Resume -and (Test-Path -LiteralPath $outputPath)) {
    foreach ($row in (Import-Csv -LiteralPath $outputPath)) {
        $results.Add($row)
        $completed[$row.instance] = $true
    }
}

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

foreach ($file in $files) {
    if ($completed.ContainsKey($file.Name)) {
        Write-Host "Skipping $($file.Name)"
        continue
    }
    $inputPath = Join-Path $inputDirectory $file.Name
    if (-not (Test-Path -LiteralPath $inputPath)) {
        Write-Warning "Skipping $($file.Name): no processed input"
        continue
    }

    $repeatCount = if ($file.Name -eq "bonus1000.txt") {
        $BonusRepetitions
    } else {
        $Repetitions
    }
    $stem = [System.IO.Path]::GetFileNameWithoutExtension($file.Name)
    $instancePath = Join-Path $workPath "$stem.instance"
    $seedPath = Join-Path $workPath "$stem.seeds"
    $resultPath = Join-Path $workPath "$stem.result"
    $lkhTemp = Join-Path $workPath "lkh"
    $depotIndex = (Get-Content -LiteralPath $inputPath | Measure-Object -Line).Lines - 1
    Write-Host "[$($results.Count + 1)/$($files.Count)] $($file.Name), repeats=$repeatCount"

    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    $seedOutput = Invoke-Checked $seedExporterPath @(
        $inputPath,
        $instancePath,
        $seedPath,
        "$repeatCount",
        "$PoolSize",
        "$Seed",
        "$MaxThreads",
        "$depotIndex",
        "$MinimumDistance"
    ) "seed export"
    $timer.Stop()
    $seedMilliseconds = $timer.ElapsedMilliseconds

    $timer.Restart()
    $maOutput = Invoke-Checked $maCetspPath @(
        "--instance_file=$instancePath",
        "--seed_file=$seedPath",
        "--result_file=$resultPath",
        "--lkh_executable=$lkhPath",
        "--lkh_temp_root=$lkhTemp",
        "--pop_size=$PoolSize",
        "--iteration=$Iterations",
        "--max_time=$MaxTime",
        "--seed=$Seed"
    ) "MA-CETSP"
    $timer.Stop()
    $maMilliseconds = $timer.ElapsedMilliseconds

    $checkOutput = Invoke-Checked $resultCheckerPath @(
        $inputPath,
        $resultPath,
        $gurobiPath
    ) "result validation"

    $benchmarkLine = Get-Content -LiteralPath $file.FullName -TotalCount 1
    $leiHao = [double]([regex]::Match(
        $benchmarkLine, 'value\s*:\s*(.+)$').Groups[1].Value)
    $prior = $priorByName[$file.Name]
    $generatedBest = [double](Read-Field $seedOutput "best_distance")
    $maDistance = [double](Read-Field $checkOutput "ma_distance")
    $finalDistance = [double](Read-Field $checkOutput "final_distance")

    $results.Add([pscustomobject]@{
        instance = $file.Name
        repetitions = $repeatCount
        pool_size = $PoolSize
        iterations = $Iterations
        generated_best = $generatedBest
        ma_cetsp = $maDistance
        final = $finalDistance
        lei_hao = $leiHao
        prior_base = if ($prior) { [double]$prior.current_base } else { [double]::NaN }
        prior_final = if ($prior) { [double]$prior.final } else { [double]::NaN }
        hybrid_gain_vs_prior_percent = if ($prior) {
            100.0 * ([double]$prior.final - $finalDistance) / [double]$prior.final
        } else { [double]::NaN }
        gap_to_lei_hao_percent = 100.0 * ($finalDistance - $leiHao) / $leiHao
        seed_ms = $seedMilliseconds
        ma_ms = $maMilliseconds
        socp_ms = [long](Read-Field $checkOutput "socp_ms")
        valid = Read-Field $checkOutput "final_valid"
    })

    $outputParent = Split-Path -Parent $outputPath
    if ($outputParent) {
        New-Item -ItemType Directory -Force -Path $outputParent | Out-Null
    }
    $results | Export-Csv -LiteralPath $outputPath -NoTypeInformation
    Write-Host ("  generated={0:F6} ma={1:F6} final={2:F6} Lei-Hao={3:F6}" -f `
        $generatedBest, $maDistance, $finalDistance, $leiHao)
}

Write-Host "Wrote $($results.Count) rows to $Output"
