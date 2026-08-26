param(
    [Parameter(Mandatory = $true)]
    [string]$BenchmarkExecutable,
    [Parameter(Mandatory = $true)]
    [string]$LkhExecutable,
    [Parameter(Mandatory = $true)]
    [string]$GurobiExecutable,
    [int]$Repetitions = 1000,
    [int]$BonusRepetitions = 10000,
    [UInt64]$Seed = 123456789,
    [int]$MaxThreads = 0,
    [int]$Rounds = 1,
    [string]$Output = "results/postprocessing/mennell_postprocess.csv"
)

$ErrorActionPreference = "Stop"
$repository = (Resolve-Path (Join-Path $PSScriptRoot "../..")).Path
$benchmarkExecutablePath = (Resolve-Path $BenchmarkExecutable).Path
$lkhExecutablePath = (Resolve-Path $LkhExecutable).Path
$gurobiExecutablePath = (Resolve-Path $GurobiExecutable).Path
$historical = Import-Csv (Join-Path $repository "results/standalone/comparison.csv")
$historicalByName = @{}
foreach ($row in $historical) {
    $historicalByName[$row.filename] = $row
}

function Read-Field([string]$Text, [string]$Name) {
    $match = [regex]::Match($Text, "(?m)^$([regex]::Escape($Name))=(.+)$")
    if (-not $match.Success) {
        throw "Benchmark output omitted '$Name'. Output:`n$Text"
    }
    return $match.Groups[1].Value.Trim()
}

$results = [System.Collections.Generic.List[object]]::new()
$benchmarkDirectory = Join-Path $repository "data/mennell/reference_tours"
$inputDirectory = Join-Path $repository "data/mennell/instances"
$files = Get-ChildItem -LiteralPath $benchmarkDirectory -Filter "*.txt" |
    Sort-Object Name

foreach ($file in $files) {
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
    Write-Host "[$($results.Count + 1)/$($files.Count)] $($file.Name), repeats=$repeatCount"

    $rawOutput = & $benchmarkExecutablePath $inputPath $repeatCount $Seed `
        $MaxThreads $lkhExecutablePath $gurobiExecutablePath $Rounds 2>&1
    $text = $rawOutput -join "`n"
    if ($LASTEXITCODE -ne 0) {
        throw "Postprocess benchmark failed for $($file.Name):`n$text"
    }

    $benchmarkLine = Get-Content -LiteralPath $file.FullName -TotalCount 1
    $leiHao = [double]([regex]::Match($benchmarkLine, 'value\s*:\s*(.+)$').Groups[1].Value)
    $base = [double](Read-Field $text "base_distance")
    $afterLkh = [double](Read-Field $text "after_lkh_distance")
    $afterSocp = [double](Read-Field $text "after_socp_distance")
    $final = [double](Read-Field $text "final_distance")
    $historicalValue = if ($historicalByName.ContainsKey($file.Name)) {
        [double]$historicalByName[$file.Name].output_value
    } else {
        [double]::NaN
    }

    $results.Add([pscustomobject]@{
        instance = $file.Name
        repetitions = $repeatCount
        lei_hao = $leiHao
        historical_cetsp = $historicalValue
        current_base = $base
        after_lkh = $afterLkh
        after_socp = $afterSocp
        final = $final
        postprocess_improvement_percent = 100.0 * ($base - $final) / $base
        final_gap_to_lei_hao_percent = 100.0 * ($final - $leiHao) / $leiHao
        solver_ms = [long](Read-Field $text "solver_ms")
        lkh_ms = [long](Read-Field $text "lkh_ms")
        socp_ms = [long](Read-Field $text "socp_ms")
        valid = Read-Field $text "valid"
    })
    $outputPath = Join-Path $repository $Output
    $results | Export-Csv -LiteralPath $outputPath -NoTypeInformation
}

Write-Host "Wrote $($results.Count) rows to $Output"
