$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectRoot
do {
    $counts = @()
    foreach ($quant in @('K','A')) {
        $csv = Join-Path $projectRoot "results/revised/long_L20_${quant}_N64_on0.5_b300000.csv"
        $stream = [System.IO.FileStream]::new($csv,[System.IO.FileMode]::Open,[System.IO.FileAccess]::Read,[System.IO.FileShare]::ReadWrite)
        $reader = [System.IO.StreamReader]::new($stream)
        try {
            $lines = 0
            while ($null -ne $reader.ReadLine()) { $lines++ }
            $counts += [Math]::Max(0,$lines - 2)
        } finally { $reader.Dispose(); $stream.Dispose() }
    }
    Write-Output "$(Get-Date -Format 'HH:mm:ss') six-cycle trajectory jobs: K $($counts[0])/256; A $($counts[1])/256"
    if (($counts | Measure-Object -Minimum).Minimum -lt 256) { Start-Sleep -Seconds 45 }
} while (($counts | Measure-Object -Minimum).Minimum -lt 256)
Start-Sleep -Seconds 2

foreach ($script in @('analysis/reanalyze.py','analysis/check_data.py','analysis/test_reanalysis.py')) {
    Write-Output "Running $script"
    $ErrorActionPreference = 'Continue'
    & .venv/Scripts/python $script 2>&1 | Tee-Object -FilePath ('results/revised/final_' + [System.IO.Path]::GetFileNameWithoutExtension($script) + '.log')
    $ErrorActionPreference = 'Stop'
    if ($LASTEXITCODE -ne 0) { throw "Failed: $script" }
}
Write-Output 'Analysis, data audit, and statistical tests complete. Visual review and interpretation still required.'
