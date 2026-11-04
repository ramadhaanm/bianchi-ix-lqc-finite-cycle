param([ValidateSet('K','A')][string]$Quant)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectRoot
$jobs = @(
    @('phase', '20', $Quant, '1000'),
    @('phase', '10', $Quant, '1000'),
    @('phase', '30', $Quant, '1000'),
    @('labels', '20', $Quant, '400'),
    @('long', '20', $Quant, '64'),
    @('phase', '20', $Quant, '200', '0.3'),
    @('phase', '20', $Quant, '200', '0.7')
)
if ($Quant -eq 'K') { $jobs += ,@('basin', '20', $Quant, '100') }
foreach ($job in $jobs) {
    $logName = ($job -join '_') + '.log'
    Write-Output "Starting $($job -join ' ')"
    $ErrorActionPreference = 'Continue' # Julia warnings on stderr are recorded, not PowerShell failures.
    $runner = if ($job[0] -eq 'long') { 'julia/run_cached.jl' } else { 'julia/reanalyze.jl' }
    & julia --threads=8 --project=julia $runner @job 2>&1 |
        Tee-Object -FilePath (Join-Path 'results/revised' $logName)
    $ErrorActionPreference = 'Stop'
    if ($LASTEXITCODE -ne 0) { throw "Run failed: $($job -join ' ')" }
}
Write-Output "ALL RUNS COMPLETE $Quant"
