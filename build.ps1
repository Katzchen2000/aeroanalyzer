<#
  build.ps1 - build AeroAnalyzer Pro with MSVC directly (no CMake required).

  This machine has the VS2022 Build Tools (cl.exe) but no CMake/Ninja, so this
  script drives cl.exe through vcvars64.bat. After 'winget install Kitware.CMake'
  you can use the CMake flow in README.md instead.

  Usage (from the project root):
    powershell -ExecutionPolicy Bypass -File build.ps1            # build app
    powershell -ExecutionPolicy Bypass -File build.ps1 -Test      # + build/run tests
    powershell -ExecutionPolicy Bypass -File build.ps1 -OpenMP    # parallel GA
    powershell -ExecutionPolicy Bypass -File build.ps1 -Run       # build + run app

  NOTE: keep this file pure ASCII. PowerShell 5.1 reads non-BOM files as ANSI,
  and a UTF-8 em-dash decodes to a curly quote that PS treats as a string
  delimiter, breaking the parse.
#>
param(
  [switch]$Test,
  [switch]$Run,
  [switch]$OpenMP,
  [switch]$Gen
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$build = Join-Path $root "build"
New-Item -ItemType Directory -Force -Path $build | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $build "build_app") | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $build "build_tests") | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $build "build_gen") | Out-Null

# ---- locate vcvars64.bat ----
$vcvars = $null
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (Test-Path $vswhere) {
  $inst = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
  if ($inst) { $vcvars = Join-Path $inst "VC\Auxiliary\Build\vcvars64.bat" }
}
if (-not $vcvars -or -not (Test-Path $vcvars)) {
  $vcvars = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
}
if (-not (Test-Path $vcvars)) {
  Write-Error "vcvars64.bat not found. Install the VS Build Tools C++ workload."
}
Write-Host "Using: $vcvars" -ForegroundColor DarkGray

$omp = ""
if ($OpenMP) { $omp = "/openmp" }
# /fp:precise (NOT /fp:fast): keep IEEE semantics so the NaN/Inf watchdog checks
# survive. See README, "Numerical strategy".
$common = "/nologo /std:c++17 /O2 /fp:precise /EHsc /W3 /utf-8 $omp /I `"$root\include`""

$installerDir = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer"

function Invoke-Cl([string]$clArgs, [string]$label) {
  # vcvars64.bat shells out to vswhere; make sure it is resolvable on PATH.
  $cmd = "set `"PATH=%PATH%;$installerDir`" && call `"$vcvars`" && cd /d `"$build`" && cl $clArgs"
  Write-Host "==> $label" -ForegroundColor Cyan
  & cmd.exe /c $cmd 2>&1 | Tee-Object -FilePath (Join-Path $build "build.log")
  if ($LASTEXITCODE -ne 0) { Write-Error "$label failed (exit $LASTEXITCODE). See build\build.log" }
}

# ---- application ----
$appArgs = "$common `"$root\src\*.cpp`" `"$root\app\main.cpp`" /Fe:aeroanalyzer.exe /Fobuild_app/"
Invoke-Cl $appArgs "Building aeroanalyzer.exe"
Write-Host "Built: $build\aeroanalyzer.exe" -ForegroundColor Green

# ---- offline surrogate generator ----
$genArgs = "$common `"$root\src\*.cpp`" `"$root\tools\build_surrogate.cpp`" /Fe:build_surrogate.exe /Fobuild_gen/"
Invoke-Cl $genArgs "Building build_surrogate.exe"
Write-Host "Built: $build\build_surrogate.exe" -ForegroundColor Green

# ---- tests ----
if ($Test) {
  $testArgs = "$common /I `"$root\tests`" `"$root\src\*.cpp`" `"$root\tests\*.cpp`" /Fe:unit_tests.exe /Fobuild_tests/"
  Invoke-Cl $testArgs "Building unit_tests.exe"
  Write-Host "Built: $build\unit_tests.exe" -ForegroundColor Green
  Write-Host "==> Running tests" -ForegroundColor Cyan
  & (Join-Path $build "unit_tests.exe")
  if ($LASTEXITCODE -ne 0) { Write-Error "TESTS FAILED (exit $LASTEXITCODE)" }
  Write-Host "Tests passed." -ForegroundColor Green
}

# ---- generate viscous surrogate ----
if ($Gen) {
  Write-Host "==> Generating viscous surrogate (config\baseline.cfg)" -ForegroundColor Cyan
  Push-Location $root
  try { & (Join-Path $build "build_surrogate.exe") "config\baseline.cfg" }
  finally { Pop-Location }
}

# ---- run app ----
if ($Run) {
  Write-Host "==> Running aeroanalyzer (config\baseline.cfg)" -ForegroundColor Cyan
  Push-Location $root
  try { & (Join-Path $build "aeroanalyzer.exe") "config\baseline.cfg" }
  finally { Pop-Location }
}
