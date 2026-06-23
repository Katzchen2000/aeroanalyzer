<#
  build_mingw.ps1 - build AeroAnalyzer Pro with the MSYS2 MinGW GNU toolchain
  (g++ / gfortran). This is the primary build going forward because it can link
  the XFOIL Fortran library and the C++ code with one toolchain (no ABI clash).
  build.ps1 (MSVC) remains as an optional XFOIL-free fallback.

  Usage (from the project root):
    powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test
    powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Run

  Keep this file pure ASCII (PowerShell 5.1 reads non-BOM files as ANSI; a UTF-8
  em-dash becomes a curly quote that breaks the parse).
#>
param(
  [switch]$Test,
  [switch]$Run,
  [switch]$Native
)

$ErrorActionPreference = "Stop"
$root  = $PSScriptRoot
$build = Join-Path $root "build"
New-Item -ItemType Directory -Force -Path $build | Out-Null

# ---- toolchain (MSYS2 MinGW64) ----
$mingwBin = "C:\msys64\mingw64\bin"
$msysBin  = "C:\msys64\usr\bin"
$gpp = Join-Path $mingwBin "g++.exe"
$gfortran = Join-Path $mingwBin "gfortran.exe"
$ar = Join-Path $mingwBin "ar.exe"
if (-not (Test-Path $gpp)) { Write-Error "g++ not found at $gpp (install MSYS2 mingw-w64 g++)" }
# Put the toolchain on PATH so built exes find libstdc++/libgomp/etc. at runtime.
$env:PATH = "$mingwBin;$msysBin;$env:PATH"
Write-Host "Using: $gpp" -ForegroundColor DarkGray

# ---- Eigen (header-only; no install needed) ----
$eigen = "C:\Users\kadan\OneDrive\c++ libraries\eigen-5.0.0"

# ---- common compile flags (NO -ffast-math: keeps NaN/Inf watchdogs alive) ----
# EIGEN_DONT_PARALLELIZE: Eigen's per-thread VLM solves run inside the OpenMP GA
# loop; let OpenMP own the parallelism (avoid nested oversubscription).
$common = @('-std=c++17','-O2','-fno-math-errno','-fopenmp',
            '-DEIGEN_DONT_PARALLELIZE',
            '-static-libgcc','-static-libstdc++',
            '-I', (Join-Path $root 'include'))
if (Test-Path (Join-Path $eigen 'Eigen\Dense')) {
  $common += @('-I', $eigen)
} else {
  Write-Error "Eigen not found at '$eigen' (Eigen/Dense missing). Fix the path in build_mingw.ps1."
}
if ($Native) { $common += '-march=native' }

$srcs = (Get-ChildItem (Join-Path $root 'src') -Filter *.cpp).FullName

function Build([string[]]$extra, [string]$outName, [string]$label) {
  Write-Host "==> $label" -ForegroundColor Cyan
  $argList = $common + $srcs + $extra + @('-o', (Join-Path $build $outName))
  & $gpp @argList
  if ($LASTEXITCODE -ne 0) { Write-Error "$label failed (exit $LASTEXITCODE)" }
  Write-Host "Built: $build\$outName" -ForegroundColor Green
}

# ---- application ----
Build @((Join-Path $root 'app\main.cpp')) "aeroanalyzer.exe" "Building aeroanalyzer.exe"

# ---- tests ----
if ($Test) {
  $tests = (Get-ChildItem (Join-Path $root 'tests') -Filter *.cpp).FullName
  Write-Host "==> Building unit_tests.exe" -ForegroundColor Cyan
  $argList = $common + @('-I', (Join-Path $root 'tests')) + $srcs + $tests + @('-o', (Join-Path $build 'unit_tests.exe'))
  & $gpp @argList
  if ($LASTEXITCODE -ne 0) { Write-Error "unit_tests build failed (exit $LASTEXITCODE)" }
  Write-Host "Built: $build\unit_tests.exe" -ForegroundColor Green
  Write-Host "==> Running tests" -ForegroundColor Cyan
  & (Join-Path $build 'unit_tests.exe')
  if ($LASTEXITCODE -ne 0) { Write-Error "TESTS FAILED (exit $LASTEXITCODE)" }
  Write-Host "Tests passed." -ForegroundColor Green
}

# ---- run optimizer ----
if ($Run) {
  Write-Host "==> Running aeroanalyzer (config\baseline.cfg)" -ForegroundColor Cyan
  Push-Location $root
  try { & (Join-Path $build 'aeroanalyzer.exe') "config\baseline.cfg" }
  finally { Pop-Location }
}
