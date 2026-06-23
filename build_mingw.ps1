<#
  build_mingw.ps1 - build AeroAnalyzer Pro with the MSYS2 MinGW GNU toolchain
  (g++ / gfortran). This is the primary build going forward because it can link
  the XFOIL Fortran library and the C++ code with one toolchain (no ABI clash).
  build.ps1 (MSVC) remains as an optional XFOIL-free fallback.

  Usage (from the project root):
    powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test
    powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Gen
    powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Run

  Keep this file pure ASCII (PowerShell 5.1 reads non-BOM files as ANSI; a UTF-8
  em-dash becomes a curly quote that breaks the parse).
#>
param(
  [switch]$Test,
  [switch]$Run,
  [switch]$Gen,
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

# ---- XFOIL Fortran static library (optional accurate surrogate engine) ----
# Compiles the real XFOIL 6.99 analysis core (tools/bin/Xfoil699src/src) plus
# our C-callable wrapper (tools/xfoil_lib/xfwrap.f) and the no-op plot/menu
# stubs into build/libxfoil.a. build_surrogate.exe links this for
# surrogate_mode = xfoil_lib. Skipped (with a warning) if gfortran is absent,
# in which case xfoil_lib mode falls back to native/synthetic at runtime.
$xfoilLib = ""
function Build-XfoilLib {
  if (-not (Test-Path $gfortran)) {
    Write-Host "gfortran not found at $gfortran -- skipping libxfoil.a (surrogate_mode=xfoil_lib unavailable)" -ForegroundColor Yellow
    return ""
  }
  $xsrc = Join-Path $root 'tools\bin\Xfoil699src\src'
  $lib  = Join-Path $root 'tools\xfoil_lib'
  if (-not (Test-Path (Join-Path $xsrc 'xfoil.f'))) {
    Write-Host "XFOIL sources not found at $xsrc -- skipping libxfoil.a" -ForegroundColor Yellow
    return ""
  }
  $obj = Join-Path $build 'xfoil_obj'
  New-Item -ItemType Directory -Force -Path $obj | Out-Null
  $fflags = @('-O2','-std=legacy','-fallow-argument-mismatch','-I', $xsrc)
  # Verified-sufficient analysis-core source set (no plot/menu/design files).
  $core = @('aread','naca','sort','spline','userio','xbl','xblsys','xgeom',
            'xoper','xpanel','xsolve','xutils')
  Write-Host "==> Building libxfoil.a (gfortran)" -ForegroundColor Cyan
  foreach ($f in $core) {
    & $gfortran @fflags -c (Join-Path $xsrc "$f.f") -o (Join-Path $obj "$f.o")
    if ($LASTEXITCODE -ne 0) { Write-Error "gfortran failed on $f.f" }
  }
  # xfoil.f defines PROGRAM XFOIL (collides with the C++ main); patch it to a
  # plain SUBROUTINE in a build-dir copy so the analysis routines still link.
  $patched = Join-Path $obj 'xfoil_nomain.f'
  (Get-Content (Join-Path $xsrc 'xfoil.f')) `
    -replace '^      PROGRAM XFOIL', '      SUBROUTINE XF_NOMAIN' |
    Set-Content -Encoding ascii $patched
  & $gfortran @fflags -c $patched              -o (Join-Path $obj 'xfoil.o')
  if ($LASTEXITCODE -ne 0) { Write-Error "gfortran failed on patched xfoil.f" }
  & $gfortran @fflags -c (Join-Path $lib 'xfwrap.f')  -o (Join-Path $obj 'xfwrap.o')
  if ($LASTEXITCODE -ne 0) { Write-Error "gfortran failed on xfwrap.f" }
  & $gfortran @fflags -c (Join-Path $lib 'xfstubs.f') -o (Join-Path $obj 'xfstubs.o')
  if ($LASTEXITCODE -ne 0) { Write-Error "gfortran failed on xfstubs.f" }
  $outLib = Join-Path $build 'libxfoil.a'
  Remove-Item -Force $outLib -ErrorAction SilentlyContinue
  $objs = (Get-ChildItem $obj -Filter *.o).FullName
  & $ar rcs $outLib @objs
  if ($LASTEXITCODE -ne 0) { Write-Error "ar failed building libxfoil.a" }
  Write-Host "Built: $outLib" -ForegroundColor Green
  return $outLib
}

function Build([string[]]$extra, [string]$outName, [string]$label) {
  Write-Host "==> $label" -ForegroundColor Cyan
  $argList = $common + $srcs + $extra + @('-o', (Join-Path $build $outName))
  & $gpp @argList
  if ($LASTEXITCODE -ne 0) { Write-Error "$label failed (exit $LASTEXITCODE)" }
  Write-Host "Built: $build\$outName" -ForegroundColor Green
}

# ---- application + offline surrogate generator ----
# The runtime optimizer never links XFOIL (stays GPL-free / dependency-free).
Build @((Join-Path $root 'app\main.cpp')) "aeroanalyzer.exe" "Building aeroanalyzer.exe"

# Only the offline surrogate generator links libxfoil.a (for xfoil_lib mode).
$xfoilLib = Build-XfoilLib
$bsExtra = @((Join-Path $root 'tools\build_surrogate.cpp'))
if ($xfoilLib -ne "") {
  # -DHAVE_XFOIL_LIB enables the xfoil_lib engine + the self-spawning worker
  # pool; the gfortran runtime libs satisfy libxfoil.a's references.
  $bsExtra = @('-DHAVE_XFOIL_LIB') + $bsExtra + @($xfoilLib, '-lgfortran', '-lquadmath')
}
Build $bsExtra "build_surrogate.exe" "Building build_surrogate.exe"

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

# ---- generate viscous surrogate ----
if ($Gen) {
  Write-Host "==> Generating viscous surrogate" -ForegroundColor Cyan
  Push-Location $root
  try { & (Join-Path $build 'build_surrogate.exe') "config\baseline.cfg" }
  finally { Pop-Location }
}

# ---- run optimizer ----
if ($Run) {
  Write-Host "==> Running aeroanalyzer (config\baseline.cfg)" -ForegroundColor Cyan
  Push-Location $root
  try { & (Join-Path $build 'aeroanalyzer.exe') "config\baseline.cfg" }
  finally { Pop-Location }
}
