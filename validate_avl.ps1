# validate_avl.ps1 -- AVL cross-check for min_drag, min_mass, knee incumbents.
# Usage: powershell -ExecutionPolicy Bypass -File validate_avl.ps1 [-AvlExe <path>] [-Run]
#
# Reads out/<name>_panel.txt (written by aeroanalyzer.exe) and out/<name>.avl, runs AVL
# at the same alpha the panel trimmed to, parses the ST dump, and prints a comparison table.
#
# Tolerances (PASS/FAIL):
#   CL  : INFO only -- panel is viscous+elevon-trimmed, AVL is inviscid clean wing;
#          expected ~25-35% gap. Printed for diagnostics but never gates PASS/FAIL.
#   e   : 10% -- span efficiency (Oswald)
#   Xnp : 3% of MAC -- neutral point location
#
# ponytail: compares CL, e, Xnp only. Add CLa/Cma derivative comparison if Xnp disagrees
# at validation (needs finite-difference alpha sweep; not worth it until there is a failure).

param(
    [string]$AvlExe = "tools/bin/avl352.exe",
    [switch]$Run
)

$ErrorActionPreference = "Stop"

$TOL_E   = 0.10   # 10%
$TOL_XNP = 0.03   # 3% of MAC

# ---- helpers ---------------------------------------------------------------

function Read-KV ($path) {
    # Returns hashtable of "key" -> double from "key = value" lines.
    $h = @{}
    foreach ($line in (Get-Content $path)) {
        if ($line -match '^\s*(\w+)\s*=\s*([0-9Ee.+\-]+)') {
            $h[$Matches[1]] = [double]$Matches[2]
        }
    }
    return $h
}

function Parse-AvlST ($path) {
    # Pull the scalars we need from AVL's "ST" dump (fixed-format text).
    # AVL prints e.g.:
    #   CLtot =  0.38471    CDtot =  0.01234    CDind =  0.00891
    #   e     =  0.9123     ...
    #   Xnp   =  0.17432    ...
    $h = @{ CLtot = $null; CDind = $null; e = $null; Xnp = $null }
    if (-not (Test-Path $path)) { return $h }
    foreach ($line in (Get-Content $path)) {
        if ($line -match 'CLtot\s*=\s*([0-9Ee.+\-]+)') { $h.CLtot = [double]$Matches[1] }
        if ($line -match 'CDind\s*=\s*([0-9Ee.+\-]+)') { $h.CDind  = [double]$Matches[1] }
        if ($line -match '\be\s*=\s*([0-9Ee.+\-]+)')   { $h.e      = [double]$Matches[1] }
        if ($line -match 'Xnp\s*=\s*([0-9Ee.+\-]+)')   { $h.Xnp    = [double]$Matches[1] }
    }
    return $h
}

function Pct-Delta ($panel, $avl) {
    if ([math]::Abs($panel) -lt 1e-12) { return 0.0 }
    return [math]::Abs($avl - $panel) / [math]::Abs($panel) * 100.0
}

function Pass-Fail ($delta_pct, $tol) {
    if ($delta_pct -le ($tol * 100.0)) { return "PASS" } else { return "FAIL" }
}

# ---- pre-flight ------------------------------------------------------------

$avlExeResolved = $AvlExe
if (-not (Test-Path $avlExeResolved)) {
    Write-Error "AVL not found at '$avlExeResolved'. Set -AvlExe or place avl352.exe at tools/bin/."
    exit 1
}

if ($Run) {
    Write-Host "Running aeroanalyzer.exe to generate incumbents..."
    & "build/aeroanalyzer.exe" "config/baseline.cfg"
    if ($LASTEXITCODE -ne 0) { Write-Error "aeroanalyzer.exe failed."; exit 1 }
}

# ---- per-incumbent loop ----------------------------------------------------

$cases    = @("min_drag", "min_mass", "knee")
$rows     = @()
$anyFail  = $false
$tmpDeck  = "$env:TEMP\avl_deck_$PID.txt"

foreach ($name in $cases) {
    $stem   = "out/$name"
    $avlFile   = "$stem.avl"
    $panelFile = "${stem}_panel.txt"
    $stFile    = "${stem}_st.txt"

    if (-not (Test-Path $avlFile) -or -not (Test-Path $panelFile)) {
        Write-Warning "${name}: missing $avlFile or $panelFile -- run aeroanalyzer.exe first (or use -Run)."
        $rows += [pscustomobject]@{
            Case="$name"; CL_panel="?"; CL_avl="?"; CL_d="?"; CL_pf="SKIP"
            e_panel="?";  e_avl="?";   e_d="?";  e_pf="SKIP"
            Xnp_panel="?";Xnp_avl="?";Xnp_d="?";Xnp_pf="SKIP"
        }
        continue
    }

    $p = Read-KV $panelFile
    $alpha = $p["alpha_deg"]
    $mac   = $p["mac"]

    # Write per-case AVL command deck.
    # Delete stale ST file first -- AVL prompts for overwrite on some versions.
    if (Test-Path $stFile) { Remove-Item $stFile -Force }

    @(
        "load $avlFile",
        "oper",
        "a a $alpha",
        "x",
        "st",
        "$stFile",
        "",
        "quit"
    ) | Set-Content $tmpDeck -Encoding ASCII

    $avlOut = "$env:TEMP\avl_stdout_$PID.txt"
    Start-Process -FilePath $avlExeResolved -RedirectStandardInput $tmpDeck `
        -RedirectStandardOutput $avlOut -Wait -NoNewWindow | Out-Null
    if (Test-Path $avlOut) { Remove-Item $avlOut -Force }

    $a = Parse-AvlST $stFile

    if ($null -eq $a.CLtot) {
        Write-Warning "${name}: AVL ST parse failed -- check $stFile manually."
        $rows += [pscustomobject]@{
            Case="$name"; CL_panel="?"; CL_avl="?"; CL_d="?"; CL_pf="ERR"
            e_panel="?";  e_avl="?";   e_d="?";  e_pf="ERR"
            Xnp_panel="?";Xnp_avl="?";Xnp_d="?";Xnp_pf="ERR"
        }
        $anyFail = $true
        continue
    }

    $cl_d  = Pct-Delta $p["CL"]   $a.CLtot
    $e_d   = Pct-Delta $p["e"]    $a.e
    # Xnp tolerance is absolute (3% MAC), converted to a %delta of MAC for comparison.
    $xnp_abs_delta = [math]::Abs($a.Xnp - $p["x_np"])
    $xnp_pct_mac   = if ($mac -gt 1e-9) { $xnp_abs_delta / $mac * 100.0 } else { 0.0 }

    $cl_pf  = "INFO"   # viscous+trimmed vs inviscid clean -- always differs, never gates
    $e_pf   = Pass-Fail $e_d   $TOL_E
    # Xnp PASS when abs error < 3% MAC.
    $xnp_pf = if ($xnp_abs_delta -le ($TOL_XNP * $mac)) { "PASS" } else { "FAIL" }

    if ($e_pf -eq "FAIL" -or $xnp_pf -eq "FAIL") { $anyFail = $true }

    $rows += [pscustomobject]@{
        Case       = $name
        CL_panel   = "{0:F4}" -f $p["CL"]
        CL_avl     = "{0:F4}" -f $a.CLtot
        CL_d       = "{0:F1}%" -f $cl_d
        CL_pf      = $cl_pf
        e_panel    = "{0:F4}" -f $p["e"]
        e_avl      = "{0:F4}" -f $a.e
        e_d        = "{0:F1}%" -f $e_d
        e_pf       = $e_pf
        Xnp_panel  = "{0:F4}" -f $p["x_np"]
        Xnp_avl    = "{0:F4}" -f $a.Xnp
        Xnp_d      = "{0:F4}m" -f $xnp_abs_delta
        Xnp_pf     = $xnp_pf
    }
}

if (Test-Path $tmpDeck) { Remove-Item $tmpDeck -Force }

# ---- print table -----------------------------------------------------------

Write-Host ""
Write-Host "  Panel vs AVL cross-check  (gates: e 10%, Xnp 3% MAC  |  CL: INFO only)"
Write-Host ("  {0,-10} {1,6} {2,6} {3,6} {4,-4}  {5,6} {6,6} {7,6} {8,-4}  {9,6} {10,6} {11,8} {12,-4}" -f `
    "case","CL_p","CL_a","dCL","","e_p","e_a","de","","Xnp_p","Xnp_a","dXnp","")
Write-Host ("  " + ("-" * 85))
foreach ($row in $rows) {
    Write-Host ("  {0,-10} {1,6} {2,6} {3,6} {4,-4}  {5,6} {6,6} {7,6} {8,-4}  {9,6} {10,6} {11,8} {12,-4}" -f `
        $row.Case, $row.CL_panel, $row.CL_avl, $row.CL_d, $row.CL_pf,
        $row.e_panel, $row.e_avl, $row.e_d, $row.e_pf,
        $row.Xnp_panel, $row.Xnp_avl, $row.Xnp_d, $row.Xnp_pf)
}
Write-Host ""

if ($anyFail) {
    Write-Host "  RESULT: FAIL -- one or more tolerances exceeded."
    exit 1
} else {
    Write-Host "  RESULT: PASS -- all cases within tolerance."
    exit 0
}
