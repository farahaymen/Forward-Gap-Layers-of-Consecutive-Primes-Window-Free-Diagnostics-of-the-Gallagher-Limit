<#
.SYNOPSIS
    Forward-Gap Partition census -- full build on Windows.

.DESCRIPTION
    The Windows equivalent of the Makefile. Nothing generated is stored in this
    repository: the census, the tables, the figures and every number quoted in
    the manuscript are produced from source by this script, so the paper cannot
    disagree with the computation.

    Run it from the repository root, in a PowerShell terminal inside VS Code.

.PARAMETER Target
    all        census, checks, analysis and PDF (the default)
    census     just the layer census
    check      the 45 offline accuracy checks
    analysis   tables, figures and macros
    paper      the PDF
    full       as 'all', but verifies all 20 checkpoints (~2 min, ~1 GB)
    lean       check the Lean proofs (needs elan; downloads Mathlib)
    clean      remove generated files, keep the census
    distclean  remove generated files and the census

.PARAMETER Verify
    Upper X for the independent reimplementation. 0 disables it.
    10000000000 checks all twenty checkpoints (what 'full' sets).

.PARAMETER NoCompiler
    Ignore any C compiler and use the pure-Python census generator.

.EXAMPLE
    .\build.ps1
.EXAMPLE
    .\build.ps1 full
.EXAMPLE
    .\build.ps1 analysis -Verify 500000000
#>

[CmdletBinding()]
param(
    [ValidateSet('all','census','check','analysis','paper','full','lean','clean','distclean','help')]
    [string]$Target = 'all',
    [long]$Verify = 2100000000,
    [switch]$NoCompiler,
    [string]$Python = ''
)

# 'Continue', not 'Stop'. In Windows PowerShell a native command that writes to
# stderr becomes a TERMINATING error when the preference is 'Stop', even with
# 2>$null in place. Every tool this script drives writes to stderr in normal
# operation: the census sieve logs progress there, latexmk chatters, and the
# Microsoft Store python stub prints its advertisement there. Success is
# therefore judged by exit code alone, which every call below checks explicitly.
$ErrorActionPreference = 'Continue'
if (Test-Path Variable:PSNativeCommandUseErrorActionPreference) {
    $PSNativeCommandUseErrorActionPreference = $false
}
Set-StrictMode -Version Latest

# 10^10 plus a margin: a prime's layer is fixed by the NEXT prime, so a sieve
# stopping exactly at 10^10 cannot classify its own last prime.
$Limit = '10000001000'
# pi(10^10), so the checkpoint schedule can be built before the run starts.
$Pi    = '455052511'

$Root    = $PSScriptRoot
$Census  = Join-Path $Root 'data\fgp_layers.csv'
$Macros  = Join-Path $Root 'paper\tables\macros.tex'
$Pdf     = Join-Path $Root 'paper\main.pdf'
$SieveC  = Join-Path $Root 'code\fgp_sieve.c'
$SieveEx = Join-Path $Root 'code\fgp_sieve.exe'

function Write-Step([string]$m) { Write-Host "==> $m" -ForegroundColor Cyan }
function Write-Ok  ([string]$m) { Write-Host $m -ForegroundColor Green }
function Fail      ([string]$m) { Write-Host $m -ForegroundColor Red; exit 1 }

# ---------------------------------------------------------------- interpreter
# Windows installs Python as 'python' or exposes the launcher as 'py'. Prefer a
# real python3 if one is on PATH, then 'python', then 'py -3'.
# Two script-scope variables rather than one array: PowerShell unwraps a
# single-element array on return, so a function returning @($exe) hands back a
# bare string and any later .Count on it fails under StrictMode.
$script:PyExe = $null
$script:PyPre = @()
$script:PyInfo = $null

# Packages the pipeline needs. sympy is used only by the self-test.
$script:NeedModules = @('numpy','scipy','pandas','matplotlib','sympy')

# Probe source. It is written to a temporary .py file and run as a script
# rather than passed with -c. Windows PowerShell 5.1 builds a single command
# line for native programs and mangles an argument that contains newlines and
# double quotes, which is exactly what an inline probe looks like; every
# interpreter then appears to fail. A file path has neither problem.
$script:ProbeFile = $null

function New-ProbeFile {
    if ($script:ProbeFile -and (Test-Path $script:ProbeFile)) { return $script:ProbeFile }
    $code = @'
import sys, importlib.util
need = ["numpy", "scipy", "pandas", "matplotlib", "sympy"]
print("%d.%d" % (sys.version_info[0], sys.version_info[1]))
print(",".join(m for m in need if importlib.util.find_spec(m) is None))
'@
    $f = Join-Path ([System.IO.Path]::GetTempPath()) ("fgp_probe_$([guid]::NewGuid().ToString('N')).py")
    Set-Content -LiteralPath $f -Value $code -Encoding ASCII
    $script:ProbeFile = $f
    return $f
}

# Every interpreter examined, with the reason it was accepted or rejected, so a
# failure to find one explains itself instead of just saying none was found.
$script:PyTried = @()

# Probes one interpreter. Returns $null if it is not a usable Python 3,
# otherwise a record carrying its version and which required packages it lacks.
function Get-PythonInfo([string]$exe, [string[]]$pre) {
    $label = (@($exe) + $pre) -join ' '
    # Windows ships stubs called python.exe and python3.exe under
    # ...\AppData\Local\Microsoft\WindowsApps that only advertise the
    # Microsoft Store. Never execute those; they are not interpreters.
    if ($exe -match '\\WindowsApps\\') {
        $script:PyTried += "$label  -- skipped, Microsoft Store stub"
        return $null
    }
    try {
        $probe = New-ProbeFile
        $out = @(& $exe @($pre + @($probe)) 2>$null)
        if ($LASTEXITCODE -ne 0) {
            $script:PyTried += "$label  -- did not run (exit code $LASTEXITCODE)"
            return $null
        }
        if ($out.Count -lt 1) {
            $script:PyTried += "$label  -- ran but printed nothing"
            return $null
        }
        $ver = "$($out[0])".Trim()
        if (-not $ver.StartsWith('3.')) {
            $script:PyTried += "$label  -- reported version '$ver', not Python 3"
            return $null
        }
        $missRaw = if ($out.Count -ge 2) { "$($out[1])".Trim() } else { '' }
        # Count it here, while it is still a real array. Assigning a
        # single-element array to a property unwraps it to a bare string, and
        # .Count on that then fails under StrictMode.
        [string[]]$miss = @()
        if ($missRaw) { $miss = [string[]]($missRaw -split ',') }
        $note = if ($miss.Count -eq 0) { 'has every required package' }
                else { "missing $($miss -join ', ')" }
        $script:PyTried += "$label  -- Python $ver, $note"
        return [pscustomobject]@{
            Exe = $exe; Pre = $pre; Version = $ver
            Missing = $miss; MissingCount = $miss.Count
        }
    } catch {
        $script:PyTried += "$label  -- could not be started ($($_.Exception.Message))"
        return $null
    }
}

function Resolve-Python {
    if ($Python) {
        $info = Get-PythonInfo $Python @()
        if (-not $info) { Fail "the interpreter given by -Python did not run as Python 3: $Python" }
        $script:PyExe = $info.Exe; $script:PyPre = $info.Pre; $script:PyInfo = $info
        return
    }

    # Gather every candidate rather than stopping at the first that starts. A
    # machine can easily carry several Pythons -- an Anaconda install, an MSYS2
    # one, a Store stub -- and only some of them have the scientific stack. The
    # one that can actually run the pipeline is the one to use, so candidates
    # are ranked by how many required packages they are missing.
    $cands = @()
    foreach ($c in @('python3','python')) {
        foreach ($p in @(Get-Command $c -All -ErrorAction SilentlyContinue)) {
            $info = Get-PythonInfo $p.Source @()
            if ($info) { $cands += $info }
        }
    }
    foreach ($p in @(Get-Command 'py' -All -ErrorAction SilentlyContinue)) {
        $info = Get-PythonInfo $p.Source @('-3')
        if ($info) { $cands += $info }
    }

    if ($cands.Count -eq 0) {
        Write-Host "Interpreters examined:" -ForegroundColor Yellow
        foreach ($t in $script:PyTried) { Write-Host "    $t" -ForegroundColor Yellow }
        Write-Host ""
        Fail @"
No working Python 3 was found.

If the terminal advertised the Microsoft Store, that is a stub rather than an
interpreter; this script skips those. Install real Python from
https://www.python.org/downloads/windows/ (tick "Add python.exe to PATH"), or
use an existing one by pointing at it directly:

    .\build.ps1 -Python "C:\Users\you\anaconda3\python.exe"
"@
    }

    $best = $cands | Sort-Object -Property MissingCount | Select-Object -First 1
    $script:PyExe = $best.Exe; $script:PyPre = $best.Pre; $script:PyInfo = $best

    if ($cands.Count -gt 1) {
        Write-Host "    found $($cands.Count) Python installs; using $($best.Exe) (Python $($best.Version))" -ForegroundColor DarkGray
    }
}

# Resolved on first use, so targets that need no Python (help, clean, lean)
# work on a machine that has none.
function Initialize-Python {
    if (-not $script:PyExe) { Resolve-Python }
}

function Invoke-Py {
    param([Parameter(ValueFromRemainingArguments=$true)][string[]]$PyArgs)
    & $script:PyExe @($script:PyPre + $PyArgs)
    if ($LASTEXITCODE -ne 0) { Fail "python exited with code $LASTEXITCODE" }
}

function Test-Deps {
    [string[]]$miss = @($script:PyInfo.Missing)
    if ($script:PyInfo.MissingCount -gt 0) {
        $pip = ($script:PyPre + @('-m','pip','install') + $miss) -join ' '
        Fail @"
The chosen Python is missing: $($miss -join ', ')

It is:
    $($script:PyExe)   (Python $($script:PyInfo.Version))

Install the packages into that same interpreter:

    & "$($script:PyExe)" $pip

If you would rather use a different Python that already has them, point at it:

    .\build.ps1 -Python "C:\path\to\python.exe"
"@
    }
}

# ---------------------------------------------------------------- census
function Build-Census {
    if (Test-Path $Census) { Write-Ok "census already present: $Census"; return }
    New-Item -ItemType Directory -Force -Path (Split-Path $Census) | Out-Null

    # fgp_sieve.c uses __builtin_ctz and unsigned __int128, which are GCC and
    # Clang extensions. MSVC (cl.exe) cannot compile it, so only gcc or clang
    # are tried; otherwise the pure-Python generator runs instead. Both produce
    # a byte-identical census.
    $cc = $null
    if (-not $NoCompiler) {
        foreach ($c in @('gcc','clang')) {
            $f = Get-Command $c -ErrorAction SilentlyContinue
            if ($f) { $cc = $f.Source; break }
        }
    }

    if ($cc) {
        Write-Step "census via the C sieve ($([System.IO.Path]::GetFileName($cc)))"
        & $cc -O3 -funroll-loops -o $SieveEx $SieveC -lm
        if ($LASTEXITCODE -ne 0) { Fail "compiling fgp_sieve.c failed" }
        # PowerShell's '>' writes UTF-16 on Windows PowerShell 5.1, which would
        # corrupt the CSV. Start-Process redirects the raw bytes instead.
        $log = Join-Path $Root 'code\sieve.log'
        $p = Start-Process -FilePath $SieveEx -ArgumentList $Limit,$Pi `
                           -NoNewWindow -Wait -PassThru `
                           -RedirectStandardOutput $Census -RedirectStandardError $log
        if ($p.ExitCode -ne 0) {
            Remove-Item $Census -ErrorAction SilentlyContinue
            Fail "fgp_sieve exited with code $($p.ExitCode); see $log"
        }
        Get-Content $log -Tail 4 | ForEach-Object { Write-Host "    $_" }
    }
    else {
        Write-Step "no gcc or clang found; census via make_census.py (slower, identical output)"
        Invoke-Py (Join-Path $Root 'code\make_census.py') `
                  '--limit' $Limit '--pi' $Pi '--out' $Census
    }
    Write-Ok "census written: $Census"
}

# ---------------------------------------------------------------- stages
function Invoke-Check {
    Build-Census
    Write-Step 'the 45 offline accuracy checks'
    Invoke-Py (Join-Path $Root 'code\selftest.py') $Census
}

function Invoke-Analysis {
    Build-Census
    Write-Step "tables, figures and macros (verify limit $Verify)"
    Invoke-Py (Join-Path $Root 'code\run_all.py') $Census `
              '--out' (Join-Path $Root 'paper') '--verify-limit' "$Verify"
}

# Runs a LaTeX tool quietly and keeps its output for the failure path.
function Invoke-Tex([string]$tool, [string[]]$toolArgs) {
    $out = & $tool @toolArgs 2>&1
    return $out
}

function Invoke-Paper {
    if (-not (Test-Path $Macros)) { Invoke-Analysis }
    if (-not (Get-Command 'pdflatex' -ErrorAction SilentlyContinue)) {
        Fail @"
pdflatex was not found on PATH.

Install MiKTeX (https://miktex.org/download) or TeX Live, then open a new
terminal. Everything else has already been built, so re-running this script
will go straight to the manuscript.
"@
    }

    # MiKTeX ships latexmk as a Perl script, and Windows has no Perl, so
    # latexmk stops with "could not find the script engine 'perl'". Rather than
    # require another install, fall back to the explicit passes latexmk would
    # have made. The document needs bibtex and then two more runs for the
    # citations and cross-references to settle; the result is identical.
    $haveLatexmk = (Get-Command 'latexmk' -ErrorAction SilentlyContinue) -and
                   (Get-Command 'perl'    -ErrorAction SilentlyContinue)

    Push-Location (Join-Path $Root 'paper')
    try {
        Remove-Item 'main.pdf' -Force -ErrorAction SilentlyContinue
        if ($haveLatexmk) {
            Write-Step 'manuscript (latexmk)'
            $log = Invoke-Tex 'latexmk' @('-pdf','-interaction=nonstopmode','main.tex')
        }
        else {
            Write-Step 'manuscript (pdflatex and bibtex; latexmk needs Perl, which is absent)'
            $log = @()
            $log += Invoke-Tex 'pdflatex' @('-interaction=nonstopmode','main.tex')
            if (Get-Command 'bibtex' -ErrorAction SilentlyContinue) {
                $log += Invoke-Tex 'bibtex' @('main')
            }
            $log += Invoke-Tex 'pdflatex' @('-interaction=nonstopmode','main.tex')
            $log += Invoke-Tex 'pdflatex' @('-interaction=nonstopmode','main.tex')
        }

        if (-not (Test-Path 'main.pdf')) {
            Write-Host 'Last lines of the LaTeX output:' -ForegroundColor Yellow
            $log | Select-Object -Last 25 | ForEach-Object { Write-Host "    $_" -ForegroundColor Yellow }
            Fail 'the manuscript did not build; the full log is paper\main.log'
        }

        # Report anything still unresolved, which usually means a pass too few.
        if (Test-Path 'main.log') {
            $l = Get-Content 'main.log' -Raw -ErrorAction SilentlyContinue
            $undef = ([regex]::Matches($l, '(Reference|Citation) .* undefined')).Count
            if ($undef -gt 0) {
                Write-Host "    warning: $undef undefined reference(s) or citation(s); run this target again" -ForegroundColor Yellow
            }
        }
    } finally { Pop-Location }
    Write-Ok "built $Pdf"
}

function Invoke-Lean {
    if (-not (Get-Command 'lake' -ErrorAction SilentlyContinue)) {
        Fail @"
lake was not found on PATH. Install the Lean toolchain first:

    curl -O --location https://elan.lean-lang.org/elan-init.ps1
    powershell -ExecutionPolicy Bypass -File .\elan-init.ps1

then open a new terminal so PATH picks it up.
"@
    }
    Push-Location (Join-Path $Root 'lean')
    try {
        Write-Step 'downloading prebuilt Mathlib (several GB on first run)'
        & lake exe cache get
        if ($LASTEXITCODE -ne 0) { Fail 'lake exe cache get failed' }
        Write-Step 'checking FGP.lean'
        & lake build
        if ($LASTEXITCODE -ne 0) { Fail 'lake build failed' }
    } finally { Pop-Location }
    Write-Ok 'Lean proofs checked'
}

function Invoke-Clean([switch]$Dist) {
    Write-Step 'cleaning'
    foreach ($g in @('paper\csv\*.csv','paper\figures\*.pdf','paper\figures\*.png',
                     'paper\tables\*.tex','paper\main.pdf','paper\main.aux',
                     'paper\main.bbl','paper\main.blg','paper\main.log',
                     'paper\main.out','paper\main.fls','paper\main.fdb_latexmk',
                     'code\sieve.log','code\fgp_sieve.exe')) {
        Remove-Item (Join-Path $Root $g) -Force -ErrorAction SilentlyContinue
    }
    foreach ($d in @('code\fgp\__pycache__','lean\.lake')) {
        Remove-Item (Join-Path $Root $d) -Recurse -Force -ErrorAction SilentlyContinue
    }
    if ($Dist) { Remove-Item $Census -Force -ErrorAction SilentlyContinue }
    Write-Ok 'clean'
}

# ---------------------------------------------------------------- dispatch
try {
switch ($Target) {
    'help'      {
@"
Forward-Gap Partition census -- Windows build.

Nothing generated is stored in this repository. The census, the tables, the
figures and every number quoted in the manuscript are produced from source by
the targets below, so the paper cannot disagree with the computation.

  .\build.ps1                 census, checks, analysis and PDF (the usual case)
  .\build.ps1 census          just the layer census
  .\build.ps1 check           the 45 offline accuracy checks
  .\build.ps1 analysis        tables, figures and macros
  .\build.ps1 paper           the PDF
  .\build.ps1 full            as the default, but verifies all 20 checkpoints
  .\build.ps1 lean            check the Lean proofs (needs elan)
  .\build.ps1 clean           remove generated files, keep the census
  .\build.ps1 distclean       remove generated files and the census

Options:
  -Verify <n>     upper X for the independent reimplementation (0 disables)
  -NoCompiler     ignore gcc/clang and use the pure-Python census generator
  -Python <path>  use a specific interpreter, e.g. your Anaconda one:
                  .\build.ps1 -Python "\$HOME\anaconda3\python.exe"

If PowerShell refuses to run the script, start it with:
  powershell -ExecutionPolicy Bypass -File .\build.ps1
"@ | Write-Host
        break }
    'clean'     { Invoke-Clean; break }
    'distclean' { Invoke-Clean -Dist; break }
    'lean'      { Invoke-Lean; break }
    'census'    { Initialize-Python; Test-Deps; Build-Census; break }
    'check'     { Initialize-Python; Test-Deps; Invoke-Check; break }
    'analysis'  { Initialize-Python; Test-Deps; Invoke-Analysis; break }
    'paper'     { Initialize-Python; Test-Deps; Invoke-Check; Invoke-Paper; break }
    'full'      { Initialize-Python; Test-Deps; $Verify = 10000000000; Invoke-Check; Invoke-Analysis; Invoke-Paper; break }
    default     { Initialize-Python; Test-Deps; Invoke-Check; Invoke-Analysis; Invoke-Paper; break }
}
}
finally {
    if ($script:ProbeFile -and (Test-Path $script:ProbeFile)) {
        Remove-Item $script:ProbeFile -Force -ErrorAction SilentlyContinue
    }
}
