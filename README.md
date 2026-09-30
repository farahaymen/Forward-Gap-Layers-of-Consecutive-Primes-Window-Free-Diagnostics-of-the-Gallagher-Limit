# Forward-Gap Partition census of consecutive prime gaps

Source only. **Nothing generated is stored here**: not the census, not the
tables, not the figures, and not the numbers quoted in the running text of the
manuscript. All of it is produced from this source by one command. That is
deliberate, and it is the reason the paper cannot disagree with the
computation: every table is an included fragment and every number in the prose
is a `\newcommand` written by the same run that produced the figures.

## Build

### Windows (PowerShell, including the VS Code terminal)

From the repository root:

```powershell
.\build.ps1
```

That compiles the sieve if a compiler is present, builds the census, runs the
45 accuracy checks, writes the tables, figures and macros, and typesets the
manuscript to `paper\main.pdf`.

If PowerShell refuses with "running scripts is disabled on this system", that
is Windows' execution policy, not a problem with the script. Either run it this
way each time, which changes no setting:

```powershell
powershell -ExecutionPolicy Bypass -File .\build.ps1
```

or allow local scripts once, for your account only:

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

**If the terminal advertises the Microsoft Store when you run it**, Windows is
shadowing Python with a stub that only opens the Store. Install real Python
from python.org, ticking "Add python.exe to PATH", and switch off both stubs
under Settings > Apps > Advanced app settings > App execution aliases. The
script skips those stubs and falls back to the `py -3` launcher on its own, so
it will work as soon as any real Python 3 is present.

Other targets:

```powershell
.\build.ps1 census       # just the layer census
.\build.ps1 check        # the 45 offline accuracy checks
.\build.ps1 analysis     # tables, figures and macros
.\build.ps1 paper        # the PDF
.\build.ps1 full         # as the default, but verifies all 20 checkpoints
.\build.ps1 lean         # check the Lean proofs
.\build.ps1 clean        # remove generated files, keep the census
.\build.ps1 distclean    # remove generated files and the census
.\build.ps1 help         # this list
```

Options: `-Verify <n>` sets the upper X for the independent reimplementation
(`0` disables it), `-NoCompiler` ignores any C compiler and uses the pure-Python
census generator, and `-Python <path>` forces a particular interpreter.

**If you have several Pythons installed** (an Anaconda one, an MSYS2 one, the
Microsoft Store stubs), the script probes each and picks the one that actually
has `numpy`, `scipy`, `pandas`, `matplotlib` and `sympy`, rather than the first
one on PATH. It prints which it chose. To override:

```powershell
.\build.ps1 -Python "$HOME\anaconda3\python.exe"
```

**On the C compiler.** `fgp_sieve.c` uses `__builtin_ctz` and
`unsigned __int128`, which are GCC and Clang extensions. MSVC (`cl.exe`) cannot
compile it. The script therefore looks for `gcc` or `clang` only, and falls
back to `make_census.py` if neither is present. The fallback is slower (about
four minutes to 10^10 rather than seventeen seconds) and produces a
byte-identical census, so nothing is lost by not having a compiler.

### macOS and Linux

```bash
make
```

Other targets:

```bash
make census      # just the layer census    (~17 s, or ~4 min with no compiler)
make check       # the 45 offline accuracy checks
make analysis    # tables, figures and macros
make paper       # the PDF
make full        # as `make`, but verifies all 20 checkpoints (~2 min, ~1 GB)
make lean        # check the Lean proofs
make clean       # remove generated files, keep the census
make distclean   # remove generated files and the census
make help        # this list
```

### By hand, on any platform

```bash
cc -O3 -funroll-loops -o code/fgp_sieve code/fgp_sieve.c -lm
./code/fgp_sieve 10000001000 455052511 > data/fgp_layers.csv
python3 code/selftest.py data/fgp_layers.csv
python3 code/run_all.py  data/fgp_layers.csv --out paper --verify-limit 2100000000
cd paper && latexmk -pdf main.tex
```

Without a C compiler, replace the first two lines with:

```bash
python3 code/make_census.py --limit 10000001000 --pi 455052511 --out data/fgp_layers.csv
```

Do not use PowerShell's `>` to redirect the sieve's output on Windows: Windows
PowerShell 5.1 writes UTF-16 there and would corrupt the census. `build.ps1`
redirects the raw bytes instead, which is one reason to prefer it over typing
the commands.

Building the manuscript before the pipeline has run stops with an error naming
the command to run, rather than a page of undefined control sequences.

## What to expect

A correct run prints `45 passed, 0 failed` from the checks, and from the
analysis a census sha256 of
`bb82712165d660a002018bedebe5bfb6e94f4839b24ef532299fe0c5c26b24c8`,
`C_fit primary: 1.2483 -> 1.1460`, `E2: 0.84606 -> 0.88815`, and
`g*: 14.13 -> 22.37`. The manuscript comes out at 20 pages with no undefined
references. With `make full` the verification line reads
`2166 layer counts compared, 0 mismatches`.

## Layout

```
Makefile             the whole build (macOS, Linux)
build.ps1            the whole build (Windows PowerShell)
code/
  fgp_sieve.c        C segmented sieve; produces the layer census
  make_census.py     pure-Python census generator, byte-identical output,
                     for machines with no C compiler
  fgp/
    core.py          census loading, H(g), effective scale
    frequency.py     Poisson GLM fit of C_fit, the 36-specification grid,
                     and the conversion to Wolf's normalisation
    local.py         differenced-interval diagnostics E2, E3, D, W2 + bootstrap
    spatial.py       R(g;N), crossover locators, the parameter-free density
                     model, and the amplitude ratio
    verify.py        independent NumPy census, cross-checks the C output exactly
    figures.py       Figures 1-6
    tables.py        LaTeX tables and the \newcommand macros used in the text
  run_all.py         single entry point
  selftest.py        offline accuracy checks (45 of them)
data/                the census lands here (ships empty)
paper/
  main.tex           the manuscript
  assumptions.tex    Appendix C, the register of assumptions
  refs.bib
  tables/            generated LaTeX fragments and macros.tex (ships empty)
  figures/           generated PDF and PNG figures (ships empty)
  csv/               generated intermediate results (ships empty)
lean/
  FGP.lean           the four formal statements (Appendix B)
  lakefile.toml      Mathlib pinned to the revision they were checked against
  lean-toolchain     Lean v4.35.0-rc3
  lake-manifest.json the nine transitive dependencies, pinned
```

## Requirements

- Python 3.10+ with `numpy`, `scipy`, `pandas`, `matplotlib`, `sympy`
  (`sympy` is used only by the self-test)
- a C compiler (gcc or clang) is **optional**: `make_census.py` produces a
  byte-identical census in pure Python, just slower (about four minutes to
  10^10 rather than seventeen seconds)
- `make` is optional too; the plain commands above do the same thing
- for the manuscript: a TeX distribution providing `pdflatex`, `bibtex` and the
  packages `natbib`, `booktabs`, `microtype`, `subcaption`, `listings`,
  `mathptmx`, `placeins`. `latexmk` is used when available but is not required:
  it is a Perl script, and MiKTeX on Windows has no Perl, so the build falls
  back to running `pdflatex`, `bibtex`, `pdflatex`, `pdflatex` itself. The
  resulting PDF is identical.

```
pip install numpy scipy pandas matplotlib sympy
```

On Windows, if `pip` is not found, use `python -m pip install ...` or
`py -3 -m pip install ...`.

## Verifying accuracy offline

```bash
make check          # or, on Windows:  .\build.ps1 check
```

Forty-five checks, no network needed. These do not re-run the pipeline and compare it
with itself; each one confronts the computation with a fact established
elsewhere:

- **Published constants** hard-coded with their sources: pi(10^k) (OEIS A006880),
  p_n for n a power of ten (A006988), twin-prime counts pi_2(10^k) (A007508),
  and the maximal prime gap record (A005250 / A002386). The g=2 layer at 10^10
  must equal 27,412,679 twin pairs exactly.
- **Internal identities**: layer counts partition p_2..p_N; counts never
  decrease as N grows; no 64-bit overflow in the prime sums; the gaps telescope
  to p_{N+1} - p_2, which recovers 10,000,000,019 from the census.
- **A second sieve** written as plainly as possible reproduces the layer counts
  and all the constants above from scratch.
- **Method checks**: H(g) against symbolic factorisation, the Poisson IRLS fit
  against Nelder-Mead on the same likelihood, the bootstrap standard error
  against the delta method, and stability across bootstrap seeds.
- **Mathematics**: the exponential moment identities by quadrature, and the
  numerical density model against the analytic first-order expansion of
  Appendix A (they agree to 3e-4 at log X = 10^4).

The suite exits non-zero if anything fails. It found real errors during
development, including a wrong arithmetic factor H(g); see Section 8.3 and
Appendix C of the paper, which lists every assumption with its status.

Two further defects were found and fixed the same way:

- The census accumulator silently discarded gaps wider than its histogram
  (GMAX = 512). Harmless below 10^10, where the maximal gap is 354, but it
  would corrupt a run past roughly 7e11 with no visible sign. Both the C sieve
  and the Python generator now abort with a diagnostic instead.
- The independent verifier sieved exactly [2, limit], so it could not classify
  its own last prime and crashed on the largest checkpoint. A prime's layer is
  fixed by the NEXT prime, which at X = 10^10 is 10,000,000,019. It now sieves
  past the limit by a margin, exactly as the census does. This is why the full
  range can now be verified rather than only the first sixteen scales.

The default run already verifies the first sixteen scales. To check every
layer count at all twenty checkpoints against the independent implementation:

```bash
make full           # or, on Windows:  .\build.ps1 full
```

That takes about two minutes and roughly 1 GB. It is the strongest check
available here, and the reported results were produced with it: 2,166 layer
counts at 20 checkpoints, zero mismatches, to X = 9,999,999,967.

## What the pipeline writes

`paper/tables/*.tex` are LaTeX `tabular` fragments included by `main.tex`.
`paper/tables/macros.tex` defines every number quoted in the running text as a
`\newcommand`, so the prose cannot drift from the computation: if the census
changes, the text changes with it.
`paper/csv/*.csv` are the same results in plain form, including the full
36-specification grid, the amplitude ratio and the layer-by-layer spatial
table.
`paper/figures/*.pdf` are vector figures drawn at the manuscript's text width,
so `\includegraphics[width=\linewidth]` scales them 1:1 and the font sizes on
the page are the ones the code sets. `.png` previews sit alongside.

None of this is under version control here; `make clean` removes all of it and
`make` puts it back.

## Verification built into the run

`run_all.py` prints, and `tab_verify.tex` records, a comparison of the C census
against an independent NumPy reimplementation that uses a different data layout
(one byte per integer rather than one bit per odd integer) and a different way
of forming gaps (array differencing rather than a running pointer). The two programs share no code, no language and no data layout. Every layer
count at every checkpoint below the verify limit must match exactly; with
`make full` that covers all twenty checkpoints, 2,166 layer counts, to
X = 9,999,999,967. The run
also checks the local Prime Number Theorem on each differenced interval, the
maximal gap below 10^10 (354, following 4,302,407,359) and pi(x) at tabulated
points.

## Lean (Appendix B)

Separate from the pipeline above. It needs the Lean toolchain rather than
Python, and the first run downloads several GB of prebuilt Mathlib.

Install the toolchain once.

Windows (PowerShell):

```powershell
curl.exe -O --location https://elan.lean-lang.org/elan-init.ps1
powershell -ExecutionPolicy Bypass -File .\elan-init.ps1
```

macOS and Linux:

```bash
curl https://elan.lean-lang.org/elan-init.sh -sSf | sh -s -- -y
export PATH="$HOME/.elan/bin:$PATH"          # add this to your shell profile
```

Then open a new terminal so `lake` is on PATH.

Then, from the repository root:

```bash
make lean           # or, on Windows:  .\build.ps1 lean
```

or by hand:

```bash
cd lean
lake exe cache get      # downloads prebuilt Mathlib, several GB, once
lake build              # checks FGP.lean
```

Those two commands are identical on Windows.

A clean run prints nothing and exits 0. Anything else is a real failure worth
reading.

`lean/lakefile.toml` pins Mathlib to revision `18dc857edd2da9c0d63672b94fbaf48a9f4810b5`
and `lean/lean-toolchain` pins Lean to `v4.35.0-rc3`; `lake-manifest.json` pins
the nine transitive dependencies. That pairing is what the file was checked
against. To move to current Mathlib:

```bash
cd lean
sed -i 's/rev = ".*"/rev = "master"/' lakefile.toml
lake update
cp .lake/packages/mathlib/lean-toolchain ./lean-toolchain   # cache needs an exact match
lake exe cache get && lake build
```

### What has and has not been checked

Every Mathlib name the file uses was verified to exist with the right
signature at the pinned revision: `Nat.nth`, `Nat.Prime`, `Finset.Icc`,
`Finset.disjoint_filter`, `Nat.exists_prime_lt_and_le_two_mul`,
`Finset.sum_range_sub`, `integral_rpow_mul_exp_neg_mul_Ioi`,
`Real.Gamma_nat_eq_factorial`, `Real.rpow_natCast` and
`MeasureTheory.integral_const_mul`. One name was wrong and is now fixed: the
constant pulled out of an integral is `integral_const_mul`, not
`integral_mul_left`.

The proof **terms** have not been machine-checked end to end, because the
Mathlib binary cache was unreachable from the machine where this was prepared
and building Mathlib from source there was not practical. Run `make lean` on
your own machine before citing Appendix B as verified. The four results are
elementary, so the remaining risk is in the tactic scripts, not the statements.

## Licence

Code released into the public domain (CC0). The census it generates is a
property of the integers and carries no licence.
