#!/usr/bin/env python3
"""Build the files to upload to Experimental Mathematics (ScholarOne).

The journal asks for (i) a PDF of the manuscript, uploaded as
"Manuscript - with author details", and (ii) the LaTeX source as a single zip,
uploaded as "LaTeX Source Files", containing everything needed to run the
paper: the .tex, the macro file, the .bib and .bbl, the .bst, and the figure
sources in EPS. This script produces both from paper/, in a flat layout:

    submission/
        latex_source/            the contents of LaTeX_Source_Files.zip
            main.tex             tables and the assumptions appendix inlined
            fgpmacros.tex        all macros, input once after the packages
            refs.bib  main.bbl  tfs.bst
            fig1_cfit.eps ...    EPS figure sources (plus PDF copies for pdflatex)
        LaTeX_Source_Files.zip
        Manuscript_with_author_details.pdf

Run it from the repository root after the paper has been built once (so that
paper/main.bbl and the figures exist):

    python3 code/make_submission.py

The manuscript PDF is compiled with pdflatex from exactly the source that goes
into the zip, so the two cannot disagree; the source is also compiled with
latex to check the EPS figures. Either check stops the script if it fails or
leaves an undefined reference. If pdflatex is not installed, the existing
paper/main.pdf is used only when it is newer than all of its inputs.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


def find_tex(engine: str):
    """Full path to a TeX engine, or None.

    Looks on PATH first, then in the standard MiKTeX and TeX Live install
    folders, because on Windows MiKTeX is often installed without being added
    to the PATH of the terminal that runs this script."""
    hit = shutil.which(engine)
    if hit:
        return hit
    if os.name == "nt":
        roots = [os.environ.get("LOCALAPPDATA", ""), os.environ.get("ProgramFiles", ""),
                 os.environ.get("ProgramFiles(x86)", "")]
        candidates = [Path(r) / "Programs" / "MiKTeX" / "miktex" / "bin" / "x64" for r in roots[:1] if r]
        candidates += [Path(r) / "MiKTeX" / "miktex" / "bin" / "x64" for r in roots[1:] if r]
        candidates += sorted(Path("C:/texlive").glob("*/bin/windows")) + sorted(Path("C:/texlive").glob("*/bin/win32"))
        for d in candidates:
            exe = d / f"{engine}.exe"
            if exe.exists():
                return str(exe)
    return None


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def flatten_main(paper: Path) -> str:
    tex = read(paper / "main.tex")

    def inline(m):
        name = m.group(1)
        if name == "fgpmacros":            # stays a separate, input-once file
            return m.group(0)
        path = paper / (name if name.endswith(".tex") else name + ".tex")
        if not path.exists():
            sys.exit(f"cannot inline \\input{{{name}}}: {path} not found")
        return read(path).strip()

    tex = re.sub(r"\\input\{([^}]+)\}", inline, tex)
    tex = tex.replace(r"\graphicspath{{figures/}}", "")
    return tex


def flatten_macros(paper: Path) -> str:
    mac = read(paper / "fgpmacros.tex")
    start = mac.index(r"\IfFileExists{tables/macros.tex}")
    end = mac.index(r"\input{tables/macros}") + len(r"\input{tables/macros}")
    numbers = read(paper / "tables" / "macros.tex").strip()
    return (mac[:start]
            + "% Numbers quoted in the text, generated from the census by code/run_all.py.\n"
            + numbers + "\n" + mac[end:])


def compile_check(src: Path, engine: str, exe: str):
    work = src.parent / f"_check_{engine}"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(src, work)
    if engine == "latex":                 # make sure the EPS sources are what gets used
        for f in work.glob("fig*.pdf"):
            f.unlink()
    r = None
    for _ in range(3):
        r = subprocess.run([exe, "-interaction=nonstopmode", "main.tex"],
                           cwd=work, capture_output=True, text=True)
    log = (work / "main.log").read_text(encoding="latin-1")
    problems = [l for l in log.splitlines() if l.startswith("!") or "undefined" in l.lower()]
    if r.returncode != 0 or problems:
        print("\n".join(problems[:20]))
        sys.exit(f"the source did not compile cleanly with {engine}; see {work}/main.log")
    out = work / ("main.pdf" if engine == "pdflatex" else "main.dvi")
    print(f"{engine}: compiled cleanly ({out.name})")
    return work


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper", default="paper", help="manuscript folder (default: paper)")
    ap.add_argument("--out", default="submission", help="output folder (default: submission)")
    ap.add_argument("--no-compile", action="store_true", help="skip the test compilations")
    a = ap.parse_args()

    paper, out = Path(a.paper), Path(a.out)
    if not (paper / "main.bbl").exists():
        sys.exit("paper/main.bbl not found: build the paper once first (make paper)")
    if out.exists():
        shutil.rmtree(out)
    src = out / "latex_source"
    src.mkdir(parents=True)

    (src / "main.tex").write_text(flatten_main(paper), encoding="utf-8")
    (src / "fgpmacros.tex").write_text(flatten_macros(paper), encoding="utf-8")
    for name in ("refs.bib", "main.bbl", "tfs.bst"):
        shutil.copy(paper / name, src / name)
    eps = sorted((paper / "figures").glob("fig*.eps"))
    pdf = sorted((paper / "figures").glob("fig*.pdf"))
    if not eps or len(eps) != len(pdf):
        sys.exit("EPS and PDF figures missing or unmatched: run the analysis first (make analysis)")
    for f in eps + pdf:
        shutil.copy(f, src / f.name)

    pdflatex, latex = find_tex("pdflatex"), find_tex("latex")
    manuscript = out / "Manuscript_with_author_details.pdf"

    if a.no_compile or pdflatex is None:
        # Without pdflatex the PDF cannot be built from the zipped source, so
        # the only option is the PDF already in paper/. Use it only if it is
        # newer than every file it is built from; otherwise the PDF and the
        # LaTeX source could disagree, which is exactly what must not happen.
        if pdflatex is None and not a.no_compile:
            print("WARNING: pdflatex was not found (looked on PATH and in the usual MiKTeX/TeX Live folders).")
        pdf = paper / "main.pdf"
        inputs = [paper / "main.tex", paper / "fgpmacros.tex", paper / "assumptions.tex",
                  paper / "main.bbl"] + list((paper / "tables").glob("*.tex")) \
                 + list((paper / "figures").glob("fig*.pdf"))
        newest = max(inputs, key=lambda f: f.stat().st_mtime)
        if not pdf.exists() or pdf.stat().st_mtime < newest.stat().st_mtime:
            sys.exit(f"paper/main.pdf is older than {newest.as_posix()}, so it may not match the LaTeX source.\n"
                     "Install MiKTeX (https://miktex.org/download) or TeX Live, open a new terminal,\n"
                     "build the paper (.\\build.ps1 paper, or make paper), then run this script again.")
        shutil.copy(pdf, manuscript)
        print("WARNING: copied paper/main.pdf without recompiling; it is newer than all of its inputs.")
    else:
        w = compile_check(src, "pdflatex", pdflatex)
        shutil.copy(w / "main.pdf", manuscript)
        shutil.rmtree(w)
        if latex is None:
            print("WARNING: latex was not found, so the EPS figures were not test-compiled.")
        else:
            shutil.rmtree(compile_check(src, "latex", latex))

    zpath = out / "LaTeX_Source_Files.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(src.iterdir()):
            z.write(f, f.name)
    print(f"wrote {out}/Manuscript_with_author_details.pdf and {zpath}")


if __name__ == "__main__":
    main()
