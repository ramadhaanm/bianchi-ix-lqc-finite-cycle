"""Build the manuscript with pdfLaTeX and BibTeX; no Perl dependency."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT/'paper/main.tex'
for command in (["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source.name],
                ["bibtex", source.stem],
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source.name],
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source.name]):
    result = subprocess.run(command, cwd=source.parent, capture_output=True, text=True,
                            errors="replace")
    if result.returncode:
        print(result.stdout)
        print(result.stderr)
        raise SystemExit(result.returncode)
    print("Completed:", " ".join(command), flush=True)
log = source.with_suffix('.log').read_text(errors='replace')
for line in log.splitlines():
    if any(word in line for word in ('Warning:', 'Overfull', 'Underfull')):
        print(line)
if 'There were undefined references' in log or 'There were undefined citations' in log:
    raise SystemExit('Unresolved references or citations after final pass.')
print('Built:', source.with_suffix('.pdf'))
