"""Create a versioned local submission archive without deleting an earlier release."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parent


def main():
    summary=json.loads((ROOT/'results/revised/summary.json').read_text())
    if summary['missing']:
        raise SystemExit('Refusing to package incomplete analysis.')
    files=[]
    patterns={
        'paper':['main.tex','main.pdf','references.bib','iopart-num-doi.bst','generated/*.tex','figs/fig[1-5]_*.pdf','figs/fig[1-5]_*.png'],
        'julia':['BianchiLQC.jl','Cycles.jl','CornishLevin.jl','reanalyze.jl','run_cached.jl','VernerDoubleCache.jl','test_validity.jl','test_baseline.jl','test_A.jl','Project.toml','Manifest.toml','uncertainty_classical.jl','basins_classical_d64.jl'],
        'derivations':['*.py','DERIVATIONS.md'],
        'analysis':['reanalyze.py','revised_figures.py','test_reanalysis.py','check_revised_data.py'],
        'results/revised':['*.csv','summary.json','RESULTS.md','*.log'],
        'results':['cl_basins_fig2_200_d64.csv','uncertainty_classical_d64.txt'],
        'review':['run_reanalysis.ps1','finish_reanalysis.ps1','rerun_fixed.sh','build_pdf.py','validate_cache.py','REVISION_NOTES.md','validation/*.log','validation/*.json','validation/VALIDATION.md'],
        '.':['REPRODUCE_REVISED.md','requirements.txt','make_revised_release.py','LICENSE','CITATION.cff','README.md'],
    }
    for sub,pats in patterns.items():
        for pat in pats:
            files.extend((ROOT/sub).glob(pat))
    files=sorted(set(p for p in files if p.is_file() and not
        (p.suffix=='.csv' and p.name.startswith('long_') and '_N64_' not in p.name)))
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    fingerprint=hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()[:12]
    dest=ROOT/f'code_release_revised_{fingerprint}.zip'
    if dest.exists():
        print('Identical version already packaged:',dest)
        return
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files:
            z.write(p,str(p.relative_to(ROOT)).replace('\\','/'))
        z.writestr('MANIFEST_SHA256.json',json.dumps(manifest,indent=2)+'\n')
    print(dest)
    print(f'{len(files)} files; {dest.stat().st_size/1e6:.1f} MB; manifest {fingerprint}')


if __name__=='__main__':
    main()
