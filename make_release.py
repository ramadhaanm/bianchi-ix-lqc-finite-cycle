"""Write MANIFEST_SHA256.json and build a versioned archive of the code, data and scripts."""
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
        'julia':['BianchiLQC.jl','Cycles.jl','CornishLevin.jl','reanalyze.jl','run_cached.jl','VernerDoubleCache.jl','test_validity.jl','test_baseline.jl','test_A.jl','Project.toml','Manifest.toml','uncertainty_classical.jl','basins_classical_d64.jl'],
        'derivations':['*.py','DERIVATIONS.md'],
        'analysis':['reanalyze.py','figures.py','test_reanalysis.py','check_data.py'],
        'results/revised':['*.csv','summary.json','RESULTS.md','*.log'],
        'results':['cl_basins_fig2_200_d64.csv','uncertainty_classical_d64.txt'],
        'tools':['run_reanalysis.ps1','finish_reanalysis.ps1','run_all.sh','validate_cache.py','validation/*.log','validation/VALIDATION.md'],
        '.':['REPRODUCE.md','requirements.txt','make_release.py','LICENSE','CITATION.cff','README.md'],
    }
    for sub,pats in patterns.items():
        for pat in pats:
            files.extend((ROOT/sub).glob(pat))
    files=sorted(set(p for p in files if p.is_file() and not
        (p.suffix=='.csv' and p.name.startswith('long_') and '_N64_' not in p.name)))
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (ROOT/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2)+chr(10))
    fingerprint=hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()[:12]
    dest=ROOT/f'code_release_{fingerprint}.zip'
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
