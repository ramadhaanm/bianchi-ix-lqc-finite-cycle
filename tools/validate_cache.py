"""Compare records from explicitly supplied reference/candidate result directories.

For a fresh cache check, generate N=2 K/A files with reanalyze.jl and run_cached.jl
in isolated, otherwise identical project copies. The verifier does not infer that
a production dataset was generated without the cache.
"""
from pathlib import Path
import argparse
import csv
import hashlib

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    with path.open() as f:
        meta = dict(kv.split('=', 1) for kv in f.readline()[2:].split())
        rows = {}
        for row in csv.DictReader(f):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'Incomplete row: {path}')
            key = (int(row['id']), row['role'])
            if key in rows:
                raise ValueError(f'Duplicate record {key}: {path}')
            rows[key] = row
    return meta, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir', type=Path, required=True)
    parser.add_argument('--candidate-dir', type=Path, required=True)
    parser.add_argument('--reference-n', type=int, default=2)
    parser.add_argument('--candidate-n', type=int, default=2)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    expected = {(k, role) for k in (1, 2) for role in ('base', 'ctrl', 'e1', 'e2')}
    lines = ['Reference: ' + str(args.reference_dir),
             'Candidate: ' + str(args.candidate_dir)]
    for quant in ('K', 'A'):
        def filename(n):
            return f'long_L20_{quant}_N{n}_on0.5_b300000.csv'
        ref_meta, ref = read(args.reference_dir / filename(args.reference_n))
        cand_meta, cand = read(args.candidate_dir / filename(args.candidate_n))
        if {k: v for k, v in ref_meta.items() if k != 'N'} != {
                k: v for k, v in cand_meta.items() if k != 'N'}:
            raise ValueError(f'Incompatible metadata or simulation sources: {quant}')
        if not (expected <= ref.keys() and expected <= cand.keys()):
            raise ValueError(f'Missing first-two-point controlled records: {quant}')
        for key in sorted(expected):
            if ref[key] != cand[key]:
                raise ValueError(f'Serialized records differ: {quant} {key}')
        lines.append(f'PASS {quant}: all fields match exactly for 8 complete controlled records.')
    lines += ['PASS: 16/16 records match; no numerical tolerance used.',
              'Execution methods must be established by the datasets generation logs.',
              'Local cache SHA256: ' + hashlib.sha256(
                  (ROOT / 'julia/VernerDoubleCache.jl').read_bytes()).hexdigest()]
    report = '\n'.join(lines) + '\n'
    if args.report:
        args.report.write_text(report, encoding='utf-8')
    print(report, end='')


if __name__ == '__main__':
    main()
