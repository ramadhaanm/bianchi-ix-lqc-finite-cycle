"""Audit structural and physical-domain invariants of the declared manuscript records."""
from pathlib import Path
import argparse
import numpy as np
from reanalyze import OUT, Run


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--core',action='store_true',help='Audit completed ordinary-precision ensembles only.')
    args=parser.parse_args()
    count = events = 0
    for path in sorted(OUT.glob('*.csv')):
        if path.name.startswith('long_') and '_N64_' not in path.name:
            continue  # exploratory/performance pilots, not production ensembles
        if args.core and path.name.startswith('long'):
            continue
        run = Run(path)
        assert run.complete, f'incomplete: {path.name}'
        for (point,role),row in run.rows.items():
            a = row['a']
            slots = list(range(len(a))) if run.meta['mode']=='phase' else [i for i in range(len(a)) if i%5!=4]
            reached = [np.all(np.isfinite(a[i])) for i in slots]
            if not all(reached):
                first = reached.index(False)
                assert not any(reached[first:]), (path.name,point,role,'event after exit')
            prev_mask = 0
            prev_min = np.inf
            prev_error = 0.
            for i,ok in zip(slots,reached):
                if not ok:
                    assert np.all(np.isnan(a[i])), (path.name,point,role,'partial event')
                    continue
                mask = int(a[i,13])
                assert mask in (0,4,6), (path.name,point,role,'invalid cutoff mask')
                assert mask | prev_mask == mask, (path.name,point,role,'cutoff re-entry')
                assert a[i,11] <= prev_min + 1e-12, (path.name,point,role,'minimum increased')
                assert a[i,12] >= prev_error, (path.name,point,role,'residual maximum decreased')
                assert 0 <= a[i,12] < 1e-3, (path.name,point,role,'invalid saved constraint')
                assert a[i,11] >= -1e-12, (path.name,point,role,'invalid saved axis')
                if i%5 in (0,2):
                    assert abs(a[i,7]-float(run.meta['onset'])) < 1e-7, (path.name,point,role,'marker event')
                if i%5 in (1,3):
                    assert abs(sum(a[i,8:11])) < 1e-8, (path.name,point,role,'volume extremum')
                if i%5==3:
                    assert a[i,7] < .01, (path.name,point,role,'recollapse marker')
                if i%5==4:
                    assert abs(a[i,6]-a[i-1,6]+5) < 1e-7, (path.name,point,role,'turnaround depth')
                prev_mask,prev_min,prev_error = mask,a[i,11],a[i,12]
                events += 1
            if row['reason']=='complete':
                assert all(reached), (path.name,point,role,'incomplete complete row')
            assert row['completed']==sum(run.valid(point,role,5*c-2) for c in range(1,run.cycles+1))
            count += 1
        print(f'PASS {path.name}: {len(run.rows)} trajectories')
    print(f'PASS: {count} trajectory records; {events} saved valid events')


if __name__=='__main__':
    main()
