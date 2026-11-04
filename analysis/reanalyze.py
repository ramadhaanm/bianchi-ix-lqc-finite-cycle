"""Analyze the pathwise-filtered raw records; generate the manuscript tables and figures.

Run after tools/run_reanalysis.ps1 finishes for both K and A:
    .venv/Scripts/python analysis/reanalyze.py
Use --partial for inspection during an unfinished run; no manuscript artifacts are then written.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/revised"
F = 14


def interval(x, seed=20261026):
    """Median, 16/84 population percentiles, percentile-bootstrap 95% CI of median."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0, median=None, p16=None, p84=None, lo=None, hi=None)
    if len(x)==1:
        value=float(x[0])
        return dict(n=1,median=value,p16=value,p84=value,lo=None,hi=None)
    rng = np.random.default_rng(seed)
    boot = np.median(x[rng.integers(len(x), size=(2000, len(x)))], axis=1)
    q = np.percentile(x, [50, 16, 84])
    ci = np.percentile(boot, [2.5, 97.5])
    return dict(n=len(x), median=float(q[0]), p16=float(q[1]), p84=float(q[2]), lo=float(ci[0]), hi=float(ci[1]))


def wilson(k, n):
    if not n:
        return [None, None]
    z = norm.ppf(.975)
    p = k/n
    center = (p+z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [0. if k==0 else max(0., float(center-half)),
            1. if k==n else min(1., float(center+half))]


class Run:
    def __init__(self, path):
        self.path = Path(path)
        with self.path.open() as f:
            self.meta = dict(kv.split('=', 1) for kv in f.readline()[2:].split())
            expected = hashlib.sha256(b''.join((ROOT/'julia'/name).read_bytes() for name in
                ('BianchiLQC.jl','Cycles.jl','CornishLevin.jl','reanalyze.jl'))).hexdigest()
            if self.meta.get('schema') != '3' or self.meta.get('source') != expected:
                raise ValueError(f"schema or simulation source mismatch: {path}")
            self.rows = {}
            self.interior = set()
            for r in csv.DictReader(f):
                if None in r or any(v is None for v in r.values()):
                    raise ValueError(f"incomplete row: {path}")
                key = (int(r['id']), r['role'])
                if 1.01 <= float(r['u']) <= 1.99 and 1.01 <= float(r['v']) <= 1.99:
                    self.interior.add(key[0])
                if key in self.rows:
                    raise ValueError(f"duplicate job {key}: {path}")
                names = list(r)
                self.rows[key] = dict(
                    a=np.array([float(r[n]) for n in names[12:]]).reshape(-1, F),
                    reason=r['reason'], completed=int(r['completed']),
                    initial=float(r['initial_lnV']), min=float(r['final_min']),
                    error=float(r['final_error']))
        self.N = int(self.meta['N'])
        self.cycles = int(self.meta['cycles'])
        self.eps = [float(x) for x in self.meta['eps'].split(';')]
        self.ids = list(range(1, self.N*self.N+1 if self.meta['mode']=='basin' else self.N+1))
        self.roles = ['base'] if self.meta['mode']=='basin' else ['base', 'ctrl'] + [f'e{i+1}' for i in range(len(self.eps))]
        expected_keys = {(k,r) for k in self.ids for r in self.roles}
        if not self.rows.keys() <= expected_keys:
            raise ValueError(f"unexpected point or role: {path}")
        self.complete = len(self.rows) == len(self.ids)*len(self.roles)

    def event(self, k, role, event):
        r = self.rows.get((k, role))
        return None if r is None else r['a'][event]

    def valid(self, k, role, event, cutoff=1):
        a = self.event(k, role, event)
        bit = {1:1, 3:2, 10:4}[cutoff]
        return bool(a is not None and np.all(np.isfinite(a)) and not (int(a[13]) & bit))

    def difference(self, k, role, event):
        a, b = self.event(k,role,event), self.event(k,'base',event)
        return (a[:6:2]-b[:6:2])+(a[1:6:2]-b[1:6:2])

    def eligible(self, k, events, cutoff=1):
        return all(self.valid(k,r,e,cutoff) for r in ('base','ctrl','e1','e2') for e in events)

    def quality(self, k, events, noise_fraction=.01, linear_fraction=.1):
        for e in events:
            small = self.difference(k,'e2',e)
            big = self.difference(k,'e1',e)*(self.eps[1]/self.eps[0])
            noise = self.difference(k,'ctrl',e)
            size = np.linalg.norm(small)
            if not 0 < size <= .1:
                return False
            if np.linalg.norm(small-big) > linear_fraction*size or np.linalg.norm(noise) > noise_fraction*size:
                return False
        return True

    def gain(self, k, start, end):
        return float(np.log10(np.linalg.norm(self.difference(k,'e2',end)) /
                              np.linalg.norm(self.difference(k,'e2',start))))

    def status(self):
        """First stop reason; for phase runs P after the last R is also attempted."""
        names = ('complete','physical','constraint','budget','integration')
        return {s:sum(r['reason']==s for (k,role),r in self.rows.items() if role=='base') for s in names}


def phase_summary(run):
    # R1 -> P1 -> On2 -> Off2 -> R2: a complete, telescoping recollapse cycle.
    endpoints = [3,4,5,7,8]
    phases = [('turnaround',3,4),('contraction',4,5),('bounce',5,7),('expansion',7,8)]
    res = dict(meta=run.meta, jobs=len(run.rows), complete=run.complete, stop=run.status(), cutoffs={})
    for cutoff in (1,3,10):
        valid = [k for k in run.ids if run.eligible(k,endpoints,cutoff)]
        good = [k for k in valid if run.quality(k,endpoints)]
        vals = {name:[run.gain(k,s,e) for k in good] for name,s,e in phases}
        total = [run.gain(k,3,8) for k in good]
        if good:
            assert np.allclose(np.sum(list(vals.values()),axis=0),total,rtol=1e-12,atol=1e-12)
        stats = {name:interval(x) for name,x in vals.items()}
        stats['total'] = interval(total)
        stats['eligible'] = len(valid)
        stats['resolved'] = len(good)
        stats['ids'] = good
        res['cutoffs'][str(cutoff)] = stats
    # Event-specific survival, separate from any later failure or quality condition.
    res['survival'] = {str(c):sum(run.valid(k,'base',5*c-2) for k in run.ids) for c in range(1,run.cycles+1)}
    depth = [run.rows[k,'base']['initial']-run.event(k,'base',1)[6]
             for k in run.ids if run.valid(k,'base',1)]
    res['volume_depth'] = interval(depth)
    res['scale_depth'] = interval(np.array(depth)/3)
    # H_i here is d ln(a_i)/dT; convert to proper time using N=gamma sqrt(Delta).
    gamma = .2375
    lapse = gamma*np.sqrt(4*np.sqrt(3)*np.pi*gamma)
    good = res['cutoffs']['1']['ids']
    cohort_depth=[run.rows[k,'base']['initial']-run.event(k,'base',1)[6] for k in good]
    res['cohort_volume_depth']=interval(cohort_depth)
    res['cohort_scale_depth']=interval(np.array(cohort_depth)/3)
    res['Hmax_On2'] = interval([np.max(np.abs(run.event(k,'base',5)[8:11]))/lapse for k in good])
    eligible = [k for k in run.ids if run.eligible(k,endpoints)]
    res['quality_sensitivity'] = {}
    for name,noise,linear in [('strict',.001,.05),('primary',.01,.1),('relaxed',.1,.2)]:
        cohort = [k for k in eligible if run.quality(k,endpoints,noise,linear)]
        res['quality_sensitivity'][name] = dict(noise_fraction=noise,linear_fraction=linear,
            n=len(cohort),bounce=interval([run.gain(k,5,7) for k in cohort]),
            total=interval([run.gain(k,3,8) for k in cohort]))
    return res


def growth_summary(run):
    out = dict(meta=run.meta, jobs=len(run.rows), complete=run.complete, stop=run.status(),
               steps={},cutoffs={str(c):{} for c in (1,3,10)})
    for ev, offset in [('B',1),('R',3)]:
        perstep = []
        for c in range(run.cycles-1):
            s,e = 5*c+offset,5*(c+1)+offset
            valid = [k for k in run.ids if run.eligible(k,[s,e])]
            good = [k for k in valid if run.quality(k,[s,e])]
            st = interval([run.gain(k,s,e) for k in good])
            st.update(eligible=len(valid),ids=good)
            perstep.append(st)
        out['steps'][ev] = perstep
        all_events = [5*c+offset for c in range(run.cycles)]
        previous=set(run.ids)
        for cutoff in (1,3,10):
            eligible=[k for k in run.ids if run.eligible(k,all_events,cutoff)]
            cohort=[k for k in eligible if run.quality(k,all_events)]
            assert set(cohort) <= previous
            previous=set(cohort)
            st=interval([run.gain(k,all_events[0],all_events[-1])/(run.cycles-1) for k in cohort])
            st.update(ids=cohort,eligible=len(eligible))
            out['cutoffs'][str(cutoff)][ev]=st
        out['steps'][ev+'_cohort']=out['cutoffs']['1'][ev]
    out['survival'] = {str(c):sum(run.valid(k,'base',5*c-2) for k in run.ids) for c in range(1,run.cycles+1)}
    return out


def matched_threshold_summary(runs):
    """Compare markers on identical point IDs, resolved at every marker's endpoints."""
    endpoints = [3,4,5,7,8]
    ids = sorted(set.intersection(*(set(r.ids[:200]) for r in runs.values())))
    cohort = [k for k in ids if all(r.eligible(k,endpoints) and r.quality(k,endpoints)
                                  for r in runs.values())]
    return {str(h):dict(ids=cohort, resolved=len(cohort),
        bounce=interval([r.gain(k,5,7) for k in cohort]),
        total=interval([r.gain(k,3,8) for k in cohort])) for h,r in runs.items()}


def label(run,k,role,event,kind):
    a = run.event(k,role,event)
    # At contracting onset the most negative directional rate is the fastest contracting axis.
    return int(np.argmin(a[8:11] if kind=='On' else a[:6:2]+a[1:6:2]))


def label_summary(run):
    res = dict(meta=run.meta, jobs=len(run.rows), complete=run.complete, stop=run.status(), events={})
    for c in range(run.cycles):
        for kind,offset in [('On',0),('B',1),('R',3)]:
            e = 5*c+offset
            # Keep base and perturbation inside the declared initial square.
            # Points are sampled inside [1.01,1.99]^2 for boundary-sensitive statistics by
            # restricting the base coordinates; epsilon <= .01 then remains inside [1,2]^2.
            # IDs are independent of radius, so the same base restriction is used for all eps.
            eligible_ids = []
            for k in run.ids:
                # Regenerate no random numbers: u,v are read below from a separate CSV pass.
                if k in run.interior:
                    eligible_ids.append(k)
            series = []
            for j,eps in enumerate([0.]+run.eps):
                role = 'ctrl' if j==0 else f'e{j}'
                valid = [k for k in eligible_ids if run.valid(k,'base',e) and run.valid(k,role,e)]
                flips = sum(label(run,k,'base',e,kind)!=label(run,k,role,e,kind) for k in valid)
                series.append(dict(eps=eps,n=len(valid),flips=flips,f=flips/len(valid) if valid else None,ci=wilson(flips,len(valid))))
            brackets = []
            arr = series[1:]
            for a,b in zip(arr,arr[1:]):
                if a['f'] is not None and b['f'] is not None and a['f']>=.1>b['f']:
                    brackets.append([b['eps'],a['eps']])
            res['events'][f'{kind}{c+1}'] = dict(series=series,brackets=brackets)
    return res


def read_run(mode,L,q,N,onset=.5,budget=300000):
    p = OUT/f'{mode}_L{L}_{q}_N{N}_on{onset}_b{budget}.csv'
    return Run(p) if p.exists() else None


def fmt(s, ci=True):
    if s['median'] is None:
        return '--'
    if ci and s['lo'] is None:
        return f"{s['median']:.2f} (CI n/a)"
    return (f"{s['median']:.2f} [{s['lo']:.2f}, {s['hi']:.2f}]" if ci else f"{s['median']:.2f}")


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--partial',action='store_true'); args=parser.parse_args()
    summary = dict(phases={},growth={},labels={},thresholds={},matched_thresholds={})
    missing=[]
    runs={}
    for q in ('K','A'):
        for L in (10,20,30):
            run=read_run('phase',L,q,1000)
            if run is None or not run.complete:
                missing.append(f'phase {L} {q}')
            if run is not None and run.rows:
                summary['phases'][f'{q}{L}']=phase_summary(run); runs[f'{q}{L}']=run
        run=read_run('long',20,q,64)
        if run is None or not run.complete: missing.append(f'long {q}')
        if run is not None and run.rows: summary['growth'][q]=growth_summary(run)
        run=read_run('labels',20,q,400)
        if run is None or not run.complete: missing.append(f'labels {q}')
        if run is not None and run.rows:
            summary['labels'][q]=label_summary(run)
        threshold_runs = {.5:runs[f'{q}20']} if f'{q}20' in runs else {}
        for onset in (.3,.7):
            run=read_run('phase',20,q,200,onset)
            if run is None or not run.complete: missing.append(f'threshold {q} {onset}')
            if run is not None and run.rows:
                summary['thresholds'][f'{q}{onset}']=phase_summary(run)
                threshold_runs[onset]=run
        if len(threshold_runs)==3 and all(r.complete for r in threshold_runs.values()):
            summary['matched_thresholds'][q]=matched_threshold_summary(threshold_runs)
    basin=read_run('basin',20,'K',100)
    if basin is None or not basin.complete: missing.append('basin')
    summary['missing']=missing
    outname='partial_summary.json' if args.partial else 'summary.json'
    if missing and not args.partial:
        raise SystemExit('Incomplete production runs: '+', '.join(missing))
    (OUT/outname).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    for name,p in summary['phases'].items():
        z=p['cutoffs']['1']
        print(name,'jobs',p['jobs'],'R survivors',p['survival'],'phase cohort',z['resolved'], '/', z['eligible'],
              'bounce',fmt(z['bounce']),'total',fmt(z['total']))
    for q,g in summary['growth'].items():
        print('long',q,g['jobs'],g['survival'],'B rate',fmt(g['steps']['B_cohort']),'R rate',fmt(g['steps']['R_cohort']))
    print('Missing:',missing)
    if not args.partial:
        make_artifacts(summary,runs,basin)


def make_artifacts(summary,runs,basin):
    # Imported here so --partial remains lightweight and cannot overwrite manuscript figures.
    from figures import make_figures, make_tables
    make_figures(summary,runs,basin)
    make_tables(summary)


if __name__=='__main__':
    main()
