#!/usr/bin/env python3
"""Final-build replay, portable-build replay, and bracket guard qualification."""
import csv,hashlib,json
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent;ROOT=Path('/private/tmp/ems-t25')
spec=importlib.util.spec_from_file_location('work',HERE/'t25-work.py');w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
def main():
    rows=[];exe=ROOT/'build/t25-expansion-map.ex'
    portable=next((ROOT/'portable-build').glob('t25-expansion-map*.ex'))
    for case,kind,program in [('final-single','single',exe),('final-late','late',exe),('portable-single','single',portable)]:
        d=ROOT/'verification'/case;d.mkdir(parents=True,exist_ok=True)
        (d/'params.txt').write_text((ROOT/'maps'/f'{kind}-N192/params.txt').read_text())
        rc=w.measure('regression',d,[program,'params.txt'],3000000000)
        assert rc==0 and 'T25_EXPANSION_MAP_COMPLETE' in (d/'run.log').read_text()
        for f in ['map.csv','angles.csv']:
            expected=(ROOT/'maps'/f'{kind}-N192'/f).read_bytes();actual=(d/f).read_bytes()
            assert actual==expected,case+' '+f
            rows.append(dict(case=case,check=f,result='BIT_IDENTICAL',bytes=len(actual),sha256=hashlib.sha256(actual).hexdigest()))
    base=(ROOT/'finders/single-N48/params.txt').read_text()
    for label,override,expected in [('valid-proof',dict(map_find_max_updates=1,map_find_thresholds='1e-7 1e-7 1e-7'),0),
        ('mismatched-proof',dict(map_finder_bracket=ROOT/'seeds/initial-h0-o+0.0000-c1.00/bracket.csv'),2)]:
        d=ROOT/'verification'/label;d.mkdir(parents=True,exist_ok=True)
        (d/'params.txt').write_text(w.replace(base,override))
        rc=w.measure('regression',d,[exe,'params.txt'],3000000000)
        assert rc==expected
        if expected==0:assert 'T25_FINDER_COMPLETE FOUND_ALL_STAGES' in (d/'run.log').read_text()
        else:assert 'bracket does not match checkpoint/time/centre/seed or signs' in (d/'run.log').read_text()
        rows.append(dict(case=label,check='bracket guard',result='PASS',bytes=0,sha256=''))
    # The bracket postprocessor never interprets an unresolved surface.
    spec=importlib.util.spec_from_file_location('b',HERE/'t25-brackets.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
    sample=list(csv.DictReader((HERE/'t25-map-single-N192.csv').open()))
    candidates=b.pointwise_barriers(sample);assert len(candidates)==1
    pair=[dict(next(r for r in sample if r['id']==candidates[0][s])) for s in ['inner','outer']]
    assert len(b.brackets(pair))==1
    pair[0]['status']='UNRESOLVED';assert not b.brackets(pair) and not b.pointwise_barriers(pair)
    rows.append(dict(case='unresolved',check='interpretation guard',result='PASS',bytes=0,sha256=''))
    with (HERE/'t25-regression.csv').open('w') as f:
        out=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');out.writeheader();out.writerows(rows)
    data=dict(status='PASS',checks=rows,portable_executable=str(portable),portable_sha256=w.hashlib.sha256(portable.read_bytes()).hexdigest())
    (HERE/'t25-regression.json').write_text(json.dumps(data,indent=2)+'\n');print('T25_REGRESSION_PASS')
if __name__=='__main__':main()
