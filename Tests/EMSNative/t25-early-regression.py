#!/usr/bin/env python3
"""Verify unchanged native map values and the explicit weaker-bracket gate."""
import csv,hashlib,importlib.util,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('early',HERE/'t25-early.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def main():
    root=e.ROOT/'regression';root.mkdir(parents=True,exist_ok=True);results=[]
    d=root/'unchanged-map';d.mkdir(exist_ok=True)
    (d/'params.txt').write_text(Path('/private/tmp/ems-t25/maps/single-N192/params.txt').read_text())
    assert e.w.measure('early-regression',d,[e.EXE,'params.txt'],3000000000)==0
    e.verified(d,'early-regression',0,'T25_EXPANSION_MAP_COMPLETE')
    comparisons=0
    for file in ('map.csv','angles.csv'):
        old=e.read(Path('/private/tmp/ems-t25/maps/single-N192')/file);new=e.read(d/file)
        assert len(old)==len(new)
        for a,b in zip(old,new):
            for key,value in a.items():assert b[key]==value,(file,key,value,b[key]);comparisons+=1
    results.append(dict(check='historical native map and angular values',result='BIT_IDENTICAL',comparisons=comparisons))
    original=e.read('/private/tmp/ems-t25/seeds/single/bracket.csv')[0]
    weaker=dict(original,status='AVERAGE_ONLY_WEAKER');e.write(root/'weaker.csv',[weaker])
    base=Path('/private/tmp/ems-t25/finders/single-N48/params.txt').read_text()
    for flag,expected in [('false',2),('true',0)]:
        d=root/('average-'+flag);d.mkdir(exist_ok=True)
        (d/'params.txt').write_text(e.w.replace(base,dict(map_finder_bracket=root/'weaker.csv',map_find_allow_average=flag,
            map_find_max_updates=1,map_find_thresholds='1e-7 1e-7 1e-7')))
        rc=e.w.measure('early-regression',d,[e.EXE,'params.txt'],3000000000)
        e.verified(d,'early-regression',expected,'T25_FINDER_COMPLETE FOUND_ALL_STAGES' if expected==0 else 'explicitly permitted weaker mean bracket')
        assert rc==expected
        results.append(dict(check='weaker seed flag '+flag,result='PASS',comparisons=0))
    e.write(HERE/'t25-early-regression.csv',results)
    e.atomic_json(HERE/'t25-early-regression.json',dict(status='PASS',checks=results,
        executable_sha256=hashlib.sha256(e.EXE.read_bytes()).hexdigest()))
    print('T25_EARLY_REGRESSION_PASS',comparisons,flush=True)
if __name__=='__main__':main()
