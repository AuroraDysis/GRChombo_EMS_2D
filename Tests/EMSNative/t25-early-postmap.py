#!/usr/bin/env python3
"""Map the numerical centres reached by completed capped flows. No searches.

This preserves the original bracket/seed/235-s receipt, rather than claiming
the recorded search began at the refined centre. Spheres and one spheroid
estimated from the numerical final contour are mapped at N96/N192. Nothing
is seeded from a static profile, an area target or a CTT surface.
"""
import csv,importlib.util,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('early',HERE/'t25-early.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def main():
    assert json.loads((HERE/'t25-early-audit.json').read_text())['all_native_receipts_verified']
    summaries=[]
    for s in e.samples():
        d=e.ROOT/s['tag'];selected=json.loads((d/'selection.json').read_text())['selections']
        if not any(x['bracket'] for x in selected):continue
        print('POSTMAP_BEGIN',s['tag'],flush=True);specs=[]
        for x in selected:
            if not x['bracket']:continue
            fd=Path(x['seed_dir'])/'finder';last=e.read(fd/'finder.csv')[-1]
            shape=e.read(fd/('finder-shape-'+last['stage']+'.csv'))
            pole=sum(float(z['radius']) for z in shape[:2]+shape[-2:])/4
            equator=sum(float(z['radius']) for z in shape[len(shape)//2-1:len(shape)//2+1])/2
            aspects=sorted(set([1.,pole/equator]));centre=float(last['centre'])
            dense=[];r=max(3*e.H,float(last['r_min'])*.9)
            while r<min(.05,float(last['r_max'])*1.1):dense.append(r);r*=1.005
            for j,q in enumerate(aspects):
                family=f'post-h{x["hole"]}-q{j}'
                for i,a in enumerate(sorted(set(e.radii()+dense))):
                    if min(a,a*q)<3*e.H or max(a,a*q)>.05:continue
                    specs.append(dict(id=f'{family}-r{i:03}',family=family,centre=centre,a=a,c=a*q))
        maps={n:e.mapping(s,'post-flow-centres',specs,n) for n in (96,192)}
        e.write(HERE/f't25-early-{s["tag"]}-post-centres.csv',maps[192])
        e.write(HERE/f't25-early-{s["tag"]}-post-barriers.csv',e.strong(maps[192]))
        e.write(HERE/f't25-early-{s["tag"]}-post-brackets.csv',e.b.brackets(maps[192]))
        roots=[];pairs=[]
        for x in selected:
            if not x['bracket']:continue
            rs=[r for r in maps[192] if r['family'].startswith(f'post-h{x["hole"]}-')]
            cs=e.strong(rs)
            if cs:chosen=min(cs,key=lambda z:float(z['a_outer'])/float(z['a_inner']))
            else:
                cs=[z for z in e.b.brackets(rs) if z['theta_inner_mean']<0 and z['theta_outer_mean']>0]
                chosen=min(cs,key=lambda z:float(z['a_outer'])/float(z['a_inner'])) if cs else None
            if not chosen:
                summaries.append(dict(time=s['time'],hole=x['hole'],status='NO_POST_BRACKET',root=None));continue
            local=[z for z in e.b.brackets(rs) if z['family']==chosen['family'] and z['theta_inner_mean']<0 and z['theta_outer_mean']>0 and
                float(z['a_inner'])>=float(chosen['a_inner']) and float(z['a_outer'])<=float(chosen['a_outer'])]
            pair=min(local,key=lambda z:float(z['a_outer'])/float(z['a_inner'])) if local else chosen
            roots.append(dict(id=f'post-h{x["hole"]}-root',family=f'post-h{x["hole"]}-root',centre=chosen['centre'],
                a=pair['mean_linear_a_root'],c=pair['mean_linear_c_root']))
            pairs.append((x,chosen,pair))
        if roots:
            results={n:e.mapping(s,'post-flow-roots',roots,n) for n in (96,192)}
            e.write(HERE/f't25-early-{s["tag"]}-post-roots.csv',results[192])
            for (x,chosen,pair),root,root96 in zip(pairs,results[192],results[96]):
                assert root['status']==root96['status']=='RESOLVED'
                summaries.append(dict(time=s['time'],hole=x['hole'],status=chosen['status'],
                    fractional_width=float(chosen['a_outer'])/float(chosen['a_inner'])-1,
                    bracket=chosen,mean_pair=pair,root=root,root_N96=root96,
                    finder_seed_was_preliminary=True,new_finder_probes=0))
        print('POSTMAP_DONE',s['tag'],[(z['hole'],z['status'],z.get('fractional_width')) for z in summaries if z['time']==s['time']],flush=True)
    e.atomic_json(HERE/'t25-early-postmap.json',dict(status='COMPLETE',rows=summaries,new_finder_probes=0))
    e.summarize();print('T25_POSTMAP_COMPLETE; no searches, no advances',flush=True)
if __name__=='__main__':main()
