#!/usr/bin/env python3
"""Derive sampled pointwise or weaker average barriers from native map CSVs.

CLI: python t25-brackets.py map.csv brackets.csv
Groups must share checkpoint, time, centre, aspect ratio and angular N. Only
consecutive, resolved, strictly nested surfaces qualify. UNRESOLVED/invalid
rows break adjacency; they are never silently skipped to manufacture a pair.
"""
import argparse,csv,math
from collections import defaultdict
from pathlib import Path
def brackets(rows):
    groups=defaultdict(list);out=[]
    for r in rows:groups[(r['checkpoint'],r['time'],r['family'],r['N'])].append(r)
    for key,rs in groups.items():
        rs.sort(key=lambda r:float(r['a']))
        for x,y in zip(rs,rs[1:]):
            assert float(x['centre'])==float(y['centre'])
            assert math.isclose(float(x['c'])/float(x['a']),float(y['c'])/float(y['a']),rel_tol=1e-12)
            assert float(x['a'])<float(y['a']) and float(x['c'])<float(y['c'])
            if x['status']!='RESOLVED' or y['status']!='RESOLVED':continue
            xm,xM,ym,yM=[float(r[c]) for r,c in [(x,'theta_min'),(x,'theta_max'),(y,'theta_min'),(y,'theta_max')]]
            mx,my=float(x['theta_mean']),float(y['theta_mean'])
            if xM<0 and ym>0:status='POINTWISE_NEGATIVE_TO_POSITIVE'
            elif xm>0 and yM<0:status='POINTWISE_POSITIVE_TO_NEGATIVE'
            elif mx*my<0:status='AVERAGE_ONLY_WEAKER'
            else:continue
            f=-mx/(my-mx);root=float(x['a'])+f*(float(y['a'])-float(x['a']))
            out.append(dict(checkpoint=key[0],time=key[1],family=key[2],N=key[3],centre=x['centre'],
                status=status,inner=x['id'],outer=y['id'],a_inner=x['a'],a_outer=y['a'],c_inner=x['c'],c_outer=y['c'],
                theta_inner_min=xm,theta_inner_max=xM,theta_outer_min=ym,theta_outer_max=yM,
                theta_inner_mean=mx,theta_outer_mean=my,mean_linear_a_root=root,
                mean_linear_c_root=root*float(x['c'])/float(x['a']),
                A_linear=float(x['A'])+f*(float(y['A'])-float(x['A'])),
                Q_linear=float(x['Q'])+f*(float(y['Q'])-float(x['Q']))))
    return out
def pointwise_barriers(rows):
    """Tight uniform-sign barriers, allowing intervening mixed-sign surfaces.

    Dense radial sampling can insert a mixed-sign surface between the two
    uniform signs. That must not be mistaken for absence of a barrier pair.
    Keep all input surfaces and report whether the pair was adjacent.
    """
    groups=defaultdict(list);out=[]
    for r in rows:groups[(r['checkpoint'],r['time'],r['family'],r['N'])].append(r)
    for rs in groups.values():
        rs.sort(key=lambda r:float(r['a']))
        for j,y in enumerate(rs):
            if y['status']!='RESOLVED' or float(y['theta_min'])<=0:continue
            # All intervening surfaces must still have valid, resolved geometry.
            xs=[]
            for i in range(j-1,-1,-1):
                x=rs[i]
                if x['status']!='RESOLVED':break
                if float(x['theta_max'])<0:xs=[(i,x)];break
            if not xs:continue
            i,x=xs[0]
            r=brackets([x,y])[0];r['adjacent_in_input']=int(j==i+1);r['intervening_mixed_surfaces']=j-i-1
            out.append(r);break # tight first negative-to-positive barrier per family
    return out
def write(path,rows):
    fields=list(rows[0]) if rows else ['checkpoint','time','family','N','centre','status','inner','outer','a_inner','a_outer','c_inner','c_outer',
        'theta_inner_min','theta_inner_max','theta_outer_min','theta_outer_max','theta_inner_mean','theta_outer_mean','mean_linear_a_root','mean_linear_c_root','A_linear','Q_linear']
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--barriers');a=p.parse_args()
    source=list(csv.DictReader(Path(a.input).open()));rows=brackets(source);write(a.output,rows)
    if a.barriers:write(a.barriers,pointwise_barriers(source))
    from collections import Counter
    print(dict(Counter(r['status'] for r in rows)))
