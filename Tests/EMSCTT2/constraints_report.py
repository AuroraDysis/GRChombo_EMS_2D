"""Summarize valid-cell composite norms, retaining all measured trends."""
from pathlib import Path
from evidence import ARCHIVE
import math,csv,json
R=Path(__file__).resolve().parent;out=R/'results'
rows=[];tables=[]
for stem in ('','near-'):
    previous={}
    lines=['| Region | N1 | fine dx/M | H RMS | order | M RMS | order | GE RMS | order | GB RMS |',
           '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for n in (80,160,320):
        text=(ARCHIVE/f'hierarchy-{stem}{n}.log').read_text()
        assert 'seconds=' in text
        for line in text.splitlines():
            if not line.startswith('NORM,'):continue
            a=line.split(',');name=a[1];count=int(a[2]);floored=int(a[3]);staging=float(a[4]);v=list(map(float,a[5:]))
            assert count and not floored and staging<=5e-13 and all(math.isfinite(x) for x in v)
            vals=[v[0],v[2],v[4],v[6]]
            orders=[math.log2(before/after) if before>0 and after>0 else None for before,after in zip(previous.get(name,vals),vals)] if name in previous else [None]*4
            previous[name]=vals
            dx=(160 if not stem else .32)/(2*n)
            row=dict(series=stem or 'global',N1=n,dx=dx,region=name,count=count,floored=floored,staging=staging)
            for k,q in enumerate(('H','M','GE','GB')):
                row.update({q+'_rms':v[2*k],q+'_max':v[2*k+1],q+'_order':orders[k]})
            rows.append(row)
            lines.append('| '+ ' | '.join([name,str(n),f'{dx:g}',*[s for k in range(3) for s in (f'{vals[k]:.5e}',f'{orders[k]:.3f}' if orders[k] is not None else '—')],f'{vals[3]:g}'])+' |')
    tables.append(('Near-hole' if stem else 'Global')+' hierarchy\n\n'+'\n'.join(lines))
with (out/'constraints.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys(),lineterminator="\n");writer.writeheader();writer.writerows(rows)
(out/'constraint-tables.md').write_text('\n\n'.join(tables)+'\n')
print('\n\n'.join(tables))
