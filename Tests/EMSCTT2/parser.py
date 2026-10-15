from pathlib import Path
from evidence import RAW, OUT
import subprocess,hashlib,json,os
R=Path(__file__).resolve().parent
out=R/'results'; source=(out/'input-snapshots/n12-final.ctt').read_text()
exe=str(next(R.glob('EMSCTT2Readback*.ex')))
profile=str(out/'input-snapshots/a0.9-e7.trumpet')
payload=source[:source.rindex('SHA256 ')]
def change(s,a,b):
    assert a in s,a
    return s.replace(a,b,1)
def seal(s):return s+'SHA256 '+hashlib.sha256(s.encode()).hexdigest()+'\n'
mutants={
 'magic':change(payload,'EMSCTT 2','EMSCTT 3'),
 'CRLF magic':change(payload,'EMSCTT 2\n','EMSCTT 2\r\n'),
 'unknown':change(payload,'# mass=1.0','# unknown=1\n# mass=1.0'),
 'duplicate':change(payload,'# mass=1.0','# mass=1.0\n# mass=1.0'),
 'missing':change(payload,'# mass=1.0\n',''),
 'geometry':change(payload,'one-unboosted-trumpet','binary'),
 'schema':change(payload,'EMSCTT/2-single-pulse','EMSCTT/1'),
 'convention':change(payload,'ems-native-draft-v1','other'),
 'free data':change(payload,'fixed-D0','fixed-E0'),
 'mass':change(payload,'# mass=1.0','# mass=2.0'),
 'sigma':change(payload,'# sigma_M=3.0','# sigma_M=0'),
 'nonfinite':change(payload,'# amplitude=0.004','# amplitude=NaN'),
 'nonconverged':change(payload,'# newton_residual=2.1760371282653068e-14','# newton_residual=1'),
 'integer':change(payload,'# panel_count=22','# panel_count=22.0'),
 'tail integer':change(payload,'# tail_L=18','# tail_L=18.5'),
 'hash syntax':change(payload,'# solver_sha256=c','# solver_sha256=x'),
 'profile binding':change(payload,'# profile_sha256=7','# profile_sha256=0'),
 'phase normalization':change(payload,'21.10140303248864','22.10140303248864'),
 'phase size':change(payload,'62.06455126251464 66','62.06455126251464 17'),
 'truncation':payload[:payload.index('VALUES')+10],
 'line-boundary truncation':payload[:payload.index('VALUES')],
 'END':change(payload,'\nEND\n','\nEND2\n'),
 'trailing':payload+'unexpected\n',
}
# Mutate actual panel rows without depending on source coefficient counts.
lines=payload.splitlines(keepends=True);panel=next(i for i,l in enumerate(lines) if l.startswith('PANEL '));row=lines[panel].split()
for name,index,value in [('map',1,'bridge'),('gap',3,str(float(row[3])*1.01)),('dimensions',4,'2'),('noninteger dimensions',4,'3.5')]:
    a=row.copy();a[index]=value
    s=lines.copy();s[panel]=' '.join(a)+'\n';mutants[name]=''.join(s)
s=lines.copy();s[panel+2]='Inf\n';mutants['nonfinite nodal']=''.join(s)
s=lines.copy();s[panel+1]='VALUES 1\n';mutants['block size']=''.join(s)
results={}
def check(name,s,expected):
    path=RAW/'parser-mutant.ctt';path.write_text(s)
    p=subprocess.run([exe,profile,str(path),'--check'],capture_output=True,text=True,env={**os.environ,'OMP_NUM_THREADS':'1'},timeout=10)
    assert p.returncode==expected,(name,p.returncode,p.stdout,p.stderr)
    results[name]=p.stdout.splitlines()[-1]
    if name=='CRLF magic': assert 'version line contains CR; use LF line endings' in results[name]
check('valid',source,0)
for name,s in mutants.items():check(name,seal(s),2)
check('checksum',source.replace('amplitude=0.004','amplitude=0.005'),2)
check('final newline',source.rstrip('\n'),2)
(RAW/'parser-mutant.ctt').unlink()
(OUT/'parser.json').write_text(json.dumps(results,indent=2)+'\n')
print(f'PASS valid + {len(results)-1} strict-parser rejection controls')

# Constructor binding and default-lapse guards are distinct from file parsing.
bindings={}
for flag in ['binary','boosted','rapidity','mass','separation','geometric-lapse','maximal-lapse','coupling','puncture','below-min','bad-panel','negative-panel']:
    p=subprocess.run([exe,profile,str(out/'input-snapshots/n12-final.ctt'),'--'+flag],text=True,capture_output=True,timeout=10,env={**os.environ,'OMP_NUM_THREADS':'2'})
    assert p.returncode!=0,(flag,p.stdout,p.stderr)
    expected=('outside represented end' if flag in ('puncture','below-min') else 'invalid panel/domain') if flag in ('puncture','below-min','bad-panel','negative-panel') else ('EMSTRUMPET coupling differs' if flag=='coupling' else ('requires an unboosted single hole' if flag=='maximal-lapse' else 'EMSCTT 2 requires'))
    assert expected in p.stdout+p.stderr,(flag,p.stdout,p.stderr)
    bindings[flag]={'exit':p.returncode,'output':(p.stdout+p.stderr).strip()}
(OUT/'bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
print('PASS eight binding rejections + four OpenMP fatal domain controls')
