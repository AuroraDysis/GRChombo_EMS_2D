#!/usr/bin/env python3
"""T13 serial bounded runner, /usr/bin/time -l and atomic done.exit."""
import argparse, ctypes, json, os, resource, signal, subprocess, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ.get('T13_OUTPUT_ROOT','/private/tmp/ems-t13'))
CAP=int(os.environ.get('T13_PROCESS_CAP_BYTES','6000000000'))
TREE_CAP=int(os.environ.get('T13_TREE_CAP_BYTES','7750000000'))
DISK_CAP=int(os.environ.get('T13_DISK_CAP_BYTES','5700000000'))
LIB=ctypes.CDLL('/usr/lib/libproc.dylib',use_errno=True)

def tree_usage(root):
    pending=[root];usage={}
    while pending:
        pid=pending.pop();a=(ctypes.c_uint64*32)()
        if LIB.proc_pid_rusage(pid,0,ctypes.byref(a))==0:
            usage[pid]=(int(a[8]),int(a[9])) # resident size, physical footprint
        children=(ctypes.c_int*128)()
        LIB.proc_listchildpids(pid,children,ctypes.sizeof(children))
        pending.extend(p for p in children if p>1 and p not in usage)
    return usage

def measure(label,command,directory):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    # Repeated analyses retain every observed peak rather than overwriting it.
    if (directory/(label+'.resources.json')).exists():
        index=0
        while (directory/(f'{label}-repeat-{index:03}.resources.json')).exists():index+=1
        for suffix in ('.resources.json','.time','.child.json'):
            old=directory/(label+suffix)
            if old.exists():old.rename(directory/(f'{label}-repeat-{index:03}'+suffix))
    env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',
             MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',
             T13_CURRENT_MEASURE=label,
             CHOMBO_HOME='/Users/auroradysis/Workspace/EMS-deps/Chombo/lib')
    stamp=directory/(label+'.time')
    child_meta=directory/(label+'.child.json')
    if child_meta.exists():child_meta.unlink()
    start=time.monotonic()
    p=subprocess.Popen(['/usr/bin/time','-l','-o',str(stamp),sys.executable,
                        str(Path(__file__).resolve()),'--child-metadata',str(child_meta),'--',*command],
                       cwd=directory,env=env,start_new_session=True)
    peaks={};tree_peak=0;reason='';disk_check=start
    while p.poll() is None:
        usage=tree_usage(p.pid)
        total=sum(max(v) for v in usage.values())+resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        tree_peak=max(tree_peak,total)
        for pid,v in usage.items():
            old=peaks.get(pid,(0,0));peaks[pid]=tuple(max(a,b) for a,b in zip(old,v))
        # The serial queue's inner job watchdog owns termination/markers. Leave
        # 250 MB for its idle ancestors and avoid racing two disk-gate killers.
        if label!='queue' and (any(max(v)>CAP for v in usage.values()) or total>TREE_CAP):
            reason=f'memory gate: {CAP} per process / {TREE_CAP} tree bytes'
        if label!='queue' and time.monotonic()-disk_check>10:
            disk_check=time.monotonic()
            if sum(q.stat().st_size for q in ROOT.rglob('*') if q.is_file())>DISK_CAP:
                reason=f'output ceiling {DISK_CAP} bytes'
        if reason:
            os.killpg(p.pid,signal.SIGTERM);break
        time.sleep(.25)
    time_rc=p.wait()
    child=json.loads(child_meta.read_text()) if child_meta.exists() else {}
    rc=child.get('returncode',time_rc)
    # A gate kill can prevent the child wrapper from writing its wait4 result.
    peak=max(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
             child.get('peak_rss_bytes',0),max((v[0] for v in peaks.values()),default=0))
    row=dict(process=label,command=json.dumps(command),directory=str(directory),
             wall_seconds=time.monotonic()-start,peak_rss_bytes=peak,
             rss_GB=peak/1e9,cap_GB=CAP/1e9,returncode=rc,gate_reason=reason,
             tree_peak_bytes=tree_peak,per_pid_peaks=peaks,
             time_returncode=time_rc,child_measurement=child,
             measurement='time -l; wait4 RUSAGE_CHILDREN fallback if sandbox sysctl denied')
    (directory/(label+'.resources.json')).write_text(json.dumps(row,indent=2)+'\n')
    if peak>=CAP: raise RuntimeError('T13 memory cap violated')
    return 125 if reason else rc

def worker(plan_path):
    plan=json.loads(Path(plan_path).read_text())
    for job in plan['jobs']:
        d=Path(job['directory']); d.mkdir(parents=True,exist_ok=True)
        if (d/'done.exit').exists(): continue
        with (d/'run.log').open('wb') as log:
            p=subprocess.run([sys.executable,str(Path(__file__).resolve()),
                  '--measure',job['name'],'--directory',str(d),'--',*job['command']],
                  stdout=log,stderr=subprocess.STDOUT)
        tmp=d/'done.exit.tmp';tmp.write_text(str(p.returncode)+'\n');tmp.replace(d/'done.exit')
        if p.returncode: return p.returncode
    return 0

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--detach');ap.add_argument('--worker')
    ap.add_argument('--child-metadata')
    ap.add_argument('--measure');ap.add_argument('--directory',default=str(ROOT))
    ap.add_argument('command',nargs=argparse.REMAINDER);a=ap.parse_args()
    if a.child_metadata:
        cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
        rc=subprocess.run(cmd).returncode
        Path(a.child_metadata).write_text(json.dumps(dict(returncode=rc,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            wrapper_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))+'\n')
        return rc
    if a.detach:
        p=Path(a.detach).resolve();d=p.parent
        with (d/'launcher.log').open('ab') as log:
            child=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--measure','queue',
                         '--directory',str(d),'--',sys.executable,str(Path(__file__).resolve()),'--worker',str(p)],
                         stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        (d/'launcher.pid').write_text(str(child.pid)+'\n');print('detached PID',child.pid,flush=True)
        return 0
    if a.worker:
        rc=worker(a.worker);d=Path(a.worker).parent
        temp=d/'done.exit.tmp';temp.write_text(str(rc)+'\n');temp.replace(d/'done.exit')
        return rc
    if a.measure:
        cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
        assert cmd,'missing command';return measure(a.measure,cmd,a.directory)
    ap.error('use --measure, --detach or --worker')

if __name__=='__main__':sys.exit(main())
