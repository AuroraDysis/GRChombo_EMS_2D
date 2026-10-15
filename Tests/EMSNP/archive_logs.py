"""Move raw logs out of the tracked tree, retaining a verified SHA-256 manifest."""
import hashlib
import json
from pathlib import Path
import shutil
from config import LOG_ROOT

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
MANIFEST=HERE/'log-archive.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records=json.loads(MANIFEST.read_text()) if MANIFEST.exists() else []
    previous={d['original']:d for d in records}
    for path in sorted(HERE.rglob('*.log')):
        relative=str(path.relative_to(ROOT));target=LOG_ROOT/relative
        sha=digest(path)
        if target.exists():assert digest(target)==sha, target
        else:
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(path,target)
        assert digest(target)==sha, target
        previous[relative]={'original':relative,'archive':relative,
            'bytes':path.stat().st_size,'lines':len(path.read_bytes().splitlines()),'sha256':sha}
        MANIFEST.write_text(json.dumps(sorted(previous.values(),key=lambda d:d['original']),indent=2)+'\n')
        path.unlink()
    records=sorted(previous.values(),key=lambda d:d['original'])
    MANIFEST.write_text(json.dumps(records,indent=2)+'\n')
    for d in records:
        path=LOG_ROOT/d['archive']
        assert path.stat().st_size==d['bytes'] and digest(path)==d['sha256'],path
    print(f'PASS: {len(records)} archived logs, {sum(d["bytes"] for d in records)} bytes, all hashes verified')


if __name__=='__main__':main()
