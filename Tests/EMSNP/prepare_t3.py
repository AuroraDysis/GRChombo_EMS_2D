"""Reuse T2 inputs, preserving its evidence directory, and reauthenticate after runs."""
import hashlib
import json
from pathlib import Path
import sys
import prepare_t2 as base
from config import input_path
base.OUT=base.HERE/'results-t3'


def main():
    if '--verify-only' in sys.argv:
        expected=json.loads((base.OUT/'inputs.json').read_text())
        for item in expected:
            path=input_path(item['path']);item['path']=str(path);digest=hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
            assert path.stat().st_size==item['bytes'] and digest.hexdigest()==item['sha256'],path
        (base.OUT/'inputs-final.json').write_text(json.dumps(expected,indent=2)+'\n')
        print('PASS: all five input hashes unchanged after controls')
        return
    base.main()
    default=base.OUT/'default/params.txt'
    default.write_text('\n'.join(line for line in default.read_text().splitlines()
                               if not line.startswith('ems_np_activate'))+'\n')
    source=(base.OUT/'rn-n256/params.txt').read_text()
    (base.OUT/'bad-activation.params.txt').write_text(base.replace(source,{'ems_radiation_activate':'false'}))


if __name__=='__main__':main()
