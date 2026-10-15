"""Prepare/verify the frozen puncture-crossing non-fatal extraction control."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from prepare_t2 import replace

HERE=Path(__file__).resolve().parent
OUT=HERE/'results-t5'


def main():
    folder=OUT/'puncture-test'
    if '--prepare' in sys.argv:
        inputs=json.loads((HERE/'results-t4/inputs.json').read_text())
        for item in inputs:
            path=Path(item['path']); digest=hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
            assert digest.hexdigest()==item['sha256'],path
        OUT.mkdir(exist_ok=True)
        (OUT/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
        source=(HERE/'results-t3/checkpoint/params.txt').read_text()
        # The remnant puncture is at the symmetry origin x=840,y=0.
        # This displaced sphere passes through it; nearby angular nodes hit the floor.
        text=replace(source,{'extraction_center':'840.0512 0','num_extraction_radii':2,
            'extraction_radii':'0.0512 50','extraction_levels':'2 2',
            'np_validation_points':'false'})
        folder.mkdir(exist_ok=True)
        for form in ('vacuum','matter','sixth-vacuum','sixth-matter'):
            (folder/form).mkdir(exist_ok=True)
        (folder/'params.txt').write_text(text)
        print('PASS input hashes and non-fatal test preparation')
        return
    for form in ('vacuum','matter'):
        with (folder/form/'ems_np.csv').open() as stream:rows=list(csv.DictReader(stream))
        invalid=[d for d in rows if float(d['r'])==.0512]
        valid=[d for d in rows if float(d['r'])==50]
        assert len(invalid)==len(valid)==7
        for d in invalid:
            assert all(math.isnan(float(v)) for k,v in d.items() if k not in ('t','r','l','m'))
        for d in valid:assert all(math.isfinite(float(v)) for v in d.values())
    logroot=Path(os.environ.get('EMS_NP_LOG_ROOT','/Users/auroradysis/Workspace/EMS/.data/fork-np-logs'))
    log=(logroot/'t5/nonfatal.log').read_text()
    warnings=[line for line in log.splitlines() if 'EMS NP warning: invalid sphere' in line]
    assert len(warnings)==2 and all('t=602' in line for line in warnings)
    result={'status':'PASS','invalid_sphere':.0512,'finite_sphere':50,
        'forms':2,'nan_rows':14,'finite_rows':14,'warnings':len(warnings),'exit':0}
    (OUT/'nonfatal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
