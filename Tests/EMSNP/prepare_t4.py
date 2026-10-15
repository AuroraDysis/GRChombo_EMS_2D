"""Authenticate the supplied checkpoint inputs and prepare the order study."""
import hashlib
import json
from pathlib import Path
import prepare_t2 as base
from config import input_path

OUT=base.HERE/'results-t4'


def authenticate():
    inputs=json.loads((base.HERE/'results-t3/inputs.json').read_text())
    for item in inputs:
        path=input_path(item['path']);item['path']=str(path);digest=hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
        assert path.stat().st_size==item['bytes'] and digest.hexdigest()==item['sha256'],path
    return inputs


def main():
    OUT.mkdir(exist_ok=True)
    (OUT/'inputs.json').write_text(json.dumps(authenticate(),indent=2)+'\n')
    source=(base.HERE/'results-t3/checkpoint/params.txt').read_text()
    for name,bounded in [('study',True),('insufficient-ghosts',False)]:
        folder=OUT/name; folder.mkdir(exist_ok=True)
        values={'restart_file':base.CHECKPOINT,'np_order_study':'true','np_order_repeats':3,
            'num_ghosts':3,'np_order_bounded':str(bounded).lower()}
        (folder/'params.txt').write_text(base.replace(source,values))
    contract={'target_order':4,'comparison_orders':[2,4,6],'safety_factor':2,
        'estimate':'abs(NP4-NP6) + abs(W4-W6)',
        'basis':'sum of within-route 4-to-6 order changes for two discretizations',
        'acceptance':'abs(NP4-W4) <= safety_factor * estimate',
        'radii':[50,75,100],'modes':[2,4,6],
        'production_cell_derivatives':'fourth order at every interpolation order',
        'ghosts':3,'production_default_ghosts':3,'theta_points':97,
        'stencil_policy':'existing Lagrange full-width search bounded to the sampled FAB',
        'negative_control':'three ghosts must reject a L6 stencil outside its FAB'}
    (OUT/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    print('PASS: five authenticated inputs; order-error contract written before study')


if __name__=='__main__':main()
