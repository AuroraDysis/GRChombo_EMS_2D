"""Prepare local frozen/t=0 controls from authenticated supplied inputs."""
from pathlib import Path
import hashlib
import json
from config import EMS, input_path

HERE=Path(__file__).resolve().parent
OUT=HERE/'results-t2'
SOURCE=EMS/'state/threads/ems-spectral-solver/runs/exp-0031/submissions/exp-0031/production/params-H15-L1680.txt'
CHECKPOINT=EMS/'.data/exp-0033/chk/EMS_001376.2d.hdf5'


def replace(text,values):
    present=set()
    lines=[]
    for line in text.splitlines():
        key=line.split('=',1)[0].strip()
        if key in values:
            line=f'{key} = {values[key]}';present.add(key)
        lines.append(line)
    lines += [f'{key} = {value}' for key,value in values.items() if key not in present]
    return '\n'.join(lines)+'\n'


def main():
    authenticated=json.loads((HERE/'results/inputs.json').read_text())
    for item in authenticated:
        path=input_path(item['path']);item['path']=str(path);digest=hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
        assert path.stat().st_size==item['bytes'] and digest.hexdigest()==item['sha256'],path
    OUT.mkdir(exist_ok=True)
    (OUT/'inputs.json').write_text(json.dumps(authenticated,indent=2)+'\n')
    common=SOURCE.read_text()
    for name,radii,np in [('checkpoint','50 75 100',True),
                         ('cavity','0.0512 0.2622 0.6069 1.0429',True),
                         ('default','50 75 100',False)]:
        folder=OUT/name
        for form in ('vacuum','matter'): (folder/form).mkdir(parents=True,exist_ok=True)
        values={'restart_file':CHECKPOINT,'ems_np_activate':str(np).lower(),
                'num_extraction_radii':len(radii.split()),'extraction_radii':radii,
                'extraction_levels':' '.join(['2']*len(radii.split()))}
        (folder/'params.txt').write_text(replace(common,values))
    template=(HERE.parent/'EMSRHFinder/results/offline/static/n96/params.txt').read_text()
    template='\n'.join(line for line in template.splitlines()
                       if not line.startswith(('restart_file','offline_')))+'\n'
    for name,N,L,radii,profile,f2,Q in [
        ('rn',256,16,'1 2 4','RN-identity',0,.7),
        ('ems',1024,4,'0.0512 0.2622 0.6069 1.0429','a0.9-e7',-.9,1.0681492364093355)]:
        for refinement in (1,2):
            folder=OUT/f'{name}-n{N*refinement}'
            folder.mkdir(parents=True,exist_ok=True)
            values={'N1':N*refinement,'N2':N*refinement//2,'L':L,'center':f'{L/2} 0',
                'star_centre':f'{L/2} 0','ems_data_path':EMS/f'artifacts/echo/{profile}.trumpet',
                'ems_f2':f2,'bh_charge':Q,'RH_activate':0,'ems_test_checkpoint':0,
                'num_plot_vars':28,'plot_vars':'chi h11 h12 h22 hww K A11 A12 A22 Aww Theta Gamma1 Gamma2 lapse shift1 shift2 B1 B2 phi Pi Lambda Bx By Bz Ex Ey Ez Xi',
                'activate_extraction':1,'ems_radiation_activate':'true','ems_np_activate':'true',
                'extraction_center':f'{L/2} 0','num_extraction_radii':len(radii.split()),
                'extraction_radii':radii,'extraction_levels':' '.join(['0']*len(radii.split())),
                'num_points_theta':97,'num_points_phi':64,'num_modes':5,'modes':'2 0 3 0 4 0 5 0 6 0',
                'sigma':1,'max_spatial_derivative_order':4,'num_ghosts':3,'verbosity':0,
                'output_path':'""','data_subpath':'.','extraction_subpath':'.','write_extraction':'false'}
            (folder/'params.txt').write_text(replace(template,values))
    assert replace('a = 1\nb = 2\n',{'a':3,'c':4})=='a = 3\nb = 2\nc = 4\n'


if __name__=='__main__':main()
