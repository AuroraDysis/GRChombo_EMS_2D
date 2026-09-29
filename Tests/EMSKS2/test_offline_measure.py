"""Small quadrature and CSV aggregation check; run with python3."""
import csv
import math
import tempfile
from pathlib import Path

from offline_measure import collect, sinc_correction


def test():
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)/'work'
        case = work/'0000-snapshot'
        for n in (48, 96):
            directory = case/f'n{n}'
            directory.mkdir(parents=True)
            correction = sinc_correction(n)
            with (directory/'surfaces.csv').open('w', newline='') as file:
                writer = csv.DictWriter(file, ['time','N_theta','stage','A','Q','area_delta','charge_delta','status'])
                writer.writeheader()
                writer.writerow(dict(time=0,N_theta=n,stage=2,A=4*math.pi/correction,
                                     Q=6/correction,area_delta=0,charge_delta=0,status='FOUND'))
            with (directory/'exterior.csv').open('w', newline='') as file:
                names = ['time','N_theta','kind','target_over_rh','R_mean','Q_S','m_mean','m_static',
                         'phi_mean','phi_static','lapse_mean','lapse_static','K_mean','K_static','Q_static']
                writer = csv.DictWriter(file, names)
                writer.writeheader()
                writer.writerow(dict(time=0,N_theta=n,kind='coordinate',target_over_rh=2,R_mean=2,
                                     Q_S=6/correction,m_mean=1,m_static=1,phi_mean=0,phi_static=0,
                                     lapse_mean=1,lapse_static=1,K_mean=0,K_static=0,Q_static=6))
        result = collect(work, Path(temp)/'out.csv')
        horizons = [row for row in result if row['kind']=='horizon']
        assert len(result)==4 and len(horizons)==2
        assert all(abs(row['A_corrected']-4*math.pi)<1e-12 and
                   abs(row['Q_corrected']-6)<1e-12 for row in horizons)
        assert all(row['angular_delta_Q_corrected']<1e-12 for row in result)
        assert all(row['coordinate_R_drift_rel']==0 for row in result if row['kind']=='exterior_coordinate')


if __name__ == '__main__':
    test()
