"""Strict output gate: complete rows, finite values, absolute tolerance 1e-12.

uv run --with numpy --with h5py python compare_perf.py reference candidate [--mpi]
Only --mpi allows checkpoint Processors datasets to differ. HDF5 object creation
timestamps are not physical data; all datasets and attributes are checked.
"""
import json
from pathlib import Path
import sys

import h5py
import numpy as np


def compare(reference, candidate, mpi=False):
    reference, candidate = Path(reference), Path(candidate)
    patterns = ('rh_surf_*.dat', 'rh_f*.dat', '*extraction*.dat', '*Extraction*.dat',
                '*mode*.dat', 'ems_radiation.csv', '*.hdf5',
                'constraint_norms.dat', 'M_ADM.dat', 'min_chi.dat')
    files = lambda root: {p.relative_to(root) for pattern in patterns
                          for p in root.rglob(pattern)}
    expected = files(reference)
    assert expected and expected == files(candidate), 'output file coverage differs'
    numeric = 0
    max_error = 0.

    def equal(a, b, where):
        nonlocal numeric, max_error
        a, b = np.asarray(a), np.asarray(b)
        assert a.shape == b.shape and a.dtype == b.dtype, where
        if a.dtype.names:
            for key in a.dtype.names:
                equal(a[key], b[key], f'{where}/{key}')
        elif np.issubdtype(a.dtype, np.number):
            assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), where
            error = float(np.max(np.abs(a.astype(float)-b.astype(float)), initial=0))
            assert error <= 1e-12, (where, error)
            numeric += a.size
            max_error = max(max_error, error)
        else:
            assert np.array_equal(a, b), where

    for name in sorted(expected):
        a, b = reference/name, candidate/name
        if a.suffix == '.hdf5':
            with h5py.File(a) as fa, h5py.File(b) as fb:
                na, nb = [], []
                fa.visit(na.append); fb.visit(nb.append)
                assert na == nb, name
                for key in [''] + na:
                    x, y = fa[key or '/'], fb[key or '/']
                    assert set(x.attrs) == set(y.attrs), (name, key)
                    for attr in x.attrs:
                        equal(x.attrs[attr], y.attrs[attr], f'{name}/{key}@{attr}')
                    if isinstance(x, h5py.Dataset):
                        if mpi and key.endswith('/Processors'):
                            continue
                        equal(x[()], y[()], f'{name}/{key}')
        else:
            # Keep every row, including repeated RH times: no dictionary dedup.
            ra, rb = a.read_text().splitlines(), b.read_text().splitlines()
            assert len(ra) == len(rb), name
            for line, (x, y) in enumerate(zip(ra, rb)):
                xx, yy = x.replace(',', ' ').split(), y.replace(',', ' ').split()
                assert len(xx) == len(yy), (name, line)
                for i, (u, v) in enumerate(zip(xx, yy)):
                    try:
                        uf, vf = float(u), float(v)
                    except ValueError:
                        assert u == v, (name, line, i)
                    else:
                        equal(uf, vf, f'{name}:{line}:{i}')
    result = dict(files=len(expected), numeric_values=numeric, max_abs_error=max_error,
                  mpi_processors_ignored=mpi)
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    compare(sys.argv[1], sys.argv[2], '--mpi' in sys.argv[3:])
