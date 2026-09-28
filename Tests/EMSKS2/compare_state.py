"""Compare all checkpointed evolved variables, including ghost cells, bitwise."""
import csv
import hashlib
import sys

import h5py
import numpy as np


def hashes(handle, level, components):
    group = handle[f"level_{level}"]
    data = group["data:datatype=0"][:]
    offsets = group["data:offsets=0"][:]
    output = [hashlib.sha256() for _ in range(components)]
    cells = 0
    for start, end in zip(offsets[:-1], offsets[1:]):
        span = int(end - start)
        assert span % components == 0
        size = span // components
        cells += size
        for component, digest in enumerate(output):
            lo = int(start) + component * size
            digest.update(data[lo:lo + size].tobytes())
    return [item.hexdigest() for item in output], cells


def compare(left_path, right_path, label):
    with h5py.File(left_path) as left, h5py.File(right_path) as right:
        components = int(left.attrs["num_components"])
        levels = int(left.attrs["num_levels"])
        assert components == int(right.attrs["num_components"])
        assert levels == int(right.attrs["num_levels"])
        assert int(left.attrs["iteration"]) == int(right.attrs["iteration"])
        assert np.array_equal(left.attrs["time"], right.attrs["time"])
        writer = csv.writer(sys.stdout)
        writer.writerow(("case", "step", "level", "variable", "cells_with_ghosts",
                         "left_sha256", "right_sha256", "bit_identical"))
        for level in range(levels):
            a, b = left[f"level_{level}"], right[f"level_{level}"]
            for key in ("boxes", "data:offsets=0"):
                assert np.array_equal(a[key][:], b[key][:]), (level, key)
            lhs, cells = hashes(left, level, components)
            rhs, other_cells = hashes(right, level, components)
            assert cells == other_cells
            for component, (h1, h2) in enumerate(zip(lhs, rhs)):
                name = left.attrs[f"component_{component}"].decode()
                assert name == right.attrs[f"component_{component}"].decode()
                writer.writerow((label, int(left.attrs["iteration"]), level, name,
                                 cells, h1, h2, h1 == h2))


if __name__ == "__main__":
    compare(*sys.argv[1:4])
