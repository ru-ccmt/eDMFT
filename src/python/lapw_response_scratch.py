"""Rank-private, ephemeral storage for immutable LAPW response intermediates."""
from dataclasses import dataclass
import os
from tempfile import TemporaryFile

import numpy as np


@dataclass(frozen=True)
class StoredMatrix:
    file: object
    offset: int
    shape: tuple
    dtype: np.dtype
    nbytes: int

    def read(self):
        # Independent of the append position and safe for concurrent readouts.
        data = os.pread(self.file.fileno(), self.nbytes, self.offset)
        if len(data) != self.nbytes:
            raise OSError("incomplete LAPW response scratch read")
        return np.frombuffer(data, self.dtype).reshape(self.shape)


class MatchingScratch:
    """One anonymous file per spin branch and MPI rank, closed with the native state."""

    def __init__(self, stack):
        self.file = stack.enter_context(TemporaryFile(mode="w+b", buffering=0))

    def store(self, value):
        array = np.ascontiguousarray(value)
        offset = self.file.tell()
        if array.nbytes and self.file.write(memoryview(array).cast("B")) != array.nbytes:
            raise OSError("incomplete LAPW response scratch write")
        return StoredMatrix(self.file, offset, array.shape, array.dtype, array.nbytes)
