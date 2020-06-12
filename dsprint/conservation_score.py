import os.path
import pickle
import gzip
import numpy as np
import bisect
from concurrent.futures import ThreadPoolExecutor, as_completed


class WigFix:

    def __init__(self, filename=None, index_file=None, dtype='float16'):
        if index_file is not None and os.path.exists(index_file):
            with open(index_file, 'rb') as f:
                self.positions, self.values = pickle.load(f)
            return

        assert filename is not None, 'Wigfix gz filename must be specified if not using an index file'

        positions = []
        last_position = 0
        values = {}

        with gzip.open(filename, 'rt') as f:
            line = f.readline()
            while line != '':
                if line.startswith('fixedStep'):
                    last_position = int(line.split(' ')[2].split('=')[1])
                    values[last_position] = []
                    positions.append(last_position)
                else:
                    values[last_position].append(float(line))
                line = f.readline()

        values = {k: np.array(v).astype(dtype) for k, v in values.items()}
        self.positions = positions
        self.values = values

        if index_file is not None:
            with open(index_file, 'wb') as f:
                pickle.dump((self.positions, self.values), f)

    def __getitem__(self, item):
        position = bisect.bisect_right(self.positions, item)
        # The return position from bisect_right is the insert position
        # This is 0 for elements < the first that we have, 1 between [<first>, <second>)
        # Subtract 1 to get the index where we can start our forward search
        if position < 1:
            return None
        start_position = self.positions[position - 1]
        return self.values[start_position][item - start_position]


class ConservationScore:

    def __init__(self, name, dtype='float16'):
        self.name = name
        self.dtype = dtype
        self.wigfixes = {}

    def _load(self, wigfix_file, index_file=None):
        print(f'Loading {wigfix_file}')
        return WigFix(wigfix_file, index_file=index_file, dtype=self.dtype)

    def load(self, chromosomes, wigfix_files, index_files=None, max_workers=2):
        assert len(chromosomes) == len(wigfix_files), 'Specify a wigfix gz file for each chromosome'
        if index_files is not None:
            assert len(chromosomes) == len(index_files), 'Specify an index file for each chromosome'
        else:
            index_files = [None] * len(chromosomes)

        with ThreadPoolExecutor(max_workers=min(max_workers, len(chromosomes))) as executor:
            futures = {executor.submit(self._load, wigfix_file, index_file): chromosome
                       for chromosome, wigfix_file, index_file in zip(chromosomes, wigfix_files, index_files)}
            for future in as_completed(futures):
                chromosome = futures[future]
                self.wigfixes[chromosome] = future.result()

    def __getitem__(self, item):
        return self.wigfixes[item]

    def mem_usage(self):
        # Return approx. mem usage of this object in MB
        return sum(x.nbytes for w in self.wigfixes for x in self[w].values.values()) / 1024. / 1024.
