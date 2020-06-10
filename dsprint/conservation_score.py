import gzip
import numpy as np
import bisect
from concurrent.futures import ThreadPoolExecutor, as_completed

from dsprint.two_stage_runner import TwoStageRunner
from dsprint.core import CHROMOSOMES


class WigFix:

    def __init__(self, filename, dtype='float16'):
        self.filename = filename
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

    def __getitem__(self, item):
        position = bisect.bisect_right(self.positions, item)
        # The return position from bisect_right is the insert position
        # This is 0 for elements < the first that we have, 1 between [<first>, <second>)
        # Subtract 1 to get the index where we can start our forward search
        if position < 1:
            return None
        start_position = self.positions[position - 1]
        return self.values[start_position][item - start_position]


class ConservationScore(TwoStageRunner):

    def __init__(self, file_pattern, dtype='float16', *args, **kwargs):
        super(ConservationScore, self).__init__(*args, **kwargs)
        self.file_pattern = file_pattern
        self.dtype = dtype
        self.wigfixes = {}

    def _load1(self, chromosome):
        print(f'Loading chromosome {chromosome}')
        file_path = self.file_pattern.format(chromosome=chromosome)
        self.wigfixes[chromosome] = WigFix(file_path, dtype=self.dtype)

    def _stage1(self, chromosomes=None, max_workers=2):
        chromosomes = chromosomes or CHROMOSOMES

        with ThreadPoolExecutor(max_workers=min(max_workers, len(chromosomes))) as executor:
            futures = {executor.submit(self._load1, chromosome): chromosome
                       for chromosome in chromosomes}
            for future in as_completed(futures):
                print(f'Finished processing chromosome {futures[future]}')

    def __getitem__(self, item):
        return self.wigfixes[item]

    def mem_usage(self):
        # Return approx. mem usage of this object in MB
        return sum(x.nbytes for w in self.wigfixes for x in self[w].values.values()) / 1024. / 1024.
