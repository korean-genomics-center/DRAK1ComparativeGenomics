import json
import os
import re

import numpy as np
from joblib import Parallel, delayed


class BlockJsonReader:
    def __init__(self, path_file, offset_start, offset_end, key_set=None):
        self.path = path_file
        self.offset_start = offset_start
        self.offset_end = offset_end
        self.key_set = key_set if key_set else set()
        self.results = {}

    def extract_matching_pairs(self):
        with open(self.path, 'rb') as f:
            f.seek(self.offset_start)
            data = f.read(self.offset_end - self.offset_start + 1024 * 1024)  # slight over-read
            text = data.decode(errors='ignore')

            # Find "key": value (value can be string, number, null, boolean)
            pattern = re.compile(r'"([^"]+)"\s*:\s*(".*?"|\d+\.?\d*|true|false|null)', re.DOTALL)
            for match in pattern.finditer(text):
                key, val_str = match.groups()
                if key in self.key_set:
                    try:
                        value = json.loads(val_str)
                        self.results[key] = value
                    except json.JSONDecodeError:
                        continue
        return self.results


class BlockJsonSplitter:
    def __init__(self, path_file, parallel=4):
        self.path = path_file
        self.parallel = parallel
        self.file_size = os.path.getsize(self.path)
        self.offsets_start = []
        self.offsets_end = []

    def split(self):
        self.offsets_start = list(map(int, np.linspace(0, self.file_size, self.parallel + 1)))[:-1]
        self.offsets_end = self.offsets_start[1:] + [self.file_size]

    def get_block_readers(self, key_set):
        return [
            BlockJsonReader(self.path, start, end, key_set)
            for start, end in zip(self.offsets_start, self.offsets_end)
        ]


def process_reader(reader):
    return reader.extract_matching_pairs()


def run_parallel_json_extract(path_json, protein_id_list, n_jobs=8):
    splitter = BlockJsonSplitter(path_json, parallel=n_jobs)
    splitter.split()

    readers = splitter.get_block_readers(set(protein_id_list))

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_reader)(reader) for reader in readers
    )

    final_result = {}
    for d in results:
        final_result.update(d)
        
    return final_result