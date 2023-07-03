import pickle
from .quarry import open_default_json
import json
import os


def benchmark_runtime(quarry, N, output_file):
    data_out = {}
    with open_default_json() as default_file:
        default_json = json.load(default_file)
    quarry.set_value(["settings"], default_json)

    for n in range(N):
        quarry.clear_cache()
        quarry.update_all(lambda s: None)
        runtimes = quarry.get_runtimes()
        for k, v in runtimes:
            if k not in data_out:
                data_out[k] = []
            data_out[k].append(v)

    with open(output_file, "wb") as f:
        pickle.dump(data_out, f)
