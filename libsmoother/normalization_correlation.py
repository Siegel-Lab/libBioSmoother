from scipy.stats.stats import pearsonr
from .quarry import open_default_json
import json
import math

def get_correlation(a, b):
    def f(xs):
        return [0 if math.isnan(x[-1]) or math.isinf(x[-1]) else x[-1] for x in a]
    return pearsonr(f(a), f(b)).statistic


def correlate(quarry, config_sample, sample_values):
    default = quarry.get_heatmap_export(lambda s: None)
    corr = []
    for sample_value in sample_values:
        config_sample(sample_value)
        sample = quarry.get_heatmap_export(lambda s: None)
        corr.append(get_correlation(sample, default))
    return corr


def correlate_icing(quarry, sample_values):
    # use cooler implementation as default
    quarry.set_value(["settings", "normalization", "normalize_by"], "cool-ice")
    def conf(sample_value):
        quarry.set_value(["settings", "normalization", "normalize_by"], "ice")
        quarry.set_value(["settings", "normalization", "num_ice_bins", "val"], sample_value)
    return correlate(quarry, conf, sample_values)


def correlate_associated_slices(quarry, sample_values):
    quarry.set_value(["settings", "normalization", "normalize_by"], "grid-seq")
    # max out num samples for the default setting
    quarry.set_value(["settings", "normalization", "grid_seq_samples", "val"], quarry.get_value([ "contigs", "genome_size" ]))
    def conf(sample_value):
        quarry.set_value(["settings", "normalization", "grid_seq_samples", "val"], sample_value)
    return correlate(quarry, conf, sample_values)


def correlate_ddd(quarry, sample_values):
    quarry.set_value(["settings", "normalization", "ddd"], True)
    # max out num samples for the default setting
    quarry.set_value(["settings", "normalization", "ddd_samples", "val_max"], 
                     quarry.get_value([ "contigs", "genome_size" ]))
    def conf(sample_value):
        quarry.set_value(["settings", "normalization", "ddd_samples", "val_max"], sample_value)
    return correlate(quarry, conf, sample_values)

SAMPLE_VALUES = [2**n for n in range(1, 14)]

def correlate_all(quarry, sample_values=SAMPLE_VALUES):
    with open_default_json() as default_file:
        default_json = json.load(default_file)
    print("name", *sample_values, sep="\t")
    for name, func in [
        ("ICE", correlate_icing),
        ("associated slices", correlate_associated_slices),
        ("DDD", correlate_ddd),
    ]:
        quarry.set_value(["settings"], default_json)
        print(name, *func(quarry, sample_values), sep="\t")