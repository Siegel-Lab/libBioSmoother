from scipy.stats.stats import pearsonr
from .quarry import open_default_json
import json
import math
import sys


def get_correlation(a, b):
    da = {x[:-1]: x[-1] for x in a}
    c = []
    d = []
    for x in b:
        k = x[:-1]
        v = x[-1]
        if not (math.isnan(v) or math.isinf(v)):
            if k in da and not (math.isnan(da[k]) or math.isinf(da[k])):
                c.append(v)
                d.append(da[k])
    return pearsonr(c, d).statistic


# @todo figure out why ice does not correlate properly...
def correlate(quarry, config_sample, sample_values, window_size, n_windows=10):
    canvas_size_x, canvas_size_y = quarry.get_canvas_size(lambda s: None)
    quarry.set_value(["area", "x_start"], 0)
    quarry.set_value(["area", "x_end"], canvas_size_x)
    quarry.set_value(["area", "y_start"], 0)
    quarry.set_value(["area", "y_end"], canvas_size_y)
    default = quarry.get_heatmap_export(lambda s: None)
    for sample_value in sample_values:
        config_sample(sample_value)
        corrs_local = []
        for x_start in range(0, canvas_size_x, window_size)[:n_windows]:
            quarry.set_value(["area", "x_start"], x_start)
            quarry.set_value(["area", "x_end"], x_start + window_size)
            for y_start in range(0, canvas_size_y, window_size)[:n_windows]:
                quarry.set_value(["area", "y_start"], y_start)
                quarry.set_value(["area", "y_end"], y_start + window_size)
                sample = quarry.get_heatmap_export(lambda s: None)
                corrs_local.append(get_correlation(sample, default))
        print(sum(corrs_local) / len(corrs_local), end="\t", flush=True)


def correlate_icing(quarry, sample_values, window_size):
    # use cooler implementation as default
    quarry.set_value(["settings", "normalization", "normalize_by"], "cool-ice")

    def conf(sample_value):
        quarry.set_value(["settings", "normalization", "normalize_by"], "ice")
        quarry.set_value(
            ["settings", "normalization", "num_ice_bins", "val"], sample_value
        )

    correlate(quarry, conf, sample_values, window_size)


def correlate_associated_slices(quarry, sample_values, window_size):
    quarry.set_value(["settings", "normalization", "normalize_by"], "grid-seq")
    # max out num samples for the default setting
    quarry.set_value(
        ["settings", "normalization", "grid_seq_samples", "val"],
        quarry.get_value(["contigs", "genome_size"]),
    )

    def conf(sample_value):
        quarry.set_value(
            ["settings", "normalization", "grid_seq_samples", "val"], sample_value
        )

    correlate(quarry, conf, sample_values, window_size)


def correlate_ddd(quarry, sample_values, window_size):
    quarry.set_value(["settings", "normalization", "ddd"], True)
    # max out num samples for the default setting
    quarry.set_value(
        ["settings", "normalization", "ddd_samples", "val_max"],
        quarry.get_value(["contigs", "genome_size"]),
    )

    def conf(sample_value):
        quarry.set_value(
            ["settings", "normalization", "ddd_samples", "val_max"], sample_value
        )

    correlate(quarry, conf, sample_values, window_size)


def __adjust_by_dividend(quarry, v):
    div = quarry.get_value(["dividend"])
    if v % div != 0:
        print("WARNING: uneven division by index dividend", file=sys.stderr)
    if v < div:
        print("WARNING: dividend larger than value", file=sys.stderr)
    return max(1, v // div)


SAMPLE_VALUES = [2**n for n in range(1, 14)]
KBP = 1000
MBP = 1000 * KBP
WINDOW_SIZES = [500 * KBP, 1 * MBP, 5 * MBP]
BIN_SIZES = [10 * KBP, 50 * KBP, 100 * KBP]


def correlate_all(
    quarry, sample_values=SAMPLE_VALUES, window_sizes=WINDOW_SIZES, bin_sizes=BIN_SIZES
):
    with open_default_json() as default_file:
        default_json = json.load(default_file)
    print("name", "bin_size", "window_size", *sample_values, sep="\t")
    for name, func in [
        ("ICE", correlate_icing),
        # ("associated slices", correlate_associated_slices),
        # ("DDD", correlate_ddd),
    ]:
        quarry.set_value(["settings"], default_json)
        quarry.set_value(["settings", "interface", "fixed_bin_size"], True)
        quarry.set_value(["settings", "interface", "add_draw_area", "val"], 0)
        quarry.set_value(
            ["contigs", "displayed_on_x"], [quarry.get_value(["contigs", "list"])[0]]
        )
        quarry.set_value(
            ["contigs", "displayed_on_y"], [quarry.get_value(["contigs", "list"])[0]]
        )
        for bin_size in bin_sizes:
            b = __adjust_by_dividend(quarry, bin_size)
            quarry.set_value(["settings", "interface", "fixed_bin_size_x", "val"], b)
            quarry.set_value(["settings", "interface", "fixed_bin_size_y", "val"], b)
            for window_size in window_sizes:
                w = __adjust_by_dividend(quarry, window_size)

                print(name, bin_size, window_size, sep="\t", end="\t", flush=True)
                func(quarry, sample_values, w)
                print()
