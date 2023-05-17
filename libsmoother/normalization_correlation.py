from .quarry import open_default_json
import json
import sys

"""
def correlate_icing(quarry, sample_values, window_size):
    # use cooler implementation as default
    quarry.set_value(["settings", "normalization", "normalize_by"], "cool-ice")

    def conf(sample_value):
        quarry.set_value(["settings", "normalization", "normalize_by"], "ice")
        quarry.set_value(
            ["settings", "normalization", "num_ice_bins", "val"], sample_value
        )
    def config_neg_ctr():
        quarry.set_value(["settings", "normalization", "normalize_by"], "dont")
    def config_pos_ctr():
        quarry.set_value(["settings", "normalization", "normalize_by"], "dont")
        #quarry.set_value(["settings", "normalization", "normalize_by"], "ice")
        #quarry.set_value(["settings", "normalization", "num_ice_bins", "val"], 50)

    return correlate(quarry, config_neg_ctr, config_pos_ctr, conf, sample_values, window_size)


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


SAMPLE_VALUES = [*range(0, 76, 5)]#[2**n for n in range(4, 14)]
KBP = 1000
MBP = 1000 * KBP
WINDOW_SIZES = [500 * KBP,
                1 * MBP,
                5 * MBP]
BIN_SIZES = [10 * KBP,
             50 * KBP,
             100 * KBP]
"""


def conf_quarry(quarry):
    #warnings.filterwarnings('ignore')
    with open_default_json() as default_file:
        default_json = json.load(default_file)
        quarry.set_value(["settings"], default_json)
    quarry.set_value(["settings", "filters", "symmetry"], "mirror")
    quarry.set_value(["settings", "filters", "cut_off_bin"], "smaller")
    quarry.set_value(["settings", "filters", "show_contig_smaller_than_bin"], True)
    quarry.set_value(["settings", "interface", "fixed_bin_size"], True)
    quarry.set_value(["settings", "interface", "add_draw_area", "val"], 0)
    quarry.set_value(["settings", "normalization", "scale"], "dont")
    quarry.set_value(["settings", "normalization", "log_base"], 0)

def set_bin_size(quarry, bin_size):
    div = quarry.get_value(["dividend"])
    if bin_size % div != 0:
        print("WARNING: uneven division by index dividend", file=sys.stderr)
    if bin_size < div:
        print("WARNING: dividend larger than value", file=sys.stderr)
    bin_size = max(1, bin_size // div)
    quarry.set_value(["settings", "interface", "fixed_bin_size_x", "val"], bin_size)
    quarry.set_value(["settings", "interface", "fixed_bin_size_y", "val"], bin_size)

def __tsv_to_ret(tsv):
    ret = [(x[:-1], x[-1]) for x in tsv]
    ret.sort()
    return [a for a, _ in ret], [b for _, b in ret]

def quarry_whole_window(quarry, window_size, n_windows):
    canvas_size_x, canvas_size_y = quarry.get_canvas_size(lambda s: None)
    quarry.set_value(["area", "x_start"], 0)
    quarry.set_value(["area", "x_end"], min(window_size * n_windows, canvas_size_x))
    quarry.set_value(["area", "y_start"], 0)
    quarry.set_value(["area", "y_end"], min(window_size * n_windows, canvas_size_y))
    return __tsv_to_ret(quarry.get_heatmap_export(lambda s: None))

def quarry_chunked_window(quarry, window_size, n_windows):
    canvas_size_x, canvas_size_y = quarry.get_canvas_size(lambda s: None)
    tsv = []
    for x_start in range(0, min(window_size * n_windows, canvas_size_x), window_size):
        quarry.set_value(["area", "x_start"], x_start)
        quarry.set_value(["area", "x_end"], x_start + window_size)
        for y_start in range(0, min(window_size * n_windows, canvas_size_y), window_size):
            quarry.set_value(["area", "y_start"], y_start)
            quarry.set_value(["area", "y_end"], y_start + window_size)
            tsv.extend(quarry.get_heatmap_export(lambda s: None))
    return __tsv_to_ret(tsv)

