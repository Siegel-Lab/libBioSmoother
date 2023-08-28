try:
    from .cooler_interface import icing

    HAS_COOLER_ICING = True
except ImportError:
    # Error handling
    HAS_COOLER_ICING = False
    pass
try:
    from bokeh.palettes import (
        Viridis256,
        Colorblind,
        Plasma256,
        Turbo256,
    )  # pyright: ignore missing import

    HAS_PALETTES = True
except ImportError:
    # Error handling
    HAS_PALETTES = False
    pass

try:
    from statsmodels.stats.multitest import (
        multipletests,
    )  # pyright: ignore missing import
    from scipy.stats import binom_test  # pyright: ignore missing import

    HAS_STATS = True
except ImportError:
    # Error handling
    HAS_STATS = False
    pass

from ._import_lib_bio_smoother_cpp import (
    PartialQuarry,
    SPS_VERSION,
    LIB_BIO_SMOOTHER_CPP_VERSION,
)

try:
    import importlib.resources as pkg_resources  # pyright: ignore missing import
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources  # pyright: ignore missing import
import json
import sys
import fileinput


def open_default_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "default.json").open("r")


class Quarry(PartialQuarry):
    def __init__(self, *args):
        PartialQuarry.__init__(self, *args)

        sps_in_index = self.get_value(["version", "lib_sps_version"])
        if sps_in_index != SPS_VERSION:
            print(
                "WARNING: the version of libSps that was used to create this index is different from the current version.",
                "This may lead to undefined behavior. Version in index:",
                sps_in_index,
                "current version:",
                SPS_VERSION,
                file=sys.stderr,
            )

        lib_bio_smoother_in_index = self.get_value(
            ["version", "lib_bio_smoother_version"]
        )
        if lib_bio_smoother_in_index != LIB_BIO_SMOOTHER_CPP_VERSION:
            print(
                "WARNING: the version of libBioSmoother that was used to create this index is different from the current version.",
                "This may lead to undefined behavior. Version in index:",
                lib_bio_smoother_in_index,
                "current version:",
                LIB_BIO_SMOOTHER_CPP_VERSION,
                file=sys.stderr,
            )

    def normalizeBinominalTestTrampoline(
        self,
        bin_values,
        num_interactions_total,
        num_bins_interacting_with,
        samples,
        max_num_interacting_with,
        p_accept,
        is_col,
        grid_height,
    ):
        fac = max_num_interacting_with / samples

        def bin_test(jdx):
            ret = []
            for idx, val in enumerate(bin_values):
                n = num_interactions_total[
                    (idx // grid_height if is_col else idx % grid_height)
                ][jdx]
                i = (
                    num_bins_interacting_with[
                        (idx // grid_height if is_col else idx % grid_height)
                    ][jdx]
                    * fac
                )
                if i > 0 and HAS_STATS:
                    p = 1 / (i + 1)
                    x = val[jdx]
                    ret.append(binom_test(x, max(x, n), p, alternative="greater"))
                else:
                    ret.append(1)
            return ret

        psx = bin_test(0)
        psy = bin_test(1)
        if len(psx) == 0 or len(psy) == 0:
            return []
        return [
            (1 if x < p_accept else 0, 1 if y < p_accept else 0)
            for x, y in zip(
                multipletests(psx, alpha=float("NaN"), method="fdr_bh")[1],
                multipletests(psy, alpha=float("NaN"), method="fdr_bh")[1],
            )
        ]

    def normalizeCoolerTrampoline(self, bin_values, axis_size):
        return icing(bin_values, axis_size)

    def __combine_hex_values(self, d):
        ## taken from: https://stackoverflow.com/questions/61488790/how-can-i-proportionally-mix-colors-in-python
        d_items = sorted(d.items())
        tot_weight = max(1, sum(d.values()))
        red = int(sum([int(k[:2], 16) * v for k, v in d_items]) / tot_weight)
        green = int(sum([int(k[2:4], 16) * v for k, v in d_items]) / tot_weight)
        blue = int(sum([int(k[4:6], 16) * v for k, v in d_items]) / tot_weight)
        zpad = lambda x: x if len(x) == 2 else "0" + x
        return "#" + zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])

    def colorPalette(self, palette_name, color_low, color_high):
        if HAS_PALETTES:
            if palette_name == "Viridis256":
                return Viridis256
            elif palette_name == "Plasma256":
                return Plasma256
            elif palette_name == "Turbo256":
                return Turbo256
            elif palette_name == "Fall":
                white = "ffffff"
                orange = "f5a623"
                red = "d0021b"
                black = "000000"
                return (
                    [
                        self.__combine_hex_values({white: 1 - x / 100, orange: x / 100})
                        for x in range(100)
                    ]
                    + [
                        self.__combine_hex_values({orange: 1 - x / 100, red: x / 100})
                        for x in range(100)
                    ]
                    + [
                        self.__combine_hex_values({red: 1 - x / 100, black: x / 100})
                        for x in range(100)
                    ]
                )
            elif palette_name == "LowToHigh":
                return [
                    self.__combine_hex_values(
                        {color_low[1:]: 1 - x / 255, color_high[1:]: x / 255}
                    )
                    for x in range(256)
                ]
            else:
                raise RuntimeError("invalid value for color_palette")
        else:
            if palette_name == "LowToHigh":
                return [
                    self.__combine_hex_values(
                        {color_low[1:]: 1 - x / 255, color_high[1:]: x / 255}
                    )
                    for x in range(256)
                ]
            else:
                return [
                    self.__combine_hex_values(
                        {"000000": 1 - x / 255, "ffffff": x / 255}
                    )
                    for x in range(256)
                ]

    def compute_biases(
        self, dataset_name, default_session, progress_print, ice_resolution=50000
    ):
        # set session as needed
        ## reset to default session
        self.set_session(default_session)
        ## set default settings
        with open_default_json() as f:
            self.set_value(["settings"], json.load(f))
        ## modify parameters as needed

        # pick the relevant dataset
        self.set_value(["replicates", "in_group_a"], [dataset_name])
        self.set_value(["replicates", "in_group_b"], [])

        # never skip a region -> if the user decides to display that region the biases might be missing
        self.set_value(["settings", "filters", "cut_off_bin"], "smaller")

        # activate the local balancing but render the whole heatmap
        self.set_value(["settings", "normalization", "normalize_by"], "ice-local")

        # render whole heatmap
        self.set_value(["settings", "export", "do_export_full"], True)

        # fix the bin size
        self.set_value(["settings", "interface", "fixed_bin_size"], True)
        div_resolution = max(1, ice_resolution // self.get_value(["dividend"]))
        self.set_value(
            ["settings", "interface", "fixed_bin_size_x", "val"], div_resolution
        )
        self.set_value(
            ["settings", "interface", "fixed_bin_size_y", "val"], div_resolution
        )

        biases_x = self.get_slice_bias(0, 0, 0, progress_print)
        biases_y = self.get_slice_bias(0, 0, 1, progress_print)

        coords_x = self.get_axis_coords(True, progress_print)
        coords_y = self.get_axis_coords(False, progress_print)

        return biases_x, coords_x, biases_y, coords_y

    def copy(self):
        # trigger the cpp copy constructor
        return Quarry(super(PartialQuarry, self))

    def set_ploidy_itr(self, ploidy_iterator):
        ploidy_map = {}
        ploidy_list = []
        ploidy_groups = {}
        curr_ploidy_group = set()
        group_count = 1
        for line in ploidy_iterator:
            line = line[:-1].strip()
            if len(line) > 0 and line[0] != "#":
                # if whole line is '-'
                if all(c == "-" for c in line):
                    group_count += 1
                    curr_ploidy_group = set()
                    continue
                chr_from, chr_to = line.split()
                if chr_to in ploidy_map:
                    print(
                        "ERROR: The target contig name",
                        chr_to,
                        "occurs multiple times in the input file. Hence, the given ploidy file is not valid and will be ignored.",
                    )
                    return
                if chr_from not in self.get_value(["contigs", "ploidy_list"]):
                    print(
                        "WARNING: The source contig name",
                        chr_from,
                        "does not occur in the dataset. It will be ignored.",
                    )
                    continue
                if chr_from in curr_ploidy_group:
                    print(
                        "WARNING: The source contig name",
                        chr_from,
                        "occurs multiple times in the same ploidy group. Is this really what you want?",
                    )
                    continue
                ploidy_map[chr_to] = chr_from
                ploidy_list.append(chr_to)
                ploidy_groups[chr_to] = group_count
                curr_ploidy_group.add(chr_from)
        self.set_value(["contigs", "list"], ploidy_list)
        self.set_value(["contigs", "displayed_on_x"], ploidy_list)
        self.set_value(["contigs", "displayed_on_y"], ploidy_list)
        self.set_value(["contigs", "ploidy_map"], ploidy_map)
        self.set_value(["contigs", "ploidy_groups"], ploidy_groups)
        self.save_session()

    def set_ploidy_list(self, ploidy_file):
        with fileinput.input(ploidy_file) as file:
            self.set_ploidy_itr(file)

    @staticmethod
    def get_libSps_version():
        return SPS_VERSION

    @staticmethod
    def get_libBioSmoother_version():
        return LIB_BIO_SMOOTHER_CPP_VERSION

    @staticmethod
    def has_cooler_icing():
        return HAS_COOLER_ICING
