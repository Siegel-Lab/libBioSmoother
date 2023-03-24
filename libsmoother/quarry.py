try:
    from .cooler_interface import icing

    HAS_COOLER_ICING = True
except ImportError:
    # Error handling
    HAS_COOLER_ICING = False
    pass
try:
    from bokeh.palettes import Viridis256, Colorblind, Plasma256, Turbo256

    HAS_PALETTES = True
except ImportError:
    # Error handling
    HAS_PALETTES = False
    pass

try:
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import binom_test

    HAS_STATS = True
except ImportError:
    # Error handling
    HAS_STATS = False
    pass

from ._import_lib_smoother_cpp import (
    PartialQuarry,
    SPS_VERSION,
    LIB_SMOOTHER_CPP_VERSION,
)

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources
import json


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
            )

        lib_smoother_in_index = self.get_value(["version", "lib_smoother_version"])
        if lib_smoother_in_index != LIB_SMOOTHER_CPP_VERSION:
            print(
                "WARNING: the version of libSmoother that was used to create this index is different from the current version.",
                "This may lead to undefined behavior. Version in index:",
                lib_smoother_in_index,
                "current version:",
                LIB_SMOOTHER_CPP_VERSION,
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
                if i > 0 and HAS_PALETTES:
                    p = 1 / (i + 1)
                    x = val[jdx]
                    ret.append(binom_test(x, max(x, n), p, alternative="greater"))
                else:
                    ret.append(1)
            return ret

        psx = bin_test(0)
        psy = bin_test(1)
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

    def compute_biases(self, dataset_name, default_session, progress_print):
        # set session as needed
        ## reset to default session
        self.set_session(default_session)
        ## set default settings
        with (pkg_resources.files("libsmoother") / "default.json").open("r") as f:
            self.set_value(["settings"], json.load(f))
        ## modify parameters as needed
        self.set_value(["replicates", "in_group_a"], [dataset_name])
        self.set_value(["replicates", "in_group_b"], [])

        self.set_value(["settings", "export", "do_export_full"], True)
        self.set_value(["settings", "interface", "fixed_number_of_bins"], True)
        genome_size = self.get_value([ "contigs", "genome_size" ])
        self.set_value(["settings", "interface", "fixed_num_bins_x", "val"], genome_size)
        self.set_value(["settings", "interface", "fixed_num_bins_y", "val"], genome_size)

        biases_x = self.get_slice_bias(0, 0, 0, progress_print)
        biases_y = self.get_slice_bias(0, 0, 1, progress_print)
        
        coords_x = self.get_axis_coords(True, progress_print)
        coords_y = self.get_axis_coords(False, progress_print)

        return biases_x, coords_x, biases_y, coords_y

    def copy(self):
        # trigger the cpp copy constructor
        return Quarry(super(PartialQuarry, self))

    @staticmethod
    def get_libSps_version():
        return SPS_VERSION

    @staticmethod
    def get_libSmoother_version():
        return LIB_SMOOTHER_CPP_VERSION

    @staticmethod
    def has_cooler_icing():
        return HAS_COOLER_ICING
