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
from ._import_lib_cm import PartialQuarry, SPS_VERSION, CM_VERSION


class Quarry(PartialQuarry):
    def __init__(self, *args):
        PartialQuarry.__init__(self, *args)
        self.print_callback = lambda s: print(s)

    def normalizeBinominalTestTrampoline(
        self, bin_values, num_interactions_total, num_bins, p_accept
    ):
        def bin_test(jdx=0):
            ret = []
            for idx, val in enumerate(bin_values):
                ret.append(binom_test(val[jdx], num_interactions_total[idx % len(num_interactions_total)][jdx], 
                                      1 / num_bins, alternative="greater"))
            return ret
        psx = bin_test(0)
        psy = bin_test(1)
        return [
            (1 if x < p_accept else 0, 1 if y < p_accept else 0)
            for x, y in zip(multipletests(psx, alpha=float("NaN"), method="fdr_bh")[1],
                            multipletests(psy, alpha=float("NaN"), method="fdr_bh")[1])
        ]

    def normalizeCoolerTrampoline(
        self, bin_values, axis_size
    ):
        return icing(bin_values, axis_size)

    def __combine_hex_values(self, d):
        ## taken from: https://stackoverflow.com/questions/61488790/how-can-i-proportionally-mix-colors-in-python
        d_items = sorted(d.items())
        tot_weight = sum(d.values())
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

    def print(self, s):
        self.print_callback(s)

    @staticmethod
    def get_libSps_version():
        return SPS_VERSION

    @staticmethod
    def get_libCm_version():
        return CM_VERSION

    @staticmethod
    def has_cooler_icing():
        return HAS_COOLER_ICING
