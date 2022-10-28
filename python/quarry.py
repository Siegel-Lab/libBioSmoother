from bokeh.palettes import Viridis256, Colorblind, Plasma256, Turbo256
from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test
from .libContactMapping import PartialQuarry, SPS_VERSION


class Quarry(PartialQuarry):
    def normalizeBinominalTestTrampoline(
        self, bin_values, num_interactions_total, num_bins, p_accept
    ):
        ps = [
            binom_test(x, num_interactions_total, 1 / num_bins, alternative="greater")
            for x in bin_values
        ]
        return [
            1 if x < p_accept else 0
            for x in multipletests(ps, alpha=float("NaN"), method="fdr_bh")[1]
        ]

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

    @staticmethod
    def get_libSps_version():
        return SPS_VERSION
