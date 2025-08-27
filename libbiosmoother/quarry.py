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
import math


def open_default_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "default.json").open("r")


def open_descriptions_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "descriptions.json").open(
        "r"
    )


def open_button_names_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "button_names.json").open(
        "r"
    )


def open_button_tabs_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "button_tabs.json").open(
        "r"
    )


def open_valid_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "valid.json").open("r")


class Quarry(PartialQuarry):
    def establish_backwards_compatibility(self):
        with open_default_json() as default_settings_file:
            default_settings = json.load(default_settings_file)

            def combine_dict(a, b, previous_keys=[]):
                r = {}
                for k in b.keys():
                    if isinstance(b[k], dict) and k in a:
                        r[k] = combine_dict(
                            a[k], b[k], previous_keys=previous_keys + [k]
                        )
                    elif isinstance(b[k], dict):
                        print(
                            "WARNING: the loaded index was missing the",
                            ".".join(previous_keys + [k]),
                            "setting, copying it from the default settings.",
                        )
                        r[k] = b[k]
                    elif k in a:
                        r[k] = a[k]
                    else:
                        print(
                            "WARNING: the loaded index was missing the",
                            ".".join(previous_keys + [k]),
                            "setting, copying it from the default settings.",
                        )
                        r[k] = b[k]
                return r

            if not self.get_value(["settings"]) is None:
                self.set_value(
                    ["settings"],
                    combine_dict(self.get_value(["settings"]), default_settings),
                )

        if not self.get_value(["coverage", "list"]) is None:
            for name in self.get_value(["coverage", "list"]):
                for entry, val in [
                    ("no_map_q", False),
                    ("no_groups", False),
                    ("no_multi_map", False),
                    ("no_category", False),
                    ("no_strand", True),
                    ("shekelyan", False),
                ]:
                    if not self.has_value(["coverage", "by_name", name, entry]):
                        self.set_value(["coverage", "by_name", name, entry], val)
                        print(
                            "WARNING: the loaded index was missing the",
                            entry,
                            "entry for the",
                            name,
                            "secondary dataset, set value to",
                            val,
                            ".",
                        )

    def __init__(self, *args):
        PartialQuarry.__init__(self, *args)

        sps_in_index = self.get_value(["version", "lib_sps_version"])
        if sps_in_index != SPS_VERSION:
            print(
                "WARNING: the version of libSps that was used to create this index is different from the current version.",
                "This may lead to undefined behavior.\nVersion in index:",
                sps_in_index,
                "\nCurrent version:",
                SPS_VERSION,
                file=sys.stderr,
            )

        lib_bio_smoother_in_index = self.get_value(
            ["version", "lib_bio_smoother_version"]
        )
        if lib_bio_smoother_in_index != LIB_BIO_SMOOTHER_CPP_VERSION:
            print(
                "WARNING: the version of libBioSmoother that was used to create this index is different from the current version.",
                "This may lead to undefined behavior.\nVersion in index:",
                lib_bio_smoother_in_index,
                "\nCurrent version:",
                LIB_BIO_SMOOTHER_CPP_VERSION,
                file=sys.stderr,
            )

        self.establish_backwards_compatibility()

    def normalizeBinominalTestTrampoline(
        self,
        bin_values,
        num_interactions_total,
        num_bins_interacting_with,
        p_accept,
        is_col,
        grid_height,
    ):
        if grid_height == 0 or len(bin_values) == 0:
            return []

        def bin_test(jdx):
            ret = []
            for idx, val in enumerate(bin_values):
                n = num_interactions_total[
                    (idx // grid_height if is_col else idx % grid_height)
                ][jdx]
                i = num_bins_interacting_with[
                    (idx // grid_height if is_col else idx % grid_height)
                ][jdx]
                x = val[jdx]
                if i > 0 and HAS_STATS:
                    p = 1 / i
                    ret.append(binom_test(x, n, p, alternative="greater"))
                else:
                    ret.append(1)
            return ret

        def split_list(l):
            grid_width = len(l) // grid_height
            assert grid_height * grid_width == len(l)
            ret_l = [[] for _ in range(grid_width if is_col else grid_height)]
            for idx, val in enumerate(l):
                ret_l[idx // grid_height if is_col else idx % grid_height].append(val)
            assert len(l) == sum(len(x) for x in ret_l)
            return ret_l

        def combine_list(l):
            ret_l = []
            if is_col:
                for val in l:
                    ret_l.extend(val)
            else:
                for val in zip(*l):
                    ret_l.extend(val)
            assert len(ret_l) == sum(len(x) for x in l)
            return ret_l

        def p_val_correction(ll):
            return [
                multipletests(l, alpha=float("NaN"), method="fdr_bh")[1] for l in ll
            ]

        def binarization(ll):
            return [[1 if x < p_accept else 0 for x in l] for l in ll]

        ret_x = combine_list(binarization(p_val_correction(split_list(bin_test(0)))))
        ret_y = combine_list(binarization(p_val_correction(split_list(bin_test(1)))))
        if len(ret_x) == 0 or len(ret_y) == 0:
            return []
        return list(zip(ret_x, ret_y))

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

    def set_ploidy_itr(self, ploidy_iterator, report_error=lambda s: None):
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
                    report_error(
                        "ERROR: The target contig name" +
                        str(chr_to) +
                        "occurs multiple times in the input file. Hence, the given ploidy file is not valid and will be ignored."
                    )
                    return
                if chr_from not in self.get_value(["contigs", "ploidy_list"]):
                    report_error(
                        "WARNING: The source contig name" +
                        str(chr_from) +
                        "does not occur in the dataset. It will be ignored."
                    )
                    continue
                if chr_from in curr_ploidy_group:
                    report_error(
                        "WARNING: The source contig name" +
                        str(chr_from) +
                        "occurs multiple times in the same ploidy group. Is this really what you want?"
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

    def set_ploidy_list(self, ploidy_file, report_error=lambda s: None):
        with fileinput.input(ploidy_file) as file:
            self.set_ploidy_itr(file, report_error=report_error)
            self.set_value(
                ["settings", "normalization", "ploidy_last_uploaded_filename"],
                ploidy_file,
            )
        return self

    def __isfloat(self, num):
        try:
            float(num)
            return True
        except ValueError:
            return False

    def __interpret_number(self, s, bot=True, report_error=lambda s: None):
        s_in = s
        if s[-1:] == "b":
            s = s[:-1]
        elif s[-2:] == "bp":
            s = s[:-2]
        fac = 1
        if len(s) > 0 and s[-1] == "g":
            fac = 1000000000
            s = s[:-1]
        if len(s) > 0 and s[-1] == "m":
            fac = 1000000
            s = s[:-1]
        if len(s) > 0 and s[-1] == "k":
            fac = 1000
            s = s[:-1]
        if len(s) > 0 and s[-1] == "h":
            fac = 100
            s = s[:-1]
        if len(s) > 1 and s[-2:] == "da":
            fac = 10
            s = s[:-2]
        s = s.replace(",", "")
        if self.__isfloat(s):
            val = (float(s) * fac) / self.get_value(["dividend"])
            if bot:
                return int(math.floor(val))
            else:
                return int(math.ceil(val))
        report_error("Could not interpret '" + str(s_in) + "' as a locus")
        return None

    def interpret_position(
        self, s, on_x_axis=True, bot=True, report_error=lambda s: None
    ):
        if s.count(":") == 0 and s.count("+-") == 1:
            x, y = s.split("+-")
            c = self.__interpret_number(y, True, report_error=report_error)
            if not c is None and bot:
                c = -c
            a = self.interpret_name(x, on_x_axis, bot, lambda x: None, report_error)
            if not a is None and not c is None:
                return [a + c]
        elif s.count(":") == 1:
            x, y = s.split(":")
            if "+-" in y:
                y1, y2 = y.split("+-")
                if len(y1) == 0:
                    b = 0
                else:
                    b = self.__interpret_number(y1, bot, report_error=report_error)
                c = self.__interpret_number(y2, True, report_error=report_error)
                if not c is None and bot:
                    c = -c
                a = self.interpret_name(
                    x,
                    on_x_axis,
                    bot if len(y1) == 0 else True,
                    lambda x: None,
                    report_error,
                )
                if not a is None and not b is None and not c is None:
                    return [a + b + c]
            else:
                b = self.__interpret_number(y, bot, report_error=report_error)
                a = self.interpret_name(
                    x, on_x_axis, True, lambda x: None, report_error
                )
                if not a is None and not b is None:
                    return [a + b]

        # try to interpret as a number
        a = self.__interpret_number(s, bot, report_error=lambda s: None)
        if not a is None:
            return [a]

        if not s is None:
            # try to interpret as a contig name
            a = self.interpret_name(s, on_x_axis, bot, lambda x: None, lambda s: None)
            if not a is None:
                return [a]

        report_error("Could not interpret '" + str(s) + "' as 'contig:locus'")
        return [None]

    def interpret_range(self, s, on_x_axis=True, report_error=lambda s: None):
        s = "".join(s.lower().split())
        if (
            s.count("..") == 1
            and s.count(":") == 1
            and s.index(":") < s.index("..")
            and s.count("[") <= 1
            and s.count("]") <= 1
        ):
            s2, z = s.split("..")
            x, y = s2.split(":")
            if x[:1] == "[":
                x = x[1:]
            if z[-1:] == "]":
                z = z[:-1]

            a = self.interpret_name(x, on_x_axis, True, lambda x: None, report_error)
            b = self.__interpret_number(y, True, report_error=report_error)
            c = self.__interpret_number(z, False, report_error=report_error)
            ret = [None, None]
            if not a is None and not b is None:
                ret[0] = a + b
            if not a is None and not c is None:
                ret[1] = a + c
        elif s.count("..") == 1 and s.count("[") <= 1 and s.count("]") <= 1:
            x, y = s.split("..")
            if x[:1] == "[":
                x = x[1:]
            if y[-1:] == "]":
                y = y[:-1]

            ret = self.interpret_position(
                x, on_x_axis, True, report_error=report_error
            ) + self.interpret_position(y, on_x_axis, False, report_error=report_error)
        else:
            if s[:1] == "[":
                s = s[1:]
            if s[-1:] == "]":
                s = s[:-1]
            ret = self.interpret_position(
                s, on_x_axis, True, report_error=report_error
            ) + self.interpret_position(
                s, on_x_axis, False, report_error=lambda x: None
            )
        if not None in ret:
            ret.sort()
            if ret[1] <= ret[0]:
                ret[1] = ret[0] + 1
        if None in ret:
            report_error("Could not interpret '" + str(s) + "' as a range")
        return ret

    def interpret_area(
        self,
        s,
        default_start_x,
        default_start_y,
        default_end_x,
        default_end_y,
        report_error=lambda s: None,
    ):
        # remove all space-like characters
        s = "".join(s.lower().split())
        if s.count(";") == 1 and s.count("x=") == 0 and s.count("y=") == 0:
            x, y = s.split(";")
            return self.interpret_range(
                x, True, report_error=report_error
            ) + self.interpret_range(y, False, report_error=report_error)

        if s.count("x=") == 1 and s[:2] == "x=" and s.count("y=") == 0:
            s = s[2:]
            return self.interpret_range(s, True, report_error=report_error) + [
                default_start_y,
                default_end_y,
            ]

        if s.count("x=") == 0 and s.count("y=") == 1 and s[:2] == "y=":
            s = s[2:]
            return [
                default_start_x,
                default_end_x,
            ] + self.interpret_range(s, True, report_error=report_error)

        if s.count("x=") == 1 and s.count("y=") == 1:
            x_pos = s.find("x=")
            y_pos = s.find("y=")
            x = s[x_pos + 2 : y_pos] if x_pos < y_pos else s[x_pos + 2 :]
            y = s[y_pos + 2 : x_pos] if y_pos < x_pos else s[y_pos + 2 :]
            return self.interpret_range(
                x, True, report_error=report_error
            ) + self.interpret_range(y, False, report_error=report_error)

        return self.interpret_range(
            s, True, report_error=report_error
        ) + self.interpret_range(s, False, report_error=lambda s: None)

    def __to_readable_pos(self, x, genome_end, contig_names, contig_starts, lcs=0):
        if len(contig_names) == 0 or len(contig_starts) == 0:
            return "n/a"
        x = int(x)
        if x < 0:
            idx = 0
        elif x >= genome_end * self.get_value(["dividend"]):
            idx = len(contig_names) - 1
            x -= contig_starts[-1] * self.get_value(["dividend"])
        else:
            idx = 0
            for idx, (start, end) in enumerate(
                zip(contig_starts, contig_starts[1:] + [genome_end])
            ):
                if x >= start * self.get_value(
                    ["dividend"]
                ) and x < end * self.get_value(["dividend"]):
                    x -= start * self.get_value(["dividend"])
                    break

        if x == 0:
            label = "0 bp"
        elif x % 1000000 == 0:
            label = "{:,}".format(x // 1000000) + " Mbp"
        elif x % 1000 == 0:
            label = "{:,}".format(x // 1000) + " kbp"
        else:
            label = "{:,}".format(x) + " bp"

        if idx >= len(contig_names):
            return "n/a"

        if lcs != 0:
            n = contig_names[idx][:-lcs]
        else:
            n = contig_names[idx]
        return n + ": " + label

    def get_readable_range(
        self, start, end, of_x_axis=True, genomic_coords=False, print=lambda x: None
    ):
        lcs = self.get_longest_common_suffix(print)
        contig_names = self.get_annotation_list(of_x_axis, print)
        if genomic_coords:
            contig_starts = self.get_contig_start_list(of_x_axis, print)
        else:
            contig_starts = self.get_tick_list(of_x_axis, print)
        if len(contig_starts) > 0:
            return (
                self.__to_readable_pos(
                    start * int(self.get_value(["dividend"])),
                    contig_starts[-1],
                    contig_names,
                    contig_starts[:-1],
                    lcs,
                )
                + " .. "
                + self.__to_readable_pos(
                    end * int(self.get_value(["dividend"])),
                    contig_starts[-1],
                    contig_names,
                    contig_starts[:-1],
                    lcs,
                )
            )
        else:
            return "n/a"

    def get_readable_area(
        self,
        start_x=None,
        start_y=None,
        end_x=None,
        end_y=None,
        genomic_coords=False,
        print=lambda x: None,
    ):
        if start_x is None:
            start_x = self.get_value(["area", "x_start"])
        if start_y is None:
            start_y = self.get_value(["area", "y_start"])
        if end_x is None:
            end_x = self.get_value(["area", "x_end"])
        if end_y is None:
            end_y = self.get_value(["area", "y_end"])
        return (
            "X=["
            + self.get_readable_range(start_x, end_x, True, genomic_coords, print)
            + "] Y=["
            + self.get_readable_range(start_y, end_y, False, genomic_coords, print)
            + "]"
        )

    @staticmethod
    def get_libSps_version():
        return SPS_VERSION

    @staticmethod
    def get_libBioSmoother_version():
        return LIB_BIO_SMOOTHER_CPP_VERSION

    @staticmethod
    def has_cooler_icing():
        return HAS_COOLER_ICING
