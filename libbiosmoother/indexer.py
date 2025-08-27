from ._import_lib_bio_smoother_cpp import (
    Index,
    SPS_VERSION,
    LIB_BIO_SMOOTHER_CPP_VERSION,
    SPS_BUILD_TIME,
    LIB_BIO_SMOOTHER_CPP_BUILD_TIME,
    MAX_NUM_FILTER_ANNOTATIONS,
)
from ._parse_and_group_reads import *
import json
import os
import copy
import time
from importlib.metadata import version
from .quarry import Quarry
from .quarry import open_default_json
import sys

try:
    from .cooler_interface import CoolerIterator

    HAS_COOLER_ITERATOR = True
except ImportError:
    # Error handling
    HAS_COOLER_ITERATOR = False
    pass

MAP_Q_MAX = 255


def touch(f_name):
    with open(f_name, "a"):  # Create file if does not exist
        pass


GENERATE_VERBOSITY = 0
PROGRESS_PRINT_TIME = 3


class Indexer:
    def __init__(self, prefix, strict=False):
        self.last_prog_print = None
        if strict:
            self.prefix = prefix + ".smoother_index"
        else:
            self.prefix = None
            for possible in [prefix, prefix + ".smoother_index"]:
                if os.path.exists(possible) and os.path.isdir(possible):
                    self.prefix = possible
                    break
        if self.prefix is None:
            raise RuntimeError("the given index", prefix, "does not exist.")
        self.progress_print("working on index", self.prefix, force_print=True)
        self.indices = None
        self.session_default = {}
        self.session = {}
        self.try_load_index()

    def try_load_index(self):
        if os.path.exists(self.prefix + "/default_session.json"):
            # self.indices = CachedSpsInterface(prefix + ".smoother_index")
            self.indices = Index(self.prefix, True)
            with open(self.prefix + "/default_session.json", "r") as f:
                self.session_default = json.load(f)
            with open(self.prefix + "/session.json", "r") as f:
                self.session = json.load(f)

    def progress_print(self, *text, force_print=False):
        t = time.time()
        if (
            self.last_prog_print is None
            or t - self.last_prog_print >= PROGRESS_PRINT_TIME
            or force_print
        ):
            print(*text)
            self.last_prog_print = t

    def save_session(self):
        with open(self.prefix + "/default_session.json", "w") as f:
            json.dump(self.session_default, f)
        with open(self.prefix + "/session.json", "w") as f:
            json.dump(self.session, f)
        del self.indices

    def set_session(self, keys, val):
        curr = self.session
        curr_def = self.session_default
        for key in keys[:-1]:
            curr = curr[key]
            curr_def = curr_def[key]
        # do not put the same object in both jsons -> it might have references to the same subobject
        # instead deep-copy the object
        curr[keys[-1]] = copy.deepcopy(val)
        curr_def[keys[-1]] = copy.deepcopy(val)

    def append_session(self, keys, val):
        curr = self.session
        curr_def = self.session_default
        for key in keys[:-1]:
            curr = curr[key]
            curr_def = curr_def[key]
        # do not put the same object in both jsons -> it might have references to the same subobject
        # instead deep-copy the object
        curr[keys[-1]].append(copy.deepcopy(val))
        curr_def[keys[-1]].append(copy.deepcopy(val))

    def create_session(
        self,
        chr_len_file_name,
        dividend,
        anno_path,
        test=False,
        filterable_annotations=[],
        map_q_thresholds=[],
    ):
        if os.path.exists(self.prefix):
            print("ERROR: The given index already exists.")
            exit()
        os.makedirs(self.prefix)
        self.set_session(
            ["version"],
            {
                "lib_sps_version": SPS_VERSION,
                "lib_sps_build_time": SPS_BUILD_TIME,
                "pip_lib_bio_smoother_version": version("libbiosmoother"),
                "lib_bio_smoother_version": LIB_BIO_SMOOTHER_CPP_VERSION,
                "lib_bio_smoother_build_time": LIB_BIO_SMOOTHER_CPP_BUILD_TIME,
            },
        )
        self.set_session(["dividend"], dividend)
        self.set_session(["previous"], None)
        self.set_session(["next"], None)
        self.set_session(["settings"], None)
        self.set_session(
            ["replicates"],
            {
                "list": [],
                "by_name": {},
                "in_group_a": [],
                "in_group_b": [],
                "in_column": [],
                "in_row": [],
                "cov_column_a": [],
                "cov_column_b": [],
                "cov_row_a": [],
                "cov_row_b": [],
            },
        )
        self.set_session(
            ["coverage"],
            {
                "list": [],
                "by_name": {},
                "in_column": [],
                "in_row": [],
                "cov_column_a": [],
                "cov_column_b": [],
                "cov_row_a": [],
                "cov_row_b": [],
            },
        )
        self.set_session(
            ["annotation"],
            {
                "list": [],
                "filterable": [],
                "by_name": {},
                "visible_x": [],
                "visible_y": [],
                "filter_present_x": [],
                "filter_present_y": [],
                "filter_absent_x": [],
                "filter_absent_y": [],
            },
        )
        self.set_session(
            ["contigs"],
            {
                "list": [],
                "ploidy_map": {},
                "ploidy_list": [],
                "ploidy_groups": {},
                "lengths": {},
                "displayed_on_x": [],
                "displayed_on_y": [],
                "displayed_on_x_ploidy": [],
                "displayed_on_y_ploidy": [],
                "genome_size": 0,
                "annotation_coordinates": "",
            },
        )
        self.set_session(["map_q_thresholds"], map_q_thresholds)
        if test:
            self.set_session(["test"], True)

        chr_len_list = []
        if chr_len_file_name != "":
            with open(chr_len_file_name, "r") as len_file:
                for line in len_file:
                    chr_name, chr_len = line.split()
                    chr_len_list.append((chr_name, chr_len))
        else:
            with open(anno_path, "r") as len_file:
                for line in len_file:
                    if line.startswith("##sequence-region "):
                        chr_name, chr_start, chr_end = line.split()[1:]
                        if chr_start != "1":
                            raise RuntimeError(
                                "Expecting sequence-regions to start at position 1, however this region does not: "
                                + line[:-1]
                            )
                        chr_len_list.append((chr_name, chr_end))
            if len(chr_len_list) == 0:
                raise RuntimeError(
                    "You are using a gff file without '##sequence-region' lines. You can use this file to load annotations but you will also have to supply a contig.sizes file using the --ctg_len parameter."
                )

        for chr_name, chr_len in chr_len_list:
            if not test or ("Chr1" in chr_name):
                self.append_session(["contigs", "list"], chr_name)
                self.append_session(["contigs", "ploidy_list"], chr_name)
                self.set_session(["contigs", "ploidy_map", chr_name], chr_name)
                self.set_session(["contigs", "ploidy_groups", chr_name], 1)
                self.set_session(
                    ["contigs", "lengths", chr_name], int(chr_len) // dividend
                )
                self.append_session(["contigs", "displayed_on_x"], chr_name)
                self.append_session(["contigs", "displayed_on_y"], chr_name)
                self.append_session(["contigs", "displayed_on_x_ploidy"], chr_name)
                self.append_session(["contigs", "displayed_on_y_ploidy"], chr_name)
                self.set_session(
                    ["contigs", "genome_size"],
                    self.session_default["contigs"]["genome_size"] + int(chr_len),
                )
        self.set_session(
            ["area"],
            {
                "x_start": 0,
                "y_start": 0,
                "x_end": self.session_default["contigs"]["genome_size"] // dividend,
                "y_end": self.session_default["contigs"]["genome_size"] // dividend,
            },
        )

        touch(self.prefix + "/sps.coords")
        touch(self.prefix + "/sps.datasets")
        touch(self.prefix + "/sps.overlays")
        touch(self.prefix + "/sps.prefix_sums")

        self.save_session()
        self.try_load_index()

        sorted_list = {}
        if anno_path != "":
            for name, chrom, start, end, info, on_forw_strnd in parse_annotations(
                anno_path
            ):
                if not chrom in self.session_default["contigs"]["list"]:
                    continue
                if name not in sorted_list:
                    sorted_list[name] = {}
                if chrom not in sorted_list[name]:
                    sorted_list[name][chrom] = []
                sorted_list[name][chrom].append((start, end, info, on_forw_strnd))
            existing_filterable_annotations = []
            for anno in filterable_annotations:
                if anno in sorted_list:
                    existing_filterable_annotations.append(anno)
            if len(existing_filterable_annotations) > MAX_NUM_FILTER_ANNOTATIONS:
                print(
                    "ERROR: The number of filterable annotations is limited to",
                    MAX_NUM_FILTER_ANNOTATIONS,
                    "but",
                    len(existing_filterable_annotations),
                    "were found.",
                )
                exit()
            self.set_session(
                ["annotation", "filterable"], existing_filterable_annotations
            )
            if len(existing_filterable_annotations) > 0:
                self.set_session(
                    ["contigs", "annotation_coordinates"],
                    existing_filterable_annotations[0],
                )
            else:
                self.set_session(["contigs", "annotation_coordinates"], "")
            for name, anno_by_chrom in sorted_list.items():
                if name not in self.session_default["annotation"]["list"]:
                    self.append_session(["annotation", "list"], name)
                    self.append_session(["annotation", "visible_x"], name)
                    self.append_session(["annotation", "visible_y"], name)

                    first_id = None
                    for chrom in self.session_default["contigs"]["list"]:
                        if chrom in anno_by_chrom:
                            annos = anno_by_chrom[chrom]
                        else:
                            annos = []
                        self.progress_print(
                            "annotating", name + "(s)", "for contig", chrom
                        )
                        idx = self.indices.anno.add_intervals(
                            annos,
                            self.session_default["dividend"],
                            verbosity=GENERATE_VERBOSITY,
                        )
                        if first_id is None:
                            first_id = idx
                    self.set_session(["annotation", "by_name", name], first_id)
                else:
                    raise RuntimeError("annotation with this name already exists")
        else:
            self.set_session(["contigs", "annotation_coordinates"], "")

        with open_default_json() as default_file:
            default_json = json.load(default_file)
        self.session["settings"] = default_json

        self.save_session()
        self.progress_print("done initializing.", force_print=True)

    def name_unique(self, name):
        return (
            name not in self.session_default["replicates"]["list"]
            and name not in self.session_default["coverage"]["list"]
        )

    def get_map_q_thresholds(self):
        t_dict = {0: 0}
        t_last = 0
        for idx, t in enumerate(self.session_default["map_q_thresholds"]):
            for i in range(t_last, t):
                t_dict[i + 1] = idx + 1
        for i in range(t, MAP_Q_MAX + 1):
            t_dict[i + 1] = idx + 2
        return t_dict

    def add_replicate(
        self,
        path,
        name,
        group="a",
        no_groups=False,
        keep_points=False,
        only_points=False,
        no_map_q=False,
        no_multi_map=False,
        no_category=False,
        no_strand=False,
        shekelyan=False,
        force_upper_triangle=False,
        columns=["chr1", "pos1", "chr2", "pos2"],
        allow_col_change=False,
    ):
        if not self.name_unique(name):
            raise RuntimeError(
                "The dataset name you provide must be unique but is not. "
                + "Use the <list> command to see all datasets."
            )

        self.progress_print("generating replicate...", force_print=True)

        self.append_session(["replicates", "list"], name)
        self.set_session(
            ["replicates", "by_name", name],
            {"ids": {}, "path": path, "ice_col": {}, "ice_row": {}},
        )
        if group in ["a", "both"]:
            self.append_session(["replicates", "in_group_a"], name)
        if group in ["b", "both"]:
            self.append_session(["replicates", "in_group_b"], name)

        contigs = self.session_default["contigs"]["list"]

        read_iterator = chr_order_heatmap(
            self.prefix,
            name,
            path,
            0,  # unused
            contigs,
            no_groups,
            "test" in self.session_default,
            force_upper_triangle,
            lambda *x: self.progress_print("loading", *x),
            columns,
            allow_col_change,
        )
        t_dict = self.get_map_q_thresholds()
        total_reads = 0

        num_itr = len(contigs) * len(contigs)
        cnt = 0
        first_id = None
        count_matrix_warning_done = False
        has_upper_triangle = False
        has_lower_triangle = False
        map_q_min = 255
        map_q_max = 0
        has_multimapper = False
        categoires_min = 2 ** (MAX_NUM_FILTER_ANNOTATIONS * 2 + 1)
        categories_max = 0
        strand_min = 4
        strand_max = 0
        for idx_x, chr_x in enumerate(contigs):
            anno_ids_x = [
                self.session_default["annotation"]["by_name"][anno] + idx_x
                for anno in self.session_default["annotation"]["filterable"]
            ]
            for idx_y, chr_y in enumerate(contigs):
                anno_ids_y = [
                    self.session_default["annotation"]["by_name"][anno] + idx_y
                    for anno in self.session_default["annotation"]["filterable"]
                ]
                cnt += 1
                self.progress_print(
                    "generating heatmap for contig-pair",
                    chr_x,
                    chr_y + ".",
                    cnt,
                    "of",
                    str(num_itr) + ":",
                    str(round(100 * cnt / num_itr, 2)) + "%",
                )
                for (
                    read_name,
                    pos_1_s,
                    pos_1_e,
                    pos_2_s,
                    pos_2_e,
                    pos_1_l,
                    pos_2_l,
                    strand_1,
                    strand_2,
                    map_q,
                    bin_cnt,
                ) in read_iterator.itr_cell(chr_x, chr_y):
                    if (
                        int(bin_cnt) > 1
                        and any(
                            int(p) % self.session_default["dividend"] != 0
                            for p in [pos_1_s, pos_1_e, pos_2_s, pos_2_e]
                        )
                        and not count_matrix_warning_done
                    ):
                        print(
                            "WARNING: The input file has a count column (i.e could be a count matrix)",
                            "and the bin position of at least one row does not match the minimal resolution",
                            "of the index.",
                            "pos1: ",
                            pos_1_s,
                            "..",
                            pos_1_e,
                            "pos2:",
                            pos_2_s,
                            "..",
                            pos_2_e,
                            "count:",
                            bin_cnt,
                            "minimal index resolution:",
                            str(self.session_default["dividend"])
                            + ". Will not show this warning again.",
                            file=sys.stderr,
                        )
                        count_matrix_warning_done = True
                    if idx_x > idx_y or (idx_x == idx_y and pos_1_s > pos_2_s):
                        has_lower_triangle = True
                    if idx_x < idx_y or (idx_x == idx_y and pos_1_s < pos_2_s):
                        has_upper_triangle = True
                    total_reads += 1
                    if no_category:
                        cat_x = [0] * MAX_NUM_FILTER_ANNOTATIONS
                        cat_y = cat_x
                    else:
                        cat_x = [
                            0 if x else 1
                            for x in self.indices.anno.get_categories(
                                [int(x) for x in pos_1_l.split(",")],
                                self.session_default["dividend"],
                                anno_ids_x,
                            )
                        ] + [0] * (MAX_NUM_FILTER_ANNOTATIONS - len(anno_ids_x))
                        cat_y = [
                            0 if x else 1
                            for x in self.indices.anno.get_categories(
                                [int(x) for x in pos_2_l.split(",")],
                                self.session_default["dividend"],
                                anno_ids_y,
                            )
                        ] + [0] * (MAX_NUM_FILTER_ANNOTATIONS - len(anno_ids_y))
                    act_pos_1_s = int(pos_2_s) // self.session_default["dividend"]
                    act_pos_2_s = int(pos_1_s) // self.session_default["dividend"]
                    if no_multi_map:
                        act_pos_1_e = act_pos_1_s
                        act_pos_2_e = act_pos_2_s
                    else:
                        act_pos_1_e = int(pos_2_e) // self.session_default["dividend"]
                        act_pos_2_e = int(pos_1_e) // self.session_default["dividend"]
                    if act_pos_1_e != act_pos_1_s or act_pos_2_e != act_pos_2_s:
                        has_multimapper = True

                    map_q = t_dict[int(map_q)]
                    if no_map_q:
                        map_q = 1
                    map_q_min = min(map_q_min, map_q)
                    map_q_max = max(map_q_max, map_q)

                    if no_strand:
                        same_strand_idx = 0
                        y_strand_idx = 0
                    else:
                        same_strand_idx = 0 if bool(strand_1) == bool(strand_2) else 1
                        y_strand_idx = 0 if bool(strand_1) else 1
                    strand_min = min(strand_min, same_strand_idx + y_strand_idx * 2)
                    strand_max = max(strand_max, same_strand_idx + y_strand_idx * 2)

                    bin_cnt = int(bin_cnt)

                    cat_pos = [
                        item for sublist in zip(cat_x, cat_y) for item in sublist
                    ]
                    cat_hash = sum(x * 2**i for i, x in enumerate(cat_pos))
                    categoires_min = min(categoires_min, cat_hash)
                    categories_max = max(categories_max, cat_hash)

                    start = [
                        act_pos_1_s,
                        act_pos_2_s,
                        MAP_Q_MAX - map_q - 1,
                        *cat_pos,
                        same_strand_idx,
                        y_strand_idx,
                    ]
                    end = [
                        act_pos_1_e,
                        act_pos_2_e,
                        MAP_Q_MAX - map_q - 1,
                        *cat_pos,
                        same_strand_idx,
                        y_strand_idx,
                    ]
                    self.indices.insert(start, end, bin_cnt)
                id = self.indices.generate(
                    fac=-2 if shekelyan else -1, verbosity=GENERATE_VERBOSITY
                )
                if first_id is None:
                    first_id = id

        if has_lower_triangle != has_upper_triangle and first_id == 0:
            print(
                "Info: Detected an exclusive",
                "lower" if has_lower_triangle else "upper",
                "triangle matrix. Setting symmetry and axis labels of index.",
            )
            self.session["settings"]["interface"]["axis_labels"] = "DNA_DNA"
            self.session["settings"]["filters"]["symmetry"] = "mirror"
            # @todo this should update default.json as well

        if total_reads == 0:
            print(
                "WARNING: the total number of reads that were added to the index is zero! Something seems off..."
            )
        if not no_map_q and map_q_min == map_q_max and total_reads > 1:
            print(
                "WARNING: Detected only a single mapping quality value (",
                map_q_min,
                "). Consider using --no_map_q if this dataset should not be filterable by mapping quality.",
            )
        if not has_multimapper and not no_multi_map and total_reads > 1:
            print(
                "WARNING: Detected no multimappers. Consider using --no_multi_map if this dataset should not be filterable by multimapping.",
            )
        if categories_max == categoires_min and not no_category and total_reads > 1:
            print(
                "WARNING: All interactions had the same annotation overlap (i.e. all interactions are overlapping some gene OR all interactions are overlapping nothing. NOT: some interactions are overlapping a gene and some others are overlapping nothing). Consider using --no_anno if this dataset should not be filterable by annotation overlap.",
            )
        if strand_min == strand_max and not no_strand and total_reads > 1:
            print(
                "WARNING: Detected only a single strand overlap type. Consider using --no_strand if this dataset should not be filterable by strand overlap.",
            )

        self.set_session(["replicates", "by_name", name, "first_dataset_id"], first_id)
        self.set_session(["replicates", "by_name", name, "num_datasets"], num_itr)
        self.set_session(["replicates", "by_name", name, "total_reads"], total_reads)
        self.set_session(["replicates", "by_name", name, "no_map_q"], no_map_q)
        self.set_session(["replicates", "by_name", name, "no_groups"], no_groups)
        self.set_session(["replicates", "by_name", name, "no_multi_map"], no_multi_map)
        self.set_session(["replicates", "by_name", name, "no_category"], no_category)
        self.set_session(["replicates", "by_name", name, "no_strand"], no_strand)
        self.set_session(["replicates", "by_name", name, "shekelyan"], shekelyan)
        self.set_session(
            ["replicates", "by_name", name, "force_upper_triangle"],
            force_upper_triangle,
        )
        self.set_session(
            ["replicates", "by_name", name, "has_upper_triangle"], has_upper_triangle
        )
        self.set_session(
            ["replicates", "by_name", name, "has_lower_triangle"], has_lower_triangle
        )

        read_iterator.cleanup()

        if not keep_points:
            self.indices.clear_points_and_desc()

        self.save_session()

        self.progress_print("done generating dataset.", force_print=True)

    def add_cool(
        self,
        path,
        name,
        group="a",
        no_category=False,
        shekelyan=False,
        force_upper_triangle=False,
    ):
        if not self.name_unique(name):
            raise RuntimeError(
                "The dataset name you provide must be unique but is not. "
                + "Use the <list> command to see all datasets."
            )

        self.progress_print(
            "generating replicate from cool (Note: this data will not be filterable on-the-fly as the required information for these filters is absent from cool files.)...",
            force_print=True,
        )

        self.append_session(["replicates", "list"], name)
        self.set_session(
            ["replicates", "by_name", name],
            {"ids": {}, "path": path, "ice_col": {}, "ice_row": {}},
        )
        if group in ["a", "both"]:
            self.append_session(["replicates", "in_group_a"], name)
        if group in ["b", "both"]:
            self.append_session(["replicates", "in_group_b"], name)

        contigs = self.session_default["contigs"]["list"]

        cool_iterator = CoolerIterator(path, self.session_default["dividend"])
        t_dict = self.get_map_q_thresholds()
        total_reads = 0

        num_itr = len(contigs) * len(contigs)
        cnt = 0
        first_id = None
        has_upper_triangle = False
        has_lower_triangle = False
        categoires_min = 2 ** (MAX_NUM_FILTER_ANNOTATIONS * 2 + 1)
        categories_max = 0
        for idx_x, chr_x in enumerate(contigs):
            anno_ids_x = [
                self.session_default["annotation"]["by_name"][anno] + idx_x
                for anno in self.session_default["annotation"]["filterable"]
            ]
            for idx_y, chr_y in enumerate(contigs):
                anno_ids_y = [
                    self.session_default["annotation"]["by_name"][anno] + idx_y
                    for anno in self.session_default["annotation"]["filterable"]
                ]
                cnt += 1
                self.progress_print(
                    "generating heatmap for contig-pair",
                    chr_x,
                    chr_y + ".",
                    cnt,
                    "of",
                    str(num_itr) + ":",
                    str(round(100 * cnt / num_itr, 2)) + "%",
                )
                for (
                    _,
                    pos_2,
                    _,
                    pos_1,
                    bin_cnt,
                ) in cool_iterator.iterate(chr_x, chr_y):
                    if idx_x > idx_y or (idx_x == idx_y and pos_1 > pos_2):
                        has_lower_triangle = True
                    if idx_x < idx_y or (idx_x == idx_y and pos_1 < pos_2):
                        has_upper_triangle = True
                    total_reads += 1
                    if no_category:
                        cat_x = [0] * MAX_NUM_FILTER_ANNOTATIONS
                        cat_y = cat_x
                    else:
                        cat_x = [
                            0 if x else 1
                            for x in self.indices.anno.get_categories(
                                [pos_1],
                                self.session_default["dividend"],
                                anno_ids_x,
                            )
                        ] + [0] * (MAX_NUM_FILTER_ANNOTATIONS - len(anno_ids_x))
                        cat_y = [
                            0 if x else 1
                            for x in self.indices.anno.get_categories(
                                [pos_2],
                                self.session_default["dividend"],
                                anno_ids_y,
                            )
                        ] + [0] * (MAX_NUM_FILTER_ANNOTATIONS - len(anno_ids_y))
                    act_pos_1 = int(pos_1) // self.session_default["dividend"]
                    act_pos_2 = int(pos_2) // self.session_default["dividend"]

                    map_q = 1
                    same_strand_idx = 0
                    y_strand_idx = 0

                    cat_pos = [
                        item for sublist in zip(cat_x, cat_y) for item in sublist
                    ]
                    cat_hash = sum(x * 2**i for i, x in enumerate(cat_pos))
                    categoires_min = min(categoires_min, cat_hash)
                    categories_max = max(categories_max, cat_hash)

                    start = [
                        act_pos_1,
                        act_pos_2,
                        MAP_Q_MAX - map_q - 1,
                        *cat_pos,
                        same_strand_idx,
                        y_strand_idx,
                    ]
                    end = [
                        act_pos_1,
                        act_pos_2,
                        MAP_Q_MAX - map_q - 1,
                        *cat_pos,
                        same_strand_idx,
                        y_strand_idx,
                    ]
                    self.indices.insert(start, end, bin_cnt)
                id = self.indices.generate(
                    fac=-2 if shekelyan else -1, verbosity=GENERATE_VERBOSITY
                )
                if first_id is None:
                    first_id = id

        if has_lower_triangle != has_upper_triangle and first_id == 0:
            print(
                "Info: Detected an exclusive",
                "lower" if has_lower_triangle else "upper",
                "triangle matrix. Setting symmetry and axis labels of index.",
            )
            self.session["settings"]["interface"]["axis_labels"] = "DNA_DNA"
            self.session["settings"]["filters"]["symmetry"] = "mirror"
            # @todo this should update default.json as well

        if total_reads == 0:
            print(
                "WARNING: the total number of reads that were added to the index is zero! Something seems off..."
            )
        if categories_max == categoires_min and not no_category and total_reads > 1:
            print(
                "WARNING: All interactions had the same annotation overlap (i.e. all interactions are overlapping some gene OR all interactions are overlapping nothing. NOT: some interactions are overlapping a gene and some others are overlapping nothing). Consider using --no_anno if this dataset should not be filterable by annotation overlap.",
            )

        self.set_session(["replicates", "by_name", name, "first_dataset_id"], first_id)
        self.set_session(["replicates", "by_name", name, "num_datasets"], num_itr)
        self.set_session(["replicates", "by_name", name, "total_reads"], total_reads)
        self.set_session(["replicates", "by_name", name, "no_map_q"], True)
        self.set_session(["replicates", "by_name", name, "no_groups"], False)
        self.set_session(["replicates", "by_name", name, "no_multi_map"], True)
        self.set_session(["replicates", "by_name", name, "no_category"], no_category)
        self.set_session(["replicates", "by_name", name, "no_strand"], True)
        self.set_session(["replicates", "by_name", name, "shekelyan"], shekelyan)
        self.set_session(
            ["replicates", "by_name", name, "force_upper_triangle"],
            force_upper_triangle,
        )
        self.set_session(
            ["replicates", "by_name", name, "has_upper_triangle"], has_upper_triangle
        )
        self.set_session(
            ["replicates", "by_name", name, "has_lower_triangle"], has_lower_triangle
        )

        self.indices.clear_points_and_desc()

        self.save_session()

        self.progress_print("done generating dataset.", force_print=True)

    def add_normalization(
        self,
        path,
        name,
        group="a",
        no_groups=False,
        keep_points=False,
        only_points=False,
        no_map_q=False,
        no_multi_map=False,
        no_category=False,
        no_strand=False,
        shekelyan=False,
        columns=["chr", "pos"],
        allow_col_change=False,
    ):
        if not self.name_unique(name):
            raise RuntimeError(
                "The track name you provide must be unique but is not. "
                + "Use the <list> command to see all tracks."
            )

        self.progress_print(
            "generating track",
            force_print=True,
        )

        self.append_session(["coverage", "list"], name)
        self.set_session(
            ["coverage", "by_name", name],
            {
                "ids": {},
                "path": path,
            },
        )
        if group in ["col", "both"]:
            self.append_session(["coverage", "in_column"], name)
        if group in ["row", "both"]:
            self.append_session(["coverage", "in_row"], name)

        contigs = self.session_default["contigs"]["list"]

        read_iterator = chr_order_coverage(
            self.prefix,
            name,
            path,
            0,  # unused
            contigs,
            no_groups,
            "test" in self.session_default,
            lambda *x: self.progress_print("loading", *x),
            columns,
            allow_col_change,
        )
        t_dict = self.get_map_q_thresholds()
        total_reads = 0

        num_itr = len(contigs)
        cnt = 0
        fist_id = None
        count_matrix_warning_done = False

        map_q_min = 255
        map_q_max = 0
        has_multimapper = False
        categoires_min = 2 ** (MAX_NUM_FILTER_ANNOTATIONS * 2 + 1)
        categories_max = 0
        strand_min = 2
        strand_max = 0

        for idx_x, chr_x in enumerate(contigs):
            anno_ids = [
                self.session_default["annotation"]["by_name"][anno] + idx_x
                for anno in self.session_default["annotation"]["filterable"]
            ]
            cnt += 1
            self.progress_print(
                "generating track for contig",
                chr_x + ".",
                cnt,
                "of",
                str(num_itr) + ":",
                str(round(100 * cnt / num_itr, 2)) + "%",
            )
            for (
                read_name,
                pos_1_s,
                pos_1_e,
                pos_1_l,
                strand_1,
                map_q,
                bin_cnt,
            ) in read_iterator.itr_cell(chr_x):
                if (
                    int(bin_cnt) > 1
                    and any(
                        int(p) % self.session_default["dividend"] != 0
                        for p in [pos_1_s, pos_1_e]
                    )
                    and not count_matrix_warning_done
                ):
                    print(
                        "WARNING: The input file has a count column (i.e could be a count matrix)",
                        "and the bin position of at least one row does not match the minimal resolution",
                        "of the index.",
                        "pos: ",
                        pos_1_s,
                        "..",
                        pos_1_e,
                        "count:",
                        bin_cnt,
                        "minimal index resolution:",
                        str(self.session_default["dividend"])
                        + ". Will not show this warning again.",
                        file=sys.stderr,
                    )
                    count_matrix_warning_done = True
                total_reads += 1
                if no_category:
                    cat = [0] * MAX_NUM_FILTER_ANNOTATIONS
                else:
                    cat = self.indices.anno.get_categories(
                        [0 if x else 1 for x in pos_1_l.split(",")],
                        self.session_default["dividend"],
                        anno_ids,
                    ) + [0] * (MAX_NUM_FILTER_ANNOTATIONS - len(anno_ids))

                act_pos_1_s = int(pos_1_s) // self.session_default["dividend"]
                if no_multi_map:
                    act_pos_1_e = act_pos_1_s
                else:
                    act_pos_1_e = int(pos_1_e) // self.session_default["dividend"]
                if act_pos_1_e != act_pos_1_s:
                    has_multimapper = True

                map_q = t_dict[int(map_q)]
                if no_map_q:
                    map_q = 1
                map_q_min = min(map_q_min, map_q)
                map_q_max = max(map_q_max, map_q)

                if no_strand:
                    strand_idx = 0
                else:
                    strand_idx = 0 if bool(strand_1) else 1
                strand_min = min(strand_min, strand_idx)
                strand_max = max(strand_max, strand_idx)

                bin_cnt = int(bin_cnt)

                cat_pos = [
                    item for sublist in zip(cat, [0] * len(cat)) for item in sublist
                ]
                cat_hash = sum(x * 2**i for i, x in enumerate(cat_pos))
                categoires_min = min(categoires_min, cat_hash)
                categories_max = max(categories_max, cat_hash)

                start = [
                    act_pos_1_s,
                    0,
                    MAP_Q_MAX - map_q - 1,
                    *cat_pos,
                    strand_idx,
                    0,
                ]
                end = [
                    act_pos_1_e,
                    0,
                    MAP_Q_MAX - map_q - 1,
                    *cat_pos,
                    strand_idx,
                    0,
                ]
                self.indices.insert(start, end, bin_cnt)
            id = self.indices.generate(
                fac=-2 if shekelyan else -1, verbosity=GENERATE_VERBOSITY
            )
            if fist_id is None:
                fist_id = id

        if total_reads == 0:
            print(
                "WARNING: the total number of reads that were added to the index is zero! Something seems off..."
            )
        if not no_map_q and map_q_min == map_q_max and not total_reads > 1:
            print(
                "WARNING: Detected only a single mapping quality value (",
                map_q_min,
                "). Consider using --no_map_q if this dataset should not be filterable by mapping quality.",
            )
        if not has_multimapper and not no_multi_map and not total_reads > 1:
            print(
                "WARNING: Detected no multimappers. Consider using --no_multi_map if this dataset should not be filterable by multimapping.",
            )
        if categories_max == categoires_min and not no_category and not total_reads > 1:
            print(
                "WARNING: Detected only a single annotation overlap type. Consider using --no_anno if this dataset should not be filterable by annotation overlap.",
            )
        if strand_min == strand_max and not no_strand and not total_reads > 1:
            print(
                "WARNING: Detected only a single strand overlap type. Consider using --no_strand if this dataset should not be filterable by strand overlap.",
            )

        self.set_session(["coverage", "by_name", name, "first_dataset_id"], fist_id)
        self.set_session(["coverage", "by_name", name, "num_datasets"], num_itr)
        self.set_session(["coverage", "by_name", name, "total_reads"], total_reads)
        self.set_session(["coverage", "by_name", name, "no_map_q"], no_map_q)
        self.set_session(["coverage", "by_name", name, "no_groups"], no_groups)
        self.set_session(["coverage", "by_name", name, "no_multi_map"], no_multi_map)
        self.set_session(["coverage", "by_name", name, "no_category"], no_category)
        self.set_session(["coverage", "by_name", name, "no_strand"], no_strand)
        self.set_session(["coverage", "by_name", name, "shekelyan"], shekelyan)

        read_iterator.cleanup()

        if not keep_points:
            self.indices.clear_points_and_desc()

        self.save_session()

        self.progress_print("done generating track.", force_print=True)
