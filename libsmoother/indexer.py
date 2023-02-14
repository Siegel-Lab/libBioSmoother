from ._import_lib_smoother_cpp import Index
from ._parse_and_group_reads import *
import json
import os
import copy
import time

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
        if self.last_prog_print is None or t - self.last_prog_print >= PROGRESS_PRINT_TIME or force_print:
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

    def create_session(self, chr_len_file_name, dividend, anno_path, annotation_order=None, test=False):
        if os.path.exists(self.prefix):
            print("ERROR: The given index already exists.")
            exit()
        os.makedirs(self.prefix)
        self.set_session(["dividend"], dividend)
        self.set_session(["previous"], None)
        self.set_session(["next"], None)
        self.set_session(["settings"], None)
        self.set_session(["replicates"], {
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
        })
        self.set_session(["coverage"], {
            "list": [],
            "by_name": {},
            "in_column": [],
            "in_row": [],
            "cov_column_a": [],
            "cov_column_b": [],
            "cov_row_a": [],
            "cov_row_b": [],
        })
        self.set_session(["annotation"], {
            "list": [],
            "by_name": {},
            "visible_x": [],
            "visible_y": [],
            "filter": None,
            "filter_row": False,
            "filter_col": True,
        })
        self.set_session(["contigs"], {
            "list": [],
            "lengths": {},
            "displayed_on_x": [],
            "displayed_on_y": [],
            "genome_size": 0,
            "column_coordinates": "full_genome",
            "row_coordinates": "full_genome",
        })
        if test:
            self.set_session(["test"], True)

        with open(chr_len_file_name, "r") as len_file:
            for line in len_file:
                chr_name, chr_len = line.split()
                if not test or ("Chr1" in chr_name):
                    self.append_session(["contigs", "list"], chr_name)
                    self.set_session(["contigs", "lengths", chr_name], int(chr_len) // dividend)
                    self.append_session(["contigs", "displayed_on_x"], chr_name)
                    self.append_session(["contigs", "displayed_on_y"], chr_name)
                    self.set_session(["contigs", "genome_size"],
                                     self.session_default["contigs"]["genome_size"] + int(chr_len))
        self.set_session(["area"],{
            "x_start": 0,
            "y_start": 0,
            "x_end": self.session_default["contigs"]["genome_size"] // dividend,
            "y_end": self.session_default["contigs"]["genome_size"] // dividend,
        })

        for o, d in [(6, 2)]:
            idx_suff = str(o) + "." + str(d)
            touch(self.prefix + "/" + idx_suff + ".coords")
            touch(self.prefix + "/" + idx_suff + ".datsets")
            touch(self.prefix + "/" + idx_suff + ".overlays")
            touch(self.prefix + "/" + idx_suff + ".prefix_sums")

        self.save_session()
        self.try_load_index()

        sorted_list = {}
        for name, chrom, start, end, info, on_forw_strnd in parse_annotations(anno_path):
            if not chrom in self.session_default["contigs"]["list"]:
                continue
            if name not in sorted_list:
                sorted_list[name] = {}
            if chrom not in sorted_list[name]:
                sorted_list[name][chrom] = []
            sorted_list[name][chrom].append((start, end, info, on_forw_strnd))
        order = []
        if annotation_order is None:
            if "gene" in sorted_list:
                order.append("gene")
        else:
            with open(annotation_order, "r") as anno_order_file:
                for line in anno_order_file:
                    if line in sorted_list:
                        order.append(line)
        for name in sorted(list(sorted_list.keys())):
            if not name in order:
                order.append(name)
        for name in order:
            chroms = sorted_list[name]
            if name not in self.session_default["annotation"]["list"]:
                self.append_session(["annotation", "list"], name)
                self.append_session(["annotation", "visible_x"], name)
                self.append_session(["annotation", "visible_y"], name)
                self.set_session(["annotation", "by_name", name], {})

                for chrom, annos in chroms.items():
                    self.progress_print("annotating", name + "(s)", "for contig", chrom)
                    self.set_session(["annotation", "by_name", name, chrom], 
                                     self.indices.anno.add_intervals(annos, self.session_default["dividend"],
                                                                     verbosity=GENERATE_VERBOSITY))
            else:
                raise RuntimeError("annotation with this name already exists")

        self.save_session()

    def name_unique(self, name):
        return (
            name not in self.session_default["replicates"]["list"]
            and name not in self.session_default["coverage"]["list"]
        )

    def __get_cat_indices_1d(self, cats):
        doubles = -1
        no_anno = True
        for idx, c in enumerate(cats):
            if c:
                no_anno = False
                idx_list.append(idx*3+1)
                doubles += 1
        if no_anno:
            return [len(cat_x)*3]
        if doubles > 0:
            for _ in range(0, doubles):
                idx_list.append(len(cat_x)*3 + 1)
        return idx_list

    def __get_cat_indices_2d(self, cat_x, cat_y):
        doubles = -1
        no_anno = True
        cat_to_idx = {
            True: {True: 1, False: 0},
            False: {True: 2}
        }
        idx_list = []
        for idx, (cx, cy) in enumerate(zip(cat_x, cat_y)):
            if cx or cy:
                no_anno = False
                idx_list.append((cat_to_idx[cx][cy]+idx*3, 1))
                doubles += 1
        if no_anno:
            return [(max(len(cat_x), len(cat_y))*3, 1)]
        if doubles > 0:
            for _ in range(0, doubles):
                idx_list.append((max(len(cat_x), len(cat_y))*3 + 1, -1))
        return idx_list


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
        shekelyan=False,
        force_upper_triangle=False,
    ):
        if not self.name_unique(name):
            raise RuntimeError(
                "The dataset name you provide must be unique but is not. "
                + "Use the <list> command to see all datasets."
            )

        self.progress_print(
            "generating replicate.",
            force_print=True
        )

        self.append_session(["replicates", "list"], name)
        self.set_session(["replicates", "by_name", name], {
            "ids": {},
            "path": path
        })
        if group in ["a", "both"]:
            self.append_session(["replicates", "in_group_a"], name)
        if group in ["b", "both"]:
            self.append_session(["replicates", "in_group_b"], name)

        read_iterator = chr_order_heatmap(
            self.prefix,
            name,
            path,
            get_filesize(path),
            self.session_default["contigs"]["list"],
            no_groups,
            "test" in self.session_default,
            force_upper_triangle,
            lambda *x: self.progress_print("loading", *x)
        )
        total_reads = 0


        num_itr = len(read_iterator)
        cnt = 0
        for chr_x in read_iterator.itr_x_axis():
            anno_ids_x = [self.session_default["annotation"]["by_name"][anno][chr_x] 
                            if chr_x in self.session_default["annotation"]["by_name"][anno]
                            else -1
                            for anno in self.session_default["annotation"]["list"]]
            for chr_y in read_iterator.itr_y_axis():
                anno_ids_y = [self.session_default["annotation"]["by_name"][anno][chr_y]
                                if chr_y in self.session_default["annotation"]["by_name"][anno]
                                else -1
                                for anno in self.session_default["annotation"]["list"]]
                cnt += 1
                self.progress_print("generating heatmap for contig-pair", chr_x, chr_y + ".", cnt, "of", num_itr,
                                    str(round(100*cnt/num_itr, 2)) + "%")
                for (
                    read_name,
                    pos_1_s,
                    pos_1_e,
                    pos_2_s,
                    pos_2_e,
                    map_q,
                ) in read_iterator.itr_cell(chr_x, chr_y):
                    total_reads += 1
                    if no_category:
                        cat_x = [False] * len(self.session_default["annotation"]["list"])
                        cat_y = cat_x
                    else:
                        cat_x = self.indices.anno.get_categories(pos_1_s, pos_1_e,
                                                                self.session_default["dividend"], anno_ids_x)
                        cat_y = self.indices.anno.get_categories(pos_2_s, pos_2_e,
                                                                self.session_default["dividend"], anno_ids_y)
                    act_pos_1_s = int(pos_2_s) // self.session_default["dividend"]
                    act_pos_2_s = int(pos_1_s) // self.session_default["dividend"]
                    if no_multi_map:
                        act_pos_1_e = act_pos_1_s
                        act_pos_2_e = act_pos_2_s
                    else:
                        act_pos_1_e = int(pos_2_e) // self.session_default["dividend"]
                        act_pos_2_e = int(pos_1_e) // self.session_default["dividend"]
                    if no_map_q:
                        map_q = 1

                    for cat_idx, val in self.__get_cat_indices_2d(cat_x, cat_y):
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - int(map_q) - 1, cat_idx]
                        end = [act_pos_1_e, act_pos_2_e, MAP_Q_MAX - int(map_q) - 1, cat_idx]
                        self.indices.insert(start, end, val)
                if (
                    chr_x
                    not in self.session_default["replicates"]["by_name"][name]["ids"]
                ):
                    self.set_session(["replicates", "by_name", name, "ids", chr_x], {})
                assert (
                    chr_y
                    not in self.session_default["replicates"]["by_name"][name]["ids"][
                        chr_x
                    ]
                )
                self.set_session(["replicates", "by_name", name, "ids", chr_x, chr_y], 
                                    self.indices.generate(fac=-2 if shekelyan else -1,
                                                          verbosity=GENERATE_VERBOSITY))

        self.set_session(["replicates", "by_name", name, "total_reads"], total_reads)
        if total_reads == 0:
            print("WARNING: the total number of reads that were added to the index is zero! Something seems off...")

        read_iterator.cleanup()

        if not keep_points:
            self.indices.clear_points_and_desc()

        self.save_session()

        print("done generating index")

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
    ):
        raise RuntimeError(
            "Categories need to be fixed first"
        )
        if not self.name_unique(name):
            raise RuntimeError(
                "The track name you provide must be unique but is not. "
                + "Use the <list> command to see all tracks."
            )

        if not (no_map_q and no_multi_map):
            self.progress_print("pre-scanning file for index parameters...", force_print=True)
            has_map_q, multi_map = has_map_q_and_multi_map(
                path,
                "test" in self.session_default,
                self.session_default["contigs"]["list"],
                lambda *x: self.progress_print("scanning", *x),
                parse_func=parse_track
            )
        if no_map_q:
            has_map_q = False
        if no_multi_map:
            multi_map = False

        self.progress_print(
            "generating track",
            "with" if has_map_q else "without",
            "mapping quality and",
            "with" if multi_map else "without",
            "multi mapping.",
            force_print=True
        )

        self.append_session(["coverage", "list"], name)
        self.set_session(["coverage", "by_name", name], {
            "ids": {},
            "has_map_q": has_map_q,
            "has_multimapping": multi_map,
            "path": path,
        })
        if group in ["col", "both"]:
            self.append_session(["coverage", "in_column"], name)
        if group in ["row", "both"]:
            self.append_session(["coverage", "in_row"], name)

        o = (1 if multi_map else 0)
        d = o + (2 if has_map_q else 1)

        read_iterator = chr_order_coverage(
            self.prefix + ".smoother_index",
            name,
            path,
            get_filesize(path),
            self.session_default["contigs"]["list"],
            no_groups,
            "test" in self.session_default,
            lambda *x: self.progress_print("loading", *x)
        )
        total_reads = 0

        cnt = 0
        num_itr = len(read_iterator)
        for chr_x in read_iterator.itr_x_axis():
            anno_ids = [self.session_default["annotation"]["by_name"][anno][chr_x] 
                            for anno in self.session_default["annotation"]["list"]
                            if chr_x in self.session_default["annotation"]["by_name"][anno]]
            cnt += 1
            self.progress_print("generating track for contig", chr_x + ".", cnt, "of", num_itr, 
                                str(round(100*cnt/num_itr, 2)) + "%")
            for (
                read_name,
                pos_1_s,
                pos_1_e,
                map_q,
            ) in read_iterator.itr_cell(chr_x):
                total_reads += 1
                if no_cat:
                    cat = []
                else:
                    cat = self.indices.anno.get_categories(pos_1_s, pos_1_e, self.session_default["dividend"], anno_ids)
                act_pos_1_s = int(pos_1_s) // self.session_default["dividend"]
                act_pos_1_e = int(pos_1_e) // self.session_default["dividend"]
                if has_map_q and multi_map:
                    start = [*cat, act_pos_1_s, MAP_Q_MAX - int(map_q) - 1]
                    end = [*cat, act_pos_1_e, MAP_Q_MAX - int(map_q) - 1]
                elif has_map_q and not multi_map:
                    start = [*cat, act_pos_1_s, MAP_Q_MAX - int(map_q) - 1]
                    end = [*cat, 0, 0]
                elif not has_map_q and multi_map:
                    start = [*cat, act_pos_1_s]
                    end = [*cat, act_pos_1_e]
                elif not has_map_q and not multi_map:
                    start = [*cat, act_pos_1_s]
                    end = [*cat, 0]
                else:
                    raise RuntimeError("this statement should never be reached")
                self.indices.insert(d + len(anno_ids), o, start, end)
            self.set_session(["coverage", "by_name", name, "ids", chr_x], 
                                self.indices.generate(d + len(anno_ids), o, verbosity=GENERATE_VERBOSITY))

        self.set_session(["coverage", "by_name", name, "total_reads"], total_reads)

        if total_reads == 0:
            print("WARNING: the total number of reads that were added to the index is zero! Something seems off...")

        read_iterator.cleanup()

        if not keep_points:
            self.indices.clear_points_and_desc()

        self.save_session()

        print("done generating track")
