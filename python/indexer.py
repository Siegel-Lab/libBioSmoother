from libContactMapping import CachedSpsInterface, DiskSpsInterface
from _parse_and_group_reads import *
import json
import os

MAP_Q_MAX = 255

def touch(f_name):
    with open(f_name, "a"):  # Create file if does not exist
        pass

GENERATE_VERBOSITY = 1

class Indexer:
    def __init__(self, prefix):
        self.prefix = prefix
        if os.path.exists(prefix + ".smoother_index/default_session.json"):
            # self.indices = CachedSpsInterface(prefix + ".smoother_index")
            self.indices = DiskSpsInterface(prefix + ".smoother_index", True)
            with open(prefix + ".smoother_index/default_session.json", "r") as f:
                self.session_default = json.load(f)
            with open(prefix + ".smoother_index/session.json", "r") as f:
                self.session = json.load(f)
        else:
            self.indices = None
            self.session_default = {}
            self.session = {}

    def progress_print(self, *text):
        print(*text)

    def save_session(self):
        with open(self.prefix + ".smoother_index/default_session.json", "w") as f:
            json.dump(self.session_default, f)
        with open(self.prefix + ".smoother_index/session.json", "w") as f:
            json.dump(self.session, f)
        del self.indices

    def set_session(self, keys, val):
        curr = self.session
        curr_def = self.session_default
        for key in keys[:-1]:
            curr = curr[key]
            curr_def = curr_def[key]
        curr[keys[-1]] = val
        curr_def[keys[-1]] = val

    def append_session(self, keys, val):
        curr = self.session
        curr_def = self.session_default
        for key in keys[:-1]:
            curr = curr[key]
            curr_def = curr_def[key]
        curr[keys[-1]].append(val)
        curr_def[keys[-1]].append(val)

    def create_session(self, chr_len_file_name, dividend, test=False):
        os.makedirs(self.prefix + ".smoother_index")
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
            "row_filter": [],
            "col_filter": [],
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
                if not test or ("Chr1_3A" in chr_name):
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

        for o, d in [(1, 0), (2, 1), (2, 0), (3, 1), (4, 2), (3, 0), (5, 2)]:
            idx_suff = str(o) + "." + str(d)
            touch(self.prefix + ".smoother_index/" + idx_suff + ".coords")
            touch(self.prefix + ".smoother_index/" + idx_suff + ".datsets")
            touch(self.prefix + ".smoother_index/" + idx_suff + ".overlays")
            touch(self.prefix + ".smoother_index/" + idx_suff + ".prefix_sums")

        self.save_session()

    def add_annotation(self, file_name):
        sorted_list = {}
        for name, chrom, start, end, info, on_forw_strnd in parse_annotations(file_name):
            if not chrom in self.session_default["contigs"]["list"]:
                continue
            if name not in sorted_list:
                sorted_list[name] = {}
            if chrom not in sorted_list[name]:
                sorted_list[name][chrom] = []
            sorted_list[name][chrom].append((start // self.session_default["dividend"], 
                                             end // self.session_default["dividend"], 
                                            info, on_forw_strnd))
        for name, chroms in sorted_list.items():
            if name not in self.session_default["annotation"]["list"]:
                self.append_session(["annotation", "list"], name)
                self.append_session(["annotation", "visible_x"], name)
                self.append_session(["annotation", "visible_y"], name)
                self.set_session(["annotation", "by_name", name], {})

                for chrom, annos in chroms.items():
                    self.progress_print("annotating", name + "(s)", "for contig", chrom)
                    self.set_session(["annotation", "by_name", name, chrom], 
                                     self.indices.anno.add_intervals(annos, verbosity=GENERATE_VERBOSITY))
            else:
                raise RuntimeError("annotation with this name already exists")

        self.save_session()

    def name_unique(self, name):
        return (
            name not in self.session_default["replicates"]["list"]
            and name not in self.session_default["coverage"]["list"]
        )

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
    ):
        if not self.name_unique(name):
            raise RuntimeError(
                "The dataset name you provide must be unique but is not. "
                + "Use the <list> command to see all datasets."
            )

        if not (no_map_q and no_multi_map):
            self.progress_print("pre-scanning file for index parameters...")
            has_map_q, multi_map = has_map_q_and_multi_map(
                path,
                "test" in self.session_default,
                self.session_default["contigs"]["list"],
            )
        if no_map_q:
            has_map_q = False
        if no_multi_map:
            multi_map = False

        self.progress_print(
            "generating index",
            "with" if has_map_q else "without",
            "mapping quality and",
            "with" if multi_map else "without",
            "multi mapping.",
        )

        self.append_session(["replicates", "list"], name)
        self.set_session(["replicates", "by_name", name], {
            "ids": {},
            "has_map_q": has_map_q,
            "has_multimapping": multi_map,
            "path": path,
        })
        if group in ["a", "both"]:
            self.append_session(["replicates", "in_group_a"], name)
        if group in ["b", "both"]:
            self.append_session(["replicates", "in_group_b"], name)

        o = 2 if multi_map else 0
        d = o + 3 if has_map_q else 2

        read_iterator = chr_order_heatmap(
            self.prefix + ".smoother_index",
            name,
            path,
            get_filesize(path),
            self.session_default["contigs"]["list"],
            no_groups,
            "test" in self.session_default,
        )
        total_reads = 0

        for chr_x in read_iterator.itr_x_axis():
            for chr_y in read_iterator.itr_y_axis():
                self.progress_print("generating heatmap for contig-pair", chr_x, chr_y)
                for (
                    read_name,
                    pos_1_s,
                    pos_1_e,
                    pos_2_s,
                    pos_2_e,
                    map_q,
                ) in read_iterator.itr_cell(chr_x, chr_y):
                    total_reads += 1
                    act_pos_1_s = int(pos_2_s) // self.session_default["dividend"]
                    act_pos_1_e = int(pos_2_e) // self.session_default["dividend"]
                    act_pos_2_s = int(pos_1_s) // self.session_default["dividend"]
                    act_pos_2_e = int(pos_1_e) // self.session_default["dividend"]
                    if has_map_q and multi_map:
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - int(map_q) - 1]
                        end = [act_pos_1_e, act_pos_2_e, MAP_Q_MAX - int(map_q) - 1]
                    elif has_map_q and not multi_map:
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - int(map_q) - 1]
                        end = [0, 0, 0]
                    elif not has_map_q and multi_map:
                        start = [act_pos_1_s, act_pos_2_s]
                        end = [act_pos_1_e, act_pos_2_e]
                    elif not has_map_q and not multi_map:
                        start = [act_pos_1_s, act_pos_2_s]
                        end = [0, 0]
                    else:
                        raise RuntimeError("this statement should never be reached")
                    self.indices.insert(d, o, start, end)
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
                                    self.indices.generate(d, o, verbosity=GENERATE_VERBOSITY))

        self.set_session(["replicates", "by_name", name, "total_reads"], total_reads)
        o = 1 if multi_map else 0
        d = o + 2 if has_map_q else 1

        for x_axis in [True, False]:
            for chr_ in (
                read_iterator.itr_x_axis() if x_axis else read_iterator.itr_y_axis()
            ):
                self.progress_print("generating tracks for contig", chr_)
                for (read_name, pos_1_s, pos_1_e, pos_2_s, pos_2_e, map_q,) in (
                    read_iterator.itr_row(chr_)
                    if x_axis
                    else read_iterator.itr_col(chr_)
                ):
                    pos_s = pos_1_s if x_axis else pos_2_s
                    pos_e = pos_1_e if x_axis else pos_2_e
                    act_pos_s = int(pos_s) // self.session_default["dividend"]
                    act_pos_e = int(pos_e) // self.session_default["dividend"]
                    if has_map_q and multi_map:
                        start = [act_pos_s, MAP_Q_MAX - int(map_q) - 1]
                        end = [act_pos_e, MAP_Q_MAX - int(map_q) - 1]
                    elif has_map_q and not multi_map:
                        start = [act_pos_s, MAP_Q_MAX - int(map_q) - 1]
                        end = [0, 0]
                    elif not has_map_q and multi_map:
                        start = [act_pos_s]
                        end = [act_pos_e]
                    elif not has_map_q and not multi_map:
                        start = [act_pos_s]
                        end = [0]
                    else:
                        raise RuntimeError("this statement should never be reached")
                    self.indices.insert(d, o, start, end)

                self.set_session(["replicates", "by_name", name, "ids", chr_, "row" if x_axis else "col"], 
                                 self.indices.generate(d, o, verbosity=GENERATE_VERBOSITY))

        read_iterator.cleanup()

        if not keep_points:
            self.indices.clear_points_and_desc()

        self.save_session()

        print("done generating index")
