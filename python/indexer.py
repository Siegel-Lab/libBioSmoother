from libContactMapping import CachedSpsInterface
import json


class Indexer:
    def __init__(self, prefix):
        self.prefix = prefix
        self.indices = CachedSpsInterface(prefix + ".smoother_index/")
        if not self.indices.loaded():
            raise RuntimeError(
                "indices with prefix " + prefix + " could not be opened."
            )
        if os.path.exists(prefix + ".smoother_index/default_session.json"):
            with open(prefix + ".smoother_index/default_session.json", "r") as f:
                self.session_default = json.load(f)
        else:
            self.session_default = {}

    def save_session(self):
        with open(self.prefix + ".smoother_index/default_session.json", "w") as f:
            json.dump(self.session_default, f)

    def create_session(
        self, chr_len_file_name, annotation_filename, dividend, test=False
    ):
        os.makedirs(self.prefix + ".smoother_index")
        self.session_default["dividend"] = dividend
        self.session_default["previous"] = None
        self.session_default["next"] = None
        self.session_default["settings"] = None  # @todo get fac default
        self.session_default["replicates"] = {
            "list": [],
            "by_name": {},
            "in_group": {},
            "coverage": {},
        }
        self.session_default["coverage"] = {
            "list": [],
            "by_name": {},
            "coverage": {},
        }
        self.session_default["contigs"] = {
            "list": [],
            "lengths": {},
            "displayed": {},
            "genome_size": 0,
        }
        self.session_default["annotation"] = {"list": [], "visible": {}, "ids": {}}

        with open(file_name, "r") as chr_len_file_name:
            for line in len_file:
                chr_name, chr_len = line.split()
                if not test or ("Chr1_3A" in chr_name):
                    self.session_default["contigs"]["list"].append(chr_name)
                    self.session_default["contigs"]["lengths"][chr_name] = int(chr_len)
                    self.session_default["contigs"]["displayed"][chr_name][
                        "displayed_on_x"
                    ] = True
                    self.session_default["contigs"]["displayed"][chr_name][
                        "displayed_on_y"
                    ] = True
                    self.session_default["contigs"]["genome_size"] += int(chr_len)

        for o, d in [(1, 0), (2, 1), (2, 0), (3, 1), (4, 2), (3, 0), (5, 2)]:
            idx_suff = str(o) + "." + str(d)
            touch(out_prefix + ".smoother_index/" + idx_suff + ".coords")
            touch(out_prefix + ".smoother_index/" + idx_suff + ".datsets")
            touch(out_prefix + ".smoother_index/" + idx_suff + ".overlays")
            touch(out_prefix + ".smoother_index/" + idx_suff + ".prefix_sums")

        self.save_session()

    def add_annotation(self, file_name):
        sorted_list = {}
        for name, chrom, start, end, info in parse_annotations(file_name):
            if not chrom in self.session_default["contigs"]["list"]:
                continue
            if name not in sorted_list:
                sorted_list[name] = {}
            if chrom not in sorted_list[name]:
                sorted_list[name][chrom] = []
            sorted_list[name][chrom].append((start, end, info))
        for name, chroms in sorted_list.items():
            if name not in self.session_default["annotation"]["list"]:
                self.session_default["annotation"]["list"].append(name)
                self.session_default["annotation"]["visible"][name] = {
                    "x": True,
                    "y": True,
                }
                self.session_default["annotation"]["ids"][name] = {}

                for chrom, annos in chroms.items():
                    for start, end, info in annos:
                        self.indices.insert(
                            1,
                            1,
                            [start // self.session_default["dividend"]],
                            [end // self.session_default["dividend"]],
                            info,
                        )
                    self.session_default["annotation"]["ids"][name][
                        chrom
                    ] = self.indices.generate(1, 1)
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
        test=False,
        cached=False,
        no_groups=False,
        without_dep_dim=True,
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
            print("pre-scanning file for index parameters...")
            has_map_q, multi_map = has_map_q_and_multi_map(
                path, test, self.session_default["contigs"]["list"]
            )
        if no_map_q:
            has_map_q = False
        if no_multi_map:
            multi_map = False
        print(
            "generating index",
            "with" if has_map_q else "without",
            "mapping quality and",
            "with" if multi_map else "without",
            "multi mapping.",
        )

        self.session_default["replicates"]["list"].append(name)
        self.session_default["replicates"]["by_name"][name] = {
            "ids": {},
            "has_map_q": has_map_q,
            "has_multimapping": multi_map,
            "path": path,
        }
        self.session_default["replicates"]["coverage"][name] = {
            "in_column": False,
            "in_row": False,
            "cov_column_a": False,
            "cov_row_a": False,
            "cov_column_b": False,
            "cov_row_b": False,
        }
        self.session_default["replicates"]["in_group"][name] = {
            "in_group_a": group in ["a", "both"],
            "in_group_b": group in ["b", "both"],
        }

        o = 2 if multi_map else 0
        d = o + 3 if has_map_q else 2

        read_iterator = chr_order_heatmap(
            self.prefix + ".smoother_index",
            name,
            path,
            get_filesize(path),
            meta.chr_sizes.chr_sizes,
            no_groups,
            test,
        )

        for chr_x in read_iterator.itr_x_axis():
            for chr_y in read_iterator.itr_y_axis():
                for (
                    read_name,
                    pos_1_s,
                    pos_1_e,
                    pos_2_s,
                    pos_2_e,
                    map_q,
                ) in read_iterator.itr_cell(chr_x, chr_y):
                    act_pos_1_s = pos_2_s // meta.dividend
                    act_pos_1_e = pos_2_e // meta.dividend
                    act_pos_2_s = pos_1_s // meta.dividend
                    act_pos_2_e = pos_1_e // meta.dividend
                    if has_map_q and multi_map:
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - map_q - 1]
                        end = [act_pos_1_e, act_pos_2_e, MAP_Q_MAX - map_q - 1]
                    elif has_map_q and not multi_map:
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - map_q - 1]
                        end = None
                    elif not has_map_q and multi_map:
                        start = [act_pos_1_s, act_pos_2_s]
                        end = [act_pos_1_e, act_pos_2_e]
                    elif not has_map_q and not multi_map:
                        start = [act_pos_1_s, act_pos_2_s]
                        end = None
                    else:
                        raise RuntimeError("this statement should never be reached")
                    self.indices.insert(d, o, start, end)
                if (
                    chr_x
                    not in self.session_default["replicates"]["by_name"][name]["ids"]
                ):
                    self.session_default["replicates"]["by_name"][name]["ids"][
                        chr_x
                    ] = {}
                assert (
                    chr_y
                    not in self.session_default["replicates"]["by_name"][name]["ids"][
                        chr_x
                    ]
                )
                self.session_default["replicates"]["by_name"][name]["ids"][chr_x][
                    chr_y
                ] = self.indices.generate(d, o)

        o = 1 if multi_map else 0
        d = o + 2 if has_map_q else 1

        for x_axis in [True, False]:
            for chr_ in (
                read_iterator.itr_x_axis() if x_axis else read_iterator.itr_y_axis()
            ):
                for (read_name, pos_1_s, pos_1_e, pos_2_s, pos_2_e, map_q,) in (
                    read_iterator.itr_row(chr_)
                    if x_axis
                    else read_iterator.itr_col(chr_)
                ):
                    act_pos_1_s = pos_2_s // meta.dividend
                    act_pos_1_e = pos_2_e // meta.dividend
                    act_pos_2_s = pos_1_s // meta.dividend
                    act_pos_2_e = pos_1_e // meta.dividend
                    if has_map_q and multi_map:
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - map_q - 1]
                        end = [act_pos_1_e, act_pos_2_e, MAP_Q_MAX - map_q - 1]
                    elif has_map_q and not multi_map:
                        start = [act_pos_1_s, act_pos_2_s, MAP_Q_MAX - map_q - 1]
                        end = None
                    elif not has_map_q and multi_map:
                        start = [act_pos_1_s, act_pos_2_s]
                        end = [act_pos_1_e, act_pos_2_e]
                    elif not has_map_q and not multi_map:
                        start = [act_pos_1_s, act_pos_2_s]
                        end = None
                    else:
                        raise RuntimeError("this statement should never be reached")
                    self.indices.insert(d, o, start, end)

                self.session_default["replicates"]["by_name"][name]["ids"][chr_][
                    "row" if x_axis else "col"
                ] = self.indices.generate(d, o)

        read_iterator.cleanup()

        if not keep_points:
            self.indices.clear_points_and_desc()

        self.save_session()

        print("done generating index")
