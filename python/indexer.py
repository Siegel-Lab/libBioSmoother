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
        with open(prefix + ".smoother_index/default_session.json", "r") as f:
            self.session_default = json.load(f)

    def save(self):
        with open(self.prefix + ".smoother_index/default_session.json", "w") as f:
            json.dump(self.session_default, f)

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
            "path": path
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

        o2 = 1 if multi_map else 0
        d2 = o + 2 if has_map_q else 1

        last_chr = None

        def generate():
            if not last_chr is None:
                if (
                    last_chr[0]
                    not in self.session_default["replicates"]["by_name"][name][
                        "ids"
                    ]
                ):
                    self.session_default["replicates"]["by_name"][name]["ids"][
                        last_chr[0]
                    ] = {}
                assert last_chr[1] not in self.session_default["replicates"]["by_name"][name]["ids"][
                    last_chr[0]
                ]
                self.session_default["replicates"]["by_name"][name]["ids"][
                    last_chr[0]
                ][last_chr[1]] = self.indices.generate(d, o)
        for (
            read_name,
            chr_1,
            pos_1_s,
            pos_1_e,
            chr_2,
            pos_2_s,
            pos_2_e,
            map_q,
        ) in group_heatmap(
            path, get_filesize(path), meta.chr_sizes.chr_sizes, no_groups, test
        ):
            if (chr_1, chr_2) != last_chr:
                generate()
                last_chr = (chr_1, chr_2)
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

        if not last_chr is None:
            generate()
        print("done generating index")

        meta.save(out_prefix + ".smoother_index/meta")
        if not keep_points:
            del index
            os.remove(out_prefix + ".smoother_index/repl" + idx_suff + ".points")
            os.remove(out_prefix + ".smoother_index/repl" + idx_suff + ".desc")
