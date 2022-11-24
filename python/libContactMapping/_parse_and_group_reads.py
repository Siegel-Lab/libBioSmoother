PRINT_MODULO = 10000
import subprocess
import errno
import os

TEST_FAC = 1000000


def simplified_filepath(path):
    if "/" in path:
        x = path[path.rindex("/") + 1 :]
    if "." in path:
        return x[: x.index(".")]
    return x


def read_xa_tag(tags):
    if tags == "notag" or len(tags) < 5:
        return []
    l = []
    for tag in tags[5:].split(";"):
        split = tag.split(",")
        if len(split) == 5:
            chrom, str_pos, _CIGAR, _NM = split
            # strand = str_pos[0]
            pos = int(str_pos[1:])
            l.append([chrom, pos])
    return l


def parse_heatmap(in_filename, test, chr_filter):
    with open(in_filename, "r") as in_file_1:
        cnt = 0
        for line in in_file_1:
            # parse file columns
            num_cols = len(line.split())
            if num_cols == 7:
                read_name, chr_1, pos_1, chr_2, pos_2, mapq_1, mapq_2 = line.split()
                tag_a = "?"
                tag_b = "?"
            elif num_cols == 9:
                (
                    read_name,
                    chr_1,
                    pos_1,
                    chr_2,
                    pos_2,
                    mapq_1,
                    mapq_2,
                    tag_a,
                    tag_b,
                ) = line.split()
            elif num_cols == 11:
                (
                    read_name,
                    _1,
                    chr_1,
                    pos_1,
                    _2,
                    _3,
                    chr_2,
                    pos_2,
                    _4,
                    mapq_1,
                    mapq_2,
                ) = line.split()
                tag_a = "?"
                tag_b = "?"
            elif num_cols == 13:
                (
                    read_name,
                    _1,
                    chr_1,
                    pos_1,
                    _2,
                    _3,
                    chr_2,
                    pos_2,
                    _4,
                    mapq_1,
                    mapq_2,
                    tag_a,
                    tag_b,
                ) = line.split()
            else:
                raise ValueError(
                    'line "'
                    + line
                    + '" has '
                    + str(num_cols)
                    + ", columns which is unexpected. There can be 7, 9, 11, or 13 columns."
                )

            if not chr_1 in chr_filter:
                continue
            if not chr_2 in chr_filter:
                continue
            if mapq_1 == "" or mapq_1 == "nomapq" or mapq_1 == "255" or mapq_1 == "*":
                mapq_1 = 0
            if mapq_2 == "" or mapq_2 == "nomapq" or mapq_2 == "255" or mapq_2 == "*":
                mapq_2 = 0
            # convert number values to ints
            pos_1, pos_2, mapq_1, mapq_2 = (
                int(x) for x in (pos_1, pos_2, mapq_1, mapq_2)
            )
            pos_1 -= 1
            pos_2 -= 1

            if cnt > TEST_FAC and test:
                break
            cnt += 1

            yield line, read_name, chr_1, int(pos_1), chr_2, int(
                pos_2
            ), mapq_1, mapq_2, tag_a, tag_b


def has_map_q_and_multi_map(in_filename, test, chr_filter):
    map_q = False
    multi_map = False
    last_map_q = None
    last_read_name = None
    for _, read_name, _, _, _, _, mapq_1, mapq_2, tag_a, tag_b in parse_heatmap(
        in_filename, test, chr_filter
    ):
        if last_read_name is None:
            last_read_name = read_name
        elif last_read_name == read_name:
            multi_map = True
        if len(tag_a) >= 5 and tag_b != "notag":
            multi_map = True
        if len(tag_b) >= 5 and tag_b != "notag":
            multi_map = True

        if last_map_q is None:
            last_map_q = mapq_1
            if mapq_1 != mapq_2:
                map_q = True
        elif last_map_q != mapq_1:
            map_q = True
        elif last_map_q != mapq_2:
            map_q = True

        if map_q and multi_map:
            break
    return map_q, multi_map


def group_heatmap(in_filename, file_size, chr_filter, no_groups=False, test=False):
    file_name = simplified_filepath(in_filename)
    curr_read_name = None
    group_1 = []
    group_2 = []

    def deal_with_group():
        nonlocal group_1
        nonlocal group_2
        do_cont = True
        chr_1_cmp = group_1[0][0]
        for chr_1, _1, _2 in group_1:
            if chr_1_cmp != chr_1:
                do_cont = False  # no reads that come from different chromosomes
        chr_2_cmp = group_2[0][0]
        for chr_2, _1, _2 in group_2:
            if chr_2_cmp != chr_2:
                do_cont = False  # no reads that come from different chromosomes
        if do_cont:
            if no_groups:
                pos_1_s = group_1[1]
                pos_1_e = group_1[1]
                pos_2_s = group_2[1]
                pos_2_e = group_2[1]
            else:
                pos_1_s = min([p for _1, p, _2 in group_1])
                pos_1_e = max([p for _1, p, _2 in group_1])
                pos_2_s = min([p for _1, p, _2 in group_2])
                pos_2_e = max([p for _1, p, _2 in group_2])
            map_q = min(
                max(x for _1, _2, x in group_1), max(x for _1, _2, x in group_2)
            )
            if len(group_1) > 1 and len(group_2) > 1:
                map_q += 1
            yield curr_read_name, chr_1, pos_1_s, pos_1_e, chr_2, pos_2_s, pos_2_e, map_q
        group_1 = []
        group_2 = []

    for idx_2, (
        _,
        read_name,
        chr_1,
        pos_1,
        chr_2,
        pos_2,
        mapq_1,
        mapq_2,
        tag_1,
        tag_2,
    ) in enumerate(parse_heatmap(in_filename, test, chr_filter)):
        if read_name != curr_read_name and len(group_1) > 0:
            yield from deal_with_group()
        curr_read_name = read_name
        if tag_1 == "notag":
            have_no_tag_1 = True
        if tag_2 == "notag":
            have_no_tag_1 = True
        group_1.append((chr_1, int(pos_1), int(mapq_1)))
        group_2.append((chr_2, int(pos_2), int(mapq_2)))
        for chr_1, pos_1 in read_xa_tag(tag_1):
            group_1.append((chr_1, int(pos_1), 0))
        for chr_2, pos_2 in read_xa_tag(tag_2):
            group_2.append((chr_2, int(pos_2), 0))

        if idx_2 % PRINT_MODULO == 0:
            print(
                "loading file",
                file_name,
                ", line",
                idx_2 + 1,
                "of",
                file_size,
                "=",
                round(100 * (idx_2 + 1) / file_size, 2),
                "%",
                end="\033[K\r",
                flush=True,
            )
    yield from deal_with_group()


class ChrOrderHeatmapIterator:
    def __init__(self, chrs, prefix):
        self.chrs = chrs
        self.prefix = prefix

    def cleanup(self):
        for chr_1 in self.chrs:
            for chr_2 in self.chrs:
                os.remove(self.prefix + "." + chr_1 + "." + chr_2)

    def itr_x_axis(self):
        for x in set(self.chrs.keys()):
            yield x

    def itr_y_axis(self):
        for y in set([chr_y for vals in self.chrs.values() for chr_y in vals]):
            yield y

    def itr_heatmap(self):
        for chr_x in self.yield_x_axis():
            for chr_y in self.yield_y_axis():
                yield chr_x, chr_y

    def itr_cell(self, chr_x, chr_y):
        with open(self.prefix + "." + chr_x + "." + chr_y, "r") as in_file:
            for line in in_file:
                yield line.split()

    def itr_row(self, chr_y):
        for chr_x in self.itr_x_axis():
            yield from self.itr_cell(chr_x, chr_y)

    def itr_col(self, chr_x):
        for chr_y in self.itr_y_axis():
            yield from self.itr_cell(chr_x, chr_y)

    def itr_diag(self):
        for x in set([*self.itr_x_axis(), *self.itr_y_axis()]):
            yield x


def chr_order_heatmap(
    index_prefix,
    dataset_name,
    in_filename,
    file_size,
    chr_filter,
    no_groups=False,
    test=False,
):
    prefix = index_prefix + "/.tmp." + dataset_name
    out_files = {}
    chrs = {}
    for (
        read_name,
        chr_1,
        pos_1_s,
        pos_1_e,
        chr_2,
        pos_2_s,
        pos_2_e,
        map_q,
    ) in group_heatmap(in_filename, file_size, chr_filter, no_groups, test):
        if chr_1 not in out_files:
            out_files[chr_1] = {}
            chrs[chr_1] = set()
        if chr_2 not in out_files[chr_1]:
            out_files[chr_1][chr_2] = open(prefix + "." + chr_1 + "." + chr_2, "x")
            chrs[chr_1].add(chr_2)

        out_files[chr_1][chr_2].write(
            "\t".join(
                [
                    read_name,
                    str(pos_1_s),
                    str(pos_1_e),
                    str(pos_2_s),
                    str(pos_2_e),
                    str(map_q),
                ]
            )
            + "\n"
        )

    for vals in out_files.values():
        for val in vals.values():
            val.close()

    return ChrOrderHeatmapIterator(chrs, prefix)


def get_filesize(path):
    if not os.path.exists(path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)
    return int(
        subprocess.run(["wc", "-l", path], stdout=subprocess.PIPE)
        .stdout.decode("utf-8")
        .split(" ")[0]
    )


def parse_annotations(annotation_file):
    with open(annotation_file, "r") as in_file_1:
        for line in in_file_1:
            if line[0] == "#":
                continue
            # parse file colum
            (
                chrom,
                db_name,
                annotation_type,
                from_pos,
                to_pos,
                _,
                strand,
                _,
                extras,
                *opt,
            ) = line.split("\t")
            yield annotation_type, chrom, int(from_pos), int(to_pos), extras.replace(
                ";", "\n"
            ).replace("%2C", ","), strand == '+'
