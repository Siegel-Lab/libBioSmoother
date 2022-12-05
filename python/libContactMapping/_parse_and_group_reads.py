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


def parse_tsv(in_filename, test, chr_filter, line_format):
    with open(in_filename, "r") as in_file_1:
        cnt = 0
        for line in in_file_1:
            # parse file columns
            num_cols = len(line.split())
            if num_cols in line_format:
                read_name, chrs, poss, mapqs, tags = line_format[num_cols](*line.split())

                cont = False
                for chr_ in chrs:
                    if not chr_ in chr_filter:
                        cont = True
                if cont:
                    continue
                mapqs = [0 if mapq in ["", "nomapq", "255", "*"] else mapq for mapq in mapqs]
                poss = [int(x)-1 for x in poss]
                mapqs = [int(x) for x in mapqs]

                if cnt > TEST_FAC and test:
                    break
                cnt += 1

                yield line, read_name, chrs, poss, mapqs, tags
            else:
                raise ValueError(
                    'line "'
                    + line
                    + '" has '
                    + str(num_cols)
                    + ", columns which is unexpected. There can be ["
                    + ", ".join(str(x) for x in line_format.keys())
                    + "] columns."
                )


def parse_heatmap(in_filename, test, chr_filter):
    yield from parse_tsv(in_filename, test, chr_filter, {
        7: lambda read_name, chr_1, pos_1, chr_2, pos_2, mapq_1, mapq_2: 
                (read_name, [chr_1, chr_2], [pos_1, pos_2], [mapq_1, mapq_2], ["?", "?"]),
        9: lambda read_name, chr_1, pos_1, chr_2, pos_2, mapq_1, mapq_2, tag_a, tag_b: 
                (read_name, [chr_1, chr_2], [pos_1, pos_2], [mapq_1, mapq_2], [tag_a, tag_b]),
        11: lambda read_name, _1, chr_1, pos_1, _2, _3, chr_2, pos_2, _4, mapq_1, mapq_2: 
                (read_name, [chr_1, chr_2], [pos_1, pos_2], [mapq_1, mapq_2], ["?", "?"]),
        13: lambda read_name, _1, chr_1, pos_1, _2, _3, chr_2, pos_2, _4, mapq_1, mapq_2, tag_a, tag_b: 
                (read_name, [chr_1, chr_2], [pos_1, pos_2], [mapq_1, mapq_2], [tag_a, tag_b]),
    })

def parse_track(in_filename, test, chr_filter):
    yield from parse_tsv(in_filename, test, chr_filter, {
        4: lambda read_name, chr_1, pos_1, mapq_1: 
                (read_name, [chr_1], [pos_1], [mapq_1], ["?"]),
        5: lambda read_name, chr_1, pos_1, mapq_1, tag_a: 
                (read_name, [chr_1], [pos_1], [mapq_1], [tag_a]),
        6: lambda read_name, _1, chr_1, pos_1, _2, mapq_1: 
                (read_name, [chr_1], [pos_1], [mapq_1], ["?"]),
        7: lambda read_name, _1, chr_1, pos_1, _2, mapq_1, tag_a: 
                (read_name, [chr_1], [pos_1], [mapq_1], [tag_a]),
    })


def has_map_q_and_multi_map(in_filename, test, chr_filter, parse_func=parse_heatmap):
    map_q = False
    multi_map = False
    last_map_q = None
    last_read_name = None
    for _, read_name, _, _, mapqs, tags in parse_func(in_filename, test, chr_filter):
        if last_read_name is None:
            last_read_name = read_name
        elif last_read_name == read_name:
            multi_map = True
        for tag in tags:
            if len(tag) >= 5 and tag != "notag":
                multi_map = True

        if last_map_q is None:
            last_map_q = mapqs[0]
        for mapq in mapqs:
            if mapq != last_map_q:
                map_q = True

        if map_q and multi_map:
            break
    return map_q, multi_map


def group_reads(in_filename, file_size, chr_filter, parse_func=parse_heatmap, no_groups=False, test=False):
    file_name = simplified_filepath(in_filename)
    curr_read_name = None
    group = []

    def deal_with_group():
        nonlocal group
        do_cont = True
        chrs = []
        for g in group:
            chr_1_cmp = g[0][0]
            for chr_, _1, _2 in g:
                if chr_1_cmp != chr_:
                    do_cont = False # no reads that come from different chromosomes
            chrs.append(chr_1_cmp)
        if do_cont:
            pos_s = []
            pos_e = []
            for g in group:
                if no_groups:
                    pos_s.append(g[0][1])
                    pos_e.append(g[0][1])
                else:
                    pos_s.append(min([p for _1, p, _2 in g]))
                    pos_e.append(max([p for _1, p, _2 in g]))
            map_q = min( [max(x for _1, _2, x in g ) for g in group] )
            if min(len(g) for g in group) > 1:
                map_q += 1
            yield curr_read_name, chrs, pos_s, pos_e, map_q
        group = []

    for idx_2, (
        _,
        read_name,
        chrs,
        poss,
        mapqs,
        tags,
    ) in enumerate(parse_func(in_filename, test, chr_filter)):
        if read_name != curr_read_name and len(group) > 0 and len(group[0]) > 0:
            yield from deal_with_group()
        curr_read_name = read_name
        for tag in tags:
            if tag == "notag":
                have_no_tag_1 = True
        for idx, (chr_, pos, mapq, tag) in enumerate(zip(chrs, poss, mapqs, tags)):
            if idx >= len(group):
                group.append([])
            group[idx].append((chr_, pos, mapq))
            for chr_1, pos_1 in read_xa_tag(tag):
                group[idx].append((chr_1, int(pos_1), 0))

        if idx_2 % PRINT_MODULO == PRINT_MODULO - 1:
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
        for chr_1 in set(self.chrs.keys()):
            for chr_2 in set(self.chrs[chr_1]):
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
        if chr_x in self.chrs and chr_y in self.chrs[chr_x]:
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
        chrs_,
        pos_s,
        pos_e,
        map_q,
    ) in group_reads(in_filename, file_size, chr_filter, parse_heatmap, no_groups, test):
        chr_1, chr_2 = chrs_
        if chr_1 not in out_files:
            out_files[chr_1] = set()
            chrs[chr_1] = set()
        if chr_2 not in out_files[chr_1]:
            out_files[chr_1].add(chr_2)
            chrs[chr_1].add(chr_2)
            ## clear file if exists
            if os.path.exists(prefix + "." + chr_1 + "." + chr_2):
                os.remove(prefix + "." + chr_1 + "." + chr_2)

        with open(prefix + "." + chr_1 + "." + chr_2, "a") as out_file:
            out_file.write(
                "\t".join(
                    [
                        read_name,
                        str(pos_s[0]),
                        str(pos_e[0]),
                        str(pos_s[1]),
                        str(pos_e[1]),
                        str(map_q),
                    ]
                )
                + "\n"
            )

    return ChrOrderHeatmapIterator(chrs, prefix)


class ChrOrderCoverageIterator:
    def __init__(self, chrs, prefix):
        self.chrs = chrs
        self.prefix = prefix

    def cleanup(self):
        for chr_ in self.chrs:
            os.remove(self.prefix + "." + chr_)

    def itr_x_axis(self):
        for x in self.chrs:
            yield x

    def itr_cell(self, chr_):
        with open(self.prefix + "." + chr_, "r") as in_file:
            for line in in_file:
                yield line.split()


def chr_order_coverage(
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
    chrs = set()
    for (
        read_name,
        chrs_,
        pos_s,
        pos_e,
        map_q,
    ) in group_reads(in_filename, file_size, chr_filter, parse_track, no_groups, test):
        if chrs_[0] not in chrs:
            ## clear file if exists
            if os.path.exists(prefix + "." + chrs_[0]):
                os.remove(prefix + "." + chrs_[0])
        chrs.add(chrs_[0])

        with open(prefix + "." + chrs_[0], "a") as out_file:
            out_file.write(
                "\t".join(
                    [
                        read_name,
                        str(pos_s[0]),
                        str(pos_e[0]),
                        str(map_q),
                    ]
                )
                + "\n"
            )

    return ChrOrderCoverageIterator(chrs, prefix)


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
