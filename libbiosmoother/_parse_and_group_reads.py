PRINT_MODULO = 1000
import errno
import os

TEST_FAC = 100000
MAX_READS_IM_MEM = 10000


def simplified_filepath(path):
    if "/" in path:
        path = path[path.rindex("/") + 1 :]
    if "." in path:
        return path[: path.index(".")]
    return path


def read_xa_tag(tags):
    if tags == "notag" or len(tags) < 5:
        return []
    l = []
    for tag in tags[5:].split(";"):
        split = tag.split(",")
        if len(split) == 5:
            chrom, str_pos, _CIGAR, _NM = split
            strand = str_pos[0]
            pos = int(str_pos[1:])
            l.append([chrom, pos, strand])
    return l


def parse_tsv(in_filename, test, chr_filter, line_format, progress_print=print):
    with open(in_filename, "r") as in_file_1:
        cnt = 0
        file_pos = 0
        file_size = get_filesize(in_filename)
        for idx_2, line in enumerate(in_file_1):
            file_pos += len(line)
            if idx_2 % PRINT_MODULO == PRINT_MODULO - 1:
                progress_print(
                    "file",
                    in_filename + ":",
                    str(round(100 * (file_pos) / file_size, 2)) + "%",
                )
            # parse file columns
            num_cols = len(line.split())
            if num_cols in line_format:
                read_name, chrs, poss, mapqs, tags, strand, bin_cnt = line_format[num_cols](
                    *line.split()
                )

                cont = False
                for chr_ in chrs:
                    if not chr_ in chr_filter:
                        cont = True
                if cont:
                    continue
                mapqs = [
                    0 if mapq in ["", "nomapq", "255", "*"] else mapq for mapq in mapqs
                ]
                poss = [max(0, int(x)) for x in poss]
                mapqs = [max(0, int(x)) for x in mapqs]
                strand = [s == "+" for s in strand]
                bin_cnt = int(bin_cnt)

                # if cnt > TEST_FAC and test:
                #    break
                cnt += 1

                yield line, read_name, chrs, poss, mapqs, tags, strand, bin_cnt
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


def parse_heatmap(in_filename, test, chr_filter, progress_print=print):
    def helper_5_col(*cols):
        chr_1, pos_1, chr_2, pos_2, cnt = cols
        if isinstance(pos_1, int) and isinstance(pos_2, int) and isinstance(cnt, int):
            return (
                None,
                [chr_1, chr_2],
                [pos_1, pos_2],
                ["", ""],
                ["?", "?"],
                ["+", "+"],
                cnt,
            )
        read_name, chr_1, pos_1, chr_2, pos_2 = cols
        return (
            read_name,
            [chr_1, chr_2],
            [pos_1, pos_2],
            ["", ""],
            ["?", "?"],
            ["+", "+"],
            1,
        )

    yield from parse_tsv(
        in_filename,
        test,
        chr_filter,
        {
            5: helper_5_col,
            5: lambda read_name, chr_1, pos_1, chr_2, pos_2: (
                read_name,
                [chr_1, chr_2],
                [pos_1, pos_2],
                ["", ""],
                ["?", "?"],
                ["+", "+"],
                1,
            ),
            7: lambda read_name, chr_1, pos_1, chr_2, pos_2, mapq_1, mapq_2: (
                read_name,
                [chr_1, chr_2],
                [pos_1, pos_2],
                [mapq_1, mapq_2],
                ["?", "?"],
                ["+", "+"],
                1,
            ),
            9: lambda read_name, chr_1, pos_1, chr_2, pos_2, mapq_1, mapq_2, tag_a, tag_b: (
                read_name,
                [chr_1, chr_2],
                [pos_1, pos_2],
                [mapq_1, mapq_2],
                [tag_a, tag_b],
                ["+", "+"],
                1,
            ),
            11: lambda read_name, str1, chr_1, pos_1, _2, str2, chr_2, pos_2, _4, mapq_1, mapq_2: (
                read_name,
                [chr_1, chr_2],
                [pos_1, pos_2],
                [mapq_1, mapq_2],
                ["?", "?"],
                [str1, str2],
                1,
            ),
            13: lambda read_name, str1, chr_1, pos_1, _2, str2, chr_2, pos_2, _4, mapq_1, mapq_2, tag_a, tag_b: (
                read_name,
                [chr_1, chr_2],
                [pos_1, pos_2],
                [mapq_1, mapq_2],
                [tag_a, tag_b],
                [str1, str2],
                1,
            ),
        },
        progress_print,
    )


def force_upper_triangle(
    in_filename, test, chr_filter, progress_print=print, parse_func=parse_heatmap
):
    for line, read_name, chrs, poss, mapqs, tags, strand, cnt in parse_func(
        in_filename, test, chr_filter, progress_print
    ):
        order = [
            (chr_filter.index(chrs[idx]), poss[idx], idx) for idx in range(len(chrs))
        ]

        chrs_out = []
        poss_out = []
        mapqs_out = []
        tags_out = []
        strand_out = []
        for _, _, idx in sorted(order):
            chrs_out.append(chrs[idx])
            poss_out.append(poss[idx])
            mapqs_out.append(mapqs[idx])
            tags_out.append(tags[idx])
            strand_out.append(strand[idx])

        yield line, read_name, chrs_out, poss_out, mapqs_out, tags_out, strand_out, cnt


def parse_track(in_filename, test, chr_filter, progress_print=print):
    yield from parse_tsv(
        in_filename,
        test,
        chr_filter,
        {
            4: lambda read_name, chr_1, pos_1, mapq_1: (
                read_name,
                [chr_1],
                [pos_1],
                [mapq_1],
                ["?"],
                ["+"],
            ),
            5: lambda read_name, chr_1, pos_1, mapq_1, tag_a: (
                read_name,
                [chr_1],
                [pos_1],
                [mapq_1],
                [tag_a],
                ["+"],
            ),
            6: lambda read_name, str1, chr_1, pos_1, _2, mapq_1: (
                read_name,
                [chr_1],
                [pos_1],
                [mapq_1],
                ["?"],
                [str1],
            ),
            7: lambda read_name, str1, chr_1, pos_1, _2, mapq_1, tag_a: (
                read_name,
                [chr_1],
                [pos_1],
                [mapq_1],
                [tag_a],
                [str1],
            ),
        },
        progress_print,
    )


def group_reads(
    in_filename,
    file_size,
    chr_filter,
    progress_print=print,
    parse_func=parse_heatmap,
    no_groups=False,
    test=False,
):
    curr_read_name = None
    curr_count = None
    group = []

    def deal_with_group():
        nonlocal group
        do_cont = True
        chrs = []
        for g in group:
            chr_1_cmp = g[0][0]
            for chr_, _1, _2, _3 in g:
                if chr_1_cmp != chr_:
                    do_cont = False  # no reads that come from different chromosomes
            chrs.append(chr_1_cmp)
        if do_cont:
            pos_s = []
            pos_e = []
            strands = []
            for g in group:
                strands.append(g[0][3])
                if no_groups:
                    pos_s.append(g[0][1])
                    pos_e.append(g[0][1])
                else:
                    pos_s.append(min([p for _1, p, _2, _3 in g]))
                    pos_e.append(max([p for _1, p, _2, _3 in g]))
            map_q = min([max(x for _1, _2, x, _3 in g) for g in group])
            if min(len(g) for g in group) > 1:
                map_q += 1
            yield curr_read_name, chrs, pos_s, pos_e, map_q, strands, curr_count
        group = []

    for (
        _,
        read_name,
        chrs,
        poss,
        mapqs,
        tags,
        strands,
        cnt,
    ) in parse_func(in_filename, test, chr_filter, progress_print):
        if (
            (curr_read_name is None or read_name != curr_read_name)
            and len(group) > 0
            and len(group[0]) > 0
        ):
            yield from deal_with_group()
        curr_read_name = read_name
        curr_count = cnt
        for idx, (chr_, pos, mapq, tag, strand) in enumerate(
            zip(chrs, poss, mapqs, tags, strands)
        ):
            if idx >= len(group):
                group.append([])
            group[idx].append((chr_, pos, mapq, strand))
            for chr_1, pos_1, str_1 in read_xa_tag(tag):
                group[idx].append((chr_1, int(pos_1), 0, str_1))
    if len(group) > 0 and len(group[0]) > 0:
        yield from deal_with_group()


class ChrOrderHeatmapIterator:
    def __init__(self, chrs, in_file, prefix):
        self.chrs = chrs
        self.in_file = in_file
        self.prefix = prefix

    def cleanup(self):
        for chr_1 in set(self.chrs.keys()):
            for chr_2 in set(self.chrs[chr_1].keys()):
                if self.in_file[chr_1][chr_2]:
                    os.remove(self.prefix + "." + chr_1 + "." + chr_2)

    def itr_x_axis(self):
        for x in set(self.chrs.keys()):
            yield x

    def itr_y_axis(self):
        for y in set([chr_y for vals in self.chrs.values() for chr_y in vals.keys()]):
            yield y

    def itr_heatmap(self):
        for chr_x in self.yield_x_axis():
            for chr_y in self.yield_y_axis():
                yield chr_x, chr_y

    def itr_cell(self, chr_x, chr_y):
        if chr_x in self.chrs and chr_y in self.chrs[chr_x]:
            if self.in_file[chr_x][chr_y]:
                with open(self.prefix + "." + chr_x + "." + chr_y, "r") as in_file:
                    for line in in_file:
                        yield line.split()
            for tup in self.chrs[chr_x][chr_y]:
                yield tup

    def itr_row(self, chr_y):
        for chr_x in self.itr_x_axis():
            yield from self.itr_cell(chr_x, chr_y)

    def itr_col(self, chr_x):
        for chr_y in self.itr_y_axis():
            yield from self.itr_cell(chr_x, chr_y)

    def itr_diag(self):
        for x in set([*self.itr_x_axis(), *self.itr_y_axis()]):
            yield x

    def __len__(self):
        x_cnt = len(list(self.itr_x_axis()))
        y_cnt = len(list(self.itr_y_axis()))
        return x_cnt * y_cnt


def chr_order_heatmap(
    index_prefix,
    dataset_name,
    in_filename,
    file_size,
    chr_filter,
    no_groups=False,
    test=False,
    do_force_upper_triangle=False,
    progress_print=print,
):
    prefix = index_prefix + "/.tmp." + dataset_name
    chrs = {}
    in_file = {}
    if do_force_upper_triangle:
        parse_func = force_upper_triangle
    else:
        parse_func = parse_heatmap
    for (
        read_name,
        chrs_,
        pos_s,
        pos_e,
        map_q,
        strands,
        cnt,
    ) in group_reads(
        in_filename, file_size, chr_filter, progress_print, parse_func, no_groups, test
    ):
        chr_1, chr_2 = chrs_
        if chr_1 not in chrs:
            chrs[chr_1] = {}
            in_file[chr_1] = {}
        if chr_2 not in chrs[chr_1]:
            chrs[chr_1][chr_2] = []
            in_file[chr_1][chr_2] = False
        chrs[chr_1][chr_2].append(
            (
                read_name,
                pos_s[0],
                pos_e[0],
                pos_s[1],
                pos_e[1],
                strands[0],
                strands[1],
                map_q,
                cnt,
            )
        )

        if len(chrs[chr_1][chr_2]) >= MAX_READS_IM_MEM:
            with open(
                prefix + "." + chr_1 + "." + chr_2,
                "a" if in_file[chr_1][chr_2] else "w",
            ) as out_file:
                for tup in chrs[chr_1][chr_2]:
                    out_file.write("\t".join([str(x) for x in tup]) + "\n")
            chrs[chr_1][chr_2] = []
            in_file[chr_1][chr_2] = True

    return ChrOrderHeatmapIterator(chrs, in_file, prefix)


class ChrOrderCoverageIterator:
    def __init__(self, chrs, in_file, prefix):
        self.chrs = chrs
        self.in_file = in_file
        self.prefix = prefix

    def cleanup(self):
        for chr_ in self.chrs.keys():
            if self.in_file[chr_]:
                os.remove(self.prefix + "." + chr_)

    def itr_x_axis(self):
        for x in self.chrs.keys():
            yield x

    def itr_cell(self, chr_):
        if self.in_file[chr_]:
            with open(self.prefix + "." + chr_, "r") as in_file:
                for line in in_file:
                    yield line.split()
        for tup in self.chrs[chr_]:
            yield tup

    def __len__(self):
        return len(self.chrs)


def chr_order_coverage(
    index_prefix,
    dataset_name,
    in_filename,
    file_size,
    chr_filter,
    no_groups=False,
    test=False,
    progress_print=print,
):
    prefix = index_prefix + "/.tmp." + dataset_name
    chrs = {}
    in_file = {}
    for (
        read_name,
        chrs_,
        pos_s,
        pos_e,
        map_q,
        strands,
    ) in group_reads(
        in_filename, file_size, chr_filter, progress_print, parse_track, no_groups, test
    ):
        if chrs_[0] not in chrs:
            chrs[chrs_[0]] = []
            in_file[chrs_[0]] = False
        chrs[chrs_[0]].append((read_name, pos_s[0], pos_e[0], strands[0], map_q))

        if len(chrs[chrs_[0]]) >= MAX_READS_IM_MEM:
            with open(
                prefix + "." + chrs_[0], "a" if in_file[chrs_[0]] else "w"
            ) as out_file:
                for tup in chrs[chrs_[0]]:
                    out_file.write("\t".join([str(x) for x in tup]) + "\n")
            chrs[chrs_[0]] = []
            in_file[chrs_[0]] = True

    return ChrOrderCoverageIterator(chrs, in_file, prefix)


def get_filesize(path):
    if not os.path.exists(path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)
    return os.path.getsize(path)
    # return int(
    #    subprocess.run(["wc", "-l", path], stdout=subprocess.PIPE)
    #    .stdout.decode("utf-8")
    #    .split(" ")[0]
    # )


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
            ).replace("%2C", ","), strand == "+"
