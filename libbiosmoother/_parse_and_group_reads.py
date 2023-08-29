PRINT_MODULO = 1000
import errno
import os
import fileinput
import sys
import gzip

TEST_FAC = 100000
MAX_READS_IM_MEM = 10000


def simplified_filepath(path):
    if "/" in path:
        path = path[path.rindex("/") + 1 :]
    if "." in path:
        return path[: path.index(".")]
    return path


def read_xa_tag(tags):
    tags.replace("XA:Z:", "")
    if len(tags) == 0 or tags == "notag":
        return []
    l = []
    for tag in tags.split(";"):
        if tag != "":
            split = tag.split(",")
            if len(split) == 4:
                chrom, str_pos, _CIGAR, _NM = split
                strand = str_pos[0]
                pos = int(str_pos[1:])
                l.append([chrom, pos, strand])
    return l


def __check_columns(columns, necessary, synonyms):
    columns = [col.lower() for col in columns]  # ignore capitalization

    # replace synonyms
    def get_synonym(col):
        if "[" in col and "]" in col:
            col = col.replace("[", "").replace("]", "")
            if col in synonyms:
                return "[" + synonyms[col] + "]"
            else:
                return "[" + col + "]"
        else:
            if col in synonyms:
                return synonyms[col]
            else:
                return col

    columns = [get_synonym(col) for col in columns]

    # check invalid optional columns
    for min_def in ["[" + col + "]" for col in necessary]:
        if min_def in columns:
            raise RuntimeError(
                ", ".join(necessary) + " cannot be given as optional columns."
            )

    # check necessary columns
    for min_def in necessary:
        if min_def not in columns:
            raise RuntimeError(
                "the given columns do not contain "
                + min_def
                + ". But this column is always necessary."
            )

    # check duplicates
    columns_non_optional = [col.replace("[", "").replace("]", "") for col in columns]
    cols_no_dot = [col for col in columns_non_optional if col != "."]
    if len(cols_no_dot) != len(set(cols_no_dot)):
        dup = []
        for idx, col in enumerate(cols_no_dot):
            if col in cols_no_dot[idx + 1 :]:
                dup.append(col)
        raise RuntimeError(
            "the given columns cannot contain duplicates. But "
            + ", ".join(dup)
            + ("is" if len(dup) == 1 else "are")
            + " given multiple times."
        )
    return columns


def setup_col_converter(
    columns, col_order, default_values, necessary_columns, synonyms
):
    have_printed_known_columns = False
    columns = __check_columns(columns, necessary_columns, synonyms)
    for col in columns:
        col = col.replace("[", "").replace("]", "")
        if col != "." and col not in col_order:
            print(
                "Warning: Unknown column name: '"
                + col
                + "'. This column will be ignored.",
                file=sys.stderr,
            )
            if not have_printed_known_columns:
                print( 
                    "Known column names are: "
                    + ", ".join("'" + c + "'" for c in col_order)
                    + ", and '.'.",
                    file=sys.stderr,
                )
                have_printed_known_columns = True
    non_opt_cols = sum(1 if "[" not in col and "]" not in col else 0 for col in columns)
    col_converter = {}
    for n in range(non_opt_cols, len(columns) + 1):
        dropped_cols = columns[:]
        idx = len(dropped_cols) - 1
        while len(dropped_cols) > n:
            assert idx >= 0
            if "[" in dropped_cols[idx] and "]" in dropped_cols[idx]:
                del dropped_cols[idx]
            idx -= 1
        dropped_cols = [col.replace("[", "").replace("]", "") for col in dropped_cols]

        col_converter[n] = [
            col_order.index(col_name)
            if col_name != "." and col_name in col_order
            else None
            for col_name in dropped_cols
        ]

    def convert(cols):
        n = min(len(cols), len(columns))
        if n not in col_converter:
            raise RuntimeError(
                "line '"
                + " ".join(cols)
                + "' does not match the expected columns:"
                + ", ".join(columns)
            )
        ret = default_values[:]
        for idx, col in zip(col_converter[n], cols):
            if not idx is None:
                ret[idx] = col
        return ret

    return convert


def parse_tsv(
    in_filename, test, chr_filter, make_line_format, default_cols, progress_print=print
):
    line_format = make_line_format(default_cols)
    with fileinput.input(in_filename) as in_file_1:
        cnt = 0
        file_pos = 0
        if in_filename != "-":
            file_size = get_filesize(in_filename)
        for idx_2, line in enumerate(in_file_1):
            file_pos += len(line)
            if idx_2 % PRINT_MODULO == PRINT_MODULO - 1:
                if in_filename != "-":
                    progress_print(
                        "file",
                        in_filename + ":",
                        str(round(100 * (file_pos) / file_size, 2)) + "%",
                    )
                else:
                    progress_print("from stdin: read " + str(idx_2) + " lines so far.")
            # ignore empty lines and comments / header lines
            if len(line) == 0:
                continue
            if line[:9] == "#columns:":
                line_format = make_line_format(line[9:].strip().split())
                # next line of code will make sure that this line of the input file is not read as actual data
            if line[0] == "#":
                continue

            # parse file columns
            read_name, chrs, poss, mapqs, tags, strand, bin_cnt = line_format(
                line.split("\t") if "\t" in line else line.split(" ")
            )

            cont = False
            for chr_ in chrs:
                if not chr_ in chr_filter:
                    cont = True
            if cont:
                continue
            mapqs = [
                0 if mapq in ["", "nomapq", "255", "*", "."] else mapq for mapq in mapqs
            ]
            poss = [max(0, int(x)) for x in poss]
            mapqs = [max(0, int(x)) for x in mapqs]
            for s in strand:
                if s not in ["+", "-"]:
                    raise RuntimeError("Invalid strand: " + s + "in line: " + line)
            strand = [s == "+" for s in strand]
            bin_cnt = int(bin_cnt)

            cnt += 1

            yield line, read_name, chrs, poss, mapqs, tags, strand, bin_cnt


def parse_heatmap(
    in_filename,
    test,
    chr_filter,
    progress_print=print,
    columns=["chr1", "pos1", "chr2", "pos2"],
):
    def make_converter(columns_in):
        col_converter = setup_col_converter(
            columns_in,
            [
                "readid",
                "chr1",
                "pos1",
                "chr2",
                "pos2",
                "strand1",
                "strand2",
                "mapq1",
                "mapq2",
                "xa1",
                "xa2",
                "cnt",
            ],
            ["-", ".", "0", ".", "0", "+", "+", "*", "*", "", "", "1"],
            ["chr1", "pos1", "chr2", "pos2"],
            {
                "chrom1": "chr1",
                "chrom2": "chr2",
                "readname": "readid",
                "count": "cnt",
            },
        )

        def convert(cols):
            (
                readid,
                chr1,
                pos1,
                chr2,
                pos2,
                strand1,
                strand2,
                mapq1,
                mapq2,
                xa1,
                xa2,
                cnt,
            ) = col_converter(cols)
            return (
                readid,
                [chr1, chr2],
                [pos1, pos2],
                [mapq1, mapq2],
                [xa1, xa2],
                [strand1, strand2],
                cnt,
            )

        return convert

    yield from parse_tsv(
        in_filename,
        test,
        chr_filter,
        make_converter,
        columns,
        progress_print=progress_print,
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


def parse_track(
    in_filename, test, chr_filter, progress_print=print, columns=["chr", "pos"]
):
    def make_converter(columns_in):
        col_converter = setup_col_converter(
            columns_in,
            ["readid", "chr", "pos", "strand", "mapq", "xa", "cnt"],
            ["-", ".", "0", "+", "*", "", "1"],
            ["chr", "pos"],
            {
                "chrom1": "chr1",
                "chrom2": "chr2",
                "readname": "readid",
                "count": "cnt",
            },
        )

        def convert(cols):
            readid, chr1, pos1, strand1, mapq1, xa1, cnt = col_converter(cols)
            return readid, [chr1], [pos1], [mapq1], [xa1], [strand1], cnt

        return convert

    yield from parse_tsv(
        in_filename,
        test,
        chr_filter,
        make_converter,
        columns,
        progress_print=progress_print,
    )


def group_reads(
    in_filename,
    file_size,
    chr_filter,
    progress_print=print,
    parse_func=parse_heatmap,
    no_groups=False,
    test=False,
    columns=["chr1", "pos1", "chr2", "pos2"],
):
    curr_read_name = None
    curr_count = None
    group = []

    def deal_with_group():
        nonlocal group
        dont_cont = False
        chrs = []
        for g in group:
            chr_1_cmp = g[0][0]
            for chr_, _1, _2, _3 in g:
                if chr_1_cmp != chr_:
                    dont_cont = True  # no reads that come from different chromosomes
            chrs.append(chr_1_cmp)

        pos_l = []
        pos_s = []
        pos_e = []
        strands = []
        for g in group:
            strands.append(g[0][3])
            if no_groups or dont_cont:
                pos_l.append([g[0][1]])
                pos_s.append(pos_l[-1][0])
                pos_e.append(pos_l[-1][0])
            else:
                pos_l.append([p for _1, p, _2, _3 in g])
                pos_s.append(min(pos_l[-1]))
                pos_e.append(max(pos_l[-1]))
        map_q = min([max(x for _1, _2, x, _3 in g) for g in group])
        if min(len(g) for g in group) > 1 and not dont_cont:
            map_q += 1
        yield curr_read_name, chrs, pos_s, pos_e, pos_l, map_q, strands, curr_count
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
    ) in parse_func(in_filename, test, chr_filter, progress_print, columns):
        if (
            (curr_read_name in [None, ".", "", "-"] or read_name != curr_read_name)
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

    def itr_cell(self, chr_x, chr_y):
        if chr_x in self.chrs and chr_y in self.chrs[chr_x]:
            if self.in_file[chr_x][chr_y]:
                with gzip.open(
                    self.prefix + "." + chr_x + "." + chr_y, "rt"
                ) as in_file:
                    for line in in_file:
                        yield line.split()
            for tup in self.chrs[chr_x][chr_y]:
                yield tup


def __make_filename_save(prefix):
    for illegal_char in "#%&{}\\<>*?/$!'\":@+`|= ":
        prefix = prefix.replace(illegal_char, "_")
    return prefix


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
    columns=["chr1", "pos1", "chr2", "pos2"],
):
    prefix = index_prefix + "/.tmp." + __make_filename_save(dataset_name)
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
        pos_l,
        map_q,
        strands,
        cnt,
    ) in group_reads(
        in_filename,
        file_size,
        chr_filter,
        progress_print,
        parse_func,
        no_groups,
        test,
        columns,
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
                ",".join(str(x) for x in pos_l[0]),
                ",".join(str(x) for x in pos_l[1]),
                strands[0],
                strands[1],
                map_q,
                cnt,
            )
        )

        if len(chrs[chr_1][chr_2]) >= MAX_READS_IM_MEM:
            with gzip.open(
                prefix + "." + chr_1 + "." + chr_2,
                "at" if in_file[chr_1][chr_2] else "wt",
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

    def itr_cell(self, chr_):
        if chr_ in self.chrs:
            if self.in_file[chr_]:
                with gzip.open(self.prefix + "." + chr_, "rt") as in_file:
                    for line in in_file:
                        yield line.split()
            for tup in self.chrs[chr_]:
                yield tup


def chr_order_coverage(
    index_prefix,
    dataset_name,
    in_filename,
    file_size,
    chr_filter,
    no_groups=False,
    test=False,
    progress_print=print,
    columns=["chr", "pos"],
):
    prefix = index_prefix + "/.tmp." + __make_filename_save(dataset_name)
    chrs = {}
    in_file = {}
    for (
        read_name,
        chrs_,
        pos_s,
        pos_e,
        pos_l,
        map_q,
        strands,
        cnt,
    ) in group_reads(
        in_filename,
        file_size,
        chr_filter,
        progress_print,
        parse_track,
        no_groups,
        test,
        columns,
    ):
        if chrs_[0] not in chrs:
            chrs[chrs_[0]] = []
            in_file[chrs_[0]] = False
        chrs[chrs_[0]].append(
            (
                read_name,
                pos_s[0],
                pos_e[0],
                ",".join(str(x) for x in pos_l[0]),
                strands[0],
                map_q,
                cnt,
            )
        )

        if len(chrs[chrs_[0]]) >= MAX_READS_IM_MEM:
            with gzip.open(
                prefix + "." + chrs_[0], "at" if in_file[chrs_[0]] else "wt"
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


def parse_annotations(annotation_file):
    with fileinput.input(annotation_file) as in_file_1:
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
