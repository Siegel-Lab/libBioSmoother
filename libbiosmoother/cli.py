from ._import_lib_bio_smoother_cpp import (
    SPS_VERSION,
    LIB_BIO_SMOOTHER_CPP_VERSION,
    COMPILER_ID,
)
from importlib.metadata import version
import argparse
from .indexer import *
from .quarry import Quarry
from .export import export_tsv, export_png, export_svg
import os
import shutil
from .parameters import list_parameters, values_for_parameter, open_valid_json
from .test import test
from .benchmark_runtime import benchmark_runtime

try:
    import cProfile

    HAS_CPROFILE = True
except ImportError:
    HAS_CPROFILE = False
from datetime import datetime


def get_path(prefix):
    for possible in [prefix, prefix + ".smoother_index"]:
        if os.path.exists(possible) and os.path.isdir(possible):
            return possible
    raise RuntimeError("the given index", prefix, "does not exist.")


def init(args):
    path = args.index_prefix.replace(".smoother_index", "")
    Indexer(path, strict=True).create_session(
        args.chr_len,
        args.dividend,
        args.anno_path,
        args.test,
        args.filterable_annotations,
        args.map_q_thresholds,
    )


def reset(args):
    possible = get_path(args.index_prefix)
    with open_default_json() as default_settings_file:
        default_settings = json.load(default_settings_file)
    with open(possible + "/default_session.json", "r") as default_session_file:
        default_session = json.load(default_session_file)
    default_session["settings"] = default_settings
    with open(possible + "/session.json", "w") as session_file:
        json.dump(default_session, session_file)


def w_perf(func, args):
    if args.perf:
        if HAS_CPROFILE:
            with cProfile.Profile() as pr:
                func()
                print("Profile:")
                pr.print_stats("tottime")
                pr.dump_stats(
                    get_path(args.index_prefix)
                    + "/profile."
                    + datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
                )
        else:
            raise RuntimeError("cProfile is not installed. Cannot write profile.")
    else:
        func()


def repl(args):
    idx = Indexer(get_path(args.index_prefix))

    def run():
        idx.add_replicate(
            args.path,
            args.name,
            args.group,
            args.no_groups,
            args.keep_points,
            args.only_points,
            args.no_map_q,
            args.no_multi_map,
            args.no_cat,
            not args.strand,
            args.shekelyan,
            args.force_upper_triangle,
            args.columns,
        )
        if args.ploidy_file is not None:
            idx.set_ploidy_list(args.ploidy_file)

    w_perf(run, args)


def norm(args):
    idx = Indexer(get_path(args.index_prefix))

    def run():
        idx.add_normalization(
            args.path,
            args.name,
            args.group,
            args.no_groups,
            args.keep_points,
            args.only_points,
            args.no_map_q,
            args.no_multi_map,
            args.no_cat,
            not args.strand,
            args.shekelyan,
            args.columns,
        )

    w_perf(run, args)


def export_smoother(args):
    session = Quarry(get_path(args.index_prefix))

    def run():
        if args.export_prefix is not None:
            session.set_value(["settings", "export", "prefix"], args.export_prefix)
        if args.export_selection is not None:
            session.set_value(
                ["settings", "export", "selection"], args.export_selection
            )
        if args.export_size is not None:
            session.set_value(["settings", "export", "size", "val"], args.export_size)

        if args.export_format is not None:
            if "tsv" in args.export_format:
                export_tsv(session)
            if "svg" in args.export_format:
                export_svg(session)
            if "png" in args.export_format:
                export_png(session)
        else:
            if session.get_value(["settings", "export", "export_format"]) == "tsv":
                export_tsv(session)
            if session.get_value(["settings", "export", "export_format"]) == "svg":
                export_svg(session)
            if session.get_value(["settings", "export", "export_format"]) == "png":
                export_png(session)

    w_perf(run, args)


def set_smoother(args):
    possible = get_path(args.index_prefix)
    with open(possible + "/session.json", "r") as in_file:
        json_file = json.load(in_file)
        tmp = json_file
        keys = args.name.split(".")
        for key in keys[:-1]:
            tmp = tmp[key]
        if isinstance(tmp[keys[-1]], bool):
            tmp[keys[-1]] = bool(args.val)
        elif isinstance(tmp[keys[-1]], float):
            tmp[keys[-1]] = float(args.val)
        elif isinstance(tmp[keys[-1]], int):
            tmp[keys[-1]] = int(args.val)
        elif isinstance(tmp[keys[-1]], str):
            tmp[keys[-1]] = str(args.val)
        else:
            print("Error: can only set string, int, bool and float values.")
    with open(possible + "/session.json", "w") as out_file:
        json.dump(json_file, out_file)


def get_smoother(args):
    possible = get_path(args.index_prefix)
    with open(possible + "/session.json", "r") as in_file:
        json_file = json.load(in_file)
        tmp = json_file
        for key in args.name.split("."):
            tmp = tmp[key]
        print(tmp)


def info_smoother(args):
    with open_valid_json() as valid_file:
        valid_json = json.load(valid_file)
        possible = get_path(args.index_prefix)
        with open(possible + "/session.json", "r") as in_file:
            json_file = json.load(in_file)
            for p in list_parameters(json_file, valid_json):
                print(
                    ".".join(p),
                    values_for_parameter(p, json_file, valid_json),
                    sep="\t",
                )


def test_smoother(args):
    test(Quarry(get_path(args.index_prefix)), args.seed, args.skip_first)


def benchmark_runtime_smoother(args):
    benchmark_runtime(
        Quarry(get_path(args.index_prefix)), args.num_experiments, args.outfile
    )


def ploidy_smoother(args):
    Quarry(get_path(args.index_prefix)).set_ploidy_list(args.ploidy_file)


def add_parsers(main_parser):
    init_parser = main_parser.add_parser("init", help="Create a new index.")
    init_parser.add_argument(
        "index_prefix",
        help="Path where the index shall be saved. Note: a folder with multiple files will be created.",
    )
    init_parser.add_argument(
        "chr_len",
        help="Path to a file that contains the length (in nucleotides) of all chromosomes. The file shall contain 2 tab seperated columns columns: The chromosome names and their size in nucleotides. The order of chromosomes in this files will be used as the display order in the viewer.",
    )
    init_parser.add_argument(
        "anno_path",
        help="Path to a file that contains the annotations",
        nargs="?",
        default="",
    )
    init_parser.add_argument(
        "-f",
        "--filterable_annotations",
        nargs="*",
        type=str,
        default=["gene"],
        help="Pick the annotations that can be used as filters in smoother. (default: %(default)s)",
    )
    init_parser.add_argument(
        "-m",
        "--map_q_thresholds",
        nargs="*",
        type=int,
        default=[3, 30],
        help="Pick several thresholds, that can then be used to filter reads by mapping quality. (default: %(default)s)",
    )
    init_parser.add_argument(
        "-d",
        "--dividend",
        type=int,
        default=10000,
        help="Divide all coordinates by this number. Larger numbers will reduce the index size and preprocessing time. However, bins with a size below this given number cannot be displayed. (default: %(default)s)",
    )
    init_parser.set_defaults(func=init)
    init_parser.add_argument("--test", help=argparse.SUPPRESS, action="store_true")

    reset_parser = main_parser.add_parser("reset", help="Reset session of an index.")
    reset_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    reset_parser.set_defaults(func=reset)

    repl_parser = main_parser.add_parser(
        "repl", help="Add a replicate to a given index."
    )
    repl_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    repl_parser.add_argument(
        "path", help="Path to the file that contains the aligned reads."
    )
    repl_parser.add_argument("name", help="Name for the new replicate.")
    repl_parser.add_argument(
        "-g",
        "--group",
        default="a",
        choices=["a", "b", "both", "neither"],
        help="Which analysis group to place the new replicate in when opening the interface. (default: %(default)s)",
    )
    repl_parser.add_argument(
        "-q",
        "--no_map_q",
        action="store_true",
        help="Do not store mapping quality information. This will make the index faster and smaller. (default: off)",
    )
    repl_parser.add_argument(
        "-m",
        "--no_multi_map",
        action="store_true",
        help="Do not multi mapping information (reads that map to multiple loci). This will make the index faster and smaller. (default: off)",
    )
    repl_parser.add_argument(
        "-c",
        "--no_cat",
        action="store_true",
        help="Do not store category information. (default: off)",
    )
    repl_parser.add_argument(
        "-s",
        "--strand",
        action="store_true",
        help="Do store strand information. (default: off)",
    )
    repl_parser.add_argument(
        "-C",
        "--columns",
        nargs="*",
        default=[
            "[readID]",
            "chr1",
            "pos1",
            "chr2",
            "pos2",
            "[strand1]",
            "[strand2]",
            "[.]",
            "[mapq1]",
            "[mapq2]",
            "[xa1]",
            "[xa2]",
            "[cnt]",
        ],
        # type=str,
        help="How columns of the input file should be interpreted. Valid column names are: 'readId', 'chr1', 'chr2', 'pos1', 'pos2', 'strand1', 'strand2', 'mapq1', 'mapq2', 'xa1', 'xa2', 'cnt', and '.'. At a minimum 'chr1', 'chr2', 'pos1' and 'pos2' need to be defined. Column names in squared brackets indicate optional columns (e.g. '[mapq1] [mapq2]'). Specify the columns as a space separated list: '-C chr1 pos1 [xa1]'. Optional columns can be omitted in the input file on a row-by-row basis. If there are multiple optional columns, they are ignored starting from the back. If '.' or '[.]' is given as a column name, that column of the input file will be ignored. If a row in the input file contains more than the given columns, the superfluous columns will be ignored. Lines in the input file that start with '#columns:' will also change this parameter for any following lines. (default: %(default)s)",
    )
    repl_parser.add_argument(
        "-u",
        "--force_upper_triangle",
        action="store_true",
        help="Mirror all interactions to the upper triangle. (default: off)",
    )
    repl_parser.set_defaults(func=repl)
    repl_parser.add_argument("--perf", help=argparse.SUPPRESS, action="store_true")
    repl_parser.add_argument(
        "--keep_points", help=argparse.SUPPRESS, action="store_true"
    )
    repl_parser.add_argument(
        "--only_points", help=argparse.SUPPRESS, action="store_true"
    )
    repl_parser.add_argument(
        "--ploidy_file", help="File that specifies the ploidy for each chromosome."
    )
    repl_parser.add_argument("--shekelyan", help=argparse.SUPPRESS, action="store_true")
    repl_parser.add_argument("--no_groups", help=argparse.SUPPRESS, action="store_true")

    norm_parser = main_parser.add_parser(
        "track",
        help="Add a normalization track to an index, using external sequencing data.",
    )
    norm_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    norm_parser.add_argument(
        "path", help="Path to the file that contains the aligned reads."
    )
    norm_parser.add_argument("name", help="Name for the new normalization track.")
    norm_parser.add_argument(
        "-g",
        "--group",
        default="neither",
        choices=["row", "col", "both", "neither"],
        help="Where to to place the new normalization track when opening the interface. (default: %(default)s)",
    )
    norm_parser.set_defaults(func=norm)
    norm_parser.add_argument(
        "--keep_points", help=argparse.SUPPRESS, action="store_true"
    )
    norm_parser.add_argument(
        "--only_points", help=argparse.SUPPRESS, action="store_true"
    )
    norm_parser.add_argument(
        "-q",
        "--no_map_q",
        action="store_true",
        help="Do not store mapping quality information. This will make the index faster and smaller. (default: off)",
    )
    norm_parser.add_argument(
        "-m",
        "--no_multi_map",
        action="store_true",
        help="Do not multi mapping information (reads that map to multiple loci). This will make the index faster and smaller. (default: off)",
    )
    norm_parser.add_argument(
        "-c",
        "--no_cat",
        action="store_true",
        help="Do not store category information. (default: off)",
    )
    norm_parser.add_argument(
        "-s",
        "--strand",
        action="store_true",
        help="Do store strand information. (default: off)",
    )
    norm_parser.add_argument(
        "-C",
        "--columns",
        nargs="*",
        default=["[readID]", "chr", "pos", "[strand]", "[mapq]", "[xa]", "[cnt]"],
        type=str,
        help="How columns of the input file should be interpreted. Valid column names are: 'readId', 'chr', 'pos', 'strand', 'mapq', 'xa', 'cnt', and '.'. At a minimum 'chr' and 'pos'. Column names in squared brackets indicate optional columns (e.g. '[mapq] [xa]'). Specify the columns as a space separated list: '-C chr1 pos1 [xa1]'. Optional columns can be omitted in the input file on a row-by-row basis. If there are multiple optional columns they are ignored starting from the back. If '.' or '[.]' is given as a column name, that column of the input file will be ignored. If a row in the input file contains more than the given columns, the superfluous columns will be ignored. Lines in the input file that start with '#columns:' will also change this parameter for any following lines. (default: %(default)s)",
    )
    norm_parser.add_argument("--shekelyan", help=argparse.SUPPRESS, action="store_true")
    norm_parser.add_argument("--no_groups", help=argparse.SUPPRESS, action="store_true")

    export_parser = main_parser.add_parser(
        "export", help="Export the current index session to a file."
    )
    export_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    export_parser.add_argument(
        "-p",
        "--export_prefix",
        help="Path where the exported file shall be saved",
    )
    export_parser.add_argument(
        "-f",
        "--export_format",
        choices=["tsv", "svg", "png"],
        nargs="*",
        help="The format which to export to.",
    )
    export_parser.add_argument(
        "-S",
        "--export_selection",
        choices=["heatmap", "col_sec_data", "row_sec_data"],
        nargs="*",
        help="Which regions shall be exported",
    )
    export_parser.add_argument(
        "-s",
        "--export_size",
        help="The size of the heatmap to be exported",
    )
    export_parser.set_defaults(func=export_smoother)

    set_parser = main_parser.add_parser("set", help="Set a parameter to a value.")
    set_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    set_parser.add_argument(
        "name",
        help="The name of the parameter to set.",
    )
    set_parser.add_argument(
        "val",
        help="The value to set.",
    )
    set_parser.set_defaults(func=set_smoother)

    get_parser = main_parser.add_parser("get", help="Get the value of a parameter.")
    get_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    get_parser.add_argument(
        "name",
        help="The name of the parameter to get.",
    )
    get_parser.set_defaults(func=get_smoother)

    info_parser = main_parser.add_parser("par_info", argument_default=argparse.SUPPRESS)
    info_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    info_parser.add_argument(
        "-n",
        "--name",
        help="The name of the parameter to get.",
    )
    info_parser.set_defaults(func=info_smoother)

    test_parser = main_parser.add_parser("test", argument_default=argparse.SUPPRESS)
    test_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    test_parser.add_argument(
        "-S", "--seed", help="Seed for random configurations.", default=None, type=int
    )
    test_parser.add_argument(
        "-s", "--skip_first", help="Skip the first x many test", default=0, type=int
    )
    test_parser.set_defaults(func=test_smoother)

    bench_parser = main_parser.add_parser(
        "benchmark", argument_default=argparse.SUPPRESS
    )
    bench_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    bench_parser.add_argument(
        "-N",
        "--num_experiments",
        help="Number of samples to take. (default: %(default)s)",
        default=100,
        type=int,
    )
    bench_parser.add_argument(
        "-o", "--outfile", help="outputfile", default="benchmark.pickle", type=str
    )
    bench_parser.set_defaults(func=benchmark_runtime_smoother)

    ploidy_parser = main_parser.add_parser("ploidy")
    ploidy_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    ploidy_parser.add_argument(
        "ploidy_file",
        help="File that specifies the ploidy for each chromosome.",
    )
    ploidy_parser.set_defaults(func=ploidy_smoother)


def make_main_parser():
    parser = argparse.ArgumentParser(description="")
    sub_parsers = parser.add_subparsers(
        help="Sub-command that shall be executed.", dest="cmd"
    )
    sub_parsers.required = True
    add_parsers(sub_parsers)
    return parser


def main():
    parser = make_main_parser()
    parser.add_argument(
        "-v", "--version", action="version", version=version("libbiosmoother")
    )
    parser.add_argument(
        "--version_smoother_cpp",
        help=argparse.SUPPRESS,
        action="version",
        version=LIB_BIO_SMOOTHER_CPP_VERSION,
    )
    parser.add_argument(
        "--version_sps", help=argparse.SUPPRESS, action="version", version=SPS_VERSION
    )
    parser.add_argument(
        "--compiler_id", action="version", help=argparse.SUPPRESS, version=COMPILER_ID
    )

    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
