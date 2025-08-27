from ._import_lib_bio_smoother_cpp import (
    SPS_VERSION,
    LIB_BIO_SMOOTHER_CPP_VERSION,
    COMPILER_ID,
    WITH_STXXL,
    UNROLL_FOR_ALL_COMBINATIONS,
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

HIDE_SUBCOMMANDS_MANUAL = "SMOOTHER_HIDE_SUBCOMMANDS_MANUAL" in os.environ

REPL_C_DEFAULT = [
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
]
NORM_C_DEFAULTS = ["[readID]", "chr", "pos", "[strand]", "[mapq]", "[xa]", "[cnt]"]


def get_path(prefix):
    for possible in [prefix, prefix + ".smoother_index"]:
        if os.path.exists(possible) and os.path.isdir(possible):
            return possible
    raise RuntimeError("the given index", prefix, "does not exist.")


def init(args):
    path = args.index_prefix.replace(".smoother_index", "")
    if len(args.ctg_len) == 0 and len(args.anno_path) == 0:
        raise RuntimeError(
            "Either a chromosome length file or an annotation file must be given."
        )
    Indexer(path, strict=True).create_session(
        args.ctg_len,
        args.dividend,
        args.anno_path,
        args.test,
        args.filterable_annotations,
        args.map_q_thresholds,
    )
    if args.ploidy_file is not None:
        Quarry(get_path(args.index_prefix)).set_ploidy_list(args.ploidy_file, print)


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

    if args.columns is None:
        cols = REPL_C_DEFAULT
        allow_col_change = True
    else:
        cols = args.columns
        allow_col_change = False

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
            args.no_anno,
            not args.strand,
            args.shekelyan,
            args.force_upper_triangle,
            cols,
            allow_col_change,
        )

    w_perf(run, args)


def cool(args):
    idx = Indexer(get_path(args.index_prefix))

    def run():
        idx.add_cool(
            args.path,
            args.name,
            args.group,
            args.no_anno,
            args.shekelyan,
            args.force_upper_triangle,
        )

    w_perf(run, args)


def norm(args):
    idx = Indexer(get_path(args.index_prefix))

    if args.columns is None:
        cols = NORM_C_DEFAULTS
        allow_col_change = True
    else:
        cols = args.columns
        allow_col_change = False

    def run():
        idx.add_normalization(
            args.path,
            args.name,
            args.axis,
            args.no_groups,
            args.keep_points,
            args.only_points,
            args.no_map_q,
            args.no_multi_map,
            args.no_anno,
            not args.strand,
            args.shekelyan,
            cols,
            allow_col_change,
        )

    w_perf(run, args)


def export_smoother(args):
    session = Quarry(get_path(args.index_prefix))

    def run():
        if args.export_prefix is not None:
            session.set_value(["settings", "export", "prefix"], args.export_prefix)
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
    if args.name == "area":
        q = Quarry(possible)
        new_area = q.interpret_area(
            args.val,
            q.get_value(["area", "x_start"]),
            q.get_value(["area", "y_start"]),
            q.get_value(["area", "x_end"]),
            q.get_value(["area", "y_end"]),
            report_error=lambda s: print(s, file=sys.stderr),
        )
        if not new_area[0] is None:
            q.set_value(["area", "x_start"], new_area[0])
        if not new_area[1] is None:
            q.set_value(["area", "x_end"], new_area[1])
        if not new_area[2] is None:
            q.set_value(["area", "y_start"], new_area[2])
        if not new_area[3] is None:
            q.set_value(["area", "y_end"], new_area[3])
        q.save_session()
    elif args.name == "settings.interface.v4c.col":
        q = Quarry(possible)
        new_range = q.interpret_range(
            args.val, True, report_error=lambda s: print(s, file=sys.stderr)
        )
        if not new_range[0] is None:
            q.set_value(["settings", "interface", "v4c", "col_from"], new_range[0])
        if not new_range[1] is None:
            q.set_value(["settings", "interface", "v4c", "col_to"], new_range[1])
        q.save_session()
    elif args.name == "settings.interface.v4c.row":
        q = Quarry(possible)
        new_range = q.interpret_range(
            args.val, False, report_error=lambda s: print(s, file=sys.stderr)
        )
        if not new_range[0] is None:
            q.set_value(["settings", "interface", "v4c", "row_from"], new_range[0])
        if not new_range[1] is None:
            q.set_value(["settings", "interface", "v4c", "row_to"], new_range[1])
        q.save_session()
    elif args.name == ".":
        with open(possible + "/session.json", "w") as out_file:
            json.dump(json_file, json.loads(args.val))
    else:
        with open(possible + "/session.json", "r") as in_file:
            json_file = json.load(in_file)
            tmp = json_file
            keys = args.name.split(".")
            for key in keys[:-1]:
                tmp = tmp[key]
            if isinstance(tmp[keys[-1]], bool):
                if args.val.lower() not in [
                    "true",
                    "1",
                    "t",
                    "y",
                    "yes",
                    "on",
                    "false",
                    "0",
                    "f",
                    "n",
                    "no",
                    "off",
                    "blub",
                ]:
                    print("Error: can only set bool values to true or false.")
                    return
                tmp[keys[-1]] = args.val.lower() in ["true", "1", "t", "y", "yes", "on"]
            elif isinstance(tmp[keys[-1]], float):
                tmp[keys[-1]] = float(args.val)
            elif isinstance(tmp[keys[-1]], int):
                tmp[keys[-1]] = int(args.val)
            elif isinstance(tmp[keys[-1]], str):
                tmp[keys[-1]] = str(args.val)
            elif isinstance(tmp[keys[-1]], list) and (
                len(tmp[keys[-1]]) == 0 or isinstance(tmp[keys[-1]][0], str)
            ):
                tmp[keys[-1]] = str(args.val).split()
            else:
                print(
                    "Error: can only set string, int, bool, float, and list of string values.",
                    key,
                    "is neither of those.",
                )
        with open(possible + "/session.json", "w") as out_file:
            json.dump(json_file, out_file)


def get_smoother(args):
    possible = get_path(args.index_prefix)
    if args.name == "area":
        q = Quarry(possible)
        print(
            q.get_readable_area(
                q.get_value(["area", "x_start"]),
                q.get_value(["area", "y_start"]),
                q.get_value(["area", "x_end"]),
                q.get_value(["area", "y_end"]),
            )
        )
    elif args.name == "settings.interface.v4c.col":
        q = Quarry(possible)
        print(
            q.get_readable_range(
                q.get_value(["settings", "interface", "v4c", "col_from"]),
                q.get_value(["settings", "interface", "v4c", "col_to"]),
                True,
            )
        )
    elif args.name == "settings.interface.v4c.row":
        q = Quarry(possible)
        print(
            q.get_readable_range(
                q.get_value(["settings", "interface", "v4c", "row_from"]),
                q.get_value(["settings", "interface", "v4c", "row_to"]),
                False,
            )
        )
    elif args.name == ".":
        with open(possible + "/session.json", "r") as in_file:
            json_file = json.load(in_file)
            print(json_file)
    else:
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
    Quarry(get_path(args.index_prefix)).set_ploidy_list(args.ploidy_file, report_error=print)


def add_parsers(main_parser):
    def fmt_defaults(defaults):
        return " (default: " + " ".join(str(x) for x in defaults) + ")"

    init_parser = main_parser.add_parser(
        "init",
        help="Generate a new index for a given reference genome.",
        description="Either --anno_path or --ctg_len must be given.",
    )
    init_parser.add_argument(
        "index_prefix",
        help="Path where the index directory will be saved. Note: a folder with multiple files will be created.",
    )
    init_parser.add_argument(
        "-c",
        "--ctg_len",
        help="Path to a 2-column tab separated file containing the contig names and their size in basepairs. The order of contigs in this file will determine the order they are displayed in. If no --anno_path is given, this parameter becomes required.",
        nargs="?",
        default="",
    )
    init_parser.add_argument(
        "-a",
        "--anno_path",
        help="Path to a GFF file containing the annotations of the reference genome. If --ctg_len is not given this parameter becomes required and the order of contigs will be inferred from the order of '##sequence-region' lines in this GFF file.",
        nargs="?",
        default="",
    )
    defaults = ["gene"]
    init_parser.add_argument(
        "-f",
        "--filterable_annotations",
        nargs="*",
        type=str,
        default=defaults,
        help="List the annotations that can be used as filters. The annotations listed must be present in the 'anno_path' GFF file."
        + fmt_defaults(defaults),
    )
    defaults = [3, 30]
    init_parser.add_argument(
        "-m",
        "--map_q_thresholds",
        nargs="*",
        type=int,
        default=defaults,
        help="List several thresholds that can be used to filter reads by mapping quality score."
        + fmt_defaults(defaults),
    )
    init_parser.add_argument(
        "-d",
        "--dividend",
        type=int,
        default=10000,
        help="Divide all coordinates by this number, this corresponds to the minimal bin size that can be displayed. Larger numbers will reduce the index size and pre-processing time. (default: %(default)s)",
    )
    init_parser.add_argument(
        "--ploidy_file", help="File that specifies the ploidy for each chromosome."
    )
    init_parser.set_defaults(func=init)
    init_parser.add_argument("--test", help=argparse.SUPPRESS, action="store_true")

    reset_parser = main_parser.add_parser(
        "reset", help="Reset a given index with to the default parameters."
    )
    reset_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    reset_parser.set_defaults(func=reset)

    repl_parser = main_parser.add_parser(
        "repl", help="Add data for a sample or replicate to an index."
    )
    repl_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    repl_parser.add_argument(
        "path", help="Path to the input pairs file containing the interactions."
    )
    repl_parser.add_argument("name", help="Name for the new replicate or sample.")
    repl_parser.add_argument(
        "-g",
        "--group",
        default="a",
        choices=["a", "b", "both", "neither"],
        help="Analysis group or condition for the new replicate. This can also be modified in the GUI. Options are: 'a', 'b', 'both', 'neither'. (default: %(default)s)",
    )
    repl_parser.add_argument(
        "-q",
        "--no_map_q",
        action="store_true",
        help="Do not store mapping quality information for the replicate/sample. This will make the index smaller. (default: off)",
    )
    repl_parser.add_argument(
        "-m",
        "--no_multi_map",
        action="store_true",
        help="Exclude multi mapping reads from the analysis. For reads mapping to multiple loci only the main mapping position will be kept, and the secondary alignments will be ignored. This will make the index smaller. Note that if multiple alignments are not marked as one primary and other secondary alignments (I.e., they are kept in different lines in the input pairs file), only the first line with a given readID will be kept (default: off)",
    )
    repl_parser.add_argument(
        "-a",
        "--no_anno",
        action="store_true",
        help="Do not store annotation information. (default: off)",
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
        default=None,
        # type=str,
        help="Define the order of columns of the input pairs file. Valid column names are: 'readId', 'chr1', 'chr2', 'pos1', 'pos2', 'strand1', 'strand2', 'mapq1', 'mapq2', 'xa1', 'xa2', 'cnt', 'pair_type', and '.'. Specify the columns as a space separated list: '-C chr1 pos1 [xa1]'. Column names in squared brackets indicate optional columns (e.g. '[mapq1] [mapq2]'). Optional columns must not appear in all lines of the input file. 'chr1', 'chr2', 'pos1' and 'pos2' must be defined and cannot be optional columns. If a row in the input file has less columns than defined, optional columns will be ignored starting from the end. If a row in the input file contains more than the given columns, the superfluous columns will be ignored starting from the end. Columns defined as '.' or '[.]' will be ignored. Lines in the input file that start with '#columns:' will set this parameter for all following lines, if the parameter has not been specified on the command line."
        + fmt_defaults(REPL_C_DEFAULT),
    )
    repl_parser.add_argument(
        "-u",
        "--force_upper_triangle",
        action="store_true",
        help="Move all interactions to the upper triangle and keep the lower triangle empty. Enable this option for symmetric data (e.g., Hi-C) with a non-redundant and non-sorted input pairs file. If the input file has sorted pairs, or the data is asymmetric (e.g., RD-SPRITE) this option is not needed. (default: off)",
    )
    repl_parser.set_defaults(func=repl)
    repl_parser.add_argument("--perf", help=argparse.SUPPRESS, action="store_true")
    repl_parser.add_argument(
        "--keep_points", help=argparse.SUPPRESS, action="store_true"
    )
    repl_parser.add_argument(
        "--only_points", help=argparse.SUPPRESS, action="store_true"
    )
    repl_parser.add_argument("--shekelyan", help=argparse.SUPPRESS, action="store_true")
    repl_parser.add_argument("--no_groups", help=argparse.SUPPRESS, action="store_true")

    cool_parser = main_parser.add_parser(
        "cool",
        help="Add data for a sample or replicate to an index from a cooler file. Note: The information necessary to apply filters (e.g. the mapping quality of reads) is not stored in cool files, therefore filters will be disabled for datasets created from cool files.",
    )
    cool_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    cool_parser.add_argument(
        "path", help="Path to the input pairs file containing the interactions."
    )
    cool_parser.add_argument("name", help="Name for the new replicate or sample.")
    cool_parser.add_argument(
        "-g",
        "--group",
        default="a",
        choices=["a", "b", "both", "neither"],
        help="Analysis group or condition for the new replicate. This can also be modified in the GUI. Options are: 'a', 'b', 'both', 'neither'. (default: %(default)s)",
    )
    cool_parser.add_argument(
        "-a",
        "--no_anno",
        action="store_true",
        help="Do not store annotation information. (default: off)",
    )
    cool_parser.add_argument(
        "-u",
        "--force_upper_triangle",
        action="store_true",
        help="Move all interactions to the upper triangle and keep the lower triangle empty. Enable this option for symmetric data (e.g., Hi-C) with a non-redundant and non-sorted input pairs file. If the input file has sorted pairs, or the data is asymmetric (e.g., RD-SPRITE) this option is not needed. (default: off)",
    )
    cool_parser.set_defaults(func=cool)
    cool_parser.add_argument("--perf", help=argparse.SUPPRESS, action="store_true")
    cool_parser.add_argument("--shekelyan", help=argparse.SUPPRESS, action="store_true")

    norm_parser = main_parser.add_parser(
        "track",
        help="Add a track of uni-dimensional sequencing data that can be used for normalisation on Smoother to an index.",
    )
    norm_parser.add_argument("--perf", help=argparse.SUPPRESS, action="store_true")
    norm_parser.add_argument(
        "index_prefix",
        help="Prefix that was used to create the index (see the init subcommand).",
    )
    norm_parser.add_argument(
        "path", help="Path to the file that contains the aligned reads."
    )
    norm_parser.add_argument("name", help="Name for the new normalization track.")
    norm_parser.add_argument(
        "-x",
        "--axis",
        default="neither",
        choices=["row", "col", "both", "neither"],
        help="Axis in which the uni-dimensional track will be displayed. This can also be modified interactively. Options are: 'row', 'col', 'both', 'neither' (default: %(default)s)",
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
        help="Do not store mapping quality information for the replicate/sample. This will make the index smaller. (default: off)",
    )
    norm_parser.add_argument(
        "-m",
        "--no_multi_map",
        action="store_true",
        help="Exclude multi mapping reads from the analysis. For reads mapping to multiple loci only the main mapping position will be kept, and all secondary alignments will be ignored. This will make the index smaller. Note that if multiple alignments are not marked as one primary and other secondary alignments and are thus kept in different lines in the input pairs file, only the first line with a given readID will be kept. (default: off)",
    )
    norm_parser.add_argument(
        "-a",
        "--no_anno",
        action="store_true",
        help="Do not store annotation information. This will make the index smaller. (default: off)",
    )
    norm_parser.add_argument(
        "-s",
        "--strand",
        action="store_true",
        help="Do store strand information. This will make the index larger. (default: off)",
    )
    norm_parser.add_argument(
        "-C",
        "--columns",
        nargs="*",
        default=None,
        type=str,
        help="Define the order of columns of the input file. Valid column names are: 'readId', 'chr', 'pos', 'strand', 'mapq', 'xa', 'cnt', and '.'. Specify the columns as a space separated list: '-C chr pos [xa]'. Column names in squared brackets indicate optional columns (e.g. '[mapq] [xa]'). Optional columns must not appear in all lines of the input file. 'chr' and 'pos' must be defined and cannot be optional columns. If a row in the input file has less columns than defined, optional columns will be ignored starting from the end. If a row in the input file contains more than the given columns, the superfluous columns will be ignored starting from the end. Columns defined as '.' or '[.]' will be ignored. Lines in the input file that start with '#columns:' will set this parameter for all following lines, if the parameter has not been specified on the command line."
        + fmt_defaults(NORM_C_DEFAULTS),
    )
    norm_parser.add_argument("--shekelyan", help=argparse.SUPPRESS, action="store_true")
    norm_parser.add_argument("--no_groups", help=argparse.SUPPRESS, action="store_true")

    export_parser = main_parser.add_parser(
        "export", help="Export the current index session to a file."
    )
    export_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    export_parser.add_argument(
        "-p",
        "--export_prefix",
        help="Path for the exported files. The appropriate file extension will be added automatically.",
    )
    export_parser.add_argument(
        "-f",
        "--export_format",
        choices=["tsv", "svg", "png"],
        nargs="*",
        help="File format to be exported. Format options are 'tsv', 'svg', 'png'. When exporting in tsv format, 3 files are saved: one for the interactome, and one for each axis. Multiple formats can be exported at once: For example, using the argument '-f tsv png' will export both the tsv text file and a png picture. (Default: tsv)",
    )
    export_parser.add_argument(
        "-s",
        "--export_size",
        type=int,
        help="Size in pixels of the heatmap. (default: 800)",
    )
    export_parser.add_argument("--perf", help=argparse.SUPPRESS, action="store_true")
    export_parser.set_defaults(func=export_smoother)

    set_parser = main_parser.add_parser(
        "set",
        help="Set the values of different parameters in an index. This can also be done on in the graphical user interface.",
    )
    set_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    set_parser.add_argument(
        "name",
        help="The name of the parameter to set. The parameter names can be found in the manual https://biosmoother.readthedocs.io",
    )
    set_parser.add_argument(
        "val",
        help="The value to set.",
    )
    set_parser.set_defaults(func=set_smoother)

    get_parser = main_parser.add_parser(
        "get",
        help="Retrieve the value of a parameter from the current session of an index. This can also be done on in the graphical user interface.",
    )
    get_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    get_parser.add_argument(
        "name",
        help="The name of the parameter to get. The parameter names can be found in the manual https://biosmoother.readthedocs.io. If a '.' is given for the name argument, the whole session will be printed. If 'settings' is given for the name argument, all settings will be printed.",
    )
    get_parser.set_defaults(func=get_smoother)

    if not HIDE_SUBCOMMANDS_MANUAL:
        info_parser = main_parser.add_parser(
            "par_info", argument_default=argparse.SUPPRESS
        )
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
            help="Path to the index directory generated with the init command.",
        )
        test_parser.add_argument(
            "-S",
            "--seed",
            help="Seed for random configurations.",
            default=None,
            type=int,
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
            help="Path to the index directory generated with the init command.",
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

    ploidy_parser = main_parser.add_parser(
        "ploidy", help="Add a ploidy file to the index."
    )
    ploidy_parser.add_argument(
        "index_prefix",
        help="Path to the index directory generated with the init command.",
    )
    ploidy_parser.add_argument(
        "ploidy_file",
        help="File that specifies the ploidy for each contig. The format of this file is described in the manual: https://biosmoother.rtdf.io",
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


def make_versioned_parser():
    parser = make_main_parser()
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=version("libbiosmoother"),
        help="Print libSmoother's version and exit.",
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
    return parser


def main():
    parser = make_versioned_parser()

    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
