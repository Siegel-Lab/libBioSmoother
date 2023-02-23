from .indexer import *
from .quarry import Quarry
from .export import export_tsv, export_png, export_svg
import argparse
import os
import shutil


def init(args):
    Indexer(args.index_prefix, strict=True).create_session(
        args.chr_len, args.dividend, args.anno_path, args.order_path, args.test
    )


def reset(args):
    for possible in [args.index_prefix, args.index_prefix + ".smoother_index"]:
        if os.path.exists(possible) and os.path.isdir(possible):
            os.remove(possible + "/session.json")
            shutil.copy(possible + "/default_session.json", possible + "/session.json")
            return
    raise RuntimeError("the given index", args.index_prefix, "does not exist.")


def repl(args):
    Indexer(args.index_prefix).add_replicate(
        args.path,
        args.name,
        args.group,
        args.no_groups,
        args.keep_points,
        args.only_points,
        args.no_map_q,
        args.no_multi_map,
        args.no_cat,
        args.no_strand,
        args.shekelyan,
        args.force_upper_triangle,
    )


def norm(args):
    Indexer(args.index_prefix).add_normalization(
        args.path,
        args.name,
        args.group,
        args.no_groups,
        args.keep_points,
        args.only_points,
        args.no_map_q,
        args.no_multi_map,
    )


def export_smoother(args):
    session = Quarry(args.index_prefix)
    if args.export_prefix is not None:
        session.set_value(["settings", "export", "prefix"], args.export_prefix)
    if args.export_selection is not None:
        session.set_value(["settings", "export", "selection"], args.export_selection)
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
        "anno_path", help="Path to a file that contains the annotations.", default=None
    )
    init_parser.add_argument(
        "--order_path",
        help="Path to a file that contains the order of annotations, default: gene first, then alphabetic for others.",
    )
    init_parser.add_argument(
        "-d",
        "--dividend",
        type=int,
        default=10000,
        help="Divide all coordinates by this number. Larger numbers will reduce the index size and preprocessing time. However, bins with a size below this given number cannot be displayed.",
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
        "--no_strand",
        action="store_true",
        help="Do not store strand information. (default: off)",
    )
    repl_parser.add_argument(
        "-u",
        "--force_upper_triangle",
        action="store_true",
        help="Mirror all interactions to the upper triangle. (default: off)",
    )
    repl_parser.set_defaults(func=repl)
    repl_parser.add_argument(
        "--keep_points", help=argparse.SUPPRESS, action="store_true"
    )
    repl_parser.add_argument(
        "--only_points", help=argparse.SUPPRESS, action="store_true"
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


def make_main_parser():
    parser = argparse.ArgumentParser(description="")
    sub_parsers = parser.add_subparsers(
        help="Sub-command that shall be executed.", dest="cmd"
    )
    sub_parsers.required = True
    add_parsers(sub_parsers)
    return parser
