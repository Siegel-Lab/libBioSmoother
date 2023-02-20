from .indexer_parser import make_main_parser
from ._import_lib_smoother_cpp import SPS_VERSION, LIB_SMOOTHER_CPP_VERSION
from importlib.metadata import version
import argparse


def main():
    parser = make_main_parser()
    parser.add_argument(
        "-v", "--version", action="version", version=version("libsmoother")
    )
    parser.add_argument(
        "--version_smoother_cpp",
        help=argparse.SUPPRESS,
        action="version",
        version=LIB_SMOOTHER_CPP_VERSION,
    )
    parser.add_argument(
        "--version_sps", help=argparse.SUPPRESS, action="version", version=SPS_VERSION
    )

    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
