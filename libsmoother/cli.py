from .indexer_parser import make_main_parser

def main():
    parser = make_main_parser()

    args = parser.parse_args()

    args.func(args)

if __name__ == "__main__":
    main()
