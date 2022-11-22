from libContactMapping import make_main_parser

if __name__ == "__main__":
    parser = make_main_parser()

    args = parser.parse_args()

    args.func(args)
