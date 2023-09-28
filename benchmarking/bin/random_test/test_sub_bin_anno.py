import fileinput
import sys

d = int(sys.argv[2])

with fileinput.input(sys.argv[1]) as anno_file:
    for line in anno_file:
        if len(line) == 0  or line[0] == '#':
            print(line, end="")
        else:
            chr1, a, b, pos1, pos2, *extra = line[:-1].split()
            pos1 = int(pos1)
            pos2 = int(pos2)
            chr1 += "_" + str(pos1 // d)
            pos1 = pos1 % d
            pos2 = pos2 % d
            print(chr1, a, b, pos1, pos2, *extra, sep="\t")