import fileinput
import sys

d = int(sys.argv[2])

with fileinput.input(sys.argv[1]) as pair_sam_in_file:
    for line in pair_sam_in_file:
        if len(line) == 0  or line[0] == '#':
            print(line, end="")
        else:
            readID, chrom1, pos1, chrom2, pos2, strand1, strand2, pair_type, mapq1, mapq2 = line.split()
            pos1 = int(pos1)
            pos2 = int(pos2)
            chrom1 += "_" + str(pos1 // d)
            pos1 = pos1 % d
            chrom2 += "_" + str(pos2 // d)
            pos2 = pos2 % d
            print(readID, chrom1, pos1, chrom2, pos2, strand1, strand2, pair_type, mapq1, mapq2, sep="\t")