import fileinput
import sys

d = int(sys.argv[2])

with fileinput.input(sys.argv[1]) as genome_file:
    for line in genome_file:
        name, l = line.split()
        l = int(l)
        for i in range(1, l // d):
            print(name + "_" + str(i), d, sep="\t")
        print(name + "_" + str(l // d), l % d, sep="\t")
