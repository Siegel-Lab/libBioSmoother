import sys
import fileinput

def count_maps(chr, xa):
    pos = [chr] + [x.split(",")[0] for x in xa.split(';') if x != '']
    #print(pos, len(pos), len(set(pos)))
    return len(pos), len(set(pos))

count_map = {}
with fileinput.input(sys.argv[1]) as in_file:
    for line in in_file:
        if line[0] == "#":
            continue
        read_name, chr1, pos1, chr2, pos2, strand1, strand2, mm, mapq1, mapq2, xa1, xa2 = line[:-1].split('\t')
        num_maps1, num_contigs1 = count_maps(chr1, xa1)
        num_maps2, num_contigs2 = count_maps(chr2, xa2)
        num_maps = max(num_maps1, num_maps2)
        num_contigs = max(num_contigs1, num_contigs2)
        if (num_maps, num_contigs) not in count_map:
            count_map[(num_maps, num_contigs)] = 0
        count_map[(num_maps, num_contigs)] += 1

print("read count", "mapping pos", "contig count", sep="\t")
for (num_maps, num_contigs), count in sorted(count_map.items()):
    print(count, num_maps, num_contigs, sep="\t")