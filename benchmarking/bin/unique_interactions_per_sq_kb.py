import fileinput
import sys

last_pos = None
pos_cnt = 0
with fileinput.input(sys.argv[1]) as pair_sam_in_file:
    for line in pair_sam_in_file:
        if len(line) == 0  or line[0] == '#':
            continue
        cols = line.split()
        pos = cols[1:5]
        if pos != last_pos:
            pos_cnt += 1
            last_pos = pos


total_genome_size = 0
with fileinput.input(sys.argv[2]) as genome_sizes_file:
    for line in genome_sizes_file:
        if len(line) == 0  or line[0] == '#':
            continue
        _, size = line.split()
        total_genome_size += int(size)


KB = 1000

ret = ((KB**2) * pos_cnt) / (total_genome_size ** 2)
print('{:.20f}'.format(ret))
