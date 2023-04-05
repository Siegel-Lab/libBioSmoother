import os

if not os.path.exists("test_data"):
    os.mkdir("test_data")

with open("test_data/contig_sizes.tsv", "w") as out_file:
    out_file.write("Chr1\t4\n")

with open("test_data/sample_data.tsv", "w") as out_file:
    read_cnt = 1
    for pos_x, pos_y, cnt in [
        (0, 0, 7),
        (1, 0, 21),
        (2, 0, 35),
        (3, 0, 14),
        
        (0, 1, 4),
        (1, 1, 24),
        (2, 1, 20),
        (3, 1, 3),

        (0, 2, 12),
        (1, 2, 54),
        (2, 2, 60),
        (3, 2, 12),

        (0, 3, 8),
        (1, 3, 48),
        (2, 3, 40),
        (3, 3, 16),
    ]:
        for _ in range(cnt):
            out_file.write("\t".join(["r" + str(read_cnt), "Chr1", str(pos_x), "Chr1", str(pos_y)]) + "\n")
            read_cnt += 1