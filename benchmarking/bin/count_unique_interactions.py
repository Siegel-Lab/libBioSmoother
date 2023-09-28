import fileinput
import sys

def cnt_unique_interactions(in_filename):
    last_pos = None
    num_interactions = 0
    num_unique_interactions = 0
    with fileinput.input(in_filename) as pair_sam_in_file:
        for line in pair_sam_in_file:
            if len(line) == 0  or line[0] == '#':
                continue
            num_interactions += 1
            cols = line.split()
            pos = cols[1:5]
            if pos != last_pos:
                num_unique_interactions += 1
                last_pos = pos

    return num_interactions, num_unique_interactions

if __name__ == "__main__":
    num_interactions, num_unique_interactions = cnt_unique_interactions(sys.argv[1])
    print(num_interactions, "number_of_interactions", sep="\t")
    print(num_unique_interactions, "number_of_unique_interactions", sep="\t")
