from bokeh.plotting import figure, output_file, save
import sys
import numpy

def read_tsv(file_name, n_keys):
    with open(file_name, "r") as in_file:
        for line in in_file.readlines():
            if line[0] == "#":
                continue
            columns = line.split("\t")
            yield "\t".join(columns[:n_keys]), float(columns[n_keys])

def plot_points(points):
    output_file(filename="compare_tsv.html")
    plot = figure(sizing_mode="stretch_both")
    plot.x(x=[p[0] for p in points], y=[p[0] for p in points])
    save(plot)


def correlate_tsv(file_a, file_b, n_keys=1):
    n_keys = int(n_keys)
    file_a_dict = dict(list(read_tsv(file_a, n_keys)))
    len_file_a = len(file_a_dict)
    points = []
    len_file_b = 0
    for key, val in read_tsv(file_b, n_keys):
        len_file_b += 1
        if key in file_a_dict:
            points.append((file_a_dict[key], val))
            del file_a_dict[key]

    print("points shared:", len(points))
    print("points only in a:", len_file_a - len(points))
    print("points only in b:", len_file_b - len(points))

    plot_points(points)

    print("correlation: ", numpy.corrcoef(points, rowvar=False))


if __name__ == "__main__":
    correlate_tsv(*sys.argv[1:])
