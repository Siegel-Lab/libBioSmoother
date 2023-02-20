from bokeh.plotting import figure, output_file, save
from bokeh.palettes import Viridis256
from bokeh.models import ColumnDataSource, LinearColorMapper
import sys
from scipy.stats.stats import pearsonr
from random import shuffle


def read_tsv(file_name, n_keys):
    with open(file_name, "r") as in_file:
        for line in in_file.readlines():
            if line[0] == "#":
                continue
            columns = line.split("\t")
            yield "\t".join(columns[:n_keys]), float(columns[n_keys])


def plot_points(x, y, hover, file_a, file_b, out_file_name):
    output_file(filename=out_file_name)
    plot = figure(sizing_mode="stretch_both", tooltips=[("desc", "@hover")])
    mapper = LinearColorMapper(palette=Viridis256)
    # plot.circle(x="x", y="y", size=10, source=ColumnDataSource(data={
    #                    "x": x, "y": y, "hover": hover, "c": [int(h.split()[1]) for h in hover]
    #            }),
    #            color={'field': 'c', 'transform': mapper})
    plot.circle(
        x="x",
        y="y",
        size=10,
        source=ColumnDataSource(data={"x": x, "y": y}),
        color="blue",
        fill_alpha=0.5,
    )
    plot.xaxis.axis_label = file_a
    plot.yaxis.axis_label = file_b
    save(plot)


def correlate_tsv(file_a, file_b, n_keys=1, out_file_name="compare_tsv.html"):
    n_keys = int(n_keys)
    file_a_dict = dict(list(read_tsv(file_a, n_keys)))
    len_file_a = len(file_a_dict)
    x = []
    y = []
    keys = []
    len_file_b = 0
    list_b = list(read_tsv(file_b, n_keys))
    # shuffle(list_b)
    for key, val in list_b:
        len_file_b += 1
        if key in file_a_dict:
            x.append(file_a_dict[key])
            y.append(val)
            keys.append(key)
            del file_a_dict[key]

    print("points shared:", len(x))
    print("points only in a:", len_file_a - len(x), len_file_a, len(x))
    print("points only in b:", len_file_b - len(x), len_file_b, len(x))

    plot_points(x, y, keys, file_a, file_b, out_file_name)

    print("correlation: ", pearsonr(x, y))


if __name__ == "__main__":
    correlate_tsv(*sys.argv[1:])
