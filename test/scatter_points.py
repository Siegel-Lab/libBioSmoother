import random
from bokeh.plotting import figure, save
import numpy as np

f = figure(x_range=[0,1], y_range=[0,1])
f.sizing_mode = "scale_both"

n = 100000

xs = [random.random() for _ in range(n)]
ys = [ x + n for x, n in zip(xs, np.random.normal(0, 0.25, n))]
f.circle(x=xs, y=ys, fill_alpha=0.2, line_color=None, fill_color="#3b528b")

save(f)