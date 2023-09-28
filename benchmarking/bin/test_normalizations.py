import libbiosmoother

import json
import sys
import shelve
import sys



def lib_sps_print(s):
    pass

def conf_quarry_basic(quarry):
    #warnings.filterwarnings('ignore')
    with libbiosmoother.open_default_json() as default_file:
        default_json = json.load(default_file)
        quarry.set_value(["settings"], default_json)
    quarry.set_value(["settings", "filters", "cut_off_bin"], "fit_chrom_smaller")
    quarry.set_value(["settings", "filters", "show_contig_smaller_than_bin"], True)
    quarry.set_value(["settings", "interface", "fixed_bin_size"], True)
    quarry.set_value(["settings", "interface", "add_draw_area", "val"], 0)
    quarry.set_value(["settings", "normalization", "scale"], "dont")
    quarry.set_value(["settings", "filters", "symmetry"], "all")


def conf_quarry_data(quarry):
    quarry.set_value(["settings", "normalization", "log_base", "val"], 0)

def conf_quarry_heatmap(quarry, max_val):
    quarry.set_value(["settings", "normalization", "color_range", "val_max"], max_val)
    quarry.set_value(["settings", "normalization", "log_base", "val"], 10)


def set_bin_size(quarry, bin_size):
    div = quarry.get_value(["dividend"])
    if bin_size % div != 0:
        print("WARNING: uneven division by index dividend", file=sys.stderr)
    if bin_size < div:
        print("WARNING: dividend larger than value", file=sys.stderr)
    bin_size = max(1, bin_size // div)
    quarry.set_value(["settings", "interface", "fixed_bin_size_x", "val"], bin_size)
    quarry.set_value(["settings", "interface", "fixed_bin_size_y", "val"], bin_size)

def tsv_to_ret(data, tsv):
    ret = [(x[:-1], x[-1]) for x in tsv]
    ret.sort()
    if False:
        data["bin_coords"] = [a for a, _ in ret]
    data["bin_vals"] = [b for _, b in ret]
    return max([b for _, b in ret])

def quarry_whole_window(data, quarry):
    canvas_size_x, canvas_size_y = quarry.get_canvas_size(lib_sps_print)
    quarry.set_value(["area"], {"x_start": 0, "x_end": canvas_size_x, "y_start": 0, "y_end": canvas_size_y})
    
    conf_quarry_data(quarry)
    max_v = tsv_to_ret(data, quarry.get_heatmap_export(lib_sps_print))
    
    conf_quarry_heatmap(quarry, max_v)
    data["heatmap"] = quarry.get_heatmap(lib_sps_print)

def print_window_amount(quarry, window_size):
    canvas_size_x, canvas_size_y = quarry.get_canvas_size(lambda s: None)
    div = quarry.get_value(["dividend"])
    print("num windows on genome:", canvas_size_x / (window_size // div), "x", canvas_size_y / (window_size // div))

def quarry_chunked_window(data, quarry, window_size):
    canvas_size_x, canvas_size_y = quarry.get_canvas_size(lambda s: None)
    div = quarry.get_value(["dividend"])
    tsv = []
    runtimes = {}
    heatmap = None
    for x_start in range(0, canvas_size_x, window_size // div):
        for y_start in range(0, canvas_size_y, window_size // div):
            quarry.set_value(["area"], {"x_start": x_start, "x_end": min(canvas_size_x, x_start + window_size // div), 
                                        "y_start": y_start, "y_end":  min(canvas_size_y, y_start + window_size // div)})
            conf_quarry_data(quarry)
            quarry.clear_cache()
            tsv.extend(quarry.get_heatmap_export(lib_sps_print))
            for k, v in quarry.get_runtimes():
                if not k in runtimes:
                    runtimes[k] = []
                runtimes[k].append(v)
            
    max_v = tsv_to_ret(data, tsv)
    for x_start in range(0, min(canvas_size_x, 3* window_size // div), window_size // div):
        for y_start in range(0, min(canvas_size_y, 3* window_size // div), window_size // div):
            quarry.set_value(["area"], {"x_start": x_start, "x_end": min(canvas_size_x, x_start + window_size // div), 
                                        "y_start": y_start, "y_end":  min(canvas_size_y, y_start + window_size // div)})
            conf_quarry_heatmap(quarry, max_v)
            heatmap_local = quarry.get_heatmap(lib_sps_print)
            if heatmap is None:
                heatmap = heatmap_local
            else:
                for key, val in heatmap_local.items():
                    heatmap[key].extend(val)
    data["heatmap"] = heatmap
    data["runtimes"] = runtimes


class ShelveTupleDict:
    def __init__(self, filename, flag='c'):
        self.data = shelve.open(filename, flag=flag)
        
    def __to_to_key(self, tup):
        return ".".join(str(x) for x in tup)
    
    def __getitem__(self, tup):
        return self.data[self.__to_to_key(tup)]
    
    def __setitem__(self, tup, val):
        self.data[self.__to_to_key(tup)] = val

    def close(self):
        self.data.close()

# usage: python3 test_normalizations.py index_file output_file bin_size window_size extra_params...

index_file = sys.argv[1]
output_file = sys.argv[2]
bin_size = int(sys.argv[3])
window_size = int(sys.argv[4])
if window_size == 0:
    whole = True
else:
    whole = False
print("index_file", index_file, "output_file", output_file, "bin_size", bin_size, 
      "window_size", window_size, "whole", whole)

with shelve.open(output_file, flag="c") as data:

    print("loading index")
    index = libbiosmoother.Quarry(index_file)
    conf_quarry_basic(index)

    set_bin_size(index, bin_size)

    extra_params = sys.argv[5:]
    for extra_param in extra_params:
        key, val = extra_param.split("=")
        if val == "True":
            val = True
        elif val == "False":
            val = False
        elif val == "None":
            val = None
        elif val.isdigit():
            val = int(val)
        elif val.replace(".", "", 1).isdigit():
            val = float(val)
        else:
            val = str(val)
        print(key.split("."), val)
        index.set_value(key.split("."), val)

    print("querying")
    if whole:
        quarry_whole_window(data, index)
    else:
        quarry_chunked_window(data, index, window_size)

print("done")
