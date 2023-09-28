import os
import sys
import pickle
import shelve

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

def read_line_from_file(filename, line_num):
    with open(filename, 'r') as f:
        return list(f.readlines())[line_num].strip()

def collect(data, file, key_prefix):
    if os.path.isfile(file + ".status"):
        if read_line_from_file(file + ".status", 0) == "OK":
            if len(sys.argv) < 2 or sys.argv[1] == "True":
                with shelve.open(file, flag="r") as in_file:
                    for key, value in in_file.items():
                        data[key_prefix + [key]] = value
            print(key_prefix, "OK")
            return True
    print(key_prefix, "is not finished")
    return False


if len(sys.argv) < 2 or sys.argv[1] == "True":
    data = ShelveTupleDict("norm_corelation.shelf", 'n')
else:
    data = None

not_done = 0

for bin_size in [50000, 100000, 500000]:
    if not collect(data, "norm_out/grid_seq_ground_truth." + str(bin_size), ["grid_seq", "GT", bin_size]):
        not_done += 1
    if not collect(data, "norm_out/radicl_seq_ground_truth." + str(bin_size), ["radicl_seq", "GT", bin_size]):
        not_done += 1
    if not collect(data, "norm_out/raw_data_radicl." + str(bin_size), ["raw_radicl", "GT", bin_size]):
        not_done += 1
    if not collect(data, "norm_out/ddd_ground_truth." + str(bin_size), ["ddd", "GT", bin_size]):
        not_done += 1
    if not collect(data, "norm_out/cooler_ground_truth." + str(bin_size), ["cooler", "GT", bin_size]):
        not_done += 1
    if not collect(data, "norm_out/ice_ground_truth." + str(bin_size), ["ice", "GT", bin_size]):
        not_done += 1
    if not collect(data, "norm_out/raw_data_hic." + str(bin_size), ["raw_hic", "GT", bin_size]):
        not_done += 1

    window_size = 5000000
    for num_samples in [1, 2, 10, 25, 100, 250, 1000, 2500, 10000]:
        if not collect(data, "norm_out/grid_seq." + str(bin_size) + "." + str(window_size) + "." + str(num_samples),
                ["grid_seq", window_size, num_samples, bin_size]):
            not_done += 1
        if not collect(data, 
                        "norm_out/radicl_seq." + str(bin_size) + "." + str(window_size) + "." + str(num_samples),
                ["radicl_seq", window_size, num_samples, bin_size]):
            not_done += 1
        if not collect(data, "norm_out/ice." + str(bin_size) + "." + str(window_size) + "." + str(num_samples), 
                ["ice", window_size, num_samples, bin_size]):
            not_done += 1
        if not collect(data, "norm_out/ddd." + str(bin_size) + "." + str(window_size) + "." + str(num_samples),
                    ["ddd", window_size, num_samples, bin_size]):
            not_done += 1

bin_size = 100000
for window_size in [1000000, 5000000, 10000000]:
    for num_samples in [1, 2, 10, 25, 100, 250, 1000, 2500, 10000]:
        if not collect(data, "norm_out/grid_seq." + str(bin_size) + "." + str(window_size) + "." + str(num_samples),
                ["grid_seq", window_size, num_samples, bin_size]):
            not_done += 1
        if not collect(data, 
                        "norm_out/radicl_seq." + str(bin_size) + "." + str(window_size) + "." + str(num_samples),
                ["radicl_seq", window_size, num_samples, bin_size]):
            not_done += 1
        if not collect(data, "norm_out/ice." + str(bin_size) + "." + str(window_size) + "." + str(num_samples), 
                ["ice", window_size, num_samples, bin_size]):
            not_done += 1
        if not collect(data, "norm_out/ddd." + str(bin_size) + "." + str(window_size) + "." + str(num_samples),
                    ["ddd", window_size, num_samples, bin_size]):
            not_done += 1

print("not done:", not_done)

if len(sys.argv) < 2 or sys.argv[1] == "True":
    data.close()