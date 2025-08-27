from pandas import DataFrame
import cooler
import os
import sys


def icing(bin_values, axis_size):
    has_non_zero = False
    for v in bin_values:
        if v != 0:
            has_non_zero = True
            break
    if not has_non_zero:
        return bin_values
    bins = DataFrame(
        data={
            "chrom": ["chr1" for _ in range(axis_size)],
            "start": [i * 10 for i in range(axis_size)],
            "end": [(i + 1) * 10 for i in range(axis_size)],
        }
    )
    pixels = DataFrame(
        data={
            "bin1_id": [i // axis_size for i in range(len(bin_values))],
            "bin2_id": [i % axis_size for i in range(len(bin_values))],
            "count": bin_values,
        }
    )
    print("cooler_interface, axis_size:", axis_size)
    cooler.create_cooler(
        ".tmp.cooler", bins, pixels, symmetric_upper=False, triucheck=False
    )
    clr = cooler.Cooler(".tmp.cooler")
    bias, stats = cooler.balance_cooler(clr)
    ret = []
    for i, v in enumerate(bin_values):
        ret.append(v * bias[i // axis_size] * bias[i % axis_size])
    os.remove(".tmp.cooler")
    return ret


class CoolerIterator:
    def __init__(self, cooler_path, bin_size=None):
        self.clr = cooler.Cooler(cooler_path)
        self.non_existant_chr_warning_delivered = False
        if not bin_size is None and self.clr.binsize != bin_size:
            raise ValueError(
                "bin size of cooler file do not match the base resolution of the smoother index."
            )

    def iterate(self, chr_x, chr_y):
        # print(chr_x, chr_y)
        if chr_x not in self.clr.chromnames or chr_y not in self.clr.chromnames:
            if not self.non_existant_chr_warning_delivered:
                if chr_x not in self.clr.chromnames:
                    print(
                        "Warning: Chromosome of index does not exist in cooler file: ",
                        chr_x,
                        "existing chr names are: ",
                        self.clr.chromnames,
                        file=sys.stderr,
                    )
                if chr_y not in self.clr.chromnames:
                    print(
                        "Warning: Chromosome of index does not exist in cooler file: ",
                        chr_y,
                        "existing chr names are: ",
                        self.clr.chromnames,
                        file=sys.stderr,
                    )
                self.non_existant_chr_warning_delivered = True
            return
        for idx, row in (
            self.clr.matrix(balance=False, as_pixels=True, join=True, sparse=True)
            .fetch(chr_x, chr_y)
            .iterrows()
        ):
            yield row["chrom1"], row["start1"], row["chrom2"], row["start2"], row[
                "count"
            ]
