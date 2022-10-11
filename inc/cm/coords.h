#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{
    std::vector<> axisCoordsHelper()
    {

    }

    void Computation::setAxisCoords( bool bX )
    {
    def bin_cols_or_rows(self, h_bin, start=0, end=None, none_for_chr_border=False, chr_filter=[], 
                         produce_smaller_bins=True, is_canceld=lambda: False):
                if end is None:
            end = self.chr_start_pos["end"]
        h_bin = max(1, h_bin)
        ret = []
        bin_start = max(0, int(start))
        if none_for_chr_border:
            ret.append(None)
        for chr_name in chr_filter:
            x_start = self.chr_start_pos[chr_name]
            x_end = x_start + self.chr_sizes[chr_name]
            x = max(int(start), x_start)
            while x <= min(end, x_end):
                if produce_smaller_bins or x + h_bin <= x_end:
                    ret.append([chr_name, x - x_start, bin_start, min(h_bin, x_end - x)])
                    bin_start += min(h_bin, x_end - x)
                if is_canceld():
                    return
                x += h_bin
            if none_for_chr_border:
                ret.append(None)
        return ret

    }
}