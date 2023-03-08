#!/usr/bin/env python

import os
import sys
from datetime import datetime

os.makedirs(os.path.dirname(sys.argv[2]), exist_ok=True)


build_time = datetime.now().strftime("%Y-%m-%d-%H:%M:%S")

with open(sys.argv[1], "r") as in_file:
    with open(sys.argv[2], "w") as out_file:
        for line in in_file.readlines():
            line = line.replace("@BUILD_TIME@", build_time)
            out_file.write(line)

# set timestamp to that of original file to hide change from cmake
in_info = os.stat(sys.argv[1])
os.utime(sys.argv[2], (in_info.st_atime, in_info.st_mtime))