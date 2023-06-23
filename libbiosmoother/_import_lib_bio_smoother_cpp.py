import os
import sys

if "PREFIX" in os.environ:
    sys.path.append(os.environ["PREFIX"] + "/lib")
if "CONDA_PREFIX" in os.environ:
    sys.path.append(os.environ["CONDA_PREFIX"] + "/lib")

from libbiosmoothercpp import *
