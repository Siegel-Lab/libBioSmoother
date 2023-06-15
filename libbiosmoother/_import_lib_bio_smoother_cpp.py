import os
import sys

if "LIB_BIO_SMOOTHER_CPP_IMPORT" in os.environ:
    sys.path.append(os.environ["LIB_BIO_SMOOTHER_CPP_IMPORT"])
if "PREFIX" in os.environ:
    sys.path.append(os.environ["PREFIX"] + "/lib")
if "CONDA_PREFIX" in os.environ:
    sys.path.append(os.environ["CONDA_PREFIX"] + "/lib")

from libbiosmoothercpp import *
