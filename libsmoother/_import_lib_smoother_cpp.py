import os
import sys

sys.path.append(os.environ["PREFIX"] + "/lib")
sys.path.append(os.environ["CONDA_PREFIX"] + "/lib")

from libsmoothercpp import *