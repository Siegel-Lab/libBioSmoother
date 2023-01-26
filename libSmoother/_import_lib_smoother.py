import os

if os.environ["smoother_import_mode"] == "debug":
    from .dbg.libSmoother import *
elif os.environ["smoother_import_mode"] == "relwdbg":
    from .rel_w_dbg.libSmoother import *
else:
    from .rel.libSmoother import *
