import os

if os.environ["smoother_import_mode"] == "debug":
    from .dbg.libContactMapping import *
elif os.environ["smoother_import_mode"] == "relwdbg":
    from .rel_w_dbg.libContactMapping import *
else:
    from .rel.libContactMapping import *
