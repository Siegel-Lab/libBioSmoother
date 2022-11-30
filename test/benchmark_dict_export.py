import sys
import os
sys.path.append(os.getcwd())
from build_rel.libContactMapping import test_py_dict, test_cpp_dict
import random
import time

py_wins = 0
cpp_wins = 0

print("a", "b", "cpp", "py", sep="\t")
for _ in range(100):
    a = random.randrange(10000)
    b = random.randrange(10000)

    t1 = time.perf_counter()
    d1 = test_py_dict(a, b)
    t2 = time.perf_counter()
    d2 = test_cpp_dict(a, b)
    t3 = time.perf_counter()

    if False:
        for key in set([*d1.keys(), *d2.keys()]):
            assert len(d1[key]) == len(d2[key])
            for x, y in zip(d1[key], d2[key]):
                assert x == y
    print(a, b, t3-t2, t2-t1, sep="\t")
    if t3-t2 < t2-t1:
        cpp_wins += 1
    else:
        py_wins += 1
print()
print("cpp won", cpp_wins, "times, while py won", py_wins, "times")