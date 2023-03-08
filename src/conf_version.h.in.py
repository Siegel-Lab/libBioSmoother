#!/usr/bin/env python

import subprocess
import sys
import os

def run_command(cmd):
    return subprocess.check_output(cmd)

git_commit_hash="-".join(run_command(["git", "log", "-1", "--format=%h-%ci"]).decode().split()[:2])
git_status="" if run_command(["git", "status", "-s"]) == "" else "D-"
cmake_version = sys.argv[2]

os.makedirs(os.path.dirname(sys.argv[3]), exist_ok=True)

# configure the new version file
out_lines = []
with open(sys.argv[1], "r") as in_file:
    for line in in_file.readlines():
        line = line.replace("@CMAKE_VERSION@", cmake_version)
        line = line.replace("@GIT_COMMIT_HASH@", git_commit_hash)
        line = line.replace("@GIT_STATUS@", git_status)
        out_lines.append(line)

file_changed = False
if not os.path.isfile(sys.argv[3]):
    file_changed=True
else:
    with open(sys.argv[3], "r") as in_file:
        lines = in_file.readlines()
        file_changed = lines == out_lines

if file_changed:
    print("writing new version file")
    with open(sys.argv[3], "w") as out_file:
        for line in out_lines:
            out_file.write(line)
