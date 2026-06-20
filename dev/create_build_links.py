#!/usr/bin/env python3
import os
import shutil
import stat
import sys


EXEC_DIRS = {"profile", "test", "tune"}
SKIP_SUFFIXES = {".o", ".obj", ".d", ".pdb", ".ilk", ".rsp"}


def is_executable_file(path):
    if not os.path.isfile(path):
        return False
    name = os.path.basename(path)
    if name.endswith(".exe"):
        return True
    if os.path.splitext(name)[1] in SKIP_SUFFIXES:
        return False
    mode = os.stat(path).st_mode
    return bool(mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH))


def link_or_copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)

    if os.path.lexists(dst):
        if os.path.isdir(dst) and not os.path.islink(dst):
            return
        os.unlink(dst)

    rel_src = os.path.relpath(src, os.path.dirname(dst))
    try:
        os.symlink(rel_src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main(argv):
    if len(argv) != 2:
        print("usage: create_build_links.py BUILDDIR", file=sys.stderr)
        return 2

    build_dir = os.path.abspath(argv[1])
    src_build_dir = os.path.join(build_dir, "src")

    if not os.path.isdir(src_build_dir):
        return 0

    for root, dirs, files in os.walk(src_build_dir):
        dirs[:] = [d for d in dirs if not d.endswith(".p")]
        rel_root = os.path.relpath(root, src_build_dir)
        parts = rel_root.split(os.sep)

        if not any(part in EXEC_DIRS for part in parts):
            continue

        for name in files:
            src = os.path.join(root, name)
            if not is_executable_file(src):
                continue

            dst = os.path.join(build_dir, rel_root, name)
            link_or_copy(src, dst)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
