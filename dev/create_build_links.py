#!/usr/bin/env python3
import os
import sys


def main(argv):
    if len(argv) != 2:
        print("usage: create_build_links.py BUILDDIR", file=sys.stderr)
        return 2

    build_dir = os.path.abspath(argv[1])
    src_build_dir = os.path.join(build_dir, "src")

    if not os.path.isdir(src_build_dir):
        return 0

    for name in os.listdir(src_build_dir):
        src = os.path.join(src_build_dir, name)
        dst = os.path.join(build_dir, name)

        if not os.path.isdir(src) or os.path.exists(dst) or os.path.islink(dst):
            continue

        rel_src = os.path.relpath(src, os.path.dirname(dst))
        try:
            os.symlink(rel_src, dst, target_is_directory=True)
        except OSError:
            # Symlinks may be unavailable on some platforms. The canonical Meson
            # paths under build/src still work, so this compatibility layer is
            # best-effort.
            pass

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
