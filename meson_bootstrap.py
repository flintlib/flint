#!/usr/bin/env python

from os.path import join, dirname, abspath
from sys import executable as python, argv
from subprocess import check_call

this_dir = dirname(abspath(__file__))
meson_build_dir = join(this_dir, '_meson_build')
generate_script = join(meson_build_dir, 'generate_meson_build.py')
check_call([python, generate_script, this_dir] + argv[1:])
