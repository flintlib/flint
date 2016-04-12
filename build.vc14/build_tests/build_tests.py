'''
Set up Visual Sudio to build a specified MPIR configuration

Copyright (C) 2011, 2012, 2013, 2014, 2015 Brian Gladman
'''

from operator import itemgetter
from os import listdir, walk, unlink, makedirs
from os.path import split, splitext, isdir, relpath, join, exists
from os.path import dirname, normpath
from copy import deepcopy
from sys import argv, exit
from filecmp import cmp
from shutil import copy
from re import compile, search
from collections import defaultdict
from uuid import uuid4
from time import sleep

from _msvccompiler import MSVCCompiler

import argparse
parser = argparse.ArgumentParser(description='Build flint tests')

# for script debugging
parser.add_argument('--debug', choices=["True", "False"], default="False")

# what to build
parser.add_argument('--platform', default="x64")
parser.add_argument('--configuration', default="Release")
parser.add_argument('--library-type', choices=["dll", "lib"], default="lib")
parser.add_argument('--interfaces-tests', choices=["True", "False"], default="True")

args = parser.parse_args()

print(args)

debug = args.debug == "True"
intd = '\\%s\\%s\\' % (args.platform, args.configuration)
library_type = args.library_type
build_interfaces_tests = args.interfaces_tests == "True"

# The path to flint, solution and project directories
script_dir = dirname(__file__)
project_name = 'flint'
build_vc = 'build.vc14'
flint_dir = normpath(join(script_dir, '../../'))
solution_dir = normpath(join(flint_dir, build_vc))

try:
  input = raw_input
except NameError:
  pass

app_type, lib_type, dll_type = 0, 1, 2
app_str = ('Application', 'StaticLibrary', 'DynamicLibrary')
app_ext = ('.exe', '.lib', '.dll')

# copy from file ipath to file opath but avoid copying if
# opath exists and is the same as ipath (this is to avoid
# triggering an unecessary rebuild).

def write_f(ipath, opath):
  if exists(ipath) and not isdir(ipath):
    if exists(opath):
      if isdir(opath) or cmp(ipath, opath):
        return
    copy(ipath, opath)

ignore_dirs = ( '.git', 'doc', 'examples', 'lib', 'exe', 'dll', 'win_hdrs')
req_extns = ( '.h', '.c', '.cc', '.cpp' )

def find_src(path):

  c, h, cx, hx, t, tx, p = [], [], [], [], [], [], []
  for root, dirs, files in walk(path):
    if 'template' in root:
      continue
    _, _t = split(root)
    if _t in ignore_dirs:
      continue
    if 'build.vc' in root:
      for di in list(dirs):
        dirs.remove(di)
    for di in list(dirs):
      if di in ignore_dirs:
        dirs.remove(di)
      if 'template' in di:
        dirs.remove(di)
    relp = relpath(root, flint_dir)
    if relp == '.':
      relp = ''
    for f in files:
      if 'template' in f:
        continue
      n, x = splitext(f)
      if x not in req_extns:
        continue
      pth, leaf = split(root)
      fp = join(relp, f)
      if leaf == 'tune':
        continue
      if leaf == 'test':
        p2, l2 = split(pth)
        l2 = '' if l2 == 'flint2' else l2
        if 'flintxx' in pth:
          tx += [(l2, fp)]
        else:
          t += [(l2, fp)]
      elif leaf == 'profile':
        p2, l2 = split(pth)
        l2 = '' if l2 == 'flint2' else l2
        p += [(l2, fp)]
      elif leaf == 'flintxx':
        cx += [fp]
      elif x == '.c':
        c += [(leaf, fp)]
      elif x == '.h':
        if n.endswith('xx'):
          hx += [fp]
        else:
          h += [fp]
  for x in (c, h, cx, hx, t, tx, p):
    x.sort()
  return (c, h, cx, hx, t, tx, p)

c, h, cx, hx, t, tx, p = find_src(flint_dir)

# def compile(self, sources,
#       output_dir=None, macros=None, include_dirs=None, debug=0,
#       extra_preargs=None, extra_postargs=None, depends=None):

# def link(self, target_desc, objects, output_filename, output_dir=None,
#       libraries=None, library_dirs=None, runtime_library_dirs=None,
#       export_symbols=None, debug=0, extra_preargs=None,
#       extra_postargs=None, build_temp=None, target_lang=None):


cc = MSVCCompiler()
error_list = []
inc_dirs = [
  '..\\',
  '..\\..\\',
  '..\\..\\..\\mpir\\' + library_type + intd,
  '..\\..\\..\\mpfr\\' + library_type + intd,
  '..\\..\\..\\pthreads\\' + library_type + intd
]
libs = [
  '..\\..\\' + library_type + intd + library_type + '_flint',
  '..\\..\\..\\mpir\\' + library_type + intd + 'mpir',
  '..\\..\\..\\mpfr\\' + library_type + intd + 'mpfr',
  '..\\..\\..\\pthreads\\' + library_type + intd + 'pthreads'
]
if (library_type == "lib"):
  macros = [('PTW32_STATIC_LIB',1)]
else:
  macros = [('PTW32_BUILD',1)]

for l2, fp in t:
  fdn, fx = splitext(fp)
  fd, fn = split(fdn)
  if (not build_interfaces_tests and "interface" in fn):
    continue
  source = [join('..\\..\\', fp)]
  p = fd.rfind('test')
  assert p >= 0
  tmp_dir = 'test\\test'
  outd = '..\\tests\\' + fd[:p] + intd
  try:
    obj = cc.compile(source, output_dir=tmp_dir, include_dirs=inc_dirs,macros=macros)
    cc.link("executable", obj, fn + '.exe', output_dir=outd, libraries=libs)
  except:
    error_list += [(l2, fp)]

print('Build Errors:')
for l2, fp in error_list:
  print('    ', l2, fp)
