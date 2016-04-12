'''
Set up Visual Sudio to build a specified MPIR configuration

Copyright (C) 2011, 2012, 2013, 2014 Brian Gladman
'''

from __future__ import print_function
from operator import itemgetter
from os import listdir, walk, unlink, makedirs, remove
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

import argparse
parser = argparse.ArgumentParser(description='Configure flint msvc build')

# for script debugging
parser.add_argument('--debug', choices=["True", "False"], default="False")

# what to build
parser.add_argument('--build-lib', choices=["True", "False"], default="False")
parser.add_argument('--build-dll', choices=["True", "False"], default="False")
parser.add_argument('--build-tests', choices=["True", "False"], default="True")
parser.add_argument('--build-profiles', choices=["True", "False"], default="True")

# add user choice
parser.add_argument('--flib-type', choices=('gc', 'reentrant', 'single'), default="single")

args = parser.parse_args()
print(args)

debug = args.debug == "True"
build_lib = args.build_lib == "True"
build_dll = args.build_dll == "True"
build_tests = args.build_tests == "True"
build_profiles = args.build_profiles == "True"

flib_type = args.flib_type

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

ignore_dirs = ( '.git', 'doc', 'examples', 'lib', 'exe', 'dll', 'win_hdrs', 'build.vc12')
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

s_sym = r'([_a-zA-Z][_a-zA-Z0-9]*)'
dec = '^(?:' + s_sym + '\s+)?' + s_sym + '\('
re_dec = compile(dec)
re_end = compile('\s*\);\s*$')

# crude parser to detect function declarations in
# FLINT header files (no longer needed but left in
# place in case this changes).

def strip_static_inline(lines, pos, lth):
  p0 = pos
  m = re_end.search(lines[pos])
  pos += 1
  if m:
    return pos
  level = 0
  while pos < lth:
    line = (lines[pos].replace('\n', ' ')).strip()
    m = re_end.search(lines[pos])
    if not level and (m or not line or pos > p0 + 5) or pos >= lth:
      return pos + 1
    level += line.count('{') - line.count('}')
    pos += 1
  return pos

def parse_hdrs(h):
  d = defaultdict(list)
  for hf in h:
    lines = open(join(flint_dir, hf), 'r').readlines()
    pos, n_lines = 0, len(lines)
    line_buf = ''
    in_dec = 0
    while pos < n_lines:
      if in_dec == 0:
        line_buf = ''
      line = (lines[pos].replace('\n', ' ')).strip()
      if not line:
        pos += 1
        continue
      if line.startswith('#define'):
        while line.endswith('\\') and pos < n_lines:
          pos += 1
          line = (lines[pos].replace('\n', ' ')).strip()
        pos += 1
        continue
      if 'INLINE' in line.upper():
        pos = strip_static_inline(lines, pos, n_lines)
        in_dec = 0
        continue
      if not in_dec:
        m = re_dec.search(line)
        if m:
          line_buf += line + ' '
          mm = re_end.search(line)
          if mm:
            pos += 1
            print(m.group(2))
            d[hf] += [(pos, m.group(2), line_buf)]
            in_dec = 0
            continue
          in_dec += 1
      else:
        line_buf += line
        mm = re_end.search(line)
        if mm:
          if in_dec < 5:
            d[hf] += [(pos - in_dec + 1, m.group(2), line_buf)]
          in_dec = 0
        else:
          in_dec += 1
      pos += 1
  return d

def write_hdrs(h):
  d = parse_hdrs(h)
  hdr_dir = join(flint_dir, 'win_hdrs')
  if not exists(hdr_dir):
    makedirs(hdr_dir)

  for hf in sorted(d.keys()):
    out_name = hf.replace('build.vc14\\', '')
    inf = open(join(flint_dir, hf), 'r')
    outf = open(join(flint_dir, join(hdr_dir, out_name)), 'w')
    lines = inf.readlines()
    for n, _p, _d in d[hf]:
      lines[n - 1] = 'FLINT_DLL ' + lines[n - 1]
    outf.writelines(lines)
    inf.close()
    outf.close()

def write_def_file(name, h):
  d = parse_hdrs(h)
  lines = ['LIBRARY ' + name + '\n', 'EXPORTS' + '\n']
  for hf in sorted(d.keys()):
    for _n, sym, _d in d[hf]:
      lines += ['    ' + sym + '\n']
  with open(join(solution_dir, name + '.def'), 'w') as outf:
    outf.writelines(lines)

# end of parser code

def filter_folders(cf_list, outf):

  f1 = r'''  <ItemGroup>
    <Filter Include="Header Files" />
    <Filter Include="Source Files" />
'''
  f2 = r'''    <Filter Include="Source Files\{0:s}" />
'''
  f3 = r'''  </ItemGroup>
'''
  c_dirs = set(i[0] for i in cf_list)
  if c_dirs:
    outf.write(f1)
    for d in sorted(c_dirs):
      if d:
        outf.write(f2.format(d))
    outf.write(f3)

def filter_headers(hdr_list, relp, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''    <ClInclude Include="{}{}">
    <Filter>Header Files</Filter>
    </ClInclude>
'''
  f3 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for h in hdr_list:
    outf.write(f2.format(relp, h))
  outf.write(f3)

def filter_csrc(cf_list, relp, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''  <ClCompile Include="{}{}">
    <Filter>Source Files</Filter>
    </ClCompile>
'''
  f3 = r'''  <ClCompile Include="{}{}">
    <Filter>Source Files\{}</Filter>
    </ClCompile>
'''
  f4 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for i in cf_list:
    if not i[0]:
      outf.write(f2.format(relp, i[1]))
    else:
      outf.write(f3.format(relp, i[1], i[0]))
  outf.write(f4)

def gen_filter(name, hf_list, cf_list):

  f1 = r'''<?xml version="1.0" encoding="utf-8"?>
<Project ToolsVersion="12.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
'''
  f2 = r'''
</Project>
'''
  fn = normpath(join(solution_dir, name))
  relp = split(relpath(flint_dir, fn))[0] + '\\'
  try:
    makedirs(split(fn)[0])
  except IOError:
    pass
  with open(fn, 'w') as outf:

    outf.write(f1)
    filter_folders(cf_list, outf)
    if hf_list:
      filter_headers(hf_list, relp, outf)
    filter_csrc(cf_list, relp, outf)
    outf.write(f2)

# generate vcxproj file

def vcx_proj_cfg(plat, outf):

  f1 = r'''  <ItemGroup Label="ProjectConfigurations">
'''
  f2 = r'''    <ProjectConfiguration Include="{1:s}|{0:s}">
    <Configuration>{1:s}</Configuration>
    <Platform>{0:s}</Platform>
    </ProjectConfiguration>
'''
  f3 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for pl in plat:
    for conf in ('Release', 'Debug'):
      outf.write(f2.format(pl, conf))
  outf.write(f3)

def vcx_globals(name, guid, outf):

  f1 = r'''  <PropertyGroup Label="Globals">
    <RootNamespace>{0:s}</RootNamespace>
    <Keyword>Win32Proj</Keyword>
    <ProjectGuid>{1:s}</ProjectGuid>
    </PropertyGroup>
'''
  outf.write(f1.format(name, guid))

def vcx_default_cpp_props(outf):

  f1 = r'''  <Import Project="$(VCTargetsPath)\Microsoft.Cpp.Default.props" />
'''
  outf.write(f1)

def vcx_library_type(plat, proj_type, outf):

  f1 = r'''  <PropertyGroup Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'" Label="Configuration">
    <ConfigurationType>{2:s}</ConfigurationType>
    <PlatformToolset>v140</PlatformToolset>
    </PropertyGroup>
'''
  for pl in plat:
    for conf in ('Release', 'Debug'):
      outf.write(f1.format(pl, conf, app_str[proj_type]))

def vcx_cpp_props(outf):

  f1 = r'''  <Import Project="$(VCTargetsPath)\Microsoft.Cpp.props" />
'''
  outf.write(f1)

def vcx_user_props(plat, outf):

  f1 = r'''  <ImportGroup Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'" Label="PropertySheets">
    <Import Project="$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props" Condition="exists('$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props')" />
    </ImportGroup>
'''
  for pl in plat:
    for conf in ('Release', 'Debug'):
      outf.write(f1.format(pl, conf))

def vcx_target_name_and_dirs(proj_name, proj_dir, plat, outf):

  f1 = r'''  <PropertyGroup>
    <_ProjectFileVersion>10.0.21006.1</_ProjectFileVersion>
'''
  f2 = r'''    <TargetName Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">{2:s}</TargetName>
    <IntDir Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">$(Platform)\$(Configuration)\</IntDir>
    <OutDir Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">$(SolutionDir){3:s}$(Platform)\$(Configuration)\</OutDir>
    <LinkIncremental>false</LinkIncremental>
'''
  f3 = r'''  </PropertyGroup>
'''
  if not proj_dir:
    proj_dir = ''
  elif not (proj_dir.endswith('\\') or proj_dir.endswith('/')):
    proj_dir += '\\'

  outf.write(f1)
  for pl in plat:
    for conf in ('Release', 'Debug'):
      outf.write(f2.format(pl, conf, proj_name, proj_dir))
  outf.write(f3)

def compiler_options(plat, proj_type, is_debug, inc_dirs, outf):

  f1 = r'''    <ClCompile>
    <Optimization>{0:s}</Optimization>
    <IntrinsicFunctions>true</IntrinsicFunctions>
    <AdditionalIncludeDirectories>{1:s}</AdditionalIncludeDirectories>
    <PreprocessorDefinitions>{2:s}%(PreprocessorDefinitions)</PreprocessorDefinitions>
    <RuntimeLibrary>MultiThreaded{3:s}</RuntimeLibrary>
    <ProgramDataBaseFileName>$(TargetDir)$(TargetName).pdb</ProgramDataBaseFileName>
    <DebugInformationFormat>{4:s}</DebugInformationFormat>
    </ClCompile>
'''

  if proj_type == app_type:
    s1 = 'DEBUG;WIN32;_CONSOLE;PTW32_STATIC_LIB;'
    s2 = ''
  if proj_type == dll_type:
    s1 = 'DEBUG;WIN32;HAVE_CONFIG_H;MSC_BUILD_DLL;PTW32_BUILD;'
    s2 = 'DLL'
  elif proj_type == lib_type:
    s1 = 'DEBUG;WIN32;_LIB;HAVE_CONFIG_H;PTW32_STATIC_LIB;'
    s2 = ''
  else:
    pass
  if flib_type == 'single':
    s1 += 'FLINT_REENTRANT=0;HAVE_TLS=1;'
  elif flib_type == 'reentrant':
    s1 += 'FLINT_REENTRANT=1;'
  else:
    pass
  if plat == 'x64':
    s1 = s1 + '_WIN64;'
  if is_debug:
    opt, defines, crt, dbf = 'Disabled', '_' + s1, 'Debug' + s2, 'ProgramDatabase'
  else:
    opt, defines, crt, dbf = 'Full', 'N' + s1, s2, 'None'
  outf.write(f1.format(opt, inc_dirs, defines, crt, dbf))

def linker_options(name, link_libs, proj_type, debug_info, outf):

  f1 = r'''    <Link>
'''
  f2 = r'''      <GenerateDebugInformation>true</GenerateDebugInformation>
'''
  f3 = r'''      <LargeAddressAware>true</LargeAddressAware>
      <AdditionalDependencies>{}%(AdditionalDependencies)</AdditionalDependencies>
'''
  f4 = '''      <ModuleDefinitionFile>$(SolutionDir){}.def</ModuleDefinitionFile>
'''
  f5 = '''    </Link>
'''
  outf.write(f1)
  if debug_info:
    outf.write(f2)
  outf.write(f3.format(link_libs))
  # no longer needed as we are using the declspec approach
  #  if proj_type == dll_type:
  #    outf.write(f4.format(name))
  outf.write(f5)

def vcx_pre_build(outf):

  f1 = r'''    <PreBuildEvent>
        <Command>..\out_copy_rename.bat ..\cpimport.h ..\..\qadic\ cpimport.h
        </Command>
    </PreBuildEvent>
'''

  outf.write(f1)

def vcx_post_build(outf, proj_type):

  f1 = r'''
  <PostBuildEvent>
      <Command>postbuild $(IntDir) {}
      </Command>
  </PostBuildEvent>
'''

  outf.write(f1.format((app_ext[proj_type][1:]).upper()))

def vcx_tool_options(name, plat, proj_type, prebuild, postbuild, inc_dirs, link_libs, debug_info, outf):

  f1 = r'''  <ItemDefinitionGroup Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">
'''
  f2 = r'''  </ItemDefinitionGroup>
'''
  for pl in plat:
    for is_debug in (False, True):
      outf.write(f1.format(pl, 'Debug' if is_debug else 'Release'))
      if prebuild:
        vcx_pre_build(outf)
      compiler_options(pl, proj_type, is_debug, inc_dirs, outf)
      if proj_type != lib_type:
        linker_options(name, link_libs, proj_type, debug_info, outf)
      if postbuild:
        vcx_post_build(outf, proj_type)
      outf.write(f2)

def vcx_hdr_items(hdr_list, relp, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''    <ClInclude Include="{}{}" />
'''
  f3 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for i in hdr_list:
    outf.write(f2.format(relp, i))
  outf.write(f3)

def vcx_c_items(cf_list, plat, relp, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''    <ClCompile Include="{}{}" />
'''
  f3 = r'''    <ClCompile Include="{}{}">
'''
  f4 = r'''        <ObjectFileName Condition="'$(Configuration)|$(Platform)'=='{0:s}|{1:s}'">$(IntDir){2:s}\</ObjectFileName>
'''
  f5 = r'''      <ExcludedFromBuild Condition="'$(Configuration)|$(Platform)'=='{0:s}|{1:s}'">true</ExcludedFromBuild>
'''
  f6 = r'''    </ClCompile>
'''
  f7 = r'''  </ItemGroup>
'''

  outf.write(f1)
  for nxd in cf_list:
    if nxd[0] in ('', build_vc):
      outf.write(f2.format(relp, nxd[1]))
    else:
      outf.write(f3.format(relp, nxd[1]))
      for cf in ('Release', 'Debug'):
        for pl in plat:
          outf.write(f4.format(cf, pl, nxd[0]))
      if nxd[0] == 'link':
        ts = split(nxd[1])[1]
        ts = ts.replace('fmpz_', '')
        ts = ts.replace('.c', '')
        if ts != flib_type:
          for cf in ('Release', 'Debug'):
            for pl in plat:
              outf.write(f5.format(cf, pl))
      outf.write(f6)
  outf.write(f7)

def gen_vcxproj(proj_name, project_dir, file_name, guid, plat, proj_type, prebuild, postbuild, debug_info, hf_list, cf_list, inc_dirs, link_libs):

  f1 = r'''<?xml version="1.0" encoding="utf-8"?>
<Project DefaultTargets="Build" ToolsVersion="14.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
'''
  f2 = r'''  <PropertyGroup Label="UserMacros" />
'''
  f3 = r'''  <Import Project="$(VCTargetsPath)\Microsoft.Cpp.targets" />
</Project>
'''

  fn = normpath(join(solution_dir, file_name))
  relp = split(relpath(flint_dir, fn))[0] + '\\'
  with open(fn, 'w') as outf:
    outf.write(f1)
    vcx_proj_cfg(plat, outf)
    vcx_globals(proj_name, guid, outf)
    vcx_default_cpp_props(outf)
    vcx_library_type(plat, proj_type, outf)
    vcx_cpp_props(outf)
    vcx_user_props(plat, outf)
    outf.write(f2)
    vcx_target_name_and_dirs(proj_name, project_dir, plat, outf)
    vcx_tool_options(proj_name, plat, proj_type, prebuild, postbuild, inc_dirs, link_libs, debug_info, outf)
    if hf_list:
      vcx_hdr_items(hf_list, relp, outf)
    vcx_c_items(cf_list, plat, relp, outf)
    outf.write(f3)

# add project files to the solution

folder_guid = "{2150E333-8FDC-42A3-9474-1A3956D46DE8}"
vcxproj_guid = "{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}"
pyproj_guid =  "{888888A0-9F3D-457C-B088-3A5042F75D52}"

s_guid = r'\s*(\{\w{8}-\w{4}-\w{4}-\w{4}-\w{12}\})\s*'
s_name = r'\s*\"([a-zA-Z][-.\\_a-zA-Z0-9]*\s*)\"\s*'
re_guid = compile(r'\s*\"\s*' + s_guid + r'\"\s*')
re_proj = compile(r'Project\s*\(\s*\"' + s_guid + r'\"\)\s*=\s*'
                  + s_name + r'\s*,\s*' + s_name + r'\s*,\s*\"' + s_guid + r'\"')
re_fmap = compile(r'\s*' + s_guid + r'\s*=\s*' + s_guid)

def read_solution_file(soln_name):
  g2fldr, g2proj = {}, {}
  gf2gpl = defaultdict(list)
  solution_path = join(solution_dir, soln_name)
  if exists(solution_path):
    lines = open(solution_path).readlines()
    for i, ln in enumerate(lines):
      m = re_proj.search(ln)
      if m:
        if m.group(1) == folder_guid and m.group(2) == m.group(3):
          # folder guid -> folder name
          g2fldr[m.group(4)] = m.group(2)
        elif (m.group(1) == vcxproj_guid and m.group(3).endswith('.vcxproj') or
              m.group(1) == pyproj_guid and m.group(3).endswith('.pyproj')):
          # project guid -> proj_type_guid, proj_name, vcx_path 
          g2proj[m.group(4)] = (m.group(1), m.group(2), m.group(3))

    for i, ln in enumerate(lines):
      m = re_fmap.search(ln)
      if m:
        if m.group(1) in g2proj and m.group(2) in g2fldr:
          gf2gpl[m.group(2)].append(m.group(1))

  for g in g2proj:
    for _, gpl in gf2gpl.items():
      if g in gpl:
        break
    else:
      gf2gpl[''].append(g)

  assert len(g2fldr.keys()) == len(gf2gpl.keys()) - (1 if '' in gf2gpl.keys() else 0)
  assert sum(len(gpl) for gf, gpl in gf2gpl.items()) == len(g2proj.keys())
  return g2fldr, g2proj, gf2gpl

sol_1 = '''Microsoft Visual Studio Solution File, Format Version 12.00
# Visual Studio 14
VisualStudioVersion = 14.0.23107.0
MinimumVisualStudioVersion = 10.0.40219.1
'''

sol_2 = '''Project("{}") = "{}", "{}", "{}"
EndProject
'''

sol_3 = '''Global
	GlobalSection(SolutionConfigurationPlatforms) = preSolution
		Debug|Win32 = Debug|Win32
		Debug|x64 = Debug|x64
		Release|Win32 = Release|Win32
		Release|x64 = Release|x64
	EndGlobalSection
	GlobalSection(SolutionProperties) = preSolution
		HideSolutionNode = FALSE
	EndGlobalSection
	GlobalSection(NestedProjects) = preSolution
'''

sol_4 = '''		{} = {}
'''

sol_5 = r'''	EndGlobalSection
EndGlobal
'''

def write_solution_file(soln_name, g2fldr, g2proj, gf2gpl):
  if len(g2fldr.keys()) > len(gf2gpl.keys()):
    for g in list(g2fldr.keys()):
      if g not in gf2gpl.keys():
        del g2fldr[g]
  assert len(g2fldr.keys()) == len(gf2gpl.keys()) - (1 if '' in gf2gpl.keys() else 0)
  assert sum(len(gpl) for gf, gpl in gf2gpl.items()) == len(g2proj.keys())

  with open(join(solution_dir, soln_name), 'w') as outf:
    outf.write(sol_1)
    for g, f in g2fldr.items():
      outf.write(sol_2.format(folder_guid, f, f, g))
    for g, (gg, f, n) in g2proj.items():
      outf.write(sol_2.format(gg, f, n, g))
    outf.write(sol_3)
    for gf, gpl in gf2gpl.items():
      if gf != '':
        for gp in gpl:
          outf.write(sol_4.format(gp, gf))
    outf.write(sol_5)

def add_proj_to_sln(soln_name, soln_folder, proj_name, file_name, p_guid, sol_inf):
  if not sol_inf:
    g2fldr, g2proj, gf2gpl = read_solution_file(soln_name)
  else:
    g2fldr, g2proj, gf2gpl = sol_inf
  if soln_folder:
    for g, f in g2fldr.items():
      if f == soln_folder:
        f_guid = g
        break
    else:
      f_guid = '{' + str(uuid4()).upper() + '}'
      g2fldr[f_guid] = soln_folder

  for g in list(g2proj.keys()):
    if g2proj[g] == (vcxproj_guid, proj_name, file_name):
      del g2proj[g]
      for _, gpl in gf2gpl.items():
        if g in gpl:
          del gpl[gpl.index(g)]
      break
  g2proj[p_guid] = (vcxproj_guid, proj_name, file_name)
  gf2gpl[f_guid if soln_folder else ''].append(p_guid)
  if not sol_inf:
    write_solution_file(soln_name, g2fldr, g2proj, gf2gpl)

c, h, cx, hx, t, tx, p = find_src(flint_dir)

# write_hdrs(h)

if not debug:
  with open('..\\..\\qadic\\CPimport.txt', 'r') as fin:
    with open('tmp.h', 'w') as fout:
      while True:
        l = fin.readline()
        if not l:
          break
        l = l.replace(' ', ',')
        l = l.replace('\n', ',\n')
        fout.writelines([l])
  write_f('tmp.h', '..\\cpimport.h')
  remove('tmp.h')

  fn = join(flint_dir, 'fmpz-conversions-{}.in'.format(flib_type))
  copy(fn , join(flint_dir, 'fmpz-conversions.h'))
  # fn = join(flint_dir, 'fft_tuning32.in')
  fn = join(flint_dir, 'fft_tuning64.in')
  copy(fn , join(flint_dir, 'fft_tuning.h'))
  sln_name = project_name + '.sln'
  # write_hdrs(h)

if build_lib:
  # set up LIB build
  guid = '{' + str(uuid4()).upper() + '}'
  proj_name = 'lib_flint'
  vcx_path = 'lib_flint\\lib_flint.vcxproj'
  gen_filter(vcx_path + '.filters', h, c)
  mode = ('Win32', 'x64')
  inc_dirs = r'..\;..\..\;..\..\..\mpir\lib\$(IntDir);..\..\..\mpfr\lib\$(IntDir);..\..\..\pthreads\lib\$(IntDir)'
  link_libs = r'..\..\..\mpir\lib\$(IntDir)mpir.lib;..\..\..\mpfr\lib\$(IntDir)mpfr.lib;..\..\..\pthreads\lib\$(IntDir)pthreads.lib'
  gen_vcxproj(proj_name, None, vcx_path, guid, mode, lib_type, True, True, True, h, c, inc_dirs, link_libs)
  add_proj_to_sln(sln_name, '', proj_name, vcx_path, guid, None)

if build_dll:
  # set up DLL build

  # no longer needed
  # write_def_file('dll_flint', h)

  guid = '{' + str(uuid4()).upper() + '}'
  proj_name = 'dll_flint'
  vcx_path = 'dll_flint\\dll_flint.vcxproj'
  gen_filter(vcx_path + '.filters', h, c)
  mode = ('Win32', 'x64')
  inc_dirs = r'..\;..\..\;..\..\..\mpir\dll\$(IntDir);..\..\..\mpfr\dll\$(IntDir);..\..\..\pthreads\dll\$(IntDir);'
  link_libs = r'..\..\..\mpir\dll\$(IntDir)mpir.lib;..\..\..\mpfr\dll\$(IntDir)mpfr.lib;..\..\..\pthreads\dll\$(IntDir)pthreads.lib;'
  gen_vcxproj(proj_name, None, vcx_path, guid, mode, dll_type, True, True, True, h, c, inc_dirs, link_libs)
  add_proj_to_sln(sln_name, '', proj_name, vcx_path, guid, None)

def gen_test(sln_name, test_name, directory, proj_dir, name, c_file):
  # set up LIB build
  guid = '{' + str(uuid4()).upper() + '}'
  proj_name = test_name[2:]
  p_dir = join(proj_dir, name)
  vcx_path = directory + '\\' + name + '_' + proj_name + '\\' + proj_name + '.vcxproj'
  gen_filter(vcx_path + '.filters', [], [('', c_file)])
  mode = ('Win32', 'x64')
  inc_dirs = r'..\..\;..\..\..\;..\..\..\..\mpir\lib\$(IntDir);..\..\..\..\mpfr\lib\$(IntDir);..\..\..\..\pthreads\lib\$(IntDir);'
  link_libs = r'..\..\..\lib\$(IntDir)lib_flint.lib;..\..\..\..\mpir\lib\$(IntDir)mpir.lib;..\..\..\..\mpfr\lib\$(IntDir)mpfr.lib;..\..\..\..\pthreads\lib\$(IntDir)pthreads.lib;'
  gen_vcxproj(test_name, p_dir, vcx_path, guid, mode, app_type, False, False, False, [], [('', c_file)], inc_dirs, link_libs)
  return vcx_path, guid

def gen_tests(sln_name, directory, proj_dir, c_files):
  sn = sln_name[:sln_name.rfind('.')]
  cnt = -1
  for name, fpath in c_files:
    cnt += 1
    soln = sn + str(cnt // 100 + 1) + '.sln'
    if cnt % 100 == 0:
      print(soln)
      g2fldr, g2proj, gf2gpl = read_solution_file(soln)
    print('  ', cnt, name, fpath)
    _, t = split(fpath)
    if t[:2] not in ('t-', 'p-'):
      continue
    p_name = t.replace('.c', '')
    vcx_name, p_guid = gen_test(sln_name, p_name, directory, proj_dir, name, fpath)
    add_proj_to_sln(soln, name, p_name, vcx_name, p_guid, (g2fldr, g2proj, gf2gpl))
    if cnt % 100 == 99:
      write_solution_file(soln, g2fldr, g2proj, gf2gpl)
  if cnt % 100:
    write_solution_file(soln, g2fldr, g2proj, gf2gpl)

if build_tests:
  gen_tests('flint-tests.sln', 'flint-tests', 'tests', t)

if build_profiles:
  gen_tests('flint-profiles.sln', 'flint-profiles', 'profiles', p)

if debug:
  print('\nC files')
  for d in c:
    print('  ', d)
  print('\nC header files')
  for d in h:
    print('  ', d)
  print('\nC++ files')
  for d in cx:
    print('  ', d)
  print('\nC++ header files')
  for d in hx:
    print('  ', d)
  print('\nC Test files')
  for d in t:
    print('  ', d)
  print('\nC++ Test files')
  for d in tx:
    print('  ', d)
  print('\nProfile files')
  for d in p:
    print('  ', d)
