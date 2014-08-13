'''
Set up Visual Sudio to build a specified MPIR configuration

Copyright (C) 2011, 2012, 2013, 2014 Brian Gladman
'''

from __future__ import print_function
from operator import itemgetter
from os import listdir, walk, unlink, makedirs
from os.path import split, splitext, isdir, relpath, join, exists, dirname
from copy import deepcopy
from sys import argv, exit
from filecmp import cmp
from shutil import copy
from re import compile, search
from collections import defaultdict
from uuid import uuid4
from time import sleep

# for script debugging
debug = False
# what to build
build_lib = False
build_dll = False
build_tests = False
build_profiles = False

# add user choice
flib_type = 'reentrant' # ('gc', 'rentrant', 'single')

# The path to flint, solution and project directories
project_name = 'flint'
build_vc = 'build.vc12'
flint_dir = '../../'
solution_dir = join(flint_dir, build_vc)
project_dir = join(flint_dir, build_vc)

try:
  input = raw_input
except NameError:
  pass

# copy from file ipath to file opath but avoid copying if
# opath exists and is the same as ipath (this is to avoid
# triggering an unecessary rebuild).

def write_f(ipath, opath):
  if exists(ipath) and not isdir(ipath):
    if exists(opath):
      if isdir(opath) or cmp(ipath, opath):
        return
    copy(ipath, opath)

ignore_dirs = ( '.git', 'doc', 'examples')

def find_src(path):
  c, h, cx, hx, t, tx, p = [], [], [], [], [], [], [] 
  for root, dirs, files in walk(path):
    for di in dirs:
      if di in ignore_dirs:
        dirs.remove(di)
    relp = relpath(root, flint_dir)
    if relp == '.':
      relp = ''
    for f in files:
      n, x = splitext(f)
      pth, leaf = split(root)
      fp = join(relp, f)
      if leaf == 'tune':
        continue
      if leaf == 'test':
        p2, l2 = split(pth)
        l2 = '' if l2 == '..' else l2
        if 'flintxx' in pth:
          tx += [(l2, fp)]
        else:
          t += [(l2, fp)]
      elif leaf == 'profile':
        p += [fp]
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
  outf.write(f1)
  for d in sorted(c_dirs):
    if d:
      outf.write(f2.format(d))
  outf.write(f3)

filter_hdr_item = r'    <ClInclude Include="..\..\{}">\n      <Filter>Header Files</Filter>\n    </ClInclude>\n'

def filter_headers(hdr_list, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''    <ClInclude Include="..\..\{}">
    <Filter>Header Files</Filter>
    </ClInclude>
'''
  f3 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for h in hdr_list:
    outf.write(f2.format(h))
  outf.write(f3)

def filter_csrc(cf_list, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''  <ClCompile Include="..\..\{}">
    <Filter>Source Files</Filter>
    </ClCompile>
'''
  f3 = r'''  <ClCompile Include="..\..\{}">
    <Filter>Source Files\{}</Filter>
    </ClCompile>
'''
  f4 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for i in cf_list:
    if not i[0]:
      outf.write(f2.format(i[1]))
    else:
      outf.write(f3.format(i[1], i[0]))
  outf.write(f4)

def gen_filter(name, hf_list, cf_list):

  f1 = r'''<?xml version="1.0" encoding="utf-8"?>
<Project ToolsVersion="12.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
'''
  f2 = r'''
</Project>
'''

  fn = join(project_dir, name)
  try:
    makedirs(split(fn)[0])
  except IOError:
    pass
  with open(fn, 'w') as outf:

    outf.write(f1)
    filter_folders(cf_list, outf)
    filter_headers(hf_list, outf)
    filter_csrc(cf_list, outf)
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

def vcx_library_type(plat, is_dll, outf):

  f1 = r'''  <PropertyGroup Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'" Label="Configuration">
    <ConfigurationType>{2:s}Library</ConfigurationType>
    <PlatformToolset>v120</PlatformToolset>
    </PropertyGroup>
'''
  for pl in plat:
    for conf in ('Release', 'Debug'):
      outf.write(f1.format(pl, conf, 'Dynamic' if is_dll else 'Static'))

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

def vcx_target_name_and_dirs(name, plat, is_dll, outf):

  f1 = r'''  <PropertyGroup>
    <_ProjectFileVersion>10.0.21006.1</_ProjectFileVersion>
'''
  f2 = r'''    <TargetName Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">{2:s}</TargetName>
    <IntDir Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">$(Platform)\$(Configuration)\</IntDir>
    <OutDir Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">$(SolutionDir)..\{3:s}\$(Platform)\$(Configuration)\</OutDir>
'''
  f3 = r'''  </PropertyGroup>
'''
  outf.write(f1)
  for pl in plat:
    for conf in ('Release', 'Debug'):
      outf.write(f2.format(pl, conf, name, ('dll' if is_dll else 'lib')))
  outf.write(f3)

def compiler_options(plat, is_dll, is_debug, inc_dirs, outf):

  f1 = r'''    <ClCompile>
    <Optimization>{0:s}</Optimization>
    <IntrinsicFunctions>true</IntrinsicFunctions>
    <AdditionalIncludeDirectories>{1:s}</AdditionalIncludeDirectories>
    <PreprocessorDefinitions>{2:s}%(PreprocessorDefinitions)</PreprocessorDefinitions>
    <RuntimeLibrary>MultiThreaded{3:s}</RuntimeLibrary>
    <ProgramDataBaseFileName>$(TargetDir)$(TargetName).pdb</ProgramDataBaseFileName>
    <DebugInformationFormat>ProgramDatabase</DebugInformationFormat>
    </ClCompile>
'''

  if is_dll:
    s1 = 'DEBUG;HAVE_CONFIG_H;MSC_BUILD_DLL;'
    s2 = 'DLL'
  else:
    s1 = 'DEBUG;_LIB;HAVE_CONFIG_H;'
    s2 = ''
  if plat == 'x64':
    s1 = s1 + '_WIN64;'
  if is_debug:
    opt, defines, crt = 'Disabled', '_' + s1, 'Debug' + s2
  else:
    opt, defines, crt = 'Full', 'N' + s1, s2
  outf.write(f1.format(opt, inc_dirs, defines, crt))

def linker_options(link_libs, outf):

  f1 = r'''    <Link>
      <GenerateDebugInformation>true</GenerateDebugInformation>
      <LargeAddressAware>true</LargeAddressAware>
      <AdditionalDependencies>{}%(AdditionalDependencies)</AdditionalDependencies>
    </Link>
'''
  outf.write(f1.format(link_libs))

def vcx_post_build(outf):

  f1 = r'''
  <PostBuildEvent>
      <Command>postbuild $(TargetPath)
      </Command>
  </PostBuildEvent>
'''

  outf.write(f1)

def vcx_tool_options(plat, is_dll, postbuild, inc_dirs, link_libs, outf):

  f1 = r'''  <ItemDefinitionGroup Condition="'$(Configuration)|$(Platform)'=='{1:s}|{0:s}'">
'''
  f2 = r'''  </ItemDefinitionGroup>
'''
  for pl in plat:
    for is_debug in (False, True):
      outf.write(f1.format(pl, 'Debug' if is_debug else 'Release'))
      compiler_options(pl, is_dll, is_debug, inc_dirs, outf)
      if is_dll:
        linker_options(link_libs, outf)
      if postbuild:
        vcx_post_build(outf)
      outf.write(f2)

def vcx_hdr_items(hdr_list, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''    <ClInclude Include="..\..\{0:s}" />
'''
  f3 = r'''  </ItemGroup>
'''
  outf.write(f1)
  for i in hdr_list:
    outf.write(f2.format(i))
  outf.write(f3)

def vcx_c_items(cf_list, plat, outf):

  f1 = r'''  <ItemGroup>
'''
  f2 = r'''    <ClCompile Include="..\..\{}" />
'''
  f3 = r'''    <ClCompile Include="..\..\{}">
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
      outf.write(f2.format(nxd[1]))
    else:
      outf.write(f3.format(nxd[1]))
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

def gen_vcxproj(proj_name, file_name, guid, plat, is_dll, postbuild, hf_list, cf_list, inc_dirs, link_libs):

  f1 = r'''<?xml version="1.0" encoding="utf-8"?>
<Project DefaultTargets="Build" ToolsVersion="4.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
'''
  f2 = r'''  <PropertyGroup Label="UserMacros" />
'''
  f3 = r'''  <Import Project="$(VCTargetsPath)\Microsoft.Cpp.targets" />
</Project>
'''

  with open(join(project_dir, file_name), 'w') as outf:
    outf.write(f1)
    vcx_proj_cfg(plat, outf)
    vcx_globals(proj_name, guid, outf)
    vcx_default_cpp_props(outf)
    vcx_library_type(plat, is_dll, outf)
    vcx_cpp_props(outf)
    vcx_user_props(plat, outf)
    outf.write(f2)
    vcx_target_name_and_dirs(proj_name, plat, is_dll, outf)
    vcx_tool_options(plat, is_dll, postbuild, inc_dirs, link_libs, outf)
    vcx_hdr_items(hf_list, outf)
    vcx_c_items(cf_list, plat, outf)
    outf.write(f3)

# add project files to the solution

s_guid = r'\s*(\{\w{8}-\w{4}-\w{4}-\w{4}-\w{12}\})\s*'
s_name = r'\s*\"([a-zA-Z][-.\\_a-zA-Z0-9]*\s*)\"\s*'
re_guid = compile(r'\s*\"\s*' + s_guid + r'\"\s*')
re_proj = compile(r'Project\s*\(\s*\"' + s_guid + r'\"\)\s*=\s*' 
                  + s_name + r'\s*,\s*' + s_name + r'\s*,\s*\"' + s_guid + r'\"')
re_fmap = compile(r'\s*' + s_guid + r'\s*=\s*' + s_guid)

def read_solution_file(soln_name):
  fd, pd, p2f = {}, {}, {}
  solution_path = join(solution_dir, soln_name)
  lines = open(solution_path).readlines()
  for i, ln in enumerate(lines):
    m = re_proj.search(ln)
    if m: 
      if m.group(2) == m.group(3):
        fd[m.group(2)] = (m.group(1), m.group(4))
      elif m.group(3).endswith('.vcxproj') or m.group(3).endswith('.pyproj'): 
        pd[m.group(2)] = (m.group(1), m.group(3), m.group(4))
    m = re_fmap.search(ln)
    if m:
      p2f[m.group(1)] = m.group(2)
  return fd, pd, p2f

sol_1 = '''Microsoft Visual Studio Solution File, Format Version 12.00
# Visual Studio 2013
VisualStudioVersion = 12.0.30626.0
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

def write_solution_file(file_name, fd, pd, p2f):
  with open(join(solution_dir, file_name), 'w') as outf:
    outf.write(sol_1)
    for f, (g1, g2) in fd.items():
      outf.write(sol_2.format(g1, f, f, g2))
    for f, (g1, pn, g2) in pd.items():
      outf.write(sol_2.format(g1, f, pn, g2))
    outf.write(sol_3)
    for f, g in p2f.items():
      outf.write(sol_4.format(f, g))
    outf.write(sol_5)

def add_proj_to_sln(soln_name, soln_folder, proj_name, file_name, guid):
  fd, pd, p2f = read_solution_file(soln_name)
  if soln_folder:
    if soln_folder in fd:
      guid1, guid2 = fd[soln_folder]
    else:
      guid1, guid2 = '{' + str(uuid4()).upper() + '}', '{' + str(uuid4()).upper() + '}'
      fd[soln_folder] = (guid1, guid2)
  pd[proj_name] = ('{' + str(uuid4()).upper() + '}', file_name, guid)
  if soln_folder:
    p2f[guid] = guid2
     
  write_solution_file(soln_name, fd, pd, p2f)

c, h, cx, hx, t, tx, p = find_src(flint_dir)

fn = join(flint_dir, 'fmpz-conversions-{}.in'.format(flib_type))
copy(fn , join(flint_dir, 'fmpz-conversions.h')) 
# fn = join(flint_dir, 'fft_tuning32.in')
fn = join(flint_dir, 'fft_tuning64.in')
copy(fn , join(flint_dir, 'fft_tuning.h')) 
sln_name = project_name + '.sln' 

if build_lib:
  # set up LIB build
  guid = '{' + str(uuid4()).upper() + '}'
  proj_name = 'lib_flint'
  vcx_path = 'lib_flint\\lib_flint.vcxproj'
  gen_filter(vcx_path + '.filters', h, c)
  mode = ('win32', 'x64')
  inc_dirs = r'.\;..\;..\..\;..\..\..\mpir\lib\$(IntDir);..\..\..\mpfr\lib\$(IntDir);..\..\..\pthreads'
  link_libs = r'..\..\..\mpir\lib\$(IntDir)mpir.lib;..\..\..\mpfr\lib\$(IntDir)mpfr.lib;..\..\..\lib\pthreads.lib'
  gen_vcxproj(proj_name, vcx_path, guid, mode, False, True, h, c, inc_dirs, link_libs)
  add_proj_to_sln(sln_name, '', proj_name, vcx_path, guid)

if build_dll:
  # set up DLL build
  guid = '{' + str(uuid4()).upper() + '}'
  proj_name = 'dll_flint'
  vcx_path = 'dll_flint\\dll_flint.vcxproj'
  gen_filter(vcx_path + '.filters', h, c)
  mode = ('win32', 'x64')
  inc_dirs = r'.\;..\;..\..\;..\..\..\mpir\dll\$(IntDir);..\..\..\mpfr\dll\$(IntDir);..\..\..\pthreads'
  link_libs = r'..\..\..\mpir\dll\$(IntDir)mpir.lib;..\..\..\mpfr\dll\$(IntDir)mpfr.lib;..\..\..\lib\pthreads.lib'
  gen_vcxproj(proj_name, vcx_path, guid, mode, True, True, h, c, inc_dirs, link_libs)
  add_proj_to_sln(sln_name, '', proj_name, vcx_path, guid)

if debug:
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
