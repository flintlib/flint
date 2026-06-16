project('FLINT', 'c', 'cpp',
  version: '3.6.0-dev',
  license: 'LGPL3+',
  meson_version: '>=1.1.0',
  default_options: [
    'buildtype=release',
    'c_std=c11',
    'cpp_std=c++11',
  ],
)

FLINT_MAJOR = meson.project_version().split('.')[0]
FLINT_MINOR = meson.project_version().split('.')[1]
FLINT_PATCH = meson.project_version().split('.')[2].split('-')[0]
FLINT_VERSION_FULL = meson.project_version()
FLINT_MAJOR_SO = '24'
FLINT_MINOR_SO = '0'
FLINT_PATCH_SO = '0'
FLINT_LIBRARY_VERSION = '.'.join([FLINT_MAJOR_SO, FLINT_MINOR_SO, FLINT_PATCH_SO])

cc = meson.get_compiler('c')
fs = import('fs')

m_dep = cc.find_library('m', required: false)

gmp_dep = dependency('gmp', version: '>= 6.2.1')
mpfr_dep = dependency('mpfr', version: '>= 4.1.0')

gmp_libdir = gmp_dep.get_variable(pkgconfig: 'libdir')
gmp_link_dep = declare_dependency(
  dependencies: gmp_dep,
  link_args: ['-L' + gmp_libdir],
)

add_project_arguments(
  '-Werror=implicit-function-declaration',
  language: 'c',
)

subdir('config')

flint_deps = [gmp_link_dep, mpfr_dep, m_dep]
if have_blas
  flint_deps += blas_dep
endif
if have_gc
  flint_deps += gc_dep
endif

add_project_arguments(flint_abi_c_args, language: ['c', 'cpp'])
add_project_link_arguments(flint_abi_c_args, language: ['c', 'cpp'])

subdir('src')

if have_assembly
  # Generate config.m4 with assembly language configuration
  config_m4 = configure_file(
    input: 'config.m4.in',
    output: 'config.m4',
    configuration: {'LSYM_PREFIX': LSYM_PREFIX},
  )

  s_files = []
  subdir('src' / assembly_dir)
  c_files_all += s_files
endif

# Build include directory after list of all source files is complete
subdir('include')

libflint = library('flint', c_files_all + headers_all,
  dependencies: flint_deps,
  include_directories: [headers_built_nodir_inc, 'src'],
  install: true,
  version: FLINT_LIBRARY_VERSION,
  soversion: FLINT_MAJOR_SO,
  c_args: [
    '-DBUILDING_FLINT',
    '-DFLINT_NOSTDIO',
    '-DFLINT_NOSTDARG',
  ] + flint_cpu_c_args
)

flint_test_dep = declare_dependency(
  link_with: libflint,
  dependencies: flint_deps,
  include_directories: [headers_built_nodir_inc, 'src'],
)

tests_opt = get_option('tests')
if tests_opt.enabled()
  foreach mod_test : mod_tests
    subdir('src' / mod_test)
  endforeach
else
  bad_test_exe = executable('bad_test', '_meson_build/bad_test.c')
  test('Need to build tests with -Dtests=enabled', bad_test_exe)
endif

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries: libflint,
  version : meson.project_version(),
  name: 'flint',
  filebase: 'flint',
  description: 'Fast Library for Number Theory'
)
