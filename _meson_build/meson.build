project('FLINT', 'c',
  version: '3.2.0-dev',
  license: 'LGPL3+',
  default_options: [
    'buildtype=release',
    'c_std=c11',
  ],
)

FLINT_MAJOR = meson.project_version().split('.')[0]
FLINT_MINOR = meson.project_version().split('.')[1]
FLINT_PATCH = meson.project_version().split('.')[2]
FLINT_VERSION_FULL = meson.project_version()

cfg_data = configuration_data()
cfg_data.set('FLINT_MAJOR', FLINT_MAJOR)
cfg_data.set('FLINT_MINOR', FLINT_MINOR)
cfg_data.set('FLINT_PATCH', FLINT_PATCH)
cfg_data.set('FLINT_VERSION_FULL', FLINT_VERSION_FULL)

cc = meson.get_compiler('c')
fs = import('fs')

m_dep = cc.find_library('m', required: false)

gmp_dep = dependency('gmp', version: '>= 6.2.1')
mpfr_dep = dependency('mpfr', version: '>= 4.1.0')

flint_deps = [gmp_dep, mpfr_dep, m_dep]

add_project_arguments(
  '-Werror=implicit-function-declaration',
  language: 'c',
)

subdir('detection')
subdir('src')

if have_assembly
  s_files = []
  subdir('src' / assembly_dir)
  c_files_all += s_files
endif

# Build include directory after list of all source files is complete
subdir('include')

install_headers(headers_all, subdir: 'flint')

libflint = library('flint', c_files_all,
  dependencies: flint_deps,
  include_directories: [headers_built_nodir_inc, 'src'],
  install: true,
  c_args: [
    '-DBUILDING_FLINT',
    '-DFLINT_NOSTDIO',
    '-DFLINT_NOSTDARG',
  ]
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
  version : '3.1.0',
  name: 'flint',
  filebase: 'flint',
  description: 'Fast Library for Number Theory'
)
