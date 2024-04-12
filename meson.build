project('FLINT', 'c',
  version: '3.1.0',
  license: 'LGPL3+',
  default_options: [
    'buildtype=release',
    'c_std=c11',
  ],
)

cc = meson.get_compiler('c')
cfg_data = configuration_data()
fs = import('fs')

m_dep = cc.find_library('m', required: false)

gmp_dep = dependency('gmp', version: '>= 6.2.1')
mpfr_dep = dependency('mpfr', version: '>= 4.1.0')

flint_deps = [gmp_dep, mpfr_dep, m_dep]

add_project_arguments(
  '-Werror=implicit-function-declaration',
  language: 'c',
)

subdir('src')

headers_all = []
foreach header : headers
  headers_all += ['src' / header]
endforeach

c_files_all = []
foreach file : c_files_all_subdir
  c_files_all += ['src' / file]
endforeach

subdir('detection')

if have_assembly
  c_files_all += s_files
endif

install_headers(headers_all, subdir: 'flint')

libflint = library('flint', c_files_all,
  dependencies: flint_deps,
  include_directories: include_directories('src'),
  install: true,
  c_args: [
    '-DBUILDING_FLINT',
    '-DFLINT_NOSTDIO',
    '-DFLINT_NOSTDARG',
  ]
)

tests_opt = get_option('tests')
if tests_opt.enabled()
  subdir('src/fmpz/test')
  subdir('src/fmpq/test')
  subdir('src/ulong_extras/test')
else
  bad_test_exe = executable('bad_test', 'bad_test.c')
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
