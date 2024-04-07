project('FLINT', 'c',
  version: '3.1.0',
  license: 'LGPL3+'
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)

gmp_dep = dependency('gmp', version: '>= 6.2.1')
mpfr_dep = dependency('mpfr', version: '>= 4.1.0')

flint_deps = [gmp_dep, mpfr_dep, m_dep]

subdir('src')

headers_all = []
foreach header : headers
  headers_all += ['src' / header]
endforeach

c_files_all = []
foreach file : c_files_all_subdir
  c_files_all += ['src' / file]
endforeach

install_headers(headers_all, subdir: 'flint')

libflint = library('flint', c_files_all,
  dependencies: flint_deps,
  include_directories: include_directories('src'),
)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries: libflint,
  version : '3.1.0',
  name: 'flint',
  filebase: 'flint',
  description: 'Fast Library for Number Theory'
)
