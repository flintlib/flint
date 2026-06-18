project('FLINT', 'c',
  version: '3.6.0-dev',
  license: 'LGPL3+',
  meson_version: '>=1.1.0',
  default_options: [
    'buildtype=release',
    'c_std=c11',
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

gmp_link_args = []
if cc.get_argument_syntax() != 'msvc'
  gmp_libdir = gmp_dep.get_variable(pkgconfig: 'libdir')
  gmp_link_args = ['-L' + gmp_libdir]
endif
gmp_link_dep = declare_dependency(
  dependencies: gmp_dep,
  link_args: gmp_link_args,
)

subdir('config')

add_project_arguments(flint_warning_c_args + flint_msvc_c_args, language: 'c')

flint_deps = [gmp_link_dep, mpfr_dep, m_dep, pthread_dep]
if have_blas
  flint_deps += blas_dep
endif
if have_gc
  flint_deps += gc_dep
endif

add_project_arguments(flint_abi_c_args, language: 'c')
add_project_link_arguments(flint_abi_c_args, language: 'c')

subdir('src')

if have_assembly
  # Generate config.m4 with assembly language configuration
  config_m4 = configure_file(
    input: 'config.m4.in',
    output: 'config.m4',
    configuration: asm_cfg,
  )

  s_files = []
  subdir('src' / assembly_dir)
  c_files_all += s_files
endif

# Build include directory after list of all source files is complete
subdir('include')

libflint_c_args = [
  '-DBUILDING_FLINT',
  '-DFLINT_NOSTDIO',
  '-DFLINT_NOSTDARG',
]

regression_c_args = libflint_c_args + flint_cpu_c_args + cc.get_supported_arguments([
  '-march=native',
  '-Wall',
  '-Wextra',
  '-Werror',
])

regression_lib = static_library('flint_regression_check',
  regression_c_files + headers_all,
  dependencies: flint_deps,
  include_directories: [headers_built_nodir_inc, 'src'],
  install: false,
  build_by_default: false,
  c_args: regression_c_args,
)

alias_target('regression-check', regression_lib)

if cc.get_argument_syntax() == 'msvc' and get_option('default_library') != 'static'
  libflint_objects = static_library('flint_objects', c_files_all + headers_all,
    dependencies: flint_deps,
    include_directories: [headers_built_nodir_inc, 'src'],
    install: false,
    c_args: libflint_c_args + ['-DFLINT_BUILD_DLL'] + flint_cpu_c_args,
  )

  libflint_object_files = libflint_objects.extract_all_objects(recursive: false)

  libflint_def = custom_target(
    'flint.def',
    output: 'flint.def',
    command: [
      python,
      files('config/generate_msvc_def.py'),
      '@OUTPUT@',
      '--object-dir',
      'libflint_objects.a.p',
    ],
    depends: libflint_objects,
  )

  libflint = shared_library('flint',
    objects: libflint_object_files,
    dependencies: flint_deps,
    include_directories: [headers_built_nodir_inc, 'src'],
    install: true,
    version: FLINT_LIBRARY_VERSION,
    soversion: FLINT_MAJOR_SO,
    vs_module_defs: libflint_def,
  )
else
  libflint = library('flint', c_files_all + headers_all,
    dependencies: flint_deps,
    include_directories: [headers_built_nodir_inc, 'src'],
    install: true,
    version: FLINT_LIBRARY_VERSION,
    soversion: FLINT_MAJOR_SO,
    c_args: libflint_c_args + flint_cpu_c_args,
  )
endif

flint_test_dep = declare_dependency(
  link_with: libflint,
  dependencies: flint_deps,
  include_directories: [headers_built_nodir_inc, 'src'],
  sources: headers_all,
  compile_args: flint_cpu_c_args,
)

examples_opt = get_option('examples')
if not examples_opt.disabled()
  subdir('examples')
endif

tests_opt = get_option('tests')
if not tests_opt.disabled()
  foreach mod_test : mod_tests
    subdir('src' / mod_test)
  endforeach
endif

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries: libflint,
  version : meson.project_version(),
  name: 'flint',
  filebase: 'flint',
  description: 'Fast Library for Number Theory'
)
