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
python3 = find_program('python3')
default_library_opt = get_option('default_library')
build_shared_library = default_library_opt in ['shared', 'both']
build_static_library = default_library_opt in ['static', 'both']

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

flint_compile_deps = []
foreach dep : flint_deps
  flint_compile_deps += dep.partial_dependency(
    compile_args: true,
    includes: true,
    sources: true,
  )
endforeach

add_project_arguments(flint_abi_c_args, language: 'c')
add_project_link_arguments(flint_abi_c_args, language: 'c')

subdir('src')

meson.add_postconf_script(
  python3,
  files('dev/create_build_links.py'),
  meson.project_build_root(),
)

if have_assembly
  # Generate config.m4 with assembly language configuration
  config_m4 = configure_file(
    input: 'config.m4.in',
    output: 'config.m4',
    configuration: asm_cfg,
  )

  s_files_static = []
  s_files_shared = []
  subdir('src' / assembly_dir)
else
  s_files_static = []
  s_files_shared = []
endif

c_files_static_all = c_files_all + s_files_static
c_files_shared_all = c_files_all + s_files_shared
c_files_static = s_files_static
c_files_shared = s_files_shared

# Build include directory after list of all source files is complete
subdir('include')

libflint_c_args = [
  '-DBUILDING_FLINT',
  '-DFLINT_NOSTDIO',
  '-DFLINT_NOSTDARG',
]

libflint_shared_link_whole = []
libflint_static_link_whole = []

if cc.get_argument_syntax() != 'msvc'
  foreach bucket : [
    ['missing_prototypes',
      c_file_buckets_missing_prototypes,
      libflint_c_args + flint_cpu_c_args + flint_library_warning_c_args],
    ['no_missing_prototypes',
      c_file_buckets_no_missing_prototypes,
      libflint_c_args + flint_cpu_c_args],
    ['unroll',
      c_file_buckets_missing_prototypes_unroll,
      libflint_c_args + flint_cpu_c_args + flint_library_warning_c_args + flint_unroll_loops_c_args],
    ['no_missing_prototypes_unroll',
      c_file_buckets_no_missing_prototypes_unroll,
      libflint_c_args + flint_cpu_c_args + flint_unroll_loops_c_args],
  ]
    bucket_name = bucket[0]
    bucket_source_groups = bucket[1]
    bucket_c_args = bucket[2]

    foreach bucket_source_group : bucket_source_groups
      bucket_source_name = bucket_source_group[0]
      bucket_sources = bucket_source_group[1]

      if bucket_sources.length() == 0
        continue
      endif

      if build_shared_library
        libflint_shared_link_whole += [static_library('flint_' + bucket_name + '_' + bucket_source_name + '_shared',
          bucket_sources + headers_all,
          dependencies: flint_compile_deps,
          include_directories: [headers_built_nodir_inc, 'src'],
          install: false,
          pic: true,
          c_args: bucket_c_args,
        )]
      endif

      if build_static_library
        libflint_static_link_whole += [static_library('flint_' + bucket_name + '_' + bucket_source_name + '_static',
          bucket_sources + headers_all,
          dependencies: flint_compile_deps,
          include_directories: [headers_built_nodir_inc, 'src'],
          install: false,
          pic: false,
          c_args: bucket_c_args,
        )]
      endif
    endforeach
  endforeach
endif

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
  override_options: ['optimization=0'],
)

alias_target('regression-check', regression_lib)

if cc.get_argument_syntax() == 'msvc' and build_shared_library
  libflint_objects = static_library('flint_objects', c_files_shared_all + headers_all,
    dependencies: flint_deps,
    include_directories: [headers_built_nodir_inc, 'src'],
    install: false,
    c_args: libflint_c_args + ['-DFLINT_BUILD_DLL'] + flint_cpu_c_args + flint_library_warning_c_args,
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

  libflint_shared = shared_library('flint',
    objects: libflint_object_files,
    dependencies: flint_deps,
    include_directories: [headers_built_nodir_inc, 'src'],
    install: true,
    version: FLINT_LIBRARY_VERSION,
    soversion: FLINT_MAJOR_SO,
    vs_module_defs: libflint_def,
  )

  if build_static_library
    libflint_static = static_library('flint', c_files_static_all + headers_all,
      dependencies: flint_deps,
      include_directories: [headers_built_nodir_inc, 'src'],
      install: true,
      c_args: libflint_c_args + flint_cpu_c_args + flint_library_warning_c_args,
    )
  endif

  libflint = libflint_shared
elif cc.get_argument_syntax() == 'msvc'
  libflint = static_library('flint', c_files_static_all + headers_all,
    dependencies: flint_deps,
    include_directories: [headers_built_nodir_inc, 'src'],
    install: true,
    c_args: libflint_c_args + flint_cpu_c_args + flint_library_warning_c_args,
  )
else
  if build_shared_library
    libflint_shared = shared_library('flint', c_files_shared + headers_all,
      dependencies: flint_deps,
      include_directories: [headers_built_nodir_inc, 'src'],
      link_whole: libflint_shared_link_whole,
      install: true,
      version: FLINT_LIBRARY_VERSION,
      soversion: FLINT_MAJOR_SO,
      c_args: libflint_c_args + flint_cpu_c_args + flint_library_warning_c_args,
    )
    libflint = libflint_shared
  endif

  if build_static_library
    libflint_static = static_library('flint', c_files_static + headers_all,
      dependencies: flint_deps,
      include_directories: [headers_built_nodir_inc, 'src'],
      link_whole: libflint_static_link_whole,
      install: true,
      c_args: libflint_c_args + flint_cpu_c_args + flint_library_warning_c_args,
    )

    if not build_shared_library
      libflint = libflint_static
    endif
  endif
endif

flint_test_dep = declare_dependency(
  link_with: libflint,
  dependencies: flint_deps,
  include_directories: [headers_built_nodir_inc, 'src'],
  sources: headers_all,
  compile_args: flint_cpu_c_args,
)

subdir('examples')

foreach mod_profile : mod_profiles
  subdir('src' / mod_profile)
endforeach
alias_target('profile', profile_exes)

foreach mod_tune : mod_tunes
  subdir('src' / mod_tune)
endforeach
alias_target('tune', tune_exes)

foreach mod_test : mod_tests
  subdir('src' / mod_test)
endforeach
alias_target('tests', test_exes)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries: libflint,
  version : meson.project_version(),
  name: 'flint',
  filebase: 'flint',
  description: 'Fast Library for Number Theory'
)
