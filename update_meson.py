#!/usr/bin/env python3

import sys
from os import listdir, makedirs
from os.path import join, dirname, abspath, isdir
from shutil import copy2
from argparse import ArgumentParser

this_dir = dirname(abspath(__file__))

files_to_copy = [
    'meson.build',
    'meson.options',
    'include/meson.build',
    'include/flint/meson.build',
    'config.m4.in',
]

# Directories in src that are not modules or that need special handling
exclude_mod_dirs = [
    'test',
    'profile',
    'interfaces',
    'python',
    'include',
]

conditional_modules = [
    'fft_small',
]

template_modules = [
    'fq_templates',
    'fq_embed_templates',
    'fq_poly_templates',
    'fq_poly_factor_templates',
    'fq_mat_templates',
    'fq_vec_templates',
]

test_only_modules = [
    'fq_zech_vec',
]

non_library_modules = template_modules + test_only_modules

unroll_modules = [
    'ulong_extras',
    'nmod',
    'nmod_vec',
    'nmod_mat',
    'nmod_poly',
    'arith',
]

# Modules that do not have a corresponding header file
mod_no_header = [
    'generic_files',
]

# Headers that do not have a corresponding module
head_no_dir = [
    'NTL-interface',
    'acb_types',
    'acf_types',
    'arb_types',
    'arf_types',
    'ca_types',
    'crt_helpers',
    'fmpq_types',
    'fmpz_mod_types',
    'fmpz_types',
    'fq_nmod_types',
    'fq_types',
    'fq_zech_types',
    'gettimeofday',
    'gr_types',
    'limb_types',
    'longlong',
    'longlong_asm_clang',
    'longlong_asm_gcc',
    'longlong_asm_gnu',
    'longlong_div_gnu',
    'longlong_msc_arm64',
    'longlong_msc_x86',
    'machine_vectors',
    'mpoly_types',
    'n_poly_types',
    'nmod_types',
    'padic_types',
    'profiler',
    'templates',
    'test_helpers',
] + non_library_modules

headers_skip = [
    'flint',
    'config',
    'flint-config',
    'flint-mparam',
    'gmpcompat',
]

mod_no_tests = [
    'ca_vec',
    'calcium',
    'fexpr_builtin',
    'fmpz_mod_vec',
    'fq_zech_mpoly',
    'fq_zech_mpoly_factor',
    'generic_files',
    'hypgeom',
]

mod_no_profiles = [
]

mod_no_tunes = [
]

regression_modules = [
    'thread_pool',
    'thread_support',
    'ulong_extras',
    'long_extras',
    'perm',
    'double_extras',
    'd_vec',
    'd_mat',
    'mpn_extras',
    'nmod',
    'nmod_vec',
    'nmod_mat',
    'nmod_poly',
    'mpn_mod',
    'fmpz',
    'fmpz_vec',
    'fmpz_mat',
    'fmpz_poly',
    'fmpz_mod',
    'fmpz_mod_vec',
    'fmpz_mod_mat',
    'fmpz_mod_poly',
    'fmpq',
    'fmpq_vec',
    'fmpq_mat',
    'fmpq_poly',
    'fq',
    'fq_vec',
    'fq_mat',
    'fq_poly',
    'fq_nmod',
    'fq_nmod_vec',
    'fq_nmod_mat',
    'fq_nmod_poly',
    'fq_zech',
    'fq_zech_mat',
    'fq_zech_poly',
    'fq_default',
    'fq_default_mat',
    'fq_default_poly',
    'fq_embed',
    'fq_nmod_embed',
    'fq_zech_embed',
    'padic',
    'padic_mat',
    'padic_poly',
    'qadic',
    'nmod_poly_factor',
    'fmpz_factor',
    'fmpz_poly_factor',
    'fmpz_mod_poly_factor',
    'fq_poly_factor',
    'fq_nmod_poly_factor',
    'fq_zech_poly_factor',
    'fq_default_poly_factor',
    'nmod_poly_mat',
    'fmpz_poly_mat',
    'mpoly',
    'nmod_mpoly',
    'fmpz_mpoly',
    'fmpz_mod_mpoly',
    'fmpq_mpoly',
    'fq_nmod_mpoly',
    'fq_zech_mpoly',
    'nmod_mpoly_factor',
    'fmpz_mpoly_factor',
    'fmpz_mod_mpoly_factor',
    'fmpq_mpoly_factor',
    'fq_nmod_mpoly_factor',
    'fq_zech_mpoly_factor',
    'fft',
    'fft_small',
]

src_meson_build = '''\
#
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

modules = [
%s
]

modules_no_header = [
%s
]

modules_no_tests = [
%s
]

headers_no_dir = [
%s
]

regression_modules = [
%s
]

unroll_modules = [
%s
]

profile_modules = [
%s
]

conditional_profile_modules = [
%s
]

tune_modules = [
%s
]

headers_all = []
c_files_all = []
c_files_missing_prototypes = []
c_files_missing_prototypes_unroll = []
c_files_no_missing_prototypes = []
c_files_no_missing_prototypes_unroll = []
c_file_buckets_missing_prototypes = []
c_file_buckets_missing_prototypes_unroll = []
c_file_buckets_no_missing_prototypes = []
c_file_buckets_no_missing_prototypes_unroll = []
regression_c_files = []
mod_tests = []
test_exes = []
mod_profiles = []
profile_exes = []
mod_tunes = []
tune_exes = []

# Select the right version of fmpz.c (see configuration)
fmpz_link_c_file = files('fmpz/link' / fmpz_c_in)
c_files_all += fmpz_link_c_file
c_files_missing_prototypes += fmpz_link_c_file
c_file_buckets_missing_prototypes += [['fmpz_link', fmpz_link_c_file]]
regression_c_files += fmpz_link_c_file

foreach mod : modules
  subdir(mod)
  if mod not in modules_no_tests
    mod_tests += [mod / 'test']
  endif
endforeach

foreach mod : modules + headers_no_dir
  if mod not in modules_no_header
    headers_all += [mod + '.h']
  endif
endforeach

# Add conditional modules

if have_fft_small
  subdir('fft_small')
  mod_tests += ['fft_small/test']
  mod_profiles += conditional_profile_modules
  headers_all += ['fft_small.h']
endif

if ntl_opt.enabled()
  mod_tests += ['interfaces/test']
endif

foreach mod : [
%s
]
  mod_tests += [mod / 'test']
endforeach

mod_profiles += profile_modules
mod_tunes += tune_modules

headers_all = files(headers_all)
'''

src_mod_meson_build = '''\
#
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

module_c_files = files(
%s
)

module_c_files_no_missing_prototypes = files(
%s
)

c_files_all += module_c_files
c_files_all += module_c_files_no_missing_prototypes

if '%s' in unroll_modules
  c_files_missing_prototypes_unroll += module_c_files
  c_files_no_missing_prototypes_unroll += module_c_files_no_missing_prototypes
  c_file_buckets_missing_prototypes_unroll += [['%s', module_c_files]]
  c_file_buckets_no_missing_prototypes_unroll += [['%s', module_c_files_no_missing_prototypes]]
else
  c_files_missing_prototypes += module_c_files
  c_files_no_missing_prototypes += module_c_files_no_missing_prototypes
  c_file_buckets_missing_prototypes += [['%s', module_c_files]]
  c_file_buckets_no_missing_prototypes += [['%s', module_c_files_no_missing_prototypes]]
endif

if '%s' in regression_modules
  regression_c_files += module_c_files
  regression_c_files += module_c_files_no_missing_prototypes
endif
'''

test_mod_meson_build = '''\
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

test_exe = executable('main',
  'main.c',
  dependencies: [flint_test_dep],
  install: false,
  build_by_default: false,
)
test_exes += test_exe

test('%s', test_exe, timeout: 0)
'''

profile_mod_meson_build = '''\
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

foreach profile_source : [
%s
]
  profile_stem = fs.stem(profile_source)
  profile_exes += executable(
    profile_stem,
    profile_source,
    dependencies: [flint_test_dep],
    include_directories: [headers_built_inc],
    install: false,
    build_by_default: false,
  )
endforeach
'''

tune_mod_meson_build = '''\
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

foreach tune_source : [
%s
]
  tune_stem = fs.stem(tune_source)
  tune_exes += executable(
    tune_stem,
    tune_source,
    dependencies: [flint_test_dep],
    include_directories: [headers_built_inc],
    install: false,
    build_by_default: false,
  )
endforeach
'''


test_mod_NTL_meson_build = '''\
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

test_exe = executable('main',
  't-NTL-interface.cpp',
  dependencies: [flint_test_dep, ntl_dep],
  install: false,
  build_by_default: false,
)
test_exes += test_exe

test('%s', test_exe, timeout: 0)
'''

examples_meson_build = '''\
#
# This file is generated automatically.
#
# To regenerate it, run:
#   python update_meson.py
#

examples = [
%s
]
example_exes = []

foreach example : examples
  example_exe = executable(
    example,
    example + '.c',
    dependencies: [flint_test_dep],
    include_directories: [headers_built_inc],
    install: false,
    build_by_default: false,
  )

  example_exes += example_exe
endforeach

alias_target('examples', example_exes)
'''

asm_submodule_arm64 = '''\

asm_deps = [
  '../asm-defs.m4',
]

if host_machine.system() == 'darwin'
  asm_deps += ['darwin.m4']
else
  asm_deps += ['arm64-defs.m4']
endif

asm_files = [
%s
]

src_dir_inc = include_directories('.')
'''

asm_submodule_x86_broadwell = '''\

asm_deps = [
    '../../asm-defs.m4',
    '../x86_64-defs.m4',
]

if host_machine.system() == 'darwin'
  asm_deps += ['darwin.m4']
endif

asm_files = [
%s
]
'''

asm_to_s_files = '''\

m4_prog = find_program('m4', native: true)

foreach asm_file: asm_files
  asm_stem = fs.stem(asm_file)

  if default_library_opt in ['static', 'both']
    s_filename = asm_stem + '.s'
    s_file = custom_target(s_filename,
      input: [asm_file, config_m4] + asm_deps,
      output: s_filename,
      command: [m4_prog, '-I', meson.project_build_root(), '@INPUT0@'],
      capture: true,
    )
    s_files_static += [s_file]
  endif

  if default_library_opt in ['shared', 'both']
    s_pic_filename = asm_stem + '_pic.s'
    s_pic_file = custom_target(s_pic_filename,
      input: [asm_file, config_m4] + asm_deps,
      output: s_pic_filename,
      command: [m4_prog, '-I', meson.project_build_root(), '-DPIC', '@INPUT0@'],
      capture: true,
    )
    s_files_shared += [s_pic_file]
  endif
endforeach
'''


asm_modules = [
    ('mpn_extras/arm64', asm_submodule_arm64),
    ('mpn_extras/x86_64/broadwell', asm_submodule_x86_broadwell),
]


def get_flint_modules(flint_root):
    """
    Scan the src directory and return a list of all modules.

    Also check for possible mismatches between subdirs and headers.
    """
    src_path = join(flint_root, 'src')
    is_mod = lambda p: isdir(join(src_path, p)) and p not in exclude_mod_dirs
    subdirs = [p for p in listdir(src_path) if is_mod(p)]
    headers = [p[:-2] for p in listdir(src_path) if p.endswith('.h')]

    # Check for mismatches between subdirs and headers apart from the
    # known exceptions (exclude_mod_dirs, mod_no_header, head_no_dir)
    extra_dirs = set(subdirs) - set(headers) - set(mod_no_header)
    extra_headers = set(headers) - set(subdirs) - set(head_no_dir) - set(headers_skip)
    if extra_dirs or extra_headers:
        print('Mismatch between subdirs and headers in src')
        print('Extra headers:\n', '\n'.join(sorted(extra_headers)))
        print('Extra subdirs:\n', '\n'.join(sorted(extra_dirs)))
        sys.exit(1)
    missing_dirs = set(mod_no_header) - set(subdirs)
    missing_headers = (set(head_no_dir) - set(headers))
    if missing_dirs or missing_headers:
        print('Missing subdirs:\n', '\n'.join(sorted(missing_dirs)))
        print('Missing headers:\n', '\n'.join(sorted(missing_headers)))
        sys.exit(1)

    return subdirs


parser = ArgumentParser(description='Generate Meson build files')
parser.add_argument('-q', '--quiet', action='store_true',
                    help='Do not print anything')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Show all steps')


def main(args):
    args = parser.parse_args(args)
    output_dir = this_dir

    for fname in files_to_copy:
        src_path = join(this_dir, fname)
        dst_path = join(output_dir, fname)
        if not args.quiet:
            print('Copying %s to %s' % (src_path, dst_path))
        copy_file(src_path, dst_path)

    modules = get_flint_modules(output_dir)

    # We will generate the meson.build file for conditional modules but the
    # main meson.build file will decide whether to include them or not.
    modules_unconditional = set(modules) - set(conditional_modules) - set(non_library_modules)

    # src/meson.build
    top_profile_dir = join(output_dir, 'src', 'profile')
    top_profile_sources = []
    if isdir(top_profile_dir):
        top_profile_sources = [
            f for f in listdir(top_profile_dir)
            if f.endswith('.c')
        ]
    profile_modules = []
    if top_profile_sources:
        profile_modules.append('profile')
    conditional_profile_modules = []
    profile_scan_modules = (
        set(modules) - set(non_library_modules) - set(conditional_modules)
    ) | set(mod_no_header)
    for mod in profile_scan_modules:
        profile_dir = join(output_dir, 'src', mod, 'profile')
        if mod not in mod_no_profiles and isdir(profile_dir):
            profile_c_files = [
                f for f in listdir(profile_dir)
                if f.endswith('.c')
            ]
            if profile_c_files:
                profile_modules.append(mod + '/profile')
    for mod in conditional_modules:
        profile_dir = join(output_dir, 'src', mod, 'profile')
        if mod not in mod_no_profiles and isdir(profile_dir):
            profile_c_files = [
                f for f in listdir(profile_dir)
                if f.endswith('.c')
            ]
            if profile_c_files:
                conditional_profile_modules.append(mod + '/profile')

    top_tune_dir = join(output_dir, 'src', 'tune')
    top_tune_sources = []
    if isdir(top_tune_dir):
        top_tune_sources = [
            f for f in listdir(top_tune_dir)
            if f.endswith('.c')
        ]
    tune_modules = []
    if top_tune_sources:
        tune_modules.append('tune')
    tune_scan_modules = (
        set(modules) - set(non_library_modules) - set(conditional_modules)
    ) | set(mod_no_header)
    for mod in tune_scan_modules:
        tune_dir = join(output_dir, 'src', mod, 'tune')
        if mod not in mod_no_tunes and isdir(tune_dir):
            tune_c_files = [
                f for f in listdir(tune_dir)
                if f.endswith('.c')
            ]
            if tune_c_files:
                tune_modules.append(mod + '/tune')

    src_meson_build_text = src_meson_build % (
        format_lines(modules_unconditional),
        format_lines(mod_no_header),
        format_lines(mod_no_tests),
        format_lines(head_no_dir),
        format_lines(regression_modules),
        format_lines(unroll_modules),
        format_lines(profile_modules),
        format_lines(conditional_profile_modules),
        format_lines(tune_modules),
        format_lines(test_only_modules),
    )
    dst_path = join(output_dir, 'src', 'meson.build')
    if not args.quiet:
        print('Writing %s' % dst_path)
    write_file(dst_path, src_meson_build_text)

    if not args.quiet:
        print('Making meson.build files in all modules')
    # src/mod/meson.build
    for mod in modules + mod_no_header:
        mod_dir = join(output_dir, 'src', mod)
        c_files = [f for f in listdir(mod_dir) if f.endswith('.c')]
        if mod == 'fmpz':
            c_files = [f for f in c_files if f != 'fmpz.c']
        c_files_no_missing_prototypes = [f for f in c_files if f == 'inlines.c']
        c_files = [f for f in c_files if f not in c_files_no_missing_prototypes]
        src_mod_meson_build_text = src_mod_meson_build % (
            format_lines(c_files),
            format_lines(c_files_no_missing_prototypes),
            mod,
            mod,
            mod,
            mod,
            mod,
            mod,
        )
        dst_path = join(mod_dir, 'meson.build')
        if args.verbose:
            print('Writing %s' % dst_path)
        write_file(dst_path, src_mod_meson_build_text)

        # src/mod/test/meson.build
        if mod not in mod_no_tests:
            test_dir = join(mod_dir, 'test')
            test_mod_meson_build_text = test_mod_meson_build % mod
            dst_path = join(test_dir, 'meson.build')
            if args.verbose:
                print('Writing %s' % dst_path)
            write_file(dst_path, test_mod_meson_build_text)

        # src/mod/profile/meson.build
        profile_dir = join(mod_dir, 'profile')
        if mod not in mod_no_profiles and isdir(profile_dir):
            profile_c_files = [f for f in listdir(profile_dir) if f.endswith('.c')]
            if profile_c_files:
                profile_mod_meson_build_text = profile_mod_meson_build % (
                    format_lines(profile_c_files),
                )
                dst_path = join(profile_dir, 'meson.build')
                if args.verbose:
                    print('Writing %s' % dst_path)
                write_file(dst_path, profile_mod_meson_build_text)

        # src/mod/tune/meson.build
        tune_dir = join(mod_dir, 'tune')
        if mod not in mod_no_tunes and isdir(tune_dir):
            tune_c_files = [f for f in listdir(tune_dir) if f.endswith('.c')]
            if tune_c_files:
                tune_mod_meson_build_text = tune_mod_meson_build % (
                    format_lines(tune_c_files),
                )
                dst_path = join(tune_dir, 'meson.build')
                if args.verbose:
                    print('Writing %s' % dst_path)
                write_file(dst_path, tune_mod_meson_build_text)

    # Build for the NTL tests
    dst_path = join(output_dir, 'src', 'interfaces', 'test', 'meson.build')
    test_mod_meson_build_text = test_mod_NTL_meson_build % 'NTL-interface'
    if args.verbose:
        print('Writing %s' % dst_path)
    write_file(dst_path, test_mod_meson_build_text)

    # examples/meson.build
    examples_dir = join(output_dir, 'examples')
    examples = [
        f[:-2] for f in listdir(examples_dir)
        if f.endswith('.c')
    ]
    examples_meson_build_text = examples_meson_build % format_lines(examples)
    dst_path = join(examples_dir, 'meson.build')
    if args.verbose:
        print('Writing %s' % dst_path)
    write_file(dst_path, examples_meson_build_text)

    # src/profile/meson.build
    if top_profile_sources:
        profile_mod_meson_build_text = profile_mod_meson_build % (
            format_lines(top_profile_sources),
        )
        dst_path = join(top_profile_dir, 'meson.build')
        if args.verbose:
            print('Writing %s' % dst_path)
        write_file(dst_path, profile_mod_meson_build_text)

    # src/tune/meson.build
    if top_tune_sources:
        tune_mod_meson_build_text = tune_mod_meson_build % (
            format_lines(top_tune_sources),
        )
        dst_path = join(top_tune_dir, 'meson.build')
        if args.verbose:
            print('Writing %s' % dst_path)
        write_file(dst_path, tune_mod_meson_build_text)

    # src/mpn_extras/*/meson.build
    for path, asm_submodule in asm_modules:
        asm_dir = join(output_dir, 'src', path)
        asm_files = [f for f in listdir(asm_dir) if f.endswith('.asm')]
        asm_submodule_text = asm_submodule % format_lines(asm_files)
        asm_submodule_text += asm_to_s_files
        dst_path = join(asm_dir, 'meson.build')
        if args.verbose:
            print('Writing %s' % dst_path)
        write_file(dst_path, asm_submodule_text)


def format_lines(lst):
    return '\n'.join(f"  '{m}'," for m in sorted(lst))


def write_file(dst_path, text):
    makedirs(dirname(dst_path), exist_ok=True)
    with open(dst_path, 'w') as fout:
        fout.write(text)


def copy_file(src_path, dst_path):
    if abspath(src_path) == abspath(dst_path):
        return
    makedirs(dirname(dst_path), exist_ok=True)
    copy2(src_path, dst_path)


def same_files(src_path, dst_path):
    with open(src_path, 'rb') as f:
        src_content = f.read()
    return same_content(src_content, dst_path)


def same_content(src_content, dst_path):
    try:
        with open(dst_path, 'rb') as f:
            dst_content = f.read()
    except FileNotFoundError:
        return False
    else:
        return src_content == dst_content


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
