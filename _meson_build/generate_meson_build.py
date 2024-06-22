import sys
from os import listdir, makedirs
from os.path import join, dirname, abspath, isdir
from shutil import copyfile
from argparse import ArgumentParser

this_dir = dirname(abspath(__file__))

files_to_copy = [
    'configure',
    'Makefile',
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
    'crt_helpers',
    'fmpq_types',
    'fmpz_mod_types',
    'fmpz_types',
    'fq_nmod_types',
    'fq_types',
    'fq_zech_types',
    'gettimeofday',
    'limb_types',
    'longlong',
    'longlong_asm_clang',
    'longlong_asm_gcc',
    'longlong_asm_gnu',
    'longlong_div_gnu',
    'longlong_msc_arm64',
    'longlong_msc_x86',
    'machine_vectors',
    'mpf-impl',
    'mpoly_types',
    'n_poly_types',
    'nmod_types',
    'padic_types',
    'profiler',
    'templates',
    'test_helpers',
]

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
    # XXX: Templates...
    'fq_templates',
    'fq_embed_templates',
    'fq_poly_templates',
    'fq_poly_factor_templates',
    'fq_mat_templates',
    'fq_vec_templates',
]

src_meson_build = '''\
#
# This file is generated automatically.
#
# To regenerate it, run:
#   python _meson_build/generate_meson_build.py
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

headers_all = []
c_files_all = []
mod_tests = []

# Select the right version of fmpz.c (see configuration)
c_files_all += files('fmpz/link' / fmpz_c_in)

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
  headers_all += ['fft_small.h']
endif

if ntl_opt.enabled()
  mod_tests += ['interfaces/test']
endif

headers_all = files(headers_all)
'''

src_mod_meson_build = '''\
#
# This file is generated automatically.
#
# To regenerate it, run:
#   python _meson_build/generate_meson_build.py
#

c_files_all += files(
%s
)
'''

test_mod_meson_build = '''\
# This file is generated automatically.
#
# To regenerate it, run:
#   python _meson_build/generate_meson_build.py
#

test_exe = executable('main',
  'main.c',
  dependencies: flint_deps,
  link_with: libflint,
  include_directories: [headers_built_nodir_inc, '../..'],
  install: false,
)

test('%s', test_exe)
'''


test_mod_NTL_meson_build = '''\
# This file is generated automatically.
#
# To regenerate it, run:
#   python _meson_build/generate_meson_build.py
#

# NTL no found by pkgconfig
cpp = meson.get_compiler('cpp')
ntl_dep = cpp.find_library('ntl')

test_exe = executable('main',
  't-NTL-interface.cpp',
  dependencies: flint_deps + [ntl_dep],
  link_with: libflint,
  include_directories: [headers_built_nodir_inc, '../..'],
  install: false,
)

test('%s', test_exe)
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
  s_filename = fs.stem(asm_file) + '.s'
  s_file = custom_target(s_filename,
    input: [asm_file, config_m4] + asm_deps,
    output: s_filename,
    command: [m4_prog, '@INPUT0@'],
    capture: true,
  )
  s_files += [s_file]
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
parser.add_argument('output_dir', default='.', help='Output directory')


def main(args):
    args = parser.parse_args(args)

    for fname in files_to_copy:
        src_path = join(this_dir, fname)
        dst_path = join(args.output_dir, fname)
        if not args.quiet:
            print('Copying %s to %s' % (src_path, dst_path))
        copy_file(src_path, dst_path)

    modules = get_flint_modules(args.output_dir)

    # We will generate the meson.build file for conditional modules but the
    # main meson.build file will decide whether to include them or not.
    modules_unconditional = set(modules) - set(conditional_modules)

    # src/meson.build
    src_meson_build_text = src_meson_build % (
        format_lines(modules_unconditional),
        format_lines(mod_no_header),
        format_lines(mod_no_tests),
        format_lines(head_no_dir),
    )
    dst_path = join(args.output_dir, 'src', 'meson.build')
    if not args.quiet:
        print('Writing %s' % dst_path)
    write_file(dst_path, src_meson_build_text)

    if not args.quiet:
        print('Making meson.build files in all modules')
    # src/mod/meson.build
    for mod in modules + mod_no_header:
        mod_dir = join(args.output_dir, 'src', mod)
        c_files = [f for f in listdir(mod_dir) if f.endswith('.c')]
        if mod == 'fmpz':
            c_files = [f for f in c_files if f != 'fmpz.c']
        src_mod_meson_build_text = src_mod_meson_build % format_lines(c_files)
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

    # Build for the NTL tests
    dst_path = join(args.output_dir, 'src', 'interfaces', 'test', 'meson.build')
    test_mod_meson_build_text = test_mod_NTL_meson_build % 'NTL-interface'
    if args.verbose:
        print('Writing %s' % dst_path)
    write_file(dst_path, test_mod_meson_build_text)

    # src/mpn_extras/*/meson.build
    for path, asm_submodule in asm_modules:
        asm_dir = join(args.output_dir, 'src', path)
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
    makedirs(dirname(dst_path), exist_ok=True)
    copyfile(src_path, dst_path)


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
