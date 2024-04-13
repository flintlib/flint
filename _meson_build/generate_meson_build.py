import sys
from os import listdir
from os.path import join, dirname, abspath, isdir
from shutil import copyfile
from argparse import ArgumentParser

this_dir = dirname(abspath(__file__))

files_to_copy = [
    'configure',
    'Makefile',
    'meson.build',
    'meson.options',
]

# Directories in src that are not modules or that need special handling
exclude_mod_dirs = [
    'test',
    'profile',
    'interfaces',
    'python',
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
    'config',
    'crt_helpers',
    'flint',
    'flint-config',
    'flint-mparam',
    'fmpq_types',
    'fmpz_mod_types',
    'fmpz_types',
    'fq_nmod_types',
    'fq_types',
    'fq_zech_types',
    'gettimeofday',
    'gmpcompat',
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

configure_file(
  input: 'flint.h.in',
  output: 'flint.h',
  configuration: cfg_data,
)

src_dir_inc = include_directories('.')

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
  include_directories: src_dir_inc,
  install: false,
)

test('%s', test_exe)
'''

asm_submodule_arm64 = '''\

asm_deps = [
  '../../../config.m4',
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
'''

asm_submodule_x86_broadwell = '''\

asm_deps = [
    '../../../config.m4',
    '../asm-defs.m4',
    'x86-defs.m4',
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
    input: [asm_file] + asm_deps,
    output: s_filename,
    command: [m4_prog, '-I..', '@INPUT0@'],
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
    extra_dirs = set(subdirs) - set(headers) != set(mod_no_header)
    extra_headers = set(headers) - set(subdirs) != set(head_no_dir)
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
parser.add_argument('--error-if-changed', action='store_true',
                    help='Exit with error code 1 if the files have changed')
parser.add_argument('output_dir', default='.', help='Output directory')


def main(args):
    args = parser.parse_args(args)

    for fname in files_to_copy:
        src_path = join(this_dir, fname)
        dst_path = join(args.output_dir, fname)
        copy_file(src_path, dst_path, args)

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
    write_file(dst_path, src_meson_build_text, args)

    # src/mod/meson.build
    for mod in modules + mod_no_header:
        mod_dir = join(args.output_dir, 'src', mod)
        c_files = [f for f in listdir(mod_dir) if f.endswith('.c')]
        src_mod_meson_build_text = src_mod_meson_build % format_lines(c_files)
        dst_path = join(mod_dir, 'meson.build')
        write_file(dst_path, src_mod_meson_build_text, args)

        # src/mod/test/meson.build
        if mod not in mod_no_tests:
            test_dir = join(mod_dir, 'test')
            test_mod_meson_build_text = test_mod_meson_build % mod
            dst_path = join(test_dir, 'meson.build')
            write_file(dst_path, test_mod_meson_build_text, args)

    # src/mpn_extras/*/meson.build
    for path, asm_submodule in asm_modules:
        asm_dir = join(args.output_dir, 'src', path)
        asm_files = [f for f in listdir(asm_dir) if f.endswith('.asm')]
        asm_submodule_text = asm_submodule % format_lines(asm_files)
        asm_submodule_text += asm_to_s_files
        dst_path = join(asm_dir, 'meson.build')
        write_file(dst_path, asm_submodule_text, args)


def format_lines(lst):
    return '\n'.join(f"  '{m}'," for m in sorted(lst))


def write_file(dst_path, text, args):
    if not args.quiet:
        print('Writing %s' % dst_path)
    if args.error_if_changed:
        if not same_content(text, dst_path):
            print('File {} has changed'.format(dst_path))
            sys.exit(1)
    with open(dst_path, 'w') as fout:
        fout.write(text)


def copy_file(src_path, dst_path, args):
    with open(src_path, 'r') as f:
        src_content = f.read()
    write_file(dst_path, src_content, args)


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
