#!/usr/bin/env python3

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, *, cwd=None, input_text=None):
    return subprocess.run(
        cmd,
        cwd=cwd,
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )


def compile_asm(compiler, c_args, asm_text):
    tmpdir = Path(tempfile.mkdtemp(prefix='flint-asm-'))
    asm_file = tmpdir / 'probe.s'
    obj_file = tmpdir / 'probe.o'
    asm_file.write_text(asm_text)
    proc = run(list(compiler) + list(c_args) + ['-c', str(asm_file), '-o', str(obj_file)])
    return proc, obj_file if proc.returncode == 0 else None


def compile_c_to_asm(compiler, c_args, c_text):
    tmpdir = Path(tempfile.mkdtemp(prefix='flint-c-'))
    c_file = tmpdir / 'probe.c'
    s_file = tmpdir / 'probe.s'
    c_file.write_text(c_text)
    proc = run(list(compiler) + list(c_args) + ['-S', str(c_file), '-o', str(s_file)])
    return proc, s_file if proc.returncode == 0 else None


def nm_output(obj_file):
    nm = shutil.which('nm')
    if nm is None:
        return None
    proc = run([nm, str(obj_file)])
    if proc.returncode != 0:
        return None
    return proc.stdout


def probe_asm_candidate(compiler, c_args, asm_text):
    proc, _ = compile_asm(compiler, c_args, asm_text)
    return proc.returncode == 0


def detect_text(compiler, c_args):
    for candidate in ['.text', '.code', '.csect .text[PR]']:
        if probe_asm_candidate(compiler, c_args, f'{candidate}\n'):
            return candidate
    return '.text'


def detect_data(compiler, c_args):
    for candidate in ['.data', '.csect .data[RW]']:
        if probe_asm_candidate(compiler, c_args, f'{candidate}\n'):
            return candidate
    return '.data'


def detect_label_suffix(compiler, c_args, text_dir):
    for candidate in ['', ':']:
        asm = f'''{text_dir}
foo{candidate}
  .byte 0
'''
        if probe_asm_candidate(compiler, c_args, asm):
            return candidate
    return ':'


def detect_gsym_prefix(compiler, c_args, host_system):
    default = '_' if host_system == 'darwin' else ''
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        c_file = tmpdir / 'probe.c'
        obj_file = tmpdir / 'probe.o'
        c_file.write_text('int gurkmacka;\n')
        proc = run(list(compiler) + list(c_args) + ['-c', str(c_file), '-o', str(obj_file)])
        if proc.returncode != 0:
            return default
        out = nm_output(obj_file)
        if out is None:
            return default
        for line in out.splitlines():
            if 'gurkmacka' not in line:
                continue
            symbol = line.split()[-1]
            if symbol == '_gurkmacka':
                return '_'
            if symbol == 'gurkmacka':
                return ''
    return default


def detect_lsym_prefix(compiler, c_args, text_dir, label_suffix, gsym_prefix):
    candidates = ['L', '.L', '$', 'L$']
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for candidate in candidates:
            asm_file = tmpdir / 'probe.s'
            obj_file = tmpdir / 'probe.o'
            asm_file.write_text(f'''{text_dir}
dummy{label_suffix}
{candidate}gurkmacka{label_suffix}
''')
            proc = run(list(compiler) + list(c_args) + ['-c', str(asm_file), '-o', str(obj_file)])
            if proc.returncode != 0:
                continue
            out = nm_output(obj_file)
            if out is None:
                continue
            if 'gurkmacka' not in out:
                return candidate
            for line in out.splitlines():
                if 'gurkmacka' not in line:
                    continue
                fields = line.split()
                if len(fields) >= 2 and fields[1] in ['t', 'T', 'n', 'N']:
                    return candidate
    return 'L' if gsym_prefix == '_' else '.L'


def detect_rodata(compiler, c_args, host_system, gsym_prefix):
    default = '    .section    __TEXT,__const' if host_system == 'darwin' else '.section .rodata'
    c_source = '''extern const int foo[];
const int foo[] = {1, 2, 3};
'''
    proc, s_file = compile_c_to_asm(compiler, c_args, c_source)
    if proc.returncode != 0 or s_file is None:
        return default

    current_section = default
    label_re = re.compile(r'^\s*' + re.escape(gsym_prefix) + r'?foo(:|\b)')
    for line in s_file.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith(('.section', '.csect', '.data', '.rdata', '.text')):
            current_section = stripped
        if label_re.search(line):
            return current_section
    return default


def detect_type_size(compiler, c_args, text_dir, globl, label_suffix, host_system):
    if host_system == 'darwin':
        return '', ''

    type_value = ''
    size_value = ''
    for marker in ['@', '#', '%']:
        asm = f'''{text_dir}
{globl} foo
.type foo,{marker}function
foo{label_suffix}
  .byte 0
'''
        proc, _ = compile_asm(compiler, c_args, asm)
        if proc.returncode == 0:
            type_value = f'.type $1,{marker}$2'
            break

    asm = f'''{text_dir}
{globl} foo
foo{label_suffix}
  .byte 0
.size foo,.-foo
'''
    proc, _ = compile_asm(compiler, c_args, asm)
    if proc.returncode == 0:
        size_value = '.size $1,$2'

    return type_value, size_value


def detect_align_logarithmic(compiler, c_args, data_dir, globl, label_suffix, gsym_prefix):
    asm = f'''{data_dir}
  .byte 1
  .align 4
{globl} foo
foo{label_suffix}
  .byte 2
'''
    proc, obj_file = compile_asm(compiler, c_args, asm)
    if proc.returncode != 0 or obj_file is None:
        return 'yes'

    out = nm_output(obj_file)
    if out is None:
        return 'yes'

    for line in out.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[-1] == f'{gsym_prefix}foo':
            try:
                addr = int(fields[0], 16)
            except ValueError:
                continue
            return 'yes' if addr in (16, 0x10) else 'no'
    return 'yes'


def detect_align_fill_0x90(compiler, c_args, text_dir):
    asm = f'''{text_dir}
  .align 4, 0x90
  .byte 0
  .align 4, 0x90
'''
    proc, _ = compile_asm(compiler, c_args, asm)
    return 'yes' if proc.returncode == 0 else 'no'


def detect_coff_type(compiler, c_args, text_dir, globl, label_suffix, gsym_prefix):
    symbol = f'{gsym_prefix}foo'
    asm = f'''{text_dir}
{globl} {symbol}
.def {symbol}
.scl 2
.type 32
.endef
{symbol}{label_suffix}
  .byte 0
'''
    proc, _ = compile_asm(compiler, c_args, asm)
    return 'yes' if proc.returncode == 0 else 'no'


def main():
    parser = argparse.ArgumentParser(description='Detect assembler syntax support')
    parser.add_argument('--host-cpu', required=True)
    parser.add_argument('--host-system', required=True)
    parser.add_argument('compiler', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    compiler = args.compiler
    if compiler and compiler[0] == '--':
        compiler = compiler[1:]
    if not compiler:
        print('no compiler command provided', file=sys.stderr)
        return 1

    c_args = []
    text_dir = detect_text(compiler, c_args)
    data_dir = detect_data(compiler, c_args)
    label_suffix = detect_label_suffix(compiler, c_args, text_dir)
    globl = '.globl'
    globl_attr = ''
    gsym_prefix = detect_gsym_prefix(compiler, c_args, args.host_system)
    rodata = detect_rodata(compiler, c_args, args.host_system, gsym_prefix)
    type_value, size_value = detect_type_size(
        compiler, c_args, text_dir, globl, label_suffix, args.host_system)
    lsym_prefix = detect_lsym_prefix(compiler, c_args, text_dir, label_suffix, gsym_prefix)
    align_logarithmic = detect_align_logarithmic(compiler, c_args, data_dir, globl, label_suffix, gsym_prefix)
    align_fill_0x90 = 'no'
    have_coff_type = 'no'
    if args.host_cpu == 'x86_64':
        align_fill_0x90 = detect_align_fill_0x90(compiler, c_args, text_dir)
        have_coff_type = detect_coff_type(compiler, c_args, text_dir, globl, label_suffix, gsym_prefix)

    values = {
        'TEXT': text_dir,
        'DATA': data_dir,
        'LABEL_SUFFIX': label_suffix,
        'GLOBL': globl,
        'GLOBL_ATTR': globl_attr,
        'GSYM_PREFIX': gsym_prefix,
        'RODATA': rodata,
        'TYPE': type_value,
        'SIZE': size_value,
        'LSYM_PREFIX': lsym_prefix,
        'ALIGN_LOGARITHMIC': align_logarithmic,
        'ALIGN_FILL_0x90': align_fill_0x90,
        'HAVE_COFF_TYPE': have_coff_type,
    }

    for key, value in values.items():
        print(f'{key}={value}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
