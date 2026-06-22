#!/usr/bin/env python3
import argparse
import json
import os
import shlex
import sys
from collections import Counter, defaultdict
from pathlib import Path


SOURCE_EXTS = ('.c', '.cc', '.cpp', '.cxx', '.S', '.s')


def canonical_generated_source(path):
    parts = Path(path).parts
    for i, part in enumerate(parts):
        if part == 'mpn_extras' and i + 2 < len(parts):
            rel = Path(*parts[i:]).as_posix()
            if rel.endswith('.s'):
                return 'generated/' + rel
    return None


def split_command(command):
    return shlex.split(command, posix=(os.name != 'nt'))


def normalise_source(path, base):
    p = Path(path)
    if not p.is_absolute():
        p = (base / p)
    try:
        source = p.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        source = p.resolve().as_posix()

    generated = canonical_generated_source(source)
    if generated is not None:
        return generated
    return source


def source_from_argv(argv, base):
    for i, arg in enumerate(argv):
        if arg == '-c' and i + 1 < len(argv):
            return normalise_source(argv[i + 1], base)

    for arg in argv:
        if arg.endswith(SOURCE_EXTS):
            return normalise_source(arg, base)

    return None


def source_is_assembly(source):
    return source.endswith('.s') or source.endswith('.S')


def is_compile_command(argv):
    return '-c' in argv and any(arg.endswith(SOURCE_EXTS) for arg in argv)


def normalise_path_arg(arg, base):
    if arg.startswith('-I') and len(arg) > 2:
        prefix, path = '-I', arg[2:]
    else:
        prefix, path = '', arg

    p = Path(path)
    if not p.is_absolute():
        p = base / p
    try:
        path = p.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        path = p.resolve().as_posix()

    return prefix + path


def normalise_relative_path(path):
    p = Path(path)
    if not p.is_absolute():
        p = Path.cwd() / p
    try:
        return p.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return p.resolve().as_posix()


def is_ignored_include(include, build_dir=None):
    if not include.startswith('-I'):
        return False
    if build_dir is None:
        return False

    path = include[2:]
    parts = Path(path).parts
    build_parts = Path(build_dir).parts
    if len(parts) < len(build_parts) + 1:
        return False

    # Meson adds a private include directory for each target's object files
    # (for example build/libflint.24.dylib.p). These are implementation
    # details of the target split and do not correspond to autotools flags.
    after_build = parts[len(build_parts):]
    if parts[:len(build_parts)] != build_parts:
        return False

    if after_build[0].startswith('libflint') and after_build[0].endswith('.p'):
        return True

    # The generated assembly subdirectory is added only to the main library
    # target. It is not relevant to C source comparison, and otherwise splits
    # the inlines.c objects into a separate difference group.
    return len(after_build) >= 3 and after_build[0] == 'src' and after_build[1] == 'mpn_extras'


def classify_args(argv, base, build_dir=None):
    includes = []
    defines = []
    flags = []

    skip_next = False
    compiler_seen = False
    options_with_arg = {'-o', '-MF', '-MQ', '-MT', '-x'}
    ignored_flags = {
        '-c',
        '-MMD',
        '-MD',
        '-MP',
        '-fdiagnostics-color=always',
    }

    for i, arg in enumerate(argv):
        if skip_next:
            skip_next = False
            continue

        if not compiler_seen:
            compiler_seen = True
            continue

        if arg in options_with_arg:
            skip_next = True
            continue

        if i > 0 and argv[i - 1] == '-c':
            continue

        if arg.endswith(SOURCE_EXTS):
            continue

        if arg in ignored_flags:
            continue

        if arg == '-I' and i + 1 < len(argv):
            include = normalise_path_arg(argv[i + 1], base)
            if not is_ignored_include(include, build_dir):
                includes.append(include)
            skip_next = True
            continue

        if arg.startswith('-I'):
            include = normalise_path_arg(arg, base)
            if not is_ignored_include(include, build_dir):
                includes.append(include)
            continue

        if arg == '-D' and i + 1 < len(argv):
            defines.append('-D' + argv[i + 1])
            skip_next = True
            continue

        if arg.startswith('-D'):
            defines.append(arg)
            continue

        flags.append(arg)

    return {
        'includes': includes,
        'defines': defines,
        'flags': flags,
    }


def counter_diff(left, right):
    left_counter = Counter(left)
    right_counter = Counter(right)
    only_left = list((left_counter - right_counter).elements())
    only_right = list((right_counter - left_counter).elements())
    return tuple(sorted(only_left)), tuple(sorted(only_right))


def command_record(command, base, build_dir=None):
    argv = split_command(command)
    if not is_compile_command(argv):
        return None

    source = source_from_argv(argv, base)
    if source is None:
        return None

    return {
        'source': source,
        'argv': argv,
        'parts': classify_args(argv, base, build_dir),
        'raw': command,
    }


def read_autotools_commands(path):
    records = {}
    base = Path.cwd()

    for line in path.read_text(errors='replace').splitlines():
        line = line.strip()
        if not line:
            continue

        try:
            record = command_record(line, base)
        except ValueError:
            continue

        if record is not None:
            records.setdefault(record['source'], record)

    return records


def meson_command_is_default_library(record):
    source = record['source']
    output = record.get('output', '')

    if not (source.startswith('src/') or source.startswith('generated/mpn_extras/')):
        return False
    if '/test/' in source or source.startswith('examples/'):
        return False
    if 'regression_check' in output:
        return False
    if '.exe.p' in output or '/test/' in output or 'examples/' in output:
        return False

    return output.startswith('libflint') or '/libflint' in output


def read_meson_commands(build_dir):
    path = build_dir / 'compile_commands.json'
    if not path.exists():
        raise SystemExit(f'{path} does not exist. Run meson setup first.')

    data = json.loads(path.read_text())
    all_records = {}
    selected = {}
    build_dir_rel = normalise_relative_path(build_dir)

    for entry in data:
        command = entry.get('command')
        if command is None and 'arguments' in entry:
            command = shlex.join(entry['arguments'])
        if command is None:
            continue

        base = Path(entry.get('directory', build_dir))
        record = command_record(command, base, build_dir_rel)
        if record is None:
            continue

        record['output'] = entry.get('output', '')
        all_records.setdefault(record['source'], record)

        if meson_command_is_default_library(record):
            selected.setdefault(record['source'], record)

    return selected, all_records


def print_list(title, items, max_items=None):
    print(title)
    if not items:
        print('  none')
        return
    printed_items = items if max_items is None else items[:max_items]
    for item in printed_items:
        print(f'  {item}')
    if max_items is not None and len(items) > max_items:
        print(f'  ... {len(items) - max_items} more')


def format_tuple(values):
    if not values:
        return 'none'
    return ' '.join(values)


def main():
    parser = argparse.ArgumentParser(
        description='Compare autotools and Meson compile commands by source file.')
    parser.add_argument('build_dir', type=Path)
    parser.add_argument('autotools_commands', type=Path)
    parser.add_argument(
        '--max-sources',
        type=int,
        default=None,
        help='limit the number of sources printed for each difference group')
    args = parser.parse_args()

    autotools = read_autotools_commands(args.autotools_commands)
    meson, meson_all = read_meson_commands(args.build_dir)

    autotools_sources = set(autotools)
    meson_sources = set(meson)
    matched = sorted(autotools_sources & meson_sources)
    only_autotools = sorted(autotools_sources - meson_sources)
    only_meson = sorted(meson_sources - autotools_sources)

    print(f'Autotools compile commands: {len(autotools)}')
    print(f'Meson default-library compile commands: {len(meson)}')
    print(f'Meson compile commands before filtering: {len(meson_all)}')
    print(f'Matched sources: {len(matched)}')
    print()

    print_list('Only in autotools:', only_autotools, args.max_sources)
    print()
    print_list('Only in Meson default library:', only_meson, args.max_sources)
    print()

    groups = defaultdict(list)
    same = 0

    for source in matched:
        if source_is_assembly(source):
            same += 1
            continue

        auto = autotools[source]['parts']
        mes = meson[source]['parts']
        diff = (
            counter_diff(auto['flags'], mes['flags']),
            counter_diff(auto['defines'], mes['defines']),
            counter_diff(auto['includes'], mes['includes']),
        )

        if diff == (((), ()), ((), ()), ((), ())):
            same += 1
        else:
            groups[diff].append(source)

    print(f'Sources with no normalized argument differences: {same}')
    print(f'Difference groups: {len(groups)}')
    print()

    for idx, (diff, sources) in enumerate(sorted(groups.items(), key=lambda x: (-len(x[1]), x[1][0])), 1):
        (flags_auto, flags_meson), (defs_auto, defs_meson), (incs_auto, incs_meson) = diff
        print(f'Group {idx}: {len(sources)} source(s)')
        print(f'  flags only in autotools: {format_tuple(flags_auto)}')
        print(f'  flags only in Meson:     {format_tuple(flags_meson)}')
        print(f'  defines only in autotools: {format_tuple(defs_auto)}')
        print(f'  defines only in Meson:     {format_tuple(defs_meson)}')
        print(f'  includes only in autotools: {format_tuple(incs_auto)}')
        print(f'  includes only in Meson:     {format_tuple(incs_meson)}')
        print('  sources:')
        printed_sources = sources
        if args.max_sources is not None:
            printed_sources = sources[:args.max_sources]
        for source in printed_sources:
            print(f'    {source}')
        if args.max_sources is not None and len(sources) > args.max_sources:
            print(f'    ... {len(sources) - args.max_sources} more')
        print()


if __name__ == '__main__':
    sys.exit(main())
