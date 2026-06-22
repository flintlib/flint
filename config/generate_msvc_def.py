#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import subprocess


DUMPBIN_SYMBOL_RE = re.compile(r"\bExternal\s+\|\s+(\S+)$")
SKIP_PREFIXES = (
    "__imp_",
    "__real@",
    "__xmm@",
    "__ymm@",
    "__local_stdio_",
)
SKIP_SYMBOLS = {
    "main",
}


def object_symbols(obj):
    proc = subprocess.run(
        ["dumpbin", "/symbols", obj],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    for line in proc.stdout.splitlines():
        if " UNDEF " in line:
            continue

        match = DUMPBIN_SYMBOL_RE.search(line)
        if not match:
            continue

        symbol = match.group(1)
        if (
            symbol in SKIP_SYMBOLS
            or symbol.startswith(SKIP_PREFIXES)
            or symbol.startswith("$")
            or symbol.startswith(".")
            or "@" in symbol
            or "?" in symbol
        ):
            continue

        yield symbol


def main():
    parser = argparse.ArgumentParser(
        description="Generate a MSVC module definition file from object files."
    )
    parser.add_argument("output")
    parser.add_argument("objects", nargs="*")
    parser.add_argument("--object-dir")
    args = parser.parse_args()

    objects = list(args.objects)
    if args.object_dir:
        objects.extend(str(obj) for obj in Path(args.object_dir).rglob("*.obj"))

    if not objects:
        raise SystemExit("no object files found")

    symbols = set()
    for obj in objects:
        symbols.update(object_symbols(obj))

    with open(args.output, "w", encoding="utf-8", newline="\n") as output:
        output.write("EXPORTS\n")
        for symbol in sorted(symbols):
            output.write(f"    {symbol}\n")


if __name__ == "__main__":
    main()
