#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)
PROTOTYPE_RE = re.compile(r"\b([A-Za-z_]\w*)\s*\([^;{}]*\)\s*;$", re.S)
SKIP_SYMBOLS = {
    "TEMPLATE",
    "_TEMPLATE",
    "TEMPLATE3",
    "TEMPLATE4",
    "main",
}
GMP_INTERNAL_SYMBOLS = {
    "mpn_add_nc",
    "mpn_addlsh1_n",
    "mpn_addlsh1_n_ip1",
    "mpn_addmul_2",
    "mpn_div_q",
    "mpn_gcd_11",
    "mpn_invert_limb",
    "mpn_modexact_1_odd",
    "mpn_rsh1add_n",
    "mpn_rsh1sub_n",
    "mpn_sub_nc",
}


def meson_headers(source_root, with_fft_small):
    src_root = source_root / "src"
    headers = set(src_root.glob("*.h"))

    for meson_build in src_root.glob("*/meson.build"):
        mod = meson_build.parent.name
        if mod == "fft_small" and not with_fft_small:
            continue

        header = src_root / f"{mod}.h"
        if header.exists():
            headers.add(header)

    return sorted(headers)


def exported_functions(header):
    text = header.read_text(encoding="utf-8", errors="ignore")
    text = COMMENT_RE.sub("", text)

    declaration = []
    in_preprocessor = False
    for line in text.splitlines():
        stripped = line.strip()

        if in_preprocessor or stripped.startswith("#"):
            in_preprocessor = stripped.endswith("\\")
            continue

        if not declaration:
            if (
                not stripped
                or line[:1].isspace()
                or stripped.startswith(("static", "typedef", "struct", "enum"))
                or "FLINT_INLINE" in stripped
                or "FLINT_FORCE_INLINE" in stripped
            ):
                continue
            declaration = [stripped]
        else:
            declaration.append(stripped)

        joined = " ".join(declaration)
        if "{" in joined or "=" in joined or "DECLSPEC_IMPORT" in joined:
            declaration = []
            continue
        if not joined.endswith(";"):
            continue

        if "static" in joined or "typedef" in joined:
            declaration = []
            continue

        match = PROTOTYPE_RE.search(joined)
        if match:
            name = match.group(1)
            if (
                name not in SKIP_SYMBOLS
                and name not in GMP_INTERNAL_SYMBOLS
                and not name.startswith("__gmp")
                and name == name.lower()
            ):
                yield name
        declaration = []


def main():
    parser = argparse.ArgumentParser(
        description="Generate a MSVC module definition file for libflint."
    )
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--with-fft-small", action="store_true")
    args = parser.parse_args()

    symbols = set()
    for header in meson_headers(args.source_root, args.with_fft_small):
        symbols.update(exported_functions(header))

    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        output.write("EXPORTS\n")
        for symbol in sorted(symbols):
            output.write(f"    {symbol}\n")


if __name__ == "__main__":
    main()
