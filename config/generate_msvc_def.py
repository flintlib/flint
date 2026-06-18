#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


C_SOURCE_RE = re.compile(r"['\"]([^'\"]+\.c)['\"]")
COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)
SIGNATURE_RE = re.compile(r"([A-Za-z_]\w*)\s*\([^;{}]*\)\s*$", re.S)
SKIP_SYMBOLS = {
    "TEMPLATE",
    "_TEMPLATE",
    "TEMPLATE3",
    "TEMPLATE4",
    "main",
}


def meson_c_sources(source_root, with_fft_small, fmpz_source):
    src_root = source_root / "src"
    sources = {src_root / "fmpz" / "link" / fmpz_source}

    for meson_build in src_root.glob("*/meson.build"):
        mod = meson_build.parent.name
        if mod == "fft_small" and not with_fft_small:
            continue

        text = meson_build.read_text(encoding="utf-8")
        for match in C_SOURCE_RE.finditer(text):
            source = meson_build.parent / match.group(1)
            if "test" in source.parts or "profile" in source.parts:
                continue
            sources.add(source)

    return sorted(sources)


def exported_functions(source):
    text = source.read_text(encoding="utf-8", errors="ignore")
    text = COMMENT_RE.sub("", text)

    signature = []
    in_preprocessor = False
    for line in text.splitlines():
        stripped = line.strip()

        if in_preprocessor or stripped.startswith("#"):
            in_preprocessor = stripped.endswith("\\")
            continue

        if not signature:
            if (
                not stripped
                or line[:1].isspace()
                or stripped.startswith(("static", "typedef"))
            ):
                continue
            signature = [stripped]
        else:
            signature.append(stripped)

        joined = " ".join(signature)
        before_brace, has_brace, _ = joined.partition("{")
        if ";" in before_brace or "=" in before_brace:
            signature = []
            continue
        if not has_brace:
            continue

        match = SIGNATURE_RE.search(before_brace.strip())
        if match:
            name = match.group(1)
            if name not in SKIP_SYMBOLS and name == name.lower():
                yield name
        signature = []


def main():
    parser = argparse.ArgumentParser(
        description="Generate a MSVC module definition file for libflint."
    )
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fmpz-source", required=True)
    parser.add_argument("--with-fft-small", action="store_true")
    args = parser.parse_args()

    symbols = set()
    for source in meson_c_sources(
        args.source_root, args.with_fft_small, args.fmpz_source
    ):
        if source.exists():
            symbols.update(exported_functions(source))

    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        output.write("EXPORTS\n")
        for symbol in sorted(symbols):
            output.write(f"    {symbol}\n")


if __name__ == "__main__":
    main()
