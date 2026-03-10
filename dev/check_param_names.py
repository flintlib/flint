#!/usr/bin/env python3
"""
Compare parameter names between FLINT header files, RST documentation,
and .c source files.

Reports mismatches where the same function has different parameter names
across these sources.

Usage:
    python3 dev/check_param_names.py [--module MODULE] [--check-src]
    python3 dev/check_param_names.py --check-src -m fmpz

Run from the FLINT root directory.
"""

import os
import re
import sys
import argparse
from collections import OrderedDict


# Decorators/macros that appear before return types in headers
DECL_PREFIXES = re.compile(
    r"^(?:"
    r"WARN_UNUSED_RESULT|"
    r"FLINT_FORCE_INLINE|"
    r"\w+_INLINE|"          # FMPZ_INLINE, GR_POLY_INLINE, etc.
    r"FLINT_DLL"
    r")\s*",
)


def find_modules(src_dir, doc_dir):
    """Find modules that have both a header and a doc file."""
    modules = []
    for f in sorted(os.listdir(src_dir)):
        if f.endswith(".h"):
            mod = f[:-2]  # strip .h
            rst = os.path.join(doc_dir, mod + ".rst")
            header = os.path.join(src_dir, f)
            if os.path.isfile(rst):
                modules.append((mod, header, rst))
    return modules


def extract_param_name(param_str):
    """
    Extract the parameter name from a C parameter declaration string.

    Examples:
        "const fmpz_t x"  -> "x"
        "slong n"          -> "n"
        "ulong * out"      -> "out"
        "FILE * file"      -> "file"
        "nn_srcptr xp"     -> "xp"
        "..."              -> None
        "void"             -> None
        "const padic_ctx_t FLINT_UNUSED(ctx)" -> "ctx"
    """
    param_str = param_str.strip()
    if not param_str or param_str == "void" or param_str == "...":
        return None

    # Handle FLINT_UNUSED(name)
    m = re.search(r"FLINT_UNUSED\((\w+)\)", param_str)
    if m:
        return m.group(1)

    # Remove array brackets: e.g. "ulong out[3]" -> name is "out"
    param_str = re.sub(r"\[.*?\]", "", param_str).strip()

    # Split on spaces and asterisks, take the last word token
    tokens = re.findall(r"\w+", param_str)
    if not tokens:
        return None

    # The last token is the name, unless the entire thing is just type(s)
    # with no name (e.g., "void" - already handled above)
    name = tokens[-1]

    # Sanity check: if the "name" looks like a type keyword, skip
    type_keywords = {
        "void", "int", "long", "short", "char", "unsigned", "signed",
        "float", "double", "size_t", "ssize_t", "FILE",
    }
    if name in type_keywords:
        return None

    # If the "name" looks like a type (ends with _t, _struct, _ptr,
    # _srcptr), it's an unnamed parameter (type only).
    type_suffixes = ("_t", "_struct", "_ptr", "_srcptr")
    if name.endswith(type_suffixes):
        return None

    return name


def split_params(params_str):
    """
    Split a parameter list string by commas, respecting nested parentheses
    and function pointer syntax.
    """
    params = []
    depth = 0
    current = []
    for ch in params_str:
        if ch in ("(", "["):
            depth += 1
            current.append(ch)
        elif ch in (")", "]"):
            depth -= 1
            current.append(ch)
        elif ch == "," and depth == 0:
            params.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
    last = "".join(current).strip()
    if last:
        params.append(last)
    return params


def parse_func_signature(sig):
    """
    Parse a C function signature string, returning (func_name, [param_names])
    or None if it can't be parsed.

    The signature should be like:
        "void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)"
    """
    # Normalize whitespace
    sig = re.sub(r"\s+", " ", sig).strip()

    # Remove trailing semicolons
    sig = sig.rstrip(";").strip()

    # Must have parentheses
    if "(" not in sig or ")" not in sig:
        return None

    # Find the first '(' at depth 0 (the parameter list opening)
    depth = 0
    paren_start = -1
    for i, ch in enumerate(sig):
        if ch == "(":
            if depth == 0:
                paren_start = i
                break
            depth += 1
        elif ch == ")":
            depth -= 1

    if paren_start < 0:
        return None

    # Find matching close paren
    depth = 1
    paren_end = -1
    for i in range(paren_start + 1, len(sig)):
        if sig[i] == "(":
            depth += 1
        elif sig[i] == ")":
            depth -= 1
            if depth == 0:
                paren_end = i
                break

    if paren_end < 0:
        return None

    before_paren = sig[:paren_start].strip()
    params_str = sig[paren_start + 1 : paren_end].strip()

    # Skip function pointer declarations: "typedef void (*name)(...)"
    if "(*" in before_paren:
        return None

    # Extract function name: last word before the paren
    name_match = re.search(r"(\w+)\s*$", before_paren)
    if not name_match:
        return None
    func_name = name_match.group(1)

    # Parse parameter names
    if not params_str or params_str == "void":
        return (func_name, [])

    param_parts = split_params(params_str)
    param_names = []
    for p in param_parts:
        # Function pointer parameters: "int (*cmp)(void *, const void *)"
        fptr_match = re.match(r".*\(\s*\*\s*(\w+)\s*\)", p)
        if fptr_match:
            param_names.append(fptr_match.group(1))
            continue

        name = extract_param_name(p)
        if name is not None:
            param_names.append(name)

    return (func_name, param_names)


def strip_c_comments_and_preprocessor(content):
    """
    Remove block comments, line comments, and preprocessor directive lines
    from C source code. Also remove extern "C" { wrappers.
    Returns cleaned text.
    """
    # Remove block comments
    content = re.sub(r"/\*.*?\*/", " ", content, flags=re.DOTALL)

    # Remove line comments (shouldn't appear in FLINT style, but just in case)
    content = re.sub(r"//[^\n]*", "", content)

    # Remove preprocessor directive lines (keep code inside #ifdef blocks)
    lines = content.split("\n")
    clean_lines = []
    in_macro_def = False
    for line in lines:
        stripped = line.strip()
        if in_macro_def:
            if not stripped.endswith("\\"):
                in_macro_def = False
            clean_lines.append("")
            continue
        if stripped.startswith("#"):
            if stripped.endswith("\\"):
                in_macro_def = True
            clean_lines.append("")
            continue
        clean_lines.append(line)

    text = "\n".join(clean_lines)

    # Remove extern "C" { wrappers
    text = re.sub(r'extern\s+"C"\s*\{', " ", text)

    return text


def extract_declarations(text):
    """
    Extract top-level declarations (ending with ';') from cleaned C text.
    Skips anything inside brace-delimited blocks (function bodies, structs).
    Returns list of (declaration_text, line_number).
    """
    declarations = []
    current_decl = []
    brace_depth = 0

    for i, ch in enumerate(text):
        if ch == "{":
            brace_depth += 1
            current_decl = []
        elif ch == "}":
            brace_depth -= 1
            if brace_depth < 0:
                brace_depth = 0
            current_decl = []
        elif brace_depth == 0:
            if ch == ";":
                decl_text = "".join(current_decl).strip()
                if decl_text and "(" in decl_text and ")" in decl_text:
                    pos_in_text = i - len(decl_text)
                    line_num = text[:pos_in_text].count("\n") + 1
                    declarations.append((decl_text, line_num))
                current_decl = []
            else:
                if not current_decl and ch == "\n":
                    pass  # skip leading newlines
                else:
                    current_decl.append(ch)

    return declarations


def extract_definitions(text):
    """
    Extract top-level function definitions from cleaned C text.
    A definition is a function signature followed by a '{' body.
    Returns list of (signature_text, line_number).
    """
    definitions = []
    brace_depth = 0
    current_sig = []
    sig_line = 1
    i = 0

    while i < len(text):
        ch = text[i]

        if ch == "{":
            if brace_depth == 0 and current_sig:
                sig_text = "".join(current_sig).strip()
                if sig_text and "(" in sig_text and ")" in sig_text:
                    definitions.append((sig_text, sig_line))
                current_sig = []
            brace_depth += 1
        elif ch == "}":
            brace_depth -= 1
            if brace_depth < 0:
                brace_depth = 0
            if brace_depth == 0:
                current_sig = []
        elif brace_depth == 0:
            if ch == ";":
                # This is a declaration, not a definition — skip it
                current_sig = []
            else:
                if not current_sig and ch == "\n":
                    pass
                else:
                    if not current_sig:
                        sig_line = text[:i].count("\n") + 1
                    current_sig.append(ch)
        i += 1

    return definitions


def parse_header_declarations(header_path):
    """
    Parse non-inline function declarations from a header file.
    Returns dict: func_name -> (param_names, line_number, raw_declaration)
    """
    with open(header_path, "r") as f:
        content = f.read()

    text = strip_c_comments_and_preprocessor(content)
    declarations = extract_declarations(text)

    functions = {}
    for decl_text, line_num in declarations:
        # Skip typedefs
        if re.match(r"\s*typedef\b", decl_text):
            continue
        # Skip extern declarations without function signatures
        if re.match(r"\s*extern\b", decl_text) and "(" not in decl_text:
            continue

        # Strip decorator prefixes
        clean = decl_text.strip()
        while True:
            m = DECL_PREFIXES.match(clean)
            if m:
                clean = clean[m.end():]
            else:
                break

        result = parse_func_signature(clean)
        if result:
            func_name, param_names = result
            # Skip ALL_CAPS names without underscores (macros)
            if func_name.isupper() and "_" not in func_name:
                continue
            if func_name not in functions:
                functions[func_name] = (param_names, line_num, decl_text.strip())

    return functions


def parse_source_definitions(src_dir, module):
    """
    Parse function definitions from .c source files for a module.
    Returns dict: func_name -> (param_names, file_path, line_number)
    """
    module_dir = os.path.join(src_dir, module)
    if not os.path.isdir(module_dir):
        return {}

    functions = {}

    for fname in sorted(os.listdir(module_dir)):
        if not fname.endswith(".c"):
            continue
        # Skip test files
        if fname.startswith("t-") or fname == "main.c":
            continue

        fpath = os.path.join(module_dir, fname)
        with open(fpath, "r") as f:
            content = f.read()

        text = strip_c_comments_and_preprocessor(content)
        definitions = extract_definitions(text)

        for sig_text, line_num in definitions:
            # Strip 'static' qualifier — we skip static functions
            clean = sig_text.strip()
            if re.match(r"\s*static\b", clean):
                continue

            # Strip decorator prefixes
            while True:
                m = DECL_PREFIXES.match(clean)
                if m:
                    clean = clean[m.end():]
                else:
                    break

            result = parse_func_signature(clean)
            if result:
                func_name, param_names = result
                if func_name not in functions:
                    functions[func_name] = (param_names, fpath, line_num)

    return functions


def parse_rst_functions(rst_path):
    """
    Parse function signatures from RST documentation.
    Returns dict: func_name -> (param_names, line_number, raw_signature)
    """
    functions = {}

    func_directive = re.compile(r"^\.\.\s+(?:c:)?function\s*::\s*(.+)$")
    continuation = re.compile(r"^\s{5,}(\S.+)$")

    with open(rst_path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        m = func_directive.match(line)
        if m:
            sigs = [(m.group(1).strip(), i + 1)]

            j = i + 1
            while j < len(lines):
                cm = continuation.match(lines[j].rstrip())
                if cm:
                    text = cm.group(1).strip()
                    if "(" in text and ")" in text:
                        sigs.append((text, j + 1))
                        j += 1
                    else:
                        break
                else:
                    break

            for sig_text, line_num in sigs:
                result = parse_func_signature(sig_text)
                if result:
                    func_name, param_names = result
                    if func_name not in functions:
                        functions[func_name] = (
                            param_names, line_num, sig_text.strip()
                        )

            i = j
        else:
            i += 1

    return functions


def compare_params(funcs_a, funcs_b):
    """
    Compare parameter names between two function dictionaries.
    Returns list of mismatch dicts.
    """
    mismatches = []
    common = set(funcs_a.keys()) & set(funcs_b.keys())

    for func_name in sorted(common):
        a_params = funcs_a[func_name][0]
        b_params = funcs_b[func_name][0]

        if len(a_params) != len(b_params):
            mismatches.append({
                "func": func_name,
                "type": "param_count",
                "a_params": a_params,
                "b_params": b_params,
                "a_info": funcs_a[func_name],
                "b_info": funcs_b[func_name],
            })
            continue

        if a_params != b_params:
            diffs = []
            for k, (a, b) in enumerate(zip(a_params, b_params)):
                if a != b:
                    diffs.append((k, a, b))
            mismatches.append({
                "func": func_name,
                "type": "param_name",
                "diffs": diffs,
                "a_params": a_params,
                "b_params": b_params,
                "a_info": funcs_a[func_name],
                "b_info": funcs_b[func_name],
            })

    return mismatches


def format_location(info, label):
    """Format a location string from function info tuple."""
    if len(info) == 3:
        # (params, line_or_path, line_or_raw)
        # Header: (params, line_num, raw_decl)
        # Doc: (params, line_num, raw_sig)
        # Source: (params, file_path, line_num)
        if isinstance(info[1], str) and os.sep in info[1]:
            return f"{label} ({info[1]}:{info[2]})"
        elif isinstance(info[1], int):
            return f"{label} (line {info[1]})"
    return label


def print_mismatches(mismatches, label_a, label_b, path_a="", path_b=""):
    """Print mismatch results."""
    for mm in mismatches:
        a_info = mm["a_info"]
        b_info = mm["b_info"]

        if mm["type"] == "param_count":
            print(f"\n  {mm['func']}: parameter count mismatch")
            loc_a = format_location(a_info, label_a)
            loc_b = format_location(b_info, label_b)
            if path_a:
                print(f"    {path_a}:{a_info[1]}: {mm['a_params']}")
            else:
                print(f"    {loc_a}: {mm['a_params']}")
            if path_b:
                print(f"    {path_b}:{b_info[1]}: {mm['b_params']}")
            else:
                print(f"    {loc_b}: {mm['b_params']}")
        else:
            print(f"\n  {mm['func']}:")
            for idx, a_name, b_name in mm["diffs"]:
                print(f"    param {idx}: {label_a} has '{a_name}'"
                      f", {label_b} has '{b_name}'")
            if path_a:
                print(f"    {label_a} ({path_a}:{a_info[1]})")
            else:
                loc = format_location(a_info, label_a)
                print(f"    {loc}")
            if path_b:
                print(f"    {label_b} ({path_b}:{b_info[1]})")
            else:
                # Source info has (params, filepath, line_num)
                if isinstance(b_info[1], str):
                    print(f"    {label_b} ({b_info[1]}:{b_info[2]})")
                else:
                    loc = format_location(b_info, label_b)
                    print(f"    {loc}")


def find_function_span(content, func_name):
    """
    Find the byte range of a function definition (signature + body) in
    C source content. Returns (sig_start, body_end) or None.
    sig_start is the index of the start of the signature.
    body_end is the index just past the closing '}'.
    """
    # Search for the function name followed by '('
    pattern = re.compile(r"\b" + re.escape(func_name) + r"\s*\(")
    for m in pattern.finditer(content):
        # Walk backward to find the start of the signature
        # (return type, possibly on previous line)
        sig_start = m.start()
        # Walk backward past whitespace and return type tokens
        j = sig_start - 1
        while j >= 0 and content[j] in " \t":
            j -= 1
        # Walk backward past the return type line
        while j >= 0 and content[j] != "\n" and content[j] != ";"\
                and content[j] != "}":
            j -= 1
        if j >= 0 and content[j] == "\n":
            # Check if previous line is the return type
            line_start = j + 1
            prev_text = content[line_start:sig_start].strip()
            if prev_text and not prev_text.startswith("/*") \
                    and not prev_text.startswith("//") \
                    and not prev_text.startswith("#"):
                sig_start = line_start

        # Walk forward from the match to find the opening '{'
        pos = m.end()
        paren_depth = 1
        while pos < len(content) and paren_depth > 0:
            if content[pos] == "(":
                paren_depth += 1
            elif content[pos] == ")":
                paren_depth -= 1
            pos += 1

        # Now pos is right after the closing ')' of the parameter list.
        # Skip whitespace to find '{' or ';'
        while pos < len(content) and content[pos] in " \t\n\r":
            pos += 1

        if pos >= len(content):
            continue

        if content[pos] == ";":
            # This is a declaration, not a definition — skip
            continue

        if content[pos] != "{":
            continue

        # Found opening brace — find matching close brace
        brace_depth = 1
        pos += 1
        while pos < len(content) and brace_depth > 0:
            if content[pos] == "{":
                brace_depth += 1
            elif content[pos] == "}":
                brace_depth -= 1
            pos += 1

        return (sig_start, pos)

    return None


def rename_param_in_range(content, start, end, old_name, new_name):
    """
    Rename a parameter within a byte range of content, avoiding
    struct member accesses (->old_name, .old_name).
    Returns modified content.
    """
    region = content[start:end]

    # Replace old_name as a word, but not when preceded by -> or .
    # (?<![>.]) excludes the second char of -> and the . of member access
    pattern = r"(?<![>.])\b" + re.escape(old_name) + r"\b"
    new_region = re.sub(pattern, new_name, region)

    return content[:start] + new_region + content[end:]


def fix_source_file(filepath, func_name, renames, dry_run=False):
    """
    Rename parameters in a function definition in a .c source file.
    renames: list of (old_name, new_name)
    Returns True if changes were made.
    """
    with open(filepath, "r") as f:
        content = f.read()

    span = find_function_span(content, func_name)
    if span is None:
        print(f"  WARNING: could not find {func_name} in {filepath}")
        return False

    start, end = span
    region = content[start:end]

    # Check for macros above the function that capture variable names.
    # If a #define before the function references old_name, renaming
    # will break macro expansions.
    preamble = content[:start]
    safe_renames = []
    for old_name, new_name in renames:
        # Check collision: new_name already exists in the function span
        new_pat = r"(?<![>.])\b" + re.escape(new_name) + r"\b"
        if re.search(new_pat, region):
            print(f"  SKIP {func_name} rename {old_name}->{new_name}"
                  f" in {filepath}: collision with existing '{new_name}'")
            continue
        # Check macro capture: old_name in a #define above.
        # Handle multi-line macros (backslash continuation).
        macro_pat = (r"#\s*define\s+\w+(?:[^\n]*\\\n)*[^\n]*\b"
                     + re.escape(old_name) + r"\b")
        if re.search(macro_pat, preamble):
            print(f"  SKIP {func_name} rename {old_name}->{new_name}"
                  f" in {filepath}: macro captures '{old_name}'")
            continue
        safe_renames.append((old_name, new_name))

    if not safe_renames:
        return False

    changed = False
    for old_name, new_name in safe_renames:
        new_content = rename_param_in_range(
            content, start, end, old_name, new_name
        )
        if new_content != content:
            end += len(new_content) - len(content)
            content = new_content
            changed = True

    if changed:
        if dry_run:
            print(f"  [dry-run] would fix {func_name} in {filepath}")
        else:
            with open(filepath, "w") as f:
                f.write(content)
            print(f"  fixed {func_name} in {filepath}")
    return changed


def fix_declaration_file(filepath, func_name, renames, dry_run=False):
    """
    Rename parameters in a function declaration in a .h header file
    or a .. function:: directive in an .rst doc file.
    Only modifies the signature line(s), not any body.
    Skips inline function definitions (only fixes pure declarations
    ending with ';').
    renames: list of (old_name, new_name)
    Returns True if changes were made.
    """
    with open(filepath, "r") as f:
        content = f.read()

    # Find the function name in the file
    pattern = re.compile(r"\b" + re.escape(func_name) + r"\s*\(")
    changed = False

    for m in pattern.finditer(content):
        # Find the full declaration: from func_name( to the closing )
        pos = m.end()
        paren_depth = 1
        while pos < len(content) and paren_depth > 0:
            if content[pos] == "(":
                paren_depth += 1
            elif content[pos] == ")":
                paren_depth -= 1
            pos += 1

        # For .h files, check if this is a declaration (;) or
        # definition ({). Skip inline definitions — they need
        # fix_source_file instead.
        if filepath.endswith(".h"):
            rest = content[pos:pos + 40].lstrip()
            if not rest.startswith(";"):
                continue  # Skip inline definitions

        param_start = m.end() - 1  # the '('
        param_end = pos              # just past ')'

        region = content[param_start:param_end]
        new_region = region
        for old_name, new_name in renames:
            pat = r"\b" + re.escape(old_name) + r"\b"
            new_region = re.sub(pat, new_name, new_region)

        if new_region != region:
            content = content[:param_start] + new_region + content[param_end:]
            changed = True
            break  # Only fix first occurrence

    if changed:
        if dry_run:
            print(f"  [dry-run] would fix {func_name} in {filepath}")
        else:
            with open(filepath, "w") as f:
                f.write(content)
            print(f"  fixed {func_name} in {filepath}")
    return changed


def collect_mismatches(modules, src_dir, check_src=False):
    """
    Collect all mismatches across modules.
    Returns list of (mod, header_path, rst_path, hdr_doc_mm, hdr_src_mm,
                      hdr_funcs, doc_funcs, src_funcs).
    """
    results = []
    for mod, header_path, rst_path in modules:
        hdr_funcs = parse_header_declarations(header_path)
        doc_funcs = parse_rst_functions(rst_path)
        src_funcs = {}

        hdr_doc_mm = compare_params(hdr_funcs, doc_funcs)

        hdr_src_mm = []
        if check_src:
            src_funcs = parse_source_definitions(src_dir, mod)
            hdr_src_mm = compare_params(hdr_funcs, src_funcs)

        results.append((mod, header_path, rst_path,
                         hdr_doc_mm, hdr_src_mm,
                         hdr_funcs, doc_funcs, src_funcs))
    return results


def apply_fixes(results, src_dir, dry_run=False):
    """
    Auto-fix mismatches. Strategy:
    - If doc and source agree but header differs: fix header
    - If doc and header agree but source differs: fix source
    - If header and source agree but doc differs: fix header+source
      (doc is source of truth for naming)
    - If all three differ: skip (needs manual review)
    """
    fixed = 0
    skipped = 0

    for (mod, header_path, rst_path,
         hdr_doc_mm, hdr_src_mm,
         hdr_funcs, doc_funcs, src_funcs) in results:

        if not hdr_doc_mm and not hdr_src_mm:
            continue

        # Build a unified view per function
        all_funcs = set()
        for mm in hdr_doc_mm:
            all_funcs.add(mm["func"])
        for mm in hdr_src_mm:
            all_funcs.add(mm["func"])

        hdr_doc_by_func = {mm["func"]: mm for mm in hdr_doc_mm}
        hdr_src_by_func = {mm["func"]: mm for mm in hdr_src_mm}

        for func_name in sorted(all_funcs):
            hd = hdr_doc_by_func.get(func_name)
            hs = hdr_src_by_func.get(func_name)

            hdr_params = hdr_funcs[func_name][0] if func_name in hdr_funcs else None
            doc_params = doc_funcs[func_name][0] if func_name in doc_funcs else None
            src_params = src_funcs[func_name][0] if func_name in src_funcs else None

            # Skip param count mismatches — need manual review
            if (hd and hd["type"] == "param_count") or \
               (hs and hs["type"] == "param_count"):
                print(f"  SKIP {func_name}: parameter count mismatch"
                      " (manual review needed)")
                skipped += 1
                continue

            # Determine which names to use
            renames_hdr = []   # (old, new) for header
            renames_src = []   # (old, new) for source file(s)

            if hd and hs:
                # Header disagrees with both doc and source
                # Check if doc and source agree
                if doc_params == src_params:
                    # doc and source agree — fix header to match
                    for idx, old, new in hd["diffs"]:
                        renames_hdr.append((old, new))
                else:
                    # All three differ — use doc as truth
                    for idx, old, new in hd["diffs"]:
                        renames_hdr.append((old, new))
                    # Also fix source to match doc
                    # Recompute source renames against doc
                    if src_params and doc_params and \
                            len(src_params) == len(doc_params):
                        for k, (s, d) in enumerate(
                                zip(src_params, doc_params)):
                            if s != d:
                                renames_src.append((s, d))
                    else:
                        print(f"  SKIP {func_name}: all three differ"
                              " (manual review needed)")
                        skipped += 1
                        continue
            elif hd and not hs:
                # Header vs doc mismatch only (no source mismatch or
                # source not checked)
                if src_params is not None:
                    # Source exists — check what it says
                    if src_params == doc_params:
                        # Source agrees with doc — fix header
                        for idx, old, new in hd["diffs"]:
                            renames_hdr.append((old, new))
                    elif src_params == hdr_params:
                        # Source agrees with header — doc is truth,
                        # fix header + source
                        for idx, old, new in hd["diffs"]:
                            renames_hdr.append((old, new))
                            renames_src.append((old, new))
                    else:
                        print(f"  SKIP {func_name}: all three differ"
                              " (manual review needed)")
                        skipped += 1
                        continue
                else:
                    # No source — just fix header to match doc
                    for idx, old, new in hd["diffs"]:
                        renames_hdr.append((old, new))
            elif hs and not hd:
                # Header vs source mismatch only
                # Header agrees with doc (or doc doesn't exist)
                # Fix source to match header
                # diffs are (idx, hdr_name, src_name) — rename
                # src_name -> hdr_name in the source file
                for idx, hdr_name, src_name in hs["diffs"]:
                    renames_src.append((src_name, hdr_name))

            # Filter out renames that would collide.
            # Case 1: new_name is already a param not being renamed.
            # Case 2: overlapping renames (swaps) where new_name
            #   equals old_name of another rename — sequential
            #   regex can't handle this correctly.
            def filter_safe_renames(renames, params, label):
                if not renames or not params:
                    return renames
                old_names = {old for old, new in renames}
                new_names = {new for old, new in renames}
                # Check for overlapping renames (any new is also an old)
                if old_names & new_names:
                    for old, new in renames:
                        print(f"  SKIP {func_name} rename {old}->{new}"
                              f" in {label}: overlapping renames")
                    return []
                # Check if new_name collides with an unchanged param
                unchanged = set(params) - old_names
                safe = []
                for old, new in renames:
                    if new in unchanged:
                        print(f"  SKIP {func_name} rename {old}->{new}"
                              f" in {label}: '{new}' already a param")
                    else:
                        safe.append((old, new))
                return safe

            renames_hdr = filter_safe_renames(
                renames_hdr, hdr_params, "header"
            )

            # Apply fixes
            if renames_hdr:
                ok = fix_declaration_file(
                    header_path, func_name, renames_hdr, dry_run
                )
                if not ok:
                    # Declaration not found (inline definition only).
                    # Try fix_source_file on the header instead.
                    ok = fix_source_file(
                        header_path, func_name, renames_hdr, dry_run
                    )
                if ok:
                    fixed += 1

            if renames_src:
                # Find the source file
                if func_name in src_funcs:
                    src_path = src_funcs[func_name][1]
                    ok = fix_source_file(
                        src_path, func_name, renames_src, dry_run
                    )
                    if ok:
                        fixed += 1

    return fixed, skipped


def main():
    parser = argparse.ArgumentParser(
        description="Check parameter name consistency between "
                    "FLINT headers, docs, and source files"
    )
    parser.add_argument(
        "--module", "-m",
        help="Only check a specific module (e.g., fmpz, acb_poly)",
    )
    parser.add_argument(
        "--check-src", action="store_true",
        help="Also check .c source files against headers",
    )
    parser.add_argument(
        "--fix", action="store_true",
        help="Auto-fix mismatches (doc is source of truth for naming)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="With --fix, show what would be changed without modifying files",
    )
    parser.add_argument(
        "--list-missing", action="store_true",
        help="Also list functions missing from source or docs",
    )
    args = parser.parse_args()

    if args.fix:
        args.check_src = True

    src_dir = "src"
    doc_dir = "doc/source"

    if not os.path.isdir(src_dir) or not os.path.isdir(doc_dir):
        print("Error: run from the FLINT root directory", file=sys.stderr)
        sys.exit(1)

    modules = find_modules(src_dir, doc_dir)
    if args.module:
        modules = [(m, h, r) for m, h, r in modules if m == args.module]
        if not modules:
            print(f"Error: module '{args.module}' not found", file=sys.stderr)
            sys.exit(1)

    results = collect_mismatches(modules, src_dir, args.check_src)

    if args.fix:
        fixed, skipped = apply_fixes(results, src_dir, args.dry_run)
        print(f"\n{'='*70}")
        print(f"Fixed: {fixed}, Skipped: {skipped}")
        if not args.dry_run and fixed > 0:
            print("Re-running check to verify...")
            results2 = collect_mismatches(modules, src_dir, True)
            remaining = sum(
                len(r[3]) + len(r[4]) for r in results2
            )
            print(f"Remaining mismatches: {remaining}")
        return

    total_hdr_doc = 0
    total_hdr_doc_checked = 0
    total_hdr_src = 0
    total_hdr_src_checked = 0

    for (mod, header_path, rst_path,
         hdr_doc_mm, hdr_src_mm,
         hdr_funcs, doc_funcs, src_funcs) in results:

        hdr_doc_common = len(set(hdr_funcs.keys()) & set(doc_funcs.keys()))
        total_hdr_doc += len(hdr_doc_mm)
        total_hdr_doc_checked += hdr_doc_common

        if args.check_src:
            hdr_src_common = len(
                set(hdr_funcs.keys()) & set(src_funcs.keys())
            )
            total_hdr_src += len(hdr_src_mm)
            total_hdr_src_checked += hdr_src_common

        if hdr_doc_mm or hdr_src_mm:
            print(f"\n{'='*70}")
            print(f"Module: {mod}")
            print(f"{'='*70}")

        if hdr_doc_mm:
            print(f"\n  --- header vs doc ---")
            print_mismatches(
                hdr_doc_mm, "header", "doc",
                path_a=header_path, path_b=rst_path,
            )

        if hdr_src_mm:
            print(f"\n  --- header vs source ---")
            print_mismatches(
                hdr_src_mm, "header", "source",
                path_a=header_path,
            )

        if args.list_missing:
            only_hdr = sorted(
                set(hdr_funcs.keys()) - set(doc_funcs.keys())
            )
            only_doc = sorted(
                set(doc_funcs.keys()) - set(hdr_funcs.keys())
            )
            if only_hdr or only_doc:
                if not hdr_doc_mm and not hdr_src_mm:
                    print(f"\n{'='*70}")
                    print(f"Module: {mod}")
                    print(f"{'='*70}")
                if only_hdr:
                    print(f"\n  In header but not in docs ({len(only_hdr)}):")
                    for fn in only_hdr:
                        print(f"    {fn}")
                if only_doc:
                    print(f"\n  In docs but not in header ({len(only_doc)}):")
                    for fn in only_doc:
                        print(f"    {fn}")

    print(f"\n{'='*70}")
    print(f"Header vs doc: {total_hdr_doc} mismatches in"
          f" {total_hdr_doc_checked} functions"
          f" across {len(modules)} modules")
    if args.check_src:
        print(f"Header vs src: {total_hdr_src} mismatches in"
              f" {total_hdr_src_checked} functions"
              f" across {len(modules)} modules")

    if total_hdr_doc > 0 or total_hdr_src > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
