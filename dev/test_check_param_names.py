#!/usr/bin/env python3
"""Tests for check_param_names.py"""

import os
import tempfile
import unittest

from check_param_names import (
    extract_param_name,
    split_params,
    parse_func_signature,
    strip_c_comments_and_preprocessor,
    extract_declarations,
    extract_definitions,
    parse_header_declarations,
    parse_source_definitions,
    parse_rst_functions,
    compare_params,
    find_function_span,
    rename_param_in_range,
    fix_source_file,
    fix_declaration_file,
    collect_mismatches,
    apply_fixes,
)


class TestExtractParamName(unittest.TestCase):
    def test_simple_types(self):
        self.assertEqual(extract_param_name("slong n"), "n")
        self.assertEqual(extract_param_name("ulong x"), "x")
        self.assertEqual(extract_param_name("int flag"), "flag")

    def test_const_qualified(self):
        self.assertEqual(extract_param_name("const fmpz_t x"), "x")
        self.assertEqual(extract_param_name("const padic_ctx_t ctx"), "ctx")

    def test_pointer_types(self):
        self.assertEqual(extract_param_name("ulong * out"), "out")
        self.assertEqual(extract_param_name("FILE * file"), "file")
        self.assertEqual(extract_param_name("nn_srcptr xp"), "xp")
        self.assertEqual(extract_param_name("char * str"), "str")

    def test_flint_unused(self):
        self.assertEqual(
            extract_param_name("const padic_ctx_t FLINT_UNUSED(ctx)"),
            "ctx",
        )
        self.assertEqual(
            extract_param_name("slong FLINT_UNUSED(n)"),
            "n",
        )

    def test_array_brackets(self):
        self.assertEqual(extract_param_name("ulong out[3]"), "out")
        self.assertEqual(extract_param_name("slong perm[10]"), "perm")

    def test_void_and_ellipsis(self):
        self.assertIsNone(extract_param_name("void"))
        self.assertIsNone(extract_param_name("..."))
        self.assertIsNone(extract_param_name(""))

    def test_type_only(self):
        """When param is just a type with no name, return None."""
        self.assertIsNone(extract_param_name("int"))
        self.assertIsNone(extract_param_name("void"))
        self.assertIsNone(extract_param_name("FILE"))

    def test_complex_types(self):
        self.assertEqual(extract_param_name("flint_bitcnt_t bits"), "bits")
        self.assertEqual(extract_param_name("gr_ctx_t ctx"), "ctx")
        self.assertEqual(extract_param_name("const fmpz_poly_t poly"), "poly")


class TestSplitParams(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(
            split_params("fmpz_t f, const fmpz_t g, const fmpz_t h"),
            ["fmpz_t f", "const fmpz_t g", "const fmpz_t h"],
        )

    def test_single_param(self):
        self.assertEqual(split_params("fmpz_t f"), ["fmpz_t f"])

    def test_empty(self):
        self.assertEqual(split_params(""), [])

    def test_void(self):
        self.assertEqual(split_params("void"), ["void"])

    def test_function_pointer(self):
        result = split_params(
            "int (*cmp)(void *, const void *, const void *), slong n"
        )
        self.assertEqual(len(result), 2)
        self.assertIn("(*cmp)", result[0])
        self.assertEqual(result[1], "slong n")

    def test_nested_parens(self):
        result = split_params("acb_calc_func_t func, void * param")
        self.assertEqual(result, ["acb_calc_func_t func", "void * param"])


class TestParseFuncSignature(unittest.TestCase):
    def test_simple(self):
        result = parse_func_signature(
            "void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)"
        )
        self.assertEqual(result, ("fmpz_add", ["f", "g", "h"]))

    def test_no_params(self):
        result = parse_func_signature("void _fmpz_cleanup(void)")
        self.assertEqual(result, ("_fmpz_cleanup", []))

    def test_return_type_pointer(self):
        result = parse_func_signature("mpz_ptr _fmpz_promote(fmpz_t f)")
        self.assertEqual(result, ("_fmpz_promote", ["f"]))

    def test_with_semicolon(self):
        result = parse_func_signature(
            "void fmpz_clear(fmpz_t f);"
        )
        self.assertEqual(result, ("fmpz_clear", ["f"]))

    def test_multiline(self):
        result = parse_func_signature(
            "void fmpz_multi_CRT_ui(fmpz_t output, "
            "nn_srcptr residues, const fmpz_comb_t comb, "
            "fmpz_comb_temp_t ctemp, int sign)"
        )
        self.assertEqual(
            result,
            ("fmpz_multi_CRT_ui", ["output", "residues", "comb",
                                     "ctemp", "sign"]),
        )

    def test_function_pointer_param(self):
        result = parse_func_signature(
            "void qsort_r(void * base, slong n, slong size, "
            "int (*cmp)(void *, const void *, const void *), void * arg)"
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[0], "qsort_r")
        self.assertIn("cmp", result[1])

    def test_no_parens(self):
        self.assertIsNone(parse_func_signature("int x"))

    def test_typedef_function_pointer(self):
        """Function pointer typedefs are filtered by the caller
        (parse_header_declarations skips typedefs), so parse_func_signature
        may return a result — it just won't be a correct function name.
        The important thing is that parse_header_declarations filters these."""
        # Direct call may parse it (incorrectly), but that's OK
        # since the caller filters typedefs before calling this function.
        pass

    def test_warn_unused_result(self):
        """Decorators should be handled by the caller, not parse_func_signature."""
        result = parse_func_signature(
            "int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)"
        )
        self.assertEqual(result, ("gr_poly_set", ["res", "src", "ctx"]))

    def test_flint_unused_param(self):
        result = parse_func_signature(
            "void foo(slong n, const nmod_t FLINT_UNUSED(mod))"
        )
        self.assertEqual(result, ("foo", ["n", "mod"]))

    def test_ellipsis(self):
        result = parse_func_signature(
            "int flint_printf(const char * fmt, ...)"
        )
        self.assertEqual(result, ("flint_printf", ["fmt"]))


class TestStripCCommentsAndPreprocessor(unittest.TestCase):
    def test_block_comment(self):
        text = strip_c_comments_and_preprocessor(
            "/* comment */\nvoid foo(int x);\n"
        )
        self.assertIn("void foo(int x);", text)
        self.assertNotIn("comment", text)

    def test_multiline_block_comment(self):
        text = strip_c_comments_and_preprocessor(
            "/* multi\nline\ncomment */\nvoid bar(void);\n"
        )
        self.assertIn("void bar(void);", text)
        self.assertNotIn("multi", text)

    def test_preprocessor_directive(self):
        text = strip_c_comments_and_preprocessor(
            "#include <stdio.h>\nvoid baz(int n);\n"
        )
        self.assertIn("void baz(int n);", text)
        self.assertNotIn("include", text)

    def test_multiline_macro(self):
        text = strip_c_comments_and_preprocessor(
            "#define FOO(x) \\\n    ((x) + 1)\nvoid quux(int a);\n"
        )
        self.assertIn("void quux(int a);", text)
        self.assertNotIn("FOO", text)

    def test_ifdef_keeps_body(self):
        """Code inside #ifdef blocks should be preserved."""
        text = strip_c_comments_and_preprocessor(
            "#ifdef HAVE_FEATURE\n"
            "void feature_func(int x);\n"
            "#endif\n"
        )
        self.assertIn("void feature_func(int x);", text)

    def test_extern_c_removed(self):
        text = strip_c_comments_and_preprocessor(
            '#ifdef __cplusplus\n'
            'extern "C" {\n'
            '#endif\n'
            'void foo(int x);\n'
            '#ifdef __cplusplus\n'
            '}\n'
            '#endif\n'
        )
        self.assertIn("void foo(int x);", text)
        self.assertNotIn('extern "C"', text)


class TestExtractDeclarations(unittest.TestCase):
    def test_simple_declaration(self):
        text = "void foo(int x);\n"
        decls = extract_declarations(text)
        self.assertEqual(len(decls), 1)
        self.assertIn("void foo(int x)", decls[0][0])

    def test_skips_function_body(self):
        text = (
            "void inline_func(int x)\n"
            "{\n"
            "    return;\n"
            "}\n"
            "void declared_func(int y);\n"
        )
        decls = extract_declarations(text)
        self.assertEqual(len(decls), 1)
        self.assertIn("declared_func", decls[0][0])

    def test_multiple_declarations(self):
        text = (
            "void foo(int a);\n"
            "int bar(slong b, ulong c);\n"
        )
        decls = extract_declarations(text)
        self.assertEqual(len(decls), 2)

    def test_multiline_declaration(self):
        text = (
            "void long_func(int a,\n"
            "    int b,\n"
            "    int c);\n"
        )
        decls = extract_declarations(text)
        self.assertEqual(len(decls), 1)
        sig = decls[0][0]
        self.assertIn("long_func", sig)

    def test_no_declarations_when_no_parens(self):
        text = "int x;\nchar * str;\n"
        decls = extract_declarations(text)
        self.assertEqual(len(decls), 0)


class TestExtractDefinitions(unittest.TestCase):
    def test_simple_definition(self):
        text = (
            "void foo(int x)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        defs = extract_definitions(text)
        self.assertEqual(len(defs), 1)
        self.assertIn("foo", defs[0][0])

    def test_skips_declarations(self):
        text = (
            "void declared_only(int x);\n"
            "void defined(int y)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        defs = extract_definitions(text)
        self.assertEqual(len(defs), 1)
        self.assertIn("defined", defs[0][0])

    def test_return_type_on_separate_line(self):
        text = (
            "void\n"
            "foo(int x, int y)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        defs = extract_definitions(text)
        self.assertEqual(len(defs), 1)
        self.assertIn("foo", defs[0][0])

    def test_multiple_definitions(self):
        text = (
            "void foo(int a)\n"
            "{\n"
            "    return;\n"
            "}\n"
            "\n"
            "int bar(int b)\n"
            "{\n"
            "    return 0;\n"
            "}\n"
        )
        defs = extract_definitions(text)
        self.assertEqual(len(defs), 2)

    def test_nested_braces(self):
        text = (
            "void foo(int x)\n"
            "{\n"
            "    if (x) {\n"
            "        return;\n"
            "    }\n"
            "}\n"
        )
        defs = extract_definitions(text)
        self.assertEqual(len(defs), 1)
        self.assertIn("foo", defs[0][0])

    def test_static_and_nonstatic(self):
        """Both static and non-static definitions should be extracted."""
        text = (
            "static void helper(int x)\n"
            "{\n"
            "    return;\n"
            "}\n"
            "\n"
            "void public_func(int y)\n"
            "{\n"
            "    helper(y);\n"
            "}\n"
        )
        defs = extract_definitions(text)
        self.assertEqual(len(defs), 2)


class TestParseHeaderDeclarations(unittest.TestCase):
    def _write_temp(self, content):
        fd, path = tempfile.mkstemp(suffix=".h")
        os.write(fd, content.encode())
        os.close(fd)
        return path

    def test_basic_header(self):
        path = self._write_temp(
            "/* copyright */\n"
            "#ifndef FOO_H\n"
            "#define FOO_H\n"
            "#ifdef __cplusplus\n"
            'extern "C" {\n'
            "#endif\n"
            "void foo_init(foo_t x);\n"
            "void foo_clear(foo_t x);\n"
            "int foo_add(foo_t r, const foo_t a, const foo_t b);\n"
            "#ifdef __cplusplus\n"
            "}\n"
            "#endif\n"
            "#endif\n"
        )
        try:
            funcs = parse_header_declarations(path)
            self.assertIn("foo_init", funcs)
            self.assertIn("foo_clear", funcs)
            self.assertIn("foo_add", funcs)
            self.assertEqual(funcs["foo_init"][0], ["x"])
            self.assertEqual(funcs["foo_add"][0], ["r", "a", "b"])
        finally:
            os.unlink(path)

    def test_skips_inline(self):
        path = self._write_temp(
            "void foo_declared(int x);\n"
            "static inline void foo_inline(int y)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        try:
            funcs = parse_header_declarations(path)
            self.assertIn("foo_declared", funcs)
            self.assertNotIn("foo_inline", funcs)
        finally:
            os.unlink(path)

    def test_skips_typedefs(self):
        path = self._write_temp(
            "typedef void (*callback_t)(int x);\n"
            "void real_func(int y);\n"
        )
        try:
            funcs = parse_header_declarations(path)
            self.assertIn("real_func", funcs)
            self.assertNotIn("callback_t", funcs)
        finally:
            os.unlink(path)

    def test_multiline_declaration(self):
        path = self._write_temp(
            "void long_func(int a,\n"
            "    int b,\n"
            "    int c);\n"
        )
        try:
            funcs = parse_header_declarations(path)
            self.assertIn("long_func", funcs)
            self.assertEqual(funcs["long_func"][0], ["a", "b", "c"])
        finally:
            os.unlink(path)

    def test_warn_unused_result(self):
        path = self._write_temp(
            "WARN_UNUSED_RESULT int foo(int x, int y);\n"
        )
        try:
            funcs = parse_header_declarations(path)
            self.assertIn("foo", funcs)
            self.assertEqual(funcs["foo"][0], ["x", "y"])
        finally:
            os.unlink(path)


class TestParseSourceDefinitions(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.src_dir = self.tmpdir
        self.mod_dir = os.path.join(self.tmpdir, "mymod")
        os.makedirs(self.mod_dir)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def _write_c(self, name, content):
        path = os.path.join(self.mod_dir, name)
        with open(path, "w") as f:
            f.write(content)

    def test_simple_definition(self):
        self._write_c("add.c",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        funcs = parse_source_definitions(self.src_dir, "mymod")
        self.assertIn("mymod_add", funcs)
        self.assertEqual(funcs["mymod_add"][0], ["r", "a", "b"])

    def test_skips_static(self):
        self._write_c("impl.c",
            "static void helper(int x)\n"
            "{\n"
            "    return;\n"
            "}\n"
            "\n"
            "void mymod_func(int y)\n"
            "{\n"
            "    helper(y);\n"
            "}\n"
        )
        funcs = parse_source_definitions(self.src_dir, "mymod")
        self.assertNotIn("helper", funcs)
        self.assertIn("mymod_func", funcs)

    def test_skips_test_files(self):
        self._write_c("t-add.c",
            "void test_add(int x)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        self._write_c("add.c",
            "void mymod_add(int y)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        funcs = parse_source_definitions(self.src_dir, "mymod")
        self.assertNotIn("test_add", funcs)
        self.assertIn("mymod_add", funcs)

    def test_return_type_on_separate_line(self):
        self._write_c("gcd.c",
            "void\n"
            "mymod_gcd(foo_t f, const foo_t g, const foo_t h)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        funcs = parse_source_definitions(self.src_dir, "mymod")
        self.assertIn("mymod_gcd", funcs)
        self.assertEqual(funcs["mymod_gcd"][0], ["f", "g", "h"])


class TestParseRstFunctions(unittest.TestCase):
    def _write_temp(self, content):
        fd, path = tempfile.mkstemp(suffix=".rst")
        os.write(fd, content.encode())
        os.close(fd)
        return path

    def test_simple_directive(self):
        path = self._write_temp(
            ".. function:: void foo_init(foo_t x)\n"
            "\n"
            "    Initializes x.\n"
        )
        try:
            funcs = parse_rst_functions(path)
            self.assertIn("foo_init", funcs)
            self.assertEqual(funcs["foo_init"][0], ["x"])
        finally:
            os.unlink(path)

    def test_continuation_lines(self):
        path = self._write_temp(
            ".. function:: void foo_init(foo_t x)\n"
            "              void foo_clear(foo_t x)\n"
            "\n"
            "    Init/clear.\n"
        )
        try:
            funcs = parse_rst_functions(path)
            self.assertIn("foo_init", funcs)
            self.assertIn("foo_clear", funcs)
        finally:
            os.unlink(path)

    def test_c_function_directive(self):
        path = self._write_temp(
            ".. c:function:: int bar(slong n, ulong k)\n"
            "\n"
            "    Computes something.\n"
        )
        try:
            funcs = parse_rst_functions(path)
            self.assertIn("bar", funcs)
            self.assertEqual(funcs["bar"][0], ["n", "k"])
        finally:
            os.unlink(path)

    def test_multiple_directives(self):
        path = self._write_temp(
            ".. function:: void alpha(int a)\n"
            "\n"
            "    Does alpha.\n"
            "\n"
            ".. function:: void beta(int b)\n"
            "\n"
            "    Does beta.\n"
        )
        try:
            funcs = parse_rst_functions(path)
            self.assertIn("alpha", funcs)
            self.assertIn("beta", funcs)
        finally:
            os.unlink(path)

    def test_line_numbers(self):
        path = self._write_temp(
            "Title\n"
            "=====\n"
            "\n"
            ".. function:: void first(int x)\n"
            "\n"
            "    First func.\n"
            "\n"
            ".. function:: void second(int y)\n"
            "\n"
            "    Second func.\n"
        )
        try:
            funcs = parse_rst_functions(path)
            self.assertEqual(funcs["first"][1], 4)
            self.assertEqual(funcs["second"][1], 8)
        finally:
            os.unlink(path)


class TestCompareParams(unittest.TestCase):
    def _make_funcs(self, **kwargs):
        """Helper: create funcs dict from name -> param_list mapping."""
        return {
            name: (params, 1, f"void {name}(...)")
            for name, params in kwargs.items()
        }

    def test_no_mismatches(self):
        a = self._make_funcs(foo=["x", "y"], bar=["a"])
        b = self._make_funcs(foo=["x", "y"], bar=["a"])
        self.assertEqual(compare_params(a, b), [])

    def test_name_mismatch(self):
        a = self._make_funcs(foo=["x", "y"])
        b = self._make_funcs(foo=["a", "b"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["func"], "foo")
        self.assertEqual(mm[0]["type"], "param_name")
        self.assertEqual(mm[0]["diffs"], [(0, "x", "a"), (1, "y", "b")])

    def test_count_mismatch(self):
        a = self._make_funcs(foo=["x", "y"])
        b = self._make_funcs(foo=["x", "y", "z"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["type"], "param_count")

    def test_only_common_compared(self):
        a = self._make_funcs(foo=["x"], only_a=["y"])
        b = self._make_funcs(foo=["x"], only_b=["z"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 0)

    def test_partial_mismatch(self):
        a = self._make_funcs(foo=["x", "y", "z"])
        b = self._make_funcs(foo=["x", "b", "z"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["diffs"], [(1, "y", "b")])

    def test_swap_detection(self):
        a = self._make_funcs(foo=["x", "val", "len", "y"])
        b = self._make_funcs(foo=["x", "len", "val", "y"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["swaps"], [(1, 2, "val", "len")])

    def test_no_swap_when_names_differ(self):
        a = self._make_funcs(foo=["a", "b"])
        b = self._make_funcs(foo=["c", "d"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["swaps"], [])

    def test_swap_among_multiple_diffs(self):
        a = self._make_funcs(foo=["w", "val", "len", "z"])
        b = self._make_funcs(foo=["x", "len", "val", "z"])
        mm = compare_params(a, b)
        self.assertEqual(len(mm), 1)
        # Only positions 1,2 are swapped; position 0 is just different
        self.assertEqual(mm[0]["swaps"], [(1, 2, "val", "len")])


class TestEndToEnd(unittest.TestCase):
    """End-to-end tests with temp header + RST + source files."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.src_dir = os.path.join(self.tmpdir, "src")
        self.doc_dir = os.path.join(self.tmpdir, "doc", "source")
        os.makedirs(self.src_dir)
        os.makedirs(self.doc_dir)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def _write(self, path, content):
        full = os.path.join(self.tmpdir, path)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "w") as f:
            f.write(content)
        return full

    def test_header_vs_doc_match(self):
        hdr = self._write("src/mymod.h",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b);\n"
        )
        rst = self._write("doc/source/mymod.rst",
            ".. function:: void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "\n"
            "    Adds a and b.\n"
        )
        h = parse_header_declarations(hdr)
        d = parse_rst_functions(rst)
        mm = compare_params(h, d)
        self.assertEqual(len(mm), 0)

    def test_header_vs_doc_mismatch(self):
        hdr = self._write("src/mymod.h",
            "void mymod_add(foo_t f, const foo_t g, const foo_t h);\n"
        )
        rst = self._write("doc/source/mymod.rst",
            ".. function:: void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "\n"
            "    Adds a and b.\n"
        )
        h = parse_header_declarations(hdr)
        d = parse_rst_functions(rst)
        mm = compare_params(h, d)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["func"], "mymod_add")

    def test_header_vs_source_match(self):
        self._write("src/mymod.h",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b);\n"
        )
        self._write("src/mymod/add.c",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        hdr = os.path.join(self.tmpdir, "src", "mymod.h")
        h = parse_header_declarations(hdr)
        s = parse_source_definitions(
            os.path.join(self.tmpdir, "src"), "mymod"
        )
        mm = compare_params(h, s)
        self.assertEqual(len(mm), 0)

    def test_header_vs_source_mismatch(self):
        self._write("src/mymod.h",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b);\n"
        )
        self._write("src/mymod/add.c",
            "void mymod_add(foo_t f, const foo_t g, const foo_t h)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        hdr = os.path.join(self.tmpdir, "src", "mymod.h")
        h = parse_header_declarations(hdr)
        s = parse_source_definitions(
            os.path.join(self.tmpdir, "src"), "mymod"
        )
        mm = compare_params(h, s)
        self.assertEqual(len(mm), 1)
        self.assertEqual(mm[0]["func"], "mymod_add")
        self.assertEqual(
            mm[0]["diffs"],
            [(0, "r", "f"), (1, "a", "g"), (2, "b", "h")],
        )

    def test_three_way_consistency(self):
        """All three sources match — no mismatches."""
        hdr = self._write("src/mymod.h",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b);\n"
            "void mymod_neg(foo_t r, const foo_t a);\n"
        )
        self._write("doc/source/mymod.rst",
            ".. function:: void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "\n"
            "    Adds.\n"
            "\n"
            ".. function:: void mymod_neg(foo_t r, const foo_t a)\n"
            "\n"
            "    Negates.\n"
        )
        self._write("src/mymod/add.c",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        self._write("src/mymod/neg.c",
            "void mymod_neg(foo_t r, const foo_t a)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        h = parse_header_declarations(hdr)
        d = parse_rst_functions(
            os.path.join(self.tmpdir, "doc", "source", "mymod.rst")
        )
        s = parse_source_definitions(
            os.path.join(self.tmpdir, "src"), "mymod"
        )
        self.assertEqual(compare_params(h, d), [])
        self.assertEqual(compare_params(h, s), [])

    def test_three_way_doc_off(self):
        """Doc has different names, but header and source agree."""
        hdr = self._write("src/mymod.h",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b);\n"
        )
        self._write("doc/source/mymod.rst",
            ".. function:: void mymod_add(foo_t res, const foo_t x, const foo_t y)\n"
            "\n"
            "    Adds.\n"
        )
        self._write("src/mymod/add.c",
            "void mymod_add(foo_t r, const foo_t a, const foo_t b)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        h = parse_header_declarations(hdr)
        d = parse_rst_functions(
            os.path.join(self.tmpdir, "doc", "source", "mymod.rst")
        )
        s = parse_source_definitions(
            os.path.join(self.tmpdir, "src"), "mymod"
        )
        hd = compare_params(h, d)
        hs = compare_params(h, s)
        self.assertEqual(len(hd), 1)  # header vs doc mismatch
        self.assertEqual(len(hs), 0)  # header vs source match


class TestRealFiles(unittest.TestCase):
    """Integration tests against actual FLINT files (skipped if not present)."""

    def setUp(self):
        if not os.path.isfile("src/fmpz.h"):
            self.skipTest("Not running from FLINT root directory")

    def test_fmpz_header_parses(self):
        funcs = parse_header_declarations("src/fmpz.h")
        self.assertGreater(len(funcs), 100)
        self.assertIn("fmpz_add", funcs)
        self.assertEqual(funcs["fmpz_add"][0], ["f", "g", "h"])

    def test_fmpz_doc_parses(self):
        funcs = parse_rst_functions("doc/source/fmpz.rst")
        self.assertGreater(len(funcs), 100)
        self.assertIn("fmpz_add", funcs)

    def test_fmpz_source_parses(self):
        funcs = parse_source_definitions("src", "fmpz")
        self.assertGreater(len(funcs), 50)
        self.assertIn("fmpz_add", funcs)

    def test_acb_header_skips_inlines(self):
        funcs = parse_header_declarations("src/acb.h")
        # acb_add is an inline — should be skipped
        self.assertNotIn("acb_add", funcs)
        # acb_mul is a declaration — should be found
        self.assertIn("acb_mul", funcs)

    def test_gr_poly_header_with_warn_unused(self):
        funcs = parse_header_declarations("src/gr_poly.h")
        self.assertIn("gr_poly_set", funcs)
        self.assertIn("gr_poly_mul", funcs)

    def test_padic_mat_header_parsed(self):
        """padic_mat header should parse correctly."""
        funcs = parse_header_declarations("src/padic_mat.h")
        self.assertIn("padic_mat_is_canonical", funcs,
                      "padic_mat_is_canonical should be in header")


class TestRenameParamInRange(unittest.TestCase):
    def test_simple_rename(self):
        content = "void foo(int x)\n{\n    return x + 1;\n}\n"
        result = rename_param_in_range(content, 0, len(content), "x", "n")
        self.assertEqual(result, "void foo(int n)\n{\n    return n + 1;\n}\n")

    def test_avoids_struct_member(self):
        """Must not rename struct->field or obj.field patterns."""
        content = (
            "void foo(int length)\n"
            "{\n"
            "    x = length + cache->length;\n"
            "}\n"
        )
        result = rename_param_in_range(
            content, 0, len(content), "length", "len"
        )
        self.assertIn("int len)", result)
        self.assertIn("x = len + cache->length", result)

    def test_avoids_dot_member(self):
        content = "void foo(int x)\n{\n    a = x + obj.x;\n}\n"
        result = rename_param_in_range(content, 0, len(content), "x", "n")
        self.assertIn("a = n + obj.x", result)

    def test_word_boundary(self):
        """Must not rename substrings."""
        content = "void foo(int x)\n{\n    int xx = x;\n}\n"
        result = rename_param_in_range(content, 0, len(content), "x", "n")
        self.assertIn("int xx = n", result)
        self.assertIn("int n)", result)


class TestFindFunctionSpan(unittest.TestCase):
    def test_simple_function(self):
        content = "void foo(int x)\n{\n    return;\n}\n"
        span = find_function_span(content, "foo")
        self.assertIsNotNone(span)
        start, end = span
        self.assertIn("foo", content[start:end])
        self.assertIn("return", content[start:end])

    def test_skips_declarations(self):
        content = (
            "void foo(int x);\n"
            "void bar(int y)\n"
            "{\n"
            "    return;\n"
            "}\n"
        )
        span = find_function_span(content, "foo")
        self.assertIsNone(span)
        span = find_function_span(content, "bar")
        self.assertIsNotNone(span)

    def test_return_type_on_separate_line(self):
        content = "void\nfoo(int x)\n{\n    return;\n}\n"
        span = find_function_span(content, "foo")
        self.assertIsNotNone(span)


class TestFixSourceFile(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def _write(self, name, content):
        path = os.path.join(self.tmpdir, name)
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_fix_simple(self):
        path = self._write("func.c",
            "void mymod_foo(foo_t f, const foo_t g)\n"
            "{\n"
            "    bar(f, g);\n"
            "}\n"
        )
        fix_source_file(path, "mymod_foo", [("f", "r"), ("g", "a")])
        with open(path) as f:
            result = f.read()
        self.assertIn("foo_t r, const foo_t a)", result)
        self.assertIn("bar(r, a);", result)

    def test_fix_avoids_struct_members(self):
        path = self._write("func.c",
            "void mymod_foo(slong length)\n"
            "{\n"
            "    x = length + cache->length;\n"
            "}\n"
        )
        fix_source_file(path, "mymod_foo", [("length", "len")])
        with open(path) as f:
            result = f.read()
        self.assertIn("slong len)", result)
        self.assertIn("x = len + cache->length;", result)


class TestFixDeclarationFile(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def _write(self, name, content):
        path = os.path.join(self.tmpdir, name)
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_fix_header_declaration(self):
        path = self._write("mod.h",
            "void mymod_foo(foo_t f, const foo_t g);\n"
            "void mymod_bar(int x);\n"
        )
        fix_declaration_file(path, "mymod_foo", [("f", "r"), ("g", "a")])
        with open(path) as f:
            result = f.read()
        self.assertIn("foo_t r, const foo_t a)", result)
        # mymod_bar should be untouched
        self.assertIn("mymod_bar(int x)", result)


class TestApplyFixes(unittest.TestCase):
    """Test the full fix workflow."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def _write(self, path, content):
        full = os.path.join(self.tmpdir, path)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "w") as f:
            f.write(content)
        return full

    def test_fix_source_to_match_header(self):
        """When header and doc agree, source gets fixed."""
        hdr_path = self._write("src/mymod.h",
            "void mymod_foo(foo_t r, const foo_t a);\n"
        )
        rst_path = self._write("doc/source/mymod.rst",
            ".. function:: void mymod_foo(foo_t r, const foo_t a)\n"
            "\n"
            "    Does foo.\n"
        )
        src_path = self._write("src/mymod/foo.c",
            "void mymod_foo(foo_t f, const foo_t g)\n"
            "{\n"
            "    bar(f, g);\n"
            "}\n"
        )
        src_dir = os.path.join(self.tmpdir, "src")
        doc_dir = os.path.join(self.tmpdir, "doc", "source")
        modules = [("mymod", hdr_path, rst_path)]
        results = collect_mismatches(modules, src_dir, check_src=True)
        fixed, skipped = apply_fixes(results, src_dir)
        self.assertEqual(fixed, 1)
        self.assertEqual(skipped, 0)
        with open(src_path) as f:
            content = f.read()
        self.assertIn("foo_t r, const foo_t a)", content)
        self.assertIn("bar(r, a);", content)

    def test_fix_header_and_source_to_match_doc(self):
        """When header and source agree but doc differs,
        both header and source get fixed to match doc."""
        hdr_path = self._write("src/mymod.h",
            "void mymod_foo(foo_t f, const foo_t g);\n"
        )
        rst_path = self._write("doc/source/mymod.rst",
            ".. function:: void mymod_foo(foo_t r, const foo_t a)\n"
            "\n"
            "    Does foo.\n"
        )
        src_path = self._write("src/mymod/foo.c",
            "void mymod_foo(foo_t f, const foo_t g)\n"
            "{\n"
            "    bar(f, g);\n"
            "}\n"
        )
        src_dir = os.path.join(self.tmpdir, "src")
        modules = [("mymod", hdr_path, rst_path)]
        results = collect_mismatches(modules, src_dir, check_src=True)
        fixed, skipped = apply_fixes(results, src_dir)
        self.assertEqual(fixed, 2)  # header + source
        with open(hdr_path) as f:
            hdr = f.read()
        self.assertIn("foo_t r, const foo_t a)", hdr)
        with open(src_path) as f:
            src = f.read()
        self.assertIn("foo_t r, const foo_t a)", src)
        self.assertIn("bar(r, a);", src)

    def test_fix_struct_member_safety(self):
        """Renaming param must not affect struct->member access."""
        hdr_path = self._write("src/mymod.h",
            "void mymod_foo(slong len);\n"
        )
        rst_path = self._write("doc/source/mymod.rst",
            ".. function:: void mymod_foo(slong len)\n"
            "\n"
            "    Does foo.\n"
        )
        src_path = self._write("src/mymod/foo.c",
            "void mymod_foo(slong length)\n"
            "{\n"
            "    x = length + cache->length;\n"
            "}\n"
        )
        src_dir = os.path.join(self.tmpdir, "src")
        modules = [("mymod", hdr_path, rst_path)]
        results = collect_mismatches(modules, src_dir, check_src=True)
        fixed, skipped = apply_fixes(results, src_dir)
        self.assertEqual(fixed, 1)
        with open(src_path) as f:
            content = f.read()
        self.assertIn("slong len)", content)
        self.assertIn("x = len + cache->length;", content)


if __name__ == "__main__":
    unittest.main()
