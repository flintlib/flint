---
title: Contributing to FLINT
summary: How to add a function to FLINT: source, header, test and documentation, the code style, how to write good test code, and how to submit changes.
order: 1
---

FLINT is developed on GitHub at https://github.com/flintlib/flint.
Contributions are submitted as pull requests there; larger changes are
best discussed first in an issue or on the `flint-devel` mailing list
(https://groups.google.com/g/flint-devel). The conventions that users of
FLINT see (types, naming, aliasing, memory management) are described in
[](guide/conventions/types-and-naming) and [](guide/conventions/memory).
The file `code_conventions.md` in the source tree has further details.

## Adding a function {#adding}

Four steps are needed to add a function `module_foo` to an existing
module:

1. Implement it in a new file `src/module/foo.c` (several closely related
   functions may share a file).
2. Add a test program `src/module/test/t-foo.c` and register it in
   `src/module/test/main.c`, with an `#include "t-foo.c"` line and an
   entry `TEST_FUNCTION(module_foo)` in the list of tests, both in
   alphabetical order.
3. Declare the function in `src/module.h`.
4. Document it in `doc/src/module/index.md`, as described in
   [](#documentation) below.

The build system finds new source and test files automatically. Helper
functions used only by `module_foo` are declared `static` in `foo.c`.
Very short functions can instead be defined in the header with the
module's `MODULE_INLINE` macro; the file `src/module/inlines.c` then
compiles them into the library. A test can be omitted for a trivial
helper that is tested indirectly through another function.

A new module also needs an entry in `Makefile.in` and in
`CMakeLists.txt`, and a nonempty `test` directory.

## Code style {#style}

- Use C11, and the GNU extension for inline assembly where applicable.
- Use FLINT's types and functions: `ulong` and `slong` with `UWORD` and
  `WORD` constants, `FLINT_BITS`, `flint_malloc` and `flint_free`,
  `flint_printf`, and `flint_throw` for errors.
- Keep names and argument order consistent with similar existing
  functions: outputs first, the random state second in functions that set
  their first argument to a random value, conversion functions named
  `module1_get_module2` and `module1_set_module2` with `module1` the more
  complex type.
- Mark input arguments passed by reference `const`, but not outputs or
  arguments passed by value.
- Do not use a leading underscore to hide a function; use `static`, or a
  descriptive name ending in `_internal` or `_helper`. A leading
  underscore means a lower-level version of a function, as described in
  [](guide/conventions/types-and-naming#underscore).
- Avoid complex macros.
- Indent by four spaces, without tabs. Put curly braces on lines of their
  own. Surround binary operators with spaces, but not the insides of
  parentheses and brackets. Write `if (cond)`, `for (...)`, `while
  (...)`, but `func(args)`. Format preprocessor conditionals with one
  space after the `#` per level (`# define`).
- Use `/* ... */` comments, not `//`.
- Keep lines to at most 79 characters unless breaking them hurts
  readability.
- Separate logically distinct chunks of code (declarations,
  initialisation, the main loop, cleanup) by a blank line.
- Put the copyright notice and license header at the top of each new
  file, as in the existing files.

The style is that of the `indent` invocation given in
`code_conventions.md`.

## Test code {#tests}

A test program starts with `TEST_FUNCTION_START(module_foo, state)`, which
provides a random state `state`, and ends with
`TEST_FUNCTION_END(state)`. The easiest way to write one is to adapt the
test of a similar function. Most tests in FLINT compute the same result
in two or more different ways and check that the results agree: with
another algorithm, through an identity (for example $(a + b) c = ac +
bc$), by comparing with GMP or another module, or by varying the
precision and other parameters. Precomputed reference values are useful
in some cases. Every test should also check aliasing: call the function
with the output equal to each input in turn.

Keep tests fast. A test program should take less than a second with the
default number of iterations, preferably less than 0.1 seconds; most
functions can be tested well in a few milliseconds. Scale the number of
iterations with `flint_test_multiplier()`, and run the test once with
`FLINT_TEST_MULTIPLIER=10` or more. Think about the corner cases and bias
the random input towards them; the `randtest` functions of the modules do
this to some extent. Check that the test fails when the function is
deliberately broken. If the input space is too large to be explored by
random tests, that is often a sign that the function should be split into
smaller parts that can be tested separately.

Test code must run cleanly under Valgrind (`make valgrind MOD=module`).
The most common cause of memory errors and leaks is a missing or
duplicated `init` or `clear` call; check that the calls match the
declarations in each block, in the test code as well. Tests should also
pass in a build configured with `--enable-assert`.

## Profiling {#profiling}

Profiling code is not needed in most cases, but benchmarks are useful
while developing a new feature. The `TIMEIT_START`, `TIMEIT_STOP` and
`TIMEIT_ONCE_START` macros and `print_memory_usage` of [](profiler) are
convenient for quick measurements. Benchmark programs that are worth
keeping go in `src/module/profile/p-foo.c` and are built with `make
profile MOD=module`; programs that compute tuning cutoffs go in
`src/module/tune`.

## Documentation {#documentation}

The documentation is written in the Markdown dialect of the `fdoc` tool
in `doc/fdoc`, described in `doc/fdoc/README.md`. Each module has an
overview page `doc/src/module/index.md`, which lists every public function
with its prototype in a `::: function` block and a short description, and
may have hand-written pages for important functions and guides. Algorithms
are described in `doc/src/algorithms`. After building FLINT, run

```
make checkdoc          # check all of the documentation and run its examples
make doc               # build the HTML pages in build/doc/html
```

or, in `doc`, for one module,

```
python3 -m fdoc check --only module
python3 -m fdoc examples --only module
```

to compare the documented prototypes with the header (and list
undocumented functions) and to compile and run the examples in the
documentation. Add an entry for your change to the release notes of the
next version in `doc/src/history/flint3/`. Document the preconditions,
the aliasing rules, and the algorithm used, with references.

## Submitting changes {#submitting}

Make sure that the tests of the affected modules pass (`make check
MOD=...`), also with `--enable-assert`, and that the documentation check
reports no problems. Then open a pull request on GitHub, describing the
change and, for performance work, giving timings before and after.
Changes are reviewed before they are merged; the continuous integration
runs the test suite on several platforms, and also builds the
documentation and runs its checks and example programs. Notable changes are recorded in
the release notes ([](history)) when a release is made.
