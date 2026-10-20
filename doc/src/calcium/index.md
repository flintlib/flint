---
title: calcium.h
header: calcium.h
description: Definitions shared by the exact real and complex number modules: triple-valued logic, output streams, function codes and hashing.
category: Exact real and complex numbers
order: 1
---

The header `calcium.h` collects small definitions that are shared by the
modules for exact real and complex numbers, formerly the separate library
Calcium: [](ca), [](ca_vec), [](ca_poly), [](ca_mat), [](ca_ext),
[](ca_field), [](qqbar) and [](fexpr). It provides output streams that
write either to a file or to a growing string, the codes that identify
the mathematical functions used to build extension numbers, a printing
helper for the triple-valued logic type `truth_t`, and a hash function
for integers. Most programs do not need to include it; it is needed for
`truth_print`, which `ca.h` does not declare.

The user guide page [](guide/types/choosing#exact) gives an overview of the exact
number types; [](ca/using) and [](qqbar/using) introduce the two main
modules.

## Triple-valued logic {#truth}

The exact number modules use two kinds of predicates.

- A predicate `int foo_is_X(const foo_t x, ...)` returns the C booleans
  1 (true) and 0 (false). It either decides a property of the
  representation, such as `ca_is_qq_elem`, or a mathematical property
  that is always decidable for the type, such as `qqbar_is_real`. A few
  predicates return 0 also when truth cannot be certified; this is stated
  in their documentation.
- A predicate `truth_t foo_check_is_X(const foo_t x, ...)` checks a
  mathematical property that may be undecidable, or too expensive to
  decide. It returns a `truth_t`: `T_TRUE`, `T_FALSE` or `T_UNKNOWN`.
  The answers `T_TRUE` and `T_FALSE` are mathematically proved;
  `T_UNKNOWN` means that the implemented methods could not decide.

The type `truth_t` is defined in `gr_types.h` together with the
combinators `truth_and`, `truth_or` and `truth_not` (see [](gr)). Since
`T_TRUE` is the enumerator 0, a `truth_t` must never be used directly as
a C condition: write `t == T_TRUE` or `t == T_FALSE`.

::: function
void truth_print(truth_t t)

Prints `T_TRUE`, `T_FALSE` or `T_UNKNOWN` to standard output, without a
newline. `truth_println` in [](gr) does the same followed by a newline.
:::

## Output streams {#streams}

The printing functions of `ca`, `fexpr` and related modules write to a
stream that holds either a file pointer or a string that grows as text is
appended. This allows the same code to print to a terminal and to build a
string in memory, as in `ca_get_str` and `fexpr_get_str`.

::: type
calcium_stream_struct
calcium_stream_t

Macros for `gr_stream_struct` and `gr_stream_t` (see [](gr)). The
structure has the fields `fp` (a file pointer, or `NULL` for a string
stream), `s` (the string), `len` (its length, excluding the terminating
null byte) and `alloc` (the allocated size of `s`). A `calcium_stream_t`
is an array of one structure, so it is passed by reference.
:::

::: function
void calcium_stream_init_file(calcium_stream_t out, FILE * fp)

Initialises `out` for writing to the file `fp`, which can be `stdout`,
`stderr` or any file opened for writing. Nothing needs to be freed
afterwards. As with all FLINT functions taking a `FILE *`, the
declaration is only visible if `stdio.h` is included before the FLINT
headers.
:::

::: function
void calcium_stream_init_str(calcium_stream_t out)

Initialises `out` for writing to a string in memory, initially empty.
When finished, the caller owns the string `out->s` and must free it
with `flint_free`.
:::

::: function
void calcium_write(calcium_stream_t out, const char * s)

Appends the null-terminated string `s` to `out`.

::: details
When writing to a string, the buffer at least doubles in size whenever
it has to grow, so that appending $n$ characters in total costs $O(n)$
time. A file stream is written with `fputs`, so `s` is written
literally.
:::
:::

::: function
void calcium_write_free(calcium_stream_t out, char * s)

Appends `s` to `out` and then frees `s` with `flint_free`. This is
convenient with functions that return a newly allocated string, such
as `arb_get_str`.
:::

::: function
void calcium_write_si(calcium_stream_t out, slong x)
void calcium_write_fmpz(calcium_stream_t out, const fmpz_t x)

Appends the decimal representation of the integer `x` to `out`.
:::

::: function
void calcium_write_arb(calcium_stream_t out, const arb_t z, slong digits, ulong flags)
void calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)

Appends the ball `z` to `out` as `arb_get_str` formats it, with `digits`
significant digits and the display style `flags` (`ARB_STR_NO_RADIUS`,
`ARB_STR_CONDENSE` and so on, see [](arb)). A complex ball is written as
`a + b*I` or `a - b*I`, or as `a` or `b*I` when the imaginary or real
part is exactly zero. A part with a NaN midpoint (an indeterminate
enclosure) is written as `[+/- inf]`.
:::

## Function codes {#func-codes}

::: type
calcium_func_code

An enumeration of the mathematical functions and constants that can
define an extension number in [](ca_ext): `CA_QQBar` (an algebraic number
given by a `qqbar_t`), the arithmetic operations `CA_Neg`, `CA_Add`,
`CA_Sub`, `CA_Mul`, `CA_Div`, the roots `CA_Sqrt`, `CA_Cbrt`, `CA_Root`,
the complex parts `CA_Floor`, `CA_Ceil`, `CA_Abs`, `CA_Sign`, `CA_Re`,
`CA_Im`, `CA_Arg`, `CA_Conjugate`, the constant `CA_Pi`, the elementary
functions `CA_Sin`, `CA_Cos`, `CA_Exp`, `CA_Log`, `CA_Pow`, `CA_Tan`,
`CA_Cot`, `CA_Cosh`, `CA_Sinh`, `CA_Tanh`, `CA_Coth`, `CA_Atan`,
`CA_Acos`, `CA_Asin`, `CA_Acot`, `CA_Atanh`, `CA_Acosh`, `CA_Asinh`,
`CA_Acoth`, Euler's constant `CA_Euler`, and the special functions
`CA_Gamma`, `CA_LogGamma`, `CA_Psi`, `CA_Erf`, `CA_Erfc`, `CA_Erfi`,
`CA_RiemannZeta` and `CA_HurwitzZeta`. The last enumerator,
`CA_FUNC_CODE_LENGTH`, is the number of codes.

::: details
The enumeration is declared in `ca_types.h`. Only some of the codes
occur in the numbers that [](ca) constructs at present (`CA_QQBar`,
`CA_Sqrt`, `CA_Pow`, `CA_Exp`, `CA_Log`, `CA_Pi`, `CA_Euler`, the complex
parts, the direct trigonometric functions, `CA_Gamma` and the error
functions); the others are reserved. The numerical order of the
enumerators is part of the structural sort order of extension numbers
(see `ca_ext_cmp_repr`).
:::
:::

::: function
const char * calcium_func_name(calcium_func_code func)

Returns a pointer to a static string with the name of the function or
constant `func`, such as `"Exp"` or `"Pi"`, as used when printing
extension numbers. The names agree with the corresponding builtin
symbols of [](fexpr_builtin). For `CA_QQBar`, for the reserved codes
`CA_Cot`, `CA_Coth`, `CA_Acot` and `CA_Acoth`, and for invalid values,
the string `"<unknown function>"` is returned.
:::

## Hashing {#hashing}

::: function
ulong calcium_fmpz_hash(const fmpz_t x)

Returns a hash value for the integer `x`, used for hashing
representations of exact numbers. The current algorithm, which may
change, returns the value of a small `x`, and the least significant limb
of a large `x` negated if `x` is negative (modulo $2^{\bits}$). Equal
integers have equal hashes.
:::
