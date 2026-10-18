---
title: fexpr.h
header: fexpr.h
description: Symbolic expressions stored in a flat-packed array of words, with construction, inspection, replacement, conversion to rational functions, plain-text and LaTeX output.
category: Exact real and complex numbers
order: 8
---

This module implements symbolic expressions. An expression is either an
*atom* or a *function call*:

- an atom is an integer (such as `0` or `-34`), a symbol (such as `x`,
  `Mul` or `SomeUserNamedSymbol`; symbols should be valid C identifiers)
  or a string (such as `"Hello, world!"`; ASCII only for now);
- a function call $e_0(e_1, \ldots, e_n)$ consists of expressions $e_0,
  \ldots, e_n$, where $e_0$ is usually, but not necessarily, a symbol.

What an expression means depends on how its symbols are interpreted. A
large set of *builtin symbols* with a standard mathematical meaning is
defined in [](fexpr_builtin): with it, `Mul(3, Add(Neg(x), y))` denotes
$3 \cdot (-x + y)$ for variables $x$ and $y$. Builtin symbols are used
for conversions to and from [](ca) numbers (`ca_get_fexpr`,
`ca_set_fexpr`), [](qqbar) numbers (`qqbar_get_fexpr_formula` and
related functions), polynomials and matrices, and for LaTeX output.
The guide [](fexpr/using) shows how to build, inspect and print
expressions.

## Design {#design}

Symbolic expressions are a data exchange format rather than a data
structure for heavy computation: functions that compute with
expressions typically convert them to a dedicated type, such as
polynomials, and convert the result back. For the same reason,
expressions cannot contain arbitrary binary objects such as `fmpz_poly_t`
or `ca_t` values as atoms. To put such an object into an expression,
represent it structurally with atoms and calls (for example a
polynomial as a list of coefficients or as an arithmetic expression),
introduce a placeholder symbol together with a separate table of values,
or, less advisably, encode it in a string.

An `fexpr_t` is stored in *flat-packed* form: a single array of words
(`ulong`) without internal pointers. The first word is a header that
gives the type of the expression and its total size in words. Small
integers, builtin symbols, and symbols and strings of at most
`FEXPR_SMALL_SYMBOL_LEN` characters fit entirely into the header; larger
atoms store their data in the following words. A function call stores
its header followed by $e_0, e_1, \ldots, e_n$ packed contiguously; for
more than four arguments, an index of word offsets follows the header
so that arguments can be located quickly.

This has several consequences:

- memory management is trivial, copying an expression copies an array,
  and testing two expressions for structural equality compares arrays;
- combining expressions into a call concatenates arrays;
- expressions can be exchanged freely between threads, and between
  machines with the same word size and byte order;
- repeated subexpressions cannot share memory (introduce local symbols
  with `Where` for this);
- extracting a subexpression for modification requires a copy, although
  read-only access is free with *views* (see `fexpr_view_arg`);
- modifying a part of an expression rebuilds the whole expression, so
  building an expression incrementally, for example a sum by repeated
  binary `Add`, costs $O(n^2)$; build balanced expressions, or a single
  call such as `Add` with many arguments, instead.

## Types and macros {#types}

::: type
fexpr_struct
fexpr_t

A symbolic expression. The structure has the fields `data` (a pointer to
the array of words) and `alloc` (the number of allocated words, 0 for a
view). An `fexpr_t` is an array of one structure, so it is passed by
reference.
:::

::: type
fexpr_ptr
fexpr_srcptr

Aliases for `fexpr_struct *` and `const fexpr_struct *`, used for arrays
of expressions.
:::

::: type
fexpr_vec_struct
fexpr_vec_t

A vector of expressions with managed length, with the fields `entries`
(an `fexpr_ptr`), `length` and `alloc`.
:::

::: macro
fexpr_vec_entry(vec, i)

Expands to a pointer to the entry with index `i` of `vec`.
:::

::: macro
FEXPR_TYPE_SMALL_INT
FEXPR_TYPE_SMALL_SYMBOL
FEXPR_TYPE_SMALL_STRING
FEXPR_TYPE_BIG_INT_POS
FEXPR_TYPE_BIG_INT_NEG
FEXPR_TYPE_BIG_SYMBOL
FEXPR_TYPE_BIG_STRING
FEXPR_TYPE_CALL0
FEXPR_TYPE_CALL1
FEXPR_TYPE_CALL2
FEXPR_TYPE_CALL3
FEXPR_TYPE_CALL4
FEXPR_TYPE_CALLN

The type codes stored in the low `FEXPR_TYPE_BITS` = 4 bits of the
header word: small atoms stored in the header, large integers (with a
sign) and large symbols and strings stored in the following words, and
function calls with 0 to 4 arguments or with an argument count stored
after the header.
:::

::: macro
FEXPR_TYPE_BITS
FEXPR_TYPE_MASK
FEXPR_TYPE(head)
FEXPR_SIZE(head)
FEXPR_HEADER_SIZE
FEXPR_COEFF_MAX
FEXPR_COEFF_MIN
FEXPR_SMALL_SYMBOL_LEN
FEXPR_BUILTIN_ID(head)

Macros for decoding header words: the type and the size in words of the
expression with header `head`, the size of a header (one word), the
range of integers stored in the header (about $\pm 2^{\bits - 5}$), the
maximum length of a symbol or string stored in the header ($\bits/8 -
1$ characters), and the index of a builtin symbol stored in the header.
:::

## Memory management {#memory}

::: function
void fexpr_init(fexpr_t expr)

Initialises `expr` and sets it to the integer 0.
:::

::: function
void fexpr_clear(fexpr_t expr)

Clears `expr`, freeing its memory.
:::

::: function
fexpr_ptr _fexpr_vec_init(slong len)
void _fexpr_vec_clear(fexpr_ptr vec, slong len)

Allocates an array of `len` expressions initialised to 0, or clears the
`len` expressions of `vec` and frees the array.
:::

::: function
void fexpr_fit_size(fexpr_t expr, slong size)

Ensures that `expr` has room for `size` words, growing the allocation to
at least twice its size if it must grow.
:::

::: function
void fexpr_set(fexpr_t res, const fexpr_t expr)

Sets `res` to a copy of `expr`.
:::

::: function
void fexpr_swap(fexpr_t a, fexpr_t b)

Swaps `a` and `b` in constant time.
:::

## Size information {#size}

::: function
slong fexpr_depth(const fexpr_t expr)

Returns the depth of `expr` as a tree: 1 for an atom, and one more than
the largest depth of the function and the arguments for a call.
:::

::: function
slong fexpr_num_leaves(const fexpr_t expr)

Returns the number of atoms in `expr`, counted with repetition.
:::

::: function
slong fexpr_size(const fexpr_t expr)
slong fexpr_size_bytes(const fexpr_t expr)
slong fexpr_allocated_bytes(const fexpr_t expr)

Returns the size of `expr` in words or in bytes, or the number of bytes
allocated for it. The counts exclude the structure itself; add
`sizeof(fexpr_struct)` for the size of the whole object.
:::

## Comparisons {#comparisons}

::: function
int fexpr_equal(const fexpr_t a, const fexpr_t b)

Returns whether `a` and `b` are the same expression (structural
equality, a comparison of the word arrays).
:::

::: function
int fexpr_equal_si(const fexpr_t expr, slong c)
int fexpr_equal_ui(const fexpr_t expr, ulong c)

Returns whether `expr` is the integer atom $c$.
:::

::: function
ulong fexpr_hash(const fexpr_t expr)

Returns a hash of `expr`; equal expressions have equal hashes.
:::

::: function
int fexpr_cmp_fast(const fexpr_t a, const fexpr_t b)

Compares `a` and `b` in a total order based on the internal
representation, returning $-1$, 0 or 1. The order is fast and suitable
for sorting and binary search, but has no mathematical meaning.
:::

## Atoms {#atoms}

::: function
int fexpr_is_integer(const fexpr_t expr)
int fexpr_is_symbol(const fexpr_t expr)
int fexpr_is_string(const fexpr_t expr)
int fexpr_is_atom(const fexpr_t expr)

Returns whether `expr` is an integer, a symbol, a string, or any atom.
:::

::: function
void fexpr_zero(fexpr_t res)
int fexpr_is_zero(const fexpr_t expr)

Sets `res` to the integer 0, or returns whether `expr` is the integer 0.
:::

::: function
int fexpr_is_neg_integer(const fexpr_t expr)

Returns whether `expr` is a negative integer atom.
:::

::: function
void fexpr_set_si(fexpr_t res, slong c)
void fexpr_set_ui(fexpr_t res, ulong c)
void fexpr_set_fmpz(fexpr_t res, const fmpz_t c)

Sets `res` to the integer atom $c$.
:::

::: function
int fexpr_get_fmpz(fmpz_t res, const fexpr_t expr)

If `expr` is an integer atom, sets `res` to its value and returns 1;
otherwise returns 0 and leaves `res` unchanged.
:::

::: function
void fexpr_set_symbol_builtin(fexpr_t res, slong id)

Sets `res` to the builtin symbol with index `id`, for example
`FEXPR_Pi` (see [](fexpr_builtin)).
:::

::: function
int fexpr_is_builtin_symbol(const fexpr_t expr, slong id)
int fexpr_is_any_builtin_symbol(const fexpr_t expr)

Returns whether `expr` is the builtin symbol with index `id`, or any
builtin symbol.
:::

::: function
void fexpr_set_symbol_str(fexpr_t res, const char * s)

Sets `res` to the symbol with name `s`. If `s` is the name of a builtin
symbol, the builtin symbol is created, so that `fexpr_set_symbol_str(e,
"Pi")` and `fexpr_set_symbol_builtin(e, FEXPR_Pi)` give the same
expression.
:::

::: function
char * fexpr_get_symbol_str(const fexpr_t expr)

Returns a newly allocated string with the name of the symbol `expr`,
which the caller must free with `flint_free`. Throws an exception if
`expr` is not a symbol.
:::

::: function
void fexpr_set_string(fexpr_t res, const char * s)
char * fexpr_get_string(const fexpr_t expr)

Sets `res` to the string atom with content `s`, or returns a newly
allocated copy of the string atom `expr`, which the caller must free
with `flint_free`.
:::

## Input and output {#io}

The plain-text form of an expression is the notation used throughout
this documentation: integers in decimal, symbols by name, strings in
double quotes, and calls as `f(e1, e2, ...)`.

::: function
void fexpr_write(calcium_stream_t stream, const fexpr_t expr)
void fexpr_print(const fexpr_t expr)
char * fexpr_get_str(const fexpr_t expr)

Writes `expr` in plain text to `stream`, prints it to standard output,
or returns it as a newly allocated string to be freed with `flint_free`.

::: details
String atoms are written between double quotes without any escaping.
:::
:::

::: function
char * fexpr_get_decimal_str(const fexpr_t expr, slong digits, ulong flags)

Returns a newly allocated string with a decimal approximation of the
value of `expr` with `digits` significant digits, or `"?"` if the
expression cannot be evaluated numerically. The expression is evaluated
in ball arithmetic with increasing precision (up to a few thousand
bits); integers, the constants `Pi`, `NumberI`, `NumberE`, `Euler`,
`CatalanConstant` and `GoldenRatio`, arithmetic, roots, elementary
functions, and some special functions such as `Gamma`, `Erf`,
`RiemannZeta`, the Bessel and Airy functions are supported. `flags` is
reserved and should be 0.
:::

## LaTeX output {#latex}

::: function
void fexpr_write_latex(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_print_latex(const fexpr_t expr, ulong flags)
char * fexpr_get_str_latex(const fexpr_t expr, ulong flags)

Writes a LaTeX rendering of `expr` to `out`, prints it to standard
output, or returns it as a newly allocated string to be freed with
`flint_free`. Builtin symbols are rendered with standard mathematical
notation (see [](fexpr_builtin)); other calls are rendered as function
applications. String atoms are not escaped.
:::

::: macro
FEXPR_LATEX_SMALL
FEXPR_LATEX_LOGIC

Flags for LaTeX output. `FEXPR_LATEX_SMALL` produces more compact
formulas, in particular fractions written inline as $p/q$ rather than
$\displaystyle \frac{p}{q}$; it is enabled automatically in subscripts,
superscripts and some other parts of formulas. `FEXPR_LATEX_LOGIC` uses
symbols for logical operators such as `Not`, `And` and `Or`, which are
rendered as words by default.
:::

## Function call structure {#calls}

::: function
slong fexpr_nargs(const fexpr_t expr)

Returns the number of arguments $n$ of the call $f(e_1, \ldots, e_n)$,
or $-1$ if `expr` is an atom.
:::

::: function
void fexpr_func(fexpr_t res, const fexpr_t expr)
void fexpr_view_func(fexpr_t res, const fexpr_t expr)

For a call $f(e_1, \ldots, e_n)$, sets `res` to a copy of $f$, or to a
view of it. A view shares the data of `expr` (it has `alloc` = 0): it
must not be initialised before nor cleared after use, must not be
modified, and is valid only as long as `expr` is not modified or
cleared.
:::

::: function
void fexpr_arg(fexpr_t res, const fexpr_t expr, slong i)
void fexpr_view_arg(fexpr_t res, const fexpr_t expr, slong i)

For a call $f(e_1, \ldots, e_n)$, sets `res` to a copy of, or a view of,
the argument $e_{i+1}$ (indices start at 0 and must satisfy $0 \le i <
n$). Throws an exception if `expr` is an atom. Access takes $O(1)$ steps
for calls with more than four arguments, thanks to the index.
:::

::: function
void fexpr_view_next(fexpr_t view)

Advances a view of $e_i$ (or of $f$) in a call to a view of the next
argument $e_{i+1}$ (or of $e_1$). It may be called on a view of the last
argument provided that the result is not used.
:::

::: function
int fexpr_is_builtin_call(const fexpr_t expr, slong id)
int fexpr_is_any_builtin_call(const fexpr_t expr)

Returns whether `expr` is a call $f(\ldots)$ where $f$ is the builtin
symbol with index `id`, or any builtin symbol.
:::

## Construction {#construction}

::: function
void fexpr_call0(fexpr_t res, const fexpr_t f)
void fexpr_call1(fexpr_t res, const fexpr_t f, const fexpr_t x1)
void fexpr_call2(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2)
void fexpr_call3(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3)
void fexpr_call4(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3, const fexpr_t x4)
void fexpr_call_vec(fexpr_t res, const fexpr_t f, fexpr_srcptr args, slong len)

Sets `res` to the call $f(x_1, \ldots, x_n)$, with the arguments given
explicitly or as the array `args` of length $n$ = `len`. The output must
not be the same object as any input.
:::

::: function
void fexpr_call_builtin1(fexpr_t res, slong f, const fexpr_t x)
void fexpr_call_builtin2(fexpr_t res, slong f, const fexpr_t x, const fexpr_t y)

Sets `res` to the call $f(x)$ or $f(x, y)$, where $f$ is the builtin
symbol with index `f`. Here `res` may be the same object as an argument.
:::

## Subexpressions and replacement {#replacement}

::: function
int fexpr_contains(const fexpr_t expr, const fexpr_t x)

Returns whether `x` occurs as a subexpression of `expr` (including
`expr` itself).
:::

::: function
int fexpr_replace(fexpr_t res, const fexpr_t expr, const fexpr_t x, const fexpr_t y)
int fexpr_replace2(fexpr_t res, const fexpr_t expr, const fexpr_t x1, const fexpr_t y1, const fexpr_t x2, const fexpr_t y2)
int fexpr_replace_vec(fexpr_t res, const fexpr_t expr, const fexpr_vec_t xs, const fexpr_vec_t ys)

Sets `res` to `expr` with every occurrence of the subexpression `x`
replaced by `y` (or of the entries of `xs` by the corresponding entries
of `ys`, which must have the same length), simultaneously and without
rescanning the replacements. Returns whether any replacement was made.
`res` may be the same object as `expr` but not as the other inputs.
:::

## Arithmetic expressions {#arithmetic}

::: function
void fexpr_set_fmpq(fexpr_t res, const fmpq_t x)

Sets `res` to the rational number $x$: an integer atom if the
denominator is 1, and otherwise `Div(p, q)`.
:::

::: function
void fexpr_set_arf(fexpr_t res, const arf_t x)
void fexpr_set_d(fexpr_t res, double x)

Sets `res` to an expression for the exact value of the floating-point
number $x = m \cdot 2^e$ ($m$ odd): an integer if $0 \le e \le 20$, a
fraction `Div(m, 2^-e)` if $-8 \le e < 0$, and otherwise `Pow(2, e)`,
`Neg(Pow(2, e))` or `Mul(m, Pow(2, e))`. Infinities are written as
`Infinity` and `Neg(Infinity)`, NaN as `Undefined`.
:::

::: function
void fexpr_set_re_im_d(fexpr_t res, double x, double y)

Sets `res` to an expression for the complex number $x + yi$, built
with `fexpr_set_d`, `NumberI`, `Mul` and `Add` (omitting zero parts).
:::

::: function
void fexpr_neg(fexpr_t res, const fexpr_t a)
void fexpr_add(fexpr_t res, const fexpr_t a, const fexpr_t b)
void fexpr_sub(fexpr_t res, const fexpr_t a, const fexpr_t b)
void fexpr_mul(fexpr_t res, const fexpr_t a, const fexpr_t b)
void fexpr_div(fexpr_t res, const fexpr_t a, const fexpr_t b)
void fexpr_pow(fexpr_t res, const fexpr_t a, const fexpr_t b)

Sets `res` to the expression `Neg(a)`, `Add(a, b)`, `Sub(a, b)`,
`Mul(a, b)`, `Div(a, b)` or `Pow(a, b)`, without any simplification.
:::

::: function
int fexpr_is_arithmetic_operation(const fexpr_t expr)

Returns whether `expr` is a call of one of the arithmetic operators
`Pos`, `Neg`, `Add`, `Sub`, `Mul` or `Div`.
:::

::: function
void fexpr_arithmetic_nodes(fexpr_vec_t nodes, const fexpr_t expr)

Sets `nodes` to the subexpressions of `expr` that are the leaves when
`expr` is viewed as an arithmetic expression: `expr` is built from
integers and the elements of `nodes` by nested applications of `Pos`,
`Neg`, `Add`, `Sub`, `Mul`, `Div` and `Pow` with an integer exponent.
Each node appears once; the nodes are not sorted.
:::

::: function
int fexpr_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fexpr_t expr, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)

Sets `res` to `expr` viewed as a rational function of the
subexpressions in `vars`, which are treated as independent
indeterminates, and returns 1. The length of `vars` must equal the
number of variables of `ctx`; `fexpr_arithmetic_nodes` gives a suitable
vector. Returns 0 if a subexpression is neither an arithmetic operation
nor an entry of `vars`, on overflow (too many terms or too large
exponents), or on division by zero.

::: details
The conversion is formal: algebraic relations between the entries of
`vars` are ignored, so that `Sub(Sqrt(2), Sqrt(2))` gives 0, but
`Sub(Pow(Sqrt(2), 2), 2)` does not, and a denominator that vanishes
only because of such relations is not detected.
:::
:::

::: function
void fexpr_set_fmpz_mpoly(fexpr_t res, const fmpz_mpoly_t poly, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)
void fexpr_set_fmpz_mpoly_q(fexpr_t res, const fmpz_mpoly_q_t frac, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)

Sets `res` to an expression for the polynomial `poly` (or the rational
function `frac`), with the entries of `vars` substituted for the
variables. The length of `vars` must equal the number of variables of
`ctx`; if `vars` is `NULL`, default symbols are used.
:::

::: function
int fexpr_expanded_normal_form(fexpr_t res, const fexpr_t expr, ulong flags)

Sets `res` to the expanded normal form of `expr` as a formal rational
function of its non-arithmetic subexpressions, and returns 1 (0 on
failure, see `fexpr_get_fmpz_mpoly_q`). The nodes are found with
`fexpr_arithmetic_nodes` and sorted with `fexpr_cmp_fast`, `expr` is
converted with `fexpr_get_fmpz_mpoly_q` and converted back with
`fexpr_set_fmpz_mpoly_q`. `flags` is reserved and should be 0.
:::

## Vectors {#vectors}

::: function
void fexpr_vec_init(fexpr_vec_t vec, slong len)
void fexpr_vec_clear(fexpr_vec_t vec)

Initialises `vec` as a vector of `len` expressions equal to 0, or clears
it.
:::

::: function
void fexpr_vec_print(const fexpr_vec_t F)

Prints the entries of `F` to standard output, comma-separated in square
brackets, without a newline.
:::

::: function
void fexpr_vec_swap(fexpr_vec_t x, fexpr_vec_t y)

Swaps `x` and `y` in constant time.
:::

::: function
void fexpr_vec_fit_length(fexpr_vec_t vec, slong len)

Ensures that `vec` has room for `len` entries, without changing its
length.
:::

::: function
void fexpr_vec_set(fexpr_vec_t dest, const fexpr_vec_t src)

Sets `dest` to a copy of `src`.
:::

::: function
void fexpr_vec_append(fexpr_vec_t vec, const fexpr_t f)

Appends a copy of `f` to `vec`.
:::

::: function
slong fexpr_vec_insert_unique(fexpr_vec_t vec, const fexpr_t f)

Returns the index of an entry of `vec` equal to `f`, appending `f` first
if there is none. The search is linear.
:::

::: function
void fexpr_vec_set_length(fexpr_vec_t vec, slong len)

Sets the length of `vec` to `len`, appending zeros or truncating as
needed.
:::

::: function
void _fexpr_vec_sort_fast(fexpr_ptr vec, slong len)

Sorts the `len` entries of `vec` with the order of `fexpr_cmp_fast`.
:::

## Internal functions {#internal}

::: function
slong _fexpr_size(const ulong * expr)
int _fexpr_is_integer(const ulong * expr)
int _fexpr_is_symbol(const ulong * expr)
int _fexpr_is_string(const ulong * expr)
int _fexpr_is_atom(const ulong * expr)

Versions of `fexpr_size` and the atom predicates that take a pointer to
the header word of an expression.
:::

The following functions write a call of a particular builtin symbol in
LaTeX. They are the writer callbacks referenced by the table of builtin
symbols (`fexpr_builtin_table`) and are called by `fexpr_write_latex`;
the argument `expr` must be a call of a symbol for which the writer is
intended, with arguments of the expected form. All have the arguments
of `fexpr_write_latex`.

::: function
void fexpr_write_latex_call(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_subscript(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_subscript_call(calcium_stream_t out, const fexpr_t expr, ulong flags)

Generic forms: a function application $f(x_1, \ldots, x_n)$, a symbol
with subscripts $f_{x_1, \ldots, x_n}$ (used for the indexed Greek and
Latin variable symbols), and a subscripted function application
$f_{x_1}(x_2, \ldots)$.
:::

::: function
void fexpr_write_latex_infix(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_mul(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_div(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_neg_pos(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_add(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_sub(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_pow(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_exp(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_factorial(calcium_stream_t out, const fexpr_t expr, ulong flags)

Operators: infix relations and operators (such as `Equal`, `Less` or
`Element`), products, quotients (as fractions, or inline with
`FEXPR_LATEX_SMALL`), unary minus and plus, sums, differences, powers,
$e^x$, and factorials, with parentheses inserted as needed.
:::

::: function
void fexpr_write_latex_sum_product(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_divsum(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_integral(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_limit(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_residue(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)

Operators binding variables through `For`: sums and products,
divisor sums, integrals, limits, residues and derivatives.
:::

::: function
void fexpr_write_latex_call1_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_call2_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_sub1_call1_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_sub1_call2_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)

Special functions with an optional derivative order as last argument:
$f(x)$ or $f^{(n)}(x)$, $f(x, y)$ or $f^{(n)}(x, y)$, and the same with
the first argument as a subscript, $f_x(y)$ or $f_x(y, z)$.
:::

::: function
void fexpr_write_latex_simple(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_simple2(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_simple2_small(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_misc_special(calcium_stream_t out, const fexpr_t expr, ulong flags)

Notation for functions with a fixed LaTeX form, such as $|x|$,
$\lfloor x \rfloor$, $\binom{n}{k}$ or $\overline{z}$ (the *simple2*
variants for two-argument functions such as binomial coefficients,
*small* writing the arguments in small style), and a few special
functions with unusual notation.
:::

::: function
void fexpr_write_latex_logic(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_setop(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_collection(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_range(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_matrix(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_cases(calcium_stream_t out, const fexpr_t expr, ulong flags)

Logical connectives and quantifiers (`Not`, `And`, `Or`, `Implies`,
`All`, `Exists`, ...), extrema and similar operators over sets
(`Maximum`, `Minimum`, `Infimum`, `ArgMax`, `Poles`, ...), collections
(`List`, `Tuple`, `Set`, also in set-builder form), integer ranges,
matrices, and piecewise definitions (`Cases`).
:::

::: function
void fexpr_write_latex_alg_structure(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_where(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_show_form(calcium_stream_t out, const fexpr_t expr, ulong flags)
void fexpr_write_latex_decimal(calcium_stream_t out, const fexpr_t expr, ulong flags)

Algebraic structures such as polynomial and power series rings
$R[x]$, $R[[x]]$; `Where` expressions with local definitions;
`ShowExpandedNormalForm`; and `Decimal` literals.
:::
