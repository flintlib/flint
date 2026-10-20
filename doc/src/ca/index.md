---
title: ca.h
header: ca.h
description: Exact real and complex numbers represented as elements of automatically constructed extension fields of the rationals, with special values, proved predicates, elementary and special functions.
category: Exact real and complex numbers
order: 2
---

A `ca_t` represents a real or complex number exactly, in a form that
supports field arithmetic and, as far as possible, proved comparisons. A
number is stored as an element of a field $K = \mathbb{Q}(a_1, \ldots,
a_n)$, where each *extension number* $a_k$ is an algebraic number given
by its minimal polynomial and an isolating ball, a constant such as
$\pi$, or a function such as $\exp$, $\log$ or $\sqrt{\cdot}$ applied to
other `ca_t` numbers. The fields are built automatically: computing
$\sqrt{2} + \pi$ creates the field $\mathbb{Q}(\sqrt{2}, \pi)$ and an
element of it. Algebraic relations between the extension numbers, such
as $a^2 = 2$ for $a = \sqrt{2}$ or $e^{x} e^{-x} = 1$, are found and
used for simplification. Besides numbers, a `ca_t` can hold the special
values *Unknown*, *Undefined*, unsigned infinity and signed infinities,
so that every operation has a result.

Arithmetic in fields with one algebraic generator uses number field
arithmetic ([](nf_elem)); other fields use multivariate rational
functions ([](fmpz_mpoly_q)) reduced modulo an ideal of relations.
Numerical values are computed in ball arithmetic ([](acb)) with
automatically increasing precision. Predicates such as
`ca_check_is_zero` return a `truth_t`: `T_TRUE` and `T_FALSE` are
proved, and `T_UNKNOWN` is returned when the implemented methods fail to
decide. Equality of algebraic numbers is always decided (within
configurable limits on degrees and sizes); equality involving
transcendental numbers is decided in many but not all cases.

The guide [](ca/using) introduces this module with examples, and
[](algorithms/exact/zero-testing#in-practice) explains what can be
decided. The methods are described in
[](algorithms/exact). For algebraic numbers alone, the canonical
representation of [](qqbar) is usually faster; `ca_t` numbers can be
converted to and from `qqbar_t`. Vectors, polynomials and matrices with
`ca_t` entries are provided by [](ca_vec), [](ca_poly) and [](ca_mat),
and the generic rings interface [](gr) offers `ca_t` contexts through
`gr_ctx_init_real_ca`, `gr_ctx_init_complex_ca` and related functions
(see [](gr_domains)).

::: summary algorithms/exact/overview :::

## Representation and conventions {#conventions}

- Every `ca_t` belongs to a context object `ca_ctx_t`, which caches the
  extension numbers and fields that have been created, together with
  numerical enclosures and the ideals of relations. All numbers that
  take part in one operation must belong to the same context, and every
  function takes the context as its last argument. Numbers are moved
  between contexts with `ca_transfer`.
- A context is mutated even by functions that only read numbers, for
  example when a new field is created or a cached enclosure is refined.
  Consequently a context and its numbers must not be used by several
  threads at the same time; use one context per thread.
- The representation of a number is not canonical in general: the same
  number can be an element of different fields, and different elements
  of one field can represent the same number if the ideal of relations
  is incomplete. Exceptions are rational numbers, which are always
  stored as elements of $\mathbb{Q}$, and elements of a field
  $\mathbb{Q}(a)$ with a single algebraic generator, which are stored
  canonically as polynomials in $a$ of degree less than $[\mathbb{Q}(a)
  : \mathbb{Q}]$. Functions with `repr` in their name compare
  representations, the `check` functions compare values.
- Functions allow any output to be the same object as any input unless
  stated otherwise.
- The results of most operations depend on the options stored in the
  context (see [](ca/options)): evaluation limits, the precision for
  numerical tests and the forms used for trigonometric functions.

### Special values {#special-values}

A `ca_t` can hold the following values besides numbers.

- *Unsigned infinity* $\tilde \infty$, the value of $1/0$ and more
  generally of meromorphic functions at their poles.
- *Signed infinity* $c \cdot \infty$, where the sign $c$ is a number
  with $|c| = 1$. The most common ones are $+\infty$, $-\infty$, $+i
  \infty$ and $-i \infty$. Signed infinities express directional limits
  and logarithmic singularities, for example $\log(0) = -\infty$.
- *Undefined*, the value of indeterminate forms such as $0/0$ and of
  functions at essential singularities, such as $\exp(\tilde \infty)$,
  where neither a number nor an infinity makes sense.
- *Unknown*, a meta-value meaning that the result could not be
  computed, either because no algorithm is implemented or because the
  computation would be too expensive. This happens for example when
  dividing by a number that can neither be proved zero nor nonzero.
  Interfaces to other languages may want to raise an exception when a
  result is *Unknown*.

A *number* is a field element and represents an element of
$\mathbb{C}$. The special values are not numbers: they cannot be
extension numbers, and the denominator of a field element always
represents a nonzero number. For any value other than *Unknown*, it is
therefore known exactly whether it is a number, unsigned infinity, a
signed infinity or *Undefined*. The functions below state how they treat
special values; a function that receives *Unknown* returns *Unknown*
unless stated otherwise, and one that receives *Undefined* usually
returns *Undefined*.

## Types {#types}

::: type
ca_struct
ca_t

A number or special value. The structure holds a word `field`, which
encodes a pointer to the field $K$ (an object owned by the context) or
one of the special values, and the element data: an inline `fmpq` for
$K = \mathbb{Q}$, an inline `nf_elem_struct` when $K = \mathbb{Q}(a)$
for a single algebraic $a$, or a pointer to a heap-allocated
`fmpz_mpoly_q_struct` for any other field. A signed infinity stores its
sign $c$ as a field element with a flag in `field`. A `ca_t` is an array
of one structure, so it is passed by reference.

::: details
Since the field is a shallow reference to an object in the context,
creating many elements of the same field is cheap. The layout is
internal; use the functions of this module (or the macros in
[](#internal)) to access the data.
:::
:::

::: type
ca_ptr
ca_srcptr

Aliases for `ca_struct *` and `const ca_struct *`, used for arrays of
numbers such as those of [](ca_vec).
:::

::: type
ca_ctx_struct
ca_ctx_t

A context object. It holds the cache of extension numbers (see
[](ca_ext)) and fields (see [](ca_field)), shortcuts to the fields
$\mathbb{Q}$ and $\mathbb{Q}(i)$, cached contexts for multivariate
polynomials in $1, 2, \ldots$ variables, and the array of options.
A `ca_ctx_t` is an array of one structure.
:::

## Context objects {#context}

::: function
void ca_ctx_init(ca_ctx_t ctx)

Initialises `ctx` and sets all options to their default values. The
fields $\mathbb{Q}$ and $\mathbb{Q}(i)$ are created immediately.
:::

::: function
void ca_ctx_clear(ca_ctx_t ctx)

Clears `ctx`, freeing all cached extension numbers and fields. All
`ca_t` numbers (and vectors, polynomials and matrices of numbers)
belonging to `ctx` must be cleared before.
:::

::: function
void ca_ctx_print(ca_ctx_t ctx)

Prints a description of `ctx` to standard output: the number of cached
extension numbers and fields, and a list of the fields with their
ideals of relations.
:::

## Context options {#options}

The member `options` of a `ca_ctx_t` is an array of `slong` values that
control evaluation limits, simplification and printing. They are read
and written with `ca_ctx_get_option` and `ca_ctx_set_option`, or
directly as `ctx->options[i]`. [](ca/options) explains the
options and their effect in detail.

Options that control evaluation should be set right after
`ca_ctx_init`: fields, relations and enclosures are cached, and a later
change does not apply retroactively to cached objects. Printing options
can be changed at any time.

::: function
slong ca_ctx_get_option(ca_ctx_t ctx, slong i)
void ca_ctx_set_option(ca_ctx_t ctx, slong i, slong value)

Returns, or sets to `value`, the option with index `i`, which must be
one of the constants `CA_OPT_VERBOSE`, ..., `CA_OPT_TRIG_FORM` below.
:::

::: macro
CA_OPT_VERBOSE
CA_OPT_PRINT_FLAGS

Whether to print debugging information about internal decisions
(default 0), and the printing style, a combination of the flags in
[](#printing) (default `CA_PRINT_DEFAULT`).
:::

::: macro
CA_OPT_MPOLY_ORD

The monomial order (`ORD_LEX`, `ORD_DEGLEX` or `ORD_DEGREVLEX`) for the
multivariate polynomials that represent elements of generic fields.
Default `ORD_LEX`. It must be set before any computation.
:::

::: macro
CA_OPT_PREC_LIMIT
CA_OPT_LOW_PREC
CA_OPT_LLL_PREC

The maximum working precision in bits for numerical evaluation (default
4096), which also bounds the size of some exact intermediate results;
the precision for quick numerical checks (default 64); and the precision
for finding integer relations with LLL (default 128).
:::

::: macro
CA_OPT_QQBAR_DEG_LIMIT
CA_OPT_POW_LIMIT
CA_OPT_SMOOTH_LIMIT

The maximum degree of algebraic numbers created internally (default
120); the largest exponent $n$ for which $x^n$ is expanded by repeated
multiplication in generic fields (default 20); and the size in bits of
the factors found by smooth integer factorisation (default 32).
:::

::: macro
CA_OPT_USE_GROEBNER
CA_OPT_GROEBNER_LENGTH_LIMIT
CA_OPT_GROEBNER_POLY_LENGTH_LIMIT
CA_OPT_GROEBNER_POLY_BITS_LIMIT

Whether to compute Gröbner bases of the ideals of relations (default 1),
and the limits for Buchberger's algorithm: the maximum number of basis
elements (default 100), number of terms of a polynomial (default 1000)
and coefficient size in bits (default 10000).
:::

::: macro
CA_OPT_VIETA_LIMIT

The maximum degree $n$ of an algebraic number for which Vieta's
formulas are added to the ideal of relations when all its conjugates are
generators of one field (default 6).
:::

::: macro
CA_OPT_TRIG_FORM
CA_TRIG_DIRECT
CA_TRIG_EXPONENTIAL
CA_TRIG_SINE_COSINE
CA_TRIG_TANGENT

The option `CA_OPT_TRIG_FORM` selects the form in which trigonometric
functions and their inverses are represented: `CA_TRIG_EXPONENTIAL`
(complex exponentials and logarithms, the default), `CA_TRIG_DIRECT`
(the functions themselves), `CA_TRIG_SINE_COSINE` (tangents as
quotients of sines and cosines) or `CA_TRIG_TANGENT` (sines and cosines
in terms of tangents of half the argument).
:::

::: macro
CA_OPT_NUM_OPTIONS

The number of options, that is, the length of the array `options`.
:::

## Memory management {#memory}

::: function
void ca_init(ca_t x, ca_ctx_t ctx)

Initialises `x` for use with the context `ctx` and sets it to the
rational number 0.
:::

::: function
void ca_clear(ca_t x, ca_ctx_t ctx)

Clears `x`, freeing any memory it uses.
:::

::: function
void ca_swap(ca_t x, ca_t y, ca_ctx_t ctx)

Swaps `x` and `y` in constant time.
:::

## Assignment and constants {#assignment}

::: function
void ca_set(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to a copy of `x`.
:::

::: function
void ca_set_si(ca_t res, slong v, ca_ctx_t ctx)
void ca_set_ui(ca_t res, ulong v, ca_ctx_t ctx)
void ca_set_fmpz(ca_t res, const fmpz_t v, ca_ctx_t ctx)
void ca_set_fmpq(ca_t res, const fmpq_t v, ca_ctx_t ctx)

Sets `res` to the integer or rational number `v`, as an element of
$\mathbb{Q}$.
:::

::: function
void ca_set_d(ca_t res, double x, ca_ctx_t ctx)
void ca_set_d_d(ca_t res, double x, double y, ca_ctx_t ctx)

Sets `res` to the exact value of the double `x`, or to $x + yi$. A
finite double becomes a dyadic rational number, $\pm$`inf` becomes
$\pm \infty$, and NaN becomes *Unknown*.
:::

::: function
void ca_transfer(ca_t res, ca_ctx_t res_ctx, const ca_t src, ca_ctx_t src_ctx)

Sets `res`, which belongs to `res_ctx`, to the value of `src`, which
belongs to `src_ctx`. The value is preserved, but the representation
may differ since it depends on the options and cached fields of the
target context.

::: details
Rational numbers are copied directly. Other values are converted to a
symbolic expression with `ca_get_fexpr` (using
`CA_FEXPR_SERIALIZATION`, so that algebraic numbers are transferred
exactly) and rebuilt with `ca_set_fexpr` in the target context. If the
contexts are the same, this is `ca_set`.
:::
:::

::: function
void ca_zero(ca_t res, ca_ctx_t ctx)
void ca_one(ca_t res, ca_ctx_t ctx)
void ca_neg_one(ca_t res, ca_ctx_t ctx)

Sets `res` to 0, 1 or $-1$, as an element of $\mathbb{Q}$.
:::

::: function
void ca_i(ca_t res, ca_ctx_t ctx)
void ca_neg_i(ca_t res, ca_ctx_t ctx)

Sets `res` to the imaginary unit $i$ or to $-i$, as an element of the
number field $\mathbb{Q}(i)$, which is always present in the context.
:::

::: function
void ca_pi(ca_t res, ca_ctx_t ctx)
void ca_pi_i(ca_t res, ca_ctx_t ctx)

Sets `res` to $\pi$, the generator of the field $\mathbb{Q}(\pi)$, or
to $\pi i$, an element of $\mathbb{Q}(\pi, i)$ (it is not an extension
number of its own).
:::

::: function
void ca_euler(ca_t res, ca_ctx_t ctx)

Sets `res` to Euler's constant $\gamma \approx 0.5772$, the generator
of the field $\mathbb{Q}(\gamma)$. It is not known whether $\gamma$ is
irrational, and the implementation assumes no relations involving it.
:::

## Special values {#special}

::: function
void ca_unknown(ca_t res, ca_ctx_t ctx)
void ca_undefined(ca_t res, ca_ctx_t ctx)
void ca_uinf(ca_t res, ca_ctx_t ctx)

Sets `res` to *Unknown*, *Undefined* or unsigned infinity $\tilde
\infty$.
:::

::: function
void ca_pos_inf(ca_t res, ca_ctx_t ctx)
void ca_neg_inf(ca_t res, ca_ctx_t ctx)
void ca_pos_i_inf(ca_t res, ca_ctx_t ctx)
void ca_neg_i_inf(ca_t res, ca_ctx_t ctx)

Sets `res` to the signed infinity $+\infty$, $-\infty$, $+i \infty$ or
$-i \infty$. Other signed infinities are obtained by multiplying one of
these by a number.
:::

## Conversions {#conversions}

::: function
void ca_set_qqbar(ca_t res, const qqbar_t x, ca_ctx_t ctx)

Sets `res` to the algebraic number `x`.

::: details
A rational `x` gives an element of $\mathbb{Q}$, and an irrational `x`
an element of a field $\mathbb{Q}(a)$ with one algebraic generator
(never of a field with several generators: $\sqrt{2} + \sqrt{3}$ is not
recognised as an element of $\mathbb{Q}(\sqrt{2}, \sqrt{3})$). A
quadratic `x` is written as an element of $\mathbb{Q}(\sqrt{N})$, where
$N$ is the discriminant of its minimal polynomial with small square
factors removed by smooth factorisation, so that for example
$\sqrt{12}$ and $1 + \sqrt{3}$ share the field $\mathbb{Q}(\sqrt{3})$.
Higher-degree numbers are made generators of their own field
$\mathbb{Q}(x)$.
:::
:::

::: function
int ca_get_fmpz(fmpz_t res, const ca_t x, ca_ctx_t ctx)
int ca_get_fmpq(fmpq_t res, const ca_t x, ca_ctx_t ctx)
int ca_get_qqbar(qqbar_t res, const ca_t x, ca_ctx_t ctx)

Attempts to evaluate `x` to an explicit integer, rational number or
algebraic number. On success, sets `res` to this number and returns 1;
otherwise returns 0 and leaves `res` unspecified.

::: details
The conversion fails if `x` is a special value or is not an integer,
rational or algebraic number, but it can also fail when `x` is one. An
element of a field with one algebraic generator is always converted. In
other fields, the extension numbers must themselves be evaluable to
algebraic numbers (see `ca_can_evaluate_qqbar`); the numerator and
denominator are then evaluated with [](qqbar) arithmetic, which is
aborted when an intermediate result, or more precisely a resultant
before factorisation, has degree larger than `CA_OPT_QQBAR_DEG_LIMIT`
or coefficients larger than about $10 \cdot$ `CA_OPT_PREC_LIMIT` bits.
These limits can be reached even if the minimal polynomial of `x` is
small.
:::
:::

::: function
int ca_can_evaluate_qqbar(const ca_t x, ca_ctx_t ctx)

Returns whether `ca_get_qqbar` has a chance to succeed, that is, whether
`x` is a number all of whose extension numbers are, recursively,
algebraic numbers, square roots, rational powers, absolute values,
signs, real and imaginary parts, conjugates, floors or ceilings of
evaluable numbers. Nothing is evaluated.
:::

## Symbolic expressions {#fexpr}

Numbers are converted to and from symbolic expressions of type `fexpr_t`
(see [](fexpr)). This is used for printing in other formats such as
LaTeX, for serialisation and for transfers between contexts.

::: function
void ca_get_fexpr(fexpr_t res, const ca_t x, ulong flags, ca_ctx_t ctx)

Sets `res` to a symbolic expression for `x`. Rational numbers and
special values give simple expressions. Other numbers give an expression
`Where(f(a_1, ..., a_n), Def(a_1, e_1), ..., Def(a_n, e_n))` that
defines symbols `a_1`, `a_2`, ... for the extension numbers involved
(recursively) and expresses $x$ as a rational function of them.

::: details
With `flags` = 0, algebraic extension numbers are written as
`Sqrt(n)`-type formulas when they are Gaussian rationals or quadratic,
and otherwise as `PolynomialRootNearest(coeffs, approx)` with a decimal
approximation; this is readable but cannot always be parsed back
unambiguously. With `flags` = `CA_FEXPR_SERIALIZATION`, algebraic
numbers are written exactly with `AlgebraicNumberSerialized` (see
`qqbar_get_fexpr_repr`), so that `ca_set_fexpr` recovers the same
number.
:::
:::

::: macro
CA_FEXPR_SERIALIZATION

Flag for `ca_get_fexpr`: write algebraic numbers in an exact form that
can be read back.
:::

::: function
int ca_set_fexpr(ca_t res, const fexpr_t expr, ca_ctx_t ctx)

Sets `res` to the value of the symbolic expression `expr` and returns
1, or returns 0 if `expr` cannot be evaluated.

::: details
The expression is evaluated by traversing it with the arithmetic and
functions of this module; no symbolic simplification takes place.
Supported are integers, decimal literals (`Decimal`), the constants
`Pi`, `NumberI`, `NumberE`, `Euler`, `GoldenRatio`, the special values
`Infinity`, `UnsignedInfinity`, `Undefined` and `Unknown`, arithmetic
(`Pos`, `Neg`, `Add`, `Sub`, `Mul`, `Div`, `Pow`), `Sqrt`, `Exp`,
`Log`, the trigonometric functions and their inverses, `Abs`, `Sign`,
`Csgn`, `Arg`, `Re`, `Im`, `Conjugate`, `Floor`, `Ceil`, `Gamma`, `Erf`,
`Erfc`, `Erfi`, the algebraic number constructors
`PolynomialRootIndexed`, `PolynomialRootNearest` and
`AlgebraicNumberSerialized`, and `Where` with `Def` for local
definitions. Every expression produced by `ca_get_fexpr` can be read.
:::
:::

## Printing {#printing}

The style of printed output is selected by the option
`CA_OPT_PRINT_FLAGS`, which is a combination of the flags below.
Integers $n$ with $|n| < 10^6$ are always printed as plain literals.

::: macro
CA_PRINT_N
CA_PRINT_DIGITS
CA_PRINT_REPR
CA_PRINT_FIELD
CA_PRINT_DEFAULT
CA_PRINT_DEBUG

Printing flags. `CA_PRINT_N` prints a decimal approximation, correct to
within one unit in the last digit, with six digits or with $d$ digits
if $d \cdot$ `CA_PRINT_DIGITS` is added to the flags. `CA_PRINT_REPR`
prints the symbolic representation, giving local names `a`, `b`, ... to
the extension numbers involved; together with `CA_PRINT_N`, field
elements are printed as `decimal {symbolic}` and extension numbers as
`decimal [definition]`. `CA_PRINT_FIELD` also prints the field of each
element with its ideal of relations, such as `QQ(a,b)/<a^2-2>`.
`CA_PRINT_DEFAULT` is `CA_PRINT_N | CA_PRINT_REPR` and `CA_PRINT_DEBUG`
is `CA_PRINT_N | CA_PRINT_REPR | CA_PRINT_FIELD`.

The decimal approximation is computed at a fixed precision of about
$3.33 d + 64$ bits, without the adaptive refinement of
`ca_get_decimal_str`, and only correct digits are printed. A number
much smaller than the numbers it is computed from is therefore printed
with fewer digits or as `0eN`: Ramanujan's constant
$e^{\pi \sqrt{163}} - 640320^3 - 744 \approx -7.5 \cdot 10^{-13}$ prints as
`0e-5 {a-262537412640768744 where ...}`. Use `ca_get_decimal_str` for a
guaranteed number of significant digits.

When the enclosure is indeterminate at this precision (for example
because the enclosure of a denominator contains zero), its real and
imaginary parts are printed as `[+/- inf]`, even for a well-defined
real number; there is no recomputation at a higher precision. For
example, the reciprocal of Ramanujan's constant prints as
`[+/- inf] + [+/- inf]*I {(1)/(a-262537412640768744) where ...}`, while
`ca_get_decimal_str` gives `-1.333462408e+12`.

::: details
For example, $-7$, $2/3$, $(\sqrt{3}+5)/2$ and $\sqrt{2} (\log(\pi) +
\pi i)$ print as follows.

```text
# CA_PRINT_DEFAULT
-7
0.666667 {2/3}
3.36603 {(a+5)/2 where a = 1.73205 [a^2-3=0]}
1.61889 + 4.44288*I {a*c+b*c*d where a = 1.14473 [Log(3.14159 {b})], b = 3.14159 [Pi], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}

# CA_PRINT_N
-7
0.666667
3.36603
1.61889 + 4.44288*I

# CA_PRINT_N | (CA_PRINT_DIGITS * 20)
-7
0.66666666666666666667
3.3660254037844386468
1.6188925298220266685 + 4.4428829381583662470*I

# CA_PRINT_REPR
-7
2/3
(a+5)/2 where a = [a^2-3=0]
a*c+b*c*d where a = Log(b), b = Pi, c = [c^2-2=0], d = [d^2+1=0]

# CA_PRINT_DEBUG
-7
0.666667 {2/3  in  QQ}
3.36603 {(a+5)/2  in  QQ(a)/<a^2-3> where a = 1.73205 [a^2-3=0]}
1.61889 + 4.44288*I {a*c+b*c*d  in  QQ(a,b,c,d)/<c^2-2, d^2+1> where a = 1.14473 [Log(3.14159 {b  in  QQ(b)})], b = 3.14159 [Pi], c = 1.41421 [c^2-2=0], d = I [d^2+1=0]}
```

A decimal approximation is computed when it is printed, which can be
expensive, and may refine the enclosures cached in the context, which
can affect later computations.
:::
:::

::: function
void ca_print(const ca_t x, ca_ctx_t ctx)
void ca_fprint(FILE * fp, const ca_t x, ca_ctx_t ctx)

Prints `x` to standard output or to the file `fp`, in the style given by
`CA_OPT_PRINT_FLAGS`, without a trailing newline.
:::

::: function
char * ca_get_str(const ca_t x, ca_ctx_t ctx)

Returns a string with the output that `ca_print` would produce. The
caller must free it with `flint_free`.
:::

::: function
void ca_printn(const ca_t x, slong n, ca_ctx_t ctx)

Prints a decimal approximation of `x` with `n` significant digits, as
`ca_get_decimal_str` computes it, to standard output.
:::

## Random generation {#random}

::: function
void ca_randtest_rational(ca_t res, flint_rand_t state, slong bits, ca_ctx_t ctx)

Sets `res` to a random rational number with numerator and denominator
of up to `bits` bits.
:::

::: function
void ca_randtest(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

Sets `res` to a random number obtained by evaluating a random
expression. At depth 0 it is a random rational number, an element of a
quadratic field with coefficients of up to `bits` bits, or a constant
such as $\pi$; for positive `depth` it may instead be an arithmetic
operation or a function applied to operands generated recursively with
depth `depth` $- 1$. The result is always a number, not a special value.
:::

::: function
void ca_randtest_special(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

Sets `res` to a random special value or, with some probability, to a
number generated by `ca_randtest`.
:::

::: function
void ca_randtest_same_nf(ca_t res, flint_rand_t state, const ca_t x, slong bits, slong den_bits, ca_ctx_t ctx)

Sets `res` to a random element of the field of `x`, with numerator
coefficients of up to `bits` bits and a denominator of up to `den_bits`
bits. The number `x` must be an element of a field $\mathbb{Q}(a)$ with
one algebraic generator.
:::

## Representation properties {#representation}

The following functions inspect the representation of a `ca_t`. They
are fast and their answers are certain, but they say nothing about
mathematical properties of the value beyond what the representation
shows. Predicates return 1 for true and 0 for false.

::: function
int ca_equal_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)

Returns whether `x` and `y` have the same representation: they belong
to the same field (with structurally identical generators) and are the
same rational function in it, or they are the same special value.

::: details
*Unknown* is treated as a value here: `ca_equal_repr` of *Unknown* and
*Unknown* is 1. If neither operand is *Unknown*, equal representations
imply equal values; the converse does not hold.
:::
:::

::: function
int ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)

Compares the representations of `x` and `y` in a total order, returning
$-1$, 0 or 1. The order is structural (it compares fields, then
coefficients) and has no meaning as an order of the numbers.
:::

::: function
ulong ca_hash_repr(const ca_t x, ca_ctx_t ctx)

Returns a hash of the representation of `x`; values with equal
representations (`ca_equal_repr`) have equal hashes. All special values
have the same hash.
:::

::: function
int ca_is_unknown(const ca_t x, ca_ctx_t ctx)
int ca_is_special(const ca_t x, ca_ctx_t ctx)

Returns whether `x` is *Unknown*, or any special value (not a number).
:::

::: function
int ca_is_qq_elem(const ca_t x, ca_ctx_t ctx)
int ca_is_qq_elem_zero(const ca_t x, ca_ctx_t ctx)
int ca_is_qq_elem_one(const ca_t x, ca_ctx_t ctx)
int ca_is_qq_elem_integer(const ca_t x, ca_ctx_t ctx)

Returns whether `x` is represented as an element of $\mathbb{Q}$, as the
element 0 or 1 of $\mathbb{Q}$, or as an integer in $\mathbb{Q}$. Since
operations normally demote rational results to $\mathbb{Q}$ (see
`ca_condense_field`), a rational value usually has this representation,
but this is not guaranteed.
:::

::: function
int ca_is_nf_elem(const ca_t x, ca_ctx_t ctx)

Returns whether `x` is represented as an element of a field
$\mathbb{Q}(a)$ with a single algebraic generator, stored as an
[](nf_elem).
:::

::: function
int ca_is_cyclotomic_nf_elem(slong * p, ulong * q, const ca_t x, ca_ctx_t ctx)

Returns whether `x` is represented as an element of a field
$\mathbb{Q}(a)$ whose generator $a$ is a root of unity. In that case,
if `p` and `q` are not `NULL`, sets them to the integers with $0 \le p <
q$ and $\gcd(p, q) = 1$ such that $a = e^{2 \pi i p / q}$.

::: details
Elements of $\mathbb{Q}$ do not count: the result is 0 for a rational
representation. The answer 0 does not prove that `x` lies in no
cyclotomic field, and $q$ need not be the conductor of the smallest
cyclotomic field containing `x`.
:::
:::

::: function
int ca_is_generic_elem(const ca_t x, ca_ctx_t ctx)

Returns whether `x` is represented as an element of a generic field,
that is, it is a number that is neither represented in $\mathbb{Q}$ nor
in a field with a single algebraic generator.
:::

## Predicates and comparisons {#predicates}

The following functions check mathematical properties, which need not be
decidable. They return a `truth_t`; `T_TRUE` and `T_FALSE` are proved
answers. [](ca/check-is-zero) describes how zero tests and
equality tests work and what to expect from them.

::: function
truth_t ca_check_is_number(const ca_t x, ca_ctx_t ctx)

Returns `T_TRUE` if `x` is a number (a field element), `T_FALSE` if it
is an infinity or *Undefined*, and `T_UNKNOWN` if it is *Unknown*.
:::

::: function
truth_t ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_one(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_neg_one(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_i(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_neg_i(const ca_t x, ca_ctx_t ctx)

Checks whether `x` is equal to the number $0$, $1$, $-1$, $i$ or $-i$.
Special values are not equal to any number; *Unknown* gives `T_UNKNOWN`.

::: details
`ca_check_is_zero` first inspects the representation, then evaluates `x`
in ball arithmetic at increasing precision up to `CA_OPT_PREC_LIMIT`,
which proves $x \ne 0$ when the ball excludes zero. To prove $x = 0$ it
tries an exact evaluation as an algebraic number, a rewriting of the
extension numbers in a normal form based on exponentials and
logarithms, and a factorisation of the numerator. The other functions
subtract the constant and call `ca_check_is_zero`. See
[](ca/check-is-zero) and [](algorithms/exact/zero-testing).
:::
:::

::: function
truth_t ca_is_zero_check_fast(const ca_t x, ca_ctx_t ctx)

The cheap first stage of `ca_check_is_zero`: decides the question from
the representation alone for elements of $\mathbb{Q}$ and of fields with
one algebraic generator, returns `T_FALSE` for special values other than
*Unknown*, and `T_UNKNOWN` otherwise. Nothing is evaluated.
:::

::: function
truth_t ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_rational(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_integer(const ca_t x, ca_ctx_t ctx)

Checks whether `x` is an algebraic number, a rational number or an
integer. Special values give `T_FALSE` (or `T_UNKNOWN` for *Unknown*).

::: details
`ca_check_is_algebraic` returns `T_TRUE` when every extension number is
algebraic by construction (see `ca_can_evaluate_qqbar`), and otherwise
`T_UNKNOWN`: it never proves transcendence. The other two functions use
numerical evaluation to disprove and exact evaluation to an algebraic
number to prove.
:::
:::

::: function
truth_t ca_check_is_real(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_negative_real(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_imaginary(const ca_t x, ca_ctx_t ctx)

Checks whether `x` is a real number, a negative real number, or a
purely imaginary number ($\operatorname{Re}(x) = 0$; zero counts as
imaginary). These functions return `T_FALSE` for all special values,
including infinities with a real or imaginary sign, and `T_UNKNOWN` for
*Unknown*.

::: details
Numerical evaluation disproves, and also proves when the ball has an
exactly zero imaginary (or real) part, which happens when the
expression is real by construction. Otherwise `ca_check_is_real` tries
to prove $x = \overline{x}$ by comparing `x` with its deep conjugate
(`ca_conj_deep`), and an exact evaluation as an algebraic number.
:::
:::

::: function
truth_t ca_check_is_undefined(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_infinity(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_uinf(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_signed_inf(const ca_t x, ca_ctx_t ctx)

Checks whether `x` is *Undefined*, any infinity, unsigned infinity, or
any signed infinity. These depend only on the kind of value, so the
answer is `T_UNKNOWN` only for *Unknown*.
:::

::: function
truth_t ca_check_is_pos_inf(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_neg_inf(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_pos_i_inf(const ca_t x, ca_ctx_t ctx)
truth_t ca_check_is_neg_i_inf(const ca_t x, ca_ctx_t ctx)

Checks whether `x` is the signed infinity $+\infty$, $-\infty$, $+i
\infty$ or $-i \infty$. For a signed infinity this tests the sign with
`ca_check_equal`, which may give `T_UNKNOWN`.
:::

::: function
truth_t ca_check_equal(const ca_t x, const ca_t y, ca_ctx_t ctx)

Checks whether $x = y$. Two numbers are compared by testing $x - y$ for
zero, after a quick check of the representations; special values are
equal if they are the same kind of value (with equal signs for signed
infinities). The result is `T_UNKNOWN` if either operand is *Unknown*,
and may be `T_UNKNOWN` if `x` and `y` cannot be distinguished
numerically and their equality cannot be proved.
:::

::: function
truth_t ca_check_lt(const ca_t x, const ca_t y, ca_ctx_t ctx)
truth_t ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx)
truth_t ca_check_gt(const ca_t x, const ca_t y, ca_ctx_t ctx)
truth_t ca_check_ge(const ca_t x, const ca_t y, ca_ctx_t ctx)

Checks whether $x < y$, $x \le y$, $x > y$ or $x \ge y$. Only real
numbers and $\pm \infty$ are comparable. If either operand is proved not
comparable (a non-real number, unsigned infinity, a signed infinity
other than $\pm \infty$, or *Undefined*) the result is `T_FALSE`, not
`T_UNKNOWN`.

::: details
Real numbers are compared by evaluating both in ball arithmetic at
increasing precision; when the balls overlap at the maximum precision,
the numbers are tested for equality, which decides `ca_check_le` and
`ca_check_ge` if they are equal and otherwise gives `T_UNKNOWN`.
:::
:::

## Arithmetic {#arithmetic}

Arithmetic on numbers is field arithmetic: the operands are first
coerced to a common field (see `ca_merge_fields`), the result is reduced
modulo the ideal of relations of that field, and a result that is
rational is demoted to $\mathbb{Q}$. In the rules for special values,
$c \infty$ and $d \infty$ denote signed infinities ($|c| = |d| = 1$) and
$z$ a number.

::: function
void ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $-x$. The negation of $c \infty$ is $(-c) \infty$; the
other special values are unchanged.
:::

::: function
void ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
void ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
void ca_add_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
void ca_add_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
void ca_add(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `res` to $x + y$.

::: details
For special values:

- $c \infty + d \infty = c \infty$ if $c = d$, and *Undefined* if $c
  \ne d$ (*Unknown* if this cannot be decided);
- $\tilde \infty + \tilde \infty$ and $\tilde \infty + c \infty$ are
  *Undefined*;
- $c \infty + z = c \infty$ and $\tilde \infty + z = \tilde \infty$;
- *Undefined* plus any value, including *Unknown*, is *Undefined*;
  otherwise *Unknown* plus any value is *Unknown*.
:::
:::

::: function
void ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
void ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
void ca_sub_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
void ca_sub_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
void ca_fmpq_sub(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx)
void ca_fmpz_sub(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx)
void ca_ui_sub(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx)
void ca_si_sub(ca_t res, slong x, const ca_t y, ca_ctx_t ctx)
void ca_sub(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `res` to $x - y$, with the same rules for special values as $x +
(-y)$.
:::

::: function
void ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
void ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
void ca_mul_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
void ca_mul_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
void ca_mul(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `res` to $x y$.

::: details
For special values:

- $c \infty \cdot d \infty = (c d) \infty$;
- $\tilde \infty \cdot \tilde \infty = \tilde \infty \cdot c \infty =
  \tilde \infty$;
- $c \infty \cdot z = (c \operatorname{sgn}(z)) \infty$ and $\tilde
  \infty \cdot z = \tilde \infty$ if $z \ne 0$;
- $c \infty \cdot 0$ and $\tilde \infty \cdot 0$ are *Undefined*;
- *Undefined* times any value is *Undefined*; otherwise *Unknown* times
  any value is *Unknown*.

When it cannot be decided whether $z = 0$, the product of $z$ and an
infinity is *Unknown*.
:::
:::

::: function
void ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $1/x$. In a field with one algebraic generator the result
has a rational denominator; in a generic field the denominator is not
rationalised.

::: details
The rules for zero and special values are $1/0 = \tilde \infty$, $1/(c
\infty) = 1/\tilde \infty = 0$, and *Undefined* and *Unknown* map to
themselves. If it cannot be decided whether $x = 0$, the result is
*Unknown*.
:::
:::

::: function
void ca_fmpq_div(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx)
void ca_fmpz_div(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx)
void ca_ui_div(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx)
void ca_si_div(ca_t res, slong x, const ca_t y, ca_ctx_t ctx)
void ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
void ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
void ca_div_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
void ca_div_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
void ca_div(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `res` to $x / y$, which is $x \cdot (1/y)$ for numbers.

::: details
For special values and division by zero:

- $z / 0 = \tilde \infty$ for $z \ne 0$, and $0/0$ is *Undefined*;
- $c \infty / z = (c / \operatorname{sgn}(z)) \infty$ for $z \ne 0$, and
  $c \infty / 0 = \tilde \infty / 0 = \tilde \infty$;
- $z / (c \infty) = z / \tilde \infty = 0$;
- a quotient of two infinities is *Undefined*;
- *Undefined* divided by or into any value is *Undefined*; otherwise
  *Unknown* gives *Unknown*.

If the case cannot be decided because a zero test fails, the result is
*Unknown*.
:::
:::

::: function
void ca_dot(ca_t res, const ca_t initial, int subtract, ca_srcptr x, slong xstep, ca_srcptr y, slong ystep, slong len, ca_ctx_t ctx)

Sets `res` to $s + (-1)^{\mathit{subtract}} \sum_{i=0}^{\mathit{len}-1}
x_{i \cdot \mathit{xstep}} \, y_{i \cdot \mathit{ystep}}$.

::: details
The initial value $s$ is `initial`, or 0 if `initial` is `NULL`.
`subtract` must be 0 or 1. A nonpositive `len` gives $s$. The steps may
be negative, to traverse arrays backwards from the given pointers. The
output may be the same object as `initial` but must not be an entry of
`x` or `y`. The products are accumulated one at a time with `ca_mul`
and `ca_add`.
:::
:::

::: function
void ca_fmpz_poly_evaluate(ca_t res, const fmpz_poly_t poly, const ca_t x, ca_ctx_t ctx)
void ca_fmpq_poly_evaluate(ca_t res, const fmpq_poly_t poly, const ca_t x, ca_ctx_t ctx)

Sets `res` to the polynomial `poly` evaluated at `x`.
:::

::: function
void ca_fmpz_mpoly_evaluate_horner(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
void ca_fmpz_mpoly_evaluate(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)

Sets `res` to the multivariate polynomial `f` (in the context `mctx`)
evaluated at the vector `x`, which must have one entry per variable.
The *horner* version uses a multivariate Horner scheme, and the default
version currently calls the *horner* version.
:::

::: function
void ca_fmpz_mpoly_q_evaluate(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)

Sets `res` to the multivariate rational function `f` evaluated at `x`,
dividing the values of the numerator and the denominator.
:::

::: function
void ca_fmpz_mpoly_q_evaluate_no_division_by_zero(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
void ca_inv_no_division_by_zero(ca_t res, const ca_t x, ca_ctx_t ctx)

Like `ca_fmpz_mpoly_q_evaluate` and `ca_inv`, but assume without
checking that the value of the denominator, respectively `x`, is
nonzero, which avoids a zero test. If the assumption is false, the
result is meaningless; `ca_inv_no_division_by_zero` throws an exception
if `x` is manifestly zero.
:::

## Powers and roots {#powers}

::: function
void ca_sqr(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $x^2$.
:::

::: function
void ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
void ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
void ca_pow_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
void ca_pow_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
void ca_pow(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `res` to $x^y$, with the principal branch $x^y = \exp(y \log x)$
for non-integer $y$.

::: details
Integer powers of rational numbers and of number field elements are
computed exactly unless the size of the result would exceed about
`CA_OPT_PREC_LIMIT` bits; in generic fields, $x^n$ is expanded by
repeated multiplication for $|n| \le$ `CA_OPT_POW_LIMIT`. Exponents
$n/2$ are computed via `ca_sqrt`, and $(z^a)^b$ is simplified to
$z^{ab}$ when this is valid. Otherwise the result is the generator of a
new field $\mathbb{Q}(x^y)$, which requires that $x$ is proved nonzero
unless $y$ is a positive rational number; if that fails, the result is
*Unknown*.

For special bases and integer exponents $n$: $\tilde \infty^n$ is
$\tilde \infty$, 0 or *Undefined* for $n > 0$, $n < 0$ and $n = 0$;
$(c \infty)^n = c^n \infty$ for $n > 0$, 0 for $n < 0$ and *Undefined*
for $n = 0$; $0^n = \tilde \infty$ for $n < 0$. Non-integer powers
involving special values give *Unknown*.
:::
:::

::: function
void ca_pow_si_arithmetic(ca_t res, const ca_t x, slong n, ca_ctx_t ctx)

Sets `res` to $x^n$, always by field arithmetic (binary exponentiation),
whereas `ca_pow_si` may create $x^n$ as a new extension number for large
$n$. Special values of `x` give *Unknown*.
:::

::: function
void ca_sqrt_inert(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_sqrt_nofactor(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_sqrt_factor(ca_t res, const ca_t x, ulong flags, ca_ctx_t ctx)
void ca_sqrt(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the principal square root of `x`, with
$\operatorname{Re}(\sqrt{x}) \ge 0$ and $\sqrt{x} = i \sqrt{|x|}$ for
negative real $x$.

::: details
The *inert* version returns the generator of the field
$\mathbb{Q}(\sqrt{x})$ without trying to simplify. The *nofactor*
version tries cheap simplifications, in particular an exact evaluation
of $\sqrt{x}$ as an algebraic number of moderate degree. The *factor*
version also writes $x = A^2 B$ by factoring the representation of $x$
(with the factorisation options `flags`, see `ca_factor`) and returns
$\pm A \sqrt{B}$ with the sign chosen to give the principal root; the
factorisation is a heuristic and $B$ need not be minimal. `ca_sqrt`
currently calls `ca_sqrt_factor` with `CA_FACTOR_POLY_SQF |
CA_FACTOR_ZZ_SMOOTH`, so that for example $\sqrt{12 \pi}$ becomes $2
\sqrt{3 \pi}$.

For special values: $\sqrt{c \infty} = \sqrt{c} \, \infty$,
$\sqrt{\tilde \infty} = \tilde \infty$, and *Undefined* and *Unknown*
map to themselves.
:::
:::

::: function
void ca_sqrt_ui(ca_t res, ulong n, ca_ctx_t ctx)

Sets `res` to $\sqrt{n}$.
:::

## Complex parts {#complex-parts}

The functions in this section return an exact algebraic number when the
argument can be evaluated as one (within `CA_OPT_QQBAR_DEG_LIMIT`), and
otherwise usually the generator of a new field such as
$\mathbb{Q}(|x|)$.

::: function
void ca_abs(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $|x|$. For special values, $|c \infty| = |\tilde \infty|
= +\infty$, and *Undefined* and *Unknown* map to themselves.
:::

::: function
void ca_sgn(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the complex sign of `x`, which is $x/|x|$ for $x \ne 0$
and 0 for $x = 0$. For special values, $\operatorname{sgn}(c \infty) =
c$, $\operatorname{sgn}(\tilde \infty)$ is *Undefined*, and *Undefined*
and *Unknown* map to themselves.
:::

::: function
void ca_csgn(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the extension of the real sign function that is 1 in the
open right half-plane, $-1$ in the open left half-plane, and the sign of
the imaginary part on the imaginary axis; equivalently
$\operatorname{csgn}(z) = z / \sqrt{z^2}$ for $z \ne 0$, and
$\operatorname{csgn}(0) = 0$. For a signed infinity
$\operatorname{csgn}(c \infty) = \operatorname{csgn}(c)$; unsigned
infinity gives *Undefined*. The result is *Unknown* when the required
sign tests fail.
:::

::: function
void ca_arg(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the argument (phase) of `x` in $(-\pi, \pi]$, with
$\arg(0) = 0$. If $\operatorname{sgn}(x)$ is a root of unity
$e^{\pi i p/q}$, the result is the rational multiple $p \pi / q$. For
special values, $\arg(c \infty) = \arg(c)$, $\arg(\tilde \infty)$ is
*Undefined*, and *Undefined* and *Unknown* map to themselves.
:::

::: function
void ca_re(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_im(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the real or imaginary part of `x`. Every infinity gives
*Undefined*, including real and imaginary infinities.
:::

::: function
void ca_conj_deep(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_conj_shallow(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_conj(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the complex conjugate $\overline{x}$.

::: details
The *shallow* version conjugates in place when this is trivial (for
example in $\mathbb{Q}(i)$ or when $x$ is real by construction) and
otherwise creates the extension number $\overline{x}$. The *deep* version
conjugates every extension number of the field of $x$, recursively
(using $\overline{\exp(z)} = \exp(\overline{z})$ and similar rules, with
care at branch cuts), and evaluates the conjugated rational function.
`ca_conj` returns $x$ or $-x$ if $x$ is proved real or imaginary, and
otherwise uses the deep version. All special values give *Unknown* at
present.
:::
:::

::: function
void ca_floor(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_ceil(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $\lfloor x \rfloor$ or $\lceil x \rceil$; for a nonreal
number, of its real part. The integer is determined numerically when the
real part is not too close to an integer; otherwise the result is the
generator of the field $\mathbb{Q}(\lfloor x \rfloor)$. All infinities
give *Undefined*.
:::

## Exponentials and logarithms {#exp-log}

::: function
void ca_exp(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $e^x$. In the generic case the result is the generator of
$\mathbb{Q}(e^x)$.

::: details
The following simplifications are made automatically: $e^0 = 1$,
$e^{\log z} = z$, $e^{(p/q) \log z} = z^{p/q}$, $e^{\pi i p / q}$ is
the algebraic root of unity when $q \le$ `CA_OPT_QQBAR_DEG_LIMIT`, and
$e^{a \log z + s} = z^a e^{s}$ when the representation of the exponent
contains a term $a \log z$ with rational $a$.

For special values, $e^{+\infty} = +\infty$, $e^{-\infty} = 0$, and
$e^{\tilde \infty}$ and $e^{\text{Undefined}}$ are *Undefined*; other
signed infinities give *Unknown*.
:::
:::

::: function
void ca_log(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the principal branch of the natural logarithm, with $-\pi
< \operatorname{Im}(\log x) \le \pi$. In the generic case the result is
the generator of $\mathbb{Q}(\log x)$.

::: details
The following simplifications are made automatically: $\log 1 = 0$;
$\log(e^z) = z + 2 \pi i k$, $\log(\sqrt{z}) = \frac{1}{2} \log z + 2
\pi i k$ and $\log(z^a) = a \log z + 2 \pi i k$ with the integer $k$
determined numerically; $\log x = \log(-x) + \pi i$ for negative real
$x$; and logarithms of roots of unity are rational multiples of $\pi
i$. Integer relations between logarithms in one field, such as $\log 8
= 3 \log 2$, are found when the field is created (see
[](algorithms/exact/field-towers)).

For special values and zero, $\log 0 = -\infty$ (*Unknown* if it cannot
be decided whether $x = 0$), the logarithm of any infinity is
$+\infty$, and *Undefined* and *Unknown* map to themselves.
:::
:::

## Trigonometric functions {#trig}

::: function
void ca_sin_cos_exponential(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
void ca_sin_cos_direct(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
void ca_sin_cos_tangent(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
void ca_sin_cos(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)

Sets `res1` to $\sin x$ and `res2` to $\cos x$. Either output may be
`NULL` if only the other function is wanted.

::: details
The *exponential* version writes $\sin x = (e^{ix} - e^{-ix})/(2i)$ and
$\cos x = (e^{ix} + e^{-ix})/2$; rational multiples of $\pi$ then give
elements of cyclotomic fields. The *direct* version returns the
generators of $\mathbb{Q}(\sin x)$ and $\mathbb{Q}(\cos x)$ (perhaps
after symmetry transformations, which may exchange the two functions),
and the *tangent* version expresses both in terms of $t = \tan(x/2)$ as
$2t/(1+t^2)$ and $(1-t^2)/(1+t^2)$. In the last two versions, very
simple algebraic values become elements of real number fields.
`ca_sin_cos` uses the form selected by `CA_OPT_TRIG_FORM`, exponential
by default, which works best for simplification and zero testing
although it introduces complex numbers in real computations;
`CA_TRIG_SINE_COSINE` selects the direct form.

For special values: $\sin(\pm i \infty) = \pm i \infty$, $\cos(\pm i
\infty) = +\infty$, the other infinities give *Undefined*, and
*Undefined* and *Unknown* map to themselves.
:::
:::

::: function
void ca_sin(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_cos(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $\sin x$ or $\cos x$, computed with `ca_sin_cos`.
:::

::: function
void ca_tan_sine_cosine(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_tan_exponential(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_tan_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_tan(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $\tan x$, as a quotient of a sine and a cosine, in terms
of complex exponentials, or directly as the generator of $\mathbb{Q}(\tan
x)$ (possibly after a transformation of the argument). `ca_tan` uses the
form selected by `CA_OPT_TRIG_FORM`. Simple algebraic values simplify to
elements of number fields.

::: details
For special values: $\tan$ has the value $\tilde \infty$ at its poles
$(n + \frac{1}{2}) \pi$; $\tan(c \infty) = i$ if
$\operatorname{Im}(c) > 0$ and $-i$ if $\operatorname{Im}(c) < 0$;
$\tan(\pm \infty)$, $\tan(\tilde \infty)$ and $\tan(\text{Undefined})$
are *Undefined*; $\tan(\text{Unknown})$ and $\cot(\text{Unknown})$ are
*Unknown*.
:::
:::

::: function
void ca_cot(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $\cot x = 1/\tan x$.
:::

::: function
void ca_atan_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_atan_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_atan(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the principal branch of $\operatorname{atan} x$.

::: details
The *logarithm* version writes $\operatorname{atan} x = \frac{i}{2}
(\log(1 - ix) - \log(1 + ix))$ (in a form that respects the branch
cuts), the *direct* version returns the generator of
$\mathbb{Q}(\operatorname{atan} x)$, possibly after transforming the
argument. Simple algebraic arguments such as $1$ or $\sqrt{3}$ give
rational multiples of $\pi$. `ca_atan` uses the logarithm form unless
`CA_OPT_TRIG_FORM` is `CA_TRIG_DIRECT`.

For special values: $\operatorname{atan}(\pm i) = \pm i \infty$,
$\operatorname{atan}(c \infty) = \operatorname{csgn}(c) \pi / 2$, and
$\operatorname{atan}(\tilde \infty)$ is *Undefined*.
:::
:::

::: function
void ca_asin_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_acos_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_asin_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_acos_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_asin(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_acos(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the principal branch of $\operatorname{asin} x$ or
$\operatorname{acos} x$.

::: details
The *logarithm* versions use $\operatorname{asin} x = -i \log(ix +
\sqrt{1 - x^2})$, the *direct* versions return the generator of
$\mathbb{Q}(\operatorname{asin} x)$, possibly after transforming the
argument; simple algebraic arguments give rational multiples of $\pi$.
The inverse cosine is computed as $\operatorname{acos} x = \pi/2 -
\operatorname{asin} x$. `ca_asin` and `ca_acos` use the logarithm forms
unless `CA_OPT_TRIG_FORM` is `CA_TRIG_DIRECT`.
:::
:::

## Special functions {#special-functions}

::: function
void ca_gamma(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to $\Gamma(x)$. Integers give factorials (or $\tilde \infty$
at the poles $0, -1, -2, \ldots$), half-integers give rational multiples
of $\sqrt{\pi}$, and otherwise the result is an element of a field
containing $\Gamma(x)$; the relation $\Gamma(x + n) = x (x+1) \cdots
(x+n-1) \Gamma(x)$ between generators that differ by integers is added
to the ideal of relations.

::: details
For special values: $\Gamma(+\infty) = +\infty$; $\Gamma(-\infty)$,
$\Gamma(\tilde \infty)$ and $\Gamma(\text{Undefined})$ are *Undefined*;
$\Gamma(\pm i \infty) = 0$; and the result is *Unknown* for *Unknown*
and for other signed infinities.
:::
:::

::: function
void ca_erf(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_erfc(ca_t res, const ca_t x, ca_ctx_t ctx)
void ca_erfi(ca_t res, const ca_t x, ca_ctx_t ctx)

Sets `res` to the error function $\operatorname{erf}(x)$, the
complementary error function $\operatorname{erfc}(x) = 1 -
\operatorname{erf}(x)$ or the imaginary error function
$\operatorname{erfi}(x) = -i \operatorname{erf}(ix)$. Odd symmetry is
used to normalise the argument, and the relations between the three
functions at related arguments are added to the ideal of relations.
For special values, $\operatorname{erf}(\pm \infty) = \pm 1$,
$\operatorname{erfc}(+\infty) = 0$, $\operatorname{erfc}(-\infty) = 2$,
$\operatorname{erf}(\pm i \infty) = \pm i \infty$, and unsigned infinity
gives *Undefined*.
:::

## Numerical evaluation {#numerical}

::: function
void ca_get_acb_raw(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)

Sets `res` to an enclosure of the value of `x`, computed by evaluating
the representation in ball arithmetic with working precision `prec`,
without any adaptive refinement. Enclosures of the extension numbers
are cached in the context. A special value gives an indeterminate ball
(see `acb_indeterminate`).
:::

::: function
void ca_get_acb(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
void ca_get_acb_accurate_parts(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)

Sets `res` to an enclosure of the value of `x`, increasing the working
precision until the result has a relative accuracy of `prec` bits, or
for the *accurate_parts* version until the real and the imaginary part
separately have `prec` accurate bits.

::: details
The precision starts at about $1.05 \cdot$ `prec` and is doubled up to
the larger of twice the initial precision and `CA_OPT_PREC_LIMIT`. At
some point an exact test is made, so that an exact zero (or a real or
imaginary number, for the *accurate_parts* version) gives an exact zero
(or zero imaginary or real part) instead of a ball that only contains
zero. Use `acb_rel_accuracy_bits` to check whether the requested
accuracy was reached. The output is not rounded to `prec` bits; call
`acb_set_round` for this.
:::
:::

::: function
char * ca_get_decimal_str(const ca_t x, slong digits, ulong flags, ca_ctx_t ctx)

Returns a decimal approximation of `x` with up to `digits` significant
digits, as computed by `ca_get_acb` (or `ca_get_acb_accurate_parts` if
`flags` is 1), correct to within one unit in the last digit. Fewer
digits are returned if the evaluation does not reach the requested
accuracy, and the string `"?"` if `x` is not a number or no finite
enclosure is obtained. The caller must free the string with
`flint_free`.
:::

## Factorisation {#factor}

::: type
ca_factor_struct
ca_factor_t

A number in factored form $b_1^{e_1} b_2^{e_2} \cdots b_n^{e_n}$, where
the bases $b_i$ and the exponents $e_i$ are `ca_t` numbers (the
exponents need not be integers). The structure has the fields `base`,
`exp` (arrays of `alloc` numbers), `length` and `alloc`.
:::

::: function
void ca_factor_init(ca_factor_t fac, ca_ctx_t ctx)
void ca_factor_clear(ca_factor_t fac, ca_ctx_t ctx)

Initialises `fac` to the empty product, which represents 1, or clears
it.
:::

::: function
void ca_factor_one(ca_factor_t fac, ca_ctx_t ctx)

Sets `fac` to the empty product.
:::

::: function
void ca_factor_print(const ca_factor_t fac, ca_ctx_t ctx)

Prints the factors of `fac` to standard output, one line `(b) ^ (e)` per
factor.
:::

::: function
void ca_factor_insert(ca_factor_t fac, const ca_t base, const ca_t exp, ca_ctx_t ctx)

Multiplies `fac` by $b^e$, where $b$ is `base` and $e$ is `exp`. If a
base with the same representation (see `ca_equal_repr`) is already
present, its exponent is increased by $e$; otherwise the factor is
appended.
:::

::: function
void ca_factor_get_ca(ca_t res, const ca_factor_t fac, ca_ctx_t ctx)

Sets `res` to the product $b_1^{e_1} \cdots b_n^{e_n}$, computed with
`ca_pow` and `ca_mul`.
:::

::: function
void ca_factor(ca_factor_t res, const ca_t x, ulong flags, ca_ctx_t ctx)

Sets `res` to a factorisation $x = b_1^{e_1} \cdots b_n^{e_n}$ of the
number `x`, which must not be a special value. `flags` combines one of
the polynomial options and one of the integer options below.

::: details
A rational number is factored over $\Z$ as selected by the integer
option. An element of a field with one algebraic generator is currently
returned as a single factor. An element $p/q$ of a generic field
$\mathbb{Q}(a_1, \ldots, a_n)$ is factored as a rational function of the
generators, that is, as a formal expression: the factors of $p$ and $q$
as polynomials in $a_1, \ldots, a_n$ (with negative exponents for $q$),
followed by the factors of the rational content. This is structural:
the bases are not checked to be irreducible as numbers.
:::
:::

::: macro
CA_FACTOR_POLY_NONE
CA_FACTOR_POLY_CONTENT
CA_FACTOR_POLY_SQF
CA_FACTOR_POLY_FULL

Polynomial factorisation in generic fields $\mathbb{Q}(a_1, \ldots,
a_n)$: none, extraction of the rational content only, squarefree
factorisation (and content), or complete factorisation into irreducible
multivariate polynomials (and content). Each option includes the
previous ones, so only one should be given.
:::

::: macro
CA_FACTOR_ZZ_NONE
CA_FACTOR_ZZ_SMOOTH
CA_FACTOR_ZZ_FULL

Factorisation of rational numbers and of the rational content: none,
extraction of small prime factors (of up to about `CA_OPT_SMOOTH_LIMIT`
bits) and perfect powers, or complete factorisation into primes, which
is prohibitively slow for general integers of more than 70 to 80 digits.
Only one should be given.
:::

## Field structure operations {#field-structure}

::: function
void ca_merge_fields(ca_t resx, ca_t resy, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `resx` and `resy` to copies of `x` and `y` represented as elements
of a common field, whose generators are the union of the generators of
both fields, sorted in the structural order used for elimination. Both
inputs must be numbers, not special values, and the outputs must be
distinct from the inputs (unless the fields are already equal or one of
them is $\mathbb{Q}$).

::: details
The common field and its ideal of relations are looked up in, or
inserted into, the field cache of the context. No generators are
eliminated even if they satisfy relations; the ideal of the merged field
records such relations instead.
:::
:::

::: function
void ca_condense_field(ca_t res, ca_ctx_t ctx)

Moves `res` to a smaller field when this is possible by inspecting its
representation: a rational value is demoted to $\mathbb{Q}$, and an
element of a generic field that only depends on one algebraic generator
becomes an element of the number field generated by it. This is done
automatically by most operations.
:::

::: function
ca_ext_ptr ca_is_gen_as_ext(const ca_t x, ca_ctx_t ctx)

If `x` is a generator $a_k$ of its field $\mathbb{Q}(a_1, \ldots, a_n)$,
returns a pointer to the extension number $a_k$; otherwise returns
`NULL`.
:::

## Internal functions {#internal}

::: macro
CA_FMPQ(x)
CA_FMPQ_NUMREF(x)
CA_FMPQ_DENREF(x)

Assuming that `x` is represented as an element of $\mathbb{Q}$, these
macros give the data as an `fmpq_t` and its numerator and denominator
as `fmpz_t`.
:::

::: macro
CA_MPOLY_Q(x)

Assuming that `x` is an element of a generic field, gives its data as an
`fmpz_mpoly_q_t`.
:::

::: macro
CA_NF_ELEM(x)

Assuming that `x` is an element of a field with one algebraic generator,
gives its data as an `nf_elem_t`.
:::

::: function
void _ca_make_field_element(ca_t x, ca_field_srcptr field, ca_ctx_t ctx)

Changes the internal representation of `x` to that of an element of
`field` (allocating or freeing data as needed). The value of `x` is
destroyed.
:::

::: function
void _ca_make_fmpq(ca_t x, ca_ctx_t ctx)

Changes the internal representation of `x` to that of an element of
$\mathbb{Q}$. The value of `x` is destroyed unless it already is one.
:::

::: function
ca_field_ptr _ca_ctx_get_field_const(ca_ctx_t ctx, calcium_func_code func)
ca_field_ptr _ca_ctx_get_field_fx(ca_ctx_t ctx, calcium_func_code func, const ca_t x)
ca_field_ptr _ca_ctx_get_field_fxy(ca_ctx_t ctx, calcium_func_code func, const ca_t x, const ca_t y)

Returns the field $\mathbb{Q}(a)$ generated by the extension number $a$
= `func`, `func(x)` or `func(x, y)`, creating and caching the extension
number and the field if necessary.
:::

::: function
void _ca_function_fx(ca_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
void _ca_function_fxy(ca_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)

Sets `res` to the generator of the field $\mathbb{Q}(f(x))$ or
$\mathbb{Q}(f(x, y))$ with $f$ = `func`, without any simplification or
checks. The arguments must be numbers.
:::

::: function
void _ca_mpoly_q_reduce_ideal(fmpz_mpoly_q_t res, ca_field_srcptr field, ca_ctx_t ctx)

Reduces the numerator and the denominator of `res`, an element of the
generic field `field`, modulo the ideal of relations of the field (by
pseudo-division with remainder), and cancels common factors.
:::

::: function
void _ca_mpoly_q_simplify_fraction_ideal(fmpz_mpoly_q_t res, ca_field_srcptr field, ca_ctx_t ctx)

Tries to simplify the fraction `res` using the ideal of relations of
`field`, currently by moving monomial factors of the denominator to the
numerator: if a relation has a leading monomial $u$ sharing a factor $g$
with the denominator, numerator and denominator are multiplied by $u/g$
and reduced, which is kept if the denominator decreases.
:::

::: macro
CA_TEST_PROPERTY(f, s, x, ctx, expected)

Used in the test suite: evaluates the predicate `f(x, ctx)` and aborts
with a message containing the string `s` if the result differs from the
`truth_t` value `expected`.
:::
