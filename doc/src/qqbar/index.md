---
title: qqbar.h
header: qqbar.h
description: Algebraic numbers represented canonically by their minimal polynomial and an isolating ball, with exact arithmetic, comparisons, roots of polynomials, trigonometric values and conversion to closed forms.
category: Algebraic numbers and number fields
order: 4
---

A `qqbar_t` represents a real or complex algebraic number $x \in
\overline{\mathbb{Q}}$ exactly, by its minimal polynomial $P \in
\Z[t]$ (irreducible, primitive, with positive leading coefficient) and
a complex ball that contains $x$ and no other root of $P$. This
representation is canonical: two `qqbar_t` values are equal if and only
if they have the same minimal polynomial and their balls overlap, so
equality, signs, comparisons and all predicates are decided exactly.
The balls are refined automatically when an operation needs more
precision.

The representation suits individual algebraic numbers of moderate
degree, up to 100 say. Arithmetic is expensive: in the general case,
$x \pm y$, $xy$ and $x/y$ for numbers of degrees $m$ and $n$ require an
annihilating polynomial of degree $mn$ that is computed from power sums
and then factored over $\Z$, followed by selecting the right factor
numerically. FLINT avoids this in many structured cases (quadratic
operands, conjugate operands, operands in a common number field or
cyclotomic field), see [](algorithms/exact/qqbar-arithmetic). For many
operations in one fixed number field it is much more efficient to
compute with `ca_t` (see [](ca)) or directly with [](nf_elem), and to
convert to `qqbar_t` only to decide predicates.

`qqbar_t` is used by [](ca) to represent the algebraic generators of
number fields and to decide predicates for algebraic numbers. The
generic rings interface provides the fields of real and complex
algebraic numbers with `qqbar_t` elements through
`gr_ctx_init_real_qqbar` and `gr_ctx_init_complex_qqbar`, with optional
limits on degrees and sizes set by `_gr_ctx_qqbar_set_limits` (see
[](gr_domains)); generic polynomial and matrix functions such as
`gr_poly_roots` and `gr_mat_eigenvalues` then work over
$\overline{\mathbb{Q}}$. The guide [](qqbar/using) introduces the
module with examples.

::: summary algorithms/exact/qqbar-representation :::

## Types and macros {#types}

::: type
qqbar_struct
qqbar_t

An algebraic number. The structure consists of an `fmpz_poly_struct`
holding the minimal polynomial and an `acb_struct` holding the
isolating ball (the *enclosure*). A `qqbar_t` is an array of one
structure, so it is passed by reference.
:::

::: type
qqbar_ptr
qqbar_srcptr

Aliases for `qqbar_struct *` and `const qqbar_struct *`, used for arrays
of algebraic numbers.
:::

::: macro
QQBAR_POLY(x)
QQBAR_COEFFS(x)
QQBAR_ENCLOSURE(x)

Give the minimal polynomial of `x` (usable as an `fmpz_poly_t`), its
coefficient array (an `fmpz *`) and the enclosure (usable as an
`acb_t`). The data must not be modified in a way that breaks the
invariants: the polynomial must remain the primitive minimal polynomial
with positive leading coefficient, and the ball must isolate the root.
:::

::: macro
QQBAR_DEFAULT_PREC

The working precision (128 bits) at which isolating balls are computed
by default and at which refinement loops start.
:::

## Memory management {#memory}

::: function
void qqbar_init(qqbar_t res)

Initialises `res` and sets it to zero.
:::

::: function
void qqbar_clear(qqbar_t res)

Clears `res`, freeing its memory.
:::

::: function
qqbar_ptr _qqbar_vec_init(slong len)
void _qqbar_vec_clear(qqbar_ptr vec, slong len)

Allocates an array of `len` algebraic numbers initialised to zero, or
clears the `len` entries of `vec` and frees the array.
:::

## Assignment {#assignment}

::: function
void qqbar_swap(qqbar_t x, qqbar_t y)

Swaps `x` and `y` in constant time.
:::

::: function
void qqbar_set(qqbar_t res, const qqbar_t x)
void qqbar_set_si(qqbar_t res, slong x)
void qqbar_set_ui(qqbar_t res, ulong x)
void qqbar_set_fmpz(qqbar_t res, const fmpz_t x)
void qqbar_set_fmpq(qqbar_t res, const fmpq_t x)

Sets `res` to `x`. A rational number $p/q$ in lowest terms with $q > 0$
has the minimal polynomial $q t - p$ and an exact enclosure.
:::

::: function
void qqbar_set_re_im(qqbar_t res, const qqbar_t x, const qqbar_t y)

Sets `res` to $x + yi$ for real `x` and `y` (not checked).
:::

::: function
int qqbar_set_d(qqbar_t res, double x)
int qqbar_set_re_im_d(qqbar_t res, double x, double y)

Sets `res` to the exact value of the double `x`, or of $x + yi$, and
returns 1. If an input is infinite or NaN, returns 0 and leaves `res`
unchanged.
:::

## Properties {#properties}

::: function
slong qqbar_degree(const qqbar_t x)

Returns the degree of `x`, that is, of its minimal polynomial.
:::

::: function
int qqbar_is_rational(const qqbar_t x)
int qqbar_is_integer(const qqbar_t x)
int qqbar_is_algebraic_integer(const qqbar_t x)

Returns whether `x` is a rational number (degree 1), an integer, or an
algebraic integer (its minimal polynomial is monic).
:::

::: function
int qqbar_is_zero(const qqbar_t x)
int qqbar_is_one(const qqbar_t x)
int qqbar_is_neg_one(const qqbar_t x)
int qqbar_is_i(const qqbar_t x)
int qqbar_is_neg_i(const qqbar_t x)

Returns whether `x` is $0$, $1$, $-1$, $i$ or $-i$. These inspect the
minimal polynomial (and the sign of the imaginary part of the enclosure
for $\pm i$) and are therefore very fast.
:::

::: function
int qqbar_is_real(const qqbar_t x)

Returns whether `x` is real. This is decided exactly, see
`qqbar_sgn_im`.
:::

::: function
void qqbar_height(fmpz_t res, const qqbar_t x)
slong qqbar_height_bits(const qqbar_t x)

Sets `res` to the height of `x`, the largest absolute value of a
coefficient of its minimal polynomial, or returns the number of bits of
the height.
:::

::: function
int qqbar_within_limits(const qqbar_t x, slong deg_limit, slong bits_limit)

Returns whether `x` has degree at most `deg_limit` and height of at most
`bits_limit` bits. A limit of 0 means no limit.
:::

::: function
int qqbar_binop_within_limits(const qqbar_t x, const qqbar_t y, slong deg_limit, slong bits_limit)

Returns whether the trivial bounds for the result of a binary operation
on `x` and `y` are within the limits: the product of the degrees must be
at most `deg_limit` and the sum of the height bits at most
`bits_limit`. The degree bound is rigorous, the height bound only a
heuristic. A limit of 0 means no limit.
:::

## Conversions {#conversions}

::: function
void _qqbar_get_fmpq(fmpz_t num, fmpz_t den, const qqbar_t x)
void qqbar_get_fmpq(fmpq_t res, const qqbar_t x)
void qqbar_get_fmpz(fmpz_t res, const qqbar_t x)

Sets the output to the rational number `x` (as a numerator and
denominator for the underscore version), or to the integer `x`. Throws
an exception if `x` is not rational, respectively not an integer.
:::

## Special values {#special}

::: function
void qqbar_zero(qqbar_t res)
void qqbar_one(qqbar_t res)
void qqbar_i(qqbar_t res)
void qqbar_phi(qqbar_t res)

Sets `res` to 0, 1, the imaginary unit $i$, or the golden ratio $\varphi
= \frac{1}{2}(1 + \sqrt{5})$.
:::

## Input and output {#io}

::: function
void qqbar_print(const qqbar_t x)

Prints `x` to standard output in a form intended for debugging: the
degree, the coefficients of the minimal polynomial (from the constant
term upwards) and the enclosure, for example `deg 2 [-1, -1, 1]
[1.61803398874989484820458683436563811772 +/- 6.00e-39]` for $\varphi$.
:::

::: function
void qqbar_printn(const qqbar_t x, slong n)
void qqbar_printnd(const qqbar_t x, slong n)

Prints a decimal approximation of `x` with `n` significant digits to
standard output, such as `1.41421 - 1.73205*I`; `qqbar_printnd` appends
the degree, as in `1.41421 - 1.73205*I (deg 4)`.
:::

::: function
char * qqbar_get_str_nd(const qqbar_t x, slong n)

Returns a string with the output of `qqbar_printnd`. The caller must
free it with `flint_free`.
:::

## Random generation {#random}

::: function
void qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits)
void qqbar_randtest_real(qqbar_t res, flint_rand_t state, slong deg, slong bits)
void qqbar_randtest_nonreal(qqbar_t res, flint_rand_t state, slong deg, slong bits)

Sets `res` to a random algebraic number, a random real algebraic number,
or a random nonreal algebraic number, with degree at most `deg` and
height of at most `bits` bits. For `qqbar_randtest_nonreal`, `deg` must
be at least 2.
:::

## Comparisons {#comparisons}

All comparisons are exact. They first compare enclosures, refining them
as needed, and fall back to exact computations (such as computing $x -
y$) only when the numbers are equal or very close.

::: function
int qqbar_equal(const qqbar_t x, const qqbar_t y)

Returns whether $x = y$: the minimal polynomials must be equal and the
enclosures must contain a common root, which is decided by refining
them until they are disjoint or their union is proved to contain a
single root.
:::

::: function
int qqbar_equal_fmpq_poly_val(const qqbar_t x, const fmpq_poly_t f, const qqbar_t y)

Returns whether $x = f(y)$. This is faster than computing $f(y)$ with
`qqbar_evaluate_fmpq_poly`: the value $f(y)$ is first compared
numerically with $x$, and only if they cannot be separated is $P(f(y))
= 0$ checked exactly in $\mathbb{Q}(y)$, where $P$ is the minimal
polynomial of $x$.
:::

::: function
int qqbar_cmp_re(const qqbar_t x, const qqbar_t y)
int qqbar_cmp_im(const qqbar_t x, const qqbar_t y)

Compares the real parts, or the imaginary parts, of `x` and `y`,
returning $-1$, 0 or 1.
:::

::: function
int qqbar_cmpabs_re(const qqbar_t x, const qqbar_t y)
int qqbar_cmpabs_im(const qqbar_t x, const qqbar_t y)
int qqbar_cmpabs(const qqbar_t x, const qqbar_t y)

Compares the absolute values of the real parts, of the imaginary parts,
or of `x` and `y` themselves, returning $-1$, 0 or 1.
:::

::: function
int qqbar_cmp_root_order(const qqbar_t x, const qqbar_t y)

Compares `x` and `y` in the order used to sort roots of polynomials,
returning $-1$, 0 or 1: real numbers come first, in decreasing order;
nonreal numbers follow, by decreasing real part, then by increasing
absolute value of the imaginary part, then with positive imaginary part
first. Complex conjugates are therefore adjacent, the one in the upper
half-plane first.
:::

::: function
ulong qqbar_hash(const qqbar_t x)

Returns a hash of `x`. Currently only the minimal polynomial is hashed,
so conjugates have equal hashes, and similar polynomials give similar
hash values; apply an integer hash function to the result if all bits
must look random.
:::

## Complex parts {#complex-parts}

::: function
void qqbar_conj(qqbar_t res, const qqbar_t x)

Sets `res` to the complex conjugate of `x` (the minimal polynomial is
unchanged).
:::

::: function
void qqbar_re(qqbar_t res, const qqbar_t x)
void qqbar_im(qqbar_t res, const qqbar_t x)
void qqbar_re_im(qqbar_t res1, qqbar_t res2, const qqbar_t x)

Sets `res` to the real part or the imaginary part of `x`, or `res1` and
`res2` to both, computed as $(x + \overline{x})/2$ and $(x -
\overline{x})/(2i)$.
:::

::: function
void qqbar_abs(qqbar_t res, const qqbar_t x)
void qqbar_abs2(qqbar_t res, const qqbar_t x)

Sets `res` to $|x|$ or $|x|^2 = x \overline{x}$.
:::

::: function
void qqbar_sgn(qqbar_t res, const qqbar_t x)

Sets `res` to the complex sign of `x`: 0 if $x = 0$, otherwise $x/|x|$.
:::

::: function
int qqbar_sgn_re(const qqbar_t x)
int qqbar_sgn_im(const qqbar_t x)

Returns the sign ($-1$, 0 or 1) of the real or imaginary part of `x`.

::: details
A sign is read off the enclosure if it excludes zero; otherwise the
enclosure is refined. Zero is certified by proving that the ball with
the real (or imaginary) part set to zero contains a root of the minimal
polynomial, which by uniqueness must be $x$. A real part can only be
zero if the minimal polynomial is even, which is checked first. See
[](algorithms/exact/qqbar-representation).
:::
:::

::: function
int qqbar_csgn(const qqbar_t x)

Returns the extension of the real sign function that is 1 in the open
right half-plane, $-1$ in the open left half-plane, the sign of the
imaginary part on the imaginary axis, and 0 at 0; equivalently
$\operatorname{csgn}(x) = x / \sqrt{x^2}$ for $x \ne 0$.
:::

## Integer parts {#integer-parts}

::: function
void qqbar_floor(fmpz_t res, const qqbar_t x)
void qqbar_ceil(fmpz_t res, const qqbar_t x)

Sets `res` to $\lfloor x \rfloor$ or $\lceil x \rceil$; for nonreal `x`,
of its real part.
:::

## Arithmetic {#arithmetic}

Arithmetic with a rational operand is an affine transformation of the
minimal polynomial (see `qqbar_scalar_op`) and is cheap. Arithmetic
with two irrational operands is described in
[](algorithms/exact/qqbar-arithmetic).

::: function
void qqbar_neg(qqbar_t res, const qqbar_t x)

Sets `res` to $-x$.
:::

::: function
void qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y)
void qqbar_add_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
void qqbar_add_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
void qqbar_add_ui(qqbar_t res, const qqbar_t x, ulong y)
void qqbar_add_si(qqbar_t res, const qqbar_t x, slong y)

Sets `res` to $x + y$.
:::

::: function
void qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y)
void qqbar_sub_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
void qqbar_sub_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
void qqbar_sub_ui(qqbar_t res, const qqbar_t x, ulong y)
void qqbar_sub_si(qqbar_t res, const qqbar_t x, slong y)
void qqbar_fmpq_sub(qqbar_t res, const fmpq_t x, const qqbar_t y)
void qqbar_fmpz_sub(qqbar_t res, const fmpz_t x, const qqbar_t y)
void qqbar_ui_sub(qqbar_t res, ulong x, const qqbar_t y)
void qqbar_si_sub(qqbar_t res, slong x, const qqbar_t y)

Sets `res` to $x - y$.
:::

::: function
void qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y)
void qqbar_mul_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
void qqbar_mul_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
void qqbar_mul_ui(qqbar_t res, const qqbar_t x, ulong y)
void qqbar_mul_si(qqbar_t res, const qqbar_t x, slong y)

Sets `res` to $x y$.
:::

::: function
void qqbar_mul_2exp_si(qqbar_t res, const qqbar_t x, slong e)

Sets `res` to $x \cdot 2^e$.
:::

::: function
void qqbar_sqr(qqbar_t res, const qqbar_t x)

Sets `res` to $x^2$.
:::

::: function
void qqbar_inv(qqbar_t res, const qqbar_t x)

Sets `res` to $1/x$, whose minimal polynomial is the reversal of that of
`x`. Throws an exception (`FLINT_DIVZERO`) if $x = 0$.
:::

::: function
void qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y)
void qqbar_div_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
void qqbar_div_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
void qqbar_div_ui(qqbar_t res, const qqbar_t x, ulong y)
void qqbar_div_si(qqbar_t res, const qqbar_t x, slong y)
void qqbar_fmpq_div(qqbar_t res, const fmpq_t x, const qqbar_t y)
void qqbar_fmpz_div(qqbar_t res, const fmpz_t x, const qqbar_t y)
void qqbar_ui_div(qqbar_t res, ulong x, const qqbar_t y)
void qqbar_si_div(qqbar_t res, slong x, const qqbar_t y)

Sets `res` to $x / y$. Throws an exception (`FLINT_DIVZERO`) if $y = 0$.
:::

::: function
void qqbar_scalar_op(qqbar_t res, const qqbar_t x, const fmpz_t a, const fmpz_t b, const fmpz_t c)

Sets `res` to $(a x + b)/c$ in one operation, for integers $a$, $b$ and
$c \ne 0$, by transforming the minimal polynomial and the enclosure.
Throws an exception if $c = 0$.
:::

## Powers and roots {#powers}

::: function
void qqbar_sqrt(qqbar_t res, const qqbar_t x)
void qqbar_sqrt_ui(qqbar_t res, ulong x)

Sets `res` to the principal square root of `x`.
:::

::: function
void qqbar_rsqrt(qqbar_t res, const qqbar_t x)

Sets `res` to $1/\sqrt{x}$. Throws an exception if $x = 0$.
:::

::: function
void qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong n)
void qqbar_pow_si(qqbar_t res, const qqbar_t x, slong n)
void qqbar_pow_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t n)
void qqbar_pow_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t n)

Sets `res` to $x^n$; for rational $n = p/q$ the principal branch
$(x^{1/q})^p$. $x^0 = 1$ for all $x$. Raising zero to a negative power
throws an exception, and so does a rational exponent whose numerator or
denominator does not fit in a word.
:::

::: function
void qqbar_root_ui(qqbar_t res, const qqbar_t x, ulong n)
void qqbar_fmpq_root_ui(qqbar_t res, const fmpq_t x, ulong n)

Sets `res` to the principal $n$-th root of `x`, for $n \ge 1$. The
minimal polynomial is found among the factors of $P(t^n)$, where $P$ is
the minimal polynomial of $x$; principal roots of positive rational
numbers and roots of roots of unity are computed directly.

::: warning
For $n = 0$, or if $n \deg(x)$ is extremely large (more than about
$2^{30}$), `qqbar_root_ui` prints a message to standard output and
returns without setting `res`, instead of throwing an exception.
:::
:::

::: function
void qqbar_fmpq_pow_si_ui(qqbar_t res, const fmpq_t x, slong m, ulong n)

Sets `res` to the principal branch of $x^{m/n}$ for a rational `x` and
$n \ge 1$. Throws an exception if $x = 0$ and $m < 0$.
:::

::: function
int qqbar_pow(qqbar_t res, const qqbar_t x, const qqbar_t y)

General exponentiation: if $x^y$ is an algebraic number, sets `res` to
it and returns 1; if $x^y$ is transcendental or undefined, returns 0. By
the Gelfond–Schneider theorem, $x^y$ for algebraic $x \notin \{0, 1\}$
and irrational algebraic $y$ is transcendental, so the result is
algebraic only if $y$ is rational, $x \in \{0, 1\}$ or $y = 0$. Unlike
the other functions, this one returns 0 instead of throwing for $0^y$
with $\operatorname{Re}(y) \le 0$, $y \ne 0$.
:::

## Numerical enclosures {#enclosures}

The following functions return polished enclosures, in which the real
and imaginary parts are both accurate to `prec` bits and are exact when
they are exactly representable (that is, when a real or imaginary part
is a sufficiently small dyadic number). Polishing can be expensive; use
`qqbar_enclosure_raw` when it is not needed, or `qqbar_cache_enclosure`
to avoid repeated work.

::: function
void qqbar_get_acb(acb_t res, const qqbar_t x, slong prec)
void qqbar_get_arb(arb_t res, const qqbar_t x, slong prec)

Sets `res` to an enclosure of `x` rounded to `prec` bits. For
`qqbar_get_arb`, if `x` is not real, `res` is set to $[\mathrm{NaN} \pm
\infty]$.
:::

::: function
void qqbar_get_arb_re(arb_t res, const qqbar_t x, slong prec)
void qqbar_get_arb_im(arb_t res, const qqbar_t x, slong prec)

Sets `res` to an enclosure of the real or imaginary part of `x` rounded
to `prec` bits.
:::

::: function
void qqbar_cache_enclosure(qqbar_t res, slong prec)

Refines the enclosure stored in `res` to somewhat more than `prec` bits
(at least `QQBAR_DEFAULT_PREC`). Operations normally refine enclosures
in temporary variables without storing the result; calling this
function first helps when `res` is used as an operand many times.
:::

## Numerator and denominator {#numerator}

::: function
void qqbar_denominator(fmpz_t res, const qqbar_t y)

Sets `res` to the leading coefficient $d$ of the minimal polynomial of
`y`, a denominator in the sense that $d y$ is an algebraic integer.
:::

::: function
void qqbar_numerator(qqbar_t res, const qqbar_t y)

Sets `res` to $d y$, where $d$ is given by `qqbar_denominator`.
:::

## Conjugates {#conjugates}

::: function
void qqbar_conjugates(qqbar_ptr res, const qqbar_t x)

Sets the entries of the array `res` to the $d$ conjugates of `x`
(including `x`), where $d$ is the degree of `x`, sorted as by
`qqbar_cmp_root_order`.
:::

## Polynomial evaluation {#evaluation}

::: function
void _qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpz * poly, const fmpz_t den, slong len, const qqbar_t x)
void qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpq_poly_t poly, const qqbar_t x)
void _qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz * poly, slong len, const qqbar_t x)
void qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz_poly_t poly, const qqbar_t x)

Sets `res` to $f(x)$ for the polynomial $f$ given by `poly` (with
denominator `den` for the underscore *fmpq* version). The polynomial is
first reduced modulo the minimal polynomial of `x` if its degree is too
large; constant and linear polynomials and rational `x` are handled
directly. In general, the minimal polynomial of $f(x)$ is the minimal
polynomial of the matrix of multiplication by $f(x)$ in
$\mathbb{Q}(x)$, see [](algorithms/exact/qqbar-arithmetic).
:::

::: function
int qqbar_evaluate_fmpz_mpoly(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the value of the multivariate polynomial `poly` at the
algebraic numbers `x` (one per variable of `ctx`) and returns 1. If a
temporary result (not necessarily the final value, in case of
cancellation) exceeds `deg_limit` or `bits_limit` in the sense of
`qqbar_binop_within_limits`, the evaluation is aborted and 0 is
returned. The evaluation uses the generic multivariate Horner scheme
`gr_fmpz_mpoly_evaluate_horner` over the $\overline{\mathbb{Q}}$
context (`gr_fmpz_mpoly_evaluate_iter` for a single term).
:::

## Polynomial roots {#roots}

::: macro
QQBAR_ROOTS_IRREDUCIBLE
QQBAR_ROOTS_UNSORTED

Flags for the root-finding functions: the polynomial is known to be
irreducible (it may still have content), so that no factorisation is
done; the roots need not be sorted.
:::

::: function
void qqbar_roots_fmpz_poly(qqbar_ptr res, const fmpz_poly_t poly, int flags)
void qqbar_roots_fmpq_poly(qqbar_ptr res, const fmpq_poly_t poly, int flags)

Sets the entries of the array `res` to the $d = \deg(\mathit{poly})$
roots of `poly`, repeated according to multiplicity. Unless
`QQBAR_ROOTS_UNSORTED` is given, the roots are sorted by
`qqbar_cmp_root_order`; copies of a repeated root are always
consecutive. See [](qqbar/roots).
:::

::: function
slong qqbar_real_roots_fmpz_poly(qqbar_ptr res, const fmpz_poly_t poly, int flags)
slong qqbar_real_roots_fmpq_poly(qqbar_ptr res, const fmpq_poly_t poly, int flags)

Sets the first entries of `res` to the real roots of `poly`, repeated
according to multiplicity, and returns their number. The array must
have room for $\deg(\mathit{poly})$ entries. Unless
`QQBAR_ROOTS_UNSORTED` is given, the roots are in decreasing order, so
that the output agrees with the first entries of the output of
`qqbar_roots_fmpz_poly`. This is faster than computing all roots, since
nonreal roots are never isolated and factors without real roots are
skipped.
:::

::: function
void qqbar_eigenvalues_fmpz_mat(qqbar_ptr res, const fmpz_mat_t mat, int flags)
void qqbar_eigenvalues_fmpq_mat(qqbar_ptr res, const fmpq_mat_t mat, int flags)

Sets the entries of `res` to the eigenvalues of the square matrix
`mat`, repeated according to algebraic multiplicity, by computing the
characteristic polynomial and calling `qqbar_roots_fmpz_poly` with the
same flags.
:::

::: function
int _qqbar_roots_poly_squarefree(qqbar_ptr roots, qqbar_srcptr coeffs, slong len, slong deg_limit, slong bits_limit)

Sets `roots` to the $d = \mathit{len} - 1$ roots of the squarefree
polynomial $f = a_0 + \cdots + a_d x^d$ with algebraic coefficients
`coeffs` ($a_d \ne 0$), and returns 1; returns 0 if a limit is
exceeded.

::: details
An annihilating polynomial with rational coefficients is the product
$g = \prod (\tilde a_0 + \tilde a_1 x + \cdots + \tilde a_d x^d)$ over
all combinations of conjugates $\tilde a_k$ of the coefficients. It is
computed numerically, as a product of polynomials with ball
coefficients by binary splitting, and multiplied by a power of the
common denominator of the coefficients so that its coefficients are
integers that can be read off once the precision is high enough. Then
$g$ is factored and the factors are matched against the numerical roots
of $f$. Since $g$ can be huge, the computation is aborted when its
degree (the number of conjugate combinations) or size exceeds
`deg_limit` or `bits_limit` (nonpositive values remove the limits, but the function can
still fail if $g$ exceeds machine size). For polynomials that may have
repeated roots, use `gr_poly_roots` with a `qqbar` context, which
performs a squarefree factorisation first; `gr_poly_roots_other`
computes algebraic roots of polynomials over other rings.
:::
:::

## Roots of unity and trigonometric functions {#trig}

The following functions take rational numbers as two words `p`, `q`
rather than as `fmpq_t`, since the minimal polynomials have degree
$O(q)$ and only small $q$ are feasible. The inputs need not be reduced
but should not be close to the limits of a word, as they may be added
or multiplied internally.

::: function
void qqbar_root_of_unity(qqbar_t res, slong p, ulong q)
void qqbar_exp_pi_i(qqbar_t res, slong p, ulong q)

Sets `res` to the root of unity $e^{2 \pi i p / q}$ or $e^{\pi i p /
q}$. The minimal polynomial is a cyclotomic polynomial.
:::

::: function
int qqbar_is_root_of_unity(slong * p, ulong * q, const qqbar_t x)

Returns whether `x` is a root of unity, that is, whether its minimal
polynomial is cyclotomic. If so and if `p` and `q` are not `NULL`, sets
them to the integers with $0 \le p < q$ and $\gcd(p, q) = 1$ such that $x
= e^{2 \pi i p / q}$.
:::

::: function
void qqbar_cos_pi(qqbar_t res, slong p, ulong q)
void qqbar_sin_pi(qqbar_t res, slong p, ulong q)
int qqbar_tan_pi(qqbar_t res, slong p, ulong q)
int qqbar_cot_pi(qqbar_t res, slong p, ulong q)
int qqbar_sec_pi(qqbar_t res, slong p, ulong q)
int qqbar_csc_pi(qqbar_t res, slong p, ulong q)

Sets `res` to $\cos(\pi x)$, $\sin(\pi x)$, $\tan(\pi x)$, $\cot(\pi x)$,
$\sec(\pi x)$ or $\csc(\pi x)$ for $x = p/q$. The functions with an
`int` return value return 1, or 0 if $x$ is a pole. The cosine uses the
minimal polynomial of $\cos(2 \pi / n)$ (`fmpz_poly_cos_minpoly`), and
the other functions are reduced to it.
:::

::: function
int qqbar_log_pi_i(slong * p, ulong * q, const qqbar_t x)

If $y = \log(x) / (\pi i)$ is algebraic, which happens exactly when $x$
is a root of unity and then $y$ is rational, sets $p/q$ to $y$ in lowest
terms with $-1 < y \le 1$ and returns 1; otherwise returns 0.
:::

::: function
int qqbar_atan_pi(slong * p, ulong * q, const qqbar_t x)
int qqbar_asin_pi(slong * p, ulong * q, const qqbar_t x)
int qqbar_acos_pi(slong * p, ulong * q, const qqbar_t x)
int qqbar_acot_pi(slong * p, ulong * q, const qqbar_t x)
int qqbar_asec_pi(slong * p, ulong * q, const qqbar_t x)
int qqbar_acsc_pi(slong * p, ulong * q, const qqbar_t x)

If $y = f(x) / \pi$ is algebraic, where $f$ is the principal branch of
the inverse trigonometric function, sets $p/q$ to $y$ (necessarily
rational) in lowest terms and returns 1; otherwise returns 0. The ranges
are $|y| < \frac{1}{2}$ for atan, $|y| \le \frac{1}{2}$ for asin and
acsc, $0 \le y \le 1$ for acos and asec, and $-\frac{1}{2} < y \le
\frac{1}{2}$ for acot.
:::

## Guessing and simplification {#guessing}

::: function
int qqbar_guess(qqbar_t res, const acb_t z, slong max_deg, slong max_bits, int flags, slong prec)

Tries to find an algebraic number of degree at most `max_deg` and height
of at most `max_bits` bits contained in the ball `z`, using LLL to find
an integer relation between $1, z, \ldots, z^{\mathit{max\_deg}}$ at
working precision `prec`. Returns 1 and sets `res` on success. The
result is only a heuristic: success does not prove that the number
behind `z` equals `res`, and failure does not prove that no such
algebraic number exists.

::: details
The precision should normally be the precision at which `z` was
computed, and of the order of `max_deg` times `max_bits` at least. If
`flags` is 1, a search with `max_deg`/4 is tried first when `max_deg`
exceeds 8. A single search is done; when the size of the answer is
unknown, call the function repeatedly with increasing parameters.
:::
:::

::: function
int qqbar_express_in_field(fmpq_poly_t res, const qqbar_t alpha, const qqbar_t x, slong max_bits, int flags, slong prec)

Tries to express `x` in the number field $\mathbb{Q}(\alpha)$: on
success, returns 1 and sets `res` to a polynomial $f$ of degree less
than $\deg \alpha$ with $f(\alpha) = x$. The identity is verified
exactly, so success is proved. Failure does not prove that $x \notin
\mathbb{Q}(\alpha)$.

::: details
A relation between $1, \alpha, \ldots, \alpha^{d-1}$ and $x$ is searched
with LLL at working precision `prec`; if $x \in \mathbb{Q}(\alpha)$, a
large enough precision always finds it. A rational `x` is handled
directly. The parameters `max_bits` and `flags` are currently ignored.
A single search is done; call the function repeatedly with increasing
`prec` if needed.
:::
:::

## Symbolic expressions and radicals {#fexpr}

The functions in this section use `fexpr_t` (see [](fexpr)) and are
declared only if `fexpr.h` is included before `qqbar.h`. See
[](qqbar/formulas) for examples.

::: function
void qqbar_get_quadratic(fmpz_t res_a, fmpz_t res_b, fmpz_t res_c, fmpz_t res_q, const qqbar_t x, int factoring)

For `x` of degree 1 or 2, computes integers $a$, $b$, $c$, $q$ with $x =
(a + b \sqrt{c})/q$, $q > 0$, $\gcd(a, b, q) = 1$ and $c$ not a perfect
square; for rational `x`, $b = c = 0$.

::: details
For a Gaussian rational, $c = -1$. Otherwise $c$ is derived from the
discriminant according to `factoring`: 0 removes only factors of two,
1 factors completely so that $c$ is squarefree (expensive for large
discriminants), and 2 removes small square factors found by smooth
factorisation, which gives nice output in most cases without the risk
of a hard factorisation.
:::
:::

::: function
int qqbar_set_fexpr(qqbar_t res, const fexpr_t expr)

Sets `res` to the algebraic number given by the expression `expr` and
returns 1, or returns 0 if `expr` is not understood.

::: details
The expression is evaluated statically with `qqbar` arithmetic; no
symbolic simplification is done. Supported are integers, `Decimal`
literals, arithmetic with algebraic operands, `Sqrt`, `Pow` with a
rational exponent, the constants `NumberI` and `GoldenRatio`,
`RootOfUnity`, `Floor`, `Ceil`, `Abs`, `Sign`, `Csgn`, `Conjugate`, `Re`,
`Im`, `Max`, `Min`, trigonometric functions of rational multiples of
`Pi`, `Exp` of rational multiples of `Pi * NumberI`, and the
constructors `AlgebraicNumberSerialized` (whose data is trusted),
`PolynomialRootIndexed` and `PolynomialRootNearest`. This covers the
output of all `qqbar_get_fexpr` functions. Expressions such as `Sub(Pi,
Pi)`, `Div(1, Infinity)` or `Sum(n, For(n, 1, 10))` are not supported,
although their values are algebraic.
:::
:::

::: function
void qqbar_get_fexpr_repr(fexpr_t res, const qqbar_t x)

Sets `res` to `AlgebraicNumberSerialized(List(coeffs), enclosure)`,
which reflects the internal representation exactly. It can be read back
with `qqbar_set_fexpr` cheaply and is the recommended format for
serialisation, although it is not meant for human readers.
:::

::: function
void qqbar_get_fexpr_root_nearest(fexpr_t res, const qqbar_t x)
void qqbar_get_fexpr_root_indexed(fexpr_t res, const qqbar_t x)

Sets `res` to `PolynomialRootNearest(List(coeffs), point)`, where the
decimal `point` is closer to `x` than to any other root, or to
`PolynomialRootIndexed(List(coeffs), index)`, where `index` (starting
from 1) is the position of `x` among the roots in the order of
`qqbar_cmp_root_order`. Both describe `x` unambiguously and readably,
and can be read back with `qqbar_set_fexpr`, at the cost of root
isolation.
:::

::: function
int qqbar_get_fexpr_formula(fexpr_t res, const qqbar_t x, ulong flags)

Tries to write `x` as a closed-form expression with arithmetic
operations, radicals and possibly exponentials or trigonometric
functions, without `PolynomialRootNearest` or `PolynomialRootIndexed`.
Returns 1 on success and 0 on failure. `flags` selects the methods (with
`flags` = 0 only rational numbers are handled).
:::

::: macro
QQBAR_FORMULA_GAUSSIANS
QQBAR_FORMULA_QUADRATICS
QQBAR_FORMULA_CYCLOTOMICS
QQBAR_FORMULA_CUBICS
QQBAR_FORMULA_QUARTICS
QQBAR_FORMULA_QUINTICS
QQBAR_FORMULA_DEPRESSION
QQBAR_FORMULA_DEFLATION
QQBAR_FORMULA_SEPARATION
QQBAR_FORMULA_ALL

Methods for `qqbar_get_fexpr_formula`: recognise Gaussian rationals $a +
bi$; write quadratic numbers as $a + b \sqrt{d}$; recognise elements of
cyclotomic fields, by trying plausible fields, finding candidate
representations with LLL and verifying them exactly (heuristic, not
guaranteed to find all); solve cubics, quartics and quintics (not yet
implemented); depression (not yet implemented); deflation, which
handles numbers $a^{1/n}$ where $a$ has a closed form; separation into
real and imaginary parts $a + bi$ or modulus and sign $m s$, tried as a
fallback when the other methods fail. `QQBAR_FORMULA_ALL` enables all
methods, which can be expensive.
:::

::: macro
QQBAR_FORMULA_EXP_FORM
QQBAR_FORMULA_TRIG_FORM
QQBAR_FORMULA_RADICAL_FORM
QQBAR_FORMULA_AUTO_FORM

The form used for elements of cyclotomic fields. The automatic form (0,
the default) gives radicals for numbers of low degree, trigonometric
functions for real numbers and complex exponentials for nonreal
numbers; the other flags (not fully implemented) request exponentials,
trigonometric functions or radicals.
:::

## Internal functions {#internal}

::: function
void qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op)

For nonconstant `A` and `B`, sets `res` to a polynomial whose roots are
$a + b$ (`op` = 0), $a - b$ (1), $a b$ (2) or $a / b$ (3) for all roots
$a$ of `A` and $b$ of `B`, with multiplicity. For `op` = 3, `B` must not
have the root zero. The polynomial is computed from power sums with
power series arithmetic over $\mathbb{Q}$, see
[](algorithms/exact/qqbar-arithmetic).
:::

::: function
void qqbar_fmpz_poly_symmetric_composed_op(fmpz_poly_t res, const fmpz_poly_t A, int op)

For `A` of degree $d \ge 2$ with roots $a_1, \ldots, a_d$, sets `res` to
a nonzero integer multiple of the polynomial of degree $d(d-1)/2$ with
roots $a_i + a_j$ (`op` = 0), $(a_i - a_j)^2$ (`op` = 1) or $a_i a_j$
(`op` = 2) for $1 \le i < j \le d$. For two distinct roots $x$, $y$ of an
irreducible `A` (such as complex conjugates), this annihilates $x + y$,
$(x - y)^2$ or $x y$ with half the degree of the polynomial from
`qqbar_fmpz_poly_composed_op`.
:::

::: function
void qqbar_binary_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)

Sets `res` to $x + y$, $x - y$, $x y$ or $x / y$ for `op` = 0, 1, 2, 3.
Trivial cases such as rational operands are not handled here (the
public arithmetic functions do that), but structure is exploited to
avoid factoring a composed polynomial of degree $\deg x \deg y$ when
possible.

::: details
In order: equal operands are handled directly; if an operand is
quadratic, closed formulas (both quadratic) or a norm of degree $2 \deg
x$ certified irreducible modulo a prime are used; for operands of equal
degree at least 4, a rational result is guessed and verified; if one
operand lies in the field generated by the other, or both in a common
cyclotomic field (detected by modular tests), the operands are
expressed in that field with LLL and exact verification and the result
is computed by linear algebra there, or the minimal polynomial of the
result is guessed with LLL and proved by a numerical evaluation at
sufficient precision, all within a budget proportional to the
estimated cost of the fallback; distinct conjugates use the symmetric
composed polynomial of half degree; otherwise the composed polynomial
is factored and the factor vanishing at the numerical value is
selected. See [](algorithms/exact/qqbar-arithmetic).
:::
:::

::: function
int _qqbar_validate_uniqueness(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)

Given a ball `z` known to contain at least one root of `poly`, tries to
prove that it contains exactly one, by one step of the interval Newton
method on a slightly enlarged ball. On success sets `res` to the
resulting (possibly smaller) enclosure of the root and returns 1;
otherwise returns 0. The working precision is chosen from the accuracy
of `z`, up to `max_prec` bits (a value of 0 selects it automatically).
Uniqueness on the enlarged ball implies uniqueness on `z` but not
existence, which must be known beforehand.
:::

::: function
int _qqbar_validate_existence_uniqueness(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)

Like `_qqbar_validate_uniqueness`, but for an arbitrary ball `z`: proves
that `z` contains exactly one root of `poly` (one interval Newton step
whose result lies in the interior of `z`). For a real (or purely
imaginary) ball the step is done on the real line (or imaginary axis),
which proves that the root lies there.
:::

::: function
void _qqbar_enclosure_raw(acb_t res, const fmpz_poly_t poly, const acb_t z, slong prec)
void qqbar_enclosure_raw(acb_t res, const qqbar_t x, slong prec)

Sets `res` to an enclosure of the root of `poly` isolated by `z` (or of
`x`) with about `prec` bits of accuracy, by interval Newton steps with
doubling precision; if a step fails to improve the accuracy enough, the
roots are recomputed from scratch at higher precision. If `z` is already
accurate enough, it is returned without rounding. The result is not
polished as by `qqbar_get_acb`.
:::

::: function
int _qqbar_acb_lindep(fmpz * rel, acb_srcptr vec, slong len, int check, slong prec)

Tries to find an integer vector `rel` with $\sum_i \mathit{rel}_i
\mathit{vec}_i \approx 0$ using LLL. The working precision is the
smaller of `prec` and the relative accuracy of `vec`; 95% of the bits
go into the LLL matrix and the remaining 5% are used to check the
relation by evaluating it in ball arithmetic. If `check` is set, the
return value says whether this check passed; otherwise it says whether
the algorithm ran normally (it fails, for example, for non-finite
input). A passed check is a heuristic filter against spurious
relations, not a proof.
:::

::: function
int _qqbar_fast_detect_simple_principal_surd(const qqbar_t x)

Returns 1 if `x` is cheaply recognised, from its minimal polynomial $q
t^d - p$ and its enclosure, as the principal root $(p/q)^{1/d}$ of a
positive rational number with $d \ge 2$, and 0 otherwise (0 does not
prove the contrary). Used to speed up `qqbar_root_ui`.
:::
