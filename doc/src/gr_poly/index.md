---
title: gr_poly.h
header: gr_poly.h
description: Dense univariate polynomials over any generic ring, with the generic algorithms for multiplication, division, GCDs, resultants, evaluation, composition, power series and factorisation over finite fields that most of FLINT's polynomial modules build on.
category: Generic rings
order: 7
---

This module implements dense univariate polynomials $f \in R[x]$ over a
ring $R$ given by a generic context `gr_ctx_t` (see [](gr)). A
`gr_poly_t` stores its coefficients contiguously, in the representation
of the elements of $R$, so the same code serves $\Z$, $\mathbb{Q}$, $\Z/n\Z$,
finite fields, real and complex balls, floating-point numbers, exact real
and complex numbers, polynomial rings and any other ring with a `gr`
context. The module contains the generic algorithms for multiplication,
division with remainder and Newton iteration, GCDs and resultants,
evaluation and interpolation, composition, power series functions,
squarefree factorisation, factorisation and irreducibility testing over
finite fields, and arithmetic modulo a fixed polynomial.

The central operations are methods of the coefficient context: a ring
can replace the generic implementation of multiplication, division,
GCD, power series inversion and a few others by its own (for example
`nmod_poly` multiplication for $\Z/n\Z$, `arb_poly` multiplication for
real balls, [](algorithms/ball-arithmetic/polynomial-multiplication)), and every algorithm in this module that multiplies calls
that method. The overview [](algorithms/polynomials/overview) says which
ring uses which algorithm.

The guide [](gr_poly/getting-started) explains contexts, status codes and
the handling of inexact coefficients with complete programs.

::: summary algorithms/polynomials/overview :::

## Types and conventions {#conventions}

::: type
gr_poly_struct
gr_poly_t

A polynomial over the ring given by a `gr_ctx_t`. The structure has the
fields `coeffs` (a pointer to `alloc` elements of the ring), `alloc`
and `length`; the coefficient of $x^i$ is the element at index $i$ of
`coeffs`. A `gr_poly_t` is an array of one `gr_poly_struct`, so it is
passed by reference.
:::

The following conventions hold throughout the module.

- **The context is the coefficient ring.** The `ctx` argument of every
  function describes $R$, not $R[x]$. A context for the polynomial ring
  itself, created with `gr_ctx_init_gr_poly`, is needed only to use
  polynomials as elements of another generic structure (matrices of
  polynomials, polynomials over polynomials) or with the generic
  element functions such as `gr_add`; its elements are `gr_poly_struct`
  objects that the functions here accept directly.
- **Status codes.** Every function returning `int` returns `GR_SUCCESS`
  or a combination of the flags `GR_DOMAIN` (the result does not exist
  in $R[x]$: division by zero, a quotient with remainder that is not
  exact, a non-invertible leading coefficient) and `GR_UNABLE` (the
  result may exist but cannot be computed: an undecidable zero test, an
  unimplemented case, an overflow). On failure the output is undefined
  but remains a valid object that can be cleared. Functions returning
  `truth_t` answer `T_TRUE`, `T_FALSE` or `T_UNKNOWN`.
- **Weak normalisation.** A polynomial is normalised by removing leading
  coefficients for which `gr_is_zero` returns `T_TRUE`. Over rings with
  undecidable equality, such as balls and exact real numbers, a leading
  coefficient that may or may not be zero (for example $[0 \pm 10^{-10}]$)
  is kept, so `length` is only an upper bound for the degree plus one.
  Functions that need the exact degree, such as division, GCD,
  resultants and factorisation, test the leading coefficient and return
  `GR_UNABLE` when it is not provably nonzero, instead of returning a
  wrong result. See [](gr_poly/getting-started#inexact) and
  [](algorithms/generic/inexact#weak-normalisation).
- **Memory management and aliasing.** Functions without a leading
  underscore resize their outputs and allow any output to be the same
  object as any input, unless stated otherwise. Several outputs of one
  function must be distinct objects.
- **Underscore functions** take coefficient arrays with explicit
  lengths, write a fixed number of coefficients (including leading
  zeros, which are not removed), and leave memory management to the
  caller. They usually require positive lengths, assume that the
  leading coefficients of the inputs are nonzero where degrees matter,
  and do not allow the output to overlap the inputs. Each block states
  its conditions.
- **Coefficient domains.** Most functions assume that $R$ is a
  commutative ring, or an approximate commutative ring such as a
  floating-point type; over noncommutative rings they can return
  meaningless results without an error flag. Functions with stronger
  requirements (integral domain, field, unique factorisation domain,
  finite field, characteristic zero) state them; most of them test the
  context predicates (`gr_ctx_is_field`, `gr_ctx_is_finite` and so on)
  and return `GR_UNABLE` or `GR_DOMAIN` when the requirement cannot be
  verified, while some only document the assumption.

The `gr_poly_t` layout is the same as that of `fmpz_poly_t`,
`fq_poly_t`, `fq_nmod_poly_t`, `fq_zech_poly_t`, `arb_poly_t`,
`acb_poly_t` and `ca_poly_t`, so a polynomial of one of these types can
be cast to a `gr_poly_t` over the matching context and passed to the
functions here. `nmod_poly_t` and `fmpz_mod_poly_t` store their modulus
in the polynomial and `fmpq_poly_t` uses a common denominator, so they
must be converted.

## Memory management {#memory}

::: function
void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx)
void gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx)

Initialises `poly` as the zero polynomial. `gr_poly_init2` also
allocates (and initialises) space for `len` coefficients.
:::

::: function
void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx)

Clears the coefficients and frees the memory of `poly`.
:::

::: function
gr_ptr gr_poly_coeff_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)
gr_srcptr gr_poly_coeff_srcptr(const gr_poly_t poly, slong i, gr_ctx_t ctx)

Returns a pointer to the coefficient of $x^i$. The index must be less
than the allocated size; coefficients at indices from `length` on are
initialised elements but not part of the polynomial.
:::

::: function
gr_ptr gr_poly_entry_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)
gr_srcptr gr_poly_entry_srcptr(const gr_poly_t poly, slong i, gr_ctx_t ctx)

Deprecated names of `gr_poly_coeff_ptr` and `gr_poly_coeff_srcptr`.
:::

::: function
slong gr_poly_length(const gr_poly_t poly, gr_ctx_t ctx)

Returns the length of `poly`, that is one more than the index of the
leading stored coefficient, or 0 for the zero polynomial. Over rings
with undecidable equality this bounds the degree plus one.
:::

::: function
void gr_poly_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)

Swaps `poly1` and `poly2` in constant time by exchanging the structures.
:::

::: function
void gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)
void _gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

Ensures that `poly` has space for at least `len` coefficients, without
changing its length or coefficients. The inline function
`gr_poly_fit_length` calls `_gr_poly_fit_length` only when the current
allocation is too small; the latter grows the allocation to at least
twice its old size, so that repeated calls cost amortised linear time.
:::

::: function
void _gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

Sets the length of `poly` to `len`, which must not exceed the
allocation, without normalising. When the polynomial is shortened, the
coefficients beyond the new length are set to zero, which frees the
memory of heap-allocated elements.
:::

::: function
void _gr_poly_set_length_normalise(gr_poly_t poly, slong len, gr_ctx_t ctx)
void _gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)

Set the length of `poly` to `len` (respectively keep the current
length) and then decrease it while the leading coefficient is provably
zero. This is the usual last step after an underscore function has
written coefficients. Both call the ring method
`GR_METHOD_POLY_SET_LENGTH_NORMALISE`, whose default
`gr_generic_poly_set_length_normalise` zeroes the dropped coefficients
like `_gr_poly_set_length`; rings with shallow elements (`nmod`,
`fmpz_mod`, `mpn_mod`, `nmod8`, `nmod32`) override it to compute only
the new length.
:::

::: function
void gr_generic_poly_set_length_normalise(gr_poly_struct * poly, slong len, gr_ctx_t ctx)

The default implementation of `_gr_poly_set_length_normalise`, declared
in `gr_generic.h`.
:::

## Vectors of polynomials {#vectors}

::: type
gr_poly_vec_struct
gr_poly_vec_t

A resizable vector of polynomials over one ring, with the fields
`entries` (an array of `alloc` objects of type `gr_poly_struct`),
`alloc` and `length`. The `ctx` argument of the functions below is the
coefficient ring of the polynomials; no context for the polynomial ring
is needed. A `gr_poly_vec_t` has the layout of a `gr_vec_t` whose
elements are polynomials, and the two can be cast to each other.
Factorisation functions return their factors in this type.
:::

::: function
void gr_poly_vec_init(gr_poly_vec_t vec, slong len, gr_ctx_t ctx)
void gr_poly_vec_clear(gr_poly_vec_t vec, gr_ctx_t ctx)

Initialises `vec` as a vector of `len` zero polynomials ($\mathit{len}
\ge 0$), respectively clears all its polynomials and frees the vector.
:::

::: function
gr_poly_struct * gr_poly_vec_entry_ptr(gr_poly_vec_t vec, slong i, gr_ctx_t ctx)
const gr_poly_struct * gr_poly_vec_entry_srcptr(const gr_poly_vec_t vec, slong i, gr_ctx_t ctx)
slong gr_poly_vec_length(const gr_poly_vec_t vec, gr_ctx_t ctx)

Return a pointer to the polynomial at index $0 \le i < \mathit{length}$,
and the length of the vector.
:::

::: function
void gr_poly_vec_fit_length(gr_poly_vec_t vec, slong len, gr_ctx_t ctx)
void gr_poly_vec_set_length(gr_poly_vec_t vec, slong len, gr_ctx_t ctx)

`gr_poly_vec_fit_length` allocates space for at least `len`
polynomials without changing the length. `gr_poly_vec_set_length`
changes the length to $\mathit{len} \ge 0$, clearing removed entries and
appending zero polynomials as needed.
:::

::: function
int gr_poly_vec_set(gr_poly_vec_t res, const gr_poly_vec_t src, gr_ctx_t ctx)

Sets `res` to a copy of `src`.
:::

::: function
int gr_poly_vec_append(gr_poly_vec_t vec, const gr_poly_t f, gr_ctx_t ctx)
void gr_poly_vec_append_swap(gr_poly_vec_t vec, gr_poly_t f, gr_ctx_t ctx)

Append `f` to the vector. `gr_poly_vec_append` copies it;
`gr_poly_vec_append_swap` swaps it into the new entry without copying,
so that `f` becomes the zero polynomial.
:::

## Assignment and conversion {#assignment}

::: function
int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

Sets `res` to a copy of `src`.
:::

::: function
int gr_poly_set_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx)
int gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx)
int gr_poly_set_gr_poly_other(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)

Set `res` to the image of an integer polynomial, a rational polynomial,
or a polynomial over the ring `x_ctx`, converting each coefficient
into $R$ with `gr_set_fmpz` or `gr_set_other`; for a rational
polynomial the integer numerators are converted and multiplied by the
inverse of the common denominator in $R$. The conversion fails (with
`GR_DOMAIN` or `GR_UNABLE`) if a coefficient has no image, for example
a fraction whose denominator is divisible by $p$ in $\Z/p\Z$.
:::

::: function
int gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx)
int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx)
int gr_poly_gen(gr_poly_t poly, gr_ctx_t ctx)

Set `poly` to $0$, $1$, $-1$ or the generator $x$. In the zero ring
the results are all the zero polynomial.
:::

::: function
int gr_poly_set_scalar(gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_set_si(gr_poly_t poly, slong c, gr_ctx_t ctx)
int gr_poly_set_ui(gr_poly_t poly, ulong c, gr_ctx_t ctx)
int gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
int gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)

Set `poly` to the constant polynomial $c$, where $c$ is an element of
$R$ or an integer or rational number converted into $R$.
:::

::: function
int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

Sets `res` to the reversal $x^{n-1} f(1/x)$ of the polynomial $f$ formed
by the first $n$ coefficients of `poly` (padded with zeros if `poly`
is shorter). The underscore version requires $\mathit{len} \le n$,
writes $n$ coefficients and allows `res` to be `poly`.
:::

::: function
int gr_poly_truncate(gr_poly_t res, const gr_poly_t poly, slong newlen, gr_ctx_t ctx)

Sets `res` to `poly` reduced modulo $x^{\mathit{newlen}}$.
:::

## Input and output {#io}

::: function
int _gr_poly_write(gr_stream_t out, gr_srcptr poly, slong len, const char * x, gr_ctx_t ctx)
int gr_poly_write(gr_stream_t out, const gr_poly_t poly, const char * x, gr_ctx_t ctx)
int _gr_poly_get_str(char ** res, gr_srcptr f, slong len, const char * x, gr_ctx_t ctx)
int gr_poly_get_str(char ** res, const gr_poly_t poly, const char * x, gr_ctx_t ctx)
int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)

Write the polynomial as an expression in the variable named `x`, in
increasing order of degree with zero terms omitted, to a stream, to a
newly allocated string `*res` (to be freed with `flint_free`), or to
standard output with the variable `x`.
:::

::: function
int _gr_poly_set_str(gr_ptr res, const char * s, const char * x, slong len, gr_ctx_t ctx)
int gr_poly_set_str(gr_poly_t res, const char * s, const char * x, gr_ctx_t ctx)

Parse a polynomial from the expression `s` in the variable named `x`.
The underscore version writes exactly `len` coefficients, zero-padding
a shorter result, and returns `GR_UNABLE` if the parsed polynomial is
longer than `len` (intermediate values may be longer). Returns
`GR_UNABLE` if the string cannot be parsed.

::: details
The parser evaluates the expression with polynomial arithmetic and is
not optimised for long inputs; it can take time quadratic in the
length.
:::
:::

::: function
int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)

Sets `poly` to a random polynomial of length at most `len`, with
coefficients from `gr_randtest`. The distribution is meant for testing:
it is not uniform and produces many special values.
:::

## Comparisons and predicates {#predicates}

::: function
truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

Returns whether the two polynomials are equal, comparing the common
coefficients and testing the extra coefficients of the longer one for
zero. The underscore version requires $\mathit{len1} \ge \mathit{len2}$.
:::

::: function
truth_t gr_poly_is_zero(const gr_poly_t poly, gr_ctx_t ctx)
truth_t gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx)
truth_t gr_poly_is_gen(const gr_poly_t poly, gr_ctx_t ctx)
truth_t gr_poly_is_scalar(const gr_poly_t poly, gr_ctx_t ctx)

Return whether `poly` is $0$, $1$, the generator $x$, or a constant
polynomial. Coefficients that are not provably zero or one make the
answer `T_UNKNOWN`.
:::

## Coefficients {#coefficients}

::: function
int gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong c, gr_ctx_t ctx)
int gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong c, gr_ctx_t ctx)
int gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t c, gr_ctx_t ctx)
int gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t c, gr_ctx_t ctx)

Set the coefficient of $x^n$ to $c$, extending `poly` with zeros if
$n \ge$ `length`, and normalise.
:::

::: function
int gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

Sets `res` to the coefficient of $x^n$, which is zero for $n \ge$
`length`.
:::

## Addition and multiplication {#arithmetic}

::: function
int gr_poly_neg(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

Sets `res` to $-f$.
:::

::: function
int _gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
int _gr_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

Set `res` to the sum or difference. The underscore versions write
$\max(\mathit{len1}, \mathit{len2})$ coefficients and allow `res` to
be either input.
:::

::: function
int _gr_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

Sets `res` to the product. The underscore version requires positive
lengths, writes $\mathit{len1} + \mathit{len2} - 1$ coefficients and
does not allow aliasing; it is `_gr_poly_mullow` with the full length,
so the algorithm is that of the ring's `GR_METHOD_POLY_MULLOW` method.
:::

::: function
int _gr_poly_mullow(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
int gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong len, gr_ctx_t ctx)

Sets `res` to the product truncated to length `len`, that is reduced
modulo $x^{\mathit{len}}$. The underscore version requires positive
lengths and $\mathit{len} \le \mathit{len1} + \mathit{len2} - 1$,
writes `len` coefficients and does not allow aliasing. It calls the
ring method `GR_METHOD_POLY_MULLOW`, which defaults to classical
multiplication and which most rings with fast multiplication override.
:::

::: function
int _gr_poly_mulmid(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
int gr_poly_mulmid(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong nlo, slong nhi, gr_ctx_t ctx)

Sets `res` to the middle product: the coefficients of $x^{\mathit{nlo}},
\ldots, x^{\mathit{nhi} - 1}$ of the product, shifted down by
$\mathit{nlo}$. The underscore version requires positive lengths and
$0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1} + \mathit{len2} -
1$, writes $\mathit{nhi} - \mathit{nlo}$ coefficients and does not allow
aliasing. It calls the ring method `GR_METHOD_POLY_MULMID`; see
[](algorithms/polynomials/middle-product).
:::

## Multiplication algorithms {#mul-algorithms}

The functions in this section implement one algorithm each, for testing
and for use in ring methods. Those that recurse (Karatsuba, Toom-3) call
`_gr_poly_mul` for their sub-products, so they recurse only if the ring
method dispatches back to them above a cutoff.

::: function
int _gr_poly_mulmid_classical(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
int gr_poly_mulmid_classical(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong nlo, slong nhi, gr_ctx_t ctx)
int _gr_poly_mullow_classical(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_mullow_classical(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)

Multiply with the classical algorithm, computing each output coefficient
as one dot product (`_gr_vec_dot_rev`), so that rings with fast dot
products (delayed reduction, exact accumulation) benefit. Squaring uses
the symmetry of the product to halve the number of multiplications.
:::

::: function
int _gr_poly_mul_karatsuba(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_mul_karatsuba(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

One level of Karatsuba multiplication [@KO1962]: three products of half
length with sub-products computed by `_gr_poly_mul`. It is not
optimised for unbalanced operands. The underscore version requires
positive lengths and does not allow aliasing.
:::

::: function
int _gr_poly_mul_toom33(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_mul_toom33(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

One level of balanced Toom-3 multiplication [@Too1963; @Coo1966]:
evaluation at $0, 1, -1, 2, \infty$ with Bodrato's sequence of
operations and five products of a third of the length, computed by
`_gr_poly_mul`. Requires a commutative ring with exact division by 2
and 3; not optimised for squaring. The underscore version requires
positive lengths and does not allow aliasing.
:::

::: function
int _gr_poly_mullow_toom_serial(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong trunc, slong num_points, slong splitting_degree, gr_ctx_t ctx)
int gr_poly_mullow_toom_serial(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, slong num_points, slong splitting_degree, gr_ctx_t ctx)

Toom-Cook multiplication of arbitrary order with low memory use: the
inputs are split into blocks of $m$ coefficients, viewed as
polynomials in $t = x^m$, and the product is computed by evaluation at
$r$ points and interpolation. Exactly one of `num_points` ($r$) and
`splitting_degree` ($m$) must be nonzero; the other is chosen
automatically.

::: details
For balanced operands, $r = 2k - 1$ points correspond to splitting both
inputs into $k$ blocks; $r = 3$ (or $m = \lceil \max(\mathit{len1},
\mathit{len2}) / 2 \rceil$) gives Karatsuba multiplication with
evaluation at $0, 1, \infty$. Even $r$ is mainly useful for unbalanced
products. The points are $0, 1, -1, 2, -2, \ldots, \infty$, and
evaluation and interpolation multiply by the Vandermonde matrix and its
inverse. In characteristic zero the inverse is computed with a common
denominator that is removed at the end by an exact division; in small
positive characteristic, where the Vandermonde matrix is singular,
`GR_UNABLE` is returned.

The evaluations are done one point at a time and the output is
accumulated incrementally, using $4m$ temporary elements for
polynomials and $2r^2$ elements or integers for the matrices. This
suits huge products under memory constraints, unlike a usual
Toom-Cook implementation that keeps many temporary buffers to save
scalar operations. The underscore version requires positive lengths
and does not allow aliasing.
:::
:::

::: function
int _gr_poly_mulmid_bivariate_KS(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
int gr_poly_mulmid_bivariate_KS(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong nlo, slong nhi, gr_ctx_t ctx)
int _gr_poly_mullow_bivariate_KS(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_mullow_bivariate_KS(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)

For coefficients that are themselves polynomials (a context created by
`gr_ctx_init_gr_poly` or `gr_ctx_init_fmpz_poly`), reduce the
bivariate product to one univariate product over the inner coefficient
ring by Kronecker substitution $y \mapsto x^k$
([](algorithms/polynomials/kronecker)). Returns `GR_UNABLE` if the
coefficients are not polynomials. The polynomial-ring contexts use this
for their `GR_METHOD_POLY_MULLOW` method when all lengths are at least 5
and the inner ring has a non-classical multiplication.
:::

::: function
int _gr_poly_mulmid_complex_reorder(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong nlo, slong nhi, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
int gr_poly_mulmid_complex_reorder(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong nlo, slong nhi, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
int _gr_poly_mullow_complex_reorder(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
int gr_poly_mullow_complex_reorder(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)

For complex coefficients stored as pairs of elements of `real_ctx`
(real part, imaginary part), compute the product from products of real
polynomials: the inputs are split into their real and imaginary parts
and the products are computed in `real_ctx`.

::: details
If one input is detected to be purely real or purely imaginary, one or
two real products suffice. Otherwise four real products are used, or
three if `karatsuba` is set (Gauss's trick, which can lose accuracy in
floating-point arithmetic). A square takes two real squarings and one
real product, or three squarings if `karatsuba` is set.
:::
:::

## Scalar arithmetic {#scalar}

::: function
int gr_poly_add_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_add_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
int gr_poly_add_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
int gr_poly_add_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
int gr_poly_add_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)
int gr_poly_sub_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_sub_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
int gr_poly_sub_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
int gr_poly_sub_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
int gr_poly_sub_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)

Set `res` to $f + c$ or $f - c$ for a constant $c$, an element of $R$ or
an integer or rational number converted into $R$.
:::

::: function
int gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_scalar_mul(gr_poly_t res, gr_srcptr c, const gr_poly_t poly, gr_ctx_t ctx)
int gr_poly_mul_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
int gr_poly_mul_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
int gr_poly_mul_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
int gr_poly_mul_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)

Set `res` to $f c$, or to $c f$ for `gr_poly_scalar_mul` (the two differ
over noncommutative rings).
:::

::: function
int gr_poly_addmul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_submul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

Set `res` to $\mathit{res} + f c$ or $\mathit{res} - f c$.
:::

::: function
int gr_poly_div_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_divexact_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

Set `res` to $f / c$, dividing each coefficient with `gr_div`
(`GR_DOMAIN` if a quotient does not exist in $R$), respectively with
`gr_divexact`, which assumes that every division is exact and may then
be faster (for example over $\Z$).
:::

## Powers {#powers}

::: function
int _gr_poly_pow_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
int gr_poly_pow_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)
int _gr_poly_pow_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
int gr_poly_pow_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)

Set `res` to $f^e$ by binary exponentiation (`pow_ui` currently calls
the `binexp` version). The underscore versions write $e(\mathit{flen} -
1) + 1$ coefficients and require $\mathit{flen} \ge 1$ and $e \ge 1$.
The non-underscore versions return `GR_UNABLE` if the length of the
result would overflow an `slong`.
:::

::: function
int gr_poly_pow_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx)

Sets `res` to $f^e$. A negative exponent is allowed only for a constant
polynomial with an invertible constant (`GR_DOMAIN` otherwise), and an
exponent that does not fit in a word gives `GR_UNABLE` unless $f$ is
constant.
:::

::: function
int _gr_poly_pow_series_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
int gr_poly_pow_series_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)
int _gr_poly_pow_series_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
int gr_poly_pow_series_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

Set `res` to $f^e \bmod x^{\mathit{len}}$ by binary exponentiation with
truncated products (`pow_series_ui` currently calls the `binexp`
version). The underscore versions write `len` coefficients.
:::

::: function
int _gr_poly_pow_series_fmpq_recurrence(gr_ptr h, gr_srcptr f, slong flen, const fmpq_t exp, slong len, int precomp, gr_ctx_t ctx)
int gr_poly_pow_series_fmpq_recurrence(gr_poly_t res, const gr_poly_t poly, const fmpq_t exp, slong len, gr_ctx_t ctx)

Set `res` to the power series $f^{p/q} \bmod x^{\mathit{len}}$ for a
rational exponent $p/q$, using the recurrence obtained from $q f h' = p
f' h$ for $h = f^{p/q}$ (often attributed to J. C. P. Miller, see
[@Knu1997, §4.7]), which costs $O(\mathit{flen} \cdot \mathit{len})$
operations. The constant term $f_0^{p/q}$ is computed with
`gr_pow_fmpq`, and the recurrence divides by $i q f_0$ for $i = 1,
\ldots, \mathit{len} - 1$, so $f_0$ and these integers must be
invertible (or the divisions exact).

::: details
The underscore version writes `len` coefficients. Its flags `precomp`
allow precomputed data: if bit 0 is set, `h[0]` already contains
$f_0^{p/q}$; if bit 1 is set, `h[i]` contains $1/i$ for $1 \le i <
\mathit{len}$ on input. Over integral domains of characteristic zero and
integer exponents, the divisions are exact divisions.
:::
:::

## Shifting {#shifting}

::: function
int _gr_poly_shift_left(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
int gr_poly_shift_left(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)
int _gr_poly_shift_right(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
int gr_poly_shift_right(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

Set `res` to $f x^n$, respectively to $\lfloor f / x^n \rfloor$ (the
coefficients of degree $n$ and higher, moved down). The underscore
versions write $\mathit{len} + n$ respectively $\mathit{len} - n$
coefficients ($n < \mathit{len}$ for the right shift) and allow
aliasing.
:::

## Division with remainder {#divrem}

Given $A, B \in R[x]$ with $B \ne 0$ and an invertible leading
coefficient $\operatorname{lc}(B)$, there are unique $Q, R$ with $A = BQ
+ R$ and $\deg R < \deg B$. If $\operatorname{lc}(B)$ is not invertible,
the quotient still exists when every leading coefficient of the
partial remainders is divisible by $\operatorname{lc}(B)$ during long
division; the *basecase* and *divconquer* algorithms then succeed and
return `GR_DOMAIN` exactly when such a division fails, while the
*newton* algorithms need $\operatorname{lc}(B)^{-1}$. Newton division
reverses the polynomials and multiplies by a power series inverse
([](algorithms/polynomials/newton#divrem)).

::: function
int _gr_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_divrem_basecase_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
int _gr_poly_divrem_basecase_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int _gr_poly_divrem_basecase(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_divrem_basecase(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_divrem_divconquer_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
int _gr_poly_divrem_divconquer_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
int _gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
int gr_poly_divrem_divconquer(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

Set $Q$ and $R$ to the quotient and remainder of $A$ by $B$. The
non-underscore functions return `GR_DOMAIN` if $B$ is zero and
`GR_UNABLE` if the leading coefficient of $B$ is not provably nonzero.
`gr_poly_divrem` calls the ring method `GR_METHOD_POLY_DIVREM`.

::: details
The *basecase* versions perform long division. The *divconquer*
versions divide blocks of $2 \mathit{lenB} - 1$ coefficients of the
dividend by $B$ recursively, with two divisions by the upper half of
$B$ and products by its lower half, and use long division for divisors
shorter than `cutoff`. The *newton* versions use Newton division. The *preinv1* versions take a precomputed
inverse `invB` of the leading coefficient, the *noinv* versions divide
by the leading coefficient with `gr_div` at each step, and the plain
versions try to invert the leading coefficient and fall back to the
*noinv* algorithm.

The default method `_gr_poly_divrem_generic` uses long division when
$\mathit{lenB} = 1$ or $\mathit{lenA} - \mathit{lenB} = 1$, and over
finite rings also when the divisor or the quotient has at most 10
coefficients; otherwise it uses Newton division, falling back to long
or divide-and-conquer division without inversion when the leading
coefficient is not invertible. Rings with fast arithmetic override the
method: `nmod`, `fmpz_mod`, `mpn_mod`, `fq`, `fq_nmod`, `fq_zech`,
`arb` and `acb` call the division of their own polynomial modules.

The underscore functions assume:

- $\mathit{lenA} \ge \mathit{lenB} \ge 1$ and a nonzero leading
  coefficient of $B$;
- $Q$ has room for $\mathit{lenA} - \mathit{lenB} + 1$ and $R$ for
  $\mathit{lenB} - 1$ coefficients;
- $Q$ is not aliased with $A$ or $B$, and $R$ is not aliased with $B$;
- $R$ may be aliased with $A$, in which case all $\mathit{lenA}$
  coefficients are used as scratch space and only the low
  $\mathit{lenB} - 1$ are meaningful on output.
:::
:::

::: function
int _gr_poly_div(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_div(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_div_basecase_preinv1(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
int _gr_poly_div_basecase_noinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int _gr_poly_div_basecase(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_div_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_div_divconquer_preinv1(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
int _gr_poly_div_divconquer_noinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
int _gr_poly_div_divconquer(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
int gr_poly_div_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
int _gr_poly_div_newton(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_div_newton(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

The same as the *divrem* functions, but compute only the quotient.
`gr_poly_div` calls the ring method
`GR_METHOD_POLY_DIV`; its default `_gr_poly_div_generic` uses long
division when the divisor or the quotient has at most 10 coefficients
and Newton division otherwise (falling back to divide-and-conquer
division without inversion).
:::

::: function
int _gr_poly_rem(gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_rem(gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

The same as the *divrem* functions, but return only the remainder (the
quotient is computed in a temporary).
:::

::: function
int _gr_poly_div_newton_n_preinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr Binv, slong lenBinv, gr_ctx_t ctx)
int gr_poly_div_newton_n_preinv(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, const gr_poly_t Binv, gr_ctx_t ctx)
int _gr_poly_divrem_newton_n_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr Binv, slong lenBinv, gr_ctx_t ctx)
int gr_poly_divrem_newton_n_preinv(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, const gr_poly_t Binv, gr_ctx_t ctx)

Division with remainder (or quotient only) given `Binv`, the inverse of
the reversal of $B$ modulo $x^{\mathit{lenB}}$, which can be computed
once with `gr_poly_inv_series` and reused for many divisions by $B$.
The quotient is the reversal of a single truncated product, and the
remainder costs one more product. The non-underscore versions return
`GR_UNABLE` if $\mathit{lenA} > 2 \mathit{lenB} - 2$.

::: details
The underscore versions require $\mathit{lenA} - \mathit{lenB} + 1 \le
\mathit{lenBinv}$ for a correct quotient; all versions return
`GR_UNABLE` when $\mathit{lenBinv} = 0$ (no inverse available), as
when the caller failed to invert $B$. `B` is not read by the quotient-only
version. For dividends of any length and for repeated reductions,
`gr_poly_preinv_divrem` with a `gr_poly_preinv_t` is more convenient.
:::
:::

## Power series inversion and division {#series-division}

For the Newton versions, the parameter `cutoff` is the length below
which the basecase algorithm is used: lengths $n < \mathit{cutoff}$
use the basecase, and $\mathit{cutoff} = n$ gives exactly one Newton
step from length $\lceil n/2 \rceil$. Calling a Newton function with $n
< \mathit{cutoff}$ can be much slower than calling the basecase
directly. The Newton iterations use middle products when the ring has
a non-generic `GR_METHOD_POLY_MULMID`, and division uses the
Karp-Markstein trick of multiplying by the dividend in the last step
[@KM1997]; see [](algorithms/polynomials/newton).

::: function
int _gr_poly_inv_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
int gr_poly_inv_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
int _gr_poly_inv_series_basecase_preinv1(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr Ainv, slong len, gr_ctx_t ctx)
int _gr_poly_inv_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
int gr_poly_inv_series_basecase(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
int _gr_poly_inv_series_newton(gr_ptr res, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
int gr_poly_inv_series_newton(gr_poly_t res, const gr_poly_t A, slong len, slong cutoff, gr_ctx_t ctx)

Set `res` to the power series inverse $A^{-1} \bmod x^{\mathit{len}}$.
Returns `GR_DOMAIN` if the constant term of $A$ is not invertible
(for example if $A = 0$, except in the zero ring).

::: details
The basecase computes the coefficients one by one from $A \cdot A^{-1}
= 1$ with dot products, in $O(\mathit{Alen} \cdot \mathit{len})$
operations; the *preinv1* version takes the inverse `Ainv` of the
constant term. The Newton version doubles the precision with $B \gets
B + B(1 - AB)$. `gr_poly_inv_series` calls the ring method
`GR_METHOD_POLY_INV_SERIES`, whose default uses the basecase when
$\mathit{Alen} \le 8$ or the ring has only classical multiplication,
and Newton iteration otherwise. The underscore versions require
$\mathit{Alen} \ge 1$, write `len` coefficients and do not allow
aliasing.
:::
:::

::: function
int _gr_poly_div_series(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int gr_poly_div_series(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_basecase_preinv1(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_srcptr Binv, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int gr_poly_div_series_basecase(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_divconquer(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, slong cutoff, gr_ctx_t ctx)
int gr_poly_div_series_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
int _gr_poly_div_series_invmul(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int gr_poly_div_series_invmul(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_newton(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, slong cutoff, gr_ctx_t ctx)
int gr_poly_div_series_newton(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)

Set `res` to the power series quotient $A / B \bmod x^{\mathit{len}}$.
Returns `GR_DOMAIN` if $B = 0$; if the constant term of $B$ is not
invertible, the quotient may still exist, and the *basecase*,
*noinv* and *divconquer* algorithms find it whenever every coefficient
division is exact.

::: details
The *basecase* computes the quotient coefficient by coefficient
(*preinv1*: with a precomputed inverse `Binv` of the constant term;
*noinv*: dividing with `gr_div` at each step); *divconquer* splits the
quotient in halves; *invmul* multiplies $A$ by $B^{-1}$; *newton* uses
Newton iteration for $B^{-1}$ to half the length and a Karp-Markstein
step. `gr_poly_div_series` calls the ring method
`GR_METHOD_POLY_DIV_SERIES`, whose default uses the basecase when
$\min(\mathit{Blen}, \mathit{len}) \le 8$ or the ring has only
classical multiplication, and Newton iteration otherwise, falling back
to divide-and-conquer if the constant term is not invertible. The
underscore versions require positive lengths, write `len` coefficients
and do not allow aliasing.
:::
:::

## Exact division {#divexact}

These functions compute a quotient $Q = A / B$ that is known to be exact
in $R[x]$ (or in $R[[x]]/(x^n)$ for series division). If the division
is not exact, they may set $Q$ to an arbitrary polynomial and return
`GR_SUCCESS`. $R$ is assumed to be an integral domain, which is not
checked.

Exact division can start from the leading coefficients (classical
division) or from the trailing ones (power series division); which is
faster depends on which coefficient is cheaper to divide by, which is
hard to predict in general. The *bidirectional* algorithms compute the
upper and the lower half of the quotient from the two ends, which halves
the work in the basecase at the cost of possibly one more coefficient
inversion. The *noinv* versions use `gr_divexact` on the coefficients
instead of inverting the leading (or trailing) coefficient.

::: function
int _gr_poly_divexact(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_divexact(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_divexact_basecase_bidirectional(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
int gr_poly_divexact_basecase_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_divexact_bidirectional(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
int gr_poly_divexact_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_divexact_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
int _gr_poly_divexact_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
int gr_poly_divexact_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

Set $Q$ to the exact quotient $A / B$. The non-underscore functions
return `GR_DOMAIN` if $B = 0$ and `GR_UNABLE` if its leading
coefficient is not provably nonzero. `gr_poly_divexact` calls the ring
method `GR_METHOD_POLY_DIVEXACT`, whose default uses the basecase for
divisors of length at most 2 and the bidirectional algorithm otherwise;
`fmpz` uses the exact division of `fmpz_poly` and `nmod` that of
`nmod_poly`. The underscore functions require $\mathit{Alen} \ge
\mathit{Blen} \ge 1$ and write $\mathit{Alen} - \mathit{Blen} + 1$
coefficients.
:::

::: function
int _gr_poly_divexact_series_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int _gr_poly_divexact_series_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int gr_poly_divexact_series_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)

Set $Q$ to the power series quotient $A / B \bmod x^{\mathit{len}}$,
assuming that it exists, by basecase division starting from the
constant terms.
:::

::: function
int gr_poly_inv(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

Sets `res` to the inverse of `poly` in $R[x]$. This succeeds only for
a constant polynomial with an invertible constant (and for the zero
polynomial over the zero ring); any polynomial of length at least 2
gives `GR_DOMAIN`.

::: warning
Over a ring with nilpotent elements a nonconstant polynomial can be a
unit (in $(\Z/4\Z)[x]$, $(1 + 2x)^2 = 1$), and `GR_DOMAIN` is then
returned incorrectly.
:::
:::

## Square roots of power series {#sqrt-series}

::: function
int _gr_poly_sqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_sqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_sqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_sqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_sqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_sqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_sqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
int gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
int _gr_poly_rsqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_rsqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_rsqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_rsqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_rsqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_rsqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_rsqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
int gr_poly_rsqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)

Set `res` to the power series square root $f^{1/2}$ or reciprocal square
root $f^{-1/2}$ modulo $x^{\mathit{len}}$. The constant term is computed
with `gr_sqrt` or `gr_rsqrt` and must exist; the computation then
divides by 2 and by the constant term.

::: details
The *basecase* versions solve for one coefficient at a time with dot
products. The *miller* versions use the recurrence of
`gr_poly_pow_series_fmpq_recurrence` with exponent $\pm 1/2$, which
costs $O(\mathit{flen} \cdot \mathit{len})$ operations but divides by
$1, \ldots, \mathit{len} - 1$. The *newton* versions use the
third-order iteration for the reciprocal square root (each step about
triples the precision), and obtain the square root from the reciprocal
square root to half the length with a final Karp-Markstein step.

The default methods `GR_METHOD_POLY_SQRT_SERIES` and
`GR_METHOD_POLY_RSQRT_SERIES` use Miller's recurrence for short inputs
($\mathit{flen} \le 8 < \mathit{len}$) over exact rings; otherwise the
square root uses Newton iteration, and the reciprocal square root the
basecase for $\mathit{len} \le 8$ or classical multiplication and
Newton iteration otherwise. A `GR_DOMAIN` failure of the square root
is turned into `GR_UNABLE` unless the ring is a field of
characteristic other than 2, since in other rings a square root may
exist that the algorithm cannot find. The underscore versions require
$\mathit{flen} \ge 1$, write `len` coefficients and do not allow
aliasing.
:::
:::

## Evaluation {#evaluation}

::: function
int _gr_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_poly_evaluate(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to $f(x)$ for $x \in R$. This calls the ring method
`GR_METHOD_POLY_EVALUATE`, whose default `gr_generic_poly_evaluate`
uses Horner's rule; `nmod`, `fq`, `fq_nmod` and `fq_zech` override it
with the evaluation of their polynomial modules.
:::

::: function
int _gr_poly_evaluate_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)
int _gr_poly_evaluate_rectangular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)
int _gr_poly_evaluate_modular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_poly_evaluate_modular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $f(x)$ with a specific algorithm. Horner's rule uses
$\mathit{len} - 1$ multiplications by $x$. The *rectangular* version
(rectangular splitting, [@PS1973]) computes $x, \ldots, x^m$ for $m
\approx \sqrt{\mathit{len}}$, evaluates the blocks of $m$ coefficients
as dot products with these powers and combines the blocks by Horner's
rule in $x^m$, so that only $O(\sqrt{\mathit{len}})$ multiplications
are full products of ring elements; this pays off when the
coefficients are cheaper than $x$ (small integers and a large $x$,
scalars and a matrix or polynomial argument). The *modular* version
also computes the powers $(x^j)^0, \ldots, (x^j)^{k-1}$ with $j \approx k
\approx \sqrt{\mathit{len}}$ and evaluates everything with $j + 1$ dot
products, which suits rings where dot products are much cheaper than
separate multiplications (delayed modular reduction).
:::

::: function
int gr_generic_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)

The default implementation of `_gr_poly_evaluate` (Horner's rule),
declared in `gr_generic.h`.
:::

::: function
int _gr_poly_evaluate_other(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_poly_evaluate_other(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int _gr_poly_evaluate_other_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_poly_evaluate_other_horner(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int _gr_poly_evaluate_other_rectangular(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_poly_evaluate_other_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

Set `res` to $f(x)$ where the coefficients of $f$ belong to `ctx` and
$x$ and `res` belong to `x_ctx`, for example an integer polynomial
evaluated at a matrix or at a real ball. The coefficients enter the
computation through `gr_mul_other` and `gr_add_other`. The default
version uses Horner's rule.
:::

## Newton basis {#newton-basis}

Let $C = (c_0, c_1, c_2, \ldots)$ be a sequence of elements of a
commutative ring $R$. The polynomials

$$
N_0 = 1, \quad N_1 = x - c_0, \quad N_2 = (x - c_0)(x - c_1), \quad
\ldots, \quad N_k = \prod_{i < k} (x - c_i)
$$

form a basis of $R[x]$, the *Newton basis* associated with $C$. A
polynomial $f = g_0 N_0 + \cdots + g_n N_n$ in this basis is represented
by the formal polynomial $g = g_0 + g_1 x + \cdots + g_n x^n$, and a
finite initial segment of $C$ by a vector `basis`. The functions below
allow aliasing between inputs and outputs, except that no output may be
aliased with `basis`.

::: function
int _gr_poly_newton_basis_from_monomial(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_ctx_t ctx)
int gr_poly_newton_basis_from_monomial(gr_poly_t res, const gr_vec_t basis, const gr_poly_t poly, gr_ctx_t ctx)
int _gr_poly_newton_basis_to_monomial(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_ctx_t ctx)
int gr_poly_newton_basis_to_monomial(gr_poly_t res, const gr_vec_t basis, const gr_poly_t poly, gr_ctx_t ctx)

Convert a polynomial of length `len` from the monomial basis to the
Newton basis or back, with $O(\mathit{len}^2)$ operations. The
underscore versions read $\mathit{len} - 1$ entries of `basis`; the
non-underscore versions return `GR_UNABLE` if `basis` is shorter.
:::

::: function
int _gr_poly_newton_basis_evaluate(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_poly_newton_basis_evaluate(gr_ptr res, const gr_vec_t basis, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to the value at $x$ of the polynomial given in the Newton
basis, by a Horner scheme in the factors $x - c_i$. The same conditions
on `basis` apply as for the conversions.
:::

::: function
int _gr_poly_newton_basis_interpolate_exact(gr_ptr res, gr_srcptr basis, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_newton_basis_interpolate_exact(gr_poly_t res, const gr_vec_t basis, const gr_vec_t ys, gr_ctx_t ctx)
int _gr_poly_newton_basis_interpolate(gr_ptr res, gr_srcptr basis, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_newton_basis_interpolate(gr_poly_t res, const gr_vec_t basis, const gr_vec_t ys, gr_ctx_t ctx)

Set `res` to the Newton basis form of the interpolating polynomial $f$
of length at most `len` with $f(c_i) = y_i$ for $0 \le i <
\mathit{len}$, computed by divided differences. The first `len` points
of `basis` must be pairwise distinct.

::: details
The ring is presumed to be an integral domain, which is not checked;
over other rings the functions may succeed without the result being an
interpolating polynomial. They return `GR_DOMAIN` or `GR_UNABLE` if
`basis` has fewer than `len` points (checked by the non-underscore
versions only), if the points are not distinct, or if $f$ has
coefficients in the fraction field of $R$ but not in $R$. The *exact*
versions assume distinct points and coefficients in $R$, use exact
divisions, and may return an arbitrary polynomial otherwise.
:::
:::

## Multipoint evaluation and interpolation {#multipoint}

The fast algorithms use a subproduct tree: the products of $(x - x_i)$
over the leaves, pairs, quadruples and so on of the points, from which
evaluation is a descent of remainders and interpolation an ascent of
linear combinations, both with $O(M(n) \log n)$ operations; see
[](algorithms/polynomials/subproduct-tree). The tree is stored as an
array of levels, each a vector of the monic products of one level.

::: function
gr_ptr * _gr_poly_tree_alloc(slong len, gr_ctx_t ctx)
void _gr_poly_tree_free(gr_ptr * tree, slong len, gr_ctx_t ctx)
int _gr_poly_tree_build(gr_ptr * tree, gr_srcptr roots, slong len, gr_ctx_t ctx)

Allocate, free and fill the subproduct tree of `len` points `roots`:
level 0 holds the linear polynomials $x - x_i$ and each higher level the
products of pairs of the level below.
:::

::: function
int _gr_poly_product_roots(gr_ptr poly, gr_srcptr xs, slong n, gr_ctx_t ctx)
int gr_poly_product_roots(gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

Set `poly` to $(x - x_0)(x - x_1) \cdots (x - x_{n-1})$. The underscore
version writes $n + 1$ coefficients.
:::

::: function
int _gr_poly_evaluate_vec_iter(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)
int gr_poly_evaluate_vec_iter(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)
int _gr_poly_evaluate_vec_fast_precomp(gr_ptr vs, gr_srcptr poly, slong plen, const gr_ptr * tree, slong len, gr_ctx_t ctx)
int _gr_poly_evaluate_vec_fast(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)
int gr_poly_evaluate_vec_fast(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

Set $y_i = f(x_i)$ for all points. The *iter* versions evaluate at each
point separately with `_gr_poly_evaluate` ($O(n \cdot \mathit{plen})$
operations); the *fast* versions reduce $f$ down the subproduct tree,
which requires the divisions by monic tree polynomials. The *precomp*
version takes a tree built with `_gr_poly_tree_build`. The
non-underscore versions resize `ys` to the length of `xs`.
:::

::: function
int _gr_poly_interpolate(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_interpolate(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
int _gr_poly_interpolate_exact(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_interpolate_exact(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
int _gr_poly_interpolate_newton(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_interpolate_newton(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
int _gr_poly_interpolate_exact_newton(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_interpolate_exact_newton(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)

Set `res` to the interpolating polynomial $f$ of length at most `len`
with $f(x_i) = y_i$, for pairwise distinct points $x_i$. The default
versions currently call the *newton* versions, which compute divided
differences in the Newton basis and convert to the monomial basis, with
$O(\mathit{len}^2)$ operations.

::: details
The ring is presumed to be an integral domain; over other rings the
functions may succeed without the result interpolating. They return
`GR_DOMAIN` or `GR_UNABLE` if `xs` and `ys` have different lengths
(checked by the non-underscore versions only), if the points are not
distinct, or if $f$ has coefficients in the fraction field of $R$ but
not in $R$ (as for the points $0, 2$ and values $0, 1$ over $\Z$). The
*exact* versions assume that the points are distinct and that $f$ has
coefficients in $R$, use exact divisions, and may return an arbitrary
polynomial otherwise.
:::
:::

::: function
int _gr_poly_interpolation_weights(gr_ptr w, const gr_ptr * tree, slong len, gr_ctx_t ctx)
int _gr_poly_interpolate_fast_precomp(gr_ptr poly, gr_srcptr ys, const gr_ptr * tree, gr_srcptr weights, slong len, gr_ctx_t ctx)
int _gr_poly_interpolate_fast(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
int gr_poly_interpolate_fast(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)

Fast interpolation with a subproduct tree. The weights are $w_i = 1 /
\prod_{j \ne i} (x_i - x_j)$, obtained by evaluating the derivative of
the root product at the points; the *precomp* version takes the tree
from `_gr_poly_tree_build` and the weights from
`_gr_poly_interpolation_weights`. This requires a field (the weights
must be invertible).
:::

## Composition and Taylor shifts {#composition}

::: function
int _gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int _gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int _gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
int _gr_poly_taylor_shift_convolution(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_poly_taylor_shift_convolution(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

Set `res` to the Taylor shift $f(x + c)$. The *horner* version uses an
optimised Horner scheme with $O(\mathit{len}^2)$ additions and
multiplications by $c$; the *divconquer* version splits $f$ in halves
and multiplies by powers of $x + c$, in $O(M(\mathit{len}) \log
\mathit{len})$ operations; the *convolution* version uses one
multiplication of length `len`, $f(x + c) = \sum_k \frac{x^k}{k!}
\sum_i (i! f_i) \frac{c^{i-k}}{(i-k)!}$, which requires division by
$1, \ldots, \mathit{len} - 1$. The underscore versions allow aliasing.

::: details
`gr_poly_taylor_shift` calls the ring method
`GR_METHOD_POLY_TAYLOR_SHIFT`, whose default uses Horner's scheme for
$\mathit{len} \le 20$ and divide-and-conquer otherwise; `fmpz` uses
`_fmpz_poly_taylor_shift`, and `arb` and `acb` the Taylor shifts of
`arb_poly` and `acb_poly`. See
[](algorithms/polynomials/composition#taylor-shift).
:::
:::

::: function
int _gr_poly_compose(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_compose(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
int _gr_poly_compose_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_compose_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
int _gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_compose_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

Set `res` to the composition $f(g(x))$, where $f$ is `poly1` and $g$ is
`poly2`, by Horner's rule or by divide and conquer (splitting $f$ in
halves and multiplying by powers of $g$ computed by repeated squaring).
The default version handles $g = a x^n + c$ by a Taylor shift by $c$
followed by a scaling and an inflation, uses Horner's rule for
$\mathit{len1} \le 7$ and divide and conquer otherwise. The underscore
versions write $(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$
coefficients and do not allow aliasing.
:::

::: function
int _gr_poly_inflate(gr_ptr poly, slong len, slong n, gr_ctx_t ctx)
int gr_poly_inflate(gr_poly_t res, const gr_poly_t poly, ulong inflation, gr_ctx_t ctx)

Set `res` to $f(x^n)$. The underscore version works in place on an
array of $(\mathit{len} - 1) n + 1$ elements whose first `len` entries
hold $f$ and whose other entries are zero.
:::

::: function
ulong gr_poly_deflation(const gr_poly_t poly, gr_ctx_t ctx)
int gr_poly_deflate(gr_poly_t res, const gr_poly_t poly, ulong deflation, gr_ctx_t ctx)

`gr_poly_deflation` returns the largest $m$ such that $f = g(x^m)$ for a
polynomial $g$: 0 for $f = 0$, 1 for a nonzero constant, and $n$ for a
monomial $c x^n$ with $n \ge 1$. A coefficient that is not provably zero
counts as nonzero, so the value is always valid but need not be maximal
over inexact rings. `gr_poly_deflate` sets `res` to $g$ given $m \ge 1$
(`GR_DOMAIN` for $m = 0$); $f$ is assumed to have this form.
:::

## Power series composition and reversion {#series-composition}

::: function
int _gr_poly_compose_series(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_compose_series(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
int _gr_poly_compose_series_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_compose_series_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
int _gr_poly_compose_series_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_compose_series_brent_kung(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
int _gr_poly_compose_series_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_compose_series_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
int _gr_poly_compose_series_kinoshita_li(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int gr_poly_compose_series_kinoshita_li(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)

Set `res` to the power series composition $f(g(x)) \bmod x^n$, where
$f$ is `poly1` and $g$ is `poly2` with $g(0) = 0$ (the non-underscore
versions return `GR_DOMAIN` if $g(0)$ is nonzero and `GR_UNABLE` if it
is not provably zero). The algorithms are Horner's rule, the baby-step
giant-step method of Brent and Kung [@BrentKung1978],
divide and conquer, and the algorithm of Kinoshita and Li [@KL2024]
with $O(M(n) \log n)$ operations.

::: details
The default version handles $n \le 3$ and a monomial $g = a x^k$
directly, uses Horner's rule when $\mathit{len1} \le 7$ or $n \le 7$,
Brent-Kung when $\mathit{len1}^2 < n$ or the full composition is shorter
than $4n$, and otherwise divide and conquer below a cutoff and
Kinoshita-Li above it. The cutoff is 500, and for rings with a real
precision $p$ (balls, floating-point numbers) it is $\max(100, 15000 /
\sqrt{\max(64, p)})$. See
[](algorithms/polynomials/composition#series). The underscore versions
require $g(0) = 0$, do not allow aliasing of the output with an input,
and write $n$ coefficients without zero-padding beyond the length of
the composition.
:::
:::

::: function
int _gr_poly_revert_series(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
int gr_poly_revert_series(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
int _gr_poly_revert_series_lagrange(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
int gr_poly_revert_series_lagrange(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
int _gr_poly_revert_series_lagrange_fast(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
int gr_poly_revert_series_lagrange_fast(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
int _gr_poly_revert_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
int gr_poly_revert_series_newton(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)

Set `res` to the compositional inverse $f^{-1}$, with $f^{-1}(f(x)) =
f(f^{-1}(x)) = x \bmod x^n$. The reversion exists when $f(0) = 0$ and
$f'(0)$ is invertible; `GR_DOMAIN` is returned otherwise (`GR_UNABLE`
if $f(0)$ is not provably zero).

::: details
The *lagrange* version uses the Lagrange inversion formula, the
*lagrange_fast* version the baby-step giant-step form of it from
[@Joh2015b], and the *newton* version Newton iteration on $f(h) = x$
with power series compositions [@BrentKung1978]. The Lagrange versions
divide by $2, \ldots, n - 1$ and return `GR_UNABLE` in too small
characteristic. The default version uses Newton iteration from five
times the composition cutoff of `_gr_poly_compose_series` on (when
Kinoshita-Li composition makes it faster), and *lagrange_fast*
otherwise, falling back to Newton iteration when the latter fails with
`GR_UNABLE`. The underscore versions do not allow aliasing.
:::
:::

## Derivatives and integrals {#derivative}

::: function
int _gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
int gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

Sets `res` to the derivative $f'$. The underscore version writes
$\mathit{len} - 1$ coefficients and allows aliasing.
:::

::: function
int _gr_poly_nth_derivative(gr_ptr res, gr_srcptr poly, ulong n, slong len, gr_ctx_t ctx)
int gr_poly_nth_derivative(gr_poly_t res, const gr_poly_t poly, ulong n, gr_ctx_t ctx)

Sets `res` to the $n$-th derivative $f^{(n)}$. The underscore version
requires $n < \mathit{len}$, writes $\mathit{len} - n$ coefficients and
allows aliasing.
:::

::: function
int _gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
int gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

Sets `res` to the integral $\int_0^x f(t)\,dt$, whose coefficient of
$x^{i+1}$ is $f_i / (i + 1)$, computed with `gr_div_ui`; the divisions
must be possible in $R$ (`GR_DOMAIN` or `GR_UNABLE` otherwise, for
example over $\Z$ or in small characteristic). The underscore version
reads $\mathit{len} - 1$ coefficients of `poly`, writes `len`
coefficients and allows aliasing.
:::

## Monic polynomials and associates {#monic}

::: function
int _gr_poly_make_monic(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
int gr_poly_make_monic(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

Sets `res` to $f / \operatorname{lc}(f)$, multiplying by the inverse
of the leading coefficient or, if it has none, dividing each
coefficient by it (which can succeed over integral domains when the
divisions are exact). Returns `GR_DOMAIN` if $f = 0$ or a division
fails, and `GR_UNABLE` if the leading coefficient is not provably
nonzero.
:::

::: function
truth_t _gr_poly_is_monic(gr_srcptr poly, slong len, gr_ctx_t ctx)
truth_t gr_poly_is_monic(const gr_poly_t res, gr_ctx_t ctx)

Return whether the leading coefficient is 1. The zero polynomial is not
monic.
:::

::: function
int gr_poly_canonical_associate(gr_poly_t res, gr_poly_t u, const gr_poly_t src, gr_ctx_t ctx)

Sets `res` to the canonical associate of $f$: $f$ multiplied by the unit
$c$ of $R$ that makes the leading coefficient canonical in the sense of
`gr_canonical_associate` (1 over a field, positive over $\Z$). If `u` is
not `NULL`, it is set to the constant polynomial $c$, so that
$\mathit{res} = u f$. The zero polynomial is its own canonical associate,
with $u = 1$. GCD functions use this to normalise their output.
:::

## Greatest common divisors {#gcd}

The GCD functions choose the algorithm from the properties of the ring:
the Euclidean algorithm and the half-gcd algorithm over fields, and the
subresultant pseudo-remainder sequence over unique factorisation
domains. The algorithms and their specialisations are compared in
[](algorithms/polynomials/gcd).

::: function
int _gr_poly_gcd(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_gcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_gcd_generic(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int _gr_poly_gcd_euclidean(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_gcd_euclidean(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_gcd_hgcd(gr_ptr G, slong * _lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
int gr_poly_gcd_hgcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
int _gr_poly_gcd_subresultant(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_gcd_subresultant(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

Set $G$ to a greatest common divisor of $A$ and $B$. The non-underscore
functions return the canonical associate (the monic GCD over a field,
the GCD with positive leading coefficient over $\Z$), the zero
polynomial if both inputs are zero, and `GR_UNABLE` if a leading
coefficient is not provably nonzero. See [](gr_poly/gcd).
:::

::: function
int _gr_poly_xgcd(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_xgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_xgcd_generic(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int _gr_poly_xgcd_euclidean(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_xgcd_euclidean(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int _gr_poly_xgcd_hgcd(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)
int gr_poly_xgcd_hgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)
int _gr_poly_xgcd_subresultant(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int gr_poly_xgcd_subresultant(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

Set $G$, $S$, $T$ with $G = SA + TB$, where $G$ is a GCD of $A$ and $B$
over a field (monic) and over other GCD domains an associate of the GCD
that may differ from it by a factor dividing the resultant. See
[](gr_poly/gcd).
:::

::: function
int gr_poly_gcd_wrapper(gr_method_poly_gcd_op gcd_impl, int canonicalise_unit, gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
int gr_poly_xgcd_wrapper(gr_method_poly_xgcd_op _xgcd_op, gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

Implement `gr_poly_gcd` and `gr_poly_xgcd` around a given underscore
implementation (with the signatures of `_gr_poly_gcd` and
`_gr_poly_xgcd`): handle zero and constant inputs, order the inputs by
length, check the leading coefficients, manage memory and aliasing, and
canonicalise the result (for the GCD only if `canonicalise_unit` is
nonzero; the cofactors of the extended GCD are scaled by the same
unit). The non-underscore algorithm-specific functions are built this
way, and rings can use them to wrap their own GCD methods.
:::

::: function
int _gr_poly_hgcd(gr_ptr r, slong * sgn, gr_ptr * M, slong * lenM, gr_ptr A, slong * lenA, gr_ptr B, slong * lenB, gr_srcptr a, slong lena, gr_srcptr b, slong lenb, slong cutoff, gr_ctx_t ctx)

Computes the half-gcd of $a$ and $b$ ([](algorithms/polynomials/half-gcd)):
a $2 \times 2$ matrix $M$ of polynomials, a sign $\sigma = \pm 1 = \det
M$ and the consecutive remainders $A$, $B$ with $(A, B)^T = \sigma
M^{-1} (a, b)^T$ and $\deg A \ge \lceil \deg a / 2 \rceil > \deg B$. The
recursion uses the Euclidean algorithm below `cutoff`.

::: details
Requires $\mathit{lena} > \mathit{lenb} > 0$ and nonzero leading
coefficients; returns `GR_DOMAIN` if a length is zero or a leading
coefficient is not invertible. $A$ and $B$ must have room for
$\mathit{lena}$ and $\mathit{lenb}$ coefficients; on return `*lenA` and
`*lenB` hold their lengths. `M[0]`, …, `M[3]` must each have room for
$\mathit{lena}$ coefficients, with the lengths written to `lenM`; if `M`
is `NULL`, only $A$ and $B$ are computed. If `r` is not `NULL`, it is
an accumulator for the resultant: it is multiplied by the factor
contributed by the remainder sequence up to $(A, B)$, as used by
`_gr_poly_resultant_hgcd`. The input pairs `(A, a)` and `(B, b)` may be
aliased.

The products of the $2 \times 2$ polynomial matrices are done by
`_gr_poly_hgcd_mat_mul`, a ring method; the default
`_gr_poly_hgcd_mat_mul_generic` uses 8 classical products when some
entry has length below 20 and Strassen's 7 products otherwise, and when
the ring provides a transformed representation (with `fft_small`) and
all entries have length at least 128, it transforms each of the 8
entries once and computes the products pointwise, in several threads
for entries of length at least 1024.
:::
:::

## Resultants {#resultant}

For nonzero polynomials $f = a_m x^m + \cdots + a_0$ and $g = b_n x^n +
\cdots + b_0$ of degrees $m$ and $n$, the resultant is

$$
\operatorname{res}(f, g) = a_m^n b_n^m \prod_{f(\alpha) = 0,\, g(\beta) = 0} (\alpha - \beta)
= \det \operatorname{Syl}(f, g),
$$

the determinant of the Sylvester matrix. The resultant is defined to be
zero if $f$ or $g$ is zero.

::: function
int _gr_poly_resultant(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_resultant(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int _gr_poly_resultant_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int _gr_poly_resultant_euclidean(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_resultant_euclidean(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int _gr_poly_resultant_hgcd(gr_ptr res, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
int gr_poly_resultant_hgcd(gr_ptr res, const gr_poly_t f, const gr_poly_t g, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
int _gr_poly_resultant_subresultant(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_resultant_subresultant(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int _gr_poly_resultant_sylvester(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_resultant_sylvester(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int _gr_poly_resultant_small(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_resultant_small(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)

Set `res` to the resultant of $f$ and $g$. The non-underscore functions
return `GR_UNABLE` if a leading coefficient is not provably nonzero.
The underscore functions require $\mathit{len1} \ge \mathit{len2} \ge
1$ and nonzero leading coefficients.

::: details
- The *euclidean* version follows the Euclidean remainder sequence and
  accumulates the leading coefficients; the *hgcd* version does the same
  with the half-gcd algorithm, with the tuning parameters
  `inner_cutoff` (passed to `_gr_poly_hgcd`) and `cutoff`. Both need
  inverses of leading coefficients and return `GR_DOMAIN` when one does
  not exist, so they are meant for fields.
- The *subresultant* version uses the subresultant pseudo-remainder
  sequence [@Col1967; @BT1971], with exact divisions only; it is valid
  over any GCD domain (in fact any integral domain with exact
  division).
- The *sylvester* version computes the determinant of the Sylvester
  matrix with `gr_mat_det`, which is division-free by default and so
  works over any commutative ring, and is also the most robust choice
  over inexact rings.
- The *small* version uses division-free straight-line formulas for
  $\mathit{len2} \le 3$ or $\mathit{len1} \le 2$ and returns `GR_UNABLE`
  for larger inputs.

`gr_poly_resultant` calls the ring method `GR_METHOD_POLY_RESULTANT`.
The default `_gr_poly_resultant_generic` tries the *small* formulas;
over finite rings and fields it then uses the half-gcd algorithm when
both lengths are at least 200 (inner cutoff 100) and the Euclidean
algorithm otherwise; if that fails or the ring is not a field, it uses
the subresultant algorithm over unique factorisation domains, and
finally the Sylvester determinant. The contexts of polynomial rings
override the method with the bivariate algorithms below.
:::
:::

::: function
int _gr_poly_resultant_multipoint(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
int gr_poly_resultant_multipoint(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int _gr_poly_resultant_modular(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, int proved, gr_ctx_t ctx)
int gr_poly_resultant_modular(gr_ptr res, const gr_poly_t f, const gr_poly_t g, int proved, gr_ctx_t ctx)

Resultants of bivariate polynomials, that is of polynomials whose
coefficients are polynomials in another variable $y$ (the ring `ctx` is
a polynomial ring), with respect to $x$. They return `GR_UNABLE` for
other rings. The algorithms work on the dense bivariate types:
`n_bpoly_mod_resultant` over a word-size prime field and
`fmpz_bpoly_resultant` over $\Z$ and $\mathbb{Q}$.

::: details
The *multipoint* version, for coefficients in $(\Z/p\Z)[y]$ with $p$ a
word-size prime, evaluates the coefficients at enough points to
determine the result (the degree bound in $y$ plus one), computes the
univariate resultants at these points and interpolates
([](algorithms/multivariate/modp-bivariate-resultant)). With an
FFT-friendly modulus the evaluation and interpolation use a DFT instead
of geometric progressions; the evaluations and the univariate
resultants are split over the thread pool when there is enough work.

The *modular* version, for coefficients in $\Z[y]$ or $\mathbb{Q}[y]$ (rational
inputs are cleared of denominators), runs the multipoint algorithm
modulo several word-size primes and combines the results by Chinese
remaindering (`fmpz_bpoly_resultant`). If `proved` is 1, enough primes are used to cover a bound
for the result, so the result is certain. Otherwise the iteration stops
once the reconstruction has not changed for a few consecutive primes,
which can be much faster but makes the result correct only with high
probability (a Monte Carlo algorithm).

The contexts created by `gr_ctx_init_gr_poly` (over `nmod`, `fmpz` or
`fmpq`) and by `gr_ctx_init_fmpz_poly` and `gr_ctx_init_fmpq_poly`
supply these algorithms through `GR_METHOD_POLY_RESULTANT`: the
multipoint algorithm when `_gr_poly_resultant_multipoint_cutoff`
accepts the sizes, and the modular algorithm with `proved` set when
the shorter input has at least 7 coefficients or the resultant has
degree at least 96 in $y$; below these sizes the subresultant
algorithm is faster.
:::
:::

::: function
int _gr_poly_resultant_multipoint_cutoff(slong len1, slong len2, slong npoints)

Returns whether the multipoint algorithm is expected to be faster than
the subresultant algorithm for bivariate inputs of the given lengths in
$x$, where `npoints` is the number of evaluation points needed (the
degree bound in $y$ plus one). This is a tuning function used by the
polynomial-ring contexts.
:::

## Squarefree factorisation {#squarefree}

::: function
int gr_poly_factor_squarefree(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp, const gr_poly_t poly, gr_ctx_t ctx)

Computes a squarefree factorisation $f = c \prod_i g_i^{e_i}$ with
pairwise coprime squarefree $g_i$ and distinct multiplicities $e_i \ge
1$, using Yun's algorithm in characteristic zero and its extension
with $p$-th roots over finite fields ([](algorithms/polynomials/squarefree)).
The factors are written to `fac`, a vector of polynomials over `ctx`,
and the multiplicities to `exp`, in no particular order.

::: details
The constant $c$ is an element of $R$. Over a field all factors are
monic and $c$ is the leading coefficient of $f$. Over other rings the
factors are primitive with canonical leading coefficients, and $c$ is a
unit multiple of the content of $f$ (which is not factored further).
The ring must be a unique factorisation domain, and in positive
characteristic a finite field (for the $p$-th roots); `GR_UNABLE` is
returned if this cannot be verified.
:::
:::

::: function
int gr_poly_squarefree_part(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

Sets `res` to the squarefree part of $f$, the canonical associate of
$\prod_i g_i$ in terms of the squarefree factorisation. In
characteristic zero this is $f / \gcd(f, f')$.
:::

::: function
truth_t gr_poly_is_squarefree(const gr_poly_t f, gr_ctx_t ctx)

Returns whether $f$ is squarefree. Unless the ring is known to be a
field, the answer is `T_UNKNOWN`. Over a field the test is $\gcd(f, f')
= 1$; in positive characteristic a nonconstant $f$ with $f' = 0$ is a
$p$-th power and hence not squarefree if the field is finite, and
`T_UNKNOWN` is returned if finiteness cannot be verified. Nonzero
constants are squarefree and zero is not.
:::

## Factorisation over finite fields {#factor-fq}

The functions in this section assume that `ctx` is a finite field
$\mathbb{F}_q$, $q = p^k$, and return `GR_UNABLE` unless
`gr_ctx_is_field`, `gr_ctx_is_finite` and
`gr_ctx_is_finite_characteristic` all return `T_TRUE` and
`gr_ctx_fq_order`, `gr_ctx_fq_prime` and `gr_ctx_fq_degree` are
implemented. This covers `nmod`, `fmpz_mod` and `mpn_mod` when the
modulus is known to be prime (see `gr_ctx_set_is_field`), and `fq`,
`fq_nmod` and `fq_zech`. For these rings `gr_factor` on elements of the
polynomial ring (`gr_ctx_init_gr_poly`) and the default
`GR_METHOD_POLY_FACTOR` use `gr_poly_factor_finite_field`. The
algorithms are explained in [](algorithms/polynomials/factoring-finite-fields)
and [](algorithms/polynomials/irreducibility).

::: function
truth_t gr_poly_is_irreducible(const gr_poly_t f, gr_ctx_t ctx)
truth_t gr_poly_is_irreducible_ben_or(const gr_poly_t f, gr_ctx_t ctx)
truth_t gr_poly_is_irreducible_ddf(const gr_poly_t f, gr_ctx_t ctx)
truth_t gr_poly_is_irreducible_rabin(const gr_poly_t f, gr_ctx_t ctx)

Return whether $f$ is irreducible over $\mathbb{F}_q$. Constants are
not irreducible.

::: details
The *ben_or* version uses Ben-Or's test [@BenOr1981], computing
$\gcd(f, x^{q^i} - x)$ for $i \le n/2$ with the Frobenius iterates
obtained by modular composition or exponentiation. The *ddf* version
runs the coarse phase of the baby-step giant-step distinct-degree
factorisation (`gr_poly_factor_distinct_deg`) and stops at the first
factor found, which is fast on inputs with small factors. The *rabin*
version uses Rabin's test [@Rab1980b]: $f$ of degree $n$ is irreducible
if and only if $x^{q^n} \equiv x \bmod f$ and $\gcd(x^{q^{n/r}} - x, f) =
1$ for every prime $r \mid n$; the Frobenius powers are computed by
modular composition from the ladder $x^{q^{2^i}}$, so that the test
costs $O((1 + \omega(n)) \log n)$ compositions after the initial $x^q$,
against about $2 \sqrt{n}$ for the distinct-degree test, but it does not
stop early on inputs with small factors.

The default version first looks for roots by evaluation when the field
is a small prime field; for degree below 600 it then calls the *ddf*
version, and from 600 on it checks squarefreeness, removes inputs with
factors of degree at most 4 with a few gcds, and calls the *rabin*
version.
:::
:::

::: function
int _gr_poly_iterated_frobenius_preinv(gr_poly_struct * h, slong n, const gr_poly_t frob, const gr_poly_preinv_t P, const fmpz_t q, gr_ctx_t ctx)
int gr_poly_iterated_frobenius_preinv(gr_poly_struct * h, slong n, const gr_poly_preinv_t P, gr_ctx_t ctx)

Set $h_i = x^{q^i} \bmod v$ for $0 \le i < n$, where `P` holds the
modulus $v$. The underscore version takes $x^q \bmod v$ as `frob`.
Depending on the size of $q$ relative to the degree of $v$, the powers
are computed by repeated exponentiation or by modular composition with
a precomputed matrix of powers of `frob`.
:::

::: function
int gr_poly_factor_distinct_deg(gr_poly_vec_t fac, fmpz_vec_t degs, const gr_poly_t poly, gr_ctx_t ctx)
int _gr_poly_factor_distinct_deg_with_frob(gr_poly_vec_t fac, fmpz_vec_t degs, const gr_poly_t poly, const gr_poly_preinv_t P, const gr_poly_t frob, const fmpz_t q, gr_ctx_t ctx)
int _gr_poly_factor_distinct_deg_with_frob_threaded(gr_poly_vec_t fac, fmpz_vec_t degs, const gr_poly_t poly, const gr_poly_preinv_t P, const gr_poly_t frob, const fmpz_t q, gr_ctx_t ctx)

Compute the distinct-degree factorisation of a squarefree polynomial:
$\operatorname{monic}(f) = \prod_i g_i$, where $g_i$ is the product of
the irreducible factors of degree $d_i$ and the $d_i$ are distinct. The
$g_i$ are written to `fac` and the $d_i$ to `degs`. The result is
undefined if $f$ is not squarefree.

::: details
This is the baby-step giant-step algorithm of Kaltofen and Shoup
[@KalSho1998], with the baby steps computed by
`_gr_poly_iterated_frobenius_preinv` and the giant steps by modular
composition with a precomputed matrix. The underscore versions assume
$f$ monic, take the precomputed modulus `P` for $f$, $x^q \bmod f$ as
`frob` and $q$, and append to `fac` and `degs`. The factorisation
functions choose the sparse representation of the modulus
automatically when it applies (`gr_poly_preinv_set`), so sparse inputs
benefit from cheap reductions.

When more than one thread is available (`flint_set_num_threads`) and the
ring is thread-safe (`gr_ctx_is_threadsafe`), the threaded variant is
used: the giant steps are computed in blocks of independent
compositions sharing one precomputed matrix, the interval polynomials
of a block and the splitting of the intervals are computed in parallel,
and the vector Brent-Kung composition of the baby steps evaluates its
polynomials in parallel. Matrix products are done outside the parallel
sections, so that `gr_mat_mul` can use threads itself. The total work
exceeds that of the serial algorithm by a few percent.
:::
:::

::: function
int gr_poly_factor_equal_deg_prob(gr_poly_t factor, flint_rand_t state, const gr_poly_t pol, slong d, gr_ctx_t ctx)

One random trial of equal-degree factorisation: for a squarefree $f$ of
degree $n = rd$, $r \ge 2$, all of whose irreducible factors have
degree $d$, tries to find a proper monic factor, written to `factor`
(set to 1 on failure). A trial succeeds with probability roughly $1/2$
(at least $4/9$). Returns `GR_DOMAIN` if $d \nmid n$ or $n < 2$.
:::

::: function
int gr_poly_factor_equal_deg(gr_poly_vec_t fac, const gr_poly_t pol, slong d, gr_ctx_t ctx)
int _gr_poly_factor_equal_deg_with_frob(gr_poly_vec_t fac, const gr_poly_t f, slong d, const gr_poly_t frob, flint_rand_t state, gr_ctx_t ctx)

Equal-degree factorisation: for a squarefree $f$ whose irreducible
factors all have degree $d$, set `fac` to the monic irreducible
factors. The result is undefined for other inputs.

::: details
This is the Cantor-Zassenhaus algorithm [@CanZas1981] in the trace
form of von zur Gathen and Shoup [@GatSho1992]: for random $a \in
\mathbb{F}_q[x]/(f)$, the trace $b = a + a^q + \cdots + a^{q^{d-1}}$
(computed by modular composition with the Frobenius $x^q$, by
exponentiation, or by a doubling scheme, depending on $q$ and $d$) lies
in the Berlekamp subalgebra, and $\gcd(f, b^{(q-1)/2} - 1)$ for odd $q$,
or $\gcd(f, b + b^2 + \cdots + b^{2^{k-1}})$ for $q = 2^k$, splits $f$
with probability about $1/2$. The splitting is repeated on the factors.
For $d = 1$ (root finding) $a$ is a random linear polynomial, which
makes the exponentiation cheaper. When several pieces remain and
threads are available over a thread-safe ring, the pieces are split in
parallel, each with its own random state. The underscore version
assumes $f$ monic, takes $x^q \bmod f$ as `frob` (unused for $d = 1$)
and a random state, and appends to `fac`.
:::
:::

::: function
int gr_poly_factor_finite_field(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp, const gr_poly_t F, gr_ctx_t ctx)
int _gr_poly_factor_finite_field(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp, const gr_poly_t F, int algorithm, gr_ctx_t ctx)
int gr_poly_factor_cantor_zassenhaus(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp, const gr_poly_t F, gr_ctx_t ctx)
int gr_poly_factor_berlekamp(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp, const gr_poly_t F, gr_ctx_t ctx)
int gr_poly_factor_kaltofen_shoup(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp, const gr_poly_t F, gr_ctx_t ctx)

Factor $F$ over $\mathbb{F}_q$ as $F = c \prod_i g_i^{e_i}$, where $c$ is
the leading coefficient, the $g_i$ are distinct monic irreducible
polynomials (written to `fac`, a vector of polynomials over `ctx`) and
the $e_i \ge 1$ are written to `exp`, in no particular order. For $F =
0$, $c$ is set to zero and the factorisation is empty.

::: details
All versions start with a squarefree factorisation and then factor
each squarefree part: with the Cantor-Zassenhaus algorithm
(distinct-degree factorisation by $\gcd(f, x^{q^i} - x)$ for $i = 1, 2,
\ldots$, then equal-degree factorisation), Berlekamp's algorithm
[@Ber1967; @Ber1970] (a basis of the kernel of the Frobenius map minus
the identity by linear algebra, then random splitting), or the
Kaltofen-Shoup algorithm (baby-step giant-step distinct-degree
factorisation, then equal-degree factorisation). The main function
first deflates inputs of the form $G(x^m)$, $m > 1$, and chooses
between Cantor-Zassenhaus and Kaltofen-Shoup from the degree and the
size of the field. The underscore version takes `algorithm`, one of
`GR_POLY_FACTOR_ALGORITHM_DEFAULT`,
`GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS`,
`GR_POLY_FACTOR_ALGORITHM_BERLEKAMP` and
`GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP`.
:::
:::

::: macro
GR_POLY_FACTOR_ALGORITHM_DEFAULT
GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS
GR_POLY_FACTOR_ALGORITHM_BERLEKAMP
GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP

The values 0 to 3 of the `algorithm` argument of
`_gr_poly_factor_finite_field`.
:::

::: function
int _gr_poly_factor_finite_field_method(gr_ptr c, gr_vec_t fac, fmpz_vec_t exp, const gr_poly_t F, int flags, gr_ctx_t ctx)

`gr_poly_factor_finite_field` with the signature of the
`GR_METHOD_POLY_FACTOR` method, for the method tables of finite field
contexts: $c$ is an element of the polynomial ring (a `gr_poly_t`) and
`fac` a vector of elements of the polynomial ring.
:::

::: function
int gr_poly_roots_finite_field(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)

`gr_poly_roots` for finite fields, usable as the `GR_METHOD_POLY_ROOTS`
method: computes the roots in $\mathbb{F}_q$ with multiplicities by a
squarefree factorisation, $\gcd(f, x^q - x)$ and equal-degree
factorisation with $d = 1$. See [](algorithms/polynomials/root-finding).
:::

## Shift equivalence and dispersion {#dispersion}

::: function
truth_t gr_poly_shift_equivalent(fmpz_t shift, const gr_poly_t p, const gr_poly_t q, gr_ctx_t ctx)

Returns whether there is an integer $n$ with $p(x + n) = q(x)$. If the
answer is `T_TRUE` and `shift` is not `NULL`, it is set to such an $n$.
:::

::: function
int gr_poly_leading_taylor_shift(gr_ptr shift, const gr_poly_t p, const gr_poly_t q, gr_ctx_t ctx)

Computes the $s$ for which $p(x + s)$ and $q(x)$ have the same two
leading coefficients: if $p$ and $q$
have degree $n$ and the same leading coefficient $c$, then $s =
(q_{n-1} - p_{n-1}) / (n c)$. Returns `GR_DOMAIN` if the degrees or the
leading coefficients differ or the division fails, and `GR_UNABLE` if a
leading coefficient is not provably nonzero. This is the candidate
tested by `gr_poly_shift_equivalent`.
:::

::: function
int gr_poly_dispersion(fmpz_t disp, gr_vec_t disp_set, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int gr_poly_dispersion_factor(fmpz_t disp, gr_vec_t disp_set, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
int gr_poly_dispersion_resultant(fmpz_t disp, gr_vec_t disp_set, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)

Compute the dispersion set of $f$ and $g$, over a unique factorisation
domain of characteristic zero: the set of integers $n \ge 0$ such that
$f(x + n)$ and $g(x)$ have a nonconstant common factor, and the
dispersion, its largest element.

::: details
Either output may be `NULL`. If the set is empty, `disp` is unchanged.
The elements of `disp_set` (of type `fmpz`) are sorted in increasing
order. The *factor* version factors $f$ and $g$ and compares the
factors as in [@ManWright1994]; the *resultant* version finds the
nonnegative integer roots of $\operatorname{res}_x(f(x + y), g(x))$ and
is mainly meant for testing.
:::
:::

::: function
int gr_poly_dispersion_from_factors(fmpz_t disp, gr_vec_t disp_set, const gr_poly_vec_t ffac, const gr_poly_vec_t gfac, gr_ctx_t ctx)

`gr_poly_dispersion_factor` for nonzero $f$ and $g$, given their
nonconstant irreducible factors (without multiplicities).
:::

::: function
int gr_poly_shiftless_decomposition(gr_ptr c, gr_poly_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_poly_t f, gr_ctx_t ctx)
int gr_poly_shiftless_decomposition_factor(gr_ptr c, gr_poly_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_poly_t f, gr_ctx_t ctx)

Compute a decomposition

$$
f = c \prod_{i=1}^{r} \prod_{j} g_i(x + h_{i,j})^{e_{i,j}}
$$

where $c$ is a constant, the $g_i$ are squarefree of degree at least 1,
$g_i(x)$ and $g_j(x + h)$ are coprime for $i \ne j$ and all $h \in \Z$,
$g_i(x)$ and $g_i(x + h)$ are coprime for all nonzero $h \in \Z$, and
the $e_{i,j} \ge 1$ and $0 = h_{i,1} < h_{i,2} < \cdots$ are integers.

::: details
`slfac` must be a vector of polynomials over `ctx`, and `slshifts` and
`slmult` vectors of vectors of `fmpz`. The *factor* version (the only
algorithm implemented) computes the irreducible factorisation and sorts
the factors into classes of shift-equivalent polynomials.
:::
:::

::: function
int _gr_poly_shiftless_decomposition_from_factors(gr_poly_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_poly_vec_t fac, const fmpz_vec_t mult, gr_ctx_t ctx)
int gr_poly_shiftless_decomposition_from_factors(gr_poly_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_poly_vec_t fac, const fmpz_vec_t mult, gr_ctx_t ctx)

The shiftless decomposition from a given irreducible factorisation
(`fac`, `mult`) of $f$, without the constant $c$. The underscore
version does not allow `slfac` to be aliased with `fac`.
:::

## Roots {#roots}

::: function
int gr_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
int gr_poly_roots_other(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, gr_ctx_t poly_ctx, int flags, gr_ctx_t ctx)

Find the roots of $f$ in the ring `ctx`, writing them without repetition
to `roots` (a vector over `ctx`) and their multiplicities to `mult`.
`gr_poly_roots_other` takes a polynomial over another ring `poly_ctx`,
for example to find the `qqbar` or `acb` roots of an integer
polynomial. The roots of the zero polynomial are undefined, and
`GR_DOMAIN` is returned for it.

::: details
Unless the ring is an algebraically closed field, the multiplicities
may add up to less than the degree: over `fmpz` only integer roots are
found. The order of the roots depends on the ring. Both functions are
ring methods (`GR_METHOD_POLY_ROOTS`, `GR_METHOD_POLY_ROOTS_OTHER`):
the default finds roots over finite fields with
`gr_poly_roots_finite_field` ([](algorithms/polynomials/root-finding))
and returns `GR_UNABLE` otherwise;
`fmpz`, `fmpq`, `qqbar`, `ca`, `arb`, `acb` and the modular rings
provide their own.
:::
:::

::: function
int _gr_poly_refine_roots_aberth(gr_ptr w, gr_srcptr f, gr_srcptr f_prime, slong deg, gr_srcptr z, int progressive, gr_ctx_t ctx)
int _gr_poly_refine_roots_wdk(gr_ptr w, gr_srcptr f, slong deg, gr_srcptr z, int progressive, gr_ctx_t ctx)

Given approximations $z_1, \ldots, z_{\mathit{deg}}$ of the complex roots
of $f = \sum_{i=0}^{\mathit{deg}} f_i x^i$, compute corrections $w_k$
such that $z_k - w_k$ approximate the roots better if the initial
approximations are close enough and the evaluation of $f$ is accurate.
These are meant to be called in a loop, over floating-point rings.

::: details
The *wdk* version performs the Weierstrass-Durand-Kerner update

$$
w_k = \frac{f(z_k)}{f_{\mathit{deg}} \prod_{j \ne k} (z_k - z_j)},
$$

and the *aberth* version the Aberth-Ehrlich update

$$
w_k = \frac{g_k}{1 - g_k \sum_{j \ne k} (z_k - z_j)^{-1}}, \qquad g_k = \frac{f(z_k)}{f'(z_k)},
$$

which needs the coefficients of $f'$ in `f_prime`. If `progressive` is
set, the already corrected values $z_j - w_j$ replace $z_j$ in the
updates (a Gauss-Seidel variant), which can speed up convergence.
:::
:::

## Elementary functions of power series {#special}

The functions in this section compute $F(h) \bmod x^n$ for a power
series $h$ and an elementary function $F$. The constant term $F(h_0)$
is computed by the corresponding scalar function (`gr_exp`, `gr_log`,
`gr_sin_cos` and so on) and must exist in $R$; the remaining
coefficients involve divisions by small integers, which must be
possible. The underscore versions write $n$ coefficients and allow
aliasing unless stated otherwise. See [](algorithms/polynomials/newton)
for the Newton iterations.

::: function
int _gr_poly_exp_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_exp_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_exp_series_generic(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int _gr_poly_exp_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_exp_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_exp_series_basecase_mul(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_exp_series_basecase_mul(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_exp_series_newton(gr_ptr f, gr_ptr g, gr_srcptr h, slong hlen, slong n, slong cutoff, gr_ctx_t ctx)
int gr_poly_exp_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx)

Set $f$ to $\exp(h) \bmod x^n$.

::: details
The *basecase* version uses the recurrence $k f_k = \sum_{j=1}^{k} j h_j
f_{k-j}$ from $f' = h' f$, with $O(\mathit{hlen} \cdot n)$ operations and
divisions by $1, \ldots, n - 1$. The *basecase_mul* version computes
the low third with the recurrence and corrects the rest with two
truncated products, using $\exp(a + b) = \exp(a)(1 + b + b^2/2 + \cdots)$
for $b$ divisible by $x^{\lceil n/3 \rceil}$. The *newton* version uses the
iteration $f \gets f(1 + h - \log f)$ with the simultaneous iteration
for $1/f$ of [@HZ2004; @Ber2008], starting from the basecase below
`cutoff`; if `g` is not `NULL`, the underscore version also writes $1/f
\bmod x^n$ to it.

`_gr_poly_exp_series` calls the ring method `GR_METHOD_POLY_EXP_SERIES`
(overridden by `nmod`), whose default `_gr_poly_exp_series_generic` uses
the basecase for $n < 10$, $\mathit{hlen} < 10$ or rings with classical
multiplication, the *basecase_mul* version for $n < 30$ when $h$ is
nearly as long as $n$, and Newton iteration with cutoff 30 otherwise.
:::
:::

::: function
int _gr_poly_log_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_log_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_log1p_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_log1p_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

Set `res` to $\log(f)$, respectively $\log(1 + f)$, modulo
$x^{\mathit{len}}$, computed as the integral of $f' / f$ (with one
power series division), with a direct series for binomials $a + b x^d$.
:::

::: function
int _gr_poly_asin_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_asin_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_asinh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_asinh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_acos_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_acos_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_acosh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_acosh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_atan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_atan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
int _gr_poly_atanh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int gr_poly_atanh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

Set `res` to the inverse trigonometric or hyperbolic function of $f$
modulo $x^{\mathit{len}}$, as the integral of $f' / \sqrt{1 - f^2}$,
$f' / \sqrt{1 + f^2}$, $-f' / \sqrt{1 - f^2}$, $f' / \sqrt{f^2 - 1}$,
$f' / (1 + f^2)$ or $f' / (1 - f^2)$ plus the value at the constant
term.
:::

::: function
int _gr_poly_sin_cos_series(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_sin_cos_series(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_sin_cos_pi_series(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_sin_cos_pi_series(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_sin_series(gr_ptr s, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_sin_series(gr_poly_t s, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_sin_pi_series(gr_ptr s, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_sin_pi_series(gr_poly_t s, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_cos_series(gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_cos_series(gr_poly_t c, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_cos_pi_series(gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_cos_pi_series(gr_poly_t c, const gr_poly_t h, slong n, gr_ctx_t ctx)

Set $s = \sin(h)$ and $c = \cos(h)$ modulo $x^n$, or $\sin(\pi h)$ and
$\cos(\pi h)$ for the `pi` versions. The default algorithm is the
basecase for $\mathit{hlen} < 20$ and below a ring-dependent cutoff (500
over finite rings, $\min(500, 20 + 100000/p)$ over rings with real
precision $p$, 25 otherwise), and Newton iteration above it.
:::

::: function
int _gr_poly_sin_cos_series_basecase(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, int times_pi, gr_ctx_t ctx)
int gr_poly_sin_cos_series_basecase(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
int _gr_poly_sin_cos_series_tangent(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, int times_pi, gr_ctx_t ctx)
int gr_poly_sin_cos_series_tangent(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
int _gr_poly_sin_cos_series_newton(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, slong cutoff, int times_pi, gr_ctx_t ctx)
int gr_poly_sin_cos_series_newton(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, slong cutoff, int times_pi, gr_ctx_t ctx)

Compute the sine and cosine of a power series with a specific
algorithm; `times_pi` multiplies the argument by $\pi$.

::: details
The *basecase* version uses the coupled recurrences from $s' = h' c$,
$c' = -h' s$, with $O(\mathit{hlen} \cdot n)$ operations. The *tangent*
version computes $t = \tan(h_1 / 2)$ for the nonconstant part $h_1$
with `_gr_poly_tan_series` and uses the half-angle formulas $\sin h_1 =
2t / (1 + t^2)$, $\cos h_1 = (1 - t^2) / (1 + t^2)$, with $O(M(n))$
operations, and then the addition theorems for the constant term
$h_0$. The *newton* version runs Newton iteration for $\exp(i h)$ with
the real and imaginary parts represented formally, so that no complex
arithmetic is needed; `cutoff` is the length computed by the basecase.
:::
:::

::: function
int _gr_poly_tan_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_tan_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_tanh_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_tanh_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_cot_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_cot_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_coth_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_coth_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_tan_pi_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_tan_pi_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
int _gr_poly_cot_pi_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
int gr_poly_cot_pi_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)

Set $f$ to $\tan$, $\tanh$, $\cot$, $\coth$, $\tan(\pi \cdot)$ or
$\cot(\pi \cdot)$ of $h$ modulo $x^n$. Over rings with real precision,
an argument with large imaginary part (for $\tan$, $\cot$) or large
real part (for $\tanh$, $\coth$) is handled by the *exponential*
algorithm for the sake of relative accuracy; otherwise the basecase is
used below a cutoff as for the sine and cosine, and Newton iteration
above it.
:::

::: function
int _gr_poly_tan_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, int func, gr_ctx_t ctx)
int gr_poly_tan_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, int func, gr_ctx_t ctx)
int _gr_poly_tan_series_newton(gr_ptr f, gr_srcptr h, slong hlen, slong n, slong cutoff, int func, gr_ctx_t ctx)
int gr_poly_tan_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, int func, gr_ctx_t ctx)
int _gr_poly_tan_series_sine_cosine(gr_ptr f, gr_srcptr h, slong hlen, slong n, int func, gr_ctx_t ctx)
int gr_poly_tan_series_sine_cosine(gr_poly_t res, const gr_poly_t h, slong len, int func, gr_ctx_t ctx)
int _gr_poly_tan_series_exponential(gr_ptr f, gr_srcptr h, slong hlen, slong n, int func, gr_ctx_t ctx)
int gr_poly_tan_series_exponential(gr_poly_t res, const gr_poly_t h, slong len, int func, gr_ctx_t ctx)

Compute a tangent-like function of a power series with a specific
algorithm. The parameter `func` selects the function: 0 for $\tan$, 1
for $\tanh$, 2 for $\cot$, 3 for $\coth$, 4 for $\tan(\pi \cdot)$, 6 for
$\cot(\pi \cdot)$ (5 and 7, the hyperbolic functions times $\pi$, are
not implemented).

::: details
The *basecase* algorithm uses an $O(n^2)$ recurrence for the
coefficients (from $f' = h'(1 \pm f^2)$). The *newton* algorithm inverts
the integral defining the inverse function by Newton iteration, with
the basecase below `cutoff`. The *sine_cosine* algorithm divides a sine
by a cosine; for the hyperbolic functions it needs the imaginary unit
in the ring. The *exponential* algorithm writes the function as a
quotient of $\pm 1$ plus an exponential, choosing the sign of the
exponent (if possible) so that the exponential is small; for the
trigonometric functions it needs the imaginary unit.
:::
:::

::: function
int _gr_poly_bessel_j_series(gr_ptr res, gr_srcptr nu, gr_srcptr z, slong zlen, slong len, gr_ctx_t ctx)
int gr_poly_bessel_j_series(gr_poly_t res, gr_srcptr nu, const gr_poly_t z, slong len, gr_ctx_t ctx)

Sets `res` to the Bessel function $J_\nu(z) \bmod x^{\mathit{len}}$ of
the first kind for a power series $z$, by computing the Taylor
coefficients of $J_\nu$ at $z(0)$ with `gr_bessel_j_jet` and composing
them with $z - z(0)$. The underscore version allows aliasing.
:::

## Precomputed moduli {#preinv}

Arithmetic modulo a fixed polynomial $f$, as in finite field extensions
and factorisation algorithms, benefits from precomputation. A
`gr_poly_preinv_t` stores precomputed data for $f$ in one of several
representations, and the functions taking it hide the choice from the
algorithms built on top:

- `GR_POLY_PREINV_NEWTON`: the inverse of the reversal of $f$ as a power
  series to length $\operatorname{len}(f)$, for Newton division with two
  multiplications per reduction ([](algorithms/polynomials/newton#divrem)).
  This is the default. The same kind without a stored inverse
  (`_gr_poly_preinv_set_plain`) means that the ring's ordinary division
  is used, which is preferable for short moduli.
- `GR_POLY_PREINV_SPARSE`: the exponents and coefficients (negated and
  divided by the leading coefficient) of the nonzero terms of $f$ below
  the leading term. A reduction costs one multiplication per nonzero
  term for each eliminated coefficient, which beats Newton division when
  $f$ has few terms, such as the trinomials and pentanomials that
  define extension fields, and does not depend on fast multiplication.
- `GR_POLY_PREINV_TRANSFORMED`: for long dense moduli, precomputed
  transforms of $f$ and of its inverse, as in Shoup's NTL. The quotient
  is the low part of $\operatorname{rev}(A_{\mathrm{hi}}) \cdot
  f^{-1}$, computed in a ring of transformed polynomials
  (`gr_ctx_init_gr_poly_transformed_repr`) with the transform of the
  inverse precomputed, and the remainder $A - Qf$ is computed modulo
  $x^L - 1$ with $L \ge \deg f$, which loses nothing since the remainder
  has degree less than $\deg f$, in a cyclic ring of transformed
  polynomials (`gr_ctx_init_gr_poly_transformed_cyclic_repr`) with the
  transform of $f \bmod (x^L - 1)$ precomputed. A reduction then costs
  about one multiplication instead of two. This needs both transformed
  representations from the base ring, currently `nmod` and `mpn_mod`
  with `fft_small` (`GR_UNABLE` otherwise). These rings select it for
  dense moduli of length at least 160, except for tiny `nmod` moduli,
  where the ordinary multiplication packs two coefficients per
  transform slot: there it is used from length 262144 on for $p \le 3$
  and from 16384 on for $4 \le p \le 15$.

Division by a precomputed modulus is the ring method
`GR_METHOD_POLY_DIVREM_PREINV`, which `nmod` overrides for the sparse
representation; the default `gr_generic_poly_divrem_preinv` handles all
representations. The functions that take $f$ and its Newton inverse as
separate polynomials (`gr_poly_mulmod_preinv`,
`gr_poly_powmod_fmpz_sliding_preinv` and so on) are wrappers that build
a shallow Newton representation.

::: type
gr_poly_preinv_struct
gr_poly_preinv_t

Precomputed data for reduction modulo a polynomial $f$. The fields
`kind` (one of the representations below), `lenf` (the length of $f$),
`f` (its coefficients) and `monic` are public; the other fields depend
on the representation. An initialised object can be used from several
threads at once: the transformed representation keeps its large scratch
elements in a pool inside the object, handed out under a lock.
:::

::: macro
GR_POLY_PREINV_NEWTON
GR_POLY_PREINV_SPARSE
GR_POLY_PREINV_TRANSFORMED

The values 0, 1 and 2 of the `kind` field of a `gr_poly_preinv_t`.
:::

::: function
void gr_poly_preinv_init(gr_poly_preinv_t P, gr_ctx_t ctx)
void gr_poly_preinv_clear(gr_poly_preinv_t P, gr_ctx_t ctx)

Initialise `P` (to an empty modulus) and clear it.
:::

::: function
int gr_poly_preinv_set(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_preinv_set(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_poly_preinv_set_newton(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_preinv_set_newton(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int _gr_poly_preinv_set_plain(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_poly_preinv_set_sparse(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_preinv_set_sparse(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_poly_preinv_set_transformed(gr_poly_preinv_t P, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_preinv_set_transformed(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_ctx_t ctx)

Set `P` to the nonzero modulus $f$, storing a copy of it, with a
representation chosen automatically, the Newton representation, the
Newton representation without an inverse (ordinary division), the
sparse representation or the transformed representation (`GR_UNABLE`
if unavailable). Return `GR_DOMAIN` if $f = 0$ and `GR_UNABLE` if its
leading coefficient cannot be inverted.

::: details
The automatic choice is the ring method `GR_METHOD_POLY_PREINV_SET`,
which can take the degree, the sparsity and the coefficient size into
account. The default `gr_generic_poly_preinv_set` chooses ordinary
division for moduli of length at most 8, the sparse representation for
at most 8 nonzero terms below the leading term, and the Newton
representation otherwise. Rings with cheap elements use larger cutoffs
for ordinary division (32 for `nmod`, with a sparse cutoff of 5
terms). `nmod`, `fmpz_mod`, `mpn_mod` and the
contexts of `gr_ctx_init_gr_poly_quotient` override the method.
:::
:::

::: function
void _gr_poly_preinv_init_newton_shallow(gr_poly_preinv_t P, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)

Initialises `P` as a Newton representation referencing, without
copying, the modulus $f$ and its precomputed inverse (which may be
empty, $\mathit{lenfinv} = 0$, to request ordinary division). The object
must not be cleared, and the referenced arrays must outlive it.
:::

::: function
int gr_generic_poly_preinv_set(gr_poly_preinv_struct * P, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_generic_poly_divrem_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, const gr_poly_preinv_struct * P, gr_ctx_t ctx)

The default implementations of the methods `GR_METHOD_POLY_PREINV_SET`
and `GR_METHOD_POLY_DIVREM_PREINV`, declared in `gr_generic.h`. A ring
that overrides one of these methods can call them for the cases it
does not handle itself.
:::

::: function
int _gr_poly_preinv_divrem(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_rem(gr_ptr R, gr_srcptr A, slong lenA, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_rem(gr_poly_t R, const gr_poly_t A, const gr_poly_preinv_t P, gr_ctx_t ctx)

Division with remainder by the precomputed modulus. Dividends of any
length are supported; Newton division is applied blockwise to long
dividends. The underscore *divrem* requires $\operatorname{len}(A) \ge
\operatorname{len}(f)$ and writes $\operatorname{len}(A) -
\operatorname{len}(f) + 1$ quotient and $\operatorname{len}(f) - 1$
remainder coefficients (not normalised); `Q` may be `NULL` when only
the remainder is wanted.
:::

::: function
int _gr_poly_preinv_mulmod(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_mulmod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_preinv_t P, gr_ctx_t ctx)

Set `res` to the product modulo $f$. The underscore version requires
reduced inputs (lengths less than $\operatorname{len}(f)$), writes
$\operatorname{len}(f) - 1$ coefficients, allows aliasing, and skips the
reduction when the product has degree less than $\deg f$.
:::

::: function
int _gr_poly_preinv_powmod_fmpz_binexp(gr_ptr res, gr_srcptr poly, slong len, const fmpz_t e, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_powmod_fmpz_binexp(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_powmod_fmpz_sliding(gr_ptr res, gr_srcptr poly, slong len, const fmpz_t e, ulong k, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_powmod_fmpz_sliding(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, ulong k, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_powmod_x_fmpz(gr_ptr res, const fmpz_t e, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_powmod_x_fmpz(gr_poly_t res, const fmpz_t e, const gr_poly_preinv_t P, gr_ctx_t ctx)

Modular exponentiation with the precomputed modulus: binary
exponentiation, sliding window exponentiation with window size `k` (0
for an automatic choice), and powers of $x$. The algorithms and the
conventions of the underscore versions are those of the corresponding
`gr_poly_powmod` functions below.
:::

::: function
int _gr_poly_preinv_compose_mod(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_compose_mod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_compose_mod_horner(gr_ptr res, gr_srcptr f, slong lenf, gr_srcptr g, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_compose_mod_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_precompute_matrix(gr_mat_t A, gr_srcptr poly1, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_precompute_matrix(gr_mat_t A, const gr_poly_t poly1, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_compose_mod_brent_kung_precomp(gr_ptr res, gr_srcptr poly1, slong len1, const gr_mat_t A, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_compose_mod_brent_kung_precomp(gr_poly_t res, const gr_poly_t poly1, const gr_mat_t A, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_compose_mod_brent_kung_vec(gr_poly_struct * res, const gr_poly_struct * polys, slong lenpolys, slong l, gr_srcptr g, slong glen, const gr_poly_preinv_t P, gr_ctx_t ctx)
int gr_poly_preinv_compose_mod_brent_kung_vec(gr_poly_struct * res, const gr_poly_struct * polys, slong len1, slong n, const gr_poly_t g, const gr_poly_preinv_t P, gr_ctx_t ctx)
int _gr_poly_preinv_mod_matrix_rows_evaluate(gr_ptr res, const gr_mat_t A, gr_srcptr h, slong n, const gr_poly_preinv_t P, gr_ctx_t ctx)

Modular composition $f(g) \bmod h$ with the precomputed modulus $h$.
The algorithms and conventions are those of the corresponding
`gr_poly_compose_mod` functions and of `gr_poly_precompute_matrix`
below, which are wrappers of these.
:::

## Modular arithmetic and composition {#modular}

The functions in this section take the modulus $f$ and, for the
*preinv* versions, the inverse `finv` of the reversal of $f$ modulo
$x^{\operatorname{len}(f)}$ as separate polynomials. They are wrappers
of the functions with a `gr_poly_preinv_t`, which should be preferred
for repeated use of one modulus. Modular composition is explained in
[](algorithms/polynomials/composition#modular).

::: function
int _gr_poly_mulmod(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_poly_mulmod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_mulmod_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
int gr_poly_mulmod_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)

Set `res` to the remainder of the product of `poly1` and `poly2` by
$f$.

::: details
`_gr_poly_mulmod` requires $\mathit{len1} + \mathit{len2} - \mathit{lenf}
> 0$, that is a product that actually needs reduction, writes
$\mathit{lenf} - 1$ coefficients, and does not allow `res` to be
aliased with `f`. `_gr_poly_mulmod_preinv` does not require reduced
inputs, but the quotient length $\mathit{len1} + \mathit{len2} -
\mathit{lenf}$ must not exceed $\mathit{lenf}$ (`GR_UNABLE` otherwise),
and `res` must not be aliased with any input.
:::
:::

::: function
int _gr_poly_powmod_ui_binexp(gr_ptr res, gr_srcptr poly, slong len, ulong e, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_poly_powmod_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong e, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_powmod_ui_binexp_preinv(gr_ptr res, gr_srcptr poly, slong len, ulong e, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
int gr_poly_powmod_ui_binexp_preinv(gr_poly_t res, const gr_poly_t poly, ulong e, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
int _gr_poly_powmod_fmpz_binexp(gr_ptr res, gr_srcptr poly, slong len, const fmpz_t e, gr_srcptr f, slong lenf, gr_ctx_t ctx)
int gr_poly_powmod_fmpz_binexp(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, const gr_poly_t f, gr_ctx_t ctx)
int _gr_poly_powmod_fmpz_binexp_preinv(gr_ptr res, gr_srcptr poly, slong len, const fmpz_t e, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
int gr_poly_powmod_fmpz_binexp_preinv(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
int _gr_poly_powmod_fmpz_sliding_preinv(gr_ptr res, gr_srcptr poly, slong len, const fmpz_t e, ulong k, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
int gr_poly_powmod_fmpz_sliding_preinv(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, ulong k, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
int _gr_poly_powmod_x_fmpz_preinv(gr_ptr res, const fmpz_t e, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
int gr_poly_powmod_x_fmpz_preinv(gr_poly_t res, const fmpz_t e, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)

Set `res` to `poly` raised to the power $e$ modulo $f$, by binary
exponentiation, by sliding window exponentiation with window size `k`
($k = 0$ chooses it from the size of $e$), or, for the *x* versions,
$x^e \bmod f$ with shifts instead of multiplications by the base.

::: details
The underscore versions with a base require a reduced base
($\mathit{len} < \mathit{lenf}$), $\mathit{lenf} \ge 2$ and $e \ge 0$,
write $\mathit{lenf} - 1$ coefficients (zero-padded, not normalised),
and allow `res` to be `poly`; the *x* version requires $\mathit{lenf}
\ge 3$. The *preinv* versions accept $\mathit{lenfinv} = 0$, in which
case ordinary division is used. The intermediate powers are kept at
their actual length, so that powers of a polynomial of small degree,
such as $x$, cost nothing until the degree reaches $\deg f$; the
iterated Frobenius $x^{p^i}$ is thus free while $p^i < \deg f$.
:::
:::

::: function
int _gr_poly_compose_mod(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_ctx_t ctx)
int gr_poly_compose_mod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)
int _gr_poly_compose_mod_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_ctx_t ctx)
int gr_poly_compose_mod_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)
int _gr_poly_compose_mod_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_ctx_t ctx)
int gr_poly_compose_mod_brent_kung(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)
int _gr_poly_compose_mod_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong inv3len, gr_ctx_t ctx)
int gr_poly_compose_mod_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)
int _gr_poly_compose_mod_horner_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong inv3len, gr_ctx_t ctx)
int gr_poly_compose_mod_horner_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)
int _gr_poly_compose_mod_brent_kung_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong inv3len, gr_ctx_t ctx)
int gr_poly_compose_mod_brent_kung_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)

Set `res` to $f(g) \bmod h$, where $f$, $g$, $h$ are `poly1`, `poly2`,
`poly3`, by Horner's rule ($\deg h$ modular multiplications) or by the
Brent-Kung baby-step giant-step algorithm [@BrentKung1978], which
computes $m \approx \sqrt{\deg h}$ powers of $g$, multiplies a matrix of
the coefficients of $f$ by the matrix of these powers, and combines the
results by Horner's rule in $g^m$, so that the cost is dominated by one
$\sqrt{n} \times n$ by $n \times \sqrt{n}$ matrix product. The *preinv*
versions take the inverse `poly3inv` of the reversal of $h$.

::: details
The default versions use Horner's rule when $\mathit{len3} < 6$ or
$\mathit{len1} \ge \mathit{len3}$ and Brent-Kung otherwise; the
Brent-Kung versions themselves return `GR_UNABLE` when
$\mathit{len1} \ge \mathit{len3}$, that is when $f$ is not reduced
modulo $h$. The non-underscore versions reduce $g$ modulo $h$ but not
$f$. The underscore functions require nonzero lengths and $g$ given with
exactly $\mathit{len3} - 1$ coefficients (zero-padded), write
$\mathit{len3} - 1$ coefficients, and do not allow aliasing.
:::
:::

::: function
int gr_poly_compose_mod_wrapper(_gr_method_compose_mod_op _compose_mod, gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)
int gr_poly_compose_mod_preinv_wrapper(_gr_method_compose_mod_preinv_op _compose_mod, gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)

Implement the non-underscore modular composition functions around an
underscore implementation: return `GR_DOMAIN` if $h = 0$, handle
constant inputs, reduce $g$ modulo $h$ (or zero-pad it to length
$\deg h$), manage memory and aliasing, and normalise the result.
:::

::: function
int _gr_poly_reduce_matrix_mod_poly(gr_mat_t A, const gr_mat_t B, const gr_poly_t f, gr_ctx_t ctx)

Sets row $i$ of $A$ to the reduction of row $i$ of $B$ (read as the
coefficients of a polynomial) modulo $f$, for $1 \le i \le
\sqrt{\deg f}$. $B$ must have at least $\sqrt{\deg f}$ rows and $\deg f$
columns, and $f$ must be nonzero.
:::

::: function
int _gr_poly_precompute_matrix(gr_mat_t A, gr_srcptr f, gr_srcptr g, slong leng, gr_srcptr ginv, slong lenginv, gr_ctx_t ctx)
int gr_poly_precompute_matrix(gr_mat_t A, const gr_poly_t f, const gr_poly_t g, const gr_poly_t ginv, gr_ctx_t ctx)

Set row $i$ of $A$ to $f^i \bmod g$ for $0 \le i < m$, $m = \lfloor
\sqrt{\deg g} \rfloor + 1$, the baby steps of Brent-Kung modular
composition with argument $f$ and modulus $g$. $A$ must have exactly
$m$ rows and $\deg g$ columns (`GR_DOMAIN` otherwise, for the
non-underscore version), $g$ must be nonzero and `ginv` must be the
inverse of the reversal of $g$.
:::

::: function
int _gr_poly_compose_mod_brent_kung_precomp_preinv(gr_ptr res, gr_srcptr f, slong lenf, const gr_mat_t A, gr_srcptr h, slong lenh, gr_srcptr hinv, slong lenhinv, gr_ctx_t ctx)
int gr_poly_compose_mod_brent_kung_precomp_preinv(gr_poly_t res, const gr_poly_t f, const gr_mat_t A, const gr_poly_t h, const gr_poly_t hinv, gr_ctx_t ctx)

Set `res` to $f(g) \bmod h$ given the matrix $A$ of powers of $g$ from
`gr_poly_precompute_matrix`, so that several compositions with the same
$g$ and $h$ share the baby steps. $h$ must be nonzero, $f$ must have
smaller degree than $h$, and `hinv` must be the inverse of the reversal
of $h$. The underscore version does not allow the output to be aliased
with an input.
:::

::: function
int _gr_poly_compose_mod_brent_kung_vec_preinv(gr_poly_struct * res, const gr_poly_struct * polys, slong lenpolys, slong l, gr_srcptr g, slong glen, gr_srcptr poly, slong len, gr_srcptr polyinv, slong leninv, gr_ctx_t ctx)
int gr_poly_compose_mod_brent_kung_vec_preinv(gr_poly_struct * res, const gr_poly_struct * polys, slong len1, slong n, const gr_poly_t g, const gr_poly_t poly, const gr_poly_t polyinv, gr_ctx_t ctx)

Set `res[i]` to `polys[i]` composed with $g$ modulo `poly` for $0 \le i
< n$ (respectively $l$), where `polys` holds `len1` (respectively
`lenpolys`) polynomials of degree less than that of `poly` and
`polyinv` is the inverse of the reversal of `poly`. Compared with $n$
separate compositions, a single matrix product with splitting parameter
$m \approx \sqrt{n \deg(\mathit{poly})}$ reduces the number of modular
multiplications by a factor $\sqrt{n}$.

::: details
The underscore version requires $g$ of length at most
$\operatorname{len}(\mathit{poly}) - 1$ and room for
$\operatorname{len}(\mathit{poly}) - 1$ coefficients in each `res[i]`.
The non-underscore version does not allow the output to be aliased
with $g$, `poly` or `polyinv`.
:::
:::

::: function
int _gr_poly_mod_matrix_rows_evaluate(gr_ptr res, const gr_mat_t A, gr_srcptr h, slong n, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong len3inv, gr_ctx_t ctx)

The last step of Brent-Kung composition: sets `res` to $\sum_i c_i h^i
\bmod \mathit{poly3}$, where the $c_i$ are the rows of $A$ and $h$ (in
Brent-Kung composition, $g^m$) and the rows have $n = \mathit{len3} - 1$
coefficients. For few rows this is Horner's rule; for many rows it uses
rectangular splitting with powers $h^0, \ldots, h^k$, $k \approx
\sqrt{r}$, and unreduced products, so that only about $r/k + k$
modular reductions are needed for $r$ rows.
:::

## Transformed representations {#transformed}

Some rings provide rings of *transformed polynomials*, in which
polynomials are stored as transforms (for example, values under a
number-theoretic transform) so that products are pointwise operations;
the contexts are constructed with
`gr_ctx_init_gr_poly_transformed_repr` and
`gr_ctx_init_gr_poly_transformed_cyclic_repr` (see [](gr)). The
half-gcd matrix products and division by a transformed precomputed
modulus use them. The functions below convert between ordinary
coefficient arrays over `base_ctx` and elements of such a ring `tctx`;
each calls a method of `tctx`.

::: function
int _gr_set_gr_poly(gr_ptr res, gr_srcptr a, slong len, gr_ctx_t base_ctx, gr_ctx_t tctx)

Sets the transformed element `res` to the polynomial with the `len`
coefficients `a`, which must not exceed the length bound of `tctx`.
:::

::: function
int _gr_get_gr_poly(gr_ptr c, slong * len, gr_srcptr x, gr_ctx_t base_ctx, gr_ctx_t tctx)
int _gr_get_gr_poly_destructive(gr_ptr c, slong * len, gr_ptr x, gr_ctx_t base_ctx, gr_ctx_t tctx)

Write the coefficients of the polynomial represented by `x` to `c` and
its length to `*len`; `c` must have room for the length bound of
`tctx`. The *destructive* version may run the inverse transform in the
storage of `x`, saving a copy, after which `x` may only be cleared or
overwritten as the destination of an assignment or a product; it falls
back to the copying conversion when the ring does not provide it.
:::

::: function
int _gr_get_gr_poly_window(gr_ptr c, gr_srcptr x, slong zl, slong zh, gr_ctx_t base_ctx, gr_ctx_t tctx)
int _gr_get_gr_poly_window_destructive(gr_ptr c, gr_ptr x, slong zl, slong zh, gr_ctx_t base_ctx, gr_ctx_t tctx)

Write the coefficients of index $\mathit{zl} \le i < \mathit{zh}$ of the
polynomial represented by `x` to `c` (zeros beyond its length), as
needed for middle products. The *destructive* version is as above.
:::

## Default method implementations {#generic}

The functions in this section are the default implementations of the
ring methods for polynomials, installed in the method table of every
context that does not provide its own (see [](gr_generic)). A ring that
overrides a method can call the default for the cases it does not
handle. The choices are described with the corresponding functions
above.

| Method | Called by | Default |
|---|---|---|
| `POLY_MULLOW` | `_gr_poly_mullow` | `_gr_poly_mullow_generic` |
| `POLY_MULMID` | `_gr_poly_mulmid` | `_gr_poly_mulmid_generic` |
| `POLY_DIVREM` | `_gr_poly_divrem` | `_gr_poly_divrem_generic` |
| `POLY_DIV` | `_gr_poly_div` | `_gr_poly_div_generic` |
| `POLY_DIVEXACT` | `_gr_poly_divexact` | `_gr_poly_divexact_generic` |
| `POLY_INV_SERIES` | `_gr_poly_inv_series` | `_gr_poly_inv_series_generic` |
| `POLY_INV_SERIES_BASECASE` | `_gr_poly_inv_series_basecase` | `_gr_poly_inv_series_basecase_generic` |
| `POLY_DIV_SERIES` | `_gr_poly_div_series` | `_gr_poly_div_series_generic` |
| `POLY_DIV_SERIES_BASECASE` | `_gr_poly_div_series_basecase` | `_gr_poly_div_series_basecase_generic` |
| `POLY_SQRT_SERIES` | `_gr_poly_sqrt_series` | `_gr_poly_sqrt_series_generic` |
| `POLY_RSQRT_SERIES` | `_gr_poly_rsqrt_series` | `_gr_poly_rsqrt_series_generic` |
| `POLY_EXP_SERIES` | `_gr_poly_exp_series` | `_gr_poly_exp_series_generic` |
| `POLY_TAYLOR_SHIFT` | `_gr_poly_taylor_shift` | `_gr_poly_taylor_shift_generic` |
| `POLY_EVALUATE` | `_gr_poly_evaluate` | `gr_generic_poly_evaluate` |
| `POLY_GCD` | `_gr_poly_gcd` | `_gr_poly_gcd_generic` |
| `POLY_XGCD` | `_gr_poly_xgcd` | `_gr_poly_xgcd_generic` |
| `POLY_HGCD_MAT_MUL` | `_gr_poly_hgcd_mat_mul` | `_gr_poly_hgcd_mat_mul_generic` |
| `POLY_RESULTANT` | `_gr_poly_resultant` | `_gr_poly_resultant_generic` |
| `POLY_SET_LENGTH_NORMALISE` | `_gr_poly_set_length_normalise` | `gr_generic_poly_set_length_normalise` |
| `POLY_PREINV_SET` | `_gr_poly_preinv_set` | `gr_generic_poly_preinv_set` |
| `POLY_DIVREM_PREINV` | `_gr_poly_preinv_divrem` | `gr_generic_poly_divrem_preinv` |
| `POLY_FACTOR` | `gr_factor` on polynomials | finite fields only |
| `POLY_ROOTS`, `POLY_ROOTS_OTHER` | `gr_poly_roots`, `gr_poly_roots_other` | finite fields only |

::: function
int _gr_poly_mullow_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
int _gr_poly_mulmid_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

The default multiplication methods. `_gr_poly_mullow_generic` is
classical multiplication. `_gr_poly_mulmid_generic` calls the ring's
`_gr_poly_mullow` for the product up to $\mathit{nhi}$ and discards the
low part when that method is not the default, and uses classical
multiplication otherwise. Generic algorithms compare the method table
entries with these functions to decide whether fast multiplication or
a native middle product is available.
:::

::: function
int _gr_poly_divrem_generic(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int _gr_poly_div_generic(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
int _gr_poly_divexact_generic(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)

The default division methods; see `gr_poly_divrem`, `gr_poly_div` and
`gr_poly_divexact`.
:::

::: function
int _gr_poly_inv_series_generic(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx)
int _gr_poly_inv_series_basecase_generic(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_generic(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int _gr_poly_div_series_basecase_generic(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
int _gr_poly_sqrt_series_generic(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int _gr_poly_rsqrt_series_generic(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
int _gr_poly_taylor_shift_generic(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)

The default methods for power series inversion and division (with
their basecases, which invert the constant term once if possible and
otherwise divide by it at each step), square roots and Taylor shifts;
see `gr_poly_inv_series`, `gr_poly_div_series`,
`gr_poly_sqrt_series` and `gr_poly_taylor_shift`.
:::

::: function
int _gr_poly_hgcd_mat_mul(gr_ptr * C, slong * lenC, gr_ptr * A, slong * lenA, gr_ptr * B, slong * lenB, gr_ptr T0, gr_ptr T1, gr_ctx_t ctx)
int _gr_poly_hgcd_mat_mul_generic(gr_ptr * C, slong * lenC, gr_ptr * A, slong * lenA, gr_ptr * B, slong * lenB, gr_ptr T0, gr_ptr T1, gr_ctx_t ctx)

Set the $2 \times 2$ matrix of polynomials $C$ to $AB$, as needed by the
half-gcd algorithm. `A`, `B` and `C` point to 4 polynomials each, in
row-major order, with their lengths in `lenA`, `lenB` and `lenC`; the
output lengths are normalised. Aliasing is not allowed. `T0` and `T1`
must have room for any product of an entry of $A$ by an entry of $B$.
`_gr_poly_hgcd_mat_mul` calls the ring method; the default is described
with `_gr_poly_hgcd`.
:::

## Test functions {#testing}

The following functions run `iters` random tests of an implementation
of a method, with polynomials of length up to `maxn`. If `ctx` is
`NULL`, a random ring is chosen for each iteration; otherwise the given
ring is used.

::: function
void _gr_poly_test_mullow(gr_method_poly_binary_trunc_op mullow_impl, gr_method_poly_binary_trunc_op mullow_ref, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void _gr_poly_test_mulmid(gr_method_poly_binary_trunc2_op mulmid_impl, gr_method_poly_binary_trunc2_op mulmid_ref, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

Test `mullow_impl` (`mulmid_impl`) as an implementation of
`_gr_poly_mullow` (`_gr_poly_mulmid`), comparing with `mullow_ref`
(`mulmid_ref`), or with classical multiplication if the reference is
`NULL`.
:::

::: function
void _gr_poly_test_approx_mulmid_pos_entrywise_accurate(gr_method_poly_binary_trunc2_op mulmid_impl, gr_method_poly_binary_trunc2_op mulmid_ref, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

Tests the accuracy of a middle product implementation over an
approximate (floating-point) ring: for polynomials with nonnegative
coefficients, the relative error of every output coefficient compared
with the reference implementation must not exceed `rel_tol`.
:::

::: function
void _gr_poly_test_divrem(gr_method_poly_binary_binary_op divrem_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void _gr_poly_test_div(gr_method_poly_binary_op div_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void _gr_poly_test_inv_series(gr_method_poly_unary_trunc_op inv_series_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void _gr_poly_test_div_series(gr_method_poly_binary_trunc_op div_series_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

Test an implementation of `_gr_poly_divrem`, `_gr_poly_div`,
`_gr_poly_inv_series` or `_gr_poly_div_series`.
:::

::: function
void _gr_poly_test_gcd_field(gr_method_poly_gcd_op gcd_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void _gr_poly_test_gcd_ufd(gr_method_poly_gcd_op gcd_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

Test an implementation of `_gr_poly_gcd` that is expected to succeed
over fields, respectively over unique factorisation domains.
:::

::: function
void _gr_poly_test_xgcd(gr_method_poly_xgcd_op xgcd_impl, flint_rand_t state, slong iters, slong maxn, int flags, gr_ctx_t ctx)

Tests an implementation of `_gr_poly_xgcd`. If `flags` is 1, it is
expected to succeed over GCD domains, otherwise only over fields.
:::
