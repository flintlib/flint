---
title: nmod_poly.h
header: nmod_poly.h
description: Dense univariate polynomials over $\Z/n\Z$ for a word-size modulus $n$, from basic arithmetic to fast multiplication, division, GCDs, evaluation, interpolation, composition and power series.
category: Univariate polynomials
order: 5
---

This module implements dense univariate polynomials over $\Z/n\Z$ for a
word-size modulus $n \ge 1$. A polynomial is an `nmod_poly_t`: an array of
`ulong` coefficients reduced modulo $n$, together with its length and an
`nmod_t` for the modulus. Multiplication switches between the classical
algorithm, several variants of
[Kronecker substitution](algorithms/polynomials/kronecker) and FFTs modulo
word-size primes ([](algorithms/fft/small-prime-fft)), so that a product of
length $m$ costs $O(m \log m)$ word operations for large $m$. Division,
power series inversion, GCDs, multipoint evaluation, interpolation and
composition are reduced to multiplication by Newton iteration, the
half-GCD algorithm and subproduct trees, and inherit its quasi-linear
complexity. Factorisation and irreducibility testing are in
`nmod_poly_factor.h`, and polynomial matrices in `nmod_poly_mat.h`.

The guides below explain how to use the module; start with
[](nmod_poly/getting-started).

::: summary algorithms/polynomials/multiplication :::

## Types and conventions {#conventions}

::: type
nmod_poly_struct
nmod_poly_t

A polynomial over $\Z/n\Z$. The structure has the fields `coeffs` (an
array of `alloc` coefficients), `alloc`, `length` and `mod` (an `nmod_t`
for $n$, see [](nmod)). An `nmod_poly_t` is an array of one
`nmod_poly_struct`, so it is passed by reference. The coefficient of $x^i$
is `coeffs[i]`.
:::

- The modulus is stored in each polynomial. All polynomials passed to one
  function must have the same modulus, which is not checked.
- A polynomial is *normalised* if its length is zero or its leading
  coefficient `coeffs[length - 1]` is nonzero. The zero polynomial has
  length 0 and degree $-1$. Coefficients are *reduced*, that is, lie in
  $[0, n)$. All functions require normalised inputs with reduced
  coefficients and produce such outputs. This is not checked.
- Functions without a leading underscore manage memory and allow every
  output to be the same object as any input, unless stated otherwise.
  Several outputs of one function (such as quotient and remainder) must be
  distinct objects.
- Functions with a leading underscore work on plain coefficient arrays
  `(poly, len)`, take the modulus as an `nmod_t`, and leave memory
  management to the caller. They usually require nonzero lengths, place
  ordering conditions on them, write a fixed number of coefficients
  (including leading zeros), and do not allow the output to overlap the
  inputs. Each function states its conditions.
- Many functions need to invert leading coefficients, or integers up to
  the length, and therefore assume that $n$ is prime: division by a
  polynomial whose leading coefficient is not a unit, GCDs, resultants,
  modular inverses, square roots, interpolation and the transcendental
  power series. This is not checked. A non-invertible element usually
  causes an exception (`FLINT_IMPINV` or `FLINT_ERROR`), but a composite
  modulus can also give wrong results without warning. The GCD and
  resultant functions always throw `FLINT_IMPINV` when they meet a
  leading coefficient that is not invertible.

The conventions for arithmetic modulo a fixed polynomial are explained in
[](nmod_poly/modular-arithmetic).

## Memory management {#memory}

::: function
void nmod_poly_init(nmod_poly_t poly, ulong n)
void nmod_poly_init_preinv(nmod_poly_t poly, ulong n, ulong ninv)

Initialise `poly` as the zero polynomial modulo $n \ge 1$, without
allocating coefficients. `nmod_poly_init_preinv` takes the precomputed
inverse `ninv = n_preinvert_limb(n)` instead of computing it.
:::

::: function
void nmod_poly_init_mod(nmod_poly_t poly, const nmod_t mod)

Initialises `poly` as the zero polynomial with the modulus `mod`, which
must have been set up with `nmod_init`. This is the cheapest way to create
polynomials with the modulus of an existing one, `nmod_poly_init_mod(a,
b->mod)`.
:::

::: function
void nmod_poly_init2(nmod_poly_t poly, ulong n, slong alloc)
void nmod_poly_init2_preinv(nmod_poly_t poly, ulong n, ulong ninv, slong alloc)

Like `nmod_poly_init` and `nmod_poly_init_preinv`, but allocate space for
`alloc` coefficients.
:::

::: function
void nmod_poly_realloc(nmod_poly_t poly, slong alloc)

Changes the allocated space of `poly` to `alloc` coefficients. If the
length exceeds `alloc`, the polynomial is truncated to length `alloc` and
normalised. For `alloc` = 0 the coefficients are freed and `poly` becomes
zero.
:::

::: function
void nmod_poly_fit_length(nmod_poly_t poly, slong alloc)

Ensures that `poly` has space for at least `alloc` coefficients. If it has
to grow, it grows to at least twice the old allocation, so that repeated
calls with increasing lengths cost amortised linear time. The length and
the coefficients are not changed.
:::

::: function
void nmod_poly_clear(nmod_poly_t poly)

Frees the memory used by `poly`, which must be initialised again before it
is used.
:::

::: function
void nmod_poly_set_mod(nmod_poly_t poly, const nmod_t mod)

Replaces the modulus of `poly` by `mod` without changing the
coefficients, which must be reduced modulo the new modulus.
:::

::: function
void _nmod_poly_set_length(nmod_poly_t poly, slong len)

Sets the length field of `poly` to `len` without zeroing or normalising
anything. `poly` must have space for `len` coefficients.
:::

::: function
void _nmod_poly_normalise(nmod_poly_t poly)

Decreases the length of `poly` until the leading coefficient is nonzero or
the length is zero. This is the usual last step after writing coefficients
with an underscore function.
:::

## Length, degree and basic properties {#properties}

::: function
slong nmod_poly_length(const nmod_poly_t poly)
slong nmod_poly_degree(const nmod_poly_t poly)

Return the length and the degree of `poly`. The zero polynomial has length
0 and degree $-1$.
:::

::: function
ulong nmod_poly_modulus(const nmod_poly_t poly)

Returns the modulus $n$ of `poly`.
:::

::: function
flint_bitcnt_t nmod_poly_max_bits(const nmod_poly_t poly)

Returns the number of bits of the largest coefficient of `poly`, or 0 for
the zero polynomial.
:::

::: function
nn_ptr nmod_poly_lead(const nmod_poly_t poly)

Returns a pointer to the leading coefficient of `poly`, or `NULL` if
`poly` is zero.
:::

::: function
slong _nmod_poly_hamming_weight(nn_srcptr a, slong len)
slong nmod_poly_hamming_weight(const nmod_poly_t A)

Return the number of nonzero coefficients.
:::

::: function
int nmod_poly_is_zero(const nmod_poly_t poly)

Returns whether `poly` is the zero polynomial.
:::

::: function
int nmod_poly_is_one(const nmod_poly_t poly)
int nmod_poly_is_gen(const nmod_poly_t poly)

Return whether `poly` is the constant 1, respectively the generator $x$.
For $n = 1$ both return 1, since every polynomial is zero modulo 1.
:::

::: function
int nmod_poly_is_unit(const nmod_poly_t poly)

Returns whether `poly` is a nonzero constant, which for prime $n$ is the
same as being a unit. For composite $n$ the result does not tell whether
the constant is invertible. For $n = 1$ it returns 1.
:::

::: function
int nmod_poly_is_monic(const nmod_poly_t poly)

Returns whether `poly` is nonzero with leading coefficient 1.
:::

## Assignment and coefficients {#assignment}

::: function
void nmod_poly_set(nmod_poly_t a, const nmod_poly_t b)

Sets `a` to a copy of `b`.
:::

::: function
void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)

Swaps `poly1` and `poly2` in constant time by swapping the structures,
including the moduli.
:::

::: function
void nmod_poly_zero(nmod_poly_t res)
void nmod_poly_one(nmod_poly_t res)

Set `res` to 0 and to 1. Modulo $n = 1$, `nmod_poly_one` sets `res` to the
zero polynomial.
:::

::: function
void nmod_poly_set_trunc(nmod_poly_t res, const nmod_poly_t poly, slong len)
void nmod_poly_truncate(nmod_poly_t poly, slong len)

Set `res` to `poly` truncated to length `len`, that is, `poly` modulo
$x^{\mathit{len}}$, and truncate `poly` in place. The result is normalised.
Nothing happens if `len` is at least the length.
:::

::: function
void _nmod_poly_reverse(nn_ptr output, nn_srcptr input, slong len, slong m)
void nmod_poly_reverse(nmod_poly_t output, const nmod_poly_t input, slong m)

Set `output` to the reversal $x^{m-1} f(1/x)$ of the input $f$ viewed as a
polynomial of length $m \ge 0$: coefficient $i$ of the output is
coefficient $m - 1 - i$ of the input, taken as 0 beyond its length. Input
coefficients from index $m$ on are ignored.

::: details
`_nmod_poly_reverse` writes exactly $m$ coefficients, including leading
zeros; `nmod_poly_reverse` normalises the result. The underscore version
allows `output == input` but no partial overlap. Reversal turns division
with remainder into power series division, which is why the
`_preinv` functions of this module take the inverse of the reversed
divisor; see [](nmod_poly/modular-arithmetic).
:::
:::

::: function
ulong nmod_poly_get_coeff_ui(const nmod_poly_t poly, slong j)

Returns the coefficient of $x^j$, or 0 if $j$ is at least the length.
Requires $j \ge 0$.
:::

::: function
void nmod_poly_set_coeff_ui(nmod_poly_t poly, slong j, ulong c)

Sets the coefficient of $x^j$ to $c \bmod n$, for any $c$. The polynomial
grows, with zero coefficients in between, or shrinks if the leading
coefficient is set to zero.
:::

## Comparison {#comparison}

::: function
int nmod_poly_equal(const nmod_poly_t a, const nmod_poly_t b)

Returns whether `a` and `b` are equal. The moduli are not compared.
:::

::: function
int nmod_poly_equal_nmod(const nmod_poly_t poly, ulong cst)
int nmod_poly_equal_ui(const nmod_poly_t poly, ulong cst)

Return whether `poly` is the constant `cst`. `nmod_poly_equal_nmod`
requires `cst` to be reduced; `nmod_poly_equal_ui` reduces it first.
:::

::: function
int nmod_poly_equal_trunc(const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)

Returns whether `poly1` and `poly2` agree modulo $x^n$. A negative $n$ is
treated as 0.
:::

## Random polynomials {#random}

The `rand` functions draw coefficients uniformly from $[0, n)$. The
`randtest` functions use `n_randtest`, which favours special values and
sparse polynomials; they are meant for test code. Functions that produce
irreducible or primitive polynomials repeat until
`nmod_poly_is_irreducible` (from `nmod_poly_factor.h`) succeeds and
require a prime modulus.

::: function
void nmod_poly_rand(nmod_poly_t poly, flint_rand_t state, slong len)
void nmod_poly_rand_monic(nmod_poly_t poly, flint_rand_t state, slong len)

Set `poly` to a random polynomial of length at most `len`, respectively a
random monic polynomial of length exactly $\mathit{len} \ge 1$.
:::

::: function
void nmod_poly_rand_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)

Sets `poly` to a random irreducible polynomial of length at most `len`.
Since `nmod_poly_is_irreducible` regards nonzero constants as irreducible,
the result can be a constant, and it need not be monic.
:::

::: function
void nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, slong len)
void nmod_poly_randtest_not_zero(nmod_poly_t poly, flint_rand_t state, slong len)
void nmod_poly_randtest_monic(nmod_poly_t poly, flint_rand_t state, slong len)

Set `poly` to a random test polynomial of length at most `len`, a nonzero
one (requires $\mathit{len} \ge 1$ and $n \ge 2$), and a monic one of
length exactly $\mathit{len} \ge 1$.
:::

::: function
void nmod_poly_randtest_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)
void nmod_poly_randtest_monic_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)

Set `poly` to a random irreducible polynomial of length at most `len`
(which may be a constant), and to a random monic irreducible polynomial
of length exactly $\mathit{len} \ge 2$.
:::

::: function
void nmod_poly_randtest_monic_primitive(nmod_poly_t poly, flint_rand_t state, slong len)

Sets `poly` to a random monic irreducible polynomial of length
$\mathit{len} \ge 2$ whose root generates the multiplicative group of the
field $\mathbb{F}_n[x]/(\mathit{poly})$. Primitivity is tested with
`fq_nmod_is_primitive`.
:::

::: function
void nmod_poly_randtest_trinomial(nmod_poly_t poly, flint_rand_t state, slong len)
void nmod_poly_randtest_pentomial(nmod_poly_t poly, flint_rand_t state, slong len)

Set `poly` to a random monic polynomial of length `len` with at most three
nonzero terms $x^{\mathit{len}-1} + a x^k + b$, $0 < k < \mathit{len} - 1$
(requires $\mathit{len} \ge 3$), respectively at most five nonzero terms
$x^{\mathit{len}-1} + a_3 x^3 + a_2 x^2 + a_1 x + a_0$ (requires
$\mathit{len} \ge 5$).
:::

::: function
int nmod_poly_randtest_trinomial_irreducible(nmod_poly_t poly, flint_rand_t state, slong len, slong max_attempts)
int nmod_poly_randtest_pentomial_irreducible(nmod_poly_t poly, flint_rand_t state, slong len, slong max_attempts)

Try up to `max_attempts` random trinomials, respectively pentomials, as
above and return 1 as soon as one is irreducible, or 0 if none is. For
`max_attempts` = 0 they try until they succeed, which does not terminate
if there is no irreducible polynomial of that shape.
:::

::: function
void nmod_poly_randtest_sparse_irreducible(nmod_poly_t poly, flint_rand_t state, slong len)

Sets `poly` to a sparse monic irreducible polynomial of length `len`: it
tries $2 \mathit{len}$ trinomials, then $2 \mathit{len}$ pentomials, then
random polynomials with a growing number of terms.
:::

## Input and output {#io}

The plain format is the length, the modulus and the coefficients from the
constant term upwards, for example `"3 7  6 0 5"` for $5x^2 + 6$ modulo 7
(two spaces before the first coefficient). The pretty format writes the
terms from the highest degree down, for example `"5*x^2+6"`.

::: function
char * nmod_poly_get_str(const nmod_poly_t poly)
char * nmod_poly_get_str_pretty(const nmod_poly_t poly, const char * x)

Return `poly` as a string in the plain format and in the pretty format
with variable name `x`. The zero polynomial is `"0 n"` and `"0"`. Free the
string with `flint_free`.
:::

::: function
int nmod_poly_set_str(nmod_poly_t poly, const char * s)

Reads `poly` from a string in the plain format and returns 1, or returns 0
if the string is not in this format. The modulus in the string is ignored
and the coefficients are not reduced: they must be smaller than the
modulus of `poly`.
:::

::: function
int nmod_poly_print(const nmod_poly_t a)
int nmod_poly_fprint(FILE * f, const nmod_poly_t poly)
int nmod_poly_print_pretty(const nmod_poly_t a, const char * x)
int nmod_poly_fprint_pretty(FILE * f, const nmod_poly_t a, const char * x)

Write the polynomial to `stdout` or to `f`, in the plain or the pretty
format. Return a positive value on success and a non-positive value on
failure.
:::

::: function
int nmod_poly_read(nmod_poly_t poly)
int nmod_poly_fread(FILE * f, nmod_poly_t poly)

Read a polynomial in the plain format from `stdin` or from `f`. Unlike
`nmod_poly_set_str`, they reinitialise `poly` with the modulus that was
read. Return 1 on success and 0 if the input is not in the plain format.
The coefficients are not reduced.
:::

## Shifting {#shifting}

::: function
void _nmod_poly_shift_left(nn_ptr res, nn_srcptr poly, slong len, slong k)
void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, slong k)

Set `res` to `poly` times $x^k$, for $k \ge 0$. The underscore version
writes $\mathit{len} + k$ coefficients and allows `res == poly`.
:::

::: function
void _nmod_poly_shift_right(nn_ptr res, nn_srcptr poly, slong len, slong k)
void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, slong k)

Set `res` to `poly` divided by $x^k$, discarding the remainder, for
$k \ge 0$. In the underscore version `len` is the length of the
*output*: it sets `res[i] = poly[i + k]` for $0 \le i < \mathit{len}$, so
`poly` must have at least $\mathit{len} + k$ coefficients. It allows
`res == poly`.
:::

## Addition and subtraction {#addition}

::: function
void _nmod_poly_add(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
void _nmod_poly_sub(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)

Set `res` to the sum and the difference of `(poly1, len1)` and
`(poly2, len2)`, writing $\max(\mathit{len1}, \mathit{len2})$
coefficients, which may include leading zeros. Any lengths, including 0,
are allowed, and `res` may be `poly1` or `poly2`.
:::

::: function
void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1)

Set `res` to $\mathit{poly1} + \mathit{poly2}$, $\mathit{poly1} -
\mathit{poly2}$ and $-\mathit{poly1}$.
:::

::: function
void nmod_poly_add_series(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)
void nmod_poly_sub_series(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)

Set `res` to the sum and the difference modulo $x^n$. A negative $n$ is
treated as 0.
:::

::: function
void nmod_poly_add_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c)
void nmod_poly_sub_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c)

Set `res` to $\mathit{poly} + c$ and $\mathit{poly} - c$, for any $c$.
:::

## Scalar multiplication {#scalar}

::: function
void nmod_poly_scalar_mul_nmod(nmod_poly_t res, const nmod_poly_t poly, ulong c)

Sets `res` to $c \cdot \mathit{poly}$, for reduced $c$. The result is
normalised, which matters for composite $n$.
:::

::: function
void nmod_poly_scalar_addmul_nmod(nmod_poly_t res, const nmod_poly_t poly, ulong c)

Sets `res` to $\mathit{res} + c \cdot \mathit{poly}$, for reduced $c$.
:::

::: function
void _nmod_poly_make_monic(nn_ptr output, nn_srcptr input, slong len, nmod_t mod)
void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input)

Set `output` to `input` divided by its leading coefficient, which must be
invertible modulo $n$ (otherwise `FLINT_IMPINV` is thrown).
`nmod_poly_make_monic` throws `FLINT_DIVZERO` for the zero polynomial; the
underscore version requires $\mathit{len} \ge 1$ and allows
`output == input`.
:::

## Multiplication {#multiplication}

`nmod_poly_mul` and the other functions without an algorithm in their
name choose the algorithm from the lengths and the number of bits of $n$:

- the classical algorithm, computing each coefficient as a dot product
  with `_nmod_vec_dot_rev` (see [](nmod_vec/dot)), for short operands;
- [Kronecker substitution](algorithms/polynomials/kronecker) (`KS`), which
  packs the coefficients into integers and multiplies those with
  `flint_mpn_mul`, and its multipoint variants `KS2` and `KS4` [@Har2009a],
  for medium lengths;
- the [small-prime FFT](algorithms/fft/small-prime-fft) of `fft_small`
  for long operands, when FLINT is built with it (the usual case on
  x86-64 and ARM64). Moduli $n \le 23$ use a repacked product of half the
  length.

The crossover points, tabulated for each bit size of $n$, are explained in
[](nmod_poly/mul). Truncated products (`mullow`), high products
(`mulhigh`) and middle products (`mulmid`) compute part of the product
for about the cost of the part; middle products are the transposes of
ordinary products and are what Newton iteration needs
([](algorithms/polynomials/middle-product)).

::: summary algorithms/polynomials/multiplication :::

### Choosing the algorithm automatically {#mul-general}

::: function
void _nmod_poly_mul(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
void nmod_poly_mul(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

Set `res` to the product $\mathit{poly1} \cdot \mathit{poly2}$. The
underscore version requires $\mathit{len1} \ge \mathit{len2} \ge 1$,
writes $\mathit{len1} + \mathit{len2} - 1$ coefficients and does not
allow `res` to overlap the inputs; `poly1 == poly2` (squaring) is allowed
and faster.
:::

::: function
void _nmod_poly_mullow(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong trunc, nmod_t mod)
void nmod_poly_mullow(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc)

Set `res` to the product truncated to length `trunc`, that is,
$\mathit{poly1} \cdot \mathit{poly2} \bmod x^{\mathit{trunc}}$. The
underscore version requires $\mathit{len1} \ge \mathit{len2} \ge 1$ and
$1 \le \mathit{trunc} \le \mathit{len1} + \mathit{len2} - 1$, writes
`trunc` coefficients and does not allow aliasing. Inputs longer than
`trunc` are allowed.
:::

::: function
void _nmod_poly_mulhigh(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong n, nmod_t mod)
void nmod_poly_mulhigh(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)

Set the coefficients of `res` from index $n$ on to those of the product.
The coefficients below $n$ are either zero or the coefficients of the
product. The underscore version requires $\mathit{len1} \ge \mathit{len2}
\ge 1$, writes $\mathit{len1} + \mathit{len2} - 1$ coefficients and does
not allow aliasing.

::: details
The classical algorithm, which sets the low coefficients to zero, is used
when $\mathit{len1} + \mathit{len2} \le 6$, or when the unreduced
coefficients fit in one word and $\mathit{len1} + \mathit{len2} < 16$.
Otherwise the coefficients from index $n$ on are computed with
`_nmod_poly_mulmid` (with the small-prime FFT for long operands when
`fft_small` is available), and the low coefficients are set to zero.
:::
:::

::: function
void _nmod_poly_mulmid(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
void nmod_poly_mulmid(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong nlo, slong nhi)

Set `res` to the coefficients $\mathit{nlo}, \ldots, \mathit{nhi} - 1$ of
the product, that is, $\lfloor \mathit{poly1} \cdot \mathit{poly2} /
x^{\mathit{nlo}} \rfloor \bmod x^{\mathit{nhi} - \mathit{nlo}}$. The
underscore version requires $\mathit{len1}, \mathit{len2} \ge 1$ in either
order and $0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1} +
\mathit{len2} - 1$, writes $\mathit{nhi} - \mathit{nlo}$ coefficients and
does not allow aliasing. The non-underscore version accepts any $0 \le
\mathit{nlo} \le \mathit{nhi}$.

::: details
Operand coefficients that do not contribute to the requested range are
dropped first. The cost is then about that of a product of lengths
$\min(\mathit{len1}, \mathit{len2})$ and $\mathit{nhi} - \mathit{nlo}$:
the classical algorithm computes only the requested coefficients,
`fft_small` computes them either as a transposed product of those
lengths or as a window of the whole product, using a cyclic convolution
when a power of two allows it, whichever has the shorter transform, and
Kronecker substitution computes the
whole product of the trimmed operands. See [](algorithms/polynomials/middle-product).
:::
:::

### Classical multiplication {#mul-classical}

::: function
void _nmod_poly_mul_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
void nmod_poly_mul_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

Set `res` to the product by the classical algorithm, with
$O(\mathit{len1} \cdot \mathit{len2})$ operations. The underscore version
requires $\mathit{len1} \ge \mathit{len2} \ge 1$ and does not allow
aliasing.

::: details
If $(n-1)^2 \min(\mathit{len1}, \mathit{len2})$ fits in one word, the
products are accumulated without reduction and every coefficient is
reduced once at the end. Otherwise each coefficient is a dot product
computed with `_nmod_vec_dot_rev`, whose method is chosen once for the
longest dot product. Squaring uses the symmetry of the product and does
about half the multiplications.
:::
:::

::: function
void _nmod_poly_mullow_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong trunc, nmod_t mod)
void nmod_poly_mullow_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc)

Set `res` to the product modulo $x^{\mathit{trunc}}$ by the classical
algorithm. The underscore version requires $\mathit{len1} \ge
\mathit{len2} \ge 1$ and $1 \le \mathit{trunc} \le \mathit{len1} +
\mathit{len2} - 1$, and does not allow aliasing.
:::

::: function
void _nmod_poly_mulhigh_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong start, nmod_t mod)
void nmod_poly_mulhigh_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong start)

Set the coefficients of `res` from index `start` on to those of the
product, by the classical algorithm, and the coefficients below `start` to
zero. The underscore version requires $\mathit{len1} \ge \mathit{len2}
\ge 1$, writes $\mathit{len1} + \mathit{len2} - 1$ coefficients and does
not allow aliasing.
:::

::: function
void _nmod_poly_mulmid_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
void nmod_poly_mulmid_classical(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong nlo, slong nhi)

Set `res` to the coefficients $\mathit{nlo}, \ldots, \mathit{nhi} - 1$ of
the product by the classical algorithm, with the same conventions as
`_nmod_poly_mulmid` and `nmod_poly_mulmid`.
:::

### Kronecker substitution {#mul-ks}

Kronecker substitution evaluates both polynomials at a power of two
$2^b$, with $b$ large enough that the coefficients of the product do not
overlap, multiplies the resulting integers and reads off the coefficients.
`KS2` evaluates at $2^b$ and $-2^b$ and `KS4` also at $2^{-b}$ and
$-2^{-b}$, which halves, respectively quarters, the size of each integer
product [@Har2009a]. All three variants read their inputs completely before
writing the output, so the output may be the same array as an input.

::: summary algorithms/polynomials/kronecker :::

::: function
void _nmod_poly_mul_KS(nn_ptr out, nn_srcptr in1, slong len1, nn_srcptr in2, slong len2, nmod_t mod)
void nmod_poly_mul_KS(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

Set `res` to the product by Kronecker substitution. The underscore version
requires $\mathit{len1}, \mathit{len2} \ge 1$ in either order.

::: details
The field width is the number of bits of $(n - 1)^2 \min(\mathit{len1},
\mathit{len2})$, a bound for the coefficients of the integer product, so
it can exceed one word; it is at most three words.
:::
:::

::: function
void _nmod_poly_mullow_KS(nn_ptr out, nn_srcptr in1, slong len1, nn_srcptr in2, slong len2, slong n, nmod_t mod)
void nmod_poly_mullow_KS(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)

Set `res` to the product modulo $x^n$ by Kronecker substitution. The
underscore version requires $\mathit{len1}, \mathit{len2} \ge 1$ and $1
\le n \le \mathit{len1} + \mathit{len2} - 1$. It multiplies the operands
truncated to length $n$ in full and unpacks the low $n$ coefficients.
:::

::: function
void _nmod_poly_mulmid_KS(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
void nmod_poly_mulmid_KS(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong nlo, slong nhi)

Set `res` to the coefficients $\mathit{nlo}, \ldots, \mathit{nhi} - 1$ of
the product by Kronecker substitution, with the same conventions as
`_nmod_poly_mulmid` and `nmod_poly_mulmid`, except that the output may
alias the inputs.
:::

::: function
void _nmod_poly_mul_KS2(nn_ptr res, nn_srcptr op1, slong n1, nn_srcptr op2, slong n2, nmod_t mod)
void nmod_poly_mul_KS2(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

Set `res` to the product by Kronecker substitution at $2^b$ and $-2^b$
[@Har2009a]. The underscore version requires $\mathit{n1} \ge \mathit{n2}
\ge 1$.

::: details
Write $f(x) = f_e(x^2) + x f_o(x^2)$ for even and odd parts. The two
integer products $h(2^b)$ and $h(-2^b)$ of length about half that of
`_nmod_poly_mul_KS` give $2 h_e(2^{2b}) = h(2^b) + h(-2^b)$ and $2^{b+1}
h_o(2^{2b}) = h(2^b) - h(-2^b)$, from which the even and odd coefficients
of the product are unpacked separately with `_nmod_poly_KS2_unpack` and
reduced with `_nmod_poly_KS2_reduce`.
:::
:::

::: function
void _nmod_poly_mul_KS4(nn_ptr res, nn_srcptr op1, slong n1, nn_srcptr op2, slong n2, nmod_t mod)
void nmod_poly_mul_KS4(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

Set `res` to the product by Kronecker substitution at $\pm 2^b$ and at the
reciprocal points $\pm 2^{-b}$ [@Har2009a]. The underscore version requires
$\mathit{n1} \ge \mathit{n2} \ge 1$.

::: details
The evaluations at $\pm 2^{-b}$ are evaluations of the reversed
polynomials. They give the high halves of the coefficients of the product,
the evaluations at $\pm 2^b$ the low halves, so each of the four integer
products needs fields of only about half the width;
`_nmod_poly_KS2_recover_reduce` combines the two halves and reduces them.
:::
:::

### FFT multiplication {#mul-fft}

These functions call the [small-prime FFT](algorithms/fft/small-prime-fft)
of `fft_small` ([](fft_small)), which multiplies modulo one to several
50-bit primes and reconstructs the coefficients modulo $n$ by the Chinese
remainder theorem ([](algorithms/fft/multimodular)). They exist only when
FLINT is built with `fft_small` (`FLINT_HAVE_FFT_SMALL`); otherwise
`_nmod_poly_mulmid_fft_small` and `_nmod_poly_mullow_fft_small` throw an
exception and the `want` functions return 0.

::: function
void _nmod_poly_mulmid_fft_small(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong znlo, slong zn, nmod_t mod)
void _nmod_poly_mullow_fft_small(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zn, nmod_t mod)

Write the coefficients $\mathit{znlo}, \ldots, \mathit{zn} - 1$
(respectively $0, \ldots, \mathit{zn} - 1$) of the product of `(a, an)`
and `(b, bn)` to `z`, using `fft_small`. Require $\mathit{an},
\mathit{bn} \ge 1$, in either order, and $0 \le \mathit{znlo} <
\mathit{zn} \le \mathit{an} + \mathit{bn} - 1$; the output must not
overlap the inputs.

::: details
`_nmod_poly_mullow_fft_small` is `_nmod_poly_mulmid_fft_small` with
$\mathit{znlo} = 0$. The coefficients of the operands that no output
coefficient uses (those of index at least $\mathit{zn}$, and the low
coefficients of each operand that meet the other only below
$\mathit{znlo}$) are dropped first. For $n \le 23$ they first try `_nmod_poly_mulmid_fft_small_repack`, which
packs two coefficients into one and multiplies modulo a single 50-bit
prime at half the length. Otherwise `_nmod_poly_mulmid_fft_small`
compares two ways of computing the window, with $s$ and $\ell$ the
lengths of the shorter and the longer operand and $m = \mathit{zn} -
\mathit{znlo}$: as a window of the full product
(`_nmod_poly_mul_mid_default_mpn_ctx`), with a cyclic convolution when a
power of two allows it, or as a transposed product of the reversed
shorter operand and the longer one, shifted and zero-padded, which costs
a product of lengths $s$ and $m$
(`_nmod_poly_mul_tr_default_mpn_ctx`). The transposed product is used if
$m < \ell$, the window does not extend beyond the full sums on both sides
($\mathit{znlo} \ge s - 1$ or $\mathit{zn} \le \ell$), and its transform is
strictly shorter. After the trimming this happens
for middle products such as those of Newton iteration, but never for low
products, whose output is at least as long as the longer operand. See
[](algorithms/polynomials/middle-product#transposed-fft).
:::
:::

::: function
int _nmod_poly_mullow_want_fft_small(slong len1, slong len2, slong n, int squaring, nmod_t mod)
int _nmod_poly_mulmid_want_fft_small(slong len2, slong m, nmod_t mod)

Return whether `fft_small` is expected to be the fastest algorithm for a
product of lengths `len1` and `len2` truncated to length $n$ (a full
product if $n = \mathit{len1} + \mathit{len2} - 1$), respectively for $m$
coefficients of a middle product whose shorter operand, after trimming,
has length `len2`.

::: details
For full products the condition is $\min(\mathit{len1}, 2 \mathit{len2})$
at least a tabulated cutoff, one for each bit size of $n$ and separate
ones for squaring, between about 80 and 1400. Moduli that are themselves
one of the primes of the default `fft_small` context need no Chinese
remaindering and use `fft_small` from $\min(\mathit{len1}, 2
\mathit{len2}) \ge 24$. For $n \le 23$ (repacking) the full product
length must exceed 370, with gaps just above the power-of-two lengths 512
and 1024 where the transform length doubles.
:::
:::

::: function
int _nmod_poly_mulmid_fft_small_repack(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong znlo, slong zn, nmod_t mod)

Internal helper of `_nmod_poly_mulmid_fft_small` for $n \le 23$: if the
coefficients are small enough, writes the coefficients $\mathit{znlo},
\ldots, \mathit{zn} - 1$ of the product to `z` and returns 1; otherwise
returns 0 without writing anything. The conditions on the arguments are
those of `_nmod_poly_mulmid_fft_small`.

::: details
With $A_2(x) = \sum_k (a_{2k} + a_{2k+1} M) x^k$ and $B_2$ likewise, for
$M = 2^{16}$ or $2^{17}$, the product $A_2 B_2$ computed modulo the prime
$p = 1108307720798209$ has coefficients $c_0 + c_1 M + c_2 M^2$ from which
the coefficients of $AB$ are read off, provided that $c_0, c_1 < M$ and
$c_2 M^2 < p$. These conditions are checked with the bound
$\min(\mathit{an}, \mathit{bn}) (n-1)^2$ or, if that fails, with the
Cauchy–Schwarz inequality applied to the even and odd parts.

::: note
In FLINT 3.6 and earlier,
the trailing coefficient of an operand of odd length entered the
Cauchy–Schwarz bound unsquared, so the check could accept a packing
that overflows: for $n = 23$ and
$a = b = 22 \sum_{i=0}^{270} x^{2i}$, `nmod_poly_mul` returned a wrong
product.
:::
:::
:::

## Powering {#powering}

::: function
void _nmod_poly_pow(nn_ptr res, nn_srcptr poly, slong len, ulong e, nmod_t mod)
void nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e)
void _nmod_poly_pow_binexp(nn_ptr res, nn_srcptr poly, slong len, ulong e, nmod_t mod)
void nmod_poly_pow_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e)

Set `res` to $\mathit{poly}^e$ by binary exponentiation; the two pairs of
functions are the same. The underscore versions require $\mathit{len} \ge
1$ and $e \ge 1$, write $(\mathit{len} - 1) e + 1$ coefficients and do
not allow aliasing. `nmod_poly_pow` returns 1 for $0^0$, but
`nmod_poly_pow_binexp` returns 0.

::: details
A factor $x^k$ of `poly` is removed first and contributes $x^{ke}$. The
squarings and multiplications by `poly` are products truncated to the
length of the result, so the cost is dominated by the last squaring, about
$M((\mathit{len} - 1) e)$.
:::
:::

::: function
void _nmod_poly_pow_trunc(nn_ptr res, nn_srcptr poly, ulong e, slong trunc, nmod_t mod)
void nmod_poly_pow_trunc(nmod_poly_t res, const nmod_poly_t poly, ulong e, slong trunc)
void _nmod_poly_pow_trunc_binexp(nn_ptr res, nn_srcptr poly, ulong e, slong trunc, nmod_t mod)
void nmod_poly_pow_trunc_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e, slong trunc)

Set `res` to $\mathit{poly}^e \bmod x^{\mathit{trunc}}$ by binary
exponentiation with truncated products; the two pairs of functions are the
same. The underscore versions read exactly `trunc` coefficients of `poly`
(zero-pad shorter inputs), require $\mathit{trunc} \ge 1$ and $e \ge 1$,
write `trunc` coefficients and do not allow aliasing.

::: details
For a power series with constant term 1 and $e$ much larger than
`trunc`, $\exp(e \log(\mathit{poly}))$ with `nmod_poly_log_series` and
`nmod_poly_exp_series` can be faster, provided the modulus is a prime
larger than `trunc`.
:::
:::

## Division {#division}

Division with remainder by $B$ requires the leading coefficient of $B$ to
be invertible modulo $n$; otherwise `FLINT_IMPINV` is thrown. Division by
the zero polynomial throws `FLINT_DIVZERO`, except modulo $n = 1$, where
`nmod_poly_divrem`, `nmod_poly_div` and `nmod_poly_divexact` return
$Q = A$ and $R = 0$. Short quotients are computed by schoolbook division;
long ones by [Newton iteration](algorithms/polynomials/newton) on the
reversed polynomials, which costs a constant number of multiplications.
[](nmod_poly/divrem) compares the functions and gives examples.

::: summary algorithms/polynomials/newton :::

::: function
void _nmod_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B)

Set $Q$ and $R$ to the quotient and remainder of $A$ by $B$, so that $A =
BQ + R$ with $\deg R < \deg B$. The underscore version requires
$\mathit{lenA} \ge \mathit{lenB} \ge 1$ and writes $\mathit{lenA} -
\mathit{lenB} + 1$ coefficients to $Q$ and $\mathit{lenB} - 1$ to $R$.
$R$ may be the same array as $A$; otherwise the outputs must not overlap
each other or the inputs.

::: details
Schoolbook division is used if the quotient has at most 80
coefficients or $\mathit{lenB}$ is small (below 400 for moduli of at most
27 bits and below 80 otherwise, in builds with AVX2; 200 and 80 without),
and Newton division otherwise. Without `fft_small`, divisors of length at
least 1024 and moduli of at least 16 bits use divide-and-conquer division
for $\mathit{lenA} \le 16384$.
:::
:::

::: function
void _nmod_poly_div(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_div(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B)

Set $Q$ to the quotient of $A$ by $B$ without computing the remainder.
The underscore version requires $\mathit{lenA} \ge \mathit{lenB} \ge 1$,
writes $\mathit{lenA} - \mathit{lenB} + 1$ coefficients and does not allow
aliasing.

::: details
Schoolbook division is used if $\mathit{lenB} \le 15$ or the quotient has
at most 16 coefficients, and Newton division of the reversed polynomials
otherwise: the quotient only depends on the top $\mathit{lenQ}$
coefficients of $A$ and $B$.
:::
:::

::: function
void _nmod_poly_rem(nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B)

Set $R$ to the remainder of $A$ by $B$. The underscore version requires
$\mathit{lenA} \ge \mathit{lenB} \ge 1$, writes $\mathit{lenB} - 1$
coefficients and does not allow aliasing. Quotients of length 2, the
typical step of the Euclidean algorithm, use a dedicated routine.
`nmod_poly_rem` throws `FLINT_DIVZERO` for $B = 0$ even modulo 1.
:::

::: function
void _nmod_poly_divexact(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_divexact(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B)

Set $Q$ to $A / B$, assuming that $B$ divides $A$; otherwise the result is
meaningless. The underscore version requires $\mathit{lenA} \ge
\mathit{lenB} \ge 1$, writes $\mathit{lenA} - \mathit{lenB} + 1$
coefficients and does not allow aliasing.

::: details
Since the division is exact, the quotient can be computed from the top
coefficients (division of the reversed polynomials) or from the bottom
coefficients (power series division), or half from each end. The
algorithm is schoolbook division from the top if $\mathit{lenB} \le 40$
or $\mathit{lenQ} \le 20$, schoolbook division from both ends if
$\mathit{lenB} \le 60$ or $\mathit{lenQ} \le 60$, and otherwise Newton
division of both halves, each of half the length.
:::
:::

::: function
ulong _nmod_poly_div_root(nn_ptr Q, nn_srcptr A, slong len, ulong c, nmod_t mod)
ulong nmod_poly_div_root(nmod_poly_t Q, const nmod_poly_t A, ulong c)

Set $Q$ to the quotient of $A$ by $x - c$ and return the remainder
$A(c)$, by Horner's rule (synthetic division). The value $c$ must be
reduced. The underscore version requires $\mathit{len} \ge 2$ (for
smaller lengths it returns 0 and writes nothing), writes $\mathit{len} -
1$ coefficients and allows `Q == A`.
:::

::: function
int _nmod_poly_divides(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
int nmod_poly_divides(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B)
int _nmod_poly_divides_classical(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
int nmod_poly_divides_classical(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B)

If $B$ divides $A$, set $Q = A / B$ and return 1; otherwise set $Q = 0$
and return 0. The zero polynomial divides only itself. The underscore
versions require $\mathit{lenA} \ge \mathit{lenB} \ge 1$, write
$\mathit{lenA} - \mathit{lenB} + 1$ coefficients and do not allow
aliasing.

::: details
Both compute the quotient $Q$ with `_nmod_poly_div` first and then check
that $BQ$ agrees with $A$ in the low $\mathit{lenB} - 1$ coefficients,
which is cheaper than computing the remainder. The classical version
checks the coefficients one at a time and stops at the first mismatch;
`_nmod_poly_divides` uses it when $\mathit{lenA} < 40$ and
$\mathit{lenB} < 20$, and otherwise checks the low coefficients in
blocks of length $\mathit{lenQ}$ by fast multiplication (or computes the
remainder if $\mathit{lenA} \ge 2 \mathit{lenB} - 1$).
:::
:::

::: function
ulong nmod_poly_remove(nmod_poly_t f, const nmod_poly_t p)

Removes the largest power $p^k$ dividing `f` from `f` and returns $k$.
Requires $\deg p \ge 1$: for a nonzero constant $p$ the function does not
terminate. For $f = 0$ it returns 0.
:::

::: function
void _nmod_poly_divrem_basecase(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B)

Division with remainder by the schoolbook algorithm, with the same
conventions as `_nmod_poly_divrem` and `nmod_poly_divrem`.

::: details
The partial remainders are kept unreduced in one, two or three words per
coefficient, depending on the number of bits of $(n - 1)^2
(\mathit{lenA} - \mathit{lenB} + 1)$, and each coefficient is reduced only
when it becomes the leading one. Quotients of length 1 and 2 are handled
by special routines.
:::
:::

::: function
void _nmod_poly_div_newton_n_preinv(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr Binv, slong lenBinv, nmod_t mod)
void nmod_poly_div_newton_n_preinv(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv)

Set $Q$ to the quotient of $A$ by $B$, given the power series inverse
`Binv` of the reversal of $B$ (see [](nmod_poly/modular-arithmetic)). The
quotient is the reversal of $\operatorname{rev}(A) \cdot \mathit{Binv}
\bmod x^{\mathit{lenQ}}$, so `Binv` must be correct modulo
$x^{\mathit{lenQ}}$, $\mathit{lenQ} = \mathit{lenA} - \mathit{lenB} + 1$.
The underscore version requires $\mathit{lenA} \ge \mathit{lenB} \ge 1$,
does not read $B$, and does not allow $Q$ to overlap `Binv`.

::: details
`Binv` is usually computed with
`nmod_poly_inv_series(Binv, rev(B), lenB)`, which is enough for
$\mathit{lenA} \le 2 \mathit{lenB} - 1$. A shorter `Binv` is allowed if
its missing coefficients are zero, as after normalisation.

::: warning Spurious message
`nmod_poly_div_newton_n_preinv` prints `Exception
(nmod_poly_div_newton_n_preinv).` to standard output if $\mathit{lenA} >
2 \mathit{lenB} - 2$, but does not stop: the quotient is correct whenever
`Binv` is precise enough, in particular for $\mathit{lenA} = 2
\mathit{lenB} - 1$ and `Binv` of length `lenB`.
:::
:::
:::

::: function
void _nmod_poly_divrem_newton_n_preinv(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr Binv, slong lenBinv, nmod_t mod)
void nmod_poly_divrem_newton_n_preinv(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t Binv)

Set $Q$ and $R$ to the quotient and remainder of $A$ by $B$, given `Binv`
as for `nmod_poly_div_newton_n_preinv`; the remainder is $A - BQ$ computed
by a short product. The underscore version requires $\mathit{lenA} \ge
\mathit{lenB} \ge 1$ and does not allow the outputs to overlap the inputs.

::: details
For quotients of length 2 the schoolbook routine is used with
`Binv[0]` $= \operatorname{lc}(B)^{-1}$. For $\mathit{lenB} > 20$ and
divisors with at most five nonzero coefficients below the leading one,
the division is done directly in $O(\mathit{lenA})$ operations, which
makes arithmetic modulo sparse polynomials, such as the trinomials and
pentanomials used to define finite fields, efficient.

::: warning Spurious message
`nmod_poly_divrem_newton_n_preinv` prints `Exception
(nmod_poly_divrem_newton_n_preinv).` to standard output if $\mathit{lenA}
> 2 \mathit{lenB} - 1$, but does not stop; the result is correct whenever
`Binv` is precise enough.
:::
:::
:::

## Power series inversion and division {#series-division}

These functions compute with power series truncated to length $n$, that
is, modulo $x^n$. The constant term of the divisor must be invertible
modulo $n$; otherwise `FLINT_IMPINV` is thrown. For long series they use
[Newton iteration](algorithms/polynomials/newton), which costs a small
constant number of multiplications of length $n$. See
[](nmod_poly/inv_series).

::: function
void _nmod_poly_inv_series(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod)
void nmod_poly_inv_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)

Set `Qinv` to the inverse of the power series $Q$ modulo $x^n$, that is,
$Q \cdot \mathit{Qinv} \equiv 1 \pmod{x^n}$. The underscore version
requires $\mathit{Qlen} \ge 1$ and $n \ge 1$, writes $n$ coefficients and
does not allow aliasing; $Q$ may be longer than $n$. The non-underscore
version throws `FLINT_DIVZERO` if $Q = 0$ or $n = 0$.
:::

::: function
void _nmod_poly_inv_series_basecase(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod)
void nmod_poly_inv_series_basecase(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)

The same, by the $O(n \min(n, \mathit{Qlen}))$ recurrence $q_0 = Q_0^{-1}$,
$q_i = -q_0 \sum_{j \ge 1} Q_j q_{i-j}$, with each sum computed as a dot
product. A binomial $Q_0 + c x^k$ is inverted directly as a geometric
series.
:::

::: function
void _nmod_poly_inv_series_newton(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod)
void nmod_poly_inv_series_newton(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)

The same as `_nmod_poly_inv_series` and `nmod_poly_inv_series`, of which
they are aliases; they use Newton iteration only when the length makes it
faster.
:::

::: function
void _nmod_poly_div_series(nn_ptr Q, nn_srcptr A, slong Alen, nn_srcptr B, slong Blen, slong n, nmod_t mod)
void nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B, slong n)

Set $Q$ to $A / B \bmod x^n$, so that $BQ \equiv A \pmod{x^n}$. The
underscore version requires $\mathit{Alen}, \mathit{Blen}, n \ge 1$,
writes $n$ coefficients and does not allow aliasing. The non-underscore
version throws `FLINT_DIVZERO` if $n = 0$, $B = 0$ or $B_0 = 0$.
:::

::: function
void _nmod_poly_div_series_basecase(nn_ptr Q, nn_srcptr A, slong Alen, nn_srcptr B, slong Blen, slong n, nmod_t mod)
void nmod_poly_div_series_basecase(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B, slong n)

The same by the $O(n \min(n, \mathit{Blen}))$ recurrence $q_i = B_0^{-1}
(A_i - \sum_{j \ge 1} B_j q_{i-j})$.
:::

## Arithmetic modulo a polynomial {#modular}

These functions compute in $(\Z/n\Z)[x]/(f)$ for a nonzero $f$ whose
leading coefficient is invertible. Elements are represented by their
remainders, polynomials of length less than $\mathit{lenf}$; the
underscore functions usually require them to be zero-padded to length
exactly $\mathit{lenf} - 1$. The `_preinv` variants take the power series
inverse `finv` of the reversal of $f$, computed once with
`nmod_poly_inv_series(finv, rev(f), lenf)`, which turns every reduction
into two products. The guide [](nmod_poly/modular-arithmetic) explains
these conventions; [](nmod_poly/powmod) covers exponentiation.

::: function
void _nmod_poly_mulmod(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nn_srcptr f, slong lenf, nmod_t mod)
void nmod_poly_mulmod(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f)

Set `res` to $\mathit{poly1} \cdot \mathit{poly2} \bmod f$, computed as a
product followed by `_nmod_poly_divrem`. The underscore version requires
$\mathit{len1}, \mathit{len2} \ge 1$ and $\mathit{len1} + \mathit{len2} -
1 \ge \mathit{lenf}$ (otherwise there is nothing to reduce; use
`_nmod_poly_mul`), writes $\mathit{lenf} - 1$ coefficients and does not
allow `res` to overlap $f$. The non-underscore version accepts inputs of
any length and throws `FLINT_ERROR` for $f = 0$.
:::

::: function
void _nmod_poly_mulmod_preinv(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod)
void nmod_poly_mulmod_preinv(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f, const nmod_poly_t finv)

Set `res` to $\mathit{poly1} \cdot \mathit{poly2} \bmod f$, reducing with
`_nmod_poly_divrem_newton_n_preinv` and the precomputed `finv`. Both
factors must be reduced modulo $f$: the non-underscore version throws
`FLINT_ERROR` if $\mathit{len1} \ge \mathit{lenf}$ or $\mathit{len2} \ge
\mathit{lenf}$. The underscore version also requires $\mathit{len1},
\mathit{len2} \ge 1$ and $\mathit{len1} + \mathit{len2} - 1 \ge
\mathit{lenf}$, writes $\mathit{lenf} - 1$ coefficients and does not
allow `res` to overlap $f$ or `finv`.
:::

::: function
int _nmod_poly_invmod(ulong * A, const ulong * B, slong lenB, const ulong * P, slong lenP, const nmod_t mod)
int nmod_poly_invmod(nmod_poly_t A, const nmod_poly_t B, const nmod_poly_t P)

If $B$ is invertible modulo $P$, set $A$ to its inverse and return 1;
otherwise return 0, leaving $A$ undefined. Requires a prime modulus. The
non-underscore version reduces $B$ modulo $P$ first, returns 0 for
$B = 0$ and throws `FLINT_ERROR` if $\mathit{lenP} < 2$. The underscore
version requires $1 \le \mathit{lenB} < \mathit{lenP}$ (zero padding of
$B$ is allowed), writes $\mathit{lenP} - 1$ coefficients and does not
allow aliasing. It calls the extended Euclidean algorithm
`_nmod_poly_gcdinv`.
:::

::: function
void _nmod_poly_powmod_ui_binexp(nn_ptr res, nn_srcptr poly, ulong e, nn_srcptr f, slong lenf, nmod_t mod)
void nmod_poly_powmod_ui_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e, const nmod_poly_t f)
void _nmod_poly_powmod_fmpz_binexp(nn_ptr res, nn_srcptr poly, fmpz_t e, nn_srcptr f, slong lenf, nmod_t mod)
void nmod_poly_powmod_fmpz_binexp(nmod_poly_t res, const nmod_poly_t poly, fmpz_t e, const nmod_poly_t f)

Set `res` to $\mathit{poly}^e \bmod f$ by left-to-right binary
exponentiation, reducing each product with `_nmod_poly_divrem`. The
non-underscore versions accept any `poly` (reducing it first) and $e \ge
0$, and throw `FLINT_ERROR` for $f = 0$. The underscore versions require
$\mathit{lenf} \ge 2$, `poly` reduced and zero-padded to length
$\mathit{lenf} - 1$, and $e \ge 1$; they write $\mathit{lenf} - 1$
coefficients and do not allow aliasing.
:::

::: function
void _nmod_poly_powmod_ui_binexp_preinv(nn_ptr res, nn_srcptr poly, ulong e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod)
void nmod_poly_powmod_ui_binexp_preinv(nmod_poly_t res, const nmod_poly_t poly, ulong e, const nmod_poly_t f, const nmod_poly_t finv)
void _nmod_poly_powmod_fmpz_binexp_preinv(nn_ptr res, nn_srcptr poly, fmpz_t e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod)
void nmod_poly_powmod_fmpz_binexp_preinv(nmod_poly_t res, const nmod_poly_t poly, fmpz_t e, const nmod_poly_t f, const nmod_poly_t finv)

The same with a precomputed `finv`, reducing with
`_nmod_poly_divrem_newton_n_preinv`. These are the functions to use for
large $f$ or many exponentiations modulo the same $f$.
:::

::: function
void _nmod_poly_powmod_x_ui_preinv(nn_ptr res, ulong e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod)
void nmod_poly_powmod_x_ui_preinv(nmod_poly_t res, ulong e, const nmod_poly_t f, const nmod_poly_t finv)
void _nmod_poly_powmod_x_fmpz_preinv(nn_ptr res, fmpz_t e, nn_srcptr f, slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod)
void nmod_poly_powmod_x_fmpz_preinv(nmod_poly_t res, fmpz_t e, const nmod_poly_t f, const nmod_poly_t finv)

Set `res` to $x^e \bmod f$, given `finv`. They use sliding-window
exponentiation in which a multiplication by $x^w$ is a shift followed by a
reduction, so each window costs one reduction instead of a full
multiplication. The underscore versions require $\mathit{lenf} \ge 3$ and
$e \ge 1$ and write $\mathit{lenf} - 1$ coefficients; the non-underscore
versions accept $e \ge 0$ and any nonzero $f$.
:::

::: function
void _nmod_poly_powers_mod_preinv_naive(nn_ptr * res, nn_srcptr f, slong flen, slong n, nn_srcptr g, slong glen, nn_srcptr ginv, slong ginvlen, const nmod_t mod)
void nmod_poly_powers_mod_naive(nmod_poly_struct * res, const nmod_poly_t f, slong n, const nmod_poly_t g)

Set `res[i]` to $f^i \bmod g$ for $0 \le i < n$, computing $f^i$ as
$f^{\lfloor i/2 \rfloor} f^{\lceil i/2 \rceil}$. The non-underscore
version takes an array of $n$ initialised polynomials, reduces $f$ first,
and throws `FLINT_ERROR` for $g = 0$; for $f = 0$ or constant $g$ it sets
`res[0]` to 1 and the others to 0. The underscore version requires $1 \le
\mathit{flen} < \mathit{glen}$, `ginv` as for `_nmod_poly_mulmod_preinv`,
and $n$ arrays of $\mathit{glen} - 1$ coefficients, all of which are
written.
:::

::: function
void _nmod_poly_powers_mod_preinv_threaded_pool(nn_ptr * res, nn_srcptr f, slong flen, slong n, nn_srcptr g, slong glen, nn_srcptr ginv, slong ginvlen, const nmod_t mod, thread_pool_handle * threads, slong num_threads)
void _nmod_poly_powers_mod_preinv_threaded(nn_ptr * res, nn_srcptr f, slong flen, slong n, nn_srcptr g, slong glen, nn_srcptr ginv, slong ginvlen, const nmod_t mod)
void nmod_poly_powers_mod_bsgs(nmod_poly_struct * res, const nmod_poly_t f, slong n, const nmod_poly_t g)

The same as the `naive` functions, computed with baby steps and giant
steps in parallel: with $k = \lfloor \sqrt{n} \rfloor$, the powers
$f^0, \ldots, f^k$ and $f^{2k}, f^{3k}, \ldots$ are computed first, and
then the threads fill in $f^{ik + j} = f^{ik} f^j$, taking blocks of $k$
powers from a shared counter. `_nmod_poly_powers_mod_preinv_threaded`
and `nmod_poly_powers_mod_bsgs` take their threads from the global pool
(see [](thread_pool)), the `_pool` version uses the `num_threads` threads
given.
:::

## Preconditioned modular multiplication {#precond}

When many products $ab \bmod d$ are needed with the same $a$ and $d$,
data depending on $a$ can be precomputed so that each product is cheaper
than `nmod_poly_mulmod_preinv`. Two methods are available. The *matrix*
method stores the $n \times n$ matrix of multiplication by $a$, whose
column $i$ is $a x^i \bmod d$, $n = \deg d$, and computes $ab \bmod d$ as
a matrix-vector product; for small moduli several entries are packed into
one word. It costs $O(n^2)$ but is very fast for $n$ up to a few hundred.
The *Shoup* method stores the quotient $\lfloor a x^n / d \rfloor$; then
the quotient of $ab$ by $d$ is the high part of a short product, and $ab
\bmod d$ costs three short products instead of four, by the polynomial
analogue of [Shoup multiplication](algorithms/arithmetic/shoup)
[@Sho1995].

::: type
nmod_poly_mulmod_precond_struct
nmod_poly_mulmod_precond_t

Precomputed data for multiplication by a fixed $a$ modulo a fixed $d$.
It stores shallow references to $a$, $d$ and the inverse of $d$, which
must stay alive and unchanged while it is used.
:::

::: macro
NMOD_POLY_MULMOD_PRECOND_NONE
NMOD_POLY_MULMOD_PRECOND_SHOUP
NMOD_POLY_MULMOD_PRECOND_MATRIX

The methods: no precomputation (`_nmod_poly_mulmod_preinv` is called),
Shoup's method, and the matrix method.
:::

::: function
void _nmod_poly_mulmod_precond_init_method(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, nn_srcptr dinv, slong lendinv, int method, nmod_t mod)
void nmod_poly_mulmod_precond_init_method(nmod_poly_mulmod_precond_t precond, const nmod_poly_t a, const nmod_poly_t d, const nmod_poly_t dinv, int method)

Initialise `precond` for multiplication by $a$ modulo $d$ with the given
method. Require $a$ reduced modulo $d$ ($\mathit{alen} < \mathit{dlen}$)
and `dinv` the power series inverse of the reversal of $d$, as for
`nmod_poly_mulmod_preinv`; the matrix method only reads `dinv[0]`.
An invalid method throws `FLINT_ERROR`.
:::

::: function
void _nmod_poly_mulmod_precond_init_num(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, nn_srcptr dinv, slong lendinv, slong num, nmod_t mod)
void nmod_poly_mulmod_precond_init_num(nmod_poly_mulmod_precond_t precond, const nmod_poly_t a, const nmod_poly_t d, const nmod_poly_t dinv, slong num)

Like the `init_method` functions, but choose the method expected to be
fastest for `num` multiplications.

::: details
With $n = \deg d$: no precomputation for $n = 1$; the matrix method for
$n \le 16$ and $\mathit{num} \ge 2$. For $d$ with at most four nonzero
coefficients below the leading one (where reduction is cheap anyway), the
matrix method for $n \le 32$ and $\mathit{num} \ge 8$ or $n \le 64$ and
$\mathit{num} \ge 16$, and otherwise no precomputation. For other $d$, the
matrix method for $n \le 64$ and $\mathit{num} \ge 4$, $n \le 128$ and
$\mathit{num} \ge 8$, $n \le 192$ and $\mathit{num} \ge 32$, or $n \le
320$, $\mathit{num} \ge 64$ and moduli below 256; otherwise Shoup's method
if $\mathit{num} \ge 8$, and no precomputation below.
:::
:::

::: function
void nmod_poly_mulmod_precond_clear(nmod_poly_mulmod_precond_t precond)

Frees the memory used by `precond`.
:::

::: function
void _nmod_poly_mulmod_precond(nn_ptr res, const nmod_poly_mulmod_precond_t precond, nn_srcptr b, slong blen, nmod_t mod)
void nmod_poly_mulmod_precond(nmod_poly_t res, const nmod_poly_mulmod_precond_t precond, const nmod_poly_t b)

Set `res` to $ab \bmod d$ for the $a$ and $d$ of `precond`. Require $b$
reduced modulo $d$. The underscore version requires $\mathit{alen},
\mathit{blen} \ge 1$, writes $\deg d$ coefficients (zero-padded) and does
not allow `res` to overlap $b$, $a$ or $d$. The non-underscore version
allows `res` to be $b$, but not $a$.
:::

## Greatest common divisor {#gcd}

All GCD functions require a prime modulus (not checked); a non-invertible
leading coefficient met by the Euclidean or the half-GCD algorithm throws
`FLINT_IMPINV`.
The non-underscore functions return the monic GCD, with $\gcd(0, 0) = 0$,
and accept their inputs in any order. Below a cutoff that depends on the
bit size of $n$ the Euclidean algorithm is used, for moduli of $\bits/2$ to
$\bits - 2$ bits and $\mathit{lenB} \ge 32$ in
[Montgomery arithmetic](algorithms/arithmetic/montgomery); above it the
[half-GCD algorithm](algorithms/polynomials/half-gcd) [@ThullYap1990]
computes the GCD in $O(M(n) \log n)$ operations. See [](nmod_poly/gcd).

::: summary algorithms/polynomials/half-gcd :::

::: function
slong _nmod_poly_gcd(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_gcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B)

Set $G$ to a greatest common divisor of $A$ and $B$. The underscore
version requires $\mathit{lenA} \ge \mathit{lenB} \ge 1$, writes at most
$\mathit{lenB}$ coefficients, returns the length of $G$, does not make $G$
monic and does not allow aliasing.

::: details
`_nmod_poly_gcd` uses the half-GCD algorithm if $\mathit{lenB}$ is at
least `nmod_poly_gcd_hgcd_cutoff(mod)`, and the Euclidean algorithm
otherwise, with Montgomery arithmetic when
`NMOD_POLY_GCD_EUCLIDEAN_USE_REDC_FAST` holds.
:::
:::

::: function
slong _nmod_poly_gcd_euclidean(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_gcd_euclidean(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B)

The GCD by the Euclidean algorithm, with the conventions of
`_nmod_poly_gcd` and `nmod_poly_gcd`. Each step is a call to
`_nmod_poly_rem`, so the cost is $O(\mathit{lenA} \cdot \mathit{lenB})$.
:::

::: function
slong _nmod_poly_gcd_euclidean_redc_fast(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_gcd_euclidean_redc_fast(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B)

The Euclidean algorithm in Montgomery arithmetic (`nmod_redc_fast`, see
[](nmod)), after a first remainder in the ordinary representation.
Require an odd prime modulus $n < 2^{\bits - 2}$, which is not checked.
:::

::: macro
NMOD_POLY_GCD_EUCLIDEAN_USE_REDC_FAST(lenB, mod)

The condition under which `_nmod_poly_gcd` and `_nmod_poly_resultant` use
the Montgomery version of the Euclidean algorithm: $\mathit{lenB} \ge 32$
and $\bits/2 \le$ `NMOD_BITS(mod)` $\le \bits - 2$. It does not check that
$n$ is odd, so an even modulus in this range gives an exception or a
crash; such moduli are not prime and therefore not supported anyway.
:::

::: function
slong _nmod_poly_gcd_hgcd(nn_ptr G, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_gcd_hgcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B)

The GCD by the half-GCD algorithm, with the conventions of
`_nmod_poly_gcd` and `nmod_poly_gcd`. Its complexity is $O(M(n) \log
n)$, where $M(n)$ is the cost of a multiplication of length $n$.

::: details
The outer loop calls `_nmod_poly_hgcd` to reduce the operands to about
half their length, and finishes with the Euclidean algorithm once they are
shorter than `nmod_poly_hgcd_outer_cutoff(mod)`.
:::
:::

::: function
slong _nmod_poly_hgcd(nn_ptr * M, slong * lenM, nn_ptr A, slong * lenA, nn_ptr B, slong * lenB, nn_srcptr a, slong lena, nn_srcptr b, slong lenb, nmod_t mod)
slong nmod_poly_hgcd(nmod_poly_t m11, nmod_poly_t m12, nmod_poly_t m21, nmod_poly_t m22, nmod_poly_t A, nmod_poly_t B, const nmod_poly_t a, const nmod_poly_t b)

The half-GCD step. For $\deg a > \deg b$, compute consecutive remainders
$A, B$ of the Euclidean remainder sequence of $a, b$ with $\deg A \ge
\deg(a)/2 > \deg B$ and the matrix
$$
M = \begin{pmatrix} m_{11} & m_{12} \\ m_{21} & m_{22} \end{pmatrix},
\qquad
\begin{pmatrix} a \\ b \end{pmatrix} = M \begin{pmatrix} A \\ B
\end{pmatrix},
$$
the product of the matrices $\left(\begin{smallmatrix} q_i & 1 \\ 1 & 0
\end{smallmatrix}\right)$ for the quotients $q_i$ of these steps. Return
$\det M = \pm 1$.

::: details
`nmod_poly_hgcd` throws `FLINT_ERROR` unless $\deg a > \deg b$; for $b =
0$ it returns the identity matrix. The underscore version requires
$\mathit{lena} > \mathit{lenb} \ge 1$; `A` and `B` and the four arrays
`M[0]`, ..., `M[3]` (for $m_{11}, m_{12}, m_{21}, m_{22}$) need space for
`lena` coefficients, and the actual lengths are written to `lenA`,
`lenB` and `lenM`. It allows `A == a` and `B == b`. Below
`nmod_poly_hgcd_iter_recursive_cutoff(mod)` the recursion is replaced by
Euclidean steps. See [](algorithms/polynomials/half-gcd).
:::
:::

::: function
slong nmod_poly_hgcd_ref(nmod_poly_t m11, nmod_poly_t m12, nmod_poly_t m21, nmod_poly_t m22, nmod_poly_t A, nmod_poly_t B, const nmod_poly_t a, const nmod_poly_t b)

A reference implementation of `nmod_poly_hgcd` by repeated division with
remainder, for testing. Same contract; no aliasing is allowed.
:::

::: function
slong nmod_poly_hgcd_iter_recursive_cutoff(nmod_t mod)
slong nmod_poly_hgcd_outer_cutoff(nmod_t mod)
slong nmod_poly_gcd_hgcd_cutoff(nmod_t mod)
slong nmod_poly_xgcd_hgcd_cutoff(nmod_t mod)

Return the tuned cutoffs for the bit size of the modulus: the length below
which the half-GCD recursion uses Euclidean steps (60 to 244), the length
below which the outer loop of the half-GCD GCD switches to the Euclidean
algorithm (310 to 1585), and the length of $B$ from which `_nmod_poly_gcd`
(470 to 2199) and `_nmod_poly_xgcd` (72 to 388) use the half-GCD algorithm.
:::

::: function
slong _nmod_poly_xgcd(nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, const nmod_poly_t A, const nmod_poly_t B)

Compute the GCD $G$ of $A$ and $B$ and cofactors $S, T$ with $SA + TB =
G$. The non-underscore version makes $G$ monic; for $B \ne 0$ the
cofactors satisfy $\deg S < \max(\deg B - \deg G, 1)$ and $\deg T <
\max(\deg A - \deg G, 1)$, and for $B = 0$ it returns $S =
\operatorname{lc}(A)^{-1}$, $T = 0$. $G$, $S$ and $T$ must be distinct
objects.

::: details
The underscore version requires $\mathit{lenA} \ge \mathit{lenB} \ge 1$
and space for $\mathit{lenB}$ coefficients in $G$, $\max(\mathit{lenB} -
1, 2)$ in $S$ and $\max(\mathit{lenA} - 1, 2)$ in $T$, all of which may be
written; it returns the length of $G$, does not make $G$ monic and does
not allow aliasing. It uses the half-GCD algorithm if $\mathit{lenB} \ge$
`nmod_poly_xgcd_hgcd_cutoff(mod)` and the Euclidean algorithm otherwise.
:::
:::

::: function
slong _nmod_poly_xgcd_euclidean(nn_ptr res, nn_ptr s, nn_ptr t, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
void nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, const nmod_poly_t A, const nmod_poly_t B)
slong _nmod_poly_xgcd_hgcd(nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
void nmod_poly_xgcd_hgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, const nmod_poly_t A, const nmod_poly_t B)

The extended GCD by the Euclidean and by the half-GCD algorithm, with the
conventions of `_nmod_poly_xgcd` and `nmod_poly_xgcd`.
:::

::: function
slong _nmod_poly_gcdinv(ulong * G, ulong * S, const ulong * A, slong lenA, const ulong * B, slong lenB, const nmod_t mod)
void nmod_poly_gcdinv(nmod_poly_t G, nmod_poly_t S, const nmod_poly_t A, const nmod_poly_t B)

Compute $G = \gcd(A, B)$ and $S$ with $G \equiv SA \pmod B$ and $\deg S <
\deg B - \deg G$, without the cofactor of $B$. If $\gcd(A, B) = 1$, then
$S = A^{-1} \bmod B$. The non-underscore version requires $\mathit{lenB}
\ge 2$ (otherwise it throws `FLINT_ERROR`), reduces $A$ modulo $B$ first,
makes $G$ monic and returns $G = S = 0$ if $B$ divides $A$. The
underscore version requires $1 \le \mathit{lenA} < \mathit{lenB}$, writes
at most `lenA` coefficients to $G$ and $\mathit{lenB} - 1$ to $S$, returns
the length of $G$ and does not make $G$ monic.
:::

## Resultant and discriminant {#resultant}

The resultant of $f = a \prod_{i=1}^{m} (x - \alpha_i)$ and $g = b
\prod_{j=1}^{k} (x - \beta_j)$ is
$$
\operatorname{res}(f, g) = a^{k} b^{m} \prod_{i,j} (\alpha_i - \beta_j)
= a^{k} \prod_{i} g(\alpha_i),
$$
the determinant of the Sylvester matrix. It is zero exactly when $f$ and
$g$ have a common factor. These functions require a prime modulus (not
checked) and return 0 if either polynomial is zero. See
[](nmod_poly/resultant).

::: function
ulong _nmod_poly_resultant(nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
ulong nmod_poly_resultant(const nmod_poly_t f, const nmod_poly_t g)

Return $\operatorname{res}(f, g)$. The underscore version requires
$\mathit{len1} \ge \mathit{len2} \ge 1$.

::: details
Small cases are handled by explicit formulas. Otherwise the Euclidean
algorithm is used below `nmod_poly_gcd_hgcd_cutoff(mod)`, in Montgomery
arithmetic when `NMOD_POLY_GCD_EUCLIDEAN_USE_REDC_FAST` holds and
$\mathit{len1} < 8 \mathit{len2}$, and the half-GCD algorithm above.
If both arguments are the same polynomial of positive degree, the result
is 0 at once.
:::
:::

::: function
ulong _nmod_poly_resultant_euclidean(nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
ulong nmod_poly_resultant_euclidean(const nmod_poly_t f, const nmod_poly_t g)

The resultant by the Euclidean algorithm, using $\operatorname{res}(f, g)
= (-1)^{\deg f \deg g} \operatorname{lc}(g)^{\deg f - \deg r}
\operatorname{res}(g, r)$ for $r = f \bmod g$.
:::

::: function
ulong _nmod_poly_resultant_hgcd(nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nmod_t mod)
ulong nmod_poly_resultant_hgcd(const nmod_poly_t A, const nmod_poly_t B)

The resultant by the half-GCD algorithm, in $O(M(n) \log n)$ operations.

::: details
The half-GCD algorithm performs the same divisions as the Euclidean
algorithm, but the remainders of a recursive call are only known in
their top coefficients until the transformation matrix is applied. The
factors $\operatorname{lc}^{\deg - \deg}$ and signs contributed by each
division are therefore recorded when the division is performed and
applied once the remainder is complete.
:::
:::

::: function
ulong _nmod_poly_discriminant(nn_srcptr poly, slong len, nmod_t mod)
ulong nmod_poly_discriminant(const nmod_poly_t f)

Return the discriminant of $f$ of degree $d$,
$$
\operatorname{disc}(f) = (-1)^{d(d-1)/2} \operatorname{lc}(f)^{-1}
\operatorname{res}(f, f') = \operatorname{lc}(f)^{2d - 2} \prod_{i < j}
(\alpha_i - \alpha_j)^2,
$$
where the power of $\operatorname{lc}(f)$ is adjusted if $\deg f' < d -
1$ (when the characteristic divides $d$). The discriminant is 0 for
constants and when $f' = 0$. The underscore version requires
$\mathit{len} \ge 2$.
:::

## Derivative and integral {#calculus}

::: function
void _nmod_poly_derivative(nn_ptr x_prime, nn_srcptr x, slong len, nmod_t mod)
void nmod_poly_derivative(nmod_poly_t x_prime, const nmod_poly_t x)

Set `x_prime` to the derivative of `x`. The underscore version requires
$\mathit{len} \ge 1$, writes $\mathit{len} - 1$ coefficients and allows
`x_prime == x`.
:::

::: function
void _nmod_poly_integral(nn_ptr x_int, nn_srcptr x, slong len, nmod_t mod)
void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x)

Set `x_int` to the integral of `x` with constant term 0. The integers
$1, \ldots, \deg x + 1$ must be invertible modulo $n$, for example $n$
prime and $n > \deg x + 1$; otherwise `FLINT_IMPINV` is thrown. Here
`len` is the length of the *output*: the underscore version reads
$\mathit{len} - 1$ coefficients, writes `len` coefficients, requires
$\mathit{len} \ge 1$ and allows `x_int == x`. It needs only one modular
inversion.
:::

## Evaluation {#evaluation}

::: function
ulong _nmod_poly_evaluate_nmod(nn_srcptr poly, slong len, ulong c, nmod_t mod)
ulong nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c)

Return $\mathit{poly}(c)$ for reduced $c$. Any length, including 0, is
allowed.

::: details
For moduli below $2^{\bits - 1}$, lengths below
`FLINT_MULMOD_SHOUP_THRESHOLD` (23 on common platforms) use Horner's rule
with `nmod_mul`, lengths up to 50 Horner's rule with
[Shoup multiplication](algorithms/arithmetic/shoup) by the fixed $c$
(with lazy reduction if $3n - 1 \le 2^{\bits}$), and longer polynomials
rectangular splitting. For larger moduli, Horner's rule is used up to
length 13 and rectangular splitting beyond.
:::
:::

::: function
ulong _nmod_poly_evaluate_nmod_horner(nn_srcptr poly, slong len, ulong c, nmod_t mod)
ulong _nmod_poly_evaluate_nmod_rectangular(nn_srcptr poly, slong len, ulong c, nmod_t mod)

Return $\mathit{poly}(c)$ by Horner's rule, and by rectangular splitting:
with $m = \min(\lfloor\sqrt{\mathit{len}}\rfloor + 1, 1024)$, the powers
$c, \ldots, c^m$ are computed once, each block of $m$ coefficients is
evaluated by a dot product (see [](nmod_vec/dot)), and the blocks are
combined by Horner's rule in $c^m$. The dot products accumulate without
intermediate reductions, which makes rectangular splitting faster than
Horner's rule for long polynomials although it does the same number of
multiplications.
:::

::: function
ulong _nmod_poly_evaluate_nmod_precomp(nn_srcptr poly, slong len, ulong c, ulong c_precomp, ulong modn)
ulong _nmod_poly_evaluate_nmod_precomp_lazy(nn_srcptr poly, slong len, ulong c, ulong c_precomp, ulong modn)

Evaluate `poly` at $c$ modulo `modn` by Horner's rule with Shoup
multiplication, given `c_precomp = n_mulmod_precomp_shoup(c, modn)`.
Require $c < \mathit{modn} < 2^{\bits - 1}$.

::: details
The `lazy` version does not reduce the additions: if all coefficients are
less than $m$ and $m \le 2^{\bits} - 2\,\mathit{modn} + 1$, it returns a
value in $[0, m + 2\,\mathit{modn} - 1)$ congruent to the result. The
coefficients need not be reduced. With reduced coefficients ($m =
\mathit{modn}$) this requires $3\,\mathit{modn} - 1 \le 2^{\bits}$, and
the result is reduced by subtracting $\mathit{modn}$ or
$2\,\mathit{modn}$.
:::
:::

::: function
void nmod_poly_evaluate_mat(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
void nmod_poly_evaluate_mat_horner(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
void nmod_poly_evaluate_mat_paterson_stockmeyer(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)

Set `dest` to the matrix $\mathit{poly}(c)$, for a square matrix $c$ with
the modulus of `poly`; `dest` must have the same dimensions and may be
`c`. The Horner
version uses $\deg(\mathit{poly})$ matrix multiplications, the
Paterson–Stockmeyer version [@PS1973] about $2\sqrt{\deg}$ matrix
multiplications and $O(\deg)$ scalar multiply-adds. `nmod_poly_evaluate_mat`
uses Horner's rule if $\mathit{len} < 5$ or $\mathit{len} \cdot r < 425$
for an $r \times r$ matrix, and Paterson–Stockmeyer otherwise.

::: details
The characteristic and minimal polynomials of $c$, which vanish at $c$,
are computed by `nmod_mat_charpoly` and
`nmod_mat_minpoly` ([](nmod_mat/charpoly)).
:::
:::

::: function
void nmod_mat_one_addmul(nmod_mat_t dest, const nmod_mat_t mat, ulong c)

Sets `dest` to $\mathit{mat} + cI$ for a square matrix `mat` and reduced
$c$, the step of Horner's rule above. `dest` may be `mat`.
:::

## Multipoint evaluation and interpolation {#multipoint}

Evaluating a polynomial of length $n$ at $n$ points, or interpolating
through $n$ points, takes $O(n^2)$ operations with Horner's rule or
Lagrange interpolation, and $O(M(n) \log n)$ operations with a
*subproduct tree*: the products $\prod (x - x_i)$ over the leaves of a
binary tree on the points. Evaluation reduces the polynomial down the
tree, interpolation combines the values up the tree. Interpolation needs
the $x_i$ to be distinct, with invertible differences (for example, $n$
prime). See [](nmod_poly/multipoint).

::: summary algorithms/polynomials/subproduct-tree :::

::: function
void _nmod_poly_evaluate_nmod_vec(nn_ptr ys, nn_srcptr coeffs, slong ilen, nn_srcptr xs, slong olen, nmod_t mod)
void nmod_poly_evaluate_nmod_vec(nn_ptr ys, const nmod_poly_t poly, nn_srcptr xs, slong olen)

Set `ys[i]` to $\mathit{poly}(\mathit{xs}[i])$ for $0 \le i < \mathit{olen}$,
for reduced points. Polynomials of length less than 32 are evaluated
point by point with `_nmod_poly_evaluate_nmod`, longer ones with a
subproduct tree.
:::

::: function
void _nmod_poly_evaluate_nmod_vec_iter(nn_ptr ys, nn_srcptr coeffs, slong ilen, nn_srcptr xs, slong olen, nmod_t mod)
void nmod_poly_evaluate_nmod_vec_iter(nn_ptr ys, const nmod_poly_t poly, nn_srcptr xs, slong olen)

The same, evaluating at each point with `_nmod_poly_evaluate_nmod`.
:::

::: function
void _nmod_poly_evaluate_nmod_vec_fast(nn_ptr ys, nn_srcptr coeffs, slong ilen, nn_srcptr xs, slong olen, nmod_t mod)
void nmod_poly_evaluate_nmod_vec_fast(nn_ptr ys, const nmod_poly_t poly, nn_srcptr xs, slong olen)

The same with a temporary subproduct tree.
:::

::: function
void _nmod_poly_evaluate_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly, slong ilen, const nn_ptr * tree, slong olen, nmod_t mod)

Evaluates `(poly, ilen)` at the `olen` points of the subproduct tree
`tree`, built with `_nmod_poly_tree_build`. The polynomial may be longer
or shorter than the number of points; it is first reduced modulo the
products at the highest useful level of the tree, and the remainders are
then reduced down the tree.
:::

::: function
nn_ptr * _nmod_poly_tree_alloc(slong len)
void _nmod_poly_tree_free(nn_ptr * tree, slong len)

Allocate and free a subproduct tree for `len` points. For `len` = 0,
`_nmod_poly_tree_alloc` returns `NULL`.
:::

::: function
void _nmod_poly_tree_build(nn_ptr * tree, nn_srcptr roots, slong len, nmod_t mod)

Builds the subproduct tree of the monic linear factors $x - r_i$, $0 \le i
< \mathit{len}$, in the space allocated by `_nmod_poly_tree_alloc`.

::: details
Level $i$, for $0 \le i < \lceil \log_2 \mathit{len} \rceil$, stores the
monic products of $2^i$ consecutive factors one after the other, each with
its $2^i + 1$ coefficients, followed by the product of the remaining
factors if `len` is not a multiple of $2^i$. The product of all factors
is not computed. For $\mathit{len} = 7$:

```
level 0:  x1 x1 x1 x1 x1 x1 x1     (x - r_i), two coefficients each
level 1:  xx1   xx1   xx1   x1
level 2:  xxxx1       xx1   x1
```
:::
:::

::: function
void _nmod_poly_interpolation_weights(nn_ptr w, const nn_ptr * tree, slong len, nmod_t mod)

Sets `w[i]` to the barycentric weight $1 / \prod_{j \ne i} (x_i - x_j) =
1 / P'(x_i)$, where $P = \prod_j (x - x_j)$ for the points of `tree`.
:::

::: function
void _nmod_poly_interpolate_nmod_vec_fast_precomp(nn_ptr poly, nn_srcptr ys, const nn_ptr * tree, nn_srcptr weights, slong len, nmod_t mod)

Sets `(poly, len)` to the polynomial of length at most `len` taking the
values `ys` at the points of `tree`, given the weights from
`_nmod_poly_interpolation_weights`, by combining the scaled values
$y_i w_i$ up the tree (fast Lagrange interpolation). Writes `len`
coefficients.
:::

::: function
void _nmod_poly_interpolate_nmod_vec(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong len, nmod_t mod)
void nmod_poly_interpolate_nmod_vec(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong len)

Set `poly` to the unique polynomial of length at most `len` with
$\mathit{poly}(x_i) = y_i$ for $0 \le i < \mathit{len}$. The points must be
distinct and all values reduced. Uses the Newton form for $\mathit{len} <
6$, the barycentric form for $\mathit{len} < 16$ and fast interpolation
otherwise. The underscore version requires $\mathit{len} \ge 1$, writes
`len` coefficients and does not allow `poly` to overlap `xs` or `ys`.
:::

::: function
void _nmod_poly_interpolate_nmod_vec_fast(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong len, nmod_t mod)
void nmod_poly_interpolate_nmod_vec_fast(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong len)

The same by fast Lagrange interpolation with a temporary subproduct tree,
in $O(M(\mathit{len}) \log \mathit{len})$ operations.
:::

::: function
void _nmod_poly_interpolate_nmod_vec_newton(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong len, nmod_t mod)
void nmod_poly_interpolate_nmod_vec_newton(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong n)

The same by divided differences, followed by a conversion from the Newton
basis to the monomial basis. It uses $O(\mathit{len}^2)$ multiplications
and $O(\mathit{len}^2)$ modular inversions.
:::

::: function
void _nmod_poly_interpolate_nmod_vec_barycentric(nn_ptr poly, nn_srcptr xs, nn_srcptr ys, slong len, nmod_t mod)
void nmod_poly_interpolate_nmod_vec_barycentric(nmod_poly_t poly, nn_srcptr xs, nn_srcptr ys, slong len)

The same by the Lagrange formula $\sum_i y_i w_i P(x)/(x - x_i)$, with
$P = \prod (x - x_i)$ and $w_i$ as above, in $O(\mathit{len}^2)$
multiplications and `len` inversions.
:::

::: function
void _nmod_poly_product_roots_nmod_vec(nn_ptr poly, nn_srcptr xs, slong n, nmod_t mod)
void nmod_poly_product_roots_nmod_vec(nmod_poly_t poly, nn_srcptr xs, slong n)

Set `poly` to $\prod_{i=0}^{n-1} (x - x_i)$, for reduced $x_i$. The
underscore version writes $n + 1$ coefficients and does not allow
aliasing. Uses the quadratic algorithm for $n < 20$ and balanced
products above.
:::

## Geometric progressions {#geometric}

Evaluation and interpolation at the points of a geometric progression $1,
q, q^2, \ldots$ cost $O(M(n))$ instead of $O(M(n) \log n)$: evaluation is
a single middle product by Bluestein's chirp transform [@Blu1970]
[@RSR1969], using $ij = \binom{i+j}{2} - \binom{i}{2} - \binom{j}{2}$,
and interpolation is two products [@BS2005]. The functions take $r$ with
$q = r^2$, so that $q^{i^2/2} = r^{i^2}$ needs no square root, and use the
points $1, r^2, r^4, \ldots$. The value $r$ must be reduced and
invertible modulo $n$. For interpolation and extrapolation the points must
be distinct with $r^{2k} - 1$ invertible, which for prime $n$ means that
the multiplicative order of $r^2$ exceeds the number of points. See
[](nmod_poly/multipoint#geometric).

::: type
nmod_geometric_progression_struct
nmod_geometric_progression_t

Precomputed data for evaluation, interpolation and extrapolation at the
points $1, r^2, \ldots, r^{2(\mathit{len} - 1)}$: scaling vectors, the
polynomials to multiply by and, optionally, their `fft_small` transforms.
The functions using it only read it, so several threads can use one
progression at the same time.
:::

::: function
void nmod_geometric_progression_init(nmod_geometric_progression_t G, ulong r, slong len, nmod_t mod)
void _nmod_geometric_progression_init_function(nmod_geometric_progression_t G, ulong r, slong len, nmod_t mod, ulong function)

Initialise `G` for the points $1, r^2, \ldots, r^{2(\mathit{len}-1)}$.
`len` must be at least the number of points and at least the length of
the polynomials used. `nmod_geometric_progression_init` prepares
evaluation, interpolation and extrapolation; the second function only the
functionalities selected by the bits of `function`: bit 0 evaluation, bit
1 interpolation, bit 2 extrapolation, and bit 3 stored transforms.

::: details
The precomputation takes $O(\mathit{len})$ operations and space for at
most $11 \mathit{len}$ coefficients: $3 \mathit{len}$ for evaluation, $3
\mathit{len}$ for interpolation and $6 \mathit{len}$ for extrapolation,
of which $\mathit{len}$ are shared with interpolation. The constraints on
$r$ are those of the section introduction; evaluation alone (`function`
= 1 or 9) only needs $r$ invertible.

If bit 3 is set, `_nmod_geometric_progression_init_fft` is called with all
three bounds equal to `len`. It stores the `fft_small` transforms of the
polynomials that the selected functionalities multiply by, for those
functionalities whose largest product would use `fft_small`, so that a
product with a stored transform costs two transforms instead of three.
The products of evaluation and extrapolation are middle products with
the fixed polynomial $v$ as the long operand; for them the stored data
are the images of $v$ through the transposed inverse transform
(`_mul_precomp_init_tr`), and each product is a transposed product
(`_nmod_poly_mul_tr_precomp`) with transforms of length $v_n$ rounded up
to a multiple of 256, where $v_n = 2 \mathit{len} - 1$ (evaluation) or
$\mathit{len} - 1$ (extrapolation), instead of a power of two; see
[](algorithms/polynomials/middle-product#transposed-fft). Interpolation
stores the forward transform of its series (`_mul_precomp_init`). A call
uses the stored data only where it is cheaper than the product without
them. The bit is ignored without `fft_small`. Each stored transform takes
$\mathit{np} \cdot 2^d$ doubles, where $\mathit{np}$ is the number of
primes `fft_small` needs for the modulus and $2^d$ is the smallest power
of two at least $2 \mathit{len} - 1$ (evaluation, interpolation) or at
least $\mathit{len} - 1$ (each of the two transforms for extrapolation).
:::
:::

::: function
void _nmod_geometric_progression_init_fft(nmod_geometric_progression_t G, slong ev_len, slong int_len, slong ext_len)

Stores the `fft_small` transforms described above for the
functionalities of `G`, with a separate bound, at most `G->len`, for each:
evaluation of polynomials of length at most `ev_len`, interpolation from
at most `int_len` points, and extrapolation from at most `ext_len`
values. A bound 0 removes the stored transforms of the functionality and
a negative bound keeps them unchanged. Smaller bounds give shorter
transforms for evaluation ($\mathit{ev\_len} + \mathit{len} - 1$ instead
of $2 \mathit{len} - 1$) and interpolation (about
$2\,\mathit{int\_len}$), and may need fewer primes; the extrapolation
transforms have length $\mathit{len} - 1$ whatever `ext_len`, which
only bounds the number of primes. Longer inputs remain possible
without the stored transforms. Must not be called while other threads use
`G`.
:::

::: function
void nmod_geometric_progression_clear(nmod_geometric_progression_t G)

Frees the memory used by `G`.
:::

::: function
void _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly, slong ilen, const nmod_geometric_progression_t G, slong olen, nmod_t mod)

Sets `vs[i]` to $\mathit{poly}(r^{2i})$ for $0 \le i < \mathit{olen}$,
using `G` prepared for evaluation. Requires $\mathit{ilen}, \mathit{olen}
\le$ `G->len`.

::: details
With $q = r^2$, $\mathit{poly}(q^j) = q^{-j^2/2} \sum_i
(p_i q^{-i^2/2}) q^{(i+j)^2/2}$, so after scaling the coefficients by
$r^{-i^2}$ the values are a middle product with the fixed polynomial
$\sum_k r^{k^2} x^k$, followed by a scaling by $r^{-j^2}$.
:::
:::

::: function
void _nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys, nn_srcptr coeffs, slong ilen, ulong r, slong olen, nmod_t mod)
void nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys, const nmod_poly_t poly, ulong r, slong olen)

The same with a temporary precomputation for evaluation.
:::

::: function
void _nmod_poly_evaluate_geometric_nmod_vec_iter(nn_ptr ys, nn_srcptr coeffs, slong ilen, ulong r, slong olen, nmod_t mod)
void nmod_poly_evaluate_geometric_nmod_vec_iter(nn_ptr ys, const nmod_poly_t poly, ulong r, slong olen)

The same by evaluating at each point with `_nmod_poly_evaluate_nmod`. No
condition on $r$.
:::

::: function
void _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nn_ptr poly, nn_srcptr v, const nmod_geometric_progression_t G, slong len, nmod_t mod)
void nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_t poly, nn_srcptr v, const nmod_geometric_progression_t G, slong len)

Set `poly` to the polynomial of length at most `len` taking the values
`v[i]` at $r^{2i}$, $0 \le i < \mathit{len}$, using `G` prepared for
interpolation. Requires $1 \le \mathit{len} \le$ `G->len`. The
underscore version writes `len` coefficients and does not allow `poly`
to overlap `v`.

::: details
Following [@BS2005, §5], the values are first converted to the Newton
basis of the progression by one short product and then to the monomial
basis by a transposed short product, with scalings by the precomputed
$\prod_{k \le i} (q^k - 1)$ and $q^{i(i-1)/2}$.
:::
:::

::: function
void nmod_poly_interpolate_geometric_nmod_vec_fast(nmod_poly_t poly, ulong r, nn_srcptr ys, slong len)

The same with a temporary precomputation for interpolation.
:::

::: function
void nmod_poly_extrapolate_geometric(nn_ptr oval, slong olen, nn_srcptr ival, slong ilen, slong offset, ulong r, nmod_t mod)
void nmod_poly_extrapolate_geometric_precomp(nn_ptr oval, slong olen, nn_srcptr ival, slong ilen, slong offset, const nmod_geometric_progression_t G)

Given the values `ival` of a polynomial $f$ of degree less than
$\mathit{ilen}$ at `ilen` consecutive points of the progression, compute
its values `oval` at `olen` consecutive points starting `offset`
positions later, without computing $f$. The points must not overlap:
$\mathit{offset} \ge \mathit{ilen}$ (forward) or $\mathit{offset} +
\mathit{olen} \le 0$ (backward).

::: details
Precisely, for any invertible $c$ and $k \ge 0$ (which the functions need
not know), if `ival[i]` $= f(c \cdot r^{2(k+i)})$ for $0 \le i <
\mathit{ilen}$, then `oval[j]` $= f(c \cdot r^{2(\ell + j)})$ for $0 \le j
< \mathit{olen}$ with $\ell = k + \mathit{offset}$. This uses one middle
product. In the forward case none of $r^2, \ldots, r^{2(\mathit{offset} +
\mathit{olen} - 1)}$ may be 1, and `G->len` must be at least
$\mathit{offset} + \mathit{olen}$; in the backward case the same holds for
$\mathit{ilen} - \mathit{offset}$. `G` must be prepared for
extrapolation; the function without `_precomp` builds a temporary `G`.
For $\mathit{ilen} \le 1$ the output is constant.
:::
:::

## Composition {#composition}

Composition $f(g)$ and the Taylor shift $f(x + c)$; see
[](nmod_poly/compose).

::: function
void _nmod_poly_compose(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
void nmod_poly_compose(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

Set `res` to the composition $\mathit{poly1}(\mathit{poly2})$. The
underscore version requires $\mathit{len1}, \mathit{len2} \ge 1$, writes
$(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$ coefficients and does not
allow aliasing.

::: details
A constant `poly2` gives an evaluation, a binomial $a x^k + c$ a
substitution of the monomial and a Taylor shift, Horner's rule is used
for $\mathit{len1} \le 7$, and otherwise the divide-and-conquer algorithm
that splits `poly1` into halves and multiplies by precomputed powers
$\mathit{poly2}^{2^k}$, in $O(M(N) \log N)$ operations for output length
$N$. See [](algorithms/polynomials/composition).
:::
:::

::: function
void _nmod_poly_compose_horner(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
void nmod_poly_compose_horner(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)

The composition by Horner's rule, with the conventions of
`_nmod_poly_compose`.
:::

::: macro
_nmod_poly_compose_divconquer
nmod_poly_compose_divconquer

Obsolete names, defined as `_nmod_poly_compose` and `nmod_poly_compose`.
:::

::: function
void _nmod_poly_taylor_shift(nn_ptr poly, ulong c, slong len, nmod_t mod)
void nmod_poly_taylor_shift(nmod_poly_t g, const nmod_poly_t f, ulong c)

Set $g(x) = f(x + c)$, for reduced $c$; the underscore version works in
place on `(poly, len)`. Uses Horner's rule for $\mathit{len} < 100$ or
$\mathit{len} > n$, and for $c = \pm 1$ and $\mathit{len} < 1000$, and
the convolution method otherwise, which requires $(\mathit{len} - 1)!$ to
be invertible modulo $n$ (for prime $n$, guaranteed by $\mathit{len} \le
n$); otherwise `FLINT_IMPINV` is thrown.
:::

::: function
void _nmod_poly_taylor_shift_horner(nn_ptr poly, ulong c, slong len, nmod_t mod)
void nmod_poly_taylor_shift_horner(nmod_poly_t g, const nmod_poly_t f, ulong c)

The Taylor shift by Horner's rule, with $O(\mathit{len}^2)$ additions
and multiplications. Works for every modulus.
:::

::: function
void _nmod_poly_taylor_shift_convolution(nn_ptr poly, ulong c, slong len, nmod_t mod)
void nmod_poly_taylor_shift_convolution(nmod_poly_t g, const nmod_poly_t f, ulong c)

The Taylor shift as one product of length `len`, using $\sum_k f_k (x +
c)^k = \sum_j \frac{x^j}{j!} \sum_k k!\, f_k \frac{c^{k-j}}{(k-j)!}$
[@ASU1975] [@GG1997], in $O(M(\mathit{len}))$ operations. Requires
$(\mathit{len} - 1)!$ to be invertible modulo $n$.
:::

## Modular composition {#compose-mod}

Modular composition computes $f(g) \bmod h$. By the algorithm of Brent
and Kung [@BrentKung1978] it costs $O(n^{1/2} M(n) + n^{(\omega + 1)/2})$
operations for $\deg h = n$, using a matrix product of size $\sqrt{n}
\times \sqrt{n}$ by $\sqrt{n} \times n$. With $m = \lfloor \sqrt n \rfloor
+ 1$, the $m \times n$ matrix $A$ holding $g^0, \ldots, g^{m-1} \bmod h$
in its rows can be precomputed and reused for many $f$. All functions
require $h \ne 0$ with invertible leading coefficient and allow `res` to be
any input (the non-underscore versions). See [](nmod_poly/compose_mod).

::: summary algorithms/polynomials/composition :::

::: function
void _nmod_poly_compose_mod(nn_ptr res, nn_srcptr f, slong lenf, nn_srcptr g, nn_srcptr h, slong lenh, nmod_t mod)
void nmod_poly_compose_mod(nmod_poly_t res, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h)

Set `res` to $f(g) \bmod h$. The non-underscore version accepts any $f$
and $g$ and throws `FLINT_DIVZERO` for $h = 0$. The underscore version
requires $\mathit{lenh} \ge 2$, $\mathit{lenf} \ge 1$ and $g$ reduced and
zero-padded to length $\mathit{lenh} - 1$, writes $\mathit{lenh} - 1$
coefficients and does not allow aliasing. It uses Horner's rule if
$\mathit{lenh} < 8$ or $\mathit{lenf} \ge \mathit{lenh}$, and the
Brent–Kung algorithm otherwise.
:::

::: function
void _nmod_poly_compose_mod_horner(nn_ptr res, nn_srcptr f, slong lenf, nn_srcptr g, nn_srcptr h, slong lenh, nmod_t mod)
void nmod_poly_compose_mod_horner(nmod_poly_t res, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h)

Modular composition by Horner's rule, with $\mathit{lenf} - 1$
multiplications modulo $h$, and the conventions of
`_nmod_poly_compose_mod`.
:::

::: function
void _nmod_poly_compose_mod_brent_kung(nn_ptr res, nn_srcptr f, slong lenf, nn_srcptr g, nn_srcptr h, slong lenh, nmod_t mod)
void nmod_poly_compose_mod_brent_kung(nmod_poly_t res, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t h)

Modular composition by the Brent–Kung algorithm, with the conventions of
`_nmod_poly_compose_mod`. Requires $\mathit{lenf} < \mathit{lenh}$
($\deg f < \deg h$); the non-underscore version throws `FLINT_ERROR`
otherwise.
:::

::: function
void _nmod_poly_compose_mod_brent_kung_preinv(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, nn_srcptr poly3, slong len3, nn_srcptr poly3inv, slong len3inv, nmod_t mod)
void nmod_poly_compose_mod_brent_kung_preinv(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t poly3, const nmod_poly_t poly3inv)

Set `res` to $\mathit{poly1}(\mathit{poly2}) \bmod \mathit{poly3}$ by the
Brent–Kung algorithm, given the inverse `poly3inv` of the reversal of
`poly3` as for `nmod_poly_mulmod_preinv`. Requires $\mathit{len1} <
\mathit{len3}$; the underscore version also requires `poly2` reduced and
zero-padded to length $\mathit{len3} - 1$.
:::

::: function
void _nmod_poly_precompute_matrix(nmod_mat_t A, nn_srcptr poly1, nn_srcptr poly2, slong len2, nn_srcptr poly2inv, slong len2inv, nmod_t mod)
void nmod_poly_precompute_matrix(nmod_mat_t A, const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t poly2inv)

Set the rows of the $m \times n$ matrix $A$ to $\mathit{poly1}^i \bmod
\mathit{poly2}$ for $0 \le i < m$, where $n = \deg \mathit{poly2}$ and $m =
\lfloor \sqrt{n} \rfloor + 1$, given the inverse `poly2inv` of the
reversal of `poly2`. $A$ must be initialised with these dimensions (the
non-underscore version throws `FLINT_ERROR` otherwise). The underscore
version requires `poly1` reduced and zero-padded to length $n$.
:::

::: function
void _nmod_poly_compose_mod_brent_kung_precomp_preinv(nn_ptr res, nn_srcptr poly1, slong len1, const nmod_mat_t A, nn_srcptr poly3, slong len3, nn_srcptr poly3inv, slong len3inv, nmod_t mod)
void nmod_poly_compose_mod_brent_kung_precomp_preinv(nmod_poly_t res, const nmod_poly_t poly1, const nmod_mat_t A, const nmod_poly_t poly3, const nmod_poly_t poly3inv)

Set `res` to $\mathit{poly1}(g) \bmod \mathit{poly3}$, where $A$ is the
matrix computed by `nmod_poly_precompute_matrix` for $g$ and `poly3`.
Requires $\mathit{len1} < \mathit{len3}$. This is the function to use
for many compositions with the same $g$ and modulus, as in the
computation of Frobenius powers in polynomial factorisation
[@GatSho1992] [@KalSho1998].
:::

::: function
void _nmod_poly_reduce_matrix_mod_poly(nmod_mat_t A, const nmod_mat_t B, const nmod_poly_t f)

Initialises $A$ as an $m \times n$ matrix, $n = \deg f$, $m = \lfloor
\sqrt n \rfloor + 1$, and sets row 0 to 1 and row $i$ to row $i$ of $B$
reduced modulo $f$ for $1 \le i < m$. $A$ must not be initialised before
the call; $B$ needs at least $m$ rows and at least $\deg f + 1$ columns.
:::

::: function
void _nmod_poly_mod_matrix_rows_evaluate(nn_ptr res, const nmod_mat_t A, nn_srcptr h, slong n, nn_srcptr poly3, slong len3, nn_srcptr poly3inv, slong len3inv, nmod_t mod)

Sets `(res, n)` to $\sum_{i} A_i h^i \bmod \mathit{poly3}$, where $A_i$ is
row $i$ of the $d \times n$ matrix $A$, `h` has length $n$ and
$\mathit{len3} = n + 1$. This is the last step of the Brent–Kung
algorithm, with $h = g^m$; it uses Horner's rule for $d \le 10$ and
rectangular splitting otherwise.
:::

## Power series composition and reversion {#series-composition}

::: function
void _nmod_poly_compose_series(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong n, nmod_t mod)
void nmod_poly_compose_series(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)

Set `res` to $\mathit{poly1}(\mathit{poly2}) \bmod x^n$, where `poly2`
must have constant term 0 (the non-underscore version throws
`FLINT_ERROR` otherwise). The underscore version requires $0 <
\mathit{len1}, \mathit{len2} \le n$, writes $n$ coefficients and does not
allow aliasing.

::: details
The algorithm is chosen by `_gr_poly_compose_series`: Horner's rule for
short inputs, a direct substitution for a monomial `poly2`, the
divide-and-conquer algorithm when the composition is not much longer than
$n$, the Brent–Kung baby-step giant-step algorithm [@BrentKung1978] for
$n < 500$, and the algorithm of Kinoshita and Li [@KL2024], with
$O(M(n) \log n)$ operations, from $n = 500$ on.
:::
:::

::: function
void _nmod_poly_revert_series(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, nmod_t mod)
void nmod_poly_revert_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)

Set `Qinv` to the compositional inverse of the power series $Q$ modulo
$x^n$, so that $Q(\mathit{Qinv}(x)) = \mathit{Qinv}(Q(x)) = x \bmod x^n$.
Requires $Q_0 = 0$ and $Q_1$ invertible; the non-underscore version
throws `FLINT_ERROR` otherwise. The underscore version writes $n$
coefficients and does not allow aliasing.

::: details
For $n < 2500$ the fast Lagrange inversion of [@Joh2015b] is used, which
divides by $1, \ldots, n - 1$; if one of them is not invertible (small
characteristic), and for larger $n$, Newton iteration with power series
composition is used instead.
:::
:::

## Square roots {#sqrt}

::: function
void _nmod_poly_sqrt_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_sqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_invsqrt_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_invsqrt_series(nmod_poly_t g, const nmod_poly_t h, slong n)

Set $g$ to $\sqrt{h}$ and $1/\sqrt{h}$ modulo $x^n$. The constant term of
$h$ must be a nonzero square modulo a prime $n \ne 2$ (it is usually 1),
and its square root is chosen by `nmod_sqrt`; otherwise an exception is
thrown. The underscore versions require $\mathit{hlen}, n \ge 1$, write
$n$ coefficients and do not allow aliasing.

::: details
For short $h$ and $n$ larger than the length the coefficients follow from
a linear recurrence (J. C. P. Miller's method for powers of power series,
[@Knu1997, §4.7]), which requires the series length to be less than the modulus. Otherwise the
basecase algorithm is used for small $n$ and Newton iteration for
$1/\sqrt{h}$ for large $n$.
:::
:::

::: function
int _nmod_poly_sqrt(nn_ptr s, nn_srcptr p, slong len, nmod_t mod)
int nmod_poly_sqrt(nmod_poly_t b, const nmod_poly_t a)

If $a$ is the square of a polynomial, set $b$ to a square root and return
1; otherwise return 0 (and set $b = 0$). Requires a prime modulus. The
underscore version requires $\mathit{len} \ge 1$, writes $\mathit{len} / 2
+ 1$ coefficients and does not allow aliasing.

::: details
For $n = 2$, $a$ is a square exactly when its odd coefficients vanish.
Otherwise the low $\mathit{len}/2 + 1$ coefficients of a candidate root
are computed with `_nmod_poly_sqrt_series` after removing an even power of
$x$ and scaling the constant term to 1, and the candidate is accepted if
its square agrees with $a$ in the remaining coefficients.
:::
:::

## Power sums {#power-sums}

The power sums of a polynomial $f$ with roots $\alpha_1, \ldots,
\alpha_d$ (with multiplicity, in an algebraic closure) are $p_k = \sum_i
\alpha_i^k$. They determine $f$ if $d!$ is invertible, by Newton's
identities $k e_k = \sum_{i=1}^{k} (-1)^{i-1} e_{k-i} p_i$, and can be
computed in $O(M(d))$ operations from the series
$\sum_{k \ge 0} p_k x^k = d - x \tilde f'(x)/\tilde f(x)$, where $\tilde
f$ is the reversal of $f$ [@Sch1982b] [@BFSS2006]. Here the power sum series
of $f$ is $\sum_{k \ge 0} p_k x^k$ with $p_0 = d$.

::: function
void _nmod_poly_power_sums(nn_ptr res, nn_srcptr poly, slong len, slong n, nmod_t mod)
void nmod_poly_power_sums(nmod_poly_t res, const nmod_poly_t poly, slong n)

Set `res` to $\sum_{k=0}^{n-1} p_k x^k$. The non-underscore version
accepts any nonzero `poly` with invertible leading coefficient (throwing
`FLINT_ERROR` for zero) and handles factors $x^k$. The underscore version
requires `poly` monic and writes $n$
coefficients. Uses Newton's identities for $10 n < \mathit{len} + 75$ and
the series formula otherwise.
:::

::: function
void _nmod_poly_power_sums_naive(nn_ptr res, nn_srcptr poly, slong len, slong n, nmod_t mod)
void nmod_poly_power_sums_naive(nmod_poly_t res, const nmod_poly_t poly, slong n)
void _nmod_poly_power_sums_schoenhage(nn_ptr res, nn_srcptr poly, slong len, slong n, nmod_t mod)
void nmod_poly_power_sums_schoenhage(nmod_poly_t res, const nmod_poly_t poly, slong n)

The power sums by Newton's identities, in $O(n \cdot \mathit{len})$
operations, and by the series formula with one power series division, in
$O(M(n))$ operations. The underscore versions require `poly` monic.
:::

::: function
void _nmod_poly_power_sums_to_poly(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod)
void nmod_poly_power_sums_to_poly(nmod_poly_t res, const nmod_poly_t Q)

Set `res` to the monic polynomial of degree $d = Q_0$ whose power sum
series is $Q$ (modulo $x^{\mathit{len}}$, the higher power sums being
determined by $d$ of them). Requires $1, \ldots, d$ invertible modulo $n$.
For $Q = 0$ the result is 1. The underscore version writes $d + 1$
coefficients. Uses Newton's identities if $n \le 12$ or $d \le 10$ and
the exponential formula $\tilde f = \exp(\int (d - Q)/x)$ otherwise.
:::

::: function
void _nmod_poly_power_sums_to_poly_naive(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod)
void nmod_poly_power_sums_to_poly_naive(nmod_poly_t res, const nmod_poly_t Q)
void _nmod_poly_power_sums_to_poly_schoenhage(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod)
void nmod_poly_power_sums_to_poly_schoenhage(nmod_poly_t res, const nmod_poly_t Q)

The same by Newton's identities, in $O(d^2)$ operations, and by the
exponential formula, in $O(M(d))$ operations.
:::

## Elementary functions of power series {#transcendental}

For a power series $h$ these functions compute $f(h) \bmod x^n$, where
$$
\exp(h) = \sum_{k \ge 0} \frac{h^k}{k!}, \quad
\log(h) = \int \frac{h'}{h}, \quad
\operatorname{atan}(h) = \int \frac{h'}{1 + h^2}, \quad
\operatorname{atanh}(h) = \int \frac{h'}{1 - h^2},
$$
$$
\operatorname{asin}(h) = \int \frac{h'}{\sqrt{1 - h^2}}, \quad
\operatorname{asinh}(h) = \int \frac{h'}{\sqrt{1 + h^2}},
$$
with integration constants 0, and the trigonometric and hyperbolic
functions are defined by the usual relations to $\exp$. The constant term
of $h$ must be 1 for $\log$ and 0 for all other functions; the
non-underscore functions throw `FLINT_ERROR` otherwise. The coefficients
involve $1/k$ for $k < n$, so the modulus must be a prime $p \ge n$, and
also $p > 2$; this is not checked. All functions cost $O(M(n))$ operations
for large $n$, typically with [Newton iteration](algorithms/polynomials/newton)
[@HZ2004] [@BZ2011], and use quadratic algorithms for small $n$. The
underscore versions require $\mathit{hlen}, n \ge 1$, write $n$
coefficients, and allow $g$ to be the same array as $h$ unless stated
otherwise.

::: function
void _nmod_poly_exp_series(nn_ptr f, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, slong n)

Set $f = \exp(h) \bmod x^n$.

::: details
The algorithm depends on the length of $h$ (not on $n$). If $h$ has
fewer than about 230 to 490 coefficients (with `fft_small`, depending on
the bit size of $p$), the recurrence $k f_k = \sum_{j = 1}^{k} j h_j f_{k -
j}$ from $f' = h' f$ is used, which costs $O(n \cdot \mathit{hlen})$; up
to about 470 to 1490 coefficients the same recurrence with the sums
computed by blocks of multiplications; beyond that Newton iteration.
:::
:::

::: function
void _nmod_poly_exp_expinv_series(nn_ptr f, nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)

Sets $f = \exp(h)$ and $g = \exp(-h) = 1/f$ modulo $x^n$, currently by
`_nmod_poly_exp_series` followed by `_nmod_poly_inv_series`. Requires
$\mathit{hlen} \le n$ and does not allow aliasing.
:::

::: function
void _nmod_poly_log_series(nn_ptr res, nn_srcptr f, slong flen, slong n, nmod_t mod)
void nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, slong n)

Set $\mathit{res} = \log(f) \bmod x^n$, for $f$ with constant term 1,
as the integral of $f'/f$.
:::

::: function
void _nmod_poly_atan_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_atanh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_atanh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_asin_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_asin_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_asinh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_asinh_series(nmod_poly_t g, const nmod_poly_t h, slong n)

Set $g$ to $\operatorname{atan}(h)$, $\operatorname{atanh}(h)$,
$\operatorname{asin}(h)$ and $\operatorname{asinh}(h)$ modulo $x^n$, from
the integrals above.
:::

::: function
void _nmod_poly_sin_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_sin_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_cos_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_cos_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_tan_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_tanh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
void nmod_poly_tanh_series(nmod_poly_t g, const nmod_poly_t h, slong n)

Set $g$ to $\sin(h)$, $\cos(h)$, $\tan(h)$ and $\tanh(h)$ modulo $x^n$.
:::

::: function
void _nmod_poly_sinh_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod)
void nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
void _nmod_poly_cosh_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod)
void nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, slong n)

Set $g$ to $\sinh(h) = (e^h - e^{-h})/2$ and $\cosh(h) = (e^h +
e^{-h})/2$ modulo $x^n$, using `_nmod_poly_exp_expinv_series`. The
underscore versions read exactly $n$ coefficients of $h$ (zero-pad
shorter inputs) and do not allow aliasing.
:::

## Roots {#roots}

Root finding and factorisation in general are in `nmod_poly_factor.h`
(`nmod_poly_roots`, `nmod_poly_factor`). The function here is a fast path
for polynomials that are expected to split into distinct linear factors,
as in sparse interpolation and the decoding of Reed–Solomon codes (see
[](nmod_poly/roots) and [](tutorials/algebra/reed-solomon)). The
algorithms for roots over finite fields are described in
[](algorithms/polynomials/root-finding).

::: function
int nmod_poly_find_distinct_nonzero_roots(ulong * roots, const nmod_poly_t P)

If the nonzero polynomial $P$ of degree $d$ has $d$ distinct nonzero roots
modulo the prime $n$ (not checked), write them to `roots[0]`, ...,
`roots[d - 1]` in no particular order and return 1; otherwise return 0, in
which case `roots` may have been partly written. Constants give 1 and no
roots. For $n = 2$ the result is 0 for $d \ge 2$.

::: details
The polynomial $\gcd(x^{(n-1)/2} - 1, P)$ is the product of the $x - r$
for the nonzero squares $r$ among the roots, and $\gcd(x^{(n-1)/2} + 1,
P)$ that for the non-squares. If their degrees do not add up to $d$, then
$P$ has a repeated root, the root 0 or an irreducible factor of degree
greater than 1, and the function returns 0. Otherwise each factor is split
further with $\gcd((x + \delta)^{(n-1)/2} - 1, \cdot)$ for random $\delta$
(`_nmod_poly_split_rabin`) [@Ber1970] [@Rab1980b]. The powers are computed
with `nmod_poly_powmod_ui_binexp_preinv`. The expected cost is $O(M(d)
\log d \log n)$ operations. The random state is local, so the function
is deterministic.
:::
:::

::: function
void _nmod_poly_split_rabin(nmod_poly_t a, nmod_poly_t b, const nmod_poly_t f, nmod_poly_t t, nmod_poly_t t2, flint_rand_t randstate)

For monic $f$ of degree at least 2 with $\deg f$ distinct nonzero roots
modulo an odd prime, find a factorisation $f = ab$ with $1 \le \deg b \le
\deg a$ by trying random $\delta$ until $\gcd((x + \delta)^{(n-1)/2} - 1,
f)$ is a proper factor. `t` and `t2` are temporaries.
:::

## Irreducible and Conway polynomials {#special}

::: function
void nmod_poly_minimal_irreducible(nmod_poly_t res, ulong n)

Sets `res` to a monic irreducible polynomial of degree $n \ge 1$ over the
prime field $\mathbb{F}_p$, $p$ the modulus of `res`, with the least
possible number of nonzero terms: a binomial $x^n + a$ if one exists,
otherwise a trinomial, a tetranomial or a pentanomial. Throws
`FLINT_ERROR` if none of these exists; it is conjectured that a
tetranomial ($p > 2$) or pentanomial ($p = 2$) always does.

::: details
An irreducible binomial of degree $n$ exists if and only if every prime
factor of $n$ divides $p - 1$, and $4 \mid p - 1$ if $4 \mid n$. Among
the polynomials with the least number of terms the function returns the
first in the following order. A trinomial $x^n + a x^k + b$ comes before
its monic reversal $x^n + a' x^{n-k} + b'$ if $k < n - k$. Writing $f =
x^n + a_1 x^{k_1} + \cdots + a_t x^{k_t}$ with $n > k_1 > \cdots > k_t$,
the tuples $(a_1, \ldots, a_t, k_1, \ldots, k_t)$ are ordered
lexicographically, which favours small coefficients (all 1 if possible)
and then small degrees of the middle terms. Sparse moduli make
reduction cheap, see `_nmod_poly_divrem_newton_n_preinv`. The choice of
defining polynomials for finite fields is discussed in
[](algorithms/finite-fields/defining-polynomials).
:::
:::

::: function
int _nmod_poly_conway(nn_ptr op, ulong prime, slong deg)

If FLINT's compressed copy of Lübeck's tables [@Lüb2004] contains the
Conway polynomial $C_{p,d}$ for $p$ = `prime` and $d$ = `deg`, writes its
$d + 1$ coefficients to `op` and returns 1; otherwise returns 0 (see
[](algorithms/finite-fields/defining-polynomials)). The
tables cover the primes below 109988, with all tabulated degrees for $p <
260$ and degree bounds that decrease with $p$ above (12 for $p < 300$,
then 9, 6, 4).
:::

::: function
ulong _nmod_poly_conway_rand(slong * degree, flint_rand_t state, int type)

Returns a random prime $p$ and sets `degree` to a random $d$ such that
`_nmod_poly_conway(op, p, d)` succeeds, for test code. For `type` 0 the
pair is any tabulated one, for 1 the degree is less than 15, for 2 the
prime is less than $2^{10}$, and for 3 both. Other values of `type` throw
`FLINT_ERROR`.
:::

## Inflation and deflation {#deflation}

::: function
void nmod_poly_inflate(nmod_poly_t result, const nmod_poly_t input, slong inflation)

Sets `result` to $p(x^k)$, where $p$ is `input` and $k$ = `inflation` $\ge
0$. For $k = 0$ the result is the constant $p(1)$.
:::

::: function
void nmod_poly_deflate(nmod_poly_t result, const nmod_poly_t input, slong deflation)

Sets `result` to $p(x^{1/k})$, where $p$ is `input` and $k$ = `deflation`
$\ge 1$ must divide the exponent of every nonzero term (not checked).
Throws `FLINT_DIVZERO` for $k = 0$.
:::

::: function
slong nmod_poly_deflation(const nmod_poly_t input)

Returns the largest $k$ such that `input` is a polynomial in $x^k$, that
is, the GCD of the exponents of its nonzero terms. Returns 0 for the zero
polynomial and 1 for constants.
:::

## Chinese remaindering {#crt}

Given pairwise coprime moduli $m_0, \ldots, m_{l-1}$ and values $v_i$,
these functions find the polynomial $r$ of least degree with $r \equiv v_i
\pmod{m_i}$. The moduli are combined pairwise along a binary tree sorted by
degree; each node stores an idempotent $I$ with $I \equiv 1$ modulo one
child and $I \equiv 0$ modulo the other, and computes $A = B + I(C - B)
\bmod M$. All polynomials must have the same prime modulus.

::: type
nmod_poly_multi_crt_struct
nmod_poly_multi_crt_t

Precomputed data for Chinese remaindering with fixed moduli: a straight
line program of instructions `_nmod_poly_multi_crt_prog_instr`.
:::

::: type
_nmod_poly_multi_crt_prog_instr

One instruction $A = B + I(C - B) \bmod M$ of the program, with the
indices of $A$, $B$ and $C$ (negative indices refer to the inputs), the
idempotent $I$ and the modulus $M$.
:::

::: function
void nmod_poly_multi_crt_init(nmod_poly_multi_crt_t CRT)
void nmod_poly_multi_crt_clear(nmod_poly_multi_crt_t CRT)

Initialise and free `CRT`.
:::

::: function
int nmod_poly_multi_crt_precompute(nmod_poly_multi_crt_t CRT, const nmod_poly_struct * moduli, slong len)
int nmod_poly_multi_crt_precompute_p(nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * moduli, slong len)

Prepare `CRT` for the $\mathit{len} \ge 1$ moduli, given as an array of
polynomials or of pointers. Return 1 on success, which happens if and only
if either $\mathit{len} = 1$ and the modulus is nonzero, or all moduli
have positive degree and are pairwise coprime. On failure (0), later
calls with `CRT` give undefined output.
:::

::: function
void nmod_poly_multi_crt_precomp(nmod_poly_t output, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * values)
void nmod_poly_multi_crt_precomp_p(nmod_poly_t output, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * values)

Set `output` to the polynomial of least degree congruent to `values[i]`
modulo the $i$-th modulus of `CRT`, for all $i$. The values need not be
reduced.
:::

::: function
int nmod_poly_multi_crt(nmod_poly_t output, const nmod_poly_struct * moduli, const nmod_poly_struct * values, slong len)

The same without a separate precomputation. Returns the value of
`nmod_poly_multi_crt_precompute`; on failure the output is undefined.
:::

::: function
slong _nmod_poly_multi_crt_local_size(const nmod_poly_multi_crt_t CRT)
void _nmod_poly_multi_crt_run(nmod_poly_struct * outputs, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * inputs)
void _nmod_poly_multi_crt_run_p(nmod_poly_struct * outputs, const nmod_poly_multi_crt_t CRT, const nmod_poly_struct * const * inputs)

Run the program of `CRT` with caller-supplied temporaries: `outputs` must
be an array of at least `_nmod_poly_multi_crt_local_size(CRT)`
initialised polynomials with the right modulus, and the result is placed
in `outputs[0]`. Since `CRT` is only read, several threads can run it
with their own `outputs`.
:::

## Berlekamp–Massey algorithm {#berlekamp-massey}

The Berlekamp–Massey algorithm [@Ber1968] [@Mas1969] finds the shortest
linear recurrence satisfied by a sequence. FLINT's implementation views it
as a partial extended Euclidean algorithm [@SKHN1975]: after $m$ points
$a_1, \ldots, a_m$ have been added, `nmod_berlekamp_massey_reduce`
computes polynomials $U, V, R$ of the extended Euclidean remainder
sequence of $x^m$ and $S = a_1 x^{m-1} + a_2 x^{m-2} + \cdots + a_m$ with
$$
U x^m + V S = R, \qquad \deg U < \deg V \le m/2, \qquad \deg R < m/2.
$$
Then $\sum_{i \ge 1} a_i x^{-i} = -U/V + R/(V x^m)$, so writing $V =
\sum_i v_i x^i$, the recurrence $\sum_i v_i a_{j+i} = 0$ holds for $1 \le
j < m - \deg R$. If $\deg R < \deg V$, the recurrence holds for the whole
sequence ($1 \le j \le m - \deg V$), and $V$ is, up to a constant
factor, the minimal polynomial (characteristic polynomial of the shortest
recurrence) of the sequence; otherwise more points are needed. Unlike
`fmpz_mod_poly_minpoly`, which returns a polynomial of least degree
annihilating all $m$ points, whose degree can be close to $m$ and which is
then not determined uniquely, $V$ has degree at most $m/2$. The modulus must be prime. See
[](tutorials/algebra/reed-solomon) for an application.

::: type
nmod_berlekamp_massey_struct
nmod_berlekamp_massey_t

The state of the algorithm: the points added so far and the polynomials
$V$ and $R$ of the last two steps of the remainder sequence.
:::

::: function
void nmod_berlekamp_massey_init(nmod_berlekamp_massey_t B, ulong p)
void nmod_berlekamp_massey_clear(nmod_berlekamp_massey_t B)

Initialise `B` with an empty sequence modulo the prime $p$, and free it.
:::

::: function
void nmod_berlekamp_massey_start_over(nmod_berlekamp_massey_t B)
void nmod_berlekamp_massey_set_prime(nmod_berlekamp_massey_t B, ulong p)

Empty the sequence of `B`; the second function also changes the prime
to $p$.
:::

::: function
void nmod_berlekamp_massey_add_points(nmod_berlekamp_massey_t B, const ulong * a, slong count)
void nmod_berlekamp_massey_add_zeros(nmod_berlekamp_massey_t B, slong count)
void nmod_berlekamp_massey_add_point(nmod_berlekamp_massey_t B, ulong a)

Append `count` reduced points, `count` zeros, or one point to the
sequence. The polynomials are updated only by
`nmod_berlekamp_massey_reduce`.
:::

::: function
int nmod_berlekamp_massey_reduce(nmod_berlekamp_massey_t B)

Brings $V$ and $R$ up to date with all points added and returns 1 if $V$
changed, 0 otherwise.

::: details
The new points are absorbed into $R$ and the remainder sequence is
continued with Euclidean steps, or with `nmod_poly_hgcd` when at least 10
more degrees have to be removed, until $\deg R < m/2$; the algorithm and
its correctness are described in
[](algorithms/multivariate/sparse-interpolation#bm-hgcd). For example, after
adding $1, 1, 2, 3$ the function returns 1 with $V = x^2 - x - 1$ (up to
a constant factor) and $\deg R < \deg V$; adding $5, 8$ and reducing
again returns 0.
:::
:::

::: function
slong nmod_berlekamp_massey_point_count(const nmod_berlekamp_massey_t B)
const ulong * nmod_berlekamp_massey_points(const nmod_berlekamp_massey_t B)

Return the number of points added and a pointer to them (which may be
`NULL` if there are none).
:::

::: function
const nmod_poly_struct * nmod_berlekamp_massey_V_poly(const nmod_berlekamp_massey_t B)
const nmod_poly_struct * nmod_berlekamp_massey_R_poly(const nmod_berlekamp_massey_t B)

Return the polynomials $V$ and $R$ as of the last call to
`nmod_berlekamp_massey_reduce`.
:::

::: function
void nmod_berlekamp_massey_print(const nmod_berlekamp_massey_t B)

Prints $V$ in pretty format with variable `#`, a comma and the points.
:::

## Bit packing {#bit-packing}

These functions convert between polynomials and integers by placing the
coefficients in fixed-width bit fields, that is, by evaluating at
$2^{\mathit{bits}}$. They are the building blocks of
[Kronecker substitution](algorithms/polynomials/kronecker).

::: function
void _nmod_poly_bit_pack(nn_ptr res, nn_srcptr poly, slong len, flint_bitcnt_t bits)

Writes $\sum_i \mathit{poly}_i 2^{i \cdot \mathit{bits}}$ to `res` as
$\lceil \mathit{len} \cdot \mathit{bits} / \bits \rceil$ limbs. Requires
$\mathit{len} \ge 1$, $1 \le \mathit{bits} < 3 \bits$ and every
coefficient less than $2^{\mathit{bits}}$.
:::

::: function
void _nmod_poly_bit_unpack(nn_ptr res, slong nlo, slong nhi, nn_srcptr mpn, flint_bitcnt_t bits, nmod_t mod)

Reads the bit fields $\mathit{nlo}, \ldots, \mathit{nhi} - 1$ of width
`bits` from `mpn`, reduces them modulo $n$ and writes them to `res[0]`,
..., `res[nhi - nlo - 1]`. The fields may hold any value. Requires $1 \le
\mathit{bits} < 3 \bits$.
:::

::: function
void nmod_poly_bit_pack(fmpz_t f, const nmod_poly_t poly, flint_bitcnt_t bit_size)

Sets $f = \mathit{poly}(2^{\mathit{bit\_size}})$. Every coefficient must
be less than $2^{\mathit{bit\_size}}$; $1 \le \mathit{bit\_size} < 3 \bits$.
For `bit_size` = 0 or the zero polynomial, $f = 0$.
:::

::: function
void nmod_poly_bit_unpack(nmod_poly_t poly, const fmpz_t f, flint_bitcnt_t bit_size)

Sets `poly` to the polynomial whose coefficients are the bit fields of
width `bit_size` of $f \ge 0$, reduced modulo $n$; the inverse of
`nmod_poly_bit_pack` for reduced coefficients. Throws `FLINT_ERROR` for
$f < 0$. For `bit_size` = 0 the result is 0.
:::

::: function
void _nmod_poly_KS2_pack1(nn_ptr res, nn_srcptr op, slong n, slong s, ulong b, ulong k, slong r)
void _nmod_poly_KS2_pack(nn_ptr res, nn_srcptr op, slong n, slong s, ulong b, ulong k, slong r)

Write $k$ zero bits followed by the $n$ values `op[0]`, `op[s]`, ...,
`op[(n-1) s]` in fields of $b$ bits to `res`, and zero-pad the output to
$r$ limbs if $r \ne 0$. Each value must be less than $2^b$. The `pack1`
version requires $b \le \bits$, the other $b \le 3 \bits$. Used by
`_nmod_poly_mul_KS2` and `_nmod_poly_mul_KS4`.
:::

::: function
void _nmod_poly_KS2_unpack(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k)
void _nmod_poly_KS2_unpack1(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k)
void _nmod_poly_KS2_unpack2(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k)
void _nmod_poly_KS2_unpack3(nn_ptr res, nn_srcptr op, slong n, ulong b, ulong k)

Skip $k$ bits of `op` and extract $n$ fields of $b$ bits, writing each
as one, two or three words to `res`. `_nmod_poly_KS2_unpack1` requires $b
\le \bits$, `unpack2` $\bits < b \le 2 \bits$, `unpack3` $2 \bits < b < 3
\bits$; `_nmod_poly_KS2_unpack` dispatches on $b$.
:::

::: function
void _nmod_poly_KS2_reduce(nn_ptr res, slong s, nn_srcptr op, slong n, ulong w, nmod_t mod)

Reduces $n$ values of $w \in \{1, 2, 3\}$ words each, stored
consecutively in `op`, modulo $n$, writing them to `res[0]`, `res[s]`,
... with stride $s$.
:::

::: function
void _nmod_poly_KS2_recover_reduce(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod)
void _nmod_poly_KS2_recover_reduce1(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod)
void _nmod_poly_KS2_recover_reduce2(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod)
void _nmod_poly_KS2_recover_reduce2b(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod)
void _nmod_poly_KS2_recover_reduce3(nn_ptr res, slong s, nn_srcptr op1, nn_srcptr op2, slong n, ulong b, nmod_t mod)

Recover $n$ coefficients of $2b$ bits from two sequences of $n + 1$ digits
of $b$ bits, unpacked from an integer product evaluated at $2^b$ (`op1`,
low digits first) and from the product of the reversed polynomials (`op2`),
propagating the borrows between them; reduce the coefficients modulo $n$
and write them with stride $s$. Used by `_nmod_poly_mul_KS4`. The variants
require $2b \le \bits$ (`1`), $\bits < 2b < 2 \bits$ (`2`), $b = \bits$
(`2b`) and $2 \bits < 2b \le 3 \bits$ (`3`); the version without suffix
dispatches on $b$.
:::

## Internal functions {#internal}

The following helpers are declared in the header for use by other
modules and by the functions above.

::: function
void _nmod_poly_divrem_basecase_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong A_len, nn_srcptr B, slong B_len, ulong invB, nmod_t mod)

`_nmod_poly_divrem_basecase` given `invB` $= \operatorname{lc}(B)^{-1}$.
Same conventions; $R$ may be the same array as $A$.
:::

::: function
void _nmod_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, ulong invL, nmod_t mod)

Division with remainder for $\mathit{lenA} = \mathit{lenB} + 1$, the
typical step of the Euclidean algorithm, given `invL` $=
\operatorname{lc}(B)^{-1}$. Writes two coefficients to $Q$ and
$\mathit{lenB} - 1$ to $R$. $R$ may be the same array as $A$, but $Q$ must
not overlap $A$. The two quotient coefficients are computed first, and each
remainder coefficient $a_i - q_1 b_{i-1} - q_0 b_i$ is then formed with a
single reduction, using Lemire, Barrett or Möller–Granlund reduction
depending on the size of $n$.
:::

::: type
nmod_poly_res_struct
nmod_poly_res_t

Bookkeeping for the resultant in a half-GCD computation (the resultant so
far, the last leading coefficient and lengths). No function of the
library uses it: the half-GCD functions above use the generic
implementation in [`gr_poly`](gr_poly).
:::

::: type
nmod_poly_matrix_precompute_arg_t

An argument structure for computing the matrix of
`nmod_poly_precompute_matrix` in a thread. It is not used by the library.
:::
