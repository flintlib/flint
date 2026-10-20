---
title: fmpzi.h
header: fmpzi.h
description: Gaussian integers $\mathbb{Z}[i]$: arithmetic, division with remainder, greatest common divisors, primality testing and polynomial multiplication.
category: Algebraic numbers and number fields
order: 1
---

This module implements the ring $\mathbb{Z}[i]$ of Gaussian integers. An
element is an `fmpzi_t`, a pair of `fmpz` integers holding the real and
imaginary parts. Besides the ring operations, the module provides
Euclidean division, which makes $\mathbb{Z}[i]$ a Euclidean domain with
respect to the norm $N(a + bi) = a^2 + b^2$, greatest common divisors
normalised to a canonical associate, removal of the prime $1 + i$,
primality testing, and fast multiplication of polynomials with Gaussian
integer coefficients. The generic ring `gr_ctx_init_fmpzi` (in `gr.h`)
is built on these functions and gives polynomials and matrices over
$\mathbb{Z}[i]$. Factorisation is not implemented.

[](fmpzi/gcd) describes the GCD functions in more detail.

::: summary algorithms/number-theory/gaussian-integers :::

## Types {#types}

::: type
fmpzi_struct
fmpzi_t

A Gaussian integer $a + bi$, with fields `a` and `b` of type `fmpz`. An
`fmpzi_t` is an array of one `fmpzi_struct`, so it is passed by
reference.
:::

::: macro
fmpzi_realref(x)
fmpzi_imagref(x)

Pointers (`fmpz *`) to the real and imaginary parts of `x`.
:::

## Basic manipulation {#basic}

::: function
void fmpzi_init(fmpzi_t x)
void fmpzi_clear(fmpzi_t x)

Initialise `x` to zero, respectively clear it.
:::

::: function
void fmpzi_swap(fmpzi_t x, fmpzi_t y)

Swaps `x` and `y` in constant time.
:::

::: function
void fmpzi_zero(fmpzi_t x)
void fmpzi_one(fmpzi_t x)

Set `x` to 0, respectively 1.
:::

::: function
void fmpzi_set(fmpzi_t res, const fmpzi_t x)

Sets `res` to `x`.
:::

::: function
void fmpzi_set_si_si(fmpzi_t res, slong a, slong b)

Sets `res` to $a + bi$.
:::

::: function
void fmpzi_print(const fmpzi_t x)

Prints `x` to `stdout` in the form `a+b*I` (for example `6+3*I` or
`30-15*I`).
:::

::: function
void fmpzi_randtest(fmpzi_t res, flint_rand_t state, flint_bitcnt_t bits)

Sets `res` to a random element whose parts are drawn by `fmpz_randtest`
with up to `bits` bits.
:::

## Properties {#properties}

::: function
int fmpzi_equal(const fmpzi_t x, const fmpzi_t y)

Returns whether $x = y$.
:::

::: function
int fmpzi_is_zero(const fmpzi_t x)
int fmpzi_is_one(const fmpzi_t x)

Return whether $x = 0$, respectively $x = 1$.
:::

::: function
int fmpzi_is_unit(const fmpzi_t x)

Returns whether $x$ is one of the units $\pm 1, \pm i$.
:::

::: function
slong fmpzi_canonical_unit_i_pow(const fmpzi_t x)

Returns $k \in \{0, 1, 2, 3\}$ such that $i^k x$ is the canonical
associate of $x$: the one with positive real part $a$ and imaginary
part $b$ in $(-a, a]$ (for $x = 0$, $k = 0$).
:::

::: function
void fmpzi_canonicalise_unit(fmpzi_t res, const fmpzi_t x)

Sets `res` to the canonical associate $i^k x$ of $x$.
:::

::: function
slong fmpzi_bits(const fmpzi_t x)

Returns the number of bits of the larger in absolute value of the real
and imaginary parts, 0 for $x = 0$.
:::

::: function
void fmpzi_norm(fmpz_t res, const fmpzi_t x)

Sets `res` to the norm $N(x) = a^2 + b^2$.
:::

## Arithmetic {#arithmetic}

::: function
void fmpzi_conj(fmpzi_t res, const fmpzi_t x)
void fmpzi_neg(fmpzi_t res, const fmpzi_t x)

Set `res` to the complex conjugate $\bar x$, respectively to $-x$.
:::

::: function
void fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
void fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

Set `res` to $x + y$, respectively $x - y$.
:::

::: function
void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
void fmpzi_sqr(fmpzi_t res, const fmpzi_t x)

Set `res` to $xy$, respectively $x^2$.

::: details
Single-word operands use two-word products. When all four parts are
multiword, the product is computed by `flint_mpn_mul_complex` (or
`flint_mpn_sqr_complex`), which uses three real multiplications in the
Karatsuba manner or a joint transform for large operands; otherwise four
`fmpz` multiplications are used (two for squaring, with special cases
for purely real or imaginary input).
:::
:::

::: function
void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp)

Sets `res` to $x^{e}$ with $x^0 = 1$, by left-to-right binary
exponentiation, with a direct formula for purely real or purely
imaginary $x$.
:::

::: function
void fmpzi_mul_i(fmpzi_t z, const fmpzi_t x)
void fmpzi_div_i(fmpzi_t z, const fmpzi_t x)

Set $z$ to $ix$, respectively $x/i = -ix$.
:::

::: function
void fmpzi_mul_i_pow_si(fmpzi_t res, const fmpzi_t z, slong k)

Sets `res` to $i^k z$ for any integer $k$.
:::

## Division {#division}

::: function
void fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y)

Sets $q$ to $x/y$, assuming that the division is exact. Throws if
$y = 0$.

::: details
Purely real and purely imaginary divisors are handled by integer exact
divisions. A quotient of fewer than 45 bits is obtained by rounding a
double-precision approximation; a quotient much shorter than the
divisor is computed from the leading parts of the operands; otherwise
the quotient is $x\bar y / N(y)$ with exact integer divisions.
:::
:::

::: function
void fmpzi_divrem(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)

Sets $q$ and $r$ to the quotient and remainder of $x$ by $y$:
$x = qy + r$, where $q$ is $x/y$ with both coordinates rounded to the
nearest integer, ties rounded upwards. Then $N(r) \le N(y)/2$, and the
pair is uniquely determined. `r` may be `NULL`. Throws if $y = 0$.
:::

::: function
void fmpzi_divrem_approx(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)

Sets $q$ and $r$ with $x = qy + r$ and $N(r) < N(y)$, choosing $q$ from
a floating-point approximation of $x/y$ when the operands have similar
sizes; the choice is implementation-defined. `r` may be `NULL`. Throws
if $y = 0$.
:::

::: function
slong fmpzi_remove_one_plus_i(fmpzi_t res, const fmpzi_t x)

Sets `res` to $x / (1 + i)^k$ for the largest $k$ such that the division
is exact, and returns $k$. For $x = 0$, sets `res` to 0 and returns 0.
:::

## GCD {#gcd}

::: function
void fmpzi_gcd_euclidean(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
void fmpzi_gcd_euclidean_improved(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
void fmpzi_gcd_binary(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
void fmpzi_gcd_shortest(fmpzi_t g, const fmpzi_t x, const fmpzi_t y)
void fmpzi_gcd(fmpzi_t g, const fmpzi_t x, const fmpzi_t y)

Set the output to the greatest common divisor of $x$ and $y$ in
canonical form (see `fmpzi_canonical_unit_i_pow`); the GCD of 0 and 0
is 0. The *euclidean* version is Euclid's algorithm with
`fmpzi_divrem`, the *euclidean_improved* version uses approximate
quotients and finishes in double-precision arithmetic, the *binary*
version is Weilert's $(1+i)$-ary algorithm [@Wei2000] (which falls back
to the Euclidean algorithm above 500 bits), and the *shortest* version
finds the GCD as a shortest vector of the lattice generated by $x$ and
$y$. `fmpzi_gcd` chooses among them.
:::

## Primality testing {#primality}

::: function
int fmpzi_is_prime(const fmpzi_t n)
int fmpzi_is_probabprime(const fmpzi_t n)

Return whether $n$ is a Gaussian prime: a unit times $1 + i$, a unit
times a rational prime $p \equiv 3 \pmod 4$, or an element whose norm
is a rational prime. `fmpzi_is_prime` proves primality of the relevant
rational integer with `fmpz_is_prime`, `fmpzi_is_probabprime` uses the
probable prime test `fmpz_is_probabprime`.
:::

## Polynomial multiplication {#poly}

These functions do not support aliasing: the output `res` must not
overlap the input polynomials.

::: function
int _fmpzi_poly_mulmid_classical_fft_small(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

Computes the coefficients $[\mathtt{nlo}, \mathtt{nhi})$ of the
product of the Gaussian integer polynomials `(poly1, len1)` and
`(poly2, len2)`, with the coefficient arithmetic done pointwise
on fft_small transforms: two transforms per coefficient (real and
imaginary parts) and four pointwise operations per complex
product, the shorter operand's transforms kept and the longer
rolling as in `_fmpz_poly_mulmid_classical_fft_small`.
Squaring shares transforms, pairs symmetric terms and computes
diagonal squares with three multiplications. The lengths may be
given in either order. Returns 1 on success and 0 when the
transformed representation is unavailable or unprofitable in which
case nothing is written and the caller should fall back to a
different routine.

:::

::: function
void _fmpzi_poly_mulmid(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

Sets `(res, nhi - nlo)` to the coefficients at indices
$[\mathtt{nlo}, \mathtt{nhi})$ of the product of the Gaussian integer
polynomials `(poly1, len1)` and `(poly2, len2)`, choosing an
algorithm. Assumes that `len1` and `len2` are positive and that
$0 \le nlo < nhi \le len1 + len2 - 1$; the lengths may be given in
either order, and `res` must not alias the inputs. This is the
multiplication used by the `fmpzi` generic ring for polynomial
multiplication, low products and middle products.

::: details
The choice depends on $b$, the average over the two operands of the
larger part size in bits, and on $q$, the number of coefficient
products in the window per operand coefficient ($n/2$ for an
$n \times n$ product, $n/3$ for the Newton middle product
$(2n - 1) \times n \to [n - 1, 2n - 1)$): the classical algorithm
for windows of at most 12 products; four real products when a part
vanishes or is much smaller than the other (they skip the vanishing
products); `_fmpzi_poly_mulmid_classical_fft_small` from 8192
bits up to a bound on $q$ increasing with $b$; real products, three
(Karatsuba) from about 256 bits, below a bound on $q$ decreasing with
$b$; and otherwise `_fmpzi_poly_mulmid_KS`, except for
`_fmpzi_poly_mulmid_SS` on squares at 384-3072 bits, short full
products at 768-1536 bits and middle products whose cyclic transform
saves a quarter of the length at 768-3072 bits (short ones) and from
24576 bits, and for huge products ($m_1 + m_2 \ge 32768$ and
$(len_1 + len_2)(m_1 + m_2) > 2.5 \cdot 10^8$ for the larger part sizes
$m_1, m_2$ of the operands), where Schoenhage-Strassen runs within
about 10% of the time in about a fifth of the memory. The boundaries were
fitted to timings on x86-64 with uniformly random coefficients.

:::
:::

::: function
void _fmpzi_poly_mulmid_KS(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

Sets `(res, nhi - nlo)` to the coefficients at indices
$[\mathtt{nlo}, \mathtt{nhi})$ of the product of the Gaussian integer
polynomials `(poly1, len1)` and `(poly2, len2)`, using Kronecker
substitution. Assumes that `len1` and `len2` are positive and
that $0 \le nlo < nhi \le len1 + len2 - 1$; the lengths may be given
in either order, and `res` must not alias the inputs.

::: details
The real and imaginary parts of each operand are packed into
integers at a common field size, and a single
`flint_mpn_mul_complex` (or `flint_mpn_sqr_complex`)
call forms both parts of the packed product. Above the transform
cutoff of those functions the four packed parts are thus
transformed once and combined pointwise (four forward and two
inverse transforms, two and two for squaring), and below it their
Karatsuba or schoolbook schemes apply. Windows are computed as
two's-complement windows of the packed parts, which unpack directly
into signed coefficients: with `flint_mpn_mullow_complex` (or
`flint_mpn_sqrlow_complex`) when the window starts at the
bottom and ends below the top of the product, and with
`flint_mpn_mulmid_complex` (or `flint_mpn_sqrmid_complex`)
when it starts above the bottom, three guard limbs below the first
field absorbing the error of the windowed products. Above the
transform cutoff the latter reconstructs only the CRT coefficients
reaching the window and, when the window allows it, uses a cyclic
transform, which is where a middle product saves transform length. When one operand is real,
the products are formed separately, without sharing the transform
of the other operand.

:::
:::

::: function
void _fmpzi_poly_mulmid_SS(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)
void _fmpzi_poly_mulmid_SS_pair(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)
void _fmpzi_poly_mulmid_SS_single(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

Sets `(res, nhi - nlo)` to the coefficients at indices
$[\mathtt{nlo}, \mathtt{nhi})$ of the product of the Gaussian integer
polynomials `(poly1, len1)` and `(poly2, len2)`, using
Schoenhage-Strassen convolution. Assumes that `len1` and `len2`
are positive and that $0 \le nlo < nhi \le len1 + len2 - 1$; the
lengths may be given in either order, and `res` must not alias the
inputs.

::: details
In the coefficient ring $\mathbb{Z}/(2^N+1)$ of the transform,
$s = 2^{N/2}$ is a square root of $-1$, and
$\mathbb{Z}[i]/(2^N+1) \cong \mathbb{Z}/(2^N+1) \times \mathbb{Z}/(2^N+1)$
via $a + bi \mapsto (a + sb, a - sb)$. The *pair* method
(`_SS_pair`) therefore computes the Gaussian product as two
ordinary convolutions, of
$u = \operatorname{re} + s \operatorname{im}$ and
$v = \operatorname{re} - s \operatorname{im}$, with a ring size
covering the components of the result, and recovered as
$\operatorname{re} = (U + V)/2$,
$\operatorname{im} = (V - U) s / 2$; all conversions are shifts and
additions. This uses six transforms and two pointwise products per
point (four transforms for squaring). The *single embedding*
(`_SS_single`) uses the homomorphism $a + bi \mapsto a + sb$ alone:
one convolution in a ring of twice the size gives
$T = \operatorname{re} + s \operatorname{im}$, from whose two's
complement representation both parts are read off, the ring being
large enough that $|T| < 2^{N-1}$. The transform work and the memory
are about the same as for the pair method, but the ring is rounded
up to a multiple of the transform length once rather than twice,
which favors it for small coefficients relative to the length; it
also halves the number of pointwise products, which favors it for
small rings, while for medium rings the pointwise products of the
larger ring cost more than twice as much. The main function chooses
between them by the rounded ring sizes, as measured on x86-64. The
convolution length is chosen as by `fft_convolution_choose_mid`, so
middle products
use the cyclic (wraparound) transform when it is shorter. When an
operand is real or purely imaginary, its nonzero part is transformed
once and multiplied by both parts of the other operand (two
convolutions sharing that transform, five transforms in all), or by
the single nonzero part of the other operand (one convolution).

:::
:::
