---
title: fmpq.h
header: fmpq.h
description: Rational numbers as fractions of multiprecision integers in canonical form, with arithmetic, rational reconstruction, continued fractions and enumeration of the rationals.
category: Rational numbers
order: 1
---

This module implements rational numbers. An `fmpq_t` is a numerator and a
denominator stored as two `fmpz` integers, so that a fraction whose
numerator and denominator fit in 62 bits takes no heap memory, and its
arithmetic uses word-size GCDs and double-word products. For larger
components the cost is that of `mpq_t` arithmetic plus a small overhead.
Besides arithmetic, the module provides rational reconstruction from a
residue, continued fractions with a subquadratic algorithm, enumeration
of the rationals, Farey neighbours, harmonic numbers and Dedekind sums.
Vectors of rationals are in `fmpq_vec`; polynomials and matrices over
$\Q$ in `fmpq_poly` and `fmpq_mat`.

Start with the guide [](fmpq/canonical-form), which explains the
canonical form that almost every function relies on and the low-level
interface on pairs of `fmpz`. The conventions common to all FLINT modules
(naming, underscore functions, aliasing, canonical forms) are described
in [](guide/conventions/types-and-naming).

## Types and conventions {#types}

::: type
fmpq

A rational number: a struct with two `fmpz` fields `num` (numerator) and
`den` (denominator). The type is defined in `flint.h`.
:::

::: type
fmpq_t

An array of length one of `fmpq`, so that an `fmpq_t` is passed by
reference, as with `mpq_t`.
:::

::: macro
fmpq_numref(x)
fmpq_denref(x)

Pointers (of type `fmpz *`) to the numerator and the denominator of the
`fmpq` that `x` points to. Defined in `flint.h`.
:::

::: function
fmpz * fmpq_numerator_ptr(fmpq_t q)
fmpz * fmpq_denominator_ptr(fmpq_t q)

Function versions of `fmpq_numref` and `fmpq_denref`, for use from
languages that cannot expand C macros.
:::

::: function
void fmpq_numerator(fmpz_t n, const fmpq_t q)
void fmpq_denominator(fmpz_t n, const fmpq_t q)

Set $n$ to a copy of the numerator or the denominator of $q$.
:::

A fraction $p/q$ is in *canonical form* if $\gcd(p, q) = 1$ and $q > 0$;
zero is $0/1$. Unless stated otherwise, every function assumes that its
`fmpq_t` inputs are canonical and returns canonical outputs. This is not
checked: a non-canonical input gives undefined results. Code that writes
the numerator and denominator directly must call `fmpq_canonicalise`
before passing the number to a library function.

Functions with a leading underscore take the numerator and denominator as
separate `fmpz_t` arguments (`rnum, rden` for the output). They follow
the same canonical form convention unless documented otherwise. Inputs may
be aliased with outputs in any combination, except that a numerator must
not be aliased with a denominator; functions without an underscore allow
any aliasing.

## Memory management {#memory}

::: function
void fmpq_init(fmpq_t x)

Initialises $x$ and sets it to zero. This does not allocate memory.
:::

::: function
void fmpq_clear(fmpq_t x)

Frees the memory used by $x$. It must be initialised again before reuse.
:::

## Canonicalisation {#canonical}

::: function
void fmpq_canonicalise(fmpq_t res)
void _fmpq_canonicalise(fmpz_t num, fmpz_t den)

Put a fraction with nonzero denominator in canonical form: divide
numerator and denominator by their GCD and make the denominator positive.
A zero numerator gives $0/1$. The result is undefined if the denominator
is zero. In the underscore version, `num` and `den` must not be aliased.
:::

::: function
int fmpq_is_canonical(const fmpq_t x)
int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den)

Return nonzero if the fraction is in canonical form, that is, the
denominator is positive and coprime to the numerator, and zero is
represented as $0/1$.
:::

## Assignment and conversion {#conversion}

::: function
void fmpq_set(fmpq_t dest, const fmpq_t src)

Sets `dest` to a copy of `src`. No canonicalisation is done, so this also
copies non-canonical fractions.
:::

::: function
void fmpq_swap(fmpq_t op1, fmpq_t op2)

Swaps `op1` and `op2` efficiently (only the structs are exchanged).
:::

::: function
void fmpq_zero(fmpq_t res)
void fmpq_one(fmpq_t res)

Set `res` to $0$ or $1$.
:::

::: function
void fmpq_set_fmpz(fmpq_t q, const fmpz_t n)

Sets $q$ to the integer $n$.
:::

::: function
void fmpq_set_si(fmpq_t res, slong p, ulong q)
void fmpq_set_ui(fmpq_t res, ulong p, ulong q)
void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q)
void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q)

Set the result to the canonical form of $p/q$. Requires $q \ne 0$. In the
underscore versions `rnum` and `rden` must not be aliased.
:::

::: function
void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)

Sets `res` to the canonical form of $p/q$, for any integers $p$ and $q \ne
0$ (not checked).
:::

::: function
void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c)

Convert from and to GMP rationals; `fmpq_get_mpz_frac` sets `a` and `b`
to the numerator and denominator of `c`. The value is copied as it is,
without canonicalisation. These functions are declared only if `gmp.h`
is included before any FLINT header.
:::

::: function
int fmpq_set_str(fmpq_t dest, const char * s, int base)

Sets `dest` to the integer or fraction `"p"` or `"p/q"` in the string
`s`, written in base `base` as accepted by `mpq_set_str`. Returns 0 on
success and $-1$ if the string is not valid, in which case `dest` is not
changed.

::: details
Like `mpq_set_str`, this function does not canonicalise: `"2/4"` gives the
non-canonical $2/4$, and `"1/0"` is accepted. Call `fmpq_canonicalise`
afterwards if the input may not be in lowest terms.
:::
:::

::: function
double fmpq_get_d(const fmpq_t f)

Returns $f$ as a `double`, rounded towards zero (by `mpq_get_d`). The
result is system dependent if $|f|$ is outside the range of `double`.
:::

::: function
int fmpq_get_mpfr(mpfr_t dest, const fmpq_t src, mpfr_rnd_t rnd)

Sets `dest` to `src` rounded in the direction `rnd` at the precision of
`dest`, and returns the ternary value of MPFR (the sign of the rounding
error). Declared only if `mpfr.h` is included before any FLINT header.
:::

::: function
char * fmpq_get_str(char * str, int b, const fmpq_t x)
char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den)

Write $x$ in base $b \in [2, 36]$ as `"p/q"`, or as `"p"` if the
denominator is 1. If `str` is `NULL`, a buffer is allocated with
`flint_malloc` (to be freed with `flint_free`) and returned; otherwise
`str` must have room for the digits of both components plus three
characters, and is returned.
:::

::: function
void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f)
void flint_mpq_clear_readonly(mpq_t z)

Make the uninitialised `mpq_t` $z$ a read-only view of $f$, without
copying large components, and release it. $f$ must not change while $z$
is in use, and $z$ must not be modified. Declared only if `gmp.h` is
included before any FLINT header.

::: details
The intended use is to pass an `fmpq_t` to a function that takes an
`mpq_t` argument:

```
fmpq_t f;
...
{
    mpq_t z;

    flint_mpq_init_set_readonly(z, f);
    foo(..., z);
    flint_mpq_clear_readonly(z);
}
```
:::
:::

::: function
void fmpq_init_set_readonly(fmpq_t f, const mpq_t z)
void fmpq_clear_readonly(fmpq_t f)

Make the uninitialised `fmpq_t` $f$ a read-only view of the `mpq_t` $z$,
and release it. $z$ must not change while $f$ is in use, and $f$ must not
be modified. `fmpq_init_set_readonly` is declared only if `gmp.h` is
included before any FLINT header.

::: details
The intended use is to call FLINT functions on GMP rationals:

```
mpq_t z;
...
{
    fmpq_t f;

    fmpq_init_set_readonly(f, z);
    foo(..., f);
    fmpq_clear_readonly(f);
}
```
:::
:::

## Input and output {#io}

::: function
int fmpq_print(const fmpq_t x)
int fmpq_fprint(FILE * file, const fmpq_t x)
int _fmpq_print(const fmpz_t num, const fmpz_t den)
int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den)

Print the fraction to `stdout` or to `file`, in decimal, as `p/q`, or as
`p` if the denominator is 1. Return a positive value on success and a
non-positive value on failure.
:::

## Random generation {#random}

::: function
void fmpq_randtest(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
void _fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)

Set the result to a random canonical fraction whose numerator and
denominator have at most `bits` bits. The distribution favours special
values ($\pm 1$, $\pm 2$, denominators 1 and 2, values with long runs of
bits) that trigger corner cases. For `bits` equal to 1 the result is an
integer. `num` and `den` must not be aliased.
:::

::: function
void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)

As `fmpq_randtest`, but the result is nonzero. Throws if `bits` is zero.
:::

::: function
void fmpq_randbits(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
void _fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)

Choose a numerator and a nonzero denominator of exactly `bits` bits each
and canonicalise, so the components of the result can have fewer bits.
`num` and `den` must not be aliased.
:::

## Comparison and properties {#comparison}

::: function
int fmpq_is_zero(const fmpq_t x)
int fmpq_is_one(const fmpq_t x)
int fmpq_is_pm1(const fmpq_t x)

Return nonzero if $x$ is $0$, $1$ or $\pm 1$ respectively.
:::

::: function
int fmpq_equal(const fmpq_t x, const fmpq_t y)
int fmpq_equal_fmpz(fmpq_t x, fmpz_t y)
int fmpq_equal_si(fmpq_t x, slong y)
int fmpq_equal_ui(fmpq_t x, ulong y)

Return nonzero if $x = y$. For canonical inputs this is a comparison of
the components.
:::

::: function
int fmpq_sgn(const fmpq_t x)

Returns $-1$, $0$ or $1$ according to the sign of $x$.
:::

::: function
int fmpq_cmp(const fmpq_t x, const fmpq_t y)
int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
int fmpq_cmp_si(const fmpq_t x, slong y)
int fmpq_cmp_ui(const fmpq_t x, ulong y)

Return a negative value, zero or a positive value according as $x < y$,
$x = y$ or $x > y$.
:::

::: function
int _fmpq_cmp(const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)
int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r)
int _fmpq_cmp_si(const fmpz_t p, const fmpz_t q, slong c)
int _fmpq_cmp_ui(const fmpz_t p, const fmpz_t q, ulong c)

Compare $p/q$ with $r/s$, $r$ or $c$, returning a negative value, zero or a
positive value. The denominators must be positive; the fractions need not
be in lowest terms.
:::

::: function
void fmpq_height(fmpz_t height, const fmpq_t x)
flint_bitcnt_t fmpq_height_bits(const fmpq_t x)

The *height* of $x = p/q$ is $\max(|p|, q)$. The first function sets
`height` to it, the second returns its number of bits.
:::

## Arithmetic {#arithmetic}

The functions in this section return canonical results for canonical
inputs. They use the structure of canonical fractions to replace a GCD of
the full result by GCDs of smaller numbers; see
[](fmpq/canonical-form).

::: function
void fmpq_neg(fmpq_t dest, const fmpq_t src)
void fmpq_abs(fmpq_t dest, const fmpq_t src)

Set `dest` to $-\mathit{src}$ or $|\mathit{src}|$.
:::

::: function
void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

Set `res` to $\mathit{op1} + \mathit{op2}$, $\mathit{op1} - \mathit{op2}$,
$\mathit{op1} \cdot \mathit{op2}$ or $\mathit{op1} / \mathit{op2}$.
`fmpq_div` throws if `op2` is zero.
:::

::: function
void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
void _fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

Set `(rnum, rden)` to the canonical sum, difference, product or quotient
of the canonical fractions `(op1num, op1den)` and `(op2num, op2den)`.
`_fmpq_div` requires a nonzero divisor (not checked).
:::

::: function
void _fmpq_add_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)
void _fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)

Set `(rnum, rden)` to the canonical sum or product of $p_1/q_1$ and
$p_2/q_2$, which must be canonical with all four values of absolute value
at most `COEFF_MAX` (that is, small `fmpz` values). The intermediate
results are computed in two words. These are the fast paths of the
general functions.
:::

::: function
void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c)
void fmpq_sub_si(fmpq_t res, const fmpq_t op1, slong c)
void fmpq_add_ui(fmpq_t res, const fmpq_t op1, ulong c)
void fmpq_sub_ui(fmpq_t res, const fmpq_t op1, ulong c)
void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)

Set `res` to $\mathit{op1} + c$ or $\mathit{op1} - c$ for an integer $c$.
:::

::: function
void _fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
void _fmpq_sub_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
void _fmpq_add_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
void _fmpq_sub_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
void _fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)
void _fmpq_sub_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)

Set `(rnum, rden)` to $p/q + r$ or $p/q - r$. No GCD is needed: $(p \pm
qr)/q$ is canonical when $p/q$ is.
:::

::: function
void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c)
void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c)
void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
void _fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
void _fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)

Set the result to the product of a fraction and an integer.
:::

::: function
void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)

Sets `res` to $\mathit{op}/x$. Throws if $x$ is zero.
:::

::: function
void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp)
void fmpq_div_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp)

Set `res` to $x \cdot 2^{\mathit{exp}}$ or $x / 2^{\mathit{exp}}$,
cancelling powers of two from the denominator or the numerator first.
:::

::: function
void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

Set the result to $\mathit{res} + \mathit{op1} \cdot \mathit{op2}$ or
$\mathit{res} - \mathit{op1} \cdot \mathit{op2}$. The product is formed in
a temporary, so any aliasing is allowed (except, in the underscore
versions, a numerator with a denominator).
:::

::: function
void fmpq_inv(fmpq_t dest, const fmpq_t src)

Sets `dest` to $1/\mathit{src}$ by swapping the components and fixing the
sign. `src` must be nonzero; this is not checked, and zero gives the
invalid value $1/0$.
:::

::: function
void fmpq_pow_si(fmpq_t res, const fmpq_t op, slong e)
void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, const fmpz_t opnum, const fmpz_t opden, slong e)

Set the result to $\mathit{op}^e$, with $0^0 = 1$. `fmpq_pow_si` throws if
$\mathit{op} = 0$ and $e < 0$; the underscore version does not check.
:::

::: function
int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)

Sets $a$ to $b^e$ and returns 1, or returns 0 (leaving $a$ unchanged) if
$e$ does not fit in a `slong` and $b \notin \{0, 1, -1\}$, in which case
the result could not be stored. Throws if $b = 0$ and $e < 0$.
:::

## Greatest common divisor {#gcd}

::: function
void fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
void _fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)

Set the result to the GCD of $p/q$ and $r/s$, defined as the
canonicalisation of $\gcd(ps, qr)/(qs)$. For canonical inputs this is
$\gcd(p, r)/\operatorname{lcm}(q, s)$, the largest rational $g$ such that
$p/q$ and $r/s$ are both integer multiples of $g$; it is nonnegative, and
zero only if both inputs are zero.

::: details
This definition is invariant under scaling numerator and denominator,
agrees with the GCD on the integers, and gives the same results as Sage
and PARI/GP. It is not the GCD of `fmpq_poly`, where the GCD of nonzero
constants is 1.

`_fmpq_gcd` does not require canonical inputs, but the denominators must
be nonzero.
:::
:::

::: function
void fmpq_gcd_cofactors(fmpq_t g, fmpz_t abar, fmpz_t bbar, const fmpq_t a, const fmpq_t b)
void _fmpq_gcd_cofactors(fmpz_t gnum, fmpz_t gden, fmpz_t abar, fmpz_t bbar, const fmpz_t anum, const fmpz_t aden, const fmpz_t bnum, const fmpz_t bden)

Set $g = \gcd(a, b)$ as in `fmpq_gcd` and the coprime integer cofactors
$\bar a = a/g$ and $\bar b = b/g$. If $a = b = 0$, then $g = 0$ and $\bar
a = \bar b = 0$. Unlike `_fmpq_gcd`, the inputs must be canonical.
:::

## Modular reduction and rational reconstruction {#reconstruction}

A rational number $n/d$ with $\gcd(d, m) = 1$ has a residue $a \equiv n
d^{-1} \pmod m$. Rational reconstruction inverts this map for fractions
with small numerator and denominator; it is how multimodular and $p$-adic
algorithms return rational results. See [](fmpq/reconstruct) for the
details and examples.

::: summary algorithms/number-theory/rational-reconstruction :::

::: function
int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)
int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod)

If the denominator $d$ of $x = n/d$ is invertible modulo $m$ = `mod`, set
`res` to the residue $a \in [0, m)$ with $n \equiv ad \pmod m$ and return
1; otherwise return 0 and leave an unspecified value in `res`. Requires
$m > 0$.
:::

::: function
int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
int _fmpq_reconstruct_fmpz_2_naive(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)

Given $m > 2$, $0 \le a < m$ and $N, D > 0$ with $2ND < m$, find the
fraction $n/d$ with $|n| \le N$, $0 < d \le D$, $\gcd(n, d) = 1$ and $n
\equiv ad \pmod m$, which is unique if it exists. Return 1 and set the
result if it exists, otherwise return 0.
:::

::: function
int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
int _fmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m)

Rational reconstruction with the balanced bounds $N = D =
\lfloor\sqrt{\lfloor (m-1)/2 \rfloor}\rfloor$, which satisfy $2ND < m$.
:::

## Enumeration of the rationals {#enumeration}

::: function
void fmpq_next_minimal(fmpq_t res, const fmpq_t x)
void _fmpq_next_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

Given a nonnegative canonical $x$, set the result to the next element of
the enumeration of the nonnegative rationals by increasing height
$\max(p, q)$: for $q = 1, 2, 3, \ldots$ the fractions $p/q$ and $q/p$ for
$1 \le p < q$ coprime to $q$, in increasing order of $p$. Starting from
zero, every nonnegative rational is generated exactly once:
$0, 1, 1/2, 2, 1/3, 3, 2/3, 3/2, 1/4, 4, 3/4, 4/3, 1/5, 5, 2/5, \ldots$

::: details
Finding the next $p$ coprime to $q$ needs GCDs, so this enumeration is
somewhat slower than `fmpq_next_calkin_wilf`.
:::
:::

::: function
void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x)
void _fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

As `fmpq_next_minimal`, with each positive value followed by its
negative: $0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots$ Starting from
zero, this generates every rational exactly once, by increasing height.
:::

::: function
void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)
void _fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

Given a nonnegative canonical $x$, set the result to the next fraction in
the breadth-first traversal of the Calkin–Wilf tree [@CalkinWilf2000]:
$0, 1, 1/2, 2, 1/3, 3/2, 2/3, 3, 1/4, 4/3, 3/5, 5/2, 2/5, \ldots$ Starting
from zero, every nonnegative rational is generated exactly once.

::: details
The successor of $x > 0$ is $1/(2\lfloor x \rfloor + 1 - x)$, which needs
no GCD, so this enumeration is cheaper than `fmpq_next_minimal`. It is not
ordered by height: some fractions of small height appear late.
:::
:::

::: function
void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x)
void _fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

As `fmpq_next_calkin_wilf`, with each positive value followed by its
negative: $0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3/2, \ldots$
:::

::: function
void fmpq_farey_neighbors(fmpq_t l, fmpq_t r, const fmpq_t x, const fmpz_t Q)

Set $l$ and $r$ to the fractions immediately below and above $x$ among
the fractions with denominator at most $Q$ (the Farey sequence of order
$Q$, extended to all of $\Q$). Throws if the denominator of $x$ exceeds
$Q$. Any aliasing is allowed.

::: details
Writing $x = p/q$, the neighbours are the fractions $p_l/q_l < x <
p_r/q_r$ with $p q_l - p_l q = 1 = p_r q - p q_r$ and $Q - q < q_l, q_r
\le Q$ [@HarWri2008, Chapter III]. The function computes $q_l$ from
$p^{-1} \bmod q$ and the rest from these identities.
:::
:::

::: function
void fmpq_simplest_between(fmpq_t x, const fmpq_t l, const fmpq_t r)
void _fmpq_simplest_between(fmpz_t x_num, fmpz_t x_den, const fmpz_t l_num, const fmpz_t l_den, const fmpz_t r_num, const fmpz_t r_den)

Set $x$ to the simplest fraction in the closed interval between $l$ and
$r$, where $a_1/b_1$ is simpler than $a_2/b_2$ if $b_1 < b_2$, or $b_1 =
b_2$ and $a_1 < a_2$. The endpoints need not be canonical, but their
denominators must be positive; the result is canonical.

::: details
`fmpq_simplest_between` accepts the endpoints in either order; the
underscore version requires $l \le r$. The simplest fraction is read off
the common prefix of the continued fraction expansions of the endpoints,
computed by `_fmpq_ball_get_cfrac`; see
[](algorithms/number-theory/continued-fractions).
:::
:::

## Continued fractions {#cfrac}

::: summary algorithms/number-theory/continued-fractions :::

::: function
slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)
slong fmpq_get_cfrac_naive(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)

Write the first $k \le n$ terms of the simple continued fraction of $x$
to $c$ and the remainder to `rem`, and return $k$. The expansion
generated is the shorter of the two. `fmpq_get_cfrac` uses a half-GCD
algorithm and runs in quasi-linear time, `fmpq_get_cfrac_naive` uses
repeated division. See [](fmpq/cfrac).
:::

::: function
void fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n)

Sets $x$ to the value of the continued fraction $[c_0; c_1, \ldots,
c_{n-1}]$, for $n > 0$ and positive $c_1, \ldots, c_{n-1}$. Uses a
balanced product of $2 \times 2$ matrices. See [](fmpq/cfrac).
:::

::: function
slong fmpq_cfrac_bound(const fmpq_t x)

Returns an upper bound for the number of terms of the continued fraction
of $x$: 1 for an integer, otherwise $\lfloor 1.4404\, b \rfloor + 2$ where
$b$ is the bit length of the denominator. The bound is not sharp.

::: details
A fraction whose shorter expansion has $k$ terms has denominator at least
the Fibonacci number $F_{k+1} \ge \varphi^{k-1}$, with equality for
$[c_0; 1, \ldots, 1, 2]$. The constant is $1/\log_2 \varphi$.
:::
:::

## Harmonic numbers {#harmonic}

::: function
void fmpq_harmonic_ui(fmpq_t x, ulong n)
void _fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n)

Set the result to the harmonic number $H_n = 1 + 1/2 + \cdots + 1/n$,
with $H_0 = 0$. Throws if $n \ge 2^{\bits - 1}$.

::: details
Values whose numerator and denominator fit in a word ($n \le 46$ on a
64-bit machine) are read from a table. For larger $n$ the function uses
$H_n = H_{\lfloor n/2 \rfloor}/2 + \sum_{k \le n,\, k \text{ odd}} 1/k$
recursively, which gives every odd $k$ a weight $(2^d - 1)/2^{d-1}$ that
depends only on the dyadic range $n/2^d < k \le n/2^{d-1}$, and sums the
weighted odd terms by [](algorithms/integers/binary-splitting),
accumulating word-size partial sums in the base case. A final GCD puts
the result in lowest terms. The cost is $O(M(n \log n) \log n)$ bit
operations.
:::
:::

## Dedekind sums {#dedekind}

The Dedekind sum $s(h, k)$ is defined for integers $h$ and $k$ by

$$
s(h,k) = \sum_{i=1}^{k-1} \left(\left(\frac{i}{k}\right)\right)
\left(\left(\frac{hi}{k}\right)\right),
\qquad
((x)) = \begin{cases} x - \lfloor x \rfloor - 1/2 & x \notin \Z, \\
0 & x \in \Z. \end{cases}
$$

It appears in the transformation formula of the Dedekind eta function and
in the exponential sums of the partition function [@Apostol1997],
[@RadGro1972]. It depends only on $h \bmod k$, satisfies $s(-h, k) =
-s(h, k)$ and $s(gh, gk) = s(h, k)$, and the empty sum gives $s(h, k) =
0$ for $k \le 1$.

::: function
void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)

Set $s$ to $s(h, k)$, for any integers $h$ and $k$. The naive version
evaluates the defining sum with $k$ terms. `fmpq_dedekind_sum` uses the
continued fraction of $k/h$ and costs about as much as an extended GCD.

::: details
By the reciprocity law $s(h,k) + s(k,h) = (h/k + k/h + 1/(hk))/12 - 1/4$
for coprime $h, k > 0$, the sum telescopes along the Euclidean algorithm
[@Hickerson1977]. With the partial quotients $a_1, \ldots, a_n$ of $k/(h
\bmod k)$, the convergent matrix $\prod_i \left(\begin{smallmatrix} a_i &
1 \\ 1 & 0 \end{smallmatrix}\right) = \left(\begin{smallmatrix} m_{11} &
m_{12} \\ m_{21} & m_{22} \end{smallmatrix}\right)$ and the alternating
sum $t = a_1 - a_2 + \cdots \pm a_n$,

$$
12\, s(h, k) = t + \frac{m_{21} - (-1)^{n} m_{12}}{m_{11}} - 3 \cdot
[n \text{ odd}].
$$

The function evaluates this in machine words when $k$ fits in a `slong`
and with the half-GCD continued fraction code otherwise. Writing $s(h,k)
= p/q$ in lowest terms, $|s| < k/12$ and $q \mid 6k$.
:::
:::

## Internal functions {#internal}

These types and functions implement the half-GCD continued fraction
algorithm behind `fmpq_get_cfrac`, `fmpq_simplest_between`,
`fmpq_dedekind_sum` and rational reconstruction. They are declared in
`fmpq.h` for use by other modules and may change without notice.

::: type
_fmpz_mat22_t
_ui_mat22_t

A $2 \times 2$ matrix of `fmpz` (fields `_11`, `_12`, `_21`, `_22`) or of
`ulong`, with a field `det` holding its determinant $\pm 1$, or 0 if
unknown. Products of the elementary matrices $\left(\begin{smallmatrix}
q & 1 \\ 1 & 0 \end{smallmatrix}\right)$ have determinant $\pm 1$.
:::

::: function
void _fmpz_mat22_init(_fmpz_mat22_t M)
void _fmpz_mat22_clear(_fmpz_mat22_t M)
void _fmpz_mat22_one(_fmpz_mat22_t M)
int _fmpz_mat22_is_one(_fmpz_mat22_t M)
flint_bitcnt_t _fmpz_mat22_bits(const _fmpz_mat22_t N)

Initialise, clear, set to the identity, test for the identity, and return
the largest bit length of an entry.
:::

::: function
void _fmpz_mat22_rmul(_fmpz_mat22_t M, const _fmpz_mat22_t N)
void _fmpz_mat22_rmul_ui(_fmpz_mat22_t M, const _ui_mat22_t N)
void _fmpz_mat22_rmul_inv_ui(_fmpz_mat22_t M, const _ui_mat22_t N)
void _fmpz_mat22_rmul_elem(_fmpz_mat22_t M, const fmpz_t q)
void _fmpz_mat22_lmul_elem(_fmpz_mat22_t M, const fmpz_t q)
void _fmpz_mat22_rmul_inv_elem(_fmpz_mat22_t M, const fmpz_t q)

Set $M$ to $MN$, $MN$, $MN^{-1}$, $ME_q$, $E_q M$ or $ME_q^{-1}$, where
$E_q = \left(\begin{smallmatrix} q & 1 \\ 1 & 0 \end{smallmatrix}\right)$
and the inverses use $\det N = \pm 1$.
:::

::: function
void _fmpz_mat22_addmul_inv_vec(fmpz_t ya, fmpz_t yb, _fmpz_mat22_t N, fmpz_t xa, fmpz_t xb)
void _fmpz_mat22_addmul_inv_mat(fmpz_t A11, fmpz_t A12, fmpz_t A21, fmpz_t A22, _fmpz_mat22_t M, fmpz_t B11, fmpz_t B12, fmpz_t B21, fmpz_t B22)

Add $N^{-1}x$ to the vector $y$, or $M^{-1}B$ to the matrix $A$. The
inputs `xa`, `xb` and $B$ may be overwritten.
:::

::: type
_fmpq_cfrac_list_t

A growable vector of continued fraction terms (`array`, `length`,
`alloc`), with a `limit` on the number of terms to generate. Setting
`length` to $-1$ discards the terms, and setting `want_alt_sum` to $\pm 1$
accumulates their alternating sum in `alt_sum` (used by
`fmpq_dedekind_sum`).
:::

::: function
void _fmpq_cfrac_list_init(_fmpq_cfrac_list_t v)
void _fmpq_cfrac_list_clear(_fmpq_cfrac_list_t v)
void _fmpq_cfrac_list_fit_length(_fmpq_cfrac_list_t v, slong len)
void _fmpq_cfrac_list_push_back(_fmpq_cfrac_list_t v, const fmpz_t a)
void _fmpq_cfrac_list_push_back_zero(_fmpq_cfrac_list_t v)
void _fmpq_cfrac_list_append_ui(_fmpq_cfrac_list_t v, const ulong * a, slong n)

Memory management and appending of terms for `_fmpq_cfrac_list_t`.
:::

::: type
_fmpq_ball_t

A closed interval $[\mathit{left\_num}/\mathit{left\_den},
\mathit{right\_num}/\mathit{right\_den}]$, or the single point
$\mathit{left\_num}/\mathit{left\_den}$ if the field `exact` is set.
:::

::: function
void _fmpq_ball_init(_fmpq_ball_t x)
void _fmpq_ball_clear(_fmpq_ball_t x)
void _fmpq_ball_swap(_fmpq_ball_t x, _fmpq_ball_t y)

Memory management for `_fmpq_ball_t`.
:::

::: function
int _fmpq_ball_gt_one(const _fmpq_ball_t x)

Returns nonzero if the endpoints are canonical and the ball lies in $(1,
\infty)$.
:::

::: function
void _fmpq_ball_get_cfrac(_fmpq_cfrac_list_t s, _fmpz_mat22_t M, int needM, _fmpq_ball_t x)

Given a ball $x > 1$, append to $s$ continued fraction terms valid for
every number in $x$, multiply $M$ on the right by the corresponding
elementary matrices, and replace $x$ by its image under the inverse
transformation. With `needM` zero and an exact ball, $M$ may be left
unchanged.
:::

::: function
void _fmpq_hgcd(_fmpq_cfrac_list_t s, _fmpz_mat22_t M, fmpz_t x_num, fmpz_t x_den)

Half-GCD step: for $a = $ `x_num` $> b = $ `x_den` $> 0$, generate the
terms valid for every real number in the open interval $(a/(b+1),
(a+1)/b)$, append them to $s$, multiply $M$ on the right by their
matrices and replace $(a, b)$ by the reduced pair, in place.
:::

::: macro
FMPQ_RECONSTRUCT_HGCD_CUTOFF

The difference in limbs between the modulus and the numerator bound from
which `_fmpq_reconstruct_fmpz_2` uses the half-GCD algorithm instead of
Lehmer's algorithm (500).
:::
