---
title: mpn_mod.h
header: mpn_mod.h
description: Arithmetic in $\Z/n\Z$ for moduli of 2 to 16 limbs, with residues stored as fixed-length limb arrays, as a generic ring with fast vector, polynomial and matrix methods.
category: Integers mod n
order: 5
---

This module implements the ring $\Z/n\Z$ for a modulus $n$ of $\ell$
limbs with $2 \le \ell \le 16$, that is $2^{\bits} \le n < 2^{16 \bits}$
($2^{64} \le n < 2^{1024}$ on 64-bit machines). An element is a residue
in $[0, n)$ stored as an array of exactly $\ell$ limbs, least significant
limb first and zero-padded. Elements have no header and own no memory:
they can live on the stack, be copied with `memcpy` and be packed
contiguously in vectors, polynomials and matrices. This fills the gap
between `nmod`, whose residues are single words, and `fmpz_mod`, whose
residues are `fmpz` integers of arbitrary size.

The module is a generic ring in the sense of `gr`: the ring is a
`gr_ctx_t` created by `gr_ctx_init_mpn_mod`, functions return `gr` status
codes, and polynomials and matrices are `gr_poly_t` and `gr_mat_t` over
this context. Every function below is installed in the ring's method
table, so `gr_add`, `gr_poly_gcd` or `gr_mat_lu` called with an `mpn_mod`
context run the code documented here. The functions can also be called
directly, which saves the method table lookup. Compared with the
`fmpz`-based ring of `gr_ctx_init_fmpz_mod`, fixed-length residues make
the basecase algorithms faster (Euclid's algorithm and Gaussian
elimination with delayed reduction), avoid temporary allocations in the
divide-and-conquer algorithms, and make additions cheap enough to favour
Karatsuba and Waksman multiplication [@Joh2025gr].

The guide [](mpn_mod/usage) shows how to set up and use the ring.

::: summary algorithms/arithmetic/multiword-overview :::

::: summary algorithms/arithmetic/multiword-mulmod :::

::: summary algorithms/arithmetic/multiword-delayed-reduction :::

## Conventions {#conventions}

- $\ell$ is the number of limbs of the modulus $n$, so that
  $2^{(\ell - 1) \bits} \le n < 2^{\ell \bits}$, and $b$ is the number of
  bits of $n$. An element is a pointer to $\ell$ limbs (`nn_ptr` for
  outputs, `nn_srcptr` for inputs) holding a value in $[0, n)$. Functions
  require reduced inputs and return reduced outputs; this is not checked.
- Vectors, polynomials and matrices store their entries contiguously,
  $\ell$ limbs per entry. Strides and lengths count entries, not limbs.
- Every function that returns `int` returns a `gr` status code:
  `GR_SUCCESS`, `GR_DOMAIN` (no result exists, as for an inverse of a
  non-unit) or `GR_UNABLE` (a result may exist but cannot be computed).
- The scalar functions allow their output to alias any input. Polynomial
  functions with a leading underscore do not allow the output to alias an
  input unless stated otherwise.
- Primality of $n$ is never tested. The context stores a flag which is
  `T_UNKNOWN` until set with `gr_ctx_set_is_field`; algorithms that need a
  field check the flag or the *pretend field* flag (see
  [](#context)).

## Limits and types {#types}

::: macro
MPN_MOD_MIN_LIMBS
MPN_MOD_MAX_LIMBS

The permitted numbers of limbs of the modulus: currently $\ell \ge 2$ and
$\ell \le 16$, that is moduli up to 1024 bits on 64-bit machines and up
to 512 bits on 32-bit machines.

::: details
Single-limb moduli are excluded because `nmod` already handles them well,
and excluding them removes special cases from the code. The upper limit
keeps elements and temporary buffers (a small multiple of $\ell$ limbs)
small enough to allocate on the stack and keeps operations such as swaps
cheap compared with a pointer to heap memory; the tuning tables also only
cover this range. For larger moduli, use `fmpz_mod` (or
`gr_ctx_init_fmpz_mod`). The upper limit might be increased in the
future.

:::
:::

::: type
_mpn_mod_ctx_struct

The data of an `mpn_mod` ring, stored behind the data pointer of the
`gr_ctx_t`. The fields are `nlimbs` ($\ell$), `d` (the limbs of $n$),
`dnormed` ($n 2^{\nu}$, where $\nu$ is the number of leading zero bits of
the top limb of $n$, so that the top bit is set), `dinv` (the
$\ell$-limb inverse of `dnormed` computed by `flint_mpn_preinvn`), `norm`
($\nu$), `is_prime` (a `truth_t`), `pretend` (whether the ring pretends
to be a field) and `factor` (a factor of $n$ recorded while pretending,
or `NULL`). The arrays have room for `MPN_MOD_MAX_LIMBS` limbs, so the
structure has a fixed size of about 400 bytes.
:::

## Context objects {#context}

::: function
int gr_ctx_init_mpn_mod(gr_ctx_t ctx, const fmpz_t n)
int _gr_ctx_init_mpn_mod(gr_ctx_t ctx, nn_srcptr n, slong nlimbs)

Initialises *ctx* to the ring $\Z/n\Z$ with elements stored as arrays of
as many limbs as $n$. Returns `GR_DOMAIN` if $n \le 0$ and `GR_UNABLE`
if $n$ does not have between `MPN_MOD_MIN_LIMBS` and `MPN_MOD_MAX_LIMBS`
limbs; in both cases *ctx* is not initialised and must not be cleared.
The underscore version takes $n$ as `nlimbs` limbs, of which the top one
must be nonzero.

::: details
The context allocates its `_mpn_mod_ctx_struct` on the heap and must be
cleared with `gr_ctx_clear`. Initialisation costs one precomputed inverse
of $n$ (`flint_mpn_preinvn`). The primality flag starts as `T_UNKNOWN`.

A caller that handles all sizes typically uses `gr_ctx_init_nmod` for one
limb, this function for 2 to 16 limbs and `gr_ctx_init_fmpz_mod` beyond.

:::
:::

::: function
void gr_ctx_init_mpn_mod_randtest(gr_ctx_t ctx, flint_rand_t state)

Initialises *ctx* to $\Z/n\Z$ for a random modulus, for testing. With
probability 1/2 the modulus is one of a few primes of the form
$2^k \pm c$ for $k$ between 64 and 1024 (flagged as prime with
probability 1/2); otherwise it is a random integer with an admissible
number of limbs.
:::

::: function
int mpn_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
void mpn_mod_ctx_clear(gr_ctx_t ctx)

The `gr` methods for writing a description of the ring
(`Integers mod n (mpn)`) and for freeing the context. Use
`gr_ctx_write` and `gr_ctx_clear`.
:::

::: macro
MPN_MOD_CTX(ctx)
MPN_MOD_CTX_NLIMBS(ctx)
MPN_MOD_CTX_MODULUS(ctx)
MPN_MOD_CTX_MODULUS_NORMED(ctx)
MPN_MOD_CTX_MODULUS_PREINV(ctx)
MPN_MOD_CTX_NORM(ctx)
MPN_MOD_CTX_IS_PRIME(ctx)
MPN_MOD_CTX_PRETEND(ctx)
MPN_MOD_CTX_FACTOR(ctx)

Access the fields of the `_mpn_mod_ctx_struct` of *ctx*: the structure
itself, $\ell$, the limbs of $n$, of $n 2^{\nu}$ and of its precomputed
inverse, the shift $\nu$, the primality flag, the pretend flag and the
recorded factor. They can be used as lvalues, but the pretend flag and
the factor should only be changed through the functions below, which take
a lock.
:::

::: macro
MPN_MOD_CTX_MODULUS_BITS(ctx)

The number of bits $b$ of $n$, that is $\ell \bits - \nu$.
:::

::: function
truth_t mpn_mod_ctx_is_field(gr_ctx_t ctx)
int mpn_mod_ctx_set_is_field(gr_ctx_t ctx, truth_t is_field)

Return or set the primality flag of the modulus. `mpn_mod_ctx_set_is_field`
always returns `GR_SUCCESS`; nothing is checked. These are the
`gr_ctx_is_field` and `gr_ctx_set_is_field` methods (the flag also
answers `gr_ctx_is_integral_domain` and
`gr_ctx_is_unique_factorization_domain`).

::: details
Setting the flag to `T_TRUE` for a prime modulus lets the generic
algorithms for fields (polynomial factorisation and root finding, square
roots, and so on) run instead of returning `GR_UNABLE`. Setting it for a
composite modulus gives undefined results.

:::
:::

::: function
truth_t mpn_mod_ctx_is_pretend_field(gr_ctx_t ctx)
int mpn_mod_ctx_set_is_pretend_field(gr_ctx_t ctx, truth_t is_pretend_field)
int mpn_mod_ctx_recover_zero_divisor(nn_ptr res, gr_ctx_t ctx)

The methods behind `gr_ctx_is_pretend_field`,
`gr_ctx_set_is_pretend_field` and `gr_ctx_recover_zero_divisor`. While
the ring pretends to be a field, an inversion of a nonzero non-unit $x$
fails with `GR_UNABLE` instead of `GR_DOMAIN` and records
$\gcd(x, n)$, a nontrivial factor of $n$, in the context (only the first
one found is kept). `mpn_mod_ctx_recover_zero_divisor` sets *res* to the
recorded factor and returns `GR_SUCCESS`, or sets *res* to zero and
returns `GR_UNABLE` if none was recorded.

::: details
`mpn_mod_ctx_is_pretend_field` returns `T_TRUE` if the ring pretends to
be a field or is flagged as prime, and `T_FALSE` otherwise. Turning the
pretense off forgets the recorded factor. The flag and the factor are
protected by a global lock, so several threads may compute with the same
context. The guide [](mpn_mod/usage) has an example.

:::
:::

## Memory management and assignment {#assignment}

::: function
void mpn_mod_init(nn_ptr x, gr_ctx_t ctx)
void mpn_mod_clear(nn_ptr x, gr_ctx_t ctx)
void mpn_mod_swap(nn_ptr x, nn_ptr y, gr_ctx_t ctx)

`mpn_mod_init` sets the $\ell$ limbs at *x* to zero, `mpn_mod_clear` does
nothing (elements own no memory) and `mpn_mod_swap` swaps the limbs of
*x* and *y*.
:::

::: function
int mpn_mod_set(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
int mpn_mod_zero(nn_ptr res, gr_ctx_t ctx)
int mpn_mod_one(nn_ptr res, gr_ctx_t ctx)
int mpn_mod_neg_one(nn_ptr res, gr_ctx_t ctx)

Set *res* to a copy of *x*, to 0, to 1 and to $n - 1$.
:::

::: function
int mpn_mod_set_ui(nn_ptr res, ulong x, gr_ctx_t ctx)
int mpn_mod_set_si(nn_ptr res, slong x, gr_ctx_t ctx)
int mpn_mod_set_fmpz(nn_ptr res, const fmpz_t x, gr_ctx_t ctx)

Set *res* to $x \bmod n$, for any integer $x$.

::: details
Since $n > 2^{\bits}$, a `ulong` needs no reduction and a negative
`slong` becomes $n - |x|$. An `fmpz` with as many limbs as $n$ or more is
reduced with `mpn_mod_set_mpn` and negated if negative.

:::
:::

::: function
int mpn_mod_set_mpn(nn_ptr res, nn_srcptr x, slong xn, gr_ctx_t ctx)

Sets *res* to $x \bmod n$, where $x \ge 0$ is given by the `xn` limbs at
*x* ($\mathit{xn} \ge 0$; the top limb may be zero). This is the reduction
used throughout the module for unreduced sums and products. *res* may be
the same pointer as *x* but must not overlap it otherwise.

::: details
Values with fewer limbs than $n$, or with $\ell$ limbs and less than $n$,
are copied. Otherwise $x$ is shifted left by $\nu$ bits and divided by
$n 2^{\nu}$, with GMP's `mpn_divrem_2` for $\ell = 2$ and with
`flint_mpn_mod_preinvn` and the precomputed inverse of the context for
$\ell \ge 3$.

:::
:::

::: function
int mpn_mod_set_other(nn_ptr res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx)

Sets *res* to the element *v* of the ring *v_ctx*, converted to
$\Z/n\Z$. Elements of another `mpn_mod` ring and of an `fmpz_mod` ring
(`gr_ctx_init_fmpz_mod`) with the same modulus are copied directly; all
other cases go through `gr_generic_set_other` (integers, rationals with
invertible denominator, and so on). Returns `GR_DOMAIN` or `GR_UNABLE`
when there is no such conversion.
:::

::: function
int mpn_mod_get_fmpz(fmpz_t res, nn_srcptr x, gr_ctx_t ctx)

Sets *res* to the integer in $[0, n)$ represented by *x*.
:::

::: function
int mpn_mod_rand(nn_ptr res, flint_rand_t state, gr_ctx_t ctx)
int mpn_mod_randtest(nn_ptr res, flint_rand_t state, gr_ctx_t ctx)

Set *res* to a random element. `mpn_mod_rand` draws uniformly from
$[0, n)$ (`fmpz_randm`); `mpn_mod_randtest` uses `fmpz_randtest_mod`,
which favours special values such as $0$, $1$ and $n - 1$ and values with
long runs of zero or one bits, for testing.
:::

::: function
int mpn_mod_write(gr_stream_t out, nn_srcptr x, gr_ctx_t ctx)

Writes *x* in decimal, as an integer in $[0, n)$.
:::

## Comparison {#comparison}

::: function
truth_t mpn_mod_is_zero(nn_srcptr x, gr_ctx_t ctx)
truth_t mpn_mod_is_one(nn_srcptr x, gr_ctx_t ctx)
truth_t mpn_mod_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
truth_t mpn_mod_equal(nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)

Return `T_TRUE` if $x = 0$, $x = 1$, $x = n - 1$ or $x = y$, and `T_FALSE`
otherwise. Since residues are reduced, these compare limbs.
:::

## Arithmetic {#arithmetic}

::: function
int mpn_mod_neg(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
int mpn_mod_add(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
int mpn_mod_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)

Set *res* to $-x$, $x + y$ and $x - y$ modulo $n$, with one conditional
correction. Two-limb moduli use the inline functions
`flint_mpn_negmod_2`, `flint_mpn_addmod_2` and `flint_mpn_submod_2`.
:::

::: function
int mpn_mod_add_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
int mpn_mod_sub_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
int mpn_mod_add_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
int mpn_mod_sub_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
int mpn_mod_add_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int mpn_mod_sub_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)

Set *res* to $x + y$ or $x - y$ modulo $n$, for reduced $x$ and any
integer $y$.
:::

::: function
int mpn_mod_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
int mpn_mod_sqr(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)

Set *res* to $xy$ and $x^2$ modulo $n$.

::: details
The product $(x 2^{\nu}) y$ of $2\ell$ limbs is reduced modulo
$n 2^{\nu}$ with the precomputed inverse of the context (a quotient
estimate from one high product and one low product, followed by at most
a few subtractions) and shifted back; see
[](algorithms/arithmetic/multiword-mulmod). Two-limb moduli use the
inline code of `flint_mpn_mulmod_preinvn_2`. When $x$ and $y$ are the
same pointer, the product is computed by squaring. This is the algorithm
of `flint_mpn_mulmod_preinvn`, with the shift applied to an input instead
of the output.

:::
:::

::: function
int mpn_mod_mul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
int mpn_mod_mul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
int mpn_mod_mul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)

Set *res* to $xy \bmod n$, for reduced $x$ and any integer $y$. The
product $x |y|$ is formed exactly and reduced with `mpn_mod_set_mpn` (an
`fmpz` with more limbs than $n$ is reduced first).
:::

::: function
int mpn_mod_addmul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
int mpn_mod_addmul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
int mpn_mod_addmul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
int mpn_mod_addmul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int mpn_mod_submul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
int mpn_mod_submul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
int mpn_mod_submul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
int mpn_mod_submul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)

Set *res* to $res + xy$ or $res - xy$ modulo $n$. The product is reduced
before it is added.
:::

::: function
int mpn_mod_fmma(nn_ptr res, nn_srcptr x1, nn_srcptr y1, nn_srcptr x2, nn_srcptr y2, gr_ctx_t ctx)

Sets *res* to $x_1 y_1 + x_2 y_2 \bmod n$ with a single reduction of the
unreduced sum (`flint_mpn_fmmamod_preinvn`), which costs little more
than one modular multiplication.
:::

## Inverses, division and powers {#inverse}

::: function
int mpn_mod_inv(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)

Sets *res* to $x^{-1} \bmod n$ and returns `GR_SUCCESS` if
$\gcd(x, n) = 1$. Otherwise returns `GR_DOMAIN`, or `GR_UNABLE` if $x$ is
nonzero and the ring pretends to be a field (recording
$\gcd(x, n)$, see [](#context)); *res* is then left unchanged.

::: details
The inverse is computed with GMP's extended gcd `mpn_gcdext`, except for
$x = \pm 1$, which are their own inverses.

:::
:::

::: function
int mpn_mod_div(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)

Sets *res* to $x y^{-1} \bmod n$, computing the inverse with
`mpn_mod_inv`. If $y$ is not invertible, sets *res* to zero and returns
the status of `mpn_mod_inv`, even if some $z$ with $yz = x$ exists.
:::

::: function
int mpn_mod_pow_ui(nn_ptr res, nn_srcptr x, ulong e, gr_ctx_t ctx)
int mpn_mod_pow_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t e, gr_ctx_t ctx)

Set *res* to $x^e \bmod n$, with $0^0 = 1$. A negative exponent requires
$x$ to be invertible and otherwise gives the status of `mpn_mod_inv`.

::: details
Exponents $e \le 64$ (`pow_ui`) or of at most 6 bits (`pow_fmpz`) use
binary exponentiation with `mpn_mod_mul` (`gr_generic_pow_ui`,
`gr_generic_pow_fmpz`); larger exponents call GMP's `mpz_powm`.
Negative exponents go through the generic code.

:::
:::

## Square roots {#sqrt}

::: function
truth_t mpn_mod_is_square(nn_srcptr x, gr_ctx_t ctx)
int mpn_mod_sqrt(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)

Test whether *x* is a square, and compute a square root of *x*.

::: details
$0$ and $1$ are squares for every modulus, and `mpn_mod_sqrt` returns
them unchanged. Otherwise both functions use the Jacobi symbol
$\left(\frac{x}{n}\right)$ for odd $n$ (`flint_mpn_is_square_mod`). A
symbol of $-1$ proves that $x$ is not a square, whether or not $n$ is
prime: `mpn_mod_is_square` returns `T_FALSE` and `mpn_mod_sqrt` sets
*res* to zero and returns `GR_DOMAIN`. A symbol of $+1$ gives `T_TRUE`
if the context is flagged as prime and `T_UNKNOWN` otherwise; for even
$n$, `mpn_mod_is_square` returns `T_UNKNOWN`.

`mpn_mod_sqrt` computes roots only for odd moduli flagged as prime
(`gr_ctx_set_is_field`), with `flint_mpn_sqrtmod_preinv` and the
precomputed inverse of the context (an exponentiation, Atkin's formula,
Tonelli–Shanks or Müller's variant of the Cipolla–Lehmer algorithm
depending on the 2-adic valuation of $n - 1$; see
[](algorithms/number-theory/sqrt-mod-p)). For even moduli, and for odd
moduli not flagged as prime when the Jacobi symbol is $+1$, it returns
`GR_UNABLE`. If the algorithm fails, which can only happen if the flag
was wrong, it also returns `GR_UNABLE`.

:::
:::

## Vector functions {#vec}

These are the vector methods of the ring, used by `_gr_vec_add`,
`_gr_vec_dot` and the like, and through them by the generic polynomial
and matrix algorithms. They work on arrays of `len` elements of $\ell$
limbs each.

::: function
int _mpn_mod_vec_zero(nn_ptr res, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_clear(nn_ptr res, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_set(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx)
void _mpn_mod_vec_swap(nn_ptr vec1, nn_ptr vec2, slong len, gr_ctx_t ctx)

Zero, clear (a no-op), copy and swap vectors, limb by limb.
:::

::: function
int _mpn_mod_vec_neg(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_add(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)

Elementwise negation, addition, subtraction and multiplication. Two-limb
moduli use inlined fixed-size code with the modulus kept in registers.
:::

::: function
int _mpn_mod_vec_mul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
int _mpn_mod_scalar_mul_vec(nn_ptr res, nn_srcptr y, nn_srcptr x, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_addmul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
int _mpn_mod_vec_submul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)

Set *res* to $x y$, $y x$, $res + x y$ and $res - x y$ elementwise, for a
vector $x$ and a scalar $y$.

::: details
Multiplication by a fixed $y$ is preconditioned when the vector is long
enough to amortise the precomputation:
`flint_mpn_mulmod_want_precond` chooses between no precomputation,
Shoup's method (a quotient precomputed for $y$, used for $\ell \le 10$)
and the matrix method (the residues $y 2^{\nu} B^i \bmod n 2^{\nu}$
precomputed, used for larger $\ell$), from 4 elements on. For two limbs,
Shoup's method is used from 4 elements on (6 if the top bit of $n$ is
set), with inlined code; see
[](algorithms/arithmetic/multiword-mulmod#fixed).

:::
:::

::: function
int _mpn_mod_vec_dot(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_dot_rev(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx)
int _mpn_mod_vec_dot_strided(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, slong stride1, nn_srcptr vec2, slong stride2, slong len, gr_ctx_t ctx)

Set *res* to $s + \sum_i a_i b_i$, or $s - \sum_i a_i b_i$ if `subtract`
is nonzero, where $s$ is *initial* (zero if *initial* is `NULL`),
$a_i$ = `vec1[i]` and $b_i$ = `vec2[i]`, `vec2[len - 1 - i]` or the
entries of the strided vectors (strides counted in elements, possibly
negative).

::: details
The sum of the full $2\ell$-limb products is accumulated exactly in
$\lceil (2 b + \lfloor \log_2 \mathit{len} \rfloor + 1) / \bits \rceil$
limbs and reduced once at the end
([delayed reduction](algorithms/arithmetic/multiword-delayed-reduction#dot)).
For $\ell \le 4$ the accumulation uses the generated fixed-size kernels
of `mpn_extras` (`flint_mpn_dot_tab` and its strided and reversed
variants); for larger $\ell$, `flint_mpn_mul_n` and `mpn_add_n` with one
carry limb. Lengths 1 and 2 use `mpn_mod_mul` and `mpn_mod_fmma`.

:::
:::

## Matrix algorithms {#mat}

All `gr_mat` functions work over this ring. The functions below are the
matrix methods of the ring, which the generic algorithms call (for
example, `gr_mat_mul` calls `mpn_mod_mat_mul`, and the block recursive
LU decomposition and triangular solving use it for their products).
Matrices are `gr_mat_t` over an `mpn_mod` context.

::: function
int mpn_mod_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int mpn_mod_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int mpn_mod_mat_mul_multi_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Set $C = AB$. `mpn_mod_mat_mul` chooses between classical multiplication
with delayed reduction (`gr_mat_mul_classical` with the dot product
above), Waksman's algorithm with delayed reduction and multimodular
multiplication via `nmod_mat_mul`. See [](mpn_mod/mat-mul).
:::

::: function
int mpn_mod_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int mpn_mod_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)

Solve $LX = B$ or $UX = B$ for a triangular matrix with invertible
diagonal (unit diagonal, not read, if `unit` is nonzero), choosing
`gr_mat_nonsingular_solve_tril_classical` (resp. `_triu_`) when $B$ has
fewer rows or columns than a cutoff depending on $\ell$ (700 for
$\ell = 2$, 200 for $\ell = 3$, down to 20 for $\ell \ge 12$), and the
block recursive versions otherwise, whose products use `mpn_mod_mat_mul`.
:::

::: function
int mpn_mod_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
int mpn_mod_mat_lu_classical_delayed(slong * res_rank, slong * P, gr_mat_t A, const gr_mat_t A_in, int rank_check, gr_ctx_t ctx)

LU decomposition with the conventions of `gr_mat_lu`: $PA = LU$ with
the permutation in *P*, $L$ unit lower triangular and $U$ in row echelon
form, both stored in the output, and the rank in `*rank`. If `rank_check`
is nonzero, the decomposition stops with rank 0 as soon as the matrix is
found to be singular. Returns `GR_DOMAIN` (or `GR_UNABLE` when pretending
to be a field) if a nonzero pivot is not invertible.

::: details
`mpn_mod_mat_lu_classical_delayed` is Gaussian elimination in which each
entry of the active submatrix is kept unreduced in $2\ell + 1$ limbs and
reduced only when it enters the pivot column or the pivot row
([](algorithms/arithmetic/multiword-delayed-reduction#lu)); the update of
a row is a single long `mpn_add_n` over its unreduced entries.
`mpn_mod_mat_lu` uses `gr_mat_lu_classical` for $\min(m, n)$ below a
first cutoff (23 for $\ell = 2$, 5 to 8 for larger $\ell$), the delayed
version below a second cutoff (200 for $\ell = 2$, decreasing to 20 for
$\ell \ge 9$), and `gr_mat_lu_recursive` above it.

:::
:::

::: function
int mpn_mod_mat_det(nn_ptr res, const gr_mat_t A, gr_ctx_t ctx)

Sets *res* to the determinant of the square matrix *A*, using cofactor
expansion up to $4 \times 4$, the division-free Berkowitz algorithm for
$5 \times 5$, and LU decomposition otherwise. If the LU decomposition
meets a non-invertible pivot (for a composite modulus), the determinant
is computed with the Berkowitz algorithm instead, so the function works
for every modulus. Returns `GR_DOMAIN` if *A* is not square.
:::

::: function
int _mpn_mod_mat_charpoly(nn_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

Sets the $m + 1$ coefficients at *res* to the characteristic polynomial
of the $m \times m$ matrix *mat*, with the Berkowitz algorithm for
$m \le 11$ and Danilevsky's algorithm (`_gr_mat_charpoly_danilevsky`)
for larger $m$. Danilevsky's algorithm divides by pivots, so for a
composite modulus it can fail; the function then falls back to the
division-free Berkowitz algorithm, and always succeeds. This is the method behind `gr_mat_charpoly`.
:::

::: function
int mpn_mod_mat_reduce_row(slong * column, gr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx)

Reduces row $m$ of *A* against the previous rows, with the conventions of
`gr_mat_reduce_row`; used by `gr_mat_minpoly_field`. The row is kept
unreduced in $2\ell + 1$ limbs per entry during the elimination. Returns
`GR_DOMAIN` (or `GR_UNABLE` when pretending) if a new pivot is not
invertible.
:::

## Polynomial algorithms {#poly}

All `gr_poly` functions work over this ring. The functions below are the
polynomial methods of the ring (multiplication, series inversion and
division, division with remainder, gcd), called by `gr_poly_mul`,
`gr_poly_divrem`, `gr_poly_gcd` and the many generic algorithms built on
them. Polynomials are arrays of coefficients, constant term first, with
$\ell$ limbs per coefficient. Lengths must be positive.

### Multiplication {#poly-mul}

A *middle product* (`mulmid`) computes the coefficients of index
$\mathit{nlo} \le i < \mathit{nhi}$ of the product and writes the
$\mathit{nhi} - \mathit{nlo}$ coefficients to *res*; it requires
$0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1} + \mathit{len2} - 1$.
A low product (`mullow`) is the middle product with $\mathit{nlo} = 0$
and $\mathit{nhi} = \mathit{len}$. All variants detect squaring
(`poly1 == poly2` and equal lengths). [](mpn_mod/poly-mul)
describes the algorithms and how they are chosen.

::: function
int _mpn_mod_poly_mullow(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
int _mpn_mod_poly_mulmid(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

Low and middle products with automatic choice of algorithm. These are the
`_gr_poly_mullow` and `_gr_poly_mulmid` methods of the ring.
:::

::: function
int _mpn_mod_poly_mullow_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
int _mpn_mod_poly_mulmid_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

Schoolbook multiplication: each output coefficient is an unreduced dot
product of fixed limb size, reduced once.
:::

::: function
int _mpn_mod_poly_mullow_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, slong cutoff, gr_ctx_t ctx)
int _mpn_mod_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, slong cutoff, gr_ctx_t ctx)

Karatsuba multiplication on unreduced integer coefficients
(`_flint_mpn_poly_mulmid_karatsuba`), recursing down to schoolbook
multiplication when the shorter operand is shorter than *cutoff*; a
*cutoff* of $-1$ selects the tuned value. The coefficients are reduced
only at the end.
:::

::: function
int _mpn_mod_poly_mulmid_mpn(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

Middle product without intermediate reductions, choosing between the
schoolbook algorithm and Karatsuba with the tuned cutoff of
`_flint_mpn_poly_mulmid`, with output sizes given by
`_flint_mpn_poly_mulmid_slimbs`. This is what `_mpn_mod_poly_mulmid`
uses for short and medium lengths.
:::

::: function
int _mpn_mod_poly_mullow_KS(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
int _mpn_mod_poly_mulmid_KS(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

Multiplication by Kronecker substitution: the coefficients are packed
into integers with
$2b + \lfloor \log_2 \min(\mathit{len1}, \mathit{len2}) \rfloor + 1$ bits
per coefficient, the integers are multiplied with `flint_mpn_mul` and the
output coefficients are extracted and reduced.
:::

::: function
int _mpn_mod_poly_mullow_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
int _mpn_mod_poly_mulmid_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

Multiplication with the small-prime FFT of `fft_small`: the coefficients
are reduced modulo up to eight 50-bit primes, multiplied by
number-theoretic transforms and reconstructed by Chinese remaindering.
Returns `GR_UNABLE`, writing nothing, if `fft_small` is not available or
the product of eight primes is too small for the coefficients of the
integer product, which limits $n$ to about 190 bits (196 bits for very
short, 189 bits for lengths around $10^5$).
:::

### Division {#poly-div}

::: function
int _mpn_mod_poly_inv_series(nn_ptr Q, nn_srcptr B, slong lenB, slong len, gr_ctx_t ctx)
int _mpn_mod_poly_div_series(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, slong len, gr_ctx_t ctx)

Power series inversion $Q = B^{-1} \bmod x^{\mathit{len}}$ and division
$Q = A / B \bmod x^{\mathit{len}}$, with the conventions of
`_gr_poly_inv_series` and `_gr_poly_div_series` (the constant term of $B$
must be invertible). They use the basecase algorithms up to a tuned
length of $B$ (between 66 and 370 depending on the size of $n$; always
for $\mathit{lenB} \le 20$) and Newton iteration above it.
:::

::: function
int _mpn_mod_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)
int _mpn_mod_poly_div(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)

Division with remainder $A = QB + R$ and the quotient alone, with the
conventions of `_gr_poly_divrem` ($\mathit{lenA} \ge \mathit{lenB} \ge 1$,
*Q* with room for $\mathit{lenA} - \mathit{lenB} + 1$ and *R* for
$\mathit{lenB} - 1$ coefficients, leading coefficient of $B$
invertible). `_mpn_mod_poly_divrem` uses the delayed-reduction basecase
when $\mathit{lenA}$, $\mathit{lenB}$ or $\mathit{lenA} - \mathit{lenB}$
is at most a cutoff (between 38 and 166 depending on the size of $n$),
and Newton division otherwise; `_mpn_mod_poly_div` uses
`_gr_poly_div_basecase` up to the series division cutoff and
`_gr_poly_div_newton` above it. Return `GR_DOMAIN` (or `GR_UNABLE` when
pretending) if the leading coefficient of $B$ is not invertible.
:::

::: function
int _mpn_mod_poly_divrem_basecase(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)
int _mpn_mod_poly_divrem_basecase_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx)

Schoolbook division with remainder. The `preinv1` version takes the
inverse `invL` of the leading coefficient of $B$; the other computes it
and returns `GR_DOMAIN` (or `GR_UNABLE`) if it does not exist.

::: details
The partial remainder is kept unreduced in
$\lceil (2b + \lfloor \log_2 (\mathit{lenA} - \mathit{lenB} + 1) \rfloor + 1) / \bits \rceil$
limbs per coefficient, so that each quotient coefficient costs
$\mathit{lenB} - 1$ plain multiplications and additions and one reduction
([](algorithms/arithmetic/multiword-delayed-reduction#divrem)).
The cases $\mathit{lenB} = 1$, $\mathit{lenA} = \mathit{lenB}$ and
$\mathit{lenA} = \mathit{lenB} + 1$ (the step of Euclid's algorithm)
have special code; a monic divisor saves the multiplications by `invL`.

:::
:::

::: function
int _mpn_mod_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx)
int _mpn_mod_poly_divrem_q1_preinv1_fmma(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx)
int _mpn_mod_poly_divrem_q1_preinv1_fmma_precond(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx)
int _mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx)

Division with remainder in the case $\mathit{lenA} = \mathit{lenB} + 1$,
where the quotient $q_1 x + q_0$ has two coefficients, given the inverse
`invL` of the leading coefficient of $B$. Require $\mathit{lenB} \ge 2$.

::: details
Each coefficient of the remainder is $a_i - q_1 b_{i-1} - q_0 b_i$. The
`fmma` version computes it with one `flint_mpn_fmmamod_preinvn`; the
`fmma_precond` version precomputes $q_0$ and $q_1$ for the matrix method
of preconditioned multiplication; the `karatsuba_precond` version
shares one product between two consecutive coefficients through
$(q_0 + q_1)(b_{i-1} + b_i)$, saving a quarter of the multiplications
([](algorithms/arithmetic/multiword-delayed-reduction#euclid)).
`_mpn_mod_poly_divrem_q1_preinv1` chooses the `fmma` version for
$\mathit{lenB} < 10$ or moduli of at most 576 bits, the Karatsuba version
for long divisors or large moduli, and the `fmma_precond` version in
between.

:::
:::

### GCD {#poly-gcd}

::: function
int _mpn_mod_poly_gcd(nn_ptr G, slong * lenG, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)
int _mpn_mod_poly_xgcd(slong * lenG, nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)

Greatest common divisor and extended gcd, with the conventions of
`_gr_poly_gcd` and `_gr_poly_xgcd`. They use the Euclidean algorithm
(`_gr_poly_gcd_euclidean`, `_gr_poly_xgcd_euclidean`), whose division
steps are `_mpn_mod_poly_divrem_q1_preinv1`, if either length is below
240, and the half-gcd algorithm (`_gr_poly_gcd_hgcd`, `_gr_poly_xgcd_hgcd`
with inner cutoff 80 and outer cutoff 240) otherwise; see
[](algorithms/polynomials/half-gcd). The modulus should be prime, or the
ring should pretend to be a field: a non-invertible leading coefficient
gives `GR_DOMAIN` (`GR_UNABLE` when pretending, recording a factor of
$n$).
:::

## Transformed polynomials {#transformed}

With `fft_small` available, the ring supports the transformed
representation of polynomials of `gr_ctx_init_gr_poly_transformed_repr`
and `gr_ctx_init_gr_poly_transformed_cyclic_repr`: polynomials of bounded
length are kept as their transforms modulo the `fft_small` primes, so
that sums of products can be accumulated pointwise and converted back
once. `gr_poly_preinv_set` uses it for dense polynomial moduli of length
at least 160. Like the `fft_small` multiplication, it is only available
for moduli of up to about 190 bits.

::: function
int _gr_mpn_mod_ctx_init_transformed_poly_repr(gr_ctx_t ctx, gr_ctx_t base, slong len_bound, slong terms_bound, const struct gr_transformed_poly_workload_struct * workload)
int _gr_mpn_mod_ctx_init_transformed_poly_cyclic_repr(gr_ctx_t ctx, gr_ctx_t base, slong * len, slong terms_bound, const struct gr_transformed_poly_workload_struct * workload)

The methods of the ring behind `gr_ctx_init_gr_poly_transformed_repr` and
`gr_ctx_init_gr_poly_transformed_cyclic_repr`. They initialise *ctx* to a
ring of transformed polynomials over the `mpn_mod` ring *base*, either
for polynomials of length up to `len_bound` (linear), or for products
modulo $x^N - 1$, where $N$ is the least power of two that is at least
`*len` and at least 256 and is written back to `*len` (cyclic). Each
element may accumulate up to `terms_bound` products of two converted
polynomials. Return `GR_UNABLE` if `fft_small` is not available, the
modulus is too large, or the workload model predicts no gain over
ordinary multiplication or exceeds the memory budget
`flint_fft_small_max_transformed_ring_size`.

::: details
The prime count is provisioned for values up to
$4 \cdot \mathit{terms\_bound} \cdot n^2$: one factor 2 for a bias that
makes the signed accumulated values nonnegative and one because the
Chinese remaindering requires values below half the product of the
primes. Unlike the `nmod` version, there is no single-prime path (the
modulus never equals an FFT prime), so the multiplicative depth is
capped at 2. The context keeps a pointer to *base*, which must outlive
it. A context is read-only after construction and may be shared by
several threads.

The truncation granularity of the transforms (whole blocks of 256 points
for long transforms, the full transform length below one block) is given
by the macros `_fft_small_len_trunc` and `_fft_small_op_trunc` of
`fft_small.h`.

:::
:::

::: function
void fft_small_fft_mpn_mod(fft_small_op_t X, nn_srcptr a, ulong an, ulong itrunc, gr_ctx_t ctx, const fft_small_plan_t P)
void fft_small_export_mpn_mod_range(nn_ptr z, const fft_small_op_t X, ulong zl, ulong zh, gr_ctx_t ctx, const fft_small_plan_t P)

Conversions between `mpn_mod` coefficient vectors and `fft_small`
operands, the analogues of `fft_small_fft_nmod` and
`fft_small_export_nmod_range`, declared only when `FLINT_HAVE_FFT_SMALL`
is set. `fft_small_fft_mpn_mod` reduces the `an` coefficients at *a*
modulo each prime of the plan *P* and transforms them with truncation
`itrunc` (at least `an`, at most the plan's `ztrunc`).
`fft_small_export_mpn_mod_range` inverse-transforms *X*, reconstructs the
coefficients of index $\mathit{zl} \le i < \mathit{zh}$ (within the
plan's window) by Chinese remaindering, reduces them modulo $n$ and
writes the $\mathit{zh} - \mathit{zl}$ results to *z*; the reconstructed
integers must be nonnegative and below half the product of the primes.
:::

## Internal functions {#internal}

::: function
int _mpn_mod_nonunit(nn_srcptr x, nn_srcptr g, slong gsize, gr_ctx_t ctx)

Returns the status of a failed inversion of the residue *x*, given
$g = \gcd(x, n)$ in `gsize` limbs: `GR_DOMAIN` if $x = 0$ or the ring
does not pretend to be a field; otherwise `GR_UNABLE`, after recording
$g$ as the zero divisor of the context if none was recorded and
$1 < g < n$.
:::

::: function
int _mpn_mod_vec_rand(nn_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx)

Sets the `len` elements at *res* to uniformly random residues with
`mpn_mod_rand`.
:::
