---
title: fmpz_mod.h
header: fmpz_mod.h
description: Arithmetic in $\Z/n\Z$ for a modulus $n$ of any size, with residues stored as `fmpz` integers, and discrete logarithms modulo primes.
category: Integers mod n
order: 3
---

This module implements arithmetic in $\Z/n\Z$ for an arbitrary modulus
$n \ge 1$. Residues are ordinary `fmpz_t` integers in the canonical range
$[0, n)$, and the modulus is held in a context object `fmpz_mod_ctx_t`,
which every function takes as its last argument. The context also selects
the fastest code for the size of $n$: `nmod` arithmetic for one-word
moduli, an inline two-limb Barrett reduction for moduli below
$2^{2\bits}$, and `fmpz` division (with a precomputed inverse for very
large moduli) beyond. The module is the scalar layer of `fmpz_mod_vec`,
`fmpz_mod_poly`, `fmpz_mod_mat`, `fmpz_mod_mpoly` and the finite fields
of `fq` with a multiword characteristic. It also provides discrete
logarithms modulo primes by the Pohlig–Hellman method.

For moduli of 2 to 16 limbs, polynomial and matrix computations are
faster with the generic ring `mpn_mod`, which stores residues as
fixed-length limb arrays; see [](algorithms/arithmetic/multiword-overview).
The generic ring with `fmpz` residues is `gr_ctx_init_fmpz_mod`. For a
single operation modulo $n$ without a context, the functions of `fmpz`
described in [](fmpz/modular) suffice.

The guide [](fmpz_mod/usage) shows the typical use.

::: summary algorithms/arithmetic/multiword-overview :::

## Conventions {#conventions}

- $n$ is the modulus of the context *ctx*. A residue is *canonical* if it
  lies in $[0, n)$.
- Unless stated otherwise, inputs must be canonical and outputs are
  canonical. Functions with an `fmpz`, `ulong` or `slong` operand
  documented as "any integer" reduce it themselves. Canonicity is not
  checked, except by `FLINT_ASSERT` in builds with `--enable-assert`.
- Canonical residues are compared with `fmpz_equal`, `fmpz_is_zero` and
  `fmpz_is_one` (except that for $n = 1$ the residue 0 is also one; see
  `fmpz_mod_is_one`).
- Outputs may alias inputs.

## Types {#types}

::: type
fmpz_mod_ctx_struct
fmpz_mod_ctx_t

The context for arithmetic modulo $n$. An `fmpz_mod_ctx_t` is an array of
one `fmpz_mod_ctx_struct`, passed by reference. The fields are `n` (the
modulus), function pointers `add_fxn`, `sub_fxn` and `mul_fxn` chosen
by `fmpz_mod_ctx_init`, `mod` (an `nmod_t` for one-word moduli),
`n_limbs` and `ninv_limbs` (the modulus and the three-limb Barrett
inverse $\lfloor 2^{4\bits}/n \rfloor$ for two-limb moduli) and
`ninv_huge` (a precomputed inverse for moduli of more than 20 limbs, or
`NULL`). Use `fmpz_mod_ctx_modulus` to read $n$.
:::

## Context {#context}

::: function
void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n)
void fmpz_mod_ctx_init_ui(fmpz_mod_ctx_t ctx, ulong n)

Initialises *ctx* for arithmetic modulo $n$. Throws an exception if
$n \le 0$.

::: details
The context records which code to use for $n$:

| Modulus | `fmpz_mod_add`, `_sub`, `_mul` use |
|---|---|
| $n < 2^{\bits}$ | `nmod` arithmetic (`_fmpz_mod_add1`, ...) |
| $n = 2^{\bits}$ | wrapping word arithmetic (`_fmpz_mod_add2s`, ...) |
| $2^{\bits} < n < 2^{2\bits}$ | two-limb arithmetic with a Barrett inverse (`_fmpz_mod_add2`, ...) |
| more than two limbs | `fmpz` arithmetic (`_fmpz_mod_addN`, ...) |

For moduli of more than 20 limbs it also precomputes an inverse for
`fmpz_fdiv_r_preinvn`, which multiplications then use instead of
`fmpz_mod`. See [](algorithms/arithmetic/multiword-mulmod#two-limb) for
the two-limb reduction.

:::
:::

::: function
void fmpz_mod_ctx_init_rand_bits(fmpz_mod_ctx_t ctx, flint_rand_t state, flint_bitcnt_t max_bits)
void fmpz_mod_ctx_init_rand_bits_prime(fmpz_mod_ctx_t ctx, flint_rand_t state, flint_bitcnt_t max_bits)

Initialise *ctx* with a random modulus, for testing:
$2 + m$ with $m$ a random integer of at most `max_bits` bits (from
`fmpz_randtest_unsigned`), respectively the least probable prime
(`fmpz_nextprime` without proof) greater than such an $m$. The modulus
can have slightly more than `max_bits` bits.
:::

::: function
void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)

Frees the memory used by *ctx*.
:::

::: function
const fmpz * fmpz_mod_ctx_modulus(const fmpz_mod_ctx_t ctx)

Returns a pointer to the modulus $n$, valid as long as the context is not
cleared or changed.

::: details
This replaces the removed function
`fmpz_mod_ctx_get_modulus_mpz_read_only`, whose name now produces a
compile-time error.

:::
:::

::: function
void fmpz_mod_ctx_set_modulus(fmpz_mod_ctx_t ctx, const fmpz_t n)
void fmpz_mod_ctx_set_modulus_ui(fmpz_mod_ctx_t ctx, ulong n)

Changes the modulus of the initialised context *ctx* to $n$, which must
be positive. Residues computed for the old modulus are not reduced.
:::

## Conversion and canonical form {#conversion}

::: function
void fmpz_mod_set_fmpz(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
void fmpz_mod_set_ui(fmpz_t a, ulong b, const fmpz_mod_ctx_t ctx)
void fmpz_mod_set_si(fmpz_t a, slong b, const fmpz_mod_ctx_t ctx)

Sets $a$ to $b \bmod n$ in $[0, n)$, for any integer $b$.
:::

::: function
int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
void fmpz_mod_assert_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)

`fmpz_mod_is_canonical` returns 1 if $0 \le a < n$ and 0 otherwise.
`fmpz_mod_assert_canonical` throws an exception if $a$ is not canonical,
in every build.
:::

## Comparison {#comparison}

::: function
int fmpz_mod_is_one(const fmpz_t a, const fmpz_mod_ctx_t ctx)

Returns 1 if $a$ is 1 modulo $n$, that is, if $a = 1$ or $n = 1$, and 0
otherwise.
:::

::: function
int fmpz_mod_equal_fmpz(const fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
int fmpz_mod_equal_si(const fmpz_t a, slong b, const fmpz_mod_ctx_t ctx)

Return 1 if the canonical residue $a$ is congruent to the integer $b$
modulo $n$ (any $b$), and 0 otherwise.
:::

## Arithmetic {#arithmetic}

::: function
void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

Set $a$ to $b + c$, $b - c$ and $-b$ modulo $n$, for canonical $b$ and
$c$.
:::

::: function
void fmpz_mod_add_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_add_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_add_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_sub_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_sub_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_sub_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx)

Set $a$ to $b + c$ or $b - c$ modulo $n$, for canonical $b$ and any
integer $c$. The sum is computed exactly and reduced with `fmpz_mod`.
:::

::: function
void fmpz_mod_fmpz_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_ui_sub(fmpz_t a, ulong b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_si_sub(fmpz_t a, slong b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

Set $a$ to $b - c$ modulo $n$, for any integer $b$ and canonical $c$.
:::

::: function
void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

Sets $a$ to $bc$ modulo $n$, for canonical $b$ and $c$, with the function
selected by the context: `nmod_mul` for one-word moduli, a two-limb
[Barrett reduction](algorithms/arithmetic/multiword-mulmod#two-limb) with
at most two corrections for two-limb moduli, and `fmpz_mul` followed by
`fmpz_mod` (or `fmpz_fdiv_r_preinvn` above 20 limbs) otherwise.
:::

::: function
void fmpz_mod_mul_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mul_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mul_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx)

Set $a$ to $bc$ modulo $n$, for canonical $b$ and any integer $c$. The
product is computed exactly and reduced.
:::

::: function
void fmpz_mod_addmul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_mod_ctx_t ctx)

Sets $a$ to $b + cd$ modulo $n$. The sum is computed exactly with `fmpz`
arithmetic and reduced once, so $b$, $c$ and $d$ may be any integers.
:::

## Inversion and division {#inversion}

::: function
int fmpz_mod_is_invertible(const fmpz_t a, const fmpz_mod_ctx_t ctx)

Returns 1 if $\gcd(a, n) = 1$, that is, if $a$ is invertible modulo $n$,
and 0 otherwise. Works for any integer $a$.
:::

::: function
void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

Sets $a$ to $b^{-1}$ modulo $n$, for canonical $b$. Throws an exception
(`FLINT_IMPINV`) if $b$ is not invertible; test first with
`fmpz_mod_is_invertible` when this is possible.
:::

::: function
int fmpz_mod_divides(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

If the congruence $ac \equiv b \pmod n$ has a solution, sets $a$ to one
solution in $[0, n)$ and returns 1; otherwise returns 0 and leaves $a$
undefined. $b$ and $c$ must be canonical.

::: details
With $g = \gcd(c, n) = xc + yn$, a solution exists if and only if $g$
divides $b$, and then $a = (b/g) x \bmod n$ is one; it is unique only if
$g = 1$. For $c = 0$ a solution exists only for $b = 0$, and $a = 0$ is
returned.

:::
:::

## Powering {#pow}

::: function
void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong e, const fmpz_mod_ctx_t ctx)

Sets $a$ to $b^e$ modulo $n$, for canonical $b$. $0^0 = 1$.
:::

::: function
int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e, const fmpz_mod_ctx_t ctx)

Sets $a$ to $b^e$ modulo $n$ and returns 1, for canonical $b$ and any
integer $e$. If $e < 0$ and $b$ is not invertible, returns 0 and leaves
$a$ undefined.

::: details
Both powering functions use `fmpz_powm` (GMP's `mpz_powm`, Montgomery
multiplication with sliding windows). When `fft_small` is available and
$n$ has 70000 bits or more, they use generic binary or sliding-window
exponentiation with `fmpz_mod_mul` instead, whose products then use the
FFT. Exponents 0, 1 and 2 and bases 0 and 1 are handled directly. A
negative exponent inverts $b$ first.

:::
:::

## Random elements {#random}

::: function
void fmpz_mod_rand(fmpz_t a, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets $a$ to a uniformly random residue in $[0, n)$.
:::

::: function
void fmpz_mod_rand_not_zero(fmpz_t a, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets $a$ to a random nonzero residue: it draws up to three uniformly
random residues and returns the first nonzero one, or 1 if all three are
zero. (For $n = 1$ there is no nonzero residue, and the result is the
non-canonical value 1.)
:::

## Discrete logarithms via Pohlig–Hellman {#discrete-log}

These functions compute discrete logarithms modulo a prime $p$ for which
the prime factors of $p - 1$ are small. A precomputation for $p$ finds a
primitive root and tables for each prime factor of $p - 1$; each
logarithm then costs a few exponentiations and baby-step giant-step
searches. See [](fmpz_mod/discrete-log) and
[](algorithms/number-theory/discrete-log).

::: type
fmpz_mod_discrete_log_pohlig_hellman_struct
fmpz_mod_discrete_log_pohlig_hellman_t
fmpz_mod_discrete_log_pohlig_hellman_entry_struct
fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct

The precomputed data: a context for $p$, $p - 1$, the primitive root
$\alpha$ and its inverse, and for each prime power $q^e$ exactly dividing
$p - 1$ an entry with the idempotent for the Chinese remainder theorem,
the element $\gamma = \alpha^{(p-1)/q}$ of order $q$ and a sorted table
of powers of $\gamma$ for the baby-step giant-step search.
:::

::: function
void fmpz_mod_discrete_log_pohlig_hellman_init(fmpz_mod_discrete_log_pohlig_hellman_t L)
void fmpz_mod_discrete_log_pohlig_hellman_clear(fmpz_mod_discrete_log_pohlig_hellman_t L)

Initialise *L*, which is not ready for computation until
`fmpz_mod_discrete_log_pohlig_hellman_precompute_prime` is called, and
free the memory used by *L*.
:::

::: function
double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t p)

Configures *L* for logarithms modulo $p$ to the base of the least
primitive root modulo $p$, and returns an estimate of the number of
multiplications modulo $p$ that one call of
`fmpz_mod_discrete_log_pohlig_hellman_run` needs. $p$ must be prime,
which is not checked.

::: details
This factors $p - 1$ (with `fmpz_factor`) and throws an exception if a
prime factor does not fit in a word. The cost of the precomputation is
dominated by this factorisation, and the tables for each prime factor $q$
have at most 100 entries.

:::
:::

::: function
const fmpz * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(fmpz_mod_discrete_log_pohlig_hellman_t L)

Returns a pointer to the primitive root $\alpha$ used as the base of the
logarithms.
:::

::: function
void fmpz_mod_discrete_log_pohlig_hellman_run(fmpz_t x, const fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t y)

Sets $x$ to the logarithm of $y$ to the base $\alpha$, the integer
$x \in [0, p - 1)$ with $\alpha^x \equiv y \pmod p$. $y$ must be reduced
modulo $p$ and nonzero (not checked); if no logarithm is found, which
can only happen for invalid input, an exception is thrown.
:::

::: function
int fmpz_next_smooth_prime(fmpz_t a, const fmpz_t b)

Sets $a$ to a prime $p > b$ for which $p - 1$ is smooth and returns 1, or
sets $a$ to 0 and returns 0 if there is none. Such primes are suitable
for `fmpz_mod_discrete_log_pohlig_hellman_precompute_prime`.

::: details
The primes come from a fixed table of primes $p$ with
$p - 1 = 2^{e_1} 5^{e_2} 7^{e_3} 11^{e_4} 13^{e_5} 17^{e_6} 19^{e_7} 23^{e_8}$:
four to eleven primes in every interval $[2^{k-1}, 2^k)$ for
$16 \le k \le 64$ (the least is 33857), and twelve each in
$[2^{127}, 2^{128})$ and $[2^{191}, 2^{192})$. The function returns the
least table entry greater than $b$, which need not be the least prime
with this property (for $2^{64} \le b < 2^{127}$ it is a prime above
$2^{127}$), and returns 0 for $b$ beyond the last entry, about
$5.8 \cdot 10^{57}$. The exact set of primes should not be relied upon.

:::
:::

## Internal functions {#internal}

The arithmetic functions selected by `fmpz_mod_ctx_init`. They have the
interface of `fmpz_mod_add`, `fmpz_mod_sub` and `fmpz_mod_mul` and
require canonical inputs and a modulus of the stated size.

::: function
void _fmpz_mod_add1(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_sub1(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_mul1(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

For $n < 2^{\bits}$, with `nmod_add`, `nmod_sub` and `nmod_mul` and the
`nmod_t` stored in the context.
:::

::: function
void _fmpz_mod_add2s(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_sub2s(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_mul2s(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

For $n = 2^{\bits}$, with wrapping word arithmetic.
:::

::: function
void _fmpz_mod_add2(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_sub2(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_mul2(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

For $2^{\bits} < n < 2^{2\bits}$, with two- and three-limb arithmetic;
`_fmpz_mod_mul2` is the Barrett reduction of
[](algorithms/arithmetic/multiword-mulmod#two-limb).
:::

::: function
void _fmpz_mod_addN(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_subN(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_mulN(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

For any modulus, with `fmpz` arithmetic; `_fmpz_mod_mulN` reduces with
`fmpz_mod`, or with `fmpz_fdiv_r_preinvn` when the context has a
precomputed inverse.
:::
