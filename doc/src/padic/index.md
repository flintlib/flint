---
title: padic.h
header: padic.h
description: Elements of the $p$-adic field $\mathbf{Q}_p$ as $p^v u$ reduced modulo $p^N$, with inversion, square roots, exponential, logarithm and Teichmüller lifts.
category: p-adic numbers
order: 1
---

This module implements arithmetic in the field $\mathbf{Q}_p$ of $p$-adic
numbers, for a prime $p$ of any size. An element is a `padic_t`, which
stores a value $x = p^v u$ with $u, v \in \Z$ together with an absolute
precision $N$; the value represents the residue class of $x$ modulo $p^N$.
The prime, a table of powers of $p$ and the printing mode live in a
separate context `padic_ctx_t` that is passed to every function. Each
function computes its result exactly from the values of its inputs and
then reduces it modulo $p^N$ for the precision $N$ of the output, much as
MPFR rounds to the precision of the output variable; the precision of the
inputs is not used, and no error bound is propagated.

The guide [](padic/precision) explains this model with an example.
Polynomials and matrices over $\mathbf{Q}_p$ are in
`padic_poly.h` and `padic_mat.h`, unramified extensions in `qadic.h`.
The newer module `padic_radix.h` implements $\mathbf{Q}_p$ for word-size
$p$ with tracked error bounds and the generic ring interface.

::: summary algorithms/p-adic/overview :::

## Representation {#representation}

::: type
padic_struct
padic_t

A $p$-adic number $p^v u + O(p^N)$. The structure has the fields `u` (an
`fmpz`), `v` and `N` (both `slong`). A `padic_t` is an array of one
`padic_struct`, so it is passed by reference.
:::

The value is in *canonical form* if $p \nmid u$, or if $u = v = 0$. It is
*reduced* if moreover $0 < u < p^{N - v}$ when $u \ne 0$, so that a
reduced element with $v \ge N$ is zero. All functions require canonical
inputs and produce reduced outputs, except where the description says
otherwise. The precision $N$ belongs to the variable: it is set when the
variable is initialised, may be changed through `padic_prec`, and is
never changed by arithmetic (`padic_swap` swaps it along with the value).
$N$ may be zero or negative.

::: summary algorithms/p-adic/precision :::

::: function
fmpz * padic_unit(const padic_t op)

Returns a pointer to the unit part $u$, which can be passed to the
functions of `fmpz.h`.
:::

::: macro
padic_val(x)
padic_prec(x)

The valuation $v$ and the precision $N$ of `x`, as lvalues of type
`slong`. Lowering the precision of a reduced element does not reduce it;
call `padic_reduce` afterwards.
:::

::: function
slong padic_get_val(const padic_t op)
slong padic_get_prec(const padic_t op)

Return the valuation $v$ and the precision $N$ of `op`.
:::

::: macro
PADIC_DEFAULT_PREC

The precision given to a variable by `padic_init`, equal to 20.
:::

## Context {#context}

::: type
padic_ctx_struct
padic_ctx_t

A context for $\mathbf{Q}_p$. The fields are `p` (an `fmpz_t`), `pinv`
(a `double` approximation of $1/p$, used when $p$ fits in a word), the
table `pow` with `pow[i]` $= p^{\mathtt{min} + i}$ for
$0 \le i < \mathtt{max} - \mathtt{min}$, the bounds `min` and `max`, and
`mode`, the printing mode, which may be changed at any time.
:::

::: type
enum padic_print_mode

The printing modes `PADIC_TERSE`, `PADIC_SERIES` and `PADIC_VAL_UNIT`.
For $x = 12 \cdot 7^{-1} \in \mathbf{Q}_7$ they print `12/7`,
`5*7^-1 + 1` and `12*7^-1`: as a rational number, as the sum of its
$p$-adic digits, and as unit times power of $p$.
:::

::: function
void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, slong min, slong max, enum padic_print_mode mode)

Initialises `ctx` for the prime $p$ with precomputed powers $p^e$ for
$\mathtt{min} \le e < \mathtt{max}$ and the given printing mode. Throws
unless $0 \le \mathtt{min} \le \mathtt{max}$. That $p$ is prime is not
checked; for composite $p$ the behaviour is undefined.

::: details
The table is used by every reduction modulo $p^{N - v}$, so it should
cover the differences $N - v$ that occur, typically from 0 to a little
more than the largest precision in use. Other powers are computed when
needed. With $\mathtt{min} = \mathtt{max}$ no table is stored.
:::
:::

::: function
void padic_ctx_clear(padic_ctx_t ctx)

Clears the context.
:::

::: function
void padic_ctx_pow_ui(fmpz_t rop, ulong e, const padic_ctx_t ctx)

Sets `rop` to $p^e$, copying it from the table when possible. Throws if
$e$ exceeds `WORD_MAX`.
:::

::: function
int _padic_ctx_pow_ui(fmpz_t rop, ulong e, const padic_ctx_t ctx)

Makes the *uninitialised* `rop` equal to $p^e$. If $p^e$ is in the
table, `rop` becomes a shallow copy of the table entry, which must be
neither modified nor cleared, and the function returns 0. Otherwise it
initialises `rop`, computes the power and returns 1, and the caller must
clear `rop`. Throws if $e$ exceeds `WORD_MAX`.
:::

## Memory management {#memory}

::: function
void padic_init(padic_t rop)
void padic_init2(padic_t rop, slong N)

Initialise `rop` to zero with precision `PADIC_DEFAULT_PREC` (that is,
20), respectively $N$.
:::

::: function
void padic_clear(padic_t rop)

Clears `rop`.
:::

::: function
void _padic_canonicalise(padic_t rop, const padic_ctx_t ctx)

Brings `rop` into canonical form by moving all factors $p$ of $u$ into
$v$, and sets $v = 0$ if $u = 0$. Does not reduce.
:::

::: function
void _padic_reduce(padic_t rop, const padic_ctx_t ctx)

Reduces the canonical `rop` modulo $p^N$: replaces $u$ by
$u \bmod p^{N - v}$ if $v < N$, and sets `rop` to zero otherwise.
:::

::: function
void padic_reduce(padic_t rop, const padic_ctx_t ctx)

Canonicalises and reduces `rop`.
:::

## Random generation {#random}

::: function
void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)

Sets `rop` to a random reduced element. Its valuation, when it is
nonzero, lies in $[-\lceil N/10 \rceil, N)$ for $N > 0$, in
$[N - \lceil -N/10 \rceil, N)$ for $N < 0$ and in $[-10, 0)$ for
$N = 0$.
:::

::: function
void padic_randtest_not_zero(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)

Like `padic_randtest`, but nonzero: after ten zero draws it sets `rop`
to $p^{N - 1}$.
:::

::: function
void padic_randtest_int(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)

Sets `rop` to a random reduced $p$-adic integer, with valuation in
$[0, N)$, or to zero if $N \le 0$.
:::

## Assignment and conversion {#conversion}

The functions that set a `padic_t` set it to the value of the input,
reduced modulo $p^N$ for the precision $N$ of the output.

::: function
void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Sets `rop` to `op`. The precision of `rop` is kept.
:::

::: function
void padic_set_si(padic_t rop, slong op, const padic_ctx_t ctx)
void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx)
void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)

Set `rop` to the integer `op`.
:::

::: function
void padic_set_fmpq(padic_t rop, const fmpq_t op, const padic_ctx_t ctx)

Sets `rop` to the canonical rational number `op`. The power of $p$ in
the denominator becomes a negative valuation, and the rest of the
denominator is inverted modulo $p^{N - v}$.
:::

::: function
void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)
void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx)

Set `rop` to the GMP integer or rational `op`. These are declared only
if `gmp.h` is included before `padic.h`.
:::

::: function
void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)

Sets `rop` to the integer $p^v u$ represented by `op`. Throws if
$v < 0$.
:::

::: function
void padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx)

Sets `rop` to the rational number $p^v u$ represented by `op`.
:::

::: function
void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx)
void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx)

GMP versions of `padic_get_fmpz` (which throws if $v < 0$) and
`padic_get_fmpq`. Declared only if `gmp.h` is included before
`padic.h`.
:::

::: function
void padic_swap(padic_t op1, padic_t op2)

Swaps `op1` and `op2`, including their precisions. This differs from an
exchange through `padic_set` when the precisions differ.
:::

::: function
void padic_zero(padic_t rop)

Sets `rop` to zero.
:::

::: function
void padic_one(padic_t rop)

Sets `rop` to one, or to zero if its precision is not positive.
:::

## Comparison {#comparison}

::: function
int padic_is_zero(const padic_t op)

Returns whether `op` is zero, that is, whether $u = 0$.
:::

::: function
int padic_is_one(const padic_t op)

Returns whether $u = 1$ and $v = 0$.
:::

::: function
int padic_equal(const padic_t op1, const padic_t op2)

Returns whether the stored values are equal, that is, whether
$u_1 = u_2$ and $v_1 = v_2$. The precisions are not compared: elements
of different precisions that agree modulo the smaller power of $p$ are
equal only if their stored values coincide.
:::

## Arithmetic {#arithmetic}

All functions in this section allow aliasing between the output and the
inputs.

::: function
void padic_add(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
void padic_sub(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)

Set `rop` to $op_1 + op_2$ and $op_1 - op_2$, reduced.
:::

::: function
void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Sets `rop` to $-op$, reduced.
:::

::: function
void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)

Sets `rop` to $op_1 \, op_2$, reduced: the units are multiplied and the
product is reduced modulo $p^{N - v_1 - v_2}$.
:::

::: function
void padic_shift(padic_t rop, const padic_t op, slong v, const padic_ctx_t ctx)

Sets `rop` to $op \cdot p^v$, reduced. Here $v$ may be negative.
:::

::: function
void padic_div(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)

Sets `rop` to $op_1 / op_2$, reduced. Throws if $op_2$ is zero.
:::

::: function
void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Sets `rop` to $1/op$, reduced. For $op = p^v u$ this is $p^{-v}$ times
the inverse of $u$ modulo $p^{N + v}$, computed by Hensel lifting; if
$N + v \le 0$ the result is zero. Throws if `op` is zero.

::: details
The inverse modulo $p$ is computed by the extended Euclidean algorithm
and lifted by the Newton iteration $z \gets 2z - u z^2$, which doubles
the number of correct digits in each step
([](algorithms/p-adic/hensel-lifting)). For $p = 2$ the iteration runs
directly on the limbs (`fmpz_invmod_2exp`).
:::
:::

::: function
int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Returns whether `op` is a square in $\mathbf{Q}_p$ and in that case sets
`rop` to a square root, reduced. If the return value is 0, the value of
`rop` is undefined.

::: details
A nonzero $x = p^v u$ is a square if and only if $v$ is even and the
unit $u$ is a square in $\mathbf{Z}_p$; for odd $p$ this holds if and
only if $u$ is a square modulo $p$, and for $p = 2$ if and only if
$u \equiv 1 \pmod 8$. The root is $p^{v/2}$ times a square root of $u$
modulo $p^{N - v/2}$, obtained from a root modulo $p$ (or the
congruence modulo 8) by Hensel lifting of the inverse square root
([](algorithms/p-adic/hensel-lifting)). The decision is made from the
exact value of `op` even when the precision of `rop` is too small to
hold a nonzero root.
:::
:::

::: function
void padic_pow_si(padic_t rop, const padic_t op, slong e, const padic_ctx_t ctx)

Sets `rop` to $op^e$, reduced, for any integer $e$, with $op^0 = 1$.
For $e < 0$ the unit is first inverted modulo $p^{N - ev}$. If `op` is
zero and $e \ne 0$ the result is zero, also for $e < 0$. Assumes that
$ev$ does not overflow a `slong`.

::: details
If the input $p^v u$ is only known modulo $p^N$, then $x^e$ is only
known modulo $p^{N + (e - 1)v}$, which is a loss of precision when
$v < 0$ ([](algorithms/p-adic/precision)). The function computes the
power of the stored value exactly and does not account for this.
:::
:::

## Exponential and logarithm {#exp-log}

[](padic/exp-log) describes the functions of this section and
the choice of algorithm in more detail.

::: function
slong _padic_exp_bound(slong v, slong N, const fmpz_t p)

Returns $n$ such that $\operatorname{ord}_p(x^j / j!) \ge N$ for all
$j \ge n$ and all $x$ with $\operatorname{ord}_p(x) = v$, so that the
exponential series can be truncated after $n$ terms. Requires
$1 \le v < N$ for odd $p$ and $2 \le v < N$ for $p = 2$.

::: details
If $p$ fits in a `slong` the result is
$\lceil ((p - 1)N - 1) / ((p - 1)v - 1) \rceil$, which follows from
$\operatorname{ord}_p(j!) \le (j - 1)/(p - 1)$; otherwise it is
$\lceil N / v \rceil$, since then $j!$ has no factor $p$ for the $j$
that occur.
:::
:::

::: function
void _padic_exp(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)
void _padic_exp_rectangular(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)
void _padic_exp_balanced(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)

Set `rop` to $\exp(p^v u)$ reduced modulo $p^N$. Require
$p^v u \ne 0$, $v < N$, and $v \ge 1$ for odd $p$, $v \ge 2$ for
$p = 2$. Allow aliasing between `rop` and $u$. `_padic_exp` uses the
rectangular algorithm for $N < 1024$ and the balanced one otherwise.
:::

::: function
int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
int padic_exp_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
int padic_exp_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Return whether the exponential series converges at `op`, that is,
whether `op` is zero or has valuation at least 1 (odd $p$) or 2
($p = 2$), and if so set `rop` to $\exp(op)$, reduced. Otherwise `rop`
is not changed. The *rectangular* version evaluates the series by
rectangular splitting, the *balanced* version by the $p$-adic analogue
of the bit-burst method with binary splitting, which is quasi-linear in
$N$ for fixed $p$.
:::

::: function
slong _padic_log_bound(slong v, slong N, const fmpz_t p)

Returns $b$ such that $iv - \operatorname{ord}_p(i) \ge N$ for all
$i \ge b$, so that the logarithm series for an argument of valuation $v$
can be truncated after $b - 1$ terms. Requires $1 \le v < N$ and
$N < 2^{\bits - 2}$; throws if $N$ is larger.
:::

::: function
void _padic_log(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
void _padic_log_rectangular(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
void _padic_log_satoh(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
void _padic_log_balanced(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)

Set $z$ to

$$
\log(1 - y) = -\sum_{i = 1}^{\infty} \frac{y^i}{i} \bmod p^N
$$

for an integer $y \ne 0$ with $v = \operatorname{ord}_p(y)$ and
$1 \le v < N$; the logarithm of $x$ is obtained with $y = 1 - x$. Do not
allow aliasing between $y$ and $z$. `_padic_log` uses the rectangular
version for $N < 512 / b$, where $b$ is the bit length of $p$, and the
balanced version otherwise. The balanced version does not use $v$.
:::

::: function
int padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx)
int padic_log_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
int padic_log_satoh(padic_t rop, const padic_t op, const padic_ctx_t ctx)
int padic_log_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Return whether the logarithm series converges at `op`, that is, whether
$\operatorname{ord}_p(op - 1) \ge 1$, and if so set `rop` to
$\log(op)$, reduced. Otherwise `rop` is not changed. The condition is
the same for every $p$, including $p = 2$. The *rectangular* version
evaluates the series by rectangular splitting, the *Satoh* version first
raises `op` to the power $p^k$ with $k = \lfloor \sqrt N \rfloor$, and
the *balanced* version writes `op` as a product of factors whose
logarithms are summed by binary splitting.
:::

## Teichmüller lifts and factorials {#special}

::: function
void _padic_teichmuller(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N)

Sets `rop` to the Teichmüller lift of the $p$-adic unit `op` modulo
$p^N$: the unique root of unity of order dividing $p - 1$ that is
congruent to `op` modulo $p$. Requires $N \ge 1$ and allows aliasing.
For $p = 2$ the result is 1.
:::

::: function
void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)

Sets `rop` to the Teichmüller lift of `op`, reduced. If `op` is zero or
divisible by $p$ the result is zero (which satisfies $t^p = t$ but is
not a root of unity). Throws if `op` has negative valuation.

::: details
The lift is the root of $t^p - t$ congruent to `op` modulo $p$, computed
by Newton iteration. Since the derivative at the root is $p - 1$, the
iteration $t \gets t - (t^p - t)/(p - 1)$, with an inverse of $p - 1$
lifted alongside, suffices ([](algorithms/p-adic/hensel-lifting)).
:::
:::

::: function
ulong padic_val_fac_ui_2(ulong n)
ulong padic_val_fac_ui(ulong n, const fmpz_t p)

Return $\operatorname{ord}_2(n!)$ and $\operatorname{ord}_p(n!)$, by
Legendre's formula $\sum_{k \ge 1} \lfloor n / p^k \rfloor$. The result
fits in a `ulong` since $\operatorname{ord}_p(n!) \le (n - 1)/(p - 1)$.
`padic_val_fac_ui` returns 0 if $p$ does not fit in a `ulong`.
:::

::: function
void padic_val_fac(fmpz_t rop, const fmpz_t op, const fmpz_t p)

Sets `rop` to $\operatorname{ord}_p(op!)$. Throws if `op` is negative.
:::

## Input and output {#io}

::: function
char * _padic_get_str(char * str, const padic_t op, const fmpz_t p, enum padic_print_mode mode)
char * padic_get_str(char * str, const padic_t op, const padic_ctx_t ctx)

Return the string representation of `op` in the given printing mode,
respectively the mode of the context. If `str` is `NULL`, a new string
is allocated with `flint_malloc`; otherwise `str` must be large enough
and is returned.
:::

::: function
int _padic_fprint(FILE * file, const fmpz_t u, slong v, const padic_ctx_t ctx)
int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx)
int _padic_print(const fmpz_t u, slong v, const padic_ctx_t ctx)
int padic_print(const padic_t op, const padic_ctx_t ctx)

Print $p^v u$, respectively `op`, to `file` or to `stdout` in the
printing mode of the context, and return 1.
:::

::: function
void padic_debug(const padic_t op)

Prints the raw fields of `op` to `stdout` in the form `(u v N)`.
:::

## Low-level functions {#low-level}

These functions work on integers rather than on `padic_t`. The lifting
functions implement the Newton iteration of
[](algorithms/p-adic/hensel-lifting), which runs through the precisions
$a_0 = N$, $a_{i+1} = \lceil a_i / 2 \rceil$ down to 1 and back up.

::: type
padic_inv_struct
padic_inv_t

Precomputed data for repeated inversion modulo $p^N$: the length `n` of
the precision sequence and an array `pow` holding the powers
$p^{a_i}$ together with scratch space.
:::

::: function
void _padic_inv_precompute(padic_inv_t S, const fmpz_t p, slong N)
void _padic_inv_clear(padic_inv_t S)

Initialise `S` for inversion modulo $p^N$, $N \ge 1$, and clear it.
:::

::: function
void _padic_inv_precomp(fmpz_t rop, const fmpz_t op, const padic_inv_t S)

Sets `rop` to the inverse of the $p$-adic unit `op` modulo $p^N$, using
the data `S`. The scratch space in `S` is overwritten, so `S` must not
be shared between threads. Allows aliasing.
:::

::: function
void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N)

Sets `rop` to the inverse of the $p$-adic unit `op` modulo $p^N$, for
$N \ge 1$. The result is in $[0, p^N)$. Allows aliasing.
:::

::: function
int _padic_sqrt(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N)

Returns whether the $p$-adic unit `op` is a square in $\mathbf{Z}_p$
and if so sets `rop` to a square root modulo $p^N$, for $N \ge 1$. For
odd $p$ this is decided modulo $p$, for $p = 2$ modulo 8. Allows
aliasing.
:::

::: function
slong * _padic_lifts_exps(slong * n, slong N)

Returns the sequence $a_0 = N$, $a_{i+1} = \lceil a_i/2 \rceil$, ending
with $a_{n-1} = 1$, and sets `*n` to its length
$\lceil \log_2 N \rceil + 1$. Requires $N \ge 1$. The array is
allocated with `flint_malloc` and must be freed with `flint_free`.
:::

::: function
void _padic_lifts_pows(fmpz * pow, const slong * a, slong n, const fmpz_t p)

Sets `pow[i]` to $p^{a_i}$ for a sequence $a$ of length $n$ as computed
by `_padic_lifts_exps`, using about one multiplication per entry.
:::
