---
title: qadic.h
header: qadic.h
description: Unramified extensions $\mathbf{Q}_q = \mathbf{Q}_p[X]/(f)$ of the $p$-adic numbers, with inversion, square roots, exponential, logarithm, the Frobenius automorphism, Teichmüller lifts, norm and trace.
category: p-adic numbers
order: 4
---

This module implements arithmetic in the unramified extension
$\mathbf{Q}_q$ of $\mathbf{Q}_p$ of degree $d$, $q = p^d$, represented
as $\mathbf{Q}_p[X]/(f)$ for a monic polynomial $f$ that is irreducible
modulo $p$. An element is a `qadic_t`, which is a `padic_poly_t` of
length at most $d$: a value $p^v g(X)$ with $g \in \Z[X]$, together with
an absolute precision $N$. As in `padic.h`, every function computes the
exact result from the stored inputs and reduces it modulo $p^N$ for the
precision $N$ of the output. The context `qadic_ctx_t` holds a
`padic_ctx_t` for $p$, the polynomial $f$ in a sparse representation
and the name of the generator.

Besides the field operations, the module provides the functions that
make unramified extensions useful in $p$-adic point counting: the
Frobenius automorphism, Teichmüller lifts, the norm and the trace, all
described in [](qadic/galois), and the exponential and logarithm.

::: summary algorithms/p-adic/unramified-extensions :::

## Types and context {#context}

::: type
qadic_struct
qadic_t

An element $p^v g(X) + O(p^N)$ of $\mathbf{Q}_q$; the same type as
`padic_poly_struct` and `padic_poly_t`. A *reduced* element has length
at most $d$, coefficients in $[0, p^{N - v})$ not all divisible by $p$
(or is zero with $v = 0$). All functions require reduced inputs and
produce reduced outputs unless stated otherwise.
:::

::: type
qadic_ctx_struct
qadic_ctx_t

A context for $\mathbf{Q}_q$. The fields are `pctx` (a
`padic_ctx_struct` for $p$, including the table of powers and the
printing mode), the arrays `a` and `j` of length `len` giving
$f(X) = \sum_{i < \mathtt{len}} a_i X^{j_i}$ with $j$ increasing and
the last term $X^d$, and `var`, the name of the generator used for
printing.
:::

::: function
void qadic_ctx_init(qadic_ctx_t ctx, const fmpz_t p, slong d, slong min, slong max, const char * var, enum padic_print_mode mode)

Initialises `ctx` for the extension of degree $d \ge 1$ of
$\mathbf{Q}_p$, with powers $p^e$ for $\mathtt{min} \le e < \mathtt{max}$
precomputed, generator name `var` and printing mode `mode`. The
defining polynomial is the Conway polynomial for $(p, d)$ if the
database contains it, otherwise for word-size $p$ the polynomial of
`nmod_poly_minimal_irreducible` (a binomial, trinomial, tetranomial or
pentanomial), and otherwise a random sparse irreducible polynomial. The
prime $p$ is not checked.
:::

::: function
void qadic_ctx_init_conway(qadic_ctx_t ctx, const fmpz_t p, slong d, slong min, slong max, const char * var, enum padic_print_mode mode)

Like `qadic_ctx_init` but always uses the Conway polynomial, and throws
if $p > 109987$ or the database has no Conway polynomial for
$(p, d)$.
:::

::: function
int _qadic_ctx_init_conway_ui(qadic_ctx_t ctx, ulong p, slong d, slong min, slong max, const char * var, enum padic_print_mode mode)

Initialises `ctx` with the Conway polynomial for $(p, d)$ and returns 1,
or returns 0 without initialising `ctx` if the database does not
contain it.
:::

::: function
void qadic_ctx_init_modulus(qadic_ctx_t ctx, const fmpz_t p, const fmpz_mod_poly_t modulus, slong min, slong max, const char * var, enum padic_print_mode mode)
void qadic_ctx_init_modulus_nmod(qadic_ctx_t ctx, ulong p, const nmod_poly_t modulus, slong min, slong max, const char * var, enum padic_print_mode mode)

Initialise `ctx` with the given defining polynomial, which must be
monic and irreducible over $\mathbf{F}_p$ (not checked). Its
coefficients, taken in $[0, p)$, define $f$; only the nonzero ones are
stored.
:::

::: function
void qadic_ctx_clear(qadic_ctx_t ctx)

Clears `ctx`.
:::

::: function
slong qadic_ctx_degree(const qadic_ctx_t ctx)

Returns the degree $d$ of the extension.
:::

::: function
void qadic_ctx_print(const qadic_ctx_t ctx)

Prints $p$, $d$ and $f$ to `stdout`.
:::

## Memory management and reduction {#memory}

::: function
void qadic_init(qadic_t rop)
void qadic_init2(qadic_t rop, slong prec)

Initialise `rop` to zero with precision `PADIC_DEFAULT_PREC`, respectively
`prec`.
:::

::: function
void qadic_clear(qadic_t rop)

Clears `rop`.
:::

::: function
void _fmpz_poly_reduce(fmpz * R, slong lenR, const fmpz * a, const slong * j, slong len)

Reduces `(R, lenR)` in place modulo the sparse monic polynomial
$f = \sum_k a_k X^{j_k}$ of degree $d = j_{\mathtt{len} - 1}$, setting
the coefficients of degree at least $d$ to zero. The array $j$ must be
increasing; `(R, lenR)` may have leading zeros.
:::

::: function
void _fmpz_mod_poly_reduce(fmpz * R, slong lenR, const fmpz * a, const slong * j, slong len, const fmpz_t p)

Like `_fmpz_poly_reduce`, followed by reduction of the remaining
coefficients modulo $p$ (typically a prime power) into $[0, p)$.
:::

::: function
void qadic_reduce(qadic_t rop, const qadic_ctx_t ctx)

Reduces `rop` modulo $f$ and $p^{N - v}$ and canonicalises it, or sets
it to zero if $v \ge N$.
:::

::: function
slong qadic_val(const qadic_t op)
slong qadic_prec(const qadic_t op)

Return the valuation $v$ and the precision $N$ of `op`.
:::

## Random generation {#random}

::: function
void qadic_randtest(qadic_t rop, flint_rand_t state, const qadic_ctx_t ctx)
void qadic_randtest_not_zero(qadic_t rop, flint_rand_t state, const qadic_ctx_t ctx)

Set `rop` to a random reduced element, respectively a random nonzero
one, with valuation drawn as for `padic_randtest`.
:::

::: function
void qadic_randtest_val(qadic_t rop, flint_rand_t state, slong v, const qadic_ctx_t ctx)

Sets `rop` to a random reduced element of valuation $v$, which is zero
if $v \ge N$.
:::

::: function
void qadic_randtest_int(qadic_t rop, flint_rand_t state, const qadic_ctx_t ctx)

Sets `rop` to a random reduced element of $\mathbf{Z}_q$, with valuation
in $[0, N)$, or to zero if $N \le 0$.
:::

## Assignment and conversion {#conversion}

::: function
void qadic_set(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Sets `rop` to `op`, reduced to the precision of `rop`.
:::

::: function
void qadic_zero(qadic_t rop)
void qadic_one(qadic_t rop)

Set `rop` to 0, respectively 1 (or 0 if the precision is not positive).
:::

::: function
void qadic_gen(qadic_t rop, const qadic_ctx_t ctx)

Sets `rop` to the generator $X$, or to zero if its precision is not
positive. For $d = 1$ it sets `rop` to the root $-a_0$ of $f$, reduced
modulo $p^N$.
:::

::: function
void qadic_set_ui(qadic_t rop, ulong op, const qadic_ctx_t ctx)

Sets `rop` to the integer `op`, reduced.
:::

::: function
void qadic_set_fmpz_poly(qadic_t rop, const fmpz_poly_t op, const qadic_ctx_t ctx)

Sets `rop` to the image of the integer polynomial `op` in
$\mathbf{Q}_q$, reduced modulo $f$ and $p^N$.
:::

::: function
int qadic_get_padic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)

If `op` lies in $\mathbf{Q}_p$, that is, has length at most 1, sets
`rop` to its value and returns 1; otherwise returns 0. The value is
canonicalised but not reduced to the precision of `rop`.
:::

## Comparison {#comparison}

::: function
int qadic_is_zero(const qadic_t op)
int qadic_is_one(const qadic_t op)

Return whether `op` is 0, respectively 1.
:::

::: function
int qadic_equal(const qadic_t op1, const qadic_t op2)

Returns whether the stored values are equal; the precisions are not
compared.
:::

## Arithmetic {#arithmetic}

All functions in this section allow aliasing between outputs and
inputs.

::: function
void qadic_add(qadic_t rop, const qadic_t op1, const qadic_t op2, const qadic_ctx_t ctx)
void qadic_sub(qadic_t rop, const qadic_t op1, const qadic_t op2, const qadic_ctx_t ctx)
void qadic_neg(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Set `rop` to $op_1 + op_2$, $op_1 - op_2$ and $-op$, reduced. These are
the corresponding `padic_poly` functions.
:::

::: function
void qadic_mul(qadic_t rop, const qadic_t op1, const qadic_t op2, const qadic_ctx_t ctx)

Sets `rop` to $op_1 \, op_2$, reduced: the integer polynomials are
multiplied in $\Z[X]$ and the product is reduced modulo $f$ and
$p^{N - v_1 - v_2}$ ([](algorithms/p-adic/unramified-extensions#mul)).
:::

::: function
void _qadic_inv(fmpz * rop, const fmpz * op, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N)

Sets `(rop, d)` to the inverse of the unit `(op, len)` modulo $f$ and
$p^N$, for $1 \le len \le d$ and $N \ge 1$, by Newton iteration from the
inverse modulo $p$. Does not allow aliasing.
:::

::: function
void qadic_inv(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Sets `rop` to $1/op$, reduced, or to zero if $N + v \le 0$. Throws if
`op` is zero.
:::

::: function
void _qadic_pow(fmpz * rop, const fmpz * op, slong len, const fmpz_t e, const fmpz * a, const slong * j, slong lena, const fmpz_t p)

Sets `(rop, 2d - 1)` to `(op, len)` raised to the power $e \ge 0$,
reduced modulo $f$ and the integer $p$ (typically a prime power); only
the first $d$ coefficients can be nonzero. Requires $1 \le len \le d$.
Does not allow aliasing.
:::

::: function
void qadic_pow(qadic_t rop, const qadic_t op, const fmpz_t e, const qadic_ctx_t ctx)

Sets `rop` to $op^e$, reduced, with $op^0 = 1$. Throws if $e < 0$.
Assumes that $e \cdot v$ fits in a `slong`.
:::

## Square roots {#sqrt}

::: function
int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Returns whether `op` is a square in $\mathbf{Q}_q$, and if so sets `rop`
to a square root, reduced. If the valuation of the root is at least the
precision of `rop`, `rop` is set to zero and squareness is decided from
the residue of the unit (modulo $p$, or modulo 8 for $p = 2$); if `op`
has no digits of precision left, 1 is returned.

::: details
For odd $p$ a square root in $\mathbf{F}_q$ is found by an exponentiation
or by Tonelli–Shanks and lifted by Newton iteration; for $p = 2$ an
Artin–Schreier equation over $\mathbf{F}_q$ is solved by linear algebra
over $\mathbf{F}_2$ first
([](algorithms/p-adic/unramified-extensions#sqrt)).
:::
:::

::: function
struct qadic2_sqrt_precomp * _qadic_char2_sqrt_precomp_init(const qadic_ctx_t ctx)
int _qadic_char2_sqrt_with_precomp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx, const struct qadic2_sqrt_precomp * data)
void _qadic_char2_sqrt_precomp_clear(struct qadic2_sqrt_precomp * data)

For $p = 2$, precompute the LUP decomposition of the matrix of the
Artin–Schreier map $Y \mapsto Y^2 + Y$ on $\mathbf{F}_q$, compute square
roots as `qadic_sqrt` using it, and free it. The structure is opaque.
The first two functions throw if $p \ne 2$.
:::

## Exponential and logarithm {#exp-log}

The algorithms are those of `padic.h`, with polynomial arithmetic modulo
$f$; see [](algorithms/p-adic/exp-log).

::: function
void _qadic_exp(fmpz * rop, const fmpz * op, slong v, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N, const fmpz_t pN)
void _qadic_exp_rectangular(fmpz * rop, const fmpz * op, slong v, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N, const fmpz_t pN)

Set `(rop, 2d - 1)` to the exponential of the nonzero element
`(op, v, len)` modulo $p^N$, where `pN` is $p^N$, assuming that the
series converges and $v < N$. Do not allow aliasing. `_qadic_exp` uses
rectangular splitting for $N < 8192/b$, where $b$ is the bit length of
$p$, and the balanced algorithm otherwise.
:::

::: function
void _qadic_exp_balanced(fmpz * rop, const fmpz * x, slong v, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N, const fmpz_t pN)

Sets `(rop, d)` to the exponential of `(x, v, len)` modulo $p^N$ by the
balanced (bit-burst) algorithm, assuming convergence and $v < N$.
Requires $1 \le len \le d$; the input may be zero. Allows aliasing.
:::

::: function
int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
int qadic_exp_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
int qadic_exp_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Return whether the exponential series converges at `op`, that is,
whether `op` is zero or has valuation at least 1 (odd $p$) or 2
($p = 2$), and if so set `rop` to $\exp(op)$, reduced.
:::

::: function
void _qadic_log(fmpz * z, const fmpz * y, slong v, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N, const fmpz_t pN)
void _qadic_log_rectangular(fmpz * z, const fmpz * y, slong v, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N, const fmpz_t pN)
void _qadic_log_balanced(fmpz * z, const fmpz * y, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N, const fmpz_t pN)

Set `(z, d)` to

$$
\log(1 - y) = -\sum_{i = 1}^{\infty} \frac{y^i}{i} \bmod p^N
$$

for a nonzero $y$ reduced modulo $p^N$, of valuation $v$ with
$1 \le v < N$. The array $z$ must have space for
$2d - 1$ coefficients. `_qadic_log_rectangular` allows aliasing between
$y$ and $z$, `_qadic_log_balanced` does not. `_qadic_log` uses
rectangular splitting for $N < 1024/b$, where $b$ is the bit length of
$p$, and the balanced algorithm otherwise.
:::

::: function
int qadic_log(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
int qadic_log_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
int qadic_log_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

If $\operatorname{ord}_p(op - 1) \ge 1$ (for every $p$, including
$p = 2$), set `rop` to $\log(op)$, reduced, and return 1; otherwise,
in particular for $op = 0$, return 0.
:::

## Frobenius, Teichmüller lifts, norm and trace {#galois}

These functions are described in more detail on [](qadic/galois).

::: function
void _qadic_frobenius_a(fmpz * rop, slong e, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N)

Sets `(rop, 2d - 1)` to $\Sigma^e(X)$ modulo $p^N$, where $\Sigma$ is
the Frobenius automorphism; only the first $d$ coefficients can be
nonzero. Requires $d \ge 2$, $N \ge 2$ and $0 < e < d$.
:::

::: function
void _qadic_frobenius(fmpz * rop, const fmpz * op, slong len, slong e, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N)

Sets `(rop, 2d - 1)` to $\Sigma^e$ applied to `(op, len)`, modulo $p^N$,
for $1 \le len \le d$ and $0 < e < d$. Does not allow aliasing.
:::

::: function
void qadic_frobenius(qadic_t rop, const qadic_t op, slong e, const qadic_ctx_t ctx)

Sets `rop` to $\Sigma^e(op)$, reduced, for any integer $e$, which is
taken modulo $d$. $\Sigma$ is the automorphism of $\mathbf{Q}_q$ that
reduces to $x \mapsto x^p$ on $\mathbf{F}_q$; it generates the Galois
group, which is cyclic of order $d$.
:::

::: function
void _qadic_teichmuller(fmpz * rop, const fmpz * op, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N)

Sets `(rop, d)` to the Teichmüller lift of the unit `(op, len)` modulo
$p^N$, for $N \ge 1$. Allows aliasing.
:::

::: function
void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Sets `rop` to the Teichmüller lift of `op`: the unique $(q - 1)$-th
root of unity congruent to `op` modulo $p$, reduced. Sets `rop` to zero
if `op` is zero or divisible by $p$, or if the precision of `rop` is not
positive. Throws if `op` has negative valuation.
:::

::: function
void _qadic_trace(fmpz_t rop, const fmpz * op, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t pN)
void qadic_trace(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Set `rop` to the trace of `op` over $\mathbf{Q}_p$, reduced modulo
`pN`, respectively to the precision of `rop`. The trace is computed from
the power sums of the roots of $f$ by Newton's identities.
:::

::: function
void _qadic_norm(fmpz_t rop, const fmpz * op, slong len, const fmpz * a, const slong * j, slong lena, const fmpz_t p, slong N)

Sets `rop` to the norm of the unit `(op, len)`, $len \ge 1$, reduced
modulo $p^N$; the result is a unit. Chooses between the analytic and the
resultant method.
:::

::: function
void qadic_norm(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
void qadic_norm_analytic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)

Set `rop` to the norm of `op` over $\mathbf{Q}_p$, reduced to the
precision of `rop`. The *analytic* version computes
$\exp(\operatorname{Tr}(\log u))$ for the unit part $u$ and throws
unless $u \equiv 1 \pmod p$ ($\bmod 4$ for $p = 2$) or $u \in
\mathbf{Q}_p$; the *resultant* version computes
$\operatorname{Res}(f, u)$ as a determinant modulo $p^N$; `qadic_norm`
chooses.
:::

## Output {#io}

::: function
int qadic_fprint_pretty(FILE * file, const qadic_t op, const qadic_ctx_t ctx)
int qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx)

Print `op` as a polynomial in the generator, in the printing mode of the
context: in `PADIC_TERSE` mode as a polynomial with rational
coefficients, in `PADIC_SERIES` mode as a sum of polynomials with
coefficients in $[0, p)$ times powers of $p$ (throws if a coefficient
is negative), and in `PADIC_VAL_UNIT` mode as `(g)*p^v`. Return 1.
:::

::: function
int qadic_debug(const qadic_t op)

Prints the raw fields of `op`, as `padic_poly_debug`. Returns 1.
:::
