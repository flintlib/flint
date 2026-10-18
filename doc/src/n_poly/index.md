---
title: n_poly.h
header: n_poly.h
description: Lightweight dense and sparse polynomial containers over Z/pZ and F_q with word-size coefficients, used internally by the multivariate polynomial code.
category: Univariate polynomials
order: 7
---

This module provides light containers for polynomials with word-size
coefficients, used internally by the multivariate polynomial modules
(`nmod_mpoly`, `fq_nmod_mpoly`, their factorisation modules and the
`fmpz_mpoly` modular algorithms). They differ from `nmod_poly_t` and
`fq_nmod_poly_t` in three ways:

- An `n_poly_t` is just an array of `ulong` with a length and an
  allocation. It does not store a modulus: functions that need one take an
  `nmod_t` (the `n_poly_mod_` functions) or an `fq_nmod_ctx_t` (the
  `n_fq_` functions) as last argument, so large arrays of polynomials cost
  no more than the coefficients.
- An element of $\mathbb{F}_q = \mathbb{F}_p[t]/(m(t))$, $d = \deg m$, is
  stored as $d$ consecutive words, the coefficients of its representative
  of degree less than $d$ (an *n_fq* element), and a polynomial over
  $\mathbb{F}_q$ (an `n_fq_poly_t`, the same structure as `n_poly_t`)
  stores $d \cdot \mathit{length}$ words. This avoids the separate
  allocation of every coefficient in `fq_nmod_poly_t`.
- There are bivariate (`n_bpoly_t`), trivariate (`n_tpoly_t`) and sparse
  (`n_polyu_t`, `n_polyun_t`) variants, and stacks of preallocated
  temporaries, which the multivariate GCD and factorisation code uses in
  its inner loops.

The functions do little argument checking (most checks are assertions
active only in debug builds), and their interfaces are designed for
internal use and may change. Most arithmetic is delegated to the
underscore functions of `nmod_poly` and `nmod_vec`. Application code
should normally use `nmod_poly`, `fq_nmod_poly`, `nmod_mpoly` or the
generic `gr_poly` instead.

## Types and conventions {#types}

::: type
n_poly_struct
n_poly_t

A dense polynomial with word-size coefficients: fields `coeffs` (an array
of `alloc` words), `alloc` and `length`. The coefficient of $x^i$ is
`coeffs[i]`.
:::

::: type
n_fq_poly_struct
n_fq_poly_t

A dense polynomial over $\mathbb{F}_q$, $q = p^d$: the same structure as
`n_poly_struct`, where `length` counts coefficients in $\mathbb{F}_q$,
coefficient $i$ occupies the words `coeffs[d*i]`, ..., `coeffs[d*i + d -
1]`, and `alloc` counts words.
:::

::: type
n_bpoly_struct
n_bpoly_t
n_fq_bpoly_struct
n_fq_bpoly_t

A dense bivariate polynomial $\sum_i A_i(y) x^i$, stored as an array
`coeffs` of `alloc` polynomials (`n_poly_struct`, or `n_fq_poly` over
$\mathbb{F}_q$) with `length` used entries. In function names, `gen0` is
the outer variable $x$ (the index into `coeffs`) and `gen1` the inner
variable $y$.
:::

::: type
n_tpoly_struct
n_tpoly_t
n_fq_tpoly_struct
n_fq_tpoly_t

A dense trivariate polynomial, stored as an array of `n_bpoly_struct`.
:::

::: type
n_polyu_struct
n_polyu_t
n_fq_polyu_struct
n_fq_polyu_t

A sparse polynomial with word-size coefficients: arrays `exps` and
`coeffs` of `alloc` entries, of which `length` are used. Several
exponents are packed into one word; the `n_polyu3` functions use three
fields of `FLINT_BITS/3` bits, the first variable in the most significant
field.
:::

::: type
n_polyun_struct
n_polyun_t
n_fq_polyun_struct
n_fq_polyun_t

A polynomial that is sparse in some variables and dense in the last one:
arrays `exps` (packed exponents, as for `n_polyu_t`) and `coeffs` (dense
`n_poly_struct` coefficients in the last variable). In the names
`n_polyu1n`, `n_polyu2n` and `n_polyu3n`, the digit is the number of
exponents packed into `exps`.
:::

::: type
n_poly_stack_struct
n_poly_stack_t
n_bpoly_stack_struct
n_bpoly_stack_t
n_polyun_stack_struct
n_polyun_stack_t
n_poly_bpoly_stack_struct
n_poly_bpoly_stack_t
n_poly_polyun_stack_struct
n_poly_polyun_stack_t

Stacks of reusable temporaries: an array of pointers to initialised
objects (`array`, `alloc`) and the index `top` of the first free one.
Objects taken from a stack keep their memory when given back, so that
loops do not allocate. The combined types hold one stack of each kind.
:::

::: type
nmod_eval_interp_struct
nmod_eval_interp_t

Precomputed data for converting between the coefficients of polynomials
of bounded degree and their values at a fixed set of points; see
[](#evaluation).
:::

- Functions taking an `nmod_t` assume coefficients reduced modulo $n$
  and produce reduced outputs; those taking an `fq_nmod_ctx_t` assume
  canonical $\mathbb{F}_q$ elements. Functions that divide need a prime
  modulus.
- A polynomial is normalised if its leading coefficient (element) is
  nonzero. The functions without a leading underscore normalise their
  output unless stated otherwise.
- Functions with a trailing underscore, such as `n_fq_poly_mul_`, take a
  stack of temporaries as an additional argument.
- Unless stated otherwise, outputs may alias inputs.


## Memory management and basic operations {#basic}

::: function
void n_poly_init(n_poly_t A)
void n_poly_init2(n_poly_t A, slong alloc)
void n_poly_clear(n_poly_t A)

Initialise `A` as the zero polynomial (with space for `alloc`
coefficients for `n_poly_init2`), respectively free it.
:::

::: function
void n_poly_realloc(n_poly_t A, slong len)
void n_poly_fit_length(n_poly_t A, slong len)

Ensure that `A` has space for at least `len` coefficients. The
allocation grows to $\max(\mathit{len}, \tfrac32 \mathit{alloc} + 1)$
and never shrinks.
:::

::: function
void n_poly_swap(n_poly_t A, n_poly_t B)

Swaps `A` and `B` efficiently.
:::

::: function
int n_poly_is_canonical(const n_poly_t A)

Returns whether `A` is normalised (its length is nonnegative and its
leading coefficient nonzero). It does not check that the coefficients are
reduced; see `n_poly_mod_is_canonical`.
:::

::: function
void _n_poly_normalise(n_poly_t A)
void n_poly_truncate(n_poly_t poly, slong len)

Normalise `A` by decreasing the length past zero leading
coefficients, respectively truncate `poly` to length `len` and
normalise.
:::

::: function
slong n_poly_degree(const n_poly_t A)
ulong n_poly_lead(const n_poly_t A)

Return the degree ($-1$ for zero) and the leading coefficient (the
polynomial must be nonzero).
:::

::: function
void n_poly_zero(n_poly_t res)
void n_poly_one(n_poly_t A)
void n_poly_set_ui(n_poly_t A, ulong c)
void n_poly_set(n_poly_t A, const n_poly_t B)

Set the polynomial to 0, 1, the constant $c$ (which must be reduced),
or a copy of `B`.
:::

::: function
int n_poly_is_zero(const n_poly_t poly)
int n_poly_is_one(const n_poly_t A)
int n_poly_equal(const n_poly_t a, const n_poly_t b)

Return whether the polynomial is 0, 1, or equal to the other one.
:::

::: function
ulong n_poly_get_coeff(const n_poly_t poly, slong j)
void n_poly_set_coeff(n_poly_t A, slong e, ulong c)
void n_poly_set_coeff_nonzero(n_poly_t A, slong j, ulong c)

Get or set the coefficient of $x^j$ (or $x^e$). `n_poly_set_coeff`
normalises the result; `n_poly_set_coeff_nonzero` requires $c \ne 0$.
The value is not reduced.
:::

::: function
void n_poly_shift_left(n_poly_t A, const n_poly_t B, slong k)
void n_poly_shift_right(n_poly_t res, const n_poly_t poly, slong k)

Multiply by $x^k$, respectively divide by $x^k$ discarding the low
coefficients.
:::

::: function
void n_poly_reverse(n_poly_t output, const n_poly_t input, slong m)

Sets `output` to the reversal of `input` viewed as a polynomial of
length $m$, normalised.
:::

::: function
void n_poly_print_pretty(const n_poly_t A, const char * x)

Prints `A` to standard output with variable name `x`.
:::

::: function
void n_poly_set_nmod_poly(n_poly_t a, const nmod_poly_t b)
void nmod_poly_set_n_poly(nmod_poly_t a, const n_poly_t b)

Copy the coefficients between an `n_poly_t` and an `nmod_poly_t`
(whose modulus is left unchanged).
:::

::: function
void n_poly_mock(n_poly_t a, const nmod_poly_t b)
void nmod_poly_mock(nmod_poly_t a, const n_poly_t b, nmod_t mod)

Make `a` a shallow view of `b`, sharing its coefficient array, so
that functions of the other type can read it without copying
(`nmod_poly_mock` also sets the modulus of the view). The view must not
be cleared or used to modify the length or allocation.
:::

## Arithmetic modulo a word-size prime {#mod}

The functions `n_poly_mod_*` treat the coefficients modulo $n$ given by
an `nmod_t` and mostly wrap the corresponding `nmod_poly` functions. The
versions with a leading underscore require the output to be distinct from
the inputs.

::: function
int n_poly_mod_is_canonical(const n_poly_t A, nmod_t mod)

Returns whether `A` is normalised with all coefficients reduced modulo
`mod`.
:::

::: function
void n_poly_mod_set_coeff_ui(n_poly_t A, slong j, ulong c, nmod_t mod)

Sets the coefficient of $x^j$ to $c$ reduced modulo `mod`.
:::

::: function
void n_poly_mod_make_monic(n_poly_t A, const n_poly_t B, nmod_t mod)

Sets `A` to `B` divided by its leading coefficient; `B` must be nonzero
with invertible leading coefficient.
:::

::: function
void n_poly_mod_taylor_shift(n_poly_t g, ulong c, nmod_t mod)

Replaces $g(x)$ by $g(x + c)$, in place.
:::

::: function
void n_poly_mod_neg(n_poly_t A, const n_poly_t B, nmod_t mod)
void n_poly_mod_add(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
void n_poly_mod_sub(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
void n_poly_mod_add_ui(n_poly_t res, const n_poly_t poly, ulong c, nmod_t ctx)

Negation, addition and subtraction, and addition of the constant $c$
(reduced first).
:::

::: function
void _n_poly_mod_scalar_mul_nmod(n_poly_t A, const n_poly_t B, ulong c, nmod_t mod)
void _n_poly_mod_scalar_mul_nmod_inplace(n_poly_t A, ulong c, nmod_t mod)
void n_poly_mod_scalar_mul_ui(n_poly_t A, const n_poly_t B, ulong c, nmod_t ctx)

Multiplication by a scalar $c$, which must be reduced for the
underscore versions; these do not normalise (so $c$ should be nonzero).
`n_poly_mod_scalar_mul_ui` reduces $c$ and normalises.
:::

::: function
void n_poly_mod_scalar_addmul_nmod(n_poly_t A, const n_poly_t B, const n_poly_t C, ulong d0, nmod_t ctx)

Sets $A = B + d_0 C$.
:::

::: function
void n_poly_mod_addmul_linear(n_poly_t A, const n_poly_t B, const n_poly_t C, ulong d1, ulong d0, nmod_t mod)

Sets $A = B + C \cdot (d_1 x + d_0)$. `A` must be distinct from `B` and
`C`.
:::

::: function
void n_poly_mod_shift_left_scalar_addmul(n_poly_t A, slong k, ulong c, nmod_t mod)

Multiplies `A` by $x^k + c$, in place. Used to build products of
linear factors.
:::

::: function
void n_poly_mod_product_roots_nmod_vec(n_poly_t A, nn_srcptr r, slong n, nmod_t mod)

Sets `A` to $\prod_{i<n} (x - r_i)$.
:::

::: function
ulong n_poly_mod_evaluate_nmod(const n_poly_t A, ulong c, nmod_t mod)

Returns the value of `A` at $c$.
:::

::: function
ulong _n_poly_eval_pow(n_poly_t P, n_poly_t alphapow, dot_params_t params, nmod_t ctx)
ulong n_poly_mod_eval_pow(n_poly_t P, n_poly_t alphapow, nmod_t ctx)

Return the value of `P` at $\alpha$, where `alphapow` holds the powers
$1, \alpha, \alpha^2, \ldots$ (at least two of them); missing powers are
appended to `alphapow`, so that one table serves many evaluations at the
same point. `_n_poly_eval_pow` takes the dot product parameters for the
length of `P`.
:::

::: function
void n_poly_mod_eval2_pow(ulong * vp, ulong * vm, const n_poly_t P, n_poly_t alphapow, nmod_t ctx)

Sets `*vp` and `*vm` to the values of `P` at $\alpha$ and $-\alpha$,
using and extending the table `alphapow` as above. This is the evaluation
step of the symmetric two-point interpolation (`2sm`) used by Brown's
GCD algorithm.
:::

::: function
ulong n_poly_mod_eval_step2(n_poly_t Acur, const n_poly_t Ainc, nmod_t mod)

One step of the evaluation of a sparse polynomial at consecutive
powers of a point: with $c_i$ = `Acur[i]` and `Ainc` holding the pairs
$(a_i, m_i)$, returns $\sum_i c_i a_i$ and replaces each $c_i$ by $c_i
m_i$. `Ainc` must have twice the length of `Acur`.
:::

::: function
ulong n_poly_mod_div_root(n_poly_t Q, const n_poly_t A, ulong c, nmod_t ctx)

Sets $Q$ to the quotient of $A$ by $x - c$ and returns the remainder
$A(c)$.
:::

::: function
void _n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t ctx)
void n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
void n_poly_mod_mullow(n_poly_t A, const n_poly_t B, const n_poly_t C, slong n, nmod_t mod)

Multiplication (`_nmod_poly_mul`), and multiplication truncated to
length $n$.
:::

::: function
void _n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)
void n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)
void _n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
void n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
void _n_poly_mod_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
void n_poly_mod_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
void _n_poly_mod_divexact(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)
void n_poly_mod_divexact(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)

Euclidean division by a nonzero $B$ with invertible leading coefficient:
quotient, remainder, both, or the quotient of an exact division.
:::

::: function
void n_poly_mod_pow(n_poly_t res, const n_poly_t poly, ulong e, nmod_t ctx)

Sets `res` to $\mathit{poly}^e$.
:::

::: function
ulong n_poly_mod_remove(n_poly_t f, const n_poly_t p, nmod_t ctx)

Divides `f` by the highest power of `p` that divides it and returns
the exponent. `p` must be nonconstant.
:::

::: function
void n_poly_mod_mulmod(n_poly_t res, const n_poly_t poly1, const n_poly_t poly2, const n_poly_t f, nmod_t mod)
void n_poly_mod_mulmod_preinv(n_poly_t A, const n_poly_t B, const n_poly_t C, const n_poly_t M, const n_poly_t Minv, nmod_t ctx)

Set the output to $B C \bmod M$. The `_preinv` version takes `Minv`,
the inverse of the reversal of $M$ (see `nmod_poly_mulmod_preinv`),
requires $B$, $C$ reduced modulo $M$, and raises an exception otherwise.
:::

::: function
void n_poly_mod_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B, nmod_t mod)
void n_poly_mod_xgcd(n_poly_t G, n_poly_t S, n_poly_t T, const n_poly_t A, const n_poly_t B, nmod_t mod)

Set $G$ to the monic greatest common divisor of $A$ and $B$, and for
`n_poly_mod_xgcd` the cofactors with $SA + TB = G$ (see `nmod_poly_xgcd`).
The modulus must be prime.
:::

::: function
int n_poly_mod_invmod(n_poly_t A, const n_poly_t B, const n_poly_t P, nmod_t mod)

If $B$ is invertible modulo $P$, sets $A$ to the inverse and returns 1;
otherwise returns 0.
:::

::: function
void n_poly_mod_inv_series(n_poly_t Qinv, const n_poly_t Q, slong n, nmod_t mod)
void n_poly_mod_div_series(n_poly_t Q, const n_poly_t A, const n_poly_t B, slong order, nmod_t ctx)

Power series inverse of $Q$ modulo $x^n$, and quotient $A/B$ modulo
$x^{\mathit{order}}$; the constant term of the divisor must be
invertible.
:::

## Elements of a finite field as word vectors {#n-fq}

In this section an element of $\mathbb{F}_q$ is a `ulong *` pointing to
$d$ words, $d$ being the degree of the field (`fq_nmod_ctx_degree`). The
underscore functions take $d$ and an `nmod_t`, or need a scratch array
`t` whose size in multiples of $d$ is given by the macros
`N_FQ_REDUCE_ITCH` (2), `N_FQ_MUL_ITCH` (4), `N_FQ_LAZY_ITCH` (6) and
`N_FQ_INV_ITCH` (1). Outputs may alias inputs only where noted for the
underscore functions.

::: function
nmod_t fq_nmod_ctx_mod(const fq_nmod_ctx_t ctx)

Returns the `nmod_t` of the prime field of `ctx`.
:::

::: function
int _n_fq_is_zero(const ulong * a, slong d)
int _n_fq_is_one(const ulong * a, slong d)
int _n_fq_is_ui(const ulong * a, slong d)
int n_fq_is_one(const ulong * a, const fq_nmod_ctx_t ctx)

Return whether the element is 0, 1, or in the prime field (all
coefficients except the constant one zero).
:::

::: function
void _n_fq_zero(ulong * a, slong d)
void _n_fq_one(ulong * a, slong d)
void _n_fq_set_nmod(ulong * a, ulong b, slong d)
void n_fq_gen(ulong * a, const fq_nmod_ctx_t ctx)
void _n_fq_set(ulong * a, const ulong * b, slong d)
void _n_fq_swap(ulong * a, ulong * b, slong d)

Set the element to 0, 1, the prime field element $b$, the generator
$t$ of the field, or a copy of `b`; or swap two elements.
:::

::: function
int _n_fq_equal(ulong * a, const ulong * b, slong d)
int n_fq_equal_fq_nmod(const ulong * a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)

Return whether two elements are equal (the second compares with an
`fq_nmod_t`).
:::

::: function
int n_fq_is_canonical(const ulong * a, const fq_nmod_ctx_t ctx)

Returns whether all $d$ words are reduced modulo $p$.
:::

::: function
void n_fq_randtest_not_zero(ulong * a, flint_rand_t state, const fq_nmod_ctx_t ctx)

Sets `a` to a random nonzero element.
:::

::: function
char * n_fq_get_str_pretty(const ulong * a, const fq_nmod_ctx_t ctx)
int n_fq_fprint_pretty(FILE * file, const ulong * a, const fq_nmod_ctx_t ctx)
void n_fq_print_pretty(const ulong * a, const fq_nmod_ctx_t ctx)

Convert the element to a string (to be freed with `flint_free`), or
print it, as a polynomial in the generator.
:::

::: function
void n_fq_get_n_poly(n_poly_t a, const ulong * b, const fq_nmod_ctx_t ctx)
void _n_fq_set_n_poly(ulong * a, const ulong * bcoeffs, slong blen, const fq_nmod_ctx_t ctx)
void n_fq_get_fq_nmod(fq_nmod_t a, const ulong * b, const fq_nmod_ctx_t ctx)
void n_fq_set_fq_nmod(ulong * a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)

Convert between an element and an `n_poly_t` (the polynomial of degree
less than $d$; `_n_fq_set_n_poly` reduces a polynomial of any length
modulo the defining polynomial) or an `fq_nmod_t`.
:::

::: function
void _n_fq_neg(ulong * a, const ulong * b, slong d, nmod_t mod)
void _n_fq_add(ulong * a, const ulong * b, const ulong * c, slong d, nmod_t mod)
void _n_fq_sub(ulong * a, const ulong * b, const ulong * c, slong d, nmod_t mod)
void n_fq_add(ulong * a, const ulong * b, const ulong * c, const fq_nmod_ctx_t ctx)
void n_fq_sub(ulong * a, const ulong * b, const ulong * c, const fq_nmod_ctx_t ctx)
void n_fq_add_si(ulong * a, const ulong * b, slong c, const fq_nmod_ctx_t ctx)
void n_fq_add_fq_nmod(ulong * a, const ulong * b, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
void n_fq_sub_fq_nmod(ulong * a, const ulong * b, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

Negation, addition and subtraction, including with a signed integer
or an `fq_nmod_t`.
:::

::: function
void _n_fq_mul_ui(ulong * a, const ulong * b, ulong c, slong d, nmod_t mod)

Multiplies an element by the prime field element $c$ (which need not
be reduced).
:::

::: function
void _n_fq_mul2(ulong * t, const ulong * b, const ulong * c, const fq_nmod_ctx_t ctx)
void _n_fq_madd2(ulong * a, const ulong * b, const ulong * c, const fq_nmod_ctx_t ctx, ulong * t)
void _n_fq_reduce(ulong * a, ulong * b, slong blen, const fq_nmod_ctx_t ctx, ulong * t)
void _n_fq_reduce2(ulong * a, ulong * b, const fq_nmod_ctx_t ctx, ulong * t)

Unreduced multiplication: `_n_fq_mul2` sets the $2d - 1$ coefficients of
the product of two elements as polynomials in $t$, `_n_fq_madd2` adds it
to $a$; `_n_fq_reduce` reduces a polynomial of length `blen` modulo the
defining polynomial to an element, and `_n_fq_reduce2` does so for length
$2d - 1$. Delaying the reduction makes sums of products cheaper.
:::

::: function
void _n_fq_mul(ulong * a, const ulong * b, const ulong * c, const fq_nmod_ctx_t ctx, ulong * t)
void _n_fq_addmul(ulong * a, const ulong * b, const ulong * c, const ulong * e, const fq_nmod_ctx_t ctx, ulong * t)
void n_fq_mul(ulong * a, const ulong * b, const ulong * c, const fq_nmod_ctx_t ctx)
void n_fq_mul_fq_nmod(ulong * a, const ulong * b, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
void n_fq_addmul(ulong * a, const ulong * b, const ulong * c, const ulong * d, const fq_nmod_ctx_t ctx)

Multiplication $a = bc$ and multiply-add $a = b + c e$ (`n_fq_addmul`:
$a = b + c d$). The underscore versions need scratch space of
`N_FQ_MUL_ITCH` times $d$ words.
:::

::: function
int _n_fq_dot_lazy_size(slong len, const fq_nmod_ctx_t ctx)

Returns the number of words (1, 2 or 3) needed to accumulate a sum of
`len` products of elements without reduction, or 0 if 3 words do not
suffice or $d > 30$.
:::

::: function
void _n_fq_mul2_lazy1(ulong * a, const ulong * b, const ulong * c, slong d)
void _n_fq_madd2_lazy1(ulong * a, const ulong * b, const ulong * c, slong d)
void _n_fq_reduce2_lazy1(ulong * a, slong d, nmod_t ctx)
void _n_fq_mul2_lazy2(ulong * a, const ulong * b, const ulong * c, slong d)
void _n_fq_madd2_lazy2(ulong * a, const ulong * b, const ulong * c, slong d)
void _n_fq_reduce2_lazy2(ulong * a, slong d, nmod_t ctx)
void _n_fq_mul2_lazy3(ulong * a, const ulong * b, const ulong * c, slong d)
void _n_fq_madd2_lazy3(ulong * a, const ulong * b, const ulong * c, slong d)
void _n_fq_reduce2_lazy3(ulong * a, slong d, nmod_t ctx)

Lazy dot products with 1, 2 or 3 words per coefficient, as chosen by
`_n_fq_dot_lazy_size`: `mul2` sets and `madd2` adds the unreduced product
of two elements in $a$ (6d words), and `reduce2` reduces the accumulated
coefficients modulo $p$, leaving $2d - 1$ coefficients to be reduced with
`_n_fq_reduce2`.
:::

::: function
void _n_fq_inv(ulong * a, const ulong * b, const fq_nmod_ctx_t ctx, ulong * t)
void n_fq_inv(ulong * a, const ulong * b, const fq_nmod_ctx_t ctx)

Inversion of a nonzero element.
:::

::: function
void _n_fq_pow_ui(ulong * a, const ulong * b, ulong e, const fq_nmod_ctx_t ctx)
void n_fq_pow_ui(ulong * a, const ulong * b, ulong e, const fq_nmod_ctx_t ctx)
void n_fq_pow_fmpz(ulong * a, const ulong * b, const fmpz_t e, const fq_nmod_ctx_t ctx)

Powering by binary exponentiation.
:::

## Polynomials over a finite field {#n-fq-poly}

Here an `n_fq_poly_t` holds `length` coefficients of $d$ words each.
The macros `n_fq_poly_init`, `n_fq_poly_clear`, `n_fq_poly_swap`,
`n_fq_poly_degree`, `n_fq_poly_zero` and `n_fq_poly_is_zero` are the
corresponding `n_poly` functions.

::: function
void n_fq_poly_init2(n_fq_poly_t A, slong alloc, const fq_nmod_ctx_t ctx)

Initialises `A` with space for `alloc` coefficients in $\mathbb{F}_q$.
:::

::: function
int n_fq_poly_is_canonical(const n_fq_poly_t a, const fq_nmod_ctx_t ctx)

Returns whether `A` is normalised with canonical coefficients.
:::

::: function
void _n_fq_poly_normalise(n_fq_poly_t A, slong d)

Normalises `A`, whose coefficients have $d$ words.
:::

::: function
void n_fq_poly_one(n_fq_poly_t A, const fq_nmod_ctx_t ctx)
void _n_fq_poly_one(n_fq_poly_t A, slong d)
void n_fq_poly_set(n_fq_poly_t A, const n_fq_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_poly_set_n_fq(n_poly_t A, const ulong * c, const fq_nmod_ctx_t ctx)
void n_fq_poly_set_fq_nmod(n_poly_t A, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

Set the polynomial to 1, a copy of `B`, or the constant $c$.
:::

::: function
int n_fq_poly_is_one(n_fq_poly_t A, const fq_nmod_ctx_t ctx)
int n_fq_poly_equal(const n_fq_poly_t A, const n_fq_poly_t B, const fq_nmod_ctx_t ctx)

Return whether the polynomial is 1, or equal to the other one.
:::

::: function
void n_fq_poly_get_coeff_n_fq(ulong * c, const n_poly_t A, slong e, const fq_nmod_ctx_t ctx)
void n_fq_poly_get_coeff_fq_nmod(fq_nmod_t c, const n_poly_t A, slong e, const fq_nmod_ctx_t ctx)
void n_fq_poly_set_coeff_n_fq(n_poly_t A, slong j, const ulong * c, const fq_nmod_ctx_t ctx)
void n_fq_poly_set_coeff_fq_nmod(n_poly_t A, slong j, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

Get or set the coefficient of $x^e$ (respectively $x^j$), as an `n_fq`
element or an `fq_nmod_t`.
:::

::: function
void n_fq_poly_get_fq_nmod_poly(fq_nmod_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_poly_set_fq_nmod_poly(n_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Convert to and from an `fq_nmod_poly_t`.
:::

::: function
void n_fq_poly_randtest(n_fq_poly_t A, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)

Sets `A` to a random polynomial of length at most `len`.
:::

::: function
void n_fq_poly_print_pretty(const n_fq_poly_t A, const char * x, const fq_nmod_ctx_t ctx)

Prints `A` with variable name `x`.
:::

::: function
void n_fq_poly_shift_left(n_poly_t A, const n_poly_t B, slong n, const fq_nmod_ctx_t ctx)
void n_fq_poly_shift_right(n_poly_t A, const n_poly_t B, slong n, const fq_nmod_ctx_t ctx)
void n_fq_poly_truncate(n_poly_t A, slong len, const fq_nmod_ctx_t ctx)

Multiply by $x^n$, divide by $x^n$ discarding low terms, or truncate
to length `len`.
:::

::: function
void n_fq_poly_neg(n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_poly_add(n_poly_t A, const n_poly_t B, const n_poly_t C, const fq_nmod_ctx_t ctx)
void n_fq_poly_sub(n_poly_t A, const n_poly_t B, const n_poly_t C, const fq_nmod_ctx_t ctx)
void n_fq_poly_add_si(n_poly_t A, const n_poly_t B, slong c, const fq_nmod_ctx_t ctx)

Negation, addition, subtraction, and addition of an integer constant.
:::

::: function
void n_fq_poly_scalar_mul_ui(n_poly_t A, const n_poly_t B, ulong c, const fq_nmod_ctx_t ctx)
void n_fq_poly_scalar_mul_n_fq(n_poly_t A, const n_poly_t B, const ulong * c, const fq_nmod_ctx_t ctx)
void n_fq_poly_scalar_addmul_n_fq(n_fq_poly_t A, const n_fq_poly_t B, const n_fq_poly_t C, const ulong * d, const fq_nmod_ctx_t ctx)

Multiplication by a scalar (a prime field element $c$ or an `n_fq`
element), and $A = B + C d$.
:::

::: function
void n_fq_poly_shift_left_scalar_submul(n_poly_t A, slong k, const ulong * c, const fq_nmod_ctx_t ctx)

Multiplies `A` by $x^k - c$, in place.
:::

::: function
void n_fq_poly_make_monic(n_fq_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)

Sets `A` to the nonzero `B` divided by its leading coefficient.
:::

::: function
void n_fq_poly_evaluate_fq_nmod(fq_nmod_t e, const n_poly_t A, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
void n_fq_poly_evaluate_n_fq(ulong * e, const n_poly_t A, const ulong * c, const fq_nmod_ctx_t ctx)

Set $e$ to the value of `A` at $c$.
:::

::: function
void n_fq_poly_eval_pow(ulong * ev, const n_fq_poly_t A, n_fq_poly_t alphapow, const fq_nmod_ctx_t ctx)

Sets `ev` to the value of `A` at $\alpha$, given the table `alphapow`
of powers $1, \alpha, \ldots$ (as an `n_fq_poly`), which is extended as
needed.
:::

::: function
void _n_fq_poly_mul_(ulong * A, const ulong * B, slong Blen, const ulong * C, slong Clen, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void n_fq_poly_mul_(n_poly_t A, const n_poly_t B, const n_poly_t C, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void n_fq_poly_mul(n_poly_t A, const n_poly_t B, const n_poly_t C, const fq_nmod_ctx_t ctx)
void n_fq_poly_mullow(n_poly_t A, const n_poly_t B, const n_poly_t C, slong order, const fq_nmod_ctx_t ctx)

Multiplication (truncated to length `order` for `n_fq_poly_mullow`).
If one factor has fewer than 4 terms the product is computed
classically with delayed reduction; otherwise by Kronecker substitution:
each coefficient is placed in a slot of $2d - 1$ words, the packed
polynomials are multiplied with `_nmod_poly_mul`, and the slots are reduced
modulo the defining polynomial. The trailing-underscore
versions take a stack of temporaries; `_n_fq_poly_mul_` works on arrays
of `Blen` and `Clen` coefficients and writes $\mathit{Blen} + \mathit{Clen}
- 1$.
:::

::: function
void n_fq_poly_pow(n_poly_t A, const n_poly_t B, ulong e, const fq_nmod_ctx_t ctx)

Sets `A` to $B^e$.
:::

::: function
ulong n_fq_poly_remove(n_poly_t f, const n_poly_t g, const fq_nmod_ctx_t ctx)

Divides `f` by the highest power of `g` that divides it and returns
the exponent.
:::

::: function
void _n_fq_poly_divrem_basecase_(ulong * Q, ulong * A, const ulong * AA, slong Alen, const ulong * B, slong Blen, const ulong * invB, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void _n_fq_poly_rem_basecase_(ulong * Q, ulong * A, const ulong * AA, slong Alen, const ulong * B, slong Blen, const ulong * invB, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void n_fq_poly_divrem_divconquer_(n_poly_t Q, n_poly_t R, const n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void n_fq_poly_divrem_(n_poly_t Q, n_poly_t R, const n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void n_fq_poly_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_poly_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)

Euclidean division by a nonzero $B$: quotient and remainder, or the
remainder only. The basecase functions work on arrays, with `invB` the
inverse of the leading coefficient of $B$, `AA` the dividend and `A`
receiving the remainder; `n_fq_poly_divrem_` uses the classical algorithm
for divisors of length below `N_FQ_POLY_DIVREM_DIVCONQUER_CUTOFF` (20)
and divide and conquer otherwise.
:::

::: function
void n_fq_poly_gcd_(n_poly_t G, const n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx, n_poly_stack_t St)
void n_fq_poly_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_poly_xgcd(n_poly_t G, n_poly_t S, n_poly_t T, const n_poly_t B, const n_poly_t C, const fq_nmod_ctx_t ctx)

Monic greatest common divisor by the Euclidean algorithm, and for
`n_fq_poly_xgcd` cofactors $S$, $T$ with $SB + TC = G$.
:::

::: function
void n_fq_poly_mulmod(n_poly_t A, const n_poly_t B, const n_poly_t C, const n_poly_t M, const fq_nmod_ctx_t ctx)

Sets $A = BC \bmod M$.
:::

::: function
void n_fq_poly_inv_series(n_poly_t A, const n_poly_t B, slong order, const fq_nmod_ctx_t ctx)

Sets `A` to the inverse of the power series `B` modulo
$x^{\mathit{order}}$; the constant term must be nonzero.
:::

::: function
void n_fq_poly_product_roots_n_fq(n_poly_t M, const ulong * H, slong length, const fq_nmod_ctx_t ctx, n_poly_stack_t St)

Sets `M` to $\prod_i (x - H_i)$ for the `length` elements of `H`.
:::

## Bivariate polynomials {#bpoly}

An `n_bpoly_t` $A = \sum_i A_i(y) x^i$ has `coeffs[i]` $= A_i$. The
macros `n_fq_bpoly_init`, `n_fq_bpoly_clear`, `n_fq_bpoly_zero`,
`n_fq_bpoly_fit_length`, `n_fq_bpoly_normalise`, `n_fq_bpoly_degree0` and
`n_fq_bpoly_degree1` are the corresponding `n_bpoly` functions.

::: function
void n_bpoly_init(n_bpoly_t A)
void n_bpoly_clear(n_bpoly_t A)
void n_bpoly_realloc(n_bpoly_t A, slong len)
void n_bpoly_fit_length(n_bpoly_t A, slong len)
void n_bpoly_swap(n_bpoly_t A, n_bpoly_t B)

Memory management: initialise as zero, free, ensure space for `len`
coefficients $A_i$ (all initialised), and swap.
:::

::: function
void n_bpoly_normalise(n_bpoly_t A)

Removes zero leading coefficients $A_i$.
:::

::: function
void n_bpoly_zero(n_bpoly_t A)
void n_bpoly_one(n_bpoly_t A)
void _n_bpoly_set(n_bpoly_t A, const n_bpoly_t B)
void n_bpoly_set(n_bpoly_t A, const n_bpoly_t B)

Set to 0, 1 or a copy of `B` (`_n_bpoly_set` requires `A` and `B` to
be distinct).
:::

::: function
int n_bpoly_is_zero(const n_bpoly_t A)
int n_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B)

Return whether the polynomial is zero, or equal to the other one.
:::

::: function
ulong n_bpoly_get_coeff(const n_bpoly_t A, slong e0, slong e1)
void n_bpoly_set_coeff(n_bpoly_t A, slong e0, slong e1, ulong c)
void n_bpoly_set_coeff_nonzero(n_bpoly_t A, slong e0, slong e1, ulong c)

Get or set the coefficient of $x^{e_0} y^{e_1}$; the `_nonzero` version
requires $c \ne 0$.
:::

::: function
slong n_bpoly_degree0(const n_bpoly_t A)
slong n_bpoly_degree1(const n_bpoly_t A)
ulong n_bpoly_bidegree(const n_bpoly_t A)

Return the degree in $x$, the degree in $y$ (the maximum over the
$A_i$), and the pair (degree in $x$, degree in $y$ of the leading
coefficient $A_{\mathrm{lead}}$) packed into a word with $x$ in the upper half.
:::

::: function
void n_bpoly_set_poly_gen0(n_bpoly_t A, const n_poly_t B)
void n_bpoly_set_poly_gen1(n_bpoly_t A, const n_poly_t B)

Set `A` to the univariate polynomial `B` in $x$, respectively in $y$.
:::

::: function
void n_bpoly_print_pretty(const n_bpoly_t A, const char * xvar, const char * yvar)

Prints `A` with the variable names `xvar` and `yvar`.
:::

::: function
int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod)

Returns whether `A` is normalised with every $A_i$ canonical modulo
`mod`.
:::

::: function
void n_bpoly_scalar_mul_nmod(n_bpoly_t A, ulong c, nmod_t ctx)

Multiplies `A` by the nonzero scalar $c$.
:::

::: function
void n_bpoly_mod_derivative_gen0(n_bpoly_t A, const n_bpoly_t B, nmod_t ctx)

Sets `A` to $\partial B/\partial x$; `A` and `B` must be distinct.
:::

::: function
void n_bpoly_mod_content_last(n_poly_t g, const n_bpoly_t A, nmod_t ctx)
void n_bpoly_mod_divexact_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx)
void n_bpoly_mod_mul_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx)

Treat `A` as a polynomial in $x$ over $\mathbb{F}_p[y]$: set $g$ to the
content $\gcd_i A_i$, or divide each $A_i$ exactly by, or multiply it
by, the polynomial $b(y)$.
:::

::: function
void n_bpoly_mod_make_primitive(n_poly_t g, n_bpoly_t A, nmod_t ctx)

Sets $g$ to the content of `A` in $\mathbb{F}_p[y]$ and divides it out,
then scales so that the leading coefficient of `A` (in $x$, then $y$) is
1, multiplying $g$ by the old leading coefficient.
:::

::: function
void n_bpoly_mod_taylor_shift_gen0(n_bpoly_t A, ulong c, nmod_t ctx)
void n_bpoly_mod_taylor_shift_gen1(n_bpoly_t A, const n_bpoly_t B, ulong c, nmod_t ctx)

Replace $x$ by $x + c$ in `A` (in place), respectively set `A` to `B`
with $y$ replaced by $y + c$.
:::

::: function
void n_bpoly_mod_add(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, nmod_t ctx)
void n_bpoly_mod_sub(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, nmod_t ctx)
void n_bpoly_mod_mul(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, nmod_t ctx)

Addition, subtraction and multiplication.
:::

::: function
int n_bpoly_mod_divides(n_bpoly_t Q, const n_bpoly_t A, const n_bpoly_t B, nmod_t ctx)

If $B$ divides $A$, sets $Q = A/B$ and returns 1, otherwise returns 0.
:::

::: function
void n_bpoly_mod_mul_series(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, slong order, nmod_t ctx)
void n_bpoly_mod_divrem_series(n_bpoly_t Q, n_bpoly_t R, const n_bpoly_t A, const n_bpoly_t B, slong order, nmod_t ctx)

Multiplication, and division with remainder in $x$, with all
coefficients $A_i(y)$ truncated modulo $y^{\mathit{order}}$ (arithmetic in
$(\mathbb{F}_p[y]/(y^{\mathit{order}}))[x]$, used in Hensel lifting). The
leading coefficient of $B$ must be invertible.
:::

::: function
void n_bpoly_mod_interp_reduce_2sm_poly(n_poly_t Ap, n_poly_t Am, const n_bpoly_t A, n_poly_t alphapow, nmod_t mod)
void n_bpoly_mod_interp_lift_2sm_poly(slong * deg1, n_bpoly_t T, const n_poly_t A, const n_poly_t B, ulong alpha, nmod_t mod)
int n_bpoly_mod_interp_crt_2sm_poly(slong * deg1, n_bpoly_t F, n_bpoly_t T, n_poly_t A, n_poly_t B, const n_poly_t modulus, n_poly_t alphapow, nmod_t mod)

Helpers for interpolation in $y$ at pairs of points $\pm\alpha$:
`reduce` evaluates `A` at $y = \alpha$ and $y = -\alpha$ (giving `Ap`,
`Am`); `lift` sets $T$ to the bivariate polynomial of degree at most 1 in
$y$ with these images; `crt` updates the interpolant $F$ modulo `modulus`
with new images (using $T$ as scratch), returns whether $F$ changed and
sets `*deg1` to its degree in $y$.
:::

::: function
int n_bpoly_mod_gcd_brown_smprime(n_bpoly_t G, n_bpoly_t Abar, n_bpoly_t Bbar, n_bpoly_t A, n_bpoly_t B, nmod_t ctx, n_poly_bpoly_stack_t Sp)

Brown's dense modular GCD algorithm [@Bro1971] for bivariate
polynomials over a word-size prime field: evaluates $y$ at pairs of
points $\pm \alpha$, computes univariate GCDs and interpolates. Sets $G =
\gcd(A, B)$, $\bar A = A/G$, $\bar B = B/G$ and returns 1, or returns 0
if it ran out of evaluation points (the prime is too small). $A$ and $B$
are clobbered. `Sp` provides the temporaries. See
[](algorithms/multivariate/overview).
:::

::: function
int n_fq_bpoly_is_canonical(const n_fq_bpoly_t A, const fq_nmod_ctx_t ctx)

Returns whether `A` is normalised with canonical coefficients.
:::

::: function
void n_fq_bpoly_one(n_fq_bpoly_t A, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_set(n_fq_bpoly_t A, const n_fq_bpoly_t B, const fq_nmod_ctx_t ctx)
int n_fq_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B, const fq_nmod_ctx_t ctx)

Set to 1, copy, or compare.
:::

::: function
void n_fq_bpoly_get_coeff_n_fq(ulong * c, const n_bpoly_t A, slong e0, slong e1, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_get_coeff_fq_nmod(fq_nmod_t c, const n_bpoly_t A, slong e0, slong e1, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_set_coeff_n_fq(n_fq_bpoly_t A, slong e0, slong e1, const ulong * c, const fq_nmod_ctx_t ctx)

Get or set the coefficient of $x^{e_0} y^{e_1}$.
:::

::: function
void n_fq_bpoly_set_fq_nmod_poly_gen0(n_bpoly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_set_n_fq_poly_gen0(n_bpoly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_set_n_fq_poly_gen1(n_bpoly_t A, const n_poly_t B, const fq_nmod_ctx_t ctx)

Set `A` to a univariate polynomial in $x$ (`gen0`) or in $y$ (`gen1`).
:::

::: function
void n_fq_bpoly_derivative_gen0(n_bpoly_t A, const n_bpoly_t B, const fq_nmod_ctx_t ctx)

Sets `A` to $\partial B/\partial x$.
:::

::: function
void n_fq_bpoly_scalar_mul_n_fq(n_fq_bpoly_t A, const ulong * c, const fq_nmod_ctx_t ctx)

Multiplies `A` by the scalar $c$.
:::

::: function
void n_fq_bpoly_taylor_shift_gen0_fq_nmod(n_bpoly_t A, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_taylor_shift_gen0_n_fq(n_fq_bpoly_t A, const ulong * alpha, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_taylor_shift_gen1_fq_nmod(n_bpoly_t A, const n_bpoly_t B, const fq_nmod_t c_, const fq_nmod_ctx_t ctx)

Replace $x$ by $x + \alpha$ in `A` (in place), or set `A` to `B` with
$y$ replaced by $y + c$.
:::

::: function
int n_fq_bpoly_gcd_brown_smprime(n_fq_bpoly_t G, n_fq_bpoly_t Abar, n_fq_bpoly_t Bbar, n_fq_bpoly_t A, n_fq_bpoly_t B, const fq_nmod_ctx_t ctx, n_poly_bpoly_stack_t Sp)

Brown's GCD algorithm over $\mathbb{F}_q$, with the same conventions as
`n_bpoly_mod_gcd_brown_smprime`.
:::

::: function
void n_fq_bpoly_print_pretty(const n_fq_bpoly_t A, const char * xvar, const char * yvar, const fq_nmod_ctx_t ctx)

Prints `A` with the given variable names.
:::

## Trivariate and sparse polynomials {#sparse}

::: function
void n_tpoly_init(n_tpoly_t A)
void n_tpoly_clear(n_tpoly_t A)
void n_tpoly_fit_length(n_tpoly_t A, slong len)
void n_tpoly_swap(n_tpoly_t A, n_tpoly_t B)

Memory management for `n_tpoly_t`.
:::

::: function
void n_polyu_init(n_polyu_t A)
void n_polyu_clear(n_polyu_t A)
void n_polyu_realloc(n_polyu_t A, slong len)
void n_polyu_fit_length(n_polyu_t A, slong len)
void n_polyu_swap(n_polyu_t A, n_polyu_t B)

Memory management for `n_polyu_t`.
:::

::: function
void n_polyu3_degrees(slong * deg0, slong * deg1, slong * deg2, const n_polyu_t A)

Sets the three degrees of the nonzero terms of `A`, whose exponents
pack three variables; $-1$ each for the zero polynomial.
:::

::: function
void n_polyu3_print_pretty(const n_polyu_t A, const char * gen0, const char * gen1, const char * var2)

Prints `A` with the three variable names.
:::

::: function
void n_polyun_init(n_polyun_t A)
void n_polyun_clear(n_polyun_t A)
void n_polyun_realloc(n_polyun_t A, slong len)
void n_polyun_fit_length(n_polyun_t A, slong len)
void n_polyun_swap(n_polyun_t A, n_polyun_t B)

Memory management for `n_polyun_t`; all `alloc` coefficient
polynomials are initialised.
:::

::: function
int n_polyun_mod_is_canonical(const n_polyun_t A, nmod_t mod)

Returns whether the terms are in strictly decreasing order of exponents
with nonzero normalised coefficients reduced modulo $n$.
:::

::: function
void n_polyun_set(n_polyun_t A, const n_polyun_t B)
void n_fq_polyun_set(n_fq_polyun_t A, const n_fq_polyun_t B, const fq_nmod_ctx_t ctx)
void n_polyun_one(n_polyun_t A)
int n_polyun_equal(const n_polyun_t A, const n_polyun_t B)

Copy, set to 1, or compare.
:::

::: function
ulong n_polyu1n_bidegree(n_polyun_t A)

Returns the exponent of the leading term and the degree of its
coefficient, packed into a word as for `n_bpoly_bidegree`.
:::

::: function
void n_polyu1n_print_pretty(const n_polyun_t A, const char * var0, const char * varlast)
void n_polyu2n_print_pretty(const n_polyun_t A, const char * gen0, const char * gen1, const char * varlast)
void n_polyu3n_print_pretty(const n_polyun_t A, const char * gen0, const char * gen1, const char * var2, const char * varlast)

Print `A`, whose exponents pack 1, 2 or 3 variables, with the given
variable names.
:::

::: function
int n_polyu1n_mod_gcd_brown_smprime(n_polyun_t G, n_polyun_t Abar, n_polyun_t Bbar, n_polyun_t A, n_polyun_t B, nmod_t ctx, n_poly_polyun_stack_t St)

Brown's GCD algorithm for polynomials stored as `n_polyun` with one
packed variable, with the same conventions as
`n_bpoly_mod_gcd_brown_smprime`.
:::

::: function
slong n_polyun_product_roots(n_polyun_t M, const n_polyun_t H, nmod_t ctx)
slong n_fq_polyun_product_roots(n_fq_polyun_t M, const n_fq_polyun_t H, const fq_nmod_ctx_t ctx, n_poly_stack_t St)

For each term of `H`, set the corresponding term of `M` (with the same
exponent) to the product $\prod_j (x - h_j)$ over the coefficients $h_j$
of that term's coefficient polynomial, and return the largest length of
these products. Used to build the master polynomials of Zippel's
interpolation.
:::

## Power caches {#pow-cache}

A power cache for $b$ holds three tables in `n_poly_t`s: `pos_direct`
($b^0, b^1, \ldots$, extended on demand up to $b^{50}$), `pos_bin` ($b,
b^2, b^3$ and then $b^{4^i}$, $2 b^{4^i}$, $3 b^{4^i}$ for base-4
exponentiation) and `neg_direct` ($b^{-1}, b^{-2}, \ldots$). It is used when
the same base is raised to many different exponents, as when evaluating
the monomials of a sparse polynomial.

::: function
void nmod_pow_cache_start(ulong b, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct)
void n_fq_pow_cache_start_n_fq(const ulong * b, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, const fq_nmod_ctx_t ctx)
void n_fq_pow_cache_start_fq_nmod(const fq_nmod_t b, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, const fq_nmod_ctx_t ctx)

Initialise the cache tables for the base $b$.
:::

::: function
ulong nmod_pow_cache_mulpow_ui(ulong a, ulong e, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, nmod_t ctx)
ulong nmod_pow_cache_mulpow_neg_ui(ulong a, ulong e, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, nmod_t ctx)
ulong nmod_pow_cache_mulpow_fmpz(ulong a, const fmpz_t e, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, nmod_t ctx)
void n_fq_pow_cache_mulpow_ui(ulong * r, const ulong * a, ulong e, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, const fq_nmod_ctx_t ctx)
void n_fq_pow_cache_mulpow_neg_ui(ulong * r, const ulong * a, ulong e, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, const fq_nmod_ctx_t ctx)
void n_fq_pow_cache_mulpow_fmpz(ulong * r, const ulong * a, const fmpz_t e, n_poly_t pos_direct, n_poly_t pos_bin, n_poly_t neg_direct, const fq_nmod_ctx_t ctx)

Return (or set $r$ to) $a b^e$, respectively $a b^{-e}$, using and
extending the cache.
:::

## Evaluation and interpolation at special points {#evaluation}

An `nmod_eval_interp_t` converts polynomials of degree at most `deg`
between coefficient form and their values at a fixed set of points, so
that products of many small polynomials can be done pointwise (in
*evaluation form*). For $p \equiv 1 \pmod 4$ the points are $0$ and
$\pm\beta_j, \pm i\beta_j$ for a square root $i$ of $-1$ and suitable
$\beta_j$; otherwise $0$ and $\pm \beta_j$. The transforms are products
with small precomputed matrices.

::: function
void nmod_eval_interp_init(nmod_eval_interp_t E)
void nmod_eval_interp_clear(nmod_eval_interp_t E)

Initialise and free `E`.
:::

::: function
int nmod_eval_interp_set_degree_modulus(nmod_eval_interp_t E, slong deg, nmod_t ctx)

Prepares `E` for polynomials of degree at most `deg` modulo the odd
prime $p$ of `ctx`. Returns 0 if this is not possible ($p$ is 2 or not
larger than `deg`, or not enough points exist), 1 otherwise.
:::

::: function
slong nmod_eval_interp_eval_length(nmod_eval_interp_t E)

Returns the number of evaluation points ($4d + 1$ or $2d + 1$).
:::

::: function
void nmod_eval_interp_from_coeffs_poly(n_poly_t v, const n_poly_t a, nmod_eval_interp_t E, nmod_t ctx)
void nmod_eval_interp_to_coeffs_poly(n_poly_t a, const n_poly_t v, nmod_eval_interp_t E, nmod_t ctx)
void nmod_eval_interp_from_coeffs_n_fq_poly(n_fq_poly_t v, const n_fq_poly_t a, nmod_eval_interp_t E, const fq_nmod_ctx_t ctx)
void nmod_eval_interp_to_coeffs_n_fq_poly(n_fq_poly_t a, const n_fq_poly_t v, nmod_eval_interp_t E, const fq_nmod_ctx_t ctx)

Convert from coefficients `a` to values `v`, and back. The `n_fq_poly`
versions transform polynomials over $\mathbb{F}_q$ coefficientwise.
:::

::: function
void nmod_evals_zero(n_poly_t a)
void nmod_evals_add_inplace(n_poly_t a, n_poly_t b, slong len, nmod_t ctx)
void nmod_evals_mul(n_poly_t a, n_poly_t b, n_poly_t c, slong len, nmod_t ctx)
void nmod_evals_addmul(n_poly_t a, n_poly_t b, n_poly_t c, slong len, nmod_t ctx)
void nmod_evals_fmma(n_poly_t a, n_poly_t b, n_poly_t c, n_poly_t d, n_poly_t e, slong len, nmod_t ctx)
void n_fq_evals_zero(n_fq_poly_t a)
void n_fq_evals_add_inplace(n_fq_poly_t a, n_fq_poly_t b, slong len, const fq_nmod_ctx_t ctx)
void n_fq_evals_mul(n_fq_poly_t a, n_fq_poly_t b, n_fq_poly_t c, slong len, const fq_nmod_ctx_t ctx)
void n_fq_evals_addmul(n_fq_poly_t a, n_fq_poly_t b, n_fq_poly_t c, slong len, const fq_nmod_ctx_t ctx)
void n_fq_evals_fmma(n_fq_poly_t a, n_fq_poly_t b, n_fq_poly_t c, n_fq_poly_t f, n_fq_poly_t e, slong len, const fq_nmod_ctx_t ctx)

Pointwise arithmetic on vectors of `len` values, where a zero vector may
be stored with length 0: $a = 0$, $a \mathrel{+}= b$, $a = bc$, $a
\mathrel{+}= bc$, and $a = bc + de$ (`fmma`).
:::

## Zippel interpolation helpers {#zippel}

In Zippel's sparse interpolation [@Zip1979] the values of a polynomial
with known monomials $m_1, \ldots, m_t$ at the powers $\alpha^1,
\alpha^2, \ldots$ of a point give a transposed Vandermonde system for the
coefficients. With the *master* polynomial $\prod_i (x - m_i)$ (from
`n_polyun_product_roots`), the system is solved in $O(t^2)$
operations.

::: function
ulong _nmod_zip_eval_step(ulong * cur, const ulong * inc, const ulong * coeffs, slong length, nmod_t ctx)
void _n_fq_zip_eval_step(ulong * res, ulong * cur, const ulong * inc, const ulong * coeffs, slong length, const fq_nmod_ctx_t ctx)
void _n_fqp_zip_eval_step(ulong * res, ulong * cur, const ulong * inc, const ulong * coeffs, slong length, slong d, nmod_t mod)

One evaluation step: return (or set `res` to) $\sum_i c_i
\mathit{cur}_i$ and multiply each `cur[i]` by `inc[i]`. The `fq` version
works with $\mathbb{F}_q$ elements throughout, the `fqp` version with
`cur` and `inc` in the prime field and `coeffs` in $\mathbb{F}_q$.
:::

::: function
int _nmod_zip_vand_solve(ulong * coeffs, const ulong * monomials, slong mlength, const ulong * evals, slong elength, const ulong * master, ulong * scratch, nmod_t ctx)
int _n_fq_zip_vand_solve(ulong * coeffs, const ulong * monomials, slong mlength, const ulong * evals, slong elength, const ulong * master, ulong * scratch, const fq_nmod_ctx_t ctx)
int _n_fqp_zip_vand_solve(ulong * coeffs, const ulong * monomials, slong mlength, const ulong * evals, slong elength, const ulong * master, ulong * scratch, const fq_nmod_ctx_t ctx)

Solve for the `mlength` coefficients given the monomial values, the
`elength` $\ge$ `mlength` evaluations and the master polynomial, and check
the extra evaluations. Return 1 on success, 0 if the extra evaluations
are inconsistent, and $-1$ if the monomial values are not distinct.
:::

## Stacks of temporaries {#stacks}

A stack hands out initialised objects that keep their memory between
uses. `fit_request(S, k)` makes sure that $k$ more objects are available
and returns a pointer to them, `request(S, k)` does the same and marks them
as taken, `take_top` takes one object (which must have been made
available with `fit_request`), `give_back(S, k)` returns the last
$k$ taken objects, and `size` returns the number of objects taken.
Objects must be given back in the reverse order.

::: function
void n_poly_stack_init(n_poly_stack_t S)
void n_poly_stack_clear(n_poly_stack_t S)
n_poly_struct ** n_poly_stack_fit_request(n_poly_stack_t S, slong k)
n_poly_struct ** n_poly_stack_request(n_poly_stack_t S, slong k)
n_poly_struct * n_poly_stack_take_top(n_poly_stack_t S)
void n_poly_stack_give_back(n_poly_stack_t S, slong k)
slong n_poly_stack_size(const n_poly_stack_t S)

The stack of `n_poly_struct`.
:::

::: function
ulong * n_poly_stack_vec_init(n_poly_stack_t S, slong len)
void n_poly_stack_vec_clear(n_poly_stack_t S)

Take an `n_poly` from the stack and return its coefficient array with
space for `len` words, respectively give it back.
:::

::: function
void n_bpoly_stack_init(n_bpoly_stack_t S)
void n_bpoly_stack_clear(n_bpoly_stack_t S)
n_bpoly_struct ** n_bpoly_stack_fit_request(n_bpoly_stack_t S, slong k)
n_bpoly_struct ** n_bpoly_stack_request(n_bpoly_stack_t S, slong k)
n_bpoly_struct * n_bpoly_stack_take_top(n_bpoly_stack_t S)
void n_bpoly_stack_give_back(n_bpoly_stack_t S, slong k)
slong n_bpoly_stack_size(const n_bpoly_stack_t S)

The stack of `n_bpoly_struct`.
:::

::: function
void n_polyun_stack_init(n_polyun_stack_t S)
void n_polyun_stack_clear(n_polyun_stack_t S)
n_polyun_struct ** n_polyun_stack_fit_request(n_polyun_stack_t S, slong k)
n_polyun_struct ** n_polyun_stack_request(n_polyun_stack_t S, slong k)
n_polyun_struct * n_polyun_stack_take_top(n_polyun_stack_t S)
void n_polyun_stack_give_back(n_polyun_stack_t S, slong k)
slong n_polyun_stack_size(const n_polyun_stack_t S)

The stack of `n_polyun_struct`.
:::
