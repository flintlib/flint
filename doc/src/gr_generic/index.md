---
title: gr_generic.h
header: gr_generic.h
description: The default implementations of methods of the generic rings interface, and generic evaluation of integer polynomials.
category: Generic rings
order: 3
---

This module contains the default implementations that
`gr_method_tab_init` installs in every method table before a ring adds
its own methods ([](algorithms/generic/method-tables#tables)). A ring
needs to implement only a few methods; the functions here supply the
rest in terms of them: mixed-type arithmetic by conversion, `addmul` by a
temporary product, powers by binary exponentiation, vector operations by
loops over elementwise methods, and trivial cases with `GR_UNABLE`
elsewhere. They are documented so that implementers know what they get
for free and what they should override, and because a ring that
overrides a method can still call the default for cases it does not
handle. Users normally call the wrappers of [`gr`](gr) instead, which
dispatch to these functions through the method table.

The defaults never call a method that they implement, so they cannot
recurse indefinitely, but they assume that the required methods exist.
Their status codes follow the conventions of [`gr`](gr): a default that
cannot handle a case returns `GR_UNABLE`. The defaults of the special
functions are in [`gr_special`](gr_special), those of `gr_poly` and
`gr_mat` methods in their modules (`_gr_poly_mullow_generic`,
`gr_mat_mul_generic` and others).

## Required methods {#required}

::: macro
gr_generic_init
gr_generic_clear
gr_generic_swap
gr_generic_randtest
gr_generic_write
gr_generic_zero
gr_generic_one
gr_generic_equal
gr_generic_set
gr_generic_set_si
gr_generic_set_ui
gr_generic_set_fmpz
gr_generic_neg
gr_generic_add
gr_generic_sub
gr_generic_mul

The methods that every ring must implement have no real default: these
names are defined as `gr_not_implemented`, which returns `GR_UNABLE`
without doing anything. If `GR_GENERIC_DEBUG_RINGS` is set to 1 when
`gr_generic.h` is compiled, they are instead functions that throw an
exception naming the missing method.
:::

::: macro
GR_GENERIC_DEBUG_RINGS

Set to 0. The value 1, which makes the placeholders above throw, helps
when developing a new ring.
:::

## Context methods {#context}

The predicates `gr_generic_ctx_predicate`, `gr_generic_ctx_predicate_true`
and `gr_generic_ctx_predicate_false`, which return `T_UNKNOWN`, `T_TRUE`
and `T_FALSE`, are documented with the method tables of
[`gr`](gr#method-tables). `gr_generic_ctx_predicate` is the default of
most context predicates, and `GR_METHOD_CTX_IS_THREADSAFE` defaults to
`T_TRUE`. Some predicates have computed defaults: `gr_ctx_is_zero_ring`
is `T_FALSE` in an integral domain and otherwise tests whether $1 = 0$;
`gr_ctx_is_rational_vector_space` and its real and complex analogues are
`T_TRUE` in the zero ring, `T_FALSE` in other rings of finite
characteristic and otherwise unknown; `gr_ctx_is_approx_commutative_ring`
is `T_TRUE` for commutative rings.


::: function
void gr_generic_ctx_clear(gr_ctx_t ctx)

Does nothing: the default for rings whose context owns no memory.
:::

::: function
truth_t gr_generic_ctx_is_pretend_field(gr_ctx_t ctx)

The default of `GR_METHOD_CTX_IS_PRETEND_FIELD`: returns `gr_ctx_is_field`, since a ring without a pretence mechanism pretends only if it is a field.
:::

::: function
slong gr_generic_ctx_ngens_0(slong * ngens, gr_ctx_t ctx)
slong gr_generic_ctx_ngens_1(slong * ngens, gr_ctx_t ctx)

Set `*ngens` to 0 (the default of `gr_ctx_ngens`) or 1, and return `GR_SUCCESS`. Declared with return type `slong` but used as methods returning `int`.
:::

::: function
gr_ptr gr_generic_ctx_base(gr_ctx_t ctx)

Returns `NULL`: the default of `gr_ctx_base` for rings without a base ring.
:::

::: function
int gr_generic_ctx_fq_degree_prime_field(slong * res, gr_ctx_t ctx)
int gr_generic_ctx_fq_order_prime_field(fmpz_t res, gr_ctx_t ctx)

For rings $\Z/n\Z$, which are finite fields exactly when $n$ is prime: set the degree to 1, respectively the order to $p$ (obtained with `gr_ctx_fq_prime`), if `gr_ctx_is_field` returns `T_TRUE`, and return `GR_UNABLE` otherwise. Installed by the `nmod`, `fmpz_mod` and `mpn_mod` contexts.
:::

## Elements and output {#elements}


::: function
void gr_generic_set_shallow(gr_ptr res, gr_srcptr x, const gr_ctx_t ctx)

Copies the `sizeof_elem` bytes of `x` to `res` with `memcpy`.
:::

::: function
int gr_generic_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx)

Ignores the number of digits and calls `gr_write`.
:::

::: function
int gr_generic_randtest_not_zero(gr_ptr x, flint_rand_t state, gr_ctx_t ctx)

Tries `gr_randtest` five times until the element is known to be nonzero, then falls back to $\pm 1$; returns `GR_DOMAIN` if even $\pm 1$ is zero (the zero ring) and `GR_UNABLE` if this cannot be decided.
:::

::: function
int gr_generic_randtest_small(gr_ptr x, flint_rand_t state, gr_ctx_t ctx)

Sets `x` to $a g + b$ with random integers $a, b \in [-3, 3]$, where $g$ is the generator of the ring if `gr_gen` succeeds and 0 otherwise, divided with probability 1/4 by a random integer in $[1, 4]$ when the division succeeds; falls back to a small integer and then to zero.
:::

::: function
int gr_generic_gens(gr_vec_t vec, gr_ctx_t ctx)
int gr_generic_gens_single(gr_vec_t vec, gr_ctx_t ctx)
int gr_generic_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)

Set `vec` to the generators of the ring: the empty vector (the default of `gr_gens`), the single generator returned by `gr_gen` (for rings with one generator), or the result of `gr_gens` (the default of `gr_gens_recursive`, for rings without base rings).
:::

::: function
truth_t gr_generic_is_zero(gr_srcptr x, gr_ctx_t ctx)
truth_t gr_generic_is_one(gr_srcptr x, gr_ctx_t ctx)
truth_t gr_generic_is_neg_one(gr_srcptr x, gr_ctx_t ctx)

Compare `x` with a newly initialised element (which is zero), with the result of `gr_one`, or with the result of `gr_neg_one`, using `gr_equal`; `T_UNKNOWN` if the constant cannot be constructed.
:::

::: function
int gr_generic_neg_one(gr_ptr res, gr_ctx_t ctx)

Sets `res` to $-1$ by `gr_one` and `gr_neg`.
:::

## Parsing {#parsing}


::: macro
GR_PARSE_BALANCE_ADDITIONS
GR_PARSE_RING_EXPONENTS

Flags for `gr_generic_set_str_expr`: reorder long sums into a balanced
tree of additions, which is much faster for huge sums in rings where
elements grow (polynomials); and parse exponents as arbitrary
expressions in the ring instead of nonnegative integer literals.
:::


::: function
int gr_generic_set_str_expr(gr_ptr res, const char * s, int flags, gr_ctx_t ctx)
int gr_generic_set_str(gr_ptr res, const char * s, gr_ctx_t ctx)
int gr_generic_set_str_balance_additions(gr_ptr res, const char * s, gr_ctx_t ctx)
int gr_generic_set_str_ring_exponents(gr_ptr res, const char * s, gr_ctx_t ctx)

Parse the expression `s` and evaluate it in the ring. `gr_generic_set_str` is `gr_generic_set_str_expr` with no flags and the default of `gr_set_str`; the other two set one flag each.

::: details
The grammar has integer literals, decimal literals with an optional exponent (`1.25e-3`, evaluated with `gr_set_fmpz_10exp_fmpz`), the binary operators `+`, `-`, `*`, `/` and `^`, unary minus, parentheses, and:

* the generators of the ring and of its base rings, by the names given by `gr_gens_recursive` and `gr_ctx_gen_name`;
* the constants `pi`, `i`, `inf` and `nan` (also `Pi`, `I`, `Inf`, `NaN`), when the ring has them;
* the functions `abs`, `sgn`, `arg`, `inv`, `sqrt`, `rsqrt`, `floor`, `ceil`, `nint`, `re`, `im`, `conj`, `fac`, `log`, `exp`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinpi`, `cospi`, `tanpi`, `gamma` and `zeta` of one argument;
* the ball notation `m +/- r`, evaluated with `gr_set_interval_mid_rad`, and `m +/- inf` with `gr_set_interval_mid_inf`, which does not require `inf` to be an element;
* `O(x^n)`, evaluated with `gr_big_o_base_fmpz` without forming $x^n$.

Exponents are nonnegative integer literals unless `GR_PARSE_RING_EXPONENTS` is set. The status of every operation is propagated, so `"1/2"` in $\Z$ returns `GR_DOMAIN`; a syntax error returns `GR_UNABLE`.
:::
:::

::: function
int gr_generic_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t expr, gr_ctx_t ctx)

Evaluates the symbolic expression `expr` in the ring, using the vectors `inputs` and `outputs` as a table of symbols and their values (the default of `gr_set_fexpr`). Declared only if `fexpr.h` is included first.
:::

::: function
int gr_generic_set_interval_mid_inf(gr_ptr res, gr_srcptr m, gr_ctx_t ctx)

Sets `res` to $m \pm \infty$ with `gr_pos_inf` and `gr_set_interval_mid_rad`; returns `GR_UNABLE` if the ring has no infinity.
:::

## Conversions {#conversions}


::: function
int gr_generic_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t xctx, gr_ctx_t ctx)

The default conversion from another ring: `gr_set` if `xctx` is `ctx`, `gr_set_fmpz` from $\Z$, `gr_set_fmpq` from $\mathbb{Q}$, `gr_set_fexpr` from symbolic expressions, and `GR_UNABLE` from any other ring. Rings override it to accept their base rings and related types.
:::

::: function
int gr_generic_set_fmpq(gr_ptr res, const fmpq_t y, gr_ctx_t ctx)

Sets `res` to $p/q$ by converting $p$ and $q$ with `gr_set_fmpz`, inverting $q$ and multiplying; returns `GR_DOMAIN` if $q$ is not invertible in the ring.
:::

::: function
int gr_generic_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_ptr x, gr_ctx_t ctx)

Set `res` to $x \cdot 2^y$ by powering 2 in the ring, respectively write `x` as $m \cdot 2^e$ with $m$ odd by converting it to a rational number with `gr_get_fmpq` (failing with `GR_DOMAIN` if its denominator is not a power of 2).
:::

::: function
int gr_generic_set_fmpz_10exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx)

Sets `res` to $x \cdot 10^y$, by multiplying or dividing by a power of 10 computed in the ring.
:::

## Arithmetic {#arithmetic}

The type variants convert the operand to the ring and call the plain
operation: `_fmpz` and `_other` variants with `gr_set_fmpz` and
`gr_set_other`, and the `_ui` and `_si` variants through an `fmpz`
(`gr_mul_ui` calls `gr_mul_fmpz`, for example). The `_fmpq` variants of
addition and multiplication convert the rational number with
`gr_set_fmpq`. A conversion failure is returned without performing the
operation.


::: function
int gr_generic_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_generic_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Defaults of the type variants of `gr_add`.
:::

::: function
int gr_generic_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_generic_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Defaults of the type variants of `gr_sub`.
:::

::: function
int gr_generic_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_generic_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Defaults of the type variants of `gr_mul`.
:::

::: function
int gr_generic_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_generic_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_generic_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

Compute the product $x y$ in a temporary with the corresponding multiplication and add it to, or subtract it from, `res`.
:::

::: function
int gr_generic_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Compute $2x$ as $x + x$ and $x^2$ as $x \cdot x$.
:::

::: function
int gr_generic_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)

Multiply by $2^y$ computed in the ring with `gr_pow_ui` (or `gr_pow_fmpz`), or divide by $2^{-y}$ for negative $y$.
:::

## Division {#division}


::: function
int gr_generic_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
truth_t gr_generic_is_invertible(gr_srcptr x, gr_ctx_t ctx)

Handle $x = 1$ and $x = -1$ (invertible, their own inverses) and $x = 0$ (not invertible, `GR_DOMAIN`), and return `GR_UNABLE` or `T_UNKNOWN` for every other element.
:::

::: function
int gr_generic_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_generic_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_generic_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_generic_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Convert the operand to the ring and call `gr_div`, returning the conversion status if it fails.
:::

::: function
int gr_generic_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

Multiplies by $1/y$ with `gr_mul_fmpq`; for $y = 0$ it converts 0 and calls `gr_div`, so that the ring decides (the zero ring allows it).
:::

::: function
int gr_generic_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Calls `gr_div`. The type variants of `gr_divexact` default to those of `gr_div`.
:::

## Powers {#powering}

These functions implement the default integer powers. The algorithms
are described in [](algorithms/generic/powering).


::: function
int gr_generic_pow_ui_binexp(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
int gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx)

Compute $x^e$ for $e \ge 0$ by left-to-right binary exponentiation, with special code for $e \le 4$; `gr_generic_pow_fmpz_binexp` returns `GR_DOMAIN` for negative $e$. The output may alias the input.
:::

::: function
int gr_generic_pow_ui_sliding(gr_ptr f, gr_srcptr g, ulong pow, gr_ctx_t ctx)
int gr_generic_pow_fmpz_sliding(gr_ptr f, gr_srcptr g, const fmpz_t pow, gr_ctx_t ctx)

Compute $x^e$ for $e \ge 0$ by the sliding window method with a window size $k$ between 1 and 11 bits chosen from the bit length of $e$, using fewer multiplications than binary exponentiation for large exponents at the cost of a table of $2^{k-1}$ odd powers. `gr_generic_pow_fmpz_sliding` returns `GR_DOMAIN` for negative $e$. The output may alias the input.
:::

::: function
int gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
int gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
int gr_generic_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t e, gr_ctx_t ctx)

The defaults of `gr_pow_ui`, `gr_pow_si` and `gr_pow_fmpz`: binary exponentiation; for a negative exponent, invert $x$ with `gr_inv` first, so that the result fails exactly when $x$ is not invertible.
:::

::: function
int gr_generic_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

The default of `gr_pow_fmpq`: an integer exponent $p/1$ uses `gr_pow_fmpz`; an exponent $p/2$ uses `gr_sqrt` (or `gr_rsqrt` for $p < 0$) followed by an integer power; any other exponent is converted to the ring and passed to `gr_pow`. Failures of the square root are reported as `GR_UNABLE`, since the result may still exist for another choice of root.
:::

::: function
int gr_generic_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Defaults of `gr_pow_other` and `gr_other_pow`: an exponent in $\Z$ or $\mathbb{Q}$ is passed to `gr_pow_fmpz` or `gr_pow_fmpq`; otherwise the operand is converted with `gr_set_other` and `gr_pow` is called.
:::

## Roots, fractions and comparisons {#other}


::: function
truth_t gr_generic_is_square(gr_srcptr x, gr_ctx_t ctx)
int gr_generic_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Handle $0$ and $1$ (and for `gr_generic_rsqrt`, compute $1/\sqrt{x}$ from `gr_sqrt` and `gr_inv` when both succeed; $1/\sqrt{0}$ is `GR_DOMAIN`), and return `GR_UNABLE` or `T_UNKNOWN` otherwise.
:::

::: function
int gr_generic_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Return $x$ as its own numerator and 1 as its denominator.
:::

::: function
int gr_generic_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_generic_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_generic_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

`gr_generic_cmp` returns `GR_UNABLE`. The others compare absolute values computed with `gr_abs` with `gr_cmp`, or convert `y` with `gr_set_other` and call `gr_cmp` or `gr_cmpabs`.
:::

## Number sequences {#sequences}


::: function
int gr_generic_bernoulli_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_bernoulli_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Defaults of `gr_bernoulli_ui` and the related functions: in $\mathbb{Q}$, the exact Bernoulli numbers from `bernoulli_fmpq_ui` and the global cache of `bernoulli`; in rings with a real precision, numerical values computed with `arb_bernoulli_ui` (or `arb_bernoulli_fmpz`) and converted; otherwise exact rational values converted with `gr_set_fmpq`, which fails if a denominator is not invertible. `gr_generic_bernoulli_fmpz` returns `GR_UNABLE` for huge $n$ in exact rings.
:::

::: function
int gr_generic_eulernum_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
int gr_generic_eulernum_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
int gr_generic_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx)

Defaults for the Euler numbers $E_n$: exact integers from `arb_fmpz_euler_number_ui`, converted to the ring, or numerical values computed with ball arithmetic in rings with a real precision.
:::

::: function
int gr_generic_stirling_s1u_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_generic_stirling_s1_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_generic_stirling_s2_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
int gr_generic_stirling_s1u_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
int gr_generic_stirling_s1_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
int gr_generic_stirling_s2_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)

Defaults for the unsigned and signed Stirling numbers of the first kind and the Stirling numbers of the second kind $S(x, y)$, and for the rows $S(x, 0), \ldots, S(x, \mathit{len} - 1)$: computed as integers with the functions of `arith` and converted to the ring.
:::

## Vector methods {#vectors}

The defaults of the vector methods of [`gr_vec`](gr_vec) loop over the
corresponding elementwise methods, which they look up once in the method
table before the loop. They return the OR of the elementwise statuses.


::: function
void gr_generic_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)
void gr_generic_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)
void gr_generic_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)
int gr_generic_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
int gr_generic_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

Defaults of `_gr_vec_init`, `_gr_vec_clear`, `_gr_vec_swap`, `_gr_vec_zero`, `_gr_vec_set` and `_gr_vec_neg`.
:::

::: function
int gr_generic_vec_normalise(slong * res, gr_srcptr vec, slong len, gr_ctx_t ctx)
slong gr_generic_vec_normalise_weak(gr_srcptr vec, slong len, gr_ctx_t ctx)

Defaults of `_gr_vec_normalise` and `_gr_vec_normalise_weak`: scan from the end with `gr_is_zero`.
:::

::: function
int gr_generic_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int gr_generic_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int gr_generic_vec_mul(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int gr_generic_vec_div(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int gr_generic_vec_divexact(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int gr_generic_vec_pow(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)

Defaults of the elementwise binary operations `_gr_vec_add` and the others.
:::

::: function
int gr_generic_vec_add_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_sub_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_mul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_div_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_divexact_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_pow_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_add_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_sub_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_mul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_div_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_divexact_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_pow_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_add_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int gr_generic_vec_sub_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int gr_generic_vec_mul_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int gr_generic_vec_div_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int gr_generic_vec_divexact_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int gr_generic_vec_pow_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int gr_generic_vec_add_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int gr_generic_vec_sub_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int gr_generic_vec_mul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int gr_generic_vec_div_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int gr_generic_vec_divexact_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int gr_generic_vec_pow_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int gr_generic_vec_add_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int gr_generic_vec_sub_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int gr_generic_vec_mul_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int gr_generic_vec_div_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int gr_generic_vec_divexact_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int gr_generic_vec_pow_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int gr_generic_scalar_add_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_sub_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_mul_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_div_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_pow_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Defaults of the operations between a vector and a scalar.
:::

::: function
int gr_generic_vec_add_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int gr_generic_vec_sub_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int gr_generic_vec_mul_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int gr_generic_vec_div_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int gr_generic_vec_divexact_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int gr_generic_vec_pow_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int gr_generic_other_add_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int gr_generic_other_sub_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int gr_generic_other_mul_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int gr_generic_other_div_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int gr_generic_other_divexact_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int gr_generic_other_pow_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int gr_generic_vec_add_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int gr_generic_vec_sub_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int gr_generic_vec_mul_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int gr_generic_vec_div_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int gr_generic_vec_divexact_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int gr_generic_vec_pow_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int gr_generic_scalar_other_add_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_other_sub_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_other_mul_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_other_div_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_other_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_scalar_other_pow_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Defaults of the operations with an operand in another ring.
:::

::: function
int gr_generic_vec_mul_scalar_2exp_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

Default of `_gr_vec_mul_scalar_2exp_si`.
:::

::: function
int gr_generic_vec_addmul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_submul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int gr_generic_vec_addmul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_submul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int gr_generic_vec_addmul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)

Defaults of `_gr_vec_addmul_scalar` and the related functions: multiply into a temporary and add or subtract.
:::

::: function
truth_t gr_generic_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)

Defaults of `_gr_vec_equal` and `_gr_vec_is_zero`: return `T_FALSE` at the first entry known to differ, otherwise `T_UNKNOWN` if some comparison was unknown and `T_TRUE` if not. `gr_generic_vec_is_zero` is declared as returning `int` but returns a `truth_t`.
:::

::: function
int gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int gr_generic_vec_dot_strided(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, slong stride1, gr_srcptr vec2, slong stride2, slong len, gr_ctx_t ctx)
int gr_generic_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
int gr_generic_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx)
int gr_generic_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)

Defaults of the dot products ([](gr_vec/dot)): multiply and add term by term in a temporary. The strided version gathers the entries by shallow copies and calls `_gr_vec_dot` when the ring has its own `VEC_DOT` method.
:::

::: function
int gr_generic_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx)

Defaults of `_gr_vec_set_powers` and `_gr_vec_reciprocals`.
:::

::: function
int gr_generic_vec_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, slong len, gr_ctx_t ctx)
int gr_generic_vec_gather(gr_ptr res, gr_srcptr vec, const slong * idx, slong len, gr_ctx_t ctx)
int gr_generic_vec_scatter(gr_ptr vec, const slong * idx, gr_srcptr src, slong len, gr_ctx_t ctx)

Defaults of `_gr_vec_set_other`, `_gr_vec_gather` and `_gr_vec_scatter`: loops over `gr_set_other` and `gr_set`.
:::

::: function
int gr_generic_vec_sqrt(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_rsqrt(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_exp(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_log(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_sin(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_cos(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int gr_generic_vec_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, slong len, gr_ctx_t ctx)

Defaults of the elementwise elementary functions of `gr_vec`: loops over `gr_sqrt`, `gr_exp` and the others.
:::

::: function
int gr_generic_get_interval_mid_rad(gr_ptr m, gr_ptr r, gr_srcptr x, gr_ctx_t ctx)
int gr_generic_vec_get_interval_mid_rad(gr_ptr m, gr_ptr r, gr_srcptr x, slong len, gr_ctx_t ctx)

Defaults of `gr_get_interval_mid_rad` (sets $m = x$ and $r = 0$, for exact rings) and of its vector version (loops over the method).
:::

## Polynomial methods {#polynomials}

The defaults of the methods `GR_METHOD_POLY_SET_LENGTH_NORMALISE`
(`gr_generic_poly_set_length_normalise`), `GR_METHOD_POLY_EVALUATE`
(`gr_generic_poly_evaluate`, Horner's rule), `GR_METHOD_POLY_PREINV_SET`
and `GR_METHOD_POLY_DIVREM_PREINV` (`gr_generic_poly_preinv_set` and
`gr_generic_poly_divrem_preinv`, division by a precomputed modulus) are
documented in [`gr_poly`](gr_poly).


::: function
int gr_generic_poly_factor(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, gr_srcptr elt, int flags, gr_ctx_t ctx)
int gr_generic_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)

Defaults of the methods `GR_METHOD_POLY_FACTOR` and `GR_METHOD_POLY_ROOTS` (factorisation and root finding of polynomials over the ring): implemented generically for finite fields (`gr_ctx_is_field` and `gr_ctx_is_finite` both `T_TRUE`), and `GR_UNABLE` for other rings.
:::

::: function
int gr_generic_poly_factor_roots(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, gr_srcptr elt, int flags, gr_ctx_t ctx)

Factors a polynomial over an integral domain that splits into linear factors, by computing its roots with `gr_poly_roots`: sets `c` to the leading coefficient (as a constant polynomial) and the factors to $x - r$ with the multiplicities of the roots. Returns `GR_DOMAIN` for the zero polynomial and `GR_UNABLE` if the leading coefficient is not known to be nonzero or if the multiplicities do not add up to the degree (not all roots found). It is the polynomial factorisation method of complex balls (`acb`) and complex floating-point numbers (`acf`), and is called by the `qqbar` and `ca` contexts for polynomials over algebraically closed fields.
:::

## Evaluating integer polynomials {#evaluation}

These functions evaluate polynomials with integer coefficients at an
element of any ring; they are not methods and cannot be overridden.


::: function
int _gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
int _gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
int _gr_fmpz_poly_evaluate(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
int gr_fmpz_poly_evaluate(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $f(x)$ for the integer polynomial $f$ (as an `fmpz` array of length `len`, or an `fmpz_poly_t`): by Horner's rule, by rectangular splitting, or (`_gr_fmpz_poly_evaluate`) by Horner's rule for length at most 6 and rectangular splitting otherwise. Declared only if `fmpz_poly.h` is included first.

::: details
Rectangular splitting [@PS1973] computes the powers $x^2, \ldots, x^m$ with $m \approx \sqrt{n}$ once and evaluates $f$ as a polynomial in $x^m$ whose coefficients are linear combinations of the precomputed powers with integer coefficients; this uses about $2\sqrt{n}$ multiplications in the ring instead of $n$, and the multiplications by the integer coefficients are cheap scalar operations (`_gr_vec_dot_fmpz`). It is preferable in rings where a multiplication is much more expensive than a multiplication by an integer, such as balls at high precision or matrices.
:::
:::

::: function
int gr_fmpz_mpoly_evaluate_iter(gr_ptr res, const fmpz_mpoly_t pol, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)
int gr_fmpz_mpoly_evaluate_horner(gr_ptr res, const fmpz_mpoly_t pol, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)
int gr_fmpz_mpoly_evaluate(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)

Set `res` to the value of the multivariate integer polynomial `f` (with context `mctx`) at the vector of values `x`, one entry for each variable: term by term, or by a multivariate Horner scheme that chooses the variable to factor out greedily; `gr_fmpz_mpoly_evaluate` uses the first method for polynomials with at most one term and word-size exponents and the second otherwise. Declared only if `fmpz_mpoly.h` is included first.
:::
