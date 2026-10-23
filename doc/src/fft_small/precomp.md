---
title: Products and divisions with a fixed operand
functions: mul_precomp_struct _mul_precomp_init _nmod_poly_mul_mid_precomp _mul_precomp_init_tr _nmod_poly_mul_tr_precomp nmod_poly_divrem_precomp_struct _nmod_poly_divrem_precomp_init _nmod_poly_divrem_precomp
---

A product costs, per prime, two forward transforms, a pointwise product
and one inverse transform. When one operand $b$ is multiplied by many
polynomials, its transforms can be computed once and kept, which saves
about a third of the work of every later product. For middle products
in which the fixed polynomial is the long operand, a transposed version
keeps its images through the transposed inverse transform instead.
`nmod_poly` uses both for evaluation, interpolation and extrapolation at
geometric progressions, where the same polynomial multiplies several
inputs. Division with remainder by a
fixed divisor builds on it: the divisor and its power series inverse are
both kept in transformed form.

## Products {#mul}

1. Call `_mul_precomp_init` with $b$ and a depth $d \ge 4$ such that
   $2^d$ is at least the lengths of the products (or windows) needed
   later. Pass `btrunc` $= 2^d$ unless all later outputs are shorter, and
   `len_bound` $= 2^d$ unless a smaller bound on the number of products
   per coefficient is known, such as $\min(a_n, b_n)$. Pass a nonzero
   `direct_len` (at least 100) to allow a transform modulo $n$ itself
   when $n$ is a suitable prime.
2. Call `_nmod_poly_mul_mid_precomp` for each $a$. It returns 0, writing
   nothing, if the window does not fit in a cyclic convolution of length
   $2^d$.
3. Free the data with `_mul_precomp_clear`.

A window $[z_l, z_h)$ of the product of lengths $a_n$ and $b_n$ fits if
$a_n + b_n - 1 \le 2^d$, or if $a_n \le 2^d$ and
$z_h \le 2^d \le a_n + b_n - 1 \le z_l + 2^d$: the coefficients of index
$i \ge 2^d$ then fold onto $i - 2^d < z_l$, outside the window. This is
the wrap-around trick for middle products [@BZ2011]: a middle product of
lengths $2m$ and $m$ fits in any transform length $2^d \ge 2m$, where the
full product would need $2^d \ge 3m - 1$.

`_nmod_poly_mul_mid_precomp` only reads `M`: it works on a copy of the
plan whose primes, Chinese remaindering data and scratch memory come from
the context `R`. Several threads can use one `M` at the same time, each
with its own context, normally `get_default_mpn_ctx()` of the thread. All
default contexts have the same primes, which is required.

## Transposed products {#transposed}

The transposed product $z_j = \sum_i r_i v_{i+j}$, $0 \le j < z_n$, is
computed by `_nmod_poly_mul_tr_mpn_ctx` as $F^T(\hat r \cdot I^T(v))$,
with $\hat r$ the forward transform of $r$ and $F^T$, $I^T$ the
transposes of the forward and inverse transforms
([](algorithms/polynomials/middle-product#transposed-fft)). For a fixed
$v$, the images $I^T(v)$ can be stored:

1. Call `_mul_precomp_init_tr` with $(v, v_n)$, a bound `len_bound` on
   the number of nonzero coefficients of the later $r$ (at most $v_n$),
   and `direct_len` as above. The plan is that of the window $[0, v_n)$
   of a product of length $v_n$, whose truncated transform length is
   $v_n$ rounded up to a multiple of 256 (for $v_n > 128$).
2. Call `_nmod_poly_mul_tr_precomp` for each $r$ with
   $r_n + z_n - 1 \le v_n$ (it returns 0 otherwise). Each call converts
   and transforms $r$, multiplies pointwise and applies $F^T$: two
   transforms instead of three.
3. Free the data with `_mul_precomp_clear`.

With $r$ the reversal of $a$, preceded by $k$ zeros, the output is the
window $[a_n + k - 1, a_n + k - 1 + z_n)$ of $av$. Storing the forward
transform of $v$ with `_mul_precomp_init` instead would need a cyclic
length $2^d$ into which that window fits; the transposed version needs
transforms of length about $v_n$ only. As for `_nmod_poly_mul_mid_precomp`, `M` is only read.

A polynomial $v$ of length 3000 prepared once, and the window
$[999, 2999)$ of its products with three polynomials of length 1000:

```c example name=fft_small_mul_tr_precomp requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    slong vn = 3000, an = 1000, zn = 2000, i, k;
    nmod_t mod;
    nn_ptr v, a, r, z, ref;
    mul_precomp_struct M[1];
    mpn_ctx_struct * R = get_default_mpn_ctx();

    nmod_init(&mod, 1000003);
    v = _nmod_vec_init(vn);
    a = _nmod_vec_init(an);
    r = _nmod_vec_init(an);
    z = _nmod_vec_init(zn + 2);
    ref = _nmod_vec_init(zn);
    for (i = 0; i < vn; i++)
        v[i] = nmod_pow_ui(3, (ulong) i * i, mod);

    /* store the images of v through the transposed inverse transform,
       for short operands of length at most an */
    if (!_mul_precomp_init_tr(M, v, vn, an, 0, mod, R))
        flint_abort();

    for (k = 1; k <= 3; k++)
    {
        for (i = 0; i < an; i++)
            a[i] = nmod_pow_ui(i + k, 3, mod);

        /* r = rev(a): the output is the window [an - 1, an - 1 + zn) of a*v */
        _nmod_poly_reverse(r, a, an, an);
        _nmod_poly_mulmid(ref, v, vn, a, an, an - 1, an - 1 + zn, mod);

        flint_printf("k = %wd: %s\n", k,
            (_nmod_poly_mul_tr_precomp(z, zn, r, an, M, mod, R)
              && _nmod_vec_equal(z, ref, zn)) ? "agrees" : "WRONG");
    }

    /* too long: an + zn - 1 > vn */
    flint_printf("window of length %wd returns %d\n", zn + 2,
        _nmod_poly_mul_tr_precomp(z, zn + 2, r, an, M, mod, R));

    _mul_precomp_clear(M);
    _nmod_vec_clear(v);
    _nmod_vec_clear(a);
    _nmod_vec_clear(r);
    _nmod_vec_clear(z);
    _nmod_vec_clear(ref);
    flint_cleanup_master();
    return 0;
}
```
```output
k = 1: agrees
k = 2: agrees
k = 3: agrees
window of length 2002 returns 0
```

## Division with remainder {#divrem}

For a divisor $b$ of length $b_n \ge 2$ with invertible leading
coefficient, `_nmod_poly_divrem_precomp_init` computes the power series
inverse $B$ of the reversal of $b$ to precision $B_n$ and stores the
transforms of $B$ and of $b$. A dividend $a$ of length
$a_n \le b_n + B_n - 1$ is then divided as follows, with
$q_n = a_n - b_n + 1 \le B_n$:

- the quotient $q$ is the window $[B_n - 1, B_n - 1 + q_n)$ of the
  product of the top $q_n$ coefficients of $a$ with $B$, which is the
  reversal of $\operatorname{rev}(a) \cdot \operatorname{rev}(b)^{-1}
  \bmod x^{q_n}$ [@GG2003, §9.1];
- the remainder is $r = a - bq$, computed modulo $x^N - 1$ for the power
  of two $N \ge \max(b_n - 1, 256)$. Since $\deg r < b_n - 1 \le N$, the
  reduction modulo $x^N - 1$ loses nothing, and only the residue of $a$
  modulo $x^N - 1$ and the cyclic product $bq$ are needed.

`_nmod_poly_divrem_mpn_ctx` performs the same steps without
precomputation.

`_nmod_poly_divrem_precomp` returns 0 without computing anything when
$q_n > B_n$, or when the quotient window does not fit in the prepared
transform.

## Examples {#examples}

A polynomial of length 1000 transformed once for length $2^{12}$ and
multiplied by three others, with full products of length 3999 and with a
middle product that uses the wrap-around:

```c example name=fft_small_mul_precomp requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    slong an = 4000, bn = 1000, i, k;
    nmod_t mod;
    nn_ptr a, b, full, z;
    mul_precomp_struct M[1];
    mpn_ctx_struct * R = get_default_mpn_ctx();

    nmod_init(&mod, 1000003);
    a = _nmod_vec_init(an);
    b = _nmod_vec_init(bn);
    full = _nmod_vec_init(an + bn - 1);
    z = _nmod_vec_init(an + bn - 1);
    for (i = 0; i < bn; i++)
        b[i] = nmod_pow_ui(2, i, mod);

    /* transform b once, for cyclic convolutions of length 2^12 */
    if (!_mul_precomp_init(M, b, bn, 4096, 12, 4096, 0, mod, R))
        flint_abort();

    for (k = 1; k <= 3; k++)
    {
        for (i = 0; i < an; i++)
            a[i] = nmod_pow_ui(i + k, 3, mod);

        /* a product of lengths 3000 and 1000 fits into 4096 coefficients */
        _nmod_poly_mul(full, a, 3000, b, bn, mod);
        flint_printf("k = %wd: full product %s", k,
            (_nmod_poly_mul_mid_precomp(z, 0, 3999, a, 3000, M, mod, R)
              && _nmod_vec_equal(z, full, 3999)) ? "agrees" : "WRONG");

        /* lengths 4000 and 1000: only windows that wrap below 999 fit */
        _nmod_poly_mul(full, a, an, b, bn, mod);
        flint_printf(", full product of length 4999 returns %d",
            _nmod_poly_mul_mid_precomp(z, 0, 4999, a, an, M, mod, R));
        flint_printf(", middle product %s\n",
            (_nmod_poly_mul_mid_precomp(z, 999, 4000, a, an, M, mod, R)
              && _nmod_vec_equal(z, full + 999, 3001)) ? "agrees" : "WRONG");
    }

    _mul_precomp_clear(M);
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
    _nmod_vec_clear(full);
    _nmod_vec_clear(z);
    flint_cleanup_master();
    return 0;
}
```
```output
k = 1: full product agrees, full product of length 4999 returns 0, middle product agrees
k = 2: full product agrees, full product of length 4999 returns 0, middle product agrees
k = 3: full product agrees, full product of length 4999 returns 0, middle product agrees
```

Division by a fixed divisor of length 500, prepared for quotients of
length up to 1500:

```c example name=fft_small_divrem_precomp requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    ulong n = 998244353;
    slong bn = 500, Bn = 1500, an, i;
    nmod_t mod;
    nmod_poly_t A, B, Q, R;
    nn_ptr q, r;
    nmod_poly_divrem_precomp_struct M[1];
    mpn_ctx_struct * ctx = get_default_mpn_ctx();

    nmod_init(&mod, n);
    nmod_poly_init(A, n);
    nmod_poly_init(B, n);
    nmod_poly_init(Q, n);
    nmod_poly_init(R, n);
    for (i = 0; i < bn; i++)
        nmod_poly_set_coeff_ui(B, i, nmod_pow_ui(7, i + 1, mod));

    /* the divisor serves dividends of length up to bn + Bn - 1 = 1999 */
    _nmod_poly_divrem_precomp_init(M, B->coeffs, bn, Bn, mod, ctx);
    q = _nmod_vec_init(Bn);
    r = _nmod_vec_init(bn - 1);

    for (an = 1000; an <= 1999; an += 333)
    {
        nmod_poly_zero(A);
        for (i = 0; i < an; i++)
            nmod_poly_set_coeff_ui(A, i, nmod_pow_ui(i + 2, 5, mod));
        nmod_poly_divrem(Q, R, A, B);

        if (!_nmod_poly_divrem_precomp(q, r, A->coeffs, an, M, mod, ctx))
            flint_printf("length %wd: transform too short\n", an);
        else
            flint_printf("length %wd: quotient %s, remainder %s\n", an,
                _nmod_vec_equal(q, Q->coeffs, an - bn + 1) ? "agrees" : "WRONG",
                (R->length <= bn - 1
                  && _nmod_vec_equal(r, R->coeffs, R->length)
                  && _nmod_vec_is_zero(r + R->length, bn - 1 - R->length))
                    ? "agrees" : "WRONG");
    }

    _nmod_poly_divrem_precomp_clear(M);
    _nmod_vec_clear(q);
    _nmod_vec_clear(r);
    nmod_poly_clear(A);
    nmod_poly_clear(B);
    nmod_poly_clear(Q);
    nmod_poly_clear(R);
    flint_cleanup_master();
    return 0;
}
```
```output
length 1000: quotient agrees, remainder agrees
length 1333: quotient agrees, remainder agrees
length 1666: quotient agrees, remainder agrees
length 1999: quotient agrees, remainder agrees
```
