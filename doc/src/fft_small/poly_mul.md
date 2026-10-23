---
title: Multiplying polynomials
functions: _nmod_poly_mul_mid_mpn_ctx _nmod_poly_mul_mid_default_mpn_ctx _nmod_poly_mul_tr_mpn_ctx _nmod_poly_mul_tr_default_mpn_ctx _fmpz_poly_mul_mid_mpn_ctx _fmpz_poly_mul_mid_default_mpn_ctx
---

These functions compute a window $[z_l, z_h)$ of the coefficients of a
product of two polynomials, over $\Z/n\Z$ for any word-size $n$ or over
$\Z$. With $z_l = 0$ and $z_h = a_n + b_n - 1$ this is the full product,
with $z_l = 0$ a low product, and with $z_l = b_n - 1$, $z_h = a_n$ the
middle product of Newton iteration. `nmod_poly_mul`, `nmod_poly_mullow`,
`nmod_poly_mulmid` and the corresponding `fmpz_poly` functions call the
`default` versions when their tuning tables say so, which for most
moduli is from lengths of about 100 to 1000 on (see [](fft_small/usage)).

## Contract {#contract}

- $a_n, b_n \ge 1$, in either order. Coefficients over $\Z/n\Z$ must be
  reduced; any $n \ge 1$ is allowed. The output must not overlap the
  inputs and receives $z_h - z_l$ coefficients; nothing is written if
  $z_l \ge z_h$, and coefficients from index $a_n + b_n - 1$ on are zero.
- The `fmpz` functions need initialised output coefficients. They return
  0 without writing anything if the coefficients are too large for eight
  primes, roughly when $\beta_a + \beta_b + \log_2 b_n \ge 397$ for
  coefficients of at most $\beta_a$ and $\beta_b$ bits; the caller then
  uses another algorithm, as `fmpz_poly_mul` does.
- Over $\Z/n\Z$ the bound on the number of products in a coefficient uses
  $\min(a_n, b_n)$ and nothing depends on the operand order. The `fmpz`
  functions use $b_n$, so a very unbalanced product may need one prime
  less with the shorter operand second.

## Number of primes {#primes}

A coefficient of the integer product of polynomials with entries in
$[0, n)$ is a sum of at most $b_n$ products less than $2^{2k}$, where $k$
is the bit length of $n$. The number of primes is the least whose product
$P$ is at least $b_n 2^{2k}$, which gives for products of two polynomials
of length $\ell$:

| bits of $n$ | $\ell = 10^2$ | $10^3$ | $10^4$ | $10^5$ | $10^6$ |
|---|---|---|---|---|---|
| 8 | 1 | 1 | 1 | 1 | 1 |
| 16 | 1 | 1 | 1 | 1 | 2 |
| 20 | 1 | 1 | 2 | 2 | 2 |
| 24 to 32 | 2 | 2 | 2 | 2 | 2 |
| 40 | 2 | 2 | 2 | 2 | 3 |
| 44 | 2 | 2 | 3 | 3 | 3 |
| 48 to 64 | 3 | 3 | 3 | 3 | 3 |

One prime suffices,
without Chinese remaindering, in two further cases:

- $n$ is one of the eight primes of the context, which have 50 bits;
- $n$ is a prime of at least 20 bits with $2^d \mid n - 1$, where $2^d$ is
  the transform length, and $b_n \ge 100$. The transform is then computed
  modulo $n$ itself, with twiddle tables set up for the call.

For `fmpz` polynomials the coefficients are bounded by
$b_n 2^{\beta_a + \beta_b}$ in absolute value and are reconstructed in
the symmetric range, which needs $P \ge b_n 2^{\beta_a + \beta_b + 1}$.

## Algorithm {#algorithm}

Each operand is reduced modulo each prime and transformed with a
truncated transform of length $2^d$, the transforms are multiplied
pointwise together with the normalisation factor, and the inverse
transform is applied. The truncation is the product length rounded up to
a multiple of 256, or the full length for transforms of 256 points or
less. If the window starts high enough that the coefficients beyond a
power of two $2^d \ge z_h$ fold onto coefficients below $z_l$, the
convolution is computed cyclically with length $2^d$. Finally each output
coefficient is reconstructed from its residues and reduced modulo $n$.
The transforms run in parallel over the primes when the second operand
has at least 1500 coefficients (for `fmpz`, 1000 with two or more primes
and 300 with four or more), and the reconstruction over ranges of
coefficients.

::: summary algorithms/fft/multimodular :::

## Transposed products {#transposed}

`_nmod_poly_mul_tr_mpn_ctx` computes
$z_j = \sum_{i < a_n} a_i v_{i+j}$ for $0 \le j < z_n$, the transpose of
the multiplication by $a$ of polynomials of length $z_n$. With $a$ the
reversal of $f$ this is the window $[a_n - 1, a_n - 1 + z_n)$ of $fv$,
the usual middle product, and every window $[z_l, z_h)$ of a product
reduces to this form after shifting and zero-padding the longer operand.
The transposed product costs one product of lengths $a_n$ and $z_n$: the
plan is that of that product, and the convolution engine replaces the
forward transform of $v$ by the transposed inverse transform and the
inverse transform by the transposed forward transform
([](fft_small/sd_fft#transposed)).

`_nmod_poly_mul_mid_mpn_ctx` computes the same window from the product of
$f$ and $v$, of length $a_n + v_n - 1$, which is up to $2 a_n + z_n - 2$.
It can use a cyclic convolution when a power of two lies between
$a_n + z_n - 1$ and that length
([](algorithms/polynomials/middle-product#lem-wrap)); otherwise its
transform is longer. For the balanced middle product, $a_n = z_n$, that
is a transform of length about $3 a_n$ against $2 a_n$.
`_nmod_poly_mulmid_fft_small`, and through it `nmod_poly_mulmid`,
compares the transform lengths of both plans and uses the transposed
product when it is strictly shorter. In timings of the window
$[n - 1, 2n - 1)$ of lengths $2n - 1$ and $n$ modulo 30- and 60-bit primes
($n$ from 300 to $10^5$), the transposed product was 1.2 to 1.9 times
faster than the window product, except where the latter wraps around at
the same transform length and the dispatcher keeps it.

The integer middle products of `_mpn_ctx_mpn_mul_range` (behind
`flint_mpn_mulmid`) and of `radix_mulmid_fft_small` choose between the
two forms in the same way ([](fft_small/mpn_mul#algorithm)). For a
fixed $v$ used with many $a$, `_mul_precomp_init_tr` stores the
transposed inverse transform of $v$ and `_nmod_poly_mul_tr_precomp`
computes each transposed product with two transforms
([](fft_small/precomp#transposed)).

## Example {#example}

A middle product modulo a 62-bit prime, and a product of integer
polynomials with large coefficients:

```c example name=fft_small_poly_mul requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    slong an = 2000, bn = 1000, i;
    nmod_t mod;
    nn_ptr a, b, full, mid;
    fmpz_poly_t f, g, h;
    fmpz * w;
    int ok;

    /* the middle product of lengths 2000 and 1000 modulo a 62-bit prime */
    nmod_init(&mod, UWORD(4611686018427387847));
    a = _nmod_vec_init(an);
    b = _nmod_vec_init(bn);
    full = _nmod_vec_init(an + bn - 1);
    mid = _nmod_vec_init(an - bn + 1);
    for (i = 0; i < an; i++)
        a[i] = nmod_pow_ui(3, i, mod);
    for (i = 0; i < bn; i++)
        b[i] = nmod_pow_ui(5, i * i, mod);

    _nmod_poly_mul(full, a, an, b, bn, mod);
    _nmod_poly_mul_mid_default_mpn_ctx(mid, bn - 1, an, a, an, b, bn, mod);
    flint_printf("middle product: %s\n",
        _nmod_vec_equal(mid, full + bn - 1, an - bn + 1) ? "agrees" : "WRONG");

    /* integer polynomials: 150-bit coefficients fit in eight primes ... */
    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(h);
    for (i = 0; i < 500; i++)
    {
        fmpz_poly_set_coeff_ui(f, i, i + 1);
        fmpz_poly_set_coeff_si(g, i, (i % 2) ? -1 : 1);
    }
    fmpz_poly_scalar_mul_2exp(f, f, 140);
    fmpz_poly_scalar_mul_2exp(g, g, 150);
    fmpz_poly_mul(h, f, g);

    w = _fmpz_vec_init(999);
    ok = _fmpz_poly_mul_mid_default_mpn_ctx(w, 0, 999, f->coeffs, 500, g->coeffs, 500);
    flint_printf("149 + 151 bits: returns %d, %s\n", ok,
        (ok && _fmpz_vec_equal(w, h->coeffs, 999)) ? "agrees" : "-");

    /* ... but 200-bit coefficients do not */
    fmpz_poly_scalar_mul_2exp(f, f, 60);
    fmpz_poly_scalar_mul_2exp(g, g, 50);
    ok = _fmpz_poly_mul_mid_default_mpn_ctx(w, 0, 999, f->coeffs, 500, g->coeffs, 500);
    flint_printf("209 + 201 bits: returns %d\n", ok);

    _fmpz_vec_clear(w, 999);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
    _nmod_vec_clear(full);
    _nmod_vec_clear(mid);
    flint_cleanup_master();
    return 0;
}
```
```output
middle product: agrees
149 + 151 bits: returns 1, agrees
209 + 201 bits: returns 0
```

The 62-bit modulus needs three primes. For the integer polynomials the
bound is $500 \cdot 2^{149 + 151 + 1} < 2^{310}$, which seven primes
(348 bits) cover; with 410 bits of coefficients even eight primes (398
bits) are too few.

The middle product of the example above as a transposed product, by the
reversal of $b$, and the same coefficients from `nmod_poly_mulmid`:

```c example name=fft_small_poly_mul_tr requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    slong an = 2000, bn = 1000, zn = an - bn + 1, i;
    nmod_t mod;
    nn_ptr a, b, rb, mid, tr;

    nmod_init(&mod, UWORD(4611686018427387847));
    a = _nmod_vec_init(an);
    b = _nmod_vec_init(bn);
    rb = _nmod_vec_init(bn);
    mid = _nmod_vec_init(zn);
    tr = _nmod_vec_init(zn);
    for (i = 0; i < an; i++)
        a[i] = nmod_pow_ui(3, i, mod);
    for (i = 0; i < bn; i++)
        b[i] = nmod_pow_ui(5, i * i, mod);

    /* coefficients [bn - 1, an) of a b */
    _nmod_poly_mulmid(mid, a, an, b, bn, bn - 1, an, mod);

    /* tr[j] = sum_i rb[i] a[i + j] with rb the reversal of b */
    _nmod_poly_reverse(rb, b, bn, bn);
    _nmod_poly_mul_tr_default_mpn_ctx(tr, zn, rb, bn, a, an, mod);
    flint_printf("transposed product: %s\n",
        _nmod_vec_equal(mid, tr, zn) ? "agrees" : "WRONG");

    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
    _nmod_vec_clear(rb);
    _nmod_vec_clear(mid);
    _nmod_vec_clear(tr);
    flint_cleanup_master();
    return 0;
}
```
```output
transposed product: agrees
```
