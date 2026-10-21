---
title: Arithmetic in transform space
summary: How to transform operands once, add several products pointwise and convert back once, with plans, operands, bounds and exports, for polynomials over $\Z/n\Z$ and for integers.
order: 2
---

A product computed by `fft_small` costs, per prime, two forward
transforms, a pointwise product, an inverse transform and a Chinese
remainder reconstruction. In a sum of products such as
$a_0 b_0 + a_1 b_1$, a dot product of polynomials or a matrix product,
the products can be added in transform space, so that only one inverse
transform and one reconstruction are needed for the whole sum, and each
operand that occurs in several products is transformed only once. The
plan and operand functions of this module expose this. They are what the
transformed rings of `gr` and the complex products of `mpn_extras` are
built on.

## Plans {#plans}

A plan fixes the primes, the transform length and the window of output
coefficients. Every operand transformed for a plan can be combined with
every other one of the same plan. The plan must be chosen for the largest
operands and for the largest number of products that will be added:

- `fft_small_plan_init_nmod` for polynomials over $\Z/n\Z$. A
  coefficient of a sum of $K$ products of polynomials of lengths $a_n$
  and $b_n$ is a sum of at most $K \min(a_n, b_n)$ products of residues,
  each less than $n^2 < 2^{2k}$ for $n$ of $k$ bits; pass these as
  `len_bound` and `prod_bits`. If operands are added in transform space
  before they are multiplied, their coefficients grow, which is accounted
  for by a larger `prod_bits` (one more bit for a sum of two operands).
- `fft_small_plan_init_nmod_cyclic` for cyclic convolutions modulo
  $x^{2^d} - 1$.
- `fft_small_plan_init_mpn` for integers, with `len_bound` the number of
  products added. Pass `is_signed` = 1 if products are subtracted.
- `fft_small_plan_init_mpn_window` if only some limbs of the result are
  needed, and `fft_small_plan_init_mpn_cyclic` for products modulo
  $B^{n} - 1$.

The plan's number of primes is the least whose product $P$ exceeds the
bound. The result is exact only if every true coefficient of the sum lies
in $[0, P)$ for the unsigned exports, or in $(-P/2, P/2]$ for the signed
ones; the bound guarantees this for sums of products of operands in the
range the plan was made for.

## Operands and their domains {#domains}

An operand created with `fft_small_op_init` holds one array of transform
values per prime. After `fft_small_fft_nmod` or `fft_small_fft_mpn` it is
in the *primal* domain. A pointwise product `fft_small_op_mul` or
`fft_small_op_sqr` puts its result in the *product* domain, which carries
the normalisation factor $(c_i 2^d)^{-1}$ of the plan: the inverse of
$2^d$ that the inverse transform needs, and the inverse of the Chinese
remainder cofactor $c_i$. The rules are:

- `mul`, `sqr`, `addmul`, `submul`, `addsqr` and `subsqr` take primal
  operands; the accumulator `Z` of the last four must be in the product
  domain.
- `add`, `sub` and `neg` take operands of the same domain and keep it.
  Primal sums are sums of the inputs; product sums are sums of products.
- A product cannot be multiplied again, since it would carry the
  normalisation twice. Every intermediate value stays in the range that
  the transforms accept, whatever sequence of these operations is
  applied.
- `fft_small_ifft` transforms a product back in place. The result is then
  converted out with `fft_small_export_nmod`, `fft_small_export_mpn` or a
  signed export.

The exports of integer results may destroy the transform data: the
two-prime reconstruction works in place. Copy an operand first if it is
needed again. For integer plans, `fft_small_plan_mpn_nslots` gives the
number of chunk coefficients (slots) that the inverse transform computes,
which is what the signed exports read, and
`fft_small_export_mpn_signed_limbs` the number of limbs that a signed
reconstruction of that many slots writes.

Operands can live in caller-provided memory (`fft_small_op_init_borrowed`)
and in pools (`fft_small_op_pool_fit`), which drivers use to avoid one
allocation per operand and to keep the memory warm.

## Example: a sum of two polynomial products {#example-nmod}

Here $a_0 b_0 + a_1 b_1$ for polynomials of length 3000 modulo a 62-bit
prime costs four forward transforms, two pointwise products and one
inverse transform per prime, against four, two and two for two separate
products, and one reconstruction instead of two.

```c example name=fft_small_op_nmod requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    slong len = 3000, zn = 2 * len - 1, i, k;
    nmod_t mod;
    nn_ptr a[2], b[2], t, ref, z;
    fft_small_plan_t P;
    fft_small_op_t A, B, Z;
    mpn_ctx_struct * R = get_default_mpn_ctx();

    nmod_init(&mod, UWORD(4611686018427387847));       /* 62 bits */
    t = _nmod_vec_init(zn);
    ref = _nmod_vec_init(zn);
    z = _nmod_vec_init(zn);
    _nmod_vec_zero(ref, zn);
    for (k = 0; k < 2; k++)
    {
        a[k] = _nmod_vec_init(len);
        b[k] = _nmod_vec_init(len);
        for (i = 0; i < len; i++)
        {
            a[k][i] = nmod_pow_ui(i + 3 + k, 11, mod);
            b[k][i] = nmod_pow_ui(i + 5 + k, 13, mod);
        }
        _nmod_poly_mul(t, a[k], len, b[k], len, mod);
        _nmod_vec_add(ref, ref, t, zn, mod);
    }

    /* the full product (window [0, zn)) of sums of two products: up to
       2 len products of (n - 1)^2 < 2^124 accumulate in a coefficient */
    if (!fft_small_plan_init_nmod(P, R, 0, zn, zn, n_round_up(len, BLK_SZ),
                                  2 * len, 2 * NMOD_BITS(mod), mod, 0))
        flint_abort();
    flint_printf("%wu primes, transform length 2^%wu\n", P->np, P->depth);

    fft_small_op_init(A, P);
    fft_small_op_init(B, P);
    fft_small_op_init(Z, P);

    for (k = 0; k < 2; k++)
    {
        fft_small_fft_nmod(A, a[k], len, n_round_up(len, BLK_SZ), mod, P);
        fft_small_fft_nmod(B, b[k], len, n_round_up(len, BLK_SZ), mod, P);
        if (k == 0)
            fft_small_op_mul(Z, A, B, P);
        else
            fft_small_op_addmul(Z, A, B, P);
    }

    fft_small_ifft(Z, P);
    fft_small_export_nmod(z, Z, mod, P);
    flint_printf("a0 b0 + a1 b1: %s\n", _nmod_vec_equal(z, ref, zn) ? "agrees" : "WRONG");

    fft_small_op_clear(A);
    fft_small_op_clear(B);
    fft_small_op_clear(Z);
    fft_small_plan_clear(P);
    for (k = 0; k < 2; k++)
    {
        _nmod_vec_clear(a[k]);
        _nmod_vec_clear(b[k]);
    }
    _nmod_vec_clear(t);
    _nmod_vec_clear(ref);
    _nmod_vec_clear(z);
    flint_cleanup_master();
    return 0;
}
```
```output
3 primes, transform length 2^13
a0 b0 + a1 b1: agrees
```

## Example: a signed sum of integer products {#example-mpn}

The difference $ab - cd$ of products of integers of about 2000 limbs,
with a signed plan and a signed export. For sums of at most two products
of operands of up to 40000 limbs the plan uses the two-prime geometry.

```c example name=fft_small_op_mpn_signed requires=fft_small
#include <gmp.h>
#include <flint/flint.h>
#include <flint/mpn_extras.h>
#include <flint/fft_small.h>

int main(void)
{
    mpz_t a, b, c, d, x;
    slong n, zn, nslots, xn;
    int sign;
    nn_ptr z;
    fft_small_plan_t P;
    fft_small_op_t A, B, Z;
    mpn_ctx_struct * R = get_default_mpn_ctx();

    mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d); mpz_init(x);
    mpz_ui_pow_ui(a, 3, 80000);
    mpz_ui_pow_ui(b, 5, 50000);
    mpz_ui_pow_ui(c, 7, 50000);
    mpz_ui_pow_ui(d, 11, 40000);
    n = FLINT_MAX(FLINT_MAX(mpz_size(a), mpz_size(b)),
                  FLINT_MAX(mpz_size(c), mpz_size(d)));

    /* x = a b - c d: two products, signed */
    if (!fft_small_plan_init_mpn(P, R, n, n, 2, 1))
        flint_abort();
    flint_printf("%wd limbs: %wu primes, %wu-bit chunks\n", n, P->np, P->bits);

    fft_small_op_init(A, P);
    fft_small_op_init(B, P);
    fft_small_op_init(Z, P);

    fft_small_fft_mpn(A, mpz_limbs_read(a), mpz_size(a), P);
    fft_small_fft_mpn(B, mpz_limbs_read(b), mpz_size(b), P);
    fft_small_op_mul(Z, A, B, P);
    fft_small_fft_mpn(A, mpz_limbs_read(c), mpz_size(c), P);
    fft_small_fft_mpn(B, mpz_limbs_read(d), mpz_size(d), P);
    fft_small_op_submul(Z, A, B, P);
    fft_small_ifft(Z, P);

    nslots = fft_small_plan_mpn_nslots(P, n, n);
    zn = fft_small_export_mpn_signed_limbs(P, nslots);
    z = flint_malloc(zn * sizeof(ulong));
    fft_small_export_mpn_signed(z, zn, &sign, Z, nslots, P);

    mpz_mul(x, a, b);
    mpz_submul(x, c, d);
    xn = mpz_size(x);
    flint_printf("sign %s, %s\n", sign ? "negative" : "positive",
        (sign == (mpz_sgn(x) < 0) && xn <= zn
          && mpn_cmp(z, mpz_limbs_read(x), xn) == 0
          && flint_mpn_zero_p(z + xn, zn - xn)) ? "agrees" : "WRONG");

    flint_free(z);
    fft_small_op_clear(A);
    fft_small_op_clear(B);
    fft_small_op_clear(Z);
    fft_small_plan_clear(P);
    mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d); mpz_clear(x);
    flint_cleanup_master();
    return 0;
}
```
```output
2194 limbs: 2 primes, 42-bit chunks
sign negative, agrees
```

With the unsigned export `fft_small_export_mpn`, a difference is correct
only if it is nonnegative, since the reconstruction returns values in
$[0, P)$.

## Transposed products {#transposed}

The transpose of the multiplication by a fixed polynomial $a$ of length
$a_n$, applied to $v$ of length at most $a_n + z_n - 1$, is
$z_j = \sum_i a_i v_{i+j}$ for $j < z_n$; with $a$ reversed it is a middle
product ([](algorithms/polynomials/middle-product#transpose)). In
transform space it is computed with the plan of the product of lengths
$a_n$ and $z_n$ (the window $[0, z_n)$ of a convolution of length
$a_n + z_n - 1$) by transposing the two transforms and keeping the
pointwise product:

- `fft_small_fft_nmod` transforms $a$ as for a product;
- `fft_small_ifft_tr_nmod` converts $v$ and applies the transposed
  inverse transform, which gives a primal operand;
- `fft_small_op_mul` multiplies them, with the normalisation of the plan;
- `fft_small_fft_tr` applies the transposed forward transform truncated
  to the output length, which leaves the result in the first slots;
- `fft_small_export_nmod` (or `fft_small_export_nmod_range`) reconstructs
  it.

The transform of $a$ can be reused for any number of vectors $v$, and
transposed products of the same plan can be added in transform space
before `fft_small_fft_tr`, since the transposed transform is linear.
`_nmod_poly_mul_tr_mpn_ctx` performs the same steps for a single product.

```c example name=fft_small_op_transposed requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fft_small.h>

int main(void)
{
    slong an = 1000, zn = 1000, vn = an + zn - 1, i;
    ulong atrunc, ztrunc;
    nmod_t mod;
    nn_ptr f, rf, v, z, ref;
    fft_small_plan_t P;
    fft_small_op_t B, X;
    mpn_ctx_struct * R = get_default_mpn_ctx();

    nmod_init(&mod, UWORD(4611686018427387847));       /* 62 bits */
    f = _nmod_vec_init(an);
    rf = _nmod_vec_init(an);
    v = _nmod_vec_init(vn);
    z = _nmod_vec_init(zn);
    ref = _nmod_vec_init(zn);
    for (i = 0; i < an; i++)
        f[i] = nmod_pow_ui(i + 3, 11, mod);
    for (i = 0; i < vn; i++)
        v[i] = nmod_pow_ui(i + 5, 13, mod);

    /* the plan of the product of lengths an and zn: window [0, zn) of a
       convolution of length an + zn - 1, at most an products of
       (n - 1)^2 < 2^124 accumulating in a coefficient */
    if (!fft_small_plan_init_nmod(P, R, 0, zn, an + zn - 1,
            n_max(_fft_small_len_trunc(an), _fft_small_len_trunc(zn)),
            an, 2 * NMOD_BITS(mod), mod, 0))
        flint_abort();
    atrunc = n_min(_fft_small_op_trunc(an, P), P->ztrunc);
    ztrunc = n_min(_fft_small_op_trunc(zn, P), P->ztrunc);
    flint_printf("%wu primes, transform length 2^%wu, truncation %wu\n",
                 P->np, P->depth, P->ztrunc);

    fft_small_op_init(B, P);
    fft_small_op_init(X, P);

    /* z[j] = sum_i rf[i] v[i + j] with rf the reversal of f */
    _nmod_poly_reverse(rf, f, an, an);
    fft_small_fft_nmod(B, rf, an, atrunc, mod, P);
    fft_small_ifft_tr_nmod(X, v, vn, mod, P);
    fft_small_op_mul(X, X, B, P);
    fft_small_fft_tr(X, ztrunc, P);
    fft_small_export_nmod(z, X, mod, P);

    /* the same as the coefficients [an - 1, an - 1 + zn) of f v */
    _nmod_poly_mulmid(ref, v, vn, f, an, an - 1, an - 1 + zn, mod);
    flint_printf("middle product: %s\n", _nmod_vec_equal(z, ref, zn) ? "agrees" : "WRONG");

    fft_small_op_clear(B);
    fft_small_op_clear(X);
    fft_small_plan_clear(P);
    _nmod_vec_clear(f);
    _nmod_vec_clear(rf);
    _nmod_vec_clear(v);
    _nmod_vec_clear(z);
    _nmod_vec_clear(ref);
    flint_cleanup_master();
    return 0;
}
```
```output
3 primes, transform length 2^11, truncation 2048
middle product: agrees
```

A window of the full product of $f$ and $v$ would need a convolution of
length $a_n + v_n - 1 = 2998$, here reduced to 2048 by wrap-around
([](algorithms/polynomials/middle-product#lem-wrap)) because a power of
two lies between 1999 and 2998; for lengths where none does, the
transposed product is the shorter one.
