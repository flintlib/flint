---
title: Computing dot products
functions: dot_params_t _nmod_vec_dot_params NMOD_VEC_DOT _nmod_vec_dot_bound_limbs
---

A dot product $\sum_{i=0}^{\mathit{len}-1} a_i b_i \bmod n$ is the inner
loop of classical polynomial multiplication, matrix multiplication,
triangular solving and power series inversion over $\Z/n\Z$. FLINT computes
it as an exact integer sum, which needs one, two or three limbs depending
on $n$ and the length, and reduces modulo $n$ only once at the end. The
best way to form the sum depends on the number of limbs, the size of $n$
and the instruction set, so it is chosen once, in advance.

::: summary algorithms/linear-algebra/dot-products :::

## Usage {#usage}

1. Call `_nmod_vec_dot_params(len, mod)` with the largest length that will
   be used. This returns a `dot_params_t`.
2. Pass it to `_nmod_vec_dot`, `_nmod_vec_dot_rev`, `_nmod_vec_dot_ptr`,
   `_nmod_vec_dot_strided` or `NMOD_VEC_DOT` for any length up to `len`
   and the same modulus.

The parameters are cheap to compute, but not free: `_nmod_vec_dot_params`
performs a few multiplications and, for some methods, a reduction. Callers
with many short dot products of varying length, such as
`_nmod_poly_mul_classical`, compute them once for the longest length.
This is always correct, because a method that is valid for a length is
valid for every shorter one.

`_nmod_vec_dot_params` is a library function, not an inline one, so the
method it selects always matches the kernels FLINT was compiled with
(AVX2, AVX-512, IFMA), whatever the compiler flags of the calling program.

All entries must be reduced modulo $n$. The method is chosen from the bound
$\mathit{len} \, (n - 1)^2$ on the unreduced sum, so unreduced entries can
make the sum overflow without any warning.

## Variants {#variants}

- `_nmod_vec_dot` multiplies `vec1[i]` by `vec2[i]`. It is the fastest
  variant and the only one besides `_nmod_vec_dot_rev` that uses the
  vectorised kernels for moduli above $2^{32}$.
- `_nmod_vec_dot_rev` multiplies `vec1[i]` by `vec2[len - 1 - i]`, as in
  the coefficients of a polynomial product.
- `_nmod_vec_dot_ptr` multiplies `vec1[i]` by `vec2[i][offset]`, a column
  of a matrix stored as an array of row pointers. It never uses the
  vectorised kernels for moduli above $2^{32}$.
- `_nmod_vec_dot_strided` reads `vec1[i * stride1]` and
  `vec2[i * stride2]`. The strides may be negative or zero. It is an inline
  function built on `NMOD_VEC_DOT`.
- `NMOD_VEC_DOT` takes two arbitrary expressions in the loop variable. It
  is expanded inline and is not explicitly vectorised.

The first three variants handle lengths up to 11 with fully unrolled inline
code when the method is a plain one-, two- or three-limb accumulation, which
is always the case for parameters computed for a length of at most 11.

## The number of limbs {#limbs}

The bound $\mathit{len} \, (n - 1)^2$ measured by
`_nmod_vec_dot_bound_limbs` is attained when all entries are $n - 1$, so
its number of limbs is the least that suffices for every input. Callers
that implement their own accumulation, for example in matrix
multiplication, use this function to choose between one-, two- and
three-limb code.
The method returned by `_nmod_vec_dot_params` encodes the same information,
which `_nmod_vec_dot_bound_limbs_from_params` reads off without any
multiplication.

On a 64-bit machine, one limb suffices for any length below $2^{64}/(n-1)^2$,
in particular for $n \le 2^{32}$ and $\mathit{len} = 1$. Two limbs suffice
for every length when $n \le 2^{32}$, and three limbs always suffice. The
exact thresholds are derived in
[How many words](algorithms/linear-algebra/dot-products#limbs).

## Examples {#examples}

Choosing the method once and checking the result against a reduction
after every term:

```c example name=nmod_vec_dot_basic
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

int main(void)
{
    nmod_t mod;
    ulong a[1000], b[1000], s, t;
    slong i, len = 1000;
    dot_params_t params;

    nmod_init(&mod, 1000003);
    for (i = 0; i < len; i++)
    {
        a[i] = (i * i) % mod.n;
        b[i] = (7 * i + 1) % mod.n;
    }

    /* choose a method once, for the longest length */
    params = _nmod_vec_dot_params(len, mod);

    s = _nmod_vec_dot(a, b, len, mod, params);

    /* the same sum, reducing after every term */
    t = 0;
    for (i = 0; i < len; i++)
        t = nmod_add(t, nmod_mul(a[i], b[i], mod), mod);

    flint_printf("dot      = %wu\n", s);
    flint_printf("check    = %wu\n", t);

    /* the same parameters serve any shorter length */
    flint_printf("dot[10]  = %wu\n", _nmod_vec_dot(a, b, 10, mod, params));
    flint_printf("rev[10]  = %wu\n", _nmod_vec_dot_rev(a, b, 10, mod, params));
    return 0;
}
```
```output
dot      = 343013
check    = 343013
dot[10]  = 14460
rev[10]  = 4065
```

A matrix-vector product with `_nmod_vec_dot`, and a column of the matrix
times a vector, once with a stride and once with `NMOD_VEC_DOT` and an
index expression:

```c example name=nmod_vec_dot_matrix
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

int main(void)
{
    nmod_t mod;
    ulong A[4 * 5], x[5], y[4], r;
    slong i, j, k;
    dot_params_t params;

    nmod_init(&mod, 65537);

    /* a 4 x 5 matrix stored by rows, and a vector x */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 5; j++)
            A[i * 5 + j] = nmod_pow_ui(i + 2, j + 10, mod);
    for (j = 0; j < 5; j++)
        x[j] = j + 1;

    /* y = A x: each entry is a dot product with a row */
    params = _nmod_vec_dot_params(5, mod);
    for (i = 0; i < 4; i++)
        y[i] = _nmod_vec_dot(A + i * 5, x, 5, mod, params);
    flint_printf("A x     = %wu %wu %wu %wu\n", y[0], y[1], y[2], y[3]);

    /* column 2 of A times y: the column has stride 5 */
    params = _nmod_vec_dot_params(4, mod);
    r = _nmod_vec_dot_strided(A + 2, 5, y, 1, 4, mod, params);
    flint_printf("col . y = %wu\n", r);

    /* the same with NMOD_VEC_DOT */
    NMOD_VEC_DOT(r, k, 4, A[k * 5 + 2], y[k], mod, params);
    flint_printf("col . y = %wu\n", r);
    return 0;
}
```
```output
A x     = 1022 55599 40049 42874
col . y = 10437
col . y = 10437
```

The number of limbs of the unreduced sum for a few moduli and lengths
(64-bit only, since some moduli exceed $2^{32}$):

```c example name=nmod_vec_dot_bound_limbs
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

int main(void)
{
    /* primes of 20, 31, 33, 50, 62 and 64 bits */
    ulong ns[6] = { 1048573, 2147483647, UWORD(8589934583),
                    UWORD(1125899906842597), UWORD(4611686018427387847),
                    UWORD(18446744073709551557) };
    slong lens[3] = { 10, 1000, 1000000 };
    int i, j;

    for (i = 0; i < 6; i++)
    {
        nmod_t mod;
        nmod_init(&mod, ns[i]);
        flint_printf("%2wu bits:", FLINT_BIT_COUNT(ns[i]));
        for (j = 0; j < 3; j++)
            flint_printf("  len %7wd -> %d", lens[j],
                         _nmod_vec_dot_bound_limbs(lens[j], mod));
        flint_printf("\n");
    }
    return 0;
}
```
```output
20 bits:  len      10 -> 1  len    1000 -> 1  len 1000000 -> 1
31 bits:  len      10 -> 2  len    1000 -> 2  len 1000000 -> 2
33 bits:  len      10 -> 2  len    1000 -> 2  len 1000000 -> 2
50 bits:  len      10 -> 2  len    1000 -> 2  len 1000000 -> 2
62 bits:  len      10 -> 2  len    1000 -> 3  len 1000000 -> 3
64 bits:  len      10 -> 3  len    1000 -> 3  len 1000000 -> 3
```
