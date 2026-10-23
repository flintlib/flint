---
title: Portable vec4d code
summary: Writing `vec4d` code that compiles on every backend and gives the same results everywhere, and the few operations whose results differ.
order: 1
---

`vec4d` is the one vector type with a large set of operations on every
backend: AVX2, AVX-512, NEON, GNU vector extensions and strict ISO C.
Code written with these operations compiles everywhere and, with the
exceptions listed under [Pitfalls](#pitfalls), gives the same results on
every backend. FLINT's `dfloat` module is written this way.

## What is portable {#portable}

The following `vec4d` operations exist on every backend:

- loads and stores (`vec4d_load`, `vec4d_store`, and their `_aligned`
  and `_unaligned` forms), `vec4d_set_d`, `vec4d_set_d4`, `vec4d_zero`,
  `vec4d_one` and `vec4d_get_index`;
- `vec4d_add`, `vec4d_sub`, `vec4d_mul`, `vec4d_div`, `vec4d_neg`,
  `vec4d_abs`, `vec4d_min`, `vec4d_max`, `vec4d_sqrt`, `vec4d_round`
  and `vec4d_floor`;
- `vec4d_fmadd`, `vec4d_fnmadd` and `vec4d_fmsub`, and the
  `_fused` variants;
- the comparisons `vec4d_cmp_*`, `vec4d_blendv`, `vec4d_movemask`, the
  bitwise operations `vec4d_bit_*` and the integer operations
  `vec4d_bits_*`;
- `vec4d_unpacklo`, `vec4d_unpackhi`, the six `vec4d_permute_*`, the two
  `vec4d_permute2_*`, `vec4d_to_index` and `vec4d_gather`.

The modular functions such as `vec4d_mulmod` are not in this list; they
need AVX2 or NEON. The second example below shows how to write them with
portable operations.

## Pitfalls {#pitfalls}

- **Fused or not.** `vec4d_fmadd`, `vec4d_fnmadd` and `vec4d_fmsub`
  round once on AVX2 and NEON but may round twice on the generic
  backends, depending on whether the compiler contracts `a * b + c`.
  Use `vec4d_fmadd_fused`, `vec4d_fmsub_fused` and `vec4d_fnmadd_fused`
  whenever the result must be exact or reproducible, as in error-free
  transformations.
- **Rounding mode.** `vec4d_round` follows the current rounding mode,
  except on NEON, which always rounds to nearest with ties to even.
- **Alignment.** `vec4d_load` and `vec4d_store` require 32-byte aligned
  addresses on AVX2 and may fault otherwise. Use `_Alignas(32)` or the
  `_unaligned` forms.
- **NaNs in min and max.** `vec4d_min` and `vec4d_max` return the lane of
  `b` for unordered lanes, except on NEON, which returns a NaN.
- **Masks are doubles.** A comparison returns a `vec4d` whose lanes are
  all-ones or all-zeros bit patterns. Such a lane is a NaN or $+0$ as a
  double, so masks should only be passed to `vec4d_blendv`,
  `vec4d_movemask` and the `vec4d_bit_*` operations.
- **ABI warnings.** On x86-64 without AVX, GCC and clang warn
  (`-Wpsabi`) that 32-byte vector arguments are passed differently than
  in older compiler versions. This is harmless for inline functions, but
  becomes an error with `-Werror`. The examples silence it with a
  pragma.

## Examples {#examples}

Evaluating a polynomial at four points with Horner's rule, then using a
comparison mask to find and clear the negative values. The lanes are
read back with `vec4d_get_index` and with a store.

```c example name=mv_vec4d_horner
/* without AVX, GCC and clang note that 32-byte vectors are passed
   differently; this is harmless for inline functions */
#if defined(__clang__)
# if __has_warning("-Wpsabi")
#  pragma clang diagnostic ignored "-Wpsabi"
# endif
#elif defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wpsabi"
#endif

#include <flint/flint.h>
#include <flint/machine_vectors.h>

int main(void)
{
    /* p(x) = 2x^3 - 3x^2 - 1, constant term first */
    double c[4] = { -1, 0, -3, 2 };
    /* vec4d_load needs 32-byte alignment on AVX2 */
    _Alignas(32) double xs[4] = { -1, 0, 1, 2 };
    double out[4];
    vec4d x, y, zero = vec4d_zero();
    int i;

    x = vec4d_load(xs);

    /* Horner's rule in four lanes at once */
    y = vec4d_set_d(c[3]);
    for (i = 2; i >= 0; i--)
        y = vec4d_fmadd(y, x, vec4d_set_d(c[i]));

    vec4d_store_unaligned(out, y);
    flint_printf("p(x)     = %wd %wd %wd %wd\n", (slong) out[0],
                 (slong) out[1], (slong) out[2], (slong) out[3]);

    /* lanes where p(x) < 0, as a bit mask */
    flint_printf("negative = %d\n", vec4d_movemask(vec4d_cmp_lt(y, zero)));

    /* replace negative lanes by zero */
    y = vec4d_blendv(y, zero, vec4d_cmp_lt(y, zero));
    flint_printf("max(p,0) = %wd %wd %wd %wd\n",
                 (slong) vec4d_get_index(y, 0), (slong) vec4d_get_index(y, 1),
                 (slong) vec4d_get_index(y, 2), (slong) vec4d_get_index(y, 3));
    return 0;
}
```
```output
p(x)     = -6 -1 -2 3
negative = 7
max(p,0) = 0 0 0 3
```

Modular multiplication of four pairs at once, with the formula of
`vec4d_mulmod` written in portable operations. The modulus is a prime
just above $2^{48}$ and the inputs are reduced, so by
[](algorithms/simd/fp-modular-arithmetic) the result lies in $(-n, n)$.
The fused operations make the splitting of the product $ab = h + \ell$
exact on every backend. The program compares 40000 products with
`n_mulmod2_preinv`.

To map the result to $[0, n)$, `vec4d_blendv` with the result itself as
the mask adds $n$ to the lanes whose sign bit is set. A zero result
comes out as $+0$, so the blend never adds $n$ to it: in round to
nearest an exact zero sum is $-0$ only if both operands are $-0$, and
$\ell$ is never $-0$ for inputs in $[0, n)$.

```c example name=mv_vec4d_mulmod
#if defined(__clang__)
# if __has_warning("-Wpsabi")
#  pragma clang diagnostic ignored "-Wpsabi"
# endif
#elif defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wpsabi"
#endif

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/machine_vectors.h>

/* the formula of vec4d_mulmod with operations every backend provides */
FLINT_FORCE_INLINE vec4d mulmod4(vec4d a, vec4d b, vec4d n, vec4d ninv)
{
    vec4d h = vec4d_mul(a, b);                    /* a*b rounded */
    vec4d q = vec4d_round(vec4d_mul(h, ninv));    /* quotient */
    vec4d l = vec4d_fmsub_fused(a, b, h);         /* a*b - h, exact */
    return vec4d_add(vec4d_fnmadd_fused(q, n, h), l);
}

int main(void)
{
    ulong n = n_nextprime(UWORD(1) << 48, 1);     /* n < 2^49 */
    ulong ninv_word = n_preinvert_limb(n);
    vec4d vn = vec4d_set_d((double) n);
    vec4d vninv = vec4d_set_d(1.0 / (double) n);
    ulong a[4], b[4], x = 1, y = 2;
    double r[4];
    slong i, j, bad = 0;

    for (i = 0; i < 10000; i++)
    {
        /* reduced inputs from two simple sequences */
        for (j = 0; j < 4; j++)
        {
            x = n_mulmod2_preinv(x, 3, n, ninv_word);
            y = n_mulmod2_preinv(y, 5, n, ninv_word);
            a[j] = x;
            b[j] = y;
        }

        vec4d va = vec4d_set_d4(a[0], a[1], a[2], a[3]);
        vec4d vb = vec4d_set_d4(b[0], b[1], b[2], b[3]);
        vec4d vr = mulmod4(va, vb, vn, vninv);    /* in (-n, n) */

        /* add n to the negative lanes to get [0, n) */
        vr = vec4d_blendv(vr, vec4d_add(vr, vn), vr);
        vec4d_store_unaligned(r, vr);

        for (j = 0; j < 4; j++)
            if ((ulong) r[j] != n_mulmod2_preinv(a[j], b[j], n, ninv_word))
                bad++;
    }

    flint_printf("n = %wu\n", n);
    flint_printf("last: %wu * %wu mod n = %wu\n", a[3], b[3], (ulong) r[3]);
    flint_printf("40000 products, %wd wrong\n", bad);
    return 0;
}
```
```output
n = 281474976710677
last: 234264994381937 * 89140169141278 mod n = 118661419208307
40000 products, 0 wrong
```
