---
title: Inverting all entries
functions: _nmod_vec_invert
---

A modular inversion costs an extended Euclidean algorithm, much more than
a multiplication. Using Montgomery's simultaneous inversion trick
[@Mon1987; @BZ2011, §2.5.1], `_nmod_vec_invert` and its `_generic` and
`_shoup` variants invert all $\mathit{len}$ entries with a single
inversion and $3(\mathit{len} - 1)$ multiplications. A first pass forms
the prefix products $v_0 v_1 \cdots v_k$ and the last of them is inverted.
A second pass, from the last entry down, multiplies this inverse by the
prefix products and the entries to obtain each $v_k^{-1}$.

## Non-invertible entries {#impinv}

If one entry is not invertible modulo $n$, the product of all entries is
not invertible either, and the `FLINT_IMPINV` exception is raised by
`n_invmod` when it inverts that product. The exception does not tell which
entry failed. To find it, or to handle non-invertible entries, invert the
entries one by one with `n_gcdinv`.

## Aliasing {#aliasing}

The prefix products are kept in temporary storage, and the second pass
reads `vec[k]` before it writes `res[k]`, so `res` may be the same array
as `vec` (but not a partially overlapping one).

## Variants {#variants}

- `_nmod_vec_invert_generic` multiplies with `nmod_mul`. It allocates
  $\mathit{len}$ limbs of temporary storage.
- `_nmod_vec_invert_shoup` stores the precomputed Shoup quotients of the
  entries and of the prefix products. Those of the prefix products come
  from `n_mulmod_and_precomp_shoup` without a division. Every
  multiplication in the second pass is then a
  [Shoup multiplication](algorithms/arithmetic/shoup). It allocates
  $3 \, \mathit{len}$ limbs.
- `_nmod_vec_invert_naive` calls `nmod_inv`, and thus `n_invmod`, on each
  entry. It is used for testing.

`_nmod_vec_invert` calls `_nmod_vec_invert_shoup` when
$n < 2^{\bits - 1}$ and `_nmod_vec_invert_generic` otherwise.

## Example {#example}

```c example name=nmod_vec_invert
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

int main(void)
{
    nmod_t mod;
    ulong a[8], inv[8];
    slong i;

    nmod_init(&mod, 101);
    for (i = 0; i < 8; i++)
        a[i] = i + 2;

    /* one modular inversion for all eight entries */
    _nmod_vec_invert(inv, a, 8, mod);

    for (i = 0; i < 8; i++)
        flint_printf("1/%wu = %wu  (check %wu)\n", a[i], inv[i],
                     nmod_mul(a[i], inv[i], mod));
    return 0;
}
```
```output
1/2 = 51  (check 1)
1/3 = 34  (check 1)
1/4 = 76  (check 1)
1/5 = 81  (check 1)
1/6 = 17  (check 1)
1/7 = 29  (check 1)
1/8 = 38  (check 1)
1/9 = 45  (check 1)
```
