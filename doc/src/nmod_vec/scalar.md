---
title: Multiplication by a fixed scalar
functions: _nmod_vec_scalar_mul_nmod _nmod_vec_scalar_addmul_nmod
---

When every entry of a vector is multiplied by the same scalar $c$, data
for multiplying by $c$ can be precomputed once and used for every entry.
`_nmod_vec_scalar_mul_nmod` and `_nmod_vec_scalar_addmul_nmod` each have
variants that differ in what they precompute and in their requirements on
$n$, and choose among them from $c$, $n$ and the length.

## The variants {#variants}

| Variant | Method | Requires |
|---|---|---|
| `_generic` | product reduced with the inverse in `mod`, as in `nmod_mul` | nothing |
| `_shoup` | [Shoup multiplication](algorithms/arithmetic/shoup) with $\lfloor c 2^{\bits} / n \rfloor$ | $n < 2^{\bits - 1}$ |
| `_redc` | Montgomery multiplication [@Mon1985] with $c$ in Montgomery form | $n$ odd |

Shoup multiplication needs one high and two low products per entry and no
reduction step, against a full product and a two-word reduction for the
generic code. It is usually the fastest choice but needs one spare bit in
$n$. Montgomery multiplication handles moduli with the top bit set, where
Shoup's method does not apply. Only `_nmod_vec_scalar_mul_nmod` has a
`_redc` variant.

::: summary algorithms/arithmetic/shoup :::

## Choice of variant {#choice}

`_nmod_vec_scalar_mul_nmod` handles $c = 0$, $c = 1$ and $c = n - 1$ by
zeroing, copying and negating. Otherwise:

- If $n$ has `FLINT_BITS` bits, it uses `_redc` for odd $n$ and length at
  least 8, and otherwise a generic loop without normalisation shifts.
- If $n$ has fewer bits, it uses `_shoup` from length
  `FLINT_MULMOD_SHOUP_THRESHOLD` on, and `_generic` below.

`_nmod_vec_scalar_addmul_nmod` handles $c = 0$, $c = 1$ and $c = n - 1$
with no operation, `_nmod_vec_add` and `_nmod_vec_sub`. Otherwise it uses
`_shoup` from length `FLINT_MULMOD_SHOUP_THRESHOLD` on when $n$ has fewer
than `FLINT_BITS` bits, and `_generic` in all other cases.

The threshold `FLINT_MULMOD_SHOUP_THRESHOLD` is set in `flint-mparam.h` for
each target; on some recent x86-64 processors it is 0 (always use Shoup).
The precomputation for Shoup's method
costs one hardware division, which a short vector does not amortise.

## Example {#example}

```c example name=nmod_vec_scalar
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

int main(void)
{
    nmod_t mod;
    ulong u[5] = { 1, 2, 3, 4, 5 };
    ulong v[5] = { 10, 20, 30, 40, 50 };
    ulong w[5], x[5], y[5];
    ulong c = 12345;
    int i;

    nmod_init(&mod, 1000003);

    /* w = c u + 3 v - v */
    _nmod_vec_scalar_mul_nmod(w, u, 5, c, mod);
    _nmod_vec_scalar_addmul_nmod(w, v, 5, 3, mod);
    _nmod_vec_sub(w, w, v, 5, mod);
    _nmod_vec_print_pretty(w, 5, mod);

    /* the variants agree; _shoup also accepts unreduced entries */
    _nmod_vec_scalar_mul_nmod_generic(x, v, 5, c, mod);
    _nmod_vec_scalar_mul_nmod_redc(y, v, 5, c, mod);     /* n is odd */
    flint_printf("generic = redc: %d\n", _nmod_vec_equal(x, y, 5));

    for (i = 0; i < 5; i++)
        y[i] = v[i] + 7 * mod.n;                         /* not reduced */
    _nmod_vec_scalar_mul_nmod_shoup(y, y, 5, c, mod);    /* n < 2^(FLINT_BITS - 1) */
    flint_printf("generic = shoup: %d\n", _nmod_vec_equal(x, y, 5));
    return 0;
}
```
```output
<length-5 integer vector mod 1000003>
[  12365   24730   37095   49460   61825]
generic = redc: 1
generic = shoup: 1
```
