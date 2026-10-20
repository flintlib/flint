---
title: Converting between radix B and binary
functions: radix_get_mpn radix_set_mpn
---

`radix_set_mpn` converts a binary limb array to radix $B = b^e$ and
`radix_get_mpn` converts back. For $B = 10^{19}$ the first is the
expensive half of printing a number in decimal, and FLINT's
`flint_mpn_get_str` (hence `fmpz_get_str`) uses it for base 10 when
`fft_small` is available; once the number is in radix $10^{19}$, writing the
digits takes linear time (`radix_get_str_decimal`).

Both directions use divide and conquer with multiplications only. To
convert from radix $B$ to binary, the input is split as
$a = a_{\mathrm{lo}} + a_{\mathrm{hi}} B^k$ with $k$ about half the length,
both halves are converted recursively, and the result is
$a_{\mathrm{lo}} + a_{\mathrm{hi}} \cdot B^k$, evaluated with binary
arithmetic and a precomputed binary power $B^k$. To convert from binary
to radix $B$, the binary input is read as a number in radix $c = 2^s$ with
$c < B$ (cutting it into $s$-bit chunks, which is free), and the same
recursion is run in reverse: the halves are converted recursively and
combined with radix-$B$ arithmetic (`radix_mul`, which uses `fft_small`)
and a precomputed radix-$B$ power $c^k$. No division by a large power is
needed, unlike in the classical subquadratic conversion.

::: summary algorithms/integers/radix-conversion :::

## Sizes and thresholds {#thresholds}

- `radix_get_mpn` uses Horner's rule (`radix_get_mpn_basecase`) below 50
  limbs, also at the leaves of the recursion.
- `radix_set_mpn` uses repeated division by $B$
  (`radix_set_mpn_basecase`) below 40 input limbs, and the recursion
  switches to it on pieces of fewer than 25 chunks.
- The tables of powers ($c^{k_0}, c^{k_1}, \ldots$ with
  $k_{i+1} = \lceil k_i / 2 \rceil$) are built by repeated squaring, with
  one exact division by $c$ whenever $k_i$ is odd. Trailing zero limbs of
  the powers (the powers of $B = 10^{19} = 2^{19} 5^{19}$ have many zero
  bits) are skipped in the multiplications. `radix_set_mpn` caches the
  most recent table per thread.
- Above 500 limbs, the two recursive halves are converted in parallel
  when threads are available (see [](guide/conventions/threading)).

## Example {#example}

Convert $3^{200}$ from binary to radix $10^{19}$, print it, and convert it
back:

```c example name=radix_conversion
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/radix.h>

int main(void)
{
    radix_t radix;
    fmpz_t x;
    nn_ptr d, t;
    slong xn, dn, tn;
    char * s;

    radix_init(radix, 10, 0);   /* B = 10^19 on 64-bit machines */

    fmpz_init(x);
    fmpz_ui_pow_ui(x, 3, 200);
    xn = fmpz_size(x);

    d = flint_malloc(radix_set_mpn_need_alloc(xn, radix) * sizeof(ulong));
    dn = radix_set_mpn(d, COEFF_TO_PTR(*x)->_mp_d, xn, radix);

    s = radix_get_str_decimal(NULL, d, dn, 0, radix);
    flint_printf("%s\n", s);
    flint_printf("%wd digits, %wd decimal limbs\n",
        (slong) _radix_size_digits(d, dn, radix), dn);

    t = flint_malloc((xn + 1) * sizeof(ulong));
    tn = radix_get_mpn(t, d, dn, radix);
    flint_printf("round trip: %d\n",
        tn == xn && mpn_cmp(t, COEFF_TO_PTR(*x)->_mp_d, xn) == 0);

    flint_free(s);
    flint_free(d);
    flint_free(t);
    fmpz_clear(x);
    radix_clear(radix);
    return 0;
}
```
```output
265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699044001
96 digits, 6 decimal limbs
round trip: 1
```
