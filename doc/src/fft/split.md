---
title: Splitting integers into coefficients
functions: fft_split_limbs fft_split_bits fft_combine_limbs fft_combine_bits
---

To multiply integers with a convolution, they are cut into chunks that
become polynomial coefficients, and the product polynomial is evaluated
at the chunk base to recover the product. These four functions do the
cutting and the evaluation; `mul_truncate_sqrt2`, `mul_mfa_truncate_sqrt2`
and `_fft_mulmod_2expp1` use them, and `fmpz_mat_mul_fft` uses them for
matrix entries.

## Splitting {#split}

`fft_split_bits(poly, x, total_limbs, bits, output_limbs)` writes chunk
$i$ of $x$, that is $\lfloor x / 2^{i \cdot \mathit{bits}} \rfloor \bmod 2^{\mathit{bits}}$,
to `poly[i]`, zero-padded to $\mathit{output\_limbs} + 1$ limbs, for
$0 \le i < \lceil \bits \cdot \mathit{total\_limbs} / \mathit{bits} \rceil$,
and returns this count. The count depends only on `total_limbs`, so
leading chunks may be zero. If `bits` is a multiple of $\bits$ it calls
`fft_split_limbs`, which copies whole limbs. The chunks need
$\lceil \mathit{bits}/\bits \rceil \le \mathit{output\_limbs}$ limbs.
Entries of `poly` beyond the count are not touched, so callers zero them
themselves when they pad the vector to the transform length. Both
functions distribute the work over the threads of the global thread pool
in blocks of 16 chunks.

## Combining {#combine}

`fft_combine_bits(res, poly, length, bits, output_limbs, total_limbs)`
adds $\sum_{i < \mathit{length}} \mathit{poly}_i \cdot 2^{i \cdot \mathit{bits}}$
to `(res, total_limbs)` and discards everything beyond `total_limbs`
limbs; `fft_combine_limbs` does the same with $B^{\mathit{coeff\_limbs}}$.
The coefficients must be less than $B^{\mathit{output\_limbs}}$, which
holds for reduced residues modulo $2^N + 1$ with
$\mathit{output\_limbs} = N/\bits$ unless a residue is $2^N$, and they must
have $\mathit{output\_limbs} + 1$ readable limbs, since `fft_combine_bits`
shifts the top limb in when `bits` is not a multiple of $\bits$.

Each coefficient is added into a window of $\mathit{output\_limbs} + 1$
limbs starting at the limb where it begins. When `res` starts at zero,
the sum of the coefficients added so far is less than
$2^{\mathit{output\_limbs} \cdot \bits + i \cdot \mathit{bits} + 1}$ by
a geometric series, which fits in the window of coefficient $i$, so no
carry leaves the window and the result is exact.

For this reason `res` must be zero on entry, which builds configured
with `--enable-assert` check. When `res` is not zero, carries out of a
window are lost, so the functions do not add to `res` in general: with
`res` $= 2^{192} - 1$ in four limbs and the constant polynomial $1$,
`fft_combine_limbs(res, poly, 1, 1, 1, 4)` would yield
$2^{192} - 2^{128}$ instead of $2^{192}$.

## Example {#example}

A 192-bit integer cut into 20-bit chunks and put together again. The top
limb holds only 5 bits, but all ten chunks of the three limbs are
written.

```c example name=fft_split_bits
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/mpn_extras.h>

int main(void)
{
    mp_limb_t x[3] = { UWORD(0x0123456789abcdef), UWORD(0xfedcba9876543210), UWORD(0x1f) };
    mp_limb_t y[3] = { 0, 0, 0 };
    mp_limb_t mem[10 * 2], * poly[10];
    slong i, len;

    for (i = 0; i < 10; i++)
        poly[i] = mem + 2 * i;           /* output_limbs + 1 = 2 limbs each */

    /* 3 limbs = 192 bits -> 10 coefficients of 20 bits */
    len = fft_split_bits(poly, x, 3, 20, 1);
    flint_printf("%wd coefficients:", len);
    for (i = 0; i < len; i++)
        flint_printf(" %wx", poly[i][0]);
    flint_printf("\n");

    /* evaluate at 2^20 into a zeroed result */
    fft_combine_bits(y, poly, len, 20, 1, 3);
    flint_printf("round trip: %d\n", mpn_cmp(x, y, 3) == 0);
    return 0;
}
```
```output
10 coefficients: bcdef 6789a 12345 32100 87654 dcba9 1ffe 0 0 0
round trip: 1
```
