---
title: Multiplying polynomials
functions: nmod_poly_mul nmod_poly_mullow nmod_poly_mulhigh nmod_poly_mulmid
---

`nmod_poly_mul` is the function to use for products in
$(\Z/n\Z)[x]$. It chooses among the classical algorithm, three variants of
Kronecker substitution and the small-prime FFT of `fft_small` according to
the lengths and the bit size of $n$. `nmod_poly_mullow`,
`nmod_poly_mulhigh` and `nmod_poly_mulmid` compute the low part, the high
part or an arbitrary range of coefficients of the product; the truncated
and middle products are the building blocks of power series arithmetic and
[Newton iteration](algorithms/polynomials/newton).

::: summary algorithms/polynomials/multiplication :::

## Algorithm selection {#selection}

`_nmod_poly_mul` uses the classical algorithm, or the rows described
below, if the shorter operand has length at most 5. Otherwise it uses the
IFMA kernel described below if it applies, then asks
`_nmod_poly_mullow_want_fft_small` whether `fft_small` is faster, then
uses the rows if they apply, and otherwise, with $c =
\min(\mathit{len1}, 2\,\mathit{len2})$ and $b$ the number of bits of $n$:

- classical if $3c < 2 \max(b, 10)$;
- [Kronecker substitution](algorithms/polynomials/kronecker) (`KS`) if
  $cb < 800$;
- `KS2`, evaluation at $\pm 2^k$, if $c(b+1)^2 < 100\,000$;
- `KS4`, evaluation at $\pm 2^k$ and $\pm 2^{-k}$, otherwise [@Har2009a].

For balanced operands of length $L$ and a few bit sizes this gives the
following ranges, measured with the installed library on x86-64 with
`fft_small`:

| Bits of $n$ | Classical | KS | KS2 | KS4 | `fft_small` |
|---|---|---|---|---|---|
| 8 | $L \le 6$ | 7–99 | 100–305 | — | $\ge 306$ |
| 16 | $\le 10$ | 11–49 | 50–98 | — | $\ge 99$ |
| 20 | $\le 13$ | 14–39 | 40–81 | — | $\ge 82$ |
| 30 | $\le 19$ | 20–26 | 27–104 | 105–184 | $\ge 185$ |
| 40 | $\le 26$ | — | 27–59 | 60–91 | $\ge 92$ |
| 50 | $\le 33$ | — | 34–38 | 39–161 | $\ge 162$ |
| 64 | $\le 42$ | — | — | 43–93 | $\ge 94$ |

With AVX512-IFMA, the IFMA kernel takes the lengths from 8 up to its bound
for $b \le 52$.

The `fft_small` cutoffs come from tables indexed by the bit size and are
somewhat higher for squaring. A modulus that is one of the 50-bit primes
of `fft_small` needs no Chinese remaindering, and `fft_small` is then used
from $c \ge 24$. For $n \le 23$, two coefficients are packed into one
before the FFT, which halves the transform length; this is used for
product lengths above 370, except just above the transform-length jumps
at 512 and 1024.

The classical algorithm computes each output coefficient as one dot
product with `_nmod_vec_dot_rev` and a single reduction (see
[](nmod_vec/dot)). Kronecker substitution packs the coefficients into
integers with fields wide enough for the unreduced coefficients of the
product, about $2b + \log_2 L$ bits, and multiplies them with
`flint_mpn_mul` ([](mpn_extras/mul)), which
is efficient when $b$ is small compared to the word size. `KS2` and
`KS4` replace one large integer product by two or four products of half
or a quarter of the field width.

For a shorter operand of length $\mathit{len2}$ from 3 to 16 (to 6 with
NEON, where longer rows lose), the products of the classical algorithm are
computed in rows: one `_nmod_vec_scalar_addmul_nmod` of the longer operand
per coefficient of the shorter one, rather than one dot product of at most
$\mathit{len2}$ terms per output coefficient, whose call, setup and final
reduction dominate for such short lengths. The rows run the vector kernels
of these scalar products ([](nmod_vec/simd)), so they are used where a
kernel exists for $b$ (`NMOD_VEC_SIMD_SCALAR_MAX_BITS`: up to 63 bits with
AVX-512, 50 with AVX2, 32 with NEON), $n$ is not a power of 2, the
classical algorithm does not accumulate its products in one word (from 30
to 32 bits on, depending on $\mathit{len2}$), the product is not a square,
and every row has at least `NMOD_VEC_SIMD_MIN_LEN` coefficients and $r
\cdot \mathit{len2}$, with $r = 2$ for the 512-bit kernels and up to 32
bits and $r = 4$ for the AVX2 kernels above 32 bits; middle products use
$2r$, since all their dot products have $\mathit{len2}$ terms. They
replace the classical algorithm and Kronecker substitution, not
`fft_small`.

When FLINT is built for a processor with AVX512-IFMA (Intel from Ice Lake
on for some machines, AMD from Zen 4 on), the products for $n \le 2^{52}$ and a shorter
operand of 8 to $L$ coefficients (4 to $L$ for 52-bit moduli) go to a
vector kernel of the classical algorithm, which computes any range of the
product, in place of the rows, the classical algorithm, Kronecker
substitution and `fft_small`. It computes the output coefficients by
blocks of 48 or 64, whose sums stay in vector registers for all their
terms: each step multiplies one coefficient of the shorter operand,
broadcast, by 8 consecutive coefficients of the longer one with
`vpmadd52luq` and `vpmadd52huq`, which accumulate the low and the high
52-bit halves of the products (the low halves alone for $n \le 2^{26}$).
With at most 511 terms the sums stay below $2^{61}$, and each output
coefficient takes one reduction at the end: a Montgomery reduction in
radix $2^{52}$ for odd $n > 2^{26}$, the shorter operand having been
multiplied by $2^{104}$ modulo $n$ beforehand, and otherwise the reduction
of the IFMA kernel of `nmod_mat_mul`. Below 8 coefficients the rows are
faster (except for 52-bit moduli, where they run the slower scalar
products of 51 to 63 bits), and above $L$, `fft_small`. The bound $L$
depends on the processor (`FLINT_NMOD_POLY_IFMA_MAX_LEN` in
`flint-mparam.h`) and follows the number of 50-bit primes that `fft_small`
needs: on Ice Lake, $L = 192$ up to 20 bits (one prime), 511 for 21 to 26
bits (two primes, and one IFMA per product), 384 for 27 to 45 bits (two
primes) and 511 above (three primes), and 96 for the primes of `fft_small`
itself; on Zen 4, where `fft_small` is faster, 192, 384, 192, 384 and 64.

For $2^{52} < n < 2^{64}$ the kernel writes each coefficient as $x_0 +
2^{52} x_1$ with $x_1 < 2^{12}$ and accumulates the products in three
limbs of weights $1$, $2^{52}$ and $2^{104}$, with seven IFMA instructions
per product of vectors. Where `vpmullq` is fast (Zen 4), from
`FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN` coefficients, it computes instead the
exact sums modulo $2^{64}$ with `vpmullq` and `vpaddq`, modulo $p_1 =
2^{52} - 47$ with two IFMA and modulo $p_2 = 2^{26} - 5$ with one, and
recovers them by the Chinese remainder theorem (the sums are below $511
n^2 < 2^{137}$, and $2^{64} p_1 p_2 > 2^{141}$). The three limbs of each
output coefficient are reduced by a two-limb Montgomery reduction for odd
$n$, the shorter operand having been multiplied by $2^{104}$ modulo $n$,
and by a mask for powers of 2; other even moduli, which would need a
scalar reduction of each output, do not go to the kernel. The kernel is
used there for a shorter operand of 8 coefficients (12 on Zen 4, 4 for
64-bit moduli, which have no rows) to $L$, $L = 256$ on Intel and 128 on
Zen 4 (511 for low products), when the product has at least 16 output
coefficients.

`nmod_poly_mullow` uses the IFMA kernel under the same conditions, up to
511 coefficients for a truncated product (256 on Zen 4 for the primes of
`fft_small`), and the rows under the same conditions, its shortest row
having $\mathit{trunc} - \mathit{len2} + 1$ coefficients, otherwise the
classical algorithm for truncation lengths below $10 + b^2/10$ (from 16
for 8-bit to 419 for 64-bit moduli), KS otherwise, and `fft_small` from
separately tuned cutoffs. `nmod_poly_mulmid` first drops the operand
coefficients that do not contribute to the requested range and then uses
the IFMA kernel, the rows, the classical algorithm, `fft_small` (which
computes a middle product at the cost of a product of the output length,
[](algorithms/polynomials/middle-product)), or KS for small moduli.
`nmod_poly_mulhigh` calls `nmod_poly_mulmid` for all but tiny inputs and
sets the low coefficients to zero.

## Examples {#examples}

The parts of one product modulo 7:

```c example name=nmod_poly_mul_parts
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t a, b, c;

    nmod_poly_init(a, 7);
    nmod_poly_init(b, 7);
    nmod_poly_init(c, 7);

    nmod_poly_set_str(a, "4 7  1 2 3 4");     /* 4x^3 + 3x^2 + 2x + 1 */
    nmod_poly_set_str(b, "3 7  5 0 6");       /* 6x^2 + 5 */

    nmod_poly_mul(c, a, b);
    flint_printf("a*b          = "); nmod_poly_print(c); flint_printf("\n");
    nmod_poly_mullow(c, a, b, 3);
    flint_printf("mullow 3     = "); nmod_poly_print(c); flint_printf("\n");
    nmod_poly_mulmid(c, a, b, 2, 5);
    flint_printf("mulmid [2,5) = "); nmod_poly_print(c); flint_printf("\n");
    nmod_poly_mulhigh(c, a, b, 4);
    flint_printf("mulhigh 4    = "); nmod_poly_print(c); flint_printf("\n");

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    return 0;
}
```
```output
a*b          = 6 7  5 3 0 4 4 3
mullow 3     = 2 7  5 3
mulmid [2,5) = 3 7  0 4 4
mulhigh 4    = 6 7  0 0 0 0 4 3
```

The middle product consists of the coefficients 2, 3 and 4 of the
product, shifted down to start at 0. The truncated product ends in the
zero coefficient of $x^2$ and is normalised to length 2. The high product
has the full length; its coefficients below index 4 are zero here because
the classical algorithm was used, but in general they are only known to be
either zero or those of the product.

All algorithms give the same product. The program compares them for
random polynomials of length 2000 modulo a 60-bit prime, where
`nmod_poly_mul` uses `fft_small`:

```c example name=nmod_poly_mul_algorithms
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>

int main(void)
{
    flint_rand_t state;
    nmod_poly_t a, b, c, d;
    ulong p = n_nextprime(UWORD(1) << 59, 1);

    flint_rand_init(state);
    nmod_poly_init(a, p);
    nmod_poly_init(b, p);
    nmod_poly_init(c, p);
    nmod_poly_init(d, p);

    nmod_poly_rand(a, state, 2000);
    nmod_poly_rand(b, state, 2000);

    nmod_poly_mul(c, a, b);
    nmod_poly_mul_classical(d, a, b);
    flint_printf("classical: %d\n", nmod_poly_equal(c, d));
    nmod_poly_mul_KS(d, a, b);
    flint_printf("KS:        %d\n", nmod_poly_equal(c, d));
    nmod_poly_mul_KS2(d, a, b);
    flint_printf("KS2:       %d\n", nmod_poly_equal(c, d));
    nmod_poly_mul_KS4(d, a, b);
    flint_printf("KS4:       %d\n", nmod_poly_equal(c, d));

    /* a truncated product is the truncation of the product */
    nmod_poly_mullow(d, a, b, 1500);
    nmod_poly_truncate(c, 1500);
    flint_printf("mullow:    %d\n", nmod_poly_equal(c, d));

    /* squaring in place */
    nmod_poly_mul(c, a, a);
    nmod_poly_mul(a, a, a);
    flint_printf("aliasing:  %d\n", nmod_poly_equal(c, a));

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(d);
    flint_rand_clear(state);
    return 0;
}
```
```output
classical: 1
KS:        1
KS2:       1
KS4:       1
mullow:    1
aliasing:  1
```

## Moduli up to 23 {#small-moduli}

For $n \le 23$ the FFT product packs two coefficients into one and checks
with a Cauchy–Schwarz bound that the packed digits cannot overflow
([](algorithms/polynomials/multiplication#repack)). In FLINT 3.6
and earlier the last even-indexed coefficient of an operand of odd length entered that bound
unsquared, and `nmod_poly_mul` could return a wrong product for inputs
whose even coefficients are nearly all $n - 1$. The following program
prints `BUG` on affected versions and `ok` on fixed ones:

```c example name=nmod_poly_mul_repack_check
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t a, c, d;
    slong i, k = 270;

    nmod_poly_init(a, 23);
    nmod_poly_init(c, 23);
    nmod_poly_init(d, 23);
    for (i = 0; i <= k; i++)              /* a = 22 (1 + x^2 + ... + x^2k) */
        nmod_poly_set_coeff_ui(a, 2 * i, 22);
    nmod_poly_mul(c, a, a);
    nmod_poly_mul_classical(d, a, a);
    flint_printf("%s\n", nmod_poly_equal(c, d) ? "ok" : "BUG");
    nmod_poly_clear(a);
    nmod_poly_clear(c);
    nmod_poly_clear(d);
    return 0;
}
```
```output
ok
```
