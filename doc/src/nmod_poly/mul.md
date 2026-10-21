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

`_nmod_poly_mul` uses the classical algorithm if the shorter operand has
length at most 5. Otherwise it asks `_nmod_poly_mullow_want_fft_small`
whether `fft_small` is faster, and if not it uses, with $c =
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

`nmod_poly_mullow` uses the classical algorithm for truncation lengths
below $10 + b^2/10$ (from 16 for 8-bit to 419 for 64-bit moduli), KS
otherwise, and `fft_small` from separately tuned cutoffs.
`nmod_poly_mulmid` first drops the operand coefficients that do not
contribute to the requested range and then uses the classical algorithm,
`fft_small` (which computes a middle product at the cost of a product of
the output length, [](algorithms/polynomials/middle-product)), or KS for
small moduli. `nmod_poly_mulhigh` calls `nmod_poly_mulmid` for all but
tiny inputs and sets the low coefficients to zero.

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
