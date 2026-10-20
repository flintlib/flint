---
title: Multiplying integers
functions: fmpz_mul fmpz_addmul fmpz_fmma
---

`fmpz_mul` is the function to use for products of integers, with
`fmpz_mul_ui` and `fmpz_mul_si` for a word-size factor. `fmpz_addmul` and
`fmpz_submul` add or subtract a product to an accumulator, which saves a
temporary and, for large values, a pass over memory; `fmpz_fmma` and
`fmpz_fmms` form $ab \pm cd$ in one call. All of them choose an algorithm
by the sizes of the operands, from a two-word product for inline values
to FFT multiplication for operands of thousands of limbs.

::: summary algorithms/integers/multiplication :::

## Fast paths {#fast-paths}

`fmpz_mul` first looks at the words of the operands:

- Two inline values have at most $\bits - 2$ bits each, so their product
  fits in two words. It is formed with one `smul_ppmm` and stored with
  `fmpz_set_signed_uiui`, inline if it fits and in an `mpz` otherwise.
  No function call into GMP is made.
- An inline value times a large one, or a large value times an `ulong`
  or `slong` (`fmpz_mul_ui`, `fmpz_mul_si`), is one pass of
  `mpn_mul_1` over the limbs.
- A zero factor clears the output.

Two large values are multiplied limb by limb. Operands of one or two
limbs each are multiplied inline (a one-limb by one-limb product has two
limbs and cannot be inline, since each factor exceeds `COEFF_MAX`). Larger
operands go to `flint_mpn_mul`, or to `flint_mpn_sqr` if both arguments
are the same object: `fmpz_mul(r, a, a)` is a squaring, which needs
roughly half of the limb products of a general product in the basecase
and one forward transform instead of two with FFTs. Two different
objects with equal values are multiplied as a general product. If the
output is the same object as an input, the input is copied first, since
the product cannot be formed in place.

## Algorithm selection {#selection}

`flint_mpn_mul` multiplies an $m$-limb by an $n$-limb operand, $m \ge n$.
On x86-64 processors with the ADX and BMI2 extensions (Intel since
Broadwell, AMD since Zen) it uses, with the thresholds of the installed
library:

| Shape | Algorithm |
|---|---|
| $m \le 16$ | FLINT's assembly basecase for each pair $(m, n)$ |
| $n = 1$ | `mpn_mul_1` |
| $n \le 230$ and $5n \ge 4m$ | Karatsuba (`flint_mpn_mul_toom22`) on FLINT's basecase |
| $n < 500$, unless $m + n > 1000$ and $3n > 500$ | GMP's `mpn_mul`: its basecase, Karatsuba and Toom–Cook 3, 4, 6.5 and 8.5, with unbalanced variants |
| otherwise | small-prime FFT of `fft_small` (`mpn_mul_default_mpn_ctx`) |

The threshold 500 is `FLINT_FFT_SMALL_MUL_THRESHOLD`, which depends on
the target processor ([](fft_small/usage#callers)). Squarings use
assembly code up to 14 limbs, then GMP's `mpn_sqr`, then `fft_small` from
`FLINT_FFT_SMALL_SQR_THRESHOLD` limbs. ARM64 also has FLINT assembly
basecases and Karatsuba code, with other thresholds
([](mpn_extras/mul)). On other processors FLINT's portable C code
handles operands of up to 7 limbs and the rest of the basecase and
Karatsuba ranges belong to GMP. Builds without `fft_small` call GMP's
`mpn_mul` (including GMP's own FFT) up to 32000 limbs and the
Schönhage–Strassen code of the [`fft`](fft) module above.
[](algorithms/integers/multiplication) describes each algorithm and its
cost.

## Accumulating products {#addmul}

`fmpz_addmul(f, g, h)` computes $f + gh$. For inline operands the product
and the sum are formed in two or three words; for a word-size factor it
is GMP's `mpz_addmul_ui`, a single fused pass; for two large factors the
product is formed in a temporary limb array and added to the limbs of $f$
in place, without allocating an `fmpz`. Use it rather than `fmpz_mul`
followed by `fmpz_add` in dot products and matrix products;
`_fmpz_vec_dot` does this for whole vectors and adds further
optimisations for small entries ([](fmpz_vec/dot)).

`fmpz_fmma` and `fmpz_fmms` compute $ab + cd$ and $ab - cd$, entirely in
words if all four operands are inline, and otherwise as a product
followed by `fmpz_addmul` or `fmpz_submul`, taking care of aliasing
between the output and the inputs.

## Examples {#examples}

Products and accumulations with Fibonacci numbers, which illustrate the
inline fast path for two inline values and the identity of Cassini,
$F_{n+1} F_{n-1} - F_n^2 = (-1)^n$:

```c example name=fmpz_mul_basic
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    fmpz_t a, b, c, d, r;

    fmpz_init(a); fmpz_init(b); fmpz_init(c); fmpz_init(d); fmpz_init(r);

    /* Cassini: F_101 F_99 - F_100^2 = (-1)^100 */
    fmpz_fib_ui(a, 101);
    fmpz_fib_ui(b, 100);
    fmpz_fib_ui(d, 99);
    fmpz_fmms(r, a, d, b, b);
    flint_printf("F_101 F_99 - F_100^2 = "); fmpz_print(r); flint_printf("\n");

    /* accumulate products without temporaries */
    fmpz_zero(r);
    fmpz_addmul(r, a, a);
    fmpz_addmul(r, b, b);
    fmpz_submul_ui(r, d, 3);
    flint_printf("F_101^2 + F_100^2 - 3 F_99 = "); fmpz_print(r); flint_printf("\n");

    /* small operands: the product is formed in two words */
    fmpz_set_si(c, -COEFF_MAX);
    fmpz_mul(r, c, c);
    flint_printf("COEFF_MAX^2 = "); fmpz_print(r); flint_printf("\n");

    fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(d); fmpz_clear(r);
    return 0;
}
```
```output
F_101 F_99 - F_100^2 = 1
F_101^2 + F_100^2 - 3 F_99 = 453973694165307953196640200709906953726748
COEFF_MAX^2 = 21267647932558653957237540927630737409
```

The program compares `fmpz_mul` with `mpz_mul` for random operands of $n$
and $2n$ limbs, from the basecase up to `fft_small`, and checks that the
squaring path agrees with `fmpz_pow_ui`:

```c example name=fmpz_mul_sizes
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    flint_rand_t state;
    fmpz_t a, b, c, d;
    mpz_t x, y, z;
    slong sizes[] = { 1, 2, 10, 100, 300, 1000, 5000 };
    int i;

    flint_rand_init(state);
    fmpz_init(a); fmpz_init(b); fmpz_init(c); fmpz_init(d);
    mpz_init(x); mpz_init(y); mpz_init(z);

    for (i = 0; i < 7; i++)
    {
        /* operands of sizes[i] and 2 sizes[i] limbs */
        fmpz_randbits(a, state, sizes[i] * FLINT_BITS);
        fmpz_randbits(b, state, 2 * sizes[i] * FLINT_BITS);

        fmpz_mul(c, a, b);              /* FLINT */
        fmpz_get_mpz(x, a);
        fmpz_get_mpz(y, b);
        mpz_mul(z, x, y);               /* GMP */
        fmpz_set_mpz(d, z);

        flint_printf("%5wd x %5wd limbs: %s", sizes[i], 2 * sizes[i],
            fmpz_equal(c, d) ? "same" : "DIFFERENT");

        fmpz_mul(c, a, a);              /* detected as a squaring */
        fmpz_pow_ui(d, a, 2);
        flint_printf(", square %s\n", fmpz_equal(c, d) ? "same" : "DIFFERENT");
    }

    fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(d);
    mpz_clear(x); mpz_clear(y); mpz_clear(z);
    flint_rand_clear(state);
    return 0;
}
```
```output
    1 x     2 limbs: same, square same
    2 x     4 limbs: same, square same
   10 x    20 limbs: same, square same
  100 x   200 limbs: same, square same
  300 x   600 limbs: same, square same
 1000 x  2000 limbs: same, square same
 5000 x 10000 limbs: same, square same
```
