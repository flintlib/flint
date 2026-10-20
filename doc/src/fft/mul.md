---
title: Multiplying large integers
functions: flint_mpn_mul_fft_main mul_truncate_sqrt2 mul_mfa_truncate_sqrt2
---

`flint_mpn_mul_fft_main` multiplies two natural numbers given as limb
arrays with the Schönhage–Strassen algorithm. It chooses a transform
length and a coefficient ring from tuning tables and then calls one of
the two lower-level functions documented below, which take these parameters
explicitly. For general use, call `flint_mpn_mul` (or `fmpz_mul`)
instead, which uses this function only in builds without `fft_small`,
from `FLINT_FFT_MUL_THRESHOLD` limbs ([](mpn_extras/mul)).

## How the product is computed {#method}

Each operand is cut into chunks of $b$ bits, which become the
coefficients of two polynomials whose product, evaluated at $2^b$, is the
product of the integers. The polynomials are multiplied by a truncated
sqrt2 convolution of length at most $4n$, $n = 2^{\mathit{depth}}$, with
coefficients modulo $2^{nw} + 1$, and the product coefficients are added
together with the right shifts. The chunk size

$$
b = \left\lfloor \frac{nw - \mathit{depth} - 1}{2} \right\rfloor
$$

guarantees that no product coefficient reaches $2^{nw}$, so the
convolution is exact; see
[Integer multiplication](algorithms/fft/schonhage-strassen#integers).
With $j_1$ and $j_2$ chunks, the convolution has $j_1 + j_2 - 1$
nonzero coefficients, which must not exceed $4n$. Only that many points
of the transform are computed (a [truncated FFT](algorithms/fft/truncated-fft)),
so the running time grows smoothly with the operand sizes.

- `mul_truncate_sqrt2` uses the recursive transforms `fft_truncate_sqrt2`
  and `ifft_truncate_sqrt2`, and multiplies the transformed coefficients
  with `flint_mpn_mulmod_2expp1_basecase`, a full product (by
  `flint_mpn_mul_n`) followed by a reduction.
- `mul_mfa_truncate_sqrt2` uses the
  [matrix Fourier algorithm](algorithms/fft/matrix-fourier) with
  $n_1 = 2^{\lfloor \mathit{depth}/2 \rfloor}$ columns, which keeps the
  working set of each partial transform in cache, runs the column and row
  transforms in threads, and multiplies the transformed coefficients with
  `fft_mulmod_2expp1`, which recurses into a negacyclic convolution for
  coefficients of more than 128 limbs. Since the pointwise products use
  `fft_mulmod_2expp1`, a limb count $nw/\bits$ above 128 must be a value
  returned by `fft_adjust_limbs`; otherwise the product is wrong. Powers
  of two and three times powers of two, the only values that
  `flint_mpn_mul_fft_main` produces, are such values.

Both functions allocate about $4n(nw/\bits + 1)$ limbs for each distinct
operand, that is, about twice the size of the product per operand when
the transform is nearly full, plus scratch space.

## Choice of parameters {#parameters}

`flint_mpn_mul_fft_main` searches for the smallest
$(\mathit{depth}, w)$ in the order $(6, 1), (6, 2), (7, 1), (7, 2), \ldots$
for which $j_1 + j_2 - 1 \le 4n$. Then:

- If $\mathit{depth} < 11$, the tuning table `FFT_TAB` gives an offset
  $k$; the transform is made $2^k$ times shorter with $4^k$ times larger
  $w$, so the coefficients become $2^k$ times larger, and $w$ is then
  decreased as far as the product still fits. The result goes to
  `mul_truncate_sqrt2`.
- Otherwise, if $j_1 + j_2 - 1 \le 3n$, the depth is decreased by one and
  $w$ multiplied by 3, which gives a transform half as long with
  coefficients 1.5 times larger and still room for the product. The
  result goes to `mul_mfa_truncate_sqrt2`.

For balanced operands of $\ell$ limbs each this gives, with the generic
tuning tables:

| $\ell$ | $(\mathit{depth}, w)$ | function | $4n$ | ring bits $nw$ | chunk bits $b$ | points used |
|---:|---|---|---:|---:|---:|---:|
| 1000 | (6, 16) | `mul_truncate_sqrt2` | 256 | 1024 | 508 | 251 |
| 10000 | (8, 10) | `mul_truncate_sqrt2` | 1024 | 2560 | 1275 | 1003 |
| 32000 | (9, 8) | `mul_truncate_sqrt2` | 2048 | 4096 | 2043 | 2005 |
| $10^5$ | (11, 2) | `mul_mfa_truncate_sqrt2` | 8192 | 4096 | 2042 | 6269 |
| $10^6$ | (13, 1) | `mul_mfa_truncate_sqrt2` | 32768 | 8192 | 4089 | 31303 |
| $10^7$ | (14, 3) | `mul_mfa_truncate_sqrt2` | 65536 | 49152 | 24568 | 52101 |

The full rules, the tables for the different processors and the
parameters of the pointwise products are given in
[Choice of parameters](algorithms/fft/schonhage-strassen#parameters).
The parameters depend on the total size, so an unbalanced product of
$10^5$ by 10 limbs is computed with the same transform as a balanced one
of about $5 \cdot 10^4$ limbs each; GMP's `mpn_mul` is much faster for
such shapes.

## Examples {#examples}

Multiplying two integers of about 600000 and 840000 bits, stored in
`mpz_t` variables. The function writes all $n_1 + n_2$ limbs of the
product; the top limb may be zero.

```c example name=fft_mul_main
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fft.h>

int main(void)
{
    mpz_t x, y, z;
    mp_size_t n1, n2, rn;
    mp_ptr r;

    mpz_init(x); mpz_init(y); mpz_init(z);
    mpz_ui_pow_ui(x, 3, 400000);
    mpz_ui_pow_ui(y, 7, 300000);
    n1 = mpz_size(x);
    n2 = mpz_size(y);

    /* the product of n1 and n2 limbs needs n1 + n2 limbs */
    r = flint_malloc((n1 + n2) * sizeof(mp_limb_t));
    flint_mpn_mul_fft_main(r, mpz_limbs_read(x), n1, mpz_limbs_read(y), n2);

    /* the top limb may be zero */
    rn = n1 + n2;
    while (r[rn - 1] == 0)
        rn--;

    mpz_mul(z, x, y);
    flint_printf("limbs: %wd x %wd -> %wd\n", n1, n2, rn);
    flint_printf("equal to mpz_mul: %d\n",
        mpz_size(z) == (size_t) rn && mpn_cmp(r, mpz_limbs_read(z), rn) == 0);

    flint_free(r);
    mpz_clear(x); mpz_clear(y); mpz_clear(z);
    return 0;
}
```
```output
limbs: 9907 x 13160 -> 23066
equal to mpz_mul: 1
```

Here `flint_mpn_mul_fft_main` uses `mul_truncate_sqrt2` with
$\mathit{depth} = 8$ and $w = 12$: a transform of length 1024 modulo
$2^{3072} + 1$ with 1531-bit chunks, of which 965 points are computed.

The same kind of product with explicit parameters: with
$\mathit{depth} = 9$ and $w = 1$, the chunks have
$b = \lfloor (512 - 10)/2 \rfloor = 251$ bits, and the program checks
that the product fits in the transform of length $4n = 2048$ before
calling both functions.

```c example name=fft_mul_params
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fft.h>

int main(void)
{
    mpz_t x, y, z;
    mp_size_t n1, n2, j1, j2;
    flint_bitcnt_t depth = 9, w = 1, b;
    mp_size_t n = WORD(1) << depth;
    mp_ptr r, s;

    mpz_init(x); mpz_init(y); mpz_init(z);
    mpz_ui_pow_ui(x, 3, 100000);
    mpz_ui_pow_ui(y, 5, 70000);
    n1 = mpz_size(x);
    n2 = mpz_size(y);
    mpz_mul(z, x, y);

    /* chunk size and number of chunks used by the convolution */
    b = (n * w - depth - 1) / 2;
    j1 = (n1 * FLINT_BITS - 1) / b + 1;
    j2 = (n2 * FLINT_BITS - 1) / b + 1;
    flint_printf("n1 = %wd, n2 = %wd, b = %wu, j1 + j2 - 1 = %wd <= 4n = %wd\n",
                 n1, n2, b, j1 + j2 - 1, 4 * n);

    r = flint_malloc((n1 + n2) * sizeof(mp_limb_t));
    s = flint_malloc((n1 + n2) * sizeof(mp_limb_t));
    mul_truncate_sqrt2(r, mpz_limbs_read(x), n1, mpz_limbs_read(y), n2, depth, w);
    mul_mfa_truncate_sqrt2(s, mpz_limbs_read(x), n1, mpz_limbs_read(y), n2, depth, w);

    flint_printf("truncate: %d, mfa: %d\n",
        mpn_cmp(r, mpz_limbs_read(z), mpz_size(z)) == 0,
        mpn_cmp(s, mpz_limbs_read(z), mpz_size(z)) == 0);

    flint_free(r); flint_free(s);
    mpz_clear(x); mpz_clear(y); mpz_clear(z);
    return 0;
}
```
```output
n1 = 2477, n2 = 2540, b = 251, j1 + j2 - 1 = 1279 <= 4n = 2048
truncate: 1, mfa: 1
```

For these sizes `flint_mpn_mul_fft_main` itself would choose
$\mathit{depth} = 7$ and $w = 10$ (636-bit chunks, 505 points of a
transform of length 512), which takes about half the time of the
parameters above.
